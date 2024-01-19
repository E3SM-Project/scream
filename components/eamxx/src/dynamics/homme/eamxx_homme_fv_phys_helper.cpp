#include "eamxx_homme_fv_phys_helper.hpp"

#include "share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/field/field.hpp"
#include "dynamics/homme/homme_dimensions.hpp"

// HOMMEXX Includes
#include "Context.hpp"
#include "FunctorsBuffersManager.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"
#include "GllFvRemap.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_utils.hpp"

// Defined in homme's gllfvremap_mod.F90
extern "C" void gfr_init_hxx();

namespace scream {

void HommeFvPhysHelper::
set_grids (const std::shared_ptr<const AbstractGrid>& dyn_grid,
           const std::shared_ptr<const AbstractGrid>& cgll_grid,
           const std::shared_ptr<const AbstractGrid>& phys_grid)
{
  m_dyn_grid = dyn_grid;
  m_cgll_grid = cgll_grid;
  m_phys_grid = phys_grid;

  const auto phys_grid_name = phys_grid->name();

  // Check the name of the physics grid, to find out if FvPhys is used
  if (phys_grid_name.size() >= 11 and
      phys_grid_name.substr(0, 10) == "Physics PG")
  {
    const auto pgN_str = phys_grid_name.substr(10, std::string::npos);
    std::istringstream ss(pgN_str);
    try {
      ss >> pgN;
    } catch (...) {
      pgN = -1;
    }
    fv_phys_active = true;
  }
}

void HommeFvPhysHelper::
requested_buffer_size_in_bytes () {
  if (not fv_phys_active) return;

  // Init FV<->GLL remap, since we'll need its
  using namespace Homme;
  auto& c = Context::singleton();
  EKAT_REQUIRE_MSG (c.has<SimulationParams>(),
      "Error! SimulationParams was not yet created in the Homme context.\n");

  auto& gfr = c.create_if_not_there<GllFvRemap>();
  gfr.reset(c.get<SimulationParams>());
  gfr_init_hxx();

  auto& fbm = c.create_if_not_there<FunctorsBuffersManager>();
  fbm.request_size(gfr.requested_buffer_size());
}

void HommeFvPhysHelper::
dyn_to_fv_phys_init (const bool restart, const util::TimeStamp& t0) {
  if (not fv_phys_active) return;

  constexpr int N = HOMMEXX_PACK_SIZE;
  using Pack = ekat::Pack<Real,N>;
  const auto nlevs = m_phys_grid->get_num_vertical_levels();
  const auto npacks = ekat::PackInfo<N>::num_packs(nlevs);
  const auto ncols = m_phys_grid->get_num_local_dofs();
  const auto qsize = m_Q_phys.get_view<const Real***>().extent(1);
  const auto FT = m_FT_phys.get_view<Pack**>();
  const auto FM = m_FM_phys.get_view<Pack***>();

  view_ND<Real,2> T;
  view_ND<Real,3> uv;
  view_ND<Real,3> q;
  if (restart) {
    constexpr int NGP = HOMMEXX_NP;
    const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
    const auto npg = pgN*pgN;

    T = view_ND<Real,2>("T_mid_tmp",nelem*npg,npacks*N);
    uv = view_ND<Real,3>("horiz_winds_tmp",nelem*npg,2,npacks*N);
    // Really need just the first tracer.
    q = view_ND<Real,3>("tracers_tmp",nelem*npg,qsize,npacks*N);

    remap_dyn_to_fv_phys(T,uv,q);
    assert(ncols == nelem*npg);

    view_ND<PackT,2> T_packed (reinterpret_cast<Pack*>(T.data()),  ncols,    npacks);
    view_ND<PackT,3> uv_packed(reinterpret_cast<Pack*>(uv.data()), ncols, 2, npacks);
    copy_prev(ncols,nlevs,T_packed,uv_packed,FT,FM);
  } else {
    T  = m_T_phys.get_view<Real**>();
    uv = m_uv_phys.get_view<Real***>();
    q  = m_Q_phys.get_view<Real***>();
    remap_dyn_to_fv_phys(T,uv,q);

    view_ND<PackT,2> T_packed (reinterpret_cast<Pack*>(T.data()),  ncols,    npacks);
    view_ND<PackT,3> uv_packed(reinterpret_cast<Pack*>(uv.data()), ncols, 2, npacks);
    copy_prev(ncols, npacks, T_packed, uv_packed, FT, FM);

    // In an initial run, the AD only reads IC for the Physics GLL fields,
    // and this class has just taken care of remapping them to the FV grid.
    // Therefore, the timestamp of the FV fields has *not* been set yet,
    // which can cause serious issues downstream. For details, see
    //   https://github.com/E3SM-Project/scream/issues/2250
    // To avoid any problem, we simply set the timestamp of FV state fields here.
    // NOTE: even if the init sequence *ever* changed, and the AD *did* read
    // IC for FV fields, this step remains safe (we're setting the same t0)
    for (auto f : {m_T_phys,m_uv_phys,m_ps_phys,m_phis_phys,m_omega_phys,m_dp_phys,m_Q_phys}) {
      f.get_header().get_tracking().update_time_stamp(t0);
    }
  }
}

void HommeFvPhysHelper::
remap_dyn_to_fv_phys () const
{
  auto T  = m_T_phys.get_view<Real**>();
  auto uv = m_uv_phys.get_view<Real***>();
  auto q  = m_Q_phys.get_view<Real***>();
  remap_dyn_to_fv_phys(T,uv,q);
}

void HommeFvPhysHelper::
remap_dyn_to_fv_phys (const view_ND<Real,2>& T_out, const view_ND<Real,3>& uv_out, const view_ND<Real,3>& Q_out) const
{
  if (not fv_phys_active) return;

  using GFR = Homme::GllFvRemap;

  const auto& c = Homme::Context::singleton();
  auto& gfr = c.get<Homme::GllFvRemap>();
  const auto time_idx = c.get<Homme::TimeLevel>().n0;
  constexpr int NGP = HOMMEXX_NP;
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const auto npg = pgN*pgN;
  const auto nlev = T_out.extent_int(1);
  const auto nq = Q_out.extent_int(1);
  assert(T_out.extent_int(0) == nelem*npg);
  assert(uv_out.extent_int(1) == 2);

  const auto ps    = GFR::Phys1T(m_ps_phys.get_view<Real*>().data(), nelem, npg);
  const auto phis  = GFR::Phys1T(m_phis_phys.get_view<Real*>().data(), nelem, npg);
  const auto T     = GFR::Phys2T(T_out.data(), nelem, npg, nlev);
  const auto omega = GFR::Phys2T(m_omega_phys.get_view<Real**>().data(),nelem, npg, nlev);
  const auto uv    = GFR::Phys3T(uv_out.data(),nelem, npg, 2, nlev);
  const auto q     = GFR::Phys3T(Q_out.data(),nelem, npg, nq, nlev);
  const auto dp    = GFR::Phys2T(m_dp_phys.get_view<Real**>().data(),nelem, npg, nlev);

  gfr.run_dyn_to_fv_phys(time_idx, ps, phis, T, omega, uv, q, &dp);
  Kokkos::fence();
}

void HommeFvPhysHelper::
remap_fv_phys_to_dyn () const {
  if (not fv_phys_active) return;

  using GFR = Homme::GllFvRemap;

  const auto& c = Homme::Context::singleton();
  auto& gfr = c.get<Homme::GllFvRemap>();
  const auto time_idx = c.get<Homme::TimeLevel>().n0;
  constexpr int NGP = HOMMEXX_NP;
  const int nelem = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const auto npg = pgN*pgN;

  const auto FT = m_FT_phys.get_view<const Real**>();
  const auto FM = m_FM_phys.get_view<const Real***>();
  const auto tracers = m_Q_phys.get_view<const Real***>();

  const auto nlev    = FT.extent_int(1);
  const auto nq      = tracers.extent_int(1);
  const auto uv_ndim = FM.extent_int(1);

  assert(FT.extent_int(0) == nelem*npg);
  assert(uv_ndim == 2);

  const auto T  = GFR::CPhys2T(FT.data(), nelem, npg, nlev);
  const auto uv = GFR::CPhys3T(FM.data(), nelem, npg, uv_ndim, nlev);
  const auto q  = GFR::CPhys3T(tracers.data(), nelem, npg, nq, nlev);

  gfr.run_fv_phys_to_dyn(time_idx, T, uv, q);
  Kokkos::fence();
  gfr.run_fv_phys_to_dyn_dss();
  Kokkos::fence();
}

void HommeFvPhysHelper::
copy_prev (const int ncols, const int nlevs,
           const view_ND<PackT,2>& T,  const view_ND<PackT,3>& uv,
           const view_ND<PackT,2>& FT, const view_ND<PackT,3>& FM)
{
  using KT = KokkosTypes<DefaultDevice>;
  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;

  const auto npacks = ekat::PackInfo<N>::num_packs(nlevs);
  const auto policy = ESU::get_default_team_policy(ncols, npacks);

  // Copy physics T,uv state to FT,M to form tendencies in next dynamics step.
  auto f = KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, npacks),
                         [&] (const int ilev) {
      FT(icol,ilev)   = T(icol,ilev);
      FM(icol,0,ilev) = uv(icol,0,ilev);
      FM(icol,1,ilev) = uv(icol,1,ilev);
    });
  };

  Kokkos::parallel_for(policy, f);

  Kokkos::fence();
}

void HommeFvPhysHelper::
clean_up ()
{
  m_phys_grid = m_cgll_grid = m_dyn_grid = nullptr;

  Field empty;

  m_FT_phys    = empty;
  m_FM_phys    = empty;
  m_T_phys     = empty;
  m_uv_phys    = empty;
  m_ps_phys    = empty;
  m_phis_phys  = empty;
  m_omega_phys = empty;
  m_dp_phys    = empty;
  m_Q_phys     = empty;
}

} // namespace scream
