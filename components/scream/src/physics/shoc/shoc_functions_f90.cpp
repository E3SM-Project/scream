#include "shoc_functions_f90.hpp"

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/scream_pack_kokkos.hpp"
#include "shoc_f90.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C interface to SHOC fortran calls. The stubs below will link to fortran definitions in shoc_iso_c.f90
//
extern "C" {

void shoc_init_c(int nlev, Real gravit, Real rair, Real rh2o, Real cpair,
                 Real zvir, Real latvap, Real latice, Real karman,
                 Real* pref_mid, int nbot_shoc, int ntop_shoc);

void shoc_grid_c(int shcol, int nlev, int nlevi, Real *zt_grid, Real *zi_grid,
                 Real *pdel, Real *dz_zt, Real *dzi_zi, Real *rho_zt);

void calc_shoc_varorcovar_c(Int shcol, Int nlev, Int nlevi,  Real tunefac,
                            Real *isotropy_zi, Real *tkh_zi, Real *dz_zi,
			    Real *invar1, Real *invar2, Real *varorcovar);

void integ_column_stability_c(Int nlev, Int shcol, Real *dz_zt, Real *pres,
			      Real *brunt, Real *brunt_int);

void compute_shr_prod_c(Int nlevi, Int nlev, Int shcol, Real *dz_zi,
                        Real *u_wind, Real *v_wind, Real *sterm);

void isotropic_ts_c(Int nlev, Int shcol, Real *brunt_int, Real *tke,
                    Real *a_diss, Real *brunt, Real *isotropy);

void adv_sgs_tke_c(Int nlev, Int shcol, Real dtime, Real *shoc_mix,
                   Real *wthv_sec, Real *sterm_zt, Real *tk,
                   Real *tke, Real *a_diss);

void eddy_diffusivities_c(Int nlev, Int shcol, Real *obklen, Real *pblh,
                          Real *zt_grid, Real *shoc_mix, Real *sterm_zt,
                          Real *isotropy, Real *tke, Real *tkh, Real *tk);

void calc_shoc_vertflux_c(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux);

void shoc_diag_second_moments_srf_c(Int shcol, Real* wthl, Real* uw, Real* vw, Real* ustar2, Real* wstar);

void shoc_diag_second_moments_ubycond_c(Int shcol, Int num_tracer, Real* thl, Real* qw, Real* wthl, Real* wqw, Real* qwthl, Real* uw, Real* vw,
      Real* wtke, Real* wtracer);

}

namespace scream {
namespace shoc {

namespace {

template <size_t N, size_t M>
void gen_random_data(const std::array<std::pair<Real, Real>, N>& ranges,
                     const std::array<Real**, M>& ptrs,
                     Real* data, Int nk)
{
  static_assert(N <= M, "Require at least as many ptrs as ranges");

  Int offset = 0;
  std::default_random_engine generator;

  for (size_t i = 0; i < N; ++i) {
    std::uniform_real_distribution<Real> data_dist(ranges[i].first, ranges[i].second);
    *ptrs[i] = data + offset;
    offset += nk;
    for(Int k = 0; k < nk; ++k) {
      (*ptrs[i])[k] = data_dist(generator);
    }
  }

  for (size_t i = N; i < M; ++i) {
    *ptrs[i] = data + offset;
    offset += nk;
  }
}

}

//
// Data struct
//

SHOCDataBase::SHOCDataBase(Int shcol_, Int nlev_, Int nlevi_,
                           const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i) :
  shcol(shcol_),
  nlev(nlev_),
  nlevi(nlevi_),
  m_total(shcol_ * nlev_),
  m_totali(shcol_ * nlevi_),
  m_ptrs(ptrs),
  m_ptrs_i(ptrs_i),
  m_data(m_ptrs.size() * m_total + m_ptrs_i.size() * m_totali, 0)
{
  init_ptrs();
}

SHOCDataBase::SHOCDataBase(const SHOCDataBase &rhs, const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i) :
  shcol(rhs.shcol),
  nlev(rhs.nlev),
  nlevi(rhs.nlevi),
  m_total(rhs.m_total),
  m_totali(rhs.m_totali),
  m_ptrs(ptrs),
  m_ptrs_i(ptrs_i),
  m_data(rhs.m_data)
{
  init_ptrs();
}

SHOCDataBase& SHOCDataBase::operator=(const SHOCDataBase& rhs)
{
  shcol    = rhs.shcol;
  nlev     = rhs.nlev;
  nlevi    = rhs.nlevi;
  m_total  = rhs.m_total;
  m_totali = rhs.m_totali;
  m_data   = rhs.m_data; // Copy

  init_ptrs();

  return *this;
}

void SHOCDataBase::init_ptrs()
{
  Int offset         = 0;
  Real *data_begin   = m_data.data();

  for (size_t i = 0; i < m_ptrs.size(); ++i) {
    *(m_ptrs[i]) = data_begin + offset;
    offset += m_total;
  }

  for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
    *(m_ptrs_i[i]) = data_begin + offset;
    offset += m_totali;
  }
}

void SHOCDataBase::randomize()
{
  std::default_random_engine generator;
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  for (size_t i = 0; i < m_ptrs.size(); ++i) {
    for (size_t j = 0; j < m_total; ++j) {
      (*(m_ptrs[i]))[j] = data_dist(generator);
    }
  }

  for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
    for (size_t j = 0; j < m_totali; ++j) {
      (*(m_ptrs_i[i]))[j] = data_dist(generator);
    }
  }
}

MomSrfData::MomSrfData(
  Int shcol_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  shcol(shcol_),
  m_nk(shcol_),
  m_data(NUM_ARRAYS*m_nk, 1.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs = {&wthl, &uw, &vw, &ustar2, &wstar};
  gen_random_data(ranges, ptrs, m_data.data(), m_nk);
}

MomSrfData::MomSrfData(const MomSrfData& rhs) :
  shcol(rhs.shcol),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  Real** ptrs[NUM_ARRAYS] = {&wthl, &uw, &vw, &ustar2, &wstar};

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
}

MomUbycondData::MomUbycondData(
  Int shcol_,
  Int num_tracer_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  shcol(shcol_),
  num_tracer(num_tracer_),
  m_nk(shcol_),
  m_data((NUM_ARRAYS+num_tracer)*m_nk, 1.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs = {&thl, &qw, &qwthl, &wthl, &wqw, &uw, &vw, &wtke};
  gen_random_data(ranges, ptrs, m_data.data(), m_nk);

  Real* data_begin = m_data.data();
  Int offset = m_nk*NUM_ARRAYS;
  const std::array< std::pair<Real, Real>, 1 > range2 = { std::make_pair(0, 1) };
  std::array<Real**, 1> ptrs_2d = {&wtracer};
  gen_random_data(range2, ptrs_2d, data_begin+offset, m_nk*num_tracer);
}

MomUbycondData::MomUbycondData(const MomUbycondData& rhs) :
  shcol(rhs.shcol),
  num_tracer(rhs.num_tracer),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  Real** ptrs[NUM_ARRAYS+1] = {&thl, &qw, &qwthl, &wthl, &wqw, &uw, &vw, &wtke, &wtracer};

  for (size_t i = 0; i < NUM_ARRAYS+1; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
}

void MomUbycondData::transpose()
{
  std::vector<Real> data(m_data.size());
  // Transpose on the wtracer
  util::transpose<util::TransposeDirection::f2c>(wtracer, data.data()+NUM_ARRAYS*shcol, shcol, num_tracer);
  m_data = data;
}

//
// Glue functions to call fortran from from C++ with the Data struct
//
// In all C++ -> Fortran bridge functions you should see shoc_init(nlev, true).
// We are provisionally following P3 here in case SHOC uses global data. The
// 'true' argument is to set shoc to use its fortran implementations instead of
// calling back to C++. We want this behavior since it doesn't make much sense
// for C++ to bridge over to fortran only to have fortran bridge back to C++.
// Anyone who wants the C++ implementation should call it directly. We need
// need to be aware of data layout since f90 is different from cxx. All these
// functions will expect incoming data to be C layout. They will transpose to f90
// before calling fortran and then back to C before returning.
//

void calc_shoc_varorcovar(SHOCVarorcovarData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  calc_shoc_varorcovar_c(d.shcol, d.nlev, d.nlevi, d.tunefac, d.isotropy_zi, d.tkh_zi,
                         d.dz_zi, d.invar1, d.invar2, d.varorcovar);
  d.transpose<util::TransposeDirection::f2c>();
}

void shoc_grid(SHOCGridData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  shoc_grid_c(d.shcol, d.nlev, d.nlevi, d.zt_grid, d.zi_grid, d.pdel, d.dz_zt,
              d.dz_zi, d.rho_zt);
  d.transpose<util::TransposeDirection::f2c>();
}

void calc_shoc_vertflux(SHOCVertfluxData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  calc_shoc_vertflux_c(d.shcol, d.nlev, d.nlevi, d.tkh_zi, d.dz_zi, d.invar,
		       d.vertflux);
  d.transpose<util::TransposeDirection::f2c>();
}

void integ_column_stability(SHOCColstabData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  integ_column_stability_c(d.nlev, d.shcol, d.dz_zt, d.pres, d.brunt, d.brunt_int);
  d.transpose<util::TransposeDirection::f2c>();
}

void compute_shr_prod(SHOCTkeshearData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  compute_shr_prod_c(d.nlevi, d.nlev, d.shcol, d.dz_zi, d.u_wind,
                       d.v_wind, d.sterm);
  d.transpose<util::TransposeDirection::f2c>();
}

void isotropic_ts(SHOCIsotropicData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  isotropic_ts_c(d.nlev, d.shcol, d.brunt_int, d.tke, d.a_diss,
                 d.brunt, d.isotropy);
  d.transpose<util::TransposeDirection::f2c>();
}

void adv_sgs_tke(SHOCAdvsgstkeData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  adv_sgs_tke_c(d.nlev, d.shcol, d.dtime, d.shoc_mix, d.wthv_sec,
                d.sterm_zt, d.tk, d.tke, d.a_diss);
  d.transpose<util::TransposeDirection::f2c>();
}

void eddy_diffusivities(SHOCEddydiffData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  eddy_diffusivities_c(d.nlev, d.shcol, d.obklen, d.pblh, d.zt_grid,
     d.shoc_mix, d.sterm_zt, d.isotropy, d.tke, d.tkh, d.tk);
  d.transpose<util::TransposeDirection::f2c>();
}

void shoc_diag_second_moments_srf(MomSrfData& d)
{
  Int nlev = 128;
  shoc_init(nlev, true);
  shoc_diag_second_moments_srf_c(d.shcol, d.wthl, d.uw, d.vw, d.ustar2, d.wstar);
}

void shoc_diag_second_moments_ubycond(MomUbycondData& d)
{
  Int nlev = 128;
  shoc_init(nlev, true);
  shoc_diag_second_moments_ubycond_c(d.shcol, d.num_tracer, d.thl, d.qw, d.wthl, d.wqw, d.qwthl, d.uw, d.vw, d.wtke, d.wtracer);
  d.transpose();
}

//
// _f function definitions. These expect data in C layout
//

void calc_shoc_vertflux_f(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 4;

  Kokkos::Array<view_2d, num_arrays> temp_d;
  Kokkos::Array<size_t, num_arrays> dim1_sizes     = {shcol,  shcol, shcol, shcol};
  Kokkos::Array<size_t, num_arrays> dim2_sizes     = {nlevi,  nlevi, nlev,  nlevi};
  Kokkos::Array<const Real*, num_arrays> ptr_array = {tkh_zi, dz_zi, invar, vertflux};

  // Sync to device
  pack::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

  view_2d
    tkh_zi_d  (temp_d[0]),
    dz_zi_d   (temp_d[1]),
    invar_d   (temp_d[2]),
    vertflux_d(temp_d[3]);

  const Int nk_pack = scream::pack::npack<Spack>(nlev);
  const auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto otkh_zi_d   = util::subview(tkh_zi_d, i);
    const auto odz_zi_d    = util::subview(dz_zi_d, i);
    const auto oinvar_d    = util::subview(invar_d, i);
    const auto overtflux_d = util::subview(vertflux_d, i);

    SHF::calc_shoc_vertflux(team, nlev, otkh_zi_d, odz_zi_d, oinvar_d, overtflux_d);
  });

  // Sync back to host
  Kokkos::Array<view_2d, 1> inout_views = {vertflux_d};
  pack::device_to_host({vertflux}, {shcol}, {nlevi}, inout_views, true);
}

void shoc_diag_second_moments_srf_f(Int shcol, Real* wthl, Real* uw, Real* vw, Real* ustar2, Real* wstar)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using view_1d    = typename SHOC::view_1d<Spack>;
  using KT         = typename SHOC::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;

  const Int nshcol_pack = scream::pack::npack<Spack>(shcol);
  Kokkos::Array<view_1d, MomSrfData::NUM_ARRAYS> mom_d;

  pack::host_to_device({wthl, uw, vw, ustar2, wstar}, shcol, mom_d);

  view_1d wthl_d(mom_d[0]),
          uw_d(mom_d[1]),
          vw_d(mom_d[2]),
          ustar2_d(mom_d[3]),
          wstar_d(mom_d[4]);

  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nshcol_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    SHOC::shoc_diag_second_moments_srf(team, shcol, wthl_d, uw_d, vw_d, ustar2_d, wstar_d);
  });

  // back to host
  Kokkos::Array<view_1d, 2> host_views = {ustar2_d, wstar_d};

  pack::device_to_host({ustar2, wstar}, shcol, host_views);

}

void shoc_diag_second_moments_ubycond_f(Int shcol, Int num_tracer, Real* thl, Real* qw, Real* wthl, Real* wqw, Real* qwthl, Real* uw, Real* vw,
      Real* wtke, Real* wtracer)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Spack      = typename SHOC::Spack;
  using view_1d    = typename SHOC::view_1d<Spack>;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using KT         = typename SHOC::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;

  Kokkos::Array<view_1d, MomUbycondData::NUM_ARRAYS> uby_1d;
  Kokkos::Array<view_2d, 1> uby_2d;

  pack::host_to_device({thl, qw, qwthl, wthl, wqw, uw, vw, wtke}, shcol, uby_1d);
  pack::host_to_device({wtracer}, shcol, num_tracer, uby_2d, true);

  view_1d thl_d(uby_1d[0]),
          qw_d(uby_1d[1]),
          qwthl_d(uby_1d[2]),
          wthl_d(uby_1d[3]),
          wqw_d(uby_1d[4]),
          uw_d(uby_1d[5]),
          vw_d(uby_1d[6]),
          wtke_d(uby_1d[7]);

  view_2d wtracer_d(uby_2d[0]);

  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, shcol);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    SHOC::shoc_diag_second_moments_ubycond(team, shcol, num_tracer, thl_d, qw_d, wthl_d, wqw_d, qwthl_d, uw_d, vw_d, wtke_d, wtracer_d);

  });

  // back to host
  Kokkos::Array<view_1d, 8> host_views = {thl_d, qw_d, qwthl_d, wthl_d, wqw_d, uw_d, vw_d, wtke_d};
  Kokkos::Array<view_2d, 1> host_2d_views = {wtracer_d};

  pack::device_to_host({thl, qw, qwthl, wthl, wqw, uw, vw, wtke}, shcol, host_views);
  pack::device_to_host({wtracer}, shcol, num_tracer, host_2d_views, true);
}

} // namespace shoc
} // namespace scream
