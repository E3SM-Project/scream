#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/rrtmgp_heating_rate.hpp"
#include "physics/rrtmgp/share/shr_orb_mod_c2f.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/util/scream_column_ops.hpp"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "YAKL/YAKL.h"
#include "ekat/ekat_assert.hpp"

namespace scream {

  using KT = KokkosTypes<DefaultDevice>;
  using ExeSpace = KT::ExeSpace;
  using MemberType = KT::MemberType;

RRTMGPRadiation::RRTMGPRadiation (const ekat::Comm& comm, const ekat::ParameterList& params) 
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}  // RRTMGPRadiation::RRTMGPRadiation

// =========================================================================================
void RRTMGPRadiation::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {

  using namespace ekat::units;

  // Gather the active gases from the rrtmgp parameter list and assign to the m_gas_names vector.
  auto active_gases = m_params.get<std::vector<std::string>>("active_gases");
  for (auto& it : active_gases) {
    // Make sure only unique names are added
    if (std::find(m_gas_names.begin(), m_gas_names.end(), it) == m_gas_names.end()) {
      m_gas_names.push_back(it);
    }
  }
  m_ngas = m_gas_names.size();

  // Declare the set of fields used by rrtmgp
  auto kgkg = kg/kg;
  kgkg.set_string("kg/kg");
  auto m3 = m * m * m;
  auto Wm2 = W / m / m;
  Wm2.set_string("Wm2");
  auto nondim = m/m;  // dummy unit for non-dimensional fields
  auto micron = m / 1000000;

  using namespace ShortFieldTagsNames;

  const auto& grid_name = m_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);
  m_ncol = grid->get_num_local_dofs();
  m_nlay = grid->get_num_vertical_levels();
  m_lat  = grid->get_geometry_data("lat");
  m_lon  = grid->get_geometry_data("lon");

  // Set up dimension layouts
  FieldLayout scalar2d_layout     { {COL   }, {m_ncol    } };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_ncol,m_nlay} };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {m_ncol,m_nlay+1} };

  constexpr int ps = SCREAM_SMALL_PACK_SIZE;
  // Set required (input) fields here
  add_field<Required>("p_mid" , scalar3d_layout_mid, Pa, grid->name(), ps);
  add_field<Required>("p_int", scalar3d_layout_int, Pa, grid->name(), ps);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid->name(), ps);
  add_field<Required>("sfc_alb_dir_vis", scalar2d_layout, nondim, grid->name());
  add_field<Required>("sfc_alb_dir_nir", scalar2d_layout, nondim, grid->name());
  add_field<Required>("sfc_alb_dif_vis", scalar2d_layout, nondim, grid->name());
  add_field<Required>("sfc_alb_dif_nir", scalar2d_layout, nondim, grid->name());
  add_field<Required>("qc", scalar3d_layout_mid, kgkg, grid->name(), ps);
  add_field<Required>("qi", scalar3d_layout_mid, kgkg, grid->name(), ps);
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid->name(), ps);
  add_field<Required>("eff_radius_qc", scalar3d_layout_mid, micron, grid->name(), ps);
  add_field<Required>("eff_radius_qi", scalar3d_layout_mid, micron, grid->name(), ps);
  add_field<Required>("qv",scalar3d_layout_mid,kgkg,grid->name(), ps);
  add_field<Required>("surf_lw_flux_up",scalar2d_layout,W/(m*m),grid->name());
  // Set of required gas concentration fields
  for (auto& it : m_gas_names) {
    if (it == "h2o") { /* Special case where water vapor is called h2o in radiation */
      // do nothing, qv has already been added.
    } else {
      add_field<Required>(it,scalar3d_layout_mid,kgkg,grid->name(), ps);
    }
  }
  // Set computed (output) fields
  add_field<Updated >("T_mid"     , scalar3d_layout_mid, K  , grid->name(), ps);
  add_field<Computed>("SW_flux_dn", scalar3d_layout_int, Wm2, grid->name(), ps);
  add_field<Computed>("SW_flux_up", scalar3d_layout_int, Wm2, grid->name(), ps);
  add_field<Computed>("SW_flux_dn_dir", scalar3d_layout_int, Wm2, grid->name(), ps);
  add_field<Computed>("LW_flux_up", scalar3d_layout_int, Wm2, grid->name(), ps);
  add_field<Computed>("LW_flux_dn", scalar3d_layout_int, Wm2, grid->name(), ps);

}  // RRTMGPRadiation::set_grids

int RRTMGPRadiation::requested_buffer_size_in_bytes() const
{
  const int interface_request = Buffer::num_1d_ncol*m_ncol*sizeof(Real) +
                                Buffer::num_2d_nlay*m_ncol*m_nlay*sizeof(Real) +
                                Buffer::num_2d_nlay_p1*m_ncol*(m_nlay+1)*sizeof(Real) +
                                Buffer::num_2d_nswbands*m_ncol*m_nswbands*sizeof(Real);

  return interface_request;
} // RRTMGPRadiation::requested_buffer_size
// =========================================================================================

void RRTMGPRadiation::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d array
  m_buffer.mu0 = decltype(m_buffer.mu0)("mu0", mem, m_ncol);
  mem += m_buffer.mu0.totElems();
  m_buffer.sfc_alb_dir_vis = decltype(m_buffer.sfc_alb_dir_vis)("sfc_alb_dir_vis", mem, m_ncol);
  mem += m_buffer.sfc_alb_dir_vis.totElems();
  m_buffer.sfc_alb_dir_nir = decltype(m_buffer.sfc_alb_dir_nir)("sfc_alb_dir_nir", mem, m_ncol);
  mem += m_buffer.sfc_alb_dir_nir.totElems();
  m_buffer.sfc_alb_dif_vis = decltype(m_buffer.sfc_alb_dif_vis)("sfc_alb_dif_vis", mem, m_ncol);
  mem += m_buffer.sfc_alb_dif_vis.totElems();
  m_buffer.sfc_alb_dif_nir = decltype(m_buffer.sfc_alb_dif_nir)("sfc_alb_dif_nir", mem, m_ncol);
  mem += m_buffer.sfc_alb_dif_nir.totElems();
  m_buffer.cosine_zenith = decltype(m_buffer.cosine_zenith)(mem, m_ncol);
  mem += m_buffer.cosine_zenith.size();

  // 2d arrays
  m_buffer.p_lay = decltype(m_buffer.p_lay)("p_lay", mem, m_ncol, m_nlay);
  mem += m_buffer.p_lay.totElems();
  m_buffer.t_lay = decltype(m_buffer.t_lay)("t_lay", mem, m_ncol, m_nlay);
  mem += m_buffer.t_lay.totElems();
  m_buffer.p_del = decltype(m_buffer.p_del)("p_del", mem, m_ncol, m_nlay);
  mem += m_buffer.p_del.totElems();
  m_buffer.qc = decltype(m_buffer.qc)("qc", mem, m_ncol, m_nlay);
  mem += m_buffer.qc.totElems();
  m_buffer.qi = decltype(m_buffer.qi)("qi", mem, m_ncol, m_nlay);
  mem += m_buffer.qi.totElems();
  m_buffer.cldfrac_tot = decltype(m_buffer.cldfrac_tot)("cldfrac_tot", mem, m_ncol, m_nlay);
  mem += m_buffer.cldfrac_tot.totElems();
  m_buffer.eff_radius_qc = decltype(m_buffer.eff_radius_qc)("eff_radius_qc", mem, m_ncol, m_nlay);
  mem += m_buffer.eff_radius_qc.totElems();
  m_buffer.eff_radius_qi = decltype(m_buffer.eff_radius_qi)("eff_radius_qi", mem, m_ncol, m_nlay);
  mem += m_buffer.eff_radius_qi.totElems();
  m_buffer.tmp2d = decltype(m_buffer.tmp2d)("tmp2d", mem, m_ncol, m_nlay);
  mem += m_buffer.tmp2d.totElems();
  m_buffer.lwp = decltype(m_buffer.lwp)("lwp", mem, m_ncol, m_nlay);
  mem += m_buffer.lwp.totElems();
  m_buffer.iwp = decltype(m_buffer.iwp)("iwp", mem, m_ncol, m_nlay);
  mem += m_buffer.iwp.totElems();
  m_buffer.sw_heating = decltype(m_buffer.sw_heating)("sw_heating", mem, m_ncol, m_nlay);
  mem += m_buffer.sw_heating.totElems();
  m_buffer.lw_heating = decltype(m_buffer.lw_heating)("lw_heating", mem, m_ncol, m_nlay);
  mem += m_buffer.lw_heating.totElems();
  m_buffer.rad_heating = decltype(m_buffer.rad_heating)("rad_heating", mem, m_ncol, m_nlay);
  mem += m_buffer.rad_heating.totElems();

  m_buffer.p_lev = decltype(m_buffer.p_lev)("p_lev", mem, m_ncol, m_nlay+1);
  mem += m_buffer.p_lev.totElems();
  m_buffer.t_lev = decltype(m_buffer.t_lev)("t_lev", mem, m_ncol, m_nlay+1);
  mem += m_buffer.t_lev.totElems();
  m_buffer.sw_flux_up = decltype(m_buffer.sw_flux_up)("sw_flux_up", mem, m_ncol, m_nlay+1);
  mem += m_buffer.sw_flux_up.totElems();
  m_buffer.sw_flux_dn = decltype(m_buffer.sw_flux_dn)("sw_flux_dn", mem, m_ncol, m_nlay+1);
  mem += m_buffer.sw_flux_dn.totElems();
  m_buffer.sw_flux_dn_dir = decltype(m_buffer.sw_flux_dn_dir)("sw_flux_dn_dir", mem, m_ncol, m_nlay+1);
  mem += m_buffer.sw_flux_dn_dir.totElems();
  m_buffer.lw_flux_up = decltype(m_buffer.lw_flux_up)("lw_flux_up", mem, m_ncol, m_nlay+1);
  mem += m_buffer.lw_flux_up.totElems();
  m_buffer.lw_flux_dn = decltype(m_buffer.lw_flux_dn)("lw_flux_dn", mem, m_ncol, m_nlay+1);
  mem += m_buffer.lw_flux_dn.totElems();

  m_buffer.sfc_alb_dir = decltype(m_buffer.sfc_alb_dir)("sfc_alb_dir", mem, m_ncol, m_nswbands);
  mem += m_buffer.sfc_alb_dir.totElems();
  m_buffer.sfc_alb_dif = decltype(m_buffer.sfc_alb_dif)("sfc_alb_dif", mem, m_ncol, m_nswbands);
  mem += m_buffer.sfc_alb_dif.totElems();

  int used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for RRTMGPRadiation.");
} // RRTMGPRadiation::init_buffers

void RRTMGPRadiation::initialize_impl() {
  using PC = scream::physics::Constants<Real>;

  // Determine orbital year. If Orbital Year is negative, use current year
  // from timestamp for orbital year; if positive, use provided orbital year
  // for duration of simulation.
  m_orbital_year = m_params.get<Int>("Orbital Year",-9999);

  // Determine whether or not we are using a fixed solar zenith angle (positive value)
  m_fixed_solar_zenith_angle = m_params.get<Real>("Fixed Solar Zenith Angle", -9999);

  // Initialize yakl
  if(!yakl::isInitialized()) { yakl::init(); }

  // Names of active gases
  auto gas_names_yakl_offset = string1d("gas_names",m_ngas);
  m_gas_mol_weights          = view_1d_real("gas_mol_weights",m_ngas);
  /* the lookup function for getting the gas mol weights doesn't work on device. */
  auto gas_mol_w_host = Kokkos::create_mirror_view(m_gas_mol_weights);
  for (int igas = 0; igas < m_ngas; igas++) {  
    /* Note: YAKL starts the index from 1 */
    gas_names_yakl_offset(igas+1)   = m_gas_names[igas];
    gas_mol_w_host[igas]            = PC::get_gas_mol_weight(m_gas_names[igas]);
  }
  Kokkos::deep_copy(m_gas_mol_weights,gas_mol_w_host);
  // Initialize GasConcs object to pass to RRTMGP initializer;
  gas_concs.init(gas_names_yakl_offset,m_ncol,m_nlay);
  rrtmgp::rrtmgp_initialize(gas_concs);

}
// =========================================================================================

void RRTMGPRadiation::run_impl (const int dt) {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using PC = scream::physics::Constants<Real>;
  using CO = scream::ColumnOps<DefaultDevice,Real>;

  // get a host copy of lat/lon 
  auto h_lat  = Kokkos::create_mirror_view(m_lat);
  auto h_lon  = Kokkos::create_mirror_view(m_lon);
  Kokkos::deep_copy(h_lat,m_lat);
  Kokkos::deep_copy(h_lon,m_lon);

  // Get data from the FieldManager
  auto d_pmid = get_field_in("p_mid").get_view<const Real**>();
  auto d_pint = get_field_in("p_int").get_view<const Real**>();
  auto d_pdel = get_field_in("pseudo_density").get_view<const Real**>();
  auto d_sfc_alb_dir_vis = get_field_in("sfc_alb_dir_vis").get_view<const Real*>();
  auto d_sfc_alb_dir_nir = get_field_in("sfc_alb_dir_nir").get_view<const Real*>();
  auto d_sfc_alb_dif_vis = get_field_in("sfc_alb_dif_vis").get_view<const Real*>();
  auto d_sfc_alb_dif_nir = get_field_in("sfc_alb_dif_nir").get_view<const Real*>();
  auto d_qv = get_field_in("qv").get_view<const Real**>();
  auto d_qc = get_field_in("qc").get_view<const Real**>();
  auto d_qi = get_field_in("qi").get_view<const Real**>();
  auto d_cldfrac_tot = get_field_in("cldfrac_tot").get_view<const Real**>();
  auto d_rel = get_field_in("eff_radius_qc").get_view<const Real**>();
  auto d_rei = get_field_in("eff_radius_qi").get_view<const Real**>();
  auto d_surf_lw_flux_up = get_field_in("surf_lw_flux_up").get_view<const Real*>();
  auto d_tmid = get_field_out("T_mid").get_view<Real**>();
  auto d_sw_flux_up = get_field_out("SW_flux_up").get_view<Real**>();
  auto d_sw_flux_dn = get_field_out("SW_flux_dn").get_view<Real**>();
  auto d_sw_flux_dn_dir = get_field_out("SW_flux_dn_dir").get_view<Real**>();
  auto d_lw_flux_up = get_field_out("LW_flux_up").get_view<Real**>();
  auto d_lw_flux_dn = get_field_out("LW_flux_dn").get_view<Real**>();

  // Create YAKL arrays. RRTMGP expects YAKL arrays with styleFortran, i.e., data has ncol
  // as the fastest index. For this reason we must copy the data.
  auto p_lay           = m_buffer.p_lay;
  auto t_lay           = m_buffer.t_lay;
  auto p_lev           = m_buffer.p_lev;
  auto p_del           = m_buffer.p_del;
  auto t_lev           = m_buffer.t_lev;
  auto mu0             = m_buffer.mu0;
  auto sfc_alb_dir     = m_buffer.sfc_alb_dir;
  auto sfc_alb_dif     = m_buffer.sfc_alb_dif;
  auto sfc_alb_dir_vis = m_buffer.sfc_alb_dir_vis;
  auto sfc_alb_dir_nir = m_buffer.sfc_alb_dir_nir;
  auto sfc_alb_dif_vis = m_buffer.sfc_alb_dif_vis;
  auto sfc_alb_dif_nir = m_buffer.sfc_alb_dif_nir;
  auto qc              = m_buffer.qc;
  auto qi              = m_buffer.qi;
  auto cldfrac_tot     = m_buffer.cldfrac_tot;
  auto rel             = m_buffer.eff_radius_qc;
  auto rei             = m_buffer.eff_radius_qi;
  auto sw_flux_up      = m_buffer.sw_flux_up;
  auto sw_flux_dn      = m_buffer.sw_flux_dn;
  auto sw_flux_dn_dir  = m_buffer.sw_flux_dn_dir;
  auto lw_flux_up      = m_buffer.lw_flux_up;
  auto lw_flux_dn      = m_buffer.lw_flux_dn;

  constexpr auto stebol = PC::stebol;
  const auto ncol = m_ncol;
  const auto nlay = m_nlay;
  // Copy data from the FieldManager to the YAKL arrays
  {
    // Determine the cosine zenith angle
    // NOTE: Since we are bridging to F90 arrays this must be done on HOST and then
    //       deep copied to a device view.
    auto d_mu0 = m_buffer.cosine_zenith;
    auto h_mu0 = Kokkos::create_mirror_view(d_mu0);
    if (m_fixed_solar_zenith_angle > 0) {
      for (int i=0; i<m_ncol; i++) {
        h_mu0(i) = m_fixed_solar_zenith_angle;
      }
    } else {
      // First gather the orbital parameters:
      double eccen, obliq, mvelp, obliqr, lambm0, mvelpp;
      auto ts = timestamp();
      auto orbital_year = m_orbital_year;
      if (orbital_year < 0) {
          orbital_year = ts.get_year();
      }
      shr_orb_params_c2f(&orbital_year, &eccen, &obliq, &mvelp, 
                         &obliqr, &lambm0, &mvelpp);
      // Use the orbital parameters to calculate the solar declination
      double delta, eccf;
      auto calday = ts.frac_of_year_in_days();
      shr_orb_decl_c2f(calday, eccen, mvelpp, lambm0,
                       obliqr, &delta, &eccf);
      // Now use solar declination to calculate zenith angle for all points
      for (int i=0;i<m_ncol;i++) {
        double lat = h_lat(i)*PC::Pi/180.0;  // Convert lat/lon to radians
        double lon = h_lon(i)*PC::Pi/180.0;
        h_mu0(i) = shr_orb_cosz_c2f(calday, lat, lon, delta, dt);
      }
    }
    Kokkos::deep_copy(d_mu0,h_mu0);

    // dz and T_int will need to be computed
    view_2d_real d_tint("T_int", m_ncol, m_nlay+1);
    view_2d_real d_dz  ("dz",    m_ncol, m_nlay);

    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();

      // Calculate dz
      const auto pseudo_density = ekat::subview(d_pdel, i);
      const auto p_mid          = ekat::subview(d_pmid, i);
      const auto T_mid          = ekat::subview(d_tmid, i);
      const auto qv             = ekat::subview(d_qv,   i);
      const auto dz             = ekat::subview(d_dz,   i);
      PF::calculate_dz<Real>(team, pseudo_density, p_mid, T_mid, qv, dz);
      team.team_barrier();

      // Calculate T_int from longwave flux up from the surface, assuming
      // blackbody emission with emissivity of 1.
      // TODO: Does land model assume something other than emissivity of 1? If so
      // we should use that here rather than assuming perfect blackbody emission.
      // NOTE: RRTMGP can accept vertical ordering surface to toa, or toa to
      // surface. The input data for the standalone test is ordered surface to
      // toa, but SCREAM in general assumes data is toa to surface. We account
      // for this here by swapping bc_top and bc_bot in the case that the input
      // data is ordered surface to toa.
      const auto T_int = ekat::subview(d_tint, i);
      const auto P_mid = ekat::subview(d_pmid, i);
      const int itop = (P_mid(0) < P_mid(nlay-1)) ? 0 : nlay-1;
      const Real bc_top = T_mid(itop);
      const Real bc_bot = sqrt(sqrt(d_surf_lw_flux_up(i)/stebol));
      if (itop == 0) {
          CO::compute_interface_values_linear(team, nlay, T_mid, dz, bc_top, bc_bot, T_int);
      } else {
          CO::compute_interface_values_linear(team, nlay, T_mid, dz, bc_bot, bc_top, T_int);
      }
      team.team_barrier();

      mu0(i+1) = d_mu0(i);
      sfc_alb_dir_vis(i+1) = d_sfc_alb_dir_vis(i);
      sfc_alb_dir_nir(i+1) = d_sfc_alb_dir_nir(i);
      sfc_alb_dif_vis(i+1) = d_sfc_alb_dif_vis(i);
      sfc_alb_dif_nir(i+1) = d_sfc_alb_dif_nir(i);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlay), [&] (const int& k) {
        p_lay(i+1,k+1)       = d_pmid(i,k);
        t_lay(i+1,k+1)       = d_tmid(i,k);
        p_del(i+1,k+1)       = d_pdel(i,k);
        qc(i+1,k+1)          = d_qc(i,k);
        qi(i+1,k+1)          = d_qi(i,k);
        cldfrac_tot(i+1,k+1) = d_cldfrac_tot(i,k);
        rel(i+1,k+1)         = d_rel(i,k);
        rei(i+1,k+1)         = d_rei(i,k);
        p_lev(i+1,k+1)       = d_pint(i,k);
        t_lev(i+1,k+1)       = d_tint(i,k);
      });

      p_lev(i+1,nlay+1) = d_pint(i,nlay);
      t_lev(i+1,nlay+1) = d_tint(i,nlay);
    });
  }
  Kokkos::fence();

  // Populate GasConcs object to pass to RRTMGP driver
  auto tmp2d = m_buffer.tmp2d;
  for (int igas = 0; igas < m_ngas; igas++) {
    auto name = m_gas_names[igas];
    auto fm_name = name=="h2o" ? "qv" : name;
    auto d_temp  = get_field_in(fm_name).get_view<const Real**>();
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_nlay, m_ncol);
    const auto gas_mol_weights = m_gas_mol_weights;

    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int k = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ncol), [&] (const int& i) {
        tmp2d(i+1,k+1) = PF::calculate_vmr_from_mmr(gas_mol_weights[igas],d_qv(i,k),d_temp(i,k)); // Note that for YAKL arrays i and k start with index 1
      });
    });
    Kokkos::fence();

    gas_concs.set_vmr(name, tmp2d);
  }

  // Compute layer cloud mass (per unit area)
  auto lwp = m_buffer.lwp;
  auto iwp = m_buffer.iwp;
  scream::rrtmgp::mixing_ratio_to_cloud_mass(qc, cldfrac_tot, p_del, lwp);
  scream::rrtmgp::mixing_ratio_to_cloud_mass(qi, cldfrac_tot, p_del, iwp);
  // Convert to g/m2 (needed by RRTMGP)
  {
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_nlay, m_ncol);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int k = team.league_rank()+1; // Note that for YAKL arrays i and k start with index 1
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ncol), [&] (const int& icol) {
      int i = icol+1;
      lwp(i,k) *= 1e3;
      iwp(i,k) *= 1e3;
    });
  });
  }
  Kokkos::fence();

  // Compute band-by-band surface_albedos. This is needed since
  // the AD passes broadband albedos, but rrtmgp require band-by-band.
  rrtmgp::compute_band_by_band_surface_albedos(
    m_ncol, m_nswbands,
    sfc_alb_dir_vis, sfc_alb_dir_nir,
    sfc_alb_dif_vis, sfc_alb_dif_nir,
    sfc_alb_dir, sfc_alb_dif);

  // Run RRTMGP driver
  rrtmgp::rrtmgp_main(
    m_ncol, m_nlay,
    p_lay, t_lay, p_lev, t_lev,
    gas_concs,
    sfc_alb_dir, sfc_alb_dif, mu0,
    lwp, iwp, rel, rei,
    sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
    lw_flux_up, lw_flux_dn, get_comm().am_i_root()
  );

  // Compute and apply heating rates
  auto sw_heating  = m_buffer.sw_heating;
  auto lw_heating  = m_buffer.lw_heating;
  auto rad_heating = m_buffer.rad_heating;
  rrtmgp::compute_heating_rate(
    sw_flux_up, sw_flux_dn, p_del, sw_heating
  );
  rrtmgp::compute_heating_rate(
    lw_flux_up, lw_flux_dn, p_del, lw_heating
  );
  {
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_nlay, m_ncol);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int k = team.league_rank()+1; // Note that for YAKL arrays i and k start with index 1
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ncol), [&] (const int& icol) {
      int i = icol+1;
      rad_heating(i,k) = sw_heating(i,k) + lw_heating(i,k);
      t_lay(i,k) = t_lay(i,k) + rad_heating(i,k) * dt;
    });
  });
  }
  Kokkos::fence();

  // Copy ouput data back to FieldManager
  {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlay+1), [&] (const int& k) {
        if (k < nlay) d_tmid(i,k) = t_lay(i+1,k+1);

        d_sw_flux_up(i,k)     = sw_flux_up(i+1,k+1);
        d_sw_flux_dn(i,k)     = sw_flux_dn(i+1,k+1);
        d_sw_flux_dn_dir(i,k) = sw_flux_dn_dir(i+1,k+1);
        d_lw_flux_up(i,k)     = lw_flux_up(i+1,k+1);
        d_lw_flux_dn(i,k)     = lw_flux_dn(i+1,k+1);
      });
    });
  }
}
// =========================================================================================

void RRTMGPRadiation::finalize_impl  () {
  gas_concs.reset();
  rrtmgp::rrtmgp_finalize();

  // Finalize YAKL
  yakl::finalize();
}
// =========================================================================================


}  // namespace scream
