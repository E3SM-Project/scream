#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/p3/atmosphere_microphysics.hpp"

#include "share/field/field_manager.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"

#include <iomanip>

namespace scream {

// Helper function to check the total water mass
Real calculate_water_mass(const std::shared_ptr<const GridsManager>& grids_mgr,const std::shared_ptr<FieldManager>& field_mgr, const bool use_precip) {
  
  using PC = scream::physics::Constants<Real>;
  constexpr Real gravit = PC::gravit;

  const auto& grid = grids_mgr->get_grid("Point Grid");

  int ncol  = grid->get_num_local_dofs();
  int nlay  = grid->get_num_vertical_levels();
  int dof   = ncol*nlay;
  REQUIRE(field_mgr->has_field("pseudo_density"));
  auto d_pdel       = field_mgr->get_field("pseudo_density").get_view<Real**>();

  // Start adding water species to mass calculation
  Real total_mass = 0.0;
  if (field_mgr->has_field("qv")) {
    auto d_tmp         = field_mgr->get_field("qv").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }
  if (field_mgr->has_field("qr")) {
    auto d_tmp         = field_mgr->get_field("qr").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }
  if (field_mgr->has_field("qc")) {
    auto d_tmp         = field_mgr->get_field("qc").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }
  if (field_mgr->has_field("qi")) {
    auto d_tmp         = field_mgr->get_field("qi").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }
  if (field_mgr->has_field("qs")) {
    auto d_tmp         = field_mgr->get_field("qs").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }

  if (use_precip) {
    auto d_tmp_liq = field_mgr->get_field("precip_liq_surf").get_view<Real*>();
    auto d_tmp_ice = field_mgr->get_field("precip_ice_surf").get_view<Real*>();
    Real result;
    Kokkos::parallel_reduce("",ncol, KOKKOS_LAMBDA(const int& icol,Real& lsum) {
      lsum += (d_tmp_liq(icol) + d_tmp_ice(icol) ) * 1000.0;
    },result);
    total_mass -= result;
  }

  return total_mass;
}

Real calculate_water_mass(const std::shared_ptr<const GridsManager>& grids_mgr,const std::shared_ptr<FieldManager>& field_mgr) {
  return calculate_water_mass(grids_mgr,field_mgr,false);
}


TEST_CASE("p3-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Time stepping parameters
  auto& ts = ad_params.sublist("Time Stepping");
  const auto dt = ts.get<int>("Time Step");
  const auto start_date = ts.get<std::vector<int>>("Start Date");
  const auto start_time = ts.get<std::vector<int>>("Start Time");
  const auto nsteps     = ts.get<int>("Number of Steps");

  util::TimeStamp t0 (start_date, start_time);
  EKAT_ASSERT_MSG (t0.is_valid(), "Error! Invalid start date.\n");

  // Need to register products in the factory *before* we create any atm process or grids manager.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("p3",&create_atmosphere_process<P3Microphysics>);
  gm_factory.register_product("Mesh Free",&create_mesh_free_grids_manager);

  // Create the driver
  AtmosphereDriver ad;

  // Init and run
  ad.initialize(atm_comm,ad_params,t0);

  // Grab views of water species and pseudo_density
  const auto& grids_mgr = ad.get_grids_manager();
  const auto& grid = grids_mgr->get_grid("Point Grid");
  const auto& field_mgr = ad.get_field_mgr(grid->name());
  

  if (atm_comm.am_i_root()) {
    printf("Start time stepping loop...       [  0%%]\n");
  }
  for (int i=0; i<nsteps; ++i) {
    const auto& wm_init = calculate_water_mass(grids_mgr,field_mgr);
    ad.run(dt);
    const auto& wm_after = calculate_water_mass(grids_mgr,field_mgr);
    Real total_precip = 0.0;
    EKAT_REQUIRE_MSG(wm_init - (wm_after + total_precip) < 1.e-12, 
       "Error in water mass change: " + std::to_string(wm_init) + " != "
       + std::to_string(wm_after) + " + " + std::to_string(total_precip));

    if (atm_comm.am_i_root()) {
      std::cout << "  - Iteration " << std::setfill(' ') << std::setw(3) << i+1 << " completed";
      std::cout << "       [" << std::setfill(' ') << std::setw(3) << 100*(i+1)/nsteps << "%]\n";
    }
  }

  // TODO: get the field repo from the driver, and go get (one of)
  //       the output(s) of P3, to check its numerical value (if possible)

  // Finalize 
  ad.finalize();

  // If we got here, we were able to run p3
  REQUIRE(true);
}

} // empty namespace
