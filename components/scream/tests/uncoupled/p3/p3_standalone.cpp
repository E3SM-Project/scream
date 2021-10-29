#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/p3/atmosphere_microphysics.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"

#include <iomanip>

namespace scream {

TEST_CASE("p3-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Load ad parameter list
  std::string fname = "input_np" + std::to_string(atm_comm.size()) + ".yaml";
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
  printf("Start time stepping loop...       [  0%%]\n");
  for (int i=0; i<nsteps; ++i) {
    ad.run(dt);
    std::cout << "  - Iteration " << std::setfill(' ') << std::setw(3) << i+1 << " completed";
    std::cout << "       [" << std::setfill(' ') << std::setw(3) << 100*(i+1)/nsteps << "%]\n";
  }

  // TODO: get the field repo from the driver, and go get (one of)
  //       the output(s) of P3, to check its numerical value (if possible)

  // Finalize 
  ad.finalize();

  // If we got here, we were able to run p3
  REQUIRE(true);
}

} // empty namespace
