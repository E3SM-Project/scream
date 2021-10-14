#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/shoc/atmosphere_macrophysics.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_parse_yaml_file.hpp"

namespace scream {

TEST_CASE("shoc-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Load ad parameter list
  std::string fname = "input_np" + std::to_string(atm_comm.size()) + ".yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // run params:
  const auto& num_iters = ad_params.get<int>("Number of Iterations",5);
  const auto& dt        = ad_params.get<Real>("dt",300.0);

  // Need to register products in the factory *before* we create any atm process or grids manager.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("shoc",&create_atmosphere_process<SHOCMacrophysics>);
  gm_factory.register_product("Mesh Free",&create_mesh_free_grids_manager);

  // Create the driver
  AtmosphereDriver ad;

  // Init and run
  util::TimeStamp time (0,0,0,0);
  ad.initialize(atm_comm,ad_params,time);

  for (int i=0; i<num_iters; ++i) {
    if (i % 10 == 0) {
      printf("  -  %5.2f%%\nRun iteration: %d, ",Real(i)/Real(num_iters)*100,i+1);
    } else {
      printf("%d, ",i+1);
    }
    ad.run(dt);
  }
  printf("\n");

  // TODO: get the field repo from the driver, and go get (one of)
  //       the output(s) of SHOC, to check its numerical value (if possible)

  // Finalize 
  ad.finalize();

  // If we got here, we were able to run shoc
  REQUIRE(true);
}

} // empty namespace
