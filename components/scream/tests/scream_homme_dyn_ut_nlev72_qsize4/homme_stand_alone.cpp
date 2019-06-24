#include "catch2/catch.hpp"

#include "control/atmosphere_driver.hpp"
#include "control/surface_coupling.hpp"
#include "dynamics/homme/atmosphere_dynamics.hpp"
#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/scream_homme_interface.hpp"

// Hommexx includes
#include "Context.hpp"
#include "SimulationParams.hpp"
#include "Types.hpp"

TEST_CASE("scream_homme_stand_alone", "scream_homme_stand_alone") {
  using namespace scream;
  using namespace scream::control;

  // Create a parameter list for inputs
  ParameterList ad_params("Atmosphere Driver");
  auto& params = ad_params.sublist("Atmosphere Processes");

  params.set("Number of Entries",1);
  params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = params.sublist("Process 0");
  p0.set<std::string>("Process Name", "dynamics");

  auto& gm_params = ad_params.sublist("Grids Manager");
  gm_params.set<std::string>("Type","Dynamics Driven");

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("dynamics",&create_atmosphere_process<HommeDynamics>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("Dynamics Driven",create_dynamics_driven_grids_manager);

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
  util::TimeStamp time (0,0,0);
  ad.initialize(atm_comm,ad_params,time);

  // Have to wait till now to get homme's parameters, cause it only gets init-ed during the driver initialization
  const auto& sp = Homme::Context::singleton().get_simulation_params();
  const int nmax = get_homme_param_value<int>("nmax");
  const int num_dyn_iters = nmax / (sp.qsplit*sp.rsplit);
  const double dt = get_homme_param_value<Real>("dt");

  for (int i=0; i<num_dyn_iters; ++i) {
    ad.run(dt);
  }
  ad.finalize();

  // If we got here, we were able to run homme
  REQUIRE(true);
}
