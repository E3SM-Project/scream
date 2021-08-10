#include "catch2/catch.hpp"

#include "control/atmosphere_driver.hpp"
#include "dynamics/homme/atmosphere_dynamics.hpp"
#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

// Hommexx includes
#include "Context.hpp"
#include "SimulationParams.hpp"
#include "Types.hpp"
#include "FunctorsBuffersManager.hpp"

#include "ekat/util/ekat_feutils.hpp"
#include "ekat/ekat_assert.hpp"

static int get_default_fpes () {
#ifdef SCREAM_FPE
  return (FE_DIVBYZERO |
          FE_INVALID   |
          FE_OVERFLOW);
#else
  return 0;
#endif
}

TEST_CASE("scream_homme_stand_alone", "scream_homme_stand_alone") {
  using namespace scream;
  using namespace scream::control;

  ekat::enable_fpes(get_default_fpes());

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("dynamics",&create_atmosphere_process<HommeDynamics>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("Dynamics Driven",create_dynamics_driven_grids_manager);

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Create the driver
  AtmosphereDriver ad;

  // Init, run, and finalize
  util::TimeStamp time (0,0,0,0);
  ad.initialize(atm_comm,ad_params,time);

  // Add checks to verify AD memory buffer and Homme FunctorsBuffersManager
  // are the same size and reference the same memory.
  auto& fbm  = Homme::Context::singleton().get<Homme::FunctorsBuffersManager>();
  auto& memory_buffer = ad.get_memory_buffer();
  EKAT_ASSERT_MSG(fbm.allocated_size()*sizeof(Real) == (long unsigned int)memory_buffer.allocated_bytes(),
                  "Error! AD memory buffer and Homme FunctorsBuffersManager have mismatched sizes.");
  EKAT_ASSERT_MSG(fbm.get_memory() == memory_buffer.get_memory(),
                  "Error! AD memory buffer and Homme FunctorsBuffersManager reference different memory.");

  // Have to wait till now to get homme's parameters, cause it only gets init-ed during the driver initialization
  const auto& sp = Homme::Context::singleton().get<Homme::SimulationParams>();
  const int nmax = get_homme_param<int>("nmax");
  const int num_dyn_iters = nmax / (sp.qsplit*sp.rsplit);
  const double dt = get_homme_param<Real>("dt");

  EKAT_ASSERT_MSG (dt>0, "Error! Time step must be positive.\n");

  for (int i=0; i<num_dyn_iters; ++i) {
    ad.run(dt);
  }
  ad.finalize();

  // If we got here, we were able to run homme
  REQUIRE(true);
}
