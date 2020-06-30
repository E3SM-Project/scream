#include <catch2/catch.hpp>
#include "ekat/scream_pack.hpp"
#include "ekat/scream_parse_yaml_file.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/se_grid.hpp"
#include "control/atmosphere_driver.hpp"

#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/p3/scream_p3_interface.hpp"
#include "physics/p3/p3_functions_f90.hpp"

namespace scream {

// === A dummy physics grids for this test === //

class DummyPhysicsGrid : public SEGrid
{
public:
  DummyPhysicsGrid (const int num_cols)
   : SEGrid("Physics",GridType::SE_NodeBased,num_cols)
  {
    // Nothing to do here
  }
  ~DummyPhysicsGrid () = default;
};

TEST_CASE("p3-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;

  constexpr int num_iters = 10;
  constexpr int num_cols  = 32;

  // Load ad parameter list
  std::string fname = "input.yaml";
  ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("P3",&create_atmosphere_process<P3Microphysics>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("User Provided",create_user_provided_grids_manager);

  // Set the dummy grid in the UserProvidedGridManager
  // Recall that this class stores *static* members, so whatever
  // we set here, will be reflected in the GM built by the factory.
  UserProvidedGridsManager upgm;
  upgm.set_grid(std::make_shared<DummyPhysicsGrid>(num_cols));
  upgm.set_reference_grid("Physics");

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  // Create the driver
  AtmosphereDriver ad;

  // Init and run (do not finalize, or you'll clear the field repo!)
  util::TimeStamp time (0,0,0,0);
  ad.initialize(atm_comm,ad_params,time);
  for (int i=0; i<num_iters; ++i) {
    ad.run(300.0);
  }

  // TODO: get the field repo from the driver, and go get (one of)
  //       the output(s) of P3, to check its numerical value (if possible)

  // Finalize 
  ad.finalize();
  upgm.clean_up();
  p3::P3GlobalForFortran::deinit();

  // If we got here, we were able to run p3
  REQUIRE(true);
}

} // empty namespace
