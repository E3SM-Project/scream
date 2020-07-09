#include <catch2/catch.hpp>
#include "ekat/scream_pack.hpp"
#include "ekat/scream_parse_yaml_file.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/se_grid.hpp"
#include "control/atmosphere_driver.hpp"

#include "physics/zm/atmosphere_macrophysics.hpp"
#include "physics/zm/scream_zm_interface.hpp"

#include <iostream>
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

protected:
};

TEST_CASE("zm-standalone", "") {
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
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("ZM",&create_atmosphere_process<ZMMacrophysics>);

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
    std :: cout << "this one's in the test cases\n";     
}

  // Finalize 
  ad.finalize();
  upgm.clean_up();

  // If we got here, we were able to run zm 
  REQUIRE(true);
}

} // empty namespace
