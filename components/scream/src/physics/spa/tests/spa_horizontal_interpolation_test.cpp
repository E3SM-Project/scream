#include "catch2/catch.hpp"

#include "spa_unit_tests_common.hpp"

#include "share/scream_types.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/spa/spa_functions.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace spa {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestHorizontalInterp {
  static void run_test()
  {
    using namespace scorpio;
    ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
    MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
    eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

    std::string remap_filename = "map_ne2np4_to_ne4np4_mono.nc";
    std::string data_filename  = "spa_file_unified_and_clipped_ne2_scream.nc";

    SPAFunctions::horizontal_interpolation(remap_filename); 
    REQUIRE(true);
  }
};

}  // namespace unit_test
}  // namespace spa
}  // namespace scream

namespace{

TEST_CASE("spa_horizontal_interpolation", "spa")
{
  using TestStruct = scream::spa::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestHorizontalInterp;

  TestStruct::run_test();
}

} // namespace
