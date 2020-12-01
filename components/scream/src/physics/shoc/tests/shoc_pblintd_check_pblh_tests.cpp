#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestPblintdCheckPblh {

  static void run_bfb()
  {
    PblintdCheckPblhData f90_data[] = {
      PblintdCheckPblhData(36,  72, 73),
      PblintdCheckPblhData(72,  72, 73),
      PblintdCheckPblhData(128, 72, 73),
      PblintdCheckPblhData(256, 72, 73),      
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(PblintdCheckPblhData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    PblintdCheckPblhData cxx_data[] = {
      PblintdCheckPblhData(f90_data[0]),
      PblintdCheckPblhData(f90_data[1]),
      PblintdCheckPblhData(f90_data[2]),
      PblintdCheckPblhData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      pblintd_check_pblh(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      pblintd_check_pblh_f(d.shcol, d.nlev, d.nlevi, d.z, d.ustar, d.check, d.pblh);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      PblintdCheckPblhData& d_f90 = f90_data[i];
      PblintdCheckPblhData& d_cxx = cxx_data[i];
      for (Int k = 0; k < d_f90.total(d_f90.pblh); ++k) {
        REQUIRE(d_f90.total(d_f90.pblh) == d_cxx.total(d_cxx.pblh));
        REQUIRE(d_f90.pblh[k] == d_cxx.pblh[k]);
      }

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("pblintd_check_pblh_bfb", "[shoc]")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdCheckPblh;

  TestStruct::run_bfb();
}

} // empty namespace
