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
struct UnitWrap::UnitTest<D>::TestAdvSgsTke {

  static void run_bfb()
  {
    AdvSgsTkeData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(AdvSgsTkeData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    AdvSgsTkeData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      adv_sgs_tke(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      adv_sgs_tke_f(d.nlev(), d.shcol(), d.dtime, d.shoc_mix, d.wthv_sec, d.sterm_zt, d.tk, d.tke, d.a_diss);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      AdvSgsTkeData& d_f90 = f90_data[i];
      AdvSgsTkeData& d_cxx = cxx_data[i];
      for (Int k = 0; k < d_f90.total(d_f90.tke); ++k) {
        REQUIRE(d_f90.total(d_f90.tke) == d_cxx.total(d_cxx.tke));
        REQUIRE(d_f90.tke[k] == d_cxx.tke[k]);
        REQUIRE(d_f90.total(d_f90.tke) == d_cxx.total(d_cxx.a_diss));
        REQUIRE(d_f90.a_diss[k] == d_cxx.a_diss[k]);
      }

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("adv_sgs_tke_bfb", "[shoc]")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestAdvSgsTke;

  TestStruct::run_bfb();
}

} // empty namespace
