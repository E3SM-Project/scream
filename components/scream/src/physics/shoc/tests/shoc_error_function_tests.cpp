#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestErrorFunction {

  static void run_property()
  {
    static constexpr int num_tests = 5;
    static constexpr Real inputs[num_tests] = {0, 100, -50, 0.25, -1.3};

    for (int i=0; i<num_tests; ++i) {
      ErrorFunctionData SDS;

      // fill in data
      SDS.input = inputs[i];

      // Call the fortran implementation
      error_function(SDS);

      REQUIRE(SDS.output >= -1);
      REQUIRE(SDS.output <= 1);
    }
  }

  static void run_bfb()
  {
    static constexpr Int num_runs = 5;

    ErrorFunctionData SDS_f90[num_runs];
    ErrorFunctionData SDS_cxx[num_runs];

    // Generate random input data
    std::default_random_engine generator;
    std::uniform_real_distribution<Real> data_dist(-1, 1);
    for (auto& d : SDS_f90) {
      d.input = data_dist(generator);
    }

    // Copy data to C++ version
    for (Int i=0; i<num_runs; ++i) {
      SDS_cxx[i].input = SDS_f90[i].input;
    }

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      error_function(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      error_function_f(d.input, &d.output);
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      ErrorFunctionData& d_f90 = SDS_f90[i];
      ErrorFunctionData& d_cxx = SDS_cxx[i];
      REQUIRE(d_f90.output == d_cxx.output);
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("error_function_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestErrorFunction;

  TestStruct::run_property();
}

TEST_CASE("error_function_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestErrorFunction;

  TestStruct::run_bfb();
}

} // namespace
