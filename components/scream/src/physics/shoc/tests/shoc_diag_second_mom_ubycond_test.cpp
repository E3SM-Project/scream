#include "catch2/catch.hpp"

//#include "share/scream_types.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <thread>

#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_arch.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/util/scream_utils.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/shoc/shoc_f90.hpp"
#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestSecondMomUbycond {

static void run_second_mom_ubycond_bfb()
{
  SHOCSecondMomentUbycondData uby_fortran[] = {
    // shcol, num_tracer, ranges
    SHOCSecondMomentUbycondData(128, 20),
    SHOCSecondMomentUbycondData(128, 20),
    SHOCSecondMomentUbycondData(128, 20),
    SHOCSecondMomentUbycondData(128, 20),
  };

  static constexpr Int num_runs = sizeof(uby_fortran) / sizeof(SHOCSecondMomentUbycondData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  SHOCSecondMomentUbycondData uby_cxx[num_runs] = {
    SHOCSecondMomentUbycondData(uby_fortran[0]),
    SHOCSecondMomentUbycondData(uby_fortran[1]),
    SHOCSecondMomentUbycondData(uby_fortran[2]),
    SHOCSecondMomentUbycondData(uby_fortran[3]),
  };

  // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    shoc_diag_second_moments_ubycond(uby_fortran[i]);
  }

  for (Int i = 0; i < num_runs; ++i) {
    SHOCSecondMomentUbycondData& d = uby_cxx[i];
    shoc_diag_second_moments_ubycond_f(d.shcol, d.num_tracer, d.thl, d.qw, d.qwthl, d.wthl, d.wqw, d.uw, d.vw,
        d.wtke, d.wtracer);
  }


  for (Int i = 0; i < num_runs; ++i) {
    Int shcol = uby_cxx[i].shcol;
    Int num_tracer = uby_cxx[i].num_tracer;
    for (Int k = 0; k < shcol; ++k) {
      REQUIRE(uby_fortran[i].thl[k]   == uby_cxx[i].thl[k]);
      REQUIRE(uby_fortran[i].qw[k]    == uby_cxx[i].qw[k]);
      REQUIRE(uby_fortran[i].qwthl[k] == uby_cxx[i].qwthl[k]);
      REQUIRE(uby_fortran[i].wthl[k]  == uby_cxx[i].wthl[k]);
      REQUIRE(uby_fortran[i].wqw[k]   == uby_cxx[i].wqw[k]);
      REQUIRE(uby_fortran[i].uw[k]    == uby_cxx[i].uw[k]);
      REQUIRE(uby_fortran[i].vw[k]    == uby_cxx[i].vw[k]);
      REQUIRE(uby_fortran[i].wtke[k]  == uby_cxx[i].wtke[k]);
      for (Int j = 0; j < num_tracer; ++j) {
        REQUIRE(uby_fortran[i].wtracer[j*k]  == uby_cxx[i].wtracer[j*k]);
      }
    }
  }
}

static void run_second_mom_ubycond_phys()
{
    // TODO
}
};

}
}
}

namespace {
TEST_CASE("second_mom_uby", "shoc") {
  using TRS = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomUbycond;
  TRS::run_second_mom_ubycond_phys();
  TRS::run_second_mom_ubycond_bfb();
}
}
