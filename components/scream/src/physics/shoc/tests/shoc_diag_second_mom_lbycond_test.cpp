#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
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
struct UnitWrap::UnitTest<D>::TestSecondMomLbycond {

static void run_second_mom_lbycond_bfb()
{
  SHOCSecondMomLbycondData lby_fortran[] = {
    // shcol
    SHOCSecondMomLbycondData(120),
    SHOCSecondMomLbycondData(120),
    SHOCSecondMomLbycondData(120),
    SHOCSecondMomLbycondData(120),
  };

  static constexpr Int num_runs = sizeof(lby_fortran) / sizeof(SHOCSecondMomLbycondData);

  for (Int i = 0; i < num_runs; ++i) {
    lby_fortran[i].randomize({ {-1, 1} });
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  SHOCSecondMomLbycondData lby_cxx[num_runs] = {
    SHOCSecondMomLbycondData(lby_fortran[0]),
    SHOCSecondMomLbycondData(lby_fortran[1]),
    SHOCSecondMomLbycondData(lby_fortran[2]),
    SHOCSecondMomLbycondData(lby_fortran[3]),
  };

  // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    shoc_diag_second_moments_lbycond(lby_fortran[i]);
  }

  for (Int i = 0; i < num_runs; ++i) {
    SHOCSecondMomLbycondData& d = lby_cxx[i];
    shoc_diag_second_moments_lbycond_f(d.shcol, d.wthl, d.wqw, d.uw, d.vw, d.ustar2, d.wstar,
                                       d.wthlo, d.wqwo, d.uwo, d.vwo, d.wtkeo, d.thlo, d.qwo, d.qwthlo);
  }

  for (Int i = 0; i < num_runs; ++i) {
    Int shcol      = lby_cxx[i].shcol;
    for (Int k = 0; k < shcol; ++k) {
      REQUIRE(lby_fortran[i].wthlo[k]  == lby_cxx[i].wthlo[k]);
      REQUIRE(lby_fortran[i].wqwo[k]   == lby_cxx[i].wqwo[k]);
      REQUIRE(lby_fortran[i].uwo[k]    == lby_cxx[i].uwo[k]);
      REQUIRE(lby_fortran[i].vwo[k]    == lby_cxx[i].vwo[k]);
      REQUIRE(lby_fortran[i].wtkeo[k]  == lby_cxx[i].wtkeo[k]);
      REQUIRE(lby_fortran[i].thlo[k]   == lby_cxx[i].thlo[k]);
      REQUIRE(lby_fortran[i].qwo[k]    == lby_cxx[i].qwo[k]);
      REQUIRE(lby_fortran[i].qwthlo[k] == lby_cxx[i].qwthlo[k]);
    }
  }
}

static void run_second_mom_lbycond_phys()
{
    // TODO
}

};

}
}
}

namespace {
TEST_CASE("second_mom_lby", "shoc") {
  using TRS = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomLbycond;

  TRS::run_second_mom_lbycond_phys();
  TRS::run_second_mom_lbycond_bfb();
}
}
