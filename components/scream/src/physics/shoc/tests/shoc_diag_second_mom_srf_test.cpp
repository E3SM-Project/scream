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
struct UnitWrap::UnitTest<D>::TestSecondMomSrf {

static void run_second_mom_srf_bfb()
{
  const std::array< std::pair<Real, Real>, MomSrfData::NUM_ARRAYS > ranges = {
    std::make_pair(-1.0, 1.0), // wthl
    std::make_pair(-1.0, 1.0), // uw
    std::make_pair(-1.0, 1.0), // vw
    std::make_pair( 0.0, 1.0), // ustar2
    std::make_pair( 0.0, 1.0), // wstar 
  };

  MomSrfData srf_fortran[] = {
    // shcol, ranges
    MomSrfData(1,  ranges),
    MomSrfData(1,  ranges),
    MomSrfData(1,  ranges),
    MomSrfData(1,  ranges),
  };

  static constexpr Int num_runs = sizeof(srf_fortran) / sizeof(MomSrfData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  MomSrfData srf_cxx[num_runs] = {
    MomSrfData(srf_fortran[0]),
    MomSrfData(srf_fortran[1]),
    MomSrfData(srf_fortran[2]),
    MomSrfData(srf_fortran[3]),
  };


  // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    shoc_diag_second_moments_srf(srf_fortran[i]);
  }


  for (Int i = 0; i < num_runs; ++i) {
    MomSrfData& d = srf_cxx[i];
    shoc_diag_second_moments_srf_f(d.shcol, d.wthl, d.uw, d.vw, d.ustar2, d.wstar);
  }


  for (Int i = 0; i < num_runs; ++i) {
    Int shcol = srf_cxx[i].shcol;
    for (Int k = 0; k < shcol; ++k) {
      REQUIRE(srf_fortran[i].wthl[k]   == srf_cxx[i].wthl[k]);
      REQUIRE(srf_fortran[i].uw[k]     == srf_cxx[i].uw[k]);
      REQUIRE(srf_fortran[i].vw[k]     == srf_cxx[i].vw[k]);
      REQUIRE(srf_fortran[i].ustar2[k] == srf_cxx[i].ustar2[k]);
      REQUIRE(srf_fortran[i].wstar[k]  == srf_cxx[i].wstar[k]);
    }
  }
}

static void run_second_mom_srf_phys()
{
    // TODO
}

};

}
}
}


namespace {
TEST_CASE("second_mom_srf", "shoc") {
  using TRS = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomSrf;

  TRS::run_second_mom_srf_phys();
  TRS::run_second_mom_srf_bfb();
}
}
