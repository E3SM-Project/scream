#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>
#include <iomanip>      // std::setprecision

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3IceMelting
{

static void ice_melting_bfb(){

  static constexpr Int max_pack_size = 16;
  REQUIRE(Spack::n <= max_pack_size);

  // make array of input data (why not pass actual variables?). Copied 1st 4 rows 4x to fill pack size.
  IceMeltingData IceMelt[max_pack_size] = {
    //rho,     t,        pres,     rhofaci,  f1pr05,   f1pr14,   xxlv,     xlf,      dv,       sc,       mu,       kap,      qv,       qitot_incld,nitot_incld
    {0.117E+01,0.299E+03,0.101E+06,0.829E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.263E-04,0.601E+00,0.185E-04,0.261E-01,0.160E-01,0.510E-02,  0.195E-12},
    {0.114E+01,0.296E+03,0.973E+05,0.842E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.268E-04,0.601E+00,0.183E-04,0.259E-01,0.149E-01,0.510E-02,  0.195E-12},
    {0.977E+00,0.287E+03,0.809E+05,0.913E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.306E-04,0.599E+00,0.179E-04,0.253E-01,0.827E-02,0.000E+00,  0.000E+00},
    {0.103E+01,0.289E+03,0.862E+05,0.887E+00,0.636E-03,0.281E-04,0.250E+07,0.334E+06,0.291E-04,0.600E+00,0.180E-04,0.254E-01,0.107E-01,0.510E-02,  0.336E+05},

    {0.117E+01,0.299E+03,0.101E+06,0.829E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.263E-04,0.601E+00,0.185E-04,0.261E-01,0.160E-01,0.510E-02,  0.195E-12},
    {0.114E+01,0.296E+03,0.973E+05,0.842E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.268E-04,0.601E+00,0.183E-04,0.259E-01,0.149E-01,0.510E-02,  0.195E-12},
    {0.977E+00,0.287E+03,0.809E+05,0.913E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.306E-04,0.599E+00,0.179E-04,0.253E-01,0.827E-02,0.000E+00,  0.000E+00},
    {0.103E+01,0.289E+03,0.862E+05,0.887E+00,0.636E-03,0.281E-04,0.250E+07,0.334E+06,0.291E-04,0.600E+00,0.180E-04,0.254E-01,0.107E-01,0.510E-02,  0.336E+05},

    {0.117E+01,0.299E+03,0.101E+06,0.829E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.263E-04,0.601E+00,0.185E-04,0.261E-01,0.160E-01,0.510E-02,  0.195E-12},
    {0.114E+01,0.296E+03,0.973E+05,0.842E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.268E-04,0.601E+00,0.183E-04,0.259E-01,0.149E-01,0.510E-02,  0.195E-12},
    {0.977E+00,0.287E+03,0.809E+05,0.913E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.306E-04,0.599E+00,0.179E-04,0.253E-01,0.827E-02,0.000E+00,  0.000E+00},
    {0.103E+01,0.289E+03,0.862E+05,0.887E+00,0.636E-03,0.281E-04,0.250E+07,0.334E+06,0.291E-04,0.600E+00,0.180E-04,0.254E-01,0.107E-01,0.510E-02,  0.336E+05},

    {0.117E+01,0.299E+03,0.101E+06,0.829E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.263E-04,0.601E+00,0.185E-04,0.261E-01,0.160E-01,0.510E-02,  0.195E-12},
    {0.114E+01,0.296E+03,0.973E+05,0.842E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.268E-04,0.601E+00,0.183E-04,0.259E-01,0.149E-01,0.510E-02,  0.195E-12},
    {0.977E+00,0.287E+03,0.809E+05,0.913E+00,0.122E+01,0.562E-01,0.250E+07,0.334E+06,0.306E-04,0.599E+00,0.179E-04,0.253E-01,0.827E-02,0.000E+00,  0.000E+00},
    {0.103E+01,0.289E+03,0.862E+05,0.887E+00,0.636E-03,0.281E-04,0.250E+07,0.334E+06,0.291E-04,0.600E+00,0.180E-04,0.254E-01,0.107E-01,0.510E-02,  0.336E+05}
  };

    //expected output
    // qimlt  ,nimlt
    //{0.717E-16,0.274E-26},
    //{0.626E-16,0.239E-26},
    //{0.000E+00,0.000E+00},
    //{0.333E-02,0.219E+05},
    
  // Sync to device
  view_1d<IceMeltingData> IceMelt_device("IceMelt", Spack::n);
  auto IceMelt_host = Kokkos::create_mirror_view(IceMelt_device);
  // This copy only copies the input variables.
  std::copy(&IceMelt[0], &IceMelt[0] + Spack::n, IceMelt_host.data());
  Kokkos::deep_copy(IceMelt_device, IceMelt_host);

  // Get data from fortran
  for (Int i = 0; i < Spack::n; ++i) {
    ice_melting(IceMelt[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
    // Init pack inputs
      Spack rho,t,pres,rhofaci,f1pr05,f1pr14,xxlv,xlf,dv,sc,mu,kap,qv,qitot_incld,nitot_incld,qimlt,nimlt;
    for (Int s = 0; s < Spack::n; ++s) {
      rho[s] = IceMelt_device(s).rho;
      t[s] = IceMelt_device(s).t;
      pres[s] = IceMelt_device(s).pres;
      rhofaci[s] = IceMelt_device(s).rhofaci;
      f1pr05[s] = IceMelt_device(s).f1pr05;
      f1pr14[s] = IceMelt_device(s).f1pr14;
      xxlv[s] = IceMelt_device(s).xxlv;
      xlf[s] = IceMelt_device(s).xlf;
      dv[s] = IceMelt_device(s).dv;
      sc[s] = IceMelt_device(s).sc;
      mu[s] = IceMelt_device(s).mu;
      kap[s] = IceMelt_device(s).kap;
      qv[s] = IceMelt_device(s).qv;
      qitot_incld[s] = IceMelt_device(s).qitot_incld;
      nitot_incld[s] = IceMelt_device(s).nitot_incld;
      qimlt[s] = IceMelt_device(s).qimlt;
      nimlt[s] = IceMelt_device(s).nimlt;
    }

    Functions::ice_melting(rho,t,pres,rhofaci,f1pr05,f1pr14,xxlv,xlf,dv,sc,mu,kap,qv,qitot_incld,nitot_incld,qimlt,nimlt);
    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      IceMelt_device(s).qimlt = qimlt[s];
      IceMelt_device(s).nimlt = nimlt[s];
    }

    });

  // Sync back to host
  Kokkos::deep_copy(IceMelt_host, IceMelt_device);

  // Validate results
  for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(IceMelt[s].qimlt == IceMelt_host(s).qimlt);
      REQUIRE(IceMelt[s].nimlt == IceMelt_host(s).nimlt);
  }
}; // TestP3IceMelting
  
}; // UnitWrap
  
} // namespace unit_test
} // namespace p3
} // namespace scream
 
namespace{

TEST_CASE("p3_ice_melting_test", "[p3_ice_melting_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3IceMelting::ice_melting_bfb();
}

} // namespace

