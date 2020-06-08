#include "catch2/catch.hpp"

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/p3/p3_f90.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestRainSed {

static void run_phys_rain_vel()
{
  // TODO
}

static void run_phys_rain_sed()
{
  // TODO
}

static void run_phys()
{
  run_phys_rain_vel();
  run_phys_rain_sed();
}

static void run_bfb_rain_vel()
{
  // Read in tables
  view_2d_table vn_table; view_2d_table vm_table; view_2d_table revap_table;
  view_1d_table mu_r_table; view_dnu_table dnu;
  Functions::init_kokkos_tables(vn_table, vm_table, revap_table, mu_r_table, dnu);

  constexpr Scalar qsmall = C::QSMALL;

  // Load some lookup inputs, need at least one per pack value
  ComputeRainFallVelocityData crfv_fortran[max_pack_size] = {
    // qr_incld,       rcldm,    rhofacr,         nr, nr_incld
    {1.1030E-04, 1.0000E+00, 1.3221E+00, 6.5282E+07, 6.2964E+05},
    {2.1437E-13, 2.0000E+00, 1.0918E+00, 1.7935E+02, 6.5337E+07},
    {5.6298E-05, 3.0000E+00, 1.1129E+00, 6.2909E+05, 1.6576E+02},
    {1.0000E-02, 4.0000E+00, 1.0774E+00, 1.6986E+02, 1.9436E+02},

    {0.0       , 1.0000E+00, 1.3221E+00, 6.5282E+07, 6.2964E+05},
    {2.1437E-13, 2.0000E+00, 1.0918E+00, 1.7935E+02, 6.5337E+07},
    {0.0       , 3.0000E+00, 1.1129E+00, 6.2909E+05, 1.6576E+02},
    {1.0000E-02, 4.0000E+00, 1.0774E+00, 1.6986E+02, 1.9436E+02},

    {1.1030E-04, 1.0000E+00, 1.3221E+00, 6.5282E+07, 6.2964E+05},
    {2.1437E-13, 2.0000E+00, 1.0918E+00, 1.7935E+02, 6.5337E+07},
    {0.0       , 3.0000E+00, 1.1129E+00, 6.2909E+05, 1.6576E+02},
    {0.0       , 4.0000E+00, 1.0774E+00, 1.6986E+02, 1.9436E+02},

    {0.0       , 1.0000E+00, 1.3221E+00, 6.5282E+07, 6.2964E+05},
    {2.1437E-13, 2.0000E+00, 1.0918E+00, 1.7935E+02, 6.5337E+07},
    {5.6298E-05, 3.0000E+00, 1.1129E+00, 6.2909E+05, 1.6576E+02},
    {1.0000E-02, 4.0000E+00, 1.0774E+00, 1.6986E+02, 1.9436E+02},
  };

  // Sync to device, needs to happen before fortran calls so that
  // inout data is in original state
  view_1d<ComputeRainFallVelocityData> crfv_device("crfv", max_pack_size);
  const auto crfv_host = Kokkos::create_mirror_view(crfv_device);
  std::copy(&crfv_fortran[0], &crfv_fortran[0] + max_pack_size, crfv_host.data());
  Kokkos::deep_copy(crfv_device, crfv_host);

  // Get data from fortran
  for (Int i = 0; i < max_pack_size; ++i) {
    if (crfv_fortran[i].qr_incld > qsmall) {
      compute_rain_fall_velocity(crfv_fortran[i]);
    }
  }

  // Calc bulk rime from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack qr_incld, rcldm, rhofacr, nr, nr_incld;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      qr_incld[s] = crfv_device(vs).qr_incld;
      rcldm[s]    = crfv_device(vs).rcldm;
      rhofacr[s]  = crfv_device(vs).rhofacr;
      nr[s]       = crfv_device(vs).nr;
      nr_incld[s] = crfv_device(vs).nr_incld;
    }

    Smask gt_small(qr_incld > qsmall);
    Spack mu_r(0), lamr(0), V_qr(0), V_nr(0);
    Functions::compute_rain_fall_velocity(
      gt_small, vn_table, vm_table, qr_incld, rcldm, rhofacr, nr, nr_incld, mu_r, lamr, V_qr, V_nr);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      crfv_device(vs).nr       = nr[s];
      crfv_device(vs).nr_incld = nr_incld[s];
      crfv_device(vs).mu_r     = mu_r[s];
      crfv_device(vs).lamr     = lamr[s];
      crfv_device(vs).V_qr     = V_qr[s];
      crfv_device(vs).V_nr     = V_nr[s];
    }
  });

  // Sync back to host
  Kokkos::deep_copy(crfv_host, crfv_device);

  // Validate results
  for (Int s = 0; s < max_pack_size; ++s) {
    REQUIRE(crfv_fortran[s].nr       == crfv_host(s).nr);
    REQUIRE(crfv_fortran[s].nr_incld == crfv_host(s).nr_incld);
    REQUIRE(crfv_fortran[s].mu_r     == crfv_host(s).mu_r);
    REQUIRE(crfv_fortran[s].lamr     == crfv_host(s).lamr);
    REQUIRE(crfv_fortran[s].V_qr     == crfv_host(s).V_qr);
    REQUIRE(crfv_fortran[s].V_nr     == crfv_host(s).V_nr);
  }
}

static void run_bfb_rain_sed()
{
  const std::array< std::pair<Real, Real>, RainSedData::NUM_ARRAYS > ranges = {
    std::make_pair(4.056E-03, 1.153E+00), // rho_range
    std::make_pair(0,         1.),        // inv_rho (ignored)
    std::make_pair(8.852E-01, 1.069E+00), // rhofacr
    std::make_pair(1.000E+00, 1.100E+00), // rcldm
    std::make_pair(2.863E-05, 8.141E-03), // inv_dzq_range
    std::make_pair(1.221E-14, 2.708E-03), // qr_incld
    std::make_pair(5.164E-10, 2.293E-03), // qr
    std::make_pair(9.558E+04, 6.596E+05), // nr
    std::make_pair(9.538E+04, 6.596E+05), // nr_incld
    std::make_pair(6.774E-15, 2.293E-03), // mu_r
    std::make_pair(7.075E-08, 2.418E-03), // lamr
    std::make_pair(7.861E-11, 4.179E-06), // qr_tend
    std::make_pair(5.164E-10, 2.733E-03), // nr_tend
    std::make_pair(4.469E-14, 2.557E-03), // rflx
  };

  RainSedData rsds_fortran[] = {
    //        kts, kte, ktop, kbot, kdir,        dt,       odt, prt_liq, ranges
    RainSedData(1,  72,   27,   72,   -1, 1.800E+03, 5.556E-04,     0.0, ranges),
    RainSedData(1,  72,   72,   27,    1, 1.800E+03, 5.556E-04,     1.0, ranges),
    RainSedData(1,  72,   27,   27,   -1, 1.800E+03, 5.556E-04,     0.0, ranges),
    RainSedData(1,  72,   27,   27,    1, 1.800E+03, 5.556E-04,     2.0, ranges),
  };

  static constexpr Int num_runs = sizeof(rsds_fortran) / sizeof(RainSedData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  RainSedData rsds_cxx[num_runs] = {
    RainSedData(rsds_fortran[0]),
    RainSedData(rsds_fortran[1]),
    RainSedData(rsds_fortran[2]),
    RainSedData(rsds_fortran[3]),
  };

    // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    rain_sedimentation(rsds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    RainSedData& d = rsds_cxx[i];
    rain_sedimentation_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                         d.qr_incld, d.rho, d.inv_rho, d.rhofacr, d.rcldm, d.inv_dzq,
                         d.dt, d.odt,
                         d.qr, d.nr, d.nr_incld, d.mu_r, d.lamr, &d.prt_liq, d.rflx,
                         d.qr_tend, d.nr_tend);
  }

  for (Int i = 0; i < num_runs; ++i) {
    // Due to pack issues, we must restrict checks to the active k space
    Int start = std::min(rsds_fortran[i].kbot, rsds_fortran[i].ktop) - 1; // 0-based indx
    Int end   = std::max(rsds_fortran[i].kbot, rsds_fortran[i].ktop);     // 0-based indx
    for (Int k = start; k < end; ++k) {
      REQUIRE(rsds_fortran[i].qr[k]       == rsds_cxx[i].qr[k]);
      REQUIRE(rsds_fortran[i].nr[k]       == rsds_cxx[i].nr[k]);
      REQUIRE(rsds_fortran[i].nr_incld[k] == rsds_cxx[i].nr_incld[k]);
      REQUIRE(rsds_fortran[i].mu_r[k]     == rsds_cxx[i].mu_r[k]);
      REQUIRE(rsds_fortran[i].lamr[k]     == rsds_cxx[i].lamr[k]);
      REQUIRE(rsds_fortran[i].rflx[k]     == rsds_cxx[i].rflx[k]);
      REQUIRE(rsds_fortran[i].qr_tend[k]  == rsds_cxx[i].qr_tend[k]);
      REQUIRE(rsds_fortran[i].nr_tend[k]  == rsds_cxx[i].nr_tend[k]);
    }
    REQUIRE(rsds_fortran[i].rflx[end]     == rsds_cxx[i].rflx[end]);
    REQUIRE(rsds_fortran[i].prt_liq == rsds_cxx[i].prt_liq);
  }
}

static void run_bfb()
{
  run_bfb_rain_vel();
  run_bfb_rain_sed();
}

};

}
}
}

namespace {

TEST_CASE("p3_rain_sed", "[p3_functions]")
{
  using TRS = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRainSed;

  scream::p3::p3_init(true); // need fortran table data

  TRS::run_phys();
  TRS::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
