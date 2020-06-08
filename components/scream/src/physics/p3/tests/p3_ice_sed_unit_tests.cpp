#include "catch2/catch.hpp"

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestIceSed {

static void run_phys_calc_bulk_rhime()
{
  // TODO
}

static void run_phys_ice_sed()
{
  // TODO
}

static void run_phys_homogeneous_freezing()
{
  // TODO
}

static void run_phys()
{
  run_phys_calc_bulk_rhime();
  run_phys_ice_sed();
  run_phys_homogeneous_freezing();
}

static void run_bfb_calc_bulk_rhime()
{
  constexpr Scalar qsmall = C::QSMALL;

  // Load some lookup inputs, need at least one per pack value
  CalcBulkRhoRimeData cbrr_fortran[max_pack_size] = {
    //     qi_tot,       qi_rim,       bi_rim
    {9.999978E-08, 9.999978E-03, 1.111108E-10},
    {0.000000E+00, 8.571428E-05, 1.000000E-02},
    {1.800685E-12, 1.818806E-13, 6.272458E-12},
    {5.164017E-10, 0.000000E+00, 0.000000E+00},

    {9.999978E-08, 0.000000E+00, 1.111108E-10},
    {5.100000E-03, 8.571428E-05, 1.000000E-02},
    {0.000000E+00, 1.818806E-13, 6.272458E-12},
    {5.164017E-10, 0.000000E+00, 0.000000E+00},

    {9.999978E-08, 9.999978E-08, 1.111108E-10},
    {5.100000E-03, 0.000000E+00, 1.000000E-02},
    {1.800685E-12, 1.818806E-13, 6.272458E-17},
    {0.000000E+00, 1.818806E-13, 0.000000E+00},

    {0.000000E+00, 9.999978E-08, 1.111108E-17},
    {5.100000E-03, 8.571428E-05, 1.000000E-02},
    {0.000000E+00, 1.818806E-13, 6.272458E-12},
    {5.164017E-10, 0.000000E+00, 0.000000E+00},
  };

  // Sync to device, needs to happen before fortran calls so that
  // inout data is in original state
  view_1d<CalcBulkRhoRimeData> cbrr_device("cbrr", max_pack_size);
  const auto cbrr_host = Kokkos::create_mirror_view(cbrr_device);
  std::copy(&cbrr_fortran[0], &cbrr_fortran[0] + max_pack_size, cbrr_host.data());
  Kokkos::deep_copy(cbrr_device, cbrr_host);

  // Get data from fortran
  for (Int i = 0; i < max_pack_size; ++i) {
    if (cbrr_fortran[i].qi_tot > qsmall) {
      calc_bulk_rho_rime(cbrr_fortran[i]);
    }
  }

  // Calc bulk rime from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack qi_tot, qi_rim, bi_rim;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      qi_tot[s] = cbrr_device(vs).qi_tot;
      qi_rim[s] = cbrr_device(vs).qi_rim;
      bi_rim[s] = cbrr_device(vs).bi_rim;
    }

    Smask gt_small(qi_tot > qsmall);
    Spack rho_rime = Functions::calc_bulk_rho_rime(gt_small, qi_tot, qi_rim, bi_rim);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      cbrr_device(vs).qi_rim   = qi_rim[s];
      cbrr_device(vs).bi_rim   = bi_rim[s];
      cbrr_device(vs).rho_rime = rho_rime[s];
    }
  });

  // Sync back to host
  Kokkos::deep_copy(cbrr_host, cbrr_device);

  // Validate results
  for (Int s = 0; s < max_pack_size; ++s) {
    REQUIRE(cbrr_fortran[s].qi_rim   == cbrr_host(s).qi_rim);
    REQUIRE(cbrr_fortran[s].bi_rim   == cbrr_host(s).bi_rim);
    REQUIRE(cbrr_fortran[s].rho_rime == cbrr_host(s).rho_rime);
  }
}

static void run_bfb_ice_sed()
{
  const std::array< std::pair<Real, Real>, IceSedData::NUM_ARRAYS > ranges = {
    std::make_pair(4.056E-03, 1.153E+00), // rho_range
    std::make_pair(0,         1.),        // inv_rho (ignored)
    std::make_pair(8.852E-01, 1.069E+00), // rhofaci
    std::make_pair(1.000E+00, 1.100E+00), // icldm
    std::make_pair(2.863E-05, 8.141E-03), // inv_dzq_range
    std::make_pair(1.221E-14, 2.708E-03), // qitot
    std::make_pair(5.164E-10, 2.293E-03), // qitot_incld
    std::make_pair(9.558E+04, 6.596E+05), // nitot
    std::make_pair(9.538E+04, 6.596E+05), // nitot_incld
    std::make_pair(6.774E-15, 2.293E-03), // qirim
    std::make_pair(7.075E-08, 2.418E-03), // qirim_incld
    std::make_pair(4.469E-14, 2.557E-03), // birim
    std::make_pair(7.861E-11, 4.179E-06), // birim_incld
    std::make_pair(5.164E-10, 2.733E-03), // qi_tend
    std::make_pair(1.370E+05, 6.582E+05), // ni_tend
  };

  IceSedData isds_fortran[] = {
    //       kts, kte, ktop, kbot, kdir,        dt,       odt, prt_sol, ranges
    IceSedData(1,  72,   27,   72,   -1, 1.800E+03, 5.556E-04,     0.0, ranges),
    IceSedData(1,  72,   72,   27,    1, 1.800E+03, 5.556E-04,     1.0, ranges),
    IceSedData(1,  72,   27,   27,   -1, 1.800E+03, 5.556E-04,     0.0, ranges),
    IceSedData(1,  72,   27,   27,    1, 1.800E+03, 5.556E-04,     2.0, ranges),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(IceSedData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  IceSedData isds_cxx[num_runs] = {
    IceSedData(isds_fortran[0]),
    IceSedData(isds_fortran[1]),
    IceSedData(isds_fortran[2]),
    IceSedData(isds_fortran[3]),
  };

    // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    ice_sedimentation(isds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    IceSedData& d = isds_cxx[i];
    ice_sedimentation_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                        d.rho, d.inv_rho, d.rhofaci, d.icldm, d.inv_dzq,
                        d.dt, d.odt,
                        d.qitot, d.qitot_incld, d.nitot, d.qirim, d.qirim_incld, d.birim, d.birim_incld,
                        d.nitot_incld, &d.prt_sol, d.qi_tend, d.ni_tend);
  }

  for (Int i = 0; i < num_runs; ++i) {
    // Due to pack issues, we must restrict checks to the active k space
    Int start = std::min(isds_fortran[i].kbot, isds_fortran[i].ktop) - 1; // 0-based indx
    Int end   = std::max(isds_fortran[i].kbot, isds_fortran[i].ktop);     // 0-based indx
    for (Int k = start; k < end; ++k) {
      REQUIRE(isds_fortran[i].qitot[k]       == isds_cxx[i].qitot[k]);
      REQUIRE(isds_fortran[i].qitot_incld[k] == isds_cxx[i].qitot_incld[k]);
      REQUIRE(isds_fortran[i].nitot[k]       == isds_cxx[i].nitot[k]);
      REQUIRE(isds_fortran[i].nitot_incld[k] == isds_cxx[i].nitot_incld[k]);
      REQUIRE(isds_fortran[i].qirim[k]       == isds_cxx[i].qirim[k]);
      REQUIRE(isds_fortran[i].qirim_incld[k] == isds_cxx[i].qirim_incld[k]);
      REQUIRE(isds_fortran[i].birim[k]       == isds_cxx[i].birim[k]);
      REQUIRE(isds_fortran[i].birim_incld[k] == isds_cxx[i].birim_incld[k]);
      REQUIRE(isds_fortran[i].qi_tend[k]     == isds_cxx[i].qi_tend[k]);
      REQUIRE(isds_fortran[i].ni_tend[k]     == isds_cxx[i].ni_tend[k]);
    }
    REQUIRE(isds_fortran[i].prt_sol == isds_cxx[i].prt_sol);
  }
}

static void run_bfb_homogeneous_freezing()
{
  const std::array< std::pair<Real, Real>, HomogeneousFreezingData::NUM_ARRAYS > ranges = {
    std::make_pair(C::homogfrze - 10, C::homogfrze + 10), // t
    std::make_pair(0.000E+00, 1.000E+00), // exner
    std::make_pair(0.000E+00, 1.000E+00), // xlf
    std::make_pair(0.000E+00, C::QSMALL*2), // qc
    std::make_pair(0.000E+00, 1.000E+00), // nc
    std::make_pair(0.000E+00, C::QSMALL*2), // qr
    std::make_pair(0.000E+00, 1.000E+00), // nr
    std::make_pair(0.000E+00, 1.000E+00), // qitot
    std::make_pair(0.000E+00, 1.000E+00), // nitot
    std::make_pair(0.000E+00, 1.000E+00), // qirim
    std::make_pair(0.000E+00, 1.000E+00), // birim
    std::make_pair(0.000E+00, 1.000E+00), // th
  };

  HomogeneousFreezingData hfds_fortran[] = {
    //                    kts, kte, ktop, kbot, kdir, ranges
    HomogeneousFreezingData(1,  72,   27,   72,   -1, ranges),
    HomogeneousFreezingData(1,  72,   72,   27,    1, ranges),
    HomogeneousFreezingData(1,  72,   27,   27,   -1, ranges),
    HomogeneousFreezingData(1,  72,   27,   27,    1, ranges),
  };

  static constexpr Int num_runs = sizeof(hfds_fortran) / sizeof(HomogeneousFreezingData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  HomogeneousFreezingData hfds_cxx[num_runs] = {
    HomogeneousFreezingData(hfds_fortran[0]),
    HomogeneousFreezingData(hfds_fortran[1]),
    HomogeneousFreezingData(hfds_fortran[2]),
    HomogeneousFreezingData(hfds_fortran[3]),
  };

    // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    homogeneous_freezing(hfds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    HomogeneousFreezingData& d = hfds_cxx[i];
    homogeneous_freezing_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                           d.t, d.exner, d.xlf,
                           d.qc, d.nc, d.qr, d.nr, d.qitot, d.nitot, d.qirim, d.birim, d.th);
  }

  for (Int i = 0; i < num_runs; ++i) {
    // Due to pack issues, we must restrict checks to the active k space
    Int start = std::min(hfds_fortran[i].kbot, hfds_fortran[i].ktop) - 1; // 0-based indx
    Int end   = std::max(hfds_fortran[i].kbot, hfds_fortran[i].ktop);     // 0-based indx
    for (Int k = start; k < end; ++k) {
      REQUIRE(hfds_fortran[i].qc[k]    == hfds_cxx[i].qc[k]);
      REQUIRE(hfds_fortran[i].nc[k]    == hfds_cxx[i].nc[k]);
      REQUIRE(hfds_fortran[i].qr[k]    == hfds_cxx[i].qr[k]);
      REQUIRE(hfds_fortran[i].nr[k]    == hfds_cxx[i].nr[k]);
      REQUIRE(hfds_fortran[i].qitot[k] == hfds_cxx[i].qitot[k]);
      REQUIRE(hfds_fortran[i].nitot[k] == hfds_cxx[i].nitot[k]);
      REQUIRE(hfds_fortran[i].qirim[k] == hfds_cxx[i].qirim[k]);
      REQUIRE(hfds_fortran[i].birim[k] == hfds_cxx[i].birim[k]);
      REQUIRE(hfds_fortran[i].th[k]    == hfds_cxx[i].th[k]);
    }
  }
}

static void run_bfb()
{
  run_bfb_calc_bulk_rhime();
  run_bfb_ice_sed();
  run_bfb_homogeneous_freezing();
}

};

}
}
}

namespace {

TEST_CASE("p3_ice_sed", "[p3_functions]")
{
  using TCS = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceSed;

  TCS::run_phys();
  TCS::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
