#include "catch2/catch.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/scream_types.hpp"

#include "p3_unit_tests_common.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestPreventLiqSupersaturation {

  static void run_bfb()
  {

    auto engine = setup_random_test();

    PreventLiqSupersaturationData f90_data[max_pack_size];

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
      d.dt = f90_data[0].dt; // Hold this fixed, this is not packed data

    }
    
    // Create copies of data for use by cxx and sync it to device. Needs to happen before 
    // fortran calls so that inout data is in original state
    view_1d<PreventLiqSupersaturationData> cxx_device("cxx_device", max_pack_size);
    const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
    std::copy(&f90_data[0], &f90_data[0] + max_pack_size, cxx_host.data());
    Kokkos::deep_copy(cxx_device, cxx_host);

    // Save copy of inout vars to check that prevent_liq_supersaturation always makes them smaller
    Spack qi2qv_sublim_tend_init,qr2qv_evap_tend_init;
    for (Int i = 0; i < max_pack_size; ++i) {
      qi2qv_sublim_tend_init[i] = cxx_host(i).qi2qv_sublim_tend;
      qr2qv_evap_tend_init[i] = cxx_host(i).qr2qv_evap_tend;
    }
      
    // Get data from fortran
    for (auto& d : f90_data) {
      prevent_liq_supersaturation(d);
    }
    
    // Get data from cxx. Run prevent_liq_supersaturation from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Scalar dt;
      Spack latent_heat_sublim, latent_heat_vapor, pres, qi2qv_sublim_tend, qidep, qinuc, qr2qv_evap_tend, qv, t_atm;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
	dt = cxx_device(vs).dt; //dt is scalar but PreventLiqSupersaturationData has diff val for each row.
        latent_heat_sublim[s] = cxx_device(vs).latent_heat_sublim;
        latent_heat_vapor[s] = cxx_device(vs).latent_heat_vapor;
        pres[s] = cxx_device(vs).pres;
        qi2qv_sublim_tend[s] = cxx_device(vs).qi2qv_sublim_tend;
        qidep[s] = cxx_device(vs).qidep;
        qinuc[s] = cxx_device(vs).qinuc;
        qr2qv_evap_tend[s] = cxx_device(vs).qr2qv_evap_tend;
        qv[s] = cxx_device(vs).qv;
        t_atm[s] = cxx_device(vs).t_atm;
      }

      Functions::prevent_liq_supersaturation(pres, t_atm, qv, latent_heat_vapor, latent_heat_sublim, dt, qidep, qinuc, qi2qv_sublim_tend, qr2qv_evap_tend);

      // Copy spacks back into cxx_device view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cxx_device(vs).qi2qv_sublim_tend = qi2qv_sublim_tend[s];
        cxx_device(vs).qr2qv_evap_tend = qr2qv_evap_tend[s];
      }

    });

    Kokkos::deep_copy(cxx_host, cxx_device);

    // Verify BFB results
    for (Int i = 0; i < max_pack_size; ++i) {
      PreventLiqSupersaturationData& d_f90 = f90_data[i];
      PreventLiqSupersaturationData& d_cxx = cxx_host[i];
      REQUIRE(d_f90.qi2qv_sublim_tend == d_cxx.qi2qv_sublim_tend);
      REQUIRE(d_f90.qr2qv_evap_tend == d_cxx.qr2qv_evap_tend);

      //Verify tendencies are always >=0:
      REQUIRE(d_cxx.qi2qv_sublim_tend>=0);
      REQUIRE(d_cxx.qr2qv_evap_tend>=0);

      //Verify function call always makes tendencies smaller
      REQUIRE(d_cxx.qi2qv_sublim_tend<=qi2qv_sublim_tend_init[i]);
      REQUIRE(d_cxx.qr2qv_evap_tend<=qr2qv_evap_tend_init[i]);

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("prevent_liq_supersaturation_bfb", "[p3]")
{
  using TestStruct = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPreventLiqSupersaturation;

  TestStruct::run_bfb();
}

} // empty namespace
