#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"

#include "p3_unit_tests_common.hpp"

//#include <thread>
//#include <array>
//#include <algorithm>
//#include <random>
//#include <iomanip>      // std::setprecision

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestEvapSublPrecip
{
  static void evaporate_rain_unit_bfb_tests(){
    
    //fortran generated data is input to the following
    //This subroutine has 12 args, only 10 are supplied here for invoking it as last 2 are intent-outs
    EvapRainData espd[max_pack_size] = {
    //qr_incld,     qc_incld,    nr_incld,    qi_incld,    cld_frac_l,  cld_frac_r,  qv,          qv_prev,     qv_sat_l,    qv_sat_i,    ab,          abi,         epsr,        epsi_tot,    t,           t_prev,      latent_heat_sublim, dqsdt,dt
      {3.447210e-03,7.568759e-03,1.282462e+07,4.559716e-03,2.181899e-01,4.996576e-01,1.355174e-03,7.689148e-03,3.109884e-03,8.479143e-03,1.323354e+00,1.801880e+00,7.912440e+02,4.717223e+02,2.788412e+02,2.304422e+02,3.376491e+05,3.185133e-03,6.303449e+02,},
      {7.713538e-03,1.421715e-03,1.500325e+07,3.424152e-04,7.246117e-01,9.761673e-02,3.144175e-03,3.480243e-03,6.210449e-03,5.644419e-03,1.020853e+00,1.123862e+00,9.815444e+02,5.866347e+02,1.005243e+02,2.664643e+02,3.382550e+05,3.459718e-03,3.346337e+02,},
      {1.108357e-03,7.919448e-03,7.522689e+05,4.003924e-03,4.391406e-01,5.353028e-01,9.714450e-03,8.163269e-03,7.232859e-03,3.540476e-04,1.613577e+00,1.365200e+00,9.686382e+02,9.874188e+01,3.080983e+02,2.462458e+02,3.331664e+05,6.260869e-03,1.129514e+02,},
      {8.763366e-03,5.775320e-03,4.081757e+07,4.872870e-03,9.093062e-01,8.637157e-01,5.111261e-06,3.510804e-03,7.919468e-03,2.164816e-03,1.111913e+00,1.457051e+00,4.440636e+02,3.745193e+02,1.495478e+02,2.390320e+02,3.272279e+05,3.685088e-03,1.564749e+03,},
      {1.772029e-03,9.813384e-04,6.163759e+05,6.511990e-03,5.961482e-01,2.418172e-01,6.505523e-03,1.566927e-03,2.162630e-03,8.318588e-03,1.742185e+00,1.773503e+00,5.219377e+01,1.346690e+02,1.225765e+02,1.770784e+02,3.640150e+05,3.694893e-03,1.739467e+03,},
      {2.120679e-03,8.847332e-03,1.831657e+07,8.027861e-03,3.950421e-03,4.655730e-01,1.258581e-03,1.156046e-03,4.378864e-03,2.513115e-03,1.042896e+00,1.513305e+00,6.661906e+02,9.245433e+02,1.651056e+02,1.242327e+02,3.620034e+05,1.663144e-03,5.727666e+02,},
      {6.818583e-03,4.402694e-03,9.635189e+07,1.215843e-03,7.442543e-01,9.732511e-01,9.959752e-03,6.339721e-03,6.630242e-03,3.272258e-03,1.929628e+00,1.219269e+00,3.573017e+02,6.716209e+02,1.173579e+02,1.438046e+02,3.220013e+05,3.378694e-03,1.248869e+03,},
      {4.682860e-04,1.241868e-05,5.672542e+07,5.023443e-03,3.341924e-01,5.611753e-01,1.540110e-03,6.376444e-03,1.069738e-03,6.207760e-03,1.531449e+00,1.745637e+00,6.867427e+02,5.414631e+02,3.453410e+02,2.398160e+02,3.701868e+05,1.279460e-03,6.330565e+02,},
      {9.565079e-03,5.340370e-03,5.688460e+07,3.615027e-03,8.522408e-01,3.928986e-01,8.402562e-03,6.102022e-03,1.424570e-03,6.383195e-03,1.524161e+00,1.578820e+00,6.580376e+02,8.870888e+02,2.226693e+02,2.003213e+02,3.757752e+05,5.097239e-03,3.384718e+02,},
      {3.250646e-03,2.675118e-03,1.815107e+07,7.853088e-03,9.806948e-01,7.486889e-01,3.255975e-03,8.011100e-03,3.417518e-03,6.616311e-04,1.456686e+00,1.901345e+00,3.438521e+02,3.756087e+02,1.249969e+02,3.129243e+02,3.792944e+05,1.431032e-03,1.056702e+03,},
      {2.368164e-03,6.522920e-03,9.013133e+06,3.347165e-03,5.839828e-01,9.744957e-01,3.529294e-03,4.658068e-03,5.006151e-03,9.455110e-03,1.152108e+00,1.033423e+00,2.548847e+02,6.425814e+02,3.434567e+02,1.980667e+02,3.794830e+05,5.189345e-03,1.224861e+03,},
      {1.195296e-03,6.032559e-03,8.642695e+06,9.974985e-03,6.349056e-01,8.434733e-01,1.116853e-04,5.921699e-03,9.604646e-03,6.112554e-03,1.994160e+00,1.517046e+00,5.212776e+01,3.082936e+02,2.951571e+02,1.349963e+02,3.409254e+05,7.670154e-03,2.924272e+02,},
      {1.679714e-03,9.615619e-03,3.118296e+07,8.024876e-03,3.137640e-01,9.696389e-01,2.687520e-03,3.426577e-03,9.123650e-03,6.113998e-03,1.505228e+00,1.080199e+00,2.957088e+01,1.899326e+02,1.774232e+02,1.056878e+02,3.379270e+05,2.822660e-03,4.617910e+02,},
      {2.503130e-03,3.617091e-03,9.913080e+06,6.901550e-03,5.781765e-01,1.999869e-01,6.985166e-03,1.190207e-03,7.204522e-03,1.357847e-03,1.459423e+00,1.821947e+00,4.741279e+02,6.710361e+02,2.176799e+02,2.688688e+02,3.622164e+05,5.888658e-03,1.427060e+03,},
      {3.627839e-03,5.863559e-03,8.680482e+07,5.928597e-03,2.893065e-01,6.883180e-01,4.446598e-03,8.259366e-03,7.609455e-03,3.984928e-03,1.015280e+00,1.838312e+00,2.369260e+02,6.717340e+02,3.355737e+02,1.851869e+02,3.638496e+05,6.984012e-03,1.017520e+03,},
      {1.222023e-04,5.561842e-03,3.750663e+07,8.378003e-03,5.483419e-01,7.604896e-01,8.352532e-03,4.889095e-03,6.985609e-03,1.998668e-03,1.552162e+00,1.105295e+00,5.962036e+02,3.634128e+02,2.942480e+02,2.131546e+02,3.507114e+05,7.215172e-03,1.589544e+03}
    };

    // Sync to device
    view_1d<EvapRainData> espd_device("espd", max_pack_size);
    auto espd_host = Kokkos::create_mirror_view(espd_device);

    // This copy only copies the input variables.
    std::copy(&espd[0], &espd[0] + max_pack_size, espd_host.data());
    Kokkos::deep_copy(espd_device, espd_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      evaporate_rain(espd[i]);
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack qr_incld,qc_incld,nr_incld,qi_incld,
	cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i,
	ab,abi,epsr,epsi_tot,t,t_prev,latent_heat_sublim,dqsdt;

      Scalar dt;
      
      // Init pack outputs
      Spack qr2qv_evap_tend, nr_evap_tend;

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qr_incld[s]    = espd_device(vs).qr_incld;
        qc_incld[s]    = espd_device(vs).qc_incld;
        nr_incld[s]    = espd_device(vs).nr_incld;
        qi_incld[s] = espd_device(vs).qi_incld;
        cld_frac_l[s]       = espd_device(vs).cld_frac_l;
        cld_frac_r[s]       = espd_device(vs).cld_frac_r;
	qv[s] = espd_device(vs).qv;
	qv_prev = espd_device(vs).qv_prev;
        qv_sat_l[s]         = espd_device(vs).qv_sat_l;
	qv_sat_i[s]    = espd_device(vs).qv_sat_i;
        ab[s]          = espd_device(vs).ab;
	abi[s]         = espd_device(vs).abi;
        epsr[s]        = espd_device(vs).epsr;
	epsi_tot[s]    = espd_device(vs).epsi_tot;
	t[s]           = espd_device(vs).t;
	t_prev[s]      = espd_device(vs).t_prev;
	latent_heat_sublim[s]=espd_device(vs).latent_heat_sublim;
        //qr2qv_evap_tend[s]       = espd_device(vs).qr2qv_evap_tend; //PMC shouldn't have to init output vars.
        //nr_evap_tend[s]       = espd_device(vs).nr_evap_tend;
      }

      Functions::evaporate_rain(qr_incld,qc_incld,nr_incld,qi_incld,
				cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i,
				ab,abi,epsr,epsi_tot,t,t_prev,latent_heat_sublim,dqsdt,dt,
				qr2qv_evap_tend,nr_evap_tend);
				
      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        espd_device(vs).qr_incld    = qr_incld[s];
        espd_device(vs).qc_incld    = qc_incld[s];
        espd_device(vs).nr_incld    = nr_incld[s];
        espd_device(vs).qi_incld = qi_incld[s];
        espd_device(vs).cld_frac_l       = cld_frac_l[s];
        espd_device(vs).cld_frac_r       = cld_frac_r[s];
        espd_device(vs).qv_sat_l         = qv_sat_l[s];
        espd_device(vs).ab          = ab[s];
        espd_device(vs).epsr        = epsr[s];
        espd_device(vs).qv          = qv[s];
        espd_device(vs).qr2qv_evap_tend       = qr2qv_evap_tend[s];
        espd_device(vs).nr_evap_tend       = nr_evap_tend[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(espd_host, espd_device);

    // Validate results
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(espd[s].qr2qv_evap_tend == espd_host(s).qr2qv_evap_tend);
      REQUIRE(espd[s].nr_evap_tend == espd_host(s).nr_evap_tend);
    }
  }

  static void run_bfb(){
    evaporate_rain_unit_bfb_tests();
  }

}; //TestEvapSublPrecip UnitWrap
  
}//namespace unit_test
}//namespace p3
}//namespace scream

namespace {

  TEST_CASE("p3_evaporate_rain_test", "[p3_unit_tests]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestEvapSublPrecip::run_bfb();
}

}// anonymous namespace
