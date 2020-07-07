#include "catch2/catch.hpp"

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/common/physics_functions.hpp"
#include "physics/common/physics_saturation_impl.hpp"

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
struct UnitWrap::UnitTest<D>::TestP3Saturation
{

  static Scalar condNum(const Scalar& svp, const Scalar& temp, const bool isIce=false){

    //compute condition number for saturation vapor pressure calc.

    //Set constant values
    //--------------------------------------
    const Scalar RV = C::RV;
    const Scalar RhoLiq = C::RhoH2O;
    const Scalar RhoIce = C::RhoIce;
    const Scalar LatVap = C::LatVap;
    const Scalar LatIce = C::LatIce;

    //==========================================================
    // Test Saturation Vapor Pressure
    //======================================================5A5A5A====
    // First calculate condition number Cond=x*f'(x)/f(x), which gives the growth in roundoff error
    // expected from the function. Here, x=temperature, f(x) = sat_X_p, and
    // f'(x) = L/T*sat_X_p/(Rv*T - sat_X_p/rho_{liq or ice}) from Clausius-Clapeyron.
    // Note this analytical solution isn't exact because it misses things like surface tension
    // so it can't replace curve fits like Flatau. But Clausius-Clapeyron is good enough for
    // getting a ballpark value, which is all we're doing here.

    auto latent = LatVap; //for liquid
    auto rho    = RhoLiq; // for liquid
    if(isIce){
      latent = LatVap+LatIce; //for ice
      rho    = RhoIce;        // for ice
    }

    return latent*svp/(RV*temp - svp/rho); //return condition number

    }

  KOKKOS_FUNCTION  static void saturation_tests(const Scalar& temperature, const Scalar& pressure,
						//Flatau (polysvp1) correct values
						const Scalar& correct_sat_ice_p, const Scalar& correct_sat_liq_p,
						const Scalar& correct_mix_ice_rp, const Scalar& correct_mix_liq_rp,
						//Muphy and Koop (MurphyKoop_svp) correct values
						const Scalar& correct_sat_ice_mk, const Scalar& correct_sat_liq_mk,
						const Scalar& correct_mix_ice_rmk, const Scalar& correct_mix_liq_rmk,
						int& errors ){

    //Allow usage of saturation functions
    using physics = scream::physics::Functions<Scalar, Device>;

    //Convert Scalar inputs to Spacks because that's what polysvp1 and qv_sat expect as inputs.
    //--------------------------------------
    const Spack temps(temperature);
    const Spack pres(pressure);

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    //Get values from polysvp1 and qv_sat (calling MurphyKoop_svp) to test against "correct" values
    //--------------------------------------
    Spack sat_ice_p  = physics::polysvp1(temps, true);
    Spack sat_liq_p  = physics::polysvp1(temps, false);
    Spack mix_ice_rp = physics::qv_sat(temps, pres, true, 0);
    Spack mix_liq_rp = physics::qv_sat(temps, pres, false,0);

    //Get values from MurphyKoop_svp and qv_sat (calling MurphyKoop_svp) to test against "correct" values
    Spack sat_ice_mk   = physics::MurphyKoop_svp(temps, true);
    Spack sat_liq_mk   = physics::MurphyKoop_svp(temps, false);
    Spack mix_ice_rmk  = physics::qv_sat(temps, pres, true, 1);
    Spack mix_liq_rmk  = physics::qv_sat(temps, pres, false,1);

    //Set error tolerances
    //--------------------------------------
    //: C::Tol is machine epsilon for single or double precision as appropriate. This will be
    //multiplied by a condition # to get the actual expected numerical uncertainty. If single precision, multiplying
    //by an additional fudge factor because the python values we're comparing against were computed via a
    //bunch of double-precision intermediate calculations which will add error.

    Scalar tol = (util::is_single_precision<Scalar>::value ) ? C::Tol*2 : C::Tol;

    //PMC note: original version looped over pack dimension, testing each entry. This isn't
    //necessary b/c packs were created by copying a scalar up to pack size. Thus just evaluating
    // 1st entry below.

    const Scalar Cond_ice_p=condNum(sat_ice_p[0],temps[0], true);
    const Scalar Cond_liq_p=condNum(sat_liq_p[0],temps[0]);

    const Scalar Cond_ice_mk=condNum(sat_ice_mk[0],temps[0], true);
    const Scalar Cond_liq_mk=condNum(sat_liq_mk[0],temps[0]);

    // Test vapor pressure against Flatau's impl of Wexler:
    // ---------------------------------------------------------
    // Now check that computed vs expected values are small enough.
    if ( std::abs(sat_ice_p[0] - correct_sat_ice_p ) > Cond_ice_p*tol ) {
      printf("esi_p for T = %f abs diff is %e but max allowed is %e\n",temperature,std::abs(sat_ice_p[0] - correct_sat_ice_p ),tol*Cond_ice_p );
      std::cout<<"sat_ice_p:"<<sat_ice_p[0]<<", "<<correct_sat_ice_p<<","<<temps[0]<<"\n";
      errors++;}
    if (std::abs(sat_liq_p[0] - correct_sat_liq_p) > Cond_liq_p*tol)  {
      printf("esl_p  for T = %f abs diff is %e but max allowed is %e\n",temperature,std::abs(sat_liq_p[0] - correct_sat_liq_p ),tol*Cond_liq_p);
      std::cout<<"sat_liq_p:"<<sat_liq_p[0]<<", "<<correct_sat_liq_p<<","<<temps[0]<<"\n";
      errors++;}

    // Test vapor pressure against Murphy and Koop:
    // ---------------------------------------------------------
    // Now check that computed vs expected values are small enough.
    if ( std::abs(sat_ice_mk[0] - correct_sat_ice_mk ) > Cond_ice_mk*tol ) {
      printf("esi_mk for T = %f abs diff is %e but max allowed is %e\n",temperature,std::abs(sat_ice_mk[0] - correct_sat_ice_mk ),tol*Cond_ice_mk );
      errors++;}
    if (std::abs(sat_liq_mk[0] - correct_sat_liq_mk) > Cond_liq_mk*tol)  {
      printf("esl_mk  for T = %f abs diff is %e but max allowed is %e\n",temperature,std::abs(sat_liq_mk[0] - correct_sat_liq_mk ),tol*Cond_liq_mk);
      std::cout<<"sat_liq_mk:"<<sat_liq_mk[0]<<", "<<correct_sat_liq_mk<<","<<temps[0]<<"\n";
      errors++;}

    //==========================================================
    // Test Saturation Mixing Ratio
    //==========================================================
    // First, compute condition # Cond=x*f'(x)/f(x), with x=temperature, f(x)=mix_X_r[0], and
    // f'(x) = L*mix_X_r[0]/(Rv*temperature**2.) from Clausius-Clapeyron. Nice cancelation leaves:

    //Set constant values
    //--------------------------------------
    const Scalar RV = C::RV;
    const Scalar LatVap = C::LatVap;
    const Scalar LatIce = C::LatIce;

    const Scalar Cond_ice_r=(LatVap+LatIce)/(RV*temps[0]);
    const Scalar Cond_liq_r=LatVap/(RV*temps[0]);

    //Test mixing-ratios against Wexler approx:
    // -------------------
    // Now check that computed vs expected values are small enough (Flatau).
    if (std::abs(mix_ice_rp[0] -  correct_mix_ice_rp) > Cond_ice_r*tol ) {
      printf("qsi_p: abs(calc-expected)=%e %e\n",std::abs(mix_ice_rp[0] -  correct_mix_ice_rp),tol*Cond_ice_r);
      errors++;}
    if (std::abs(mix_liq_rp[0] -  correct_mix_liq_rp) > Cond_liq_r*tol ) {
      printf("qsl_p: abs(calc-expected)=%e %e\n",std::abs(mix_liq_rp[0] -  correct_mix_liq_rp),tol*Cond_liq_r);
      std::cout<<"qsl_p:"<<mix_liq_rp[0]<<", "<<correct_mix_liq_rp<<","<<temps[0]<<"\n";
      errors++;}

    // Now check that computed vs expected values are small enough (Murphy and Koop).
    if (std::abs(mix_ice_rmk[0] -  correct_mix_ice_rmk) > Cond_ice_r*tol ) {
      printf("qsi_mk: abs(calc-expected)=%e %e\n",std::abs(mix_ice_rmk[0] -  correct_mix_ice_rmk),tol*Cond_ice_r);
      errors++;}
    if (std::abs(mix_liq_rmk[0] -  correct_mix_liq_rmk) > Cond_liq_r*tol ) {
      printf("qsl_mk: abs(calc-expected)=%e %e\n",std::abs(mix_liq_rmk[0] -  correct_mix_liq_rmk),tol*Cond_liq_r);
      errors++;}

  }

  static void run()
  {
    /*Originally written by Kyle Pressel, updated by Peter Caldwell on 4/5/20.
     *This code tests polysvp1 and qv_sat at 0 degrees C, at a very cold T, and at a very hot T
     *to make sure our impl gets the same answer as Flatau et al 1992:
     *(https://journals.ametsoc.org/doi/abs/10.1175/1520-0450%281992%29031%3C1507%3APFTSVP%3E2.0.CO%3B2
     *For 0 degrees C, polysvp values can be read directly from Flatau. For other cases, I independently
     *coded up the Flatau scheme in python and used it to derive the expected values. My python code is
     *in https://github.com/E3SM-Project/scream-docs.git analysis-scripts/test_qv_sat.py
     */


    int nerr = 0;
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("TestTableIce::run", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {

      errors = 0;
      const auto tmelt = C::Tmelt;

      // Just Freezing Case: Test values @ 273.15K @ 1e5 Pa
      //---------------------------------------
      // This is the nicest test b/c polysvp1 is a polynomial fit around (T-273.15)
      // so T=273.15 collapses back to the intercept coefficient which can be read
      // directly from the RHS of table 4 of Flatau et al 1992.
      // Note that ice values are identical to liquid values b/c C++ uses liq val for T>=0 C.

      saturation_tests(tmelt, 1e5, 611.23992100000009, 611.23992100000009,
		       0.0038251131382843278, 0.0038251131382843278,
		       611.2126978267946, 611.2126978267946,
		       0.0038249417291628678, 0.0038249417291628678, errors);

      // Cold Case: Test values @ 243.15K @ 1e5 Pa
      //---------------------------------------
      saturation_tests(243.15, 1e5, 38.024844602056795, 51.032583257624964,
		       0.00023659331311441935, 0.00031756972127516819,
		       38.01217277745647, 50.93561537896607,
		       0.00023651443812988484, 0.0003169659941390894, errors);

      //Warm Case: Test values @ 303.15 @ 1e5 Pa
      //---------------------------------------
      saturation_tests(303.15, 1e5, 4245.1933273786717, 4245.1933273786717,
		       0.027574442204332306, 0.027574442204332306,
		       4246.814076877233, 4246.814076877233,
		       0.027585436614272162, 0.027585436614272162, errors);

      //Following values are picked from Murphy and Koop (2005)
      //Table C1 titled: "VALUES RECOMMENDED FOR CHECKING COMPUTER CODES"

      saturation_tests(150, 1e5, 0.0565113360640801, 0.17827185346988017,
		       3.514840868181163e-07, 1.1088004687434155e-06,
		       6.1061006509816675e-06, 1.5621037177920032e-05,
		       3.79781500154674e-11, 9.715825652191167e-11, errors);

      saturation_tests(180, 1e5, 0.0565113360640801, 0.17827185346988017,
		       3.514840868181163e-07, 1.1088004687434155e-06,
		       0.005397500125274297, 0.01123923029036248,
		       3.3570864981545485e-08, 6.990471437859482e-08, errors);

      saturation_tests(210, 1e5, 0.7021857060894199, 1.2688182238880685,
		       4.3674192198021676e-06, 7.891776277281936e-06,
		       0.7020234713180218, 1.2335424085746476,
		       4.366410153074117e-06, 7.672365591576349e-06, errors);

      saturation_tests(240, 1e5, 27.280908658710246, 37.77676490183603,
		       0.00016972553017335565, 0.0002350491600776575,
		       27.272365420780556, 37.66700070557609,
		       0.00016967236474679822, 0.0002343659437158037, errors);

      saturation_tests(273.16, 1e5, 611.6840516537769, 611.6840516537769,
		       0.003827909594290528, 0.003827909594290528,
		       611.6570436443282, 611.6570436443282,
		       0.0038277395384149105, 0.0038277395384149105, errors);

      saturation_tests(300, 1e5, 3535.4066341569387, 3535.4066341569387,
		       0.022795088436007804, 0.022795088436007804,
		       3536.7644130514645, 3536.7644130514645,
		       0.022804163906259393, 0.022804163906259393, errors);

      //Change Pressure Case: Test values @ 243.15 @ 500 mb
      //---------------------------------------
      saturation_tests(243.15, 5e4, 38.024844602056795, 51.032583257624964,
		       0.00047336669164733106, 0.00063546390177500586,
		       38.01217277745647, 50.93561537896607,
		       0.00047320882161578, 0.0006342552147122389, errors);
    }, nerr);

    Kokkos::fence();
    REQUIRE(nerr == 0);
  }
}; //end of TestP3Saturation struct

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace{

TEST_CASE("p3_saturation_test", "[p3_saturation_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Saturation::run();

 } // TEST_CASE

} // namespace
