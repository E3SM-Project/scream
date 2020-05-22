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
struct UnitWrap::UnitTest<D>::TestP3SubgridVarianceScaling
{

  //-----------------------------------------------------------------
  static void run_bfb_tests(){
    //test that C++ and F90 implementations are BFB

    //First create a Spack spanning the gammut of possible relvar vals
    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);
    Real relvars_info[max_pack_size] = {
      0.1,0.5,1.0,2.0,
      3.0,4.0,5.0,6.0,
      6.5,7.0,8.0,9.0,
      9.1,9.5,9.8,10.};
    
    // Sync relvars to device
    view_1d<Spack> relvars_device("relvars",1);
    auto relvars_host = Kokkos::create_mirror_view(relvars_device);    
    for (Int s = 0; s < Spack::n; ++s) {
      relvars_host(0)[s] = relvars_info[s];}

    Kokkos::deep_copy(relvars_device, relvars_host);
    
    //Get view of output from C++ fn call. 
    view_1d<Spack> scaling_device("c_scaling", 1);
    auto scaling_host = Kokkos::create_mirror_view(scaling_device);
    //note no copy here - just need to get output *back* from device.

    Spack f_scaling; //function output from F90 on CPU
    //input data for F90 needs to be a struct defined in p3_functions_f90.hpp
    SubgridVarianceScalingData f_data;
    
    //Loop over exponent values and make sure BFB is maintained in each case.
    Real expons[3] = {1.0,2.47,0.1};
    for (Int i = 0; i < 3; ++i) {

      // Call Fortran version
      for (Int j = 0; j < Spack::n; ++j) {
	f_data.relvar=relvars_host(0)[j];
	f_data.expon =expons[i];
        f_scaling[j] = subgrid_variance_scaling(f_data);
      }

      // Get expon value onto device
      view_1d<Scalar> expon_device("expon_device", 1);
      auto expon_host = Kokkos::create_mirror_view(expon_device);
      expon_host(0)=expons[i]; //this line causes build errors
      Kokkos::deep_copy(expon_device, expon_host);
      
      // Run the lookup from a kernel and copy results back to host
      Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& j) {
	  
	  // Call the function on C++
	  scaling_device(0) = Functions::subgrid_variance_scaling(
		              relvars_device(0),expon_device(0) );

	  //For debugging
	  for (Int s = 0; s < Spack::n; ++s) {	  
	    printf("relvar=%f, expon=%f, C++=%e\n",relvars_device(0)[s],
		   expon_device(0),scaling_device(0)[s]);
	  }
	  
	}); //end of parallel for
      
      // Copy results back to host
      Kokkos::deep_copy(scaling_host, scaling_device);
      
      // Validate results
      //Scalar c_scaling;
      for (Int s = 0; s < Spack::n; ++s) {

	//For debugging
	//printf("relvar=%f, expon=%f, C++=%e, F90=%e, diff=%e\n",relvars_host(0)[s],
	//       expon_host(0),scaling_host(0)[s],f_scaling[s],scaling_host(0)[s] - f_scaling[s]);
	
	REQUIRE(f_scaling[s] == scaling_host(0)[s] );
      }
    } //end loop over expons[i]

  } //end function run_bfb_tests

  //-----------------------------------------------------------------
  KOKKOS_FUNCTION static void subgrid_variance_scaling_linearity_test(const Scalar& relvar,
    int& errors){
    //If expon=1, subgrid_variance_scaling should be 1 

    Scalar tol = (util::is_single_precision<Scalar>::value ) ? C::Tol*2 : C::Tol;
    
    //Get value from C++ code
    const Spack relvars(relvar);
    Spack c_scaling = Functions::subgrid_variance_scaling(relvars,1.0);
    
    if ( std::abs(c_scaling[0] -  1) > tol ){
      printf("subgrid_variance_scaling should be 1 for expon=1, but is %e\n",c_scaling[0]);
	errors++;}
  }
  
  //-----------------------------------------------------------------
  KOKKOS_FUNCTION static void subgrid_variance_scaling_relvar1_test(int& errors){
    //If relvar=1, subgrid_variance_scaling should be factorial(expon)  
    
    Scalar tol = (util::is_single_precision<Scalar>::value ) ? C::Tol*2 : C::Tol;
    
    //Get value from C++ code
    const Spack ones(1);
    Spack c_scaling = Functions::subgrid_variance_scaling(ones,4.0);
    
    Real fact = std::tgamma(5.0); //factorial(n) = gamma(n+1) 
    
    if ( std::abs(c_scaling[0] -  fact) > tol ){ 
      printf("subgrid_variance_scaling should be factorial(expon) when relvar=1.\n");
      printf("For expon=4, should be %f but is=%f\n",fact,c_scaling[0]);
      errors++;}
  }

  //-----------------------------------------------------------------
  KOKKOS_FUNCTION static void subgrid_variance_scaling_relvar3_test(int& errors){
  //If expon=3, subgrid variance scaling should be relvar^3+3*relvar^2+2*relvar/relvar^3

  Scalar tol = (util::is_single_precision<Scalar>::value ) ? C::Tol*2 : C::Tol;
  
  static constexpr Int max_pack_size = 16;
  //tested that pack size is at least this big outside this fn because exception handling
  //not allowed on device.
  Real relvar_info[max_pack_size] = {0.1,0.5,1.0,2.0,
				     3.0,4.0,5.0,6.0,
				     6.5,7.0,8.0,9.0,
				     9.1,9.5,9.8,10.};

  //workaround b/c can't assign directly to Spack
  Spack relvars;
  for (Int s = 0; s < Spack::n; ++s) {
    relvars[s] = relvar_info[s];
  }
  
  //Get value from C++ code
  Spack c_scaling = Functions::subgrid_variance_scaling(relvars,3.0);

  Spack targ=1+3/relvars + 2/pack::pow(relvars,2.0);

  //Expected discrepancy is condition # * tolerance
  //For expon=3, expected val is 1+3/relvar + 2/relvar**2.
  //Condition number is x*f'(x)/f(x) = (3*relvar + 4)/(relvar**2. + 3*relvar+2)
  const Spack cond_num = (3.*relvars + 4.)/(pack::pow(relvars,2.0) +3*relvars+2.0);
  Spack max_tol = 1e2*tol*cond_num; //small # which is empirically big enough to pass.
  
  for (Int s = 0; s < Spack::n; ++s) {
    if ( std::abs(targ[s] - c_scaling[s])>max_tol[s]){
      printf("When expon=3, subgrid_variance_scaling doesn't match analytic expectation\n");
      printf("val = %e, expected = %e, diff = %e, tol = %e, cond_num = %e\n",c_scaling[s],c_scaling[s],
	     std::abs(targ[s] - c_scaling[s]),tol,cond_num[s]);
      errors++;
    } // end if
  }   //end for
  }   //end relvar3_test
  
  //-----------------------------------------------------------------
  static void run_property_tests(){
    /*This function executes all the SGS variance scaling tests by looping
     *over a bunch of test and summing their return statuses.
     *If that sum is zero, no errors have occurred. Otherwise you have errors.
     *We do that loop in the parallel reduce below.
     */
    int nerr = 0;

    //functions below use Spack size <16 but can't deal w/ exceptions on GPU, so do it here.
    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("SGSvarScaling::run", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {
	errors = 0;
		
	//If expon=1, subgrid_variance_scaling should be 1
	//                            args = relvar,return error count
	subgrid_variance_scaling_linearity_test(10.,errors); 
	subgrid_variance_scaling_linearity_test(0.1,errors); 
	
	//If relvar=1, subgrid_variance_scaling should be factorial(expon)
	//                            args = return error count
	subgrid_variance_scaling_relvar1_test(errors);
	
	//If expon=3, subgrid variance scaling should be relvar^3+3*relvar^2+2*relvar/relvar^3
	//                          args = return error count
	subgrid_variance_scaling_relvar3_test(errors);
	
      }, nerr);
      
    Kokkos::fence();
    REQUIRE(nerr == 0);
  } //end of TestP3SubgridVarianceScaling struct
  
}; // UnitWrap
  
} // namespace unit_test
} // namespace p3
} // namespace scream
 
namespace{

TEST_CASE("p3_subgrid_variance_scaling_test", "[p3_subgrid_variance_scaling_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3SubgridVarianceScaling::run_bfb_tests();
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3SubgridVarianceScaling::run_property_tests();
}

} // namespace

