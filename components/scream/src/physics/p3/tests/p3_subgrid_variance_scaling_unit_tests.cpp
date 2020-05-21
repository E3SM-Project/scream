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
  
  KOKKOS_FUNCTION static void subgrid_variance_scaling_bfb_tests(const Scalar& relvar,
    const Scalar& expon, int& errors){
    //test that C++ and F90 implementations are BFB
    
    //Get value from F90 code
    Scalar f_scaling = subgrid_variance_scaling(relvar,expon);
    
    //Get value from C++ code
    const Spack relvars(relvar);
    Spack c_scaling = Functions::subgrid_variance_scaling(relvars,expon);
    
    if ( f_scaling != c_scaling[0] ){ //c_scaling is identical in each index, so [0] is sufficient
      printf("F90 and C++ subgrid_variance_scaling not BFB for relvar,expon,F90 val, C++val=",relvar,expon,f_scaling,c_scaling[0]);
      errors++;}
    
  } //end of bfb_tests
  
  KOKKOS_FUNCTION static void subgrid_variance_scaling_linearity_test(const Scalar& relvar,
    int& errors){
    //If expon=1, subgrid_variance_scaling should be 1 

    Scalar tol = (util::is_single_precision<Scalar>::value ) ? C::Tol*2 : C::Tol;
    
    //Get value from C++ code
    const Spack relvars(relvar);
    Spack c_scaling = Functions::subgrid_variance_scaling(relvars,1.0);
    
    if ( std::abs(c_scaling[0] -  1) > tol ){
      printf("subgrid_variance_scaling !=1 for exponent=1. relvar,scaling=",relvar,c_scaling)
	error++;}
  }
  
  KOKKOS_FUNCTION static void subgrid_variance_scaling_relvar1_test(int& errors){
    //If relvar=1, subgrid_variance_scaling should be factorial(expon)  
    
    Scalar tol = (util::is_single_precision<Scalar>::value ) ? C::Tol*2 : C::Tol;
    
    //Get value from C++ code
    const Spack relvars(relvar);
    Spack c_scaling = Functions::subgrid_variance_scaling(1.0,4.0);
    
    fact = std::factorial(4);
    
    if ( std::abs(c_scaling[0] -  fact) > tol ){
      printf("subgrid_variance_scaling not factorial(expon) for relvar=1. For expon=4, should be,is=",fact,c_scaling);
      error++;}
  }
  
  
  KOKKOS_FUNCTION static void subgrid_variance_scaling_relvar3_test(int& errors);
  //If expon=3, subgrid variance scaling should be relvar^3+3*relvar^2+2*relvar/relvar^3

  Scalar tol = (util::is_single_precision<Scalar>::value ) ? C::Tol*2 : C::Tol;

  static constexpr Int max_pack_size = 16;
  REQUIRE(Spack::n <= max_pack_size);
  Spack relvars = {
    0.1,0.5,1.0,2.0,
    3.0,4.0,5.0,6.0,
    6.5,7.0,8.0,9.0,
    9.1,9.5,9.8,10.}

  //Get value from C++ code
  Spack c_scaling = Functions::subgrid_variance_scaling(relvars,3.0);

  for (Int s = 0; s < Spack::n; ++s) {
    targ=1+3/relvars[s] + 2/relvars[s]**2;
    if ( std::abs(targ - c_scaling[s])>tol){
      printf("When expon=3, subgrid_variance_scaling doesn't match analytic expectation");
      errors++;
    }
  }
  
  static void run_subgrid_variance_scaling(){
    /*This function executes all the SGS variance scaling tests by looping
     *over a bunch of test and summing their return statuses.
     *If that sum is zero, no errors have occurred. Otherwise you have errors.
     *We do that loop in the parallel reduce below.
     */
    int nerr = 0;
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("SGSvarScaling::run", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {
	errors = 0;
	
	//First test that C++ and F90 implementations are BFB
	subgrid_variance_scaling_bfb_tests(1.0,1.0,errors);
	subgrid_variance_scaling_bfb_tests(1.0,2.0,errors);
	subgrid_variance_scaling_bfb_tests(1.0,-2.0,errors);
	subgrid_variance_scaling_bfb_tests(2.5,1.0,errors);
	subgrid_variance_scaling_bfb_tests(0.1,2.47,errors);
	subgrid_variance_scaling_bfb_tests(10.0,-1.1,errors);
	subgrid_variance_scaling_bfb_tests(10.0,0.1,errors);
	
	//If expon=1, subgrid_variance_scaling should be 1 
	subgrid_variance_scaling_linearity_test(10.,errors); 
	subgrid_variance_scaling_linearity_test(0.1,errors); 
	
	//If relvar=1, subgrid_variance_scaling should be factorial(expon)  
	subgrid_variance_scaling_relvar1_test(errors);
	
	//If expon=3, subgrid variance scaling should be relvar^3+3*relvar^2+2*relvar/relvar^3
	subgrid_variance_scaling_relvar3_test(2.0,errors);
	subgrid_variance_scaling_relvar3_test(3.5,errors);
	
      }, nerr);
      
    Kokkos::fence();
    REQUIRE(nerr == 0);
  }
}; //end of TestP3SubgridVarianceScaling struct
  
