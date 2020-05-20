#ifndef P3_FUNCTIONS_SUBGRID_VARIANCE_SCALING_IMPL_HPP
#define P3_FUNCTIONS_SUBGRID_VARIANCE_SCALING_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::subgrid_variance_scaling(const Spack& relvar, const Scalar& expon)
{
  /* We assume subgrid variations in qc follow a gamma distribution with inverse 
     relative variance relvar = 1/(var(qc)/qc**2). In this case, if the tendency
     for a given process is of the form A*qc**expon for a local value of qc, then
     the average process rate over the PDF is 
     gamma(relvar+expon)/[gamma(relvar)*relvar**expon]*A*average(qc)**expon. This 
     function calculates the local process rate => cell-average process rate scaling
     factor gamma(relvar+expon)/[gamma(relvar)*relvar**expon]. See Morrison and 
     Gettelman (2008; JCLI) eq 9 for details.
  */

  Spack result;
  Spack exponent=Spack(expon);
  result=pack::tgamma( relvar+exponent )/( pack::tgamma(relvar)*pack::pow(relvar,exponent) );
  
  return result;
				
}

} // namespace p3
} // namespace scream

#endif
