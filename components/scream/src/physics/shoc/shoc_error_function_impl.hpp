#ifndef SHOC_ERROR_FUNCTION_IMPL_HPP
#define SHOC_ERROR_FUNCTION_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc error_function. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::error_function(const Scalar &input,
                                    Scalar       &output)
{
  output = std::erf(input);
}

} // namespace shoc
} // namespace scream

#endif
