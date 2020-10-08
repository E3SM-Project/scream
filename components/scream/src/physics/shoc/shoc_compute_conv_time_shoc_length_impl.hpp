#ifndef SHOC_COMPUTE_CONV_TIME_SHOC_LENGTH_IMPL_HPP
#define SHOC_COMPUTE_CONV_TIME_SHOC_LENGTH_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc compute_conv_time_shoc_length. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_conv_time_shoc_length(
  const Scalar& pblh,
  Scalar&       conv_vel,
  Scalar&       tscale)
{
  // Take the cubed root and avoid negative values
  conv_vel = std::pow((conv_vel>0 ? conv_vel : 0), C::THIRD);

  // Compute eddy turnover timescale. If convective
  // velocity scale is zero then set to a minimum
  // threshold (100).
  if (conv_vel > 0) {
      tscale = pblh/conv_vel;
  } else {
    tscale = 100;
  }
}

} // namespace shoc
} // namespace scream

#endif
