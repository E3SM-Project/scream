#ifndef P3_CALCULATE_MASS_IN_PACK_HPP
#define P3_CALCULATE_MASS_IN_PACK_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
S Functions<S,D>
::calculate_mass_of_pack(
    const Spack& qv, const Spack& qc, const Spack& qr, const Spack&qi, const Spack& rho, const Smask& context
  )
{
  Spack tw;
  tw.set(context, (qv+qc+qr+qi)*rho);
  return ekat::reduce_sum(tw);
}

} // namespace p3
} // namespace scream

#endif // P3_CALCULATE_MASS_IN_PACK_HPP

