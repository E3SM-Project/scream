#ifndef P3_FUNCTIONS_PREVENT_ICE_OVERDEPLETION_IMPL_HPP
#define P3_FUNCTIONS_PREVENT_ICE_OVERDEPLETION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 ice overdepletion prevention.
 * Clients should NOT #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::prevent_ice_overdepletion(
  const Spack& pres, const Spack& t, const Spack& qv, const Spack& xxls, const Scalar& odt,
  Spack& qidep, Spack& qisub)
{
  constexpr Scalar cp = C::CP;
  constexpr Scalar rv = C::RH2O;

  const auto dumqvi = qv_sat(t,pres,true);
  const auto qdep_satadj = (qv-dumqvi) /
    (1 + square(xxls) * dumqvi / (cp * rv * square(t))) * odt;
  qidep *= min(1, max(0,  qdep_satadj) / max(qidep, sp(1.e-20)));
  qisub *= min(1, max(0, -qdep_satadj) / max(qisub, sp(1.e-20)));
}

} // namespace p3
} // namespace scream

#endif
