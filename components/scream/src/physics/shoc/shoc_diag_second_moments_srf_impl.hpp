
#pragma once
// Code placed here is included only once per translation unit

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_diag_second_moments_srf(
    const MemberType& team, const Int& shcol, 
    const uview_1d<const Spack>& wthl_sfc, const uview_1d<const Spack>& uw_sfc, const uview_1d<const Spack>& vw_sfc,
    const view_1d<Spack>& ustar2, const view_1d<Spack>& wstar)
{
 // Purpose of this subroutine is to diagnose surface
 // properties needed for the the lower
 // boundary condition for the second order moments needed 
 // for the SHOC parameterization.  

 const auto one      = C::ONE;
 const auto zero     = C::ZERO;
 const auto third    = C::THIRD;
 const auto ggr      = C::gravit;
 const auto basetemp = C::basetemp;

 Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, shcol), [&] (Int k) {

    ustar2(k) = pack::sqrt(uw_sfc(k)*uw_sfc(k)+vw_sfc(k)*vw_sfc(k));

    const auto is_wthl_ge_zero = wthl_sfc(k) >= zero;

    wstar(k).set(is_wthl_ge_zero,
              pack::pow(one/basetemp*ggr*wthl_sfc(k)*one, third));

    wstar(k).set(!is_wthl_ge_zero, zero);
   }
 );

}

} // namespace shoc
} // namespace scream

