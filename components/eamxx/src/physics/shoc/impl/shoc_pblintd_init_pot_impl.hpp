
#ifndef SHOC_PBLINTD_INIT_POT_IMPL_HPP
#define SHOC_PBLINTD_INIT_POT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
template<typename TempViewType>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_pblintd_init_pot(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& thl,
    const uview_1d<const Spack>& ql,
    const TempViewType&          q,
    const uview_1d<Spack>&       thv)
{
   // Compute virtual potential temperature
   const auto lcond = C::LatVap;
   const auto cp    = C::Cpair;
   const auto eps   = C::ZVIR;
   const auto one   = C::ONE;

   const Int nlev_pack = ekat::npack<Spack>(nlev);

   Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const int& k) {
     const auto th = thl(k) + (lcond/cp)*ql(k);
     thv(k) = th * (one + eps*q(k) - ql(k));
   });
}

} // namespace shoc
} // namespace scream

#endif
