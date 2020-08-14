
#pragma once
// Code placed below is included only once per translation unit

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_diag_second_moments_ubycond(
    const MemberType& team, const Int& shcol, const Int& num_tracer,
    const view_1d<Spack>& thl_sec, const view_1d<Spack>& qw_sec, const view_1d<Spack>& wthl_sec, const view_1d<Spack>& wqw_sec,
    const view_1d<Spack>& qwthl_sec, const view_1d<Spack>& uw_sec, const view_1d<Spack>& vw_sec,const view_1d<Spack>& wtke_sec, 
    const view_2d<Spack>& wtracer_sec)
{
  // Purpose of this subroutine is to diagnose the upper
  //  boundary condition for the second order moments
  //  needed for the SHOC parameterization.  Currently
  //  set all to zero.

 const auto zero = C::ZERO;

 const Int nshcol_pack = scream::pack::npack<Spack>(shcol);
 const Int num_pack = scream::pack::npack<Spack>(num_tracer);

 Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nshcol_pack), [&] (Int& k) {

      wthl_sec(k)  = zero;
      wqw_sec(k)   = zero;
      uw_sec(k)    = zero;
      vw_sec(k)    = zero;
      wtke_sec(k)  = zero;
      thl_sec(k)   = zero;
      qw_sec(k)    = zero;
      qwthl_sec(k) = zero;

    }
  );

 Kokkos::parallel_for(Kokkos::TeamThreadRange(team, shcol), [&] (int& k) {
   Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, num_pack),[&] (int& p) {

      wtracer_sec(k,p) = zero;

    });
   });

}

} // namespace shoc
} // namespace scream

