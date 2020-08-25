#ifndef SHOC_INTEG_COLUMN_STABILITY_IMPL_HPP
#define SHOC_INTEG_COLUMN_STABILITY_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::integ_column_stability(
  const MemberType& team,
  const Int& nlev,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& brunt,
  Scalar& brunt_int)
{
  const Int nlev_pack = scream::pack::npack<Spack>(nlev);
  brunt_int = 0;
  static constexpr auto troppres = 80000;

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {
      //auto range_pack1 = scream::pack::range<IntSmallPack>(k*Spack::n);

      //Find if pressure if greater that troposheric pressure
      auto press_gt_troppress = (pres(k) > troppres);
      if(press_gt_troppress.any()){
	//brunt_int = brunt_int + dz_zt(k) * brunt(k);
      }
    });


  /* brunt_int(1:shcol) = 0._rtype
     do k = 1, nlev
     do i = 1, shcol
        if (pres(i,k) .gt. troppres) then
	   brunt_int(i) = brunt_int(i) + dz_zt(i,k)*brunt(i,k)
        endif
     enddo
     enddo*/

}

} // namespace shoc
} // namespace scream

#endif
