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
  const Int nlev_pack = ekat::pack::npack<Spack>(nlev);
  brunt_int = 0;
  static constexpr auto troppres = 80000;

  //Use parallel_reduce to compute reduction in brunt_int
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k, Scalar& val) {
      //Find if pressure is greater than tropospheric pressure
      auto press_gt_troppress = (pres(k) > troppres);

      Spack my_result(0);
      my_result.set(press_gt_troppress, dz_zt(k) * brunt(k));
      for (int s = 0; s < Spack::n; ++s) {
	val += my_result[s];
      }
    }, brunt_int);

}

} // namespace shoc
} // namespace scream

#endif
