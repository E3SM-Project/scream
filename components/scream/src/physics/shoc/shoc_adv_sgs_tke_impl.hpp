#ifndef SHOC_ADV_SGS_TKE_IMPL_HPP
#define SHOC_ADV_SGS_TKE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc adv_sgs_tke. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::adv_sgs_tke(const Int& nlev, const Int& shcol, const Spack& dtime, const uview_1d<const Spack>& shoc_mix, const uview_1d<const Spack>& wthv_sec, const uview_1d<const Spack>& sterm_zt, const uview_1d<const Spack>& tk, const uview_1d<Spack>& tke, const uview_1d<Spack>& a_diss)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
