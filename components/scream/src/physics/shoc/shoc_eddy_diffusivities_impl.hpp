#ifndef SHOC_EDDY_DIFFUSIVITIES_IMPL_HPP
#define SHOC_EDDY_DIFFUSIVITIES_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc eddy_diffusivities. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::eddy_diffusivities(const Int& nlev, const Int& shcol, const uview_1d<const Spack>& obklen, const uview_1d<const Spack>& pblh, const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& shoc_mix, const uview_1d<const Spack>& sterm_zt, const uview_1d<const Spack>& isotropy, const uview_1d<const Spack>& tke, const uview_1d<Spack>& tkh, const uview_1d<Spack>& tk)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
