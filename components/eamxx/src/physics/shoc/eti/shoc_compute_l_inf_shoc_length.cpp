#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::compute_l_inf_shoc_length(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::scratch_view_1d<Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>&  tke,
  Func::Scalar&                             l_inf);

} // namespace shoc
} // namespace scream

#include "shoc_compute_l_inf_shoc_length_impl.hpp"
