#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::compute_shoc_vapor(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::uview_1d<const Func::Spack>&  qw,
  const Func::uview_1d<const Func::Spack>&  ql,
  const Func::scratch_view_1d<Func::Spack>& qv);

} // namespace shoc
} // namespace scream

#include "shoc_compute_shoc_vapor_impl.hpp"
