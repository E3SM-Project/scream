#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::compute_shoc_temperature(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::uview_1d<const Func::Spack>&  thetal,
  const Func::uview_1d<const Func::Spack>&  ql,
  const Func::uview_1d<const Func::Spack>&  inv_exner,
  const Func::scratch_view_1d<Func::Spack>& tabs);

template void Func::compute_shoc_temperature(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::uview_1d<const Func::Spack>&  thetal,
  const Func::uview_1d<const Func::Spack>&  ql,
  const Func::uview_1d<const Func::Spack>&  inv_exner,
  const Func::uview_1d<Func::Spack>& tabs);

} // namespace shoc
} // namespace scream

#include "shoc_compute_shoc_temperature_impl.hpp"
