#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::compute_brunt_shoc_length(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::scratch_view_1d<Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>&  thv,
  const Func::uview_1d<const Func::Spack>&  thv_zi,
  const Func::uview_1d<Func::Spack>&        brunt);

template void Func::compute_brunt_shoc_length(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::uview_1d<Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>&  thv,
  const Func::uview_1d<const Func::Spack>&  thv_zi,
  const Func::uview_1d<Func::Spack>&        brunt);

template void Func::compute_brunt_shoc_length(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::uview_1d<const Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>&  thv,
  const Func::uview_1d<const Func::Spack>&  thv_zi,
  const Func::uview_1d<Func::Spack>&        brunt);

} // namespace shoc
} // namespace scream

#include "shoc_compute_brunt_shoc_length_impl.hpp"
