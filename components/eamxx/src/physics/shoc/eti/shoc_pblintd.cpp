#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing pblintd on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::pblintd(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Int&                                npbl,
  const Func::uview_1d<const Func::Spack>&  z,
  const Func::uview_1d<const Func::Spack>&  zi,
  const Func::uview_1d<const Func::Spack>&  thl,
  const Func::uview_1d<const Func::Spack>&  ql,
  const Func::scratch_view_1d<Func::Spack>& q,
  const Func::uview_1d<const Func::Spack>&  u,
  const Func::uview_1d<const Func::Spack>&  v,
  const Func::Scalar&                       ustar,
  const Func::Scalar&                       obklen,
  const Func::Scalar&                       kbfs,
  const Func::uview_1d<const Func::Spack>&  cldn,
  const Func::Workspace&                    workspace,
  Func::Scalar&                             pblh);

template void Func::pblintd(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Int&                                npbl,
  const Func::uview_1d<const Func::Spack>&  z,
  const Func::uview_1d<const Func::Spack>&  zi,
  const Func::uview_1d<const Func::Spack>&  thl,
  const Func::uview_1d<const Func::Spack>&  ql,
  const Func::uview_1d<Func::Spack>& q,
  const Func::uview_1d<const Func::Spack>&  u,
  const Func::uview_1d<const Func::Spack>&  v,
  const Func::Scalar&                       ustar,
  const Func::Scalar&                       obklen,
  const Func::Scalar&                       kbfs,
  const Func::uview_1d<const Func::Spack>&  cldn,
  const Func::Workspace&                    workspace,
  Func::Scalar&                             pblh);

template void Func::pblintd(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Int&                                npbl,
  const Func::uview_1d<const Func::Spack>&  z,
  const Func::uview_1d<const Func::Spack>&  zi,
  const Func::uview_1d<const Func::Spack>&  thl,
  const Func::uview_1d<const Func::Spack>&  ql,
  const Func::uview_1d<const Func::Spack>& q,
  const Func::uview_1d<const Func::Spack>&  u,
  const Func::uview_1d<const Func::Spack>&  v,
  const Func::Scalar&                       ustar,
  const Func::Scalar&                       obklen,
  const Func::Scalar&                       kbfs,
  const Func::uview_1d<const Func::Spack>&  cldn,
  const Func::Workspace&                    workspace,
  Func::Scalar&                             pblh);

} // namespace shoc
} // namespace scream

#include "shoc_pblintd_impl.hpp"
