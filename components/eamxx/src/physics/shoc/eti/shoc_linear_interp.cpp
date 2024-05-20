#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing linear_interp on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::linear_interp(
  const Func::MemberType&                   team,
  const Func::uview_1d<const Func::Spack>&  x1,
  const Func::uview_1d<const Func::Spack>&  x2,
  const Func::uview_1d<const Func::Spack>&  y1,
  const Func::scratch_view_1d<Func::Spack>& y2,
  const Int&                                km1,
  const Int&                                km2,
  const Func::Scalar&                       minthresh);

template void Func::linear_interp(
  const Func::MemberType&                   team,
  const Func::uview_1d<const Func::Spack>&  x1,
  const Func::uview_1d<const Func::Spack>&  x2,
  const Func::scratch_view_1d<Func::Spack>& y1,
  const Func::uview_1d<Func::Spack>&        y2,
  const Int&                                km1,
  const Int&                                km2,
  const Func::Scalar&                       minthresh);

template void Func::linear_interp(
  const Func::MemberType&            team,
  const Func::uview_1d<Func::Spack>& x1,
  const Func::uview_1d<Func::Spack>& x2,
  const Func::uview_1d<Func::Spack>& y1,
  const Func::uview_1d<Func::Spack>& y2,
  const Int&                         km1,
  const Int&                         km2,
  const Func::Scalar&                minthresh);

template void Func::linear_interp(
  const Func::MemberType&                  team,
  const Func::uview_1d<const Func::Spack>& x1,
  const Func::uview_1d<const Func::Spack>& x2,
  const Func::uview_1d<Func::Spack>&       y1,
  const Func::uview_1d<Func::Spack>&       y2,
  const Int&                               km1,
  const Int&                               km2,
  const Func::Scalar&                      minthresh);

template void Func::linear_interp(
  const Func::MemberType&                  team,
  const Func::uview_1d<const Func::Spack>& x1,
  const Func::uview_1d<const Func::Spack>& x2,
  const Func::uview_1d<const Func::Spack>& y1,
  const Func::uview_1d<Func::Spack>&       y2,
  const Int&                               km1,
  const Int&                               km2,
  const Func::Scalar&                      minthresh);

} // namespace shoc
} // namespace scream

#include "shoc_linear_interp_impl.hpp"
