#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::vd_shoc_solve(
  const Func::MemberType&       team,
  const Func::uview_1d<Scalar>& du,
  const Func::uview_1d<Scalar>& dl,
  const Func::uview_1d<Scalar>& d,
  const Func::scratch_view_2d<Spack>&  var);

template void Func::vd_shoc_solve(
  const Func::MemberType&       team,
  const Func::uview_1d<Scalar>& du,
  const Func::uview_1d<Scalar>& dl,
  const Func::uview_1d<Scalar>& d,
  const Func::view_2d<Spack>&  var);

template void Func::vd_shoc_solve(
  const Func::MemberType&       team,
  const Func::uview_1d<Scalar>& du,
  const Func::uview_1d<Scalar>& dl,
  const Func::uview_1d<Scalar>& d,
  const Func::uview_2d<Spack>&  var);

} // namespace shoc
} // namespace scream

#include "shoc_tridiag_solver_impl.hpp"
