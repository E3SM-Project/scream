#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::integ_column_stability(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::scratch_view_1d<Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>&  pres,
  const Func::uview_1d<const Func::Spack>&  brunt,
  Func::Scalar&                             brunt_int);

} // namespace shoc
} // namespace scream

#include "shoc_integ_column_stability_impl.hpp"
