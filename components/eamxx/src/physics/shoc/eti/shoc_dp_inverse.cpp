#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::dp_inverse(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::scratch_view_1d<Func::Spack>& rho_zt,
  const Func::scratch_view_1d<Func::Spack>& dz_zt,
  const Func::uview_1d<Func::Spack>&        rdp_zt);

} // namespace shoc
} // namespace scream

#include "shoc_dp_inverse_impl.hpp"
