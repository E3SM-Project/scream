#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing compute_shr_prod on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::compute_shr_prod(
  const Func::MemberType&                   team,
  const Int&                                nlevi,
  const Int&                                nlev,
  const Func::scratch_view_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  u_wind,
  const Func::uview_1d<const Func::Spack>&  v_wind,
  const Func::uview_1d<Func::Spack>&        sterm);

} // namespace shoc
} // namespace scream

#include "shoc_compute_shr_prod_impl.hpp"
