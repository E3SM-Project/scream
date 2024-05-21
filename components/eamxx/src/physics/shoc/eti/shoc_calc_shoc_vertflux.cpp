#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing vertflux on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::calc_shoc_vertflux(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::uview_1d<const Func::Spack>&  tkh_zi,
  const Func::scratch_view_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  invar,
  const Func::uview_1d<Func::Spack>&        vertflux);
  
} // namespace shoc
} // namespace scream

#include "shoc_calc_shoc_vertflux_impl.hpp"
