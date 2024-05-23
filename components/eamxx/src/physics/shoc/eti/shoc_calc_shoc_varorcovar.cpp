#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::calc_shoc_varorcovar(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::Scalar&                       tunefac,
  const Func::uview_1d<const Func::Spack>&  isotropy_zi,
  const Func::uview_1d<const Func::Spack>&  tkh_zi,
  const Func::scratch_view_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  invar1,
  const Func::uview_1d<const Func::Spack>&  invar2,
  const Func::uview_1d<Func::Spack>&        varorcovar);

template void Func::calc_shoc_varorcovar(
  const Func::MemberType&                  team,
  const Int&                               nlev,
  const Func::Scalar&                      tunefac,
  const Func::uview_1d<const Func::Spack>& isotropy_zi,
  const Func::uview_1d<const Func::Spack>& tkh_zi,
  const Func::uview_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>& invar1,
  const Func::uview_1d<const Func::Spack>& invar2,
  const Func::uview_1d<Func::Spack>&       varorcovar);

template void Func::calc_shoc_varorcovar(
  const Func::MemberType&                  team,
  const Int&                               nlev,
  const Func::Scalar&                      tunefac,
  const Func::uview_1d<const Func::Spack>& isotropy_zi,
  const Func::uview_1d<const Func::Spack>& tkh_zi,
  const Func::uview_1d<const Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>& invar1,
  const Func::uview_1d<const Func::Spack>& invar2,
  const Func::uview_1d<Func::Spack>&       varorcovar);

} // namespace shoc
} // namespace scream

#include "shoc_calc_shoc_varorcovar_impl.hpp"
