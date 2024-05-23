#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::eddy_diffusivities(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::Scalar&                       Ckh,
  const Func::Scalar&                       Ckm,
  const Func::Scalar&                       pblh,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::scratch_view_1d<Func::Spack>& tabs,
  const Func::uview_1d<const Func::Spack>&  shoc_mix,
  const Func::uview_1d<const Func::Spack>&  sterm_zt,
  const Func::uview_1d<const Func::Spack>&  isotropy,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<Func::Spack>&        tkh,
  const Func::uview_1d<Func::Spack>&        tk);

template void Func::eddy_diffusivities(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::Scalar&                       Ckh,
  const Func::Scalar&                       Ckm,
  const Func::Scalar&                       pblh,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<Func::Spack>& tabs,
  const Func::uview_1d<const Func::Spack>&  shoc_mix,
  const Func::uview_1d<const Func::Spack>&  sterm_zt,
  const Func::uview_1d<const Func::Spack>&  isotropy,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<Func::Spack>&        tkh,
  const Func::uview_1d<Func::Spack>&        tk);

template void Func::eddy_diffusivities(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::Scalar&                       Ckh,
  const Func::Scalar&                       Ckm,
  const Func::Scalar&                       pblh,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>& tabs,
  const Func::uview_1d<const Func::Spack>&  shoc_mix,
  const Func::uview_1d<const Func::Spack>&  sterm_zt,
  const Func::uview_1d<const Func::Spack>&  isotropy,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<Func::Spack>&        tkh,
  const Func::uview_1d<Func::Spack>&        tk);

} // namespace shoc
} // namespace scream

#include "shoc_eddy_diffusivities_impl.hpp"
