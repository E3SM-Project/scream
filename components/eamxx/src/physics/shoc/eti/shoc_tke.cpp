#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::shoc_tke(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       dtime,
  const Func::Scalar&                       lambda_low,
  const Func::Scalar&                       lambda_high,
  const Func::Scalar&                       lambda_slope,
  const Func::Scalar&                       lambda_thresh,
  const Func::Scalar&                       Ckh,
  const Func::Scalar&                       Ckm,
  const Func::uview_1d<const Func::Spack>&  wthv_sec,
  const Func::uview_1d<const Func::Spack>&  shoc_mix,
  const Func::scratch_view_1d<Func::Spack>& dz_zi,
  const Func::scratch_view_1d<Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>&  pres,
  const Func::scratch_view_1d<Func::Spack>& tabs,
  const Func::uview_1d<const Func::Spack>&  u_wind,
  const Func::uview_1d<const Func::Spack>&  v_wind,
  const Func::uview_1d<const Func::Spack>&  brunt,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::Scalar&                       pblh,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        tke,
  const Func::uview_1d<Func::Spack>&        tk,
  const Func::scratch_view_1d<Func::Spack>& tkh,
  const Func::uview_1d<Func::Spack>&        isotropy);

template void Func::shoc_tke(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       dtime,
  const Func::Scalar&                       lambda_low,
  const Func::Scalar&                       lambda_high,
  const Func::Scalar&                       lambda_slope,
  const Func::Scalar&                       lambda_thresh,
  const Func::Scalar&                       Ckh,
  const Func::Scalar&                       Ckm,
  const Func::uview_1d<const Func::Spack>&  wthv_sec,
  const Func::uview_1d<const Func::Spack>&  shoc_mix,
  const Func::uview_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>&  pres,
  const Func::uview_1d<Func::Spack>& tabs,
  const Func::uview_1d<const Func::Spack>&  u_wind,
  const Func::uview_1d<const Func::Spack>&  v_wind,
  const Func::uview_1d<const Func::Spack>&  brunt,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::Scalar&                       pblh,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        tke,
  const Func::uview_1d<Func::Spack>&        tk,
  const Func::uview_1d<Func::Spack>&        tkh,
  const Func::uview_1d<Func::Spack>&        isotropy);

template void Func::shoc_tke(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       dtime,
  const Func::Scalar&                       lambda_low,
  const Func::Scalar&                       lambda_high,
  const Func::Scalar&                       lambda_slope,
  const Func::Scalar&                       lambda_thresh,
  const Func::Scalar&                       Ckh,
  const Func::Scalar&                       Ckm,
  const Func::uview_1d<const Func::Spack>&  wthv_sec,
  const Func::uview_1d<const Func::Spack>&  shoc_mix,
  const Func::uview_1d<const Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>&  pres,
  const Func::uview_1d<const Func::Spack>& tabs,
  const Func::uview_1d<const Func::Spack>&  u_wind,
  const Func::uview_1d<const Func::Spack>&  v_wind,
  const Func::uview_1d<const Func::Spack>&  brunt,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::Scalar&                       pblh,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        tke,
  const Func::uview_1d<Func::Spack>&        tk,
  const Func::uview_1d<Func::Spack>&        tkh,
  const Func::uview_1d<Func::Spack>&        isotropy);

} // namespace shoc
} // namespace scream

#include "shoc_tke_impl.hpp"
