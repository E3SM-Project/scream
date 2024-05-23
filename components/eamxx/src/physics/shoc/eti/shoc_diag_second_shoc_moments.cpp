#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing diag_second_shoc_moments on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::diag_second_shoc_moments(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       thl2tune,
  const Func::Scalar&                       qw2tune,
  const Func::Scalar&                       qwthl2tune,
  const Func::Scalar&                       w2tune,
  const Func::uview_1d<const Func::Spack>&  thetal,
  const Func::uview_1d<const Func::Spack>&  qw,
  const Func::uview_1d<const Func::Spack>&  u_wind,
  const Func::uview_1d<const Func::Spack>&  v_wind,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<const Func::Spack>&  isotropy,
  const Func::scratch_view_1d<Func::Spack>& tkh,
  const Func::uview_1d<const Func::Spack>&  tk,
  const Func::scratch_view_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::uview_1d<const Func::Spack>&  shoc_mix,
  const Func::Scalar&                       wthl_sfc,
  const Func::Scalar&                       wqw_sfc,
  const Func::Scalar&                       uw_sfc,
  const Func::Scalar&                       vw_sfc,
  Func::Scalar&                             ustar2,
  Func::Scalar&                             wstar,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        thl_sec,
  const Func::uview_1d<Func::Spack>&        qw_sec,
  const Func::uview_1d<Func::Spack>&        wthl_sec,
  const Func::uview_1d<Func::Spack>&        wqw_sec,
  const Func::uview_1d<Func::Spack>&        qwthl_sec,
  const Func::uview_1d<Func::Spack>&        uw_sec,
  const Func::uview_1d<Func::Spack>&        vw_sec,
  const Func::uview_1d<Func::Spack>&        wtke_sec,
  const Func::uview_1d<Func::Spack>&        w_sec);

template void Func::diag_second_shoc_moments(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       thl2tune,
  const Func::Scalar&                       qw2tune,
  const Func::Scalar&                       qwthl2tune,
  const Func::Scalar&                       w2tune,
  const Func::uview_1d<const Func::Spack>&  thetal,
  const Func::uview_1d<const Func::Spack>&  qw,
  const Func::uview_1d<const Func::Spack>&  u_wind,
  const Func::uview_1d<const Func::Spack>&  v_wind,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<const Func::Spack>&  isotropy,
  const Func::uview_1d<Func::Spack>&  tkh,
  const Func::uview_1d<const Func::Spack>&  tk,
  const Func::uview_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::uview_1d<const Func::Spack>&  shoc_mix,
  const Func::Scalar&                       wthl_sfc,
  const Func::Scalar&                       wqw_sfc,
  const Func::Scalar&                       uw_sfc,
  const Func::Scalar&                       vw_sfc,
  Func::Scalar&                             ustar2,
  Func::Scalar&                             wstar,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        thl_sec,
  const Func::uview_1d<Func::Spack>&        qw_sec,
  const Func::uview_1d<Func::Spack>&        wthl_sec,
  const Func::uview_1d<Func::Spack>&        wqw_sec,
  const Func::uview_1d<Func::Spack>&        qwthl_sec,
  const Func::uview_1d<Func::Spack>&        uw_sec,
  const Func::uview_1d<Func::Spack>&        vw_sec,
  const Func::uview_1d<Func::Spack>&        wtke_sec,
  const Func::uview_1d<Func::Spack>&        w_sec);

template void Func::diag_second_shoc_moments(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       thl2tune,
  const Func::Scalar&                       qw2tune,
  const Func::Scalar&                       qwthl2tune,
  const Func::Scalar&                       w2tune,
  const Func::uview_1d<const Func::Spack>&  thetal,
  const Func::uview_1d<const Func::Spack>&  qw,
  const Func::uview_1d<const Func::Spack>&  u_wind,
  const Func::uview_1d<const Func::Spack>&  v_wind,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<const Func::Spack>&  isotropy,
  const Func::uview_1d<const Func::Spack>&  tkh,
  const Func::uview_1d<const Func::Spack>&  tk,
  const Func::uview_1d<const Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::uview_1d<const Func::Spack>&  shoc_mix,
  const Func::Scalar&                       wthl_sfc,
  const Func::Scalar&                       wqw_sfc,
  const Func::Scalar&                       uw_sfc,
  const Func::Scalar&                       vw_sfc,
  Func::Scalar&                             ustar2,
  Func::Scalar&                             wstar,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        thl_sec,
  const Func::uview_1d<Func::Spack>&        qw_sec,
  const Func::uview_1d<Func::Spack>&        wthl_sec,
  const Func::uview_1d<Func::Spack>&        wqw_sec,
  const Func::uview_1d<Func::Spack>&        qwthl_sec,
  const Func::uview_1d<Func::Spack>&        uw_sec,
  const Func::uview_1d<Func::Spack>&        vw_sec,
  const Func::uview_1d<Func::Spack>&        wtke_sec,
  const Func::uview_1d<Func::Spack>&        w_sec);

} // namespace shoc
} // namespace scream

#include "shoc_diag_second_shoc_moments_impl.hpp"
