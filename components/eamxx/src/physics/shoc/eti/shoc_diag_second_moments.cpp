#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing diag_second_moments on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::diag_second_moments(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Real&                               thl2tune,
  const Real&                               qw2tune,
  const Real&                               qwthl2tune,
  const Real&                               w2tune,
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
  const Func::uview_1d<Func::Spack>&        isotropy_zi,
  const Func::uview_1d<Func::Spack>&        tkh_zi,
  const Func::uview_1d<Func::Spack>&        tk_zi,
  const Func::uview_1d<Func::Spack>&        thl_sec,
  const Func::uview_1d<Func::Spack>&        qw_sec,
  const Func::uview_1d<Func::Spack>&        wthl_sec,
  const Func::uview_1d<Func::Spack>&        wqw_sec,
  const Func::uview_1d<Func::Spack>&        qwthl_sec,
  const Func::uview_1d<Func::Spack>&        uw_sec,
  const Func::uview_1d<Func::Spack>&        vw_sec,
  const Func::uview_1d<Func::Spack>&        wtke_sec,
  const Func::uview_1d<Func::Spack>&        w_sec);

template void Func::diag_second_moments(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Real&                               thl2tune,
  const Real&                               qw2tune,
  const Real&                               qwthl2tune,
  const Real&                               w2tune,
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
  const Func::uview_1d<Func::Spack>&        isotropy_zi,
  const Func::uview_1d<Func::Spack>&        tkh_zi,
  const Func::uview_1d<Func::Spack>&        tk_zi,
  const Func::uview_1d<Func::Spack>&        thl_sec,
  const Func::uview_1d<Func::Spack>&        qw_sec,
  const Func::uview_1d<Func::Spack>&        wthl_sec,
  const Func::uview_1d<Func::Spack>&        wqw_sec,
  const Func::uview_1d<Func::Spack>&        qwthl_sec,
  const Func::uview_1d<Func::Spack>&        uw_sec,
  const Func::uview_1d<Func::Spack>&        vw_sec,
  const Func::uview_1d<Func::Spack>&        wtke_sec,
  const Func::uview_1d<Func::Spack>&        w_sec);

template void Func::diag_second_moments(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Real&                               thl2tune,
  const Real&                               qw2tune,
  const Real&                               qwthl2tune,
  const Real&                               w2tune,
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
  const Func::uview_1d<Func::Spack>&        isotropy_zi,
  const Func::uview_1d<Func::Spack>&        tkh_zi,
  const Func::uview_1d<Func::Spack>&        tk_zi,
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

#include "shoc_diag_second_moments_impl.hpp"
