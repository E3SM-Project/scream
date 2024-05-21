#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::update_prognostics_implicit(
  const Func::MemberType&             team,
  const Int&                          nlev,
  const Int&                          nlevi,
  const Int&                          num_qtracers,
  const Func::Scalar&                 dtime,
  const Func::scratch_view_1d<Spack>& dz_zt,
  const Func::scratch_view_1d<Spack>& dz_zi,
  const Func::scratch_view_1d<Spack>& rho_zt,
  const Func::uview_1d<const Spack>&  zt_grid,
  const Func::uview_1d<const Spack>&  zi_grid,
  const Func::uview_1d<const Spack>&  tk,
  const Func::uview_1d<const Spack>&  tkh,
  const Func::Scalar&                 uw_sfc,
  const Func::Scalar&                 vw_sfc,
  const Func::Scalar&                 wthl_sfc,
  const Func::Scalar&                 wqw_sfc,
  const Func::uview_1d<const Spack>&  wtracer_sfc,
  const Func::Workspace&              workspace,
  const Func::uview_1d<Spack>&        thetal,
  const Func::uview_1d<Spack>&        qw,
  const Func::uview_2d<Spack>&        qtracers,
  const Func::uview_1d<Spack>&        tke,
  const Func::uview_1d<Spack>&        u_wind,
  const Func::uview_1d<Spack>&        v_wind);

} // namespace shoc
} // namespace scream

#include "shoc_update_prognostics_implicit_impl.hpp"
