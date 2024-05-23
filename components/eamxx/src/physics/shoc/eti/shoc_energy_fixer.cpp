#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::shoc_energy_fixer(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       dtime,
  const Int&                                nadv,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::Scalar&                       se_b,
  const Func::Scalar&                       ke_b,
  const Func::Scalar&                       wv_b,
  const Func::Scalar&                       wl_b,
  const Func::Scalar&                       se_a,
  const Func::Scalar&                       ke_a,
  const Func::Scalar&                       wv_a,
  const Func::Scalar&                       wl_a,
  const Func::Scalar&                       wthl_sfc,
  const Func::Scalar&                       wqw_sfc,
  const Func::scratch_view_1d<Func::Spack>& rho_zt,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<const Func::Spack>&  pint,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        host_dse);

template void Func::shoc_energy_fixer(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       dtime,
  const Int&                                nadv,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::Scalar&                       se_b,
  const Func::Scalar&                       ke_b,
  const Func::Scalar&                       wv_b,
  const Func::Scalar&                       wl_b,
  const Func::Scalar&                       se_a,
  const Func::Scalar&                       ke_a,
  const Func::Scalar&                       wv_a,
  const Func::Scalar&                       wl_a,
  const Func::Scalar&                       wthl_sfc,
  const Func::Scalar&                       wqw_sfc,
  const Func::uview_1d<Func::Spack>& rho_zt,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<const Func::Spack>&  pint,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        host_dse);

template void Func::shoc_energy_fixer(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       dtime,
  const Int&                                nadv,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::Scalar&                       se_b,
  const Func::Scalar&                       ke_b,
  const Func::Scalar&                       wv_b,
  const Func::Scalar&                       wl_b,
  const Func::Scalar&                       se_a,
  const Func::Scalar&                       ke_a,
  const Func::Scalar&                       wv_a,
  const Func::Scalar&                       wl_a,
  const Func::Scalar&                       wthl_sfc,
  const Func::Scalar&                       wqw_sfc,
  const Func::uview_1d<const Func::Spack>& rho_zt,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<const Func::Spack>&  pint,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        host_dse);

} // namespace shoc
} // namespace scream

#include "shoc_energy_fixer_impl.hpp"
