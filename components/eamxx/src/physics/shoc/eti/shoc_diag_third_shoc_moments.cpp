#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::diag_third_shoc_moments(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       c_diag_3rd_mom,
  const Func::uview_1d<const Func::Spack>&  w_sec,
  const Func::uview_1d<const Func::Spack>&  thl_sec,
  const Func::uview_1d<const Func::Spack>&  wthl_sec,
  const Func::uview_1d<const Func::Spack>&  isotropy,
  const Func::uview_1d<const Func::Spack>&  brunt,
  const Func::uview_1d<const Func::Spack>&  thetal,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::scratch_view_1d<Func::Spack>& dz_zt,
  const Func::scratch_view_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        w3);

template void Func::diag_third_shoc_moments(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       c_diag_3rd_mom,
  const Func::uview_1d<const Func::Spack>&  w_sec,
  const Func::uview_1d<const Func::Spack>&  thl_sec,
  const Func::uview_1d<const Func::Spack>&  wthl_sec,
  const Func::uview_1d<const Func::Spack>&  isotropy,
  const Func::uview_1d<const Func::Spack>&  brunt,
  const Func::uview_1d<const Func::Spack>&  thetal,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<Func::Spack>& dz_zt,
  const Func::uview_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        w3);

template void Func::diag_third_shoc_moments(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       c_diag_3rd_mom,
  const Func::uview_1d<const Func::Spack>&  w_sec,
  const Func::uview_1d<const Func::Spack>&  thl_sec,
  const Func::uview_1d<const Func::Spack>&  wthl_sec,
  const Func::uview_1d<const Func::Spack>&  isotropy,
  const Func::uview_1d<const Func::Spack>&  brunt,
  const Func::uview_1d<const Func::Spack>&  thetal,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<const Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        w3);

} // namespace shoc
} // namespace scream

#include "shoc_diag_third_shoc_moments_impl.hpp"
