#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::compute_diag_third_shoc_moment(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       c_diag_3rd_mom,
  const Func::uview_1d<const Func::Spack>&  w_sec,
  const Func::uview_1d<const Func::Spack>&  thl_sec,
  const Func::uview_1d<const Func::Spack>&  wthl_sec,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::scratch_view_1d<Func::Spack>& dz_zt,
  const Func::scratch_view_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<const Func::Spack>&  isotropy_zi,
  const Func::uview_1d<const Func::Spack>&  brunt_zi,
  const Func::uview_1d<const Func::Spack>&  w_sec_zi,
  const Func::uview_1d<const Func::Spack>&  thetal_zi,
  const Func::uview_1d<Func::Spack>&        w3);
  
} // namespace shoc
} // namespace scream

#include "shoc_compute_diag_third_shoc_moment_impl.hpp"
