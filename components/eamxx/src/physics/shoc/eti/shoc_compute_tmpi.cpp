#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::compute_tmpi(
  const Func::MemberType&                   team,
  const Int&                                nlevi,
  const Func::Scalar&                       dtime,
  const Func::uview_1d<const Func::Spack>&  rho_zi,
  const Func::scratch_view_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<Func::Spack>&        tmpi);

template void Func::compute_tmpi(
  const Func::MemberType&                   team,
  const Int&                                nlevi,
  const Func::Scalar&                       dtime,
  const Func::uview_1d<const Func::Spack>&  rho_zi,
  const Func::uview_1d<Func::Spack>& dz_zi,
  const Func::uview_1d<Func::Spack>&        tmpi);

template void Func::compute_tmpi(
  const Func::MemberType&                   team,
  const Int&                                nlevi,
  const Func::Scalar&                       dtime,
  const Func::uview_1d<const Func::Spack>&  rho_zi,
  const Func::uview_1d<const Func::Spack>& dz_zi,
  const Func::uview_1d<Func::Spack>&        tmpi);

} // namespace shoc
} // namespace scream

#include "shoc_compute_tmpi_impl.hpp"
