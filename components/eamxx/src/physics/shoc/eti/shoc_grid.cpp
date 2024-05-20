#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::shoc_grid(
  const Func::MemberType&                  team,
  const Int&                               nlev,
  const Int&                               nlevi,
  const Func::uview_1d<const Func::Spack>& zt_grid,
  const Func::uview_1d<const Func::Spack>& zi_grid,
  const Func::uview_1d<const Func::Spack>& pdel,
  const Func::uview_1d<Func::Spack>&       dz_zt,
  const Func::uview_1d<Func::Spack>&       dz_zi,
  const Func::scratch_view_1d<Spack>&      rho_zt);

template void Func::shoc_grid(
  const Func::MemberType&                  team,
  const Int&                               nlev,
  const Int&                               nlevi,
  const Func::uview_1d<const Func::Spack>& zt_grid,
  const Func::uview_1d<const Func::Spack>& zi_grid,
  const Func::uview_1d<const Func::Spack>& pdel,
  const Func::uview_1d<Func::Spack>&       dz_zt,
  const Func::uview_1d<Func::Spack>&       dz_zi,
  const Func::uview_1d<Func::Spack>&       rho_zt);

} // namespace shoc
} // namespace scream

#include "shoc_grid_impl.hpp"
