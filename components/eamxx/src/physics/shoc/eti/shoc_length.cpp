#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template  void Func::shoc_length(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Int&                                nlevi,
  const Func::Scalar&                       length_fac,
  const Func::Scalar&                       dx,
  const Func::Scalar&                       dy,
  const Func::uview_1d<const Func::Spack>&  zt_grid,
  const Func::uview_1d<const Func::Spack>&  zi_grid,
  const Func::scratch_view_1d<Func::Spack>& dz_zt,
  const Func::uview_1d<const Func::Spack>&  tke,
  const Func::uview_1d<const Func::Spack>&  thv,
  const Func::Workspace&                    workspace,
  const Func::uview_1d<Func::Spack>&        brunt,
  const Func::uview_1d<Func::Spack>&        shoc_mix);

} // namespace shoc
} // namespace scream

#include "shoc_length_impl.hpp"
