#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

template struct Functions<Real,DefaultDevice>;

using Func = Functions<Real,DefaultDevice>;

template void Func::shoc_pblintd_init_pot(
  const Func::MemberType&                   team,
  const Int&                                nlev,
  const Func::uview_1d<const Func::Spack>&  thl,
  const Func::uview_1d<const Func::Spack>&  ql,
  const Func::scratch_view_1d<Func::Spack>& q,
  const Func::uview_1d<Func::Spack>&        thv);

} // namespace shoc
} // namespace scream

#include "shoc_pblintd_init_pot_impl.hpp"
