#ifndef P3_CALCULATE_MASS_OF_COLUMN_HPP
#define P3_CALCULATE_MASS_OF_COLUMN_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calculate_mass_of_column(
    const MemberType& team, const Int& nlev, const uview_1d<const Spack>& qv, const uview_1d<const Spack>& qc, const uview_1d<const Spack>& qr, const uview_1d<const Spack>&qi, const uview_1d<const Spack>& rho, Scalar& w_mass
  )
{
  w_mass = 0.0;
  using ExeSpaceUtils = ekat::ExeSpaceUtils<typename KT::ExeSpace>;


  ExeSpaceUtils::view_reduction(team,0,nlev,
                                [&] (const int k) -> Spack {

    const auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    const auto range_mask = range_pack < nlev;

    Spack return_val(0); //initialize return value for brunt_int
    return_val.set(range_mask, (qv(k) + qc(k) + qi(k) + qr(k)) * rho(k));// compute brunt_int for each column

    return return_val ;

  }, w_mass);
}

} // namespace p3
} // namespace scream

#endif // P3_CALCULATE_MASS_OF_COLUMN_HPP

