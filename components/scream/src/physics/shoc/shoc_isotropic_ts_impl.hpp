#ifndef SHOC_ISOTROPIC_TS_IMPL_HPP
#define SHOC_ISOTROPIC_TS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
  namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::isotropic_ts(
  const MemberType&            team,
  const Int&                   nlev,
  const Scalar&                brunt_int,
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& a_diss,
  const uview_1d<const Spack>& brunt,
  const uview_1d<Spack>&       isotropy)
{

  /*-------------------------------------------------
   * Compute the return to isotropic timescale as per
   * Canuto et al. 2004.  This is used to define the
   * eddy coefficients as well as to diagnose higher
   * moments in SHOC
   *-------------------------------------------------*/

  //shared constant
  static constexpr auto ggr = C::gravit;

  //Local constants
  static constexpr Scalar lambda_low   = sp(0.001);
  static constexpr Scalar lambda_high  = sp(0.04);
  static constexpr Scalar lambda_slope = sp(0.65);
  static constexpr Scalar brunt_low    = sp(0.02);
  static constexpr Scalar maxiso       = 20000; // Return to isotropic timescale [s]


  const Int nlev_pack = ekat::pack::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {

      //define the time scale
      const Spack tscale = (2*tke(k))/a_diss(k);

      //define a damping term "lambda" based on column stability
      Scalar lambda = lambda_low + ((brunt_int/ggr)-brunt_low)*lambda_slope;

      const auto lambda_gt_high = lambda > lambda_high;
      if(lambda_gt_high){
        lambda = lambda_high; // take min(lambda_high,lambda)
      }

      const auto lambda_lt_low = lambda < lambda_low;
      if(lambda_lt_low){
        lambda = lambda_low;  //take max(lambda_low,lambda)
      }

      Spack buoy_sgs_save = brunt(k);

      buoy_sgs_save
      if(brunt(k) <=0){
        lambda = 0;
      }

      //Compute the return to isotropic timescale
      isotropy(k) = tscale/(1 + lambda * buoy_sgs_save * tscale * tscale);
      isotropy(k).set(isotropy(k) < maxiso, maxiso); //take min(isotropy,maxiso)
    });
}

} // namespace shoc
} // namespace scream

#endif
