#ifndef SHOC_EDDY_DIFFUSIVITIES_IMPL_HPP
#define SHOC_EDDY_DIFFUSIVITIES_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc eddy_diffusivities. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::eddy_diffusivities(
  const MemberType&            team,
  const Int&                   nlev,
  const Scalar&                obklen,
  const Scalar&                pblh,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& shoc_mix,
  const uview_1d<const Spack>& sterm_zt,
  const uview_1d<const Spack>& isotropy,
  const uview_1d<const Spack>& tke,
  const uview_1d<Spack>&       tkh,
  const uview_1d<Spack>&       tk)
{
  // Parameters

  // Critical value of Monin-Obukhov length [m]
  const Scalar obk_crit = 0;
  // Transition depth [m] above PBL top to allow
  // stability diffusivities
  const Int pbl_trans = 200;
  // Turbulent coefficients
  const Scalar Ckh = 0.1;
  const Scalar Ckm = 0.1;
  // Eddy coefficients for stable PBL diffusivities
  const Scalar Ckh_s = 0.1;
  const Scalar Ckm_s = 0.1;
  // Maximum allowable value for PBL stability correction
  const Scalar eddycorr_max = 0.1;

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {

    // Default definition of eddy diffusivity for heat and momentum
    tkh(k) = Ckh*isotropy(k)*tke(k);
    tk(k) = Ckm*isotropy(k)*tke(k);

    // If surface layer is stable, based on near surface dimensionless Monin-Obukov
    //  add a correction background diffusivity of which the stength is based
    //  primarily on shear production and SHOC length scale to promote mixing
    //  within the PBL
    const Smask condition = (zt_grid(k) < pblh+pbl_trans) && (obklen > obk_crit);
    tkh(k).set(condition,tkh(k)+ekat::min(Ckh_s*
           ekat::square(shoc_mix(k))*ekat::sqrt(sterm_zt(k)),eddycorr_max));
    tk(k).set(condition,tk(k)+ekat::min(Ckm_s*
           ekat::square(shoc_mix(k))*ekat::sqrt(sterm_zt(k)),eddycorr_max));

  });
}

} // namespace shoc
} // namespace scream

#endif
