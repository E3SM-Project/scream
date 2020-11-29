#ifndef SHOC_ADV_SGS_TKE_IMPL_HPP
#define SHOC_ADV_SGS_TKE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc adv_sgs_tke. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::adv_sgs_tke(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   shcol,
  const Real&                  dtime,
  const uview_1d<const Spack>& shoc_mix,
  const uview_1d<const Spack>& wthv_sec,
  const uview_1d<const Spack>& sterm_zt,
  const uview_1d<const Spack>& tk,
  const uview_1d<Spack>&       tke,
  const uview_1d<Spack>&       a_diss)
{

  //constants shared between physics
  static constexpr Scalar ggr      = C::gravit;
  static constexpr Scalar basetemp = C::basetemp;
  static constexpr Scalar mintke   = scream::shoc::Constants<Real>::mintke;
  static constexpr Scalar maxtke   = scream::shoc::Constants<Real>::maxtke;

  //declare some constants
  static constexpr Scalar Cs  = 0.15;
  static constexpr Scalar Ck  = 0.1;
  Spack Ce  = ekat::cube(Spack(Ck))/std::pow(Cs,4); //BALLI
  Spack Ce1 = Ce/sp(0.7)*sp(0.19);
  Spack Ce2 = Ce/sp(0.7)*sp(0.51);
  Spack Cee = Ce1 + Ce2;


  const Int nlev_pack = ekat::npack<Spack>(nlev);

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {

      // Compute buoyant production term
      const Spack a_prod_bu = (ggr/basetemp)*wthv_sec(k);

      tke(k)    = ekat::max(0,tke(k));

      // Shear production term, use diffusivity from previous timestep
      const Spack a_prod_sh = tk(k)*sterm_zt(k);

      // Dissipation term
      a_diss(k)=Cee/shoc_mix(k)*tke(k)*1.5;//balli it is **1.5

      // March equation forward one timestep
      tke(k)=ekat::max(mintke,tke(k)+dtime*(ekat::max(0,a_prod_sh+a_prod_bu)-a_diss(k)));

      tke(k)=ekat::min(tke(k),maxtke);
    });

}

} // namespace shoc
} // namespace scream

#endif
