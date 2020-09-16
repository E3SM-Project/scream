#ifndef SHOC_ADV_SGS_TKE_IMPL_HPP
#define SHOC_ADV_SGS_TKE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "shoc_constants.hpp" // for ETI only but harmless for GPU

namespace scream {
  namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::adv_sgs_tke(
  const MemberType& team,
  const Int& nlev,
  const Real& dtime,
  const uview_1d<const Spack>& shoc_mix,
  const uview_1d<const Spack>& wthv_sec,
  const uview_1d<const Spack>& sterm_zt,
  const uview_1d<const Spack>& tk,
  const uview_1d<Spack>& tke,
  const uview_1d<Spack>& a_diss)
{

  //constants shared between physics
  static constexpr auto ggr      = C::gravit;
  static constexpr auto basetemp = C::basetemp;

  //constants shared within SHOC only
  static constexpr auto mintke   = SC::mintke;
  static constexpr auto maxtke   = SC::maxtke;

  //Local constants
  static constexpr Scalar Cs = sp(0.15);
  static constexpr Scalar Ck = sp(0.1);
  static constexpr Scalar Ce = pow(Ck,3)/pow(Cs,4);

  static constexpr Scalar Ce1 = Ce/sp(0.7)*sp(0.19);
  static constexpr Scalar Ce2 = Ce/sp(0.7)*sp(0.51);
  static constexpr Scalar Cee = Ce1 + Ce2;

  const Int nlev_pack = ekat::pack::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {

    // Find max(0,tke)
    tke(k).set(tke(k)<0, 0);

    // Compute buoyant production term
    const auto a_prod_bu = (ggr/basetemp) * wthv_sec(k);

    // Shear production term, use diffusivity from previous timestep
    const auto a_prod_sh = tk(k) * sterm_zt(k);

    // Dissipation term
    a_diss(k) = Cee/shoc_mix(k)*pow(tke(k),sp(1.5));

    //compute total production and take max(0, total production)
    auto prodTotal = a_prod_sh + a_prod_bu;
    prodTotal.set(prodTotal < 0, 0); //take max(0, prodTotal)


    //-----------------------------------------------------
    // March equation forward one timestep
    //-----------------------------------------------------
    tke(k) = tke(k) + dtime * (prodTotal - a_diss(k));

    //Take max(mintke, tke_tmp1)
    tke(k).set(tke(k) < mintke, mintke);

    //Take min(tke,maxtke)
    tke(k).set(tke(k)>maxtke, maxtke);

  });
}

} // namespace shoc
} // namespace scream

#endif
