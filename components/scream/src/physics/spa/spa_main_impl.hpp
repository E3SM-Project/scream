#ifndef SPA_MAIN_IMPL_HPP
#define SPA_MAIN_IMPL_HPP

#include "physics/spa/spa_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace spa {

/*-----------------------------------------------------------------*/
template <typename S, typename D>
void SPAFunctions<S,D>
::main(
  const Int nj,
  const Int nk,
  const Int c_month,
  const Real days_from_c_month,
  const MonthlyGHG& GHG,
  const MonthlyCCN& CCN,
  const PrescribedAero& prescribed_aero)
{
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);
  Kokkos::parallel_for(
    "spa main loop",
    policy,
    KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();
    const auto c_data = ekat::subview(CCN.CCN,i,c_month);
    const auto n_data = ekat::subview(CCN.CCN,i,(c_month+1) % 12);
    const auto output = ekat::subview(prescribed_aero.CCN,i);
    
    printf("ASD spa_main: %d, %f\n",i,days_from_c_month);
 
    // Perform the spatial interpolation
    //  - horizontal interpolation
    //    interpolate from coarse data to fine resolution.
    //  - vertical interpolation given that we may
    //    a) have a different number of levels in the simulation then from the data, and
    //    b) we have floating pressure levels.
    //  NOTE: The spatial interpolation infrastructure has not been implemented yet.  TODO: implement when possible.

    // Interpolate data given the current and next month of data.
    temporal_interpolation(team, nk, days_from_c_month, c_data, n_data, output);

  });
  Kokkos::fence();
} // main
/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_FUNCTION
void SPAFunctions<S,D>
::temporal_interpolation(
  const MemberType& team,
  const Int& nk,
  const Real& c_time_in_days,
  const uview_1d<const Spack>& c_month_data,
  const uview_1d<const Spack>& n_month_data,
  const uview_1d<Spack>&       interp_data)
{
  team.team_barrier();
  const Int nk_pack = ekat::npack<Spack>(nk);
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
      printf("ASD spa_temporal: %d\n",k);
      auto slp = (n_month_data(k) - c_month_data(k))/31.0;
      interp_data(k) = c_month_data(k) + slp*c_time_in_days;
  }); // Kokkos_parallel_for nk_pack
  team.team_barrier();
} // calc_icefrac
/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_FUNCTION
void SPAFunctions<S,D>
::spatial_interpolation(
  const MemberType& team,
  const Int& nk,
  const uview_1d<const Spack>& c_month_data,
  const uview_1d<const Spack>& n_month_data,
  const uview_1d<Spack>& c_month_interp_data,
  const uview_1d<Spack>& n_month_interp_data)
{
  team.team_barrier();
  const Int nk_pack = ekat::npack<Spack>(nk);
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
      //Do Nothing
  }); // Kokkos_parallel_for nk_pack
  team.team_barrier();
} // calc_totalfrac
/*-----------------------------------------------------------------*/

} // namespace spa
} // namespace scream

#endif // SPA_MAIN_IMPL_HPP
