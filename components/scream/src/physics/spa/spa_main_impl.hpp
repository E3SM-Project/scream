#ifndef SPA_MAIN_IMPL_HPP
#define SPA_MAIN_IMPL_HPP

#include "physics/spa/spa_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace spa {

/*-----------------------------------------------------------------*/
template <typename S, typename D>
void SPAFunctions<S,D>
::init()
{
  // Translation of pseudo-code written up by @PeterCaldwell.
  // TODO: Delete once full implementation of SPA is complete.
  //
  // At initialization we need to calculate the set of weights and support
  // points for each column in the model simulation.
  //
  // This looks something like:
  // 1. Get first 2 months of data - which includes
  //   - time for each data point, likely the beginning of a month
  //   - ps (surface pressure).  This is needed to construct the
  //     vertical pressure profile at each interpolation time.
  //   - All of the actual data
  // 2. Get weights and support point for horizontal interpolation.  We have two options here:
  //   a) Calculate the wieghts online by going through the lat/lon locations
  //      for each model column and determining which 4 data columns construct
  //      a quadrilateral around it, these are the support points.  Then calculate 
  //      weights based on those 4 points.
  //   b) Use a previously constructred horizontal remap file and read the weights
  //      directly.  This would most likely be an NCO remap file.
  // 
  // We will store these weights and support points for the remainder of the simulation.

}
/*-----------------------------------------------------------------*/
template <typename S, typename D>
void SPAFunctions<S,D>
::main(
  const Int nj,
  const Int nk,
  const Int c_month,
  const InterpolationData InterpData,
  const MonthlyGHG& GHG,
  const MonthlyCCN& CCN,
  const PrescribedAero& prescribed_aero)
{
  // Translation of pseudo-code written up by @PeterCaldwell.
  // TODO: Delete once full implementation of SPA is complete.
  //
  // We minimize the size of the data needed by using the following
  // order to interpolation:
  // 1. Temporal interpolation - Interpolate data at all data columns
  //    to the current model time using linear interpolation.
  // 2. Vertical interpolation/extrapolation - Project the vertical data
  //    from the set of pressure levels according to the data to the current
  //    set of pressure levels.  Note, we will have used temporal interpolation
  //    above to interpolate the surface pressure for the data.
  // 3. Horizontal interpolation - Using the vertically/temporally interpolated
  //    values at each data column interpolate onto the model columns using
  //    previously calculated weights and support points.

  // First check if we have entered a new month.  If so, shift the 2 months of data we
  // use, most likely reading the next months data from file.

  // Loop over all packs in parallel for temporal and vertical interpolation - since these
  // actions are horizontally independent.
  
  // Note, the weights for temporal interpolation only need to be calculated once per timestep:
  Real t_fact_0 = InterpData.current_time_in_days/InterpData.length_of_month_in_days;
  Real t_fact_1 = 1.0 - t_fact_0;

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
    
    // Temporal interpolation of data given the current and next month of data.
    const Int nk_pack = ekat::npack<Spack>(nk);
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
        output(k) = t_fact_0*n_data(k) + t_fact_1*c_data(k);
    }); // Kokkos_parallel_for nk_pack
    team.team_barrier();

  });
  Kokkos::fence();

  // Finish with application of horizontal interpolation.  Note, we need to finish
  // the temporal and vertical interpolation since each model column relies on all
  // columns having been projected into this time/vertical space.

} // main
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
