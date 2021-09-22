#ifndef SPA_FUNCTIONS_IMPL_HPP
#define SPA_FUNCTIONS_IMPL_HPP

#include "physics/share/physics_constants.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/grid/point_grid.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "ekat/util/ekat_lin_interp.hpp"

/*-----------------------------------------------------------------
 * The main SPA routines used to convert SPA data into a format that
 * is usable by the rest of the atmosphere processes.
 *
 * SPA or Simple Prescribed Aerosols provides a way to prescribe
 * aerosols for an atmospheric simulation using pre-computed data.
 * 
 * The data is typically provided at a frequency of monthly, and
 * does not necessarily have to be on the same horizontal or vertical
 * domain as the atmospheric simulation.
 *
 * In order to accomodate coarse temporal resolution and a potentially
 * different spatial resolution it is necessary to perform a series
 * of interpolations, which make up the main body of the SPA routines.
 *
 * The interpolations can be broken into three categories.
 * 1. Horizontal Interpolation: TODO - not done yet.
 * The SPA data set does not have to be provided on the same grid as
 * the atmospheric simulation.  Whenever SPA data is loaded, it is
 * interpolated horizontally onto the simulation grid to provide
 * forcing at every location.  This can be done by,
 *   a) Preloaded remapping wieghts which are applied at every
 *      horizontal column.
 *   b) Online calculation of remapping weights given a set of lat/lon
 *      for the source data and comparing it with the lat/lon of each
 *      column in the simulation.  TODO: THIS HAS NOT BEEN IMPLEMENTED YET
 *
 * 2. Temporal Interpolation:
 * As noted above, the SPA data is provided at some fixed frequency.  Typically
 * as monthly data.  As a result, the data must be interpolated to the current
 * time of the simulation at each time step.  Temporal interpolation follows
 * a basic linear interpolation and is performed for all SPA data at all columns
 * and levels.
 * Note: There is also a temporal interpolation of the surface pressure for the SPA
 * data, which is used in the vertical reconstruction of the pressure profile.
 *
 * 3. Vertical Interpolation:
 * Given that the SPA data has been generated elsewhere it is very likely that
 * the vertical pressure profiles of the data won't match the simulation pressure
 * profiles.  The vertical SPA data structure must be remapped onto the simulation
 * pressure profile.
 * This is done using the EKAT linear interpolation code, see /externals/ekat/util/ekat_lin_interp.hpp
 * The SPA pressure profiles are calculated using the surface pressure which was
 * temporally interpolated in the last step and the set of hybrid coordinates (hyam and hybm)
 * that are used in EAM to construct the physics pressure profiles.
 * The SPA data is then projected onto the simulation pressure profile (pmid)
 * using EKAT linear interpolation. 
-----------------------------------------------------------------*/

namespace scream {
namespace spa {

// Helper function
template<typename ScalarT,typename ScalarS>
KOKKOS_INLINE_FUNCTION
ScalarT linear_interp(const ScalarT& x0, const ScalarT& x1, const ScalarS& t_norm);
/*-----------------------------------------------------------------*/
// The main SPA routine which handles projecting SPA data onto the
// horizontal columns and vertical pressure profiles of the atmospheric
// state.
// Inputs:
//   time_state: A structure defined in spa_functions.hpp which handles
//     the current temporal state of the simulation.
//   pressure_state: A structure defined in spa_functions.hpp which handles
//     the vertical pressure profile for the atmospheric simulation state, and
//     all of the data needed to reconstruct the vertical pressure profile for
//     the SPA data.  See hybrid coordinate (hyam,hybm) and surface pressure (PS)
//   data_beg: A structure defined in spa_functions.hpp which handles the full
//     set of SPA data for the beginning of the month.
//   data_end: Similar to data_beg, but for SPA data for the end of the month.
//   data_out: A structure defined in spa_functions.hpp which handles the full
//     set of SPA data projected onto the pressure profile of the current atmosphere
//     state.  This is the data that will be passed to other processes.
//   ncols_atm, nlevs_atm: The number of columns and levels in the simulation grid.
//     (not to be confused with the number of columns and levels used for the SPA data, 
//      which can be different.)
//   nswbands, nlwbands: The number of shortwave (sw) and longwave (lw) aerosol bands 
//     for the data that will be passed to radiation.
template <typename S, typename D>
void SPAFunctions<S,D>
::spa_main(
  const SPATimeState& time_state,
  const SPAPressureState& pressure_state,
  const SPAData&   data_beg,
  const SPAData&   data_end,
  const SPAOutput& data_out,
  const Int ncols_atm,
  const Int nlevs_atm,
  const Int nswbands,
  const Int nlwbands)
{
  // Gather time stamp info
  auto& t_now = time_state.t_now;
  auto& t_beg = time_state.t_beg_month;
  auto& t_len = time_state.days_this_month;

  // For now we require that the Data in and the Data out have the same number of columns.
  EKAT_REQUIRE(ncols_atm==pressure_state.ncols);

  // Set up temporary arrays that will be used for the spa interpolation.
  view_2d<Spack> p_src("p_mid_src",ncols_atm,pressure_state.nlevs), 
                 ccn3_src("ccn3_src",ncols_atm,pressure_state.nlevs);
  view_3d<Spack> aer_g_sw_src("aer_g_sw_src",ncols_atm,nswbands,pressure_state.nlevs),
                 aer_ssa_sw_src("aer_ssa_sw_src",ncols_atm,nswbands,pressure_state.nlevs),
                 aer_tau_sw_src("aer_tau_sw_src",ncols_atm,nswbands,pressure_state.nlevs),
                 aer_tau_lw_src("aer_tau_lw_src",ncols_atm,nlwbands,pressure_state.nlevs);

  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nlevs_atm);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncols_atm, nk_pack);
  // SPA Main loop
  // Parallel loop order:
  // 1. Loop over all horizontal columns (i index)
  // 2. Loop over all aerosol bands (n index) - where applicable
  // 3. Loop over all vertical packs (k index)
  Kokkos::parallel_for(
    "spa main loop",
    policy,
    KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();  // SCREAM column index

    // Get single-column subviews of all 2D inputs, i.e. those that don't have aerosol bands
    const auto& ps_beg_sub                 = pressure_state.ps_this_month(i);
    const auto& ps_end_sub                 = pressure_state.ps_next_month(i);
    const auto& pmid_sub                   = ekat::subview(pressure_state.pmid, i);
    const auto& p_src_sub                  = ekat::subview(p_src, i);

    const auto& ccn3_beg_sub               = ekat::subview(data_beg.CCN3, i);
    const auto& ccn3_end_sub               = ekat::subview(data_end.CCN3, i);
    const auto& ccn3_src_sub               = ekat::subview(ccn3_src, i);

    // First Step: Horizontal Interpolation if needed - Skip for Now
  
    // Second Step: Temporal Interpolation
    // Use basic linear interpolation function y = b + mx
    auto t_norm = (t_now-t_beg)/t_len;
    /* Determine PS for the source data at this time */
    auto ps_src = linear_interp(ps_beg_sub,ps_end_sub,t_norm);
    /* Reconstruct the vertical pressure profile for the data and time interpolation
     * of the data.
     * Note: CCN3 has the same dimensions as pressure so we handle that time interpolation
     *       in this loop as well.  */
    using C = scream::physics::Constants<Real>;
    static constexpr auto P0 = C::P0;
    const Int nk_pack = ekat::npack<Spack>(pressure_state.nlevs);
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
        // Reconstruct vertical temperature profile using the hybrid coordinate system.
        p_src_sub(k)    = ps_src * pressure_state.hybm(k) + P0 * pressure_state.hyam(k);
        // Time interpolation for CCN3
        ccn3_src_sub(k) = linear_interp(ccn3_beg_sub(k),ccn3_end_sub(k),t_norm);
    });
    team.team_barrier();
    /* Loop over all SW variables with nswbands */
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nswbands), [&] (int n) {
      const auto& aer_g_sw_beg_sub           = ekat::subview(data_beg.AER_G_SW, i, n);
      const auto& aer_g_sw_end_sub           = ekat::subview(data_end.AER_G_SW, i, n);
      const auto& aer_g_sw_src_sub           = ekat::subview(aer_g_sw_src, i, n);

      const auto& aer_ssa_sw_beg_sub         = ekat::subview(data_beg.AER_SSA_SW, i, n);
      const auto& aer_ssa_sw_end_sub         = ekat::subview(data_end.AER_SSA_SW, i, n);
      const auto& aer_ssa_sw_src_sub         = ekat::subview(aer_ssa_sw_src, i, n);
  
      const auto& aer_tau_sw_beg_sub         = ekat::subview(data_beg.AER_TAU_SW, i, n);
      const auto& aer_tau_sw_end_sub         = ekat::subview(data_end.AER_TAU_SW, i, n);
      const auto& aer_tau_sw_src_sub         = ekat::subview(aer_tau_sw_src, i, n);

      /* Now loop over fastest index, the number of vertical packs */
      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(team, nk_pack), [&] (Int k) {
        aer_g_sw_src_sub(k)   = linear_interp(aer_g_sw_beg_sub(k),   aer_g_sw_end_sub(k),   t_norm);
        aer_ssa_sw_src_sub(k) = linear_interp(aer_ssa_sw_beg_sub(k), aer_ssa_sw_end_sub(k), t_norm);
        aer_tau_sw_src_sub(k) = linear_interp(aer_tau_sw_beg_sub(k), aer_tau_sw_end_sub(k), t_norm);
      });
    });
    team.team_barrier();
    /* Loop over all LW variables with nlwbands */
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nlwbands), [&] (int n) {
      const auto& aer_tau_lw_beg_sub         = ekat::subview(data_beg.AER_TAU_LW, i, n);
      const auto& aer_tau_lw_end_sub         = ekat::subview(data_end.AER_TAU_LW, i, n);
      const auto& aer_tau_lw_src_sub         = ekat::subview(aer_tau_lw_src, i, n);

      /* Now loop over fastest index, the number of vertical packs */
      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(team, nk_pack), [&] (Int k) {
        aer_tau_lw_src_sub(k) = linear_interp(aer_tau_lw_beg_sub(k), aer_tau_lw_end_sub(k), t_norm);
      });
    });
    team.team_barrier();
  });
  Kokkos::fence();

  // Third Step: Vertical interpolation, project the SPA data onto the pressure profile for this simulation.
  // This is done using the EKAT linear interpolation routine, see /externals/ekat/util/ekat_lin_interp.hpp
  // for more details. 
  using LIV = ekat::LinInterp<Real,Spack::n>;
  Real minthreshold = 0.0;  // Hard-code a minimum value for aerosol concentration to zero.

  LIV VertInterp(ncols_atm,pressure_state.nlevs,nlevs_atm,minthreshold);
  /* Parallel loop strategy:
   * 1. Loop over all simulation columns (i index)
   * 2. Where applicable, loop over all aerosol bands (n index)
   */
  const Int most_bands = std::max(nlwbands, nswbands);
  typename LIV::TeamPolicy band_policy(ncols_atm, ekat::OnGpu<typename LIV::ExeSpace>::value ? most_bands : 1, VertInterp.km2_pack());
  Kokkos::parallel_for("vertical-interp-spa",
    band_policy,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    const int i = team.league_rank();
    /* Setup the linear interpolater for this column. */
    if (team.team_rank()==0) {
      const auto tvr = Kokkos::ThreadVectorRange(team, VertInterp.km2_pack());
    
      VertInterp.setup(team,tvr,
                       ekat::subview(p_src,i),
                       ekat::subview(pressure_state.pmid,i));
    }
    team.team_barrier();
    /* Conduct vertical interpolation for the 2D variable CCN3 */
    if (team.team_rank()==0) {
      const auto tvr = Kokkos::ThreadVectorRange(team, VertInterp.km2_pack());
    
      VertInterp.lin_interp(team,tvr,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(ccn3_src,i),
                            ekat::subview(data_out.CCN3,i));
    }
    /* Conduct vertical interpolation for the LW banded data - nlwbands (n index) */
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nlwbands), [&] (int n) {
      const auto& tvr = Kokkos::ThreadVectorRange(team, VertInterp.km2_pack());
      VertInterp.lin_interp(team,
                            tvr,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_tau_lw_src,i,n),
                            ekat::subview(data_out.AER_TAU_LW,i,n));
    });
    /* Conduct vertical interpolation for the SW banded data - nswbands (n index) */
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nswbands), [&] (int n) {
      const auto& tvr = Kokkos::ThreadVectorRange(team, VertInterp.km2_pack());
      VertInterp.lin_interp(team,
                            tvr,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_g_sw_src,i,n),
                            ekat::subview(data_out.AER_G_SW,i,n));
      VertInterp.lin_interp(team,
                            tvr,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_ssa_sw_src,i,n),
                            ekat::subview(data_out.AER_SSA_SW,i,n));
      VertInterp.lin_interp(team,
                            tvr,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_tau_sw_src,i,n),
                            ekat::subview(data_out.AER_TAU_SW,i,n));
    });
  });
  Kokkos::fence();

}
/*-----------------------------------------------------------------*/
// WARNING:  We currently only support remap weights from an ncremap file
//           We will probably want to make this more general, perhaps with
//           an enum that switches between valid remap file types.
// TODO: Fix this so that we can use different types of remap file.
template <typename S, typename D>
void SPAFunctions<S,D>
::read_remap_weights_from_file(
  const ekat::Comm& comm,
  const std::string& remap_file,
  const SPAInterp& horiz_weights 
)
{
  // We need to create a local abstract grid object for the input interface to use.  
  // In the case of SPA remap data the current strategy is to have all mpi ranks pull 
  // in all of the remap data and then later parse what is needed.
  auto grid = std::make_shared<PointGrid>("Remap Point",horiz_weights.length,horiz_weights.length,1);
  PointGrid::dofs_list_type dofs_gids ("phys dofs", horiz_weights.length);
  grid->set_dofs(dofs_gids);

  // Note: We need to load both real data and integer data from the remap file.
  //       To accomodate the two different types we go through the same set of
  //       actions in scope for real, and then again for integer.
  using namespace ShortFieldTagsNames;
  FieldLayout layout ({RMP_N_S},{horiz_weights.length});
  {
    using view_h = typename view_1d<Real>::HostMirror;
    std::vector<std::string>          fnames;
    std::map<std::string,view_h>      host_views;
    std::map<std::string,FieldLayout> layouts;

    // The names of the fields to load from file
    fnames.push_back("S");    // weights
    host_views["S"]   = Kokkos::create_mirror_view(horiz_weights.weights);
    layouts.emplace("S",layout);

    // Pointers to the data views to populate
    // Parameter list to control input reading
    ekat::ParameterList weights_params;
    weights_params.set("Fields",fnames);
    weights_params.set("Filename",remap_file);

    AtmosphereInput weights_reader(comm,weights_params);
    // Read data
    weights_reader.init(grid,host_views,layouts);
    weights_reader.read_variables();
    weights_reader.finalize();

    // Map back to device views  
    Kokkos::deep_copy(horiz_weights.weights,host_views["S"]);
  }
  // Currently the scorpio IO interface does not support integer
  // input.  To get around this we manually go through all of the
  // steps usally covered by running the init, read_variables and
  // finalize routines in an io class object.
  // TODO: When the scorpio io interface is able to handle integer
  //       input, the following should be simplified by calling
  //       the io class object directly similar to the above code
  //       used for the weights.
  {
    using view_h = typename view_1d<Int>::HostMirror;
    std::vector<std::string>          fnames;
    std::map<std::string,view_h>      host_views;
    std::map<std::string,FieldLayout> layouts;

    // The names of the fields to load from file
    fnames.push_back("row");  // destination grid index
    host_views["row"] = Kokkos::create_mirror_view(horiz_weights.dst_grid_loc);
    layouts.emplace("row",layout);

    fnames.push_back("col");  // source grid index
    host_views["col"] = Kokkos::create_mirror_view(horiz_weights.src_grid_loc);
    layouts.emplace("col",layout);

    ekat::ParameterList indices_params;
    indices_params.set("Fields",fnames);
    indices_params.set("Filename",remap_file);

    AtmosphereInput indices_reader(comm,indices_params);
    // Init
    //   set_grid
    //   m_layouts
    //   m_host_views
    //   init_scorpio_structures
    //     register_infile
    //     register_variables
    //     set_degrees_of_freedom
    //     set_decomp
    indices_reader.set_grid(grid);
    scorpio::register_infile(remap_file);
    indices_reader.register_variables(remap_file, fnames, layouts, 4); // PIO type INT = 4, so we hard-code here.
    indices_reader.set_degrees_of_freedom(remap_file, fnames, layouts);
    scorpio::set_decomp(remap_file);
    indices_reader.set_m_is_inited(true);
    // TODO: Maybe all we need is to make `init_scorpio_structures` a callable function and don't need the 4 lines above? 
    
    // Read Variables
    //     grid_read_data_array
    //     assign the data read to the actual field, if needed.
    for (auto const& name : fnames) {
      scorpio::grid_read_data_array(remap_file,name,horiz_weights.length,host_views.at(name).data());
    }

    indices_reader.finalize();
    
//    // Map back to device views  
    Kokkos::deep_copy(horiz_weights.src_grid_loc,host_views["col"]);
    Kokkos::deep_copy(horiz_weights.dst_grid_loc,host_views["row"]);
  }
  printf("\n");
  for (int i = 0; i<10; i++) {
    printf(" (%2d) = %16.12f, %d, %d\n",i,horiz_weights.weights(i),horiz_weights.src_grid_loc(i),horiz_weights.dst_grid_loc(i));
  }

}
/*-----------------------------------------------------------------*/
template <typename S, typename D>
void SPAFunctions<S,D>
::spa_update_monthly_data(
  const ekat::Comm& comm,
  const std::string& spa_data_file,
  const util::TimeStamp& ts,
  const SPAInterp& horiz_weights, 
        SPATimeState& time_state,
        SPAPressureState& pressure_state,
        SPAData&   data_beg,
        SPAData&   data_end,
  const Int ncols_atm,
  const Int nlevs_atm,
  const Int nswbands,
  const Int nlwbands)
{
  // We always want to update the current time in the time_state.
  time_state.t_now = ts.get_julian_day();
  // Now we check if we have to update the data that changes monthly
  if (ts.get_months() != time_state.current_month) {
    // Update the SPA time state information
    time_state.current_month = ts.get_months();
    time_state.t_beg_month = util::julian_day(ts.get_years(),ts.get_months(),0,0);
    time_state.days_this_month = (Real)ts.get_dpm();

    // Update the SPA forcing data for this month and next month
    // Start by copying next months data to this months data structure.  
    // NOTE: If the timestep is bigger than monthly this could cause the wrong values
    //       to be assigned.  A timestep greater than a month is very unlikely so we
    //       will proceed.
    data_beg = data_end;
    pressure_state.ps_this_month = pressure_state.ps_next_month;
    // Next we load next months data into the data_end structure and ps_next_month, and apply the horizontal weights.
    auto grid = std::make_shared<PointGrid>("Physics",horiz_weights.src_grid_ncols * horiz_weights.src_grid_nlevs,
                                                             horiz_weights.src_grid_ncols, horiz_weights.src_grid_nlevs);
    PointGrid::dofs_list_type dofs_gids ("phys dofs", horiz_weights.src_grid_ncols);
    grid->set_dofs(dofs_gids);

    using namespace ShortFieldTagsNames;
    FieldLayout scalar2d_layout_mid { {COL}, {horiz_weights.src_grid_ncols} };
    FieldLayout scalar3d_layout_mid { {COL,LEV}, {horiz_weights.src_grid_ncols, horiz_weights.src_grid_nlevs} };
    FieldLayout scalar3d_swband_layout { {COL,SWBND, LEV}, {horiz_weights.src_grid_ncols, nswbands, horiz_weights.src_grid_nlevs} }; 
    FieldLayout scalar3d_lwband_layout { {COL,LWBND, LEV}, {horiz_weights.src_grid_ncols, nlwbands, horiz_weights.src_grid_nlevs} }; 
    using view_h = typename view_1d<Real>::HostMirror;
    std::vector<std::string>          fnames;
    std::map<std::string,view_h>      host_views;
    std::map<std::string,FieldLayout> layouts;
    // Define each input variable we need
    fnames.push_back("PS");
    host_views["PS"] = view_h("",horiz_weights.src_grid_ncols);
    layouts.emplace("PS", scalar2d_layout_mid);
    //
    fnames.push_back("CCN3");
    host_views["CCN3"] = view_h("",horiz_weights.src_grid_ncols*horiz_weights.src_grid_nlevs);
    layouts.emplace("CCN3",scalar3d_layout_mid);
    //
    fnames.push_back("AER_G_SW");
    host_views["AER_G_SW"] = view_h("AER_G_SW",horiz_weights.src_grid_ncols*horiz_weights.src_grid_nlevs*nswbands);
    layouts.emplace("AER_G_SW",scalar3d_swband_layout);
    //
    fnames.push_back("AER_SSA_SW");
    host_views["AER_SSA_SW"] = view_h("AER_SSA_SW",horiz_weights.src_grid_ncols*horiz_weights.src_grid_nlevs*nswbands);
    layouts.emplace("AER_SSA_SW",scalar3d_swband_layout);
    //
    fnames.push_back("AER_TAU_SW");
    host_views["AER_TAU_SW"] = view_h("AER_TAU_SW",horiz_weights.src_grid_ncols*horiz_weights.src_grid_nlevs*nswbands);
    layouts.emplace("AER_TAU_SW",scalar3d_swband_layout);
    //
    fnames.push_back("AER_TAU_LW");
    host_views["AER_TAU_LW"] = view_h("AER_TAU_LW",horiz_weights.src_grid_ncols*horiz_weights.src_grid_nlevs*nlwbands);
    layouts.emplace("AER_TAU_LW",scalar3d_lwband_layout);
    //
    ekat::ParameterList spa_data_in_params;
    spa_data_in_params.set("Fields",fnames);
    spa_data_in_params.set("Filename",spa_data_file);
    AtmosphereInput spa_data_in(comm,spa_data_in_params);
    spa_data_in.init(grid,host_views,layouts);
    spa_data_in.read_variables(time_state.current_month+1);
    spa_data_in.finalize();

  } // if (ts.get_months() != time_state.current_month) 
}
/*-----------------------------------------------------------------*/
// A helper function to manage basic linear interpolation in time.
// The inputs of x0 and x1 represent the data to interpolate from at
// times t0 and t1, respectively.  To keep the signature of the function 
// simple we use
//    t_norm = (t-t0)/(t1-t0).
template<typename ScalarT, typename ScalarS>
inline ScalarT linear_interp(
  const ScalarT& x0,
  const ScalarT& x1,
  const ScalarS& t_norm)
{
  return (1.0 - t_norm) * x0 + t_norm * x1;
}
/*-----------------------------------------------------------------*/

} // namespace spa
} // namespace scream

#endif // SPA_FUNCTIONS_IMPL_HPP
