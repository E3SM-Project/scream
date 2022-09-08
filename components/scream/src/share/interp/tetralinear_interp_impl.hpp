#ifndef TETRALINEAR_INTERP_IMPL_HPP
#define TETRALINEAR_INTERP_IMPL_HPP

#include <share/io/scream_scorpio_interface.hpp>

#include <set>
#include <vector>

// This file contains the implementation of TetralinearInterp. Include it
// after any type-specific specializations for TetralinearInterpTraits.

namespace scream {
namespace interpolators {

template <typename DataSet>
void TetralinearInterp<DataSet>::init_from_file_(const std::string& data_file,
                                                 const HCoordView&  tgt_lons,
                                                 const HCoordView&  tgt_lats) {
  // Create a coarse SE grid.
  CoarseGrid coarse_grid(data_file);

  // Map the given (lon, lat) columns to interpolation weights for the given
  // grid.
  std::vector<int> global_columns;
  compute_tetralinear_interp_weights(coarse_grid, tgt_lons, tgt_lats,
                                     weights_, global_columns);

  // Read in all relevant data from source datasets.
  times_ = Traits::dataset_times(data_file);
  EKAT_ASSERT(times_.size() >= 2); // there must be at least 1 dataset!
  EKAT_ASSERT(times_[0] == 0.0);   // the first dataset must have a zero time
  data_.resize(times_.size()-1);   // number of datasets
  for (size_t i = 0; i < times_.size(); ++i) {
    Traits::read_dataset(data_file, int(i), global_columns, data_[i]);
  }
}

template <typename DataSet>
void TetralinearInterp<DataSet>::do_time_interpolation_(Real time,
                                                        DataSet& data) const {
  // If the time falls outside of the set of times covered by the datasets,
  // map it back into that range assuming periodicity.
  Real period = times_.back();
  if (time < 0.0) {        // less than zero!
    time += period;
  } else if (time > period) {  // greater than the period!
    time -= period;
  }

  // Find the bounding times.
  auto time_iter = std::lower_bound(times_.begin(), times_.end(), time);
  EKAT_ASSERT(time_iter != times_.end());
  size_t t1, t2 = std::distance(times_.begin(), time_iter);
  EKAT_ASSERT(t2 > 0); // because the first dataset is defined at t == 0
  if (t2 == times_.size()) {
    // The time indices t1 and t2 refer to the last and first datasets,
    // respectively, because we assume time is periodic.
    t2 = 0;
    t1 = times_.size() - 1;
  } else {
    t1 = t2 - 1;
  }

  // Interpolate!
  Real dt = times_[t2] - times_[t1];
  Real frac_dt = time - times_[t1];
  int num_cols = data.extent(0), num_levels = data.extent(1);
  auto team_policy = ThreadTeamPolicy(num_cols, num_levels);
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
    int i = team.league_rank();
    Traits::linear_combination(team,
                               1.0 - frac_dt, data_[t1], frac_dt, data_[t2],
                               data_);
    });
}

template <typename DataSet>
void TetralinearInterp<DataSet>::
do_horizontal_interpolation_(const DataSet& src_data, DataSet& tgt_data) const {
  EKAT_ASSERT(Traits::num_vertical_levels(src_data) ==
              Traits::num_vertical_levels(tgt_data));
  int num_levels = Traits::num_vertical_levels(src_data);
  auto team_policy = ThreadTeamPolicy(weights_.capacity(), num_levels);
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
    int k = team.league_rank();
    if( weights_.valid_at(k) ) {
      auto col = weights_.key_at(k);
      auto wts = weights_.value_at(k);
      Traits::interpolate_horizontally(team, wts, src_data, col, tgt_data);
    }
  });
}

template <typename DataSet>
void TetralinearInterp<DataSet>::
do_vertical_interpolation_(const VCoordView& src_vcoords,
                           const DataSet& src_data,
                           const VCoordView& tgt_vcoords,
                           DataSet& tgt_data) const {
  int num_cols = tgt_vcoords.extent(0),
      num_src_levels = src_vcoords.extent(1),
      num_tgt_levels = tgt_vcoords.extent(1);
  ekat::LinInterp<Real, Pack::n> vert_interp(num_cols, num_src_levels, num_tgt_levels);

  // Set up the vertical interpolator.
  auto team_col_policy = ThreadTeamPolicy(num_cols, num_tgt_levels);
  Kokkos::parallel_for(team_col_policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
    int i = team.league_rank();
    vert_interp.setup(team,
                      ekat::subview(src_vcoords, i),
                      ekat::subview(tgt_vcoords, i));
  });

  // Now do the interpolation over all columns and variables in parallel.
  int num_vars = Traits::num_variables(tgt_data);
  auto team_col_var_policy = ThreadTeamPolicy(num_cols*num_vars, num_tgt_levels);
  Kokkos::parallel_for(team_col_var_policy,
    KOKKOS_LAMBDA(const ThreadTeam& team) {
      int icol = team.league_rank() / num_vars;
      int ivar = team.league_rank() % num_vars;
      VCoordColumnView src_col_vcoords = ekat::subview(src_vcoords, icol);
      VCoordColumnView tgt_col_vcoords = ekat::subview(tgt_vcoords, icol);
      Traits::interpolate_vertically(team, vert_interp, icol, ivar,
                                     src_vcoords, src_data,
                                     tgt_vcoords, tgt_data);
  });
}

template <typename DataSet>
void TetralinearInterp<DataSet>::
compute_vertical_coords_(const DataSet& data,
                         VCoordView& vcoords) const {
  EKAT_ASSERT(Traits::num_columns(data) == vcoords.extent(0));
  EKAT_ASSERT(Traits::num_vertical_levels(data) == vcoords.extent(1));
  auto team_policy = ThreadTeamPolicy(vcoords.extent(0), vcoords.extent(1));
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
    int i = team.league_rank();
    Traits::compute_vertical_coords(team, i, data,
                                    ekat::subview(vcoords, i));
  });
}

template <typename DataSet>
void TetralinearInterp<DataSet>::
interpolate_(Real time,
             const VCoordView& tgt_vcoords,
             DataSet& tgt_data) const {
  // Perform time interpolation on the source data.
  int num_src_cols = Traits::num_columns(data_[0]);
  int num_src_levels = Traits::num_vertical_levels(data_[0]);
  DataSet data_t = Traits::allocate(num_src_cols, num_src_levels);
  do_time_interpolation_(time, data_t);

  // Perform horizontal interpolation on the time-interpolated data.
  int num_tgt_cols = Traits::num_columns(data_[0]);
  DataSet data_th = Traits::allocate(num_tgt_cols, num_src_levels);
  do_horizontal_interpolation_(data_t, data_th);

  // Compute vertical coordinates on all columns for the temporally (t) and
  // horizontally (h) interpolated dataset.
  VCoordView src_vcoords("src_vcoords", num_tgt_cols, num_src_levels);
  compute_vertical_coords_(data_th, src_vcoords);

  // Perform vertical interpolation.
  // NOTE: this assumes that the number of vertical levels is the same in the
  // NOTE: source and target data.
  do_vertical_interpolation_(src_vcoords, data_th, tgt_vcoords, tgt_data);
}

} // namespace interpolators
} // namespace scream

#endif // TETRALINEAR_INTERP_IMPL_HPP
