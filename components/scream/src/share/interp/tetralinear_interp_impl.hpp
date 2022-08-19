#ifndef TETRALINEAR_INTERP_IMPL_HPP
#define TETRALINEAR_INTERP_IMPL_HPP

#include <share/io/scream_scorpio_interface.hpp>

#include <set>
#include <vector>

// This file contains the implementation of TetralinearInterp. Include it
// after any type-specific specializations for TetralinearInterpTraits.

namespace scream {
namespace interpolators {

template <typename Data>
void TetralinearInterp<Data>::init_from_file_(const std::string& data_file,
                                              const HCoordView&  tgt_lats,
                                              const HCoordView&  tgt_lons) {
  // Create a coarse SE grid.
  auto coarse_grid = CoarseGrid(data_file);

  // Map the given (lat, lon) columns to interpolation weights for the given
  // grid.
  std::vector<int> global_columns;
  compute_tetralinear_interp_weights(coarse_grid, tgt_lats, tgt_lons,
                                     weights_, global_columns);

  // Read in all relevant data from source datasets.
  int n_times = 12; // FIXME
  data_.resize(n_times);
  for (int i = 0; i < n_times; ++i) {
    Traits::read_from_file(data_file, i, global_columns, data_[i]);
  }

}

template <typename Data>
void TetralinearInterp<Data>::do_time_interpolation_(Real time, Data& data) const {
  // Find the bounding times.
  auto time_iter = std::lower_bound(times_.begin(), times_.end(), time);
  size_t t1, t2;
  if ((time_iter == times_.end()) || (time_iter == times_.begin())) {
    // The times t1 and t2 are the last and first time, respectively, because
    // we assume time is periodic.
    t1 = times_.size() - 1;
    t2 = 0;
  } else {
    t2 = std::distance(times_.begin(), time_iter);
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

template <typename Data>
void TetralinearInterp<Data>::
do_horizontal_interpolation_(const Data& src_data, Data& tgt_data) const {
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

template <typename Data>
void TetralinearInterp<Data>::
do_vertical_interpolation_(const VCoordView& src_vcoords,
                           const Data& src_data,
                           const VCoordView& tgt_vcoords,
                           Data& tgt_data) const {
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

template <typename Data>
void TetralinearInterp<Data>::interpolate_(Real time,
                                           const VCoordView& src_vcoords,
                                           const VCoordView& tgt_vcoords,
                                           Data& tgt_data) const {
  // Perform time interpolation on the source data.
  int num_src_cols = Traits::num_columns(data_[0]);
  int num_src_levels = Traits::num_vertical_levels(data_[0]);
  Data data_t = Traits::allocate(num_src_cols, num_src_levels);
  do_time_interpolation_(time, data_t);

  // Perform horizontal interpolation on the time-interpolated data.
  int num_tgt_cols = Traits::num_columns(data_[0]);
  Data data_th = Traits::allocate(num_tgt_cols, num_src_levels);
  do_horizontal_interpolation(data_t, data_th);

  // Perform vertical interpolation.
  // NOTE: this assumes that the number of vertical levels is the same in the
  // NOTE: source and target data.
  do_vertical_interpolation(src_vcoords, data_th, tgt_vcoords, tgt_data);
}

} // namespace interpolators
} // namespace scream

#endif // TETRALINEAR_INTERP_IMPL_HPP
