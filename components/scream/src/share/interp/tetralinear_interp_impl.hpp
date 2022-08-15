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
                                     h_weights_, global_columns);

  // Read in all relevant data from source datasets.
  int n_times = 12; // FIXME
  data_.resize(n_times);
  for (int i = 0; i < n_times; ++i) {
    Traits::read_from_file(data_file, i, global_columns, data_[i]);
  }

}

template <typename Data>
void TetralinearInterp<Data>::do_time_interpolation_(Real time, Data& data) {
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
  Traits::linear_combination(1.0 - frac_dt, data_[t1], frac_dt, data_[t2],
                             data_);
}

template <typename Data>
void TetralinearInterp<Data>::
do_vertical_interpolation_(const Data& src_data,
                           const VCoordView& tgt_vcoords,
                           Data& data) {
  VCoordView src_vcoords;
  Traits::compute_vertical_coords(src_data, src_vcoords);
  Traits::interpolate_vertically(src_vcoords, src_data, tgt_vcoords, data);
}

template <typename Data>
void TetralinearInterp<Data>::interpolate_(Real time,
                                           const VCoordView& vcoords,
                                           Data& data) {
  // Perform time interpolation.
  Data data_t = Traits::allocate(data_[0]);
  do_time_interpolation_(time, data_t);

  // Perform vertical interpolation.
  // NOTE: this assumes that the number of vertical levels is the same in the
  // NOTE: source and target data.
  Data data_tv = Traits::allocate(data_t);
  do_vertical_interpolation(data_t, vcoords, data_tv);

  // Perform horizontal interpolation by applying weights to data_tv to obtain
  // data.
  Traits::apply_column_weights(h_weights_, data_tv, data);
}

} // namespace interpolators
} // namespace scream

#endif // TETRALINEAR_INTERP_IMPL_HPP
