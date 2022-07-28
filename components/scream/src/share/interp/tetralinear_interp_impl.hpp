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
TetralinearInterp<Data>::init_from_file_(const std::string& data_file,
                                         const HCoordView&  tgt_lats,
                                         const HCoordView&  lgt_lons) {
  // Read ne, np from the file.
  int ne, np;

  // Create a coarse SE grid.
  auto coarse_grid = CoarseSEGrid(ne, np);

  // Map the given (lat, lon) columns to elements in the given grid.
  h_weights_ = compute_tetralinear_interp_weights(coarse_grid, tgt_lats,
                                                  tgt_lons);

  // Fetch all source grid dofs associated with our target points.
  std::vector<int> local_columns;
  {
    std::set<int> unique_local_cols;
    for (auto iter = h_weights_.begin(); iter != h_weights_.end(); ++iter) {
      const auto& col_weights = iter->second;
      unique_local_cols.insert(col_weights.grid_dofs[0]);
      unique_local_cols.insert(col_weights.grid_dofs[1]);
      unique_local_cols.insert(col_weights.grid_dofs[2]);
      unique_local_cols.insert(col_weights.grid_dofs[3]);
    }
    local_columns.resize(unique_local_cols.size());
    std::copy(unique_local_cols.begin(), unique_local_cols.end(),
              local_columns.begin());
  }

  // Read in all relevant data from source datasets.
  data_.resize(n_times);
  for (int i = 0; i < n_times; ++i) {
    Traits::read_from_file(data_file, i, local_columns, data_[i]);
  }

  // Reindex the horizontal weights using local source column indices so we can
  // apply the weights properly to the locally stored source data.
  std::unordered_map<int, int> g2l_columns;
  for (size_t i = 0; i < local_columns.size(); ++i) {
    g2l_columns[local_columns[i]] = i;
  }
  for (auto iter = h_weights_.begin(); iter != h_weights_.end(); ++iter) {
    for (int i = 0; i < 4; ++i) {
      iter->second.columns[i] = g2l_columns[iter->second.columns[i]];
    }
  }
}

template <typename Data>
void TetralinearInterp<Data>::
do_time_interpolation_(Real time, Data& data) {
  // Find the bounding times.
  auto time_iter = std::lower_bound(times_.begin(), times_.end(), time);
  size_t t1, t2;
  if ((time_iter == times_.end()) || (time_iter == times_.begin())) {
    // The times t1 and t2 are the last and first time, respectively, because
    // we assume time is periodic.
    t1 = times.size() - 1;
    t2 = 0;
  } else {
    t2 = std::distance(time_iter);
    t1 = t2 - 1;
  }

  // Interpolate!
  Real dt = times_[t2] - times_[t1];
  Real frac_dt = time - times_[t1];
  Traits::linear_combination(1.0 - frac_dt, data_[t1], frac_dt, data_[t2],
                             data_t);
}

template <typename Data>
void TetralinearInterp<Data>::
do_vertical_interpolation_(const Data& src_data,
                           const VCoordView& tgt_vcoords,
                           Data& data) {
  VCoordView src_vcoords;
  Traits::compute_vertical_coords(src_data, src_vcoords);
  Traits::interpolate_vertically(src_vcoords, src_data, vcoords, data);
}

template <typename Data>
void TetralinearInterp<Data>::interpolate_(Real time,
                                           const VCoordView& vcoords,
                                           Data& data) {
  // Perform time interpolation.
  Data data_t = Traits::allocate(data_[0]);
  do_time_interpolation(time, data_t);

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
