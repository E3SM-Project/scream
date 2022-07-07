#ifndef SCREAM_INTERPOLATOR_TRAITS_HPP
#define SCREAM_INTERPOLATOR_TRAITS_HPP

// This file defines some types needed to support the linear interpolator
// (LatLonLinInterpolator) in scream_interpolators.hpp. In particular, it
// defines LatLonLinInterpolatorTraits, which associates a set of functions
// with a type so that type can be interpolated linearly from one set of
// columns at fixed latitude/longitude to another.
//
// If you want to interpolate your own data type, you can add support for it
// in either of two ways:
//
// 1. You can use the generic definition of LatLonLinInterpolatorTraits below
//    and define the methods it calls on your type.
// 2. If you'd rather not modify your type to support linear interpolation, you
//    can include this file and define a specialization of LatLonLinInterpolator
//    for your type. Make sure you define this specialization (or at least
//    declare it) before you include scream_interpolators.hpp, where the traits
//    class is used in the implementation of the interpolator. Why?
//    "Because C++."

#include <share/scream_types.hpp>

#include <ekat/kokkos/ekat_kokkos_types.hpp>

namespace scream {

//------------------------------------------------------------------------
//                        Lat-Lon-Linear Interpolation
//------------------------------------------------------------------------

// This type represents a set of 4 indices associated with another index in a
// mapping.
struct LatLonLinColumnWeights {
  int columns[4];
  Real weights[4];
};

// This type represents a mapping from a target column index to 4 source column
// indices, and is used to store horizontal interpolation weights.
using LatLonLinColumnWeightMap =
  Kokkos::UnorderedMap<int, LatLonLinColumnWeights>;

// LatLonLinInterpolatorTraits -- traits to support linear interpolation from
// a set of source columns at fixed latitudes and longitudes to target data
// on another set of columns. Specialize this type for your data class, or
// implement the methods used in this generic implementation.
template <typename Data>
struct LatLonLinInterpolatorTraits {
  using KokkosTypes = ekat::KokkosTypes;
  using Pack = ekat::Pack<Real, SCREAM_SMALL_PACK_SIZE>;

  // A 1D view storing structured horizontal coordinate data
  using HCoordView = KokkosTypes::view_1d<Real>;

  // A 2D view storing vertical (packed) coordinate data. The first index is
  // the column index, and the second is the vertical pack index.
  using VCoordView = KokkosTypes::view_2d<Pack>;

  // Returns a new Data with storage identical to the given Data.
  static Data allocate(const Data& prototype) {
    return prototype.allocate();
  }

  // Reads data from the given file at the given time index, storing the
  // result in data.
  static void read_from_file(const std::string& data_file, int time_index,
                             Data& data) {
    data.read_from_file(data_file, time_index);
  }

  // Forms the linear combination y = a*x1 + b*x2 where a and b are scalar
  // weights and x1, x2, and y are Data. This method is used for time
  // interpolation.
  static void linear_combination(Real a, Real b,
                                 const Data& x1, const Data& x2,
                                 Data& y) {
    // This generic implementation just calls the linear_combination method to
    // form y = a*x1 + b*x2 + c*x3 + d*x4 in place.
    y.linear_combination(a, b, x1, x2);
  }

  // Computes vertical coordinates for the given Data, populating vcoords.
  static void compute_vertical_coords(const Data& data,
                                      VCoordView& vcoords) {
    // This generic implementation just calls the compute_vertical_coords method
    // on the Data type.
    data.compute_vertical_coords(policy, vcoords);
  }

  // Performs vertical interpolation from a source set of vertical coordinates
  // to a target set.
  static void interpolate_vertically(const VCoordView& src_vcoords,
                                     const Data& src_data,
                                     const VCoordView& tgt_vcoords,
                                     Data& tgt_data) {
    // This generic implementation just calls the interpolate_vertically method
    // on the Data type.
    tgt_data.interpolate_vertically(src_vcoords, src_data, tgt_vcoords);
  }

  // Applies the given set of (horizontal) column weights to column data in x,
  // generating a weighted linear combination of x's data to compute y.
  // This method is used for horizontal interpolation.
  static void apply_column_weights(const LatLonColumnWeightMap& W,
                                   const Data& x, Data& y) {
    y.apply_column_weights(W, x);
  }
};

} // namespace scream

#endif // SCREAM_INTERPOLATOR_TRAITS_HPP
