#ifndef TETRALINEAR_INTERP_TRAITS_HPP
#define TETRALINEAR_INTERP_TRAITS_HPP

// This file defines some types needed to support the TetralinearInterp type.
// In particular, it defines TetralinearInterpTraits, which associates a set of
// functions with a type so that type can be interpolated linearly from a coarse
// cubed-sphere grid to a set of latitude/longitude target points.
//
// If you want to interpolate your own data type, you can add support for it
// in either of two ways:
//
// 1. You can use the generic definition of TetralinearInterpTraits below
//    and define the methods it calls on your type.
// 2. If you'd rather not modify your type to support linear interpolation, you
//    can include this file and define a specialization of
//    TetralinearInterpTraits for your type. Make sure you define this
//    specialization (or at least declare it) before you include
//    tetralinear_interp.hpp, where the traits class is used in the
//    implementation of the interpolator. Why? "Because C++."

#include <share/scream_types.hpp>

#include <ekat/ekat_pack.hpp>
#include <ekat/kokkos/ekat_kokkos_types.hpp>
#include <ekat/kokkos/ekat_subview_utils.hpp>
#include <ekat/util/ekat_lin_interp.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include <type_traits>

namespace scream {
namespace interpolators {

using KokkosTypes = ekat::KokkosTypes<ekat::DefaultDevice>;
using KokkosHostTypes = ekat::KokkosTypes<ekat::HostDevice>;
using Pack = ekat::Pack<Real, SCREAM_SMALL_PACK_SIZE>;

// A 1D view storing structured horizontal coordinate data (on device)
using HCoordView = KokkosTypes::view_1d<Real>;

// A 1D view storing structured horizontal coordinate data (on host)
using HostHCoordView = KokkosHostTypes::view_1d<Real>;

// A 2D view storing vertical (packed) coordinate data (on device). The first
// index is the column index, and the second is the vertical pack index.
using VCoordView = KokkosTypes::view_2d<Pack>;

// A 2D view storing vertical (packed) coordinate data (on host). The first
// index is the column index, and the second is the vertical pack index.
using HostVCoordView = KokkosHostTypes::view_2d<Pack>;


/* -- C++17 stuff!
// A 1D view storing a single column of vertical coordinate data, equivalent
// to the return type of ekat::subview applied to VCoordView.
using VCoordColumnView = std::invoke_result_t<ekat::subview, VCoordView, int>;
 */

/* Replace this with the above C++17 stuff when we switch over. */
using VCoordColumnView = KokkosTypes::view_1d<Pack>;

// Policy and thread team type for doing parallel dispatch on interpolated data.
using ThreadTeamPolicy = typename KokkosTypes::TeamPolicy;
using ThreadTeam = typename KokkosTypes::MemberType;

//------------------------------------------------------------------------
//                        Tetralinear Interpolation
//------------------------------------------------------------------------

// This type represents a set of indices and weights associated with a the
// support for a target latitude/longitude pair.
struct TetralinearInterpWeights {
  int  indices[4]; // global grid element vertex indices
  Real weights[4]; // interpolation weights for each vertex
};

// This type represents an on-device mapping from a target column index to 4
// source column indices, and is used to store horizontal interpolation weights.
using TetralinearInterpWeightMap =
  Kokkos::UnorderedMap<int, TetralinearInterpWeights>;

// This is the host version of LatLonLinInterpWeightMap.
using HostTetralinearInterpWeightMap =
  Kokkos::UnorderedMap<int, TetralinearInterpWeights, ekat::HostDevice>;

// TetralinearInterpoTraits -- traits to support tetralinear interpolation
// from a set of source columns at fixed latitudes and longitudes to target data
// on another set of columns. Specialize this type for your data class, or
// implement the methods used in this generic implementation.
template <typename Data>
struct TetralinearInterpTraits {

  // Reads data from the given file at the given time index for the given
  // columns, storing the result in data. Call on HOST only.
  // NOTE: the columns vector defines a local-to-global mapping of required
  // NOTE: source column data, allowing you to store only locally-relevant data.
  // NOTE: (i.e. columns[local_column_index] == global_column_index)
  static void read_from_file(const std::string& filename, int time_index,
                             const std::vector<int>& columns, Data& data) {
    data.read_from_file(filename, time_index, columns);
  }

  // Returns a new Data with storage ample for the given number of columns and
  // vertical levels. Call on HOST only.
  static Data allocate(int num_columns, int num_vertical_levels) {
    return Data::allocate(num_columns, num_vertical_levels);
  }

  // Returns the number of horizontal columns represented in the data.
  static int num_columns(const Data& data) {
    return data.num_columns();
  }

  // Returns the number of vertical levels represented in the data.
  static int num_vertical_levels(const Data& data) {
    return data.num_vertical_levels();
  }

  // Forms the linear combination y = a*x1 + b*x2 on a single column, where a
  // and b are scalar weights and x1, x2, and y are Data. This method forms the
  // linear combination in parallel across all vertical levels of Data using the
  // given thread team, and is used for time interpolation.
  // Column:                    team.league_rank()
  // Number of vertical levels: team.team_size()
  static KOKKOS_INLINE_FUNCTION
  void linear_combination(const ThreadTeam& team,
                          Real a, Real b, const Data& x1, const Data& x2,
                          Data& y) {
    y.linear_combination(team, a, b, x1, x2);
  }

  // Performs horizontal interpolation by applying weights to source data
  // at a specific target column, across a range of vertical levels in that
  // column indicated by the given thread team.
  // Number of vertical levels: team.team_size()
  static KOKKOS_INLINE_FUNCTION
  void interpolate_horizontally(const ThreadTeam& team,
                                const TetralinearInterpWeights& src_weights,
                                const Data& src_data,
                                int tgt_column, Data& tgt_data) {
    tgt_data.interpolate_horizontally(team, src_weights, src_data, tgt_column);
  }

  // Performs linear vertical interpolation from data for a variable with the
  // given index on a given column from a set of source vertical coordinates to
  // a target set. Interpolation is performed in parallel over vertical levels
  // using a thread range defined by the given thread team, using the given
  // vertical interpolator.
  static KOKKOS_INLINE_FUNCTION
  void interpolate_vertically(const ThreadTeam& team,
                              const ekat::LinInterp<Real, Pack::n>& vert_interp,
                              int column_index, int variable_index,
                              const VCoordColumnView& src_col_vcoords,
                              const Data& src_data,
                              const VCoordColumnView& tgt_col_vcoords,
                              Data& tgt_data) {
    // This generic implementation just calls the interpolate_vertically method
    // on the Data type.
    tgt_data.interpolate_vertically(team, vert_interp,
                                    variable_index, column_index,
                                    src_col_vcoords, src_data,
                                    tgt_col_vcoords);
  }
};

} // namespace interpolators
} // namespace scream

#endif // TETRALINEAR_INTERP_TRAITS_HPP
