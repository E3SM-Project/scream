#ifndef TETRALINEAR_INTERP_HPP
#define TETRALINEAR_INTERP_HPP

#include <share/interp/coarse_grid.hpp>
#include <share/interp/tetralinear_interp_traits.hpp>

namespace scream {
namespace interpolators {

// Problems encountered with interpolators produce exceptions of this type.
class InterpException: public std::exception {
public:
  /// Constructs an exception containing the given descriptive message.
  InterpException(const std::string& message) : _message(message) {}

  const char* what() const throw() { return _message.c_str(); }

private:
  std::string _message;
};

// This class implements 4D space-time interpolation from one set of columns at
// fixed latitudes/longitudes to another set. The details of the interpolation
// process are implemented by the traits class TetralinearInterpTraits<DataSet>,
// which can either be a generic implementation relying on methods defined on
// types, or a specialization for a type without those methods.
template <typename DataSet>
class TetralinearInterp {
public:
  using Traits     = TetralinearInterpTraits<DataSet>;
  using HCoordView = typename Traits::HCoordView;
  using VCoordView = typename Traits::VCoordView;

  // Constructs a tetralinear (4D) lat-lon-vert-time interpolator from data in
  // the given file for the purpose of interpolating field data from a coarse
  // cubed-sphere grid to a set of target points defined by the given latitudes
  // and longitudes. The data file is a NetCDF4 file containing field data
  // defined on a cubed-sphere grid with 2 gauss-legendre-lobatto points per
  // element per horizontal dimension (neXnp2). This constructor throws an
  // exception if any error is encountered.
  TetralinearInterp(const std::string& data_file,
                    const HCoordView&  latitudes,
                    const HCoordView&  longitudes) {
    init_from_file_(data_file, latitudes, longitudes);
  }

  TetralinearInterp(const TetralinearInterp&) = default;
  ~TetralinearInterp() = default;

  // Call operator: interpolates field source data from the data file to
  // the given data at the given time, with source and target vertical
  // coordinates provided by an atmospheric state.
  void operator()(Real time,
                  const VCoordView& src_vcoords,
                  const VCoordView& tgt_vcoords,
                  DataSet& tgt_data) const {
    interpolate_(time, src_vcoords, tgt_vcoords, tgt_data);
  }

  // Unsupported operations
  TetralinearInterp() = delete;
  TetralinearInterp& operator=(const TetralinearInterp&) = delete;

private:

  // Reads data from the file with the given name, populating data fields below
  void init_from_file_(const std::string& data_file,
                       const HCoordView& lats, const HCoordView& lons);

  // Interpolate the source data at the given time, placing the result in data.
  void do_time_interpolation_(Real time, DataSet& data) const;

  // Interpolate source data horizontally by applying precomputed interpolation
  // weights to src_data, placing the result in tgt_data.
  void do_horizontal_interpolation_(const DataSet& src_data,
                                    DataSet& tgt_data) const;

  // Computes vertical coordinates for the given dataset across all columns.
  void compute_vertical_coords_(const DataSet& data, VCoordView& vcoords) const;

  // Interpolate the source data from the source vertical coordinates to the
  // given target vertical coordinates, placing the result in data.
  void do_vertical_interpolation_(const VCoordView& src_vcoords,
                                  const DataSet& src_data,
                                  const VCoordView& tgt_vcoords,
                                  DataSet& tgt_data) const;

  // Interpolates source data at the given time from the given source vertical
  // coordinate profile to the target profile, storing the result in data.
  void interpolate_(Real time,
                    const VCoordView& tgt_vcoords,
                    DataSet& data) const;

  // Data stored at times in a periodic sequence
  std::vector<Real>    times_;
  std::vector<DataSet> data_;

  // Horizontal interpolation weights (fixed in time), encoded in a sparse map
  // representing a linear operator.
  TetralinearInterpWeightMap weights_;
};

// This function maps each target lat/lon pair to its corresponding support in
// the given spectral element grid. The support for a lat/lon pair is defined by
// the set of 4 vertices for the quadrilateral cell/subcell bounding that pair
// (see tetralinear_interp_traits.hpp for details). The mapping is stored in
// weights, and the list of relevant global vertex indices for the grid is
// stored in global_vertex_indices.
void
compute_tetralinear_interp_weights(const CoarseGrid& coarse_grid,
                                   const HostHCoordView& tgt_lons,
                                   const HostHCoordView& tgt_lats,
                                   TetralinearInterpWeightMap& weights,
                                   std::vector<int>& global_vertex_indices);

} // namespace interpolators

} // namespace scream

#endif // TETRALINEAR_INTERP_HPP

#include <share/interp/tetralinear_interp_impl.hpp>
