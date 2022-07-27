#ifndef SCREAM_INTERPOLATORS_HPP
#define SCREAM_INTERPOLATORS_HPP

#include <share/grid/point_grid.hpp>
#include <share/grid/se_grid.hpp>
#include <share/util/scream_interpolator_traits.hpp>

#include <ekat/mpi/ekat_comm.hpp>

namespace scream {
namespace interpolators {

// Problems encountered with interpolators produce exceptions of this type.
class InterpolationException: public std::exception {
public:
  /// Constructs an exception containing the given descriptive message.
  InterpolationException(const std::string& message) : _message(message) {}

  const char* what() const throw() { return _message.c_str(); }

private:
  std::string _message;
};

// This class implements 4D space-time interpolation from one set of columns at
// fixed latitudes/longitudes to another set. The details of the interpolation
// process are implemented by the traits class
// TetralinearInterpolatorTraits<Data>, which can either be a generic
// implementation relying on methods defined on types, or a specialization for
// a type without those methods.
template <typename Data>
class TetralinearInterpolator {
public:
  using Traits     = TetralinearInterpolatorTraits<Data>;
  using HCoordView = Traits::HCoordView;
  using VCoordView = Traits::VCoordView;

  // Constructs a tetralinear (4D) lat-lon interpolator from data in the given
  // file for the purpose of interpolating field data from the to the points
  // defined by the given set of latitudes and longitudes. The data file should
  // be a NetCDF4 file containing field data defined on a cubed-sphere grid
  // (neXnpY). This constructor throws an exception if any error is encountered.
  TetralinearInterpolator(const std::string& data_file,
                          const HCoordView&  latitudes,
                          const HCoordView&  longitudes) {
    init_from_file_(comm, data_file, latitudes, longitudes);
  }

  TetralinearInterpolator(const TetralinearInterpolator&) = default;
  ~TetralinearInterpolator() = default;

  // Call operator: interpolates field source data from the data file to
  // the given data at the given time and vertical coordinates.
  void operator()(Real time, const VCoordView& vcoords, Data& data) {
    interpolate_(time, vcoords, data);
  }

  // Unsupported operations
  TetralinearInterpolator() = delete;
  TetralinearInterpolator& operator=(const TetralinearInterpolator&) = delete;

private:

  // Reads data from the file with the given name, populating data fields below
  void init_from_file_(const std::string& data_file,
                       const HCoordView& lats, const HCoordView& lons);

  // Interpolate the source data at the given time, placing the result in data.
  void do_time_interpolation_(const std::vector<Real>& times_,
                              const std::vector<Data>& src_data,
                              Real time, Data& data);

  // Interpolate the source data from its own vertical coordinates to the given
  // target vertical coordinates, placing the result in data.
  void do_vertical_interpolation_(const Data& src_data,
                                  const VCoordView& tgt_vcoords,
                                  Data& data);

  // Interpolates source data at the given time, storing the result in data.
  void interpolate_(Real time, const VCoordView& vcoords, Data& data);

  // Data stored at times in a periodic sequence
  std::vector<Real> times_;
  std::vector<Data> data_;

  // Horizontal interpolation weights (fixed in time), encoded in a sparse map
  // representing a linear operator.
  TetralinearInterpWeightMap h_weights_;
};

// This type represents a coarse grid from which data is interpolated.
struct CoarseGrid;

// This function reads coarse grid data from the file with the given name.
std::unique_ptr<CoarseGrid> read_coarse_grid(const std::string& filename);

// This function maps each target lat/lon pair to its corresponding support in
// the given spectral element grid. The support for a lat/lon pair is defined by
// the set of 4 vertices for the quadrilateral cell/subcell bounding that pair
// (see scream_interpolator_traits.hpp for details).
TetralinearInterpWeightMap
compute_tetralinear_interp_weights(const CoarseGrid& coarse_grid,
                                   const HostHCoordView& tgt_lats,
                                   const HostHCoordView& tgt_lons);

} // namespace interpolators

} // namespace scream

#endif // SCREAM_INTERPOLATORS_HPP

#include <share/util/scream_interpolators_impl.hpp>
