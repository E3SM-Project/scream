#ifndef SCREAM_INTERPOLATORS_HPP
#define SCREAM_INTERPOLATORS_HPP

#include <share/util/scream_interpolator_traits.hpp>

#include <ekat/mpi/ekat_comm.hpp>

namespace scream {
namespace interpolators {

// This type represents a source data file used by interpolators, accessed using
// the Scorpio/NetCDF interface. It provides only basic functionality needed
// by interpolators.
class SourceDataFile final {
  int pio_id_;   // Scorpio I/O system ID
  int file_id_;  // Scorpio file ID
 public:

  // Opens the file with the given name on rank 0 of the given communicator.
  SourceDataFile(const ekat::Comm& comm, const std::string& filename);

  // Closes the file.
  ~SourceDataFile();

  // Fetches times from the source file into the given array.
  void get_times(std::vector<Real>& times) const;

  // Fetches n values from the 1d array variable with the given name, placing
  // them into the given array.
  void get_array(const std::string& name, std::vector<Real>& values) const;

  // Disallowed stuff
  SourceDataFile() = delete;
  SourceDataFile(const SourceDataFile&) = delete;
  SourceDataFile& operator=(const SourceDataFile&) = delete;
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

  // Constructs a linear lat-lon interpolator from data in the given file for
  // the purpose of interpolating field data from the to the points defined by the given
  // set of latitudes and longitudes. The data file should be a NetCDF4 file
  // containing field data defined on a cubed-sphere grid (neXnpY).
  TetralinearInterpolator(const ekat::Comm& comm,
                          const std::string& data_file,
                          const HCoordView&  latitudes,
                          const HCoordView&  longitudes): comm_(comm) {
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

  // MPI communicator
  ekat::Comm comm_;

  // Data stored at times in a periodic sequence
  std::vector<Real> times_;
  std::vector<Data> data_;

  // Horizontal interpolation weights (fixed in time), encoded in a sparse map
  // representing a linear operator.
  TetralinearColumnWeightMap h_weights_;
};

// This type defines a bounding box in (latitude, longitude) coordinates.
struct BoundingBox {
  Real min_latitude, max_latitude;
  Real min_longitude, max_longitude;
};

// This function reads interpolation data from the data file with the given
// name, populating the given vectors on each MPI process according to the
// given bounding box.
void read_source_coordinates(const ekat::Comm& comm,
                             const SourceDataFile& data_file,
                             const BoundingBox& bbox,
                             std::vector<Real>& src_times,
                             std::vector<Real>& src_lats,
                             std::vector<Real>& src_lons);

// This function produces an unordered map whose keys are indices into the
// lat/lon arrays and whose values are TetralinearColumnWeights
// (see scream_interpolator_traits.hpp for details).
TetralinearColumnWeightMap
compute_tetralinear_column_weights(const std::vector<Real>& src_lats,
                                   const std::vector<Real>& src_lons,
                                   const HCoordView& latitudes,
                                   const HCoordView& longitudes);

} // namespace interpolators

} // namespace scream

#endif // SCREAM_INTERPOLATORS_HPP

#include <share/util/scream_interpolators_impl.hpp>
