#ifndef SCREAM_INTERPOLATORS_HPP
#define SCREAM_INTERPOLATORS_HPP

#include <share/io/scream_scorpio_interface.hpp>
#include <share/util/scream_interpolator_traits.hpp>

namespace scream {

// This class implements 4D space-time interpolation from one set of columns at
// fixed latitudes/longitudes to another set. The details of the interpolation
// process are implemented by the traits class
// LatLonLinInterpolatorTraits<Data>, which can either be a generic
// implementation relying on methods defined on types, or a specialization for
// a type without those methods.
template <typename Data>
class LatLonLinInterpolator {
public:
  using Traits = LatLonLinInterpolatorTraits<Data>;

  // Constructs a linear lat-lon interpolator from data in the given file for
  // the purpose of interpolating field data from the to the points defined by the given
  // set of latitudes and longitudes. The data file should be a NetCDF4 file
  // containing field data defined on a cubed-sphere grid (neXnpY).
  LatLonLinInterpolator(const std::string& data_file,
                        const HCoordView&  latitudes,
                        const HCoordView&  longitudes) {
    init_from_file_(data_file, latitudes, longitudes);
  }

  LatLonLinInterpolator(const LatLonLinInterpolator&) = default;
  ~LatLonLinInterpolator() = default;

  // Call operator: interpolates field source data from the data file to
  // the given data at the given time and vertical coordinates.
  void operator()(Real time, const VCoordView& vcoords, Data& data) {
    interpolate_(time, vcoords, data);
  }

  // Unsupported operations
  LatLonLinInterpolator() = delete;
  LatLonLinInterpolator& operator=(const LatLonLinInterpolator&) = delete;

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
  LatLonLinColumnWeightMap h_weights_;
};

// This function reads data from the given data file and produces an unordered
// map whose keys are indices into the lat/lon arrays and whose values are
// LatLonLinColumnWeights (see scream_interpolator_traits.hpp for details).
LatLonLinColumnWeightMap
compute_latlonlin_column_weights(const std::string& data_file,
                                 const HCoordView& latitudes,
                                 const HCoordView& longitudes);

} // namespace scream

#endif // SCREAM_INTERPOLATORS_HPP

#include <share/util/scream_interpolators_impl.hpp>
