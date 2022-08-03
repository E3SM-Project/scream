#include <share/io/scream_scorpio_interface.hpp>
#include <share/interp/coarse_grid.hpp>
#include <share/interp/tetralinear_interp.hpp>

#include <unordered_map>

namespace scream {
namespace interpolators {

namespace {

// a teensy tinsy kd-tree
struct KDTree {
  // Construct a kd-tree containing a set of points defined by
  // (longitude, latitude) pairs.
  KDTree(const HostHCoordView& lons, const HostHCoordView& lats) {
  }

  // returns a set of points bounded by the given rectangle whose corners are
  // defined by low and high longitude and latitude coordinates.
  std::vector<int> find_in_rect(Real lon_lo, Real lat_lo,
                                Real lon_hi, Real lat_hi) const {
    std::vector<int> pts;
    return std::move(pts);
  }
};

Real min(Real a, Real b, Real c, Real d) {
  return std::min(a, std::min(b, std::min(c, d)));
}

Real max(Real a, Real b, Real c, Real d) {
  return std::max(a, std::max(b, std::max(c, d)));
}

}

TetralinearInterpWeightMap
compute_tetralinear_interp_weights(const CoarѕeGrid& coarse_grid,
                                   const HostHCoordView& tgt_lons,
                                   const HostHCoordView& tgt_lats) {
  EKAT_ASSERT(tgt_lats.extent(0) == tgt_lons.extent(0));

  // Construct a kd-tree containing the target (lat, lon) points.
  KDTree tree(tgt_lats, tgt_lons);

  HostTetralinearInterpWeightMap host_mapping;

  // Loop over the elements of the coarse grid and find the target points that
  // fall within each one.
  for (int e = 0; e < coarse_grid.num_elems(); ++e) {
    // Create a bounding box that captures all points within the element.
    const ElemData& elem_data = coarse_grid.elem_data[e];
    Real lon_lo = min(coarse_grid.longitudes[elem_data.vertices[0]],
                      coarse_grid.longitudes[elem_data.vertices[1]],
                      coarse_grid.longitudes[elem_data.vertices[2]],
                      coarse_grid.longitudes[elem_data.vertices[3]]);
    Real lat_lo = min(coarse_grid.latitudes[elem_data.vertices[0]],
                      coarse_grid.latitudes[elem_data.vertices[1]],
                      coarse_grid.latitudes[elem_data.vertices[2]],
                      coarse_grid.latitudes[elem_data.vertices[3]]);
    Real lon_hi = max(coarse_grid.longitudes[elem_data.vertices[0]],
                      coarse_grid.longitudes[elem_data.vertices[1]],
                      coarse_grid.longitudes[elem_data.vertices[2]],
                      coarse_grid.longitudes[elem_data.vertices[3]]);
    Real lat_hi = max(coarse_grid.latitudes[elem_data.vertices[0]],
                      coarse_grid.latitudes[elem_data.vertices[1]],
                      coarse_grid.latitudes[elem_data.vertices[2]],
                      coarse_grid.latitudes[elem_data.vertices[3]]);

    // Find the target points in this bounding box.
    auto pts_in_elem = tree.find_in_rect(lon_lo, lat_lo, lon_hi, lat_hi);

    // Compute horizontal interp weights for all points actually within the
    // element.
    for (int pt_index: pts_in_elem) {
      Real lon = tgt_lons[pt_index];
      Real lat = tgt_lats[pt_index];

      // Compute the reference coordinates (a, b) for this point within the
      // subelement. The reference coordinates are in [0, 1] within the
      // subelement defined by the GLL points bounding the (lat, lon) point.
      Real a, b;
      std::tie(a, b) = coarse_grid.ll_to_ref(e, lat, lon);

      // Compute the interpolation weights.
      host_mapping.insert(pt_index, TetralinearInterpWeights{
        elem_data.vertices, // support indices
        { // interpolation weights
          (1.0 - a) * (1.0 - b),
          (1.0 - a) * b,
          a * b,
          a * (1.0 - b)
        }
      });
    }
  }

  // Copy the mapping from host to device.
  TetralinearInterpWeightMap device_mapping;
  kokkos::deep_copy(device_mapping, host_mapping);
  return device_mapping;
}

} // namespace interpolators
} // namespace scream


