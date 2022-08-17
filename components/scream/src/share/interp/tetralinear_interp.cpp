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
    return pts;
  }
};

Real min(Real a, Real b, Real c, Real d) {
  return std::min(a, std::min(b, std::min(c, d)));
}

Real max(Real a, Real b, Real c, Real d) {
  return std::max(a, std::max(b, std::max(c, d)));
}

}

void
compute_tetralinear_interp_weights(const CoarseGrid& coarse_grid,
                                   const HostHCoordView& tgt_lons,
                                   const HostHCoordView& tgt_lats,
                                   TetralinearInterpWeightMap& weights,
                                   std::vector<int>& global_vertex_indices) {
  EKAT_ASSERT(tgt_lats.extent(0) == tgt_lons.extent(0));

  // Construct a kd-tree containing the target (lat, lon) points.
  KDTree tree(tgt_lats, tgt_lons);

  HostTetralinearInterpWeightMap host_mapping;

  // Loop over the elements of the coarse grid and find the target points that
  // fall within each one.
  std::set<int> unique_vertices;
  for (int e = 0; e < coarse_grid.num_elems(); ++e) {
    const CoarseGrid::Element& elem = coarse_grid.elements[e];

    // Add the element's vertices to our unique list of global vertices.
    for (int v = 0; v < 4; ++v) {
      unique_vertices.insert(elem.vertices[v]);
    }

    // Create a bounding box that captures all points within the element.
    Real lon_lo = min(coarse_grid.longitudes[elem.vertices[0]],
                      coarse_grid.longitudes[elem.vertices[1]],
                      coarse_grid.longitudes[elem.vertices[2]],
                      coarse_grid.longitudes[elem.vertices[3]]);
    Real lat_lo = min(coarse_grid.latitudes[elem.vertices[0]],
                      coarse_grid.latitudes[elem.vertices[1]],
                      coarse_grid.latitudes[elem.vertices[2]],
                      coarse_grid.latitudes[elem.vertices[3]]);
    Real lon_hi = max(coarse_grid.longitudes[elem.vertices[0]],
                      coarse_grid.longitudes[elem.vertices[1]],
                      coarse_grid.longitudes[elem.vertices[2]],
                      coarse_grid.longitudes[elem.vertices[3]]);
    Real lat_hi = max(coarse_grid.latitudes[elem.vertices[0]],
                      coarse_grid.latitudes[elem.vertices[1]],
                      coarse_grid.latitudes[elem.vertices[2]],
                      coarse_grid.latitudes[elem.vertices[3]]);

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
      auto coords = coarse_grid.ll_to_ref(e, lat, lon);
      Real a = coords.first, b = coords.second;

      // Compute the interpolation weights.
      host_mapping.insert(pt_index, TetralinearInterpWeights{
        { // support indices
          elem.vertices[0],
          elem.vertices[1],
          elem.vertices[2],
          elem.vertices[3]
        },
        { // interpolation weights
          (1.0 - a) * (1.0 - b),
          (1.0 - a) * b,
          a * b,
          a * (1.0 - b)
        }
      });
    }
  }

  // Now reindex the vertices within the weights to reflect the fact that we
  // only load data on the relevant elements.
  global_vertex_indices.resize(unique_vertices.size());
  std::copy(unique_vertices.begin(), unique_vertices.end(),
            global_vertex_indices.begin());
  std::unordered_map<int, int> g2l;
  for (size_t i = 0; i < global_vertex_indices.size(); ++i) {
    g2l[global_vertex_indices[i]] = i;
  }
  for (uint64_t i = 0; i < host_mapping.capacity(); ++i) {
    if (host_mapping.valid_at(i)) {
      for (int v = 0; v < 4; ++v) {
        auto& mapping = host_mapping.value_at(i);
        mapping.indices[v] = g2l[mapping.indices[v]];
      }
    }
  }

  // Copy the mapping from host to device.
  Kokkos::deep_copy(weights, host_mapping);
}

} // namespace interpolators
} // namespace scream


