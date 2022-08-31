#include <share/io/scream_scorpio_interface.hpp>
#include <share/interp/coarse_grid.hpp>
#include <share/interp/tetralinear_interp.hpp>

#include <unordered_map>

namespace scream {
namespace interpolators {

namespace {

// a teensy tinsy 2D kd-tree
class KdTree {
 public:
  // Constructs a kd-tree containing a set of points defined by
  // (longitude, latitude) pairs stored in the given (host) views.
  KdTree(const HostHCoordView& lons, const HostHCoordView& lats)
    : root_(nullptr), num_points_(0) {
    EKAT_ASSERT(lons.size() == lats.size());
    for (size_t i = 0; i < lons.size(); ++i) {
      insert_(lons(i), lats(i));
    }
  }

  // Finds a set of points bounded by the given rectangle whose corners are
  // defined by low and high longitude and latitude coordinates, storing their
  // indices in the given result vector.
  void find_in_rect(Real lon_lo, Real lat_lo, Real lon_hi, Real lat_hi,
                    std::vector<int>& result) const {
    result.clear();
    Real low[2] = {lon_lo, lat_lo}, high[2] = {lon_hi, lat_hi};
    find_in_rect_(low, high, root_, 0, result);
  }

 private:

  struct Node {
    int   index;
    Real  coords[2]; // (lon, lat)
    Node *left;
    Node *right;

    Node(int index_, Real coords_[2]):
      index(index_), left(nullptr), right(nullptr) {
      coords[0] = coords_[0];
      coords[1] = coords_[1];
    }
  };

  Node *root_;
  int num_points_;

  void insert_(Real lon, Real lat) {
    Real coords[2] = {lon, lat};
    insert_(coords, root_, 0);
  }

  void insert_(Real coords[2], Node*& node, int dim) {
    if (node == nullptr) {
      node = new Node(num_points_++, coords);
    } else if(coords[dim] < node->coords[dim]) {
      insert_(coords, node->left, 1-dim);
    } else {
      insert_(coords, node->right, 1-dim);
    }
  }

  void find_in_rect_(Real low[2], Real high[2], Node* n, int dim,
                     std::vector<int>& result) const {
    if(n != nullptr) {
      if(low[0] <= n->coords[0] && high[0] >= n->coords[0] &&
         low[1] <= n->coords[1] && high[1] >= n->coords[1]) { // inside bounds
        result.push_back(n->index);
      }

      if(low[dim] <= n->coords[dim]) {
        find_in_rect_(low, high, n->left, 1-dim, result);
      }
      if(high[dim] >= n->coords[dim]) {
        find_in_rect_(low, high, n->right, 1-dim, result);
      }
    }
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

  // Construct a kd-tree containing the target (lon, lat) points.
  KdTree tree(tgt_lons, tgt_lats);

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
    std::vector<int> pts_in_elem;
    tree.find_in_rect(lon_lo, lat_lo, lon_hi, lat_hi, pts_in_elem);

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


