#include <share/io/scream_scorpio_interface.hpp>
#include <share/interp/coarse_se_grid.hpp>
#include <share/interp/tetralinear_interp.hpp>

#include <unordered_map>

namespace scream {
namespace interpolators {

TetralinearInterpWeightMap
compute_tetralinear_interp_weights(const CoarѕeSEGrid& coarse_grid,
                                   const HostHCoordView& tgt_lats,
                                   const HostHCoordView& tgt_lons) {
  EKAT_ASSERT(tgt_lats.extent(0) == tgt_lons.extent(0));

  // Construct a kd-tree containing the target (lat, lon) points.
  KDTree tree(tgt_lats, tgt_lons);

  HostTetralinearInterpWeightMap host_mapping;

  // Loop over the elements of the coarse grid and find the target points that
  // fall within each one.
  for (int e = 0; e < coarse_grid.num_elems(); ++e) {
    // Create a bounding box that captures all points within the element.
    Real lat_lo = min(coarse_grid.elem_data.latitudes[0][0],
                      coarse_grid.elem_data.latitudes[0][NP]);
    Real lat_hi = max(coarse_grid.elem_data.latitudes[0][0],
                      coarse_grid.elem_data.latitudes[0][NP]);
    Real lon_lo = min(coarse_grid.elem_data.latitudes[0][0],
                      coarse_grid.elem_data.latitudes[NP][0]);
    Real lon_hi = max(coarse_grid.elem_data.latitudes[0][0],
                      coarse_grid.elem_data.latitudes[NP][0]);

    // Find the target points in this bounding box.
    auto pts_in_elem = tree.find_in_rect(lat_lo, lon_lo, lat_hi, lon_hi);

    // Compute horizontal interp weights for all points actually within the
    // element.
    for (int pt_index: pts_in_elem) {
      Real lat = tgt_lats[pt_index];
      Real lon = tgt_lons[pt_index];

      // Find the subelement (ei, ej) within e that contains the given
      // (lat, lon) pair. A subelement index of -1 indicates that the element e
      // doesn't actually contain this point.
      int ei, ej;
      std::tie(ei, ej) = coarse_grid.subelement(e, lat, lon);
      if ((ei != -1) && (ej != -1)) {
        // Compute the reference coordinates (a, b) for this point within the
        // subelement. The reference coordinates are in [0, 1] within the
        // subelement defined by the GLL points bounding the (lat, lon) point.
        Real a, b;
        std::tie(a, b) = coarse_grid.ll_to_ref(e, ei, ej, lat, lon);

        // Compute the interpolation weights.
        host_mapping.insert(pt_index, TetralinearInterpWeights{
          { // support indices
            coarse_grid.elem_data[e].gll_pts[ei][ej],
            coarse_grid.elem_data[e].gll_pts[ei][ej+1],
            coarse_grid.elem_data[e].gll_pts[ei+1][ej+1],
            coarse_grid.elem_data[e].gll_pts[ei+1][ej],
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
  }

  // Copy the mapping from host to device.
  TetralinearInterpWeightMap device_mapping;
  kokkos::deep_copy(device_mapping, host_mapping);
  return device_mapping;
}

} // namespace interpolators
} // namespace scream


