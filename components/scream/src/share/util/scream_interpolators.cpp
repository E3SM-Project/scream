#include <share/io/scream_scorpio_interface.hpp>
#include <share/util/scream_interpolators.hpp>

#include <unordered_map>

namespace scream {
namespace interpolators {

// This type represents a coarse grid from which data is interpolated. This grid
// consists of a set of quadrilateral elements and their vertices, as well as
// any gauss-legendre-lobatto points (if this grid is interpreted as a spectral
// element grid). It is loaded into memory on all ranks of the atm MPI
// communicator so a correspondence can be established between the elements of
// the coarse grid and target points to which data is interpolated.
struct CoarseGrid {
  // Latitudes and longitudes of element vertices, sorted in ascending order.
  // This includes coordinates of ALL gauss-legendre-lobatto points, not just
  // those belonging to element vertices.
  std::vector<Real> latitudes, longitudes;

  // This type contains the indices of latitudes and longitudes for each of an
  // element's (shared) vertices (corners) within the above latitudes and
  // longitudes vectors, and indices identifying degrees of freedom for field
  // data.
  struct ElemData {
    int latitudes[2];      // lower and upper latitude indices
    int longitudes[2];     // lower and upper longitude indices
    int dof_indices[4][4]; // dof indices for each gll point in the element
                           // (up to np=4; only the first np slots are used)
  };

  // This mapping associates an element with its data.
  std::vector<ElemData> elem_data;

  // number of gauss-legendre lobatto points within an element in each spatial
  // dimension (2 points indicates there are no such points interior to each
  // element, so a latitude-longitude coarse grid has 2 "gauss-legendre-lobatto"
  // points)
  int num_gll_pts;

  // Constructs a grid from a given number of elements and vertices, with the
  // option to specify a number of gauss-legendre-lobatto points for spectral
  // element grids. Here we distinguish between a vertex, which identifies a
  // corner shared between elements, and the set of all gauss-legendre-lobatto
  // points including vertices and points interior to elements.
  CoarseGrid(int num_elems, int num_vertices, int num_gll_points = 2)
    : latitudes(num_vertices*num_gll_points/2),
      longitudes(num_vertices*num_gll_points/2),
      elem_data(num_elems), num_gll_pts(num_gll_points) {
  }

  // destructor
  ~CoarseGrid() = default;

  // number of elements in the grid
  int num_elems()    const { return elem_data.extent(0); }

  // number of (shared) element vertices in the grid
  int num_vertices() const { return 2*latitudes.extent(0)/num_gll_pts; }

  // disallowed stuff
  CoarseGrid() = delete;
  CoarseGrid(const CoarseGrid&) = delete;
  CoarseGrid& operator=(const CoarseGrid&) = delete;
};

std::unique_ptr<CoarseGrid> read_coarse_grid(const std::string& filename) {
  comm_self = ekat::Comm(MPI_COMM_SELF);

  int ne = scorpio::get_int_attribute_c2f(filename.c_str(), "ne");
  int np = scorpio::get_int_attribute_c2f(filename.c_str(), "np");
  int num_elem = 6*ne*ne;
  int nlev = scorpio::get_dimlen_c2f(filename.c_str(), "lev");

  std::unique_ptr<CoarseGrid> grid(new CoarseGrid(num_elem, num_vertices, np));

  // FIXME: What do we expect in the file???

  return grid;
}

namespace {

// Given a number of gauss-legendre-lobatto points, compute their locations.
void compute_gll_pts(int num_gll_pts, Real *pts) {
  EKAT_ASSERT((num_gll_pts == 2) || (num_gll_pts == 4));
  if (num_gll_pts == 2) {
    pts[0] = 0.0;
    pts[1] = 1.0;
  } else { // (num_gll_pts == 4)
    pts[0] = 0.0;
    // FIXME
    pts[3] = 1.0;
  }
}

// Given a set of element vertex indices, a number of gauss-legendre-lobatto
// points, and a set of and a pair of element reference coordinates, generate
// the corresponding interpolation weights.
TetralinearInterpWeight ref_to_weights(const CoarseGrid::ElemData& elem_data,
                                       Real* gll_pts, int num_gll_pts,
                                       Real a, Real b) {
  EKAT_ASSERT((num_gll_pts == 2) || (num_gll_pts == 4));
  if (num_gll_pts == 2) {
    Real W0_lat = (1.0 - (a - elem_data.latitudes[0]));
    Real W1_lat = elem_data.latitudes[1] - a;
    Real W0_lon = (1.0 - (b - elem_data.longitudes[0]));
    Real W1_lon = elem_data.longitudes[1] - b;
    return TetralinearInterpWeight{
      { // indices
        elem_data.dof_indices[0][0],
        elem_data.dof_indices[0][1],
        elem_data.dof_indices[1][0],
        elem_data.dof_indices[1][1]
      },
      { // weights
        W0_lat * W0_lon,
        W0_lat * W1_lon,
        W1_lat * W0_lon,
        W1_lat * W1_lon
      }
    };
  } else { // (num_gll_pts == 4)
    // Figure out the subelement containing (a, b).
    // FIXME
  }
}

}

TetralinearInterpWeightMap
compute_tetralinear_interp_weights(const CoarѕeGrid& coarse_grid,
                                   const HostHCoordView& tgt_lats,
                                   const HostHCoordView& tgt_lons,
                                   std::function<void(Real[2], Real[2],
                                                      Real, Real,
                                                      Real&, Real&)> lat_lon_to_ref) {
  EKAT_ASSERT(tgt_lats.extent(0) == tgt_lons.extent(0));

  // Precompute gauss-legendre-lobatto points in element reference coordinates.
  Real gll_pts[4];
  compute_gll_pts(coarse_grid.num_gll_pts, gll_pts);

  // Generate a mapping from coarse latitude/longitude indices to the index of
  // the element containing them.
  std::unordered_map<std::pair<int, int>, int> ll_to_elem;
  for (int e = 0; e < coarse_grid.num_elems(); ++e) {
    int lo_lat = coarѕe_grid.elem_data[e].latitudes[0];
    int lo_lon = coarѕe_grid.elem_data[e].longitudes[0];
    for (int igp = 0; igp < coarse_grid.num_gll_pts; ++igp) {
      int ilat = lo_lat + igp;
      for (int jgp = 0; jgp < coarse_grid.num_gll_pts; ++jgp) {
        int jlon = lo_lon + jgp;
        auto key = std::pair<int, int>(ilat, jlon);
        ll_to_elem[key] = e;
      }
    }
  }

  // Now map the target points to their coarse elements and compute
  // interpolation weights.
  HostTetralinearInterpWeightMap host_mapping;
  int num_tgt_pts = tgt_lats.extent(0);
  for (int i = 0; i < num_tgt_pts; ++i) {
    // Find the coarse latitude/longitude point closest to this one using
    // binary search.
    Real tgt_lat = tgt_lats(i), tgt_lon = tgt_lons(i);
    auto lat_iter = std::lower_bound(coarse_grid.latitudes.begin(),
                                     coarse_grid.latitudes.end(),
                                     tgt_lat);
    int coarse_lat = static_cast<int>(
      std::distance(coarse_grid.latitudes.begin(), lat_iter);
    auto lon_iter = std::lower_bound(coarse_grid.longitudes.begin(),
                                     coarse_grid.longitudes.end(),
                                     tgt_lon);
    int coarse_lon = static_cast<int>(
      std::distance(coarse_grid.longitudes.begin(), lon_iter);

    // Get the element for this coarse lat/lon pair.
    auto key = std::pair<int, int>(coarse_lat, coarse_lon);
    int elem = ll_to_elem[key];

    // Extract the coordinates of the element's corners and compute the
    // reference coordinates (a, b) within the element for the target point.
    const CoarseGrid::ElemData& elem_data = coarse_grid.elem_data[elem];
    Real corner_lats[2] = {
      coarse_grid.latitudes[elem_data.latitudes[0]],
      coarse_grid.latitudes[elem_data.latitudes[1]]
    };
    Real corner_lons[2] = {
      coarse_grid.longitudes[elem_data.longitudes[0]],
      coarse_grid.longitudes[elem_data.longitudes[1]]
    };
    Real a, b;
    lat_lon_to_ref(corner_lats, corner_lons, tgt_lat, tgt_lon, a, b);

    // Compute interpolation indices and weights from the reference coordinates.
    host_mapping[i] = ref_to_weights(elem_data,
                                     gll_pts, coarse_grid.num_gll_pts,
                                     a, b);
  }

  // Copy the mapping from host to device.
  TetralinearInterpWeightMap device_mapping;
  kokkos::deep_copy(device_mapping, host_mapping);
  return host_mapping;
}

} // namespace interpolators
} // namespace scream


