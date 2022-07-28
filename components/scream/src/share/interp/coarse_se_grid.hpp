#ifndef SCREAM_INTERP_COARSE_SE_GRID_HPP
#define SCREAM_INTERP_COARSE_SE_GRID_HPP

#include <vector>

namespace scream {
namespace interpolators {

// This type represents a coarse spectral element grid from which data is
// interpolated. This grid consists of a set of quadrilateral elements and their
// vertices, as well as any gauss-legendre-lobatto (GLL) points (if this grid is
// interpreted as a spectral element grid). It is loaded into memory on all
// ranks of the atm MPI communicator so a correspondence can be established
// between the elements of the coarse grid and target points to which data is
// interpolated.
//
// NOTE: this class only works on the host, NOT the device.
struct CoarseSEGrid {
  // maximum allowable value of np
  constexpr int np_max = 4;

  // Latitudes and longitudes of all GLL points.
  std::vector<Real> latitudes, longitudes;

  // This type contains the indices of latitudes and longitudes for each of the
  // element's GLL points. Here are the GLL points for an element in an np=4
  // grid:
  //
  // (3,0)               (3,3)
  //   +-----+-------+-----+
  //   |     |       |     |
  //   |     |       |     |
  //   +-----+-------+-----+
  //   |     |       |     |
  //   |     |       |     |
  //   |     |       |     |
  //   +-----+-------+-----+
  //   |     |       |     |
  //   |     |       |     |
  //   +-----+-------+-----+
  // (0,0)               (0,3)
  //
  // In general, element corner coordinate indices are (moving clockwise
  // around the quadrilateral element e, starting at the lower left corner):
  //
  // 4 +-------+ 3
  //   |       |
  //   |   e   |
  //   |       |
  // 1 +-------+ 2
  //
  // 1. e.gll_pts[0][0]
  // 2. e.gll_pts[0][np-1]
  // 3. e.gll_pts[np-1][np-1]
  // 4. e.gll_pts[np-1][0]
  struct ElemData {
    int gll_pts[np_max][np_max];
  };

  // This mapping associates an element with its GLL point indices.
  std::vector<ElemData> elem_data;

  // number of elements on a side of a cubed-sphere panel
  int ne;

  // number of gauss-legendre lobatto points per element per dimension
  int np;

  // Constructs a spectral element grid with the given specification.
  CoarseSEGrid(int ne_, int np_);

  // destructor
  ~CoarseSEGrid() = default;

  // total number of elements in the grid
  int num_elems() const { return 6*ne*ne; }

  // total number of GLL points in the grid
  int num_gll_points() const { return 6*ne*ne*(np-1)*(np-1) + 2; }

  // Returns a 2-tuple that identifies the "subelement" within the element e,
  // in which the given latitude/longitude point appears. (-1, -1) indicates
  // that the given point falls outside the element e.
  //
  // For np=2, an element contains only one subelement (0, 0): itself. Here are
  // the subelements for a grid with np=4:
  //
  //   +-----+-------+-----+
  //   | 2,0 |  2,1  | 2,2 |
  //   |     |       |     |
  //   +-----+-------+-----+
  //   |     |       |     |
  //   | 1,0 |  1,1  | 1,2 |
  //   |     |       |     |
  //   +-----+-------+-----+
  //   | 0,0 |  0,1  | 0,2 |
  //   |     |       |     |
  //   +-----+-------+-----+
  std::tuple<int, int> subelement(int e, Real lat, Real lon) const {
    const auto& elem = elem_data[e];
    for (int i = 0; i < np-1; ++i) {
      for (int j = 0; j < np-1; ++j) {
        Real lat_lo = min(latitudes[elem.gll_pts[i][j]],
                          latitudes[elem.gll_pts[i][j+1]]);
        Real lat_hi = max(latitudes[elem.gll_pts[i+1][j]];
                          latitudes[elem.gll_pts[i+1][j+1]]);
        Real lon_lo = min(longitudes[elem.gll_pts[i][j]],
                          longitudes[elem.gll_pts[i+1][j]]);
        Real lon_hi = max(longitudes[elem.gll_pts[i][j+1]],
                          longitudes[elem.gll_pts[i+1][j+1]]);
        if ((lat >= lat_lo) && (lat < lat_hi) &&
            (lon >= lon_lo) && (lon < lon_hi)) {
          return std::tuple<int, int>(i, j);
        }
      }
    }
    return std::tuple<int, int>(-1, -1);
  }

  // disallowed stuff
  CoarseSEGrid() = delete;
  CoarseSEGrid(const CoarseSEGrid&) = delete;
  CoarseSEGrid& operator=(const CoarseSEGrid&) = delete;
};

} // namespace interpolators
} // namespace scream

#endif // SCREAM_INTERP_COARSE_SE_GRID_HPP
