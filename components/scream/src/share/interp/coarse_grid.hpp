#ifndef SCREAM_INTERP_COARSE_GRID_HPP
#define SCREAM_INTERP_COARSE_GRID_HPP

#include <share/scream_types.hpp>

#include <vector>

namespace scream {
namespace interpolators {

// This type represents a coarse quadrilateral grid from which data is
// interpolated. This grid consists of a set of quad elements and their
// vertices, along with latitude and longitude coordinates for those vertices.
// It is loaded into memory on all ranks of the atm MPI communicator so a
// correspondence can be established between the elements of the coarse grid and
// target points to which data is interpolated.
//
// NOTE: this class only works on the host, NOT the device.
struct CoarseGrid {
  // Latitudes and longitudes of element vertices.
  std::vector<Real> latitudes, longitudes;

  // This type contains the global indices of latitude and longitude
  // coordinates for each of the element's vertices (corners).
  //
  // In general, element corner coordinate indices are (moving clockwise
  // around the quadrilateral element e, starting at the lower left corner):
  //
  // 3 +-------+ 2
  //   |       |
  //   |   e   |
  //   |       |
  // 0 +-------+ 1
  //
  struct Element {
    int vertices[4];
  };

  // vector of elements mapping element global IDs to vertices
  std::vector<Element> elements;

  // number of elements on a side of a cubed-sphere panel
  int ne;

  // Reads a coarse grid from the file with the given name.
  explicit CoarseGrid(const std::string& filename);

  // destructor
  ~CoarseGrid() = default;

  // total number of elements in the grid
  int num_elems() const { return 6*ne*ne; }

  // Given (longitude, latitude) coordinates for a point that falls within an
  // element e, return that point's reference coordinates (a, b) within the
  // element.
  std::pair<Real, Real> ll_to_ref(int e, Real lon, Real lat,
                                  Real tol = 1e-12, int max_iter = 10) const;

  // disallowed stuff
  CoarseGrid() = delete;
  CoarseGrid(const CoarseGrid&) = delete;
  CoarseGrid& operator=(const CoarseGrid&) = delete;
};

} // namespace interpolators
} // namespace scream

#endif // SCREAM_INTERP_COARSE_GRID_HPP
