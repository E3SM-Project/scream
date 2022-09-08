#include <catch2/catch.hpp>

#include <share/interp/coarse_grid.hpp>
#include <share/interp/tetralinear_interp.hpp>

namespace {

using scream::Real;
using namespace scream::interpolators;

TEST_CASE("coarse_grid", "") {

  // Read an ne4 grid from a file.
  CoarseGrid grid("ne4.g");
  REQUIRE(grid.ne == 4);
  int num_elems = 6*4*4;
  REQUIRE(grid.num_elems() == int(num_elems));
  REQUIRE(grid.elements.size() == size_t(num_elems));
  int num_vertices = num_elems*num_elems + 2;
  REQUIRE(grid.latitudes.size() == size_t(num_vertices));
  REQUIRE(grid.longitudes.size() == size_t(num_vertices));

  // Test the global-to-local coordinate mapping in each element.
  for (int e = 0; e < num_elems; ++e) {
    Real lons[4], lats[4];
    for (int v = 0; v < 4; ++v) {
      lons[v] = grid.longitudes[grid.elements[e].vertices[v]];
      lats[v] = grid.latitudes[grid.elements[e].vertices[v]];
    }
    // test element corners
    auto ab = grid.ll_to_ref(e, lons[0], lats[0]);
    REQUIRE(std::abs(ab.first - 0.0) < 1e-12);
    REQUIRE(std::abs(ab.second - 0.0) < 1e-12);

    ab = grid.ll_to_ref(e, lons[1], lats[1]);
    REQUIRE(std::abs(ab.first - 1.0) < 1e-12);
    REQUIRE(std::abs(ab.second - 0.0) < 1e-12);

    ab = grid.ll_to_ref(e, lons[2], lats[2]);
    REQUIRE(std::abs(ab.first - 1.0) < 1e-12);
    REQUIRE(std::abs(ab.second - 1.0) < 1e-12);

    ab = grid.ll_to_ref(e, lons[3], lats[3]);
    REQUIRE(std::abs(ab.first - 0.0) < 1e-12);
    REQUIRE(std::abs(ab.second - 1.0) < 1e-12);

    // test center
    Real center_lon = 0.25 * (lons[0] + lons[1] + lons[2] + lons[3]);
    Real center_lat = 0.25 * (lats[0] + lats[1] + lats[2] + lats[3]);
    ab = grid.ll_to_ref(e, center_lon, center_lat);
    REQUIRE(std::abs(ab.first - 0.5) < 1e-12);
    REQUIRE(std::abs(ab.second - 0.5) < 1e-12);
  }
}

// This data type stores a scalar field--specifically, the "cosine bell"
// initial condition for the linear advection equation problem proposed in
//
// D. L. Williamson, J. B. Drake, J. J. Hack, R. Jakob, and P. N. Swarztrauber.
// A standard test set for numerical approximations to the shallow water
// equations in spherical geometry. Journal of Computational Physics,
// 102:211–224, 1992.
struct CosineBell {
  // Solution data (col, level)
  KokkosTypes::view_2d<Pack> data;

  // member functions used by traits (see TetralinearInterpTraits for details)

  static CosineBell allocate(int num_cols, int num_levels) {
    return CosineBell{
      KokkosTypes::view_2d<Pack>("solution", num_cols, num_levels)
    };
  }

  void read_from_file(const std::string& filename,
                      int time_index,
                      const std::vector<int>& columns) {
    // Instead of reading data from a file, we just populate the ne4 grid
    // "manually" with the cosine bell solution
    //
    //               {(h0/2)(1 + cos(pi*r/R)), r <  R
    // h(lon, lat) = {
    //               {          0            , r >= R
    //
    // where h0(z) = z and r is the great circle distance between (lon, lat)
    // and the "center point" (lon_c, lat_c) = (3*pi/2, 0):
    //
    // r = arccos[sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c)]
    //
    // The z coordinate is expressed in meters and runs from 0 to 1000, meaning
    // that 0 < h0(z) < 1000.

    // Read horizontal coordinate data from the file.
    CoarseGrid grid(filename);
    int num_cols = int(grid.latitudes.size());

    // We select 128 equally-spaced vertical levels.
    int num_levels = 128;
    Real dz = 1000.0 / num_levels;

    // The center position of the cosine bell moves along the equator, making
    // a complete revolution over 12 months.
    Real omega = 2*M_PI / 12; // longitudinal velocity
    Real lon_c = 1.5*M_PI + time_index * omega, lat_c = 0.0;

    // Construct the solution on the coordinates at its proper location at
    // the given time.
    KokkosTypes::view_2d<Pack> h_data("h", num_cols, num_levels);
    Real R = 1.0/3.0;
    for (int i = 0; i < num_cols; ++i) {
      Real lon = grid.longitudes[i], lat = grid.latitudes[i];
      Real r = acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c));
      for (int k = 0; k < num_levels; ++k) {
        Real z = (k + 0.5) * dz;
        if (r < R) {
          h_data(i, k) = 0.5 * z * (1.0 + cos(M_PI*r/R));
        } else {
          h_data(i, k) = 0.0;
        }
      }
    }
    Kokkos::deep_copy(data, h_data);
  }

  KOKKOS_INLINE_FUNCTION
  void linear_combination(const ThreadTeam& team,
                          Real a, Real b,
                          const CosineBell& x1, const CosineBell& x2) {
    int i = team.league_rank();
    int num_levels = team.team_size();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levels),
      KOKKOS_LAMBDA(int j) {
        data(i, j) = a * x1.data(i, j) + b * x2.data(i, j);
    });
  }

  KOKKOS_INLINE_FUNCTION
  void interpolate_horizontally(const ThreadTeam& team,
                                const TetralinearInterpWeights& weights,
                                const CosineBell& src_data, int tgt_column) {
    int i = tgt_column;
    int i1 = weights.indices[0], i2 = weights.indices[1],
        i3 = weights.indices[2], i4 = weights.indices[3];

    Real W1 = weights.weights[0], W2 = weights.weights[1],
         W3 = weights.weights[2], W4 = weights.weights[3];

    // Apply the interpolation weights at each level.
    int num_levels = team.team_size();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levels),
      KOKKOS_LAMBDA(int k) {
        data(i, k) = W1 * src_data.data(i1, k) + W2 * src_data.data(i2, k) +
                     W3 * src_data.data(i3, k) + W4 * src_data.data(i4, k);
    });
  }

  KOKKOS_INLINE_FUNCTION
  void compute_vertical_coords(const ThreadTeam& team,
                               int column_index,
                               VCoordColumnView& col_vcoords) const {
    // All columns have the same vertical coordinates.
    int num_levels = team.team_size();
    Real dz = 1000.0 / num_levels;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levels),
      KOKKOS_LAMBDA(int k) { col_vcoords(k) = (k + 0.5) * dz; });
  }

  KOKKOS_INLINE_FUNCTION
  void interpolate_vertically(const ThreadTeam& team,
                              const ekat::LinInterp<Real, Pack::n>& vert_interp,
                              int variable_index, int column_index,
                              const VCoordColumnView& src_col_vcoords,
                              const CosineBell& src_data,
                              const VCoordColumnView& tgt_col_vcoords) {
    EKAT_ASSERT(variable_index == 0);
    VCoordColumnView src_col_var = ekat::subview(src_data.data, column_index);
    VCoordColumnView tgt_col_var = ekat::subview(data, column_index);
    vert_interp.lin_interp(team, src_col_vcoords, tgt_col_vcoords,
                           src_col_var, tgt_col_var, column_index);
  }
};

TEST_CASE("tetralinear_interp_cosine_bell", "") {

  // Construct an interpolator that interpolates our cosine bell profile to a
  // selected set of (lon, lat) points.
  HostHCoordView lons;
  HostHCoordView lats;
  TetralinearInterp<CosineBell> interp("ne4.g", lons, lats);

  // Interpolate weekly over the course of a year.
//  for (size_t i = 0; i < 52; ++i) {
//  }
//  interp(t, src_vcoords, tgt_vcoords, tgt_data);
}

} // anonymous namespace
