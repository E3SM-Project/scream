#include <catch2/catch.hpp>

#include <share/interp/tetralinear_interp.hpp>

namespace {

using scream::Real;
using namespace scream::interpolators;

TEST_CASE("coarse_grid", "") {

  // Read an ne4 grid from a file.
  CoarseGrid grid("ne4.g");
  REQUIRE(grid.ne == 4);
  int num_elems = 6*4*4;
  REQUIRE(grid.num_elems() == num_elems);
  REQUIRE(grid.elements.size() == num_elems);
  int num_vertices = num_elems*num_elems + 2;
  REQUIRE(grid.latitudes.size() == num_vertices);
  REQUIRE(grid.longitudes.size() == num_vertices);

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

// Here's a data type that stores a 3D coordinate field on a unit sphere,
// defined on vertices of quadrilateral prism elements. We use it to test the
// tetralinear interpolator.
struct CoordField {
  // 3D (2D horizontal + 1D vertical) coordinate data
  KokkosTypes::view_2d<Real> data;

  // member functions used by traits (see TetralinearInterpTraits for details)

  static CoordField allocate(int num_cols, int num_levels) const {
    return CoordField{
      KokkosTypes::view_2d<Real>("coords", num_cols, num_levels)
    };
  }

  KOKKOS_INLINE_FUNCTION
  void linear_combination(const ThreadTeam& team,
                          Real a, Real b,
                          const CoordField& x1, const CoordField& x2) {
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
                                const CoordField& src_data, int tgt_column) {
    int i = tgt_column;
    int i1 = weights.indices[0], i2 = weights.indices[1],
        i3 = weights.indices[2], i4 = weights.indices[3];

    Real W1 = weights.weights[0], W2 = weights.weights[1],
         W3 = weights.weights[2], W4 = weights.weights[3];

    // apply the interpolation weights at each level.
    int num_levels = team.team_size();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levels),
      KOKKOS_LAMBDA(int j) {
        data(i, j) = W1 * src_data.data(i1, j) + W2 * src_data.data(i2, j) +
                     W3 * src_data.data(i3, j) + W4 * src_data.data(i4, j);
    });
  }

  KOKKOS_INLINE_FUNCTION
  void interpolate_vertically(const ThreadTeam& team,
                              const ekat::LinInterp<Real, Pack::n>& vert_interp,
                              int variable_index, int column_index,
                              const VCoordColumnView& src_col_vcoords,
                              const CoordField& src_data,
                              const VCoordColumnView& tgt_col_vcoords) {
    VCoordColumnView src_col_var = ekat::subview(src_data.data,
                                                 column_index,
                                                 variable_index);
    VCoordColumnView tgt_col_var = ekat::subview(data,
                                                 column_index,
                                                 variable_index);
    vert_interp.lin_interp(team, src_col_vcoords, tgt_col_vcoords,
                           src_col_var, tgt_col_var, column_index);
  }
};

TEST_CASE("tetralinear_interp", "") {
  /*
  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_local_elems = 10;
  const int num_gp = 4;
  const int num_levels = 72;

  auto gm = create_mesh_free_grids_manager(comm,num_local_elems,num_gp,num_levels,0);
  gm->build_grids(std::set<std::string>{"SE Grid"});

  // SE grid
  auto se_grid = gm->get_grid("SE Grid");

  REQUIRE(se_grid->type() == GridType::SE);
  REQUIRE(se_grid->name() == "SE Grid");
  REQUIRE(se_grid->get_num_vertical_levels() == num_levels);
  REQUIRE(se_grid->get_num_local_dofs() == num_local_elems*num_gp*num_gp);

  auto layout = se_grid->get_2d_scalar_layout();
  REQUIRE(layout.tags().size() == 3);
  REQUIRE(layout.tag(0) == EL);
  REQUIRE(layout.tag(1) == GP);
  REQUIRE(layout.tag(2) == GP);

  REQUIRE (se_grid->is_unique());

  const auto max_gid = se_grid->get_global_max_dof_gid();
  const auto min_gid = se_grid->get_global_min_dof_gid();
  REQUIRE( (max_gid-min_gid+1)==se_grid->get_num_global_dofs() );
  */
}

} // anonymous namespace
