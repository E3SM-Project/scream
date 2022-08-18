#include <catch2/catch.hpp>

#include <share/interp/tetralinear_interp.hpp>

namespace {

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

  CoordField allocate() const {
    return CoordField{
      KokkosTypes::view_2d<Real>("coords", data.extent(0), data.extent(1))
    };
  }

  void linear_combination(Real a, Real b,
                          const CoordField& x1, const CoordField& x2) {
    int num_cols = data.extent(0), num_levels = data.extent(1);
    auto team_policy = TeamPolicy(num_cols, num_levels);
    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const TeamType& team) {
      int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levels),
        KOKKOS_LAMBDA(int j) {
          data(i, j) = a * x1.data(i, j) + b * x2.data(i, j);
        });
    });
  }

  void compute_vertical_coords(VCoordView& vcoords) {
    int num_cols = data.extent(0), num_levels = data.extent(1);
    auto team_policy = TeamPolicy(num_cols, num_levels);
    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const TeamType& team) {
      int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levels),
        KOKKOS_LAMBDA(int j) {
          vcoords(i, j) = data(i, j);
        });
    });
  }

  void interpolate_vertically(const VCoordView& src_vcoords,
                              const Data& src_data,
                              const VCoordView& tgt_vcoords) {
    // FIXME: Seems like this could be better decomposed.
    EKAT_REQUIRE(src_vcoords.extent(0)==tgt_vcoords.extent(0));
    EKAT_REQUIRE(src_vcoords.extent(1)==tgt_vcoords.extent(1));
  int num_cols = tgt_vcoords.extent(0), num_levels = tgt_vcoords.extent(1);

  using LIV = ekat::LinInterp<Real,Spack::n>;

  const int nlevs_src = input.nlevs;
  const int nlevs_tgt = output.nlevs;

  LIV vert_interp(ncols,nlevs_src,nlevs_tgt);

  // We can ||ize over columns as well as over variables and bands
  const int num_vars = 1+input.nswbands*3+input.nlwbands;
  const int num_vert_packs = ekat::PackInfo<Spack::n>::num_packs(nlevs_tgt);
  const auto policy_setup = ESU::get_default_team_policy(ncols, num_vert_packs);

  // Setup the linear interpolation object
  Kokkos::parallel_for("spa_vert_interp_setup_loop", policy_setup,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {

    const int icol = team.league_rank();

    // Setup
    vert_interp.setup(team, ekat::subview(p_src,icol),
                            ekat::subview(p_tgt,icol));
  });
  Kokkos::fence();

  // Now use the interpolation object in || over all variables.
  const int outer_iters = ncols*num_vars;
  const auto policy_interp = ESU::get_default_team_policy(outer_iters, num_vert_packs);
  Kokkos::parallel_for("spa_vert_interp_loop", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {

    const int icol = team.league_rank() / num_vars;
    const int ivar = team.league_rank() % num_vars;

    const auto x1 = ekat::subview(p_src,icol);
    const auto x2 = ekat::subview(p_tgt,icol);

    const auto y1 = get_var_column(input, icol,ivar);
    const auto y2 = get_var_column(output,icol,ivar);

    vert_interp.lin_interp(team, x1, x2, y1, y2, icol);
  });
  Kokkos::fence();
  }

  void apply_interp_weights(const TetralinearInterpWeightMap& W,
                            const Data& x) {
    int num_levels = data.extent(1);
    auto team_policy = TeamPolicy(W.capacity(), num_levels);
    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const TeamType& team) {
      int k = team.league_rank();
      if( W.valid_at(k) ) {
        // fetch mapped weights for this column.
        auto i   = W.key_at(k);
        auto wts = W.value_at(k);

        int i1 = wts.indices[0],
            i2 = wts.indices[1],
            i3 = wts.indices[2],
            i4 = wts.indices[3];

        Real W1 = wts.weights[0],
             W2 = wts.weights[1],
             W3 = wts.weights[2],
             W4 = wts.weights[3];

        // apply the interpolation weights at each level.
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levels),
          KOKKOS_LAMBDA(int j) {
            data(i, j) = W1 * x.data(i1, j) + W2 * x.data(i2, j) +
                         W3 * x.data(i3, j) + W4 * x.data(i4, j);
          });
      }
    });
  }
};

TEST_CASE("tetralinear_interp", "") {
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
}

} // anonymous namespace
