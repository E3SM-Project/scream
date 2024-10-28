#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_universal_constants.hpp"

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ncols,
                                        const int nlevs) {
  const int num_global_cols = ncols * comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names", vos_t{"Point Grid"});
  auto &pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type", "point_grid");
  pl.set("aliases", vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm, gm_params);
  gm->build_grids();

  return gm;
}

TEST_CASE("horiz_avg") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;
  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 3;
  constexpr int dim3  = 4;
  constexpr int dim4  = 2;
  const int ngcols    = 6 * comm.size();

  auto gm1   = create_gm(comm, ngcols, 1);
  auto gm2   = create_gm(comm, ngcols, nlevs);
  auto grid1 = gm1->get_grid("Physics");
  auto grid2 = gm2->get_grid("Physics");

  // Input (randomized) qc
  FieldLayout scalar1d_layout{{COL}, {ngcols}};
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldLayout scalar3d_layout{{COL, CMP, LEV}, {ngcols, dim3, nlevs}};
  FieldLayout scalar4d_layout{{COL, CMP, CMP, LEV},
                              {ngcols, dim3, dim4, nlevs}};

  FieldIdentifier qc1_fid("qc", scalar1d_layout, kg / kg, grid1->name());
  FieldIdentifier qc2_fid("qc", scalar2d_layout, kg / kg, grid2->name());
  FieldIdentifier qc3_fid("qc", scalar3d_layout, kg / kg, grid2->name());
  FieldIdentifier qc4_fid("qc", scalar4d_layout, kg / kg, grid2->name());

  Field qc1(qc1_fid);
  Field qc2(qc2_fid);
  Field qc3(qc3_fid);
  Field qc4(qc4_fid);

  qc1.allocate_view();
  qc2.allocate_view();
  qc3.allocate_view();
  qc4.allocate_view();

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.0, 200.0);

  auto engine = scream::setup_random_test();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  // REQUIRE_THROWS(diag_factory.create("HorizAvgDiag", comm,
  //                                    params));  // No 'field_name' parameter

  // Set time for qc and randomize its values
  qc1.get_header().get_tracking().update_time_stamp(t0);
  qc2.get_header().get_tracking().update_time_stamp(t0);
  qc3.get_header().get_tracking().update_time_stamp(t0);
  qc4.get_header().get_tracking().update_time_stamp(t0);
  randomize(qc1, engine, pdf);
  randomize(qc2, engine, pdf);
  randomize(qc3, engine, pdf);
  randomize(qc4, engine, pdf);

  // Create and set up the diagnostic
  params.set("grid_name", grid1->name());
  params.set<std::string>("field_name", "qc");
  auto diag1 = diag_factory.create("HorizAvgDiag", comm, params);
  auto diag2 = diag_factory.create("HorizAvgDiag", comm, params);
  auto diag3 = diag_factory.create("HorizAvgDiag", comm, params);
  auto diag4 = diag_factory.create("HorizAvgDiag", comm, params);
  diag1->set_grids(gm1);
  diag2->set_grids(gm2);
  diag3->set_grids(gm2);
  diag4->set_grids(gm2);

  auto area = grid1->get_geometry_data("area");

  diag1->set_required_field(qc1);
  diag1->initialize(t0, RunType::Initial);

  diag1->compute_diagnostic();
  auto diag1_f = diag1->get_diagnostic();

  FieldIdentifier diag0_fid("qc_horiz_avg_manual",
                            scalar1d_layout.clone().strip_dim(COL), kg / kg,
                            grid1->name());
  Field diag0(diag0_fid);
  diag0.allocate_view();
  auto diag0_v = diag0.get_view<Real>();

  auto qc1_v  = qc1.get_view<Real *>();
  auto area_v = area.get_view<const Real *>();

  // calculate total area
  Real atot = 0.0;
  Kokkos::parallel_reduce(
      "HorizAvgDiag::compute_diagnostic_impl::total_area", ngcols,
      KOKKOS_LAMBDA(const int icol, Real &local_atot) {
        local_atot += area_v[icol];
      },
      atot);
  // calculate weighted avg
  Real wavg = 0.0;
  Kokkos::parallel_reduce(
      "HorizAvgDiag::compute_diagnostic_impl::weighted_sum", ngcols,
      KOKKOS_LAMBDA(const int icol, Real &local_wavg) {
        local_wavg += (area_v[icol] / atot) * qc1_v[icol];
      },
      wavg);
  Kokkos::deep_copy(diag0_v, wavg);

  auto diag1_v = diag1_f.get_view<Real>();

  bool result = false;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<>(0, 1),
      KOKKOS_LAMBDA(const int, bool &local_result) {
        local_result = (diag1_v() == diag0_v());
      },
      result);
  REQUIRE(result);

  // Try known cases
  // Set qc1_v to 1.0 to get weighted average of 1.0
  Kokkos::deep_copy(qc1_v, 1.0);
  Kokkos::deep_copy(diag0_v, 1.0);
  diag1->compute_diagnostic();
  auto diag1_v2 = diag1_f.get_view<Real>();
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<>(0, 1),
      KOKKOS_LAMBDA(const int, bool &local_result) {
        local_result = (diag1_v2() == diag0_v());
      },
      result);
  REQUIRE(result);

  // other diags
  auto qc2_v = qc2.get_view<Real **>();
  Kokkos::deep_copy(qc2_v, 5.0);
  FieldIdentifier diag2_fid("qc_horiz_avg_manual",
                            scalar2d_layout.clone().strip_dim(COL), kg / kg,
                            grid2->name());
  Field diag2_manual(diag2_fid);
  diag2_manual.allocate_view();
  auto diag2_manual_v = diag2_manual.get_view<Real *>();
  Kokkos::deep_copy(diag2_manual_v, 5.0);

  diag2->set_required_field(qc2);
  diag2->initialize(t0, RunType::Initial);
  diag2->compute_diagnostic();
  auto diag2_f = diag2->get_diagnostic();

  REQUIRE(views_are_equal(diag2_f, diag2_manual));

  auto qc3_v = qc3.get_view<Real ***>();
  FieldIdentifier diag3_manual_fid("qc_horiz_avg_manual",
                                   scalar3d_layout.clone().strip_dim(COL),
                                   kg / kg, grid2->name());
  Field diag3_manual(diag3_manual_fid);
  diag3_manual.allocate_view();
  auto diag3_manual_v = diag3_manual.get_view<Real **>();
  // calculate diag3_manual by hand
  TeamPolicy p(dim3 * nlevs, ngcols);
  Kokkos::parallel_for(
      "HorizAvgDiag::compute_diagnostic_impl::manual_diag3", p,
      KOKKOS_LAMBDA(const TeamMember &m) {
        const int idx = m.league_rank();
        const int j   = idx / nlevs;
        const int k   = idx % nlevs;
        Real sum      = 0.0;
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(m, ngcols),
            [&](const int icol, Real &accum) {
              accum += (area_v(icol) / atot) * qc3_v(icol, j, k);
            },
            sum);
        Kokkos::single(Kokkos::PerTeam(m),
                       [&]() { diag3_manual_v(j, k) = sum; });
      });
  diag3->set_required_field(qc3);
  diag3->initialize(t0, RunType::Initial);
  diag3->compute_diagnostic();
  auto diag3_f = diag3->get_diagnostic();
  REQUIRE(views_are_equal(diag3_f, diag3_manual));

  // TODO: add a different flavor of testing
  // TODO: how to test the MPI part of this rigorously?
  // TODO: how to test a different type of grid (especially to test the area
  // weighting)
}

}  // namespace scream
