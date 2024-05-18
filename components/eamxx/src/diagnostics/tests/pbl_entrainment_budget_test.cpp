#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/scream_common_physics_functions.hpp"
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

TEST_CASE("entrainment_budget") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using PC = scream::physics::Constants<Real>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // time stamps
  util::TimeStamp t0({2022, 1, 1}, {0, 1, 0});
  util::TimeStamp t1({2022, 1, 1}, {0, 2, 0});
  util::TimeStamp t2({2022, 1, 1}, {0, 3, 0});

  // Create a grids manager - single column for these tests
  const int ngcols = 1 * comm.size();
  const int nlevs = 20;
  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier qv_fid("qv", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier tm_fid("T_mid", scalar2d_layout, K, grid->name());
  FieldIdentifier pm_fid("p_mid", scalar2d_layout, Pa, grid->name());
  FieldIdentifier pd_fid("pseudo_density", scalar2d_layout, Pa, grid->name());
  FieldIdentifier sd_fid("SW_flux_dn", scalar2d_layout, W / m * m,
                         grid->name());
  FieldIdentifier su_fid("SW_flux_up", scalar2d_layout, W / m * m,
                         grid->name());
  FieldIdentifier ld_fid("LW_flux_dn", scalar2d_layout, W / m * m,
                         grid->name());
  FieldIdentifier lu_fid("LW_flux_up", scalar2d_layout, W / m * m,
                         grid->name());

  Field qc(qc_fid);
  qc.allocate_view();
  qc.get_header().get_tracking().update_time_stamp(t0);
  Field qv(qv_fid);
  qv.allocate_view();
  qv.get_header().get_tracking().update_time_stamp(t0);
  Field tm(tm_fid);
  tm.allocate_view();
  tm.get_header().get_tracking().update_time_stamp(t0);
  Field pm(pm_fid);
  pm.allocate_view();
  pm.get_header().get_tracking().update_time_stamp(t0);
  Field pd(pd_fid);
  pd.allocate_view();
  pd.get_header().get_tracking().update_time_stamp(t0);
  Field sd(sd_fid);
  sd.allocate_view();
  sd.get_header().get_tracking().update_time_stamp(t0);
  Field su(su_fid);
  su.allocate_view();
  su.get_header().get_tracking().update_time_stamp(t0);
  Field ld(ld_fid);
  ld.allocate_view();
  ld.get_header().get_tracking().update_time_stamp(t0);
  Field lu(lu_fid);
  lu.allocate_view();
  lu.get_header().get_tracking().update_time_stamp(t0);

  // gravity
  constexpr Real g = PC::gravit;

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0, 0.05);
  auto engine = scream::setup_random_test();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  constexpr int ntests = 1;
  for(int itest = 0; itest < ntests; ++itest) {
    // Randomize
    randomize(qc, engine, pdf);
    randomize(qv, engine, pdf);
    randomize(tm, engine, pdf);
    randomize(pm, engine, pdf);
    randomize(pd, engine, pdf);
    randomize(sd, engine, pdf);
    randomize(su, engine, pdf);
    randomize(ld, engine, pdf);
    randomize(lu, engine, pdf);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    auto diag = diag_factory.create("PBLEntrainmentBudget", comm, params);
    diag->set_grids(gm);
    diag->set_required_field(qc);
    diag->set_required_field(qv);
    diag->set_required_field(tm);
    diag->set_required_field(pm);
    diag->set_required_field(pd);
    diag->set_required_field(sd);
    diag->set_required_field(su);
    diag->set_required_field(ld);
    diag->set_required_field(lu);

    diag->initialize(t0, RunType::Initial);

    qc.get_header().get_tracking().update_time_stamp(t1);
    qv.get_header().get_tracking().update_time_stamp(t1);
    tm.get_header().get_tracking().update_time_stamp(t1);
    pm.get_header().get_tracking().update_time_stamp(t1);
    pd.get_header().get_tracking().update_time_stamp(t1);
    sd.get_header().get_tracking().update_time_stamp(t1);
    su.get_header().get_tracking().update_time_stamp(t1);
    ld.get_header().get_tracking().update_time_stamp(t1);
    lu.get_header().get_tracking().update_time_stamp(t1);

    qc.deep_copy(0.001);
    qv.deep_copy(0.002);
    tm.deep_copy(290.0);
    pm.deep_copy(1e5);
    pd.deep_copy(0.0);
    sd.deep_copy(0.0);
    su.deep_copy(0.0);
    ld.deep_copy(0.0);
    lu.deep_copy(0.0);

    // Run diag
    diag->init_timestep(t1);

    randomize(qc, engine, pdf);
    randomize(qv, engine, pdf);
    randomize(tm, engine, pdf);
    randomize(pm, engine, pdf);
    randomize(pd, engine, pdf);
    randomize(sd, engine, pdf);
    randomize(su, engine, pdf);
    randomize(ld, engine, pdf);
    randomize(lu, engine, pdf);

    qc.get_header().get_tracking().update_time_stamp(t2);
    qv.get_header().get_tracking().update_time_stamp(t2);
    tm.get_header().get_tracking().update_time_stamp(t2);
    pm.get_header().get_tracking().update_time_stamp(t2);
    pd.get_header().get_tracking().update_time_stamp(t2);
    sd.get_header().get_tracking().update_time_stamp(t2);
    su.get_header().get_tracking().update_time_stamp(t2);
    ld.get_header().get_tracking().update_time_stamp(t2);
    lu.get_header().get_tracking().update_time_stamp(t2);

    diag->compute_diagnostic();

    diag->get_diagnostic().sync_to_host();

    const auto out_hf = diag->get_diagnostic();
    auto out_hv       = out_hf.get_view<Real **, Host>();

    for(int idx = 0; idx < 8; ++idx) {
      std::cout << "idx = " << idx << " ---> " << out_hv(0, idx) << std::endl;
    }
  }
}  // TEST_CASE("entrainment_budget")

}  // namespace scream
