#include "catch2/catch.hpp"

#include "diagnostics/simple_op.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_setup_random_test.hpp"

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,num_global_cols);
  gm->build_grids();

  return gm;
}

TEST_CASE("simpl_op")
{
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;

  constexpr int packsize = SCREAM_PACK_SIZE;

  ekat::Comm comm(MPI_COMM_WORLD);

  // Stuff for random test
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(-1.0,1.0);
  auto engine = scream::setup_random_test(&comm);

  // Create a grids manager
  const int ncols = 3;
  const int nlevs = packsize*2 + 1;
  auto gm = create_gm(comm,ncols,nlevs);
  auto grid = gm->get_grid("Point Grid");

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create random input fields
  const auto units = ekat::units::Units::invalid();

  FieldIdentifier fid_0d ("f1_0d",FL({           },{             }),units,grid->name());
  FieldIdentifier fid_1d ("f1_1d",FL({COL        },{ncols        }),units,grid->name());
  FieldIdentifier fid_2d ("f1_2d",FL({COL,LEV    },{ncols,  nlevs}),units,grid->name());
  FieldIdentifier fid_3d ("f1_3d",FL({COL,CMP,LEV},{ncols,2,nlevs}),units,grid->name());

  Field f1_0d (fid_0d,true);
  Field f1_1d (fid_1d,true);
  Field f1_2d (fid_2d,true);
  Field f1_3d (fid_3d,true);
  Field f2_0d = f1_0d.clone("f2_0d");
  Field f2_1d = f1_1d.clone("f2_1d");
  Field f2_2d = f1_2d.clone("f2_2d");
  Field f2_3d = f1_3d.clone("f2_3d");

  for (auto f : {&f1_0d,&f1_1d,&f1_2d,&f1_3d,&f2_0d,&f2_1d,&f2_2d,&f2_3d}) {
    f->get_header().get_tracking().update_time_stamp(t0);
  }

  // Helper
  auto get_f2 = [&](int rank) {
    return rank==0 ? f2_0d
                   : rank==1 ? f2_1d
                             : rank==2 ? f2_2d
                                       : f2_3d;
  };

  ekat::ParameterList params;
  auto create_diag = [&] () {
    return std::make_shared<SimpleOp>(comm,params);
  };

  printf(" -> Testing exceptions handling ...\n");
  params.set<std::string>("op","foobar");
  REQUIRE_THROWS (create_diag()); // Invalid op
  params.set<std::string>("op","add");
  REQUIRE_THROWS (create_diag()); // Missing f1
  params.set("f1",f1_1d);
  params.set<std::string>("f2_type","foobar");
  REQUIRE_THROWS (create_diag()); // Invalid f2_type
  params.set<std::string>("f2_type","field");
  REQUIRE_THROWS (create_diag()); // Missing f2
  printf(" -> Testing exceptions handling ... OK!\n");

  printf(" -> Testing op  neg ...............\n");
  params.set<std::string>("op","neg");
  for (auto f1 : {f1_0d, f1_1d, f1_2d, f1_3d}) {
    randomize(f1,engine,pdf);
    params.set("f1",f1);
    auto d = create_diag();
    d->set_required_field(f1);
    d->initialize(t0,RunType::Initial);
    d->compute_diagnostic();
    f1.scale(-1);
    if (!views_are_equal(f1,d->get_diagnostic())) {
      auto ddata = d->get_diagnostic().get_internal_view_data<Real>();
      auto fdata = f1.get_internal_view_data<Real>();
      int s = f1.get_header().get_identifier().get_layout().size();
      std::cout << " diag:";
      for (int i=0; i<s; ++i) { std::cout << " " << ddata[i]; } std::cout << "\n";
      std::cout << " expected:";
      for (int i=0; i<s; ++i) { std::cout << " " << fdata[i]; } std::cout << "\n";
    }
    REQUIRE (views_are_equal(f1,d->get_diagnostic()));
  }
  printf(" -> Testing op  neg ............... OK!\n");

  for (const std::string& op : {"add","sub","prod","div"}) {
    printf (" -> Testing op %4s ...............\n",op.c_str());
    params.set<std::string>("op",op);
    Real scale1 = pdf(engine);
    params.set("f1_scale",scale1);

    // f2_type=value
    params.set<std::string>("f2_type","value");
    Real val = pdf(engine);
    params.set("f2_value",val);
    for (auto f1 : {f1_0d, f1_1d, f1_2d, f1_3d}) {
      randomize(f1,engine,pdf);
      params.set("f1",f1);
      auto d = create_diag();
      d->set_required_field(f1);
      d->initialize(t0,RunType::Initial);
      d->compute_diagnostic();

      // Manually re-create scale1*f1 op val
      auto data1 = f1.get_internal_view_data<Real,Host>();
      for (int i=0; i<f1.get_header().get_identifier().get_layout().size(); ++i) {
        if (op=="add") {
          data1[i] = scale1*data1[i] + val;
        } else if (op=="sub") {
          data1[i] = scale1*data1[i] - val;
        } else if (op=="prod") {
          data1[i] = (scale1*val)*data1[i];
        } else {
          data1[i] = (scale1/val)*data1[i];
        }
      }
      f1.sync_to_dev();
      if (not views_are_equal(f1,d->get_diagnostic())) {
        auto datad = d->get_diagnostic().get_internal_view_data<Real,Host>();
        int s = f1.get_header().get_identifier().get_layout().size();
        std::cout << "f2type: value\n";
        std::cout << " diag:";
        for (int i=0; i<s; ++i) { std::cout << " " << datad[i]; } std::cout << "\n";
        std::cout << " exp :";
        for (int i=0; i<s; ++i) { std::cout << " " << data1[i]; } std::cout << "\n";
      }
      REQUIRE (views_are_equal(f1,d->get_diagnostic()));
    }

    // f2_type=field
    params.set<std::string>("f2_type","field");
    Real scale2 = pdf(engine);
    params.set("f2_scale",scale2);
    for (auto f1 : {f1_1d, f1_2d, f1_3d}) {
      auto f2 = get_f2(f1.rank());
      randomize(f1,engine,pdf);
      randomize(f2,engine,pdf);
      params.set("f1",f1);
      params.set("f2",f2);
      auto d = create_diag();
      d->set_required_field(f1);
      d->set_required_field(f2);
      d->initialize(t0,RunType::Initial);
      d->compute_diagnostic();

      // Manually re-create scale1*f1 op scale2*f2
      auto data1 = f1.get_internal_view_data<Real,Host>();
      auto data2 = f2.get_internal_view_data<Real,Host>();
      for (int i=0; i<f1.get_header().get_identifier().get_layout().size(); ++i) {
        if (op=="add") {
          data1[i] = scale1*data1[i] + scale2*data2[i];
        } else if (op=="sub") {
          data1[i] = scale1*data1[i] - scale2*data2[i];
        } else if (op=="prod") {
          data1[i] = scale1*data1[i] * scale2*data2[i];
        } else {
          data1[i] = scale1*data1[i] / scale2*data2[i];
        }
      }
      f1.sync_to_dev();
      if (not views_are_equal(f1,d->get_diagnostic())) {
        auto datad = d->get_diagnostic().get_internal_view_data<Real,Host>();
        int s = f1.get_header().get_identifier().get_layout().size();
        std::cout << "f2type: field\n";
        std::cout << " diag:";
        for (int i=0; i<s; ++i) { std::cout << " " << datad[i]; } std::cout << "\n";
        std::cout << " exp :";
        for (int i=0; i<s; ++i) { std::cout << " " << data1[i]; } std::cout << "\n";
      }
      REQUIRE (views_are_equal(f1,d->get_diagnostic()));
    }
    printf (" -> Testing op %4s ............... OK!\n",op.c_str());
  }
}

} // namespace scream
