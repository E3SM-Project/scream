#include <catch2/catch.hpp>

#include "share/grid/remap/horizontal_remap_utility.hpp"

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_output_manager.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "share/field/field_manager.hpp"

#include "share/util/scream_setup_random_test.hpp"

#include "ekat/util/ekat_units.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
namespace {

using namespace scream;

std::shared_ptr<GridsManager> get_test_gm(const ekat::Comm& comm, const Int num_gcols, const Int num_levs);

std::shared_ptr<FieldManager> get_test_fm(std::shared_ptr<const AbstractGrid> grid);

AtmosphereInput setup_data_input(std::shared_ptr<const FieldManager> field_manager);

ekat::ParameterList setup_output_manager_params();

GSMap get_test_map(const ekat::Comm& comm, const view_1d<gid_type>& dofs_gids, const gid_type min_dof);

void run(std::mt19937_64& engine, const ekat::Comm& comm)
{
  // Step 1: Generate a set of random remapping (source col, target col, weight).
  //         Store these in local 1D views and create a remap file.
} // end function run

//===============================================================================
TEST_CASE("horizontal_remap_test", "[horizontal_remap_test]"){

  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  ekat::Comm comm (MPI_COMM_WORLD);
  auto engine = scream::setup_random_test();

// Run tests of vertically dimensioned-functions for both Real and Pack,
// and for (potentially) different pack sizes
  if (comm.am_i_root()) {
    printf(" -> Number of randomized runs: %d\n\n", num_runs);
    printf(" -> Testing horizontal remapping...");
  }
  for (int irun=0; irun<num_runs; ++irun) {
    run(enginei, comm);
  }
  if (comm.am_i_root()) {
    printf("ok!\n");
  }

} // TEST_CASE
//===============================================================================
TEST_CASE("horizontal_remap_from_file","") {

{

  ekat::Comm comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  Int num_gcols = 218;
  Int num_levs  = 72;
  Int num_dcols = 866;
  Int num_dlevs = 72;
  MPI_Fint fcomm = MPI_Comm_c2f(comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // Construct a timestamp
  util::TimeStamp t0 ({2000,1,1},{0,0,0});

  // First set up a field manager and grids manager to interact with the inputs
  auto gm_d            = get_test_gm(comm,num_dcols,num_dlevs);
  auto grid_d          = gm_d->get_grid("Point Grid");
  auto field_manager_d = get_test_fm(grid_d);
  field_manager_d->init_fields_time_stamp(t0);
  // Get input parameters and setup input
  auto data_input = setup_data_input(field_manager_d);
  // Get view of PS for data
  auto ps_d_f   = field_manager_d->get_field("PS");
  auto ps_d_v   = ps_d_f.get_view<Real*>();
  auto ps_d_v_h = ps_d_f.get_view<Real*,Host>();

  // Now set up a field_manager and grids manager for the test
  auto gm_t            = get_test_gm(comm,num_gcols,num_levs);
  auto grid_t          = gm_t->get_grid("Point Grid");
  auto field_manager_t = get_test_fm(grid_t);
  field_manager_t->init_fields_time_stamp(t0);
  // Get view of PS for remapping
  auto ps_t_f   = field_manager_t->get_field("PS");
  auto ps_t_v   = ps_t_f.get_view<Real*>();
  auto ps_t_v_h = ps_t_f.get_view<Real*,Host>();

  // Setup an output manager for this data
  auto om_params = setup_output_manager_params();
  OutputManager om;
  om.setup(comm,om_params,field_manager_t,gm_t,t0,t0,false);
  comm.barrier();
  return;

  // Set up GSMap
  auto dofs_gids = grid_t->get_dofs_gids();
  auto min_dof   = grid_t->get_global_min_dof_gid();
  auto test_map  = get_test_map(comm,dofs_gids,min_dof);
  // Time loop
  util::TimeStamp time = t0;
  for (int ii=0;ii<12;ii++) {
    data_input.read_variables(ii+1);
    ps_d_f.sync_to_host();
    test_map.apply_remap(ps_d_v,ps_t_v);
    ps_t_f.sync_to_host();
    time += 1;
    om.run(time);
  }
  om.finalize();

  // The rest of the test is in a CPRNC comparison which is called after this test has finished.

}

} // end TEST_CASE horizontal_remap_from_file

TEST_CASE("horizontal_remap_units", "") {

  using namespace scream;
  using IPDF = std::uniform_int_distribution<int>;
  using RPDF = std::uniform_real_distribution<Real>;
  auto engine = scream::setup_random_test();

  // Test the remap_segment structure
  {
    IPDF pdf_seg_len(2,100);
    RPDF pdf_seg_wgt(0,1);
    RemapSegment test_segment;
    Int len = pdf_seg_len(engine);
    test_segment.m_length = len-1;
    test_segment.source_dofs = view_1d<gid_type>("",len);
    test_segment.weights     = view_1d<Real>("",len);
    ekat::genRandArray(test_segment.weights,engine,pdf_seg_wgt);
    Real wgt_normalize = 0;
    auto weights_h = Kokkos::create_mirror_view(test_segment.weights);
    Kokkos::deep_copy(weights_h,test_segment.weights);
    for (int ii=0; ii<len; ii++) {
      wgt_normalize += weights_h(ii);
    }
    for (int ii=0; ii<len; ii++) {
      weights_h(ii) /= wgt_normalize;
    };
    Kokkos::deep_copy(test_segment.weights,weights_h);
    // Checking this segment should fail because the weight view is a different length than
    // the reported length for the segment.
    REQUIRE_THROWS(test_segment.check());
    // Now switch to passing by bumping up the length
    test_segment.m_length = len;
    test_segment.check();
    // Make it fail again by adjusting one of the weights so that they no longer add to 1.0
    weights_h(0) += 0.1;
    Kokkos::deep_copy(test_segment.weights,weights_h);
    REQUIRE_THROWS(test_segment.check());
    weights_h(0) -= 0.2;
    Kokkos::deep_copy(test_segment.weights,weights_h);
    REQUIRE_THROWS(test_segment.check());
  }


  // Create a map with two remap segments for some simple map structure testing.
  {
    // Create a GSMap to test
    GSMap test_map;
    // Create a remap segment
    RemapSegment test_seg;
    IPDF pdf_seg_len(2,100);
    RPDF pdf_seg_wgt(0,1);
    Int len = pdf_seg_len(engine);
    test_seg.m_length = len;
    test_seg.m_dof    = 1;
    test_seg.weights     = view_1d<Real>("",len);
    test_seg.source_dofs = view_1d<gid_type>("",len);
    ekat::genRandArray(test_seg.weights,engine,pdf_seg_wgt);
    Real wgt_normalize = 0;
    auto weights_h = Kokkos::create_mirror_view(test_seg.weights);
    Kokkos::deep_copy(weights_h,test_seg.weights);
    for (int ii=0; ii<len; ii++) {
      wgt_normalize += weights_h(ii);
    }
    for (int ii=0; ii<len; ii++) {
      weights_h(ii) /= wgt_normalize;
    };
    Kokkos::deep_copy(test_seg.weights,weights_h);
    // Add segment to map
    test_map.add_segment(test_seg);
    test_map.check();
    REQUIRE(test_map.get_num_of_segs()==1);
    auto check_seg = test_map.get_segment(1);
    REQUIRE(check_seg.m_length==len);
    weights_h = Kokkos::create_mirror_view(check_seg.weights);
    Kokkos::deep_copy(weights_h,check_seg.weights);
    weights_h(0) = weights_h(0)/2.0;
    Real new_wgt = weights_h(0);
    Kokkos::deep_copy(check_seg.weights,weights_h);
    REQUIRE_THROWS(test_map.check());
    REQUIRE_THROWS(test_map.get_segment(2));
    // Create a new segment that has the same DOF as the first segment
    RemapSegment new_seg;
    new_seg.m_length = 1;
    new_seg.m_dof    = 1;
    new_seg.weights     = view_1d<Real>("",1);
    new_seg.source_dofs = view_1d<gid_type>("",1);
    weights_h = Kokkos::create_mirror_view(new_seg.weights);
    weights_h(0)  = new_wgt;
    Kokkos::deep_copy(new_seg.weights,weights_h);
    test_map.add_segment(new_seg);
    // Because the new segment is the same DOF, add_segment should have
    // appended it to the first segment.
    REQUIRE(test_map.get_num_of_segs()==1);
    test_map.check();
    check_seg = test_map.get_segment(1);
    REQUIRE(check_seg.m_length==(len+1));
    // Create a completely new segment for a different DOF and check
    new_seg.m_dof = 2;
    weights_h(0) = 1.0;
    Kokkos::deep_copy(new_seg.weights,weights_h);
    test_map.add_segment(new_seg);
    REQUIRE(test_map.get_num_of_segs()==2);
    test_map.check();
    check_seg = test_map.get_segment(2);
    REQUIRE(check_seg.m_length==(1));
  }

  // Test retrieving unique columns from a full map 
  {
    // Create a GSMap to test
    GSMap test_map;
    // Create a set of segments each with overlapping source dofs
    IPDF pdf_seg_len(4,10),
         pdf_num_seg(5,15);
    Int num_of_segs = pdf_num_seg(engine);
    gid_type seg_start = 0;
    Int seg_len;
    for (int iseg=0; iseg<num_of_segs; iseg++) {
      seg_len = pdf_seg_len(engine);
      RemapSegment seg(iseg,seg_len);
      Kokkos::parallel_for("", seg_len, KOKKOS_LAMBDA (const int& ii) {
        seg.source_dofs(ii) = seg_start + ii;
      });
      test_map.add_segment(seg);
      seg_start += seg_len-2;
    }
    // All of the segments combined should overlap by 2 on either side.  The set of unique
    // columns for all segments should be 0,1,...,seg_start+seg_len-1.  So we check that.
    const gid_type seg_end = seg_start + 2;
    test_map.set_unique_source_dofs();
    const auto& unique_dofs_vec = test_map.get_unique_dofs();
    view_1d<gid_type> unique_dofs("",unique_dofs_vec.size());
    auto unique_dofs_h = Kokkos::create_mirror_view(unique_dofs);
    Kokkos::deep_copy(unique_dofs_h,unique_dofs);
    REQUIRE(unique_dofs_h.extent(0) == seg_end);
    for (int ii=0; ii<unique_dofs.size(); ii++) {
      unique_dofs_h(ii) = unique_dofs_vec[ii];
    }
    Kokkos::deep_copy(unique_dofs,unique_dofs_h);
    Kokkos::parallel_for("", seg_end, KOKKOS_LAMBDA (const int& ii) {
      EKAT_KERNEL_REQUIRE(unique_dofs(ii) == ii);
    });
    
    

  }
} // TEST_CASE horizontal remap
/*===================================================================================================*/
std::shared_ptr<FieldManager> get_test_fm(std::shared_ptr<const AbstractGrid> grid)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;
  using FR = FieldRequest;

  // Create a fm
  auto fm = std::make_shared<FieldManager>(grid);

  const int num_lcols = grid->get_num_local_dofs();
  const int num_levs = grid->get_num_vertical_levels();

  // Create some fields for this fm
  std::vector<FieldTag> tag_h  = {COL};
  std::vector<FieldTag> tag_v  = {LEV};
  std::vector<FieldTag> tag_2d = {COL,LEV};

  std::vector<Int>     dims_h  = {num_lcols};
  std::vector<Int>     dims_v  = {num_levs};
  std::vector<Int>     dims_2d = {num_lcols,num_levs};

  const std::string& gn = grid->name();

  FieldIdentifier fid_ps("PS",FL{tag_h,dims_h},Pa,gn);

  // Register fields with fm
  // Make sure packsize isn't bigger than the packsize for this machine, but not so big that we end up with only 1 pack.
  fm->registration_begins();
  fm->register_field(FR{fid_ps});
  fm->registration_ends();

  // Initialize these fields
  auto f1 = fm->get_field(fid_ps);
  auto f1_host = f1.get_view<Real*,Host>(); 

  // Update timestamp
  util::TimeStamp time ({2000,1,1},{0,0,0});
  fm->init_fields_time_stamp(time);
  // Sync back to device
  f1.sync_to_dev();

  return fm;
} // end get_test_fm
/*==========================================================================================================*/
std::shared_ptr<GridsManager> get_test_gm(const ekat::Comm& comm, const Int num_gcols, const Int num_levs)
{
  ekat::ParameterList gm_params;
  gm_params.set("Number of Global Columns",num_gcols);
  gm_params.set("Number of Vertical Levels",num_levs);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids(std::set<std::string>{"Point Grid"});
  return gm;
} // end get_test_gm
/*==========================================================================================================*/
AtmosphereInput setup_data_input(std::shared_ptr<const FieldManager> field_manager)
{
  std::string fname = "horiz_remap_from_file_control.yaml";
  ekat::ParameterList params("horiz_remap_test");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,params) );
  auto& data_file = params.get<std::string>("Data File");
  using vos_type = std::vector<std::string>;
  ekat::ParameterList in_params("Input Parameters");
  in_params.set<std::string>("Filename",data_file);
  in_params.set<vos_type>("Field Names",{"PS"});
  AtmosphereInput data_input(in_params,field_manager);
  return data_input;
}
/*==========================================================================================================*/
ekat::ParameterList setup_output_manager_params()
{
  using vos_type = std::vector<std::string>;
  ekat::ParameterList om_params;
  om_params.set<std::string>("Casename","horizontal_remap_output");
  om_params.set<std::string>("Averaging Type","Instant");
  om_params.set<Int>("Max Snapshots Per File",12);
  om_params.set<vos_type>("Field Names",{"PS"});
  auto& om_control = om_params.sublist("Output Control");
  om_control.set<Int>("Frequency",1);
  om_control.set<std::string>("Frequency Units","Steps");
  om_control.set<bool>("Timestamp in Filename", false);
  om_control.set<bool>("MPI Ranks in Filename", false);
  om_control.set<bool>("AVG Type in Filename", false);
  om_control.set<bool>("Frequency in Filename", false);
  return om_params;
}
/*==========================================================================================================*/
GSMap get_test_map(const ekat::Comm& comm, const view_1d<gid_type>& dofs_gids, const gid_type min_dof)
{

  GSMap test_map;
  test_map.set_dofs_gids(dofs_gids);
  std::string fname = "horiz_remap_from_file_control.yaml";
  ekat::ParameterList params("horiz_remap_test");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,params) );
  auto& remap_file = params.get<std::string>("Remap File");
  test_map.set_segments_from_file(remap_file,comm,dofs_gids,min_dof);
  printf("ASD - (%d) - %d, %d\n",comm.rank(),test_map.get_num_of_segs(),test_map.get_num_of_dofs());
//  test_map.check();
//  test_map.set_unique_source_dofs();
  return test_map;

}
/*==========================================================================================================*/

} // end namespace
