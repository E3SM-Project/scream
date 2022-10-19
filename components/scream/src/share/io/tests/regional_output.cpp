// A set of tests for regional output
// 1. Test one_to_one regional output where
//     OUT(n) = IN(n) for all n = 0...ncols-1
// 2. Test one_to_one_shifted whre
//     OUT(n) = IN(ncols-1-n) for n = 0...ncols-1
// 3. Test subset output where
//     OUT(n) = IN(n) for n = 0...floor(ncols/2)
// 4. Test two_to_one coarsening output where
//     OUT(n) = 0.5*( IN(n) + IN(ncols-1-n) ) for n = 0...floor(ncols/2)


// Make a routine that writes a remap file given (remap_file_name,n_a, n_b, n_s, col, row, S)

// Create a routine to set col, row, n_s, S from (n_a, n_b, type)

// Create a routine to create the parameter list for output stream given (remap_file_name)

// Create a routine to provide src grid and src field manager (per io.cpp)

// Create a routine that returns source data given (icol,jlev)

// Create test which will
//   -  Create source grid and field manager
//   -  Populate a set of fields w/ different layouts, note we may want to add something when a field doesn't have columns, maybe no remapping, just writing.
//   -  Create each remap file using routines above
//   -  Run one timestep and create output
//   -  Read input, compare w/ what the true value should be.
#include <catch2/catch.hpp>

#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "share/util/scream_time_stamp.hpp"

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_units.hpp"

namespace {

using namespace scream;

using KT = KokkosTypes<DefaultDevice>;
template <typename S>
using view_1d = typename KT::template view_1d<S>;

std::shared_ptr<FieldManager>
get_test_fm(std::shared_ptr<const AbstractGrid> grid);

std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs);

void create_remap_file(const ekat::Comm& comm, const std::string& filename, std::shared_ptr<const AbstractGrid> grid,
                       const int n_a, const int n_b,
                       const view_1d<int>& col, const view_1d<int>& row, const view_1d<Real>& S);

Real generate_data_2d(const int col);
Real generate_data_2d(const int col,const int cmp);
Real generate_data_3d(const int col, const int lev);
Real generate_data_3d(const int col, const int lev,const int cmp);

/*===================================================================================================*/
void run()
{
  // Get a comm group for this test, set up basic grid
  ekat::Comm comm(MPI_COMM_WORLD);
  MPI_Fint fcomm = MPI_Comm_c2f(comm.mpi_comm());  
  scorpio::eam_init_pio_subsystem(fcomm);

  int num_gcols = 20*comm.size();
  int num_levs  = 2*SCREAM_SMALL_PACK_SIZE + 1;
  auto gm = get_test_gm(comm,num_gcols,num_levs);
  auto grid = gm->get_grid("Point Grid");
  // Create a test field manager, initialize the fields.
  auto field_manager = get_test_fm(grid);

  // ------------------------------------------------------------------------------------------- //
  // Create the remap file for each regional output test case
  int n_a, n_b, n_s;
  view_1d<int>  col, row;
  view_1d<Real> S;
  // Case 1: Test 1-1 mapping, i.e. no actual regional output, but we will run as
  //         though we are remapping.
  // The map will follow the pattern:
  //   Y(n) = 1.0 * X(n) for n = 0,...,ncol-1
  std::string case1_filename = "one_to_one_remap_np" + std::to_string(comm.size()) + ".nc";
  n_a = num_gcols;
  n_b = num_gcols;
  n_s = num_gcols;
  col = view_1d<int>("",n_s);
  row = view_1d<int>("",n_s);
  S   = view_1d<Real>("",n_s);
  Kokkos::parallel_for("", n_a, KOKKOS_LAMBDA (const int& ii) {
    col(ii) = ii+1;
    row(ii) = ii+1;
    S(ii)   = 1.0;
  });
  create_remap_file(comm,case1_filename,grid,n_a,n_b,col,row,S);
  
  // ------------------------------------------------------------------------------------------- //
  // Setup and write output using our regional output cases.
  // First we create a generic parameter list for the output manager.  We will simply change the
  // one parameter that controls the remap file name and Casename for each case.
  ekat::ParameterList params;
  params.set<std::string>("Casename","dummy");
  params.set<std::string>("remap_file","dummy");
  params.set<std::string>("Averaging Type","Instant");
  params.set<int>("Max Snapshots Per File",1);
  using vos_type = std::vector<std::string>;
  params.set<vos_type>("Field Names",
      {"field_2d", "field_2d_vector", "field_3d", "field_3d_vector"});
  auto& psub = params.sublist("output_control");
  psub.set<bool>("MPI Ranks in Filename",true);
  psub.set<int>("Frequency",1);
  psub.set<std::string>("frequency_units","nsteps");
  params.print();
  // Case 1:
  util::TimeStamp time ({2000,1,1},{0,0,0});
  params.set<std::string>("Casename","regional_one_to_one");
  params.set<std::string>("remap_file",case1_filename);
  OutputManager om_case_1;
  om_case_1.setup(comm,params,field_manager,gm,time,time,false);
  comm.barrier();
  om_case_1.run(time);
  om_case_1.finalize();

  

  // All Done 
  scorpio::eam_pio_finalize();
} // end function run
/*==========================================================================================================*/
std::shared_ptr<GridsManager> get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs)
{
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",num_gcols);
  gm_params.set("number_of_vertical_levels",num_levs);
  auto gm = create_mesh_free_grids_manager(io_comm,gm_params);
  gm->build_grids();
  return gm;
}
/*===================================================================================================*/
std::shared_ptr<FieldManager> get_test_fm(std::shared_ptr<const AbstractGrid> grid)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using FL = FieldLayout;
  using FR = FieldRequest;

  Units nondim(0,0,0,0,0,0,0);

  auto fm = std::make_shared<FieldManager>(grid);

  const int num_lcols = grid->get_num_local_dofs();
  const int num_levs  = grid->get_num_vertical_levels();
  const auto dof_gids = grid->get_dofs_gids();
  const int num_cmp   = 2;

  // TODO: Do we want to also include a field that is only levels, such a field wouldn't be remapped at all.
  std::vector<FieldTag> tag_2d        = {COL};
  std::vector<FieldTag> tag_2d_vector = {COL,CMP};
  std::vector<FieldTag> tag_3d        = {COL,LEV};
  std::vector<FieldTag> tag_3d_vector = {COL,CMP,LEV};

  std::vector<Int>      dims_2d        = {num_lcols};
  std::vector<Int>      dims_2d_vector = {num_lcols,num_cmp};
  std::vector<Int>      dims_3d        = {num_lcols,num_levs};
  std::vector<Int>      dims_3d_vector = {num_lcols,num_cmp,num_levs};

  const std::string& gn = grid->name();

  FieldIdentifier fid2d("field_2d",FL{tag_2d,dims_2d},nondim,gn);
  FieldIdentifier fid2d_vec("field_2d_vector",FL{tag_2d_vector,dims_2d_vector},nondim,gn);
  FieldIdentifier fid3d("field_3d",FL{tag_3d,dims_3d},nondim,gn);
  FieldIdentifier fid3d_vec("field_3d_vector",FL{tag_3d_vector,dims_3d_vector},nondim,gn);

  fm->registration_begins();
  fm->register_field(FR{fid2d});
  fm->register_field(FR{fid2d_vec});
  fm->register_field(FR{fid3d});
  fm->register_field(FR{fid3d_vec});
  fm->registration_ends();

  auto f2d      = fm->get_field(fid2d);
  auto f2dvec   = fm->get_field(fid2d_vec);
  auto f3d      = fm->get_field(fid3d);
  auto f3dvec   = fm->get_field(fid3d_vec);
  auto f2d_h    = f2d.get_view<Real*,Host>();
  auto f2dvec_h = f2dvec.get_view<Real**,Host>();
  auto f3d_h    = f3d.get_view<Real**,Host>();
  auto f3dvec_h = f3dvec.get_view<Real***,Host>();

  // Initialize these fields
  for (int col=0;col<num_lcols;col++) {
     f2d_h(col) = generate_data_2d(col);
   for (int cmp=0;cmp<num_cmp;cmp++) {
     f2dvec_h(col,cmp) = generate_data_2d(col,cmp);
     for (int lev=0;lev<num_levs;lev++) {
       f3dvec_h(col,cmp,lev) = generate_data_3d(col,lev,cmp);
     }
   } 
   for (int lev=0;lev<num_levs;lev++) {
     f3d_h(col,lev) = generate_data_3d(col,lev);
   }
  }
  // Update timestamp
  util::TimeStamp time ({2000,1,1},{0,0,0});
  fm->init_fields_time_stamp(time);

  f2d.sync_to_dev();
  f2dvec.sync_to_dev();
  f3d.sync_to_dev();
  f3dvec.sync_to_dev();

  return fm;
}
/*===================================================================================================*/
// A set of functions that will generate basic data to be used for the test.  The data follows the
// equation:
//   f = cmp*10000.0 + col + 100*lev
// The function is built on calling the full expression but using 0 in cases where a certain dimension
// isn't present.
// Note, because indices are usually 0-based we add 1 to all indices in the actual expression and pass
// -1 in the case of missing dimensions.
Real generate_data_2d(const int col) {
  return generate_data_3d(col,-1,-1);
}

Real generate_data_2d(const int col, const int cmp) {
  return generate_data_3d(col,-1,cmp);
}

Real generate_data_3d(const int col, const int lev) {
  return generate_data_3d(col,lev,-1);
}

Real generate_data_3d(const int col, const int lev, const int cmp) {
  return (col+1) + 100*(lev+1) + 10000*(cmp+1);
}
/*===================================================================================================*/
//  Function to simplify creation of remap files on the fly
void create_remap_file(const ekat::Comm& comm, const std::string& filename, std::shared_ptr<const AbstractGrid> grid,
                       const int n_a, const int n_b,
                       const view_1d<int>& col, const view_1d<int>& row, const view_1d<Real>& S)
{
  // Gather data from all ranks
  const int n_s_local = S.extent(0);
  int  n_s;
  int* n_s_all_ranks = new Int[comm.size()]; 
  comm.all_gather(&n_s_local,n_s_all_ranks,1);
  comm.all_reduce(&n_s_local,&n_s,1,MPI_SUM);
  Int n_s_offset = 0; // The offset for this rank, used for writing the remap to output
  for (int ii=0; ii<comm.rank(); ii++) {
    n_s_offset += n_s_all_ranks[ii];
  }
  scorpio::register_file(filename,scorpio::Write);
  //   - Dimensions
  scorpio::register_dimension(filename,"n_s","n_s",n_s);
  scorpio::register_dimension(filename,"n_a","n_a",n_a);
  scorpio::register_dimension(filename,"n_b","n_b",n_b);
  //   - Variables
  std::string remap_decomp_tag_r = "n_s_real";
  std::string remap_decomp_tag_i = "n_s_int";
  std::vector<std::string> vec_of_remap_dims = {"n_s"};
  scorpio::register_variable(filename,"col","col","unitless",vec_of_remap_dims,"int","int",remap_decomp_tag_i);
  scorpio::register_variable(filename,"row","row","unitless",vec_of_remap_dims,"int","int",remap_decomp_tag_i);
  scorpio::register_variable(filename,"S","S",    "unitless",vec_of_remap_dims,"real","real",remap_decomp_tag_r);
  //   - DOFs
  std::vector<int64_t> var_dof(n_s_local);
  std::iota(var_dof.begin(),var_dof.end(),n_s_offset);
  scorpio::set_dof(filename,"col",var_dof.size(),var_dof.data());
  scorpio::set_dof(filename,"row",var_dof.size(),var_dof.data());
  scorpio::set_dof(filename,"S",var_dof.size(),var_dof.data());
  //   - End definition of file
  scorpio::eam_pio_enddef(filename);
  //   - Write the data to file, note, need to copy to host
  auto col_h = Kokkos::create_mirror_view(col);
  auto row_h = Kokkos::create_mirror_view(row);
  auto S_h   = Kokkos::create_mirror_view(S);
  Kokkos::deep_copy(col_h, col);
  Kokkos::deep_copy(row_h, row);
  Kokkos::deep_copy(S_h,   S);
  scorpio::grid_write_data_array(filename,"col",col_h.data(),n_s_local);
  scorpio::grid_write_data_array(filename,"row",row_h.data(),n_s_local);
  scorpio::grid_write_data_array(filename,"S",S_h.data(),n_s_local);
  // All done writing 
  scorpio::eam_pio_closefile(filename);
}
/*===================================================================================================*/
TEST_CASE("regional_output","io")
{
  run();
} // TEST_CASE regional_output

} // namespace
