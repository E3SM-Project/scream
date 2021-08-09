#include "share/io/scorpio_input.hpp"

#include <numeric>
#include <string>

namespace scream
{

// ====================== IMPLEMENTATION ===================== //
/*
 *  pull_input calls the three typical init, run and finalize routines in order
 *  and makes the input immediately availalbe. Note, init, run and finalize are
 *  also outward facing but do not need those respective sections of the AD, 
 *  i.e. during ad.init, ad.run and ad.finalize.
 */
/* ---------------------------------------------------------- */
void AtmosphereInput::
pull_input(const std::string& filename, const std::string& var_name,
           const std::vector<std::string>& var_dims,
           const bool has_columns, const std::vector<int>& dim_lens,
           const int padding, const grid_ptr_type& grid, Real* data)
{
  using namespace scream::scorpio;
  
  auto gids_dev = m_grid_mgr->get_grid(m_grid_name)->get_dofs_gids();
  m_gids_host = Kokkos::create_mirror_view( gids_dev );
  Kokkos::deep_copy(m_gids_host,gids_dev);

  // Open the file in PIO
  register_infile(filename);

  // Register the variable information in PIO
  std::string io_decomp_tag = get_io_decomp(var_dims);
  get_variable(filename, var_name, var_name, var_dims.size(), var_dims, PIO_REAL, io_decomp_tag);

  // Determine the degree's of freedom for this variable on this rank, and register with PIO
  int data_length = 1;
  for (auto& dim : dim_lens) {
    data_length *= dim;
  }

  std::vector<Int> var_dof = get_var_dof_offsets(data_length, has_columns);
  set_dof(filename,var_name,var_dof.size(),var_dof.data());
  set_decomp(filename);
  grid_read_data_array(filename,var_name,dim_lens,data_length,padding,data);
  eam_pio_closefile(filename);  
}

/* ---------------------------------------------------------- */
AtmosphereInput::view_1d_host AtmosphereInput::pull_input(const std::string& name)
{
/*  Run through the sequence of opening the file, reading input and then closing the file.  
 *  Overloaded case to deal with just one output and to not put output into field manager.
 *  This is an inefficient way to handle this special case because the same file would need
 *  to be opened and closed multiple times if there is more than one variable.  TODO: Make this
 *  a lot more efficient, most likely will need to overhaul the "pull_input" paradigrid_mgr.
 */
  using namespace scream::scorpio;
  if (name=="avg_count") {
    view_1d_host l_view("",1);
    grid_read_data_array(m_filename,name,m_dofs_sizes.at(name),l_view.data());
    return l_view;
  } else {
    // Read into the host view of the field
    read_input(name);

    // Create a 1d  view of the host view in the field
    auto field = m_field_mgr->get_field(name);
    const auto& fl  = field.get_header().get_identifier().get_layout();
    const auto& fap = field.get_header().get_alloc_properties();
    const int last_dim = fap.get_last_extent();
    view_1d_host l_view("",fap.get_num_scalars());
    switch (fl.rank()) {
      case 1:
        Kokkos::deep_copy(l_view,field.get_view<Real*,Host>());
        break;
      case 2:
        Kokkos::deep_copy(view_ND_host<2>(l_view.data(),fl.dim(0),last_dim),field.get_view<Real**,Host>());
        break;
      case 3:
        Kokkos::deep_copy(view_ND_host<3>(l_view.data(),fl.dim(0),fl.dim(1),last_dim),field.get_view<Real***,Host>());
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
    }

    return l_view;
  }
} 

/* ---------------------------------------------------------- */
void AtmosphereInput::pull_input()
{
/*  Run through the sequence of opening the file, reading input and then closing the file.  */
  using namespace scream::scorpio;

  init();

  for (auto const& name : m_fields_names) {
    // Read from file, into the host copy of the field
    read_input(name);
  }
  finalize();
} 

/* ---------------------------------------------------------- */
void AtmosphereInput::init() 
{
/* 
 * Call all of the necessary SCORPIO routines to open the file, gather the variables and dimensions,
 * and set the degrees-of-freedom for reading.
 * Organize local structures responsible for writing over field manager data.
 */
  using namespace scorpio;
  using namespace ShortFieldTagsNames;

  m_filename = m_params.get<std::string>("FILENAME");

  // Parse the parameters that controls this input instance.
  m_grid_name = m_params.get<std::string>("GRID");

  // If this is RHIST type make sure its noted
  if (m_params.isParameter("RHIST")) {
    m_is_rhist = m_params.get<bool>("RHIST");
  }

  // Check setup against information from grid manager:
  if (!m_grid_set) {
    m_grid = m_grid_mgr->get_grid(m_grid_name);
    m_grid_set = true;
  }
  if (!m_grid_set) {
  EKAT_REQUIRE_MSG(m_grid->get_2d_scalar_layout().tags().front()==COL,
      "Error with input grid! scorpio_input.hpp class only supports input on a Physics based grid for now.\n");
    m_grid = m_grid_mgr->get_grid(m_grid_name);
    m_grid_set = true;
  }
  EKAT_REQUIRE_MSG(m_grid->get_2d_scalar_layout().tags().front()==COL,
      "Error with input grid! scorpio_input.hpp class only supports input on a Physics based grid for now.\n");

  auto gids_dev = m_grid->get_dofs_gids();
  m_gids_host = Kokkos::create_mirror_view( gids_dev );
  Kokkos::deep_copy(m_gids_host,gids_dev); 

  // Note, only the total number of columns is distributed over MPI ranks, need to sum over all procs this size to properly register COL dimension.
  int total_dofs;
  int local_dofs = m_gids_host.size();
  MPI_Allreduce(&local_dofs, &total_dofs, 1, MPI_INT, MPI_SUM, m_comm.mpi_comm());
  EKAT_REQUIRE_MSG(m_comm.size()<=total_dofs,"Error, PIO interface only allows for the IO comm group size to be less than or equal to the total # of columns in grid.  Consider decreasing size of IO comm group.\n");

  // Create map of fields in this input with the field_identifier in the field manager.
  auto& var_params = m_params.sublist("FIELDS");
  for (int var_i=0; var_i<var_params.get<Int>("Number of Fields");++var_i) {
    /* Determine the variable name */
    std::string var_name = var_params.get<std::string>(ekat::strint("field",var_i+1));
    m_fields_names.push_back(var_name);
  }

  // Register new netCDF file for input.
  register_infile(m_filename);

  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables();
  set_degrees_of_freedom();

  // Finish the definition phase for this file.
  set_decomp  (m_filename); 
} // init

/* ---------------------------------------------------------- */
void AtmosphereInput::finalize() 
{
/* Cleanup by closing the input file */
  scorpio::eam_pio_closefile(m_filename);
} // finalize

/* ---------------------------------------------------------- */
void AtmosphereInput::register_variables()
{
/* Register each variable in IO stream with the SCORPIO interface.  See scream_scorpio_interface.* for details.
 * This is necessary for the SCORPIO routines to be able to lookup variables in the io file with the appropriate
 * degrees of freedom assigned to each core and the correct io decomposition.
 */
  using namespace scorpio;

  // Cycle through all fields and register.
  for (auto const& name : m_fields_names) {
    auto field = m_field_mgr->get_field(name);
    auto& fid  = field.get_header().get_identifier();

    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::vector<std::string> vec_of_dims = get_vec_of_dims(fid.get_layout());
    std::string io_decomp_tag           = get_io_decomp(vec_of_dims);
    get_variable(m_filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);
    // TODO  Need to change dtype to allow for other variables. 
    //  Currently the field_manager only stores Real variables so it is not an issue,
    //  but in the future if non-Real variables are added we will want to accomodate that.
    //TODO: Should be able to simply inquire fromt he netCDF the dimensions for each variable.
  }
  if (m_is_rhist) {
    get_variable(m_filename,"avg_count","avg_count",1,{"cnt"}, PIO_INT, "cnt");
  }
} // register_variables
/* ---------------------------------------------------------- */
void AtmosphereInput::
read_input(const std::string& name)
{
  using namespace scorpio;

  auto field = m_field_mgr->get_field(name);
  const auto& fh  = field.get_header();
  const auto& fl  = fh.get_identifier().get_layout();
  const auto& fap = fh.get_alloc_properties();

  // Get all the info for this field.
  auto l_dims = fl.dims();
  const auto padding = fap.get_padding();

  // Strategy: create a temp contiguous view (except possibly for paddin),
  // and use it for reading from file. Then, deep copy back to the input view.
  using field_type = decltype(field);
  using RT         = typename field_type::RT;

  // Use a 1d view of correct size for scorpio reading
  view_1d_host temp_view("",fap.get_num_scalars());
  grid_read_data_array(m_filename,name,l_dims,m_dofs_sizes.at(name),
                       padding,temp_view.data());

  // Get the host view of the field properly reshaped, and deep copy
  // from temp_view (properly reshaped as well)
  auto rank   = fl.rank();
  switch (rank) {
    case 1:
      {
        // Easy: can deep copy from the 1d view directly
        Kokkos::deep_copy(field.get_view<RT*,Host>(),temp_view);
        break;
      }
    case 2:
      {
        // Reshape temp_view to a 2d view, then copy
        Kokkos::deep_copy(field.get_view<RT**,Host>(),
                          view_ND_host<2>(temp_view.data(),fl.dim(0),fl.dim(1)+padding));
        break;
      }
    case 3:
      {
        // Reshape temp_view to a 3d view, then copy
        Kokkos::deep_copy(field.get_view<RT***,Host>(),
                          view_ND_host<3>(temp_view.data(),fl.dim(0),fl.dim(1),fl.dim(2)+padding));
        break;
      }
    default:
      EKAT_ERROR_MSG (
          "Error! Rank-" + std::to_string(rank) + " field not yet supported in AtmosphereInput.\n");
  }

  // Sync to device
  field.sync_to_dev();
}

/* ---------------------------------------------------------- */
std::vector<std::string>
AtmosphereInput::get_vec_of_dims(const FieldLayout& layout)
{
  // Given a set of dimensions in field tags, extract a vector of strings
  // for those dimensions to be used with IO
  std::vector<std::string> dims_names;
  dims_names.reserve(layout.tags().size());
  for (int i=0; i<layout.rank(); ++i) {
    dims_names.push_back(scorpio::get_nc_tag_name(layout.tag(i),layout.dim(i)));
  }

  // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO,
  // may need to delete this line when switching to fully C++/C implementation.
  std::reverse(dims_names.begin(),dims_names.end());

  return dims_names;
}

/* ---------------------------------------------------------- */
std::string AtmosphereInput::get_io_decomp(const std::vector<std::string>& dims_names)
{
/* Given a vector of field dimensions, create a unique decomp string to register with I/O/
 * Note: We are hard-coding for only REAL input here.  TODO: would be to allow for other dtypes
 */
  std::string io_decomp_tag = "Real";
  for (const auto& dim : dims_names) {
    io_decomp_tag += "-" + dim;
  }

  return io_decomp_tag;
}

/* ---------------------------------------------------------- */
void AtmosphereInput::set_degrees_of_freedom()
{
/* 
 * Use information from the grids manager to determine which indices for each field are associated
 * with this computational core.
 */
  using namespace scorpio;
  using namespace ShortFieldTagsNames;

  // Cycle through all fields and set dof.
  for (auto const& name : m_fields_names) {
    auto field = m_field_mgr->get_field(name);
    auto& fid  = field.get_header().get_identifier();

    // Given dof_len and n_dim_len it should be possible to create an integer array
    // of "global input indices" for this field and this rank. For every column (i.e. gid)
    // the PIO indices would be (gid * n_dim_len),...,( (gid+1)*n_dim_len - 1).
    const bool has_col_tag = fid.get_layout().has_tag(COL);
    std::vector<Int> var_dof = get_var_dof_offsets(fid.get_layout().size(), has_col_tag);
    set_dof(m_filename,name,var_dof.size(),var_dof.data());
    m_dofs_sizes.emplace(std::make_pair(name,var_dof.size()));
  }

  if (m_is_rhist) {
    Int var_dof[1] = {0};
    set_dof(m_filename,"avg_count",1,var_dof);
    m_dofs_sizes.emplace(std::make_pair("avg_count",1));
  }
  // TODO: Gather DOF info directly from grid manager
} // set_degrees_of_freedom

/* ---------------------------------------------------------- */
std::vector<Int> AtmosphereInput::get_var_dof_offsets(const int dof_len, const bool has_cols)
{
  std::vector<Int> var_dof(dof_len);

  // Gather the offsets of the dofs of this variable w.r.t. the *global* array.
  // These are not the dofs global ids (which are just labels, and can be whatever,
  // and in fact are not even contiguous when Homme generates the dof gids).
  // So, if the returned vector is {2,3,4,5}, it means that the 4 dofs on this rank
  // correspond to the 3rd,4th,5th, and 6th dofs globally.
  if (has_cols) {
    const int num_cols = m_gids_host.size();

    // Note: col_size might be *larger* than the number of vertical levels, or even smalle.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    Int col_size = dof_len/num_cols;

    // Compute the number of columns owned by all previous ranks.
    Int offset = 0;
    m_comm.scan_sum(&num_cols,&offset,1);

    // Compute offsets of all my dofs
    std::iota(var_dof.begin(), var_dof.end(), offset*col_size);
  } else {
    // This field is *not* defined over columns, so it is not partitioned.
    std::iota(var_dof.begin(),var_dof.end(),0);
  } 

  return var_dof; 
}
/* ---------------------------------------------------------- */

} // namespace scream
