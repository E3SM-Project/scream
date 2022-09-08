#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/scream_array_utils.hpp"

#include "ekat/util/ekat_units.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"

#include <numeric>
#include <fstream>

namespace scream
{

// This helper function updates the current output val with a new one,
// according to the "averaging" type, and according to the number of
// model time steps since the last output step.
KOKKOS_INLINE_FUNCTION
void combine (const Real& new_val, Real& curr_val, const OutputAvgType avg_type)
{
  switch (avg_type) {
    case OutputAvgType::Instant:
      curr_val = new_val;
      break;
    case OutputAvgType::Max:
      curr_val = ekat::impl::max(curr_val,new_val);
      break;
    case OutputAvgType::Min:
      curr_val = ekat::impl::min(curr_val,new_val);
      break;
    case OutputAvgType::Average:
      curr_val += new_val;
      break;
    default:
      EKAT_KERNEL_ERROR_MSG ("Unexpected value for m_avg_type. Please, contact developers.\n");
  }
}

AtmosphereOutput::
AtmosphereOutput (const ekat::Comm& comm, const ekat::ParameterList& params,
                  const std::shared_ptr<const fm_type>& field_mgr,
                  const std::shared_ptr<const gm_type>& grids_mgr)
 : m_comm      (comm)
{
  using vos_t = std::vector<std::string>;

  // Figure out what kind of averaging is requested
  auto avg_type = params.get<std::string>("Averaging Type");
  m_avg_type = str2avg(avg_type);
  EKAT_REQUIRE_MSG (m_avg_type!=OutputAvgType::Invalid,
      "Error! Unsupported averaging type '" + avg_type + "'.\n"
      "       Valid options: Instant, Max, Min, Average. Case insensitive.\n");

  set_field_manager (field_mgr);

  // By default, IO is done directly on the field mgr grid
  m_grids_manager = grids_mgr;
  std::shared_ptr<const grid_type> fm_grid, io_grid;
  io_grid = fm_grid = field_mgr->get_grid();
  if (params.isParameter("Field Names")) {
    // This simple parameter list option does *not* allow to remap fields
    // to an io grid different from that of the field manager. In order to
    // use that functionality, you need the full syntax
    m_fields_names = params.get<vos_t>("Field Names");
  } else if (params.isSublist("Fields")){
    const auto& f_pl = params.sublist("Fields");
    const auto& grid_name = io_grid->name(); 
    if (f_pl.isSublist(grid_name)) {
      const auto& pl = f_pl.sublist(grid_name);
      m_fields_names = pl.get<vos_t>("Field Names");

      // Check if the user wants to remap fields on a different grid first
      if (pl.isParameter("IO Grid Name")) {
        io_grid = grids_mgr->get_grid(pl.get<std::string>("IO Grid Name"));
      }
    }
  }

  // Try to set the IO grid (checks will be performed)
  set_grid (io_grid);

  if (io_grid->name()!=fm_grid->name()) {
    // We build a remapper, to remap fields from the fm grid to the io grid
    m_remapper = grids_mgr->create_remapper(fm_grid,io_grid);

    // Register all output fields in the remapper.
    m_remapper->registration_begins();
    for (const auto& fname : m_fields_names) {
      auto f = get_field(fname);
      const auto& src_fid = f.get_header().get_identifier();
      EKAT_REQUIRE_MSG(src_fid.data_type()==DataType::RealType,
          "Error! I/O supports only Real data, for now.\n");
      m_remapper->register_field_from_src(src_fid);
    }
    m_remapper->registration_ends();

    // Now create a new FM on io grid, and create copies of output fields from FM.
    auto io_fm = std::make_shared<fm_type>(io_grid);
    io_fm->registration_begins();
    for (int i=0; i<m_remapper->get_num_fields(); ++i) {
      const auto& tgt_fid = m_remapper->get_tgt_field_id(i);
      io_fm->register_field(FieldRequest(tgt_fid));
    }
    io_fm->registration_ends();

    // Now that fields have been allocated on the io grid, we can bind them in the remapper
    for (const auto& fname : m_fields_names) {
      auto src = get_field(fname);
      auto tgt = io_fm->get_field(fname);
      m_remapper->bind_field(src,tgt);
    }

    // This should never fail, but just in case
    EKAT_REQUIRE_MSG (m_remapper->get_num_fields()==m_remapper->get_num_bound_fields(),
        "Error! Something went wrong while building the scorpio input remapper.\n");

    // Reset the field manager
    set_field_manager(io_fm);
  }

  // Setup I/O structures
  init ();
}

/* ---------------------------------------------------------- */
void AtmosphereOutput::restart (const std::string& filename)
{
  // Create an input stream on the fly, and init averaging data
  ekat::ParameterList res_params("Input Parameters");
  res_params.set<std::string>("Filename",filename);
  res_params.set("Field Names",m_fields_names);

  AtmosphereInput hist_restart (res_params,m_io_grid,m_host_views_1d,m_layouts);
  hist_restart.read_variables();
  hist_restart.finalize();
  for (auto& it : m_host_views_1d) {
    const auto& name = it.first;
    const auto& host = it.second;
    const auto& dev  = m_dev_views_1d.at(name);
    Kokkos::deep_copy(dev,host);
  }
}

void AtmosphereOutput::init()
{
  // Register any diagnostics needed by this output stream
  set_diagnostics();

  for (const auto& var_name : m_fields_names) {
    register_dimensions(var_name);
  }

  // Now that the fields have been gathered register the local views which will be used to determine output data to be written.
  register_views();


} // init
/*-----*/
void AtmosphereOutput::run (const std::string& filename, const bool is_write_step, const int nsteps_since_last_output)
{
  // If we do INSTANT output, but this is not an write step,
  // we can immediately return
  if (not is_write_step and m_avg_type==OutputAvgType::Instant) {
    return;
  }

  using namespace scream::scorpio;

  // If needed, remap fields from their grid to the unique grid, for I/O
  if (m_remapper) {
    m_remapper->remap(true);

    for (int i=0; i<m_remapper->get_num_fields(); ++i) {
      // Need to update the time stamp of the fields on the IO grid,
      // to avoid throwing an exception later
      auto src = m_remapper->get_src_field(i);
      auto tgt = m_remapper->get_tgt_field(i);

      auto src_t = src.get_header().get_tracking().get_time_stamp();
      tgt.get_header().get_tracking().update_time_stamp(src_t);
    }
  }

  // Take care of updating and possibly writing fields.
  for (auto const& name : m_fields_names) {
    // Get all the info for this field.
    const auto  field = get_field(name,true); // If diagnostic, must evaluate it
    const auto& layout = m_layouts.at(name);
    const auto& dims = layout.dims();
    const auto  rank = layout.rank();

    // Safety check: make sure that the field was written at least once before using it.
    EKAT_REQUIRE_MSG (field.get_header().get_tracking().get_time_stamp().is_valid(),
        "Error! Output field '" + name + "' has not been initialized yet\n.");

    const bool is_diagnostic = (m_diagnostics.find(name) != m_diagnostics.end());
    const bool is_aliasing_field_view =
        m_avg_type==OutputAvgType::Instant &&
        field.get_header().get_alloc_properties().get_padding()==0 &&
        field.get_header().get_parent().expired() &&
        not is_diagnostic;

    // Manually update the 'running-tally' views with data from the field,
    // by combining new data with current avg values.
    // NOTE: this is skipped for instant output, if IO view is aliasing Field view.
    auto view_dev = m_dev_views_1d.at(name);
    auto data = view_dev.data();
    KT::RangePolicy policy(0,layout.size());
    const auto extents = layout.extents();

    auto avg_type = m_avg_type;
    // If the dev_view_1d is aliasing the field device view (must be Instant output),
    // then there's no point in copying from the field's view to dev_view
    if (not is_aliasing_field_view) {
      switch (rank) {
        case 1:
        {
          auto new_view_1d = field.get_view<const Real*,Device>();
          auto avg_view_1d = view_Nd_dev<1>(data,dims[0]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int i) {
            combine(new_view_1d(i), avg_view_1d(i),avg_type);
          });
          break;
        }
        case 2:
        {
          auto new_view_2d = field.get_view<const Real**,Device>();
          auto avg_view_2d = view_Nd_dev<2>(data,dims[0],dims[1]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j;
            unflatten_idx(idx,extents,i,j);
            combine(new_view_2d(i,j), avg_view_2d(i,j),avg_type);
          });
          break;
        }
        case 3:
        {
          auto new_view_3d = field.get_view<const Real***,Device>();
          auto avg_view_3d = view_Nd_dev<3>(data,dims[0],dims[1],dims[2]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j,k;
            unflatten_idx(idx,extents,i,j,k);
            combine(new_view_3d(i,j,k), avg_view_3d(i,j,k),avg_type);
          });
          break;
        }
        case 4:
        {
          auto new_view_4d = field.get_view<const Real****,Device>();
          auto avg_view_4d = view_Nd_dev<4>(data,dims[0],dims[1],dims[2],dims[3]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l;
            unflatten_idx(idx,extents,i,j,k,l);
            combine(new_view_4d(i,j,k,l), avg_view_4d(i,j,k,l),avg_type);
          });
          break;
        }
        case 5:
        {
          auto new_view_5d = field.get_view<const Real*****,Device>();
          auto avg_view_5d = view_Nd_dev<5>(data,dims[0],dims[1],dims[2],dims[3],dims[4]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l,m;
            unflatten_idx(idx,extents,i,j,k,l,m);
            combine(new_view_5d(i,j,k,l,m), avg_view_5d(i,j,k,l,m),avg_type);
          });
          break;
        }
        case 6:
        {
          auto new_view_6d = field.get_view<const Real******,Device>();
          auto avg_view_6d = view_Nd_dev<6>(data,dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l,m,n;
            unflatten_idx(idx,extents,i,j,k,l,m,n);
            combine(new_view_6d(i,j,k,l,m,n), avg_view_6d(i,j,k,l,m,n),avg_type);
          });
          break;
        }
        default:
          EKAT_ERROR_MSG ("Error! Field rank (" + std::to_string(rank) + ") not supported by AtmosphereOutput.\n");
      }
    }

    if (is_write_step) {
      if (avg_type==OutputAvgType::Average) {
        // Divide by steps count only when the summation is complete
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int i) {
          data[i] /= nsteps_since_last_output;
        });
      }
      // Bring data to host
      auto view_host = m_host_views_1d.at(name);
      Kokkos::deep_copy (view_host,view_dev);
      grid_write_data_array(filename,name,view_host.data());
    }
  }
} // run

long long AtmosphereOutput::
res_dep_memory_footprint () const {
  long long rdmf = 0;
  if (m_remapper) {
    // The IO is done on a different grid. The FM stored here is
    // not shared with anyone else, so we can safely add its footprint
    for (const auto& it : *m_field_mgr) {
      const auto& fap = it.second->get_header().get_alloc_properties();
      if (fap.is_subfield()) {
        continue;
      }
      rdmf += fap.get_alloc_size();
    }
  }

  for (const auto& fn : m_fields_names) {
    bool is_diagnostic = (m_diagnostics.find(fn) != m_diagnostics.end());
    bool can_alias_field_view =
        m_avg_type==OutputAvgType::Instant && not is_diagnostic &&
        m_field_mgr->get_field(fn).get_header().get_alloc_properties().get_padding()==0 &&
        m_field_mgr->get_field(fn).get_header().get_parent().expired();

    if (not can_alias_field_view) {
      rdmf += m_dev_views_1d.size()*sizeof(Real);
    }
  }

  return rdmf;
}

/* ---------------------------------------------------------- */

void AtmosphereOutput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");
  EKAT_REQUIRE_MSG (field_mgr->get_grid(), "Error! Field manager stores an invalid grid pointer.\n");

  // All good, store it
  m_field_mgr = field_mgr;
}

void AtmosphereOutput::
set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (grid, "Error! Input grid pointer is invalid.\n");
  EKAT_REQUIRE_MSG (grid->is_unique(),
      "Error! I/O only supports grids which are 'unique', meaning that the\n"
      "       map dof_gid->proc_id is well defined.\n");
  EKAT_REQUIRE_MSG (
      (grid->get_global_max_dof_gid()-grid->get_global_min_dof_gid()+1)==grid->get_num_global_dofs(),
      "Error! In order for IO to work, the grid must (globally) have dof gids in interval [gid_0,gid_0+num_global_dofs).\n");

  EKAT_REQUIRE_MSG(m_comm.size()<=grid->get_num_global_dofs(),
      "Error! PIO interface requires the size of the IO MPI group to be\n"
      "       no greater than the global number of columns.\n"
      "       Consider decreasing the size of IO MPI group.\n");

  // The grid is good. Store it.
  m_io_grid = grid;
}

void AtmosphereOutput::register_dimensions(const std::string& name)
{
/* 
 * Checks that the dimensions associated with a specific variable will be registered with IO file.
 * INPUT:
 *   field_manager: is a pointer to the field_manager for this simulation.
 *   name: is a string name of the variable who is to be added to the list of variables in this IO stream.
 */
  using namespace scorpio;

  // Store the field layout
  const auto& fid = get_field(name).get_header().get_identifier();
  const auto& layout = fid.get_layout();
  m_layouts.emplace(name,layout);

  // Now check taht all the dims of this field are already set to be registered.
  for (int i=0; i<layout.rank(); ++i) {
    // check tag against m_dims map.  If not in there, then add it.
    const auto& tags = layout.tags();
    const auto& dims = layout.dims();
    const auto tag_name = get_nc_tag_name(tags[i],dims[i]);
    auto tag_loc = m_dims.find(tag_name);
    auto is_partitioned = m_io_grid->get_partitioned_dim_tag()==tags[i];
    if (tag_loc == m_dims.end()) {
      int tag_len = 0;
      if(tags[i] == m_io_grid->get_partitioned_dim_tag()) {
        // This is the dimension that is partitioned across ranks.
        tag_len = m_io_grid->get_partitioned_dim_global_size();
      } else {
        tag_len = layout.dim(i);
      }
      m_dims.emplace(std::make_pair(get_nc_tag_name(tags[i],dims[i]),tag_len));
    } else {  
      EKAT_REQUIRE_MSG(m_dims.at(tag_name)==dims[i] or is_partitioned,
        "Error! Dimension " + tag_name + " on field " + name + " has conflicting lengths");
    }
  }
} // register_dimensions
/* ---------------------------------------------------------- */
void AtmosphereOutput::register_views()
{
  // Cycle through all fields and register.
  for (auto const& name : m_fields_names) {
    auto field = get_field(name);
    bool is_diagnostic = (m_diagnostics.find(name) != m_diagnostics.end());

    // These local views are really only needed if the averaging time is not 'Instant',
    // to store running tallies for the average operation. However, we create them
    // also for Instant avg_type, for simplicity later on.

    // If we have an 'Instant' avg type, we can alias the 1d views with the
    // views of the field, provided that the field does not have padding,
    // and that it is not a subfield of another field (or else the view
    // would be strided).
    //
    // We also don't want to alias to a diagnostic output since it could share memory
    // with another diagnostic.
    bool can_alias_field_view =
        m_avg_type==OutputAvgType::Instant &&
        field.get_header().get_alloc_properties().get_padding()==0 &&
        field.get_header().get_parent().expired() &&
        not is_diagnostic;

    const auto size = m_layouts.at(name).size();
    if (can_alias_field_view) {
      // Alias field's data, to save storage.
      m_dev_views_1d.emplace(name,view_1d_dev(field.get_internal_view_data<Real,Device>(),size));
      m_host_views_1d.emplace(name,view_1d_host(field.get_internal_view_data<Real,Host>(),size));
    } else {
      // Create a local view.
      m_dev_views_1d.emplace(name,view_1d_dev("",size));
      m_host_views_1d.emplace(name,Kokkos::create_mirror(m_dev_views_1d[name]));

      // Init dev view with an "identity" for avg_type
      switch (m_avg_type) {
        case OutputAvgType::Instant:
          // No averaging
          break;
        case OutputAvgType::Max:
          Kokkos::deep_copy(m_dev_views_1d[name],-std::numeric_limits<Real>::infinity());
          break;
        case OutputAvgType::Min:
          Kokkos::deep_copy(m_dev_views_1d[name],std::numeric_limits<Real>::infinity());
          break;
        case OutputAvgType::Average:
          Kokkos::deep_copy(m_dev_views_1d[name],0);
          break;
        default:
          EKAT_ERROR_MSG ("Unrecognized averaging type.\n");
      }
    }
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::register_variables(const std::string& filename)
{
  using namespace scorpio;

  // Cycle through all fields and register.
  for (auto const& name : m_fields_names) {
    auto field = get_field(name);
    auto& fid  = field.get_header().get_identifier();
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::string io_decomp_tag = "Real";  // Note, for now we only assume REAL variables.  This may change in the future.
    std::vector<std::string> vec_of_dims;
    const auto& layout = fid.get_layout();
    std::string units = to_string(fid.get_units());
    for (int i=0; i<fid.get_layout().rank(); ++i) {
      const auto tag_name = get_nc_tag_name(layout.tag(i), layout.dim(i));
      // Concatenate the dimension string to the io-decomp string
      io_decomp_tag += "-" + tag_name;
      // If tag==CMP, we already attached the length to the tag name
      if (layout.tag(i)!=ShortFieldTagsNames::CMP) {
        io_decomp_tag += "_" + std::to_string(layout.dim(i));
      }
      vec_of_dims.push_back(tag_name); // Add dimensions string to vector of dims.
    }
    // TODO: Do we expect all vars to have a time dimension?  If not then how to trigger? 
    // Should we register dimension variables (such as ncol and lat/lon) elsewhere
    // in the dimension registration?  These won't have time. 
    io_decomp_tag += "-time";
    // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO,
    // may need to delete this line when switching to fully C++/C implementation.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end());
    vec_of_dims.push_back("time");  //TODO: See the above comment on time.

     // TODO  Need to change dtype to allow for other variables.
    // Currently the field_manager only stores Real variables so it is not an issue,
    // but in the future if non-Real variables are added we will want to accomodate that.
    register_variable(filename, name, name, units, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);
  }
} // register_variables
/* ---------------------------------------------------------- */
std::vector<scorpio::offset_t>
AtmosphereOutput::get_var_dof_offsets(const FieldLayout& layout)
{
  std::vector<scorpio::offset_t> var_dof(layout.size());

  // Gather the offsets of the dofs of this variable w.r.t. the *global* array.
  // Since we order the global array based on dof gid, and we *assume* (we actually
  // check this during set_grid) that the grid global gids are in the interval
  // [gid_0, gid_0+num_global_dofs), the offset is simply given by
  // (dof_gid-gid_0)*column_size (for partitioned arrays).
  // NOTE: we allow gid_0!=0, so that we don't have to worry about 1-based numbering
  //       vs 0-based numbering. The key feature is that the global gids are a
  //       contiguous array. The starting point doesn't matter.
  // NOTE: a "dof" in the grid object is not the same as a "dof" in scorpio.
  //       For a SEGrid 3d vector field with (MPI local) layout (nelem,2,np,np,nlev),
  //       scorpio sees nelem*2*np*np*nlev dofs, while the SE grid sees nelem*np*np dofs.
  //       All we need to do in this routine is to compute the offset of all the entries
  //       of the MPI-local array w.r.t. the global array. So long as the offsets are in
  //       the same order as the corresponding entry in the data to be read/written, we're good.
  if (layout.has_tag(ShortFieldTagsNames::COL)) {
    const int num_cols = m_io_grid->get_num_local_dofs();

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

    auto dofs = m_io_grid->get_dofs_gids();
    auto dofs_h = Kokkos::create_mirror_view(dofs);
    Kokkos::deep_copy(dofs_h,dofs);

    // Precompute this *before* the loop, since it involves expensive collectives.
    // Besides, the loop might have different length on different ranks, so
    // computing it inside might cause deadlocks.
    auto min_gid = m_io_grid->get_global_min_dof_gid();
    for (int icol=0; icol<num_cols; ++icol) {
      // Get chunk of var_dof to fill
      auto start = var_dof.begin()+icol*col_size;
      auto end   = start+col_size;

      // Compute start of the column offset, then fill column adding 1 to each entry
      auto gid = dofs_h(icol);
      auto offset = (gid-min_gid)*col_size;
      std::iota(start,end,offset);
    }
  } else if (layout.has_tag(ShortFieldTagsNames::EL)) {
    auto layout2d = m_io_grid->get_2d_scalar_layout();
    const int num_my_elems = layout2d.dim(0);
    const int ngp = layout2d.dim(1);
    const int num_cols = num_my_elems*ngp*ngp;

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

    auto dofs = m_io_grid->get_dofs_gids();
    auto dofs_h = Kokkos::create_mirror_view(dofs);
    Kokkos::deep_copy(dofs_h,dofs);

    // Precompute this *before* the loop, since it involves expensive collectives.
    // Besides, the loop might have different length on different ranks, so
    // computing it inside might cause deadlocks.
    auto min_gid = m_io_grid->get_global_min_dof_gid();
    for (int ie=0,icol=0; ie<num_my_elems; ++ie) {
      for (int igp=0; igp<ngp; ++igp) {
        for (int jgp=0; jgp<ngp; ++jgp,++icol) {
          // Get chunk of var_dof to fill
          auto start = var_dof.begin()+icol*col_size;
          auto end   = start+col_size;

          // Compute start of the column offset, then fill column adding 1 to each entry
          auto gid = dofs_h(icol);
          auto offset = (gid-min_gid)*col_size;
          std::iota(start,end,offset);
    }}}
  } else {
    // This field is *not* defined over columns, so it is not partitioned.
    std::iota(var_dof.begin(),var_dof.end(),0);
  } 

  return var_dof; 
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::set_degrees_of_freedom(const std::string& filename)
{
  using namespace scorpio;
  using namespace ShortFieldTagsNames;

  // Cycle through all fields and set dof.
  for (auto const& name : m_fields_names) {
    auto field = get_field(name);
    const auto& fid  = field.get_header().get_identifier();
    auto var_dof = get_var_dof_offsets(fid.get_layout());
    set_dof(filename,name,var_dof.size(),var_dof.data());
    m_dofs.emplace(std::make_pair(name,var_dof.size()));
  }

  /* TODO: 
   * Gather DOF info directly from grid manager
  */
} // set_degrees_of_freedom
/* ---------------------------------------------------------- */
void AtmosphereOutput::setup_output_file(const std::string& filename)
{
  using namespace scream::scorpio;

  // Register dimensions with netCDF file.
  for (auto it : m_dims) {
    register_dimension(filename,it.first,it.first,it.second);
  }

  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables(filename);

  // Set the offsets of the local dofs in the global vector.
  set_degrees_of_freedom(filename);
}
/* ---------------------------------------------------------- */
// General get_field routine for output.
// This routine will first check if a field is in the local field
// manager.  If not it will next check to see if it is in the list
// of available diagnostics.  If neither of these two options it
// will throw an error.
Field AtmosphereOutput::get_field(const std::string& name, const bool eval_diagnostic)
{
  if (m_field_mgr->has_field(name)) {
    return m_field_mgr->get_field(name);
  } else if (m_diagnostics.find(name) != m_diagnostics.end()) {
    const auto& diag = m_diagnostics[name];
    if (eval_diagnostic) {
      diag->compute_diagnostic();
    }
    return diag->get_diagnostic();
  } else {
    EKAT_ERROR_MSG ("Field " + name + " not found in output field manager or diagnostics list");
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::set_diagnostics()
{
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  for (const auto& fname : m_fields_names) {
    if (!m_field_mgr->has_field(fname)) {
      // Construct a diagnostic by this name
      ekat::ParameterList params;
      params.set<std::string>("Diagnostic Name", fname);
      const auto grid_name = m_io_grid->name(); 
      params.set<std::string>("Grid", grid_name);
      auto diag = diag_factory.create(fname,m_comm,params);
      diag->set_grids(m_grids_manager);
      m_diagnostics.emplace(fname,diag);
    }
  }
  // Set required fields for all diagnostics
  for (const auto& dd : m_diagnostics) {
    const auto& diag = dd.second;
    for (const auto& req : diag->get_required_field_requests()) {
      // Any required fields should be in the field manager, so we can gather
      // the TimeStamp for initialization from the first of these.
      const auto& req_field = get_field(req.fid.name());
      diag->set_required_field(req_field.get_const());
    }

    // Note: this inits with an invalid timestamp. If by any chance we try to
    //       output the diagnostic without computing it, we'll get an error.
    diag->initialize(util::TimeStamp(),RunType::Initial);
  }
}

} // namespace scream
