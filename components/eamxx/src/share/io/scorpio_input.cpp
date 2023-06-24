#include "share/io/scorpio_input.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "share/io/scream_io_utils.hpp"

#include <ekat/util/ekat_string_utils.hpp>

#include <memory>
#include <numeric>

namespace scream
{

AtmosphereInput::
AtmosphereInput (const ekat::ParameterList& params,
                const std::shared_ptr<const fm_type>& field_mgr)
{
  init(params,field_mgr);
}

AtmosphereInput::
AtmosphereInput (const ekat::ParameterList& params,
                 const std::shared_ptr<const grid_type>& grid,
                 const std::map<std::string,view_1d_host>& host_views_1d,
                 const std::map<std::string,FieldLayout>&  layouts)
{
  init (params,grid,host_views_1d,layouts);
}

AtmosphereInput::
AtmosphereInput (const std::string& filename,
                 const std::shared_ptr<const grid_type>& grid,
                 const std::vector<Field>& fields)
{
  // Create param list and field manager on the fly
  ekat::ParameterList params;
  params.set("Filename",filename);
  auto& names = params.get<std::vector<std::string>>("Field Names",{});

  auto fm = std::make_shared<fm_type>(grid);
  for (auto& f : fields) {
    fm->add_field(f);
    names.push_back(f.name());
  }
  init(params,fm);
}

AtmosphereInput::
~AtmosphereInput ()
{
  // In practice, this should always be true, but since we have a do-nothing default ctor,
  // it is possible to create an instance without ever using it. Since finalize would
  // attempt to close the pio file, we need to call it only if init happened.
  if (m_inited_with_views || m_inited_with_fields) {
    finalize();
  }
}

void AtmosphereInput::
init (const ekat::ParameterList& params,
      const std::shared_ptr<const fm_type>& field_mgr)
{
  EKAT_REQUIRE_MSG (not m_inited_with_views,
      "Error! Input class was already inited (with user-provided views).\n");
  EKAT_REQUIRE_MSG (not m_inited_with_fields,
      "Error! Input class was already inited (with fields).\n");

  m_params = params;
  m_fields_names = m_params.get<decltype(m_fields_names)>("Field Names");
  m_filename = m_params.get<std::string>("Filename");

  // Sets the internal field mgr, and possibly sets up the remapper
  set_field_manager(field_mgr);

  // Init scorpio internal structures
  init_scorpio_structures ();

  m_inited_with_fields = true;
}

void AtmosphereInput::
init (const ekat::ParameterList& params,
      const std::shared_ptr<const grid_type>& grid,
      const std::map<std::string,view_1d_host>& host_views_1d,
      const std::map<std::string,FieldLayout>&  layouts)
{
  EKAT_REQUIRE_MSG (not m_inited_with_views,
      "Error! Input class was already inited (with user-provided views).\n");
  EKAT_REQUIRE_MSG (not m_inited_with_fields,
      "Error! Input class was already inited (with fields).\n");

  m_params = params;
  m_filename = m_params.get<std::string>("Filename");

  // Set the grid associated with the input views
  set_grid(grid);

  EKAT_REQUIRE_MSG (host_views_1d.size()==layouts.size(),
      "Error! Input host views and layouts maps has different sizes.\n"
      "       host_views_1d size: " + std::to_string(host_views_1d.size()) + "\n"
      "       layouts size: " + std::to_string(layouts.size()) + "\n");

  m_layouts = layouts;
  m_host_views_1d = host_views_1d;

  // Loop over one of the two maps, store key in m_fields_names,
  // and check that the two maps have the same keys
  for (const auto& it : m_layouts) {
    m_fields_names.push_back(it.first);
    EKAT_REQUIRE_MSG (m_host_views_1d.count(it.first)==1,
        "Error! Input layouts and views maps do not store the same keys.\n");
  }

  // Init scorpio internal structures
  init_scorpio_structures ();

  m_inited_with_views = true;
}

/* ---------------------------------------------------------- */

void AtmosphereInput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");
  EKAT_REQUIRE_MSG (field_mgr->get_grid(), "Error! Field manager stores an invalid grid pointer.\n");

  // If resetting a field manager we want to check that the layouts of all fields are the same.
  if (m_field_mgr) {
    for (auto felem = m_field_mgr->begin(); felem != m_field_mgr->end(); felem++) {
      auto name = felem->first;
      auto field_curr = m_field_mgr->get_field(name);
      auto field_new  = field_mgr->get_field(name);
      // Check Layouts
      auto lay_curr   = field_curr.get_header().get_identifier().get_layout();
      auto lay_new    = field_new.get_header().get_identifier().get_layout();
      EKAT_REQUIRE_MSG(lay_curr==lay_new,"ERROR!! AtmosphereInput::set_field_manager - setting new field manager which has different layout for field " << name <<"\n"
		      << "    Old Layout: " << to_string(lay_curr) << "\n"
		      << "    New Layout: " << to_string(lay_new) << "\n");
    }
  }

  m_field_mgr = field_mgr;

  // Store grid and fm
  set_grid(m_field_mgr->get_grid());

  // Init fields specs
  for (auto const& name : m_fields_names) {
    auto f = m_field_mgr->get_field(name);
    const auto& fh  = f.get_header();
    const auto& fap = fh.get_alloc_properties();
    const auto& fid = fh.get_identifier();
    const auto& fl  = fid.get_layout();

    // Store the layout
    m_layouts.emplace(name,fl);

    // If we can alias the field's host view, do it.
    // Otherwise, create a temporary.
    bool can_alias_field_view = fh.get_parent().expired() && fap.get_padding()==0;
    if (can_alias_field_view) {
      auto data = f.get_internal_view_data<Real,Host>();
      m_host_views_1d[name] = view_1d_host(data,fl.size());
    } else {
      // We have padding, or the field is a subfield (or both).
      // Either way, we need a temporary view.
      m_host_views_1d[name] = view_1d_host("",fl.size());
    }
  }
}


/* ---------------------------------------------------------- */
void AtmosphereInput::
set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (grid, "Error! Input grid pointer is invalid.\n");
  const bool skip_grid_chk = m_params.get<bool>("Skip_Grid_Checks",false);
  if (!skip_grid_chk) {
    EKAT_REQUIRE_MSG (grid->is_unique(),
        "Error! I/O only supports grids which are 'unique', meaning that the\n"
        "       map dof_gid->proc_id is well defined.\n");
    EKAT_REQUIRE_MSG (
        (grid->get_global_max_dof_gid()-grid->get_global_min_dof_gid()+1)==grid->get_num_global_dofs(),
        "Error! IO requires DOF gids to (globally)  be in interval [gid_0,gid_0+num_global_dofs).\n"
        "   - global min GID : " + std::to_string(grid->get_global_min_dof_gid()) + "\n"
        "   - global max GID : " + std::to_string(grid->get_global_max_dof_gid()) + "\n"
        "   - num global dofs: " + std::to_string(grid->get_num_global_dofs()) + "\n");
  }

  EKAT_REQUIRE_MSG(grid->get_comm().size()<=grid->get_num_global_dofs(),
      "Error! PIO interface requires the size of the IO MPI group to be\n"
      "       no greater than the global number of columns.\n"
      "       Consider decreasing the size of IO MPI group.\n");

  // The grid is good. Store it.
  m_io_grid = grid;
}

/* ---------------------------------------------------------- */
// Note: The (zero-based) time_index argument provides a way to control which
//       time step to read input from in the file.  If a negative number is
//       provided the routine will read input at the last time level set by
//       running eam_update_timesnap.
void AtmosphereInput::read_variables (const int time_index)
{
  EKAT_REQUIRE_MSG (m_inited_with_views || m_inited_with_fields,
      "Error! Scorpio structures not inited yet. Did you forget to call 'init(..)'?\n");

  for (auto const& name : m_fields_names) {

    // Read the data
    auto v1d = m_host_views_1d.at(name);
    scorpio::read_var(m_filename,name,v1d.data(),time_index);

    // If we have a field manager, make sure the data is correctly
    // synced to both host and device views of the field.
    if (m_field_mgr) {

      auto f = m_field_mgr->get_field(name);
      const auto& fh  = f.get_header();
      const auto& fl  = fh.get_identifier().get_layout();
      const auto& fap = fh.get_alloc_properties();

      // Check if the stored 1d view is sharing the data ptr with the field
      const bool can_alias_field_view = fh.get_parent().expired() && fap.get_padding()==0;

      // If the 1d view is a simple reshape of the field's Host view data,
      // then we're already done. Otherwise, we need to manually copy.
      if (not can_alias_field_view) {
        // Get the host view of the field properly reshaped, and deep copy
        // from temp_view (properly reshaped as well).
        auto rank = fl.rank();
        auto view_1d = m_host_views_1d.at(name);
        switch (rank) {
          case 1:
            {
              // No reshape needed, simply copy
              auto dst = f.get_view<Real*,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                dst(i) = view_1d(i);
              }
              break;
            }
          case 2:
            {
              // Reshape temp_view to a 2d view, then copy
              auto dst = f.get_view<Real**,Host>();
              auto src = view_Nd_host<2>(view_1d.data(),fl.dim(0),fl.dim(1));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  dst(i,j) = src(i,j);
              }}
              break;
            }
          case 3:
            {
              // Reshape temp_view to a 3d view, then copy
              auto dst = f.get_view<Real***,Host>();
              auto src = view_Nd_host<3>(view_1d.data(),fl.dim(0),fl.dim(1),fl.dim(2));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  for (int k=0; k<fl.dim(2); ++k) {
                    dst(i,j,k) = src(i,j,k);
              }}}
              break;
            }
          case 4:
            {
              // Reshape temp_view to a 4d view, then copy
              auto dst = f.get_view<Real****,Host>();
              auto src = view_Nd_host<4>(view_1d.data(),fl.dim(0),fl.dim(1),fl.dim(2),fl.dim(3));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  for (int k=0; k<fl.dim(2); ++k) {
                    for (int l=0; l<fl.dim(3); ++l) {
                      dst(i,j,k,l) = src(i,j,k,l);
              }}}}
              break;
            }
          case 5:
            {
              // Reshape temp_view to a 5d view, then copy
              auto dst = f.get_view<Real*****,Host>();
              auto src = view_Nd_host<5>(view_1d.data(),fl.dim(0),fl.dim(1),fl.dim(2),fl.dim(3),fl.dim(4));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  for (int k=0; k<fl.dim(2); ++k) {
                    for (int l=0; l<fl.dim(3); ++l) {
                      for (int m=0; m<fl.dim(4); ++m) {
                        dst(i,j,k,l,m) = src(i,j,k,l,m);
              }}}}}
              break;
            }
          case 6:
            {
              // Reshape temp_view to a 6d view, then copy
              auto dst = f.get_view<Real******,Host>();
              auto src = view_Nd_host<6>(view_1d.data(),fl.dim(0),fl.dim(1),fl.dim(2),fl.dim(3),fl.dim(4),fl.dim(5));
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  for (int k=0; k<fl.dim(2); ++k) {
                    for (int l=0; l<fl.dim(3); ++l) {
                      for (int m=0; m<fl.dim(4); ++m) {
                        for (int n=0; n<fl.dim(5); ++n) {
                          dst(i,j,k,l,m,n) = src(i,j,k,l,m,n);
              }}}}}}
              break;
            }
          default:
            EKAT_ERROR_MSG ("Error! Unexpected field rank (" + std::to_string(rank) + ").\n");
        }
      }

      // Sync to device
      f.sync_to_dev();
    }
  }
} 

/* ---------------------------------------------------------- */
void AtmosphereInput::finalize() 
{
  // Protect in case of double call to finalize
  if (m_inited_with_views or m_inited_with_fields) {
    scorpio::release_file(m_filename);

    m_field_mgr = nullptr;
    m_io_grid   = nullptr;

    m_host_views_1d.clear();
    m_layouts.clear();

    m_inited_with_views = false;
    m_inited_with_fields = false;
  }
} // finalize

/* ---------------------------------------------------------- */
void AtmosphereInput::init_scorpio_structures() 
{
  scorpio::register_file(m_filename,scorpio::Read);

  // Register variables with netCDF file.
  register_variables();
  set_degrees_of_freedom();
}

/* ---------------------------------------------------------- */
void AtmosphereInput::register_variables()
{
  // Register each variable in IO stream with the SCORPIO interface.
  // This allows SCORPIO to lookup vars in the nc file with the correct
  // dof decomposition across different ranks.

  // Cycle through all fields
  for (auto const& name : m_fields_names) {
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    auto vec_of_dims   = get_vec_of_dims(m_layouts.at(name));

    // Register the variable
    // TODO  Need to change dtype to allow for other variables. 
    //  Currently the field_manager only stores Real variables so it is not an issue,
    //  but in the future if non-Real variables are added we will want to accomodate that.
    EKAT_REQUIRE_MSG (scorpio::has_variable(m_filename,name,vec_of_dims),
        "Error! Input file does not store a required variable.\n"
        " - filename: " + m_filename + "\n"
        " - varname : " + name + "\n");
  }
}

/* ---------------------------------------------------------- */
std::vector<std::string>
AtmosphereInput::get_vec_of_dims(const FieldLayout& layout)
{
  // Given a set of dimensions in field tags, extract a vector of strings
  // for those dimensions to be used with IO
  using namespace ShortFieldTagsNames;
  std::vector<std::string> dims_names;
  dims_names.reserve(layout.rank());
  for (int i=0; i<layout.rank(); ++i) {
    const FieldTag t = layout.tag(i);
    if (t==CMP) {
      dims_names.push_back("dim" + std::to_string(layout.dim(i)));
    } else {
      dims_names.push_back(m_io_grid->get_dim_name(t));
    }
  }

  return dims_names;
}

/* ---------------------------------------------------------- */
void AtmosphereInput::set_degrees_of_freedom()
{
  using namespace ShortFieldTagsNames;

  // First, set the offsets for the decomposed dim
  const auto decomp_tag  = m_io_grid->get_partitioned_dim_tag();

  bool has_decomposed_layouts = false;
  for (const auto& it : m_layouts) {
    if (it.second.has_tag(decomp_tag)) {
      has_decomposed_layouts = true;
      break;
    }
  }

  // If none of the input vars are decomposed on this grid,
  // then there's nothing to do here
  if (not has_decomposed_layouts) {
    return;
  }

  const auto decomp_name = m_io_grid->name();
  const int local_dim = m_io_grid->get_partitioned_dim_local_size();
  auto decomp_dim = m_io_grid->get_dim_name(decomp_tag);
  if (decomp_tag==COL) {
    // For PointGrid, we can use the dofs to compute the offsets on the partition dim
    std::vector<scorpio::offset_t> offsets(local_dim);
    auto dofs = m_io_grid->get_dofs_gids().get_view<const gid_t*,Host>();
    auto min_dof = m_io_grid->get_global_min_dof_gid();
    for (int idof=0; idof<local_dim; ++idof) {
      offsets[idof] = dofs[idof] - min_dof;
    }
    scorpio::set_dim_decomp(m_filename,decomp_dim,offsets,decomp_name);
  } else if (decomp_tag==EL) {
    // For SEGrid, the dofs are not on the partitioned dim (dofs are nelem*np*np, while
    // the partitioned dim is just the elements). To assign offsets along the decomp dim,
    // simply assing an elem gid in a contiguous way across MPI ranks.
    int start;
    m_io_grid->get_comm().scan(&local_dim,&start,1,MPI_SUM);
    start -= local_dim;

    scorpio::set_dim_decomp(m_filename,decomp_dim,start,local_dim,decomp_name);
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported decomposed dimension tag.\n"
        " - io grid  : " << m_io_grid->name() + "\n"
        " - field tag: " << e2str(decomp_tag) + "\n");
  }

  // Next, loop over all vars, and if the layout include decomp_tag,
  // set set the decomposition in PIO. Notice that, for non-decomposed tags,
  // there is nothing to do, since read_var will just call PIOc_get_var.
  for (const auto& it : m_layouts) {
    const auto& varname = it.first;

    // NOTE: if var is not decomposed, this function doesn't do much
    scorpio::set_var_decomp(m_filename,varname,decomp_dim,decomp_name);
  }
} // set_degrees_of_freedom

} // namespace scream
