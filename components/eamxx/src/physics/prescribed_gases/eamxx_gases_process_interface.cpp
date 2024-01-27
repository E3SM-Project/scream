#include "eamxx_gases_process_interface.hpp"
#include "share/util/scream_universal_constants.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/grid/remap/do_nothing_remapper.hpp"

namespace scream
{

// =========================================================================================
PrescribedGases::PrescribedGases (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  m_datafiles  = m_params.get<std::vector<std::string>>("gases_filename");
  m_fields_gases = m_params.get<std::vector<std::string>>("gases_fields");
  // If we are doing horizontal refine-remapping, we need to get the mapfile from user
  m_refine_remap_file = m_params.get<std::string>(
      "gases_refine_remap_mapfile", "no-file-given");
  auto src_pres_type = m_params.get<std::string>("source_pressure_type","TIME_DEPENDENT_3D_PROFILE");
  if (src_pres_type=="TIME_DEPENDENT_3D_PROFILE") {
    m_src_pres_type = TIME_DEPENDENT_3D_PROFILE;
  } else if (src_pres_type=="STATIC_1D_VERTICAL_PROFILE") {
    m_src_pres_type = STATIC_1D_VERTICAL_PROFILE;
    // Check for a designated source pressure file, default to first data source if not given.
    m_static_vertical_pressure_file = m_params.get<std::string>("source_pressure_file",m_datafiles[0]);
  } else {
    EKAT_ERROR_MSG("ERROR! PrescribedGases::parameter_list - unsupported source_pressure_type provided.  Current options are [TIME_DEPENDENT_3D_PROFILE,STATIC_1D_VERTICAL_PROFILE].  Please check");
  }

  // TODO: Add some warning messages here.
  // 1. if m_fields_gases is empty or =NONE then we will skip altogether.
  // 2. Populate with constants for fields not read from file?
}

// =========================================================================================
void PrescribedGases::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };

  constexpr int ps = 1;
  auto Q = kg/kg;
  auto molmol = mol/mol;
  molmol.set_string("mol/mol");
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);

  // Add gases to FM
  for (auto name : m_fields_gases) {
      add_field<Computed>(name, scalar3d_layout_mid, molmol, grid_name, ps);
  }

  //Now need to read in the file
  if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
    m_num_src_levs = scorpio::get_dimlen(m_datafiles[0],"lev");
  } else {
    m_num_src_levs = scorpio::get_dimlen(m_static_vertical_pressure_file,"lev");
  }

  /* Check for consistency between gas file, map file, and remapper */

  // Number of columns globally
  auto m_num_cols_global = m_grid->get_num_global_dofs(); 

  // Get the information from the first data file
  int num_cols_src = scorpio::get_dimlen(m_datafiles[0],"ncol");

  if (num_cols_src != m_num_cols_global) {
    // If differing cols, check if remap file is provided
    EKAT_REQUIRE_MSG(m_refine_remap_file != "no-file-given",
                     "Error! PrescribedGases::set_grids - the number of columns in the data file "
                     << std::to_string(num_cols_src) << " does not match the number of columns in the "
                     << "model grid " << std::to_string(m_num_cols_global) << ".  Please check the "
                     << "gases data file and/or the model grid.");
    // If remap file is provided, check if it is consistent with the gases data file
    // First get the data from the mapfile
    int num_cols_remap_a = scorpio::get_dimlen(m_refine_remap_file,"n_a");
    int num_cols_remap_b = scorpio::get_dimlen(m_refine_remap_file,"n_b");
    // Then, check if n_a (source) and n_b (target) are consistent
    EKAT_REQUIRE_MSG(num_cols_remap_a == num_cols_src,
                     "Error! PrescribedGases::set_grids - the number of columns in the data file "
                     << std::to_string(num_cols_src) << " does not match the number of columns in the "
                     << "mapfile " << std::to_string(num_cols_remap_a) << ".  Please check the "
                     << "gases data file and/or the mapfile.");
    EKAT_REQUIRE_MSG(num_cols_remap_b == m_num_cols_global,
                     "Error! PrescribedGases::set_grids - the number of columns in the model grid "
                     << std::to_string(m_num_cols_global) << " does not match the number of columns in the "
                     << "mapfile " << std::to_string(num_cols_remap_b) << ".  Please check the "
                     << "model grid and/or the mapfile.");
    // If we get here, we are good to go!
    m_refine_remap = true;
  } else {
    // If the number of columns is the same, we don't need to do any remapping,
    // but print a warning if the user provided a mapfile
    if (m_refine_remap_file != "no-file-given") {
      std::cout << "Warning! PrescribedGases::set_grids - the number of columns in the data file "
                << std::to_string(num_cols_src) << " matches the number of columns in the "
                << "model grid " << std::to_string(m_num_cols_global) << ".  The mapfile "
                << m_refine_remap_file << " will NOT be used.  Please check the "
                << "gases data file and/or the model grid." << std::endl;
    }
    // Set m_refine_remap to false
    m_refine_remap = false;
  }
}

// =============================================================================================================
void PrescribedGases::initialize_impl (const RunType /* run_type */)
{
  using namespace ShortFieldTagsNames;
  // Set up pointers for grids
  // external grid: from source data
  std::shared_ptr<scream::AbstractGrid> grid_ext;
  // temporary grid: after vertical interpolation
  std::shared_ptr<scream::AbstractGrid> grid_tmp;
  // internal grid: after horizontal interpolation
  std::shared_ptr<scream::AbstractGrid> grid_int;

  // Initialize the refining remapper stuff at the outset,
  // because we need to know the grid information;
  // for now, we are doing the horizontal interpolation last,
  // so we use the m_grid (model physics) as the target grid
  grid_int = m_grid->clone(m_grid->name(), true);
  if (m_refine_remap) {
    // P2P remapper
    m_refine_remapper =
        std::make_shared<RefiningRemapperP2P>(grid_int, m_refine_remap_file);
    // Get grid from remapper, and clone it
    auto grid_ext_const = m_refine_remapper->get_src_grid();
    grid_ext = grid_ext_const->clone(grid_ext_const->name(), true);
    // Finally, grid_ext may have different levels
    grid_ext->reset_num_vertical_lev(m_num_src_levs);
  } else {
    // We set up a DoNothingRemapper, which will do nothing
    m_refine_remapper = std::make_shared<DoNothingRemapper>(grid_int, grid_int);
    // We clone physics grid, but maybe we have different levels
    grid_ext = m_grid->clone(m_grid->name(), true);
    grid_ext->reset_num_vertical_lev(m_num_src_levs);
  }
  // The temporary grid is the external grid, but with
  // the same number of levels as the internal (physics) grid
  grid_tmp = grid_ext->clone(grid_ext->name(), true);
  grid_tmp->reset_num_vertical_lev(m_num_levs);

  // Declare the layouts for the helper fields (int: internal)
  FieldLayout scalar2d_layout_mid { {LEV}, {m_num_levs} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };

  // Get the number of external cols on current rank
  auto m_num_cols_ext = grid_ext->get_num_local_dofs();
  // Declare the layouts for the helper fields (tmp: temporary))
  FieldLayout scalar2d_layout_mid_tmp { {LEV}, {m_num_levs}};
  FieldLayout scalar3d_layout_mid_tmp { {COL,LEV}, {m_num_cols_ext, m_num_levs} };
  // Declare the layouts for the helper fields (ext: external)
  FieldLayout scalar2d_layout_mid_ext { {LEV}, {m_num_src_levs}};
  FieldLayout scalar3d_layout_mid_ext { {COL,LEV}, {m_num_cols_ext, m_num_src_levs} };

  // Initialize the time interpolator
  m_time_interp = util::TimeInterpolation(grid_ext, m_datafiles);

  constexpr int ps = SCREAM_PACK_SIZE;
  // To be extra careful, this should be the ext_grid
  const auto& grid_ext_name = grid_ext->name();
  if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
    auto pmid_ext = create_helper_field("p_mid_ext", scalar3d_layout_mid_ext, grid_ext_name, ps);
    m_time_interp.add_field(pmid_ext.alias("p_mid"),true);
  } else if (m_src_pres_type == STATIC_1D_VERTICAL_PROFILE) {
    // Load p_levs from source data file
    ekat::ParameterList in_params;
    in_params.set("Filename",m_static_vertical_pressure_file);
    in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
    std::map<std::string,view_1d_host<Real>> host_views;
    std::map<std::string,FieldLayout>  layouts;
    auto pmid_ext = create_helper_field("p_mid_ext", scalar2d_layout_mid_ext, grid_ext_name, ps);
    auto pmid_ext_v = pmid_ext.get_view<Real*,Host>();
    in_params.set<std::vector<std::string>>("Field Names",{"p_levs"});
    host_views["p_levs"] = pmid_ext_v;
    layouts.emplace("p_levs",scalar2d_layout_mid_ext);
    AtmosphereInput src_input(in_params,grid_ext,host_views,layouts);
    src_input.read_variables(-1);
    src_input.finalize();
    pmid_ext.sync_to_dev();
  }

  // Open the registration!
  m_refine_remapper->registration_begins();

  // To create helper fields for later; we do both tmp and ext...
  for (auto name : m_fields_gases) {
    std::string name_ext = name + "_ext";
    std::string name_tmp = name + "_tmp";
    // Helper fields that will temporarily store the target state
    auto grid_int_name = grid_int->name();
    auto grid_ext_name = grid_ext->name();
    auto grid_tmp_name = grid_tmp->name();
    auto field  = get_field_out(name);
    auto layout = field.get_header().get_identifier().get_layout();
    auto field_ext = create_helper_field(name_ext, scalar3d_layout_mid_ext, grid_ext_name, ps);
    auto field_tmp = create_helper_field(name_tmp, scalar3d_layout_mid_tmp, grid_tmp_name, ps);
    Field field_int;
    if (m_refine_remap) {
      field_int = create_helper_field(name, scalar3d_layout_mid, grid_int_name, ps);
    } else {
      field_int             = field_tmp.alias(name);
      m_helper_fields[name] = field_int;
    }

    // Register the fields with the remapper
    m_refine_remapper->register_field(field_tmp, field_int);
    // Add the fields to the time interpolator
    m_time_interp.add_field(field_ext.alias(name), true);
  }
  m_time_interp.initialize_data_from_files();

  // Close the registration!
  m_refine_remapper->registration_ends();

}

// =========================================================================================
void PrescribedGases::run_impl (const double dt)
{
  using namespace scream::vinterp;

  // Have to add dt because first time iteration is at 0 seconds where you will
  // not have any data from the field. The timestamp is only iterated at the
  // end of the full step in scream.
  auto ts = timestamp()+dt;

  // Perform time interpolation
  m_time_interp.perform_time_interpolation(ts);

  // Process data
  const auto& p_mid_v = get_field_in("p_mid").get_view<const mPack**>();
  view_Nd<mPack,2> p_mid_ext_p;
  view_Nd<mPack,1> p_mid_ext_1d;
  if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
    p_mid_ext_p = get_helper_field("p_mid_ext").get_view<mPack**>();
  } else if (m_src_pres_type == STATIC_1D_VERTICAL_PROFILE) {
    p_mid_ext_1d = get_helper_field("p_mid_ext").get_view<mPack*>();
  }

  for (auto name : m_fields_gases) {
    auto atm_state_field = get_field_out(name);           // int horiz, int vert
    auto int_state_field = get_helper_field(name);        // int horiz, int vert
    auto ext_state_field = get_helper_field(name+"_ext"); // ext horiz, ext vert
    auto tmp_state_field = get_helper_field(name+"_tmp"); // ext horiz, int vert
    auto ext_state_view  = ext_state_field.get_view<mPack**>();
    auto tmp_state_view  = tmp_state_field.get_view<mPack**>();
    auto atm_state_view  = atm_state_field.get_view<mPack**>();  // TODO: Right now assume whatever field is defined on COLxLEV
    auto int_state_view  = int_state_field.get_view<mPack**>();
    auto int_mask_view = m_buffer.int_mask_view;
    // Masked values in the source data can lead to strange behavior in the vertical interpolation.
    // We pre-process the data and map any masked values (sometimes called "filled" values) to the
    // nearest un-masked value.
    // Here we are updating the ext_state_view, which is the time interpolated but unmapped data.
    Real var_fill_value = constants::DefaultFillValue<Real>().value;
    // Query the helper field for the fill value, if not present use default
    if (ext_state_field.get_header().has_extra_data("mask_value")) {
      var_fill_value = ext_state_field.get_header().get_extra_data<float>("mask_value");
    }
    const int num_cols           = ext_state_view.extent(0);
    const int num_vert_packs     = ext_state_view.extent(1);
    const int num_src_levs       = m_num_src_levs;
    const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);
    Kokkos::parallel_for("correct_for_masked_values", policy,
       	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto ext_state_view_1d = ekat::subview(ext_state_view,icol);
      Real fill_value;
      int  fill_idx = -1;
      // Scan top to surf and backfill all values near TOM that are masked.
      for (int kk=0; kk<num_src_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	    const auto iidx  = kk % mPack::n;
        // Check if this index is masked
        if (ext_state_view_1d(ipack)[iidx]!=var_fill_value) {
          fill_value = ext_state_view_1d(ipack)[iidx];
          fill_idx = kk;
          for (int jj=0; jj<kk; ++jj) {
                const auto jpack = jj / mPack::n;
            const auto jidx  = jj % mPack::n;
            ext_state_view_1d(jpack)[jidx] = fill_value;
          }
          break;
        }
      }
      // Now fill the rest, the fill_idx should be non-negative.  If it isn't that means
      // we have a column that is fully masked 
      for (int kk=fill_idx+1; kk<num_src_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	    const auto iidx  = kk % mPack::n;
        // Check if this index is masked
        if (!(ext_state_view_1d(ipack)[iidx]==var_fill_value)) {
          fill_value = ext_state_view_1d(ipack)[iidx];
        } else {
          ext_state_view_1d(ipack)[iidx] = fill_value;
        }
      }
    });

    // Vertical Interpolation onto atmosphere state pressure levels
    // Note that we are going from ext to tmp here
    if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
      perform_vertical_interpolation<Real,1,2>(p_mid_ext_p,
                                               p_mid_v,
                                               ext_state_view,
                                               tmp_state_view,
                                               int_mask_view,
                                               m_num_src_levs,
                                               m_num_levs);
    } else if (m_src_pres_type == STATIC_1D_VERTICAL_PROFILE) {
      perform_vertical_interpolation<Real,1,2>(p_mid_ext_1d,
                                               p_mid_v,
                                               ext_state_view,
                                               tmp_state_view,
                                               int_mask_view,
                                               m_num_src_levs,
                                               m_num_levs);
    }

    // Check that none of the targets are masked, if they are, set value to
    // nearest unmasked value above.
    // NOTE: We use an algorithm whichs scans from TOM to the surface.
    //       If TOM is masked we keep scanning until we hit an unmasked value,
    //       we then set all masked values above to the unmasked value.
    //       We continue scanning towards the surface until we hit an unmasked value, we
    //       then assign that masked value the most recent unmasked value, until we hit the
    //       surface.
    // Here we change the int_state_view which represents the vertically interpolated fields onto
    // the simulation grid.
    const int num_levs = m_num_levs;
    Kokkos::parallel_for("correct_for_masked_values", policy,
       	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto int_mask_view_1d  = ekat::subview(int_mask_view,icol);
      auto tmp_state_view_1d = ekat::subview(tmp_state_view,icol);
      Real fill_value;
      int  fill_idx = -1;
      // Scan top to surf and backfill all values near TOM that are masked.
      for (int kk=0; kk<num_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	    const auto iidx  = kk % mPack::n;
        // Check if this index is masked
        if (!int_mask_view_1d(ipack)[iidx]) {
          fill_value = tmp_state_view_1d(ipack)[iidx];
          fill_idx = kk;
          for (int jj=0; jj<fill_idx; ++jj) {
                const auto jpack = jj / mPack::n;
            const auto jidx  = jj % mPack::n;
            tmp_state_view_1d(jpack)[jidx] = fill_value;
          }
          break;
        }
      }
      // Now fill the rest, the fill_idx should be non-negative.  If it isn't that means
      // we have a column that is fully masked
      for (int kk=fill_idx+1; kk<num_levs; ++kk) {
        const auto ipack = kk / mPack::n;
        const auto iidx  = kk % mPack::n;
        // Check if this index is masked
        if (!int_mask_view_1d(ipack)[iidx]) {
          fill_value = tmp_state_view_1d(ipack)[iidx];
        } else {
          tmp_state_view_1d(ipack)[iidx] = fill_value;
        }
      }
    });

  }

  // Refine-remap onto target atmosphere state horiz grid (int);
  // note that if the data comes from the same grid as the model,
  // this remap step is a no-op; otherwise, we refine-remap from tmp to int
  m_refine_remapper->remap(true);
  for (auto name : m_fields_gases) {
    auto atm_state_field = get_field_out(name);
    auto int_state_field = get_helper_field(name);
    auto atm_state_view  = atm_state_field.get_view<mPack**>();
    auto int_state_view  = int_state_field.get_view<mPack**>();
    Kokkos::deep_copy(atm_state_view,int_state_view);
  }
}

// =========================================================================================
void PrescribedGases::finalize_impl()
{
  m_time_interp.finalize();
}
// =========================================================================================
Field PrescribedGases::create_helper_field (const std::string& name,
                                  const FieldLayout& layout,
                                  const std::string& grid_name,
                                  const int ps)
{
  using namespace ekat::units;
  // For helper fields we don't bother w/ units, so we set them to non-dimensional
  FieldIdentifier id(name,layout,Units::nondimensional(),grid_name);

  // Create the field. Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  if (ps>=0) {
    f.get_header().get_alloc_properties().request_allocation(ps);
  } else {
    f.get_header().get_alloc_properties().request_allocation();
  }
  f.allocate_view();
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[name] = f;
  return m_helper_fields[name];
}

// =========================================================================================
size_t PrescribedGases::requested_buffer_size_in_bytes() const {
  return m_buffer.num_2d_midpoint_mask_views*m_num_cols*m_num_levs*sizeof(mMask);
}

// =========================================================================================
void PrescribedGases::init_buffers(const ATMBufferManager& buffer_manager) {
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error, PrescribedGases::init_buffers! Buffers size not sufficient.\n");
  mMask* mem = reinterpret_cast<mMask*>(buffer_manager.get_memory());

  m_buffer.int_mask_view = decltype(m_buffer.int_mask_view)(mem,m_num_cols,m_num_levs);
  mem += m_buffer.int_mask_view.size();

  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error: PrescribedGases::init_buffers! Used memory != requested memory.");
}
// =========================================================================================

} // namespace scream
