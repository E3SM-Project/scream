#include "eamxx_nudging_process_interface.hpp"

namespace scream
{

  //using namespace spa;
// =========================================================================================
Nudging::Nudging (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  m_datafiles  = m_params.get<std::vector<std::string>>("nudging_filename");
  m_timescale = m_params.get<int>("nudging_timescale",0);
  m_fields_nudge = m_params.get<std::vector<std::string>>("nudging_fields");
  auto src_pres_type = m_params.get<std::string>("source_pressure_type","DYNAMIC");
  if (src_pres_type=="DYNAMIC") {
    m_src_pres_type = DYNAMIC;
  } else if (src_pres_type=="STATIC") {
    m_src_pres_type = STATIC;
  } else {
    EKAT_ERROR_MSG("ERROR! Nudging::parameter_list - unsupported source_pressure_type provided.  Current options are [DYNAMICS,STATIC].  Please check");
  }
  // TODO: Add some warning messages here.
  // 1. if m_timescale is <= 0 we will do direct replacement.
  // 2. if m_fields_nudge is empty or =NONE then we will skip nudging altogether.
}

// =========================================================================================
void Nudging::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
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
  Q.set_string("kg/kg");
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);

  /* ----------------------- WARNING --------------------------------*/
  /* The following is a HACK to get things moving, we don't want to
   * add all fields as "updated" long-term.  A separate stream of work
   * is adapting the infrastructure to allow for a generic "add_field" call
   * to be used here which we can then setup using the m_fields_nudge variable
   */
  add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
  add_field<Updated>("qv",    scalar3d_layout_mid, Q, grid_name, "tracers", ps);
  add_field<Updated>("u",     scalar3d_layout_mid, m/s, grid_name, ps);
  add_field<Updated>("v",     scalar3d_layout_mid, m/s, grid_name, ps);
  /* ----------------------- WARNING --------------------------------*/

  //Now need to read in the file
  scorpio::register_file(m_datafiles[0],scorpio::Read);
  m_num_src_levs = scorpio::get_dimlen(m_datafiles[0],"lev");
  scorpio::eam_pio_closefile(m_datafiles[0]);
}
// =========================================================================================
void Nudging::apply_tendency(Field& base, const Field& next, const int dt)
{
  // Calculate the weight to apply the tendency
  const Real dtend = Real(dt)/Real(m_timescale);
  EKAT_REQUIRE_MSG(dtend>=0,"Error! Nudging::apply_tendency - timescale tendency of " << std::to_string(dt) 
		  << " / " << std::to_string(m_timescale) << " = " << std::to_string(dtend) 
		  << " is invalid.  Please check the timescale and/or dt");
  // Now apply the tendency.
  Field tend = base.clone();
  // Use update internal to set tendency, will be (1.0*next - 1.0*base), note tend=base at this point.
  tend.update(next,Real(1.0),Real(-1.0));
  base.update(tend,dtend,Real(1.0));
}
// =========================================================================================
void Nudging::initialize_impl (const RunType /* run_type */)
{
  using namespace ShortFieldTagsNames;

  // Initialize the time interpolator
  auto grid_ext = m_grid->clone("Point Grid", false);
  grid_ext->reset_num_vertical_lev(m_num_src_levs);
  FieldLayout scalar2d_layout_mid { {LEV}, {m_num_src_levs} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_src_levs} };
  m_time_interp = util::TimeInterpolation(grid_ext, m_datafiles);

  constexpr int ps = 1;  // TODO: I think this could be the regular packsize, right?
  const auto& grid_name = m_grid->name();
  if (m_src_pres_type == DYNAMIC) {
    create_helper_field("p_mid_ext", scalar3d_layout_mid, grid_name, ps);
    auto pmid_ext = get_helper_field("p_mid_ext");
    m_time_interp.add_field(pmid_ext,"p_mid",true);
  } else if (m_src_pres_type == STATIC) {
    // Load p_levs from source data file
    create_helper_field("p_mid_ext", scalar2d_layout_mid, grid_name, ps);
    auto pmid_ext = get_helper_field("p_mid_ext");
    auto pmid_ext_v = pmid_ext.get_view<Real*,Host>();
    // Extract p_levs information from file.  The default is to assume that data is in the nudging source data, but we give the
    // user the runtime option to specify a different file.
    auto plevs_file = m_params.get<std::string>("static_pressure_levels_file",m_datafiles[0]);
    ekat::ParameterList in_params;
    in_params.set<std::vector<std::string>>("Field Names",{"p_levs"});
    in_params.set("Filename",plevs_file);
    in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
    std::map<std::string,view_1d_host<Real>> host_views;
    std::map<std::string,FieldLayout>  layouts;
    host_views["p_levs"] = pmid_ext_v;
    layouts.emplace("p_levs",scalar2d_layout_mid);
    AtmosphereInput src_input(in_params,grid_ext,host_views,layouts);
    src_input.read_variables(-1);
    src_input.finalize();
  }
  for (auto name : m_fields_nudge) {
    std::string name_ext = name + "_ext";
    // Helper fields that will temporarily store the target state, which can then
    // be used to back out a nudging tendency
    auto field  = get_field_out(name);
    auto layout = field.get_header().get_identifier().get_layout();
    create_helper_field(name,     layout,              grid_name, ps);
    create_helper_field(name_ext, scalar3d_layout_mid, grid_name, ps);
    auto field_ext = get_helper_field(name_ext);
    m_time_interp.add_field(field_ext,name,true);
  }
  m_time_interp.initialize_data_from_files();

}

// =========================================================================================
void Nudging::run_impl (const double dt)
{
  using namespace scream::vinterp;

  // Have to add dt because first time iteration is at 0 seconds where you will
  // not have any data from the field. The timestamp is only iterated at the
  // end of the full step in scream.
  auto ts = timestamp()+dt;

  // Perform time interpolation
  m_time_interp.perform_time_interpolation(ts);
  auto current_file = m_time_interp.get_current_file();

  // Process data and nudge the atmosphere state
  const auto& p_mid_v     = get_field_in("p_mid").get_view<const mPack**>();
  view_Nd<mPack,2> p_mid_ext_p;
  view_Nd<mPack,1> p_mid_ext_1d;
  if (m_src_pres_type == DYNAMIC) {
    const auto& p_mid_ext   = get_helper_field("p_mid_ext").get_view<mPack**>();
    p_mid_ext_p = view_Nd<mPack,2>(reinterpret_cast<mPack*>(p_mid_ext.data()),
  				     m_num_cols,m_num_src_levs);
  } else if (m_src_pres_type == STATIC) {
    p_mid_ext_1d   = get_helper_field("p_mid_ext").get_view<mPack*>();
  }
  for (auto name : m_fields_nudge) {
    auto atm_state_field = get_field_out(name);
    auto int_state_field = get_helper_field(name);
    auto ext_state_field = get_helper_field(name+"_ext").get_view<mPack**>();
    auto atm_state_view  = atm_state_field.get_view<mPack**>();  // TODO: Right now assume whatever field is defined on COLxLEV
    auto int_state_view  = int_state_field.get_view<mPack**>();
    view_Nd<mMask,2> int_mask_view("",m_num_cols,m_num_levs); // Track mask of interpolated values
    const view_Nd<mPack,2> ext_state_view(reinterpret_cast<mPack*>(ext_state_field.data()),
                                          m_num_cols,m_num_src_levs);
    // Masked values in the source data can lead to strange behavior in the vertical interpolation.
    // We pre-process the data and map any masked values (sometimes called "filled" values) to the
    // nearest un-masked value.
    Real var_fill_value;
    scorpio::get_variable_metadata(current_file,name,"_FillValue",var_fill_value);
    const int num_cols           = int_state_view.extent(0);
    const int num_vert_packs     = int_state_view.extent(1);
    const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);
    Kokkos::parallel_for("correct_for_masked_values", policy,
       	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto int_mask_view_1d  = ekat::subview(int_mask_view,icol);
      auto int_state_view_1d = ekat::subview(int_state_view,icol);
      Real fill_value;
      int  fill_idx = -1;
      // Scan top to surf and backfill all values near TOM that are masked.
      for (int kk=0; kk<m_num_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	const auto iidx  = kk % mPack::n;
        // Check if this index is masked
	if (!(int_state_view_1d(ipack)[iidx]==var_fill_value)) {
	  fill_value = int_state_view_1d(ipack)[iidx];
	  fill_idx = kk;
	  for (int jj=0; jj<kk; ++jj) {
            const auto jpack = jj / mPack::n;
	    const auto jidx  = jj % mPack::n;
	    int_state_view_1d(jpack)[jidx] = fill_value;
	  }
	  break;
	}
      }
      // Now fill the rest, the fill_idx should be non-negative.  If it isn't that means
      // we have a column that is fully masked 
      for (int kk=fill_idx+1; kk<m_num_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	const auto iidx  = kk % mPack::n;
        // Check if this index is masked
	if (!(int_state_view_1d(ipack)[iidx]==var_fill_value)) {
	  fill_value = int_state_view_1d(ipack)[iidx];
	} else {
	  int_state_view_1d(ipack)[iidx] = fill_value;
	}
      }
    });


    // Vertical Interpolation onto atmosphere state pressure levels
    if (m_src_pres_type == DYNAMIC) {
      perform_vertical_interpolation<Real,1,2>(p_mid_ext_p,
                                               p_mid_v,
                                               ext_state_view,
                                               int_state_view,
                                               int_mask_view,
                                               m_num_src_levs,
                                               m_num_levs);
    } else if (m_src_pres_type == STATIC) {
      perform_vertical_interpolation<Real,1,2>(p_mid_ext_1d,
                                               p_mid_v,
                                               ext_state_view,
                                               int_state_view,
                                               int_mask_view,
                                               m_num_src_levs,
                                               m_num_levs);
    }

    // Check that none of the nudging targets are masked, if they are, set value to
    // nearest unmasked value above.
    // NOTE: We use an algorithm whichs scans from TOM to the surface.
    //       If TOM is masked we keep scanning until we hit an unmasked value,
    //       we then set all masked values above to the unmasked value.
    //       We continue scanning towards the surface until we hit an unmasked value, we
    //       then assign that masked value the most recent unmasked value, until we hit the
    //       surface.
    Kokkos::parallel_for("correct_for_masked_values", policy,
       	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto int_mask_view_1d  = ekat::subview(int_mask_view,icol);
      auto int_state_view_1d = ekat::subview(int_state_view,icol);
      Real fill_value;
      int  fill_idx = -1;
      // Scan top to surf and backfill all values near TOM that are masked.
      for (int kk=0; kk<m_num_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	const auto iidx  = kk % mPack::n;
        // Check if this index is masked
	if (!int_mask_view_1d(ipack)[iidx]) {
	  fill_value = int_state_view_1d(ipack)[iidx];
	  fill_idx = kk;
	  for (int jj=0; jj<kk; ++jj) {
            const auto jpack = jj / mPack::n;
	    const auto jidx  = jj % mPack::n;
	    int_state_view_1d(jpack)[jidx] = fill_value;
	  }
	  break;
	}
      }
      // Now fill the rest, the fill_idx should be non-negative.  If it isn't that means
      // we have a column that is fully masked 
      for (int kk=fill_idx+1; kk<m_num_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	const auto iidx  = kk % mPack::n;
        // Check if this index is masked
	if (!int_mask_view_1d(ipack)[iidx]) {
	  fill_value = int_state_view_1d(ipack)[iidx];
	} else {
	  int_state_view_1d(ipack)[iidx] = fill_value;
	}
      }
    });

    // Apply the nudging tendencies to the ATM state
    if (m_timescale <= 0) {
      // We do direct replacement
      Kokkos::deep_copy(atm_state_view,int_state_view);
    } else {
      // Back out a tendency and apply it.
      apply_tendency(atm_state_field, int_state_field, dt);
    }

  }
}

// =========================================================================================
void Nudging::finalize_impl()
{
  m_time_interp.finalize();
}

} // namespace scream
