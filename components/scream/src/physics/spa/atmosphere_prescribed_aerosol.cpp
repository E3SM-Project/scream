#include "atmosphere_prescribed_aerosol.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
  using namespace spa;
// =========================================================================================
SPA::SPA (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

// =========================================================================================
void SPA::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  auto nondim = Units::nondimensional();

  const auto& grid_name = m_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column
  m_dofs_gids = grid->get_dofs_gids();
  m_total_global_dofs = grid->get_num_global_dofs();
  m_min_global_dof    = grid->get_global_min_dof_gid();

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  FieldLayout scalar2d_layout     { {COL},     {m_num_cols} };
  FieldLayout scalar1d_layout_mid { {LEV},     {m_num_levs} };
  // Use VAR field tag for gases for now; consider adding a tag?
  FieldLayout scalar3d_swband_layout { {COL,SWBND, LEV}, {m_num_cols, m_nswbands, m_num_levs} }; 
  FieldLayout scalar3d_lwband_layout { {COL,LWBND, LEV}, {m_num_cols, m_nlwbands, m_num_levs} }; 

  // Set of fields used strictly as input
  constexpr int ps = Pack::n;
  add_field<Required>("p_mid"      , scalar3d_layout_mid, Pa,     grid_name, ps);
  add_field<Required>("hyam"       , scalar1d_layout_mid, nondim, grid_name, ps); // TODO: These fields should  be loaded from file and not registered with the field manager.
  add_field<Required>("hybm"       , scalar1d_layout_mid, nondim, grid_name, ps); // TODO: These fields should  be loaded from file and not registered with the field manager.

  // Set of fields used strictly as output
  add_field<Computed>("nc_activated",   scalar3d_layout_mid,    1/kg,   grid_name,ps);
  add_field<Computed>("aero_g_sw",      scalar3d_swband_layout, nondim, grid_name,ps);
  add_field<Computed>("aero_ssa_sw",    scalar3d_swband_layout, nondim, grid_name,ps);
  add_field<Computed>("aero_tau_sw",    scalar3d_swband_layout, nondim, grid_name,ps);
  add_field<Computed>("aero_tau_lw",    scalar3d_lwband_layout, nondim, grid_name,ps);

  // Init output data structure
  SPAData_out.init(m_num_cols,m_num_levs,m_nswbands,m_nlwbands,false);
}
// =========================================================================================
size_t SPA::requested_buffer_size_in_bytes() const
{
  using PackInfo = ekat::PackInfo<Spack::n>;

  // Recall: the quantities in spa_temp defined over vlevs have 1 Real of
  //         padding in each column (at beginning and end).
  //         That's why we have m_num_levs+2
  const int nlevs = m_num_levs+2;
  const int num_mid_packs = PackInfo::num_packs(nlevs);
  const int nlevs_alloc = num_mid_packs*Spack::n;

  // We have
  //  - one (ncols) view (spa_temp's ps)
  //  - two (ncols,nlevs) mid view (p_mid_src, spa_temp's ccn)
  //  - three (ncols,nswbands,nlevs) views (spa_temp's aer_g_sw, aer_ssa_sw, aer_tau_sw)
  //  - one (ncols,nlwbands,nlevs) view (aer_tau_lw)
  const int num_reals = m_num_cols*(1+nlevs_alloc*(2 + 3*m_nswbands + m_nlwbands));

  return num_reals*sizeof(Real);
}

// =========================================================================================
void SPA::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Buffers size not sufficient.\n");

  using PackInfo = ekat::PackInfo<Spack::n>;

  // Short names make following rows fit on text editor screen
  // Recall: the quantities in spa_temp defined over vlevs have 1 Real of
  //         padding in each column (at beginning and end).
  //         That's why we have m_num_levs+2
  const int nlevs  = m_num_levs+2;
  const int npacks = PackInfo::num_packs(nlevs);
  const int ncols  = m_num_cols;
  const int nswb   = m_nswbands;
  const int nlwb   = m_nlwbands;

  Spack* mem = reinterpret_cast<Spack*>(buffer_manager.get_memory());

  // Source pressure levels
  m_buffer.p_mid_src = decltype(m_buffer.p_mid_src)(mem, ncols, npacks);
  mem += m_buffer.p_mid_src.size();

  // SPA temporaries
  auto& spa_data = m_buffer.spa_temp.data;
  spa_data.init(ncols,nlevs,nswb,nlwb,false);

  spa_data.CCN3 = decltype(spa_data.CCN3)(mem, ncols, npacks);
  mem += spa_data.CCN3.size();

  spa_data.AER_G_SW = decltype(spa_data.AER_G_SW)(mem, ncols, nswb, npacks);
  mem += spa_data.AER_G_SW.size();
  spa_data.AER_SSA_SW = decltype(spa_data.AER_SSA_SW)(mem, ncols, nswb, npacks);
  mem += spa_data.AER_SSA_SW.size();
  spa_data.AER_TAU_SW = decltype(spa_data.AER_TAU_SW)(mem, ncols, nswb, npacks);
  mem += spa_data.AER_TAU_SW.size();

  spa_data.AER_TAU_LW = decltype(spa_data.AER_TAU_LW)(mem, ncols, nlwb, npacks);
  mem += spa_data.AER_TAU_LW.size();

  Real* r_mem = reinterpret_cast<Real*>(mem);
  m_buffer.spa_temp.PS = decltype(m_buffer.spa_temp.PS)(r_mem,ncols);
  r_mem += m_buffer.spa_temp.PS.size();

  size_t used_mem = (r_mem - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(),
      "Error! Used memory != requested memory for SPA.\n"
      "   - used mem     : " + std::to_string(used_mem) + "\n"
      "   - requested mem: " + std::to_string(requested_buffer_size_in_bytes()) + "\n");
}

// =========================================================================================
void SPA::initialize_impl (const RunType /* run_type */)
{
  // Initialize SPA pressure state stucture and set pointers for the SPA output data to
  // field managed variables.
  SPAData_out.CCN3               = get_field_out("nc_activated").get_view<Pack**>();
  SPAData_out.AER_G_SW           = get_field_out("aero_g_sw").get_view<Pack***>();
  SPAData_out.AER_SSA_SW         = get_field_out("aero_ssa_sw").get_view<Pack***>();
  SPAData_out.AER_TAU_SW         = get_field_out("aero_tau_sw").get_view<Pack***>();
  SPAData_out.AER_TAU_LW         = get_field_out("aero_tau_lw").get_view<Pack***>();

  // Retrieve the remap and data file locations from the parameter list:
  EKAT_REQUIRE_MSG(m_params.isParameter("SPA Remap File"),"ERROR: SPA Remap File is missing from SPA parameter list.");
  EKAT_REQUIRE_MSG(m_params.isParameter("SPA Data File"),"ERROR: SPA Data File is missing from SPA parameter list.");
  m_spa_remap_file = m_params.get<std::string>("SPA Remap File");
  m_spa_data_file = m_params.get<std::string>("SPA Data File");

  // Set the SPA remap weights.  
  // TODO: We may want to provide an option to calculate weights on-the-fly. 
  //       If so, then the EKAT_REQUIRE_MSG above will need to be removed and 
  //       we can have a default m_spa_data_file option that is online calculation.
  using ci_string = ekat::CaseInsensitiveString;
  ci_string no_filename = "none";
  if (m_spa_remap_file == no_filename) {
    printf("WARNING: SPA Remap File has been set to 'NONE', assuming that SPA data and simulation are on the same grid - skipping horizontal interpolation");
    SPAFunc::set_remap_weights_one_to_one(m_total_global_dofs,m_min_global_dof,m_dofs_gids,SPAHorizInterp);
  } else {
    SPAFunc::get_remap_weights_from_file(m_spa_remap_file,m_total_global_dofs,m_min_global_dof,m_dofs_gids,SPAHorizInterp);
  }
  // Note: only the number of levels associated with this data haven't been set.  We can
  //       take this information directly from the spa data file.
  scorpio::register_file(m_spa_data_file,scorpio::Read);
  const int source_data_nlevs = scorpio::get_dimlen_c2f(m_spa_data_file.c_str(),"lev")+2; // Add 2 for padding
  SPAHorizInterp.m_comm = m_comm;

  // Initialize the size of the SPAData structures:
  SPAData_start = SPAFunc::SPAInput(m_dofs_gids.size(), source_data_nlevs, m_nswbands, m_nlwbands);
  SPAData_end   = SPAFunc::SPAInput(m_dofs_gids.size(), source_data_nlevs, m_nswbands, m_nlwbands);

  // Update the local time state information and load the first set of SPA data for interpolation:
  auto ts = timestamp();
  SPATimeState.inited = false;
  SPATimeState.current_month = ts.get_month();
  SPAFunc::update_spa_timestate(m_spa_data_file,m_nswbands,m_nlwbands,ts,SPAHorizInterp,SPATimeState,SPAData_start,SPAData_end);

  // NOTE: we *assume* hybrid v coordinates don't change with time.
  //       IF this ever ceases to be the case, you need to remove these
  //       lines, and have spa_main interpolate those during the call
  //       to performe_time_interpolation.
  m_buffer.spa_temp.hyam = SPAData_start.hyam;
  m_buffer.spa_temp.hybm = SPAData_start.hybm;

  // Set property checks for fields in this process
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("nc_activated"),0.0,1.0e6,false);
  // upper bound set to 1.01 as max(g_sw)=1.00757 in current ne4 data assumingly due to remapping
  // add an epslon to max possible upper bound of aero_ssa_sw

  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("aero_g_sw"),0.0,1.0,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("aero_ssa_sw"),0.0,1.0,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("aero_tau_sw"),0.0,1.0,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("aero_tau_lw"),0.0,1.0,true);
}

// =========================================================================================
void SPA::run_impl (const int  dt )
{
  /* Gather time and state information for interpolation */
  auto ts = timestamp();
  // We always want to update the current time in the time_state.
  SPATimeState.t_now = ts.frac_of_year_in_days() + dt/86400.;
  
  
  /* Update time state and if the month has changed, update the data.*/
  SPAFunc::update_spa_timestate(m_spa_data_file,m_nswbands,m_nlwbands,ts,SPAHorizInterp,SPATimeState,SPAData_start,SPAData_end);

  // Call the main SPA routine to get interpolated aerosol forcings.
  const auto& pmid_tgt = get_field_in("p_mid").get_view<const Pack**>();
  SPAFunc::spa_main(SPATimeState, pmid_tgt, m_buffer.p_mid_src,
                    SPAData_start,SPAData_end,m_buffer.spa_temp,SPAData_out);
}

// =========================================================================================
void SPA::finalize_impl()
{
  // Do nothing
}

} // namespace scream
