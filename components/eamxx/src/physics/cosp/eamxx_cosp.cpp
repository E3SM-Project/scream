#include "eamxx_cosp.hpp"
#include "cosp_functions.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
// =========================================================================================
Cosp::Cosp (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

// =========================================================================================
void Cosp::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  auto nondim = Units::nondimensional();

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar2d_layout     { {COL},      {m_num_cols}              };
  FieldLayout scalar3d_layout_mid { {COL,LEV},  {m_num_cols,m_num_levs}   };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {m_num_cols,m_num_levs+1} };
  //FieldLayout scalar4d_ctptau     { {COL,TAU,CTP}, {m_num_cols,m_num_tau,m_num_ctp} };

  constexpr int ps = Pack::n;

  // Set of fields used strictly as input
  //                  Name in AD     Layout               Units   Grid       Group
  add_field<Required>("surf_radiative_T", scalar2d_layout    , K,      grid_name, ps);
  //add_field<Required>("surfelev",    scalar2d_layout    , m,      grid_name, ps);
  //add_field<Required>("landmask",    scalar2d_layout    , nondim, grid_name, ps);
  //add_field<Required>("horiz_wind",  scalar3d_layout_mid, m/s,    grid_name, ps);
  add_field<Required>("sunlit",           scalar2d_layout    , nondim, grid_name, ps);
  add_field<Required>("p_mid",             scalar3d_layout_mid, Pa,     grid_name, ps);
  add_field<Required>("p_int",             scalar3d_layout_int, Pa,     grid_name, ps);
  //add_field<Required>("height_mid",  scalar3d_layout_mid, m,      grid_name, ps);
  //add_field<Required>("height_int",  scalar3d_layout_int, m,      grid_name, ps);
  add_field<Required>("T_mid",            scalar3d_layout_mid, K,      grid_name, ps);
  add_field<Required>("qv",               scalar3d_layout_mid, Q,      grid_name, "tracers", ps);
  add_field<Required>("qc",               scalar3d_layout_mid, Q,      grid_name, "tracers", ps);
  add_field<Required>("qi",               scalar3d_layout_mid, Q,      grid_name, "tracers", ps);
  add_field<Required>("qr",               scalar3d_layout_mid, Q,      grid_name, "tracers", ps);
  add_field<Required>("qm",               scalar3d_layout_mid, Q,      grid_name, "tracers", ps);
  add_field<Required>("cldfrac_tot_for_analysis", scalar3d_layout_mid, nondim, grid_name, ps);
  // Optical properties, should be computed in radiation interface
  add_field<Required>("dtau067",     scalar3d_layout_mid, nondim, grid_name, ps); // 0.67 micron optical depth
  add_field<Required>("dtau105",     scalar3d_layout_mid, nondim, grid_name, ps); // 10.5 micron optical depth
  // Effective radii, should be computed in either microphysics or radiation interface
  add_field<Required>("reff_qc",     scalar3d_layout_mid, m,      grid_name, ps);
  add_field<Required>("reff_qi",     scalar3d_layout_mid, m,      grid_name, ps);
  //add_field<Required>("reff_qr",     scalar3d_layout_mid, m,      grid_name, ps);
  //add_field<Required>("reff_qm",     scalar3d_layout_mid, m,      grid_name, ps);

  // Set of fields used strictly as output
  //add_field<Computed>("isccp_ctptau_hist", scalar4d_layout_ctptau, nondim, grid_name, ps);
  add_field<Computed>("isccp_cldtot", scalar2d_layout, nondim, grid_name, ps);
}

// =========================================================================================
void Cosp::initialize_impl (const RunType /* run_type */)
{
  // Set property checks for fields in this process
  using Interval = FieldWithinIntervalCheck;
  //add_postcondition_check<Interval>(get_field_out("cldfrac_tot"),m_grid,0.0,1.0,false);
  //add_postcondition_check<Interval>(get_field_out("cldfrac_tot_for_analysis"),m_grid,0.0,1.0,false);
  std::cout << "BRHDEBUG: call CospFunc::initialize()" << std::endl;
  CospFunc::initialize(m_num_cols, m_num_subcols, m_num_levs);
}

// =========================================================================================
void Cosp::run_impl (const double /* dt */)
{
  // Get fields from field manager
  // These should maybe get switched to Packs for convenience with other routines that expect packs
  auto qv      = get_field_in("qv").get_view<const Pack**>();
  auto qc      = get_field_in("qc").get_view<const Pack**>();
  auto qi      = get_field_in("qi").get_view<const Pack**>();
  auto qr      = get_field_in("qr").get_view<const Pack**>();
  auto qm      = get_field_in("qm").get_view<const Pack**>();
  auto sunlit  = get_field_in("sunlit").get_view<const Pack*>();  // Grab incoming shortwave and construct sunlit
  auto skt     = get_field_in("surf_radiative_T").get_view<const Pack*>();
  auto T_mid   = get_field_in("T_mid").get_view<const Pack**>();
  auto p_mid   = get_field_in("p_mid").get_view<const Pack**>();
  auto p_int   = get_field_in("p_int").get_view<const Pack**>();
  auto cldfrac = get_field_in("cldfrac_tot_for_analysis").get_view<const Pack**>();
  auto reff_qc = get_field_in("reff_qc").get_view<const Pack**>();
  auto reff_qi = get_field_in("reff_qi").get_view<const Pack**>();
  auto dtau067 = get_field_in("dtau067").get_view<const Pack**>();
  auto dtau105 = get_field_in("dtau105").get_view<const Pack**>();

  auto cldfrac_tot_for_analysis = get_field_in("cldfrac_tot_for_analysis").get_view<const Pack**>();
  //auto isccp_ctptau_hist = get_field_out("isccp_ctptau_hist").get_view<Pack**>();
  auto isccp_cldtot = get_field_out("isccp_cldtot").get_view<Pack*>();

  // Call COSP wrapper routines
  Real emsfc_lw = 0.99;
  CospFunc::main(
          m_num_cols, m_num_subcols, m_num_levs,
          emsfc_lw, sunlit, skt, T_mid, p_mid, p_int, qv,
          cldfrac, reff_qc, reff_qi, dtau067, dtau105,
          isccp_cldtot
  );
}

// =========================================================================================
void Cosp::finalize_impl()
{
  // Finalize COSP wrappers
    CospFunc::finalize();
}
// =========================================================================================

} // namespace scream