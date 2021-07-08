#include "atmosphere_prescribed_aero.hpp"

#include "share/io/scorpio_input.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
//  using namespace cld_fraction;
// =========================================================================================
SPA::SPA (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_spa_comm (comm)
 , m_spa_params (params)
{
  // Nothing to do here
}

// =========================================================================================
void SPA::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto& grid_name = m_spa_params.get<std::string>("Grid");
//  auto grid  = grids_manager->get_grid(grid_name);
  m_grid     = grids_manager->get_grid(grid_name);
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  constexpr int ps = SCREAM_SMALL_PACK_SIZE;
  // Set of fields used as input and output
  add_field<Updated>("nc_activated",       scalar3d_layout_mid, 1/kg,     grid_name, ps);
  
  // Collect the spa data filename from the parameter list
  m_input_data_filename = m_spa_params.get<std::string>("Input Data Filename");
  printf("ASD Init - %s\n",m_input_data_filename.c_str());

  // DUMMY values for now, TODO: get real data from file.
  MonthlyCCN.CCN = SPAF::view_3d<Spack>("CCN",m_num_cols,12,m_num_levs);
}

// =========================================================================================
void SPA::initialize_impl (const util::TimeStamp& /* t0 */)
{
  // Initialize tracking of which month the simulation is in.
  auto ts = timestamp();
  m_current_month = ts.get_months() + 1;
  m_next_month    = m_current_month + 1;

  // Set the SPA structures with pointers to FM data
  PrescribedAero.CCN = m_spa_fields_out["nc_activated"].get_reshaped_view<Pack**>();

  // Get data from input file
  using input_type = AtmosphereInput;
  input_type spa_input(m_spa_comm,"Physics",m_grid);
  std::vector<std::string> var_dims = {"ncol","lev"};
  std::vector<int> dim_lens = {m_num_cols,m_num_levs};
//  spa_input.pull_input(m_input_data_filename,"CCN3",var_dims, true, dim_lens, nc_activated.data());
}

// =========================================================================================
void SPA::run_impl (const Real dt)
{

  // Get a copy of the current timestamp (at the beginning of the step)
  // Needed for interpolation of SPA data.
  auto ts = timestamp();
  {
  std::string time_str = ts.to_string();
  printf("ASD - %s\n",time_str.c_str());
  }

  // Gather information from the timestamp which is needed for temporal interpolation of
  // SPA data.
  //   mm_0 is the month we are currently in.
  //   mm_1 is next month.
  // combined these two will be used with the SPA data to calculate the slope and intercept
  // for the temporal interpolation.
  //   sec_from_mm_0 is the distance into the month (in sec) that we are currently at.
  // This will be used with the slope/intercept form to calculate the SPA value at this point
  // in time.
  int mm_0 = ts.get_months() + 1;
  if (mm_0 != m_current_month) {
    m_current_month = mm_0;
    m_next_month    = (mm_0 % 12) + 1;
  }

  Real days_from_mm_0 = ts.get_days() + ts.get_seconds()/86400.0;
  printf("   - mm_0 = %d, mm_1 = %d, sec = %f\n",m_current_month,m_next_month,days_from_mm_0);

  SPAF::main(m_num_cols,m_num_levs,ts.get_months(),days_from_mm_0,MonthlyGHG,MonthlyCCN,PrescribedAero);

  // Advance the timestamp
  ts += dt;
  for (auto& f : m_spa_fields_out) {
    printf("ASD update timestamp for %s\n",f.first.c_str());
    f.second.get_header().get_tracking().update_time_stamp(ts);
  }
}

// =========================================================================================
void SPA::finalize_impl()
{
  // Do nothing
}

void SPA::set_required_field_impl (const Field<const Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_spa_fields_in.emplace(name,f);

  // Add myself as customer to the field
  add_me_as_customer(f);
}

void SPA::set_computed_field_impl (const Field<      Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_spa_fields_out.emplace(name,f);

  // Add myself as provider for the field
  add_me_as_provider(f);
}

} // namespace scream
