#include "atmosphere_prescribed_aero.hpp"

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
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  constexpr int ps = SCREAM_SMALL_PACK_SIZE;
  // Set of fields used as input and output
  add_field<Updated>("nc_activated",       scalar3d_layout_mid, 1/kg,     grid_name, ps);
  
  // Collect the spa data filename from the parameter list
  m_input_data_filename = m_spa_params.get<std::string>("Input Data Filename");
  printf("ASD Init - %s\n",m_input_data_filename.c_str());
}

// =========================================================================================
void SPA::initialize_impl (const util::TimeStamp& /* t0 */)
{
  // Do nothing for now
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
  

//  SPAFunc::main(m_num_cols,m_num_levs,qi,liq_cld_frac,ice_cld_frac,tot_cld_frac);

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
