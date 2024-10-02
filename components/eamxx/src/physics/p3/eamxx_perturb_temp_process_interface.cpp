#include "physics/p3/eamxx_perturb_temp_process_interface.hpp"

namespace scream {

// =========================================================================================
PerturbTemp::PerturbTemp(const ekat::Comm &comm,
                         const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}

// =========================================================================================
void PerturbTemp::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;
  
  auto q_unit = kg / kg;
  q_unit.set_string("kg/kg");

  auto n_unit = 1 / kg;  // units of number mixing ratios of tracers
  n_unit.set_string("#/kg");


  grid_ = grids_manager->get_grid("Physics");

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column

  const auto &grid_name = grid_->name();

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  const FieldLayout scalar3d_layout_mid{{COL, LEV}, {ncol_, nlev_}};
  const FieldLayout scalar3d_layout_int{{COL, ILEV}, {ncol_, nlev_ + 1}};

  add_field<Required>("T_mid", scalar3d_layout_mid, K,
                      grid_name);  // temperature [K]

  add_field<Required>("p_mid", scalar3d_layout_mid, Pa,
                      grid_name);  // pressure at mid points in [Pa]
  add_field<Required>("p_int", scalar3d_layout_int, Pa,
                      grid_name);  // total pressure
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,
                      grid_name);  // pseudo density in [Pa]
  add_field<Required>("qv", scalar3d_layout_mid, q_unit, grid_name,
                      "tracers");  // specific humidity
  add_field<Required>("qc", scalar3d_layout_mid, q_unit, grid_name,
                      "tracers");  // liquid cloud water [kg/kg] wet
  add_field<Required>("qi", scalar3d_layout_mid, q_unit, grid_name,
                      "tracers");  // ice cloud water [kg/kg] wet
  add_field<Updated>("nc", scalar3d_layout_mid, n_unit, grid_name,
                     "tracers");  // cloud liquid wet number mixing ratio
  add_field<Required>("ni", scalar3d_layout_mid, n_unit, grid_name,
                      "tracers");  // ice number mixing ratio
}

// =========================================================================================
void PerturbTemp::initialize_impl(const RunType run_type) {
  // Gather runtime options
  //(e.g.) runtime_options.lambda_low    = m_params.get<double>("lambda_low");
}

// =========================================================================================
void PerturbTemp::run_impl(const double dt) {
  std::cout << "End of Perturb run" << std::endl;
}

// =========================================================================================
}  // namespace scream