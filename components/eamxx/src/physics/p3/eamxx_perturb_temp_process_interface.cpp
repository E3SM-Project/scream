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

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto q_unit = kg / kg;
  q_unit.set_string("kg/kg");

  auto n_unit = 1 / kg;  // units of number mixing ratios of tracers
  n_unit.set_string("#/kg");

  auto m3 = m * m * m;  // meter cubed

  grid_                 = grids_manager->get_grid("Physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column
}

// =========================================================================================
void PerturbTemp::initialize_impl(const RunType run_type) {
  // Gather runtime options
  //(e.g.) runtime_options.lambda_low    = m_params.get<double>("lambda_low");
}

// =========================================================================================
void PerturbTemp::run_impl(const double dt) {
  std::cout << "End of derydep run" << std::endl;
}

// =========================================================================================
}  // namespace scream