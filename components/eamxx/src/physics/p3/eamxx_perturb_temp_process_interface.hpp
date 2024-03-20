#ifndef EAMXX_PERTURB_TEMP_HPP
#define EAMXX_PERTURB_TEMP_HPP

// For declaring a class derived from atm process class
#include <share/atm_process/atmosphere_process.hpp>

// For component name
#include <string>

namespace scream {

class PerturbTemp final : public scream::AtmosphereProcess {
  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;

 public:
  // Constructor
  PerturbTemp(const ekat::Comm &comm, const ekat::ParameterList &params);

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // The type of subcomponent
  AtmosphereProcessType type() const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name() const { return "PerturbTemp"; }

  // grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  // Initialize variables
  void initialize_impl(const RunType run_type) override;

  // Run the process by one time step
  void run_impl(const double dt) override;

  // Finalize
  void finalize_impl(){/*Do nothing*/};

};  // PerturbTemp

}  // namespace scream

#endif  // EAMXX_PERTURB_TEMP_HPP