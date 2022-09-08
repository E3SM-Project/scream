#ifndef EAMXX_POTENTIAL_TEMP_DIAGNOSTIC_HPP
#define EAMXX_POTENTIAL_TEMP_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class PotentialTemperatureDiagnostic : public AtmosphereDiagnostic
{
public:
  using Pack          = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using PF            = scream::PhysicsFunctions<DefaultDevice>;

  // Constructors
  PotentialTemperatureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Set type to diagnostic
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // The name of the diagnostic
  std::string name () const { return "Potential Temperature"; } 

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:

  void finalize_impl   () { /* Nothing to do */ }

  // Keep track of field dimensions
  Int m_num_cols; 
  Int m_num_levs;

}; // class PotentialTemperatureDiagnostic

} //namespace scream

#endif // EAMXX_POTENTIAL_TEMP_DIAGNOSTIC_HPP
