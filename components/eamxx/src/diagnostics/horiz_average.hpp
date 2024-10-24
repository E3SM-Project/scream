#ifndef EAMXX_HORIZ_AVERAGE_HPP
#define EAMXX_HORIZ_AVERAGE_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will sum entries of the field across the COL tag dimension,
 * producing an N-1 dimensional field
 */

class HorizAverage : public AtmosphereDiagnostic
{
public:
  // Constructors
  HorizAverage (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic
  std::string name () const { return m_diag_name; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:
  void initialize_impl (const RunType /*run_type*/);

  std::string m_diag_name;
};

} // namespace scream

#endif // EAMXX_HORIZ_AVERAGE_HPP
