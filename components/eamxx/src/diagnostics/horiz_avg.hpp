#ifndef EAMXX_HORIZ_AVERAGE_HPP
#define EAMXX_HORIZ_AVERAGE_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will area-weighted average entries of the field across the COL tag dimension,
 * producing an N-1 dimensional field that is the area-weighted average of the input field.
 */

class HorizAvgDiag : public AtmosphereDiagnostic {
  using KT = ekat::KokkosTypes<DefaultDevice>;
  using const_view_1d = typename KT::template view_1d<const Real>;
 public:
  // Constructors
  HorizAvgDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic
  std::string name() const { return m_diag_name; }

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();

 protected:
  void initialize_impl(const RunType /*run_type*/);

  // Name of each field (because the diagnostic impl is generic)
  std::string m_diag_name;

  // Need grid for area field
  const_view_1d m_area;
};

}  // namespace scream

#endif  // EAMXX_HORIZ_AVERAGE_HPP
