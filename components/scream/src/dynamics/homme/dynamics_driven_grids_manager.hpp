#ifndef SCREAM_DYNAMICS_DRIVEN_GRIDS_MANAGER_HPP
#define SCREAM_DYNAMICS_DRIVEN_GRIDS_MANAGER_HPP

#include "share/grid/grids_manager.hpp"

namespace scream
{

// This class is meant to be used for small unit tests, where we want to
// test the Atmosphere Driver (AD) capabilities, without bothering too much
// about grids-related features. This manager lets you set pre-built grids
// in it rather than building them inside the manager.
class DynamicsDrivenGridsManager : public GridsManager
{
public:

  DynamicsDrivenGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p);

  ~DynamicsDrivenGridsManager ();

  std::string name () const { return "Dynamics Driven Grids Manager"; }

  void build_grids (const std::set<std::string>& grid_names);

  std::set<std::string> supported_grids () const { return m_valid_grid_names; }

  const grid_repo_type& get_repo () const { return m_grids; }

#ifndef KOKKOS_ENABLE_CUDA
protected:
#endif

  void build_dynamics_grid (const bool DG);
  void build_physics_grid  (const std::string& name);

protected:

  std::string get_reference_grid_name () const {
    return m_ref_grid_name;
  }

  remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const;

  void build_grid_codes ();

  ekat::Comm      m_comm;

  grid_repo_type  m_grids;

  std::string     m_ref_grid_name;

  // Admissible grid names
  std::set<std::string> m_valid_grid_names;

  // For each admissible grid name, store an integer code
  std::map<std::string, int> m_grid_codes;
};

inline std::shared_ptr<GridsManager>
create_dynamics_driven_grids_manager (const ekat::Comm& comm, const ekat::ParameterList& p) {
  return std::make_shared<DynamicsDrivenGridsManager>(comm,p);
}

} // namespace scream

#endif // SCREAM_DYNAMICS_DRIVEN_GRIDS_MANAGER_HPP
