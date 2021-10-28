#ifndef SCREAM_CLD_FRACTION_HPP
#define SCREAM_CLD_FRACTION_HPP

#include "physics/cld_fraction/cld_fraction_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the calculation of the subgrid cloud fractions
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class CldFraction : public AtmosphereProcess
{
public:
  using field_type       = Field<      Real>;
  using const_field_type = Field<const Real>;

  using CldFractionFunc = cld_fraction::CldFractionFunctions<Real, DefaultDevice>;
  using Spack           = CldFractionFunc::Spack;
  using Smask           = CldFractionFunc::Smask;
  using Pack            = ekat::Pack<Real,Spack::n>;

  // Constructors
  CldFraction (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "CldFraction"; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:

  // The three main overrides for the subcomponent
  void initialize_impl ();
  void run_impl        (const int dt);
  void finalize_impl   ();

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_levs;

}; // class CldFraction

} // namespace scream

#endif // SCREAM_CLD_FRACTION_HPP
