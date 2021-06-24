#ifndef SCREAM_PRESCRIBED_AERO_HPP
#define SCREAM_PRESCRIBED_AERO_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the Simple Prescribed Aerosols
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class SPA : public AtmosphereProcess
{
public:
  using field_type       = Field<      Real>;
  using const_field_type = Field<const Real>;

//  using CldFractionFunc = cld_fraction::CldFractionFunctions<Real, DefaultDevice>;
//  using Spack           = CldFractionFunc::Spack;
//  using Smask           = CldFractionFunc::Smask;
//  using Pack            = ekat::Pack<Real,Spack::n>;

  // Constructors
  SPA (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "SPA"; }

  // The communicator used by subcomponent
  const ekat::Comm& get_comm () const { return m_spa_comm; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_spa_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const util::TimeStamp& t0);
  void run_impl        (const Real dt);
  void finalize_impl   ();

  // Setting the fields in the atmospheric process
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);

  std::map<std::string,const_field_type>  m_spa_fields_in;
  std::map<std::string,field_type>        m_spa_fields_out;

  ekat::Comm          m_spa_comm;
  ekat::ParameterList m_spa_params;

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_levs;

  // The name of the SPA data file
  std::string m_input_data_filename; 

}; // class SPA

} // namespace scream

#endif // SCREAM_PRESCRIBED_AERO_HPP
