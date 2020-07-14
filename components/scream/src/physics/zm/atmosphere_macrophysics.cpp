#include "ekat/scream_assert.hpp"
#include "physics/zm/scream_zm_interface.hpp"
#include "physics/zm/atmosphere_macrophysics.hpp"
#include <iostream>
namespace scream

{
ZMMacrophysics::ZMMacrophysics (const Comm& comm,const ParameterList& /* params */)
  : m_zm_comm (comm)
{
}
void ZMMacrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
}

// =========================================================================================
void ZMMacrophysics::initialize (const util::TimeStamp& t0)
{
//  std :: cout << "In ZMMacrophysics::initialize \n"<< std::endl;
  //zm_init_f90 (q_ptr);
  zm_init_f90 ();

}

// =========================================================================================
void ZMMacrophysics::run (const Real dt)
{
 // std :: cout << "In ZMMAcrophysics::run\n"<< std::endl;
  zm_main_f90 ();
  //zm_main_f90 (dt,q_ptr,fq_ptr);
}
// =========================================================================================
void ZMMacrophysics::finalize()
{
 // std :: cout << "In ZMMacrophysics :: finalize()HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE \n"<< std::endl;
  zm_finalize_f90 ();
}
// =========================================================================================
void ZMMacrophysics::register_fields (FieldRepository<Real, device_type>& field_repo) const {
     std :: cout << ".......\n";
     for (auto& fid : m_required_fields) {
     field_repo.register_field(fid);
   }
   for (auto& fid : m_computed_fields) {
     field_repo.register_field(fid);
   }
 }

void ZMMacrophysics::set_required_field_impl (const Field<const Real, device_type>& f) {
   // Store a copy of the field. We need this in order to do some tracking checks
   // at the beginning of the run call. Other than that, there would be really
   // no need to store a scream field here; we could simply set the view ptr
   // in the Homme's view, and be done with it.
   m_zm_fields_in.emplace(f.get_header().get_identifier().name(),f);
 
   // Add myself as customer to the field
   f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
 }
void ZMMacrophysics::set_computed_field_impl (const Field<      Real, device_type>& f) {
   // Store a copy of the field. We need this in order to do some tracking updates
   // at the end of the run call. Other than that, there would be really
   // no need to store a scream field here; we could simply set the view ptr
   // in the Homme's view, and be done with it.
   m_zm_fields_out.emplace(f.get_header().get_identifier().name(),f);
 
   // Add myself as provider for the field
   f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}
} // namespace scream
