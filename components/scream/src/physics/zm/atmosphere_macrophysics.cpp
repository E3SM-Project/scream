#include "ekat/scream_assert.hpp"
#include "physics/zm/scream_zm_interface.hpp"
#include "physics/zm/atmosphere_macrophysics.hpp"
#include <iostream>
namespace scream

{


const Real& lchnk = 0;
const Real& ncol = 0;
Real* t;
Real* qh;
Real* prec;
Real* jctop; 
Real* jcbot; 
Real* pblh;
Real *zm; 
Real* geos; 
Real* zi;
			
Real* qtnd; 
Real* heat; 
Real* pap; 
Real* paph; 
Real* dpp; 
const Real &delt = 0;

Real* mcon;
Real* cme;
Real* cape;
Real* tpert;
Real* dlf;
Real* plfx;
Real* zdu;
Real* rprd; 
Real* mu;
Real* md; 
Real* du;
Real* eu; 
Real* ed;
Real* dp; 
Real* dsubcld; 
Real* jt; 
Real* maxg; 
Real* ideep;
const Real& lengath = 0; 
Real* ql; 
Real* rliq; 
Real* landfrac; 
Real* hu_nm1;
Real* cnv_nm1; 
Real* tm1; 
Real* qm1; 
Real** t_star; 
Real** q_star;
Real *dcape; 
Real* q; 
Real** tend_s; 
Real** tend_q; 
Real** cld; 
Real* snow; 
Real* ntprprd; 
Real* ntsnprd; 
Real** flxprec; 
Real** flxsnow;
const Real& ztodt = 0; 
Real* pguall; 
Real* pgdall; 
Real* icwu; 
const Real& ncnst = 0; 
const Real& limcnv_in = 0;
const bool& no_deep_pbl_in = true;


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
  zm_init_f90 (limcnv_in, no_deep_pbl_in);

}

// =========================================================================================
void ZMMacrophysics::run (const Real dt)
{
   
//  zm_main_f90	(lchnk, ncol, t, qh, prec,
//			jctop, jcbot, pblh, zm, geos, zi,
//			qtnd, heat, pap, paph, dpp, delt,
//			mcon, cme, cape, tpert, dlf, plfx,
//			zdu, rprd, mu, md, du, eu, 
//			ed, dp, dsubcld, jt, maxg, ideep,
//			lengath, ql, rliq, landfrac, hu_nm1,
//			cnv_nm1, tm1, qm1, t_star, q_star,
//			dcape, q, tend_s, tend_q, cld, 
//			snow, ntprprd, ntsnprd, flxprec,flxsnow,
//			ztodt, pguall, pgdall, icwu, ncnst );
   zm_main_f90(lchnk, ncol, t, qh, prec, jctop, jcbot, pblh, zm, geos, zi, qtnd, 
		heat, pap, paph, dpp, delt,
		mcon, cme, cape, tpert, dlf, plfx,
		zdu, rprd, mu, md, du, eu, 
		ed, dp, dsubcld, jt, maxg, ideep,
		lengath, ql, rliq, landfrac, hu_nm1, cnv_nm1,
		tm1, qm1, t_star, q_star, dcape, q, tend_s, tend_q, cld,
		snow, ntprprd, ntsnprd, flxprec, flxsnow 
//		tm1, qm1, t_star, q_star
); 
//		snow, ntprprd, ntsnprd, flxprec,flxsnow,
//		ztodt, pguall, pgdall, icwu, ncnst);
   std :: cout << "In ZMMAcrophysics::run\n"<< std::endl;
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
