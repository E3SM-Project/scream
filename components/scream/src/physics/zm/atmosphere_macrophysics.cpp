#include "ekat/scream_assert.hpp"
#include "physics/zm/zm_inputs_initializer.hpp"
#include "physics/zm/scream_zm_interface.hpp"
#include "physics/zm/atmosphere_macrophysics.hpp"
#include <iostream>
namespace scream

{
std::vector<std::string> zm_inputs = {"t"};


const Int& lchnk = 0;
const Int& ncol = 0;
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

Real*** fracis;


ZMMacrophysics::ZMMacrophysics (const Comm& comm,const ParameterList& /* params */)
  : m_zm_comm (comm)
{
  m_initializer = create_field_initializer<ZMInputsInitializer>();
    
}
void ZMMacrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{


  using namespace units;
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  
  constexpr int NVL = 72;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */
  constexpr int QSZ =  35;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */
  auto grid = grids_manager->get_grid("Physics");
  const int num_dofs = grid->get_num_local_dofs();
  const int nc = num_dofs;


  auto VL = FieldTag::VerticalLevel;
  auto CO = FieldTag::Column;
  auto VR = FieldTag::Variable;
  
  FieldLayout scalar3d_layout { {CO,VL}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout { {CO,VR,VL}, {nc,QSZ,NVL} };
  FieldLayout q_forcing_layout  { {CO,VR,VL}, {nc,QSZ,NVL} };
  auto nondim = m/m;

  m_required_fields.emplace("t",  vector3d_layout, Q, grid->name()); // in/out??
  m_computed_fields.emplace("t",  vector3d_layout, Q, grid->name()); // in/out??

  // We may have to init some fields from within P3. This can be the case in a P3 standalone run.
  // Some options:
  //  - we can tell P3 it can init all inputs or specify which ones it can init. We call the
  //    resulting list of inputs the 'initializaable' (or initable) inputs. The default is
  //    that no inputs can be inited.
  //  - we can request that P3 either inits no inputs or all of the initable ones (as specified
  //    at the previous point). The default is that P3 must be in charge of init ing ALL or NONE
  //    of its initable inputs.
  // Recall that:
  //  - initable fields may not need initialization (e.g., some other atm proc that
  //    appears earlier in the atm dag might provide them).

  using strvec = std::vector<std::string>;
//  const strvec& allowed_to_init = m_zm_params.get<strvec>("Initializable Inputs",strvec(0));
//  const bool can_init_all = m_zm_params.get<bool>("Can Initialize All Inputs", false);
//  const bool init_all_or_none = m_zm_params.get<bool>("Must Init All Inputs Or None", true);
//  
//  const strvec& initable = can_init_all ? zm_inputs : allowed_to_init;
  if (zm_inputs.size()>0) {
    bool all_inited = true, all_uninited = true;
    for (const auto& name : zm_inputs) {
      const auto& f = m_zm_fields_in.at(name);
      const auto& track = f.get_header().get_tracking();
      if (track.get_init_type()==InitType::None) {
        // Nobody claimed to init this field. P3InputsInitializer will take care of it
        m_initializer->add_me_as_initializer(f);
        all_uninited &= true;
        all_inited &= false;
      } else {
        all_uninited &= false;
        all_inited &= true;
      }
    }

    // In order to gurantee some consistency between inputs, it is best if P3
    // initializes either none or all of the inputs.
   // scream_require_msg (!init_all_or_none || all_inited || all_uninited,
   //                     "Error! Some zm inputs were marked to be inited by P3, while others weren't.\n"
   //                     "       P3 was requested to init either all or none of the inputs.\n");
  }//

}

//void add_field_ptr(string input){
//  auto dev = m_zm_fields_out.at(input).get_view();
//  auto host = Kokkos::create_mirror_view(q_dev);
//  Kokkos::deep_copy(q_host,q_dev);
//  auto ptr = q_host.data();
//  return ptr
//
//}


// =========================================================================================
void ZMMacrophysics::initialize (const util::TimeStamp& t0)
{
  std::vector<std::string> zm_inputs = {"t"};
  
  zm_init_f90 (limcnv_in, no_deep_pbl_in);

}



// =========================================================================================
void ZMMacrophysics::run (const Real dt)
{
   
   zm_main_f90(lchnk, ncol, t, qh, prec, jctop, jcbot, pblh, zm, geos, zi, qtnd, 
		heat, pap, paph, dpp, delt,
		mcon, cme, cape, tpert, dlf, plfx,
		zdu, rprd, mu, md, du, eu, ed, dp, dsubcld, jt, maxg, ideep,
		lengath, ql, rliq, landfrac, hu_nm1, cnv_nm1,
		tm1, qm1, t_star, q_star, dcape, q, tend_s, tend_q, cld,
		snow, ntprprd, ntsnprd, flxprec, flxsnow,
		ztodt, pguall, pgdall, icwu, ncnst, fracis); 
}
// =========================================================================================
void ZMMacrophysics::finalize()
{
  zm_finalize_f90 ();
}
// =========================================================================================
void ZMMacrophysics::register_fields (FieldRepository<Real, device_type>& field_repo) const {
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
  const auto& name = f.get_header().get_identifier().name();
  m_zm_fields_in.emplace(name,f);

  // Add myself as customer to the field
  // @Meredith: Diff between add_me_as_a_customer and get_tracking().add_customer? 
  add_me_as_customer(f);
}

void ZMMacrophysics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_zm_fields_out.emplace(name,f);
  // Add myself as provider for the field
  add_me_as_provider(f);
}
} // namespace scream
