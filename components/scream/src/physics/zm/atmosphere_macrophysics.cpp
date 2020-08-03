#include "ekat/scream_assert.hpp"
#include "physics/zm/scream_zm_interface.hpp"
#include "physics/zm/atmosphere_macrophysics.hpp"
#include "physics/zm/zm_inputs_initializer.cpp"
#include <iostream>
namespace scream

{
//std::vector<std::string> zm_inputs = {"t", "qh," "prec", "jctop", "jcbot", "pblh", "zm",
//				"geos", "zi", "qtnd", "heat", "pap", "paph", "dpp", "delt", 
//				"mcon", "cme", "cape", "tpert", "dlf", "pflx", "zdu", "rprd",
//				"mu", "md", "du", "eu", "ed", "dp", "dsubcld", "jt", "maxg",
//				"ideep", "lengath", "ql", "rliq", "landfrac", "hu_nm1",
//      				"cnv_nm1", "tm1", "qm1", "t_star", "q_star", "dcape", "q",
//				"tend_s", "tend_q", "cld", "snow", "ntprprd", "ntsnprd",
//				"flxprec", "flxsnow", "ztodt", "pguall", "pgdall", "icwu",
//				 "ncnst", "fracis"};
//
//

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

  using namespace ShortFieldTagsNames;
 
  FieldLayout scalar3d_layout_mid { {COL,VL}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout scalar3d_layout_int { {COL,VL}, {nc,NVL+1} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout_mid{ {COL,CMP,VL}, {nc,QSZ,NVL} };
  FieldLayout tracers_layout { {COL,VAR,VL}, {nc,QSZ,NVL} };
  
  auto nondim = m/m;

  m_required_fields.emplace("t",  scalar3d_layout_int, Q, grid->name()); // in/out??
  m_computed_fields.emplace("t",  scalar3d_layout_int, Q, grid->name()); // in/out??

}

// =========================================================================================
void ZMMacrophysics::initialize (const util::TimeStamp& t0)
{
  
  zm_init_f90 (limcnv_in, no_deep_pbl_in);
//  using strvec = std::vector<std::string>;
//  const strvec& initable = zm_inputs;
//  for (const auto& name : initable) {
//    const auto& f = m_zm_fields_in.at(name);
//    const auto& track = f.get_header().get_tracking();
//    if (track.get_init_type()==InitType::None) {
//      // Nobody claimed to init this field. P3InputsInitializer will take care of it
//      m_initializer->add_me_as_initializer(f);
//    }
//   }
}

// =========================================================================================
void ZMMacrophysics::run (const Real dt)
{
  std::vector<const Real*> in;
  std::vector<Real*> out;

  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_zm_fields_in) {
    Kokkos::deep_copy(m_zm_host_views_in.at(it.first),it.second.get_view());
  }
  for (auto& it : m_zm_fields_out) {
    Kokkos::deep_copy(m_zm_host_views_out.at(it.first),it.second.get_view());
  }

   
   zm_main_f90(lchnk, ncol, m_raw_ptrs_out["t"], qh, prec, jctop, jcbot, 
                pblh, zm, geos, zi, qtnd, heat, pap, paph, dpp, delt,
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
  // @Meredith: Diff between add_me_as_a_customer and get_tracking().add_customer? 
  
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_zm_fields_in.emplace(name,f);
  m_zm_host_views_in[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_in[name] = m_zm_host_views_in[name].data();

  // Add myself as customer to the field
  add_me_as_customer(f);

}

void ZMMacrophysics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_zm_fields_out.emplace(name,f);
  m_zm_host_views_out[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_out[name] = m_zm_host_views_out[name].data();

  // Add myself as provider for the field
  add_me_as_provider(f);
  }
} // namespace scream
