#include "ekat/scream_assert.hpp"
#include "physics/zm/scream_zm_interface.hpp"
#include "physics/zm/atmosphere_macrophysics.hpp"
#include "physics/zm/zm_inputs_initializer.cpp"
#include <iostream>
#include <unordered_map>
#include <typeinfo>

namespace scream

{


const Int& lchnk = 0;
const Int& ncol = 0;
const Real &delt = 0;

const Real& lengath = 0; 
Real** t_star; 
Real** q_star;
Real** tend_s; 
Real** tend_q; 

Real** cld; 
Real** flxprec; 
Real** flxsnow;
const Real& ztodt = 0; 
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
  using namespace std;


  using namespace units;
  auto Q = kg/kg;
  auto nondim = m/m;
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
  
  unordered_map<string, GridOpts> map;

  	
  GridOpts t;
  t.name = "t";
  t.isOut = true;  
  map.insert({t.name, t});

  GridOpts qh;
  qh.name = "qh";
  qh.isOut = true;  
  map.insert({qh.name, qh});

  GridOpts prec;
  prec.name = "prec";
  prec.isOut = true;  
  map.insert({prec.name, prec});

  GridOpts jctop;
  jctop.name = "jctop";
  jctop.isOut = true;  
  map.insert({jctop.name, jctop});

  GridOpts jcbot;
  jcbot.name = "jcbot";
  jcbot.isOut = true;  
  map.insert({jcbot.name, jcbot});

  GridOpts pblh;
  pblh.name = "pblh";
  pblh.isOut = true;  
  map.insert({pblh.name, pblh});

  GridOpts zm;
  zm.name = "zm";
  zm.isOut = true;  
  map.insert({zm.name, zm});

  GridOpts geos;
  geos.name = "geos";
  geos.isOut = true;  
  map.insert({geos.name, geos});

  GridOpts zi;
  zi.name = "zi";
  zi.isOut = true;  
  map.insert({zi.name, zi});

  GridOpts qtnd;
  qtnd.name = "qtnd";
  qtnd.isOut = true;  
  map.insert({qtnd.name, qtnd});

  GridOpts heat;
  heat.name = "heat";
  heat.isOut = true;  
  map.insert({heat.name, heat});

  GridOpts pap;
  pap.name = "pap";
  pap.isOut = true;  
  map.insert({pap.name, pap});

  GridOpts paph;
  paph.name = "paph";
  paph.isOut = true;  
  map.insert({paph.name, paph});

  GridOpts dpp;
  dpp.name = "dpp";
  dpp.isOut = true;  
  map.insert({dpp.name, dpp});

  GridOpts mcon;
  mcon.name = "mcon";
  mcon.isOut = true;  
  map.insert({mcon.name, mcon});

  GridOpts cme;
  cme.name = "cme";
  cme.isOut = true;  
  map.insert({cme.name, cme});

  GridOpts cape;
  cape.name = "cape";
  cape.isOut = true;  
  map.insert({cape.name, cape});

  GridOpts tpert;
  tpert.name = "tpert";
  tpert.isOut = true;  
  map.insert({tpert.name, tpert});

  GridOpts dlf;
  dlf.name = "dlf";
  dlf.isOut = true;  
  map.insert({dlf.name, dlf});

  GridOpts pflx;
  pflx.name = "pflx";
  pflx.isOut = true;  
  map.insert({pflx.name, pflx});

  GridOpts zdu;
  zdu.name = "zdu";
  zdu.isOut = true;  
  map.insert({zdu.name, zdu});

  GridOpts rprd;
  rprd.name = "rprd";
  rprd.isOut = true;  
  map.insert({rprd.name, rprd});
  
  GridOpts mu;
  mu.name = "mu";
  mu.isOut = true;  
  map.insert({mu.name, mu});
  
  GridOpts md;
  md.name = "md";
  md.isOut = true;  
  map.insert({md.name, md});
 // 
  GridOpts du;
  du.name = "du";
  du.isOut = true;  
  map.insert({du.name, du});
  
  GridOpts eu;
  eu.name = "eu";
  eu.isOut = true;  
  map.insert({eu.name, eu});
  
  
  GridOpts ed;
  ed.name = "ed";
  ed.isOut = true;  
  map.insert({ed.name, ed});
  
  GridOpts dp;
  dp.name = "dp";
  dp.isOut = true;  
  map.insert({dp.name, dp});
  
  GridOpts dsubcld;
  dsubcld.name = "dsubcld";
  dsubcld.isOut = true;  
  map.insert({dsubcld.name, dsubcld});
  
  GridOpts jt;
  jt.name = "jt";
  jt.isOut = true;  
  map.insert({jt.name, jt});
  
  GridOpts maxg;
  maxg.name = "maxg";
  maxg.isOut = true;  
  map.insert({maxg.name, maxg});
  
  GridOpts ideep;
  ideep.name = "ideep";
  ideep.isOut = true;  
  map.insert({ideep.name, ideep});
  
  GridOpts ql;
  ql.name = "ql";
  ql.isOut = true;  
  map.insert({ql.name, ql});
  
  GridOpts rliq;
  rliq.name = "rliq";
  rliq.isOut = true;  
  map.insert({rliq.name, rliq});
  
  GridOpts landfrac;
  landfrac.name = "landfrac";
  landfrac.isOut = true;  
  map.insert({landfrac.name, landfrac});
  
  GridOpts hu_nm1;
  hu_nm1.name = "hu_nm1";
  hu_nm1.isOut = true;  
  map.insert({hu_nm1.name, hu_nm1});
  
  GridOpts cnv_nm1;
  cnv_nm1.name = "cnv_nm1";
  cnv_nm1.isOut = true;  
  map.insert({cnv_nm1.name, cnv_nm1});
  
  GridOpts tm1;
  tm1.name = "tm1";
  tm1.isOut = true;  
  map.insert({tm1.name, tm1});
  
  GridOpts qm1;
  qm1.name = "qm1";
  qm1.isOut = true;  
  map.insert({qm1.name, qm1});
  
  GridOpts dcape;
  dcape.name = "dcape";
  dcape.isOut = true;  
  map.insert({dcape.name, dcape});
  
  GridOpts q;
  q.name = "q";
  q.isOut = true;  
  map.insert({q.name, q});
  
  GridOpts snow;
  snow.name = "snow";
  snow.isOut = true;  
  map.insert({snow.name, snow});
  
  GridOpts ntprprd;
  ntprprd.name = "ntprprd";
  ntprprd.isOut = true;  
  map.insert({ntprprd.name, ntprprd});
  
  GridOpts ntsnprd;
  ntsnprd.name = "ntsnprd";
  ntsnprd.isOut = true;  
  map.insert({ntsnprd.name, ntsnprd});
  
  GridOpts pguall;
  pguall.name = "pguall";
  pguall.isOut = true;  
  map.insert({pguall.name, pguall});
  
  GridOpts pgdall;
  pgdall.name = "pgdall";
  pgdall.isOut = true;  
  map.insert({pgdall.name, pgdall});
  
  GridOpts icwu;
  icwu.name = "icwu";
  icwu.isOut = true;  
  map.insert({icwu.name, icwu});
  
  for ( auto i = map.begin(); i != map.end(); ++i){
    cout << "the name is " << (i->second).name << "\n";
    m_required_fields.emplace((i->second).name, scalar3d_layout_int, Q, grid->name());
    if ( (i->second).isOut == true ) {
      cout << "the name is " << (i->second).name << "\n";
      m_computed_fields.emplace((i->second).name, scalar3d_layout_int, Q, grid->name());
    }
  }
 

  
}

// =========================================================================================
void ZMMacrophysics::initialize (const util::TimeStamp& t0)
{
  m_current_ts = t0;
  
  zm_init_f90 (limcnv_in, no_deep_pbl_in);
  
  using strvec = std::vector<std::string>;
  const strvec& initable = zm_inputs;
  for (const auto& name : initable) {
    const auto& f = m_zm_fields_in.at(name);
    const auto& track = f.get_header().get_tracking();
    if (track.get_init_type()==InitType::None) {
      // Nobody claimed to init this field. P3InputsInitializer will take care of it
      m_initializer->add_me_as_initializer(f);
    }
   }
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

   
  zm_main_f90(lchnk, ncol, m_raw_ptrs_out["t"], m_raw_ptrs_out["qh"], m_raw_ptrs_out["prec"],
              m_raw_ptrs_out["jctop"], m_raw_ptrs_out["jcbot"], m_raw_ptrs_out["pblh"], 
              m_raw_ptrs_out["zm"], m_raw_ptrs_out["geos"], m_raw_ptrs_out["zi"], 
	      m_raw_ptrs_out["qtnd"], m_raw_ptrs_out["heat"], m_raw_ptrs_out["pap"], 
	      m_raw_ptrs_out["paph"], m_raw_ptrs_out["dpp"], delt, m_raw_ptrs_out["mcon"], 
	      m_raw_ptrs_out["cme"], m_raw_ptrs_out["cape"], m_raw_ptrs_out["tpert"], 
              m_raw_ptrs_out["dlf"], m_raw_ptrs_out["plfx"], m_raw_ptrs_out["zdu"], 
	      m_raw_ptrs_out["rprd"], m_raw_ptrs_out["mu"], m_raw_ptrs_out["md"], 
              m_raw_ptrs_out["du"], m_raw_ptrs_out["eu"], m_raw_ptrs_out["ed"], 
	      m_raw_ptrs_out["dp"], m_raw_ptrs_out["dsubcld"], m_raw_ptrs_out["jt"], 
	      m_raw_ptrs_out["maxg"], m_raw_ptrs_out["ideep"], lengath, m_raw_ptrs_out["ql"],
              m_raw_ptrs_out["rliq"], m_raw_ptrs_out["landfrac"], m_raw_ptrs_out["hu_nm1"], 
              m_raw_ptrs_out["cnv_nm1"], m_raw_ptrs_out["tm1"], m_raw_ptrs_out["qm1"], 
              t_star, q_star, m_raw_ptrs_out["dcape"], m_raw_ptrs_out["q"], 
              tend_s, tend_q, cld, m_raw_ptrs_out["snow"], m_raw_ptrs_out["ntprprd"], 
              m_raw_ptrs_out["ntsnprd"], flxprec, flxsnow, ztodt, m_raw_ptrs_out["pguall"], 
	      m_raw_ptrs_out["pgdall"], m_raw_ptrs_out["icwu"], ncnst, fracis); 
  m_current_ts += dt;
  for (int i = 0; i < zm_inputs.size(); i++){
  	m_zm_fields_out.at(zm_inputs[i]).get_header().get_tracking().update_time_stamp(m_current_ts);
  }
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
