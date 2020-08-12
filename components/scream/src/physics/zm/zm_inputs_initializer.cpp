#include "physics/zm/zm_inputs_initializer.hpp"
#include "physics/zm/scream_zm_interface.hpp"

#include <array>
#include <string>
#include <iostream>
#include <typeinfo>

namespace scream
{

using namespace std;

const int INPUT_SIZE = 63;

vector<string> zm_inputs = {"limcnv_in", "no_deep_pbl_in","lchnk", "ncol", "t", "qh", "prec", 
				"jctop", "jcbot", "pblh", "zm", "geos", "zi", "qtnd", "heat", 
				"pap", "paph", "dpp", "delt", "mcon", "cme", "cape", "tpert",
				"dlf", "pflx", "zdu", "rprd", "mu", "md", "du", "eu", "ed", 
				"dp", "dsubcld", "jt", "maxg", "ideep", "lengath", "ql", 
				"rliq", "landfrac", "hu_nm1", "cnv_nm1", "tm1", "qm1", "t_star", 
				"q_star", "dcape", "q", "tend_s", "tend_q", "cld", "snow", 
				"ntprprd", "ztodt", "ntsnprd", "flxprec", "flxsnow", "pguall", 
				"pgdall", "icwu", "ncnst", "fracis"};


void ZMInputsInitializer::add_field (const field_type &f)
{
  const auto& id = f.get_header().get_identifier();
  
  m_fields.emplace(id.name(),f);
  m_fields_id.insert(id);
}



void ZMInputsInitializer :: initialize_fields(){

  int count = 0;
  for (int j = 0; j < zm_inputs.size(); j++){
    count += m_fields.count(zm_inputs[j]);
    if (count==0){
      return;
    }
  } 
  scream_require_msg (count==INPUT_SIZE,
    "Error! ZMInputsInitializer is expected to init _______________.\n"
    "Only " + to_string(count) + " of those have been found.\n"
    "Please, check the atmosphere processes you are using,"
    "and make sure they agree on who's initializing each field.\n");



  for (int i = 0; i < zm_inputs.size(); i++){
    //Get and store device view using input name
    Kokkos::View<Real*, Kokkos::LayoutRight, HostDevice> d_v = m_fields.at(zm_inputs[i]).get_view();

    //Create and store host mirrors using device views
    Kokkos::View<scream::Real*, Kokkos::LayoutRight, HostDevice> h_m = Kokkos::create_mirror_view(d_v);
    //Create and store host mirrors raw pointers
    Real* r_p = h_m.data();
    Kokkos::deep_copy(d_v, h_m);
  }
}


} // namespace scream


