#include "physics/zm/zm_inputs_initializer.hpp"
#include "physics/zm/scream_zm_interface.hpp"

#include <array>
#include <string>
#include <iostream>
#include <typeinfo>
#include <unordered_map>

namespace scream
{

using namespace std;

const int INPUT_SIZE = 22;

std::vector<std::string> zm_inputs = {"t", "qh", "prec", "jctop", "jcbot", "pblh", "zm",
				"geos", "zi", "qtnd", "heat", "pap", "paph", "dpp",// ,"delt"
				"mcon", "cme", "cape", "tpert", "dlf", "pflx", "zdu", "rprd"};
//				"mu", "md", "du", "eu", "ed", "dp", "dsubcld", "jt", "maxg",
//				"ideep", "lengath", "ql", "rliq", "landfrac", "hu_nm1",
//      				"cnv_nm1", "tm1", "qm1", "t_star", "q_star", "dcape", "q",
//				"tend_s", "tend_q", "cld", "snow", "ntprprd", "ntsnprd",
//				"flxprec", "flxsnow", "ztodt", "pguall", "pgdall", "icwu",
//				 "ncnst", "fracis"};

struct FieldInput{
  string name;
  Kokkos::View<Real*, Kokkos::LayoutRight, HostDevice> device_view;
  Real *raw_ptr;
  Kokkos::View<scream::Real*, Kokkos::LayoutRight, HostDevice> mirror;
  bool is_output; 
};


void ZMInputsInitializer::add_field (const field_type &f)
{
  const auto& id = f.get_header().get_identifier();
  
  m_fields.emplace(id.name(),f);
  m_fields_id.insert(id);
}

unordered_map<std::string, FieldInput> map;

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
    "Only " + std::to_string(count) + " of those have been found.\n"
    "Please, check the atmosphere processes you are using,"
    "and make sure they agree on who's initializing each field.\n");



  for (int i = 0; i < zm_inputs.size(); i++){
    FieldInput field;
    field.name = zm_inputs[i];
    //Get and store device view using input name
    Kokkos::View<Real*, Kokkos::LayoutRight, HostDevice> d_v = m_fields.at(field.name).get_view();

    //Create and store host mirrors using device views
    Kokkos::View<scream::Real*, Kokkos::LayoutRight, HostDevice> h_m = Kokkos::create_mirror_view(d_v);

    //Create and store host mirrors raw pointers
    Real* r_p = field.mirror.data();
    Kokkos::deep_copy(d_v, h_m);
    map[zm_inputs[i]] = field;
  }
}


} // namespace scream




//  for (int i = 0; i < zm_inputs.size(); i++){
//    FieldInput field;
//    field.name = zm_inputs[i];
//    //Get and store device view using input name
//    Kokkos::View<Real*, Kokkos::LayoutRight, HostDevice> d_v = m_fields.at(field.name).get_view();
//   field.device_view = d_v;
//
//    //Create and store host mirrors using device views
//    Kokkos::View<scream::Real*, Kokkos::LayoutRight, HostDevice> h_m = Kokkos::create_mirror_view(d_v);
//   field.mirror = h_m;
//
//    //Create and store host mirrors raw pointers
//    Real* r_p = field.mirror.data();
//   field.raw_ptr = r_p;
//    Kokkos::deep_copy(d_v, h_m);
//    map[zm_inputs[i]] = field;
//  }
//}


