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

const int INPUT_SIZE = 1;

std::vector<std::string> zm_inputs = {"t"};

struct FieldInput{
  string name;
  Kokkos::View<Real*, Kokkos::LayoutRight, HostDevice> device_view;
  Real *raw_ptr;
  Kokkos::View<scream::Real*, Kokkos::LayoutRight, HostDevice> mirror;
  bool is_output; 
};


void ZMInputsInitializer::add_field (const field_type &f)
{
  std :: cout << "here i am!\n";
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


//std::vector<std::string> zm_inputs = {"t", "qh," "prec", "jctop", "jcbot", "pblh", "zm",
//				"geos", "zi", "qtnd", "heat", "pap", "paph", "dpp", "delt", 
//				"mcon", "cme", "cape", "tpert", "dlf", "pflx", "zdu", "rprd",
//				"mu", "md", "du", "eu", "ed", "dp", "dsubcld", "jt", "maxg",
//				"ideep", "lengath", "ql", "rliq", "landfrac", "hu_nm1",
//      				"cnv_nm1", "tm1", "qm1", "t_star", "q_star", "dcape", "q",
//				"tend_s", "tend_q", "cld", "snow", "ntprprd", "ntsnprd",
//				"flxprec", "flxsnow", "ztodt", "pguall", "pgdall", "icwu",
//				 "ncnst", "fracis"};










// =========================================================================================
//void ZMInputsInitializer::initialize_fields ()
//{
//  // Safety check: if we're askedfd to init anything at all,
//  // then we should have been asked to init 7 fields.
//  int count = 0;
//  count += m_fields.count("q");
//  count += m_fields.count("T");
//  count += m_fields.count("ast");
//  count += m_fields.count("naai");
//  count += m_fields.count("ncnuc");
//  count += m_fields.count("pmid");
//  count += m_fields.count("dp");
//  count += m_fields.count("zi");
//
//
//  for (int i = 0 ; i < zm_inputs.size(); i++ ){
//	count += m_fields.count(zm_inputs[i]);
//  }
//  std::cout << "count is " + std :: to_string(count);
//  if (count==0) {
//    return;
//  }
//
//  scream_require_msg (count==8,
//    "Error! ZMInputsInitializer is expected to init 'q','T','ast','naai','ncnuc','pmid','dp','zi'.\n"
//    "       Only " + std::to_string(count) + " of those have been found.\n"
//    "       Please, check the atmosphere processes you are using,"
//    "       and make sure they agree on who's initializing each field.\n");
//
//
////  for (int i = 0 ; i < zm_inputs.size(); i++ ){
////    Kokkos
//  // Get device views
//  auto d_q     = m_fields.at("q").get_view();
//  auto d_T     = m_fields.at("T").get_view();
//  auto d_ast   = m_fields.at("ast").get_view();
//  auto d_naai  = m_fields.at("naai").get_view();
//  auto d_ncnuc = m_fields.at("ncnuc").get_view();
//  auto d_pmid  = m_fields.at("pmid").get_view();
//  auto d_pdel  = m_fields.at("dp").get_view();
//  auto d_zi    = m_fields.at("zi").get_view();
//
//
//  // Create host mirrors
//  auto h_q     = Kokkos::create_mirror_view(d_q);
//  auto h_T     = Kokkos::create_mirror_view(d_T);
//  auto h_ast   = Kokkos::create_mirror_view(d_ast);
//  auto h_naai  = Kokkos::create_mirror_view(d_naai);
//  auto h_ncnuc = Kokkos::create_mirror_view(d_ncnuc);
//  auto h_pmid  = Kokkos::create_mirror_view(d_pmid);
//  auto h_pdel  = Kokkos::create_mirror_view(d_pdel); auto h_zi    = Kokkos::create_mirror_view(d_zi);
//  auto h_zi    = Kokkos::create_mirror_view(d_zi);
//
//  // Get host mirros' raw pointers
//  auto q     = h_q.data();
//  auto T     = h_T.data();
//  auto ast   = h_ast.data();
//  auto naai  = h_naai.data();
//  auto ncnuc = h_ncnuc.data();
//  auto pmid  = h_pmid.data();
//  auto pdel  = h_pdel.data();
//  auto zi    = h_zi.data();
//
//  // Call f90 routine
////  zm_standalone_init_f90 (q, T, zi, pmid, pdel, ast, naai, ncnuc);
//  
//  // Deep copy back to device
//  Kokkos::deep_copy(d_q,h_q);
//  Kokkos::deep_copy(d_T,h_T);
//  Kokkos::deep_copy(d_ast,h_ast);
//  Kokkos::deep_copy(d_naai,h_naai);
//  Kokkos::deep_copy(d_ncnuc,h_ncnuc);
//  Kokkos::deep_copy(d_pmid,h_pmid);
//  Kokkos::deep_copy(d_pdel,h_pdel);
//  Kokkos::deep_copy(d_zi,h_zi);
//
//  // If we are in charge of init-ing FQ as well, init it to 0.
//  if (m_fields.count("FQ")==1) {
//    // Init FQ to 0
//    auto d_FQ = m_fields.at("FQ").get_view();
//    Kokkos::deep_copy(d_FQ,Real(0));
//  }
//}
//  for (int i = 0; i < zm_inputs.size(); i++){
//    FieldInput field;
//    field.name = zm_inputs[i];
//    //Get and store device view using input name
//    Kokkos::View<Real*, Kokkos::LayoutRight, HostDevice> d_v = m_fields.at(field.name).get_view();
// //   field.device_view = d_v;
//
//    //Create and store host mirrors using device views
//    Kokkos::View<scream::Real*, Kokkos::LayoutRight, HostDevice> h_m = Kokkos::create_mirror_view(d_v);
// //   field.mirror = h_m;
//
//    //Create and store host mirrors raw pointers
//    Real* r_p = field.mirror.data();
// //   field.raw_ptr = r_p;
//    Kokkos::deep_copy(d_v, h_m);
//    map[zm_inputs[i]] = field;
//  }
//}


} // namespace scream
