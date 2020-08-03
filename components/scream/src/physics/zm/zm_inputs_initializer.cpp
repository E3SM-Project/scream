#include "physics/zm/zm_inputs_initializer.hpp"
#include "physics/zm/scream_zm_interface.hpp"

#include <array>
#include <string>
#include <iostream>
#include <unordered_map>

namespace scream
{

using namespace std;


std::vector<std::string> zm_inputs = {"t", "qh," "prec", "jctop", "jcbot", "pblh", "zm",
				"geos", "zi", "qtnd", "heat", "pap", "paph", "dpp", "delt", 
				"mcon", "cme", "cape", "tpert", "dlf", "pflx", "zdu", "rprd",
				"mu", "md", "du", "eu", "ed", "dp", "dsubcld", "jt", "maxg",
				"ideep", "lengath", "ql", "rliq", "landfrac", "hu_nm1",
      				"cnv_nm1", "tm1", "qm1", "t_star", "q_star", "dcape", "q",
				"tend_s", "tend_q", "cld", "snow", "ntprprd", "ntsnprd",
				"flxprec", "flxsnow", "ztodt", "pguall", "pgdall", "icwu",
				 "ncnst", "fracis"};

//template<typename T>


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

//unordered_map<std::string, FieldInput map;

void ZMInputsInitializer :: initialize_fields(){
  std :: cout << "here i am!\n";
  int count = 0;
  for (int i = 0 ; i < zm_inputs.size(); i++ ){
    count += m_fields.count(zm_inputs[i]);
  }
  std :: cout << "count is " + count;
  vector<FieldInput> fields;
    FieldInput t;
//  FieldInput qh;
//  FieldInput prec;
//  FieldInput jctop;
//  FieldInput jcbot;
//  FieldInput pblh;
    fields.push_back(t);
//  fields.push_back(qh);
//  fields.push_back(prec);
//  fields.push_back(jctop);
//  fields.push_back(jcbot);
//  fields.push_back(pblh);
  for (int i = 0; i < fields.size(); i++){
    //Get and store device view using input name
    Kokkos::View<Real*, Kokkos::LayoutRight, HostDevice> d_v = m_fields.at(fields[i].name).get_view();
    fields[i].device_view = d_v;
    //Create and store host mirrors using device views
    Kokkos::View<scream::Real*, Kokkos::LayoutRight, HostDevice> m = Kokkos::create_mirror_view(d_v);
    fields[i].mirror = m;
    //Create and store host mirrors raw pointers
    Real* r_p = fields[i].mirror.data();
    fields[i].raw_ptr = r_p;
    //Kokkos::deep_copy(d_v, r_p);
}


}

vector<FieldInput> create_fields(){
  vector<FieldInput> fields;
  FieldInput t;
  FieldInput qh;
  FieldInput prec;
  FieldInput jctop;
  FieldInput jcbot;
  FieldInput pblh;
  fields.push_back(t);
  fields.push_back(qh);
  fields.push_back(prec);
  fields.push_back(jctop);
  fields.push_back(jcbot);
  fields.push_back(pblh);
  return fields;
}

//void do_everything(){
//  vector<FieldInput> fields = create_fields();
//  for (int i = 0; i < fields.size(); i++){
//    auto d_v = m_fields.at(fields[i].name).get_view();
//    auto p = fields[i].
//    auto m = Kokkos::create_mirror_view(fields[i].device_view);
//    fields[i].name = zm_inputs[i];
//    fields[i].device_view = d_v;
//    fields[i].mirror = m;
//    fields[i].raw_ptr = i;
//    Kokkos::deep_copy(fields[i].device_view, fields[i].raw_ptr);
//  return fields; 
//}
//
//void assignNames(vector<FieldInput> fields){
//  for (int i = 0; i < fields.size(); i++){
//    fields[i].name = zm_inputs[i];
//  }
//}
//
//void get_device_views(vector<FieldInput> fields){
//  for (int i = 0; i < fields.size(); i++){
//    auto d_v = m_fields.at(fields[i].name.get_view());
//    fields[i].device_view = d_v;
//  }
//}
//
//void make_mirror(vector<FieldInput> fields){
//  for (int i = 0; i < fields.size(); i++){
//    auto m = Kokkos::create_mirror_view(fields[i].device_view);
//    fields[i].mirror = m;
//  }
//}
//
//// Get host mirros' raw pointers
//void get_ptrs(vector<FieldInput> fields){
//  for (int i = 0; i < fields.size(); i++){
//    auto p = fields[i].
//    fields[i].raw_ptr = i;
//  }
//}
//
//// Deep copy back to device
//Kokkos::deep_copy(d_q,h_q);
//void copy_to_device(vector<FieldInput> fields){
//  for (int i = 0; i < fields.size(); i++){
//    Kokkos::deep_copy(fields[i].device_view, fields[i].raw_ptr);
//  }
//}
//



//Make structs and put them in an array

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

} // namespace scream
