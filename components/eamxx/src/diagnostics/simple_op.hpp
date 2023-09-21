#ifndef EAMXX_SIMPLE_OP_DIAGNOSTIC_HPP
#define EAMXX_SIMPLE_OP_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will perform elementary unary/binary op on 1-2 input fields
 */

class SimpleOp : public AtmosphereDiagnostic
{
public:
  // Constructors
  SimpleOp (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic
  std::string name () const { return m_name; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> /* grids_manager */) {}

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:
  void initialize_impl (const RunType /*run_type*/);

  enum Op { ADD, SUB, PROD, DIV, NEG };
  enum F2Type { FIELD, VALUE};

  static std::string op_str(const Op op) {
    std::string s;
    switch (op) {
      case ADD:  s = "add";  break;
      case SUB:  s = "sub";  break;
      case PROD: s = "prod"; break;
      case DIV:  s = "div";  break;
      case NEG:  s = "neg";  break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported value for SimpleOp::Op enum.\n");
    }
    return s;
  }

  static std::string f2type_str(const F2Type t) {
    std::string s;
    switch (t) {
      case VALUE: s = "value";  break;
      case FIELD: s = "field";  break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported value for SimpleOp::F2Type enum.\n");
    }
    return s;
  }

  Op            m_op;
  std::string   m_name;

  F2Type        m_f2_type;
  Real          m_f2_value;
  Real          m_f2_scale;
  Real          m_f1_scale;
};

} //namespace scream

#endif // EAMXX_SIMPLE_OP_DIAGNOSTIC_HPP
