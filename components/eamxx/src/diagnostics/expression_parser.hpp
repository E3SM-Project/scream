#ifndef EAMXX_EXPRESSION_PARSER_DIAGNOSTIC_HPP
#define EAMXX_EXPRESSION_PARSER_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class ExprParser : public AtmosphereDiagnostic
{
public:
  using Pack          = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using PF            = scream::PhysicsFunctions<DefaultDevice>;
  using KT            = KokkosTypes<DefaultDevice>;
  using MemberType    = typename KT::MemberType;
  using view_1d       = typename KT::template view_1d<Pack>;

  // Constructors
  ExprParser (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic
  std::string name () const override { return m_expr; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> /* grids_manager */) override {}

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl () override;
protected:

  void validate_expression () const;

  std::string m_expr;

}; // class IceWaterPathDiagnostic

} //namespace scream

#endif // EAMXX_EXPRESSION_PARSER_DIAGNOSTIC_HPP
