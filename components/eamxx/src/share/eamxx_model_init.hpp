#ifndef EAMXX_MODEL_INIT_HPP
#define EAMXX_MODEL_INIT_HPP

#include "share/field/field.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"

#include <ekat/logging/ekat_logger.hpp>
#include <ekat/util/ekat_factory.hpp>

#include <vector>

namespace scream
{

// Interface class for initializing fields at model init time
// Derived classes can implement whatever strategy is needed
// For instance, they could simply look for them in the IC file
// An alternative (which is the case for PG2) is to read a version
// of the fields on a different grid, and then remap to the actual
// fields that were requested

class ModelInit {
public:
  using strvec_t = std::vector<std::string>;
  template<typename T>
  using strmap_t = std::map<std::string,T>;

  ModelInit (const strmap_t<Field>& eamxx_inputs,
             const std::shared_ptr<const GridsManager>& gm,
             const ekat::ParameterList& ic_pl,
             const util::TimeStamp& run_t0);

  virtual ~ModelInit () = default;

  // Note: if all derived classes simply call set_constant_field->read_ic_file->read_topo_file,
  //       you may want to put that impl in the base class
  virtual void set_initial_conditions (const std::shared_ptr<ekat::logger::LoggerBase>& logger) = 0;

protected:

  // Separate inputs into topo fields and IC fields
  void separate_inputs ();

  void set_constant_fields (const std::shared_ptr<ekat::logger::LoggerBase>& logger);
  virtual void read_ic_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger) = 0;
  virtual void read_topo_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger) = 0;

  // ------------------------ Attributes ------------------------ //

  // Parameters
  ekat::ParameterList   m_params;

  // Grids manager
  std::shared_ptr<const GridsManager> m_gm;

  // Run start date/time
  util::TimeStamp     m_run_t0;

  // Map name->field for all eamxx input fields (passed at construction time)
  strmap_t<Field>   m_eamxx_inputs;

  // Subset of EAMxx input fields that are const-inited
  std::set<Field>   m_const_fields;

  // Subset of EAMxx input fields that are topography-related
  // NOTE: when set_constant_fields is called, we remove const fields from these lists
  strvec_t          m_topo_fields_names;
  std::set<Field>   m_topo_fields;

  // Subset of EAMxx input fields that are not topography-related
  // NOTE: when set_constant_fields is called, we remove const fields from these lists
  strvec_t          m_ic_fields_names;
  std::set<Field>   m_ic_fields;
};

// The driver will use this factory to agnostically build the concrete model init
using ModelInitFactory =
    ekat::Factory<ModelInit,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<ModelInit>,
                  const std::map<std::string,Field>&,
                  const std::shared_ptr<const GridsManager>&,
                  const ekat::ParameterList&,
                  const util::TimeStamp&>;

} // namespace scream

#endif // EAMXX_MODEL_INIT_HPP
