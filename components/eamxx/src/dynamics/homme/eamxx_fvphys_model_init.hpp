#ifndef EAMXX_FVPHYS_MODEL_INIT_HPP
#define EAMXX_FVPHYS_MODEL_INIT_HPP

#include "share/eamxx_model_init.hpp"

namespace scream
{

class FvPhysModelInit : public ModelInit
{
public:
  FvPhysModelInit (const strmap_t<Field>& eamxx_inputs,
                   const std::shared_ptr<const GridsManager>& gm,
                   const ekat::ParameterList& ic_pl,
                   const util::TimeStamp& run_t0);

  ~FvPhysModelInit () = default;

  void set_initial_conditions (const std::shared_ptr<ekat::logger::LoggerBase>& logger);

protected:

  void read_ic_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger) override;
  void read_topo_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger) override;

  void remap_fields ();

  strmap_t<Field>   m_gll_fields;
};

inline std::shared_ptr<ModelInit>
create_fvphys_model_init (const std::map<std::string,Field>& eamxx_inputs,
                          const std::shared_ptr<const GridsManager>& gm,
                          const ekat::ParameterList& params,
                          const util::TimeStamp& run_t0)
{
  return std::make_shared<FvPhysModelInit>(eamxx_inputs,gm,params,run_t0);
}

} // namespace scream

#endif // EAMXX_FVPHYS_MODEL_INIT_HPP
