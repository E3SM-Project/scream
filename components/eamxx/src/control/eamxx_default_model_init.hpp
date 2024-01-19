#ifndef EAMXX_DEFAULT_MODEL_INIT_HPP
#define EAMXX_DEFAULT_MODEL_INIT_HPP

#include "share/eamxx_model_init.hpp"

namespace scream
{

// Concrete implementation of ModelInit where we simply read in
// the eamxx input fields from input file, using the 'Physics'
// grid from the grids manager.
class DefaultModelInit : public ModelInit {
public:

  DefaultModelInit (const strmap_t<Field>& eamxx_inputs,
                    const std::shared_ptr<const GridsManager>& gm,
                    const ekat::ParameterList& ic_pl,
                    const util::TimeStamp& run_t0);

  ~DefaultModelInit () = default;

  void set_initial_conditions (const std::shared_ptr<ekat::logger::LoggerBase>& logger);

protected:
  void read_ic_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger) override;
  void read_topo_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger) override;
};

inline std::shared_ptr<ModelInit>
create_default_model_init (const std::map<std::string,Field>& eamxx_inputs,
                           const std::shared_ptr<const GridsManager>& gm,
                           const ekat::ParameterList& params,
                           const util::TimeStamp& run_t0)
{
  return std::make_shared<DefaultModelInit>(eamxx_inputs,gm,params,run_t0);
}

inline void register_default_model_init () {
  ModelInitFactory::instance().register_product("default",&create_default_model_init);
}

} // namespace scream

#endif // EAMXX_DEFAULT_MODEL_INIT_HPP
