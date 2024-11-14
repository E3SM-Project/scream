#ifndef EAMXX_DATA_INTERPOLATION_HPP
#define EAMXX_DATA_INTERPOLATION_HPP

#include "share/grid/abstract_grid.hpp"

#include "share/util/scream_time_stamp.hpp"

#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

#include "share/io/scorpio_input.hpp"

namespace scream{

class DataInterpolation
{
public:
  // Constructor(s) & Destructor
  DataInterpolation (const std::shared_ptr<const AbstractGrid>& grid,
                     const ekat::ParameterList& params);
  ~DataInterpolation () = default;

  void run (const util::TimeStamp& ts);

  void set_field (const Field& f);

  void complete_setup (const util::TimeStamp& t0);

protected:

  using strvec_t = std::vector<std::string>;

  void update_end_fields ();

  struct TimeState {
    util::TimeStamp beg;
    util::TimeStamp end;
    int file_index;
  };

  enum Phase {
    AfterRead = 0,
    AfterVInterp,
    AfterHInterp
  };

  std::vector<Field>& get_fields (Phase phase, bool beg);

  TimeState m_time_state;

  std::map<Phase,std::vector<Field>> m_fields_beg;
  std::map<Phase,std::vector<Field>> m_fields_end;
  std::vector<Field>  m_tgt_fields;

  strvec_t m_input_files;

  ekat::ParameterList       m_params;
};

} // namespace scream

#endif // EAMXX_DATA_INTERPOLATION_HPP
