#ifndef SCREAM_PIO2_IO_MANAGER_HPP
#define SCREAM_PIO2_IO_MANAGER_HPP

#include "control/io/scream_io_manager.hpp"
#include <pio.h>

namespace scream {

template<typename ScalarType, typename DeviceType>
class PIO2_IOManager : public IOManager<ScalarType,DeviceType> {
public:
  using base_type   = IOManager<ScalarType,DeviceType>;
  using scalar_type = typename base_type::scalar_type;
  using device_type = typename base_type::device_type;
  using field_type  = typename base_type::field_type;

  // Constructor(s) & Destructor
  PIO2_IOManager (const Comm& comm, const ParameterList& params);
  ~PIO2_IOManager ();

  // A string identifying the backend library used for I/O operations 
  std::string backend_name () const override { return "PIO2"; }

  // The state of the io repo. To be used to check if registration_[begins|ends] have been called.
  RepoState get_io_repo_state () const override { return m_repo_state; }

  // Get the stamp of next restart/output step after (or equal to) the given time
  util::TimeStamp get_next_restart_step (const util::TimeStamp& current_time) const override;
  util::TimeStamp get_next_output_step (const util::TimeStamp& current_time) const override;

  // Read/write from/to file all restart/input/output fields at given time step.
  // NOTE: derived classes are free to choose their strategy if the exact time was not
  //       found in the input file (interpolate, closest, ...). Either way, they must
  //       return whether the requested time stamp was found in the file.
  //       Restart files are required to match the exact time stamp.
  void read_restart_fields  (const util::TimeStamp& time) override;
  void write_restart_fields (const util::TimeStamp& time) override;
  bool read_input_fields    (const util::TimeStamp& time, const bool must_match) override;
  void write_output_fields  (const util::TimeStamp& time) override;

  const std::vector<std::string>& get_output_fields_names  () const = 0;
  const std::vector<std::string>& get_input_fields_names   () const = 0;
  const std::vector<std::string>& get_restart_fields_names () const = 0;

  bool has_restart_data () const { return m_has_restart; }

private:

  // Begin/ends registration phase
  void do_registration_begins () override ;
  void do_registration_ends   () override ;

  // Register a field as input (to be read from file)
  void register_input_field  (const field_type& f) override;
  // Register a field as output (to be written to file)
  void register_output_field (const field_type& f) override;
  // Register a field as both input and output (combination of the two above)
  void register_io_field     (const field_type& f) override;

  Comm              m_atm_comm;

  RepoState         m_repo_state;

  ParameterList     m_io_params;

  bool              m_has_restart;

  // Use pimpl idiom
  struct PIO2_Data;

  std::unique_ptr<PIO2_Data> m_pio2_struct;
};

// A creator function for the PIO2_IOManager class
template<typename ScalarType, typename DeviceType>
inline std::shared_ptr<IOManager<ScalarType,DeviceType>>
create_pio2_io_manager (const Comm& comm, const ParameterList& params) {
  return std::make_shared<PIO2_IOManager<ScalarType,DeviceType>>(comm,params);
}

} // namespace scream

#endif // SCREAM_PIO2_IO_MANAGER_HPP
