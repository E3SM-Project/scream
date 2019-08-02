#ifndef SCREAM_IO_MANAGER_HPP
#define SCREAM_IO_MANAGER_HPP

#include "share/field/field.hpp"
#include "share/field/field_repository.hpp"
#include "share/grid/abstract_grid.hpp"

#include "share/util/time_stamp.hpp"
#include "share/util/factory.hpp"
#include "share/mpi/scream_comm.hpp"
#include "share/parameter_list.hpp"

namespace scream {

/*
 * A class responsible to handle I/O operations from/to file
 *
 * This class is the interface that the driver uses to read/write fields
 * from/to file. Here are a few notes on the interface and its indended behavior:
 *
 *  - The fields need to be "registered" in the IOManager during a registration
 *    phase (delimited by the 'registration_begins' and 'registration_ends' calls).
 *    During the registration phase, the user can pass a number of field repositories.
 *    The IOManager will grab all fields belonging to the 'Input', 'Output', or
 *    'Restart' groups. Additionally, it will also grab output/input files that
 *    are requested from the parameter list (if any of these fields were not yet
 *    marked as input/output (as it is most likely the case), the IOManager will
 *    proceed to add the fields to those groups).
 *  - The IOManager has three sets of fields: input, output and restart. 'Restart'
 *    fields are those needed to restart a simulation. Atmosphere processes are
 *    responsible to mark fields as 'Restart', by adding them to the 'Restart' group.
 *    'Output' fields are fields not necessarily needed for restart, but that
 *    the user is interested in saving (perhaps for post-processing purposes).
 *    Finally, 'input' fields are those that the user wants to load but that are
 *    not necessarily needed for restart (e.g., a 'measurement' field, that we
 *    want to compare the results against).
 *  - The file stream used for restart is *guaranteed* to be different from the
 *    stream used for output. If the user accidentally requests a file name for
 *    the output that coincides with one of the restart file name, the IOManager
 *    will error out. The restart file name has the time stamp in the file name,
 *    so that each restart file is guaranteed to not class with other restart files.
 */

template<typename ScalarType, typename DeviceType>
class IOManager {
public:
  using scalar_type     = ScalarType;
  using device_type     = DeviceType;
  using field_type      = Field<scalar_type,device_type>;
  using field_repo_type = FieldRepository<scalar_type,device_type>;
  using grid_type       = AbstractGrid;

  virtual ~IOManager () = default;

  // A string identifying the backend library used for I/O operations
  virtual std::string backend_name () const = 0;

  // The state of the io repo. To be used to check if registration_[begins|ends] have been called.
  virtual RepoState get_io_repo_state () const = 0;

  // Set/get the grid on which i/o has to be performed
  virtual void set_io_grid (const std::shared_ptr<grid_type>& grid) = 0;
  virtual const std::shared_ptr<grid_type>& get_io_grid () const = 0;

  // Begin/ends registration phase
  void registration_begins ();
  void registration_ends   ();

  // Scans the field repo for the fields that needs I/O (either b/c of their group or b/c requested from file)
  void register_fields (const field_repo_type& field_repo);

  // Get the stamp of next restart/output step after (or equal to) the given time
  virtual util::TimeStamp get_next_restart_step (const util::TimeStamp& current_time) const = 0;
  virtual util::TimeStamp get_next_output_step (const util::TimeStamp& current_time) const = 0;

  // Read/write from/to file all restart/input/output fields at given time step.
  // NOTE: if must_match=false, derived classes are free to choose their strategy
  //       if the exact time was not found in the input file (interpolate, closest, ...).
  // NOTE: the IOManager shold *not* error out if there are no input/output fields.
  //       Instead, it should simply return.
  virtual void read_restart_fields  (const util::TimeStamp& time) = 0;
  virtual void write_restart_fields (const util::TimeStamp& time) = 0;
  virtual void read_input_fields    (const util::TimeStamp& time, const bool must_match) = 0;
  virtual void write_output_fields  (const util::TimeStamp& time) = 0;

  // The names of the restart/input/output fields
  virtual const std::vector<std::string>& get_restart_fields_names () const = 0;
  virtual const std::vector<std::string>& get_input_fields_names  () const = 0;
  virtual const std::vector<std::string>& get_output_fields_names () const = 0;

  virtual bool has_restart_data () const = 0;

protected:

  virtual void do_registration_begins () = 0;
  virtual void do_registration_ends   () = 0;

  // Register a field as input (to be read from file)
  virtual void register_input_field  (const field_type& f) = 0;
  // Register a field as output (to be written to file)
  virtual void register_output_field (const field_type& f) = 0;
  // Register a field as both input and output (combination of the two above)
  virtual void register_io_field     (const field_type& f) = 0;
};

// A short name for the factory for atmosphere processes
template<typename ScalarType, typename DeviceType>
using IOFactory =
    util::Factory<IOManager<ScalarType,DeviceType>,
                  util::CaseInsensitiveString,
                  std::shared_ptr<IOManager<ScalarType,DeviceType>>,
                  const Comm&,const ParameterList&>;

// ======================= IMPLEMENTATION ======================= //

template<typename ScalarType, typename DeviceType>
void IOManager<ScalarType,DeviceType>::registration_begins ()
{
  scream_require_msg (static_cast<bool>(this->get_io_grid()),
                      "Error! You must set the I/O grid *before* registering fields in the IOManager.\n");

  this->do_registration_begins();
}

template<typename ScalarType, typename DeviceType>
void IOManager<ScalarType,DeviceType>::
register_fields(const field_repo_type& field_repo) {

  // Get restart/input/output group names from the repo, to avoid mispelling
  const auto& restart_group_name = field_repo.get_restart_group_name();
  const auto& input_group_name   = field_repo.get_input_group_name();
  const auto& output_group_name  = field_repo.get_output_group_name();

  // Before registering anything, for each field that was requested as
  // input/output from the paremeter list, we add "Input" and "Output"
  // to the field's groups.
  for (const auto& name : this->get_input_fields_names()) {
    const auto& aliases = field_repo.get_field_aliases(name);
    // Look for the field defined on the reference grid
    for (const auto& it : aliases) {
      if (it.first.get_grid_name()==get_io_grid()->name()) {
        it.second.get_header_ptr()->get_tracking().add_to_group(input_group_name);
        break;
      }
    }
  }
  for (const auto& name : this->get_output_fields_names()) {
    const auto& aliases = field_repo.get_field_aliases(name);
    // Look for the field defined on the reference grid
    for (const auto& it : aliases) {
      if (it.first.get_grid_name()==get_io_grid()->name()) {
        it.second.get_header_ptr()->get_tracking().add_to_group(output_group_name);
        break;
      }
    }
  }

  // Register all 'Restart' fields as input-output fields in the IOManager
  for (const auto& name : field_repo.get_field_groups().at(restart_group_name)) {
    const auto& aliases = field_repo.get_field_aliases(name);
    // Look for the field defined on the reference grid
    for (const auto& it : aliases) {
      if (it.first.get_grid_name()==this->get_io_grid()->name()) {
        this->register_io_field(it.second);
        break;
      }
    }
  }

  // Register all 'Output' fields as output fields in the IOManager
  for (const auto& name : field_repo.get_field_groups().at(output_group_name)) {
    const auto& aliases = field_repo.get_field_aliases(name);
    // Look for the field defined on the reference grid
    for (const auto& it : aliases) {
      if (it.first.get_grid_name()==get_io_grid()->name()) {
        this->register_output_field(it.second);
        break;
      }
    }
  }

  // Register all 'Input' fields as input fields in the IOManager
  for (const auto& name : field_repo.get_field_groups().at(input_group_name)) {
    const auto& aliases = field_repo.get_field_aliases(name);
    // Look for the field defined on the reference grid
    for (const auto& it : aliases) {
      if (it.first.get_grid_name()==get_io_grid()->name()) {
        this->register_input_field(it.second);
        break;
      }
    }
  }
}

template<typename ScalarType, typename DeviceType>
void IOManager<ScalarType,DeviceType>::registration_ends ()
{
  // TODO: add any type of check here

  this->do_registration_ends();
}


} // namespace scream

#endif // SCREAM_IO_MANAGER_HPP
