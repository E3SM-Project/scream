#ifndef SCREAM_OUTPUT_MANAGER_HPP
#define SCREAM_OUTPUT_MANAGER_HPP

#include "share/io/scorpio_output.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/io/scream_io_utils.hpp"

#include "share/field/field_manager.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

namespace scream
{

/*
 * The OutputManager Class handles all individual output streams.
 * Rather than have the SCREAM-AD or unit tests call the scorpio_output
 * class objects directly, the OutputManager stores these in a vector,
 * m_output_streams.
 *
 * Similar to a typical atmospheric process object, the OutputManager
 * has an init, run and finalize routine which is called during the AD
 * during those respective steps.
 *
 * PROPER USAGE:
 * The output manager requires a communication group, a parameter list of control
 * variables, a grid manager and a field manager.
 * Each of these four things are set using the setter function 'set_X' where
 *   X = comm, for the EKAT comm group.
 *     = params, for the parameter list.  In typical SCREAM runs this parameter
 *       list is a sublist of the scream control yaml file.
 *     = grids, for the grids mananger
 *     = fm, for the field manager.
 * The setup of the output manager in the SCREAM-AD is one of the last steps to
 * ensure that all four of the above objects have already been constructed.
 * see /control/atmospheric_driver.cpp for an example.
 *
 * For UNIT TESTS:
 * The output manager does require a comm group, list of parameters, grids manager
 * and field manager to work.  If output is desired in a unit test then these
 * must be established.  There are examples in /src/share/io/tests of how to
 * establish a simple grids manager and field manager.  As well as how to
 * locally create a parameter list.
 *
 * Adding output streams mid-simulation:
 * It is possible to add an output stream after init has been called by calling
 * the internal function 'add_output_stream' which takes an EKAT parameter list as input.
 * See comments in add_output_stream below for more details.
 *
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 *  (2021-08-19) Luca Bertagna (SNL)
 */
class OutputManager
{
public:
  using fm_type = FieldManager<Real>;
  using gm_type = GridsManager;

  // Constructor(s) & Destructor
  OutputManager () = default;
  virtual ~OutputManager () = default;

  // Set up the manager, creating all output streams. Inputs:
  //  - params: the parameter list with file/fields info, as well as method of output options
  //  - model_restart_output: whether this output stream is to write a model restart file
  void setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
              const std::shared_ptr<fm_type>& field_mgr,
              const std::shared_ptr<const gm_type>& grids_mgr,
              const util::TimeStamp& t0,
              const bool is_model_restart_output,
              const bool is_restarted_run);

  void setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
              const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs,
              const std::shared_ptr<const gm_type>& grids_mgr,
              const util::TimeStamp& t0,
              const bool is_model_restart_output,
              const bool is_restarted_run);
  void run(util::TimeStamp& current_ts);
  void finalize();

protected:

  std::string compute_filename_root (const IOControl& control, const IOFileSpecs& file_specs) const;

  // Craft the restart parameter list
  void set_params (const ekat::ParameterList& params,
                   const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs);

  using output_type     = AtmosphereOutput;
  using output_ptr_type = std::shared_ptr<output_type>;

  std::vector<output_ptr_type>   m_output_streams;
  ekat::Comm                     m_io_comm;
  ekat::ParameterList            m_params;

  // The output filename root
  std::string       m_casename;

  // How to combine multiple snapshots in the output: Instant, Max, Min, Average
  std::string       m_avg_type;

  // Whether this OutputManager handles a model restart file, or normal model output.
  bool m_is_model_restart_output;

  // Frequency of output and checkpointing
  // See scream_io_utils.hpp for details.
  IOControl m_output_control;
  IOControl m_checkpoint_control;

  // The specs (name, capacity, size) of output and checkpoint file.
  // See scream_io_utils.hpp for details.
  IOFileSpecs m_output_file_specs;
  IOFileSpecs m_checkpoint_file_specs;

  // Whether this run is the restart of a previous run, in which case
  // we might have to load an output checkpoint file (depending on avg type)
  bool m_is_restarted_run;

  // The simulation start date/time. We use this to produce a 'time'
  // var in the output file, corresponding to seconds_since_start_of_simulation.
  util::TimeStamp   m_t0;
};

} // namespace scream

#endif // SCREAM_OUTPUT_MANAGER_HPP
