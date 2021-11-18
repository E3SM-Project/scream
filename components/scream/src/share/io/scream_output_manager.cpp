#include "scream_output_manager.hpp"

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <fstream>
#include <memory>

namespace scream
{

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::shared_ptr<fm_type>& field_mgr,
       const std::shared_ptr<const gm_type>& grids_mgr,
       const util::TimeStamp& t0,
       const bool is_model_restart_output,
       const bool is_restarted_run)
{
  using map_t = std::map<std::string,std::shared_ptr<fm_type>>;
  map_t fms;
  fms[field_mgr->get_grid()->name()] = field_mgr;
  setup(io_comm,params,fms,grids_mgr,t0,is_model_restart_output,is_restarted_run);
}

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs,
       const std::shared_ptr<const gm_type>& grids_mgr,
       const util::TimeStamp& t0,
       const bool is_model_restart_output,
       const bool is_restarted_run)
{
  m_io_comm = io_comm;
  m_t0      = t0;
  m_is_restarted_run = is_restarted_run;
  m_is_model_restart_output = is_model_restart_output;

  // Check for model restart output
  set_params(params,field_mgrs);

  // Output control
  auto& out_control_pl = m_params.sublist("Output Control");
  m_output_control.frequency  = out_control_pl.get<int>("Frequency");
  m_output_control.frequency_units = out_control_pl.get<std::string>("Frequency Units");
  m_output_control.nsteps_since_last_write = 0;

  // File specs
  m_output_file_specs.max_snapshots_in_file = m_params.get<int>("Max Snapshots Per File");
  m_output_file_specs.num_snapshots_in_file = 0;
  m_output_file_specs.filename_with_time_string = out_control_pl.get("Timestamp in Filename",true);
  m_output_file_specs.filename_with_mpiranks    = out_control_pl.get("MPI Ranks in Filename",false);

  // For each grid, create a separate output stream.
  if (field_mgrs.size()==1) {
    auto output = std::make_shared<output_type>(m_io_comm,m_params,field_mgrs.begin()->second,grids_mgr);
    m_output_streams.push_back(output);
  } else {
    const auto& fields_pl = m_params.sublist("Fields");
    for (auto it=fields_pl.sublists_names_cbegin(); it!=fields_pl.sublists_names_cend(); ++it) {
      const auto& gname = *it;
      EKAT_REQUIRE_MSG (grids_mgr->has_grid(gname),
          "Error! Output requested on grid '" + gname + "', but the grids manager does not store such grid.\n");

      EKAT_REQUIRE_MSG (field_mgrs.find(gname)!=field_mgrs.end(),
          "Error! Output requested on grid '" + gname + "', but no field manager is available for such grid.\n");

      auto output = std::make_shared<output_type>(m_io_comm,m_params,field_mgrs.at(gname),grids_mgr);
      m_output_streams.push_back(output);
    }
  }

  const auto has_restart_data = (m_avg_type!=OutputAvgType::Instant || m_output_control.frequency>1);
  if (has_restart_data) {
    if (m_params.isSublist("Checkpoint Control")) {
      // Output control
      auto& pl = m_params.sublist("Checkpoint Control");
      m_checkpoint_control.frequency  = pl.get<int>("Frequency");
      m_checkpoint_control.frequency_units = pl.get<std::string>("Frequency Units");
      m_checkpoint_control.nsteps_since_last_write = 0;

      // File specs
      m_checkpoint_file_specs.max_snapshots_in_file = 1;
      m_checkpoint_file_specs.num_snapshots_in_file = 0;
      m_checkpoint_file_specs.filename_with_time_string = pl.get("Timestamp in Filename",true);
      m_checkpoint_file_specs.filename_with_mpiranks    = pl.get("MPI Ranks in Filename",false);
    }
  }

  // A restarted run needs at the very least to figure out the start date/time of
  // the original simulation.
  // On top of that, if this is normal output (not the model restart output),
  // and if the output specs require it, we might have to restart the output history.
  // E.g., we might save 30-day avg value for field F, but due to job size
  // break the run into three 10-day runs. We then need to save the state of
  // our averaging in a "restart" file (e.g., the current avg and avg_count).
  // Note: the user might decide *not* to restart the output, so give the option
  //       of disabling the restart. Also, the user might want to change the
  //       casename, so allow to specify a different casename for the restart file.
  if (m_is_restarted_run) {

    // Allow to skip history restart, or to specify a casename for the restart file
    // that is different from the casename of the current output.
    auto& restart_pl = m_params.sublist("Restart");
    bool perform_history_restart = has_restart_data && restart_pl.get("Perform Restart",true);
    auto hist_restart_casename = restart_pl.get("Casename",m_casename);

    if (m_is_model_restart_output) {
      // If this is the model restart output, We look for the original start date/time
      // in the model restart file. For "normal" output, we fully trust the t0 that was
      // passed to the setup method. The reasons are the following:
      //  - history restart files may not have been written, if the last model restart
      //    step coincided with an output step.
      //  - reading the model restart file would require this OutputManager to be provided
      //    with the casename of the model restart output file.
      // It is simpler to just have the AD handle the retrieval of the original simulation
      // start date/time, and pass it to the setup method. This can be easily achieved by
      // setting up the OM for the model restart *before* those of the normal output.

      // Find out model restart file name
      auto model_restart_filename = find_filename_in_rpointer(hist_restart_casename,".r.nc");

      // Recover the original simulation start date/time
      if (m_io_comm.am_i_root()) {
        std::cout << "Restarting simulation from '" << model_restart_filename << "' nc file.\n";
      }
      ekat::ParameterList res_params("Input Parameters");
      res_params.set<std::string>("Filename",model_restart_filename);
      AtmosphereInput model_restart (m_io_comm,res_params);
      int start_date = model_restart.read_int_scalar("start_date"); // YYYYMMDD
      int start_time = model_restart.read_int_scalar("start_time"); // HHMMSS
      model_restart.finalize();

      std::vector<int> date = {start_date/10000, (start_date/100) % 100, start_date % 100};
      std::vector<int> time = {start_time/10000, (start_time/100) % 100, start_time % 100};
      m_t0 = util::TimeStamp(date,time);
    } else {
      // Not a model restart output. Trust what the input t0 is.
      m_t0 = t0;
    }

    if (perform_history_restart) {
      auto output_restart_filename = find_filename_in_rpointer(hist_restart_casename,".rhist.nc");

      ekat::ParameterList res_params("Input Parameters");
      res_params.set<std::string>("Filename",output_restart_filename);
      AtmosphereInput output_restart (m_io_comm,res_params);
      // Also restart each stream
      for (auto stream : m_output_streams) {
        stream->restart(output_restart_filename);
      }
      // And restart the avg count.
      m_output_control.nsteps_since_last_write = output_restart.read_int_scalar("avg_count");
      output_restart.finalize();
    }
  }

  if (m_params.get("Save Initial State", false)) {
    this->run(t0);
  }
}
/*===============================================================================================*/
void OutputManager::run(const util::TimeStamp& timestamp)
{
  using namespace scorpio;

  // Check if we need to open a new file
  ++m_output_control.nsteps_since_last_write;
  ++m_checkpoint_control.nsteps_since_last_write;

  // Check if this is a write step (and what kind)
  const bool is_output_step     = m_output_control.is_write_step();
  const bool is_checkpoint_step = m_checkpoint_control.is_write_step() && not is_output_step;
  const bool is_write_step      = is_output_step || is_checkpoint_step;

  // If neither output or checkpoint, these won't be used anyways,
  // so no need to check if is_write_step == true.
  auto& control   = is_checkpoint_step ? m_checkpoint_control : m_output_control;
  auto& filespecs = is_checkpoint_step ? m_checkpoint_file_specs : m_output_file_specs;
  auto& filename  = filespecs.filename;

  std::cout << "timestamp: " << timestamp.to_string() << "\n"
            << "  output: " << (is_output_step ? "yes\n" : "no\n")
            << "  checkpoint: " << (is_checkpoint_step ? "yes\n" : "no\n")
            << "  casename: " << m_casename << "\n"
            << "  output freq: " << m_output_control.frequency << "\n"
            << "  checkpt freq: " << m_checkpoint_control.frequency << "\n"
            << "  steps since last output: " << m_output_control.nsteps_since_last_write << "\n"
            << "  steps since last checkpoint: " << m_checkpoint_control.nsteps_since_last_write << "\n";

  // Compute filename (if write step)
  if (is_write_step) {
    // Check if we need to open a new file
    if (not filespecs.is_open) {
      // Compute new file name
      filename = compute_filename_root(control,filespecs);
      if (filespecs.filename_with_time_string) {
        filename += "." + timestamp.to_string();
      }
      if (is_output_step) {
        std::cout << "It's output. Model restart: " << (m_is_model_restart_output ? "yes\n" : "no\n");
        filename += m_is_model_restart_output ? ".r.nc" : ".nc";
      } else if (is_checkpoint_step) {
        filename += ".rhist.nc";
      } else {
        filename += ".nc";
      }

      // Register new netCDF file for output. First, check no other output managers
      // are trying to write on the same file
      EKAT_REQUIRE_MSG (not is_file_open_c2f(filename.c_str(),Write),
          "Error! File '" + filename + "' is currently open for write. Cannot share with other output managers.\n");
      register_file(filename,Write);

      // Note: time has an unknown length. Setting its "length" to 0 tells the scorpio to
      // set this dimension as having an 'unlimited' length, thus allowing us to write
      // as many timesnaps to file as we desire.
      register_dimension(filename,"time","time",0);

      // Register time as a variable.
      register_variable(filename,"time","time",1,{"time"},  PIO_REAL,"time");

      // Make all output streams register their dims/vars
      for (auto& it : m_output_streams) {
        it->setup_output_file(filename);
      }

      // Set degree of freedom for "time"
      int time_dof[1] = {0};
      set_dof(filename,"time",0,time_dof);

      // Finish the definition phase for this file.
      eam_pio_enddef (filename); 
      if (is_checkpoint_step) { 
        set_int_attribute_c2f (filename.c_str(),"avg_count",m_output_control.nsteps_since_last_write);
      }
      auto t0_date = m_t0.get_date()[0]*10000 + m_t0.get_date()[1]*100 + m_t0.get_date()[2];
      auto t0_time = m_t0.get_time()[0]*10000 + m_t0.get_time()[1]*100 + m_t0.get_time()[2];
      set_int_attribute_c2f(filename.c_str(),"start_date",t0_date);
      set_int_attribute_c2f(filename.c_str(),"start_time",t0_time);
      filespecs.is_open = true;
    }

    // If we are going to write an output checkpoint file, or a model restart file,
    // we need to append to the filename ".rhist" or ".r" respectively, and add
    // the filename to the rpointer.atm file.
    if (m_is_model_restart_output || is_checkpoint_step) {
      if (m_io_comm.am_i_root()) {
        std::ofstream rpointer;
        rpointer.open("rpointer.atm",std::ofstream::app);  // Open rpointer file and append to it
        rpointer << filename << std::endl;
      }
    }

    // Update time and nsteps in the output file
    pio_update_time(filename,timestamp.seconds_from(m_t0));
    if (m_is_model_restart_output) {
      // Only write nsteps on model restart
      set_int_attribute_c2f(filename.c_str(),"nsteps",timestamp.get_num_steps());
    }
  }

  // Run the output streams
  std::cout << "number of streams: " << m_output_streams.size() << "\n";
  for (auto& it : m_output_streams) {
    // Note: filename might referencing an invalid string, but it's only used
    //       in case is_write_step=true, in which case it will *for sure* contain
    //       a valid file name.
    printf("running stream, is write: %s\n",is_write_step ? "yes" : "no");
    printf("time: %s\n", timestamp.to_string().c_str());
    printf("filename: %s\n", filename.c_str());
    it->run(filename,is_write_step,m_output_control.nsteps_since_last_write);
  }

  if (is_write_step) {
    // We're adding one snapshot to the file
    ++filespecs.num_snapshots_in_file;

    // Finish up any updates to output file
    sync_outfile(filename);

    // Check if we need to close the output file
    if (filespecs.file_is_full()) {
      eam_pio_closefile(filename);
      filespecs.num_snapshots_in_file = 0;
      filespecs.is_open = false;
      control.nsteps_since_last_write = 0;
    }

    // Whether we wrote an output or a checkpoint, the checkpoint counter needs to be reset
    m_checkpoint_control.nsteps_since_last_write = 0;
  }
}
/*===============================================================================================*/
void OutputManager::finalize()
{
  // Swapping with an empty mgr is the easiest way to cleanup.
  OutputManager other;
  std::swap(*this,other);
}

std::string OutputManager::
compute_filename_root (const IOControl& control, const IOFileSpecs& file_specs) const
{
  return m_casename + "." +
         e2str(m_avg_type) + "." +
         control.frequency_units+ "_x" +
         std::to_string(control.frequency) +
         (file_specs.filename_with_mpiranks ? ".np" + std::to_string(m_io_comm.size()) : "");
}

std::string OutputManager::
find_filename_in_rpointer (const std::string& casename, const std::string& suffix) const
{
  std::string filename;
  bool found = false;
  std::string content,line;
  if (m_io_comm.am_i_root()) {
    std::ifstream rpointer_file;
    rpointer_file.open("rpointer.atm");
    // Note: keep swallowing line, even after the first match, since we never wipe
    //       rpointer.atm, so it might contain multiple matches, and we want to pick
    //       the last one (which is the last restart file that was written).
    while (rpointer_file >> line) {
      content += line + "\n";
      if (line.find(casename) != std::string::npos &&
          line.find(suffix) != std::string::npos) {
        found = true;
        filename = line;
      }
    }
  }
  m_io_comm.broadcast(&found,1,0);

  // If the history restart file is not found, it must be because the last
  // model restart step coincided with a model output step, in which case
  // a restart history file is not written.
  // If that's the case, *disable* output restart, by setting
  //   'Restart'->'Perform Restart' = false
  // in the input parameter list
  EKAT_REQUIRE_MSG (found,
      "Error! Output restart requested, but the no history restart file found in 'rpointer.atm'.\n"
      "   restart file name root: " + casename + "\n"
      "   rpointer content:\n" + content);

  // Have the root rank communicate the nc filename
  broadcast_string(filename,m_io_comm,m_io_comm.root_rank());

  return filename;
}

void OutputManager::
set_params (const ekat::ParameterList& params,
            const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs)
{
  m_params = params;
  if (m_is_model_restart_output) {
    using vos_t = std::vector<std::string>;

    // We build some restart parameters internally
    auto avg_type = m_params.get<std::string>("Averaging Type","INSTANT");
    m_avg_type = str2avg(avg_type);
    EKAT_REQUIRE_MSG (m_avg_type==OutputAvgType::Instant,
        "Error! For restart output, the averaging type must be 'Instant'.\n"
        "   Note: you don't have to specify this parameter for restart output.\n");
    m_output_file_specs.max_snapshots_in_file = m_params.get("Max Snapshots Per File",1);
    EKAT_REQUIRE_MSG (m_output_file_specs.max_snapshots_in_file==1,
        "Error! For restart output, max snapshots per file must be 1.\n"
        "   Note: you don't have to specify this parameter for restart output.\n");

    auto& fields_pl = m_params.sublist("Fields");
    for (const auto& it : field_mgrs) {
      const auto& fm = it.second;
      vos_t fnames;
      // There may be no RESTART group on this grid
      if (fm->has_group("RESTART")) {
        auto restart_group = fm->get_groups_info().at("RESTART");
        EKAT_REQUIRE_MSG (not fields_pl.isParameter(it.first),
          "Error! For restart output, don't specify the fields names. We will create this info internally.\n");
        for (const auto& n : restart_group->m_fields_names) {
          fnames.push_back(n);
        }
      }
      fields_pl.sublist(it.first).set("Fields Names",fnames);
    }
    m_casename = m_params.get<std::string>("Casename", "scream_restart");
  } else {
    auto avg_type = m_params.get<std::string>("Averaging Type");
    m_avg_type = str2avg(avg_type);
    EKAT_REQUIRE_MSG (m_avg_type!=OutputAvgType::Invalid,
        "Error! Unsupported averaging type '" + avg_type + "'.\n"
        "       Valid options: Instant, Max, Min, Average. Case insensitive.\n");

    m_output_file_specs.max_snapshots_in_file = m_params.get<int>("Max Snapshots Per File");
    m_casename = m_params.get<std::string>("Casename");
  }
}
/*===============================================================================================*/
} // namespace scream
