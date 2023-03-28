#include "scream_output_manager.hpp"

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/util/scream_timing.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#if __cplusplus > 201700L
#include <filesystem>
namespace filesystem_ns = std::filesystem;
#else
#include <experimental/filesystem>
namespace filesystem_ns = std::experimental::filesystem;
#endif

#include <fstream>
#include <memory>

namespace scream
{

// Local helper functions:
void set_file_header(const std::string& filename);

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::shared_ptr<fm_type>& field_mgr,
       const std::shared_ptr<const gm_type>& grids_mgr,
       const util::TimeStamp& run_t0,
       const util::TimeStamp& case_t0,
       const bool is_model_restart_output)
{
  using map_t = std::map<std::string,std::shared_ptr<fm_type>>;
  map_t fms;
  fms[field_mgr->get_grid()->name()] = field_mgr;
  setup(io_comm,params,fms,grids_mgr,run_t0,case_t0,is_model_restart_output);
}

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs,
       const std::shared_ptr<const gm_type>& grids_mgr,
       const util::TimeStamp& run_t0,
       const util::TimeStamp& case_t0,
       const bool is_model_restart_output)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (run_t0.is_valid(),
      "Error! Invalid case_t0 timestamp: " + case_t0.to_string() + "\n");
  EKAT_REQUIRE_MSG (run_t0.is_valid(),
      "Error! Invalid run_t0 timestamp: " + case_t0.to_string() + "\n");
  EKAT_REQUIRE_MSG (case_t0<=run_t0,
      "Error! The case_t0 timestamp must precede run_t0.\n"
      "   run_t0 : " + run_t0.to_string() + "\n"
      "   case_t0: " + case_t0.to_string() + "\n");

  m_io_comm = io_comm;
  m_run_t0 = run_t0;
  m_case_t0 = case_t0;
  m_is_restarted_run = (case_t0<run_t0);
  m_is_model_restart_output = is_model_restart_output;

  // Check for model restart output
  set_params(params,field_mgrs);

  // Output control
  EKAT_REQUIRE_MSG(m_params.isSublist("output_control"),
      "Error! The output control YAML file for " + m_casename + " is missing the sublist 'output_control'");
  auto& out_control_pl = m_params.sublist("output_control");
  // Determine which timestamp to use a reference for output frequency.  Two options:
  // 	1. use_case_as_start_reference: TRUE  - implies we want to calculate frequency from the beginning of the whole simulation, even if this is a restarted run.
  // 	2. use_case_as_start_reference: FALSE - implies we want to base the frequency of output on when this particular simulation started.
  // Note, (2) is needed for restarts since the restart frequency in CIME assumes a reference of when this run began.
  const bool start_ref = out_control_pl.get<bool>("use_case_as_start_reference",!m_is_model_restart_output);
  m_output_control.frequency_units = out_control_pl.get<std::string>("frequency_units");
  // In case output is disabled, no point in doing anything else
  if (m_output_control.frequency_units=="none" || m_output_control.frequency_units=="never") {
    m_output_disabled = true;
    return;
  }
  m_output_control.frequency = out_control_pl.get<int>("Frequency");
  EKAT_REQUIRE_MSG (m_output_control.frequency>0,
      "Error! Invalid frequency (" + std::to_string(m_output_control.frequency) + ") in Output Control. Please, use positive number.\n");

  m_output_control.timestamp_of_last_write = start_ref ? m_case_t0 : m_run_t0;

  // File specs
  m_output_file_specs.max_snapshots_in_file = m_params.get<int>("Max Snapshots Per File",-1);
  m_output_file_specs.filename_with_mpiranks    = out_control_pl.get("MPI Ranks in Filename",false);
  m_output_file_specs.save_grid_data            = out_control_pl.get("save_grid_data",!m_is_model_restart_output);

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

  // For normal output, setup the geometry data streams, which we used to write the
  // geo data in the output file when we create it.
  if (m_output_file_specs.save_grid_data) {
    std::set<std::shared_ptr<const AbstractGrid>> grids;
    for (auto& it : m_output_streams) {
      grids.insert(it->get_io_grid());
    }

    // If 2+ grids are present, we mandate suffix on all geo_data fields,
    // to avoid clashes of names.
    bool use_suffix = grids.size()>1;
    for (auto grid : grids) {
      std::vector<Field> fields;
      for (const auto& fn : grid->get_geometry_data_names()) {
        const auto& f = grid->get_geometry_data(fn);
        if (use_suffix) {
          fields.push_back(f.clone(f.name()+"_"+grid->m_short_name));
        } else {
          fields.push_back(f.clone());
        }
      }
      auto output = std::make_shared<output_type>(m_io_comm,fields,grid);
      m_geo_data_streams.push_back(output);
    }
  }

  // Check if this kind of output needs history restart files (in general)
  const auto has_restart_data = m_avg_type!=OutputAvgType::Instant
    && (m_output_control.frequency_units!="nsteps" || m_output_control.frequency>1);

  if (has_restart_data && m_params.isSublist("Checkpoint Control")) {
    // Output control
    // TODO: It would be great if there was an option where, if Checkpoint Control was not a sublist, we
    //       could query the restart control information and just use that. 
    auto& pl = m_params.sublist("Checkpoint Control");
    m_checkpoint_control.frequency_units           = pl.get<std::string>("frequency_units");

    if (m_checkpoint_control.frequency_units!="never" && m_checkpoint_control.frequency_units!="none") {
      m_checkpoint_control.timestamp_of_last_write   = run_t0;
      m_checkpoint_control.frequency = pl.get<int>("Frequency");
      EKAT_REQUIRE_MSG (m_output_control.frequency>0,
          "Error! Invalid frequency (" + std::to_string(m_checkpoint_control.frequency) + ") in Checkpoint Control. Please, use positive number.\n");

      // File specs
      m_checkpoint_file_specs.max_snapshots_in_file = 1;
      m_checkpoint_file_specs.filename_with_mpiranks    = pl.get("MPI Ranks in Filename",false);
      m_checkpoint_file_specs.save_grid_data = false;
    }
  } else {
    // If there is no restart data or there is but no checkpoint control sublist then we initialize
    // the checkpoint control so that it never writes checkpoints.
    m_checkpoint_control.frequency_units = "never";
  }

  // If this is normal output (not the model restart output) and the output specs
  // require it, we need to restart the output history.
  // E.g., we might save 30-day avg value for field F, but due to job size
  // break the run into three 10-day runs. We then need to save the state of
  // our averaging in a "restart" file (e.g., the current avg).
  // Note: the user might decide *not* to restart the output, so give the option
  //       of disabling the restart. Also, the user might want to change the
  //       filename_prefix, so allow to specify a different filename_prefix for the restart file.
  if (m_is_restarted_run) {
    // Allow to skip history restart, or to specify a filename_prefix for the restart file
    // that is different from the filename_prefix of the current output.
    auto& restart_pl = m_params.sublist("Restart");
    bool perform_history_restart = restart_pl.get("Perform Restart",true);
    auto hist_restart_casename = restart_pl.get("filename_prefix",m_casename);

    if (perform_history_restart) {
      if (has_restart_data) {
        // This is non-instantaneous output, with freq>1 or freq=1 and freq_units!=nsteps
        // We need to read
        using namespace scorpio;
        auto fn = find_filename_in_rpointer(hist_restart_casename,false,m_io_comm,m_run_t0);

        // We will use this to see if the last output file still has room.
        m_output_control.nsamples_since_last_write = get_attribute<int>(fn, "num_samples_since_last_write");

        // If the type/freq of output needs restart data, we need to restart the streams
        if (m_output_control.nsamples_since_last_write>0) {
          for (auto stream : m_output_streams) {
            stream->restart(fn);
          }
        }
      }

      // Whether we do have restart data or not, let's find the last output file that was
      // created, to a) check if there's still room in it, and b) ensure we are not
      // changing the control settings for output.
      auto cwd = filesystem_ns::current_path();
      auto match = m_casename
                 + "." + e2str(m_avg_type)
                 + "." + m_output_control.frequency_units
                 + "_x" + std::to_string(m_output_control.frequency);

      util::TimeStamp last_file_ts;
      for (auto entry : filesystem_ns::directory_iterator(cwd)) {
        std::string fn = entry.path().filename();
        if (fn.substr(0,match.size())==match) {
          auto t_str = fn.substr(fn.size()-19,16);
          auto t = util::str_to_time_stamp(t_str);
          EKAT_REQUIRE_MSG (t.is_valid(),
              "Error! Something went wrong while extracting time string from filename.\n"
              " - filename: " + fn + "\n"
              " - time str: " + t_str + "\n");
          if (not last_file_ts.is_valid() || last_file_ts<t) {
            last_file_ts = t;
          }
        }
      }

      // If no output file was found from previous runs, last_write_ts remains m_case_t0.
      // This is correct, since INSTANT output would have written at least t=0, so this MUST
      // be some average type. And since no snap was written, we are still in the first
      // averaging window, which starts precisely at m_case_t0.
      if (last_file_ts.is_valid()) {
        auto last_output_fname = compute_filename(m_output_control,m_output_file_specs,false,last_file_ts);
        // There was at least one snapshot written in the previous run, so there was a file.

        // Check if we need to resume filling the output file
        scorpio::register_file(last_output_fname,scorpio::Read);
        int num_snaps = scorpio::get_dimlen_c2f(last_output_fname.c_str(),"time");
        scorpio::eam_pio_closefile(last_output_fname);
        m_resume_output_file = num_snaps<m_output_file_specs.max_snapshots_in_file;

        // Check consistency of output specs across restart
        auto old_freq = scorpio::get_attribute<int>(last_output_fname,"frequency");
        auto old_freq_units = scorpio::get_attribute<std::string>(last_output_fname,"frequency_units");
        auto old_avg_type = scorpio::get_attribute<std::string>(last_output_fname,"avg_type");
        EKAT_REQUIRE_MSG (old_freq == m_output_control.frequency,
            "Error! Cannot change frequency when performing history restart.\n"
            "  - old freq: " << old_freq << "\n"
            "  - new freq: " << m_output_control.frequency << "\n");
        EKAT_REQUIRE_MSG (old_freq_units == m_output_control.frequency_units,
            "Error! Cannot change frequency units when performing history restart.\n"
            "  - old freq units: " << old_freq_units << "\n"
            "  - new freq units: " << m_output_control.frequency_units << "\n");
        EKAT_REQUIRE_MSG (old_avg_type == e2str(m_avg_type),
            "Error! Cannot change avg type when performing history restart.\n"
            "  - old avg type: " << old_avg_type + "\n"
            "  - new avg type: " << e2str(m_avg_type) << "\n");

        // We can also check the time of the last write
        scorpio::register_file(last_output_fname,scorpio::Read);
        auto time = scorpio::read_curr_time_c2f(last_output_fname.c_str());
        scorpio::eam_pio_closefile(last_output_fname);

        m_output_control.timestamp_of_last_write = m_case_t0 + std::round(time*86400);
      }

      // If we need to resume output file, let's open the file immediately, so the run method remains the same
      if (m_resume_output_file) {
        setup_file(m_output_file_specs,m_output_control,last_file_ts);
      }
    }
  }


  if (m_avg_type!=OutputAvgType::Instant) {
    // Init the left hand point of time_bnds based on run/case t0.
    m_time_bnds.resize(2);
    m_time_bnds[0] = m_run_t0.days_from(m_case_t0);
  } else if (m_output_control.output_enabled() && m_run_t0==m_case_t0 && !m_is_model_restart_output) {
    this->run(m_run_t0);
  }
}

void OutputManager::
add_global (const std::string& name, const ekat::any& global) {
  EKAT_REQUIRE_MSG (m_globals.find(name)==m_globals.end(),
      "Error! Global attribute was already set in this output manager.\n"
      "  - global att name: " + name + "\n");

  m_globals[name] = global;
}

/*===============================================================================================*/
void OutputManager::run(const util::TimeStamp& timestamp)
{
  // In case output is disabled, no point in doing anything else
  if (m_output_disabled) {
    return;
  }
  using namespace scorpio;

  std::string timer_root = m_is_model_restart_output ? "EAMxx::IO::restart" : "EAMxx::IO::standard";
  start_timer(timer_root); 
  // Check if we need to open a new file
  ++m_output_control.nsamples_since_last_write;
  ++m_checkpoint_control.nsamples_since_last_write;

  // Check if this is a write step (and what kind)
  const bool is_t0_output       = timestamp==m_case_t0;
  const bool is_output_step     = m_output_control.is_write_step(timestamp) || is_t0_output;
  const bool is_checkpoint_step = m_checkpoint_control.is_write_step(timestamp) && not is_output_step;
  const bool is_write_step      = is_output_step || is_checkpoint_step;

  // If neither output or checkpoint, these won't be used anyways,
  // so no need to check if is_write_step == true.
  auto& control   = is_checkpoint_step ? m_checkpoint_control : m_output_control;
  auto& filespecs = is_checkpoint_step ? m_checkpoint_file_specs : m_output_file_specs;
  auto& filename  = filespecs.filename;

  // Compute filename (if write step)
  start_timer(timer_root+"::get_new_file"); 
  if (is_write_step) {
    // Check if we need to open a new file
    if (not filespecs.is_open) {
      // Register all dims/vars, write geometry data (e.g. lat/lon/hyam/hybm)
      setup_file(filespecs,control,timestamp);
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
    pio_update_time(filename,timestamp.days_from(m_case_t0));
    if (m_is_model_restart_output) {
      // Only write nsteps on model restart
      set_attribute(filename,"nsteps",timestamp.get_num_steps());
    } else if (is_checkpoint_step) {
      // Update the date of last write
      scorpio::write_timestamp (filename,"last_write",control.timestamp_of_last_write);
    }
  }
  stop_timer(timer_root+"::get_new_file"); 

  // Run the output streams
  start_timer(timer_root+"::run_output_streams"); 
  for (auto& it : m_output_streams) {
    // Note: filename might reference an invalid string, but it's only used
    //       in case is_write_step=true, in which case it will *for sure* contain
    //       a valid file name.
    it->run(filename,is_write_step,m_output_control.nsamples_since_last_write,is_t0_output);
  }
  stop_timer(timer_root+"::run_output_streams"); 

  if (is_write_step) {
    for (const auto& it : m_globals) {
      const auto& name = it.first;
      const auto& any = it.second;
      set_any_attribute(filename,name,any);
    }

    start_timer(timer_root+"::update_snapshot_tally"); 
    // We're adding one snapshot to the file
    ++filespecs.num_snapshots_in_file;

    if (m_time_bnds.size()>0) {
      m_time_bnds[1] = timestamp.days_from(m_case_t0);
      scorpio::grid_write_data_array(filename, "time_bnds", m_time_bnds.data(), 2);
      m_time_bnds[0] = m_time_bnds[1];
    }

    if (is_checkpoint_step) {
      set_attribute (filename,"num_samples_since_last_write",control.nsamples_since_last_write);
    }

    // Since we wrote to file we need to reset the nsamples_since_last_write, the timestamp ...
    control.nsamples_since_last_write = 0;
    control.timestamp_of_last_write = timestamp;
    // ... and the local views - unless it is a checkpoint step, then we keep local views.
    if (!is_checkpoint_step) {
      for (auto& it : m_output_streams) {
        it->reset_dev_views();
      }
    }

    // Check if we need to close the output file
    if (filespecs.file_is_full()) {
      eam_pio_closefile(filename);
      filespecs.num_snapshots_in_file = 0;
      filespecs.is_open = false;
    }

    // Whether we wrote an output or a checkpoint, the checkpoint counter needs to be reset
    m_checkpoint_control.nsamples_since_last_write = 0;
    m_checkpoint_control.timestamp_of_last_write = timestamp;

    stop_timer(timer_root+"::update_snapshot_tally"); 
  }
  stop_timer(timer_root); 
}
/*===============================================================================================*/
void OutputManager::finalize()
{
  // Close any output file still open
  if (m_output_file_specs.is_open) {
    scorpio::eam_pio_closefile (m_output_file_specs.filename);
  }
  if (m_checkpoint_file_specs.is_open) {
    scorpio::eam_pio_closefile (m_checkpoint_file_specs.filename);
  }

  // Swapping with an empty mgr is the easiest way to cleanup.
  OutputManager other;
  std::swap(*this,other);
}

long long OutputManager::res_dep_memory_footprint () const {
  long long mf = 0;
  for (const auto& os : m_output_streams) {
    mf += os->res_dep_memory_footprint();
  }

  return mf;
}

std::string OutputManager::
compute_filename (const IOControl& control,
                  const IOFileSpecs& file_specs,
                  const bool is_checkpoint_step,
                  const util::TimeStamp& timestamp) const
{
  std::string suffix =
    is_checkpoint_step ? ".rhist"
                       : (m_is_model_restart_output ? ".r" : "");
  auto filename = m_casename + suffix;

  // Always add avg type and frequency info
  filename += "." + e2str(m_avg_type);
  filename += "." + control.frequency_units+ "_x" + std::to_string(control.frequency);

  // Optionally, add number of mpi ranks (useful mostly in unit tests, to run multiple MPI configs in parallel)
  // NOTE: we do *not* allow this for checkpoints, since it would be risky if it gets somehow enabled
  //       inside an ERP cime test.
  if (not is_checkpoint_step && file_specs.filename_with_mpiranks) {
    filename += ".np" + std::to_string(m_io_comm.size());
  }

  // Always add a time stamp
  if (m_avg_type==OutputAvgType::Instant || is_checkpoint_step) {
    filename += "." + timestamp.to_string();
  } else {
    filename += "." + control.timestamp_of_last_write.to_string();
  }

  return filename + ".nc";
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
      fields_pl.sublist(it.first).set("Field Names",fnames);
    }
    m_casename = m_params.get<std::string>("filename_prefix");
    // Match precision of Fields
    m_params.set<std::string>("Floating Point Precision","real");
  } else {
    auto avg_type = m_params.get<std::string>("Averaging Type");
    m_avg_type = str2avg(avg_type);
    EKAT_REQUIRE_MSG (m_avg_type!=OutputAvgType::Invalid,
        "Error! Unsupported averaging type '" + avg_type + "'.\n"
        "       Valid options: Instant, Max, Min, Average. Case insensitive.\n");

    m_output_file_specs.max_snapshots_in_file = m_params.get<int>("Max Snapshots Per File",-1);
    m_casename = m_params.get<std::string>("filename_prefix");

    // Allow user to ask for higher precision for normal model output,
    // but default to single to save on storage
    const auto& prec = m_params.get<std::string>("Floating Point Precision", "single");
    EKAT_REQUIRE_MSG (prec=="single" || prec=="double" || prec=="real",
        "Error! Invalid floating point precision '" + prec + "'.\n");
  }
}
/*===============================================================================================*/
void OutputManager::
setup_file (      IOFileSpecs& filespecs,
            const IOControl& control,
            const util::TimeStamp& timestamp)
{
  using namespace scorpio;

  const bool is_checkpoint_step = &control==&m_checkpoint_control;
  auto& filename = filespecs.filename;

  // Compute new file name
  // If this is normal output, with some sort of average, and we're not resuming an existing
  // output file, then the timestamp should be the one of the last write, since that's when
  // the current avg window started.
  auto file_ts = m_resume_output_file || m_avg_type==OutputAvgType::Instant || is_checkpoint_step
               ? timestamp : control.timestamp_of_last_write;

  filename = compute_filename (control,filespecs,is_checkpoint_step,file_ts);

  // Register new netCDF file for output. Check if we need to append to an existing file
  auto mode = m_resume_output_file ? Append : Write;
  register_file(filename,mode);

  if (m_resume_output_file) {
    // No need to add time dims/vars to the nc file, simply add them to our online scream-scorpio metadata
    get_variable(filename,"time","time",{"time"}, "double","time");
    if (m_avg_type!=OutputAvgType::Instant) {
      get_variable(filename,"time_bnds","time_bnds",{"dim2","time"},"double","time-dim2");
      scorpio::offset_t time_bnds_dofs[2] = {0,1};
      set_dof(filename,"time_bnds",2,time_bnds_dofs);
    }
  } else {
    // Note: length=0 is how scorpio recognizes that this is an 'unlimited' dimension, which
    // allows to write as many timesnaps as we desire.
    register_dimension(filename,"time","time",0,false);

    // Register time as a variable.
    auto time_units="days since " + m_case_t0.get_date_string() + " " + m_case_t0.get_time_string();
    register_variable(filename,"time","time",time_units,{"time"}, "double", "double","time");
#ifdef SCREAM_HAS_LEAP_YEAR
    set_variable_metadata (filename,"time","calendar","gregorian");
#else
    set_variable_metadata (filename,"time","calendar","noleap");
#endif
    if (m_avg_type!=OutputAvgType::Instant) {
      // First, ensure a 'dim2' dimension with len=2 is registered.
      register_dimension(filename,"dim2","dim2",2,false);
      
      // Register time_bnds var, with its dofs
      register_variable(filename,"time_bnds","time_bnds",time_units,{"dim2","time"},"double","double","time-dim2");
      scorpio::offset_t time_bnds_dofs[2] = {0,1};
      set_dof(filename,"time_bnds",2,time_bnds_dofs);

      // Make it clear how the time_bnds should be interpreted
      set_variable_metadata(filename,"time_bnds","note","right endpoint accummulation");

      // I'm not sure what's the point of this, but CF conventions seem to require it
      set_variable_metadata (filename,"time","bounds","time_bnds");
    }

    // Add info regarding freq, freq_units, and avg_type
    set_attribute(filename,"frequency_units",control.frequency_units);
    set_attribute(filename,"frequency",control.frequency);
    set_attribute(filename,"avg_type",e2str(m_avg_type));
  }

  // Set degree of freedom for "time"
  scorpio::offset_t time_dof[1] = {0};
  set_dof(filename,"time",0,time_dof);

  std::string fp_precision = is_checkpoint_step
                           ? "real"
                           : m_params.get<std::string>("Floating Point Precision");

  // Make all output streams register their dims/vars
  for (auto& it : m_output_streams) {
    it->setup_output_file(filename,fp_precision,m_resume_output_file);
  }

  // If grid data is needed,  also register geo data fields. Skip if file is resumed,
  // since grid data was written in the previous run
  if (filespecs.save_grid_data and not m_resume_output_file) {
    for (auto& it : m_geo_data_streams) {
      it->setup_output_file(filename,fp_precision,false);
    }
  }

  // Finish the definition phase for this file.
  auto t0_date = m_case_t0.get_date()[0]*10000 + m_case_t0.get_date()[1]*100 + m_case_t0.get_date()[2];
  auto t0_time = m_case_t0.get_time()[0]*10000 + m_case_t0.get_time()[1]*100 + m_case_t0.get_time()[2];

  set_attribute(filename,"start_date",t0_date);
  set_attribute(filename,"start_time",t0_time);
  set_attribute(filename,"averaging_type",e2str(m_avg_type));
  set_attribute(filename,"averaging_frequency_units",m_output_control.frequency_units);
  set_attribute(filename,"averaging_frequency",m_output_control.frequency);
  set_attribute(filename,"max_snapshots_per_file",m_output_file_specs.max_snapshots_in_file);
  set_file_header(filename);
  eam_pio_enddef (filename); 

  if (m_avg_type!=OutputAvgType::Instant) {
    // Unfortunately, attributes cannot be set in define mode (why?), so this could
    // not be done while we were setting the time_bnds
    set_attribute(filename,"sample_size",control.frequency);
  }

  if (filespecs.save_grid_data) {
    // Immediately run the geo data streams
    for (const auto& it : m_geo_data_streams) {
      it->run(filename,true,0);
    }
  }

  filespecs.is_open = true;

  m_resume_output_file = false;
}
/*===============================================================================================*/
void set_file_header(const std::string& filename)
{
  using namespace scorpio;

  // TODO: All attributes marked TODO below need to be set.  Hopefully by a universal value that reflects
  // what the attribute is.  For example, git-hash should be the git-hash associated with this version of
  // the code at build time for this executable.

  // TODO: probably want to make sure that new versions are reflected here.
  set_attribute<std::string>(filename,"source","E3SM Atmosphere Model Version 4 (EAMxx)");
  set_attribute<std::string>(filename,"case","");  // TODO
  set_attribute<std::string>(filename,"title","EAMxx History File");
  set_attribute<std::string>(filename,"compset","");  // TODO
  set_attribute<std::string>(filename,"git_hash","");  // TODO
  set_attribute<std::string>(filename,"host","");  // TODO
  set_attribute<std::string>(filename,"version","");  // TODO
  set_attribute<std::string>(filename,"initial_file","");  // TODO
  set_attribute<std::string>(filename,"topography_file","");  // TODO
  set_attribute<std::string>(filename,"contact","");  // TODO
  set_attribute<std::string>(filename,"institution_id","");  // TODO
  set_attribute<std::string>(filename,"product","");  // TODO
  set_attribute<std::string>(filename,"component","ATM");
  set_attribute<std::string>(filename,"conventions","");  // TODO
}

} // namespace scream
