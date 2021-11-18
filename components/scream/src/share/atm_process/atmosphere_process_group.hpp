#ifndef SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP
#define SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP

#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_parameter_list.hpp"

#include <string>
#include <list>

namespace scream
{

/*
 *  A class representing a group of atmosphere processes as a single process.
 *
 *  This class allows to create nested lists of processes, by representing a group
 *  of processes as a single atmosphere process.
 *  All the calls to setup/run methods are simply forwarded to the stored list of
 *  atm processes, and the stored list of required/computed fields is simply a
 *  concatenation of the correspong lists in the underlying atm processes.
 *  The only caveat is required fields in sequential scheduling: if an atm proc
 *  requires a field that is computed by a previous atm proc in the group,
 *  that field is not exposed as a required field of the group.
 */

class AtmosphereProcessGroup : public AtmosphereProcess
{
public:
  using atm_proc_type     = AtmosphereProcess;

  // Constructor(s)
  explicit AtmosphereProcessGroup (const ekat::Comm& comm, const ekat::ParameterList& params);

  virtual ~AtmosphereProcessGroup () = default;

  // The type of the block (e.g., dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Group; }

  // The type of grids on which the process is defined
  std::set<std::string> get_required_grids () const { return m_required_grids; }

  // The name of the block
  std::string name () const { return m_group_name; }

  // Grab the proper grid from the grids manager
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  void final_setup ();

  // --- Methods specific to AtmosphereProcessGroup --- //
  int get_num_processes () const { return m_atm_processes.size(); }

  std::shared_ptr<const atm_proc_type> get_process (const int i) const {
    return m_atm_processes.at(i);
  }

  ScheduleType get_schedule_type () const { return m_group_schedule_type; }

  // Initialize memory buffer for each process
  void initialize_atm_memory_buffer (ATMBufferManager& memory_buffer);

  // The APG class needs to perform special checks before establishing whether
  // a required group/field is indeed a required group for this APG
  void set_required_field (const Field<const Real>& field);
  void set_required_group (const FieldGroup<const Real>& group);

  // Gather internal fields from all processes in the group
  void gather_internal_fields ();

protected:

  // Adds fid to the list of required/computed fields of the group (as a whole).
  void process_required_field (const FieldRequest& req);
  void process_required_group (const GroupRequest& req);

  // The initialization, run, and finalization methods
  void initialize_impl(const RunType run_type);
  void initialize_impl ();
  void run_impl        (const int dt);
  void finalize_impl   (/* what inputs? */);

  void run_sequential (const Real dt);
  void run_parallel   (const Real dt);

  // The methods to set the fields/groups in the right processes of the group
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);
  void set_required_group_impl (const FieldGroup<const Real>& group);
  void set_computed_group_impl (const FieldGroup<      Real>& group);

  // The name of the group. This is usually a concatenation of the names of the individual processes
  std::string       m_group_name;
  int               m_group_size;

  // The list of atm processes in this group
  std::vector<std::shared_ptr<atm_proc_type>>  m_atm_processes;

  // The grids required by this process
  std::set<std::string>  m_required_grids;

  // The schedule type: Parallel vs Sequential
  ScheduleType   m_group_schedule_type;
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP
