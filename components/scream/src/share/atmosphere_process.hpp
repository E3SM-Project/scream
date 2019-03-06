#ifndef SCREAM_ATMOSPHERE_PROCESS_HPP
#define SCREAM_ATMOSPHERE_PROCESS_HPP

#include <string>
#include <set>

#include "share/scream_assert.hpp"
#include "share/mpi/scream_comm.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field_repository.hpp"
#include "share/field/field.hpp"

namespace scream
{

enum class AtmosphereProcessType {
  Coupling,   // Process responsible of interfacing with the component coupler
  Dynamics,   // Process responsible of handling the dynamics
  Physics,    // Process handling a physics parametrization
  Group       // Process that groups a bunch of processes (so they look as a single process)
};

/*
 *  The abstract interface of a process of the atmosphere 
 *
 *  The process will handle a particular part of the atmosphere component.
 *  This includes both physics (i.e., parametrizations) and dynamics.
 *  The atmosphere driver will take care of calling init/run/finalize
 *  methods of each process, in an order that the driver
 *  establishes. A process must provide a list of fields
 *  that it needs as input, together with a list of fields that
 *  are computed.
 */

class AtmosphereProcess
{
public:
  using device_type      = DefaultDevice; // may need to template class on this
  using host_device_type = HostDevice;

  virtual ~AtmosphereProcess () = default;

  // The type of the block (e.g., dynamics or physics)
  virtual AtmosphereProcessType type () const = 0;

  // The name of the block 
  virtual std::string name () const = 0;

  // The communicator associated with this atm process
  virtual const Comm& get_comm () const = 0;

  // These are the three main interfaces:
  //   - the initialize method sets up all the stuff the process needs to run,
  //     including arrays/views, parameters, and precomputed stuff.
  //   - the run method time-advances the process by one time step.
  //     We could decide whether we want to assume that other process may
  //     be running at the same time, or whether each process can assume
  //     that no other process is currently inside a call to 'run'.
  //   - the finalize method makes sure, if necessary, that all resources are freed.
  // The initialize/finalize method should be called just once per simulation (should
  // we enforce that? It depends on what resources we init/free, and how), while the
  // run method can (and usually will) be called multiple times.
  // We should put asserts to verify that the process has been init-ed, when
  // run/finalize is called.
  virtual void initialize (const Comm& comm) = 0;
  virtual void run        (/* what inputs? */) = 0;
  virtual void finalize   (/* what inputs? */) = 0;

  // These methods set fields in the atm process. Fields live on device and they are all 1d.
  // If the process *needs* to store the field as n-dimensional field, use the
  // template functio 'reinterpret_field' (see field.hpp for details).
  void set_required_field (const Field<const Real, device_type>& f) {
    error::runtime_check(requires(f.get_header().get_identifier()),
                         "Error! This atmosphere process does not require this field. "
                         "Something is wrong up the call stack. Please, contact developers.\n");
    set_required_field_impl (f);
  }
  void set_computed_field (const Field<Real, device_type>& f) {
    error::runtime_check(computes(f.get_header().get_identifier()),
                         "Error! This atmosphere process does not compute this field. "
                         "Something is wrong up the call stack. Please, contact developers.\n");
    set_computed_field_impl (f);
  }

  // Register required/computed fields in the field repo
  virtual void register_fields (FieldRepository<Real, device_type>& field_repo) const = 0;

  // These two methods allow the driver to figure out what process need
  // a given field and what process updates a given field.
  virtual const std::set<FieldIdentifier>& get_required_fields () const = 0;
  virtual const std::set<FieldIdentifier>& get_computed_fields () const = 0;

  // NOTE: C++20 will introduce the method 'contains' for std::set. Till then, use find and check result.
  bool requires (const FieldIdentifier& id) const { return get_required_fields().find(id)!= get_computed_fields().end(); }
  bool computes (const FieldIdentifier& id) const { return get_computed_fields().find(id)!= get_computed_fields().end(); }

protected:
  virtual void set_required_field_impl (const Field<const Real, device_type>& f) = 0;
  virtual void set_computed_field_impl (const Field<      Real, device_type>& f) = 0;
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_HPP
