#ifndef SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP
#define SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_process.hpp"

namespace scream
{

/*
 * A class representing an atmosphere diagnostic.
 * Diagnostics are similar to atmosphere processes in that
 * they take field(s) as input and produce some output.
 * 
 * This subclass definition makes a few distinctions that are
 * unique to diagnostics.
 * 
 * 1) Fields from the field manager can only be accessed as "Required".
 *    Meaning that an atmosphere diagnostic can not compute or update
 *    a field in the field manager.
 *    - note, if you want to update something in the field manager it
 *            should happen in an atmosphere process.
 * 2) Diagnostics produce a single field output.
 *    Typically the field output will be to buffered memory, meaning it
 *    could be overwritten later in the simulation.
 *    A diagnostic output is meant to be used for OUTPUT or a PROPERTY CHECK.
 */

// TODO: inheriting from AtmosphereProcess is conceptually wrong. It was done
//       out of convenience, but we should revisit that choice.

class AtmosphereDiagnostic : public AtmosphereProcess
{
public:

  // Constructor(s)
  AtmosphereDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  virtual ~AtmosphereDiagnostic () = default;

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // Getting the diagnostic output
  Field get_diagnostic () const;

  // Avoid shadowing the base class method, but allow calling without passing a dummy dt
  using AtmosphereProcess::run;
  void run () { this->run_impl(0); } // Pass null dt, but should not be used anyways.

  void set_computed_field (const Field& f) final;
  void set_computed_group (const FieldGroup& group) final;
protected:

  // Diagnostics are meant to return a field
  Field m_diagnostic_output;
};

// A short name for the factory for atmosphere diagnostics
// WARNING: you do not need to write your own creator function to register your atmosphere diagnostic in the factory.
//          You could, but there's no need. You can simply register the common one right below, using your
//          atmosphere diagnostic subclass name as templated argument. If you roll your own creator function, you
//          *MUST* ensure that it correctly sets up the self pointer after creating the shared_ptr.
//          This is *necessary* until we can safely switch to std::enable_shared_from_this.
//          For more details, see the comments at the top of ekat_std_enable_shared_from_this.hpp.
using AtmosphereDiagnosticFactory =
    ekat::Factory<AtmosphereDiagnostic,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<AtmosphereDiagnostic>,
                  const ekat::Comm&,const ekat::ParameterList&>;

// Create an atmosphere process, and correctly set up the (weak) pointer to self.
template <typename AtmDiagType>
inline std::shared_ptr<AtmosphereDiagnostic>
create_atmosphere_diagnostic (const ekat::Comm& comm, const ekat::ParameterList& p) {
  auto ptr = std::make_shared<AtmDiagType>(comm,p);
  ptr->setSelfPointer(ptr);
  return ptr;
}
} //namespace scream

#endif // SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP
