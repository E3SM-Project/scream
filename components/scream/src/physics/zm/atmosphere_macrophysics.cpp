#include "ekat/scream_assert.hpp"
#include "physics/zm/scream_zm_interface.hpp"
#include "physics/zm/atmosphere_macrophysics.hpp"

namespace scream
{


// =========================================================================================
void ZMMacrophysics::initialize (const util::TimeStamp& t0)
{
  cout << "In ZMMacrophysics::initialize\n";
  zm_init_f90 (q_ptr);

}

// =========================================================================================
void ZMMacrophysics::run (const Real dt)
{
  count << "In ZMMActophysics::run\n";
  zm_main_f90 (dt,q_ptr,fq_ptr);

}
// =========================================================================================
void SHOCMacrophysics::finalize()
{
  count << "In ZMMacrophysics :: finalize() \n";
  zm_finalize_f90 ();
}
// =========================================================================================


} // namespace scream
