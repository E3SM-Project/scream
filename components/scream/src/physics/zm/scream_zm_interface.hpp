#ifndef SCREAM_SHOC_INTERFACE_HPP
#define SCREAM_SHOC_INTERFACE_HPP

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"

// Put everything into a scream namespace
namespace scream {

extern "C"
{

// Fortran routines to be called from C
void zm_init_f90     ();
void zm_main_f90     ();
void zm_finalize_f90 ();

} // extern "C"

} // namespace scream

#endif // SCREAM_SHOC_INTERFACE_HPP
