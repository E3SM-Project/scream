#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

#include "scream_scorpio_types.hpp"

#include <ekat/mpi/ekat_comm.hpp>
#include <ekat/std_meta/ekat_std_any.hpp>

#include <string>
#include <vector>
#include <memory>

namespace scream {
namespace scorpio {

inline std::string default_time_name () { return "time"; }

// =================== Global operations ================= //

void init_pio_subsystem(const ekat::Comm& comm, const int atm_id = 0);
bool is_pio_subsystem_inited ();
void finalize_pio_subsystem ();

// =================== File operations ================= //

// Opens a file, returns const handle to it (useful for Read mode, to get dims/vars)
void register_file (const std::string& filename, const FileMode mode);

// Release a file (if in Write mode, sync and close the file);
void release_file  (const std::string& filename);

// Check if file is open. If mode!=Unset, also checks that it's open with given mode
bool is_file_open (const std::string& filename, const FileMode mode = Unset);

// Force a flush to file (for Write mode only)
void sync_file (const std::string &filename);

// Reopen/ends the definition phase
void redef (const std::string &filename);
void enddef (const std::string &filename);

// =================== Dimensions operations ======================= //

// Define dim on output file (cannot call on Read/Append files)
void define_dim (const std::string& filename, const std::string& dimname, const int length);

// Check that the given dimension is in the file. If length>0, also check that the length is as expected.
bool has_dimension (const std::string& filename, const std::string& dimname, const int length = -1);

int get_dimlen (const std::string& filename, const std::string& dimname);

// =================== Decompositions operations ==================== //

// Create a decomposition along a particular dimension
// Notes:
// - we declare a decomposition along a single dimension.
// - specifying a decomp_name allows to have have multiple different
//   decompositions for the same dim (used at least once, by CoarseningRemapper)
// - when we call set_var_decomp, we will pick the corresponding decomposition
//   along the specified dimension
// - set_dim_decomp requires *offsets* in the global array, not global indices
//   (for 0-based indices, they're the same, but for 1-based indices they're not)
// - the decomp name is optional, since in most cases we only decompose a dimension
//   in one way; but if you know you will (or might) have multiple decompositions for
//   the same dim, then you should give each decomp a unique name

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const std::vector<offset_t>& my_offsets,
                     const std::string& decomp_name = "DEFAULT");

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const offset_t start, const offset_t count,
                     const std::string& decomp_name = "DEFAULT");

// Initializes (or recycles) a decomposition for the given variable in the given file
// Notes:
//  - dimname specifies the decomposed dimension, while dim_decomp_name specifies
//    which decomposition along that dim should be used (in case set_dim_decomp
//    was called multiple times for that dimension during the lifetime of the PIOFile)
//  - if a decomp is already present, this function throws if the new decomp is different.
//    Use reset_var_decomp (below) to allow a change of decomp
//  - throw_if_var_does_not_have_decomp_dim=false allows to call these method on all vars,
//    without checking first which ones are indeed decomposed
//  - the actual PIO decomp objects are created here
//  - for efficiency reasons, decompositions are kept in memory even after a file is closed,
//    so a call to this function may be quick (if we can recycle an existing decomp)
//  - if you wish to clear any currently unused decomposition, use free_unused_decomps()
void set_var_decomp (const std::string& filename,
                     const std::string& varname,
                     const std::string& dimname,
                     const std::string& dim_decomp_name = "DEFAULT",
                     const bool throw_if_var_does_not_have_decomp_dim = false);

void set_vars_decomp (const std::string& filename,
                      const std::vector<std::string>& varnames,
                      const std::string& dimname,
                      const std::string& dim_decomp_name = "DEFAULT",
                      const bool throw_if_var_does_not_have_decomp_dim = false);

void reset_var_decomp (const std::string& filename,
                       const std::string& varname,
                       const std::string& dimname,
                       const std::string& dim_decomp_name = "DEFAULT",
                       const bool throw_if_var_does_not_have_decomp_dim = false);

void reset_vars_decomp (const std::string& filename,
                        const std::vector<std::string>& varnames,
                        const std::string& dimname,
                        const std::string& dim_decomp_name = "DEFAULT",
                        const bool throw_if_var_does_not_have_decomp_dim = false);
// Clean up any currently unused decompositions (meaning the decomp is associated to *no* variables)
void free_unused_decomps ();

// ================== Variable operations ================== //

// Define var on output file (cannot call on Read/Append files)
void define_var (const std::string& filename, const std::string& varname,
                 const std::string& units, const std::vector<std::string>& dimensions,
                 const std::string& dtype, const std::string& nc_dtype,
                 const bool time_dependent = false);

// Shortcut when units are not used, and dtype==nc_dtype
void define_var (const std::string& filename, const std::string& varname,
                 const std::vector<std::string>& dimensions,
                 const std::string& dtype,
                 const bool time_dependent = false);

void change_var_dtype (const std::string& filename, const std::string& varname,
                       const std::string& dtype);

// Check that the given variable is in the file.
// Notes:
//  - if dims!={}: check also that the dimensions are as expected
//  - if units!="": also check the units attribute
//  - if time_dep is
//    - "yes" or "true" (case insensitive): check that var is time dep
//    - "no" or "false" (case insensitive): check that var is *not* time dep
//    - anything else: don't check time dep flag
bool has_variable (const std::string& filename, const std::string& varname,
                   const std::vector<std::string>& dims = {},
                   const std::string& units = "",
                   const std::string& time_dep = "");

std::vector<std::string> get_vardims (const std::string& filename,
                                      const std::string& varname);

// Defines both a time dimension and a time variable
void define_time (const std::string& filename, const std::string& units,
                  const std::string& time_name = default_time_name());

// Update value of time variable, increasing time dim length
void update_time(const std::string &filename, const double time);

// Retrieves the time variable value(s)
double get_time (const std::string& filename, const int time_index = -1);
std::vector<double> get_all_times (const std::string& filename);

// Read variable into user provided buffer.
// If time dim is present, read given time slice (time_index=-1 means "read last record).
// If time dim is not present, time_index must be -1 (error out otherwise)
// NOTE: ETI in the cpp file for int, float, double.
template<typename T>
void read_var (const std::string &filename, const std::string &varname, T* buf, const int time_index = -1);

// Write data from user provided buffer into the requested variable
// NOTE: ETI in the cpp file for int, float, double.
template<typename T>
void write_var (const std::string &filename, const std::string &varname, const T* buf, const T* fillValue = nullptr);

// =============== Attributes operations ================== //

// To specify GLOBAL attributes, pass "GLOBAL" as varname
ekat::any get_any_attribute (const std::string& filename,
                             const std::string& varname,
                             const std::string& attname);
void set_any_attribute (const std::string& filename,
                        const std::string& varname,
                        const std::string& attname,
                        const ekat::any& att);

template<typename T>
T get_attribute (const std::string& filename,
                 const std::string& varname,
                 const std::string& attname)
{
  auto att = get_any_attribute(filename,varname,attname);
  return ekat::any_cast<T>(att);
}

template<typename T>
T get_global_attribute (const std::string& filename,
                        const std::string& attname)
{
  return get_attribute<T>(filename,"GLOBAL",attname);
}

template<typename T>
void set_attribute (const std::string& filename,
                    const std::string& varname,
                    const std::string& attname,
                    const T& att)
{
  ekat::any a(att);
  set_any_attribute(filename,varname,attname,a);
}

template<typename T>
void set_global_attribute (const std::string& filename,
                           const std::string& attname,
                           const T& att)
{
  set_attribute(filename,"GLOBAL",attname,att);
}

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
