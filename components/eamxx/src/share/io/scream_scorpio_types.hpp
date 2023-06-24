#ifndef SCREAM_SCORPIO_TYPES_HPP
#define SCREAM_SCORPIO_TYPES_HPP

#include <string>
#include <map>
#include <vector>
#include <memory>

namespace scream {
namespace scorpio {

// Enum denoting how we open a file
enum FileMode {
  Unset = 0,
  Read = 1,
  Write = 2,
  Append = Read | Write
};

std::string e2str (const FileMode mode);

// The type used by PIOc for offsets
using offset_t = std::int64_t;

/*
 * The following PIOxyz types each represent an entity that
 * is associated in scorpio with an id. For each of them, we add some
 * additional metadata, to avoid having to call PIOc_inq_xyz every time
 * we need such information.
 * While these types are defined in this public header, the customers
 * of scream_io will never be able to get any object of these types out
 * of the internal database. In particular, all data is stored in a
 * ScorpioSession singleton class, whose declaration is hidden inside
 * scream_scorpio_interface.cpp.
 */

// The basic common data of any PIO entity
struct PIOEntity {
  int ncid = -1;            // PIO access all data via their id
  std::string name;         // In EAMxx, we prefer to use a name
};

// A dimension
struct PIODim : public PIOEntity {
  int length       = -1;
  bool unlimited   = false;

  // A map decomp_name->decomp, for all decompositions
  // that this dimension is given. In general, we only use one decomp,
  // but there is one instance where we first read a file with a simple
  // decomp (linear), and later we read the same file in a scattered decomp.
  // See share/grid/remap/coarsening_remapper.cpp for an example.
  std::map<std::string,std::vector<offset_t>> decomps;
};

// A decomposition
struct PIODecomp : public PIOEntity {
  std::vector<offset_t> offsets;
};

// A variable
struct PIOVar : public PIOEntity {
  // Note: if time_dep=true, we will add it to the list of dims passed
  // to scorpio, but the time dim will not appear in this list.
  std::vector<std::shared_ptr<const PIODim>> dims;

  std::string dtype;
  std::string nc_dtype;
  std::string units;

  bool time_dep = false;
  
  // Extra safety measure: for time_dep vars, use this to check that
  // we are not writing more slices than the current time dim length.
  int num_records = 0;

  std::string decomp_dim;
  std::string dim_decomp_name;
  std::shared_ptr<const PIODecomp> decomp;

  // Used only if a) var is not decomposed, and b) dtype!=nc_dtype
  int size = -1; // Product of all dims
  std::vector<char> buf;
};

// A file, which is basically a container for dims and vars
struct PIOFile : public PIOEntity {
  std::map<std::string,std::shared_ptr<PIODim>>   dims;
  std::map<std::string,std::shared_ptr<PIOVar>>   vars;

  std::shared_ptr<PIODim> time_dim;
  FileMode mode;
  bool enddef = false;

  // We keep track of how many places are currently using this file, so that we
  // can close it only when they are all done.
  int num_customers = 0;
};

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
