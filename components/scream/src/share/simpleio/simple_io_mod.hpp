#ifndef SIMPLE_IO_MOD_HPP
#define SIMPLE_IO_MOD_HPP

#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"
#include <netcdf.h>

namespace scream {
namespace simpleio {

#ifdef SCREAM_DOUBLE_PRECISION
const int NC_REAL = NC_DOUBLE;
#else
const int NC_REAL = NC_FLOAT;
#endif

// TODO overload functions so that nc_field is one of the only inputs needed.
struct nc_field {
  std::string fieldname;  // Field Name
  int ftype;              // Data type for field
  std::string units;      // Field Units
  int ndims;              // Dimensionality of field
  std::vector<std::string> dimensions; // List of dimension names  
};

// Set of netcdf wrappers
using scream::Real;
using scream::Int;
int init_output1(std::string filename, int ndims, std::vector<std::string> &dimnames, int* dimrng);
void init_output2(int ncid);
void regfield(int ncid,std::string field_name,int field_type,int ndim,std::vector<std::string> &field_dim,std::string units);
void writefield(int ncid,std::string field_name, Real &data, int time_dim);
void readfield(int ncid,std::string field_name, Real &data, int time_dim);
int init_input(std::string filename);

void finalize_io(int ncid);

} // namespace simpleio
} // namespace scream

#endif
