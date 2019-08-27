#ifndef SIMPLE_IO_MOD_HPP
#define SIMPLE_IO_MOD_HPP

#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace simpleio {

// TODO, compile with standard netCDF, at which point these will already be defined.
#ifdef SCREAM_DOUBLE_PRECISION
const int NC_REAL = 6;
#else
  int NC_REAL = 5;
#endif
const int NC_INT = 4;


int hello_world(int myint, double myreal, char** mychar,int len_planet, std::string (&planets)[len_planet]);
int init_output1(char** filename, int ndims, std::string (&dimnames)[ndims], int* dimrng);
int init_output2();
int finalize_output();
int regfield(char** field_name,int field_type,int ndim,std::string (&field_dim)[ndim],char** units);

} // namespace simpleio
} // namespace scream

#endif
