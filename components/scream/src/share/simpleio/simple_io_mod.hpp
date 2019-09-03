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

typedef std::vector<double> dbl_vector;
typedef std::vector<dbl_vector> dbl_2dmatrix;
typedef std::vector<dbl_2dmatrix> dbl_3dmatrix;

int hello_world(int myint, double myreal, char** mychar,int len_planet, std::string (&planets)[len_planet]);
int init_output1(char** filename, int ndims, std::string (&dimnames)[ndims], int* dimrng);
int init_output2();
int finalize_output();
int init_input(char** filename);
int finalize_input();
int regfield(char** field_name,int field_type,int ndim,std::string (&field_dim)[ndim],char** units);
// Overload write functions to accomodate varying size of arrays
int writefield(char** field_name, int time_dim, int  data_dim, double* field_data);
int writefield(char** field_name, int time_dim, int* data_dim, double  field_data[data_dim[0]][data_dim[1]]);
int writefield(char** field_name, int time_dim, int* data_dim, double  field_data[data_dim[0]][data_dim[1]][data_dim[2]]);
// Overload read functions to accomodate varying size of arrays
int readfield(char** field_name, int tim_dim, dbl_vector& res);
int readfield(char** field_name, int tim_dim, std::vector< std::vector<double> >& res);
int readfield(char** field_name, int tim_dim, std::vector< std::vector< std::vector <double> > >& res);

} // namespace simpleio
} // namespace scream

#endif
