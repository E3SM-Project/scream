#ifndef SIMPLE_IO_MOD_HPP
#define SIMPLE_IO_MOD_HPP

#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"
#include <netcdf.h>

namespace scream {
namespace simpleio {

// TODO, compile with standard netCDF, at which point these will already be defined.
#ifdef SCREAM_DOUBLE_PRECISION
const int NC_REAL = NC_DOUBLE;
#else
  int NC_REAL = NC_FLOAT;
#endif

// TODO overload functions so that nc_field is one of the only inputs needed.
struct nc_field {
  std::string fieldname;  // Field Name
  int ftype;              // Data type for field
  std::string units;      // Field Units
  int ndims;              // Dimensionality of field
  std::string dimensions[10]; // List of dimension names  //TODO: find a better way to do this.  Currently we assume that no field will have more than 10 dimensions.
};

typedef std::vector<double> dbl_vector;
typedef std::vector<dbl_vector> dbl_2dmatrix;
typedef std::vector<dbl_2dmatrix> dbl_3dmatrix;

void init_output2_c(int ncid);
void finalize_output_c(int ncid);
void regfield_c(int ncid,std::string field_name,int field_type,int ndim,std::string (&field_dim)[ndim],std::string units);
void writefield_c(int ncid,std::string field_name, double &data, int time_dim);

int init_output1(char** filename, int ndims, std::string (&dimnames)[ndims], int* dimrng);
int init_output1_c(std::string filename, int ndims, std::string (&dimnames)[ndims], int* dimrng);
int init_output2();
int finalize_output();
int init_input(char** filename);
int finalize_input();
int regfield(char** field_name,int field_type,int ndim,std::string (&field_dim)[ndim],char** units);
// Overload write functions to accomodate types TODO: ADD an int version
int writefield(char** field_name, int tim_dim, int dlen[], double res[]);
// Overload read functions to accomodate types TODO: ADD an int version
int readfield(char** field_name, int tim_dim, int dlen[], double res[]);

int simple_cpp_netcdf();

} // namespace simpleio
} // namespace scream

#endif
