#include "simple_io_mod.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"

using scream::Real;
using scream::Int;
extern "C" {


  void simple_out_init1(char** filename,Int ndims,std::string (&dimnames)[ndims],Int* dimrng);
  void simple_out_regfield(char** field_name,Int field_type,Int ndim,std::string (&field_dim)[ndim],char** units);
  void simple_out_init2();
  void simple_out_finalize();

  void simple_in_init(char** field_name);
  void simple_in_finalize();

  void simple_io_get_field_size(char** field_name, Int ncid, Int ndims,  Int* ldims, Int* nspatial_dims);

  void io_readfield_real(char** field_name, Int time_dim[2], Int data_dim,  Real field_data[data_dim]);
  void io_writefield_real(char** field_name, Int time_dim[2], Int data_dim,  Real field_data[data_dim]);

} // extern C

namespace scream {
namespace simpleio {
/* ----------------------------------------------------------------- */
int init_output1(char** filename, int ndims, std::string (&dimnames)[ndims],
    int* dimrng)
{
  simple_out_init1(filename,ndims,dimnames,dimrng);

  return 0;
}
/* ----------------------------------------------------------------- */
int writefield(char** field_name, int time_dim, int dlen, double field_data[])
{

  int ndims = 10;
  int nspatial_dims=0;
  int len[ndims];
  int time_in[] = {1,1};
  simple_io_get_field_size(field_name,1,ndims,len,&nspatial_dims);
  for (int ii=0;ii<nspatial_dims;ii++) {
    if (len[ii]==-999) {
      time_in[0] = ii+1; // Convert to location for Fortran
      time_in[1] = time_dim;
      nspatial_dims = nspatial_dims-1;
      break;
    }
  }

  int flen = 1;
  for (int i=0;i<nspatial_dims;i++) {
    flen *=len[i];
  }
//  scream::scream_require_msg(flen=dlen,
//                     "Error! Inconsistency in field lengths when loading data. \n"
//                     "Check field info for: " + field_name + "\n");
  
  io_writefield_real(field_name,time_in,flen,field_data);

  return 0;
}
/* ----------------------------------------------------------------- */
int readfield(char** field_name, int time_dim, int dlen, double res[])
{

  int ndims = 10;
  int nspatial_dims=0;
  int len[ndims];
  int time_in[] = {1,1};
  simple_io_get_field_size(field_name,0,ndims,len,&nspatial_dims);
  for (int ii=0;ii<nspatial_dims;ii++) {
    if (len[ii]==-999) {
      time_in[0] = ii+1; // Convert to location for Fortran
      time_in[1] = time_dim;
      nspatial_dims = nspatial_dims-1;
      break;
    }
  }

  int flen = 1;
  for (int i=0;i<nspatial_dims;i++) {
    flen *=len[i];
  }
//  scream::scream_require_msg(flen==dlen,
//                     "Error! Inconsistency in field lengths when loading data. \n"
//                     "Check field info for: " + field_name + "\n");
  
  double field_data[flen];

  io_readfield_real(field_name,time_in,flen,field_data);

  // unwrap fortran field data to match C++ format
  int src_ind = 0;
  int tgt_ind = 0;
  for (int i=0; i<flen; i++) {
    src_ind = i;
    res[tgt_ind] = field_data[src_ind];
    tgt_ind +=1;
  }

  return 0;
}
/* ----------------------------------------------------------------- */
int init_output2()
{
  simple_out_init2();

  return 0;
}
/* ----------------------------------------------------------------- */
int finalize_output()
{
  simple_out_finalize();

  return 0;
}
/* ----------------------------------------------------------------- */
int init_input(char** filename)
{
  simple_in_init(filename);

  return 0;
}
/* ----------------------------------------------------------------- */
int finalize_input()
{
  simple_in_finalize();

  return 0;
}
/* ----------------------------------------------------------------- */
int regfield(char** field_name,int field_type,int ndim,std::string (&field_dim)[ndim],char** units)
{
  simple_out_regfield(field_name,field_type,ndim,field_dim,units);

  return 0;
}
/* ----------------------------------------------------------------- */


} // namespace simpleio
} // namespace scream
