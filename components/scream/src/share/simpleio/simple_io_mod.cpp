#include "simple_io_mod.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"

using scream::Real;
using scream::Int;
extern "C" {

  void hello_world_f90(Int myint, Real myreal, char** mychar,Int len_planet, std::string (&planets)[len_planet]);
  void simple_out_init1(char** filename,Int ndims,std::string (&dimnames)[ndims],Int* dimrng);
  void simple_out_init2();
  void simple_out_finalize();
  void simple_out_regfield(char** field_name,Int field_type,Int ndim,std::string (&field_dim)[ndim],char** units);
  void io_writefield_1d_real(char** field_name, Int time_dim, Int data_dim,  Real field_data[data_dim]);
  void io_writefield_2d_real(char** field_name, Int time_dim, Int* data_dim, Real field_data[data_dim[0]][data_dim[1]]);
  void io_writefield_3d_real(char** field_name, Int time_dim, Int* data_dim, Real field_data[data_dim[0]][data_dim[1]][data_dim[2]]);
  void simple_in_init(char** field_name);
  void simple_in_finalize();
  void simple_io_get_field_size(char** field_name, Int ndims, Int* ldims);
  void io_readfield_real(char** field_name, Int time_dim[2], Int data_dim,  Real field_data[data_dim]);
  void io_readfield_1d_real(char** field_name, Int time_dim, Int data_dim,  Real field_data[data_dim]);
  //void io_readfield_2d_real(char** field_name, Int time_dim, Int* data_dim, Real field_data[data_dim[0]][data_dim[1]]);
  void io_readfield_2d_real(char** field_name, Int time_dim, Int* data_dim, Real** field_data);
  void io_readfield_3d_real(char** field_name, Int time_dim, Int* data_dim, Real*** field_data);
 // void io_readfield_3d_real(char** field_name, Int time_dim, Int* data_dim, Real field_data[data_dim[0]][data_dim[1]][data_dim[2]]);

} // extern C

namespace scream {
namespace simpleio {
  typedef std::vector<double> dbl_vector;
  typedef std::vector<dbl_vector> dbl_2dmatrix;
  typedef std::vector<dbl_2dmatrix> dbl_3dmatrix;
/* ----------------------------------------------------------------- */
int hello_world(int myint, double myreal, char** mychar,int len_planet, 
    std::string (&planets)[len_planet]) 
{
  
  hello_world_f90(myint, myreal,mychar,len_planet,planets);
  return 0;
}
/* ----------------------------------------------------------------- */
int init_output1(char** filename, int ndims, std::string (&dimnames)[ndims],
    int* dimrng)
{
  simple_out_init1(filename,ndims,dimnames,dimrng);

  return 0;
}
/* ----------------------------------------------------------------- */
/* Overload read field function to handle different array sizes     */
/* ----------------------------------------------------------------- */
//int readfield(char** field_name, int time_dim, dbl_vector& res)
//{
//
//  int ndims = 1;
//  int len[1];
//  simple_io_get_field_size(field_name,ndims,len);
//
//  double *field_data = NULL;
//  field_data = new double[len[0]];
//  io_readfield_1d_real(field_name,time_dim,len[0],field_data);
//
//  for (int i=0; i<len[0]; i++) {
//    res.push_back( field_data[i] );
//  }
//  return 0;
//}
/* ----------------------------------------------------------------- */
int readfield(char** field_name, int time_dim, dbl_3dmatrix& res)
{

  int ndims = 4;
  int len[]={1,1,1,1};
  int time_in[] = {1,1};
  simple_io_get_field_size(field_name,ndims,len);
  if (len[ndims-1]>1) {
    time_in[0] = ndims;
    time_in[1] = time_dim;
  }
  int flen = 1;
  for (const auto& e:len) {
    flen *=e;
  }
  double field_data[flen];

  io_readfield_real(field_name,time_in,flen,field_data);

  int ind = 0;
  for (int i=0; i<len[2]; i++) {
    dbl_2dmatrix mat;
    for (int j=0; j<len[1]; j++) {
      dbl_vector row;
      for (int k=0; k<len[0]; k++) {
        row.push_back( field_data[ind] );
        ind +=1;
      }
      mat.push_back(row);
    }
    res.push_back( mat );
  }

  return 0;
}
/* ----------------------------------------------------------------- */
int readfield(char** field_name, int time_dim, dbl_2dmatrix& res)
{

  int ndims = 3;
  int len[]={1,1,1};
  int time_in[] = {1,1};
  simple_io_get_field_size(field_name,ndims,len);
  if (len[ndims-1]>1) {
    time_in[0] = ndims;
    time_in[1] = time_dim;
  }
  int flen = 1;
  for (const auto& e:len) {
    flen *=e;
  }
  double field_data[flen];

  io_readfield_real(field_name,time_in,flen,field_data);

  int ind = 0;
  for (int j=0; j<len[1]; j++) {
    dbl_vector row;
    for (int k=0; k<len[0]; k++) {
      row.push_back( field_data[ind] );
      ind +=1;
    }
    res.push_back(row);
  }

  return 0;
}
/* ----------------------------------------------------------------- */
int readfield(char** field_name, int time_dim, dbl_vector& res)
{

  int ndims = 2;
  int len[]={1,1};
  int time_in[] = {1,1};
  simple_io_get_field_size(field_name,ndims,len);
  if (len[ndims-1]>1) {
    time_in[0] = ndims;
    time_in[1] = time_dim;
  }
  int flen = 1;
  for (const auto& e:len) {
    flen *=e;
  }
  double field_data[flen];

  io_readfield_real(field_name,time_in,flen,field_data);

  int ind = 0;
  for (int k=0; k<len[0]; k++) {
    res.push_back( field_data[ind] );
    ind +=1;
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
/* Overload write field function to handle different array sizes     */
/* ----------------------------------------------------------------- */
int writefield(char** field_name, int time_dim, int data_dim, double* field_data)
{

 io_writefield_1d_real(field_name,time_dim,data_dim,field_data); 

  return 0;
}
/* ----------------------------------------------------------------- */
int writefield(char** field_name, int time_dim, int* data_dim, double field_data[data_dim[0]][data_dim[1]])
{
  io_writefield_2d_real(field_name,time_dim,data_dim,field_data); 

  return 0;
}
/* ----------------------------------------------------------------- */
int writefield(char** field_name, int time_dim, int* data_dim, double field_data[data_dim[0]][data_dim[1]][data_dim[2]])
{
  io_writefield_3d_real(field_name,time_dim,data_dim,field_data); 

  return 0;
}
/* ----------------------------------------------------------------- */


} // namespace simpleio
} // namespace scream
