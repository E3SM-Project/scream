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
int indice_convert_c_to_fortran(const int c_ind, const int ndims, const int cdims[ndims], const int fdims_map[ndims])
{
  int f_ind=0;
  int nhat;
  /*
   * Step 1: back out the multi-dimensional index in the c++ array from the n location.
   * This follows the recursive definition that,
   *   For A[n] -> A[d0][d1]...[dN], with dimension ranges D0, D1,..., DN
   *   let nh = n
   *   d0 = floor(nh/prod(Dj,j=1,N)),
   *   nh = nh - d0*prod(Dj,j=1,N),
   *   d(1) = floor(nh/prod(Dj,j=2,N)),
   *   and so on.
   */
  nhat = c_ind;
  int loc[ndims];
  int Dprod=1;
  for (int ii=1;ii<ndims;ii++) {
    Dprod = Dprod * cdims[ii];
  }
  for (int di = 0;di<ndims-1;++di) {
    loc[di] = nhat/Dprod; // Taking advantage of the fact the integer quotient always rounds down.
    nhat -= loc[di]*Dprod;
    Dprod /= cdims[di+1];
  }
  loc[ndims-1] = nhat;
  /*
   * Step 2: Now that we have the indices in c++ array space determine the index in fortran array space
   * Use the following formula to determine the location for the C++ value,
   *     Given an array A with dimensions [D0][D1]...[DN],
   *     the element, A[d0][d1][d2]...[dN] is stored at location n such that,
   *     n = d0 + SUM(di*prod(Dj,j=0,i-1),i=1,N)
   * for example given the array with bounds A[2][4][6], location A[1][3][5] is stored at A[n] where,
   *   n = 1 + 3*2 + 5*4*2 = 47
   * thus A[1][3][5] is equivalent to *A[47]
   * fdims_map has been defined such that each c_dim has it's corresponding product for construction of n-out
   */
  for (int di=0;di<ndims;di++) {
    f_ind += loc[di]*fdims_map[di];  
  }

  return f_ind;
}
/* ----------------------------------------------------------------- */
int indice_convert_fortran_to_c(const int f_ind, const int ndims, const int fdims[ndims], const int cdims_map[ndims])
{
  int c_ind=0;
  int nhat;
  /*
   * Step 1: back out the multi-dimensional index in the fortran array from the n location.
   * This follows the recursive definition that,
   *   For A[n] -> A[d0][d1]...[dN], with dimension ranges D0, D1,..., DN
   *   let nh = n
   *   dN = floor(nh/prod(Dj,j=0,N-1)),
   *   nh = nh - dN*prod(Dj,j=0,N-1),
   *   d(N-1) = floor(nh/prod(Dj,j=0,N-2)),
   *   and so on.
   */
  nhat = f_ind;
  int loc[ndims];
  int Dprod=1;
  for (int ii=0;ii<ndims-1;ii++) {
    Dprod = Dprod * fdims[ii];
  }
  for (int di = ndims-1;di>0;--di) {
    loc[di] = nhat/Dprod; // Taking advantage of the fact the integer quotient always rounds down.
    nhat -= loc[di]*Dprod;
    Dprod /= fdims[di-1];
  }
  loc[0] = nhat;
  /*
   * Step 2: Now that we have the indices in fortran array space determine the index in C++ array space
   * Use the following formula to determine the location for the C++ value,
   *     Given an array A with dimensions [D0][D1]...[DN],
   *     the element, A[d0][d1][d2]...[dN] is stored at location n such that,
   *     n = SUM(di*prod(Dj,j=i+1,N),i=0,N-1) + dN
   * for example given the array with bounds A[2][4][6], location A[1][3][5] is stored at A[n] where,
   *   n = 1*(4*6) + 3*6 + 5 = 47
   * thus A[1][3][5] is equivalent to *A[47]
   * cdims_map has been defined such that each f_dim has it's corresponding product for construction of n-out
   */
  for (int di=0;di<ndims;di++) {
    c_ind += loc[di]*cdims_map[di];  
  }

  return c_ind;
}
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
int init_output1(char** filename, int ndims, std::string (&dimnames)[ndims],
    int* dimrng)
{
  simple_out_init1(filename,ndims,dimnames,dimrng);

  return 0;
}
/* ----------------------------------------------------------------- */
int writefield(char** field_name, int time_dim, int dlen[], double field_data[])
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

  // unwrap c++ field data to match fortran format
  // we assume that dlen and len have the same dimensions but possibly in a different order
  // we assume that no two dimnsions have the same range.  TODO: fix this to be more sophisticated in case two dimensions do have the same range.
  int fdims_map[nspatial_dims];
  for (int ii=0; ii<nspatial_dims;ii++) {
    for (int jj=0; jj<nspatial_dims;jj++) {
      if (len[ii]==dlen[jj]) {
        fdims_map[ii]=1;
        for (int kk=0;kk<jj;kk++) {
          fdims_map[ii] *= dlen[kk];
        }
        break;
      }
    }
  }
  /* 
   * C++ stores arrays in row-column form,
   * Fortran stores arrays in column-row form.
   * Thus we can follow the c++ array via the len vector to count in row-column through the field_data vector.
   * And map those dimensions to the appropriate dimensional location in fortran using dlen.
   */
  double write_data[flen];
  int f_ind = 0;
  for (int c_ind=0; c_ind<flen; c_ind++) {
    f_ind = indice_convert_c_to_fortran(c_ind, nspatial_dims, len, fdims_map);
    write_data[f_ind] = field_data[c_ind];
  }
  
  io_writefield_real(field_name,time_in,flen,write_data);

  return 0;
}
/* ----------------------------------------------------------------- */
int readfield(char** field_name, int time_dim, int dlen[], double res[])
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
  int len_chk[2] = {1,0};
  for (int i=0;i<nspatial_dims;i++) {
    flen *=len[i];
    len_chk[0] *= dlen[i];
    len_chk[1] += len[i];  // len_chk[0] should equal flen if dlen and len have the same set of dimensions
    len_chk[1] -= dlen[i]; // len_chk[1] should equal zero if dlen and len have the same set of dimensions
  }
//  scream::scream_require_msg(flen==dlen,
//                     "Error! Inconsistency in field lengths when loading data. \n"
//                     "Check field info for: " + field_name + "\n");
  
  double field_data[flen];

  io_readfield_real(field_name,time_in,flen,field_data);

  // unwrap fortran field data to match C++ format
  // we assume that dlen and len have the same dimensions but possibly in a different order
  // we assume that no two dimnsions have the same range.  TODO: fix this to be more sophisticated in case two dimensions do have the same range.
  int cdims_map[nspatial_dims];
  for (int ii=0; ii<nspatial_dims;ii++) {
    for (int jj=0; jj<nspatial_dims;jj++) {
      if (len[ii]==dlen[jj]) {
        cdims_map[ii]=1;
        for (int kk=jj+1;kk<nspatial_dims;kk++) {
          cdims_map[ii] *= dlen[kk];
        }
        break;
      }
    }
  }
  /* 
   * C++ stores arrays in row-column form,
   * Fortran stores arrays in column-row form.
   * Thus we can follow the Fortran array via the len vector to count in column-row through the field_data vector output from readfield.
   * And map those dimensions to the appropriate dimensional location in C++ using dlen.
   */
  int c_ind = 0;
  for (int f_ind=0; f_ind<flen; f_ind++) {
    c_ind = indice_convert_fortran_to_c(f_ind, nspatial_dims, len, cdims_map);
    res[c_ind] = field_data[f_ind];
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

} // namespace simpleio
} // namespace scream
