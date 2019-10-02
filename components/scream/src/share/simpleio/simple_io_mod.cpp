#include "simple_io_mod.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"
#include <netcdf.h>

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
void convert_string_to_char(std::string str_in,char str_out[str_in.size()+1]) {

  str_in.copy(str_out,str_in.size()+1);
  str_out[str_in.size()] = '\0';

}
/* ----------------------------------------------------------------- */
int init_input_c(std::string filename) {


  int nc_err, ncid;
  char filename_in[filename.size()+1];

  convert_string_to_char(filename,filename_in);
  nc_err = nc_open(filename_in, NC_NOWRITE, &ncid);
  scream_require_msg(nc_err==0,"NC Create Error! \n");

  return ncid;
}
/* ----------------------------------------------------------------- */
int init_output1_c(std::string filename, int ndims, std::string (&dimnames)[ndims],
    int* dimrng)
{

  int nc_err, ncid;
  int dimids[ndims];

  char filename_in[filename.size()+1];
  convert_string_to_char(filename,filename_in);
  nc_err = nc_create(filename_in, NC_CLOBBER, &ncid);
  scream_require_msg(nc_err==0,"NC Create Error! \n");

  char dimname_in[256];
  for (int i=0;i<ndims;i++) {
    convert_string_to_char(dimnames[i],dimname_in);
    nc_err = nc_def_dim(ncid,dimname_in, dimrng[i], &dimids[i]);
    scream_require_msg(nc_err==0,"NC Define Dim Error! \n");
  }

  return ncid;

}
/* ----------------------------------------------------------------- */
void init_output2_c(int ncid)
{

  int nc_err = nc_enddef(ncid);
  scream_require_msg(nc_err==0,"NC Error! Finish NC metadata \n");

}
/* ----------------------------------------------------------------- */
int init_output1(char** filename, int ndims, std::string (&dimnames)[ndims],
    int* dimrng)
{

  simple_out_init1(filename,ndims,dimnames,dimrng);

  return 0;
}
/* ----------------------------------------------------------------- */
void regfield_c(int ncid,std::string field_name,int field_type,int ndim,std::string (&field_dim)[ndim],std::string units)
{

  char field_in[field_name.size()+1];
  char dimname[256];
  char units_in[256];
  int dimids[ndim];
  int varid;
  int nc_err;

  // Extract dimension id dependent on the dimension names
  for (int ii=0;ii<ndim;ii++) {
    convert_string_to_char(field_dim[ii],dimname);
    nc_err = nc_inq_dimid(ncid,dimname,&dimids[ii]);
    scream_require_msg(nc_err==0,"NC Get dimid Error! \n");
  }
 
  // Register the field 
  convert_string_to_char(field_name,field_in);
  nc_err = nc_def_var(ncid,field_in,field_type,ndim,dimids,&varid);
  scream_require_msg(nc_err==0,"NC Var Create Error! \n");

  // Add the units attribute
  convert_string_to_char(units,units_in);
  nc_err = nc_put_att(ncid,varid,"units",NC_CHAR,units.size()+1,units_in);
  scream_require_msg(nc_err==0,"NC Var Units Error! \n");
}
/* ----------------------------------------------------------------- */
void field_io_c(int ncid,std::string field_name, double &data, int time_dim, const std::string flag)
{
  char fieldname[256];
  int varid;
  int nc_err;
  int ndims;

  convert_string_to_char(field_name,fieldname);

  nc_err = nc_inq_varid(ncid,fieldname,&varid);
  scream_require_msg(nc_err==0,"NC Var Inq ID Error! \n");
  nc_err = nc_inq_varndims(ncid,varid,&ndims);
  size_t start[ndims];
  size_t count[ndims];
  int dimids[ndims];
  nc_err = nc_inq_vardimid(ncid,varid,dimids);
  for (int ii=0;ii<ndims;ii++) {
    nc_err = nc_inq_dimlen(ncid,dimids[ii],&count[ii]);
    start[ii] = 0;
  }
  if (time_dim >= 0) {
    count[0] = 1;
    start[0] = time_dim;
  }

  if (flag=="read") {
    nc_err = nc_get_vara(ncid,varid, start, count, &data);
    scream_require_msg(nc_err==0,"NC Write Data Error! \n");
  } else if (flag=="write") {  
    nc_err = nc_put_vara(ncid,varid, start, count, &data);
    scream_require_msg(nc_err==0,"NC Write Data Error! \n");
  } else {
    scream_require_msg(0==1,"Error in IO, incorrect flag: " << flag);
  }

}
/* ----------------------------------------------------------------- */
void writefield_c(int ncid,std::string field_name, double &data, int time_dim)
{
  field_io_c(ncid, field_name, data, time_dim, "write");
}
/* ----------------------------------------------------------------- */
void readfield_c(int ncid,std::string field_name, double &data, int time_dim)
{
  field_io_c(ncid, field_name, data, time_dim, "read");
}
/* ----------------------------------------------------------------- */
void finalize_io_c(int ncid)
{

  int nc_err = nc_close(ncid);
  scream_require_msg(nc_err==0,"NC Error! Close NC file \n");

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
  int dlen_prod = 1;
  int len_chk = 0;
  for (int i=0;i<nspatial_dims;i++) {
    flen *=len[i];
    dlen_prod *= dlen[i];
    len_chk += len[i];
    len_chk -= dlen[i];
  }
  scream_require_msg(flen==dlen_prod,
                     "I/O Error! Inconsistency in field lengths for writefield, CHK1. \n");
  scream_require_msg(len_chk==0,
                     "I/O Error! Inconsistency in field lengths for writefield, CHK2. \n");

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
  int dlen_prod = 1;
  int len_chk = 0;
  for (int i=0;i<nspatial_dims;i++) {
    flen *=len[i];
    dlen_prod *= dlen[i];
    len_chk += len[i]; 
    len_chk -= dlen[i];
  }
  scream_require_msg(flen==dlen_prod,
                     "I/O Error! Inconsistency in field lengths for readfield, CHK1. \n");
  scream_require_msg(len_chk==0,
                     "I/O Error! Inconsistency in field lengths for readfield, CHK2. \n");
  
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
/* ----------------------------------------------------------------- */
int  simple_cpp_netcdf() {

  static const int NDIMS = 2;
  static const int NX = 6;
  static const int NY = 12;

  std::cout << "netcdf c is running!\n";

  int ncid, retval;
  retval = nc_create("simple_xy.nc", NC_CLOBBER, &ncid);

  int xDim, yDim;
  retval = nc_def_dim(ncid,"x", NX, &xDim);
  retval = nc_def_dim(ncid,"y", NY, &yDim);

  int data;
  int dimids[2] = {xDim, yDim};
  retval = nc_def_var(ncid,"data", NC_INT, 2, dimids, &data);

  retval = nc_enddef(ncid);

  int dataOut[NX][NY];
  for(int i = 0; i < NX; i++)
      for(int j = 0; j < NY; j++)
         dataOut[i][j] = i * NY + j;

  retval = nc_put_var_int(ncid,data,&dataOut[0][0]);
  retval = nc_close(ncid);

  return retval;

}

} // namespace simpleio
} // namespace scream
