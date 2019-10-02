#include "simple_io_mod.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"
#include <netcdf.h>

using scream::Real;
using scream::Int;

namespace scream {
namespace simpleio {

/* ----------------------------------------------------------------- */
void convert_string_to_char(std::string str_in,char str_out[str_in.size()+1]) {

  str_in.copy(str_out,str_in.size()+1);
  str_out[str_in.size()] = '\0';

}
/* ----------------------------------------------------------------- */
int init_input(std::string filename) {


  int nc_err, ncid;
  char filename_in[filename.size()+1];

  convert_string_to_char(filename,filename_in);
  nc_err = nc_open(filename_in, NC_NOWRITE, &ncid);
  scream_require_msg(nc_err==0,"NC Create Error! \n");

  return ncid;
}
/* ----------------------------------------------------------------- */
int init_output1(std::string filename, int ndims, std::string (&dimnames)[ndims],
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
void init_output2(int ncid)
{

  int nc_err = nc_enddef(ncid);
  scream_require_msg(nc_err==0,"NC Error! Finish NC metadata \n");

}
/* ----------------------------------------------------------------- */
void regfield(int ncid,std::string field_name,int field_type,int ndim,std::string (&field_dim)[ndim],std::string units)
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
void field_io(int ncid,std::string field_name, double &data, int time_dim, const std::string flag)
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
void writefield(int ncid,std::string field_name, double &data, int time_dim)
{
  field_io(ncid, field_name, data, time_dim, "write");
}
/* ----------------------------------------------------------------- */
void readfield(int ncid,std::string field_name, double &data, int time_dim)
{
  field_io(ncid, field_name, data, time_dim, "read");
}
/* ----------------------------------------------------------------- */
void finalize_io(int ncid)
{

  int nc_err = nc_close(ncid);
  scream_require_msg(nc_err==0,"NC Error! Close NC file \n");

}


} // namespace simpleio
} // namespace scream
