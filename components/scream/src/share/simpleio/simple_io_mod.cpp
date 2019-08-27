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
  //void simple_out_writefield(char** field_name, Real* field_data);
  void io_writefield_1d_real(char** field_name, Int fsize, Real* field_data);

} // extern C

namespace scream {
namespace simpleio {
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
int regfield(char** field_name,int field_type,int ndim,std::string (&field_dim)[ndim],char** units)
{
  simple_out_regfield(field_name,field_type,ndim,field_dim,units);

  return 0;
}
/* ----------------------------------------------------------------- */


} // namespace simpleio
} // namespace scream
