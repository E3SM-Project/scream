#include <catch2/catch.hpp>

#include "share/scream_config.hpp"
#include "share/scream_pack.hpp"
#include "share/util/scream_utils.hpp"
#include "share/simpleio/simple_io_mod.hpp"


namespace {

TEST_CASE("simple_out_mod", "test_simple_out") {

  SECTION("Hello World") {
/* TODO: Delete
 * Simple hello_world to test input/output and c-binding
 */
  int myint = 1;
  double myreal = 3.14;
  char* mychar = "Aaron";
  std::string planets[] = {"Mercury", "Venus", "Earth", "Mars", "Saturn", "Jupiter", "Neptune", "Uranus"};
  const int len_planet = sizeof(planets)/sizeof(planets[0]);

  // Simple Hello World - TODO delete
  int nerr = scream::simpleio::hello_world(myint, myreal, &mychar,len_planet,planets);
  REQUIRE (nerr==0);
  } // Hello world
  SECTION("Test Output") {
/* ====================================================================
 * Create sample output for simple_io test 
 * ==================================================================== */
  // Create sample output:
  double lats[6]={0,0,0,0,0,0};
  double lons[12];
  double pres[2][6][12], temp[2][6][12];
  for (int n=0; n<6; n++) {
    lats[n] = 25.0 + (n) * 5.0;
  }
  for (int n=0; n<12; n++) {
    lons[n] = -125.0 + (n) * 5.0;
  }
  int n = 0;
  for (int k=0; k<2;k++) {
    for (int i=0; i<6; i++) {
      for (int j=0; j<12; j++) {
        pres[k][i][j] = 900.0 + float(n);
        temp[k][i][j] = 9.0   + float(n);
        n = n+1;
      }
    }
  }

/* ====================================================================
 * Open netCDF output file and add all meta-data, including
 * All dimensions for any output field
 * Register all output fields
 * ==================================================================== */
  // Variables for creating the output netcdf file:
  char* filename = "pres_temp_4D.nc";
  std::string dimnames[] = {"longitude","latitude","level","time"};
  const int ndims = sizeof(dimnames)/sizeof(dimnames[0]);
  int dimrng[] = {12,6,2,0};
  int dim3d[]  = {12,6,2};

  // Create the netcdf file
  int init1_err = scream::simpleio::init_output1(&filename, ndims, dimnames, dimrng);
  REQUIRE (init1_err==0);
  // Register all of the fields in the file
  int reg = 0;
  char *fieldname, *units;
  int ftype, reg1;
  std::string vec_lat[] = {"latitude"};
  std::string vec_lon[] = {"longitude"};
  std::string vec_3d_state[] = {"longitude","latitude","level","time"};
  int ndim_1d = 1;
  int ndim_3dtime = 4;

  // LAT
  fieldname = "latitude";
  units = "degrees_north";
  ftype = scream::simpleio::NC_REAL;
  reg1 = scream::simpleio::regfield(&fieldname,ftype,ndim_1d,vec_lat,&units);
  reg = std::max(reg,reg1);
  // LON
  fieldname = "longitude";
  units = "degrees_east";
  ftype = scream::simpleio::NC_REAL;
  reg1 = scream::simpleio::regfield(&fieldname,ftype,ndim_1d,vec_lon,&units);
  reg = std::max(reg,reg1);
  // PRESSURE
  fieldname = "pressure";
  units = "hPa";
  ftype = scream::simpleio::NC_REAL;
  reg1 = scream::simpleio::regfield(&fieldname,ftype,ndim_3dtime,vec_3d_state,&units);
  reg = std::max(reg,reg1);
  // TEMP
  fieldname = "temperature";
  units = "celcius";
  reg1 = scream::simpleio::regfield(&fieldname,ftype,ndim_3dtime,vec_3d_state,&units);
  reg = std::max(reg,reg1);
  // Finished with fields
  REQUIRE(reg==0);

  // Finish netCDF metadata:
  int init2_err = scream::simpleio::init_output2();
  REQUIRE (init2_err==0);

/* ====================================================================
 * Write sample output to netCDF output file.
 * ==================================================================== */
 int wrt=0, wrt1;
 int timelev = 1;
 fieldname = "latitude";
 wrt1 = scream::simpleio::writefield(&fieldname,timelev,dimrng[1],lats);
 wrt = std::max(wrt,wrt1);
 fieldname = "longitude";
 wrt1 = scream::simpleio::writefield(&fieldname,timelev,dimrng[0],lons);
 wrt = std::max(wrt,wrt1);
 fieldname = "pressure";
 wrt1 = scream::simpleio::writefield(&fieldname,timelev,dim3d,pres);
 wrt = std::max(wrt,wrt1);
 fieldname = "temperature";
 wrt1 = scream::simpleio::writefield(&fieldname,timelev,dim3d,temp);
 wrt = std::max(wrt,wrt1);
 timelev++;
 fieldname = "pressure";
 wrt1 = scream::simpleio::writefield(&fieldname,timelev,dim3d,pres);
 wrt = std::max(wrt,wrt1);
 fieldname = "temperature";
 wrt1 = scream::simpleio::writefield(&fieldname,timelev,dim3d,temp);
 wrt = std::max(wrt,wrt1);

 REQUIRE(wrt==0);

/* ====================================================================
 * Finished with test, close netCDF file
 * ==================================================================== */
  // All done, not close the file
  int final_err = scream::simpleio::finalize_output();
  REQUIRE (final_err==0);
  } // Test output
/* ====================================================================
 * Test Input
 * ==================================================================== */
  SECTION("Test Input") {

  // Open the test file
  typedef std::vector<double> dbl_vector;
  typedef std::vector<dbl_vector> dbl_2dmatrix;
  typedef std::vector<dbl_2dmatrix> dbl_3dmatrix;
  char* filename = "pres_temp_4D.nc";
  char* fieldname;
  int rd;

  rd = scream::simpleio::init_input(&filename);

  int time_dim = 1;

  fieldname  = "latitude";
  dbl_vector lat;
  rd = scream::simpleio::readfield(&fieldname,time_dim,lat);
  for (double n : lat){
    std::cout << n << "\n";
  }

  fieldname = "pressure";
  dbl_3dmatrix pres;
  rd = scream::simpleio::readfield(&fieldname,time_dim,pres);
  for (dbl_2dmatrix n1 : pres){
    for (dbl_vector n2 : n1) {
      for (double n3 : n2) {
        std::cout << n3 << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

  rd = scream::simpleio::finalize_input();
  REQUIRE(rd==0);

  } // Test Input

} // TEST_CASE
} // empty namespace
