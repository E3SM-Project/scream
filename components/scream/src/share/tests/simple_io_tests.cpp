#include <catch2/catch.hpp>

#include "share/scream_config.hpp"
#include "share/scream_pack.hpp"
#include "share/util/scream_utils.hpp"
#include "share/simpleio/simple_io_mod.hpp"


namespace {

TEST_CASE("simple_io_mod", "test_simple_io") {

  SECTION("Test Output") {
/* ====================================================================
 * Create sample output for simple_io test 
 *
 * Open netCDF output file and add all meta-data, including
 * All dimensions for any output field
 * Register all output fields
 * ==================================================================== */
  std::string filename = "pres_temp_4D.nc";
  std::string dimnames[] = {"time","longitude","latitude","level"};
  const int ndims = sizeof(dimnames)/sizeof(dimnames[0]);
  int dimrng[] = {0,12,6,2};
  int ncid = scream::simpleio::init_output1(filename, ndims, dimnames, dimrng);

  // Construct an array of fields
  scream::simpleio::nc_field *myfields = new scream::simpleio::nc_field[4];
  /* Latitude */
  myfields[0].fieldname = "latitude";
  myfields[0].ftype = scream::simpleio::NC_REAL;
  myfields[0].units = std::string("Degrees North");
  myfields[0].ndims = 1;
  myfields[0].dimensions[0] = {dimnames[2]}; 
  /* Longitude */
  myfields[1].fieldname = "longitude";
  myfields[1].ftype = scream::simpleio::NC_REAL;
  myfields[1].units = "Degrees East";
  myfields[1].ndims = 1;
  myfields[1].dimensions[0] = {dimnames[1]}; 
  /* Pressure */
  myfields[2].fieldname = "pressure";
  myfields[2].ftype = scream::simpleio::NC_REAL;
  myfields[2].units = "hPa";
  myfields[2].ndims = 4;
  for (int ii=0;ii<4;ii++) {myfields[2].dimensions[ii] = dimnames[ii];} 
  /* Temperature */
  myfields[3].fieldname = "temperature";
  myfields[3].ftype = scream::simpleio::NC_REAL;
  myfields[3].units = "K";
  myfields[3].ndims = 4;
  for (int ii=0;ii<4;ii++) {myfields[3].dimensions[ii] = dimnames[ii];} 
  /* Register Fields */
  for ( int ii=0;ii<4;ii++ ) {
    scream::simpleio::regfield(ncid,myfields[ii].fieldname,myfields[ii].ftype,myfields[ii].ndims,myfields[ii].dimensions,myfields[ii].units);
  }
  scream::simpleio::init_output2(ncid);

  // Now write some data
/* ====================================================================
 * Create sample output for simple_io test 
 * ==================================================================== */
  // Create sample output:
  double lats[6]={0,0,0,0,0,0};
  double lons[12];
  double pres[12][6][2], temp[12][6][2];
  for (int n=0; n<6; n++) {
    lats[n] = 25.0 + (n) * 5.0;
  }
  for (int n=0; n<12; n++) {
    lons[n] = -125.0 + (n) * 5.0;
  }
  for (int k=0; k<2;k++) {
    for (int i=0; i<6; i++) {
      for (int j=0; j<12; j++) {
        pres[j][i][k] = j*100.+i*10.+k;
        temp[j][i][k] = j+i+k;
      }
    }
  }

/* ====================================================================
 * Write sample output to netCDF output file.
 * ==================================================================== */
  scream::simpleio::writefield(ncid,myfields[0].fieldname,lats[0],-1);
  scream::simpleio::writefield(ncid,myfields[1].fieldname,lons[0],-1);
  scream::simpleio::writefield(ncid,myfields[2].fieldname,pres[0][0][0],0);
  scream::simpleio::writefield(ncid,myfields[3].fieldname,temp[0][0][0],0);
  scream::simpleio::writefield(ncid,myfields[2].fieldname,pres[0][0][0],1);
  scream::simpleio::writefield(ncid,myfields[3].fieldname,temp[0][0][0],1);


/* ====================================================================
 * Finished with test, close netCDF file
 * ==================================================================== */
  scream::simpleio::finalize_io(ncid);
  delete [] myfields;
  }

 /* ==================================================================== */
  SECTION("Test Input") {

  // Open the test file
  std::string filename = "pres_temp_4D.nc";
  int ncid_in;

  ncid_in = scream::simpleio::init_input(filename);
  std::string fieldname;

  fieldname  = "latitude";
  double lat[6];
  scream::simpleio::readfield(ncid_in,fieldname,lat[0],-1);
  std::cout << "LAT READ\n";
  for (double n : lat){
    std::cout << n << "\n";
  }
//
  fieldname = "pressure";
  double pres[12][6][2];
  scream::simpleio::readfield(ncid_in,fieldname,pres[0][0][0],0);
  std::cout << "PRES READ\n";
  for (int k=0;k<2;k++) {
    for (int j=0;j<6;j++) {
      for (int i=0;i<12;i++) {
        std::cout << pres[i][j][k] << "\t";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  scream::simpleio::finalize_io(ncid_in);
  } // Test Input

} // TEST_CASE
} // empty namespace
