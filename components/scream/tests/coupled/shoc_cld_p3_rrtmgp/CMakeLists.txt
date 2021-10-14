INCLUDE (ScreamUtils)

# Create the test
SET (NEED_LIBS shoc cld_fraction p3 scream_rrtmgp scream_control scream_share physics_share)
CreateUnitTest(shoc_cld_p3_rrtmgp shoc_cld_p3_rrtmgp.cpp "${NEED_LIBS}" LABELS "shoc;cld;p3;rrtmgp;physics"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestProfile(NUM_STEPS 2 8 48)  # 1h 4h 24h

## Copy input files to run directory
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
set (RRTMGP_EXT_DIR ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shoc_cld_p3_rrtmgp_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_cld_p3_rrtmgp_init_ne2np4.nc COPYONLY)

configure_file(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)

configure_file(${RRTMGP_EXT_DIR}/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${RRTMGP_EXT_DIR}/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${RRTMGP_EXT_DIR}/examples/all-sky/rrtmgp-allsky.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${RRTMGP_EXT_DIR}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${RRTMGP_EXT_DIR}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_np${MPI_RANKS}.yaml)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/shoc_cld_p3_rrtmgp_standalone_output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/shoc_cld_p3_rrtmgp_standalone_output_np${MPI_RANKS}.yaml)
  
endforeach()
## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if max mpi ranks is >1
# This test requires CPRNC
include (BuildCprnc)
BuildCprnc()
SET (BASE_TEST_NAME "shoc_cld_p3_rrtmgp")
foreach (MPI_RANKS RANGE 2 ${SCREAM_TEST_MAX_RANKS})
  set (SRC_FILE "${BASE_TEST_NAME}_standalone_output_np${MPI_RANKS}.INSTANT.Steps_x5.0000-01-01.002500.nc")
  set (TGT_FILE "${BASE_TEST_NAME}_standalone_output_np1.INSTANT.Steps_x5.0000-01-01.002500.nc")
  configure_file (${SCREAM_BASE_DIR}/cmake/CprncTest.cmake
                  ${CMAKE_CURRENT_BINARY_DIR}/CprncTest_np${MPI_RANKS}.cmake @ONLY)
  add_test (NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb"
            COMMAND cmake -P CprncTest_${MPI_RANKS}.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

endforeach()
