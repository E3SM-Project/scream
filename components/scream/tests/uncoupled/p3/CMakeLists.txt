INCLUDE (ScreamUtils)

SET (NEED_LIBS p3 scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(p3_standalone "p3_standalone.cpp" "${NEED_LIBS}" LABELS "p3;physics"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_np${MPI_RANKS}.yaml)
  CONFIGURE_FILE(p3_standalone_output.yaml p3_standalone_output_np${MPI_RANKS}.yaml)
  
endforeach()

# Copy p3 lookup tables to local data directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_universal_constants.inp
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(p3_init_ne2np4.nc p3_init_ne2np4.nc COPYONLY)

## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if max mpi ranks is >1
# This test requires CPRNC
include (BuildCprnc)
BuildCprnc()
SET (BASE_TEST_NAME "p3")
foreach (MPI_RANKS RANGE 2 ${SCREAM_TEST_MAX_RANKS})
  set (SRC_FILE "${BASE_TEST_NAME}_standalone_output_np${MPI_RANKS}.INSTANT.Steps_x5.0000-01-01.002500.nc")
  set (TGT_FILE "${BASE_TEST_NAME}_standalone_output_np1.INSTANT.Steps_x5.0000-01-01.002500.nc")
  configure_file (${SCREAM_BASE_DIR}/cmake/CprncTest.cmake
                  ${CMAKE_CURRENT_BINARY_DIR}/CprncTest_np${MPI_RANKS}.cmake @ONLY)
  add_test (NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb"
            COMMAND cmake -P CprncTest_${MPI_RANKS}.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

endforeach()
