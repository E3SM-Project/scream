include (ScreamUtils)

<<<<<<< HEAD
# Create the test
set (NEED_LIBS p3 scream_control scream_share)
CreateUnitTest(p3_standalone "p3_standalone.cpp" "${NEED_LIBS}" LABELS "p3;physics"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestProfile(NUM_STEPS 2 8 48)  # 1h 4h 24h

## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_np${MPI_RANKS}.yaml)
  CONFIGURE_FILE(p3_standalone_output.yaml p3_standalone_output_np${MPI_RANKS}.yaml)
  
endforeach()

# Copy p3 lookup tables to local data directory
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
configure_file(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${SCREAM_DATA_DIR}/p3_universal_constants.inp
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(p3_init_ne2np4.nc p3_init_ne2np4.nc COPYONLY)
