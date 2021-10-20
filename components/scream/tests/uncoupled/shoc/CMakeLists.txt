include (ScreamUtils)

# Create the test
set (NEED_LIBS shoc scream_control scream_share)
CreateUnitTest(shoc_standalone "shoc_standalone.cpp" "${NEED_LIBS}" LABELS "shoc;physics"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

# Set AD configurable options
SetVarDependingOnTestProfile(NUM_STEPS 2 24 48)
set (ATM_TIME_STEP 1800)

# Copy input files to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shoc_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_init_ne2np4.nc COPYONLY)

## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_np${MPI_RANKS}.yaml)
  CONFIGURE_FILE(shoc_standalone_output.yaml shoc_standalone_output_np${MPI_RANKS}.yaml)
  
endforeach()
