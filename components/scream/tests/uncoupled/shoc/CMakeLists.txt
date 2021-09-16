INCLUDE (ScreamUtils)

SET (NEED_LIBS shoc scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(shoc_standalone "shoc_standalone.cpp" "${NEED_LIBS}" LABELS "shoc;physics"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/shoc_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_init_ne2np4.nc COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/shoc_standalone_output.yaml 
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_standalone_output.yaml COPYONLY)

## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_np${MPI_RANKS}.yaml)
  CONFIGURE_FILE(shoc_standalone_output.yaml shoc_standalone_output_np${MPI_RANKS}.yaml)
  
endforeach()
