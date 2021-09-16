INCLUDE (ScreamUtils)

SET (NEED_LIBS spa scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(spa_stand_alone "spa_stand_alone.cpp" "${NEED_LIBS}" LABELS "spa;physics"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_np${MPI_RANKS}.yaml)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_standalone_output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/spa_standalone_output_np${MPI_RANKS}.yaml)

endforeach()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/spa_init_ne2np4.nc COPYONLY)
