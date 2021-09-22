INCLUDE (ScreamUtils)

SET (NEED_LIBS spa scream_control scream_share)

# TODO - SPA_DATA_PATH should be set depending on which machine is doing the build
#        We will have to do this before this can be merged to master.
SET (SPA_DATA_PATH "/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/chem/spa")

# Test atmosphere processes
CreateUnitTest(spa_stand_alone "spa_stand_alone.cpp" "${NEED_LIBS}" LABELS "spa;physics")

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/spa_init_ne2np4.nc COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/spa_output.yaml COPYONLY)
