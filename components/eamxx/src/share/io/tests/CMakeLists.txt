#include(ScreamUtils)
include(ScreamUtils)

# The restart test requires cprnc
include (BuildCprnc)
BuildCprnc()

## Test basic output (no packs, no diags, all avg types, all freq units)
CreateUnitTest(io_basic "io_basic.cpp" "scream_io" LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test packed I/O
CreateUnitTest(io_packed "io_packed.cpp" "scream_io" LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test diagnostic output
CreateUnitTest(io_diags "io_diags.cpp" "scream_io" LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

# Test output on SE grid
configure_file(io_test_se_grid.yaml io_test_se_grid.yaml)
CreateUnitTest(io_test_se_grid "io_se_grid.cpp" scream_io LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test output restart
# NOTE: Each restart test is a "setup" for the restart_check test,
# and cannot run in parallel with other restart tests,
# due to contention of the rpointer file
foreach (AVG_TYPE IN ITEMS INSTANT AVERAGE)
  configure_file(io_test_restart.yaml io_test_restart_${AVG_TYPE}.yaml)
  CreateUnitTest(output_restart_${AVG_TYPE}_test "output_restart.cpp" scream_io LABELS "io"
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    PROPERTIES RESOURCE_LOCK rpointer_file FIXTURES_SETUP restart_setup_${AVG_TYPE}
    EXEC_ARGS --ekat-test-params ifile=io_test_restart_${AVG_TYPE}.yaml
  )
  foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
    set (SRC_FILE monolithic.${AVG_TYPE}.nsteps_x10.np${MPI_RANKS}.2000-01-01-00000.nc)
    set (TGT_FILE restarted.${AVG_TYPE}.nsteps_x10.np${MPI_RANKS}.2000-01-01-00000.nc)
    add_test (NAME io_test_restart_check_np${MPI_RANKS}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_property(TEST io_test_restart_check_np${MPI_RANKS}
                 PROPERTY FIXTURES_REQUIRED restart_setup_${AVG_TYPE})
  endforeach()
endforeach()

## Test remap output
CreateUnitTest(io_remap_test "io_remap_test.cpp" "scream_io;diagnostics" LABELS "io,remap"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

