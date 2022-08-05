# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  include(ScreamUtils)

  # Test utils
  CreateUnitTest(utils "utils_tests.cpp" scream_share)

  # Test column ops
  CreateUnitTest(column_ops "column_ops.cpp" scream_share)

  # Test fields
  CreateUnitTest(field "field_tests.cpp" scream_share)

  # Test field utils
  CreateUnitTest(field_utils "field_utils.cpp" scream_share
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test property checks
  CreateUnitTest(property_checks "property_checks.cpp" scream_share)

  # Test grids
  CreateUnitTest(grid "grid_tests.cpp" scream_share
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test common physics functions
  CreateUnitTest(common_physics "common_physics_functions_tests.cpp" scream_share)

  # Test atmosphere processes
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/atm_process_tests_parse_list.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/atm_process_tests_parse_list.yaml COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/atm_process_tests_named_procs.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/atm_process_tests_named_procs.yaml COPYONLY)
  CreateUnitTest(atm_proc "atm_process_tests.cpp" scream_share)

  # Test horizontal remapping utility
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/horizontal_remap_baseline.nc
                 ${CMAKE_CURRENT_BINARY_DIR}/horizontal_remap_baseline.nc COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/horiz_remap_from_file_control.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/horiz_remap_from_file_control.yaml)
  CreateUnitTest(horizontal_remap "horizontal_remap_test.cpp" "scream_share;scream_io"
    LABELS "horiz_remap"
    PROPERTIES FIXTURES_SETUP horiz_remap
  )
  include (BuildCprnc)
  BuildCprnc()
  add_test (NAME "verify_horizontal_remap"
            COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake horizontal_remap_output.nc horizontal_remap_baseline.nc
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties("verify_horizontal_remap" PROPERTIES LABELS "horiz_remap"
            RESOURCE_LOCK "verify_horizontal_remap"
            FIXTURES_REQUIRED horiz_remap)
 

endif()
