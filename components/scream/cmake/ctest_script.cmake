cmake_minimum_required(VERSION 3.9)

set(CTEST_BUILD_NAME "scream_unit_tests_${BUILD_NAME_MOD}")

get_filename_component(SCREAM_ROOT ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)
set(CTEST_SOURCE_DIRECTORY "${SCREAM_ROOT}")
set(CTEST_BINARY_DIRECTORY "${BUILD_WORK_DIR}")

if(NOT DEFINED dashboard_model)
  set(dashboard_model Experimental)
endif()
if(NOT DEFINED dashboard_track)
  set(dashboard_track E3SM_SCREAM)
endif()

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

ctest_start(${dashboard_model} TRACK ${dashboard_track})

if (NOT SKIP_CONFIG_STEP)
  # Config step is not to be skipped
  separate_arguments(OPTIONS_LIST UNIX_COMMAND "${CMAKE_COMMAND}")
  ctest_configure(OPTIONS "${OPTIONS_LIST}" RETURN_VALUE CONFIG_ERROR_CODE)
endif()

if (CONFIG_ERROR_CODE)
  set (TEST_FAILS TRUE)
else ()
  if (DEFINED ENV{SCREAM_BUILD_PARALLEL_LEVEL})
    ctest_build(FLAGS "-j$ENV{SCREAM_BUILD_PARALLEL_LEVEL}" RETURN_VALUE BUILD_ERROR_CODE)
  else()
    ctest_build(FLAGS "-j4" RETURN_VALUE BUILD_ERROR_CODE)
  endif()

  # Need this code so that build errors don't get buried
  if (BUILD_ERROR_CODE)
    set(TEST_FAILS TRUE)
    file(GLOB MATCHES "${CTEST_BINARY_DIRECTORY}/Testing/Temporary/LastBuild*.log")
      if (MATCHES)
        foreach (MATCH IN LISTS MATCHES)
          file(READ ${MATCH} BUILD_OUTPUT)
          message("Build failed with output:")
          message("${BUILD_OUTPUT}")
        endforeach()
      endif()
  else()
    if (NOT BUILD_ONLY)

      if (DEFINED ENV{SCREAM_FORCE_RUN_DIFF})
        set(INCLUDE_REGEX "p3_run_and_cmp")
      elseif (DEFINED ENV{SCREAM_FORCE_RUN_DIFF_BFB_UNIT})
        set(INCLUDE_REGEX "p3_tests_ut")
      else()
        set(INCLUDE_REGEX ".") # all
      endif()

      if (DEFINED ENV{CTEST_PARALLEL_LEVEL})
        ctest_test(PARALLEL_LEVEL $ENV{CTEST_PARALLEL_LEVEL} RETURN_VALUE TEST_ERROR_CODE INCLUDE ${INCLUDE_REGEX})
      else()
        ctest_test(RETURN_VALUE TEST_ERROR_CODE INCLUDE ${INCLUDE_REGEX})
      endif()
      if (TEST_ERROR_CODE)
        set(TEST_FAILS TRUE)
      endif()
      if (DO_COVERAGE)
        ctest_coverage()
      endif()
    endif()
  endif()
endif ()

if (NOT NO_SUBMIT)
  ctest_submit(RETRY_COUNT 10 RETRY_DELAY 60)
endif()

if (TEST_FAILS)
  message(FATAL_ERROR "Test had fails")
endif()
