INCLUDE (ScreamUtils)

# Test atmosphere processes
if (NOT ${SCREAM_BASELINES_ONLY})

    # Required libraries
    find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATHS}/lib)
    set (NEED_LIBS scream_rrtmgp rrtmgp ${NETCDF_C} scream_control scream_share physics_share yakl)

    CreateUnitTest(
        rrtmgp_standalone_unit rrtmgp_standalone_unit.cpp "${NEED_LIBS}" LABELS "rrtmgp;physics"
        EXE_ARGS "--ekat-test-params rrtmgp_inputfile=${CMAKE_CURRENT_BINARY_DIR}/data/rrtmgp-allsky.nc,rrtmgp_baseline=${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc"
    )
    target_include_directories(rrtmgp_standalone_unit PUBLIC
          ${SCREAM_BASE_DIR}/src/physics/rrtmgp
          ${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/kernels
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/extensions/cloud_optics
    )

    # Copy yaml input file to run directory
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

    # Copy RRTMGP initial condition to local data directory
    file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)

    # Copy RRTMGP absorption coefficient lookup tables to local data directory
    configure_file(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
    configure_file(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
    configure_file(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/examples/all-sky/rrtmgp-allsky.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
    configure_file(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
    configure_file(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)

    # The RRTMGP stand-alone test that runs multi-step
    # Set AD configurable options
    SetVarDependingOnTestProfile(NUM_STEPS 2 24 48)
    set (ATM_TIME_STEP 1800)

    CreateUnitTest(rrtmgp_standalone "rrtmgp_standalone.cpp" "${NEED_LIBS}" LABELS "rrtmgp;physics"
      MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    )
    ## Copy (and configure) yaml files needed by tests
    foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
      
      CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                     ${CMAKE_CURRENT_BINARY_DIR}/input_np${MPI_RANKS}.yaml)
      CONFIGURE_FILE(rrtmgp_standalone_output.yaml rrtmgp_standalone_output_np${MPI_RANKS}.yaml)
      
    endforeach()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rrtmgp_init_ne2np4.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/rrtmgp_init_ne2np4.nc COPYONLY)

  ## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
  ## Only if max mpi ranks is >1
  # This test requires CPRNC
  include (BuildCprnc)
  BuildCprnc()
  SET (BASE_TEST_NAME "rrtmgp")
  foreach (MPI_RANKS RANGE 2 ${SCREAM_TEST_MAX_RANKS})
    set (SRC_FILE "${BASE_TEST_NAME}_standalone_output_np${MPI_RANKS}.INSTANT.Steps_x5.0000-01-01.002500.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_standalone_output_np1.INSTANT.Steps_x5.0000-01-01.002500.nc")
    configure_file (${SCREAM_BASE_DIR}/cmake/CprncTest.cmake
                    ${CMAKE_CURRENT_BINARY_DIR}/CprncTest_np${MPI_RANKS}.cmake @ONLY)
    add_test (NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb"
              COMMAND cmake -P CprncTest_${MPI_RANKS}.cmake
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  
  endforeach()
endif()
