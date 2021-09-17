INCLUDE (ScreamUtils)

# Test atmosphere processes
if (NOT ${SCREAM_BASELINES_ONLY})

    # Required libraries
    find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATHS}/lib)
    set (NEED_LIBS scream_rrtmgp rrtmgp ${NETCDF_C} scream_control scream_share physics_share yakl)

    ## Create the unit test for rrtmgp that compares against a rrtmgp-allsky data file.
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

    ## Create free running rrtmgp stand alone test that runs from an initial condition file.
    CreateUnitTest(
        rrtmgp_standalone "rrtmgp_standalone.cpp" "${NEED_LIBS}" LABELS "rrtmgp;physics"
        MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    )

    # Copy yaml input file to run directory
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input_unit.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/input_unit.yaml COPYONLY)
    ## Copy (and configure) yaml files needed by tests
    foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
      
      CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                     ${CMAKE_CURRENT_BINARY_DIR}/input_np${MPI_RANKS}.yaml)
      CONFIGURE_FILE(rrtmgp_standalone_output.yaml rrtmgp_standalone_output_np${MPI_RANKS}.yaml)
      
    endforeach()

    # Copy RRTMGP initial condition to local data directory
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/rrtmgp_init_ne2np4.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/rrtmgp_init_ne2np4.nc COPYONLY)

    FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
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
endif()
