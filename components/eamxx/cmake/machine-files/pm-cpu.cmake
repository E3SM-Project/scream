include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/amd-zen3.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

#message(STATUS "pm-cpu CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()

set(PYTHON_EXECUTABLE "/opt/cray/pe/python/3.9.13.1/bin/python3" CACHE STRING "" FORCE)
set(PYTHON_LIBRARIES "/opt/cray/pe/python/3.9.13.1/lib/libpython3.9.so.1.0" CACHE STRING "" FORCE)
option (SCREAM_ENABLE_ML_CORRECTION "Whether to enable ML correction parametrization" ON)
set(HDF5_DISABLE_VERSION_CHECK 1 CACHE STRING "" FORCE)
execute_process(COMMAND source /global/cfs/cdirs/e3sm/eamxx-ml/python_venv/3.9.13/screamML/bin/activate)