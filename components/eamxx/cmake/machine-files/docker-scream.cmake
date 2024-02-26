include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# include (${EKAT_MACH_FILES_PATH}/kokkos/amd-zen3.cmake)
# include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-a100.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)

set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

# set(BLAS_LIBRARIES /opt/conda/lib/libblas.so CACHE STRING "")
# set(LAPACK_LIBRARIES /opt/conda/lib/liblapack.so CACHE STRING "")
# set(SCREAM_INPUT_ROOT "/storage/inputdata/" CACHE STRING "")
set(PYBIND11_PYTHON_VERSION 3.9 CACHE STRING "")
set(Python_EXECUTABLE /global/homes/a/andrep/.conda/envs/fv3net/bin/python CACHE STRING "" FORCE)
set(PYTHON_EXECUTABLE /global/homes/a/andrep/.conda/envs/fv3net/bin/python CACHE STRING "" FORCE)
set(Python_INCLUDE_DIRS /global/homes/a/andrep/.conda/envs/fv3net/include/python3.9 CACHE STRING "" FORCE)
option (SCREAM_ENABLE_ML_CORRECTION "Whether to enable ML correction parametrization" ON)

if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()
