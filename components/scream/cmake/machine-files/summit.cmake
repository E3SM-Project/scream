# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-v100.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)

set(SCREAM_MPIRUN_EXE "jsrun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-n" CACHE STRING "")
