# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# Remove this if you are using a resource manager (slurm etc)
set (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES True CACHE BOOL "")
