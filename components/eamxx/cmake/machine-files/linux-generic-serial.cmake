# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/generic.cmake)

# Additional settings
set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION TRUE CACHE BOOL "")
set(Kokkos_ENABLE_OPENMP FALSE CACHE BOOL "")

# Remove this if you are using a resource manager (slurm etc)
set (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES True CACHE BOOL "")
