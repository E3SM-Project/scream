
################ these are CPU settings

# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")

#not clear how SCREAM_MACHINE is used. Sometimes it is ENV variable.
set(SCREAM_MACHINE "crusher" CACHE STRING "")


################ these are GPU settings
#this is in cache, e3sm var
#COMPILER:UNINITIALIZED=crayclanggpu
if(${COMPILER} STREQUAL "crayclanggpu")

include (${EKAT_MACH_FILES_PATH}/kokkos/mi250.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/hip.cmake)

endif()


