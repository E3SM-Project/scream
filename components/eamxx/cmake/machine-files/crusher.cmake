
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

#set(SCREAM_MACHINE "crusher-gpu" CACHE STRING "")

#is this needed here?
SET(MPICH_DIR "/opt/cray/pe/mpich/8.1.16/ofi/crayclang/10.0" CACHE STRING "")

set(CMAKE_CXX_FLAGS "--amdgpu-target=gfx90a -fno-gpu-rdc  -I${MPICH_DIR}/include" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "-I${MPICH_DIR}/include" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "-I${MPICH_DIR}/include" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.16/gtl/lib -lmpi_gtl_hsa" CACHE STRING "" FORCE)

endif()


