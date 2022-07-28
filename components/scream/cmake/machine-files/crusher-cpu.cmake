# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MACHINE "crusher-cpu" CACHE STRING "")

#SET(MPICH_DIR "/opt/cray/pe/mpich/8.1.12/ofi/crayclang/10.0" CACHE STRING "")

#set(CMAKE_CXX_FLAGS "-I${MPICH_DIR}/include" CACHE STRING "" FORCE)
#set(CMAKE_C_FLAGS "-I${MPICH_DIR}/include" CACHE STRING "" FORCE)
#set(CMAKE_Fortran_FLAGS "-I${MPICH_DIR}/include" CACHE STRING "" FORCE)
#set(CMAKE_EXE_LINKER_FLAGS "-L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.12/gtl/lib -lmpi_gtl_hsa" CACHE STRING "" FORCE)

















