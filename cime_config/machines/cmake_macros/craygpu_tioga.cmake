string(APPEND SLIBS " -L$ENV{PNETCDF_PATH}/lib -lpnetcdf")
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
#set(PIO_FILESYSTEM_HINTS "gpfs")
string(APPEND CXX_LIBS " -lstdc++")

SET(CMAKE_C_COMPILER "cc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "")
#SET(CMAKE_CXX_COMPILER "CC" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "hipcc" CACHE STRING "")

#string(APPEND CPPDEFS " -I${MPICH_DIR}/include")
#string(APPEND CXXFLAGS " -I${MPICH_DIR}/include ") #-I/opt/cray/pe/mpich/8.1.24/ucx/cray/10.0/include")
string(APPEND CXXFLAGS " -I/opt/cray/pe/mpich/8.1.24/ucx/cray/10.0/include")
#string(APPEND CFLAGS " -I${MPICH_DIR}/include")
#string(APPEND LDFLAGS " -L${MPICH_DIR}/lib -lmpi")
#ndk string(APPEND LDFLAGS " -lmpi_gtl_cuda")
#string(APPEND LDFLAGS " -L/opt/cray/pe/mpich/default/gtl/lib -lmpi_gtl_hsa")
#string(APPEND LDFLAGS " -L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.16/gtl/lib -lmpi_gtl_hsa")

# For YAKL's -lroctx64 -lrocfft; the rocm module doesn't set this.
string(APPEND LDFLAGS " -L$ENV{ROCM_PATH}/lib")

#string(APPEND CUDA_FLAGS " -ccbin CC -O2 -arch sm_80 --use_fast_math")
#string(APPEND CUDA_FLAGS " -ccbin CC")

#string(APPEND LDFLAGS " --gcc-toolchain=/opt/rh/gcc-toolset-10/root/usr")
#string(APPEND LDFLAGS " --gcc-toolchain=/opt/rh/gcc-toolset-12/root/usr")

#string(APPEND CFLAGS " -hnoacc -hfp0 -hipa0")
string(APPEND FFLAGS " -hnoacc -hfp0 -hipa0")
  
if (NOT DEBUG)
  #this resolves a crash in mct in docn init
  #string(APPEND CFLAGS " -O2 -hnoacc -hfp0 -hipa0")
  #string(APPEND FFLAGS " -O2 -hnoacc -hfp0 -hipa0")
  #string(APPEND CFLAGS " -hnoacc -hfp0 -hipa0")
  #string(APPEND FFLAGS " -hnoacc -hfp0 -hipa0")
  string(REPLACE ' -O3' ' ' CXX_FLAGS ${CMAKE_CXX_FLAGS})
  string(REPLACE ' -O3' ' ' CXXFLAGS ${CXXFLAGS})
  string(REPLACE ' -O2' ' ' CXX_FLAGS ${CMAKE_CXX_FLAGS})
  string(REPLACE ' -O2' ' ' CXXFLAGS ${CXXFLAGS})
  string(REPLACE ' -O3' ' ' FFLAGS ${FFLAGS})
  string(REPLACE ' -O2' ' ' FFLAGS ${FFLAGS})
endif()


#3: /home/jenkins/crayftn/craylibs/google-perftools/src/tcmalloc.cc:647] Attempt to free invalid pointer: 0xa15d810
#This is a signature for a known problem with tcmalloc, that the CCE Fortran compiler uses by default for the releases currently on Tioga.
#They next craype release should make tcmalloc non-default, I believe.
#In the meantime, could you rebuild adding "-hsystem_alloc" to the fortran compiler flags? (and link flags!)
string(APPEND FFLAGS " -hsystem_alloc")
string(APPEND LDFLAGS " -hsystem_alloc")
 
if (DEBUG)
  #string(APPEND FFLAGS " --save-temps")
  #string(APPEND FFLAGS " -mllvm -print-after-all") # lots of output
  # commented out to try -g with cce15 ?
  string(REPLACE ' -g' ' ' CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  string(REPLACE ' -g' ' ' CXXFLAGS ${CXXFLAGS})
endif()

set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")


#set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "See SCREAM issue #2080.")
