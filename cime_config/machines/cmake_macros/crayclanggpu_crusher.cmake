#the only file that will be loaded for compiler crayclanggpu and machine crusher

if (compile_threaded)
  string(APPEND FFLAGS   " -fopenmp")
  string(APPEND CFLAGS   " -fopenmp")
  string(APPEND CXXFLAGS " -fopenmp")
  string(APPEND LDFLAGS  " -fopenmp")
endif()
if (DEBUG)
  string(APPEND CFLAGS   " -O0 -g")
  string(APPEND FFLAGS   " -O0 -g")
  string(APPEND CXXFLAGS " -O0 -g")
  string(APPEND CPPDEFS " -DYAKL_DEBUG")
endif()

string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRCRAY")
string(APPEND FC_AUTO_R8 " -s real64")
string(APPEND FFLAGS " -f free -N 255 -h byteswapio -em")

if (NOT compile_threaded)
	string(APPEND FFLAGS " -M1077")
endif()

string(APPEND FFLAGS_NOOPT " -O0")
set(HAS_F2008_CONTIGUOUS "TRUE")
string(APPEND LDFLAGS " -Wl,--allow-multiple-definition")

set(SUPPORTS_CXX "TRUE")
set(CXX_LINKER "FORTRAN")

#this differs from scream CMAKE_CXX_COMPILER...
set(MPICC "cc")
set(MPICXX "hipcc")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "hipcc")
set(SFC "ftn")

if (NOT DEBUG)
  string(APPEND CFLAGS   " -O2")
  string(APPEND CXXFLAGS " -O2")
  string(APPEND FFLAGS   " -O2")
endif()

if (COMP_NAME STREQUAL elm)
  string(APPEND FFLAGS " -hfp0")
endif()
string(APPEND FFLAGS " -hipa0 -hnoacc")

string(APPEND SLIBS " -L$ENV{PNETCDF_PATH}/lib -lpnetcdf")
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
set(PIO_FILESYSTEM_HINTS "gpfs")
string(APPEND CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")

#this is in config_machines
#set(PIO_FILESYSTEM_HINTS "romio_cb_read=disable")

string(APPEND CXX_LIBS " -lstdc++")

string(APPEND CXXFLAGS " -I$ENV{MPICH_DIR}/include")
string(APPEND SLIBS    " -L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa")
string(APPEND SLIBS    " -L$ENV{ROCM_PATH}/lib")

set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "See SCREAM issue #2080.")

#scream builds its kokkos from ekat
#string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On")

