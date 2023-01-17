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

#do we need -DFORTRANUNDERSCORE -DNO_R16 ?
string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRCRAY")

#do we need this?
string(APPEND FC_AUTO_R8 " -s real64")

#do we need this?
string(APPEND FFLAGS " -f free -N 255 -h byteswapio -em")

#do we need this?
if (NOT compile_threaded)
	string(APPEND FFLAGS " -M1077")
endif()

string(APPEND FFLAGS_NOOPT " -O0")

#do we need this?
set(HAS_F2008_CONTIGUOUS "TRUE")

#do we need this?
string(APPEND LDFLAGS " -Wl,--allow-multiple-definition")

#do we need these two?
set(SUPPORTS_CXX "TRUE")
set(CXX_LINKER "FORTRAN")

#same, why not CMAKE... 
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

#this is in file2
if (COMP_NAME STREQUAL elm)
  string(APPEND FFLAGS " -hfp0")
endif()
#-hipa0 -hzero are in file 2
#do we need -hzero em ef?
string(APPEND FFLAGS " -hipa0 -hzero -em -ef -hnoacc")

#this is in file2
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")

#this is in config_machines
set(PIO_FILESYSTEM_HINTS "romio_cb_read=disable")

string(APPEND CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")
string(APPEND CXX_LIBS " -lstdc++")

string(APPEND CXXFLAGS " -I$ENV{MPICH_DIR}/include --amdgpu-target=gfx90a")
string(APPEND SLIBS    " -L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa")

#ZEN3? scream builds its own
string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On")

