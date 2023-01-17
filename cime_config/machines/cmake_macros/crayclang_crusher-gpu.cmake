#crayclang = file1
#crayclang-crusher = file2
#crayclanggpu-crusher = file3

#present in gpu file3 and in file1, file2
if (compile_threaded)
  string(APPEND CFLAGS " -fopenmp")
  string(APPEND FFLAGS " -fopenmp")
  string(APPEND CXXFLAGS " -fopenmp")
  string(APPEND LDFLAGS " -fopenmp")
endif()

#present in crayclang-crusher file2
string(APPEND SLIBS " -L$ENV{PNETCDF_PATH}/lib -lpnetcdf")

#present in file3 and in file2, why in file 3 then?
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
#present in gpu file, same
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")

#not present anywhere, instead in gpu file3
#Az: set(PIO_FILESYSTEM_HINTS "romio_cb_read=disable") -- 
#whcih is already set in config_machines
set(PIO_FILESYSTEM_HINTS "gpfs")
#present in gpu file3 and in file2, why twice?
string(APPEND CXX_LIBS " -lstdc++")

#not present 
#instead in files 1,2,3 vars MPICC, SCC, etc. set
SET(CMAKE_C_COMPILER "cc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "hipcc" CACHE STRING "")

#present in file3 with  --amdgpu-target=gfx90a
#did i not use it cause kokkos puts it in?
string(APPEND CXXFLAGS " -I${MPICH_DIR}/include")
#present in file3 as APPEND_SLIBS
string(APPEND LDFLAGS " -L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.16/gtl/lib -lmpi_gtl_hsa")

#not present
# For YAKL's -lroctx64 -lrocfft; the rocm module doesn't set this.
string(APPEND LDFLAGS " -L$ENV{ROCM_PATH}/lib")

# -O2 present in file3
#-hnoacc present in file 3
#-hfp0 present in file 3 for ELM only
#-hipa0 present in file 3
#this resolves a crash in mct in docn init
if (NOT DEBUG)
  string(APPEND CFLAGS " -O2 -hnoacc -hfp0 -hipa0")
  string(APPEND FFLAGS " -O2 -hnoacc -hfp0 -hipa0")
endif()

#present in file3
string(APPEND CPPDEFS " -DCPRCRAY")

#add this to the file3
set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "See SCREAM issue #2080.")
