#interactive job
#bsub -W 2:00 -nnodes 1 -P cli115 -Is /bin/bash


#cmake -C ~/acme-fork-lb/components/homme/cmake/machineFiles/summit.cmake -DHOMMEXX_MPI_ON_DEVICE=FALSE ~/acme-fork-lb/components/homme/

#cmake -C ~/acme-MASTER-GB/components/homme/cmake/machineFiles/crusher-gpumpi.cmake -DE3SM_KOKKOS_PATH=/ccs/home/onguba/kokkos-crusher-june2022/bld-hipcc ~/acme-MASTER-GB/components/homme/

#SET (HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")
SET (HOMMEXX_CUDA_MAX_WARP_PER_TEAM "16" CACHE STRING  "")

SET (NETCDF_DIR $ENV{CRAY_NETCDF_HDF5PARALLEL_PREFIX} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{CRAY_HDF5_PARALLEL_PREFIX} CACHE FILEPATH "")

SET(BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")

SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")

SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

SET(USE_TRILINOS OFF CACHE BOOL "")

#CUDA_BUILD is set in SetCompilersFlags, after findPackage(Cuda)
#i haven't extend it to hip, set it here instead
SET(HIP_BUILD TRUE CACHE BOOL "")

#uncomment this if using internal kokkos build
SET(Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
####SET(CMAKE_CXX_STANDARD "14" CACHE STRING "")
#SET(Kokkos_ENABLE_DEBUG OFF CACHE BOOL "")
SET(Kokkos_ARCH_VEGA90A ON CACHE BOOL "")
SET(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
SET(Kokkos_ENABLE_HIP ON CACHE BOOL "")
####SET(Kokkos_ENABLE_CUDA_LAMBDA OFF CACHE BOOL "")
SET(Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")

SET(CMAKE_C_COMPILER "cc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS "--gcc-toolchain=$ENV{MEMBERWORK}/cli115/workaround" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "hipcc" CACHE STRING "")

#SET(MPICH_DIR "/opt/cray/pe/mpich/8.1.12/ofi/crayclang/10.0" CACHE STRING "")

SET(Extrae_LIBRARY "-I$ENV{CRAY_MPICH_DIR}/include -L$ENV{CRAY_MPICH_DIR}/lib -lmpi $ENV{PE_MPICH_GTL_DIR_amd_gfx90a} $ENV{PE_MPICH_GTL_LIBS_amd_gfx90a}" CACHE STRING "")

SET(ADD_Fortran_FLAGS "--gcc-toolchain=$ENV{MEMBERWORK}/cli115/workaround -h flex_mp=intolerant -h thread0 -G2 ${Extrae_LIBRARY}" CACHE STRING "")
SET(ADD_C_FLAGS "-g -O -ffp-model=strict ${Extrae_LIBRARY}" CACHE STRING "")
SET(ADD_CXX_FLAGS "-g -std=c++14 -O -ffp-model=strict -munsafe-fp-atomics --offload-arch=gfx90a -fno-gpu-rdc -Wno-unused-command-line-argument -Wno-unsupported-floating-point-opt -Wno-#pragma-messages ${Extrae_LIBRARY}" CACHE STRING "")
SET(ADD_LINKER_FLAGS "-g -O ${Extrae_LIBRARY}" CACHE STRING "")


set (ENABLE_OPENMP OFF CACHE BOOL "")
set (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

set (HOMME_TESTING_PROFILE "short" CACHE STRING "")
SET (BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")

set (USE_NUM_PROCS 8 CACHE STRING "")

SET (USE_MPI_OPTIONS "-c7 --gpu-bind=closest --gpus-per-task=1" CACHE FILEPATH "")
