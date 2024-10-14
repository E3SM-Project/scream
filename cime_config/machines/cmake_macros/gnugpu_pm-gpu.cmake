string(APPEND CONFIG_ARGS " --host=cray")
set(USE_CUDA "TRUE")
string(APPEND CPPDEFS " -DGPU")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")
string(APPEND CMAKE_CUDA_FLAGS " -ccbin CC -O2 -arch sm_80 --use_fast_math")
#string(APPEND CMAKE_CUDA_FLAGS " -ccbin CC -O2 -arch sm_80 --use_fast_math -lineinfo") # ndk, only beneficial when using ncu profiling *and* using "ncu --set full"
string(APPEND KOKKOS_OPTIONS " -DKokkos_ARCH_AMPERE80=On -DKokkos_ENABLE_CUDA=On -DKokkos_ENABLE_CUDA_LAMBDA=On -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=Off -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=Off")
set(CMAKE_CUDA_ARCHITECTURES "80")

#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -g -lineinfo")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -lineinfo --use_fast_math")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -lineinfo --maxrregcount 64")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -lineinfo")
string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2 -g")
set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")

# ndk need this to allow nvtx calls in fortran. needs to only be at final link, hence only for comp_name==cpl
if (COMP_NAME STREQUAL cpl)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " /opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda/12.2/targets/x86_64-linux/lib/libnvToolsExt.so")
endif()

# experimenting with YAKL flags
#string(APPEND CPPDEFS  " -DYAKL_AUTO_PROFILE -DYAKL_PROFILE")

#string(APPEND CPPDEFS  "-DNVCC_WRAPPER_SHOW_COMMANDS_BEING_RUN")

# don't actually need MAM, plus it auto adds openmp flags which we do not want
set(SCREAM_ENABLE_MAM OFF) # ndk
