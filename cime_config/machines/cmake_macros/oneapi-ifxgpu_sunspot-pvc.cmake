
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -fsycl-device-code-split=per_kernel -fsycl-max-parallel-link-jobs=16")
if (compile_threaded)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fiopenmp -fopenmp-targets=spir64")
endif()

if (DEBUG)
#undefined reference to `__msan....
#https://community.intel.com/t5/Intel-Fortran-Compiler/Linking-errors-when-using-memory-sanitizer-in-fortran-project/m-p/1521476
#When you compile with -check uninit (or -check all) you also need to link with that compiler option.
#  string(APPEND CMAKE_EXE_LINKER_FLAGS " -check uninit")
endif()

string(APPEND KOKKOS_OPTIONS " -DCMAKE_CXX_STANDARD=17 -DKokkos_ENABLE_SERIAL=On -DKokkos_ARCH_INTEL_PVC=On -DKokkos_ENABLE_SYCL=On -DKokkos_ENABLE_EXPLICIT_INSTANTIATION=Off -DCMAKE_CXX_FLAGS='-fsycl-device-code-split=per_kernel'")
string(APPEND SYCL_FLAGS " -\-intel -fsycl -fsycl-targets=spir64_gen -mlong-double-64 -Xsycl-target-backend \"-device 12.60.7\"")

set(SCREAM_MPI_ON_DEVICE ON CACHE STRING "")



