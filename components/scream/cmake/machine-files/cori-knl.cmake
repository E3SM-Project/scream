set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

# Load knl arch and openmp backend for kokkos
include (${EKAT_MACH_FILES_PATH}/kokkos/intel-knl.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
       set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if(DEBUG)
    #set(CMAKE_CXX_FLAGS "-DRRTMGP_EXPENSIVE_CHECKS -Werror=vla -mcmodel=large"  CACHE STRING "" FORCE)
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10) # only works with gnu v10 and above
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch -g -Wall -fbacktrace -fcheck=bounds -ffpe-trap=invalid,zero,overflow"  CACHE STRING "" FORCE)
    endif()
  else() # OPT
    #set(CMAKE_CXX_FLAGS "-g1"  CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS "-DYAKL_PROFILE -DHAVE_MPI"  CACHE STRING "" FORCE)
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)  # only works with gnu v10 and above
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch -g -Wall -fbacktrace"  CACHE STRING "" FORCE)
    endif()
  endif()
else()   # intel
  if(DEBUG)
    set(CMAKE_CXX_FLAGS "-traceback"  CACHE STRING "" FORCE)
    set(CMAKE_Fortran_FLAGS "-fpe0 -O0 -g -check bounds -ftz -traceback"  CACHE STRING "" FORCE)
  else() # OPT
    set(CMAKE_CXX_FLAGS "-g -traceback -DYAKL_PROFILE -DHAVE_MPI"  CACHE STRING "" FORCE)
    set(CMAKE_Fortran_FLAGS "-O2 -g -ftz -traceback"  CACHE STRING "" FORCE)
  endif()
endif()

if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if (compile_threaded)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      string(APPEND CMAKE_C_FLAGS " -fopenmp")
      string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
      string(APPEND CMAKE_Fortran_FLAGS " -fopenmp")
      string(APPEND CMAKE_EXE_LINKER_FLAGS " -fopenmp")
    else() #intel
      string(APPEND CMAKE_C_FLAGS " -qopenmp")
      string(APPEND CMAKE_CXX_FLAGS " -qopenmp")
      string(APPEND CMAKE_Fortran_FLAGS " -qopenmp")
      string(APPEND CMAKE_EXE_LINKER_FLAGS " -qopenmp")
    endif()
  endif()
endif()

# Fixes some openmpi link problems we observed on cori. This hack is
# not necessary if CRAYPE_LINK_TYPE=dynamic is in the environment.
set(SCREAM_CORI_HACK True CACHE BOOL "")
set(SCREAM_MACHINE "cori-knl" CACHE STRING "")
