# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

#EAMXX_ENABLE_GPU: FALSE
#CUDA_BUILD: FALSE
#HIP_BUILD: FALSE
# how best to control if I set gpu build or not? for now, if GNU, assume CPU testing, otherwise (Cray, AMD) using GPU?
message(STATUS "tioga CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(EAMXX_ENABLE_GPU "FALSE" CACHE STRING "" FORCE)
  set(HIP_BUILD "FALSE" CACHE STRING "" FORCE)
  set(USE_HIP "FALSE" CACHE STRING "" FORCE)
else()
  set(EAMXX_ENABLE_GPU "TRUE" CACHE STRING "" FORCE)
  set(HIP_BUILD "TRUE" CACHE STRING "" FORCE)
  set(USE_HIP "TRUE" CACHE STRING "" FORCE)
endif()

#message(STATUS "tioga USE_KOKKOS=${USE_KOKKOS}")
#ndk i think this is for external kokkos -- we use internal ekat kokkos, so it should be FALSE
#set(USE_KOKKOS "TRUE" CACHE STRING "" FORCE)

message(STATUS "tioga PROJECT_NAME=${PROJECT_NAME} USE_HIP=${USE_HIP} USE_CUDA=${USE_CUDA} EAMXX_ENABLE_GPU=${EAMXX_ENABLE_GPU} USE_KOKKOS=${USE_KOKKOS}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if (USE_HIP)
    include (${EKAT_MACH_FILES_PATH}/kokkos/mi250.cmake)
    include (${EKAT_MACH_FILES_PATH}/kokkos/hip.cmake)
    #set(HIP_BUILD TRUE CACHE BOOL "" FORCE)
  else()
    include (${EKAT_MACH_FILES_PATH}/kokkos/amd-zen3.cmake)
    include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
    #include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
  endif()
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/mi250.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/hip.cmake)
endif()

#option(Kokkos_ARCH_AMPERE80 "" ON)
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

message(STATUS "tioga CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-n" CACHE STRING "")
set(SCREAM_MACHINE "tioga" CACHE STRING "")
