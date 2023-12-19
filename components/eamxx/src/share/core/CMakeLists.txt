add_library(eamxx_core
  scream_config.cpp
  scream_session.cpp
)

# SCREAM_SRC_DIR allows all downstream targets
# to include files using paths relative to
# components/eamxx/src
target_include_directories(eamxx_core PUBLIC
  ${SCREAM_SRC_DIR}
  ${SCREAM_BIN_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_options(eamxx_core PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>
)
target_link_libraries(eamxx_core PUBLIC ekat)

# This used to be in `scream_config.h/scream_config.f`, but we did accidentally remove
# the includes at least twice. To avoid this from happening again, we add this crucial
# macro as a compile definition of the eamxx_core, which is linked to all other
# scream libraries.
if (SCREAM_CIME_BUILD)
  target_compile_definitions (eamxx_core PUBLIC SCREAM_CIME_BUILD)
endif()

# We have some issues with RDC and repeated libraries in the link line.
# It's a known issue, and nvcc_wrapper has a flag for handling this.
if (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
  if(CMAKE_VERSION VERSION_LESS "3.13.0")
    target_link_libraries(eamxx_core PUBLIC "--remove-duplicate-link-files")
  else()
    target_link_options(eamxx_core PUBLIC "--remove-duplicate-link-files")
  endif()
endif()
include (EkatSetCompilerFlags)
SetOmpFlags(eamxx_core Fortran)
