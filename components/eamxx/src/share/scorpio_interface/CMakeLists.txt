add_library(eamxx_scorpio_interface
  scream_scorpio_interface.F90
  scream_scorpio_interface.cpp
  scream_scorpio_interface_iso_c2f.F90
)

target_link_libraries(eamxx_scorpio_interface PUBLIC eamxx_utils ekat piof pioc)
set_target_properties(eamxx_scorpio_interface PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(eamxx_scorpio_interface PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules)

if (DEFINED ENV{ADIOS2_ROOT})
  target_include_directories(eamxx_scorpio_interface PRIVATE $ENV{ADIOS2_ROOT}/include)
endif ()

if (SCREAM_CIME_BUILD)
  # To get atm pio subsystem info
  target_link_libraries(eamxx_io PUBLIC csm_share)
endif()

