#########################################
#       SCORPIO interface library       #
#########################################

# This small lib contains some interfaces to scorpio, which are formulated
# in terms of names rather than IDs. Internally, we store a list of
# small structs, holding handlers to PIO data and extra metadata, which
# allows to perform sanity checks and print more helpful messages

set(SCREAM_SCORPIO_INTERFACE_SRCS
  scream_scorpio_types.cpp
  scream_scorpio_interface.cpp
)

add_library(scream_scorpio_interface ${SCREAM_SCORPIO_INTERFACE_SRCS})
target_link_libraries(scream_scorpio_interface PUBLIC pioc ekat)
target_include_directories(scream_scorpio_interface PUBLIC
  ${SCREAM_BIN_DIR}/src   # For scream_config.h
)

if (DEFINED ENV{ADIOS2_DIR})
  target_include_directories(scream_scorpio_interface PRIVATE $ENV{ADIOS2_DIR}/include)
endif ()

if (SCREAM_CIME_BUILD)
  # Add interface to E3SM shr lib (to retrieve PIO subsystem info)
  target_sources (scream_scorpio_interface PRIVATE
    scream_shr_interface_c2f.F90
  )
  target_compile_definitions (scream_scorpio_interface PRIVATE SCREAM_CIM_BUILD)
  set_target_properties(scream_scorpio_interface PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
  )
  target_include_directories(scream_scorpio_interface PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/modules)
  target_link_libraries(scream_scorpio_interface PUBLIC csm_share)
endif()

##################################
#       EAMxx I/O library        #
##################################

# This library allows for a simple interaction between EAMxx
# data structures (such as grids, fields, remappers,...) and
# the scorpio interface library

set(SCREAM_SCORPIO_SRCS
  scream_output_manager.cpp
  scorpio_input.cpp
  scorpio_output.cpp
  scream_io_utils.cpp
)

# Create io lib
add_library(scream_io ${SCREAM_SCORPIO_SRCS})

target_link_libraries(scream_io PUBLIC scream_share scream_scorpio_interface)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
