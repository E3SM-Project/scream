# In this folder we create several small libs, implementing
# several aspects of the eamxx infrastructure

# The infrastructure (except for io) used to all be compiled
# as a single library. We still offer that lib, for bwd compatibility
# but also for ease of use: most of the code needs several parts
# of the infrastructure, so they can simply link against eamxx_share
# However, this gives us more granularity in how we link different
# parts of eamxx together, avoiding circular dependencies,
# and making it easier to debug

# Contains basic definitions and config settings
add_subdirectory (core)

# Common utils that are needed across the board
add_subdirectory (util)

# Field stack classes
add_subdirectory (field)

# Grid classes
add_subdirectory (grid)

# Hooks to scorpio functions
add_subdirectory (scorpio_interface)

# Remapping between grids interfaces
add_subdirectory (remap)

# Classes that perform checks on fields
add_subdirectory (property_checks)

# The grids/field manager classes link together a few of the above
# libs, so handle them separately
add_subdirectory (managers)

# Atmosphere processes interfaces
add_subdirectory (atm_process)

# IO library
add_subdirectory(io)

# Interpolation routines (vertical and time)
add_subdirectory(interpolation)

# Some hacks needed to deal with fv physics grid
# If/when they are no longer necessary, we can remove them
add_subdirectory (fvphyshack)

# Classes to compute diagnostic quantities
add_subdirectory (diagnostics)

# Wrap all the above mini-libs into a single lib, for convenience
add_library (eamxx_share INTERFACE)
target_link_libraries (eamxx_share INTERFACE
  eamxx_core
  eamxx_utils
  eamxx_field
  eamxx_grid
  eamxx_scorpio_interface
  eamxx_remap
  eamxx_property_checks
  eamxx_atm_process
  eamxx_interpolation
  eamxx_io
  eamxx_managers
)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory (unit_test_support)
endif()

if (NOT SCREAM_CIME_BUILD)
  # The "build_cf_dictionary" target downloads an XML file containing valid field
  # names and aliases based on the CF conventions, and transforms it into a YAML
  # file that we use to validate our own fields against it. This target is only
  # available if CMake finds wget on your system.
  find_program(WGET wget)
  if (NOT WGET STREQUAL "WGET-NOTFOUND")
    set(CF_CONVENTIONS_XML_FILE cf-standard-name-table.xml)
    set(CF_CONVENTIONS_URL "http://cfconventions.org/Data/cf-standard-names/76/src/${CF_CONVENTIONS_XML_FILE}")
    add_custom_target(build_cf_dictionary
      ${WGET} ${CF_CONVENTIONS_URL}
      COMMAND ${PROJECT_SOURCE_DIR}/scripts/cf-xml-to-yaml ${CF_CONVENTIONS_XML_FILE}
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data
      COMMENT "Downloading CF conventions XML file and building YAML dictionary...")
  endif()
endif()
