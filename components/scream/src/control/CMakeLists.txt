MACRO (GatherControlSources)
  SET (CONTROL_DIR ${SCREAM_SRC_DIR}/control)

  set(SCREAM_CONTROL_SOURCES
    ${CONTROL_DIR}/cxx_f90_interface.cpp
    ${CONTROL_DIR}/atmosphere_driver.cpp
    ${CONTROL_DIR}/atmosphere_process_factory.cpp
    ${CONTROL_DIR}/atmosphere_process_group.cpp
    ${CONTROL_DIR}/surface_coupling.cpp
  )

  set(SCREAM_CONTROL_HEADERS
    atmosphere_driver.hpp
    atmosphere_process_factory.hpp
    atmosphere_process_group.hpp
    surface_coupling.hpp
  )
ENDMACRO(GatherControlSources)
