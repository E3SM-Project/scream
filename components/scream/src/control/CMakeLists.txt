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
# link_directories(${SCREAM_LIBRARY_DIRS} ${SCREAM_TPL_LIBRARY_DIRS})

# add_library(scream_control ${SCREAM_CONTROL_SOURCES})
# target_include_directories(scream_control PUBLIC ${SCREAM_SRC_DIR})
# target_link_libraries(scream_control scream_dynamics scream_share ${SCREAM_TPL_LIBRARIES})
