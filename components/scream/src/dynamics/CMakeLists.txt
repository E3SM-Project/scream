MACRO(GatherDynamicsSources)
  SET (DYNAMICS_DIR ${SCREAM_SRC_DIR}/dynamics)

  set(SCREAM_DYNAMICS_SOURCES
    ${DYNAMICS_DIR}/atmosphere_dynamics.cpp
  )

  set(SCREAM_DYNAMICS_HEADERS
    ${DYNAMICS_DIR}/atmosphere_dynamics.hpp
  )

ENDMACRO(GatherDynamicsSources)
