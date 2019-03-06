INCLUDE (ScreamUtils)

SET (num_mpi_ranks 1)
SET (config_defines "")

GatherDynamicsSources()
GatherControlSources()

SET (SOURCES
     ${SCREAM_CONTROL_SOURCES}
     ${SCREAM_DYNAMICS_SOURCES}
     control_tests.cpp
)

# Test control folder
CreateUnitTest(control "${SOURCES}" "scream_share" num_mpi_ranks config_defines)
