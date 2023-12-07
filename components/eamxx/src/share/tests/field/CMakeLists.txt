# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  include(ScreamUtils)

  # Test fields
  CreateUnitTest(field "field_tests.cpp")

  # Test field utils
  CreateUnitTest(field_utils "field_utils.cpp"
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
endif()
