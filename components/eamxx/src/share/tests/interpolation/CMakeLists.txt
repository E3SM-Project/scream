# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  include(ScreamUtils)

  # Test vertical interpolation
  CreateUnitTest(vertical_interp "vertical_interp_tests.cpp")

  # Test vertical remap
  CreateUnitTest(time_interpolation "eamxx_time_interpolation_tests.cpp"
    LIBS eamxx_io
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
endif()
