if (NOT SCREAM_BASELINES_ONLY)
  include(ScreamUtils)

  CreateUnitTest(nudging_tests "nudging_tests.cpp"
      LIBS nudging eamxx_io
      LABELS "physics_nudging" )

endif()
