# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  include(ScreamUtils)

  # Test coarsening remap
  CreateUnitTest(coarsening_remapper "coarsening_remapper_tests.cpp"
    LIBS eamxx_io
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  if (EAMXX_ENABLE_EXPERIMENTAL_CODE)
    # Test refining remap (RMA version)
    CreateUnitTest(refining_remapper_rma "refining_remapper_rma_tests.cpp"
      LIBS eamxx_io
      MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
  endif()

  # Test refining remap (P2P version)
  CreateUnitTest(refining_remapper_p2p "refining_remapper_p2p_tests.cpp"
    LIBS eamxx_io
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test vertical remap
  CreateUnitTest(vertical_remapper "vertical_remapper_tests.cpp"
    LIBS eamxx_io
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test horizontal remapping utility
  CreateUnitTest(horizontal_remap "horizontal_remap_test.cpp"
    LIBS eamxx_io
    LABELS horiz_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )
endif()
