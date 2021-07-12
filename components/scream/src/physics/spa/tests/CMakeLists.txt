include(ScreamUtils)

set(NEED_LIBS spa physics_share scream_share)
set(SPA_TESTS_SRCS
  spa_horizontal_interpolation_test.cpp
  spa_unit_tests.cpp
    ) # SPA_TESTS_SRCS

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  CreateUnitTest(spa_tests "${SPA_TESTS_SRCS}" "${NEED_LIBS}"
                 THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
                 LABELS "spa;physics")
endif()
