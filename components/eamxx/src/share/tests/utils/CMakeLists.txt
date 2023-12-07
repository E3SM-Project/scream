# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  include(ScreamUtils)

  # Test utils
  CreateUnitTest(utils "utils_tests.cpp")

  # Test column ops
  CreateUnitTest(column_ops "column_ops.cpp")

  # Test common physics functions
  CreateUnitTest(common_physics "common_physics_functions_tests.cpp")
endif()
