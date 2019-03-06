
set(SHOC_TEST_SRCS
  shoc_tests.cpp
)

add_executable(shoc_tests ${SHOC_TEST_SRCS})

target_include_directories(shoc_tests PUBLIC ${CATCH_INCLUDE_DIR} ${SCREAM_SRC_DIR})

target_link_libraries(shoc_tests shoc)

add_test(run_shoc_tests shoc_tests)
