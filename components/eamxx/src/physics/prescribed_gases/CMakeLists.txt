add_library(prescribed_gases eamxx_gases_process_interface.cpp)
target_compile_definitions(prescribed_gases PUBLIC EAMXX_HAS_PRESCRIBED_GASES)
target_link_libraries(prescribed_gases physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
    #add_subdirectory(tests)
endif()
