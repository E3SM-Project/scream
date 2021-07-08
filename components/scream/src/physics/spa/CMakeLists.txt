set(SPA_SRCS
  atmosphere_prescribed_aero.cpp
  spa.cpp
)

set(SPA_HEADERS
  atmosphere_prescribed_aero.hpp
  spa_functions.hpp
  spa_main_impl.hpp
)

add_library(spa ${SPA_SRCS})
target_include_directories(spa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(spa physics_share scream_share)
target_compile_options(spa PUBLIC)
