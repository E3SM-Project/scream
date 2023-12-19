add_library(eamxx_utils
  scream_time_stamp.cpp
  scream_timing.cpp
  scream_utils.cpp
  scream_bfbhash.cpp
)

target_link_libraries(eamxx_utils PUBLIC eamxx_core gptl)
