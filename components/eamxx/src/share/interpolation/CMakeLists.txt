add_library(eamxx_interpolation
  eamxx_time_interpolation.cpp
)

target_link_libraries (eamxx_interpolation PUBLIC eamxx_field)
