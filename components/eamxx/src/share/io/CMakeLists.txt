add_library(eamxx_io
  scream_output_manager.cpp
  scorpio_input.cpp
  scorpio_output.cpp
  scream_io_utils.cpp
  eamxx_time_interpolation.cpp
)
target_link_libraries (eamxx_io PUBLIC
  eamxx_scorpio_interface
  eamxx_field
  eamxx_grid
  eamxx_remap
  eamxx_utils
  eamxx_atm_process
)
