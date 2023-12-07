# Create test_support lib
add_library(scream_test_support
  mesh_free_grids_manager.cpp
  scream_test_session.cpp
)
target_link_libraries (scream_test_support PUBLIC eamxx_grid eamxx_remap ekat)
