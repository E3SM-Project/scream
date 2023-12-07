add_library(eamxx_remap
  abstract_remapper.cpp
  coarsening_remapper.cpp
  horiz_interp_remapper_base.cpp
  refining_remapper_p2p.cpp
  vertical_remapper.cpp
  horizontal_remap_utility.cpp
)

if (EAMXX_ENABLE_EXPERIMENTAL_CODE)
  target_sources (eamxx_remap PUBLIC 
    refining_remapper_rma.cpp
  )
endif()

target_link_libraries(eamxx_remap PUBLIC eamxx_grid eamxx_scorpio_interface)
