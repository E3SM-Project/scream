add_library(eamxx_fvphyshack
  eamxx_fv_phys_rrtmgp_active_gases_workaround.cpp
)

target_link_libraries (eamxx_fvphyshack PUBLIC eamxx_remap ekat)
