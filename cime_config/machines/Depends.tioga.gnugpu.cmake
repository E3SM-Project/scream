list(APPEND SLIST
  eamxx/src/physics/p3/atmosphere_microphysics_run.cpp
)

if (DEBUG)
  foreach(ITEM IN LISTS SLIST)
    #e3sm_add_flags("${ITEM}" " -fno-default-inline -finline-limit=512")
    #e3sm_add_flags("${ITEM}" " -mllvm -amdgpu-early-inline-all=true -mllvm -amdgpu-function-calls=false")
  endforeach()
endif()
