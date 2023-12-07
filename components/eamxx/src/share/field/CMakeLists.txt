add_library (eamxx_field
  field_alloc_prop.cpp
  field_identifier.cpp
  field_header.cpp
  field_layout.cpp
  field_tracking.cpp
  field.cpp
  field_group.cpp
)

target_link_libraries (eamxx_field PUBLIC eamxx_core eamxx_utils)
