string(APPEND CPPDEFS " -DSYSDARWIN")
if (COMP_CLASS STREQUAL cpl)
  string(APPEND LDFLAGS " -all_load")
endif()
