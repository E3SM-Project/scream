if (NOT DEBUG)
  string(APPEND CFLAGS " -O2")
  string(APPEND FFLAGS " -O2")
endif()
if (DEBUG)
  string(APPEND CXXFLAGS " -O0")
endif()
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_SLASHPROC")
endif()
string(APPEND SLIBS " -L$ENV{ESSL_PATH}/lib64 -lessl")
string(APPEND CXX_LIBS " -lstdc++")
set(MPICXX "mpiCC")
set(PIO_FILESYSTEM_HINTS "gpfs")