
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" CXX11_SUPPORTED)
if(${CXX11_SUPPORTED})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(FATAL_ERROR, "The C++ compiler does not support C++11")
endif()
