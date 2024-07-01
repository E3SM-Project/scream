include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE)

set(pybind11_ROOT /home/e3sm-jenkins/.local/lib/python3.8/site-packages CACHE PATH "Path to pybind11 installation folder")
