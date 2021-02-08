cmake_minimum_required (VERSION 3.13.0)
set(CMAKE_CXX_COMPILER "icc")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -g  -O2 -qopenmp -funroll-loops -ffast-math  -fPIC -mkl=sequential")

project (blockblas)
add_subdirectory(src)
add_subdirectory(tests)
