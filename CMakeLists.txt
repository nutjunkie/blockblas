cmake_minimum_required (VERSION 3.13.0)

if (MYMPI)
   message("Configuring MPI version")

   find_package(MPI REQUIRED)
   include_directories(${MPI_INCLUDE_PATH})
   add_compile_definitions(MYMPI)

   set(CMAKE_CXX_COMPILER "mpicxx")
   set(CMAKE_FC_COMPILER "mpiifort")
   set(CMAKE_C_COMPILER "mpiicc")

   #set(CMAKE_Fortran_FLAGS "-g -O3 -qopenmp -cpp -fPIC -DMKL -DMPI")
   set(CMAKE_Fortran_FLAGS "-g -O0 -debug all -qopenmp -cpp -fPIC -DMKL -DMPI")
   set(CMAKE_CXX_FLAGS  "-g -O0 -debug all -std=c++11 -O2 -fopenmp -fPIC -mkl")
   set(CMAKE_C_FLAGS  "-g -O0 -debug all -qopenmp -fPIC  -mkl")

else()
   message("Configuring OpenMP only version")
   set(CMAKE_CXX_COMPILER "icc")
   set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -g  -O0 -qopenmp -funroll-loops -ffast-math  -fPIC -mkl=sequential")
endif()


project(blockblas)
add_subdirectory(src)
add_subdirectory(feast)
add_subdirectory(tests)
