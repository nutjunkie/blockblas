add_executable(feast feast.C util.C)
target_include_directories(feast PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))
target_link_libraries(feast LINK_PUBLIC blockblas)

add_executable(tile_test tile_test.C util.C)
target_include_directories(tile_test PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))
target_link_libraries(tile_test LINK_PUBLIC blockblas)
