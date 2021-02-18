add_executable(spherium spherium.C util.C)

target_include_directories(spherium PUBLIC 
   $(CMAKE_CURRENT_SOURCE_DIR) 
   $(CMAKE_CURRENT_SOURCE_DIR)/../feast
)

target_link_libraries(spherium
   feast
   blockblas
   irc
   ifcore
   ifcoremt
)

#if(MPI_COMPILE_FLAGS)
#   set_target_properties(spherium PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#endif()
#
#if(MPI_LINK_FLAGS)
#   set_target_properties(spherium PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
#endif()

#add_executable(zfeast PCdense_zfeast_syev.c)
#target_include_directories(zfeast PUBLIC 
#   $(CMAKE_CURRENT_SOURCE_DIR)/../feast
#)
#
#target_link_libraries(zfeast
#   feast
#   irc
#   ifcore
#   ifcoremt
#)


add_executable(tile_test tile_test.C util.C)
target_include_directories(tile_test PUBLIC $(CMAKE_CURRENT_SOURCE_DIR) )
target_link_libraries(tile_test LINK_PUBLIC blockblas)

