project(conics)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(conics source/conic.cpp)
target_include_directories(conics PUBLIC include/)

target_include_directories(conics PUBLIC ../common/include/)

add_subdirectory(../misc/ build_misc)
target_link_libraries(conics misc)

add_subdirectory(../misc_math/ build_misc_math)
target_link_libraries(conics misc_math)

add_subdirectory(../affinities/ build_affinities)
target_link_libraries(conics affinities)


