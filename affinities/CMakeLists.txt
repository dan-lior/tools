project(affinities)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(affinities source/dummy.cpp)
target_include_directories(affinities PUBLIC include/)

add_subdirectory(../common/ build_common)
target_link_libraries(affinities common)