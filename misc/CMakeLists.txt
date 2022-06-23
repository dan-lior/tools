project(misc)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(misc source/misc.cpp)
target_include_directories(misc PUBLIC include/)

add_subdirectory(../common/ build_common)
target_link_libraries(misc common)

target_include_directories(misc PUBLIC "../../json/single_include/")
