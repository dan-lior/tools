project(polylines)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(disjoint_unions source/disjoint_union.cpp)
target_include_directories(disjoint_unions PUBLIC include/)

target_include_directories(disjoint_unions PUBLIC ../common/include/)
target_include_directories(disjoint_unions PUBLIC "../../../eigen-3.4.0/")




