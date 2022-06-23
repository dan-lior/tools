project(affinities)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(affinities source/dummy.cpp)
target_include_directories(affinities PUBLIC include/)

target_include_directories(affinities PUBLIC ../common/include/) 
target_include_directories(affinities PUBLIC "../../eigen-3.4.0/")
