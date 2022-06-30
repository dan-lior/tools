project(affinities)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(affinities INTERFACE)
target_include_directories(affinities INTERFACE include/)
