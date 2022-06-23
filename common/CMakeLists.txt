project(common)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(common source/dummy.cpp)
target_include_directories(common PUBLIC include/)
target_include_directories(common PUBLIC "../../eigen-3.4.0/")
