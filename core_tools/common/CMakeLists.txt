project(common)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(common INTERFACE)
target_include_directories(common INTERFACE include/)
target_include_directories(common INTERFACE "../../../eigen-3.4.0/")
