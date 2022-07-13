project(misc)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(misc source/misc.cpp)
target_include_directories(misc PUBLIC include/)

target_include_directories(misc PUBLIC "../common/include/") 
target_include_directories(misc PUBLIC "../misc_math/include/") 
target_include_directories(misc PUBLIC "../../../eigen-3.4.0/")
target_include_directories(misc PUBLIC "../../../json/single_include/")

