project(misc_math)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(misc_math source/misc_math.cpp)
target_include_directories(misc_math PUBLIC include/)

target_include_directories(misc_math PUBLIC ../common/include/) 
target_include_directories(misc_math PUBLIC "../../eigen-3.4.0/")
