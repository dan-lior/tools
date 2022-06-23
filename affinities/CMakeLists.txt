project(affinities)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wunused") 

add_library(affinities INTERFACE)
target_include_directories(affinities INTERFACE include/)

add_subdirectory(../common/ build_common)
target_include_directories(affinities INTERFACE ../common/include/)
