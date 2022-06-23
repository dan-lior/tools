project(common)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wunused") 

INCLUDE_DIRECTORIES("/home/dan/git-repos/eigen-3.4.0")

add_library(common INTERFACE)
target_include_directories(common INTERFACE include/)
