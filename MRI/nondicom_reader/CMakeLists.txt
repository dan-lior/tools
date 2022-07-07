project(nondicom_reader)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)
#set (CMAKE_CXX_FLAGS "-Wall")
#set (CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wunused")

add_subdirectory(../../core_tools build_core_tools)

add_library(nondicom_reader source/read_nondicom.cpp)
target_include_directories(nondicom_reader PUBLIC "include/")

target_include_directories(nondicom_reader PUBLIC "../../../eigen-3.4.0")
target_include_directories(nondicom_reader PUBLIC "/home/dan/git-repos/json/single_include/")

target_include_directories(nondicom_reader PUBLIC "../../core_tools/common/include")
target_include_directories(nondicom_reader PUBLIC "../../core_tools/voxels/include")
target_include_directories(nondicom_reader PUBLIC "../../core_tools/affinities/include")
target_include_directories(nondicom_reader PUBLIC "../../core_tools/misc_math/include")
target_include_directories(nondicom_reader PUBLIC "../../core_tools/misc/include")

target_link_libraries(nondicom_reader PUBLIC misc)


# just for testing...
add_executable(dummy main.cpp)
target_link_libraries(dummy nondicom_reader)
target_link_libraries(dummy voxels)