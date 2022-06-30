project(nondicom_reader)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_subdirectory(../../core_tools build_core_tools)

add_executable(nondicom_reader source/main.cpp)
target_include_directories(nondicom_reader PUBLIC "../../../eigen-3.4.0")
target_include_directories(nondicom_reader PUBLIC "/home/dan/git-repos/json/single_include/")

target_include_directories(nondicom_reader PUBLIC "../../core_tools/common/include")
target_include_directories(nondicom_reader PUBLIC "../../core_tools/voxels/include")
target_include_directories(nondicom_reader PUBLIC "../../core_tools/affinities/include")
target_include_directories(nondicom_reader PUBLIC "../../core_tools/misc_math/include")


