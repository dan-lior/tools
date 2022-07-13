project(voxels)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(voxels "source/voxel_box.cpp")
target_include_directories(voxels PUBLIC "include/")

target_include_directories(voxels PUBLIC "../common/include/") 
target_include_directories(voxels PUBLIC "../../../eigen-3.4.0/")
target_include_directories(voxels PUBLIC "../misc/include/") 
target_include_directories(voxels PUBLIC "../misc_math/include/") 
target_include_directories(voxels PUBLIC "../affinities/include/") 

target_link_libraries(voxels misc)
target_link_libraries(voxels misc_math)

