project(polylines)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(polylines source/polyline.cpp)
target_include_directories(polylines PUBLIC include/)

target_include_directories(polylines PUBLIC ../common/include/) 
target_include_directories(polylines PUBLIC "../../../eigen-3.4.0/")

target_link_libraries(polylines misc)
target_link_libraries(polylines misc_math)
target_link_libraries(polylines affinities)
