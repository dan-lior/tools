project(triangulations)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(triangulations source/triangulation.cpp)
target_include_directories(triangulations PUBLIC include/)

target_link_libraries(triangulations disjoint_unions)
target_link_libraries(triangulations polylines)

# these next three are already provided by polylines:
#target_link_libraries(triangulations misc)
#target_link_libraries(triangulations misc_math)
#target_link_libraries(triangulations affinities)

