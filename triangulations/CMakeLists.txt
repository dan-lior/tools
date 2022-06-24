project(triangulations)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(triangulations source/triangulation.cpp)
target_include_directories(triangulations PUBLIC include/)

target_include_directories(triangulations PUBLIC ../common/include/)

add_subdirectory(../disjoint_unions/ build_disjoint_unions)
target_link_libraries(triangulations disjoint_unions)

add_subdirectory(../polylines/ build_polylines)
target_link_libraries(triangulations polylines)

# these next three are already provided by polylines:

# add_subdirectory(../misc/ build_misc)
# target_link_libraries(triangulations misc)

# add_subdirectory(../misc_math/ build_misc_math)
# target_link_libraries(triangulations misc_math)

# add_subdirectory(../affinities/ build_affinities)
# target_link_libraries(triangulations affinities)

