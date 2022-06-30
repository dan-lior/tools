project(conics)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(conics source/conic.cpp)
target_include_directories(conics PUBLIC include/)
target_include_directories(conics PUBLIC ../common/include/)

target_link_libraries(conics misc)
target_link_libraries(conics misc_math)
target_link_libraries(conics affinities)


