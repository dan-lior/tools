project(coretools)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

add_library(coretools dummy.cpp)

add_subdirectory(common build_common)
target_include_directories(coretools PUBLIC common/include/) 

add_subdirectory(affinities build_affinities)
target_include_directories(coretools PUBLIC affinities/include/) 

add_subdirectory(voxels build_voxels)
target_link_libraries(coretools voxels) 

add_subdirectory(misc build_misc)
target_link_libraries(coretools misc)

add_subdirectory(misc_math build_misc_math)
target_link_libraries(coretools misc_math)

add_subdirectory(polylines build_polylines)
target_link_libraries(coretools polylines)

add_subdirectory(conics build_conics)
target_link_libraries(coretools conics)

add_subdirectory(disjoint_unions build_disjoint_unions)
target_link_libraries(coretools disjoint_unions)

add_subdirectory(triangulations build_triangulations)
target_link_libraries(coretools triangulations)

#add_subdirectory(quadrics build_quadrics)
#target_link_libraries(coretools quadrics)
