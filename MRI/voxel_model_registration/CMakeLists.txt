cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)
project(dans_utils)


#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wunused") 
#set(CMAKE_CXX_FLAGS "-Wall ") 


set(quadric_dir "/home/dan/git-repos/AAOCA-workflow-tools/quadric_fit/")

#add_library(allquadrics ${quadric_dir}src/quadricfitting.cpp ${quadric_dir}src/quadricfitting.h ${quadric_dir}src/roots3and4.cpp ${quadric_dir}src/roots3and4.h)    
add_library(allquadrics ${quadric_dir}src/quadricfitting.cpp ${quadric_dir}src/roots3and4.cpp) 
target_include_directories(allquadrics PRIVATE ${quadric_dir}include)
target_include_directories(allquadrics PRIVATE ${quadric_dir}lib/jama125)
target_include_directories(allquadrics PRIVATE ${quadric_dir}lib/tnt_126)
target_include_directories(allquadrics PRIVATE ${quadric_dir}lib/CLAPACK/include)
target_link_libraries(allquadrics PRIVATE lapack blas)

include_directories (voxel_model_registration PUBLIC "/home/dan/git-repos/eigen/" )

add_executable(voxel_model_registration source/voxel_model_registration.cpp)
target_include_directories(voxel_model_registration PRIVATE ${quadric_dir}include)
target_include_directories(voxel_model_registration PRIVATE ${quadric_dir}src)
target_include_directories(voxel_model_registration PRIVATE ../include)
target_link_libraries(voxel_model_registration PRIVATE allquadrics)
