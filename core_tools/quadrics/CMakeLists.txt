add_subdirectory(quadric_fit)

add_library(quadrics source/quadric.cpp)
target_include_directories(quadrics PUBLIC include/)
target_include_directories(quadrics PUBLIC ../affinities/include/)

target_link_libraries(quadrics misc)
target_link_libraries(quadrics misc_math)
target_link_libraries(quadrics conics)
target_link_libraries(quadrics common)


# #target_link_libraries(quadrics affinity)
# target_link_libraries(quadrics quadric_fit)



# add_library(allquadrics ../src/quadricfitting.cpp ../src/quadricfitting.h ../src/roots3and4.cpp ../src/roots3and4.h)    
# target_include_directories(allquadrics PRIVATE include)
# target_include_directories(allquadrics PRIVATE lib/jama125)
# target_include_directories(allquadrics PRIVATE lib/tnt_126)
# target_include_directories(allquadrics PRIVATE lib/CLAPACK/include)
# target_link_libraries(allquadrics PRIVATE lapack blas)

# add_executable(dans_quadric_fit ../src/dans_quadric_fit.cpp)
# target_include_directories(dans_quadric_fit PRIVATE include)
# target_link_libraries(dans_quadric_fit PRIVATE allquadrics)