cmake_minimum_required(VERSION 2.6)
PROJECT (GraphApproach)

set (CMAKE_CXX_STANDARD 14)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories("/usr/local/include/igraph")
#include_directories("/usr/local/include/eigen3")
include_directories("src")



add_executable(GraphApproach "src/GraphApproach.cpp")

target_link_libraries(GraphApproach "/usr/local/lib/libigraph.a"  ${ITK_LIBRARIES}  ${VTK_LIBRARIES})
#target_link_libraries(MultilevelImageSegmentation  ${ITK_LIBRARIES}  ${VTK_LIBRARIES})


