cmake_minimum_required(VERSION 3.0x)
project (csg)
set (CMAKE_CXX_STANDARD 11)

set( EIGEN3_INCLUDE_DIR /usr/include/eigen3 )

# the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*.cpp")


# Generate the shared library
#
add_library(csg SHARED
   src/aabb.cpp
   src/triangle.cpp
   src/trimesh.cpp
   src/libcsg.cpp
)

# Set the location for library installation
install(TARGETS csg DESTINATION lib)

set ( PROJECT_LINK_LIBS csg)
link_directories( build lib )
include_directories( include ${EIGEN3_INCLUDE_DIR} )

# CSG_TEST
add_executable(csg_test src/csg_test.cpp)
target_link_libraries(csg_test ${PROJECT_LINK_LIBS})

# AABB_TEST
add_executable(aabb_test src/aabb_test.cpp)
target_link_libraries(aabb_test ${PROJECT_LINK_LIBS})

