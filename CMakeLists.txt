cmake_minimum_required(VERSION 3.0x)
project (csg)
set (CMAKE_CXX_STANDARD 11)

set(EIGEN3_INCLUDE_DIR /usr/local/include/eigen3)


# Generate the shared library
#
add_library(csg SHARED
   src/aabb.cpp
   src/trimesh.cpp
   src/libcsg.cpp
   src/triangle.c
)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g")

# Set the location for library installation
#install(TARGETS csg DESTINATION lib)

set ( PROJECT_LINK_LIBS csg)
link_directories( build lib )
include_directories( include ${EIGEN3_INCLUDE_DIR} )


# CSG_TEST
add_executable(csg_test src/csg_test.cpp)
target_link_libraries(csg_test ${PROJECT_LINK_LIBS})

# AABB_TEST
add_executable(aabb_test src/aabb_test.cpp)
target_link_libraries(aabb_test ${PROJECT_LINK_LIBS})


## SWIG for Python bindings (pycsg)
##
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(swig/pycsg.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(swig/pycsg.i PROPERTIES SWIG_FLAGS "")
swig_add_library(pycsg
                 TYPE SHARED
                 LANGUAGE python
                 SOURCES swig/pycsg.i src/trimesh.cpp src/aabb.cpp src/triangle.c src/libcsg.cpp)
swig_link_libraries(pycsg ${PYTHON_LIBRARIES})


# Install SWIG generated DSOs for Python integration
install(TARGETS ${SWIG_MODULE_pycsg_REAL_NAME}
  DESTINATION $ENV{HOME}/.local/lib/python3.6/site-packages)

install(FILES build/pycsg.py
  DESTINATION  $ENV{HOME}/.local/lib/python3.6/site-packages)
