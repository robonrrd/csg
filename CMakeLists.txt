cmake_minimum_required(VERSION 3.14)

project (csg)
set (CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED NO_MODULE)


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
target_link_libraries(csg_test ${PROJECT_LINK_LIBS} Eigen3::Eigen)
# AABB_TEST
add_executable(aabb_test src/aabb_test.cpp)
target_link_libraries(aabb_test ${PROJECT_LINK_LIBS} Eigen3::Eigen)


enable_testing()
add_subdirectory(tests)


## SWIG for Python bindings (pycsg) — optional; skipped if SWIG is not found
##
find_package(SWIG)
if(SWIG_FOUND)
   find_package(Python COMPONENTS Interpreter Development)
   if(Python_FOUND)
      include(${SWIG_USE_FILE})

      # Find Python user site-packages directory
      execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import site; print(site.getusersitepackages())"
        OUTPUT_VARIABLE PYTHON_USER_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )

      set(CMAKE_SWIG_FLAGS "")
      set_source_files_properties(swig/pycsg.i PROPERTIES CPLUSPLUS ON)
      set_source_files_properties(swig/pycsg.i PROPERTIES SWIG_FLAGS "")
      swig_add_library(pycsg
                       TYPE SHARED
                       LANGUAGE python
                       SOURCES swig/pycsg.i src/trimesh.cpp src/aabb.cpp src/triangle.c src/libcsg.cpp)
      target_link_libraries(pycsg PRIVATE Python::Module)
      set_target_properties(pycsg PROPERTIES SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)

      # Install SWIG generated DSOs for Python integration
      install(TARGETS ${SWIG_MODULE_pycsg_REAL_NAME}
        DESTINATION ${PYTHON_USER_SITE_PACKAGES})

      install(FILES ${CMAKE_BINARY_DIR}/pycsg.py
        DESTINATION ${PYTHON_USER_SITE_PACKAGES})
   else()
      message(STATUS "Python development headers not found — Python bindings (pycsg) will not be built")
   endif()
else()
   message(STATUS "SWIG not found — Python bindings (pycsg) will not be built")
endif()
