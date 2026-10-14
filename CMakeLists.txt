cmake_minimum_required(VERSION 3.14)

project (csg)
include(GNUInstallDirs)
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

      # Python install location. Defaults to the invoking user's site-packages
      # (resolved at configure time, so run cmake as the user, not via sudo);
      # override with -DPYCSG_INSTALL_DIR=...
      execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import site; print(site.getusersitepackages())"
        OUTPUT_VARIABLE _PYCSG_USER_SITE
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
      set(PYCSG_INSTALL_DIR "${_PYCSG_USER_SITE}" CACHE PATH
          "Directory where the pycsg Python module is installed")
      set(PYTHON_USER_SITE_PACKAGES ${PYCSG_INSTALL_DIR})

      # Place the generated pycsg.py in the build root so the install path below
      # is predictable regardless of CMake version.
      set(CMAKE_SWIG_OUTDIR ${CMAKE_BINARY_DIR})

      set(CMAKE_SWIG_FLAGS "")
      set_source_files_properties(swig/pycsg.i PROPERTIES CPLUSPLUS ON)
      set_source_files_properties(swig/pycsg.i PROPERTIES SWIG_FLAGS "")

      # Only wrap pycsg.i. link against the already-built CSG library
      swig_add_library(pycsg
                       TYPE SHARED
                       LANGUAGE python
                       SOURCES swig/pycsg.i)

      # Provide include paths so SWIG can parse the CSG and Eigen headers when
      # it processes pycsg.i (SWIG_USE_TARGET_INCLUDE_DIRECTORIES uses these).
      target_include_directories(pycsg PRIVATE include ${EIGEN3_INCLUDE_DIR})
      target_link_libraries(pycsg PRIVATE csg Python::Module Eigen3::Eigen)
      set_target_properties(pycsg PROPERTIES SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)
      # Let the installed module find the installed libcsg.so without ldconfig
      set_target_properties(pycsg PROPERTIES
         INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")

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


## C++ install (library + headers). Declared last so that "make install" installs
## the Python module first: if this step fails for lack of root, Python is
## already installed.
##
# Install the shared library and public headers for use from C++
# (to /usr/local by default; override with -DCMAKE_INSTALL_PREFIX=...)
install(TARGETS csg
   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES include/libcsg.h include/trimesh.h include/aabb.h
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
