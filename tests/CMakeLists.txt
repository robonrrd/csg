# All of the tests are in one program, csg_tests (csg_tests.cpp, csg_tests.h).
# It also has two command-line tools: `csg_tests cut CLAY KNIFE` and
# `csg_tests aabb OBJ_FILE`.
add_executable(csg_tests csg_tests.cpp csg_tests.h)

target_link_libraries(csg_tests ${PROJECT_LINK_LIBS} Eigen3::Eigen)

# DATA_DIR points at the source tree's data/ folder so tests find OBJ files
# regardless of the build directory.
target_compile_definitions(csg_tests
    PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data"
)

# The tests link against libcsg.so, and LD_LIBRARY_PATH (e.g. /usr/local/lib)
# takes precedence over the binary's build-tree RUNPATH, so an older installed
# libcsg.so would be tested instead of the one just built.  Put the freshly built
# library first.

# CTest integration: `ctest` or `make test` will pick this up.
add_test(
    NAME csg_tests
    COMMAND csg_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
# (Some of the tests are for bugs that hang; they give up by themselves, but this is a
# limit in case something else does.)
set_tests_properties(csg_tests PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:csg>:$ENV{LD_LIBRARY_PATH}"
    TIMEOUT 300
)

# `make tests` builds and immediately runs the tests, with their output visible.
add_custom_target(tests
    COMMAND ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:csg>:$ENV{LD_LIBRARY_PATH}" $<TARGET_FILE:csg_tests>
    DEPENDS csg_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running CSG tests"
    VERBATIM
)


# ---------------------------------------------------------------------------
# Checks of how the sources build.  A test with WILL_FAIL documents a bug that this
# code still has: when the bug is fixed the test fails, and WILL_FAIL should be
# removed.
# ---------------------------------------------------------------------------

# The library must not fall off the end of a non-void function. With assertions
# enabled (Debug) `assert(false)` hides that, so check a build without them.
# (Fails: the free function otherMesh() in libcsg.cpp ends in assert(false).)
add_test(
    NAME csg_no_missing_return
    COMMAND ${CMAKE_CXX_COMPILER} -std=c++11 -DNDEBUG -Wreturn-type -Werror=return-type
            -I${CMAKE_SOURCE_DIR}/include -I${EIGEN3_INCLUDE_DIR}
            -c ${CMAKE_SOURCE_DIR}/src/libcsg.cpp -o ${CMAKE_BINARY_DIR}/csg_no_missing_return.o
)
set_tests_properties(csg_no_missing_return PROPERTIES WILL_FAIL TRUE)

# trimesh.cpp and aabb.cpp must link on their own in an optimized build: a missing Eigen
# header only shows up there (the template instantiations come from other files otherwise).
# (Fails: trimesh.cpp uses cross() without including <Eigen/Geometry>.)
file(WRITE ${CMAKE_BINARY_DIR}/csg_link_trimesh.cpp
"// Links trimesh.cpp and aabb.cpp on their own, optimized.
#include <iostream>
#include \"trimesh.h\"

int main()
{
   const std::vector<double> v = {0, 0, 0,  1, 0, 0,  0, 1, 0};
   const std::vector<unsigned int> f = {0, 1, 2};
   CSG::TriMesh mesh(v, f);
   const Eigen::Vector3d n = mesh.faceNormal(0);
   std::cout << mesh << n.transpose() << std::endl;
   return (n - Eigen::Vector3d::UnitZ()).norm() < 1e-12 ? 0 : 1;
}
")
add_test(
    NAME csg_trimesh_links_standalone
    COMMAND ${CMAKE_CXX_COMPILER} -std=c++11 -O2 -DNDEBUG
            -I${CMAKE_SOURCE_DIR}/include -I${EIGEN3_INCLUDE_DIR}
            ${CMAKE_BINARY_DIR}/csg_link_trimesh.cpp
            ${CMAKE_SOURCE_DIR}/src/trimesh.cpp ${CMAKE_SOURCE_DIR}/src/aabb.cpp
            -o ${CMAKE_BINARY_DIR}/csg_link_trimesh
)
set_tests_properties(csg_trimesh_links_standalone PROPERTIES WILL_FAIL TRUE)

# The headers must have include guards: a source file that reaches one twice must compile
file(WRITE ${CMAKE_BINARY_DIR}/csg_header_twice.cpp
"#include \"libcsg.h\"\n#include \"libcsg.h\"\n#include \"trimesh.h\"\n#include \"aabb.h\"\nint main() { return 0; }\n")
add_test(
    NAME csg_headers_have_include_guards
    COMMAND ${CMAKE_CXX_COMPILER} -std=c++11 -fsyntax-only
            -I${CMAKE_SOURCE_DIR}/include -I${EIGEN3_INCLUDE_DIR}
            ${CMAKE_BINARY_DIR}/csg_header_twice.cpp
)
