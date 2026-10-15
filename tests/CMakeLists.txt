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
set_tests_properties(csg_tests PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:csg>:$ENV{LD_LIBRARY_PATH}"
)

# `make tests` builds and immediately runs the tests, with their output visible.
add_custom_target(tests
    COMMAND ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:csg>:$ENV{LD_LIBRARY_PATH}" $<TARGET_FILE:csg_tests>
    DEPENDS csg_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running CSG tests"
    VERBATIM
)
