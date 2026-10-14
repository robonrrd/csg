add_executable(csg_mesh_tests csg_mesh_tests.cpp)

target_link_libraries(csg_mesh_tests ${PROJECT_LINK_LIBS} Eigen3::Eigen)

# DATA_DIR points at the source tree's data/ folder so tests find OBJ files
# regardless of the build directory.
target_compile_definitions(csg_mesh_tests
    PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data"
)

# The tests link against libcsg.so, and LD_LIBRARY_PATH (e.g. /usr/local/lib)
# takes precedence over the binaries' build-tree RUNPATH, so an older installed
# libcsg.so would be tested instead of the one just built.  Put the freshly built
# library first.

# CTest integration — `ctest` or `cmake --test .` will pick this up.
add_test(
    NAME csg_mesh_tests
    COMMAND csg_mesh_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(csg_mesh_tests PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:csg>:$ENV{LD_LIBRARY_PATH}"
)

# `make tests` builds and immediately runs the test binary.
add_custom_target(tests
    COMMAND ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:csg>:$ENV{LD_LIBRARY_PATH}" $<TARGET_FILE:csg_mesh_tests>
    DEPENDS csg_mesh_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running CSG mesh tests"
    VERBATIM
)

# Regression tests for the bugs tracked in TODO.md (data files in data/bugs/).
add_executable(csg_bug_tests csg_bug_tests.cpp)
target_link_libraries(csg_bug_tests ${PROJECT_LINK_LIBS} Eigen3::Eigen)
target_compile_definitions(csg_bug_tests
    PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data"
)
add_test(
    NAME csg_bug_tests
    COMMAND csg_bug_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(csg_bug_tests PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:csg>:$ENV{LD_LIBRARY_PATH}"
)
add_custom_target(bug_tests
    COMMAND ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:csg>:$ENV{LD_LIBRARY_PATH}" $<TARGET_FILE:csg_bug_tests>
    DEPENDS csg_bug_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running CSG bug regression tests"
    VERBATIM
)
