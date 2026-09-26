add_executable(csg_mesh_tests csg_mesh_tests.cpp)

target_link_libraries(csg_mesh_tests ${PROJECT_LINK_LIBS} Eigen3::Eigen)

# DATA_DIR points at the source tree's data/ folder so tests find OBJ files
# regardless of the build directory.
target_compile_definitions(csg_mesh_tests
    PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data"
)

# CTest integration — `ctest` or `cmake --test .` will pick this up.
add_test(
    NAME csg_mesh_tests
    COMMAND csg_mesh_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# `make tests` builds and immediately runs the test binary.
add_custom_target(tests
    COMMAND csg_mesh_tests
    DEPENDS csg_mesh_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running CSG mesh tests"
    VERBATIM
)
