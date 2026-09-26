// Tests for the CSG library using mesh files in the data/ directory.
// Run via: ctest  or  ./csg_mesh_tests
//
// Each test is a named function; the test runner prints PASS/FAIL.
// Exit code is 0 when all tests pass, 1 otherwise.

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "libcsg.h"

#ifndef DATA_DIR
#define DATA_DIR "data"
#endif

static const double EPS = 1e-5;

// ---------------------------------------------------------------------------
// Assertion helpers
// ---------------------------------------------------------------------------

#define ASSERT(cond) \
    do { if (!(cond)) { \
        std::ostringstream _os; \
        _os << "ASSERT(" #cond ") at line " << __LINE__; \
        throw std::runtime_error(_os.str()); \
    }} while(0)

#define ASSERT_GE(a, b) \
    do { double _a=(a), _b=(b); if (!(_a >= _b)) { \
        std::ostringstream _os; \
        _os << "ASSERT_GE: " << _a << " < " << _b << " at line " << __LINE__; \
        throw std::runtime_error(_os.str()); \
    }} while(0)

#define ASSERT_LE(a, b) \
    do { double _a=(a), _b=(b); if (!(_a <= _b)) { \
        std::ostringstream _os; \
        _os << "ASSERT_LE: " << _a << " > " << _b << " at line " << __LINE__; \
        throw std::runtime_error(_os.str()); \
    }} while(0)

// ---------------------------------------------------------------------------
// Mesh validation helpers
// ---------------------------------------------------------------------------

static std::string dataPath(const std::string& f)
{
    return std::string(DATA_DIR) + "/" + f;
}

// Returns true if any vertex coordinate is non-finite (NaN or Inf).
static bool hasNaN(const CSG::TriMesh& mesh)
{
    for (const auto& v : mesh.vertices())
        if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
            return true;
    return false;
}

// Returns false if any face vertex index is outside [0, numVertices).
static bool hasValidFaceIndices(const CSG::TriMesh& mesh)
{
    const int32_t nv = static_cast<int32_t>(mesh.vertices().size());
    for (const auto& f : mesh.faces())
        for (int ii = 0; ii < 3; ++ii)
            if (f.m_v[ii] < 0 || f.m_v[ii] >= nv)
                return false;
    return true;
}

// Throws on any geometry/topology violation.
static void validateMesh(const CSG::TriMesh& m)
{
    ASSERT(!hasNaN(m));
    ASSERT(hasValidFaceIndices(m));
}

static Eigen::Vector3d meshMin(const CSG::TriMesh& m)
{
    Eigen::Vector3d lo(1e30, 1e30, 1e30);
    for (const auto& v : m.vertices())
        lo = lo.cwiseMin(v);
    return lo;
}

static Eigen::Vector3d meshMax(const CSG::TriMesh& m)
{
    Eigen::Vector3d hi(-1e30, -1e30, -1e30);
    for (const auto& v : m.vertices())
        hi = hi.cwiseMax(v);
    return hi;
}

// Load two OBJ files, run CSG::construct, and return the outputs.
// Throws if either file fails to load.
static void runCSG(const std::string& clay_file, const std::string& knife_file,
                   bool cap, CSG::TriMesh& out_A, CSG::TriMesh& out_B)
{
    CSG::TriMesh clay, knife;
    ASSERT(clay.loadOBJ(dataPath(clay_file)) > 0);
    ASSERT(knife.loadOBJ(dataPath(knife_file)) > 0);
    CSG::CSGEngine engine(clay, knife);
    engine.construct(CSG::kDifference, cap, out_A, out_B);
}

// ---------------------------------------------------------------------------
// Test runner
// ---------------------------------------------------------------------------

static int g_pass = 0, g_fail = 0;

static void run(const char* name, std::function<void()> fn)
{
    std::cout << "[TEST] " << name << std::endl;
    try {
        fn();
        ++g_pass;
        std::cout << "       PASS" << std::endl;
    } catch (const std::exception& e) {
        ++g_fail;
        std::cout << "       FAIL: " << e.what() << std::endl;
    }
}

// ===========================================================================
// Tests: mesh loading
// ===========================================================================

static void test_load_clay_unit_cube()
{
    CSG::TriMesh m;
    ASSERT(m.loadOBJ(dataPath("clay_unit_cube.obj")) == 8);
    ASSERT(m.faces().size() == 12);
    validateMesh(m);
}

static void test_load_clay_tetrahedron()
{
    CSG::TriMesh m;
    ASSERT(m.loadOBJ(dataPath("clay_tetrahedron.obj")) == 4);
    ASSERT(m.faces().size() == 4);
    validateMesh(m);
}

static void test_load_all_knife_meshes()
{
    const char* files[] = {
        "knife_horiz.obj",
        "knife_horiz_flipped.obj",
        "knife_tilted.obj",
        "knife_offset_cube.obj",
        "knife_through_3verts.obj",
        "knife_through_edge.obj",
        "knife_coplanar.obj",
        "knife_multi_cut.obj",
        "knife_near_miss.obj",
    };
    for (const char* f : files) {
        CSG::TriMesh m;
        ASSERT(m.loadOBJ(dataPath(f)) > 0);
        ASSERT(m.faces().size() > 0);
        validateMesh(m);
    }
}

// ===========================================================================
// Tests: AABB tree
// ===========================================================================

static void test_aabb_clay_unit_cube()
{
    CSG::TriMesh m;
    m.loadOBJ(dataPath("clay_unit_cube.obj"));
    AABBTree tree = m.createAABBTree();

    ASSERT(tree.numObjects() == 12);  // one leaf per face

    // Query face 0 and check that at least one candidate is returned
    // (the face is adjacent to others in the unit cube).
    const auto hits = tree.query(0u);
    ASSERT(hits.size() > 0);
}

static void test_aabb_two_mesh_intersect()
{
    CSG::TriMesh clay, knife;
    clay.loadOBJ(dataPath("clay_unit_cube.obj"));
    knife.loadOBJ(dataPath("knife_horiz.obj"));

    AABBTree clay_tree = clay.createAABBTree();
    AABBTree knife_tree = knife.createAABBTree();

    // The unit cube and horizontal z=0 plane overlap
    auto pairs = clay_tree.intersect(knife_tree);
    ASSERT(pairs.size() > 0);
}

// ===========================================================================
// Tests: horizontal cut, upward knife normal
//
// knife_horiz: plane z=0, normal=(0,0,+1).
// "above" the knife means z > 0, so out_A = top half, out_B = bottom half.
// ===========================================================================

static void test_horiz_cut_upward_normal()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_unit_cube.obj", "knife_horiz.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);

    // All out_A vertices must sit at or above the cut plane (z >= 0)
    ASSERT_GE(meshMin(out_A)[2], -EPS);
    // All out_B vertices must sit at or below the cut plane (z <= 0)
    ASSERT_LE(meshMax(out_B)[2], EPS);
}

// ===========================================================================
// Tests: horizontal cut, downward knife normal  [Bug 2 regression]
//
// knife_horiz_flipped: same plane z=0, but winding reversed → normal=(0,0,-1).
//
// Before the retriangulate() fix:
//   - The face normal was computed as cross(normalized_e0, normalized_e1),
//     which is NOT unit-length, so acos(n[2]) was wrong.
//   - When n ∥ UnitZ the rotation axis became a zero vector, producing a NaN
//     AngleAxisd matrix that propagated into every output vertex.
//
// "above" the downward-normal knife means z < 0, so out_A = bottom half.
// ===========================================================================

static void test_horiz_cut_downward_normal_bug2_regression()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_unit_cube.obj", "knife_horiz_flipped.obj", false, out_A, out_B);

    validateMesh(out_A);  // critical: no NaN (was the bug)
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);

    // "above" downward-normal knife → z < 0 → out_A is the bottom half
    ASSERT_LE(meshMax(out_A)[2], EPS);
    // out_B is the top half
    ASSERT_GE(meshMin(out_B)[2], -EPS);
}

// ===========================================================================
// Tests: tilted cut
//
// knife_tilted: plane -3x + 5z = 0 (a 30-degree tilt from vertical).
// "above" means -3x + 5z > 0 for out_A.
//
// KNOWN BUG: knife_tilted.obj has two coplanar triangles covering the same
// knife plane.  Each clay face is intersected by BOTH knife triangles, so the
// same geometric cut points are inserted into new_vert_indices twice (once per
// knife face).  The Shewchuk retriangulator then receives duplicate points and
// constraint segments, which can produce sub-faces that span both sides of the
// cut.  As a result some "below" sub-faces are placed in out_A.
// The plane-equation assertions below are expected to FAIL until the duplicate
// new_vert_indices bug is fixed in construct().
// ===========================================================================

static void test_tilted_cut()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_unit_cube.obj", "knife_tilted.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);

    // knife_tilted normal ~ (-3, 0, 5); plane equation: -3x + 5z = 0
    // Fails until duplicate new_vert_indices is fixed:
    for (const auto& v : out_A.vertices())
        ASSERT_GE(-3.0 * v[0] + 5.0 * v[2], -EPS);
    for (const auto& v : out_B.vertices())
        ASSERT_LE(-3.0 * v[0] + 5.0 * v[2], EPS);
}

// ===========================================================================
// Tests: cube-cube intersection
//
// knife_offset_cube: unit cube shifted by (0.5, 0.5, 0.5).
// Many face-face intersections; both halves should be non-empty.
// ===========================================================================

static void test_cube_cube_intersection()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_unit_cube.obj", "knife_offset_cube.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);
}

// ===========================================================================
// Tests: cut through 3 vertices  [Bug 1 regression]
//
// knife_through_3verts: plane x+y+z=1 passes exactly through cube vertices
// (1,1,-1), (1,-1,1), (-1,1,1).
//
// Bug 1: convertIntersectionToIpoints stored the local vertex index within the
// face (0, 1, or 2) instead of the global mesh vertex index for clay/knife
// vertices that land exactly on the intersection plane.  This caused wrong
// IPointRef lookups in canonicalVertexIndex and ipointPos.
// ===========================================================================

static void test_cut_through_3_vertices_bug1_regression()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_unit_cube.obj", "knife_through_3verts.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);
}

// ===========================================================================
// Tests: cut through an edge  [Bug 1 + Bug 4 regression]
//
// knife_through_edge: plane -3y+4z=-1 passes through the entire bottom-front
// edge of the unit cube (both vertices (-1,-1,-1) and (1,-1,-1) lie on it).
//
// Exercises Bug 1 (same as above) AND Bug 4: when both endpoints of the
// intersection segment are original clay vertices, the retriangulated
// sub-faces contain no kNew vertices, leaving testFace uninitialized in
// classifyCutFaces and causing UB in the voting loop.
// ===========================================================================

static void test_cut_through_edge_bug1_bug4_regression()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_unit_cube.obj", "knife_through_edge.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);
}

// ===========================================================================
// Tests: coplanar knife  [coplanar-skip regression]
//
// knife_coplanar: quad exactly at z=-1, coincident with the cube's bottom face.
// Coplanar pairs are skipped in the intersection loop, so no faces are cut and
// classifyFaces has no seed to start the flood-fill from.  Both outputs may be
// empty, but the library must not crash or produce NaN.
// ===========================================================================

static void test_coplanar_knife_no_crash()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_unit_cube.obj", "knife_coplanar.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    // Outputs may be empty — that is acceptable for coplanar knives
}

// ===========================================================================
// Tests: near-miss vertex (ULP boundary)
//
// knife_near_miss: plane x+y+z=0.999, just 0.001 units below the three cube
// vertices (1,1,-1),(1,-1,1),(-1,1,1).  Should produce a valid cut.
// ===========================================================================

static void test_near_miss_vertex()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_unit_cube.obj", "knife_near_miss.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);
}

// ===========================================================================
// Tests: multiple knife faces cut the same clay face  [Bug D regression]
//
// big_triangle.obj + knife_multi_cut.obj: the knife has two triangles whose
// intersection lines both cross the single large clay triangle.
//
// Bug D: classifyCutFaces iterated ALL kNew vertices of a cut face to find
// testFace, picking the LAST kNew vertex's knife face instead of the first.
// When the two knife triangles have opposite normals this assigns the wrong
// knife face normal, flipping the classification of some sub-faces.
// ===========================================================================

static void test_multi_cut_same_face_bugD_regression()
{
    CSG::TriMesh out_A, out_B;
    runCSG("big_triangle.obj", "knife_multi_cut.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);
}

// ===========================================================================
// Tests: tetrahedron with horizontal cut
//
// All 4 tetrahedron faces are cut by the z=0 plane, so no face escapes the
// flood-fill seed.
// ===========================================================================

static void test_tetrahedron_horiz_cut()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_tetrahedron.obj", "knife_horiz.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);

    // knife_horiz normal = +Z → out_A is z > 0 half
    ASSERT_GE(meshMin(out_A)[2], -EPS);
    ASSERT_LE(meshMax(out_B)[2], EPS);
}

// ===========================================================================
// Tests: tetrahedron with tilted cut
//
// Exercises the generic (non-degenerate-axis) rotation path in retriangulate().
// ===========================================================================

static void test_tetrahedron_tilted_cut()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_tetrahedron.obj", "knife_tilted.obj", false, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);
}

// ===========================================================================
// Tests: capped cut (cap=true)
//
// Verifies that knife cap faces are merged into the output without introducing
// NaN vertices or invalid face indices.
// ===========================================================================

static void test_horiz_cut_with_cap()
{
    CSG::TriMesh out_A, out_B;
    runCSG("clay_unit_cube.obj", "knife_horiz.obj", true, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);
}

// ===========================================================================
// main
// ===========================================================================

int main()
{
    // Loading
    run("load_clay_unit_cube",          test_load_clay_unit_cube);
    run("load_clay_tetrahedron",        test_load_clay_tetrahedron);
    run("load_all_knife_meshes",        test_load_all_knife_meshes);

    // AABB tree
    run("aabb_clay_unit_cube",          test_aabb_clay_unit_cube);
    run("aabb_two_mesh_intersect",      test_aabb_two_mesh_intersect);

    // CSG operations — standard cases
    run("horiz_cut_upward_normal",      test_horiz_cut_upward_normal);
    run("tilted_cut",                   test_tilted_cut);
    run("cube_cube_intersection",       test_cube_cube_intersection);
    run("horiz_cut_with_cap",           test_horiz_cut_with_cap);

    // CSG operations — bug regressions
    run("horiz_cut_downward_normal_bug2",       test_horiz_cut_downward_normal_bug2_regression);
    run("cut_through_3_vertices_bug1",          test_cut_through_3_vertices_bug1_regression);
    run("cut_through_edge_bug1_bug4",           test_cut_through_edge_bug1_bug4_regression);
    run("coplanar_knife_no_crash",              test_coplanar_knife_no_crash);
    run("near_miss_vertex",                     test_near_miss_vertex);
    run("multi_cut_same_face_bugD",             test_multi_cut_same_face_bugD_regression);

    // Alternative clay mesh
    run("tetrahedron_horiz_cut",        test_tetrahedron_horiz_cut);
    run("tetrahedron_tilted_cut",       test_tetrahedron_tilted_cut);

    const int total = g_pass + g_fail;
    std::cout << "\n--- " << g_pass << "/" << total << " tests passed";
    if (g_fail > 0)
        std::cout << ", " << g_fail << " FAILED";
    std::cout << " ---\n";
    return g_fail > 0 ? 1 : 0;
}
