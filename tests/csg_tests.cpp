// Tests for the CSG library.
//
//     csg_tests [NAME ...]          run the tests (only those whose names contain
//                                   one of the NAMEs, if any are given)
//     csg_tests cut CLAY KNIFE [nocap]
//                                   cut CLAY (an OBJ file) with KNIFE and write
//                                   the two halves to A.obj and B.obj
//     csg_tests aabb OBJ_FILE       build an AABB tree of a mesh and print stats
//
// Run via ctest or `make tests`.  Every test is a named function; the runner
// prints PASS/FAIL, and the exit code is 0 when every test behaves as expected.
// The runner, assertions and mesh-validation helpers are in csg_tests.h.
//
// Bugs that have been fixed have regression tests here.  Tests for bugs that
// have not been fixed yet are registered with runKnownBug(): they are expected to
// fail, and the runner reports an unexpected pass so the registration can be
// switched to run() as soon as the bug is fixed.

#include <chrono>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <random>
#include <set>

#include "csg_tests.h"

using meshcheck::analyze;
using meshcheck::Report;


// ===========================================================================
// Helpers
// ===========================================================================

// All the faces of A and of B in one mesh.  Without a cap the two halves of a cut
// clay mesh are the clay surface split along the cut, so together they must form
// exactly the closed clay surface again.
static CSG::TriMesh combine(const CSG::TriMesh& A, const CSG::TriMesh& B)
{
   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> verts;
   std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> faces;
   for (const CSG::TriMesh* m : {&A, &B})
   {
      const int base = static_cast<int>(verts.size());
      verts.insert(verts.end(), m->vertices().begin(), m->vertices().end());
      for (const auto& f : m->faces())
         faces.push_back(Eigen::Vector3i(base + f.m_v[0], base + f.m_v[1], base + f.m_v[2]));
   }
   return CSG::TriMesh(verts, faces);
}

static bool hasVertexAt(const CSG::TriMesh& m, const Eigen::Vector3d& p)
{
   for (const auto& v : m.vertices())
      if ((v - p).norm() < 1e-9)
         return true;
   return false;
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

static void test_aabb_tree_is_valid_for_a_large_mesh()
{
    CSG::TriMesh m;
    ASSERT(m.loadOBJ(dataPath("bunny.obj")) > 0);
    AABBTree tree = m.createAABBTree();

    ASSERT(tree.numObjects() == m.faces().size());
    tree.validate();  // throws if the tree's structure or bounds are wrong
}

// Random points, ten of which have a near-duplicate (closer than the sphere
// radius) among the others; querying a point must find exactly its duplicate.
static void test_aabb_point_welding()
{
    const uint32_t num_points = 1000;
    const uint32_t num_duplicates = 10;
    const double extent = 10.0;
    const double radius = 1e-06;

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uniform(0.0, extent);

    AABBTree tree(0.05, num_points);
    std::vector<Eigen::Vector3d> originals;
    for (uint32_t ii = 0; ii < num_points - num_duplicates; ++ii)
    {
        const Eigen::Vector3d pt(uniform(rng), uniform(rng), uniform(rng));
        tree.addSphere(ii, pt, radius);
        originals.push_back(pt);
    }
    const uint32_t first_duplicate = num_points - num_duplicates;
    for (uint32_t ii = 0; ii < num_duplicates; ++ii)
    {
        // within half a radius of point ii
        const Eigen::Vector3d pt = originals[ii] + Eigen::Vector3d(0.5 * radius, 0, 0);
        tree.addSphere(first_duplicate + ii, pt, radius);
    }

    tree.validate();
    ASSERT(tree.numObjects() == num_points);

    // A point with no duplicate finds nothing...
    ASSERT(tree.query(200u).empty());

    // ...and one with a duplicate finds it (and the duplicate finds the point)
    for (uint32_t ii = 0; ii < num_duplicates; ++ii)
    {
        const auto hits = tree.query(first_duplicate + ii);
        ASSERT(hits.size() == 1);
        ASSERT(hits[0] == ii);

        const auto back = tree.query(ii);
        ASSERT(back.size() == 1);
        ASSERT(back[0] == first_duplicate + ii);
    }
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
    constructFiles("clay_unit_cube.obj", "knife_horiz.obj", false, out_A, out_B);

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
    constructFiles("clay_unit_cube.obj", "knife_horiz_flipped.obj", false, out_A, out_B);

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
    constructFiles("clay_unit_cube.obj", "knife_tilted.obj", false, out_A, out_B);

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
    constructFiles("clay_unit_cube.obj", "knife_offset_cube.obj", false, out_A, out_B);

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
    constructFiles("clay_unit_cube.obj", "knife_through_3verts.obj", false, out_A, out_B);

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
    constructFiles("clay_unit_cube.obj", "knife_through_edge.obj", false, out_A, out_B);

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
    constructFiles("clay_unit_cube.obj", "knife_coplanar.obj", false, out_A, out_B);

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
    constructFiles("clay_unit_cube.obj", "knife_near_miss.obj", false, out_A, out_B);

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
    constructFiles("big_triangle.obj", "knife_multi_cut.obj", false, out_A, out_B);

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
    constructFiles("clay_tetrahedron.obj", "knife_horiz.obj", false, out_A, out_B);

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
    constructFiles("clay_tetrahedron.obj", "knife_tilted.obj", false, out_A, out_B);

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
    constructFiles("clay_unit_cube.obj", "knife_horiz.obj", true, out_A, out_B);

    validateMesh(out_A);
    validateMesh(out_B);
    ASSERT(out_A.faces().size() > 0);
    ASSERT(out_B.faces().size() > 0);
}


// ===========================================================================
// BUG-9: cut faces with no new vertex are dropped from both outputs
// (fixed: faces are now classified by flood-filling regions bounded by the cuts)
//
// The knife plane x+y+z=1 passes exactly through the cube vertices (1,0,0),
// (0,1,0) and (0,0,1).  The three cube faces that only touch the plane at one
// of those vertices (x=1, y=1, z=1) must end up in the output, each in full.
// ===========================================================================

static void test_bug09_touching_faces_are_kept()
{
   CSG::TriMesh A, B;
   constructFiles("bugs/bug09_clay.obj", "bugs/bug09_knife.obj", false, A, B);

   const Eigen::Vector3d axes[3] = {Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY(),
                                    Eigen::Vector3d::UnitZ()};
   for (int ax = 0; ax < 3; ++ax)
   {
      // Each face of the unit cube has area 1 and is split between (or wholly
      // inside) A and B; none of it may vanish.
      for (double d : {0.0, 1.0})
      {
         const double total = meshcheck::areaOnPlane(A, axes[ax], d) +
                              meshcheck::areaOnPlane(B, axes[ax], d);
         if (std::fabs(total - 1.0) > 1e-9)
         {
            std::ostringstream os;
            os << "cube face " << "xyz"[ax] << "=" << d << " has area " << total
               << " in the outputs, expected 1";
            throw std::runtime_error(os.str());
         }
      }
   }
   // Together the two pieces cover the whole cube surface.
   ASSERT_NEAR(analyze(A).area + analyze(B).area, 6.0, 1e-9);
}

// ===========================================================================
// BUG-10: classifyCutFaces() looked up the cutting face with the wrong index
//
// The knife is a V-shaped fold (two triangles, normals pointing up) whose
// crease lies along the bottom-left edge of the unit cube, so some intersection
// points are existing cube vertices.  Those points shift the indices of
// m_newPoints relative to m_newPointPositions, and the wrong knife triangle was
// used to decide which side of the knife a sub-face of the top face is on.
//
// Surface: z = max(-2x, 2x).  On the top face (z=1) the region above the knife
// is |x| < 0.5, i.e. exactly half of the face.
// ===========================================================================

static void test_bug10_cutting_face_lookup()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeTriangles({
      {0, -6, 0,   0, 6, 0,   -6, 0, 12},   // x < 0 half of the fold
      {0, -6, 0,   6, 0, 12,  0, 6, 0}});   // x > 0 half of the fold
   CSG::TriMesh A, B;
   construct(clay, knife, false, A, B);

   const double topA = meshcheck::areaOnPlane(A, Eigen::Vector3d::UnitZ(), 1.0);
   const double topB = meshcheck::areaOnPlane(B, Eigen::Vector3d::UnitZ(), 1.0);
   ASSERT_NEAR(topA, 0.5, 1e-9);
   ASSERT_NEAR(topB, 0.5, 1e-9);
}

// ===========================================================================
// Mesh components that the knife never touches must not vanish
//
// Three unit boxes stacked in z (z in [-3,-2], [0,1] and [2,3]); the knife is
// the plane z=0.5 covering only the middle box.  The middle box is split in
// half; the top box lies above the knife (so belongs in A), the bottom box
// below it (B).  Nothing cuts those two, so the fill cannot reach them from a
// cut and they have to be classified against the knife directly.
// ===========================================================================

static void test_untouched_components_are_kept()
{
   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> verts;
   std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> faces;
   for (double z : {-3.0, 0.0, 2.0})
   {
      const CSG::TriMesh box =
         meshcheck::makeBox(Eigen::Vector3d(0, 0, z), Eigen::Vector3d(1, 1, z + 1));
      const int base = static_cast<int>(verts.size());
      verts.insert(verts.end(), box.vertices().begin(), box.vertices().end());
      for (const auto& f : box.faces())
         faces.push_back(Eigen::Vector3i(base + f.m_v[0], base + f.m_v[1], base + f.m_v[2]));
   }
   const CSG::TriMesh clay(verts, faces);
   const CSG::TriMesh knife = meshcheck::makeTriangles({
      {-1, -1, 0.5,   2, -1, 0.5,   2, 2, 0.5},
      {-1, -1, 0.5,   2, 2, 0.5,   -1, 2, 0.5}});

   CSG::TriMesh A, B;
   construct(clay, knife, false, A, B);

   // top box (6) + the upper half of the middle box (top 1 + sides 2)
   ASSERT_NEAR(analyze(A).area, 9.0, 1e-9);
   ASSERT_NEAR(analyze(B).area, 9.0, 1e-9);
}

// ===========================================================================
// BUG-53: points that are a vertex of both meshes
//
// The knife is a single triangle whose three vertices are the cube vertices
// (1,0,0), (0,1,0) and (0,0,1), so every point where it meets the cube is a vertex
// of both meshes.  Such a point is stored as the clay's vertex, but the knife's
// own faces refer to it by the knife's index, so classifying (and assembling)
// the knife faces for the cap used to lose the connection between them.
//
// The triangle's normal is (1,1,1), so 'above' is the large piece (volume 5/6)
// and 'below' is the corner (volume 1/6); flipping the triangle swaps them.  With
// a cap both pieces must be closed, consistently oriented solids.
// ===========================================================================

static void checkCornerCut(const CSG::TriMesh& A, const CSG::TriMesh& B,
                           double expected_volume_A, double expected_volume_B)
{
   const Report ra = analyze(A);
   const Report rb = analyze(B);
   ASSERT(ra.closed());
   ASSERT(rb.closed());
   ASSERT_NEAR(ra.volume, expected_volume_A, 1e-9);
   ASSERT_NEAR(rb.volume, expected_volume_B, 1e-9);
}

static void test_bug53_knife_vertices_on_clay_vertices()
{
   const CSG::TriMesh clay =
      meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeTriangles({{1, 0, 0,   0, 1, 0,   0, 0, 1}});
   const CSG::TriMesh knife_flipped = meshcheck::makeTriangles({{1, 0, 0,   0, 0, 1,   0, 1, 0}});

   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);
   checkCornerCut(A, B, 5.0 / 6.0, 1.0 / 6.0);

   construct(clay, knife_flipped, true, A, B);
   checkCornerCut(A, B, 1.0 / 6.0, 5.0 / 6.0);

   // Without a cap the pieces are the clay surface, split.
   construct(clay, knife, false, A, B);
   const Report both = analyze(combine(A, B));
   ASSERT(both.closed());
   ASSERT_NEAR(both.area, 6.0, 1e-9);
}

// ===========================================================================
// BUG-14: constraint segments that cross inside a face
//
// Two knife triangles, in the planes x=0.7 and y=0.3, cross each other along the
// vertical line through (0.7, 0.3).  On the top face of the cube the two cuts are
// therefore two segments that cross at (0.7, 0.3, 1), a point that no triangle
// intersection produced: the triangulator has to create it.  (The point is
// chosen away from the diagonal that splits the top face into two triangles; on
// the diagonal the cuts would only meet at an edge of the triangulation.)  Its
// extra vertex used to be indexed past the end of the vertex table, and the cuts
// have to meet at it so that the faces fit together.
// ===========================================================================

static void test_bug14_crossing_cuts_in_a_face()
{
   const CSG::TriMesh clay =
      meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeTriangles({
      {0.7, -10, -5,   0.7, 10, -5,   0.7, 0, 10},    // x = 0.7
      {-10, 0.3, -5,   10, 0.3, -5,   0, 0.3, 10}});  // y = 0.3
   CSG::TriMesh A, B;
   construct(clay, knife, false, A, B);

   // The crossing point is a vertex of the output
   ASSERT(hasVertexAt(A, Eigen::Vector3d(0.7, 0.3, 1.0)) ||
         hasVertexAt(B, Eigen::Vector3d(0.7, 0.3, 1.0)));

   // Nothing is lost, and every crack is closed: A and B together are the
   // whole, watertight, surface of the cube
   const Report both = analyze(combine(A, B));
   ASSERT(both.closed());
   ASSERT_NEAR(both.area, 6.0, 1e-9);
   ASSERT_NEAR(both.volume, 1.0, 1e-9);
}

// ===========================================================================
// BUG-54: a knife that only touches the clay at a point
//
// The apex of the knife triangle lies exactly on the interior of the cube's top
// face and the rest of the triangle is above the cube, so the two triangles
// intersect in a segment of zero length.  That is not a cut: all of the clay must
// come out, intact, on one side of the knife.
// ===========================================================================

static void test_bug54_point_contact()
{
   const CSG::TriMesh clay =
      meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeTriangles({{0.5, 0.5, 1,   -2, -2, 3,   3, -2, 3}});
   CSG::TriMesh A, B;
   construct(clay, knife, false, A, B);

   const Report both = analyze(combine(A, B));
   ASSERT(both.closed());
   ASSERT_NEAR(both.area, 6.0, 1e-9);
   ASSERT_NEAR(both.volume, 1.0, 1e-9);

   // ...all on one side of the knife (the clay is entirely below it)
   ASSERT(A.faces().empty() != B.faces().empty());
}



// ===========================================================================
// More bugs, from the project's bug tracker
//
// The tests above are for the first bugs found here.  The ones below were written
// for a long list of bugs, BUG-5 to BUG-49, numbered as in the bug tracker,
// TODO.md, which describes every bug and says whether it has been fixed.  The
// data files data/bugs/bugNN_*.obj are named by the same numbers.  BUG-53 and
// BUG-54 above were not on the original list.  Tests that pass here are
// registered with run().  The others document bugs that this code still has: some
// of them crash or hang, so they are registered with runKnownBug(), which runs
// them in a child process.
// ===========================================================================

// ===========================================================================
// BUG-20: the face used to classify a sub-face was an arbitrary one
//
// A valley-shaped knife z = 0.5 + |x - 0.5| (two triangles, normals up) cuts the
// unit cube so that the sub-faces of the front and back faces are bounded by
// intersection lines from BOTH knife triangles.  Whichever triangle's plane is
// used for a sub-face, the answer must agree with "above the valley".  Above the
// valley on the y=0 face is the triangle (0,1)-(1,1)-(0.5,0.5), area 0.25.
// ===========================================================================

static void test_bug20_two_cutting_faces()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeTriangles({
      {0.5, -6, 0.5,   0.5, 6, 0.5,   -5.5, 0, 6.5},   // x < 0.5 half
      {0.5, -6, 0.5,   6.5, 0, 6.5,   0.5, 6, 0.5}});  // x > 0.5 half
   CSG::TriMesh A, B;
   construct(clay, knife, false, A, B);

   struct Expect { Eigen::Vector3d n; double d; double above; };
   const Expect expect[] = {
      {Eigen::Vector3d::UnitY(), 0.0, 0.25}, {Eigen::Vector3d::UnitY(), 1.0, 0.25},
      {Eigen::Vector3d::UnitZ(), 1.0, 1.0},  {Eigen::Vector3d::UnitZ(), 0.0, 0.0},
      {Eigen::Vector3d::UnitX(), 0.0, 0.0},  {Eigen::Vector3d::UnitX(), 1.0, 0.0}};
   for (const Expect& e : expect)
   {
      const double a = meshcheck::areaOnPlane(A, e.n, e.d);
      const double b = meshcheck::areaOnPlane(B, e.n, e.d);
      if (std::fabs(a - e.above) > 1e-9 || std::fabs(b - (1.0 - e.above)) > 1e-9)
      {
         std::ostringstream os;
         os << "cube face n=(" << e.n.transpose() << ") d=" << e.d << ": A area " << a
            << ", B area " << b << "; expected " << e.above << " / " << 1.0 - e.above;
         throw std::runtime_error(os.str());
      }
   }
}

// ===========================================================================
// BUG-12: mergeMeshes() computed a vertex weld map and never used it
//
// With cap=true the cap's vertices on the cut line coincide with the vertices of
// the clay piece, but they stayed separate vertices, leaving the seam open as far
// as the vertex indices are concerned (a viewer or exporter sees a crack, and
// anything that works on shared vertices breaks).  The cut cube pieces must have
// exactly the vertices they geometrically have.
// ===========================================================================

static void test_bug12_cap_seam_is_welded()
{
   CSG::TriMesh A, B;
   constructFiles("clay_unit_cube.obj", "knife_horiz.obj", true, A, B);

   // Raw indices, not a geometric weld: every output vertex must be unique.
   for (const CSG::TriMesh* m : {&A, &B})
   {
      const Report r = analyze(*m);
      ASSERT(r.closed());
      ASSERT(r.duplicate_vertices == 0);
      ASSERT(r.vertices == 12);  // 4 cube corners + 4 points on the vertical edges + 4 face centres
   }
}

// ===========================================================================
// BUG-13: assembleMesh() mixed two vertex numberings
//
// The weld map is built with the surface being assembled first (so for the knife
// cap the knife vertices come first) but was looked up with canonicalVertexIndex(),
// which always numbers clay vertices first.  In addition, uncut faces used the raw
// mesh index while cut faces used the welded canonical index, so output vertices
// with different origins could be given the same slot in used_map, producing
// wrong cap triangles when the knife has more vertices than the clay.
//
// The clay is a unit cube (8 vertices); the knife is a finely subdivided
// horizontal grid at z = 0.5 (many more vertices), so plenty of uncut knife
// triangles lie inside the cube.  Both halves must be closed with volume 0.5.
// ===========================================================================

static CSG::TriMesh makeGrid(int n, double lo, double hi, double z)
{
   std::vector<double> v;
   std::vector<unsigned int> f;
   for (int j = 0; j <= n; ++j)
      for (int i = 0; i <= n; ++i)
      {
         v.push_back(lo + (hi - lo) * i / n);
         v.push_back(lo + (hi - lo) * j / n);
         v.push_back(z);
      }
   for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
      {
         const unsigned int a = j * (n + 1) + i, b = a + 1, c = a + (n + 1), d = c + 1;
         f.insert(f.end(), {a, b, d,  a, d, c});  // counter-clockwise: normal +z
      }
   return CSG::TriMesh(v, f);
}

static void test_bug13_cap_assembly_indices()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = makeGrid(7, -0.8, 1.9, 0.5);
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   for (const CSG::TriMesh* m : {&A, &B})
   {
      const Report r = analyze(*m);
      ASSERT(r.closed());
      ASSERT(r.duplicate_vertices == 0);
      ASSERT_NEAR(r.volume, 0.5, 1e-9);
   }
}

// ===========================================================================
// BUG-21: AABBTree could not be built or queried when empty
//
// The constructor wrote nodes[-1] for zero objects and query() read
// nodes[NULL_NODE] on an empty tree.  This is reached in normal use: when the
// knife does not touch the clay one of the two output meshes is empty and
// mergeMeshes()/assembleMesh() build trees from zero vertices.
// ===========================================================================

static void test_bug21_empty_aabb_tree()
{
   const CSG::TriMesh empty;
   AABBTree tree = empty.createAABBTree();
   ASSERT(tree.numObjects() == 0);
   ASSERT(tree.query(AABB(Eigen::Vector3d(-1, -1, -1), Eigen::Vector3d(1, 1, 1))).empty());

   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   AABBTree clay_tree = clay.createAABBTree();
   ASSERT(tree.intersect(clay_tree).empty());
   ASSERT(clay_tree.intersect(tree).empty());
   tree.validate();
}

static void test_bug21_knife_far_from_clay_does_not_crash()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeTriangles({{10, 10, 10,   12, 10, 10,   10, 12, 10}});
   for (bool cap : {false, true})
   {
      CSG::TriMesh A, B;
      construct(clay, knife, cap, A, B);
      ASSERT(A.vertices().size() == 0 || analyze(A).faces > 0);
      ASSERT(B.vertices().size() == 0 || analyze(B).faces > 0);
   }
}

// ===========================================================================
// BUG-15: tolerances were absolute / ULP based instead of scale aware
//
// Points were snapped to vertices with a 2 ULP relative test, welded with an
// absolute 1e-6 radius and triangles with an altitude below an absolute 1e-6 were
// thrown away.  Scaling the same scene must give the same result: here a unit cube
// cut by a plane through a vertex at an irrational orientation (so intersection
// points are never exactly on mesh vertices), at scales from 1e-8 to 1e6.
// ===========================================================================

static void test_bug15_scale_invariance()
{
   // Cuts the cube [0,1]^3 into two pieces; the plane passes through the origin
   // vertex with normal (0.3, 0.5, -0.7) / |n|.
   const double n[3] = {0.3, 0.5, -0.7};
   const double nlen = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
   const Eigen::Vector3d normal(n[0] / nlen, n[1] / nlen, n[2] / nlen);
   Eigen::Vector3d u = normal.cross(Eigen::Vector3d::UnitX()).normalized();
   const Eigen::Vector3d w = normal.cross(u);

   // Volume of the part of the cube with n.x > 0 is 1 - volume of the corner cut
   // off by the plane; computed once for scale 1 and scaled by s^3.
   double expected_a = -1.0;

   for (double s : {1e-8, 1e-4, 1.0, 1e3, 1e6})
   {
      const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d::Zero(), Eigen::Vector3d::Constant(s));
      // One big triangle in the plane through the origin containing it, normal n
      const double big = 20.0 * s;
      const Eigen::Vector3d p0 = big * (std::cos(0.3) * u + std::sin(0.3) * w);
      const Eigen::Vector3d p1 = big * (std::cos(2.4) * u + std::sin(2.4) * w);
      const Eigen::Vector3d p2 = big * (std::cos(4.5) * u + std::sin(4.5) * w);
      const Eigen::Vector3d orient = (p1 - p0).cross(p2 - p0);
      const bool flip = orient.dot(normal) < 0;
      const CSG::TriMesh knife = meshcheck::makeTriangles({
         {p0[0], p0[1], p0[2],
          (flip ? p2 : p1)[0], (flip ? p2 : p1)[1], (flip ? p2 : p1)[2],
          (flip ? p1 : p2)[0], (flip ? p1 : p2)[1], (flip ? p1 : p2)[2]}});

      CSG::TriMesh A, B;
      construct(clay, knife, true, A, B);

      const Report ra = analyze(A, 1e-6 * s), rb = analyze(B, 1e-6 * s);
      std::ostringstream where;
      where << "scale " << s << ": ";
      if (!ra.closed() || !rb.closed())
         throw std::runtime_error(where.str() + "a half is not a closed mesh");
      if (ra.faces < 4 || rb.faces < 4)
         throw std::runtime_error(where.str() + "a half is missing");

      const double a = ra.volume / (s * s * s), b = rb.volume / (s * s * s);
      if (std::fabs(a + b - 1.0) > 1e-6)
         throw std::runtime_error(where.str() + "volumes do not add up to the cube's");
      if (expected_a < 0)
         expected_a = a;
      else if (std::fabs(a - expected_a) > 1e-6)
         throw std::runtime_error(where.str() + "volume differs from the unit-scale result");
   }
}

// ===========================================================================
// BUG-17: welding used a box test that reached further than the tolerance and
// never removed the faces that collapsed
//
// Points were welded when their AABBs overlapped: spheres of radius r fattened by
// 5% on each side overlap up to 2.2 r apart per axis, i.e. further than the
// tolerance itself.  The clay is a unit cube whose corner at the origin is
// chamfered by a tiny triangle with legs of 1.5 tol (its three vertices are further
// apart than the tolerance, so they are distinct points), cut by the plane z = 0.5
// far away.  The chamfer must survive in the lower half; its vertices were welded
// together, leaving a face with repeated vertices in the output.
// ===========================================================================

static void test_bug17_weld_does_not_reach_too_far()
{
   // CSGEngine's tolerance is 1e-9 of the diagonal of the bounding box of both meshes.
   const double big = 3.0;
   const double zc = 0.5;
   const std::vector<double> knife_tri = {-big, -big, zc,   2 * big, -big, zc,   -big, 2 * big, zc};
   Eigen::Vector3d lo(-big, -big, 0), hi(2 * big, 2 * big, 1);
   const double tol = 1e-9 * (hi - lo).norm();
   const double L = 1.5 * tol;

   // Unit cube with the corner at the origin cut off by the plane x+y+z = L:
   // vertex 0 = a = (L,0,0), 8 = b = (0,L,0), 9 = c = (0,0,L)
   std::vector<double> v = {
      L, 0, 0,   1, 0, 0,   1, 1, 0,   0, 1, 0,
      0, 0, 1,   1, 0, 1,   1, 1, 1,   0, 1, 1,
      0, L, 0,   0, 0, L};
   std::vector<unsigned int> f = {
      8, 3, 2,  8, 2, 1,  8, 1, 0,     // z = 0 (fan from b over the pentagon b,3,2,1,a)
      0, 1, 5,  0, 5, 4,  0, 4, 9,     // y = 0 (a,1,5,4,c)
      9, 4, 7,  9, 7, 3,  9, 3, 8,     // x = 0 (c,4,7,3,b)
      4, 5, 6,  4, 6, 7,               // z = 1
      3, 7, 6,  3, 6, 2,               // y = 1
      1, 2, 6,  1, 6, 5,               // x = 1
      0, 9, 8};                        // the chamfer
   const CSG::TriMesh clay(v, f);
   const CSG::TriMesh knife = meshcheck::makeTriangles({knife_tri});

   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   // With a weld tolerance below the chamfer size nothing may collapse
   for (const CSG::TriMesh* m : {&A, &B})
   {
      const Report r = analyze(*m, 0.5 * tol);
      ASSERT(r.degenerate_faces == 0);
      ASSERT(r.closed());
   }
   // B (below the knife) is the lower half with the chamfered corner
   ASSERT_NEAR(analyze(B, 0.5 * tol).volume, 0.5 - L * L * L / 6.0, 1e-9);
   ASSERT_NEAR(analyze(A, 0.5 * tol).volume, 0.5, 1e-9);
}

// ===========================================================================
// BUG-18: faces the flood fill cannot reach were silently dropped from both outputs
//
// Uncut faces get their side by flood fill from classified cut faces.  A part of
// the clay that is not connected to anything the knife touches (another shell of
// the clay, or the whole clay when the knife misses it) was never reached and
// ended up in neither output.  Such a part must be classified by where it is
// relative to the knife.
// ===========================================================================

static CSG::TriMesh makeBoxes(const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& boxes)
{
   std::vector<double> v;
   std::vector<unsigned int> f;
   for (const auto& b : boxes)
   {
      const CSG::TriMesh box = meshcheck::makeBox(b.first, b.second);
      const unsigned int base = static_cast<unsigned int>(v.size() / 3);
      for (const auto& p : box.vertices())
         v.insert(v.end(), {p[0], p[1], p[2]});
      for (const auto& tri : box.faces())
         f.insert(f.end(), {base + tri.m_v[0], base + tri.m_v[1], base + tri.m_v[2]});
   }
   return CSG::TriMesh(v, f);
}

static void test_bug18_unreached_shells_are_classified()
{
   // Box 1 straddles the knife (z = 0.5); box 2 is entirely above it, box 3 below it
   const CSG::TriMesh clay = makeBoxes({
      {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1)},
      {Eigen::Vector3d(0, 0, 2), Eigen::Vector3d(1, 1, 3)},
      {Eigen::Vector3d(0, 0, -3), Eigen::Vector3d(1, 1, -2)}});
   const CSG::TriMesh knife = meshcheck::makeTriangles({{-20, -20, 0.5,   40, -20, 0.5,   -20, 40, 0.5}});

   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   const Report ra = analyze(A), rb = analyze(B);
   ASSERT(ra.closed());
   ASSERT(rb.closed());
   ASSERT_NEAR(ra.volume, 1.5, 1e-9);  // upper half of box 1 + box 2
   ASSERT_NEAR(rb.volume, 1.5, 1e-9);  // lower half of box 1 + box 3
}

static void test_bug18_knife_misses_clay()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));

   // Knife below the cube, normal +z: the whole cube is above it
   {
      const CSG::TriMesh knife = meshcheck::makeTriangles({{-20, -20, -5,   40, -20, -5,   -20, 40, -5}});
      CSG::TriMesh A, B;
      construct(clay, knife, true, A, B);
      ASSERT_NEAR(analyze(A).volume, 1.0, 1e-9);
      ASSERT(analyze(A).closed());
      ASSERT(B.faces().size() == 0);
   }
   // Knife above the cube: the whole cube is below it
   {
      const CSG::TriMesh knife = meshcheck::makeTriangles({{-20, -20, 5,   40, -20, 5,   -20, 40, 5}});
      CSG::TriMesh A, B;
      construct(clay, knife, true, A, B);
      ASSERT(A.faces().size() == 0);
      ASSERT_NEAR(analyze(B).volume, 1.0, 1e-9);
      ASSERT(analyze(B).closed());
   }
}

// ===========================================================================
// BUG-19: nearly coplanar triangles took the generic intersection path
//
// Coplanarity was only detected when all three signed distances were exactly 0.
// A knife that is coplanar with a cube face to within rounding (here a grid in the
// plane z = 1, the top face, with 1e-16 of noise in z) was handled by the generic
// code with a denominator of 0, which produced NaN points; Triangle then never
// returned.  The knife is within tolerance of the top face: no intersection
// geometry can be constructed, but the library must terminate with finite output,
// and the cube, which lies entirely below the knife, ends up in B.
// ===========================================================================

static void test_bug19_nearly_coplanar_knife()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = loadMesh("bugs/bug19_knife.obj");
   CSG::TriMesh A, B;
   construct(clay, knife, false, A, B);

   for (const CSG::TriMesh* m : {&A, &B})
      for (const auto& v : m->vertices())
         ASSERT(std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]));

   ASSERT(A.faces().size() == 0);
   const Report rb = analyze(B);
   ASSERT_NEAR(rb.area, 6.0, 1e-6);  // the whole cube surface, nothing lost
}

// ===========================================================================
// BUG-27: faces coplanar with the knife were not cut
//
// A pair of triangles that is coplanar was skipped ("TODO: will handle this
// later"), so a clay face that lies in the knife's plane was not subdivided where
// the knife's edges cross it, while its neighbours (which the knife does touch) were
// subdivided along the shared edge.  The result had the right surface but was not
// watertight: T-junctions along those edges.
//
// Same scene as BUG-19: the cube's top face is coplanar with the knife grid; the
// cube is below the knife, so B should be the closed cube.
// ===========================================================================

static void test_bug27_coplanar_knife_is_watertight()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = loadMesh("bugs/bug19_knife.obj");
   CSG::TriMesh A, B;
   construct(clay, knife, false, A, B);
   ASSERT(analyze(B).closed());
}

// Coplanar faces of two solids: a box that shares the top face plane of the unit cube
// (same orientation: both outward normals are +z) is cut out of the cube; A (clay minus
// knife) and B (clay intersected with the knife) must both be closed solids.
static void test_bug27_coplanar_same_orientation()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeBox(Eigen::Vector3d(0.25, 0.25, 0.5), Eigen::Vector3d(0.75, 0.75, 1));
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   const Report ra = analyze(A), rb = analyze(B);
   ASSERT(ra.closed());
   ASSERT(rb.closed());
   ASSERT_NEAR(ra.volume, 1.0 - 0.125, 1e-9);
   ASSERT_NEAR(rb.volume, 0.125, 1e-9);
}

// A box that only touches the top face of the cube (opposite orientations: the box's
// bottom faces down, the cube's top faces up) removes nothing from the cube and the
// intersection is empty.
static void test_bug27_coplanar_opposite_orientation()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeBox(Eigen::Vector3d(0.25, 0.25, 1), Eigen::Vector3d(0.75, 0.75, 1.5));
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   const Report ra = analyze(A);
   ASSERT(ra.closed());
   ASSERT_NEAR(ra.volume, 1.0, 1e-9);
   ASSERT_NEAR(analyze(B).volume, 0.0, 1e-9);
}

// ===========================================================================
// BUG-23: CSGEngine::construct() kept state from previous calls
//
// The points created by the intersections and the lists of cutting faces were
// members that were never cleared, so every further call added to them.  No wrong
// output was found (the stale entries were simply unused), but the state grew
// with every call; calling construct() repeatedly, with and without caps, must
// give identical results each time.
// ===========================================================================

static void test_bug23_construct_is_repeatable()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = loadMesh("bugs/bug20_knife.obj");
   CSG::CSGEngine engine(clay, knife);

   Report first[2];
   for (int round = 0; round < 4; ++round)
   {
      const bool cap = (round % 2) == 1;
      CSG::TriMesh A, B;
      engine.construct(CSG::kDifference, cap, A, B);
      const Report ra = analyze(A), rb = analyze(B);
      if (round < 2)
      {
         first[round] = ra;  // remember what the first call with / without cap gave
         continue;
      }
      const Report& want = first[round - 2];
      ASSERT(ra.faces == want.faces && ra.vertices == want.vertices);
      ASSERT_NEAR(ra.area, want.area, 1e-12);
      ASSERT_NEAR(ra.volume, want.volume, 1e-12);
      ASSERT(rb.faces > 0);
   }
}

// ===========================================================================
// BUG-28: CSGEngine kept references to the caller's meshes
//
// The constructor took const references and stored them, so an engine built from
// meshes that are destroyed (or temporaries, which is how Python tends to use it)
// before construct() was called read freed memory.  The engine must keep its own
// copy.
// ===========================================================================

static std::unique_ptr<CSG::CSGEngine> makeEngineFromTemporaries()
{
   CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   CSG::TriMesh knife = meshcheck::makeTriangles({{-5, -5, 0.5,   9, -5, 0.5,   -5, 9, 0.5}});
   return std::unique_ptr<CSG::CSGEngine>(new CSG::CSGEngine(clay, knife));
   // clay and knife are destroyed here
}

static void test_bug28_engine_owns_its_meshes()
{
   std::unique_ptr<CSG::CSGEngine> engine = makeEngineFromTemporaries();

   // Reuse (and scribble over) the memory the meshes occupied
   std::vector<std::vector<double>> churn;
   for (int ii = 0; ii < 200; ++ii)
      churn.push_back(std::vector<double>(64, std::nan("")));

   CSG::TriMesh A, B;
   engine->construct(CSG::kDifference, true, A, B);
   ASSERT_NEAR(analyze(A).volume, 0.5, 1e-9);
   ASSERT_NEAR(analyze(B).volume, 0.5, 1e-9);
}

// ===========================================================================
// BUG-29: TriMesh::faceNormal() was not a unit vector
//
// It took the cross product of two *normalized* edges, whose length is the sine of
// the angle between them.  It must return the unit normal (and zero, not NaN, for a
// degenerate triangle).
// ===========================================================================

static void test_bug29_face_normal_is_unit()
{
   const CSG::TriMesh m = meshcheck::makeTriangles({
      {0, 0, 0,   1, 0, 0,   0.2, 1, 0.7},          // skewed
      {0, 0, 0,   10, 0, 0,  10, 0.1, 0.05},        // thin
      {0, 0, 0,   0, 1, 0,   1, 1, 0},              // flat, facing down
      {0, 0, 0,   1, 1, 1,   2, 2, 2}});            // collinear: no normal
   for (uint32_t f = 0; f < 3; ++f)
   {
      const Eigen::Vector3d n = m.faceNormal(f);
      const auto& v = m.vertices();
      const auto& t = m.faces()[f];
      const Eigen::Vector3d expected = (v[t.m_v[1]] - v[t.m_v[0]]).cross(v[t.m_v[2]] - v[t.m_v[0]]).normalized();
      ASSERT_NEAR(n.norm(), 1.0, 1e-12);
      ASSERT_NEAR((n - expected).norm(), 0.0, 1e-12);
   }
   const Eigen::Vector3d degenerate = m.faceNormal(3);
   ASSERT(std::isfinite(degenerate[0]) && std::isfinite(degenerate[1]) && std::isfinite(degenerate[2]));
   ASSERT_NEAR(degenerate.norm(), 0.0, 0.0);
}

// ===========================================================================
// BUG-30: the OBJ writer lost precision
//
// operator<<(TriMesh) wrote coordinates with the stream default of 6 significant
// digits, so a mesh written out and read back was a different mesh (and the output
// of a cut at a small scale was garbage).  Coordinates must round-trip exactly.
// ===========================================================================

static void test_bug30_obj_round_trip()
{
   const std::vector<double> v = {
      0.1 + 0.2, 1.0 / 3.0, 12345678.123456789,
      1.23456789012345e-8, -2.0 / 3.0e5, 1e-300,
      1.0, 0.0, -0.0};
   const std::vector<unsigned int> f = {0, 1, 2};
   const CSG::TriMesh mesh(v, f);

   const std::string path = "csg_bug30_round_trip.obj";
   {
      std::ofstream out(path);
      out << mesh;
   }
   CSG::TriMesh back;
   ASSERT(back.loadOBJ(path) == 3);
   std::remove(path.c_str());

   ASSERT(back.faces().size() == 1);
   for (int ii = 0; ii < 3; ++ii)
      for (int jj = 0; jj < 3; ++jj)
         ASSERT(back.vertices()[ii][jj] == mesh.vertices()[ii][jj]);
}

// ===========================================================================
// BUG-31: the OBJ loader mishandled valid and invalid input
//
//  * "v x y z r g b" (vertex colours) and "v x y z w" lines were ignored, shifting
//    every later vertex index;
//  * polygons with more than three corners were silently dropped;
//  * negative (relative) indices were not resolved and produced negative vertex
//    indices;
//  * faces that refer to a vertex that does not exist were accepted, so later code
//    indexed out of range.
// ===========================================================================

static void test_bug31_obj_loader()
{
   const std::string path = "csg_bug31_probe.obj";
   {
      std::ofstream out(path);
      out << "# vertex colours, a quad, relative indices and a bad face\n"
          << "v 0 0 0 1 0 0\n"
          << "v 1 0 0 0 1 0\n"
          << "v 1 1 0 0 0 1\n"
          << "v 0 1 0 1 1 1\n"
          << "vt 0 0\nvt 1 0\nvt 1 1\nvt 0 1\n"
          << "f 1/1 2/2 3/3 4/4\n"       // a quad: two triangles
          << "f -4 -3 -2\n"              // relative: vertices 1 2 3
          << "f 1 2 99\n"                // no such vertex: skipped
          << "f 0 1 2\n";                // index 0 is not valid: skipped
   }
   CSG::TriMesh m;
   const uint32_t num_vertices = m.loadOBJ(path);
   std::remove(path.c_str());

   ASSERT(num_vertices == 4);
   ASSERT(m.faces().size() == 3);
   // quad fan: (1 2 3) (1 3 4); relative: (1 2 3)
   const int32_t expected[3][3] = {{0, 1, 2}, {0, 2, 3}, {0, 1, 2}};
   for (int f = 0; f < 3; ++f)
      for (int c = 0; c < 3; ++c)
         ASSERT(m.faces()[f].m_v[c] == expected[f][c]);
   for (const auto& face : m.faces())
      for (int c = 0; c < 3; ++c)
         ASSERT(face.m_v[c] >= 0 && face.m_v[c] < static_cast<int32_t>(m.vertices().size()));
}

// ===========================================================================
// BUG-32: invalid meshes were accepted
//
// TriMesh(vertices, faces) took lists whose length is not a multiple of 3 (silently
// truncating them) and faces that refer to vertices that do not exist, and
// CSGEngine accepted such meshes, which are then indexed out of range.  Both must
// refuse with std::invalid_argument.
// ===========================================================================

template <class Exception>
static bool throwsWhen(std::function<void()> fn)
{
   try { fn(); } catch (const Exception&) { return true; }
   return false;
}

static void test_bug32_invalid_meshes_are_rejected()
{
   const std::vector<double> v = {0, 0, 0,  1, 0, 0,  0, 1, 0};

   // Valid input still works
   ASSERT(!throwsWhen<std::invalid_argument>([&] { CSG::TriMesh m(v, std::vector<unsigned int>{0, 1, 2}); }));
   // Face index out of range
   ASSERT(throwsWhen<std::invalid_argument>([&] { CSG::TriMesh m(v, std::vector<unsigned int>{0, 1, 3}); }));
   // Lists that are not whole triples
   ASSERT(throwsWhen<std::invalid_argument>([&] { CSG::TriMesh m(v, std::vector<unsigned int>{0, 1}); }));
   ASSERT(throwsWhen<std::invalid_argument>([&] {
      CSG::TriMesh m(std::vector<double>{0, 0, 0, 1, 0}, std::vector<unsigned int>{0, 0, 0}); }));

   // A mesh edited into an invalid state through the accessors is refused by the engine
   CSG::TriMesh bad(v, std::vector<unsigned int>{0, 1, 2});
   bad.faces()[0].m_v[2] = 7;
   const CSG::TriMesh good = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   ASSERT(throwsWhen<std::invalid_argument>([&] { CSG::CSGEngine e(good, bad); }));
   ASSERT(throwsWhen<std::invalid_argument>([&] { CSG::CSGEngine e(bad, good); }));
   bad.faces()[0].m_v[2] = -1;
   ASSERT(throwsWhen<std::invalid_argument>([&] { CSG::CSGEngine e(good, bad); }));
}

// ===========================================================================
// BUG-33: the hash of a pair of face indices collided heavily
//
// std::hash<std::pair<uint32_t,uint32_t>> was h1 ^ (h2 << 1).  For the pairs of
// indices that the broad phase produces (small numbers) that has only a few
// hundred distinct values, so the unordered_set of candidate pairs degenerated into
// long chains.  A 100 x 100 grid of pairs must hash to (almost) 10000 different
// values.
// ===========================================================================

static void test_bug33_pair_hash_distribution()
{
   std::hash<std::pair<uint32_t, uint32_t>> hasher;
   std::set<size_t> values;
   for (uint32_t a = 0; a < 100; ++a)
      for (uint32_t b = 0; b < 100; ++b)
         values.insert(hasher(std::make_pair(a, b)));
   ASSERT(values.size() >= 9990);
   ASSERT(hasher(std::make_pair(1u, 2u)) != hasher(std::make_pair(2u, 1u)));
}

// ===========================================================================
// BUG-35: the 'operation' argument of construct() was silently ignored
//
// kUnion, kIntersection and even kBadValue all did a difference.  Only kDifference
// is implemented, so anything else must be refused instead of returning the wrong
// meshes.
// ===========================================================================

static void test_bug35_unimplemented_operations_are_refused()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeTriangles({{-5, -5, 0.5,   9, -5, 0.5,   -5, 9, 0.5}});
   CSG::CSGEngine engine(clay, knife);
   CSG::TriMesh A, B;

   for (CSG::CSGOperation op : {CSG::kUnion, CSG::kIntersection, CSG::kBadValue})
      ASSERT(throwsWhen<std::invalid_argument>([&] { engine.construct(op, false, A, B); }));
   ASSERT(!throwsWhen<std::invalid_argument>([&] { engine.construct(CSG::kDifference, false, A, B); }));
}

// ===========================================================================
// BUG-36: T-junctions where a cut ends on an edge shared with a neighbouring face
//
// The knife box [0.75,1.25] x [0,2] x [-0.75,0.75] cuts the unit cube so that cut lines
// of the cube's y = 1 face meet exactly at (0.75, 1, 0.75), a point on the diagonal
// that separates the face's two triangles.  The triangle on one side of the diagonal
// is cut and gets that point; the other one only touches the knife in that point, and
// the exact triangle-triangle test reported no intersection for it, so it was not
// subdivided and the two meshes met in a T-junction (an open edge).  Neighbouring faces
// must agree on the points on their shared edge.
// ===========================================================================

static void test_bug36_t_junction_on_shared_edge()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeBox(Eigen::Vector3d(0.75, 0, -0.75), Eigen::Vector3d(1.25, 2, 0.75));
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   const Report ra = analyze(A), rb = analyze(B);
   ASSERT(ra.closed());
   ASSERT(rb.closed());
   // clay minus knife and clay intersected with knife: the knife covers x in [0.75,1],
   // z in [0,0.75] of the cube (all y)
   ASSERT_NEAR(rb.volume, 0.25 * 0.75, 1e-9);
   ASSERT_NEAR(ra.volume, 1.0 - 0.25 * 0.75, 1e-9);
}

// ===========================================================================
// BUG-37: faces without area in the input produced NaN points and hung Triangle
//
// A clay face with a repeated vertex (or three collinear ones) has no plane: the
// intersection construction and closestPointOnTriangle() divided by zero, the NaN
// points reached the retriangulation of the face (neighbours of a cut face are given
// its edge points too, BUG-36) and Triangle never returned.  The degenerate faces
// are just left alone: construct() must finish, give finite coordinates, and cut
// the rest of the mesh as usual.
// ===========================================================================

static void test_bug37_degenerate_input_faces()
{
   const std::vector<double> cube_vertices = {
      0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0,  0, 0, 1,  1, 0, 1,  1, 1, 1,  0, 1, 1,
      0, 0, 0.25};  // 8: a point on the edge 0-4, which the knife crosses
   const std::vector<unsigned int> cube_faces = {
      0, 2, 1,  0, 3, 2,  4, 5, 6,  4, 6, 7,  0, 1, 5,  0, 5, 4,
      3, 7, 6,  3, 6, 2,  0, 4, 7,  0, 7, 3,  1, 2, 6,  1, 6, 5};
   const CSG::TriMesh plane = meshcheck::makeTriangles({{-3, -3, 0.5,   9, -3, 0.5,   -3, 9, 0.5}});

   // Extra faces, all on edges that the knife crosses (so that they get cut points from
   // their neighbours): a repeated vertex (on the edge 1-5), three collinear vertices
   // (0, 8 and 4 on the edge 0-4), and a face that is a point
   const std::vector<std::vector<unsigned int>> extras = {{1, 1, 5}, {0, 8, 4}, {5, 5, 5}};
   for (const auto& extra : extras)
   {
      std::vector<unsigned int> faces = cube_faces;
      faces.insert(faces.end(), extra.begin(), extra.end());
      const CSG::TriMesh clay(cube_vertices, faces);

      for (bool cap : {false, true})
      {
         CSG::TriMesh A, B;
         construct(clay, plane, cap, A, B);
         for (const CSG::TriMesh* m : {&A, &B})
            for (const auto& v : m->vertices())
               ASSERT(std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]));
         if (cap)
         {
            // The degenerate face has no volume: the halves still add up to the cube
            ASSERT_NEAR(analyze(A).volume, 0.5, 1e-9);
            ASSERT_NEAR(analyze(B).volume, 0.5, 1e-9);
         }
      }
   }

   // A degenerate knife triangle next to the real one
   {
      const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
      const CSG::TriMesh knife = meshcheck::makeTriangles({
         {-3, -3, 0.5,   9, -3, 0.5,   -3, 9, 0.5},
         {0.2, 0.2, 0.5,   0.2, 0.2, 0.5,   0.7, 0.3, 0.5}});
      CSG::TriMesh A, B;
      construct(clay, knife, true, A, B);
      ASSERT_NEAR(analyze(A).volume, 0.5, 1e-9);
      ASSERT_NEAR(analyze(B).volume, 0.5, 1e-9);
   }
}

// ===========================================================================
// BUG-38: faces were classified against the triangles the intersection test reported
//
// A sub-face of a cut face was classified by its position relative to the nearest of
// the knife triangles that the (exact) intersection test had reported for the face.
// Where two surfaces just touch, that list is not reliable: of the four knife triangles
// around a pole vertex that coincides with a clay vertex, only two were reported for one
// of the clay faces around it, the nearest triangle among those was a different feature
// than the real nearest one, and the face landed on the wrong side of the knife (an open
// edge).  The nearest triangle must be found among all triangles of the other surface.
//
// Two UV spheres of radius 1 with the same centre and different tessellations (their
// poles coincide, the surfaces touch and cross everywhere): both halves must be closed
// and their volumes add up to the clay's.
// ===========================================================================

static void check_sphere_pair(int clay_segments, int clay_rings, int knife_segments, int knife_rings)
{
   const CSG::TriMesh clay = meshcheck::makeUVSphere(Eigen::Vector3d::Zero(), 1.0, clay_segments, clay_rings);
   const CSG::TriMesh knife = meshcheck::makeUVSphere(Eigen::Vector3d::Zero(), 1.0, knife_segments, knife_rings);
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   const Report ra = analyze(A), rb = analyze(B), rc = analyze(clay);
   std::ostringstream where;
   where << "spheres (" << clay_segments << "," << clay_rings << ") / ("
         << knife_segments << "," << knife_rings << "): ";
   if (!ra.closed() || !rb.closed())
      throw std::runtime_error(where.str() + "a half is not a closed mesh");
   if (std::fabs(ra.volume + rb.volume - rc.volume) > 1e-9)
      throw std::runtime_error(where.str() + "volumes do not add up");
}

static void test_bug38_touching_spheres()
{
   check_sphere_pair(6, 4, 4, 3);
   check_sphere_pair(4, 4, 6, 3);
   check_sphere_pair(6, 6, 4, 4);
   check_sphere_pair(24, 12, 16, 8);
}

// ===========================================================================
// BUG-39: the library wrote progress messages to stdout unconditionally
//
// "Building AABB trees..", "N AABB collisions found", "Final assembled mesh: ..." and
// "N potential intersections found" were printed by every construct() call, which
// pollutes the stdout of any program (or the Blender / Python console) that uses
// the library.  Nothing may be printed unless setVerbose(true) is called.  (This
// library has no setVerbose(); only the quiet half is tested.)
// ===========================================================================

static std::string capturedStdout(std::function<void()> fn)
{
   std::ostringstream captured;
   std::streambuf* old = std::cout.rdbuf(captured.rdbuf());
   try { fn(); } catch (...) { std::cout.rdbuf(old); throw; }
   std::cout.rdbuf(old);
   return captured.str();
}

static void test_bug39_quiet_by_default()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = meshcheck::makeTriangles({{-5, -5, 0.5,   9, -5, 0.5,   -5, 9, 0.5}});

   const std::string quiet = capturedStdout([&] {
      CSG::TriMesh A, B;
      construct(clay, knife, true, A, B);
   });
   ASSERT(quiet.empty());

}

// ===========================================================================
// BUG-40: NaN and infinite vertex coordinates were accepted
//
// They turn into NaN distances everywhere (the scene tolerance, the AABBs, the
// intersection points) and garbage or a hang later.  The engine must refuse them.
// ===========================================================================

static void test_bug40_non_finite_coordinates_are_rejected()
{
   const CSG::TriMesh good = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   for (double bad_value : {std::nan(""), std::numeric_limits<double>::infinity()})
   {
      CSG::TriMesh bad = good;
      bad.vertices()[3][1] = bad_value;
      ASSERT(throwsWhen<std::invalid_argument>([&] { CSG::CSGEngine e(bad, good); }));
      ASSERT(throwsWhen<std::invalid_argument>([&] { CSG::CSGEngine e(good, bad); }));
   }
   ASSERT(!throwsWhen<std::invalid_argument>([&] { CSG::CSGEngine e(good, good); }));
}

// ===========================================================================
// BUG-42: retriangulate() compared every pair of cut points of a face
//
// A clay triangle crossed by a fine knife has a cut point per crossing; welding them
// by comparing every pair is quadratic in their number.  The accordion below (a knife
// folded N times across one big clay triangle) must take about 4 times as long for
// 4 times the panels, not 16 times.  (The quadratic loop is only visible in an
// unoptimized build, where it made 2000 panels take 15 s; with optimization Triangle's
// own time dominates, so this test passes either way there.)
// ===========================================================================

static double accordionSeconds(int panels)
{
   const CSG::TriMesh clay = meshcheck::makeTriangles({{-1, -1, 0,   3, -1, 0,   -1, 3, 0}});
   std::vector<double> v;
   std::vector<unsigned int> f;
   for (int i = 0; i <= panels; ++i)
   {
      const double x = static_cast<double>(i) / panels;
      const double z = (i % 2 == 0) ? -1.0 : 1.0;
      v.insert(v.end(), {x, -2.0, z,  x, 4.0, z});
   }
   for (int i = 0; i < panels; ++i)
   {
      const unsigned int a = 2 * i, b = 2 * i + 1, c = 2 * i + 2, d = 2 * i + 3;
      f.insert(f.end(), {a, b, d,  a, d, c});
   }
   const CSG::TriMesh knife(v, f);

   const auto start = std::chrono::steady_clock::now();
   CSG::TriMesh A, B;
   construct(clay, knife, false, A, B);
   const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
   ASSERT(A.faces().size() == static_cast<size_t>(2 * panels));  // every panel cuts the triangle in two
   return seconds;
}

static void test_bug42_retriangulation_scales_linearly()
{
   accordionSeconds(50);  // warm up
   const double small = accordionSeconds(200);
   const double large = accordionSeconds(800);
   std::ostringstream os;
   os << "200 panels took " << small << " s, 800 panels " << large << " s";
   if (large > 10.0 * std::max(small, 1e-3))
      throw std::runtime_error(os.str() + ": much more than 4x");
}

// ===========================================================================
// BUG-43: Triangle called exit() when it found an internal inconsistency
//
// triexit() in triangle.c ends the program.  For a cube and a box whose faces and
// edges coincide to within rounding at coordinates around 1e8 (so the cut segments of
// a face are numerically inconsistent) Triangle reports "Topological inconsistency
// after splitting a segment" and exited: the program that uses the library (Blender, a
// Python session) went away with it.  The failure must be survivable: the face is left
// uncut and construct() returns.
// ===========================================================================

static void test_bug43_triangle_failure_does_not_exit()
{
   const CSG::TriMesh clay = loadMesh("bugs/bug43_clay.obj");
   const CSG::TriMesh knife = loadMesh("bugs/bug43_knife.obj");
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);   // used to exit(1)

   for (const CSG::TriMesh* m : {&A, &B})
      for (const auto& v : m->vertices())
         ASSERT(std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]));
   // (Here the box contains the cube, so A may well be empty; what must hold is that the two
   // halves still make up the whole cube.)
   ASSERT_NEAR(analyze(A, 1e-6).volume + analyze(B, 1e-6).volume, analyze(clay, 1e-6).volume, 1e-6);
}

// ===========================================================================
// BUG-44: the tolerance could be smaller than the resolution of the coordinates
//
// The tolerance is 1e-9 of the size of the scene.  Far from the origin neighbouring
// doubles are further apart than that (1.5e-8 at 6.7e7), so a point computed as the
// intersection of two edges is rounded to a grid coarser than the tolerance and no
// longer matches the mesh vertex it coincides with exactly.  The tolerance must not be
// below a few units in the last place of the largest coordinate.
//
// The scene: a unit cube cut by a box, rotated and moved to about (2^26, 2^26, 2^26)
// (data/bugs/bug44_*.obj, found by tools/fuzz_solid.py --nice --rotate --offset 67108864).
// ===========================================================================

static void test_bug44_tolerance_respects_coordinate_resolution()
{
   const CSG::TriMesh clay = loadMesh("bugs/bug44_clay.obj");
   const CSG::TriMesh knife = loadMesh("bugs/bug44_knife.obj");
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   const double tol = 1e-6;
   const Report ra = analyze(A, tol), rb = analyze(B, tol), rc = analyze(clay, tol);
   ASSERT(ra.closed());
   ASSERT(rb.closed());
   ASSERT_NEAR(rc.volume, 1.0, 1e-6);
   ASSERT_NEAR(ra.volume + rb.volume, 1.0, 1e-6);
}

// ===========================================================================
// BUG-45: Triangle's internal errors could not be recovered from, and were printed
//
// A random soup of triangles on a coarse grid (coincident vertices, collinear and
// overlapping segments, tools/fuzz_soup.py --nice --seed 5 --only 245) makes Triangle
// fail while inserting crossing segments: "Internal error in segmentintersection(): ...
// Please report this bug", which it printed to stdout before ending the program, and
// with the exception hook of BUG-43 it left the face uncut.  The segments are now split at
// the points where they cross, touch or overlap before Triangle sees them, so the
// face is retriangulated normally, and nothing is printed on stdout or stderr.
// ===========================================================================

static std::string capturedStderr(std::function<void()> fn)
{
   std::ostringstream captured;
   std::streambuf* old = std::cerr.rdbuf(captured.rdbuf());
   try { fn(); } catch (...) { std::cerr.rdbuf(old); throw; }
   std::cerr.rdbuf(old);
   return captured.str();
}

static void test_bug45_triangle_failure_is_recovered()
{
   const CSG::TriMesh clay = loadMesh("bugs/bug45_clay.obj");
   const CSG::TriMesh knife = loadMesh("bugs/bug45_knife.obj");
   CSG::TriMesh A, B;

   std::string err;
   const std::string out = capturedStdout([&] {
      err = capturedStderr([&] { construct(clay, knife, true, A, B); });
   });
   ASSERT(out.empty());   // Triangle's own messages are not passed on
   ASSERT(err.empty());   // and the face was not given up on ("could not retriangulate")
   for (const CSG::TriMesh* m : {&A, &B})
      for (const auto& v : m->vertices())
         ASSERT(std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]));
}

// ===========================================================================
// BUG-5: clay_tetrahedron.obj had inverted face winding
// BUG-6: the output arrays of Triangle were never freed
// BUG-8: AABBTree::intersect() assumed that the keys are 0, 1, 2, ...
//
// (Fixed before the second review of the code; the tests are here so that every bug in
// TODO.md has one.)
// ===========================================================================

static void test_bug05_tetrahedron_faces_point_outwards()
{
   const CSG::TriMesh tet = loadMesh("clay_tetrahedron.obj");
   const Report r = analyze(tet);
   ASSERT(r.closed());
   ASSERT(r.volume > 0.0);                       // outward facing: positive signed volume
   ASSERT_NEAR(r.volume, 8.0 / 3.0, 1e-9);       // the regular tetrahedron with edges 2*sqrt(2)
}

static long residentKilobytes()
{
   std::ifstream statm("/proc/self/statm");
   long pages = 0, resident = -1;
   if (statm >> pages >> resident)
      return resident * 4;   // pages of 4 KiB
   return -1;
}

static void test_bug06_construct_does_not_leak()
{
   const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   const CSG::TriMesh knife = loadMesh("bugs/bug20_knife.obj");
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);   // warm up

#if defined(__SANITIZE_ADDRESS__)
   return;  // AddressSanitizer keeps freed memory in quarantine, so resident memory says nothing
#endif
   const long before = residentKilobytes();
   if (before < 0)
      return;  // no /proc: cannot measure
   for (int ii = 0; ii < 600; ++ii)
      construct(clay, knife, true, A, B);
   const long grown = residentKilobytes() - before;
   if (grown > 512)   // each leaked call used to cost a few kilobytes
   {
      std::ostringstream os;
      os << "resident memory grew by " << grown << " KiB over 600 calls";
      throw std::runtime_error(os.str());
   }
}

static void test_bug08_aabb_intersect_with_arbitrary_keys()
{
   // Two trees whose keys are not 0..n-1
   AABBTree a(0.0, 4), b(0.0, 4);
   a.insertAABB(10, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
   a.insertAABB(77, Eigen::Vector3d(5, 5, 5), Eigen::Vector3d(6, 6, 6));
   b.insertAABB(1000, Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Vector3d(2, 2, 2));
   b.insertAABB(2000, Eigen::Vector3d(9, 9, 9), Eigen::Vector3d(10, 10, 10));

   const auto hits = a.intersect(b);   // used to throw std::invalid_argument for key 0
   ASSERT(hits.size() == 1);
   ASSERT(hits.count(std::make_pair(10u, 1000u)) == 1);
}

// ===========================================================================
// BUG-46: a cut point on an edge could be dropped from one of the two faces sharing it
//
// A point on an edge of a face is on it only up to rounding: in the face's 2D frame it can be
// a hair outside the triangle.  Triangle's convex hull then bulges out to include it; the
// sliver between the hull and the edge lies outside of the face, and dropping it (all slivers
// are dropped) removes the point from the boundary of that face, while the neighbouring face,
// where the point came out a hair inside, keeps it: a T-junction, here an open triangle in
// the intersection of a cube and the box that contains it (data/bugs/bug43_*.obj, where a
// cube edge and a box edge coincide).
// ===========================================================================

static void test_bug46_boundary_points_stay_inside_the_face()
{
   const CSG::TriMesh clay = loadMesh("bugs/bug43_clay.obj");
   const CSG::TriMesh knife = loadMesh("bugs/bug43_knife.obj");
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   // The box contains the cube: B is the whole cube, closed
   const Report rb = analyze(B, 1e-6);
   ASSERT(rb.closed());
   ASSERT_NEAR(rb.volume, analyze(clay, 1e-6).volume, 1e-6);
   ASSERT(A.faces().size() == 0);
}

// ===========================================================================
// BUG-47: neighbouring faces could still disagree about the vertices on their shared edge
//
// Several cut points on one edge of a face (here a cube edge that a box edge crosses at four
// points, at coordinates around 1e4 so that nothing is exactly collinear) are triangulated
// by Triangle as a thin strip along the edge; which of the points stay on the face's boundary
// once the slivers of the strip are dropped is arbitrary and can differ between the two faces
// that share the edge: a T-junction.  The output is therefore conformed: every triangle edge is
// split at the vertices that lie on it (data/bugs/bug47_*.obj, found by
// tools/fuzz_solid.py --nice --rotate --offset 1e4).
// ===========================================================================

static void test_bug47_output_has_no_t_junctions()
{
   const CSG::TriMesh clay = loadMesh("bugs/bug47_clay.obj");
   const CSG::TriMesh knife = loadMesh("bugs/bug47_knife.obj");
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);

   const Report ra = analyze(A), rb = analyze(B), rc = analyze(clay);
   ASSERT(ra.closed());
   ASSERT(rb.closed());
   ASSERT_NEAR(ra.volume + rb.volume, rc.volume, 1e-6);
}

// ===========================================================================
// BUG-48: the T-junction conforming pass never ended on nearly degenerate triangles
//
// conformTJunctions() (BUG-47) split triangles at vertices found within the tolerance of
// their edges and looked at the pieces' edges again with the same test: for a random
// triangle soup on a coarse grid (data/bugs/bug48_*.obj, tools/fuzz_soup.py --nice
// --seed 71 --tris 25) a piece again had a vertex within the tolerance of one of its edges,
// over and over.  The vertices on the edges are now found once, in the original triangle,
// and everything after that is combinatorial.  (A hang makes this test fail: ctest gives
// csg_bug_tests a time limit.)
// ===========================================================================

static void test_bug48_conforming_terminates()
{
   const CSG::TriMesh clay = loadMesh("bugs/bug48_clay.obj");
   const CSG::TriMesh knife = loadMesh("bugs/bug48_knife.obj");
   CSG::TriMesh A, B;
   construct(clay, knife, true, A, B);
   for (const CSG::TriMesh* m : {&A, &B})
      for (const auto& v : m->vertices())
         ASSERT(std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]));
}

// ===========================================================================
// BUG-49: Triangle could loop forever instead of failing
//
// Triangle inserts the segments of a PSLG one after the other and, where two of them cross,
// a vertex that it computes in floating point; for nearly collinear, overlapping segments
// through other vertices (a random triangle soup on a coarse grid, data/bugs/bug49_*.obj,
// tools/fuzz_soup.py --nice --seed 71 --tris 25: 35 points and 19 segments in one face,
// data/bugs/bug49_pslg.txt, drawn by tools/plot_pslg.py) that vertex can lie on the wrong side
// of a segment, which leaves the triangles around a vertex inconsistent.  finddirection() then
// turned around that vertex for ever.  The turns are now counted, so that this is a failure
// like the others (BUG-45) and the face is triangulated after noding its segments.  (A hang
// fails this test through the ctest time limit.)
// ===========================================================================

static void test_bug49_triangle_does_not_loop()
{
   const CSG::TriMesh clay = loadMesh("bugs/bug49_clay.obj");
   const CSG::TriMesh knife = loadMesh("bugs/bug49_knife.obj");
   CSG::TriMesh A, B;

   const std::string err = capturedStderr([&] { construct(clay, knife, true, A, B); });
   ASSERT(err.empty());   // no face was given up on ("could not retriangulate")
   for (const CSG::TriMesh* m : {&A, &B})
      for (const auto& v : m->vertices())
         ASSERT(std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]));
}


// ===========================================================================
// BUG-58: the weld radius was an absolute 1e-6, and a face's rotation into 2D had
// a threshold
//
// Points that are the same point are welded; the radius was 1e-6 whatever the size
// of the scene (with a margin of 5% on top of it), so a scene smaller than that
// collapsed into one point.  It is now a number of units in the last place of the
// scene's coordinates.  (The same bug, for the other tolerances of the code, is
// BUG-15.)  A cube of any size cut in half must give two closed halves.
//
// Likewise, a face whose normal was within 1e-10 of Z was projected into 2D as it
// was, and one within a few 1e-6 of it with a rotation computed from acos().  A
// cube turned by tiny angles, so that the faces that the knife cuts are nearly
// perpendicular to Z, must give halves of half its volume (the knife goes through
// the cube's centre, which divides a cube in two equal parts at any angle).
// ===========================================================================

static void test_bug58_weld_is_relative_to_the_scene_scale()
{
    for (double s : {1e-9, 1e-6, 1.0, 1e6})
    {
        const CSG::TriMesh clay = meshcheck::makeBox(Eigen::Vector3d::Zero(), Eigen::Vector3d::Constant(s));
        const CSG::TriMesh knife = meshcheck::makeTriangles({
            {-3 * s, -3 * s, 0.5 * s,   6 * s, -3 * s, 0.5 * s,   -3 * s, 6 * s, 0.5 * s}});
        CSG::TriMesh A, B;
        construct(clay, knife, true, A, B);

        std::ostringstream where;
        where << "box of size " << s << ": ";
        for (const CSG::TriMesh* half : {&A, &B})
        {
            const Report r = analyze(*half, 1e-9 * s);
            if (!r.closed())
                throw std::runtime_error(where.str() + "a half is not a closed mesh");
            if (std::fabs(r.volume / (s * s * s) - 0.5) > 1e-6)
                throw std::runtime_error(where.str() + "a half does not have half the volume");
        }
    }
}

static void test_bug58_nearly_axis_aligned_faces_are_cut()
{
    const Eigen::Vector3d centre = Eigen::Vector3d::Constant(0.5);
    const CSG::TriMesh unit = meshcheck::makeBox(Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones());
    // The knife: x = 0.5, facing +x, big enough to cover the cube in any orientation
    const CSG::TriMesh knife = meshcheck::makeTriangles({
        {0.5, -5, -5,   0.5, 5, -5,   0.5, 0, 8}});

    for (double flip : {0.0, M_PI})   // a face's normal may be near +Z or near -Z
    {
        for (double tilt : {0.0, 1e-12, 1e-10, 1e-9, 1e-6, 1e-3})
        {
            const Eigen::Matrix3d rot =
                (Eigen::AngleAxisd(flip, Eigen::Vector3d::UnitX()) *
                 Eigen::AngleAxisd(tilt, Eigen::Vector3d(1, 2, 0).normalized())).matrix();
            std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> v;
            for (const auto& p : unit.vertices())
                v.push_back(rot * (p - centre) + centre);
            std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> f;
            for (const auto& tri : unit.faces())
                f.push_back(Eigen::Vector3i(tri.m_v[0], tri.m_v[1], tri.m_v[2]));
            const CSG::TriMesh clay(v, f);

            CSG::TriMesh A, B;
            construct(clay, knife, true, A, B);

            std::ostringstream where;
            where << "cube turned by " << flip << " + " << tilt << ": ";
            for (const CSG::TriMesh* half : {&A, &B})
            {
                const Report r = analyze(*half);
                if (!r.closed())
                    throw std::runtime_error(where.str() + "a half is not a closed mesh");
                if (std::fabs(r.volume - 0.5) > 1e-9)
                    throw std::runtime_error(where.str() + "a half does not have half the volume");
            }
        }
    }
}


// ===========================================================================
// Command-line tools
// ===========================================================================

// csg_tests cut CLAY KNIFE [nocap]: for looking at the results of a cut
static int cutTool(int argc, char* argv[])
{
   if (argc != 2 && argc != 3)
   {
      std::cout << "Usage: csg_tests cut CLAY_OBJ_FILE KNIFE_OBJ_FILE [nocap]" << std::endl;
      return 1;
   }

   CSG::TriMesh clay, knife;
   clay.loadOBJ(argv[0]);
   knife.loadOBJ(argv[1]);

   const bool cap = !(argc == 3 && std::string(argv[2]) == "nocap");  // end-caps?
   CSG::TriMesh A, B;
   construct(clay, knife, cap, A, B);

   std::ofstream("A.obj") << A << std::endl;
   std::ofstream("B.obj") << B << std::endl;
   return 0;
}

// csg_tests aabb OBJ_FILE: builds and validates an AABB tree of a mesh
static int aabbTool(int argc, char* argv[])
{
   if (argc != 1)
   {
      std::cout << "Usage: csg_tests aabb OBJ_FILE" << std::endl;
      return 1;
   }

   std::cout << "Reading " << argv[0] << ".." << std::endl;
   CSG::TriMesh mesh;
   mesh.loadOBJ(argv[0]);

   AABBTree tree = mesh.createAABBTree();
   tree.validate();
   std::cout << "numObjects = " << tree.numObjects() << "  "
             << "height = " << tree.getHeight() << "  "
             << "node count = " << tree.getNodeCount() << std::endl;
   return 0;
}


// ===========================================================================
// main
// ===========================================================================

int main(int argc, char* argv[])
{
   if (argc >= 2 && std::string(argv[1]) == "cut")
      return cutTool(argc - 2, argv + 2);
   if (argc >= 2 && std::string(argv[1]) == "aabb")
      return aabbTool(argc - 2, argv + 2);

   for (int ii = 1; ii < argc; ++ii)
      g_filters.push_back(argv[ii]);

   // Loading
   run("load_clay_unit_cube",          test_load_clay_unit_cube);
   run("load_clay_tetrahedron",        test_load_clay_tetrahedron);
   run("load_all_knife_meshes",        test_load_all_knife_meshes);

   // AABB tree
   run("aabb_clay_unit_cube",          test_aabb_clay_unit_cube);
   run("aabb_two_mesh_intersect",      test_aabb_two_mesh_intersect);
   run("aabb_large_mesh",              test_aabb_tree_is_valid_for_a_large_mesh);
   run("aabb_point_welding",           test_aabb_point_welding);

   // CSG operations: standard cases
   run("horiz_cut_upward_normal",      test_horiz_cut_upward_normal);
   run("tilted_cut",                   test_tilted_cut);
   run("cube_cube_intersection",       test_cube_cube_intersection);
   run("horiz_cut_with_cap",           test_horiz_cut_with_cap);

   // CSG operations: early bug regressions
   run("horiz_cut_downward_normal_bug2",       test_horiz_cut_downward_normal_bug2_regression);
   run("cut_through_3_vertices_bug1",          test_cut_through_3_vertices_bug1_regression);
   run("cut_through_edge_bug1_bug4",           test_cut_through_edge_bug1_bug4_regression);
   run("coplanar_knife_no_crash",              test_coplanar_knife_no_crash);
   run("near_miss_vertex",                     test_near_miss_vertex);
   run("multi_cut_same_face_bugD",             test_multi_cut_same_face_bugD_regression);

   // Alternative clay mesh
   run("tetrahedron_horiz_cut",        test_tetrahedron_horiz_cut);
   run("tetrahedron_tilted_cut",       test_tetrahedron_tilted_cut);

   // Regressions for the bugs found while reworking the face classification
   run("bug09_touching_faces_are_kept",          test_bug09_touching_faces_are_kept);
   run("bug10_cutting_face_lookup",              test_bug10_cutting_face_lookup);
   run("untouched_components_are_kept",          test_untouched_components_are_kept);
   run("bug53_knife_vertices_on_clay_vertices",  test_bug53_knife_vertices_on_clay_vertices);
   run("bug14_crossing_cuts_in_a_face",          test_bug14_crossing_cuts_in_a_face);
   run("bug54_point_contact",                    test_bug54_point_contact);
   run("bug58_weld_is_relative_to_the_scene_scale",  test_bug58_weld_is_relative_to_the_scene_scale);
   run("bug58_nearly_axis_aligned_faces_are_cut",    test_bug58_nearly_axis_aligned_faces_are_cut);

   // More bugs from the tracker: the fixed ones
   run("bug05_tetrahedron_faces_point_outwards",  test_bug05_tetrahedron_faces_point_outwards);
   run("bug06_construct_does_not_leak",           test_bug06_construct_does_not_leak);
   run("bug08_aabb_intersect_with_arbitrary_keys", test_bug08_aabb_intersect_with_arbitrary_keys);
   run("bug12_cap_seam_is_welded",                test_bug12_cap_seam_is_welded);
   run("bug13_cap_assembly_indices",              test_bug13_cap_assembly_indices);
   run("bug18_unreached_shells_are_classified",   test_bug18_unreached_shells_are_classified);
   run("bug20_two_cutting_faces",                 test_bug20_two_cutting_faces);
   run("bug23_construct_is_repeatable",           test_bug23_construct_is_repeatable);
   run("bug27_coplanar_same_orientation",         test_bug27_coplanar_same_orientation);
   run("bug29_face_normal_is_unit",               test_bug29_face_normal_is_unit);
   run("bug37_degenerate_input_faces",            test_bug37_degenerate_input_faces);
   run("bug42_retriangulation_scales_linearly",   test_bug42_retriangulation_scales_linearly);

   // More bugs from the tracker: the ones not fixed here (crashes and hangs included)
   runKnownBug("bug15_scale_invariance", "BUG-15", test_bug15_scale_invariance);
   run("bug17_weld_does_not_reach_too_far", test_bug17_weld_does_not_reach_too_far);
   run("bug18_knife_misses_clay", test_bug18_knife_misses_clay);
   runKnownBug("bug19_nearly_coplanar_knife", "BUG-19", test_bug19_nearly_coplanar_knife, 2);
   run("bug21_empty_aabb_tree", test_bug21_empty_aabb_tree);
   run("bug21_knife_far_from_clay_does_not_crash", test_bug21_knife_far_from_clay_does_not_crash);
   runKnownBug("bug27_coplanar_knife_is_watertight", "BUG-27", test_bug27_coplanar_knife_is_watertight, 2);
   runKnownBug("bug27_coplanar_opposite_orientation", "BUG-27", test_bug27_coplanar_opposite_orientation);
   runKnownBug("bug28_engine_owns_its_meshes", "BUG-28", test_bug28_engine_owns_its_meshes, 2);
   runKnownBug("bug30_obj_round_trip", "BUG-30", test_bug30_obj_round_trip);
   runKnownBug("bug31_obj_loader", "BUG-31", test_bug31_obj_loader);
   runKnownBug("bug32_invalid_meshes_are_rejected", "BUG-32", test_bug32_invalid_meshes_are_rejected);
   runKnownBug("bug33_pair_hash_distribution", "BUG-33", test_bug33_pair_hash_distribution);
   runKnownBug("bug35_unimplemented_operations_are_refused", "BUG-35", test_bug35_unimplemented_operations_are_refused);
   runKnownBug("bug36_t_junction_on_shared_edge", "BUG-36", test_bug36_t_junction_on_shared_edge);
   runKnownBug("bug38_touching_spheres", "BUG-38", test_bug38_touching_spheres);
   runKnownBug("bug39_quiet_by_default", "BUG-39", test_bug39_quiet_by_default);
   runKnownBug("bug40_non_finite_coordinates_are_rejected", "BUG-40", test_bug40_non_finite_coordinates_are_rejected);
   runKnownBug("bug43_triangle_failure_does_not_exit", "BUG-43", test_bug43_triangle_failure_does_not_exit);
   runKnownBug("bug44_tolerance_respects_coordinate_resolution", "BUG-44", test_bug44_tolerance_respects_coordinate_resolution);
   runKnownBug("bug45_triangle_failure_is_recovered", "BUG-45", test_bug45_triangle_failure_is_recovered);
   runKnownBug("bug46_boundary_points_stay_inside_the_face", "BUG-46", test_bug46_boundary_points_stay_inside_the_face);
   runKnownBug("bug47_output_has_no_t_junctions", "BUG-47", test_bug47_output_has_no_t_junctions);
   runKnownBug("bug48_conforming_terminates", "BUG-48", test_bug48_conforming_terminates);
   runKnownBug("bug49_triangle_does_not_loop", "BUG-49", test_bug49_triangle_does_not_loop, 2);

   return summarize();
}
