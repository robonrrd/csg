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

#include <fstream>
#include <random>

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
// BUG-8: points that are a vertex of both meshes
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

static void test_bug08_knife_vertices_on_clay_vertices()
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
// BUG-13: constraint segments that cross inside a face
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

static void test_bug13_crossing_cuts_in_a_face()
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
// BUG-14: a knife that only touches the clay at a point
//
// The apex of the knife triangle lies exactly on the interior of the cube's top
// face and the rest of the triangle is above the cube, so the two triangles
// intersect in a segment of zero length.  That is not a cut: all of the clay must
// come out, intact, on one side of the knife.
// ===========================================================================

static void test_bug14_point_contact()
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
   run("bug08_knife_vertices_on_clay_vertices",  test_bug08_knife_vertices_on_clay_vertices);
   run("bug13_crossing_cuts_in_a_face",          test_bug13_crossing_cuts_in_a_face);
   run("bug14_point_contact",                    test_bug14_point_contact);

   return summarize();
}
