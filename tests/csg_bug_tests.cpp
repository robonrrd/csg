// Regression tests for the bugs tracked in TODO.md.
//
// Every bug (fixed or not) has at least one test here.  Tests for bugs that
// have not been fixed yet are registered with runKnownBug(): they are expected
// to fail, and the runner reports an unexpected pass so the registration can
// be switched to run() as soon as the bug is fixed.
//
// Exit code is 0 when every test behaves as expected.

#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "libcsg.h"
#include "mesh_check.h"

#ifndef DATA_DIR
#define DATA_DIR "data"
#endif

using meshcheck::analyze;
using meshcheck::Report;

#define CHECK(cond) \
   do { if (!(cond)) { \
      std::ostringstream _os; \
      _os << "CHECK(" #cond ") failed at line " << __LINE__; \
      throw std::runtime_error(_os.str()); \
   }} while (0)

#define CHECK_NEAR(a, b, tol) \
   do { double _a = (a), _b = (b); if (!(std::fabs(_a - _b) <= (tol))) { \
      std::ostringstream _os; \
      _os << "CHECK_NEAR(" #a ", " #b "): " << _a << " vs " << _b << " at line " << __LINE__; \
      throw std::runtime_error(_os.str()); \
   }} while (0)

static std::string dataPath(const std::string& f) { return std::string(DATA_DIR) + "/" + f; }

static CSG::TriMesh loadMesh(const std::string& file)
{
   CSG::TriMesh m;
   CHECK(m.loadOBJ(dataPath(file)) > 0);
   return m;
}

static void construct(const CSG::TriMesh& clay, const CSG::TriMesh& knife, bool cap,
                      CSG::TriMesh& A, CSG::TriMesh& B)
{
   CSG::CSGEngine engine(clay, knife);
   engine.construct(CSG::kDifference, cap, A, B);
}

static void constructFiles(const std::string& clay, const std::string& knife, bool cap,
                           CSG::TriMesh& A, CSG::TriMesh& B)
{
   const CSG::TriMesh c = loadMesh(clay);
   const CSG::TriMesh k = loadMesh(knife);
   construct(c, k, cap, A, B);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

static int g_pass = 0, g_fail = 0, g_xfail = 0, g_xpass = 0;

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

// A test that documents a bug that is not fixed yet.  It must fail; if it passes
// the bug has been fixed and the test should be registered with run().
static void runKnownBug(const char* name, const char* bug, std::function<void()> fn)
{
   std::cout << "[TEST] " << name << "  (known bug " << bug << ")" << std::endl;
   try {
      fn();
      ++g_xpass;
      std::cout << "       UNEXPECTED PASS: " << bug << " seems fixed; use run() instead"
                << std::endl;
   } catch (const std::exception& e) {
      ++g_xfail;
      std::cout << "       XFAIL (" << bug << "): " << e.what() << std::endl;
   }
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
   CHECK_NEAR(analyze(A).area + analyze(B).area, 6.0, 1e-9);
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
   CHECK_NEAR(topA, 0.5, 1e-9);
   CHECK_NEAR(topB, 0.5, 1e-9);
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
   CHECK_NEAR(analyze(A).area, 9.0, 1e-9);
   CHECK_NEAR(analyze(B).area, 9.0, 1e-9);
}

// ===========================================================================
// Helpers for the tests below
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
   CHECK(ra.closed());
   CHECK(rb.closed());
   CHECK_NEAR(ra.volume, expected_volume_A, 1e-9);
   CHECK_NEAR(rb.volume, expected_volume_B, 1e-9);
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
   CHECK(both.closed());
   CHECK_NEAR(both.area, 6.0, 1e-9);
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
   CHECK(hasVertexAt(A, Eigen::Vector3d(0.7, 0.3, 1.0)) ||
         hasVertexAt(B, Eigen::Vector3d(0.7, 0.3, 1.0)));

   // Nothing is lost, and every crack is closed: A and B together are the
   // whole, watertight, surface of the cube
   const Report both = analyze(combine(A, B));
   CHECK(both.closed());
   CHECK_NEAR(both.area, 6.0, 1e-9);
   CHECK_NEAR(both.volume, 1.0, 1e-9);
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
   CHECK(both.closed());
   CHECK_NEAR(both.area, 6.0, 1e-9);
   CHECK_NEAR(both.volume, 1.0, 1e-9);

   // ...all on one side of the knife (the clay is entirely below it)
   CHECK(A.faces().empty() != B.faces().empty());
}

// ===========================================================================

int main()
{
   run("bug10_cutting_face_lookup", test_bug10_cutting_face_lookup);
   run("bug09_touching_faces_are_kept", test_bug09_touching_faces_are_kept);
   run("untouched_components_are_kept", test_untouched_components_are_kept);
   run("bug08_knife_vertices_on_clay_vertices", test_bug08_knife_vertices_on_clay_vertices);
   run("bug13_crossing_cuts_in_a_face", test_bug13_crossing_cuts_in_a_face);
   run("bug14_point_contact", test_bug14_point_contact);

   const int total = g_pass + g_fail + g_xfail + g_xpass;
   std::cout << "\n--- " << total << " tests: " << g_pass << " passed, " << g_xfail
             << " expected failures, " << g_fail << " FAILED, " << g_xpass
             << " unexpected passes ---\n";
   return (g_fail > 0 || g_xpass > 0) ? 1 : 0;
}
