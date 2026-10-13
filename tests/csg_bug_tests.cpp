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

int main()
{
   runKnownBug("bug09_touching_faces_are_kept", "BUG-9", test_bug09_touching_faces_are_kept);

   const int total = g_pass + g_fail + g_xfail + g_xpass;
   std::cout << "\n--- " << total << " tests: " << g_pass << " passed, " << g_xfail
             << " expected failures, " << g_fail << " FAILED, " << g_xpass
             << " unexpected passes ---\n";
   return (g_fail > 0 || g_xpass > 0) ? 1 : 0;
}
