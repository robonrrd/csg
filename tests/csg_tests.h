// Shared infrastructure for the CSG tests (see csg_tests.cpp).
//
//  - assertion macros and a small test runner
//  - helpers to load meshes and run the CSG engine on them
//  - meshcheck: mesh validation.  The CSG library produces triangle soups that
//    are only meaningful if they are topologically sound, so the tests check more
//    than "the output is non-empty": vertices are welded by position, then every
//    edge is counted so that cracks, non-manifold edges and inconsistent winding
//    are all visible.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Geometry>  // cross()

#include "libcsg.h"

#ifndef DATA_DIR
#define DATA_DIR "data"
#endif


// ---------------------------------------------------------------------------
// Assertions.  A failed assertion throws; the runner reports it as a failure.
// ---------------------------------------------------------------------------

static const double EPS = 1e-5;

#define ASSERT(cond) \
   do { if (!(cond)) { \
      std::ostringstream _os; \
      _os << "ASSERT(" #cond ") at line " << __LINE__; \
      throw std::runtime_error(_os.str()); \
   }} while (0)

#define ASSERT_GE(a, b) \
   do { double _a = (a), _b = (b); if (!(_a >= _b)) { \
      std::ostringstream _os; \
      _os << "ASSERT_GE: " << _a << " < " << _b << " at line " << __LINE__; \
      throw std::runtime_error(_os.str()); \
   }} while (0)

#define ASSERT_LE(a, b) \
   do { double _a = (a), _b = (b); if (!(_a <= _b)) { \
      std::ostringstream _os; \
      _os << "ASSERT_LE: " << _a << " > " << _b << " at line " << __LINE__; \
      throw std::runtime_error(_os.str()); \
   }} while (0)

#define ASSERT_NEAR(a, b, tol) \
   do { double _a = (a), _b = (b); if (!(std::fabs(_a - _b) <= (tol))) { \
      std::ostringstream _os; \
      _os << "ASSERT_NEAR(" #a ", " #b "): " << _a << " vs " << _b \
          << " at line " << __LINE__; \
      throw std::runtime_error(_os.str()); \
   }} while (0)


// ---------------------------------------------------------------------------
// Test runner
//
// Tests named on the command line (as substrings) are the only ones run;
// with no names, all of them are.
// ---------------------------------------------------------------------------

static int g_pass = 0, g_fail = 0, g_xfail = 0, g_xpass = 0;
static std::vector<std::string> g_filters;

inline bool testSelected(const std::string& name)
{
   if (g_filters.empty())
      return true;
   for (const std::string& f : g_filters)
      if (name.find(f) != std::string::npos)
         return true;
   return false;
}

inline void run(const char* name, std::function<void()> fn)
{
   if (!testSelected(name))
      return;
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
inline void runKnownBug(const char* name, const char* bug, std::function<void()> fn)
{
   if (!testSelected(name))
      return;
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

// Prints the totals; returns the process exit code (0 when every test behaved as
// expected)
inline int summarize()
{
   std::cout << "\n--- " << (g_pass + g_fail + g_xfail + g_xpass) << " tests: " << g_pass
             << " passed, " << g_xfail << " expected failures, " << g_fail << " FAILED, "
             << g_xpass << " unexpected passes ---\n";
   return (g_fail > 0 || g_xpass > 0) ? 1 : 0;
}


// ---------------------------------------------------------------------------
// Loading meshes and running the CSG engine
// ---------------------------------------------------------------------------

inline std::string dataPath(const std::string& f)
{
   return std::string(DATA_DIR) + "/" + f;
}

inline CSG::TriMesh loadMesh(const std::string& file)
{
   CSG::TriMesh m;
   ASSERT(m.loadOBJ(dataPath(file)) > 0);
   return m;
}

inline void construct(const CSG::TriMesh& clay, const CSG::TriMesh& knife, bool cap,
                      CSG::TriMesh& A, CSG::TriMesh& B)
{
   CSG::CSGEngine engine(clay, knife);
   engine.construct(CSG::kDifference, cap, A, B);
}

// Loads two OBJ files from the data directory and runs CSG::construct on them.
// Throws if either file fails to load.
inline void constructFiles(const std::string& clay, const std::string& knife, bool cap,
                           CSG::TriMesh& A, CSG::TriMesh& B)
{
   const CSG::TriMesh c = loadMesh(clay);
   const CSG::TriMesh k = loadMesh(knife);
   construct(c, k, cap, A, B);
}


// ---------------------------------------------------------------------------
// Basic mesh sanity
// ---------------------------------------------------------------------------

// Returns true if any vertex coordinate is non-finite (NaN or Inf).
inline bool hasNaN(const CSG::TriMesh& mesh)
{
   for (const auto& v : mesh.vertices())
      if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
         return true;
   return false;
}

// Returns false if any face vertex index is outside [0, numVertices).
inline bool hasValidFaceIndices(const CSG::TriMesh& mesh)
{
   const int32_t nv = static_cast<int32_t>(mesh.vertices().size());
   for (const auto& f : mesh.faces())
      for (int ii = 0; ii < 3; ++ii)
         if (f.m_v[ii] < 0 || f.m_v[ii] >= nv)
            return false;
   return true;
}

// Throws on any geometry/topology violation.
inline void validateMesh(const CSG::TriMesh& m)
{
   ASSERT(!hasNaN(m));
   ASSERT(hasValidFaceIndices(m));
}

inline Eigen::Vector3d meshMin(const CSG::TriMesh& m)
{
   Eigen::Vector3d lo(1e30, 1e30, 1e30);
   for (const auto& v : m.vertices())
      lo = lo.cwiseMin(v);
   return lo;
}

inline Eigen::Vector3d meshMax(const CSG::TriMesh& m)
{
   Eigen::Vector3d hi(-1e30, -1e30, -1e30);
   for (const auto& v : m.vertices())
      hi = hi.cwiseMax(v);
   return hi;
}


// ---------------------------------------------------------------------------
// Topological mesh validation
// ---------------------------------------------------------------------------

namespace meshcheck
{

struct Report
{
   size_t vertices = 0;
   size_t faces = 0;
   size_t duplicate_vertices = 0;  // vertices that coincide with an earlier vertex
   size_t degenerate_faces = 0;    // repeated (welded) vertex or ~zero area
   size_t boundary_edges = 0;      // edges used by exactly one face (cracks / holes)
   size_t nonmanifold_edges = 0;   // edges used by more than two faces
   size_t inconsistent_edges = 0;  // edges traversed in the same direction twice
   double volume = 0.0;            // signed volume (meaningful for closed meshes)
   double area = 0.0;

   // A closed, consistently oriented 2-manifold.
   bool closed() const
   {
      return faces > 0 && boundary_edges == 0 && nonmanifold_edges == 0 &&
             inconsistent_edges == 0 && degenerate_faces == 0;
   }
};

// Maps each vertex to the lowest-numbered vertex within 'tol' (per axis).
inline std::vector<uint32_t> weldVertices(const std::vector<Eigen::Vector3d,
                                          Eigen::aligned_allocator<Eigen::Vector3d>>& v,
                                          double tol)
{
   const size_t n = v.size();
   std::vector<uint32_t> order(n);
   for (size_t i = 0; i < n; ++i)
      order[i] = static_cast<uint32_t>(i);
   std::sort(order.begin(), order.end(),
             [&](uint32_t a, uint32_t b) { return v[a][0] < v[b][0]; });

   std::vector<uint32_t> rep(n);
   for (size_t i = 0; i < n; ++i)
      rep[i] = static_cast<uint32_t>(i);
   for (size_t ii = 0; ii < n; ++ii)
   {
      for (size_t jj = ii + 1; jj < n && v[order[jj]][0] - v[order[ii]][0] <= tol; ++jj)
      {
         const uint32_t a = order[ii], b = order[jj];
         if ((v[a] - v[b]).cwiseAbs().maxCoeff() <= tol)
         {
            // union by smallest index
            uint32_t ra = rep[a], rb = rep[b];
            while (rep[ra] != ra) ra = rep[ra];
            while (rep[rb] != rb) rb = rep[rb];
            if (ra < rb) rep[rb] = ra; else rep[ra] = rb;
         }
      }
   }
   for (size_t i = 0; i < n; ++i)
   {
      uint32_t r = rep[i];
      while (rep[r] != r) r = rep[r];
      rep[i] = r;
   }
   return rep;
}

inline double triangleArea(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                           const Eigen::Vector3d& c)
{
   return 0.5 * (b - a).cross(c - a).norm();
}

inline Report analyze(const CSG::TriMesh& mesh, double tol = 1e-6)
{
   Report r;
   const auto& verts = mesh.vertices();
   const auto& faces = mesh.faces();
   r.vertices = verts.size();
   r.faces = faces.size();

   const std::vector<uint32_t> rep = weldVertices(verts, tol);
   for (size_t i = 0; i < rep.size(); ++i)
      if (rep[i] != i)
         ++r.duplicate_vertices;

   std::map<std::pair<uint32_t, uint32_t>, int> directed;
   for (const auto& f : faces)
   {
      const uint32_t id[3] = {rep[f.m_v[0]], rep[f.m_v[1]], rep[f.m_v[2]]};
      const Eigen::Vector3d& p0 = verts[f.m_v[0]];
      const Eigen::Vector3d& p1 = verts[f.m_v[1]];
      const Eigen::Vector3d& p2 = verts[f.m_v[2]];

      const double area = triangleArea(p0, p1, p2);
      r.area += area;
      r.volume += p0.dot(p1.cross(p2)) / 6.0;

      if (id[0] == id[1] || id[1] == id[2] || id[0] == id[2] || area < 1e-12)
      {
         ++r.degenerate_faces;
         continue;  // a collapsed face contributes no meaningful edges
      }
      for (int e = 0; e < 3; ++e)
         ++directed[{id[e], id[(e + 1) % 3]}];
   }

   for (const auto& kv : directed)
   {
      const uint32_t a = kv.first.first, b = kv.first.second;
      if (a > b && directed.count({b, a}))
         continue;  // counted when we visit (b, a)
      const auto rev = directed.find({b, a});
      const int fwd = kv.second;
      const int bwd = (rev == directed.end()) ? 0 : rev->second;
      const int total = fwd + bwd;
      if (total == 1)
         ++r.boundary_edges;
      else if (total > 2)
         ++r.nonmanifold_edges;
      else if (fwd == 2 || bwd == 2)
         ++r.inconsistent_edges;
   }
   return r;
}

// Total area of faces whose vertices all lie within 'tol' of the plane n.x = d.
inline double areaOnPlane(const CSG::TriMesh& mesh, const Eigen::Vector3d& n, double d,
                          double tol = 1e-6)
{
   double total = 0.0;
   for (const auto& f : mesh.faces())
   {
      const Eigen::Vector3d& a = mesh.vertices()[f.m_v[0]];
      const Eigen::Vector3d& b = mesh.vertices()[f.m_v[1]];
      const Eigen::Vector3d& c = mesh.vertices()[f.m_v[2]];
      if (std::fabs(n.dot(a) - d) <= tol && std::fabs(n.dot(b) - d) <= tol &&
          std::fabs(n.dot(c) - d) <= tol)
         total += triangleArea(a, b, c);
   }
   return total;
}

// Builds an axis-aligned box [lo, hi] with outward-facing triangles.
inline CSG::TriMesh makeBox(const Eigen::Vector3d& lo, const Eigen::Vector3d& hi)
{
   std::vector<double> v = {
      lo[0], lo[1], lo[2],   hi[0], lo[1], lo[2],   hi[0], hi[1], lo[2],   lo[0], hi[1], lo[2],
      lo[0], lo[1], hi[2],   hi[0], lo[1], hi[2],   hi[0], hi[1], hi[2],   lo[0], hi[1], hi[2]};
   std::vector<unsigned int> f = {
      0, 2, 1,  0, 3, 2,   // z = lo
      4, 5, 6,  4, 6, 7,   // z = hi
      0, 1, 5,  0, 5, 4,   // y = lo
      3, 7, 6,  3, 6, 2,   // y = hi
      0, 4, 7,  0, 7, 3,   // x = lo
      1, 2, 6,  1, 6, 5};  // x = hi
   return CSG::TriMesh(v, f);
}

// Builds a mesh from an explicit list of triangles (9 doubles each).
inline CSG::TriMesh makeTriangles(const std::vector<std::vector<double>>& tris)
{
   std::vector<double> v;
   std::vector<unsigned int> f;
   for (const auto& t : tris)
   {
      for (int i = 0; i < 9; ++i)
         v.push_back(t[i]);
      const unsigned int base = static_cast<unsigned int>(f.size());
      f.push_back(base);
      f.push_back(base + 1);
      f.push_back(base + 2);
   }
   return CSG::TriMesh(v, f);
}

}  // namespace meshcheck
