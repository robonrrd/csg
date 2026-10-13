// Mesh validation helpers shared by the CSG test programs.
//
// The CSG library produces triangle soups that are only meaningful if they are
// topologically sound, so the tests check more than "the output is non-empty":
// vertices are welded by position, then every edge is counted so that cracks,
// non-manifold edges and inconsistent winding are all visible.
#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include <Eigen/Geometry>  // cross()

#include "libcsg.h"

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
