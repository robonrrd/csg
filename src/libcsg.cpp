//
//
//
//
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <functional>
#include <iomanip>      // std::setprecision
#include <limits>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <Eigen/Dense>
#include "libcsg.h"
extern "C"
{
#include "triangle.h"
}


//#define DEBUG

namespace CSG
{

// Magic numbers and tolerances.
// We seek to minimize the number of special constants as they create edge cases
// and problems. All magic numbers are expressed in terms of 'units in the last
// place'
constexpr uint32_t POINT_ULP = 2; // the threshhold for point equality

// Two points of the scene that are this many units in the last place (of the largest
// coordinate of the scene) apart, or less, are the same point and are welded.
constexpr uint32_t WELD_ULP = 1024;


//
// Helper functions
//

// A simple disjoint-set (union-find) structure, used to weld clusters
// of approximately-coincident vertices into a single canonical
// representative.  My previous "map each point to its nearest
// already-seen match" scheme was not transitive: for a chain of
// points A~B~C (A close to B, B close to C, but A NOT close enough to
// C), that scheme can leave A and C resolved to two different final
// vertices instead of one. This contributes to triangle sliver problems
// further down the line.
// Union-find (with path compression, and union-by-smallest-index so
// that e.g. original mesh vertices are always preferred as the
// representative) guarantees every point in a connected cluster
// resolves to the same root, regardless of the order in which matches
// are discovered.
class UnionFind
{
 public:
   explicit UnionFind(size_t n) : m_parent(n)
   {
      for (size_t ii = 0; ii < n; ++ii)
         m_parent[ii] = static_cast<uint32_t>(ii);
   }

   uint32_t find(uint32_t x)
   {
      while (m_parent[x] != x)
      {
         m_parent[x] = m_parent[m_parent[x]];
         x = m_parent[x];
      }
      return x;
   }

   void unite(uint32_t a, uint32_t b)
   {
      const uint32_t ra = find(a);
      const uint32_t rb = find(b);
      if (ra == rb)
         return;
      // Union by smallest index, so lower-indexed (e.g. original mesh)
      // vertices are always preferred as the canonical representative.
      if (ra < rb)
         m_parent[rb] = ra;
      else
         m_parent[ra] = rb;
   }

 private:
   std::vector<uint32_t> m_parent;
};


// Unites the points of 'points' that are within 'tolerance' of each other
// (in every coordinate).  Only the points listed in 'queries' look for neighbours, and
// they unite with a neighbour only if 'wanted(point, neighbour)'.
static void uniteClosePoints(const std::vector<VECTOR3D>& points, double tolerance,
                             const std::vector<uint32_t>& queries,
                             const std::function<bool(uint32_t, uint32_t)>& wanted,
                             UnionFind& uf)
{
   if (points.empty())
      return;  // nothing to weld (and an AABBTree cannot be built from nothing)

   // The tree finds the points whose boxes (of half the tolerance around each point)
   // overlap; no margin is added to the boxes, and the distance is checked below
   AABBTree tree(0.0, points.size());
   for (uint32_t ii = 0; ii < points.size(); ++ii)
      tree.addSphere(ii, points[ii], 0.5 * tolerance);

   for (uint32_t ii : queries)
   {
      for (uint32_t jj : tree.query(ii))
      {
         if (wanted(ii, jj) && ((points[ii] - points[jj]).cwiseAbs().maxCoeff() <= tolerance))
            uf.unite(ii, jj);
      }
   }
}


template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, uint32_t ulp)
{
   // The machine epsilon has to be scaled to the magnitude of the values used
   // and multiplied by the desired precision in ULPs (units in the last place)
   return std::fabs(x - y) <= std::numeric_limits<T>::epsilon() * std::fabs(x + y) * ulp
          // unless the result is subnormal
          || std::fabs(x - y) < std::numeric_limits<T>::min();
}


bool point_almost_equal(const Eigen::Vector3d& p, const Eigen::Vector3d& q)
{
   return almost_equal(p[0], q[0], POINT_ULP) && almost_equal(p[1], q[1], POINT_ULP) &&
          almost_equal(p[2], q[2], POINT_ULP);
}


inline double orient2d(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c)
{
   return ((a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]));
}

//
// Triangle-triangle intersection functions
//
bool intersectionTestVertex(const Eigen::Vector2d& P1, const Eigen::Vector2d& Q1,
                            const Eigen::Vector2d& R1, const Eigen::Vector2d& P2,
                            const Eigen::Vector2d& Q2, const Eigen::Vector2d& R2)
{
   if (orient2d(R2, P2, Q1) >= 0.0)
      if (orient2d(R2, Q2, Q1) <= 0.0)
         if (orient2d(P1, P2, Q1) > 0.0)
         {
            if (orient2d(P1, Q2, Q1) <= 0.0)
               return true;
            else
               return false;
         }
         else
         {
            if (orient2d(P1, P2, R1) >= 0.0)
               if (orient2d(Q1, R1, P2) >= 0.0)
                  return true;
               else
                  return false;
            else
               return false;
         }
      else if (orient2d(P1, Q2, Q1) <= 0.0)
         if (orient2d(R2, Q2, R1) <= 0.0)
            if (orient2d(Q1, R1, Q2) >= 0.0)
               return true;
            else
               return false;
         else
            return false;
      else
         return false;
   else if (orient2d(R2, P2, R1) >= 0.0)
      if (orient2d(Q1, R1, R2) >= 0.0)
         if (orient2d(P1, P2, R1) >= 0.0)
            return true;
         else
            return false;
      else if (orient2d(Q1, R1, Q2) >= 0.0)
      {
         if (orient2d(R2, R1, Q2) >= 0.0)
            return true;
         else
            return false;
      }
      else
         return false;
   else
      return false;
}


bool intersectionTestEdge(const Eigen::Vector2d& P1, const Eigen::Vector2d& Q1,
                          const Eigen::Vector2d& R1, const Eigen::Vector2d& P2,
                          const Eigen::Vector2d& Q2, const Eigen::Vector2d& R2)
{
   if (orient2d(R2, P2, Q1) >= 0.0)
   {
      if (orient2d(P1, P2, Q1) >= 0.0)
      {
         if (orient2d(P1, Q1, R2) >= 0.0)
            return true;
         else
            return false;
      }
      else
      {
         if (orient2d(Q1, R1, P2) >= 0.0)
         {
            if (orient2d(R1, P1, P2) >= 0.0)
               return true;
            else
               return false;
         }
         else
            return false;
      }
   }
   else
   {
      if (orient2d(R2, P2, R1) >= 0.0)
      {
         if (orient2d(P1, P2, R1) >= 0.0)
         {
            if (orient2d(P1, R1, R2) >= 0.0)
               return true;
            else
            {
               if (orient2d(Q1, R1, R2) >= 0.0)
                  return true;
               else
                  return false;
            }
         }
         else
            return false;
      }
      else
         return false;
   }
}


bool ccwTriTriIntersection2d(const Eigen::Vector2d& p1, const Eigen::Vector2d& q1,
                             const Eigen::Vector2d& r1, const Eigen::Vector2d& p2,
                             const Eigen::Vector2d& q2, const Eigen::Vector2d& r2)
{
   if (orient2d(p2, q2, p1) >= 0.0)
   {
      if (orient2d(q2, r2, p1) >= 0.0)
      {
         if (orient2d(r2, p2, p1) >= 0.0)
            return true;
         else
            return intersectionTestEdge(p1, q1, r1, p2, q2, r2);
      }
      else
      {
         if (orient2d(r2, p2, p1) >= 0.0)
            return intersectionTestEdge(p1, q1, r1, r2, p2, q2);
         else
            return intersectionTestVertex(p1, q1, r1, p2, q2, r2);
      }
   }
   else
   {
      if (orient2d(q2, r2, p1) >= 0.0)
      {
         if (orient2d(r2, p2, p1) >= 0.0)
            return intersectionTestEdge(p1, q1, r1, q2, r2, p2);
         else
            return intersectionTestVertex(p1, q1, r1, q2, r2, p2);
      }
      else
         return intersectionTestVertex(p1, q1, r1, r2, p2, q2);
   }
}


bool triTriOverlapTest2d(const Eigen::Vector2d& p1, const Eigen::Vector2d& q1,
                         const Eigen::Vector2d& r1, const Eigen::Vector2d& p2,
                         const Eigen::Vector2d& q2, const Eigen::Vector2d& r2)
{
   if (orient2d(p1, q1, r1) < 0.0)
      if (orient2d(p2, q2, r2) < 0.0)
         return ccwTriTriIntersection2d(p1, r1, q1, p2, r2, q2);
      else
         return ccwTriTriIntersection2d(p1, r1, q1, p2, q2, r2);
   else if (orient2d(p2, q2, r2) < 0.0)
      return ccwTriTriIntersection2d(p1, q1, r1, p2, r2, q2);
   else
      return ccwTriTriIntersection2d(p1, q1, r1, p2, q2, r2);
}


bool coplanarTriTri3d(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1,
                      const Eigen::Vector3d& r1, const Eigen::Vector3d& p2,
                      const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
                      const Eigen::Vector3d& normal_1, const Eigen::Vector3d& normal_2)
{
   Eigen::Vector2d P1, Q1, R1;
   Eigen::Vector2d P2, Q2, R2;

   double n_x, n_y, n_z;
   n_x = ((normal_1[0] < 0) ? -normal_1[0] : normal_1[0]);
   n_y = ((normal_1[1] < 0) ? -normal_1[1] : normal_1[1]);
   n_z = ((normal_1[2] < 0) ? -normal_1[2] : normal_1[2]);


   /* Projection of the triangles in 3D onto 2D such that the area of
      the projection is maximized. */

   if ((n_x > n_z) && (n_x >= n_y))
   {
      // Project onto plane YZ
      P1[0] = q1[2];
      P1[1] = q1[1];
      Q1[0] = p1[2];
      Q1[1] = p1[1];
      R1[0] = r1[2];
      R1[1] = r1[1];

      P2[0] = q2[2];
      P2[1] = q2[1];
      Q2[0] = p2[2];
      Q2[1] = p2[1];
      R2[0] = r2[2];
      R2[1] = r2[1];
   }
   else if ((n_y > n_z) && (n_y >= n_x))
   {
      // Project onto plane XZ
      P1[0] = q1[0];
      P1[1] = q1[2];
      Q1[0] = p1[0];
      Q1[1] = p1[2];
      R1[0] = r1[0];
      R1[1] = r1[2];

      P2[0] = q2[0];
      P2[1] = q2[2];
      Q2[0] = p2[0];
      Q2[1] = p2[2];
      R2[0] = r2[0];
      R2[1] = r2[2];
   }
   else
   {
      // Project onto plane XY
      P1[0] = p1[0];
      P1[1] = p1[1];
      Q1[0] = q1[0];
      Q1[1] = q1[1];
      R1[0] = r1[0];
      R1[1] = r1[1];

      P2[0] = p2[0];
      P2[1] = p2[1];
      Q2[0] = q2[0];
      Q2[1] = q2[1];
      R2[0] = r2[0];
      R2[1] = r2[1];
   }

   return triTriOverlapTest2d(P1, Q1, R1, P2, Q2, R2);
};


// This function is called when the triangles are known to intersect. It
// constructs the segment of intersection of the two triangles if they are not
// coplanar.
bool constructIntersection(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1,
                           const Eigen::Vector3d& r1, const Eigen::Vector3d& p2,
                           const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
                           const Eigen::Vector3d& N1, const Eigen::Vector3d& N2,
                           Eigen::Vector3d& source, Eigen::Vector3d& target,
                           double& alpha, double& beta)
{
   Eigen::Vector3d v1 = q1 - p1;
   Eigen::Vector3d v2 = r2 - p1;
   Eigen::Vector3d N = v1.cross(v2);
   Eigen::Vector3d v = p2 - p1;

   if (v.dot(N) > 0.0)
   {
      v1 = r1 - p1;
      N = v1.cross(v2);
      if (v.dot(N) <= 0.0)
      {
         v2 = q2 - p1;
         N = v1.cross(v2);
         if (v.dot(N) > 0.0)
         {
            v1 = p1 - p2;
            v2 = p1 - r1;
            alpha = v1.dot(N2) / v2.dot(N2);
            v1 = alpha * v2;
            source = p1 - v1;
            v1 = p2 - p1;
            v2 = p2 - r2;
            beta = v1.dot(N1) / v2.dot(N1);
            v1 = beta * v2;
            target = p2 - v1;
            return true;
         }
         else
         {
            v1 = p2 - p1;
            v2 = p2 - q2;
            alpha = v1.dot(N1) / v2.dot(N1);
            v1 = alpha * v2;
            source = p2 - v1;
            v1 = p2 - p1;
            v2 = p2 - r2;
            beta = v1.dot(N1) / v2.dot(N1);
            v1 = beta * v2;
            target = p2 - v1;
            return true;
         }
      }
      else
      {
         return false;
      }
   }
   else
   {
      v2 = q2 - p1;
      N = v1.cross(v2);
      if (v.dot(N) < 0.0)
      {
         return false;
      }
      else
      {
         v1 = r1 - p1;
         N = v1.cross(v2);
         if (v.dot(N) >= 0.0)
         {
            v1 = p1 - p2;
            v2 = p1 - r1;
            alpha = v1.dot(N2) / v2.dot(N2);
            v1 = alpha * v2;
            source = p1 - v1;
            v1 = p1 - p2;
            v2 = p1 - q1;
            beta = v1.dot(N2) / v2.dot(N2);
            v1 = beta * v2;
            target = p1 - v1;
            return true;
         }
         else
         {
            v1 = p2 - p1;
            v2 = p2 - q2;
            alpha = v1.dot(N1) / v2.dot(N1);
            v1 = alpha * v2;
            source = p2 - v1;
            v1 = p1 - p2;
            v2 = p1 - q1;
            beta = v1.dot(N2) / v2.dot(N2);
            v1 = beta * v2;
            target = p1 - v1;
            return true;
         }
      }
   }
}


bool triTriInter3d(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1, const Eigen::Vector3d& r1,
                   const Eigen::Vector3d& p2, const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
                   const Eigen::Vector3d& N1, const Eigen::Vector3d& N2,
                   Eigen::Vector3d& source, Eigen::Vector3d& target, double& alpha, double& beta,
                   const double dp2, const double dq2, const double dr2, bool& coplanar)
{
   alpha = -1.0;
   beta = -1.0;

   if (dp2 > 0.0)
   {
      if (dq2 > 0.0)
         return constructIntersection(p1, r1, q1, r2, p2, q2, N1, N2, source, target, alpha, beta);
      else if (dr2 > 0.0)
         return constructIntersection(p1, r1, q1, q2, r2, p2, N1, N2, source, target, alpha, beta);
      else
         return constructIntersection(p1, q1, r1, p2, q2, r2, N1, N2, source, target, alpha, beta);
   }
   else if (dp2 < 0.0)
   {
      if (dq2 < 0.0)
         return constructIntersection(p1, q1, r1, r2, p2, q2, N1, N2, source, target, alpha, beta);
      else if (dr2 < 0.0)
         return constructIntersection(p1, q1, r1, q2, r2, p2, N1, N2, source, target, alpha, beta);
      else
         return constructIntersection(p1, r1, q1, p2, q2, r2, N1, N2, source, target, alpha, beta);
   }
   else
   {
      if (dq2 < 0.0)
      {
         if (dr2 >= 0.0)
            return constructIntersection(p1, r1, q1, q2, r2, p2, N1, N2, source, target,
                                         alpha, beta);
         else
            return constructIntersection(p1, q1, r1, p2, q2, r2, N1, N2, source, target,
                                         alpha, beta);
      }
      else if (dq2 > 0.0)
      {
         if (dr2 > 0.0)
            return constructIntersection(p1, r1, q1, p2, q2, r2, N1, N2, source, target,
                                         alpha, beta);
         else
            return constructIntersection(p1, q1, r1, q2, r2, p2, N1, N2, source, target,
                                         alpha, beta);
      }
      else
      {
         if (dr2 > 0.0)
            return constructIntersection(p1, q1, r1, r2, p2, q2, N1, N2, source, target,
                                         alpha, beta);
         else if (dr2 < 0.0)
            return constructIntersection(p1, r1, q1, r2, p2, q2, N1, N2, source, target,
                                         alpha, beta);
         else
         {
            coplanar = true;
            return coplanarTriTri3d(p1, q1, r1, p2, q2, r2, N1, N2);
         }
      }
   }
}


//   The following version computes the segment of intersection of the two
//   triangles if it exists.  coplanar returns whether the triangles are
//   coplanar source and target are the endpoints of the line segment of
//   intersection
//
bool triTriIntersectionTest3d(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1,
                              const Eigen::Vector3d& r1, const Eigen::Vector3d& p2,
                              const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
                              bool& coplanar, Eigen::Vector3d& source, Eigen::Vector3d& target,
                              double& alpha, double& beta)
{
   double dp1, dq1, dr1, dp2, dq2, dr2;
   Eigen::Vector3d v1, v2, v;
   Eigen::Vector3d N1, N2, N;

   coplanar = false;
   alpha = -1.0;
   beta = -1.0;

   // Compute distance signs  of p1, q1 and r1 to the plane of triangle(p2,q2,r2)
   v1 = p2 - r2;
   v2 = q2 - r2;
   N2 = v1.cross(v2);

   v1 = p1 - r2;
   dp1 = v1.dot(N2);
   v1 = q1 - r2;
   dq1 = v1.dot(N2);
   v1 = r1 - r2;
   dr1 = v1.dot(N2);

   if (((dp1 * dq1) > 0.0) && ((dp1 * dr1) > 0.0))
      return false;

   // Compute distance signs  of p2, q2 and r2 to the plane of triangle(p1,q1,r1)
   v1 = q1 - p1;
   v2 = r1 - p1;
   N1 = v1.cross(v2);

   v1 = p2 - r1;
   dp2 = v1.dot(N1);
   v1 = q2 - r1;
   dq2 = v1.dot(N1);
   v1 = r2 - r1;
   dr2 = v1.dot(N1);

   if (((dp2 * dq2) > 0.0) && ((dp2 * dr2) > 0.0))
      return false;

   // Permutation in a canonical form of T1's vertices
   if (dp1 > 0.0)
   {
      if (dq1 > 0.0)
         return triTriInter3d(r1, p1, q1, p2, r2, q2, N1, N2, source, target, alpha, beta,
                              dp2, dr2, dq2, coplanar);
      else if (dr1 > 0.0)
         return triTriInter3d(q1, r1, p1, p2, r2, q2, N1, N2, source, target, alpha, beta,
                              dp2, dr2, dq2, coplanar);
      else
         return triTriInter3d(p1, q1, r1, p2, q2, r2, N1, N2, source, target, alpha, beta,
                              dp2, dq2, dr2, coplanar);
   }
   else if (dp1 < 0.0)
   {
      if (dq1 < 0.0)
         return triTriInter3d(r1, p1, q1, p2, q2, r2, N1, N2, source, target, alpha, beta,
                              dp2, dq2, dr2, coplanar);
      else if (dr1 < 0.0)
         return triTriInter3d(q1, r1, p1, p2, q2, r2, N1, N2, source, target, alpha, beta,
                              dp2, dq2, dr2, coplanar);
      else
         return triTriInter3d(p1, q1, r1, p2, r2, q2, N1, N2, source, target, alpha, beta,
                              dp2, dr2, dq2, coplanar);
   }
   else
   {
      if (dq1 < 0.0)
      {
         if (dr1 >= 0.0)
            return triTriInter3d(q1, r1, p1, p2, r2, q2, N1, N2, source, target, alpha, beta,
                                 dp2, dr2, dq2, coplanar);
         else
            return triTriInter3d(p1, q1, r1, p2, q2, r2, N1, N2, source, target, alpha, beta,
                                 dp2, dq2, dr2, coplanar);
      }
      else if (dq1 > 0.0)
      {
         if (dr1 > 0.0)
            return triTriInter3d(p1, q1, r1, p2, r2, q2, N1, N2, source, target, alpha, beta,
                                 dp2, dr2, dq2, coplanar);
         else
            return triTriInter3d(q1, r1, p1, p2, q2, r2, N1, N2, source, target, alpha, beta,
                                 dp2, dq2, dr2, coplanar);
      }
      else
      {
         if (dr1 > 0.0)
            return triTriInter3d(r1, p1, q1, p2, q2, r2, N1, N2, source, target, alpha, beta,
                                 dp2, dq2, dr2, coplanar);
         else if (dr1 < 0.0)
            return triTriInter3d(r1, p1, q1, p2, r2, q2, N1, N2, source, target, alpha, beta,
                                 dp2, dr2, dq2, coplanar);
         else
         {
            // triangles are co-planar
            coplanar = true;
            return coplanarTriTri3d(p1, q1, r1, p2, q2, r2, N1, N2);
         }
      }
   }
}


// Convenience wrapper around triTriIntersectionTest3d
TriangleIntersection intersect(
    const Triangle& tri_a, const std::vector<VECTOR3D>& verts_a,
    const Triangle& tri_b, const std::vector<VECTOR3D>& verts_b)
{
   TriangleIntersection retval;

   retval.intersect = triTriIntersectionTest3d(
       verts_a[tri_a.m_v[0]], verts_a[tri_a.m_v[1]], verts_a[tri_a.m_v[2]],
       verts_b[tri_b.m_v[0]], verts_b[tri_b.m_v[1]], verts_b[tri_b.m_v[2]],
       retval.coplanar, retval.p[0], retval.p[1], retval.alpha, retval.beta);

   return retval;
}


// Using the 3d position of an intersection point, determine if it is exactly
// a vertex on one (or both) of the intersecting triangles. If so, we track this
// information
std::vector<IPoint> CSGEngine::convertIntersectionToIpoints(const TriangleIntersection& ix,
                                                            uint32_t cidx, uint32_t kidx)
{
   const Triangle& ct = m_clay.faces()[cidx];
   const Triangle& kt = m_knife.faces()[kidx];

   // Fill out non-parent fields
   std::vector<IPoint> pts(2);
   for (uint32_t ii = 0; ii < 2; ++ii)
   {
      pts[ii].cidx = cidx;
      pts[ii].kidx = kidx;
   }


   // Determine the 'parenthood': is this point an exact duplicate of a point on
   // the clay mesh, the knife mesh or neither?
   for (uint32_t jj = 0; jj < 2; ++jj)
   {
      // Is it a point on the clay face?
      uint32_t clay_v_idx = 4;
      for (uint32_t ii = 0; ii < 3; ++ii)
      {
         if (point_almost_equal(ix.p[jj], m_clay.vertices()[ct.m_v[ii]]))
         {
            clay_v_idx = ii;
            break;
         }
      }

      // Is it a point on the knife face?
      uint32_t knife_v_idx = 4;
      for (uint32_t ii = 0; ii < 3; ++ii)
      {
         if (point_almost_equal(ix.p[jj], m_knife.vertices()[kt.m_v[ii]]))
         {
            knife_v_idx = ii;
            break;
         }
      }

      if ((clay_v_idx < 4) && (knife_v_idx == 4))
      {
         // This point is a clay mesh vertex
         pts[jj].ref.parent = kClay;
         pts[jj].ref.idx = ct.m_v[clay_v_idx];
      }
      else if ((clay_v_idx == 4) && (knife_v_idx < 4))
      {
         // This point is a knife mesh vertex
         pts[jj].ref.parent = kKnife;
         pts[jj].ref.idx = kt.m_v[knife_v_idx];
      }
      else if ((clay_v_idx < 4) && (knife_v_idx < 4))
      {
         // This point is a vertex on both clay and knife; store the clay global index
         pts[jj].ref.parent = kBoth;
         pts[jj].ref.idx = ct.m_v[clay_v_idx];
      }
      else
      {
         // A new point
         m_newPointPositions.push_back(ix.p[jj]);
         pts[jj].ref.parent = kNew;
         pts[jj].ref.idx = m_newPointPositions.size()-1;
      }
   }

   return pts;
}


// Returns a point index that spans clay vertices, knife vertices, and
// new vertices:
// 0...N   - Clay vertices
// N+1...M - Knife vertices
// M+1...L - New vertices
uint32_t CSGEngine::canonicalVertexIndex(const IPointRef& ref) const
{
   uint32_t idx = ref.idx;
   if ((ref.parent == kClay) || (ref.parent == kBoth))
      return idx;
   idx += m_clay.vertices().size();
   if (ref.parent == kKnife)
      return idx;

   return idx + m_knife.vertices().size();
}


// A key identifying the undirected edge between two vertex indices
inline uint64_t edgeKey(uint32_t a, uint32_t b)
{
   if (a > b)
      std::swap(a, b);
   return (static_cast<uint64_t>(a) << 32) | b;
}


void initializeTriangulateio(struct triangulateio& io)
{
   io.pointlist = 0;
   io.pointattributelist = 0;
   io.pointmarkerlist = 0;
   io.numberofpoints = 0;
   io.numberofpointattributes = 0;
   io.trianglelist = 0;
   io.triangleattributelist = 0;
   io.trianglearealist = 0;
   io.neighborlist = 0;
   io.numberoftriangles = 0;
   io.numberofcorners = 0;
   io.numberoftriangleattributes = 0;
   io.segmentlist = 0;
   io.segmentmarkerlist = 0;
   io.numberofsegments = 0;
   io.holelist = 0;
   io.numberofholes = 0;
   io.regionlist = 0;
   io.numberofregions = 0;
   io.edgelist = 0;
   io.edgemarkerlist = 0;
   io.normlist = 0;
   io.numberofedges = 0;
}


std::unordered_set<std::pair<uint32_t, uint32_t>> intersectAABBs(const TriMesh& clay,
                                                                 const TriMesh& knife)
{
   std::cout << "Building AABB trees.." << std::endl;
   AABBTree clay_tree = clay.createAABBTree();
   AABBTree knife_tree = knife.createAABBTree();
   std::cout << " ..done" << std::endl;

   // Calculate potential triangle intersections
   std::unordered_set<std::pair<uint32_t, uint32_t>> ix = clay_tree.intersect(knife_tree);
   std::cout << ix.size() << " AABB collisions found" << std::endl;

   return ix;
}


IParent otherMesh(IParent which)
{
   if (which == kClay)
      return kKnife;
   if (which == kKnife)
      return kClay;

   assert(false);
}


const Eigen::Vector3d centroid(const TriMesh& mesh, uint32_t ff)
{
   const Triangle& face = mesh.faces()[ff];
   return (mesh.vertices()[face.m_v[0]] + mesh.vertices()[face.m_v[1]] +
           mesh.vertices()[face.m_v[2]])/3.0;
}


// The shortest altitude of a triangle.  It is the altitude onto the longest
// edge, 2*area/longest, with the area taken from the cross product's norm.
// (Heron's formula is not usable here: for sliver triangles it subtracts nearly
// equal numbers, can go slightly negative, and then yields NaN.)
double triangleAltitude(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1,
                        const Eigen::Vector3d& p2)
{
   const double longest = std::sqrt(std::max(std::max((p0-p1).squaredNorm(),
                                                       (p1-p2).squaredNorm()),
                                             (p2-p0).squaredNorm()));
   if (longest == 0.0)
      return 0.0;

   return (p1-p0).cross(p2-p0).norm() / longest;
}

// Ratio between the longest edge and the shortest altitude. A well-shaped
// triangle has a ratio of order one; a sliver (or a collapsed triangle, for
// which the ratio is infinite) has a large one.
double triangleAspectRatio(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1,
                           const Eigen::Vector3d& p2)
{
   const double longest = std::sqrt(std::max(std::max((p0-p1).squaredNorm(),
                                                       (p1-p2).squaredNorm()),
                                             (p2-p0).squaredNorm()));
   const double altitude = triangleAltitude(p0, p1, p2);
   if (altitude == 0.0)
      return std::numeric_limits<double>::infinity();

   return longest / altitude;
}


// CSGEngine member functions
//

CSGEngine::CSGEngine(const TriMesh& in_clay, const TriMesh& in_knife)
   : m_clay(in_clay)
   , m_knife(in_knife)
{
   // empty
}


const Eigen::Vector3d& CSGEngine::ipointPos(const IPointRef& ref) const
{
   if ((ref.parent == kClay) || (ref.parent == kBoth))
      return m_clay.vertices()[ref.idx];
   else if (ref.parent == kKnife)
      return m_knife.vertices()[ref.idx];
   else
      return m_newPointPositions[ref.idx];
}

bool CSGEngine::degenerateTriangle(const IFace& face) const
{
   // Check for duplicate vertex indices (only useful after we've
   // done point welding or duplicate-finding)
   if ((face.v[0] == face.v[1]) || (face.v[0] == face.v[2]) || (face.v[1] == face.v[2]))
      return true;

   const Eigen::Vector3d& p0 = ipointPos(face.v[0]);
   const Eigen::Vector3d& p1 = ipointPos(face.v[1]);
   const Eigen::Vector3d& p2 = ipointPos(face.v[2]);

   // Both tests below are in units in the last place, like point equality: a
   // triangle is degenerate when its shape cannot be told from a line or a
   // point at the precision the coordinates carry.
   constexpr double ulp = std::numeric_limits<double>::epsilon() * POINT_ULP;

   // Squashed triangle: the shortest altitude is below the resolution of the
   // coordinates themselves, so the apex is within point-equality distance of
   // the line through the opposite edge.  (The comparison is written so that a
   // NaN altitude, from non-finite coordinates, also counts as degenerate.)
   const double coord_magnitude = std::max(p0.cwiseAbs().maxCoeff(),
                                  std::max(p1.cwiseAbs().maxCoeff(), p2.cwiseAbs().maxCoeff()));
   const double min_altitude = triangleAltitude(p0, p1, p2);
   if (!(min_altitude > ulp * coord_magnitude))
      return true;

   // Bad aspect ratio: the shortest altitude is below the resolution of the
   // triangle's own longest edge. This catches slivers far from the origin
   // or tiny triangles that the test above would let through.
   const double aspect_ratio = triangleAspectRatio(p0, p1, p2);
   if (!(aspect_ratio < 1.0 / ulp))
      return true;

   return false;
}

std::vector<IFace> CSGEngine::retriangulate(const TriMesh& mesh, IParent which_mesh, uint32_t fidx,
                                            const std::vector<uint32_t>& new_vert_indices)
{
   const uint32_t numPoints = 3 + new_vert_indices.size();
   const uint32_t numSegments = new_vert_indices.size() / 2;

#ifdef DEBUG
   std::cout << "numPoints=" << numPoints << "  numSegments=" << numSegments << std::endl;
#endif

   // Rotate the triangle into the XY plane ('triangle' is 2D only)
   const Eigen::Vector3d& p0 = mesh.vertices()[mesh.faces()[fidx].m_v[0]];
   const Eigen::Vector3d& p1 = mesh.vertices()[mesh.faces()[fidx].m_v[1]];
   const Eigen::Vector3d& p2 = mesh.vertices()[mesh.faces()[fidx].m_v[2]];

   const Eigen::Vector3d n = (p1 - p0).cross(p2 - p0).normalized();  // normal

   // A face without area has no plane to cut in
   if (n.squaredNorm() == 0.0)
      return std::vector<IFace>();

   // Rotate the triangle into the XY plane: the frame of the face has its x axis along
   // the face's first edge and its z axis along its normal, so that the normal goes to
   // +Z (and the winding is kept).  Unlike a rotation about the axis n x Z, this has no
   // special cases for a normal along Z.
   const Eigen::Vector3d x_axis = (p1 - p0).normalized();
   const Eigen::Vector3d y_axis = n.cross(x_axis).normalized();
   Eigen::Matrix3d rot;
   rot.row(0) = x_axis;
   rot.row(1) = y_axis;
   rot.row(2) = n;

   std::vector<VECTOR3D> pts_2d(numPoints);
   pts_2d[0] = rot * p0;
   pts_2d[1] = rot * p1;
   pts_2d[2] = rot * p2;
   for (uint32_t ii = 0; ii < new_vert_indices.size(); ++ii)
      pts_2d[3 + ii] = rot * ipointPos(m_newPoints[new_vert_indices[ii]].ref);


#ifdef DEBUG
   std::cout << "Retriangulating face " << fidx << std::endl
             << " Triangle points:" << std::endl
             << "  p0: " << p0.transpose() << std::endl
             << "  p1: " << p1.transpose() << std::endl
             << "  p2: " << p2.transpose() << std::endl
             << " Intersection points:" << std::endl;
   for (uint32_t ii = 0; ii < new_vert_indices.size(); ++ii)
      std::cout << "  i" << ii << ": "
                << ipointPos(m_newPoints[new_vert_indices[ii]].ref).transpose()
                << std::endl;

   std::cout << " 2d points" << std::endl;
   for (const auto& pt : pts_2d)
      std::cout << "  " << pt.transpose() << std::endl;
   std::cout << std::endl;
#endif

   // Weld points that are the same so that triangle doesn't see near-duplicates.
   // Points 0,1,2 are the triangle's primal vertices; points 3.. are the new
   // (intersection) points.

   // Union-find (with the lowest index as the representative, so primal
   // vertices win) makes chains of nearly-coincident points all resolve to one
   // point.
   UnionFind uf(numPoints);

   for (size_t ii=0; ii<new_vert_indices.size(); ++ii)
   {
      const size_t v_idx = ii+3;
      // Test against primal points
      const Eigen::Vector3d pt = ipointPos(m_newPoints[new_vert_indices[ii]].ref);
      if (point_almost_equal(p0, pt))
         uf.unite(v_idx, 0);
      if (point_almost_equal(p1, pt))
         uf.unite(v_idx, 1);
      if (point_almost_equal(p2, pt))
         uf.unite(v_idx, 2);

      for (size_t jj=0; jj<ii; ++jj)
      {
         if (point_almost_equal(ipointPos(m_newPoints[new_vert_indices[jj]].ref), pt))
            uf.unite(v_idx, jj+3);
      }
   }

   // triangle is only given the representatives
   std::vector<int> compact_index(numPoints, -1);
   std::vector<uint32_t> compact_to_local;
   for (uint32_t ii = 0; ii < numPoints; ++ii)
   {
      if (uf.find(ii) == ii)
      {
         compact_index[ii] = compact_to_local.size();
         compact_to_local.push_back(ii);
      }
   }
   const uint32_t numCompact = compact_to_local.size();

   struct triangulateio in;
   initializeTriangulateio(in);

   in.numberofpoints = numCompact;
   in.pointlist = (double*)malloc(numCompact * 2 * sizeof(double));
   for (uint32_t ii = 0; ii < numCompact; ++ii)
   {
      in.pointlist[2 * ii] = pts_2d[compact_to_local[ii]][0];
      in.pointlist[2 * ii + 1] = pts_2d[compact_to_local[ii]][1];
   }

   // The constraint segments are the intersection lines along which the face is
   // cut. Each is marked with (2 + the index of the face on the other mesh that
   // made the cut), so that the sub-faces on either side of it can be told
   // apart later. (Triangle uses markers 0 and 1 for 'not a segment' and 'edge
   // of the triangle').
   std::vector<int> segments;
   std::vector<int> segment_markers;
   for (uint32_t ii = 0; ii < numSegments; ++ii)
   {
      const uint32_t a = uf.find(3 + 2 * ii);
      const uint32_t b = uf.find(4 + 2 * ii);
      if (a == b)
         continue;  // zero length: the faces only touch at a point

      const IPoint& ip = m_newPoints[new_vert_indices[2 * ii]];
      segments.push_back(compact_index[a]);
      segments.push_back(compact_index[b]);
      segment_markers.push_back(2 + ((which_mesh == kClay) ? ip.kidx : ip.cidx));
   }

   in.numberofsegments = segment_markers.size();
   in.segmentlist = (int*)malloc((segments.size() + 1) * sizeof(int));
   in.segmentmarkerlist = (int*)malloc((segment_markers.size() + 1) * sizeof(int));
   std::copy(segments.begin(), segments.end(), in.segmentlist);
   std::copy(segment_markers.begin(), segment_markers.end(), in.segmentmarkerlist);

   struct triangulateio out;
   initializeTriangulateio(out);

   // p: use the input segments
   // c: keep the convex hull (the original triangle)
   // e: output edges, with the markers of the segments they lie on
   // z: zero-indexed
   char flags[] = "pQcze";
   triangulate(flags, &in, &out, 0);

   // Resolve the output points to references. The first numCompact are our
   // input points; any others were created by 'triangle' where two segments
   // cross, and are new points.
   std::vector<IPointRef> out_refs(out.numberofpoints);
   for (uint32_t ii = 0; ii < static_cast<uint32_t>(out.numberofpoints); ++ii)
   {
      if (ii < numCompact)
      {
         const uint32_t local = compact_to_local[ii];
         if (local < 3)
            out_refs[ii] = IPointRef(which_mesh, mesh.faces()[fidx].m_v[local]);
         else
            out_refs[ii] = m_newPoints[new_vert_indices[local - 3]].ref;
      }
      else
      {
         const Eigen::Vector3d flat(out.pointlist[2 * ii], out.pointlist[2 * ii + 1],
                                    pts_2d[0][2]);
         m_newPointPositions.push_back(rot.transpose() * flat);
         out_refs[ii] = IPointRef(kNew, m_newPointPositions.size() - 1);
      }
   }

   // Which output edges lie along a cut
   std::unordered_map<uint64_t, int32_t> edge_cutter;
   for (int ii = 0; ii < out.numberofedges; ++ii)
   {
      if (out.edgemarkerlist[ii] >= 2)
         edge_cutter[edgeKey(out.edgelist[2 * ii], out.edgelist[2 * ii + 1])] =
            out.edgemarkerlist[ii] - 2;
   }

   std::vector<IFace> new_faces;
   for (int ii = 0; ii < out.numberoftriangles; ii++)
   {
      IFace new_face;
      new_face.orig = fidx;
      for (int jj = 0; jj < 3; jj++)
      {
         const int a = out.trianglelist[ii * 3 + jj];
         const int b = out.trianglelist[ii * 3 + (jj + 1) % 3];
         new_face.v[jj] = out_refs[a];

         auto itr = edge_cutter.find(edgeKey(a, b));
         if (itr != edge_cutter.end())
            new_face.cutter[jj] = itr->second;
      }

      if (!degenerateTriangle(new_face))
      {
         new_faces.push_back(new_face);
      }
#ifdef DEBUG
      else
      {
         std::cout << "Skipping face: " << new_face.v[0].idx << "," << new_face.v[1].idx
                   << "," << new_face.v[2].idx << std::endl;
      }
#endif
   }

   free(in.pointlist);
   free(in.segmentlist);
   free(in.segmentmarkerlist);
   if (out.pointlist)         trifree((int*)out.pointlist);
   if (out.pointmarkerlist)   trifree(out.pointmarkerlist);
   if (out.trianglelist)      trifree(out.trianglelist);
   if (out.segmentlist)       trifree(out.segmentlist);
   if (out.segmentmarkerlist) trifree(out.segmentmarkerlist);
   if (out.edgelist)          trifree(out.edgelist);
   if (out.edgemarkerlist)    trifree(out.edgemarkerlist);

   return new_faces;
}


// Finds the closest point to 'p' on the triangle (a, b, c)
// (From Ericson, "Real-Time Collision Detection")
static Eigen::Vector3d closestPointOnTriangle(const Eigen::Vector3d& p, const Eigen::Vector3d& a,
                                              const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
   const Eigen::Vector3d ab = b - a;
   const Eigen::Vector3d ac = c - a;
   const Eigen::Vector3d ap = p - a;
   const double d1 = ab.dot(ap);
   const double d2 = ac.dot(ap);
   if (d1 <= 0.0 && d2 <= 0.0)
      return a;  // vertex a

   const Eigen::Vector3d bp = p - b;
   const double d3 = ab.dot(bp);
   const double d4 = ac.dot(bp);
   if (d3 >= 0.0 && d4 <= d3)
      return b;  // vertex b

   const double vc = d1 * d4 - d3 * d2;
   if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
      return a + (d1 / (d1 - d3)) * ab;  // edge ab

   const Eigen::Vector3d cp = p - c;
   const double d5 = ab.dot(cp);
   const double d6 = ac.dot(cp);
   if (d6 >= 0.0 && d5 <= d6)
      return c;  // vertex c

   const double vb = d5 * d2 - d1 * d6;
   if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
      return a + (d2 / (d2 - d6)) * ac;  // edge ac

   const double va = d3 * d6 - d5 * d4;
   if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
      return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b);  // edge bc

   const double denom = 1.0 / (va + vb + vc);
   return a + ab * (vb * denom) + ac * (vc * denom);  // inside the face
}


void CSGEngine::computeWeldTolerance()
{
   double scale = 0.0;
   for (const TriMesh* mesh : {&m_clay, &m_knife})
   {
      for (const auto& v : mesh->vertices())
      {
         const double largest = v.cwiseAbs().maxCoeff();
         if (std::isfinite(largest))
            scale = std::max(scale, largest);
      }
   }
   m_weldTolerance = WELD_ULP * std::numeric_limits<double>::epsilon() * scale;
}


void CSGEngine::buildVertexMap()
{
   const uint32_t num_clay_points = m_clay.vertices().size();
   const uint32_t num_knife_points = m_knife.vertices().size();
   const uint32_t num_new_points = m_newPointPositions.size();
   const uint32_t num_total_points = num_clay_points + num_knife_points + num_new_points;
   const uint32_t first_new_point = num_clay_points + num_knife_points;

   // All the points, numbered as in canonicalVertexIndex()
   std::vector<VECTOR3D> points;
   points.reserve(num_total_points);
   points.insert(points.end(), m_clay.vertices().begin(), m_clay.vertices().end());
   points.insert(points.end(), m_knife.vertices().begin(), m_knife.vertices().end());
   points.insert(points.end(), m_newPointPositions.begin(), m_newPointPositions.end());

   // We use union-find (rather than mapping each new point to a single arbitrary
   // "first" match) so that clusters of 3+ mutually-close points all resolve
   // to one canonical vertex instead of splitting across two or more
   // nearly-coincident vertices. Union-by-smallest-index means original mesh
   // vertices (which have the lowest indices) are preferred as the representative.
   UnionFind uf(num_total_points);

   // Every new point is welded to whatever is at the same place.  Clay and knife
   // vertices at the same place (e.g. a clay vertex lying on a knife vertex, where the
   // point is stored only once as a clay vertex) are the same vertex too.
   std::vector<uint32_t> queries;
   for (uint32_t ii = first_new_point; ii < num_total_points; ++ii)
      queries.push_back(ii);
   for (uint32_t ii = 0; ii < num_clay_points; ++ii)
      queries.push_back(ii);

   auto wanted = [=](uint32_t from, uint32_t to) {
      return (from >= first_new_point) || ((to >= num_clay_points) && (to < first_new_point));
   };
   uniteClosePoints(points, m_weldTolerance, queries, wanted, uf);

   m_vertexMap.resize(num_total_points);
   for (uint32_t ii = 0; ii < num_total_points; ++ii)
      m_vertexMap[ii] = uf.find(ii);
}


int CSGEngine::sideOfSurface(const Eigen::Vector3d& point, const TriMesh& surface) const
{
   constexpr double ulp = std::numeric_limits<double>::epsilon() * POINT_ULP;

   double best_dist2 = std::numeric_limits<double>::infinity();
   double best_alignment = -1.0;
   int best_side = 0;

   for (uint32_t ff=0; ff<surface.faces().size(); ++ff)
   {
      const Eigen::Vector3d normal = surface.faceNormal(ff);
      if (normal.squaredNorm() == 0.0)
         continue;  // a degenerate face has no sides

      const Triangle& face = surface.faces()[ff];
      const Eigen::Vector3d closest = closestPointOnTriangle(
         point, surface.vertices()[face.m_v[0]], surface.vertices()[face.m_v[1]],
         surface.vertices()[face.m_v[2]]);

      const Eigen::Vector3d offset = point - closest;
      const double dist2 = offset.squaredNorm();
      const double signed_dist = offset.dot(normal);
      const double alignment = (dist2 > 0.0) ? std::fabs(signed_dist) / std::sqrt(dist2) : 0.0;

      // The nearest triangle decides.  Where the closest point is on an edge or
      // vertex several triangles tie, and we take the one the point is most
      // directly in front of (or behind).
      const bool nearer = dist2 < best_dist2 * (1.0 - ulp);
      const bool tied = !nearer && (dist2 <= best_dist2 * (1.0 + ulp));
      if (nearer || (tied && alignment > best_alignment))
      {
         best_dist2 = dist2;
         best_alignment = alignment;
         best_side = (signed_dist > 0.0) ? 1 : ((signed_dist < 0.0) ? -1 : 0);
      }
   }

   return best_side;
}


void CSGEngine::classifyFaces(IParent which_surface, const std::vector<IFace>& new_faces,
                              const std::vector<bool>& is_face_cut,
                              std::vector<char>& cut_face_status,
                              std::vector<char>& uncut_face_status)
{
   // Classify every face of a surface as above (1) or below (-1) the other
   // surface, w.r.t. the normals of the other surface's faces.
   //
   // 1. Group the faces (the new faces that replaced the cut faces, and the
   //    original faces that were not cut) into regions by flood-filling across
   //    shared edges.  The fill is not allowed to cross the edges along which
   //    the faces were cut: those are the borders between above and below.
   // 2. Every new face that lies along a cut votes for its region's side, by
   //    which side of the cutting face's plane it is on.
   // 3. A region that no cut borders, or with no clear majority, is tested
   //    against the other surface directly.
   const TriMesh& original_mesh = getMesh(which_surface);
   const TriMesh& cutter_mesh = otherMesh(which_surface);

   const size_t num_new = new_faces.size();
   const size_t num_orig = original_mesh.faces().size();
   const size_t num_nodes = num_new + num_orig;  // new faces first, then original faces

   auto vertexId = [this](const IPointRef& ref) { return m_vertexMap[canonicalVertexIndex(ref)]; };

   // Map from (welded) edge to the faces on it; and the set of edges that are cuts
   std::unordered_map<uint64_t, std::vector<uint32_t>> edge_faces;
   std::unordered_set<uint64_t> cut_edges;

   for (uint32_t ff=0; ff<num_new; ++ff)
   {
      for (uint32_t kk=0; kk<3; ++kk)
      {
         const uint32_t a = vertexId(new_faces[ff].v[kk]);
         const uint32_t b = vertexId(new_faces[ff].v[(kk+1)%3]);
         if (a == b)
            continue;

         edge_faces[edgeKey(a, b)].push_back(ff);
         if (new_faces[ff].cutter[kk] >= 0)
            cut_edges.insert(edgeKey(a, b));
      }
   }

   // Original faces that were replaced by new faces are not nodes
   for (uint32_t ff=0; ff<num_orig; ++ff)
   {
      if (is_face_cut[ff])
         continue;

      const Triangle& face = original_mesh.faces()[ff];
      for (uint32_t kk=0; kk<3; ++kk)
      {
         const uint32_t a = vertexId(IPointRef(which_surface, face.m_v[kk]));
         const uint32_t b = vertexId(IPointRef(which_surface, face.m_v[(kk+1)%3]));
         if (a == b)
            continue;

         edge_faces[edgeKey(a, b)].push_back(num_new + ff);
      }
   }

   UnionFind regions(num_nodes);
   for (const auto& edge : edge_faces)
   {
      if (cut_edges.count(edge.first))
         continue;  // the border between regions

      for (size_t ii=1; ii<edge.second.size(); ++ii)
         regions.unite(edge.second[0], edge.second[ii]);
   }

   // Faces along a cut vote.  The face is above or below the cutting face's plane
   // according to its vertex opposite the cut edge (that vertex can't be on the
   // plane, or the face would be degenerate).
   constexpr double ulp = std::numeric_limits<double>::epsilon() * POINT_ULP;
   std::vector<int> votes(num_nodes, 0);
   for (uint32_t ff=0; ff<num_new; ++ff)
   {
      for (uint32_t kk=0; kk<3; ++kk)
      {
         if (new_faces[ff].cutter[kk] < 0)
            continue;

         const uint32_t cutter = new_faces[ff].cutter[kk];
         const Eigen::Vector3d normal = cutter_mesh.faceNormal(cutter);
         const Eigen::Vector3d origin = cutter_mesh.vertices()[cutter_mesh.faces()[cutter].m_v[0]];
         const Eigen::Vector3d apex = ipointPos(new_faces[ff].v[(kk+2)%3]);

         const double height = (apex - origin).dot(normal);
         const double scale = std::max(apex.cwiseAbs().maxCoeff(), origin.cwiseAbs().maxCoeff());
         if (!(std::fabs(height) > ulp * scale))
            continue;

         votes[regions.find(ff)] += (height > 0.0) ? 1 : -1;
      }
   }

   // Position of a point on a face (node), for testing it against the other surface
   auto facePoint = [&](uint32_t node) -> Eigen::Vector3d
   {
      if (node < num_new)
      {
         return (ipointPos(new_faces[node].v[0]) + ipointPos(new_faces[node].v[1]) +
                 ipointPos(new_faces[node].v[2])) / 3.0;
      }
      return centroid(original_mesh, node - num_new);
   };

   // The regions with no (or no clear) vote are classified by testing their faces
   // against the other surface: the first face that is clearly on one side decides.
   // (Not just any face: one might lie right in the other surface.)
   std::unordered_map<uint32_t, std::vector<uint32_t>> unvoted_regions;
   for (uint32_t node=0; node<num_nodes; ++node)
   {
      const bool active = (node < num_new) || !is_face_cut[node - num_new];
      if (active && (votes[regions.find(node)] == 0))
         unvoted_regions[regions.find(node)].push_back(node);
   }

   std::vector<char> region_status(num_nodes, 0);
   for (const auto& region : unvoted_regions)
   {
      for (uint32_t node : region.second)
      {
         const int side = sideOfSurface(facePoint(node), cutter_mesh);
         if (side != 0)
         {
            region_status[region.first] = side;
            break;
         }
      }
   }

   auto statusOf = [&](uint32_t node) -> char
   {
      const uint32_t root = regions.find(node);
      if (votes[root] != 0)
         return (votes[root] > 0) ? 1 : -1;
      return region_status[root];
   };

   cut_face_status.assign(num_new, 0);
   for (uint32_t ff=0; ff<num_new; ++ff)
      cut_face_status[ff] = statusOf(ff);

   uncut_face_status.assign(num_orig, 0);
   for (uint32_t ff=0; ff<num_orig; ++ff)
   {
      if (!is_face_cut[ff])
         uncut_face_status[ff] = statusOf(num_new + ff);
   }

#ifdef DEBUG
   std::cout << "classification results, uncut (original) faces: " << std::endl;
   for (size_t ii=0; ii<uncut_face_status.size(); ++ii)
   {
      if (is_face_cut[ii])
         continue;
      std::cout << " " << ii << " : " << int(uncut_face_status[ii]) << std::endl;
   }
#endif
}


TriMesh CSGEngine::assembleMesh(IParent which_surface, char side,
                                const std::vector<IFace>& new_faces,
                                const std::vector<char>& cut_face_status,
                                const std::vector<char>& uncut_face_status)
{
   const TriMesh& original_mesh = getMesh(which_surface);
   const size_t num_total_points = m_vertexMap.size();
   const std::vector<uint32_t>& vertex_map = m_vertexMap;

   std::vector<VECTOR3D> vertices;
   std::vector<VECTOR3I> faces;

   std::vector<uint32_t> used_map(num_total_points, UINT32_MAX);

   //  original faces that were not cut
   uint32_t count = 0;
   for (uint32_t ii=0; ii<original_mesh.faces().size(); ++ii)
   {
      if (uncut_face_status[ii] == side)
      {
         Eigen::Vector3i face;
         for (uint32_t jj=0; jj<3; ++jj)
         {
            // Key by the welded canonical index (the same key the cut faces
            // below use), so vertices shared between cut and uncut faces
            // are emitted once, and knife-local indices can't alias clay ones.
            const uint32_t local_idx = original_mesh.faces()[ii].m_v[jj];
            const uint32_t idx =
               vertex_map[canonicalVertexIndex(IPointRef(which_surface, local_idx))];
            if (used_map[idx] == UINT32_MAX)
            {
               vertices.push_back( original_mesh.vertices()[local_idx] );
               face[jj] = count;
               used_map[idx] = count;
               count++;
            }
            else
            {
               face[jj] = used_map[idx];
            }
         }
         faces.push_back(face);
      }
   }

   // new faces, resulting from cuts
   for (uint32_t ii=0; ii<new_faces.size(); ++ii)
   {
      if (cut_face_status[ii] == side)
      {
         Eigen::Vector3i face;
         for (uint32_t jj=0; jj<3; ++jj)
         {
            const auto& pt = new_faces[ii].v[jj];
            uint32_t idx = vertex_map[canonicalVertexIndex(pt)];
            if (used_map[idx] == UINT32_MAX)
            {
               vertices.push_back(ipointPos(pt));
               face[jj] = count;
               used_map[idx] = count;
               count++;
            }
            else
            {
               face[jj] = used_map[idx];
            }
         }
         faces.push_back(face);
      }
   }

   std::cout << "Final assembled mesh: " << vertices.size() << " vertices" << std::endl;
   std::cout << "                      " << faces.size() << " triangles" << std::endl;
   return TriMesh(vertices, faces);
}


TriMesh CSGEngine::mergeMeshes(const TriMesh& exterior, const TriMesh& cap)
{
   //
   // Weld the vertices that are at the same place.  'vertex_map' maps from the
   // unwelded vertex numbering to the welded numbering. Union-find guarantees
   // clusters of 3+ mutually-close vertices all resolve to a single canonical vertex;
   // union-by-smallest-index means exterior vertices (lower indices) are preferred as
   // the representative.
   //
   const uint32_t num_ext_vertices = exterior.vertices().size();
   const uint32_t num_cap_vertices = cap.vertices().size();
   const uint32_t num_total_vertices = num_ext_vertices + num_cap_vertices;

   std::vector<VECTOR3D> points;
   points.reserve(num_total_vertices);
   points.insert(points.end(), exterior.vertices().begin(), exterior.vertices().end());
   points.insert(points.end(), cap.vertices().begin(), cap.vertices().end());

   std::vector<uint32_t> queries(num_total_vertices);
   for (uint32_t ii = 0; ii < num_total_vertices; ++ii)
      queries[ii] = ii;

   UnionFind uf(num_total_vertices);
   uniteClosePoints(points, m_weldTolerance, queries,
                    [](uint32_t, uint32_t) { return true; }, uf);

   std::vector<uint32_t> vertex_map(num_total_vertices);
   for (uint32_t ii=0; ii<num_total_vertices; ++ii)
      vertex_map[ii] = uf.find(ii);

   // Construct the merged mesh
   TriMesh output;
   std::vector<VECTOR3D>& vertices = output.vertices();
   std::vector<VECTOR3I> faces;
   std::vector<uint32_t> used_map(num_total_vertices, UINT32_MAX);

   uint32_t count = 0;
   const size_t num_exterior_faces = exterior.faces().size();
   for (uint32_t ii=0; ii<num_exterior_faces; ++ii)
   {
      Eigen::Vector3i face;
      for (uint32_t jj=0; jj<3; ++jj)
      {
         const uint32_t local_idx = exterior.faces()[ii].m_v[jj];
         const uint32_t idx = vertex_map[local_idx];
         if (used_map[idx] == UINT32_MAX)
         {
            vertices.push_back( exterior.vertices()[local_idx] );
            face[jj] = count;
            used_map[idx] = count;
            count++;
         }
         else
         {
            face[jj] = used_map[idx];
         }
      }
      faces.push_back(face);
   }
   const size_t num_cap_faces = cap.faces().size();
   for (uint32_t ii=0; ii<num_cap_faces; ++ii)
   {
      Eigen::Vector3i face;
      for (uint32_t jj=0; jj<3; ++jj)
      {
         const uint32_t cfvi = cap.faces()[ii].m_v[jj];
         const uint32_t idx = vertex_map[cfvi + num_ext_vertices];
         if (used_map[idx] == UINT32_MAX)
         {
            vertices.push_back( cap.vertices()[cfvi] );
            face[jj] = count;
            used_map[idx] = count;
            count++;
         }
         else
         {
            face[jj] = used_map[idx];
         }
      }
      faces.push_back(face);
   }
   output.setFaces(faces);

   return output;
}


// Returns (number of edges where the cap's winding is consistent with the
// exterior's) - (number where it is inconsistent), among the edges they share
// in a merged (welded) mesh whose first 'num_exterior_faces' faces are the
// exterior.  A positive result means the cap is oriented correctly; negative
// means it should be flipped; zero means there is no evidence either way.
static int capWindingAgreement(const TriMesh& merged, size_t num_exterior_faces)
{
   std::set<std::pair<int32_t, int32_t>> exterior_edges;  // directed
   const auto& faces = merged.faces();
   for (size_t ff = 0; ff < num_exterior_faces; ++ff)
      for (int ii = 0; ii < 3; ++ii)
         exterior_edges.insert({faces[ff].m_v[ii], faces[ff].m_v[(ii + 1) % 3]});

   int agreement = 0;
   for (size_t ff = num_exterior_faces; ff < faces.size(); ++ff)
   {
      for (int ii = 0; ii < 3; ++ii)
      {
         const int32_t a = faces[ff].m_v[ii];
         const int32_t b = faces[ff].m_v[(ii + 1) % 3];
         if (exterior_edges.count({b, a}))
            ++agreement;  // opposite direction: consistent
         if (exterior_edges.count({a, b}))
            --agreement;  // same direction: inconsistent
      }
   }
   return agreement;
}


TriMesh CSGEngine::mergeWithCap(const TriMesh& exterior, const TriMesh& cap,
                                const TriMesh& cap_flipped)
{
   TriMesh merged = mergeMeshes(exterior, cap);
   if (capWindingAgreement(merged, exterior.faces().size()) < 0)
      return mergeMeshes(exterior, cap_flipped);
   return merged;
}


void CSGEngine::construct(CSGOperation operation, bool cap, TriMesh& out_A, TriMesh& out_B)
{
   computeWeldTolerance();

   // Create and intersect AABB trees
   std::unordered_set<std::pair<uint32_t, uint32_t>> ix = intersectAABBs(m_clay, m_knife);

   // Calculate actual triangle intersections, with intersection data

   // Sets of indices of faces that were cut
   std::unordered_set<uint32_t> cut_clay_faces;
   std::unordered_set<uint32_t> cut_knife_faces;
   // Vectors of bools, indicating whether the face was cut; faster than seeking
   // within an unordered_set
   std::vector<bool> is_clay_face_cut(m_clay.faces().size(), false);
   std::vector<bool> is_knife_face_cut(m_knife.faces().size(), false);


   // We track vertices on faces created by the intersection process in
   //  m_newPoints, which contains IPoints which refer to either clay or knife
   //  mesh vertices, or m_newPointPositions

   // Indices into the m_newPoints vector
   std::unordered_map<uint32_t, std::vector<uint32_t>> clay_face_new_verts;
   std::unordered_map<uint32_t, std::vector<uint32_t>> knife_face_new_verts;

   for (auto ix_itr=ix.begin(); ix_itr != ix.end(); ++ix_itr)
   {
      const uint32_t c_idx = ix_itr->first;
      const uint32_t k_idx = ix_itr->second;
      const Triangle& ct = m_clay.faces()[c_idx];
      const Triangle& kt = m_knife.faces()[k_idx];

      TriangleIntersection trix = intersect(ct, m_clay.vertices(), kt, m_knife.vertices());
      if (!trix.intersect) // just because AABBs intersect doesn't mean the faces do
         continue;

      if (trix.coplanar)
         continue;  // TODO: will handle this later

      // Record which faces have been cut, so we can replace them later
      // with the diced up versions
      cut_clay_faces.insert(c_idx);
      cut_knife_faces.insert(k_idx);
      is_clay_face_cut[c_idx] = true;
      is_knife_face_cut[k_idx] = true;

      // Convert the raw positions to indices into one of three arrays: existing clay mesh
      // vertices, existing knife mesh vertices, or new vertices (which may be duplicates
      // of other new vertices; we'll merge and clean these up later)
      std::vector<IPoint> pts = convertIntersectionToIpoints(trix, c_idx, k_idx);
      m_newPoints.insert(m_newPoints.end(), pts.begin(), pts.end());

      // Record which elements of 'm_newPoints' correspond to cuts in each face
      auto itr = clay_face_new_verts.find(c_idx);
      if (itr == clay_face_new_verts.end())
         clay_face_new_verts.insert({c_idx, std::vector<uint32_t>()});
      clay_face_new_verts[c_idx].push_back(m_newPoints.size() - 1);
      clay_face_new_verts[c_idx].push_back(m_newPoints.size() - 2);

      itr = knife_face_new_verts.find(k_idx);
      if (itr == knife_face_new_verts.end())
          knife_face_new_verts.insert({k_idx, std::vector<uint32_t>()});
      knife_face_new_verts[k_idx].push_back(m_newPoints.size() - 1);
      knife_face_new_verts[k_idx].push_back(m_newPoints.size() - 2);
   }

   // Step through the cut faces and retriangulate them. New triangles are stored
   // as triplets of point references
   std::vector<IFace> new_clay_faces;
   for (auto itr = cut_clay_faces.begin(); itr != cut_clay_faces.end(); ++itr)
   {
      auto result = retriangulate(m_clay, kClay, *itr, clay_face_new_verts[*itr]);
      new_clay_faces.insert(new_clay_faces.end(), result.begin(), result.end());
   }

   std::vector<IFace> new_knife_faces;
   if (cap)
   {
       for (auto itr = cut_knife_faces.begin(); itr != cut_knife_faces.end(); ++itr)
       {
           auto result = retriangulate(m_knife, kKnife, *itr, knife_face_new_verts[*itr]);
           new_knife_faces.insert(new_knife_faces.end(), result.begin(), result.end());
       }
   }

   // All points now exist; weld the coincident ones
   buildVertexMap();

   // Classify all faces as either above the triangle that intersected them (w.r.t. the
   // intersecting triangle's face normal) or below.  Faces are grouped into regions
   // that are connected across shared edges, but not across the edges where the
   // cuts were made, and each region is classified as a whole (see classifyFaces).
   //
   std::vector<char> clay_cut_face_status;
   std::vector<char> clay_uncut_face_status;
   classifyFaces(kClay, new_clay_faces, is_clay_face_cut, clay_cut_face_status,
                 clay_uncut_face_status);

   // We only need to partition the knife mesh's faces if we're going to put
   // caps on the two halves of the cut clay mesh
   std::vector<char> knife_cut_face_status;
   std::vector<char> knife_uncut_face_status;
   if (cap)
   {
      classifyFaces(kKnife, new_knife_faces, is_knife_face_cut, knife_cut_face_status,
                    knife_uncut_face_status);
   }

   constexpr char ABOVE = 1;
   constexpr char BELOW = -1;
   TriMesh clay_above = assembleMesh(kClay, ABOVE, new_clay_faces, clay_cut_face_status,
                                     clay_uncut_face_status);
   TriMesh clay_below = assembleMesh(kClay, BELOW, new_clay_faces, clay_cut_face_status,
                                     clay_uncut_face_status);

   if (cap)
   {
       // We only care about the knife faces _below_ the clay mesh, because those
       // are the faces _inside_ the clay.
       TriMesh knife_below = assembleMesh(kKnife, BELOW, new_knife_faces, knife_cut_face_status,
                                          knife_uncut_face_status);
       TriMesh knife_below_flipped = knife_below;
       for (auto& face : knife_below_flipped.faces())
           std::swap(face.m_v[1], face.m_v[2]);

       // Which orientation of the cap closes each half depends on the knife's
       // normal relative to the clay, so pick whichever winding is consistent
       // with the exterior along their shared (cut) edges.
       out_A = mergeWithCap(clay_above, knife_below, knife_below_flipped);
       out_B = mergeWithCap(clay_below, knife_below, knife_below_flipped);
   }
   else
   {
       out_A = clay_above;
       out_B = clay_below;
   }
}

}  // namespace CSG
