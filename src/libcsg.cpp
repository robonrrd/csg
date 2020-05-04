//
//
//
//
#include <assert.h>
#include <iomanip>      // std::setprecision
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


//
// Helper functions
//
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
         pts[jj].ref.idx = clay_v_idx;
      }
      else if ((clay_v_idx == 4) && (knife_v_idx < 4))
      {
         // This point is a knife mesh vertex
         pts[jj].ref.parent = kKnife;
         pts[jj].ref.idx = knife_v_idx;
      }
      else if ((clay_v_idx < 4) && (knife_v_idx < 4))
      {
         // This point is a vertex on both clay and knife
         pts[jj].ref.parent = kBoth;
         // TODO: should  we preserve both clay and knife vert indices here?
         pts[jj].ref.idx = clay_v_idx;
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


// For debugging
double triangleAspectRatio(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1,
                           const Eigen::Vector3d& p2)
{
   // Ratio between longest and shortest leg
   const double d01 = (p0-p1).norm();
   const double d12 = (p1-p2).norm();
   const double d20 = (p2-p0).norm();

   const double mmax = std::max(std::max(d01, d12), d20);
   const double mmin = std::min(std::min(d01, d12), d20);

   return mmax/mmin;
}

double triangleAltitude(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1,
                           const Eigen::Vector3d& p2)
{
   const double a = (p0-p1).norm();
   const double b = (p1-p2).norm();
   const double c = (p2-p0).norm();
   const double s = (a+b+c)*0.5;

   const double h_a = 2*sqrt(s*(s-a)*(s-b)*(s-c))/a;
   const double h_b = 2*sqrt(s*(s-a)*(s-b)*(s-c))/b;
   const double h_c = 2*sqrt(s*(s-a)*(s-b)*(s-c))/c;

   return std::min(std::min(h_a, h_b), h_c);
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
   if (ref.parent == kClay)
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

   // Check for bad aspect ratio
   const double BAD_VALUE = 1e+10; // TODO
   const double aspect_ratio = triangleAspectRatio(ipointPos(face.v[0]), ipointPos(face.v[1]),
                                                   ipointPos(face.v[2]));
   if (aspect_ratio > BAD_VALUE)
      return true;

   // Check for bad altitude (a squashed triangle)
   const double BAD_ALTITUDE = 1e-6; // TODO
   const double min_altitude = triangleAltitude(ipointPos(face.v[0]), ipointPos(face.v[1]),
                                                 ipointPos(face.v[2]));

   if (min_altitude < BAD_ALTITUDE)
      return true;

   return false;
}

std::vector<IFace> CSGEngine::retriangulate(const TriMesh& mesh, IParent which_mesh, uint32_t fidx,
                                            const std::vector<uint32_t>& new_vert_indices) const
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

   const Eigen::Vector3d n = (p1 - p0).normalized().cross((p2 - p0).normalized());  // normal
   const Eigen::Vector3d axis = n.cross(Eigen::Vector3d::UnitZ()).normalized();
   const double angle = acos(n[2]);
   const Eigen::Matrix3d rot = Eigen::AngleAxisd(angle, axis).matrix();

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

   struct triangulateio in;
   initializeTriangulateio(in);

   in.numberofpoints = numPoints;
   in.pointlist = (double*)malloc(numPoints * 2 * sizeof(double));

   for (uint32_t ii = 0; ii < numPoints; ++ii)
   {
      in.pointlist[2 * ii] = pts_2d[ii][0];
      in.pointlist[2 * ii + 1] = pts_2d[ii][1];
   }

   in.numberofsegments = numSegments;
   in.segmentlist = (int*)malloc(numSegments * 2 * sizeof(int));
   for (uint32_t ii = 0; ii < numSegments; ++ii)
   {
      in.segmentlist[2 * ii] = 3 + ii;
      in.segmentlist[2 * ii + 1] = 4 + ii;
   }

   struct triangulateio out;
   initializeTriangulateio(out);

   char flags[] = "Qcz";  // pcze
   triangulate(flags, &in, &out, 0);

   // triangulate is extremely precise and may create vertices which are the
   // same for any reasonable geometric purpose. Here, we weld all vertices that
   // are closer than some epsilon. We do this in a brute-force way because
   // we'll have very few vertices in this function
   const size_t num_initial_verts = 3 + new_vert_indices.size();
   std::vector<uint32_t> vertex_mapping(num_initial_verts);
   for (size_t ii=0; ii<num_initial_verts; ++ii)
   {
      vertex_mapping[ii] = ii;
   }

   for (size_t ii=0; ii<new_vert_indices.size(); ++ii)
   {
      const size_t v_idx = ii+3;
      // Test against primal points
      const Eigen::Vector3d pt = ipointPos(m_newPoints[new_vert_indices[ii]].ref);
      if (point_almost_equal(p0, pt))
      {
         vertex_mapping[v_idx] = 0;
         continue;
      }
      if (point_almost_equal(p1, pt))
      {
         vertex_mapping[v_idx] = 1;
         continue;
      }
      if (point_almost_equal(p2, pt))
      {
         vertex_mapping[v_idx] = 2;
         continue;
      }
      for (size_t jj=0; jj<new_vert_indices.size(); ++jj)
      {
         if (jj >= ii)
            continue;
         if (point_almost_equal(ipointPos(m_newPoints[new_vert_indices[jj]].ref), pt))
         {
            vertex_mapping[v_idx] = jj+3;
            continue;
         }
      }
   }

#ifdef DEBUG
   std::cout << "Vertex mapping:"  << std::endl;
   for (size_t ii=0; ii<num_initial_verts; ++ii)
   {
      std::cout << "  " << ii << " -> " << vertex_mapping[ii] << std::endl;
   }
#endif

   // Convert the local triangle vertex indices into global IPointRefs
   std::vector<IFace> new_faces;
   for (uint32_t ii = 0; ii < out.numberoftriangles; ii++)
   {
      IFace new_face;
      new_face.orig = fidx;
      for (uint32_t jj = 0; jj < 3; jj++)
      {
         uint32_t idx = vertex_mapping[out.trianglelist[ii * 3 + jj]];
         if (idx < 3) // original vertices
         {
            new_face.v[jj].parent = which_mesh;
            new_face.v[jj].idx = mesh.faces()[fidx].m_v[idx];
         }
         else
         {
            new_face.v[jj] = m_newPoints[new_vert_indices[idx-3]].ref;
         }
      }

      if (!degenerateTriangle(new_face))
      {
         new_faces.push_back(new_face);
      }
#if DEBUG
      else
      {
         std::cout << "Skipping face: " << new_face.v[0].idx << "," << new_face.v[1].idx
                   << "," << new_face.v[2].idx << std::endl;
      }
#endif
   }

   free(in.pointlist);
   free(in.segmentlist);

   return new_faces;

}


std::vector<char> CSGEngine::classifyCutFaces(const std::vector<IFace>& in_faces,
                                              IParent which_surface)
{
   // status:
   // -1 - below
   //  0 - unknown
   //  1 - above
   const size_t num_faces = in_faces.size();
   std::vector<char> status(num_faces, 0);

   for (size_t fidx = 0; fidx < num_faces; ++fidx)
   {
      const IFace& ff = in_faces[fidx];

      // Find a face on the opposite surface (the one that cut this one) to
      // compare against for above/belowness
      IPointRef testFace; // repurposing this structure to track a face
      for (size_t ii=0; ii<3; ++ii)
      {
         if (ff.v[ii].parent == kNew)
         {
            if (which_surface == kClay)
            {
               testFace.parent = kKnife;
               testFace.idx = m_newPoints[ff.v[ii].idx].kidx;
            }
            else
            {
               testFace.parent = kClay;
               testFace.idx = m_newPoints[ff.v[ii].idx].cidx;
            }
         }
      }

      // TODO: precompute or lazily compute face normals
      Eigen::Vector3d normal;
      Eigen::Vector3d center;
      if (testFace.parent == kClay)
      {
         normal = m_clay.faceNormal(testFace.idx);
         center = centroid(m_clay, testFace.idx);
      }
      else
      {
         normal = m_knife.faceNormal(testFace.idx);
         center = centroid(m_knife, testFace.idx);
      }

      int32_t vote_above = 0;
      int32_t vote_below = 0;
      double height_above = 0;
      double height_below = 0;
      // Test each vertex
      Eigen::Vector3d face_center(0,0,0);
      for (size_t ii=0; ii<3; ++ii)
      {
         const Eigen::Vector3d vpos = ipointPos(ff.v[ii]);
         face_center += vpos; // accumulate the face center for use later

         if (ff.v[ii].parent != kNew) // new points are on the border, so can't be tested
         {
            const double nn = (vpos - center).dot(normal);
            if (nn > 0)
            {
               vote_above++;
               height_above += nn;
            }
            else if (nn < 0)
            {
               vote_below++;
               height_below -= nn;
            }
         }
      }

      // Test triangle center
      face_center /= 3.0;
      const double nn = (face_center - center).dot(normal);
      if (nn > 0)
      {
         vote_above++;
         height_above += nn;
      }
      else if (nn < 0)
      {
         vote_below++;
         height_below -= nn;
      }

      if (vote_above > vote_below)
      {
         status[fidx] = 1;
      }
      else if (vote_below > vote_above)
      {
         status[fidx] = -1;
      }
      else
      {
         if (height_above > height_below)
         {
            status[fidx] = 1;
         }
         else if (height_below > height_above)
         {
            status[fidx] = -1;
         }
         else
         {
            status[fidx] = 0;
         }
      }
   }

   return status;
}


void CSGEngine::classifyFaces(IParent which_surface, const std::vector<IFace>& new_faces,
                              const std::vector<bool>& is_face_cut,
                              std::vector<char>& cut_face_status,
                              std::vector<char>& uncut_face_status)
{
   const TriMesh& original_mesh = ( (which_surface == kClay) ? m_clay : m_knife );

   // Classify cut (new) faces into 'above' and 'below' sets, depending on if
   // they are above or below the faces that cut them (with respect to that
   // face's normal)
   const size_t numNewFaces = new_faces.size();
   cut_face_status = classifyCutFaces(new_faces, which_surface);

   // Flood-fill uncut faces:
   //  1. create adjacency map between uncut faces and cut faces, using
   //     shared 'original' (not new ones created by cutting) vertices
   //  2. Select a classified cut face and begin filling
   //  3. All unlabeled faces belong to the other set

   // Map from a canonical vertex index to a canonical face index
   std::unordered_map<uint32_t, std::vector<uint32_t>> adjacency_map;

   // First, add all the cut faces (these have been classified)
   for (size_t ff=0; ff<numNewFaces; ++ff)
   {
      if (cut_face_status[ff] == 0)
      {
         //std::cout << "Skipping unclassified cut face " << ff << "!" << std::endl;
         continue;
      }

      const auto& face = new_faces[ff];
      for (size_t ii=0; ii<3; ++ii)
      {
         if (face.v[ii].parent == kNew) // skip new vertices
            continue;

         uint32_t idx = canonicalVertexIndex(face.v[ii]);
         auto itr = adjacency_map.find(idx);
         if (itr == adjacency_map.end())
         {
            adjacency_map.insert({idx, std::vector<uint32_t>(1,ff)});
         }
         else
         {
            itr->second.push_back(ff);
         }
      }
   }

   // Then add the uncut primal faces (original faces of the clay won't be
   // included in any result, they have been replaced by the new faces above)
   for(size_t ff=0; ff<original_mesh.faces().size(); ++ff)
   {
      if (is_face_cut[ff])
         continue;

      const auto& face = original_mesh.faces()[ff];
      uint32_t f_idx = ff + numNewFaces; // offset indexing
      for (size_t vv=0; vv<3; ++vv)
      {
         const uint32_t idx = canonicalVertexIndex(IPointRef(which_surface, face.m_v[vv]));
         auto itr = adjacency_map.find(idx);
         if (itr == adjacency_map.end())
         {
            adjacency_map.insert({idx, std::vector<uint32_t>(1,f_idx)});
         }
         else
         {
            itr->second.push_back(f_idx);
         }
      }
   }

#if DEBUG
   {
      // Print out adjacency map
      std::cout << "Adjacency map for " << ( (which_surface == kClay) ? "clay" : "knife" )
                << std::endl;
      std::cout << " (canonical vertex index -> canonical face index)" << std::endl;
      std::cout << " original face offset = " << numNewFaces << std::endl;
      for (auto itr=adjacency_map.cbegin(); itr!=adjacency_map.cend(); ++itr)
      {
         std::cout << itr->first << " -> ";
         for (size_t ii=0; ii<itr->second.size(); ++ii)
         {
            std::cout << itr->second[ii] << " ";
         }
         std::cout << std::endl;
      }
      std::cout << std::endl;
   }
#endif

   std::stack<uint32_t> faces_to_classify;
   std::vector<bool> uncut_face_classified(original_mesh.faces().size(), false);
   uncut_face_status.clear();
   uncut_face_status.resize(original_mesh.faces().size(), 0);
   for (size_t ff=0; ff<numNewFaces; ++ff)
   {
      // Find a successfully classified face
      const char cur_status = cut_face_status[ff];
      if (cur_status != 0)
      {
         //std::cout << "Seeding with cut face " << ff << std::endl;
         for (size_t vv=0; vv<3; ++vv)
         {
            uint32_t idx = canonicalVertexIndex(new_faces[ff].v[vv]);
            auto itr = adjacency_map.find(idx);
            if (itr != adjacency_map.end())
            {
               for (auto itr2=itr->second.begin(); itr2!=itr->second.end(); ++itr2)
               {
                  if (*itr2 < numNewFaces)
                  {
                     // don't add cut (i.e. classified) faces
                     continue;
                  }
                  if (uncut_face_classified[*itr2-numNewFaces])
                  {
                     // don't add already-classified uncut faces
                     continue;
                  }
                  faces_to_classify.push(*itr2);
                  //std::cout << " adding face " << *itr2 << " to stack" << std::endl;
               }
            }
         }
         //std::cout << std::endl;

         // go through the stack until we're finished
         while (!faces_to_classify.empty())
         {
            const uint32_t f_idx = faces_to_classify.top();
            faces_to_classify.pop();

            //std::cout << "Examining face " << f_idx << std::endl;
            if (f_idx < numNewFaces)  // cut faces
            {
               if (cut_face_status[f_idx] != 0)
               {
                  // don't examine classified cut faces
                  continue;
               }
               else
               {
                  cut_face_status[f_idx] = cur_status;
                  //std::cout << "..newly classified cut face" << std::endl;
               }
            }
            else // original (uncut) faces
            {
               const size_t original_index = f_idx - numNewFaces;
               if (uncut_face_classified[original_index])
               {
                  // don't add already-classified uncut faces
                  continue;
               }

               uncut_face_status[original_index] = cur_status;
               uncut_face_classified[original_index] = true;
               //std::cout << "..newly classified uncut face" << std::endl;

               // add this face's unclassified neigbors to the stack
               const auto& face = original_mesh.faces()[original_index];
               for (size_t vv=0; vv<3; ++vv)
               {
                  IPointRef rpt(which_surface, face.m_v[vv]);
                  const uint32_t idx = canonicalVertexIndex(rpt);
                  //std::cout << "Looking for vertex " << idx << " in adj map" << std::endl;
                  auto itr = adjacency_map.find(idx);
                  if (itr == adjacency_map.end())
                     continue;

                  for (size_t kk=0; kk<itr->second.size(); ++kk)
                  {
                     const uint32_t cf = itr->second[kk];
                     //std::cout << "..examining face " << cf << " from adj map" << std::endl;
                     if (cf < numNewFaces)
                     {
                        continue; // skip cut faces
                     }
                     if (!uncut_face_classified[cf - numNewFaces])
                     {
                        //std::cout << "...adding attached face " << cf << std::endl;
                        faces_to_classify.push(cf);
                     }
                  }
               }
            }
         }
      }
   }

#ifdef DEBUG
   std::cout << "classification results, uncut (original) faces: " << std::endl;
   for (int ii=0; ii<uncut_face_status.size(); ++ii)
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
   const TriMesh& opposite_mesh = otherMesh(which_surface);

   const size_t num_original_points = original_mesh.vertices().size();
   const size_t num_opposite_points = opposite_mesh.vertices().size();
   const size_t num_new_points = m_newPoints.size();
   const size_t num_total_points = num_original_points + num_opposite_points + num_new_points;
   const double radius = 1e-06; // TODO: express in ULP

   //
   // Weld vertices using an AABB tree
   //

   // 'vertex_map' maps from the unwelded vertex numbering to the welded
   // numbering
   std::vector<uint32_t> vertex_map(num_total_points, UINT32_MAX);
   AABBTree tree(0.05, num_total_points);

   // Add original mesh vertices
   for (uint32_t ii=0; ii<num_original_points; ++ii)
   {
      tree.addSphere(ii, original_mesh.vertices()[ii], radius);
      // we declare that 'original_mesh' vertices come first and therefore are always
      // used if a duplicate point shows up
      vertex_map[ii] = ii;
   }

   // Add opposite mesh vertices
   for (uint32_t ii=0; ii<num_opposite_points; ++ii)
   {
      const uint32_t idx = ii+num_original_points;
      tree.addSphere(idx, opposite_mesh.vertices()[ii], radius);
      vertex_map[idx] = idx;
   }

   // Add new points created by intersections
   for (uint32_t ii=0; ii<num_new_points; ++ii)
   {
      const uint32_t idx = ii+num_original_points+num_opposite_points;
      tree.addSphere(idx, m_newPointPositions[ii], radius);
   }

   // Query the tree with every new point, to detect duplicates
   for (uint32_t ii=0; ii<num_new_points; ++ii)
   {
      const size_t idx = ii+num_original_points+num_opposite_points;
      auto dupes = tree.query(idx);
      if (dupes.size() == 0)
      {
         vertex_map[idx] = idx;
      }
      else
      {
         vertex_map[idx] = dupes[0]; // use the first
      }
   }

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
            uint32_t idx = original_mesh.faces()[ii].m_v[jj];
            if (used_map[idx] == UINT32_MAX)
            {
               vertices.push_back( original_mesh.vertices()[idx] );
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
   // Weld vertices using an AABB tree
   //
   const size_t num_ext_vertices = exterior.vertices().size();
   const size_t num_cap_vertices = cap.vertices().size();
   const size_t num_total_vertices = num_ext_vertices + num_cap_vertices;
   const double radius = 1e-06; // TODO: express in ULP

   AABBTree tree(0.05, num_total_vertices);

   // Add exterior mesh vertices
   for (size_t ii=0; ii<num_ext_vertices; ++ii)
   {
      tree.addSphere(ii, exterior.vertices()[ii], radius);
   }

   // Add cap mesh vertices
   for (uint32_t ii=0; ii<num_cap_vertices; ++ii)
   {
      const uint32_t idx = ii+num_ext_vertices;
      tree.addSphere(idx, cap.vertices()[ii], radius);
   }

   // Query the tree with every vertex to detect duplicates
   // 'vertex_map' maps from the unwelded vertex numbering to the welded
   // numbering
   std::vector<uint32_t> vertex_map(num_total_vertices, UINT32_MAX);
   for (size_t ii=0; ii<num_ext_vertices; ++ii)
   {
      auto dupes = tree.query(ii);
      if (dupes.size() == 0)
      {
         vertex_map[ii] = ii;
      }
      else
      {
         std::sort(dupes.begin(), dupes.end(), std::less<uint32_t>());
         vertex_map[ii] = dupes[0]; // use the first (smallest)
      }
   }
   for (size_t ii=0; ii<num_cap_vertices; ++ii)
   {
      const uint32_t idx = ii+num_ext_vertices;
      auto dupes = tree.query(idx);
      if (dupes.size() == 0)
      {
         vertex_map[idx] = idx;
      }
      else
      {
         std::sort(dupes.begin(), dupes.end(), std::less<uint32_t>());
         vertex_map[idx] = dupes[0]; // use the first (smallest)
      }
   }

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
         const uint32_t idx = exterior.faces()[ii].m_v[jj];
         if (used_map[idx] == UINT32_MAX)
         {
            vertices.push_back( exterior.vertices()[idx] );
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
         const uint32_t idx = cfvi + num_ext_vertices;
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


void CSGEngine::construct(CSGOperation operation, bool cap, TriMesh& out_A, TriMesh& out_B)
{
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

   // Classify all faces as either above the triangle that intersected them (w.r.t. the
   // intersecting triangle's face normal) or below.
   // 1. Classify cut faces based on the face that cut them
   // 2. Construct a vertex adjacency map of all faces that share a primal
   //    vertex (i.e. verticies on the original mesh; not ones created by
   //    the retriangulation process)
   // 3. Do a breadth-first fill of all unclassified faces
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
       out_A = mergeMeshes(clay_above, knife_below);
       out_B = mergeMeshes(clay_below, knife_below);
   }
   else
   {
       out_A = clay_above;
       out_B = clay_below;
   }
}

}  // namespace CSG
