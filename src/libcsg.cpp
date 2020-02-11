//
//
//
//
#include <assert.h>
#include <set>
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
// Magic numbers and tolerances. We seek to minimize the number of
// special constants as they create edge case and problems but some times
// they're necessary.
// All magic numbers are expressed in terms of 'units in the last place'
constexpr uint32_t ALPHA_ULP = 2;
constexpr uint32_t POINT_ULP = 2;


template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, uint32_t ulp)
{
   // the machine epsilon has to be scaled to the magnitude of the values used
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


bool checkMinMax(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1, const Eigen::Vector3d& r1,
                 const Eigen::Vector3d& p2, const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
                 Eigen::Vector3d& n1)
{
   Eigen::Vector3d v1 = p2 - q1;
   Eigen::Vector3d v2 = p1 - q1;
   n1 = v1.cross(v2);

   v1 = q2 - q1;
   if (v1.dot(n1) > 0.0)
      return false;

   v1 = p2 - p1;
   v2 = r1 - p1;
   n1 = v1.cross(v2);
   v1 = r2 - p1;
   if (v1.dot(n1) > 0.0)
      return false;

   return true;
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


/* Permutation in a canonical form of T2's vertices */
bool triTri3d(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1, const Eigen::Vector3d& r1,
              const Eigen::Vector3d& p2, const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
              Eigen::Vector3d& n1, Eigen::Vector3d& n2,
              const double dp2, const double dq2, const double dr2)
{
   if (dp2 > 0.0)
   {
      if (dq2 > 0.0)
         return checkMinMax(p1, r1, q1, r2, p2, q2, n1);
      else if (dr2 > 0.0)
         return checkMinMax(p1, r1, q1, q2, r2, p2, n1);
      else
         return checkMinMax(p1, q1, r1, p2, q2, r2, n1);
   }
   else if (dp2 < 0.0)
   {
      if (dq2 < 0.0)
         return checkMinMax(p1, q1, r1, r2, p2, q2, n1);
      else if (dr2 < 0.0)
         return checkMinMax(p1, q1, r1, q2, r2, p2, n1);
      else
         return checkMinMax(p1, r1, q1, p2, q2, r2, n1);
   }
   else
   {
      if (dq2 < 0.0)
      {
         if (dr2 >= 0.0)
            return checkMinMax(p1, r1, q1, q2, r2, p2, n1);
         else
            return checkMinMax(p1, q1, r1, p2, q2, r2, n1);
      }
      else if (dq2 > 0.0)
      {
         if (dr2 > 0.0)
            return checkMinMax(p1, r1, q1, p2, q2, r2, n1);
         else
            return checkMinMax(p1, q1, r1, q2, r2, p2, n1);
      }
      else
      {
         if (dr2 > 0.0)
            return checkMinMax(p1, q1, r1, r2, p2, q2, n1);
         else if (dr2 < 0.0)
            return checkMinMax(p1, r1, q1, r2, p2, q2, n1);
         else
            return coplanarTriTri3d(p1, q1, r1, p2, q2, r2, n1, n2);
      }
   }
}


/*
*
*  Three-dimensional Triangle-Triangle Overlap Test
*
*/


bool triTriOverlapTest3d(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1,
                         const Eigen::Vector3d& r1, const Eigen::Vector3d& p2,
                         const Eigen::Vector3d& q2, const Eigen::Vector3d& r2)
{
   double dp1, dq1, dr1, dp2, dq2, dr2;
   Eigen::Vector3d v1, v2;
   Eigen::Vector3d n1, n2;

   /* Compute distance signs  of p1, q1 and r1 to the plane of
      triangle(p2,q2,r2) */
   v1 = p2 - r2;
   v2 = q2 - r2;
   n2 = v1.cross(v2);

   v1 = p1 - r2;
   dp1 = v1.dot(n2);
   v1 = q1 - r2;
   dq1 = v1.dot(n2);
   v1 = r1 - r2;
   dr1 = v1.dot(n2);

   if (((dp1 * dq1) > 0.0) && ((dp1 * dr1) > 0.0))
      return false;

   /* Compute distance signs  of p2, q2 and r2 to the plane of
     triangle(p1,q1,r1) */
   v1 = q1 - p1;
   v2 = r1 - p1;
   n1 = v1.cross(v2);

   v1 = p2 - r1;
   dp2 = v1.dot(n1);
   v1 = q2 - r1;
   dq2 = v1.dot(n1);
   v1 = r2 - r1;
   dr2 = v1.dot(n1);

   if (((dp2 * dq2) > 0.0) && ((dp2 * dr2) > 0.0))
      return false;

   /* Permutation in a canonical form of T1's vertices */
   if (dp1 > 0.0)
   {
      if (dq1 > 0.0)
         return triTri3d(r1, p1, q1, p2, r2, q2, n1, n2, dp2, dr2, dq2);
      else if (dr1 > 0.0)
         return triTri3d(q1, r1, p1, p2, r2, q2, n1, n2, dp2, dr2, dq2);
      else
         return triTri3d(p1, q1, r1, p2, q2, r2, n1, n2, dp2, dq2, dr2);
   }
   else if (dp1 < 0.0)
   {
      if (dq1 < 0.0)
         return triTri3d(r1, p1, q1, p2, q2, r2, n1, n2, dp2, dq2, dr2);
      else if (dr1 < 0.0)
         return triTri3d(q1, r1, p1, p2, q2, r2, n1, n2, dp2, dq2, dr2);
      else
         return triTri3d(p1, q1, r1, p2, r2, q2, n1, n2, dp2, dr2, dq2);
   }
   else
   {
      if (dq1 < 0.0)
      {
         if (dr1 >= 0.0)
            return triTri3d(q1, r1, p1, p2, r2, q2, n1, n2, dp2, dr2, dq2);
         else
            return triTri3d(p1, q1, r1, p2, q2, r2, n1, n2, dp2, dq2, dr2);
      }
      else if (dq1 > 0.0)
      {
         if (dr1 > 0.0)
            return triTri3d(p1, q1, r1, p2, r2, q2, n1, n2, dp2, dr2, dq2);
         else
            return triTri3d(q1, r1, p1, p2, q2, r2, n1, n2, dp2, dq2, dr2);
      }
      else
      {
         if (dr1 > 0.0)
            return triTri3d(r1, p1, q1, p2, q2, r2, n1, n2, dp2, dq2, dr2);
         else if (dr1 < 0.0)
            return triTri3d(r1, p1, q1, p2, r2, q2, n1, n2, dp2, dr2, dq2);
         else
            return coplanarTriTri3d(p1, q1, r1, p2, q2, r2, n1, n2);
      }
   }
}


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
    const Triangle& tri_a,
    const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& verts_a,
    const Triangle& tri_b,
    const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& verts_b)
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

#ifdef DEBUG
   std::cout << "alpha=" << ix.alpha << "  beta=" << ix.beta << std::endl;
   if (almost_equal(ix.alpha, 1.0, ALPHA_ULP))
   {
      std::cout << "alpha is almost 1.0" << std::endl;
   }
   else if (almost_equal(ix.alpha, 0.0, ALPHA_ULP))
   {
      std::cout << "alpha is almost 0.0" << std::endl;
   }
#endif

   return pts;
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

#if 0
std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> retriangulate(
    const TriMesh& mesh, uint32_t fidx,
    const std::vector<IPoint>& new_verts,
    const std::vector<uint32_t>& new_vert_indices)
{
   const uint32_t numPoints = 3 + new_vert_indices.size();
   const uint32_t numSegments = new_vert_indices.size() / 2;

#ifdef DEBUG
   std::cout << "numPoints=" << numPoints << "  numSegments=" << numSegments << std::endl;
#endif

   // Rotate the triangle into the XY plane ('triangle' is2D only)
   const Eigen::Vector3d& p0 = mesh.vertices()[mesh.faces()[fidx].m_v[0]];
   const Eigen::Vector3d& p1 = mesh.vertices()[mesh.faces()[fidx].m_v[1]];
   const Eigen::Vector3d& p2 = mesh.vertices()[mesh.faces()[fidx].m_v[2]];

   const Eigen::Vector3d n = (p1 - p0).normalized().cross((p2 - p0).normalized());  // normal
   const Eigen::Vector3d axis = n.cross(Eigen::Vector3d::UnitZ()).normalized();
   const double angle = acos(n[2]);
   const Eigen::Matrix3d rot = Eigen::AngleAxisd(angle, axis).matrix();

   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> pts_2d(numPoints);
   pts_2d[0] = rot * p0;
   pts_2d[1] = rot * p1;
   pts_2d[2] = rot * p2;
   for (uint32_t ii = 0; ii < new_vert_indices.size(); ++ii)
      pts_2d[3 + ii] = rot * new_verts[new_vert_indices[ii]].pos;


#ifdef DEBUG
   std::cout << "Retriangulating face " << fidx << std::endl
             << " Triangle points:" << std::endl
             << "  p0: " << p0.transpose() << std::endl
             << "  p1: " << p1.transpose() << std::endl
             << "  p2: " << p2.transpose() << std::endl
             << " Intersection points:" << std::endl;
   for (uint32_t ii = 0; ii < new_vert_indices.size(); ++ii)
      std::cout << "  i" << ii << ": " << new_verts[new_vert_indices[ii]].pos.transpose()
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

   char flags[] = "cz";  // pcze
   triangulate(flags, &in, &out, 0);

   std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>>
       new_tris(out.numberoftriangles);

   for (uint32_t ii = 0; ii < out.numberoftriangles; ii++)
      for (uint32_t jj = 0; jj < 3; jj++)
         new_tris[ii][jj] = out.trianglelist[ii * 3 + jj];

   free(in.pointlist);
   free(in.segmentlist);

   return new_tris;
}
#endif


std::vector<std::pair<uint32_t, uint32_t>> intersectAABBs(const TriMesh& clay, const TriMesh& knife)
{
   std::cout << "Building AABB trees.." << std::endl;
   AABBTree clay_tree = clay.createAABBTree();
   AABBTree knife_tree = knife.createAABBTree();
   std::cout << " ..done" << std::endl;

   // Calculate potential triangle intersections
   std::vector<std::pair<uint32_t, uint32_t>> ix = clay_tree.intersect(knife_tree);
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


std::vector<IFace> CSGEngine::retriangulate(const TriMesh& mesh, IParent which_mesh, uint32_t fidx,
                                            const std::vector<uint32_t>& new_vert_indices) const
{
   const uint32_t numPoints = 3 + new_vert_indices.size();
   const uint32_t numSegments = new_vert_indices.size() / 2;

#ifdef DEBUG
   std::cout << "numPoints=" << numPoints << "  numSegments=" << numSegments << std::endl;
#endif

   // Rotate the triangle into the XY plane ('triangle' is2D only)
   const Eigen::Vector3d& p0 = mesh.vertices()[mesh.faces()[fidx].m_v[0]];
   const Eigen::Vector3d& p1 = mesh.vertices()[mesh.faces()[fidx].m_v[1]];
   const Eigen::Vector3d& p2 = mesh.vertices()[mesh.faces()[fidx].m_v[2]];

   const Eigen::Vector3d n = (p1 - p0).normalized().cross((p2 - p0).normalized());  // normal
   const Eigen::Vector3d axis = n.cross(Eigen::Vector3d::UnitZ()).normalized();
   const double angle = acos(n[2]);
   const Eigen::Matrix3d rot = Eigen::AngleAxisd(angle, axis).matrix();

   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> pts_2d(numPoints);
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

   char flags[] = "cz";  // pcze
   triangulate(flags, &in, &out, 0);

   // Convert the locally-indices triangle vertex indices into IPointRefs
   std::vector<IFace> new_faces(out.numberoftriangles);
   for (uint32_t ii = 0; ii < out.numberoftriangles; ii++)
   {
      new_faces[ii].orig = fidx;
      for (uint32_t jj = 0; jj < 3; jj++)
      {
         uint32_t idx = out.trianglelist[ii * 3 + jj];
         if (idx < 3) // original vertices
         {
            new_faces[ii].v[jj].parent = which_mesh;
            new_faces[ii].v[jj].idx = mesh.faces()[fidx].m_v[idx];
         }
         else
         {
            new_faces[ii].v[jj] = m_newPoints[new_vert_indices[idx-3]].ref;
         }
      }
   }

#if 0
   // dump raw output
   std::cout << "# raw output from triangulate" << std::endl;
   std::cout << "v " << p0.transpose() << std::endl;
   std::cout << "v " << p1.transpose() << std::endl;
   std::cout << "v " << p2.transpose() << std::endl;

   for (uint32_t ii = 0; ii < new_vert_indices.size(); ++ii)
       std::cout << "v " << ipointPos(m_newPoints[new_vert_indices[ii]].ref).transpose()
                 << std::endl;

   std::cout << "# faces" << std::endl;
   for (uint32_t ii = 0; ii < out.numberoftriangles; ii++)
   {
       std::cout << "f";
       for (uint32_t jj = 0; jj < 3; jj++)
       {
           std::cout << " " << out.trianglelist[ii * 3 + jj]+1;
       }
       std::cout << std::endl;
   }
   std::cout << std::endl << std::endl;
#endif

   free(in.pointlist);
   free(in.segmentlist);

   return new_faces;

}


const Eigen::Vector3d centroid(const TriMesh& mesh, uint32_t ff)
{
   const Triangle& face = mesh.faces()[ff];
   return (mesh.vertices()[face.m_v[0]] + mesh.vertices()[face.m_v[1]] +
           mesh.vertices()[face.m_v[2]])/3.0;
}


std::vector<char>  CSGEngine::classifyCutFaces(const std::vector<IFace>& in_faces,
                                               IParent which_surface)
{
   // status:
   // -1 - below
   //  0 - unknown
   //  1 - above
   const size_t num_faces = in_faces.size();
   std::vector<char> status(num_faces, 0);
   uint32_t num_unknown = 0;
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

      int vote_above = 0;
      int vote_below = 0;
      // Test each vertex
      Eigen::Vector3d face_center(0,0,0);
      for (size_t ii=0; ii<3; ++ii)
      {
         const Eigen::Vector3d vpos = ipointPos(ff.v[ii]);
         face_center += vpos;

         if (ff.v[ii].parent == kNew)
         {
            Eigen::Vector3d vv = vpos - center;
            if (vv.dot(normal) > 0)
               vote_above++;
            else if (vv.dot(normal) < 0)
               vote_below++;
         }
      }

      // Test triangle center
      face_center /= 3.0f;
      Eigen::Vector3d vv = face_center - center;
      if (vv.dot(normal) > 0)
         vote_above++;
      else if (vv.dot(normal) < 0)
         vote_below++;

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
         status[fidx] = 0;
         num_unknown++;
      }
   }

   std::cout << num_unknown << " unclassified cut faces remain" << std::endl;
   return status;
}


void CSGEngine::construct(CSGOperation operation, bool cap, TriMesh& out_A, TriMesh& out_B)
{
   // Create and intersect AABB trees:
   std::vector<std::pair<uint32_t, uint32_t>> ix = intersectAABBs(m_clay, m_knife);
   const uint32_t ix_sz = ix.size();

   // Calculate actual triangle intersections, with intersection data

   // Sets of indices of faces that were cut
   std::unordered_set<uint32_t> cut_clay_faces;
   std::unordered_set<uint32_t> cut_knife_faces;
   // Vectors of bools, indicating whether the face was cut; faster than seeking
   // within an unordered_set
   std::vector<bool> is_clay_face_cut(m_clay.faces().size(), false);
   std::vector<bool> is_knife_face_cut(m_knife.faces().size(), false);

   // Indices into the m_newPoints vector, which contains IPoints which refer to
   // either clay or knife mesh verticies, or m_newPointPositions
   std::unordered_map<uint32_t, std::vector<uint32_t>> clay_face_new_verts;
   std::unordered_map<uint32_t, std::vector<uint32_t>> knife_face_new_verts;

   for (uint32_t ii = 0; ii < ix_sz; ++ii)
   {
      const uint32_t c_idx = ix[ii].first;
      const uint32_t k_idx = ix[ii].second;

      const Triangle& ct = m_clay.faces()[c_idx];
      const Triangle& kt = m_knife.faces()[k_idx];
      TriangleIntersection trix = intersect(ct, m_clay.vertices(), kt, m_knife.vertices());
      if (!trix.intersect) // just because AABBs intersect doesn't mean the faces do
         continue;

      if (trix.coplanar)
         continue;  // TODO: will handle this later

      // Record which faces have been cut, so we can ignore replace them later
      // with the diced up versions
      cut_clay_faces.insert(c_idx);
      cut_knife_faces.insert(k_idx);
      is_clay_face_cut[c_idx] = true;
      is_knife_face_cut[k_idx] = true;

      // Convert the raw positions to indices into one of three arrays: existing clay mesh
      // vertices, existing knife mesh vertices, or new vertices (which may be duplicates
      // of other new vertices; we'll clean this up later)
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

#if 0
   // Dump new faces to obj
   std::cout << "# Clay vertices" << std::endl;
   for (uint32_t ii=0; ii<m_clay.vertices().size(); ++ii)
   {
      std::cout << "v " << m_clay.vertices()[ii].transpose() << std::endl;
   }
   std::cout << "# New vertices" << std::endl;
   for (uint32_t ii=0; ii<m_newPointPositions.size(); ++ii)
   {
      std::cout << "v " << m_newPointPositions[ii].transpose() << std::endl;
   }

   for (uint32_t ii=0; ii<new_clay_faces.size(); ++ii)
   {
       std::cout << "f ";
       for (uint32_t jj=0; jj<3; ++jj)
       {
           const auto& pt = new_clay_faces[ii].v[jj];
           if (pt.parent == kClay)
               std::cout << pt.idx+1 << " ";
           else if (pt.parent == kNew)
               std::cout << pt.idx+m_clay.vertices().size()+1 << " ";
       }
       std::cout << std::endl;
   }
#endif

   std::vector<IFace> new_knife_faces;
   if (cap)
   {
       for (auto itr = cut_knife_faces.begin(); itr != cut_knife_faces.end(); ++itr)
       {
           auto result = retriangulate(m_knife, kKnife, *itr, knife_face_new_verts[*itr]);
           new_knife_faces.insert(new_knife_faces.end(), result.begin(), result.end());
       }
   }

   // Classify cut faces into 'above' and 'below' sets, depending on if they are
   // above or below the faces that cut them (with respect to that face's normal)
   std::vector<char> cut_face_status = classifyCutFaces(new_clay_faces, kClay);
//   std::vector<IFace> clay_above, clay_below;

}


}  // namespace CSG
