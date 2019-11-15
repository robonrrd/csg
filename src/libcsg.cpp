#include <bitset>
#include <unordered_set>

#include <Eigen/Dense>

#include "libcsg.h"
extern "C"
{
#include "triangle.h"
}

#define DEBUG

namespace CSG
{
// Magic numbers and tolerances. We seek to minimize the number of
// special constants, and they create edge case and problems, but some times
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
std::vector<IPoint> convertToIPoints(
    const TriangleIntersection& ix,
    const Triangle& ct, uint32_t cidx,
    const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& c_verts,
    const Triangle& kt, uint32_t kidx,
    const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& k_verts)
{
   std::vector<IPoint> pts;

   // Fill out non-parent fields
   IPoint p[2];
   for (uint32_t ii = 0; ii < 2; ++ii)
   {
      p[ii].pos = ix.p[ii];
      p[ii].cidx = cidx;
      p[ii].kidx = kidx;
   }

   // Determine the 'parenthood': is this point an exact duplicate of a point on
   // the clay mesh, the knife mesh or neither?
   for (uint32_t jj = 0; jj < 2; ++jj)
   {
      // Is it a point on the clay face?
      uint32_t clay_v_idx = 4;
      for (uint32_t ii = 0; ii < 3; ++ii)
      {
         if (point_almost_equal(ix.p[jj], c_verts[ct.m_v[ii]]))
         {
            clay_v_idx = ii;
            break;
         }
      }

      // Is it a point on the knife face?
      uint32_t knife_v_idx = 4;
      for (uint32_t ii = 0; ii < 3; ++ii)
      {
         if (point_almost_equal(ix.p[jj], k_verts[kt.m_v[ii]]))
         {
            knife_v_idx = ii;
            break;
         }
      }

      if ((clay_v_idx < 4) && (knife_v_idx == 4))
      {
         p[jj].parent = kClay;
         p[jj].p_idx = clay_v_idx;
      }
      else if ((clay_v_idx == 4) && (knife_v_idx < 4))
      {
         p[jj].parent = kKnife;
         p[jj].p_idx = knife_v_idx;
      }
      else if ((clay_v_idx < 4) && (knife_v_idx < 4))
      {
         p[jj].parent = kBoth;
         // TODO: should  we preserve both clay and knife vert indices here?
         p[jj].p_idx = clay_v_idx;
      }
      else
      {
         p[jj].parent = kNew;
         p[jj].p_idx = 0;  // don't care, but don't leave it uninitialized
      }
      pts.push_back(p[jj]);
   }

   std::cout << "alpha=" << ix.alpha << "  beta=" << ix.beta << std::endl;
   if (almost_equal(ix.alpha, 1.0, ALPHA_ULP))
   {
      std::cout << "alpha is almost 1.0" << std::endl;
   }
   else if (almost_equal(ix.alpha, 0.0, ALPHA_ULP))
   {
      std::cout << "alpha is almost 0.0" << std::endl;
   }

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


void CSG(const TriMesh& clay, const TriMesh& knife, CSGOperation operation, TriMesh& A, TriMesh& B)
{
   std::cout << "Building AABB trees.." << std::endl;
   AABBTree clay_tree = clay.createAABBTree();
   AABBTree knife_tree = knife.createAABBTree();
   std::cout << "..done" << std::endl;

   // Calculate potential triangle intersections
   std::vector<std::pair<uint32_t, uint32_t>> ix = clay_tree.intersect(knife_tree);
   std::cout << ix.size() << " AABB collisions found" << std::endl;

   // Vector to store new vertices. Both clay and knife will share these vertices
   std::vector<IPoint> new_verts;

   // Calculate actual triangle intersections, with intersection data
   const uint32_t ix_sz = ix.size();
   std::unordered_set<uint32_t> cut_clay_faces;
   std::unordered_set<uint32_t> cut_knife_faces;
   std::vector<bool> is_clay_face_cut(clay.faces().size(), false);
   std::vector<bool> is_knife_face_cut(knife.faces().size(), false);
   std::unordered_map<uint32_t, std::vector<uint32_t>> clay_face_new_verts;
   for (uint32_t ii = 0; ii < ix_sz; ++ii)
   {
      const uint32_t c_idx = ix[ii].first;
      const uint32_t k_idx = ix[ii].second;

      cut_clay_faces.insert(c_idx);
      cut_knife_faces.insert(k_idx);
      is_clay_face_cut[c_idx] = true;
      is_knife_face_cut[k_idx] = true;

      const Triangle& ct = clay.faces()[c_idx];
      const Triangle& kt = knife.faces()[k_idx];
      TriangleIntersection trix = intersect(ct, clay.vertices(), kt, knife.vertices());
      if (!trix.intersect)
         continue;

      if (trix.coplanar)
         continue;  // TODO: will handle this later

      std::vector<IPoint> pts = convertToIPoints(trix, ct, c_idx, clay.vertices(), kt, k_idx,
                                                 knife.vertices());
      new_verts.insert(new_verts.end(), pts.begin(), pts.end());

      auto itr = clay_face_new_verts.find(c_idx);
      if (itr == clay_face_new_verts.end())
         clay_face_new_verts.insert({c_idx, std::vector<uint32_t>()});
      clay_face_new_verts[c_idx].push_back(new_verts.size() - 1);
      clay_face_new_verts[c_idx].push_back(new_verts.size() - 2);
   }

   // Step through the cut faces and retriangulate them
   std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> new_clay_faces;
   for (auto itr = cut_clay_faces.begin(); itr != cut_clay_faces.end(); ++itr)
   {
      auto result = retriangulate(clay, *itr, new_verts, clay_face_new_verts[*itr]);
      new_clay_faces.insert(new_clay_faces.end(), result.begin(), result.end());
   }

}

}  // namespace CSG
