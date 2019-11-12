#include <string>
#include <Eigen/Dense>

#include "libcsg.h"




namespace CSG
{

// Triangle
std::ostream& operator<<(std::ostream& os, const Triangle& tri)
{
   os << "v: " << tri.m_v[0] << ", " << tri.m_v[1] << ", " << tri.m_v[2] << "  "
      << "uv: " << tri.m_uv[0] << ", " << tri.m_uv[1] << ", " << tri.m_uv[2]  << "  "
      << "n: "  << tri.m_n[0] << ", " << tri.m_n[1] << ", " << tri.m_n[2];

   return os;
}

/////


inline double orient2d(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c)
{
   return ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]));
}



bool intersectionTestVertex(const Eigen::Vector2d& P1, const Eigen::Vector2d& Q1,
                            const Eigen::Vector2d& R1, const Eigen::Vector2d& P2,
                            const Eigen::Vector2d& Q2, const Eigen::Vector2d& R2)
{
   if (orient2d(R2,P2,Q1) >= 0.0)
      if (orient2d(R2,Q2,Q1) <= 0.0)
         if (orient2d(P1,P2,Q1) > 0.0)
         {
            if (orient2d(P1,Q2,Q1) <= 0.0)
               return true;
            else
               return false;
         }
         else
         {
            if (orient2d(P1,P2,R1) >= 0.0)
               if (orient2d(Q1,R1,P2) >= 0.0)
                  return true;
               else
                  return false;
            else
               return false;
         }
      else
         if (orient2d(P1,Q2,Q1) <= 0.0)
            if (orient2d(R2,Q2,R1) <= 0.0)
               if (orient2d(Q1,R1,Q2) >= 0.0)
                  return true;
               else
                  return false;
            else
               return false;
         else
            return false;
   else
      if (orient2d(R2,P2,R1) >= 0.0)
         if (orient2d(Q1,R1,R2) >= 0.0)
            if (orient2d(P1,P2,R1) >= 0.0)
               return true;
            else
               return false;
         else
            if (orient2d(Q1,R1,Q2) >= 0.0)
            {
               if (orient2d(R2,R1,Q2) >= 0.0)
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
   if (orient2d(R2,P2,Q1) >= 0.0)
   {
      if (orient2d(P1,P2,Q1) >= 0.0)
      {
         if (orient2d(P1,Q1,R2) >= 0.0)
            return true;
         else
            return false;
      }
      else
      {
         if (orient2d(Q1,R1,P2) >= 0.0)
         {
            if (orient2d(R1,P1,P2) >= 0.0)
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
      if (orient2d(R2,P2,R1) >= 0.0)
      {
         if (orient2d(P1,P2,R1) >= 0.0)
         {
            if (orient2d(P1,R1,R2) >= 0.0)
               return true;
            else
            {
               if (orient2d(Q1,R1,R2) >= 0.0)
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
   if (orient2d(p2,q2,p1) >= 0.0)
   {
      if (orient2d(q2,r2,p1) >= 0.0)
      {
         if (orient2d(r2,p2,p1) >= 0.0)
            return true;
         else
            return intersectionTestEdge(p1,q1,r1,p2,q2,r2);
      }
      else
      {
         if ( orient2d(r2,p2,p1) >= 0.0 )
            return intersectionTestEdge(p1,q1,r1,r2,p2,q2);
         else
            return intersectionTestVertex(p1,q1,r1,p2,q2,r2);
      }
   }
   else
   {
      if ( orient2d(q2,r2,p1) >= 0.0 )
      {
         if ( orient2d(r2,p2,p1) >= 0.0 )
            return intersectionTestEdge(p1,q1,r1,q2,r2,p2);
         else
            return intersectionTestVertex(p1,q1,r1,q2,r2,p2);
      }
      else
         return intersectionTestVertex(p1,q1,r1,r2,p2,q2);
   }
}

bool triTriOverlapTest2d(const Eigen::Vector2d& p1, const Eigen::Vector2d& q1,
                             const Eigen::Vector2d& r1, const Eigen::Vector2d& p2,
                             const Eigen::Vector2d& q2, const Eigen::Vector2d& r2)
{
   if ( orient2d(p1,q1,r1) < 0.0 )
      if ( orient2d(p2,q2,r2) < 0.0 )
         return ccwTriTriIntersection2d(p1,r1,q1,p2,r2,q2);
      else
         return ccwTriTriIntersection2d(p1,r1,q1,p2,q2,r2);
  else
     if ( orient2d(p2,q2,r2) < 0.0f )
        return ccwTriTriIntersection2d(p1,q1,r1,p2,r2,q2);
     else
        return ccwTriTriIntersection2d(p1,q1,r1,p2,q2,r2);
}

bool coplanarTriTri3d(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1,
                      const Eigen::Vector3d& r1, const Eigen::Vector3d& p2,
                      const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
                      const Eigen::Vector3d& normal_1, const Eigen::Vector3d& normal_2)
{
   Eigen::Vector2d P1, Q1, R1;
   Eigen::Vector2d P2, Q2, R2;

   double n_x, n_y, n_z;
   n_x = ((normal_1[0]<0)?-normal_1[0]:normal_1[0]);
   n_y = ((normal_1[1]<0)?-normal_1[1]:normal_1[1]);
   n_z = ((normal_1[2]<0)?-normal_1[2]:normal_1[2]);


   /* Projection of the triangles in 3D onto 2D such that the area of
      the projection is maximized. */

   if (( n_x > n_z ) && ( n_x >= n_y ))
   {
      // Project onto plane YZ
      P1[0] = q1[2]; P1[1] = q1[1];
      Q1[0] = p1[2]; Q1[1] = p1[1];
      R1[0] = r1[2]; R1[1] = r1[1];

      P2[0] = q2[2]; P2[1] = q2[1];
      Q2[0] = p2[2]; Q2[1] = p2[1];
      R2[0] = r2[2]; R2[1] = r2[1];
   }
   else if (( n_y > n_z ) && ( n_y >= n_x ))
   {
      // Project onto plane XZ
      P1[0] = q1[0]; P1[1] = q1[2];
      Q1[0] = p1[0]; Q1[1] = p1[2];
      R1[0] = r1[0]; R1[1] = r1[2];

      P2[0] = q2[0]; P2[1] = q2[2];
      Q2[0] = p2[0]; Q2[1] = p2[2];
      R2[0] = r2[0]; R2[1] = r2[2];
   }
   else
   {
      // Project onto plane XY
      P1[0] = p1[0]; P1[1] = p1[1];
      Q1[0] = q1[0]; Q1[1] = q1[1];
      R1[0] = r1[0]; R1[1] = r1[1];

      P2[0] = p2[0]; P2[1] = p2[1];
      Q2[0] = q2[0]; Q2[1] = q2[1];
      R2[0] = r2[0]; R2[1] = r2[1];
   }

   return triTriOverlapTest2d(P1,Q1,R1,P2,Q2,R2);
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
         return checkMinMax(p1,r1,q1,r2,p2,q2, n1);
      else if (dr2 > 0.0)
         return checkMinMax(p1,r1,q1,q2,r2,p2, n1);
      else
         return checkMinMax(p1,q1,r1,p2,q2,r2, n1);
   }
  else if (dp2 < 0.0)
  {
     if (dq2 < 0.0)
        return checkMinMax(p1,q1,r1,r2,p2,q2, n1);
     else if (dr2 < 0.0)
        return checkMinMax(p1,q1,r1,q2,r2,p2, n1);
     else
        return checkMinMax(p1,r1,q1,p2,q2,r2, n1);
  }
  else
  {
     if (dq2 < 0.0)
     {
        if (dr2 >= 0.0)
           return checkMinMax(p1,r1,q1,q2,r2,p2, n1);
        else
           return checkMinMax(p1,q1,r1,p2,q2,r2, n1);
     }
     else if (dq2 > 0.0)
     {
        if (dr2 > 0.0)
           return checkMinMax(p1,r1,q1,p2,q2,r2, n1);
        else
           return checkMinMax(p1,q1,r1,q2,r2,p2, n1);
     }
     else
     {
        if (dr2 > 0.0)
           return checkMinMax(p1,q1,r1,r2,p2,q2, n1);
        else if (dr2 < 0.0)
           return checkMinMax(p1,r1,q1,r2,p2,q2, n1);
        else
           return coplanarTriTri3d(p1,q1,r1,p2,q2,r2, n1, n2);
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
        return triTri3d(r1,p1,q1,p2,r2,q2,n1,n2,dp2,dr2,dq2);
     else if (dr1 > 0.0)
        return triTri3d(q1,r1,p1,p2,r2,q2,n1,n2,dp2,dr2,dq2);
     else
        return triTri3d(p1,q1,r1,p2,q2,r2,n1,n2,dp2,dq2,dr2);
  }
  else if (dp1 < 0.0)
  {
     if (dq1 < 0.0)
        return triTri3d(r1,p1,q1,p2,q2,r2,n1,n2,dp2,dq2,dr2);
     else if (dr1 < 0.0)
        return triTri3d(q1,r1,p1,p2,q2,r2,n1,n2,dp2,dq2,dr2);
     else
        return triTri3d(p1,q1,r1,p2,r2,q2,n1,n2,dp2,dr2,dq2);
  }
  else
  {
     if (dq1 < 0.0)
     {
        if (dr1 >= 0.0)
           return triTri3d(q1,r1,p1,p2,r2,q2,n1,n2,dp2,dr2,dq2);
        else
           return triTri3d(p1,q1,r1,p2,q2,r2,n1,n2,dp2,dq2,dr2);
     }
     else if (dq1 > 0.0)
     {
        if (dr1 > 0.0)
           return triTri3d(p1,q1,r1,p2,r2,q2,n1,n2,dp2,dr2,dq2);
        else
           return triTri3d(q1,r1,p1,p2,q2,r2,n1,n2,dp2,dq2,dr2);
     }
     else
     {
        if (dr1 > 0.0)
           return triTri3d(r1,p1,q1,p2,q2,r2,n1,n2,dp2,dq2,dr2);
        else if (dr1 < 0.0)
           return triTri3d(r1,p1,q1,p2,r2,q2,n1,n2,dp2,dr2,dq2);
        else
           return coplanarTriTri3d(p1,q1,r1,p2,q2,r2,n1,n2);
     }
  }
}






/*
 *
 *  Three-dimensional Triangle-Triangle Intersection
 *
 */

/*
   This macro is called when the triangles surely intersect
   It constructs the segment of intersection of the two triangles
   if they are not coplanar.
*/

bool constructIntersection(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1,
                           const Eigen::Vector3d& r1, const Eigen::Vector3d& p2,
                           const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
                           const Eigen::Vector3d& N1, const Eigen::Vector3d& N2,
                           Eigen::Vector3d& source, Eigen::Vector3d& target)
{
   double alpha;
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
            v1 = alpha*v2;
            source = p1 - v1;
            v1 = p2 - p1;
            v2 = p2 - r2;
            alpha = v1.dot(N1) / v2.dot(N1);
            v1 = alpha*v2;
            target = p2 -v1;
            return true;
         }
         else
         {
            v1 = p2 - p1;
            v2 = p2 - q2;
            alpha = v1.dot(N1) / v2.dot(N1);
            v1 = alpha*v2;
            source = p2 - v1;
            v1 = p2 - p1;
            v2 = p2 - r2;
            alpha = v1.dot(N1) / v2.dot(N1);
            v1 = alpha*v2;
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
            v1 = alpha*v2;
            source = p1 - v1;
            v1 = p1 - p2;
            v2 = p1 - q1;
            alpha = v1.dot(N2) / v2.dot(N2);
            v1 = alpha*v2;
            target = p1 - v1;
            return true;
         }
         else
         {
            v1 = p2 - p1;
            v2 = p2 - q2;
            alpha = v1.dot(N1) / v2.dot(N1);
            v1 = alpha*v2;
            source = p2 - v1;
            v1 = p1 - p2;
            v2 = p1 - q1;
            alpha = v1.dot(N2) / v2.dot(N2);
            v1 = alpha*v2;
            target = p1 - v1;
            return true;
         }
      }
   }
}


bool triTriInter3d(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1, const Eigen::Vector3d& r1,
                   const Eigen::Vector3d& p2, const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
                   const Eigen::Vector3d& N1, const Eigen::Vector3d& N2,
                   Eigen::Vector3d& source, Eigen::Vector3d& target,
                   const double dp2, const double dq2, const double dr2, bool& coplanar)
{
    if (dp2 > 0.0)
    {
        if (dq2 > 0.0)
           return constructIntersection(p1,r1,q1,r2,p2,q2,N1,N2,source,target);
        else if (dr2 > 0.0)
            return constructIntersection(p1,r1,q1,q2,r2,p2,N1,N2,source,target);
        else
            return constructIntersection(p1,q1,r1,p2,q2,r2,N1,N2,source,target);
    }
    else if (dp2 < 0.0)
    {
        if (dq2 < 0.0)
            return constructIntersection(p1,q1,r1,r2,p2,q2,N1,N2,source,target);
        else if (dr2 < 0.0)
            return constructIntersection(p1,q1,r1,q2,r2,p2,N1,N2,source,target);
        else
            return constructIntersection(p1,r1,q1,p2,q2,r2,N1,N2,source,target);
    }
    else
    {
        if (dq2 < 0.0)
        {
            if (dr2 >= 0.0)
                return constructIntersection(p1,r1,q1,q2,r2,p2,N1,N2,source,target);
            else
                return constructIntersection(p1,q1,r1,p2,q2,r2,N1,N2,source,target);
        }
        else if (dq2 > 0.0)
        {
            if (dr2 > 0.0)
                return constructIntersection(p1,r1,q1,p2,q2,r2,N1,N2,source,target);
            else
                return constructIntersection(p1,q1,r1,q2,r2,p2,N1,N2,source,target);
        }
        else
        {
            if (dr2 > 0.0)
                return constructIntersection(p1,q1,r1,r2,p2,q2,N1,N2,source,target);
            else if (dr2 < 0.0)
                return constructIntersection(p1,r1,q1,r2,p2,q2,N1,N2,source,target);
            else
            {
                coplanar = true;
                return coplanarTriTri3d(p1,q1,r1,p2,q2,r2,N1,N2);
            }
        }
    }
}


/*
   The following version computes the segment of intersection of the
   two triangles if it exists.
   coplanar returns whether the triangles are coplanar
   source and target are the endpoints of the line segment of intersection
*/

bool triTriIntersectionTest3d(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1,
                              const Eigen::Vector3d& r1, const Eigen::Vector3d& p2,
                              const Eigen::Vector3d& q2, const Eigen::Vector3d& r2,
                              bool& coplanar, Eigen::Vector3d& source, Eigen::Vector3d& target)
{
   double dp1, dq1, dr1, dp2, dq2, dr2;
   Eigen::Vector3d v1, v2, v;
   Eigen::Vector3d N1, N2, N;

   coplanar = false;

   // Compute distance signs  of p1, q1 and r1
   // to the plane of triangle(p2,q2,r2)
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
         return triTriInter3d(r1,p1,q1,p2,r2,q2,N1,N2,source,target,dp2,dr2,dq2,coplanar);
      else if (dr1 > 0.0)
         return triTriInter3d(q1,r1,p1,p2,r2,q2,N1,N2,source,target,dp2,dr2,dq2,coplanar);
      else
         return triTriInter3d(p1,q1,r1,p2,q2,r2,N1,N2,source,target,dp2,dq2,dr2,coplanar);
   }
   else if (dp1 < 0.0)
   {
      if (dq1 < 0.0)
         return triTriInter3d(r1,p1,q1,p2,q2,r2,N1,N2,source,target,dp2,dq2,dr2,coplanar);
      else if (dr1 < 0.0)
         return triTriInter3d(q1,r1,p1,p2,q2,r2,N1,N2,source,target,dp2,dq2,dr2,coplanar);
      else
         return triTriInter3d(p1,q1,r1,p2,r2,q2,N1,N2,source,target,dp2,dr2,dq2,coplanar);
   }
   else
   {
      if (dq1 < 0.0)
      {
         if (dr1 >= 0.0)
            return triTriInter3d(q1,r1,p1,p2,r2,q2,N1,N2,source,target,dp2,dr2,dq2,coplanar);
         else
            return triTriInter3d(p1,q1,r1,p2,q2,r2,N1,N2,source,target,dp2,dq2,dr2,coplanar);
      }
      else if (dq1 > 0.0)
      {
         if (dr1 > 0.0)
            return triTriInter3d(p1,q1,r1,p2,r2,q2,N1,N2,source,target,dp2,dr2,dq2,coplanar);
         else
            return triTriInter3d(q1,r1,p1,p2,q2,r2,N1,N2,source,target,dp2,dq2,dr2,coplanar);
      }
      else
      {
         if (dr1 > 0.0)
            return triTriInter3d(r1,p1,q1,p2,q2,r2,N1,N2,source,target,dp2,dq2,dr2,coplanar);
         else if (dr1 < 0.0)
            return triTriInter3d(r1,p1,q1,p2,r2,q2,N1,N2,source,target,dp2,dr2,dq2,coplanar);
         else
         {
            // triangles are co-planar
            coplanar = true;
            return coplanarTriTri3d(p1,q1,r1,p2,q2,r2,N1,N2);
         }
      }
   }
}





/*
*
*  Two dimensional Triangle-Triangle Overlap Test
*
*/


TriangleIntersection Triangle::intersect(
   const Triangle& tri,
   const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& vertices )
{
   TriangleIntersection retval;

   retval.intersect =  triTriIntersectionTest3d(
      vertices[m_v[0]], vertices[m_v[1]], vertices[m_v[2]],
      vertices[tri.m_v[0]], vertices[tri.m_v[1]], vertices[tri.m_v[2]],
      retval.coplanar, retval.p0, retval.p1);

   return retval;
}


} // namespace CSG
