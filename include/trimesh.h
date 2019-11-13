#include <Eigen/Core>
#include <Eigen/StdVector>

#include <string>
#include <iostream>

#include "aabb.h"

#pragma once

namespace CSG
{


class TriangleIntersection
{
public:
   friend std::ostream& operator<<(std::ostream& os, const TriangleIntersection& ix)
   {
      if (!ix.intersect)
      {
         os << "no intersection";
      }
      else
      {
         os << "intersection. ";
         if (ix.coplanar)
            os << "coplanar. ";
         os << "p0:" << ix.p0[0] << ", " << ix.p0[1] << ", " << ix.p0[2] << "  "
            << "p1:" << ix.p1[0] << ", " << ix.p1[1] << ", " << ix.p1[2] << "  ";
      }
      return os;
   }

   bool intersect;     // is there an intersection at all?
   bool coplanar;      // are the triangles coplanar
   Eigen::Vector3d p0; // the two points defining the line of intersection
   Eigen::Vector3d p1;
};


//> Triangle holds information for a single triangle, including indices into
//> its owning TriMesh
class Triangle
{
public:
   static constexpr int32_t NO_ENTRY = -1;
   friend std::ostream& operator<<(std::ostream& os, const Triangle& tri);

   TriangleIntersection intersect(
      const Triangle& tri,
      const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& vertices);

   int32_t m_v[3];   //> vertices
   int32_t m_n[3];   //> normals
   int32_t m_uv[3];  //> UVs
};


//> TriMesh holds information defining a triangluated surface (actually, a triangle
//> soup, as no connectivity is checked or assumed)
class TriMesh
{
 public:
   // Load OBJ file and return the number of vertices in the new mesh.
   // Returning '0' means the operation failed
   uint32_t loadOBJ(const std::string& filename);

   const std::vector<Triangle>& faces() const { return m_faces; }

   const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& vertices() const
   { return m_vertices; }

   const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& normals() const
   { return m_normals; }

   const std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>& uvs() const
   { return m_uvs; }

   AABBTree createAABBTree() const;

 private:
   // Primal data
   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> m_vertices;
   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> m_normals;
   std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> m_uvs;
   // Faces index into the vectors above
   std::vector<Triangle> m_faces;
};


}  //namespace CSG
