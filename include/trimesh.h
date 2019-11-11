#include <Eigen/Core>
#include <Eigen/StdVector>

#include <string>
#include <iostream>

#include "aabb.h"

#pragma once

namespace CSG
{
//> Triangle holds information for a single triangle, including indices into
//> its owning TriMesh
class Triangle
{
 public:
   static constexpr int32_t NO_ENTRY = -1;
   friend std::ostream& operator<<(std::ostream& os, const Triangle& tri);

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

   const std::vector<Triangle>& faces();

   const std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>>& vertices()
   { return m_vertices; }

   const std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>>& normals()
   { return m_normals; }

   const std::vector<Eigen::Vector2f, Eigen::aligned_allocator<Eigen::Vector2f>>& uvs()
   { return m_uvs; }

   AABBTree createAABBTree() const;

 private:
   // Primal data
   std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> m_vertices;
   std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> m_normals;
   std::vector<Eigen::Vector2f, Eigen::aligned_allocator<Eigen::Vector2f>> m_uvs;
   // Faces index into the vectors above
   std::vector<Triangle> m_faces;
};


}  //namespace CSG
