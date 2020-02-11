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
   Triangle();

#ifndef SWIG
   friend std::ostream& operator<<(std::ostream& os, const Triangle& tri);
#endif

   int32_t m_v[3];   //> vertices
   int32_t m_n[3];   //> normals
   int32_t m_uv[3];  //> UVs
};


//> TriMesh holds information defining a triangluated surface (actually, a triangle
//> soup, as no connectivity is checked or assumed)
class TriMesh
{
 public:
   TriMesh();
   // Construct triangle mesh from flat vectors
   TriMesh(const std::vector<double>& vertices, const std::vector<unsigned int>& faces);

   // Load OBJ file and return the number of vertices in the new mesh.
   // Returning '0' means the operation failed
   uint32_t loadOBJ(const std::string& filename);

   const std::vector<Triangle>& faces() const { return m_faces; }

#ifdef SWIG
   // Don't expose Eigen types to Python
   std::vector<double>& vertices() const;
   std::vector<double>& normals() const;// vertex normals
   std::vector<double>& uvs() const;
//   std::vector<double> faceNormal(uint32_t f) const; // face normals
#else
   const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& vertices() const
   { return m_vertices; }

   const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& normals() const
   { return m_normals; }

   const std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>& uvs() const
   { return m_uvs; }

   Eigen::Vector3d faceNormal(uint32_t f) const;
#endif

   AABBTree createAABBTree() const;

 private:
   Eigen::Vector3d faceNormalImpl(uint32_t f) const;

   // Primal data
   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> m_vertices;
   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> m_normals;
   std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> m_uvs;
   // Faces index into the vectors above
   std::vector<Triangle> m_faces;
};


}  //namespace CSG
