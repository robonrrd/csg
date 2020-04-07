#include <Eigen/Core>
#include <Eigen/StdVector>
#include <string>
#include <iostream>
#include "aabb.h"

#pragma once

// Some macros to make the aligned std::vector code cleaner
#define VECTOR3D Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>
#define VECTOR2D Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>
#define VECTOR3I Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>


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
#ifndef SWIG
   TriMesh(const std::vector<VECTOR3D>& vertices,
           const std::vector<VECTOR3I>& faces);
#endif
   // Load OBJ file and return the number of vertices in the new mesh.
   // Returning '0' means the operation failed
   uint32_t loadOBJ(const std::string& filename);
#ifndef SWIG
   friend std::ostream &operator<<( std::ostream &output, const TriMesh &M );
#endif

   // Non-Eigen types for Python binding
   std::vector<double> mesh_vertices() const;
   std::vector<double> mesh_normals() const;// vertex normals
   std::vector<double> mesh_uvs() const;
   std::vector<unsigned int> mesh_faces() const;

#ifndef SWIG
   // Accessors that we don't wish to expose to Python
   const std::vector<VECTOR3D>& vertices() const      { return m_vertices; }
   std::vector<VECTOR3D>& vertices()                  { return m_vertices; }

   const std::vector<VECTOR3D>& normals() const       { return m_normals; }
   std::vector<VECTOR3D>& normals()                   { return m_normals; }

   const std::vector<VECTOR2D>& uvs() const           { return m_uvs; }
   std::vector<VECTOR2D>& uvs()                       { return m_uvs; }

   const std::vector<Triangle>& faces() const         { return m_faces; }
   std::vector<Triangle>& faces()                     { return m_faces; }
   void setFaces(const std::vector<VECTOR3I>& in);

   Eigen::Vector3d faceNormal(uint32_t f) const;
#endif

   AABBTree createAABBTree() const;

 private:
   Eigen::Vector3d faceNormalImpl(uint32_t f) const;

   // Primal data
   std::vector<VECTOR3D> m_vertices;
   std::vector<VECTOR3D> m_normals;
   std::vector<VECTOR2D> m_uvs;
   // Faces index into the vectors above
   std::vector<Triangle> m_faces;
};


}  //namespace CSG
