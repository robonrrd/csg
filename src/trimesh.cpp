#include "libcsg.h"

#include <fstream>
#include <iterator>
#include <string>
#include <stdlib.h>

static constexpr int32_t NO_ENTRY = -1;


// Helper functions for OBJ parsing
bool triplet_valid(const std::vector<uint32_t>& idx, uint32_t ii)
{
   return ((idx[ii] > 0) && (idx[ii + 1] > 0) && (idx[ii + 2] > 0));
}

const std::string WHITESPACE = " \n\r\t\f\v";
std::string ltrim(const std::string& s)
{
   size_t start = s.find_first_not_of(WHITESPACE);
   return (start == std::string::npos) ? "" : s.substr(start);
}

template<char delimiter>
class WordDelimitedBy : public std::string
{};

template<char T>
std::istream& operator>>(std::istream& is, WordDelimitedBy<T>& output)
{
   std::getline(is, output, T);
   return is;
}


namespace CSG
{


// Triangle class member functions
//
#ifndef SWIG
std::ostream& operator<<(std::ostream& os, const Triangle& tri)
{
   os << "v: " << tri.m_v[0] << ", " << tri.m_v[1] << ", " << tri.m_v[2] << "  "
      << "uv: " << tri.m_uv[0] << ", " << tri.m_uv[1] << ", " << tri.m_uv[2]  << "  "
      << "n: "  << tri.m_n[0] << ", " << tri.m_n[1] << ", " << tri.m_n[2];

   return os;
}
#endif


Triangle::Triangle()
{
   for (uint32_t ii=0; ii<3; ++ii)
   {
      m_v[ii]  = NO_ENTRY;
      m_n[ii]  = NO_ENTRY;
      m_uv[ii] = NO_ENTRY;
   }
}


// TriMesh class member functions
//
TriMesh::TriMesh()
{
   // empty
}


TriMesh::TriMesh(const std::vector<double>& in_vertices, const std::vector<unsigned int>& in_faces)
{
   const uint32_t numVertices = in_vertices.size()/3;
   m_vertices.resize(numVertices);
   for (uint32_t ii=0; ii<numVertices; ++ii)
      m_vertices[ii] = Eigen::Vector3d(in_vertices[3*ii], in_vertices[3*ii+1], in_vertices[3*ii+2]);

   const uint32_t numFaces = in_faces.size()/3;
   m_faces.resize(numFaces);
   for (uint32_t ii=0; ii<numFaces; ++ii)
      for (uint32_t jj=0; jj<3; ++jj)
      {
         m_faces[ii].m_v[jj] = in_faces[3*ii+jj];
         m_faces[ii].m_n[jj] = NO_ENTRY;
         m_faces[ii].m_uv[jj] = NO_ENTRY;
      }
}

#ifndef SWIG
TriMesh::TriMesh(
    const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& in_vertices,
    const std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>>& in_faces)
    : m_vertices(in_vertices)
{
    const uint32_t numFaces = in_faces.size();
    m_faces.resize(numFaces);
    for (uint32_t ii=0; ii<numFaces; ++ii)
    {
        for (uint32_t jj=0; jj<3; ++jj)
        {
            m_faces[ii].m_v[jj] = in_faces[ii][jj];
            m_faces[ii].m_n[jj] = NO_ENTRY;
            m_faces[ii].m_uv[jj] = NO_ENTRY;
        }
    }
}


void
TriMesh::setFaces(const std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>>& in)
{
    const size_t numFaces = in.size();
    m_faces = std::vector<Triangle>(numFaces);
    for (size_t ii=0; ii<numFaces; ++ii)
    {
        for (size_t jj=0; jj<3; ++jj)
        {
            m_faces[ii].m_v[jj] = in[ii][jj];
        }
    }
}
#endif

uint32_t TriMesh::loadOBJ(const std::string& path)
{
   m_vertices.clear();
   m_faces.clear();

   std::ifstream file (path);
   if (!file.is_open())
   {
      printf("Impossible to open the file !\n");
      return 0;
   }

   std::vector<uint32_t> vertexIndices, uvIndices, normalIndices;
   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> temp_vertices;
   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> temp_normals;
   std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> temp_uvs;

   std::string in_line;
   while (std::getline (file, in_line))
   {
      std::istringstream iss(ltrim(in_line));
      std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                       std::istream_iterator<std::string>());

      if (results.size() == 0) // skip blank lines
         continue;

      if ((results[0] == "v") && results.size() == 4)
      {
         Eigen::Vector3d vertex;
         for (int ii=0; ii<3; ++ii)
            vertex[ii] = strtof(results[ii+1].c_str(), NULL);
         temp_vertices.push_back(vertex);
      }

      else if ((results[0] == "vt") && results.size() == 3)
      {
         Eigen::Vector2d uv;
         for (int ii=0; ii<2; ++ii)
            uv[ii] = strtof(results[ii+1].c_str(), NULL);
         temp_uvs.push_back(uv);
      }

      else if ((results[0] == "vn") && results.size() == 4)
      {
         Eigen::Vector3d normal;
         for (int ii=0; ii<3; ++ii)
            normal[ii] = strtof(results[ii+1].c_str(), NULL);
         temp_normals.push_back(normal);
      }

      else if ((results[0] == "f") && results.size() == 4)
      {
         // A face can be specified (correctly) in three different ways:
         //  "f  1 2 3"          # simple vertex indices
         //  "f  1/10 2/11 3/13  # vertex indices plus UV indices
         //  "f  1/10/30  2/11/31  3/13/33 # vertex, UV and normal indices
         uint32_t vertexIndex[3] = {0, 0, 0};
         uint32_t uvIndex[3] = {0, 0, 0};
         uint32_t normalIndex[3] = {0, 0, 0};

         for (int ii=0; ii<3; ++ii)
         {
            std::istringstream riss(results[ii+1]);
            std::vector<std::string> elems((std::istream_iterator<WordDelimitedBy<'/'>>(riss)),
                                           std::istream_iterator<WordDelimitedBy<'/'>>());

            if (elems.size() > 0)
               vertexIndex[ii] = strtol(elems[0].c_str(), NULL, 10);
            if (elems.size() > 1)
               uvIndex[ii] = strtol(elems[1].c_str(), NULL, 10);
            if (elems.size() == 3)
               normalIndex[ii] = strtol(elems[2].c_str(), NULL, 10);
            if ((elems.size() == 0) || (elems.size() > 3))
            {
               std::cerr << "Error parsing face: '" << in_line << std::endl;
            }
         }
         for (int ii=0; ii<3; ++ii)
         {
            vertexIndices.push_back(vertexIndex[ii]);
            uvIndices.push_back(uvIndex[ii]);
            normalIndices.push_back(normalIndex[ii]);
         }
      }
   }
   file.close();

   m_vertices = temp_vertices;
   for (int32_t ii = 0; ii < vertexIndices.size();)
   {
      Triangle tri;

      // vertices
      for (int32_t jj = 0; jj < 3; ++jj)
         tri.m_v[jj] = vertexIndices[ii + jj] - 1;

      // normals
      if (triplet_valid(normalIndices, ii))
         for (int32_t jj = 0; jj < 3; ++jj)
            tri.m_n[jj] = normalIndices[ii + jj] - 1;
      else
         for (int32_t jj = 0; jj < 3; ++jj)
            tri.m_n[jj] = NO_ENTRY;

      // UVs
      if (triplet_valid(uvIndices, ii))
         for (int32_t jj = 0; jj < 3; ++jj)
            tri.m_uv[jj] = uvIndices[ii + jj] - 1;
      else
         for (int32_t jj = 0; jj < 3; ++jj)
            tri.m_uv[jj] = NO_ENTRY;

      m_faces.push_back(tri);
      ii += 3;
   }

   return m_vertices.size();
}

std::ostream &operator<<( std::ostream &output, const TriMesh &M )
{
    const auto& vertices = M.vertices();
    const auto& faces = M.faces();

    output << "# libcsg trimesh" << std::endl
           << "# " << vertices.size() << " vertices" << std::endl;
    for (const auto& v : vertices)
    {
        output << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    output << std::endl << "# " << faces.size() << " triangles" << std::endl;
    for (const auto& f : faces)
    {
        output << "f " << f.m_v[0]+1 << " " << f.m_v[1]+1 << " " << f.m_v[2]+1 << std::endl;
    }

    return output;
}


AABBTree TriMesh::createAABBTree() const
{
    AABBTree tree(0.05, m_faces.size());
    for (uint32_t idx=0; idx<m_faces.size(); ++idx)
    {
        Eigen::Vector3d lower;
        Eigen::Vector3d upper;
        for (uint32_t ii=0; ii<3; ++ii)
        {
            lower[ii] = std::min( std::min(m_vertices[m_faces[idx].m_v[0]][ii],
                                           m_vertices[m_faces[idx].m_v[1]][ii]),
                                  m_vertices[m_faces[idx].m_v[2]][ii]);

            upper[ii] = std::max( std::max(m_vertices[m_faces[idx].m_v[0]][ii],
                                           m_vertices[m_faces[idx].m_v[1]][ii]),
                                  m_vertices[m_faces[idx].m_v[2]][ii]);
        }

        tree.insertAABB(idx, lower, upper);
    }

    return tree;
}


Eigen::Vector3d TriMesh::faceNormalImpl(uint32_t f) const
{
    // TODO: lazily evaluate face normals and store them

    const Triangle& face = faces()[f];
    const Eigen::Vector3d p01 = (m_vertices[face.m_v[1]] - m_vertices[face.m_v[0]]).normalized();
    const Eigen::Vector3d p12 = (m_vertices[face.m_v[2]] - m_vertices[face.m_v[1]]).normalized();

    return p01.cross(p12);
}


// Non-Eigen accessors for Python wrapping
std::vector<double> TriMesh::mesh_vertices() const
{
   std::vector<double> out(m_vertices.size()*3);
   for (uint32_t ii=0; ii<m_vertices.size(); ++ii)
      for (uint32_t jj=0; jj<3; ++jj)
         out[ii*3+jj] = m_vertices[ii][jj];
   return out;
}

std::vector<double> TriMesh::mesh_normals() const
{
   std::vector<double> out(m_normals.size()*3);
   for (uint32_t ii=0; ii<m_normals.size(); ++ii)
      for (uint32_t jj=0; jj<3; ++jj)
         out[ii*3+jj] = m_normals[ii][jj];
   return out;
}

std::vector<double> TriMesh::mesh_uvs() const
{
   std::vector<double> out(m_uvs.size()*2);
   for (uint32_t ii=0; ii<m_uvs.size(); ++ii)
      for (uint32_t jj=0; jj<2; ++jj)
         out[ii*2+jj] = m_uvs[ii][jj];
   return out;
}

std::vector<unsigned int> TriMesh::mesh_faces() const
{
   const size_t num_faces = m_faces.size();
   std::vector<unsigned int> out(num_faces*3);
   for (uint32_t ii=0; ii<num_faces; ++ii)
      for (uint32_t jj=0; jj<3; ++jj)
         out[ii*3+jj] = m_faces[ii].m_v[jj];
   return out;
}

#ifndef SWIG
Eigen::Vector3d TriMesh::faceNormal(uint32_t f) const
{
    return faceNormalImpl(f);
}
#endif

}  // namespace CSG
