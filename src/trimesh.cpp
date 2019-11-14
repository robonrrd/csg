#include "libcsg.h"

#include <fstream>
#include <iterator>
#include <string>
#include <stdlib.h>


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
std::ostream& operator<<(std::ostream& os, const Triangle& tri)
{
   os << "v: " << tri.m_v[0] << ", " << tri.m_v[1] << ", " << tri.m_v[2] << "  "
      << "uv: " << tri.m_uv[0] << ", " << tri.m_uv[1] << ", " << tri.m_uv[2]  << "  "
      << "n: "  << tri.m_n[0] << ", " << tri.m_n[1] << ", " << tri.m_n[2];

   return os;
}


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
            tri.m_n[jj] = Triangle::NO_ENTRY;

      // UVs
      if (triplet_valid(uvIndices, ii))
         for (int32_t jj = 0; jj < 3; ++jj)
            tri.m_uv[jj] = uvIndices[ii + jj] - 1;
      else
         for (int32_t jj = 0; jj < 3; ++jj)
            tri.m_uv[jj] = Triangle::NO_ENTRY;

      m_faces.push_back(tri);
      ii += 3;
   }

   return m_vertices.size();
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


}  // namespace CSG
