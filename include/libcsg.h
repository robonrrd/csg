// Public API for CSG library
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <string>

#include "trimesh.h"

namespace CSG
{

enum CSGOperation
{
   kBadValue = 0,
   kIntersection = 1,
   kUnion = 2,
   kDifference = 3
};

enum IParent
{
   kNone,
   kClay,
   kKnife,
   kBoth,
   kNew
};

// A reference to a point (vertex) in some mesh or list of vertices
class IPointRef
{
public:
   IPointRef() : parent(kNone), idx(0) {};
   IPointRef(IParent p, uint32_t i) : parent(p), idx(i) {};

   bool operator==(const IPointRef& other) const
   {
      return (parent == other.parent) && (idx == other.idx);
   }
   bool operator!=(const IPointRef& other) const
   {
      return !(*this == other);
   }

   IParent parent;
   uint32_t idx;  // index into parent face's mesh (if applicable)
};

// A structure to contain a reference to a point created during the intersection
// of two triangles.  The point is usually a new point, stored in 'm_newPoints'
// below, but can also be an original point from the clay or knife mesh.
class IPoint
{
 public:
   IPointRef ref;

   // Indicies of triangles that created this point
   uint32_t cidx;
   uint32_t kidx;
};

class IFace
{
public:
   IPointRef v[3];
   uint32_t  orig; // original face
};

enum TriTriIntersectionType
{
   kPointPoint,
   kPointEdge,
   kEdgeEdge
};


class TriangleIntersection
{
 public:
#ifndef SWIG
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
         os << "p0: [" << ix.p[0][0] << ", " << ix.p[0][1] << ", " << ix.p[0][2] << "]  "
            << "p1: [" << ix.p[1][0] << ", " << ix.p[1][1] << ", " << ix.p[1][2] << "]  ";
         os << "alpha:" << ix.alpha << "  beta:" << ix.beta;
      }
      return os;
   }
#endif

   bool intersect;        // is there an intersection at all?
   bool coplanar;         // are the triangles coplanar
   Eigen::Vector3d p[2];  // the two points defining the line of intersection
   double alpha;          // clay points
   double beta;           // knife points
};


class CSGEngine
{
 public:
   CSGEngine(const TriMesh& clay, const TriMesh& knife);

   void construct(CSGOperation operation, bool cap, TriMesh& output_A, TriMesh& output_B);

 private:
   // Member functions
   const Eigen::Vector3d& ipointPos(const IPointRef& ref) const;

   std::vector<IPoint> convertIntersectionToIpoints(const TriangleIntersection& ix,
                                                    uint32_t clay_face_idx,
                                                    uint32_t knife_face_idx);

   uint32_t canonicalVertexIndex(const IPointRef& ref) const;

   bool degenerateTriangle(const IFace& triangle) const;

   std::vector<IFace> retriangulate(const TriMesh& mesh, IParent which_mesh, uint32_t fidx,
                                    const std::vector<uint32_t>& new_vert_indices) const;

   std::vector<char> classifyCutFaces(const std::vector<IFace>& in_faces, IParent which);

   void classifyFaces(IParent which_surface, const std::vector<IFace>& new_faces,
                      const std::vector<bool>& is_face_cut,
                      std::vector<char>& cut_face_status,
                      std::vector<char>& uncut_face_status);

   TriMesh assembleMesh(IParent which_surface, char side,
                        const std::vector<IFace>& new_faces,
                        const std::vector<char>& cut_face_status,
                        const std::vector<char>& uncut_face_status);

   TriMesh mergeMeshes(const TriMesh& exterior, const TriMesh& cap);

   // Get a reference to the mesh specified by IParent
   const TriMesh& getMesh(IParent which_mesh) { return (which_mesh == kClay) ? m_clay : m_knife; }
   const TriMesh& otherMesh(IParent which_mesh) { return (which_mesh == kClay) ? m_knife : m_clay; }

   // Data members
   const TriMesh& m_clay;
   const TriMesh& m_knife;

   // Storage for points "created" by intersections of triangles.  The points may
   // be existing points on one of the two CSG surfaces, or (most likely) it will
   // be new.  Both result meshes will share these vertices
   std::vector<IPoint> m_newPoints;

   // The positions of the genuinely new points created by triangle intersections
   std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> m_newPointPositions;
};


}  //namespace CSG
