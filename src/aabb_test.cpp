//
//

#include <iostream>
#include <cstdlib>
#include "trimesh.h"
#include "aabb.h"


void
constructionTest(const CSG::TriMesh& input)
{
    std::cout << "Construction test:" << std::endl;

    std::cout << " Building AABB tree" << std::endl;
    AABBTree tree = input.createAABBTree();
    std::cout << "  ..done" << std::endl;

    tree.validate();

    std::cout << " numObjects = " << tree.numObjects() << "  "
              << " height = " << tree.getHeight() << "  "
              << " node count = " << tree.getNodeCount() << "  "
              << std::endl;

    // for (size_t ii=0; ii<tree.numObjects(); ++ii)
    // {
    //     auto aabb = tree.getAABB(ii);
    //     std::cout << "AABB " << ii << "  "
    //               << "c: " << aabb.center().transpose() << "  "
    //               << "l: " << aabb.lowerBound.transpose() << "  "
    //               << "u: " << aabb.upperBound.transpose() << "  "
    //               << "s: " << aabb.surfaceArea() << "  "
    //               << std::endl;
    // }
    std::cout << std::endl;
}


void
pointWeldingTest(uint32_t nPoints, uint32_t nDuplicates)
{
    std::cout << "Point welding test:" << std::endl;
    const float extent = 10.0f;
    const double radius = 1e-06;

    std::vector<Eigen::Vector3d> duplicate_points;
    AABBTree tree(0.05, nPoints);
    for (uint32_t ii=0; ii<(nPoints-nDuplicates); ++ii)
    {
        Eigen::Vector3d pt(extent*drand48(), extent*drand48(), extent*drand48());
        tree.addSphere(ii, pt, radius);

        if (ii < nDuplicates)
        {
            duplicate_points.push_back(pt);
        }
    }
    for (uint32_t ii=0; ii<nDuplicates; ++ii)
    {
        Eigen::Vector3d pt = duplicate_points[ii];
        double rr = radius*0.5;
        pt += Eigen::Vector3d(rr*drand48(), 0,0);//r*drand48(), rr*drand48());
        tree.addSphere((nPoints-nDuplicates)+ii, duplicate_points[ii], radius);
    }

    tree.validate();

    std::cout << "numObjects = " << tree.numObjects() << "  "
             << "height = " << tree.getHeight() << "  "
             << "node count = " << tree.getNodeCount() << "  "
             << std::endl;

    std::cout << std::endl;

    std::cout << "No duplicates." << std::endl;
    std::vector<uint32_t> qq = tree.query(200);
    for (auto& q: qq)
        std::cout << " " << q;
    std::cout << std::endl;

    std::cout << "One duplicate." << std::endl;
    qq = tree.query(991);
    for (auto& q: qq)
        std::cout << " " << q;
    std::cout << std::endl;
}


int
main(int argc, char *argv[])
{
   using namespace CSG;

   if (argc != 2)
   {
       std::cout << "Usage: aabb_test OBJ_FILE" << std::endl;
       return 0;
   }

   std::cout << "Reading " << argv[1] << ".." << std::endl;
   TriMesh mesh;
   const std::string path(argv[1]);
   mesh.loadOBJ( path );
   std::cout << " ..done" << std::endl;

   constructionTest(mesh);

   pointWeldingTest( 1000, 10 );
}
