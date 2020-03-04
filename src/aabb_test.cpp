//
//

#include <iostream>
#include "trimesh.h"
#include "aabb.h"

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
   TriMesh input;
   const std::string path(argv[1]);
   input.loadOBJ( path );
   std::cout << " ..done" << std::endl;

   std::cout << "Building AABB tree" << std::endl;
   AABBTree tree = input.createAABBTree();
   std::cout << " ..done" << std::endl;

   tree.validate();

   //std::vector<std::pair<uint32_t, uint32_t>> ix = tree.intersect(tree);
   std::cout << "numObjects = " << tree.numObjects() << "  "
             << "height = " << tree.getHeight() << "  "
             << "node count = " << tree.getNodeCount() << "  "
             << std::endl;


   for (size_t ii=0; ii<tree.numObjects(); ++ii)
   {
       auto aabb = tree.getAABB(ii);
       std::cout << "AABB " << ii << "  "
                 << "c: " << aabb.center().transpose() << "  "
                 << "l: " << aabb.lowerBound.transpose() << "  "
                 << "u: " << aabb.upperBound.transpose() << "  "
                 << "s: " << aabb.surfaceArea() << "  "
                 << std::endl;
   }
}
