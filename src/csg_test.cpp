//
//

#include <iostream>
#include "libcsg.h"

int
main(int argc, char *argv[])
{
   using namespace CSG;


   if (argc != 2)
   {
       std::cout << "Usage: csg_test OBJ_FILE" << std::endl;
       return 0;
   }

   TriMesh input;
   const std::string path(argv[1]);
   input.loadOBJ( path );

   std::cout << "Building AABB tree" << std::endl;
   AABBTree tree = input.createAABBTree();
   std::cout << " ..done" << std::endl;
}
