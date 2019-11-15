//
//

#include <iostream>
#include "libcsg.h"

int
main(int argc, char *argv[])
{
   using namespace CSG;


   if (argc != 3)
   {
       std::cout << "Usage: csg_test CLAY_OBJ_FILE  KNIFE_OBJ_FILE" << std::endl;
       return 0;
   }

   TriMesh clay;
   {
       const std::string path(argv[1]);
       clay.loadOBJ( path );
   }
   TriMesh knife;
   {
       const std::string path(argv[2]);
       knife.loadOBJ( path );
   }

   CSGOperation operation = CSGOperation::kDifference;
   TriMesh A, B;
   CSG::CSG(clay, knife, operation, A, B);
}
