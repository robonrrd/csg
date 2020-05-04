//
//

#include <fstream>
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


   CSG::CSGEngine engine(clay, knife);

   CSGOperation operation = CSGOperation::kDifference;
   bool cap = true; // create end-caps for solid output?
   TriMesh A, B;
   engine.construct(operation, cap, A, B);

   //std::cout << "output A:" << std::endl
   //          << A
   //          << std::endl;
   //std::cout << "output B:" << std::endl
   //          << B
   //          << std::endl;

   std::ofstream outfile;
   outfile.open("A.obj");
   outfile << A << std::endl;
   outfile.close();

   outfile.open("B.obj");
   outfile << B << std::endl;
   outfile.close();
}
