//
//

#include <iostream>
#include <stdlib.h>
#include "libcsg.h"

int
main(int argc, char *argv[])
{
   using namespace CSG;


   if (argc != 2)
   {
       std::cout << "Usage: intersection_test RANDOM_SEED" << std::endl;
       return 0;
   }

   // intersecting seeds: 14, 17

   uint32_t seed = std::stol(argv[1]);
   std::cout << "Starting random seed: " << seed << std::endl;
   while (true)
   {
       srand48(seed++);

       std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertices;
       for (uint32_t ii=0; ii<6; ++ii)
           vertices.push_back(Eigen::Vector3d(drand48(), drand48(), drand48()));

       Triangle a, b;
       for (uint32_t ii=0; ii<3; ++ii)
       {
           a.m_v[ii] = ii;
           b.m_v[ii] = ii+3;
       }

       TriangleIntersection ix = a.intersect(b, vertices);
       std::cout << "  " << ix << std::endl;
       if (ix.intersect)
       {
           std::cout << "# seed " << seed << std::endl;
           // dump the two triangles as an obj file
           for (uint32_t ii=0; ii<6; ++ii)
           {
               std::cout << "v " << vertices[ii][0] << " "<< vertices[ii][1] << " "<< vertices[ii][2]
                         << std::endl;
           }
           std::cout << "# intersection points" << std::endl
                     << "v " << ix.p0[0] << " " << ix.p0[1] << " " << ix.p0[2] << std::endl
                     << "v " << ix.p1[0] << " " << ix.p1[1] << " " << ix.p1[2] << std::endl
                     << std::endl;

           std::cout << "f 1 2 3" << std::endl
                     << "f 4 5 6" << std::endl;

           break;
       }
   }
}
