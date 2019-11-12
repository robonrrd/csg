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

   uint32_t seed = std::stol(argv[1]);
   std::cout << "Random seed: " << seed << std::endl;
   srand48(seed);

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
   std::cout << ix << std::endl;
}
