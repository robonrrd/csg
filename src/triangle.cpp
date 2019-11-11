#include "libcsg.h"

#include <string>


namespace CSG
{

// Triangle
std::ostream& operator<<(std::ostream& os, const Triangle& tri)
{
   os << "v: " << tri.m_v[0] << ", " << tri.m_v[1] << ", " << tri.m_v[2] << "  "
      << "uv: " << tri.m_uv[0] << ", " << tri.m_uv[1] << ", " << tri.m_uv[2]  << "  "
      << "n: "  << tri.m_n[0] << ", " << tri.m_n[1] << ", " << tri.m_n[2];

   return os;
}

} // namespace CSG
