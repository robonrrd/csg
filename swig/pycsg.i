%module pycsg

%include "std_vector.i"

%{
#include "Eigen/Core"
#include "libcsg.h"
extern "C"
{
#include "triangle.h"
}

%}


namespace std {
  %template(DoubleVector) vector<double>;
  %template(UnsignedIntVector) vector<unsigned int>;
}

// triangle.h is intentionally NOT %include-d here: it is an internal
// implementation detail (Shewchuk's triangulator) and is not part of the
// public Python API.  It is included above in the %{ %} block so that the
// compiled wrapper code can see it, but SWIG does not parse it.
%include "trimesh.h"
%include "aabb.h"
%include "libcsg.h"

