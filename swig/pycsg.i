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

%include "triangle.h"
%include "trimesh.h"
%include "aabb.h"
%include "libcsg.h"

