# libCSG - A constructive "solid" geometry library #

This is a clean re-implementation of work I did years ago at a visual effects company.  This library implements the three core CSG operations (intersection, union, and difference) on triangulated meshes.  Unlike typical CSG implementations, the meshes to 
be operated on do not need to be solid (*i.e.* watertight) or even manifold.
