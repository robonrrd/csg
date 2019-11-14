# libCSG - A constructive "solid" geometry library #

This is a clean re-implementation of work I did years ago at a visual effects company.  This library implements the three core CSG operations (intersection, union, and difference) on triangulated meshes.  Unlike typical CSG implementations, the meshes to 
be operated on do not need to be solid (*i.e.* watertight) or even manifold.

## Dependencies ##
Requires [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download).


## To Build ##
In the top-level ```csg``` directory:
```
mkdir build
cd build
cmake ..
make
```

## Explanation of the Algorithm ##
Essentially, the CSG algorithm finds all intersections of the triangles in mesh "A" and
mesh "B", retriangulated the intersected faces to preserve these new edges, and then
seperates the cut mesh (what we called the "clay") into two new meshes.

In more detail:
1. Create AABB tree of the two meshes (the "clay" and the "knife").
2. Using the AABB tree, do a broadphase intersection step, to find pairs of
bounding boxes that may contain intersecting triangles
3. Determine exact triangle-triangle intersections, yeilding potentially new
points that describe new edges.
4. Retriangulate the cut triangles with any added points and edges from step 3.
5. Categorize the cut triangles into two new surfaces, based on whether they
are above or below the face that cut them. We flood-fill the membership among
the uncut triangles.

