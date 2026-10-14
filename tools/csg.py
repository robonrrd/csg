#!/usr/bin/env python3
"""
Run a CSG operation on two OBJ meshes and show the two resulting pieces in an
interactive Open3D window, each piece in its own colour.

--clay and --knife are OBJ files

--operation is accepted but currently unused (the library is always asked for
      kDifference)


--cap is a boolean (defaults to false) that controls whether the cut clay meshes are
      capped with the cut knife surface


The spacebar moves the two pieces apart and back again on a second press.
"""

import argparse
import os
import sys

import numpy as np
import open3d as o3d
import pycsg

COLOR_A = (0.85, 0.25, 0.20)   # red-ish
COLOR_B = (0.20, 0.45, 0.85)   # blue-ish

def to_trimesh(path):
    m = o3d.io.read_triangle_mesh(path)   # polygons are triangulated on load
    v = np.asarray(m.vertices, dtype=float)
    f = np.asarray(m.triangles)
    if len(v) == 0 or len(f) == 0:
        sys.exit('%s: no geometry found' % path)
    return pycsg.TriMesh(v.ravel().tolist(), f.ravel().astype(int).tolist())


def to_o3d(mesh, color):
    v = np.array(mesh.mesh_vertices(), dtype=float).reshape(-1, 3)
    f = np.array(mesh.mesh_faces(), dtype=np.int32).reshape(-1, 3)
    m = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(v),
                                  o3d.utility.Vector3iVector(f))
    m.compute_vertex_normals()
    m.paint_uniform_color(color)
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--clay', help='OBJ file used as the clay')
    ap.add_argument('--knife', help='OBJ file used as the knife')
    ap.add_argument('--operation', nargs='?', default='difference',
                    help='unused, reserved for future use')
    ap.add_argument('--cap', action="store_true", help='cap the cut')

    args = ap.parse_args()

    clay, knife = to_trimesh(args.clay), to_trimesh(args.knife)
    A, B = pycsg.TriMesh(), pycsg.TriMesh()
    pycsg.CSGEngine(clay, knife).construct(pycsg.kDifference, args.cap, A, B)
    print('A: %d vertices, %d faces' % (len(A.mesh_vertices()) // 3, len(A.mesh_faces()) // 3))
    print('B: %d vertices, %d faces' % (len(B.mesh_vertices()) // 3, len(B.mesh_faces()) // 3))

    geoms = [to_o3d(m, c) for m, c in ((A, COLOR_A), (B, COLOR_B))
             if len(m.mesh_faces())]
    if not geoms:
        sys.exit('both results are empty')

    # Spacebar toggles an "exploded" view: the two pieces move apart along the line
    # between their centroids by half the largest extent of the original clay.
    # With only one piece there is nothing to separate, so the key does nothing.
    cv = np.array(clay.mesh_vertices(), dtype=float).reshape(-1, 3)
    distance = 0.5 * float((cv.max(axis=0) - cv.min(axis=0)).max())
    offsets = []
    if len(geoms) == 2:
        axis = geoms[0].get_center() - geoms[1].get_center()
        norm = np.linalg.norm(axis)
        axis = axis / norm if norm > 1e-12 else np.array([1.0, 0.0, 0.0])
        offsets = [axis * distance, -axis * distance]
    exploded = [False]

    def toggle_explode(vis):
        if not offsets:
            return False
        sign = -1.0 if exploded[0] else 1.0
        for g, off in zip(geoms, offsets):
            g.translate(sign * off, relative=True)
            vis.update_geometry(g)
        exploded[0] = not exploded[0]
        return False

    vis = o3d.visualization.VisualizerWithKeyCallback()
    vis.create_window(window_name='CSG: A red, B blue (space: separate pieces)')
    for g in geoms:
        vis.add_geometry(g)
    vis.get_render_option().mesh_show_back_face = True
    vis.register_key_callback(ord(' '), toggle_explode)
    vis.run()
    vis.destroy_window()

if __name__ == '__main__':
    main()
