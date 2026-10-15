#!/usr/bin/env python3
"""Tests for the Python bindings (pycsg). Run with PYTHONPATH set to the build directory.

The bugs (BUG-n) are those of tests/csg_tests.cpp.  A test marked as an expected failure
documents a bug that this code still has; when it is fixed, the test reports an
"unexpected success": remove the decorator.
"""
import sys
import unittest

import pycsg

CUBE_VERTS = [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
CUBE_FACES = [0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4,
              3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5]
# A big triangle in the plane z = 0.5, normal +z
PLANE_VERTS = [-5, -5, 0.5, 9, -5, 0.5, -5, 9, 0.5]
PLANE_FACES = [0, 1, 2]


def signed_volume(mesh):
    v = mesh.mesh_vertices()
    f = mesh.mesh_faces()
    total = 0.0
    for i in range(0, len(f), 3):
        a, b, c = (v[3 * j:3 * j + 3] for j in f[i:i + 3])
        total += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                  + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
    return total


class PyCsgTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_invalid_mesh_raises(self):
        """BUG-32: invalid input is an exception in Python, not a crash or silent garbage."""
        with self.assertRaises(ValueError):
            pycsg.TriMesh([0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 1, 5])      # no such vertex
        with self.assertRaises(ValueError):
            pycsg.TriMesh([0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 1])          # not a triple
        with self.assertRaises(ValueError):
            pycsg.TriMesh([0, 0, 0, 1, 0], [0, 0, 0])                   # not a triple

    def test_cut_cube(self):
        clay = pycsg.TriMesh(CUBE_VERTS, CUBE_FACES)
        knife = pycsg.TriMesh(PLANE_VERTS, PLANE_FACES)
        a, b = pycsg.TriMesh(), pycsg.TriMesh()
        pycsg.CSGEngine(clay, knife).construct(pycsg.kDifference, True, a, b)
        self.assertAlmostEqual(signed_volume(a), 0.5, places=9)
        self.assertAlmostEqual(signed_volume(b), 0.5, places=9)

    @unittest.skip("BUG-28: CSGEngine keeps references to the caller's meshes, so this "
                   "reads freed memory (the result is undefined, it may crash)")
    def test_engine_from_temporaries(self):
        """BUG-28: the engine must not depend on the caller keeping its meshes alive."""
        engine = pycsg.CSGEngine(pycsg.TriMesh(CUBE_VERTS, CUBE_FACES),
                                 pycsg.TriMesh(PLANE_VERTS, PLANE_FACES))
        junk = [pycsg.TriMesh(CUBE_VERTS, CUBE_FACES) for _ in range(50)]  # reuse freed memory
        a, b = pycsg.TriMesh(), pycsg.TriMesh()
        engine.construct(pycsg.kDifference, True, a, b)
        self.assertAlmostEqual(signed_volume(a), 0.5, places=9)
        self.assertAlmostEqual(signed_volume(b), 0.5, places=9)
        self.assertEqual(len(junk), 50)


if __name__ == '__main__':
    unittest.main()
