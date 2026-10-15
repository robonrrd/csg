#!/usr/bin/env python3
"""Tests for blender/csg.py using stand-ins for Blender's bpy / bmesh modules.

Run with PYTHONPATH set to the build directory (for pycsg).

The bugs (BUG-n) are those of tests/csg_tests.cpp.  The tests marked as expected failures
document bugs that the add-on still has; when one is fixed the test reports an
"unexpected success": remove the decorator.
"""
import importlib.util
import os
import sys
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ADDON = os.path.join(HERE, '..', 'blender', 'csg.py')


class FakeOperator:
    """bpy.types.Operator: the add-on reports problems through self.report()."""
    def __init__(self):
        self.reports = []

    def report(self, kinds, message):
        self.reports.append((kinds, message))


class FakeMatrix:
    def __mul__(self, co):
        return list(co)

    __matmul__ = __mul__


class FakePoly:
    def __init__(self, index, loop_start, loop_total):
        self.index, self.loop_start, self.loop_total = index, loop_start, loop_total


class FakeLoop:
    def __init__(self, vertex_index):
        self.vertex_index = vertex_index


class FakeVertex:
    def __init__(self, co):
        self.co = co


class FakeMesh:
    def __init__(self, coords, polygons):
        """polygons: list of vertex index lists (any length)."""
        self.vertices = [FakeVertex(c) for c in coords]
        self.polygons, self.loops = [], []
        for poly in polygons:
            self.polygons.append(FakePoly(len(self.polygons), len(self.loops), len(poly)))
            self.loops.extend(FakeLoop(i) for i in poly)


class FakeObject:
    type = 'MESH'

    def __init__(self, name, coords, polygons):
        self.name = name
        self.data = FakeMesh(coords, polygons)
        self.matrix_world = FakeMatrix()


def load_addon(selected, active):
    bpy = types.ModuleType('bpy')
    bpy.types = types.SimpleNamespace(Operator=FakeOperator)
    bpy.utils = types.SimpleNamespace(register_class=lambda c: None, unregister_class=lambda c: None)
    bpy.context = types.SimpleNamespace(selected_objects=selected, active_object=active,
                                        scene=types.SimpleNamespace())
    sys.modules['bpy'] = bpy
    sys.modules['bmesh'] = types.ModuleType('bmesh')
    spec = importlib.util.spec_from_file_location('csg_addon', ADDON)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


QUAD = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]


class AddonTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_wrong_selection_is_reported(self):
        """BUG-34: this used to raise NameError (`Error` is not defined)."""
        only = FakeObject('only', QUAD, [[0, 1, 2]])
        addon = load_addon([only], only)
        op = addon.CSG()
        result = op.execute(bpy_context(addon))
        self.assertEqual(result, {'CANCELLED'})
        self.assertEqual(len(op.reports), 1)
        self.assertIn('exactly two', op.reports[0][1])

    @unittest.expectedFailure
    def test_non_triangles_are_triangulated(self):
        """BUG-34: the loops of a quad were appended as four indices, corrupting the face list."""
        quad = FakeObject('quad', QUAD, [[0, 1, 2, 3]])
        addon = load_addon([quad], quad)
        mesh = addon.build_trimesh_from_obj(quad)
        faces = list(mesh.mesh_faces())
        self.assertEqual(len(faces), 6)               # two triangles
        self.assertEqual(faces, [0, 1, 2, 0, 2, 3])   # a fan


def bpy_context(addon):
    return sys.modules['bpy'].context


if __name__ == '__main__':
    unittest.main()
