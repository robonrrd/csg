#!/usr/bin/env python3
"""Fuzz the CSG engine with a folded (V-shaped or ridge-shaped) knife.

The clay is the unit cube; the knife is two large triangles meeting at a crease
along the y axis, i.e. the graph of z = f(x) with f piecewise linear.  For such
a knife the correct answer is known exactly (clip each cube face against the two
half-spaces), so each clay face's area in output A ("above" the knife) and B
("below") can be compared against an oracle.  Crease position and slopes are
drawn from small sets of "nice" numbers so that the knife regularly passes
exactly through cube vertices and edges, which is where the coincident-point
bugs live.

    PYTHONPATH=build tools/fuzz_fold.py [--n 500] [--seed 1] [--save DIR]

Needs the Python bindings (pycsg) built by CMake.
"""
import argparse, itertools, math, os, random, sys

import pycsg

CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
CUBE_FACES = [  # outward-facing triangles (same as tests/mesh_check.h makeBox)
    (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
    (3, 7, 6), (3, 6, 2), (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5)]
CUBE_QUADS = {  # plane label -> quad
    'z0': (0, 3, 2, 1), 'z1': (4, 5, 6, 7), 'y0': (0, 1, 5, 4),
    'y1': (3, 7, 6, 2), 'x0': (0, 4, 7, 3), 'x1': (1, 2, 6, 5)}


def clip(poly, plane, keep_above):
    """Clip a convex polygon against n.x > d (or <)."""
    n, d = plane
    out = []
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]
        da = sum(n[k] * a[k] for k in range(3)) - d
        db = sum(n[k] * b[k] for k in range(3)) - d
        if not keep_above:
            da, db = -da, -db
        if da >= 0:
            out.append(a)
        if (da > 0 > db) or (db > 0 > da):
            t = da / (da - db)
            out.append(tuple(a[k] + t * (b[k] - a[k]) for k in range(3)))
    return out


def poly_area(poly):
    s = [0.0, 0.0, 0.0]
    for i in range(1, len(poly) - 1):
        u = [poly[i][k] - poly[0][k] for k in range(3)]
        w = [poly[i + 1][k] - poly[0][k] for k in range(3)]
        s[0] += u[1] * w[2] - u[2] * w[1]
        s[1] += u[2] * w[0] - u[0] * w[2]
        s[2] += u[0] * w[1] - u[1] * w[0]
    return 0.5 * math.sqrt(sum(c * c for c in s))


def oracle(cx, cz, s1, s2, quad):
    """Area of the quad above / below the fold z = cz + f(x-cx)."""
    pts = [CUBE[i] for i in quad]
    p1 = ((-s1, 0.0, 1.0), cz - s1 * cx)  # z - s1*(x-cx) > cz   (plane 1, x<cx side)
    p2 = ((-s2, 0.0, 1.0), cz - s2 * cx)
    if s1 < s2:   # valley: above = both half-spaces
        a = clip(clip(pts, p1, True), p2, True) if pts else []
        above = poly_area(a) if len(a) >= 3 else 0.0
    else:         # ridge: above = either half-space
        b = clip(clip(pts, p1, False), p2, False)
        below = poly_area(b) if len(b) >= 3 else 0.0
        above = poly_area(pts) - below
    return above, poly_area(pts) - above


def make_knife(cx, cz, s1, s2, big=6.0):
    v1 = [(cx, -big, cz), (cx, big, cz), (cx - big, 0.0, cz - s1 * big)]
    v2 = [(cx, -big, cz), (cx, big, cz), (cx + big, 0.0, cz + s2 * big)]
    tris = []
    for t in (v1, v2):
        u = [t[1][k] - t[0][k] for k in range(3)]
        w = [t[2][k] - t[0][k] for k in range(3)]
        nz = u[0] * w[1] - u[1] * w[0]
        tris.append(t if nz > 0 else [t[0], t[2], t[1]])  # normals point up (+z)
    verts, faces = [], []
    for t in tris:
        base = len(verts) // 3
        for p in t:
            verts.extend(p)
        faces.extend([base, base + 1, base + 2])
    return pycsg.TriMesh(verts, faces)


def face_areas(mesh):
    v = mesh.mesh_vertices()
    f = mesh.mesh_faces()
    areas = {k: 0.0 for k in CUBE_QUADS}
    for i in range(0, len(f), 3):
        P = [tuple(v[3 * j:3 * j + 3]) for j in f[i:i + 3]]
        for lab, quad in CUBE_QUADS.items():
            axis = {'x': 0, 'y': 1, 'z': 2}[lab[0]]
            val = float(lab[1])
            if all(abs(p[axis] - val) < 1e-7 for p in P):
                areas[lab] += poly_area(P)
    return areas


def run_case(cx, cz, s1, s2):
    verts = [c for p in CUBE for c in p]
    faces = [i for t in CUBE_FACES for i in t]
    clay = pycsg.TriMesh(verts, faces)
    knife = make_knife(cx, cz, s1, s2)
    A, B = pycsg.TriMesh(), pycsg.TriMesh()
    pycsg.CSGEngine(clay, knife).construct(pycsg.kDifference, False, A, B)
    return face_areas(A), face_areas(B)


def check(cx, cz, s1, s2, tol=1e-7):
    gotA, gotB = run_case(cx, cz, s1, s2)
    errs = []
    for lab, quad in CUBE_QUADS.items():
        wantA, wantB = oracle(cx, cz, s1, s2, quad)
        if abs(gotA[lab] - wantA) > tol or abs(gotB[lab] - wantB) > tol:
            errs.append(f'{lab}: A {gotA[lab]:.4f} (want {wantA:.4f}), B {gotB[lab]:.4f} (want {wantB:.4f})')
    return errs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=300)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--exhaustive', action='store_true')
    args = ap.parse_args()

    positions = [0.0, 0.25, 0.5, 0.75, 1.0]
    heights = [0.0, 0.25, 0.5, 0.75, 1.0]
    slopes = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
    cases = [c for c in itertools.product(positions, heights, slopes, slopes) if c[2] != c[3]]
    if not args.exhaustive:
        random.Random(args.seed).shuffle(cases)
        cases = cases[:args.n]
    bad = 0
    for c in cases:
        errs = check(*c)
        if errs:
            bad += 1
            print('FAIL cx=%g cz=%g s1=%g s2=%g' % c)
            for e in errs[:3]:
                print('   ', e)
    print(f'{bad} / {len(cases)} cases failed')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
