"""Tests of `ddg_tools.simplify`, `build_hierarchy` and `qem_simplify` (surface simplification
with intrinsic error metrics, crate ddg-intrinsic-simplification), without Blender.

Run through `scripts/py-test.sh` like test_ddg_tools.py. Standard library only; numpy and scipy
are used when importable.
"""

import array
import math
import sys
import unittest

import ddg_tools


def icosphere(level=2):
    """A unit icosphere: (positions, triangles) as array.array."""
    t = (1.0 + math.sqrt(5.0)) / 2.0
    pts = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t),
           (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    pts = [tuple(c / math.sqrt(x * x + y * y + z * z) for c in (x, y, z)) for x, y, z in pts]
    tris = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4),
            (11, 10, 2), (10, 7, 6), (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8),
            (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    for _ in range(level):
        cache, new = {}, []

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                m = [(pts[a][i] + pts[b][i]) / 2 for i in range(3)]
                n = math.sqrt(sum(c * c for c in m))
                pts.append(tuple(c / n for c in m))
                cache[key] = len(pts) - 1
            return cache[key]

        for a, b, c in tris:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        tris = new
    return (array.array('d', [c for p in pts for c in p]),
            array.array('i', [v for t in tris for v in t]))


def bumpy_grid(n=14):
    """An open grid with a bump (a disk topologically): (positions, triangles)."""
    positions, triangles = [], []
    for j in range(n + 1):
        for i in range(n + 1):
            x, y = i / n, j / n
            positions += [x, y, 0.2 * math.sin(3 * x) * math.cos(2 * y)]
    for j in range(n):
        for i in range(n):
            a, b = j * (n + 1) + i, j * (n + 1) + i + 1
            c, d = a + n + 1, b + n + 1
            triangles += [a, b, d, a, d, c]
    return array.array('d', positions), array.array('i', triangles)


def centroid_points(triangles):
    nf = len(triangles) // 3
    return array.array('i', range(nf)), array.array('d', [1 / 3] * (3 * nf))


class SimplifyTest(unittest.TestCase):
    def test_sphere_ratio_and_arrays(self):
        positions, triangles = icosphere(3)
        n = len(positions) // 3
        s = ddg_tools.simplify(positions, triangles, ratio=0.2)
        self.assertEqual(s.n_input_vertices, n)
        self.assertLessEqual(s.n_vertices, round(0.2 * n) + 1)
        self.assertEqual(s.report["stop"], "targetreached")
        self.assertEqual(len(s.vertices), s.n_vertices)
        self.assertEqual(len(s.faces), 3 * s.n_faces)
        self.assertEqual(len(s.faces_compact), 3 * s.n_faces)
        self.assertEqual(len(s.lengths), 3 * s.n_faces)
        # A closed genus-0 mesh with V vertices has 2V - 4 faces.
        self.assertEqual(s.n_faces, 2 * s.n_vertices - 4)
        verts = list(s.vertices)
        for k, f in enumerate(s.faces):
            self.assertEqual(f, verts[s.faces_compact[k]])
        self.assertEqual(len(s.point_face), n)
        self.assertEqual(len(s.point_bary), 3 * n)
        self.assertTrue(all(0 <= f < s.n_faces for f in s.point_face))
        for p in range(n):
            self.assertAlmostEqual(sum(s.point_bary[3 * p:3 * p + 3]), 1.0, places=9)
        self.assertGreaterEqual(min(s.lengths), 0.0)

    def test_target_and_errors(self):
        positions, triangles = icosphere(2)
        s = ddg_tools.simplify(positions, triangles, target=30)
        self.assertEqual(s.n_vertices, 30)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^target or ratio"):
            ddg_tools.simplify(positions, triangles)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^target and ratio"):
            ddg_tools.simplify(positions, triangles, target=10, ratio=0.5)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^ratio"):
            ddg_tools.simplify(positions, triangles, ratio=2.0)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^weights"):
            ddg_tools.simplify(positions, triangles, ratio=0.5, weights=1.5)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^representation"):
            ddg_tools.simplify(positions, triangles, ratio=0.5, representation="floats")
        with self.assertRaisesRegex(ddg_tools.InputError, r"^masses"):
            ddg_tools.simplify(positions, triangles, ratio=0.5, masses=array.array('d', [1.0]))
        with self.assertRaises(ddg_tools.DdgError):
            ddg_tools.simplify(positions, triangles, ratio=0.5, fixed=array.array('i', [10 ** 6]))
        # Three triangles on one edge: not a manifold.
        pos = array.array('d', [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1])
        tri = array.array('i', [0, 1, 2, 0, 1, 3, 0, 1, 4])
        with self.assertRaises(ddg_tools.MeshError):
            ddg_tools.simplify(pos, tri, ratio=0.5)

    def test_determinism_and_memory(self):
        positions, triangles = icosphere(3)
        a = ddg_tools.simplify(positions, triangles, ratio=0.15)
        b = ddg_tools.simplify(positions, triangles, ratio=0.15)
        self.assertEqual(list(a.vertices), list(b.vertices))
        self.assertEqual(list(a.faces), list(b.faces))
        self.assertEqual(list(a.face_colours), list(b.face_colours))
        c = ddg_tools.simplify(positions, triangles, ratio=0.15, memory=False)
        self.assertEqual(c.n_vertices, a.n_vertices)
        d = ddg_tools.simplify(positions, triangles, ratio=0.15, weights=1.0)
        self.assertEqual(d.n_vertices, a.n_vertices)

    def test_projective_sign(self):
        positions, triangles = icosphere(3)
        faces, bary = centroid_points(triangles)
        a = ddg_tools.simplify(positions, triangles, ratio=0.1, points=(faces, bary))
        c = ddg_tools.simplify(positions, triangles, ratio=0.1, points=(faces, bary),
                               projective_sign="negative")
        p = ddg_tools.simplify(positions, triangles, ratio=0.1, points=(faces, bary),
                               projective_sign="positive")
        # The default is "negative"; the coarse mesh does not depend on the sign, the map does.
        self.assertEqual(list(a.point_bary), list(c.point_bary))
        self.assertEqual(list(c.vertices), list(p.vertices))
        self.assertEqual(list(c.faces), list(p.faces))
        self.assertNotEqual(list(c.point_bary), list(p.point_bary))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.simplify(positions, triangles, ratio=0.1, projective_sign="plus")

    def test_fixed_vertices_stay(self):
        positions, triangles = icosphere(3)
        fixed = array.array('i', [0, 3, 40, 100])
        s = ddg_tools.simplify(positions, triangles, ratio=0.1, fixed=fixed)
        for v in fixed:
            self.assertIn(v, list(s.vertices))

    def test_extra_points_and_colouring(self):
        positions, triangles = icosphere(3)
        n = len(positions) // 3
        faces, bary = centroid_points(triangles)
        s = ddg_tools.simplify(positions, triangles, ratio=0.1, points=(faces, bary))
        total = n + len(faces)
        self.assertEqual(len(s.point_face), total)
        self.assertTrue(all(f >= 0 for f in s.point_face))
        self.assertEqual(len(s.point_colours), total)
        colours = list(s.face_colours)
        self.assertEqual(len(colours), s.n_faces)
        self.assertEqual(max(colours) + 1, s.n_colours)
        for p in range(total):
            self.assertEqual(s.point_colours[p], colours[s.point_face[p]])
        # Faces that share a vertex differ in colour.
        fc = s.faces_compact
        at = {}
        for f in range(s.n_faces):
            for v in set(fc[3 * f:3 * f + 3]):
                at.setdefault(v, []).append(f)
        for fs in at.values():
            self.assertEqual(len({colours[f] for f in fs}), len(fs))
        # Every coarse face is the face of some centroid on a closed mesh with many fine faces.
        used = {s.point_face[n + f] for f in range(len(faces))}
        self.assertGreater(len(used), s.n_faces // 2)

    def test_point_errors(self):
        positions, triangles = icosphere(1)
        with self.assertRaises(ddg_tools.DdgError):
            ddg_tools.simplify(positions, triangles, ratio=0.5,
                               points=(array.array('i', [10 ** 5]), array.array('d', [.3, .3, .4])))
        with self.assertRaises(ddg_tools.DdgError):
            ddg_tools.simplify(positions, triangles, ratio=0.5,
                               points=(array.array('i', [0]), array.array('d', [.5, .5, .5])))
        with self.assertRaises(TypeError):
            ddg_tools.simplify(positions, triangles, ratio=0.5, points=array.array('i', [0]))

    def test_prolongation(self):
        positions, triangles = icosphere(3)
        n = len(positions) // 3
        s = ddg_tools.simplify(positions, triangles, ratio=0.2)
        rows, cols, values = s.prolongation()
        self.assertEqual(rows.typecode, 'i')
        self.assertEqual(values.typecode, 'd')
        sums = [0.0] * n
        for r, v in zip(rows, values):
            sums[r] += v
        self.assertTrue(all(abs(x - 1.0) < 1e-9 for x in sums))
        self.assertTrue(all(0 <= c < s.n_vertices for c in cols))
        # A kept vertex maps to itself.
        for k, v in enumerate(s.vertices):
            entries = [(c, val) for r, c, val in zip(rows, cols, values) if r == v]
            self.assertEqual(entries, [(k, 1.0)])
            if k > 5:
                break
        pr, pc, pv = s.point_prolongation()
        self.assertEqual(max(pr) + 1, n)

    def test_vector_prolongation(self):
        positions, triangles = icosphere(3)
        n = len(positions) // 3
        s = ddg_tools.simplify(positions, triangles, ratio=0.2)
        v = s.vector_prolongation()
        self.assertEqual(v["shape"], (n, s.n_vertices))
        self.assertEqual(len(v["rows"]), len(v["re"]))
        self.assertEqual(len(set(v["rows"])), n)  # no empty row
        self.assertEqual(v["report"]["n_default"], 0)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^vector prolongation"):
            ddg_tools.simplify(positions, triangles, ratio=0.2,
                               refine_after=True).vector_prolongation()
        with self.assertRaisesRegex(ddg_tools.InputError, r"^n"):
            s.vector_prolongation(0)

    def test_scipy_helpers(self):
        positions, triangles = icosphere(2)
        s = ddg_tools.simplify(positions, triangles, ratio=0.3)
        try:
            import scipy.sparse  # noqa: F401
            import numpy  # noqa: F401
        except ImportError:
            with self.assertRaises(ImportError):
                s.prolongation_scipy()
            return
        import numpy as np
        p = s.prolongation_scipy()
        self.assertEqual(p.shape, (len(positions) // 3, s.n_vertices))
        self.assertTrue(np.allclose(np.asarray(p.sum(axis=1)).ravel(), 1.0))
        pv = s.vector_prolongation_scipy()
        self.assertEqual(pv.dtype, np.complex128)

    def test_refinement(self):
        positions, triangles = bumpy_grid(10)
        n = len(positions) // 3
        s = ddg_tools.simplify(positions, triangles, ratio=0.3, refine_before=True,
                               refine_after=True)
        self.assertIsNotNone(s.report["pre_refinement"])
        self.assertIsNotNone(s.report["post_refinement"])
        self.assertGreaterEqual(s.report["post_refinement"]["min_angle"], math.radians(24.0))
        self.assertEqual(len(s.point_face), n)
        self.assertEqual(len(s.vertices), s.n_vertices)
        for rep in ("integer", "signposts"):
            t = ddg_tools.simplify(positions, triangles, ratio=0.3, refine_before=True,
                                   representation=rep)
            self.assertEqual(t.report["pre_refinement"]["representation"], rep)

    def test_anisotropy(self):
        positions, triangles = icosphere(3)
        s = ddg_tools.simplify(positions, triangles, ratio=0.3, anisotropy=0.2,
                               anisotropy_clamp=True)
        self.assertIsNotNone(s.report["anisotropy"])
        self.assertLessEqual(s.report["anisotropy"]["tau_used"], 0.2 + 1e-12)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^anisotropy_direction"):
            ddg_tools.simplify(positions, triangles, ratio=0.3, anisotropy=0.1,
                               anisotropy_direction="sideways")

    def test_boundary_mesh(self):
        positions, triangles = bumpy_grid(12)
        n = len(positions) // 3
        s = ddg_tools.simplify(positions, triangles, ratio=0.25)
        self.assertLessEqual(s.n_vertices, round(0.25 * n) + 1)
        self.assertGreater(s.report["n_boundary"] + s.report["n_interior"], 0)

    def test_gil_released(self):
        import threading
        positions, triangles = icosphere(4)
        results = []

        def run():
            results.append(ddg_tools.simplify(positions, triangles, ratio=0.05).n_vertices)

        threads = [threading.Thread(target=run) for _ in range(2)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0], results[1])


class ApplicationsTest(unittest.TestCase):
    def setUp(self):
        self.positions, self.triangles = icosphere(3)
        self.n = len(self.positions) // 3
        self.s = ddg_tools.simplify(self.positions, self.triangles, ratio=0.25,
                                    fixed=array.array('i', [0, 1]))

    def test_poisson_mean_zero_and_dirichlet(self):
        f = array.array('d', [self.positions[3 * v + 2] for v in range(self.n)])
        r = self.s.poisson(f)
        self.assertEqual(len(r["fine"]), self.n)
        self.assertEqual(len(r["coarse"]), self.s.n_vertices)
        self.assertEqual(r["report"]["n_fixed"], 0)
        # z is a degree-1 spherical harmonic: -Lap z = 2 z, so the solution is about z/2 (sign
        # per the crate's convention); only check it is a multiple of z.
        ratio = [r["fine"][v] / f[v] for v in range(self.n) if abs(f[v]) > 0.5]
        self.assertLess(max(ratio) - min(ratio), 0.5 * abs(sum(ratio) / len(ratio)))
        d = self.s.poisson(f, fixed=array.array('i', [0, 1]), values=array.array('d', [0.0, 1.0]))
        self.assertEqual(d["report"]["n_fixed"], 2)
        self.assertAlmostEqual(d["fine"][1], 1.0, places=9)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^fixed"):
            self.s.poisson(f, fixed=array.array('i', [10 ** 6]))
        with self.assertRaisesRegex(ddg_tools.InputError, r"^f"):
            self.s.poisson(array.array('d', [1.0]))

    def test_geodesics(self):
        exact = self.s.geodesic(0, method="exact")
        heat = self.s.geodesic(0, method="heat")
        self.assertAlmostEqual(exact["fine"][0], 0.0, places=9)
        # The antipode of a unit sphere's vertex is at distance about pi.
        far = max(range(self.n), key=lambda v: exact["fine"][v])
        self.assertAlmostEqual(exact["fine"][far], math.pi, delta=0.25)
        err = sum(abs(a - b) for a, b in zip(exact["fine"], heat["fine"])) / self.n
        self.assertLess(err, 0.4)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^method"):
            self.s.geodesic(0, method="fast")
        with self.assertRaises(ddg_tools.DdgError):
            self.s.geodesic(10 ** 6)

    def test_voronoi(self):
        v = self.s.voronoi(array.array('i', [0, 1]))
        self.assertEqual(len(v["fine_labels"]), self.n)
        self.assertEqual(set(v["fine_labels"]), {0, 1})
        self.assertEqual(v["fine_labels"][0], 0)
        self.assertEqual(v["fine_labels"][1], 1)
        with self.assertRaises(ddg_tools.DdgError):
            self.s.voronoi(array.array('i', []))


class HierarchyTest(unittest.TestCase):
    def test_levels_and_solves(self):
        positions, triangles = icosphere(4)
        n = len(positions) // 3
        h = ddg_tools.build_hierarchy(positions, triangles, ratio=0.25, n_min=60)
        self.assertGreaterEqual(h.n_levels, 3)
        self.assertEqual(h.level_sizes[0], n)
        self.assertLessEqual(h.level_sizes[-1], 60)
        self.assertEqual(h.report["counts"], h.level_sizes)
        rows, cols, values = h.prolongation(0)
        self.assertEqual(max(rows) + 1, n)
        self.assertEqual(max(cols) + 1, h.level_sizes[1])
        self.assertEqual(len(h.mass(0)), n)
        self.assertEqual(len(h.level_vertices(1)), h.level_sizes[1])
        self.assertEqual(len(h.faces(1)), 3 * (2 * h.level_sizes[1] - 4))
        lr, lc, lv = h.laplacian(0)
        self.assertEqual(max(lr) + 1, n)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^level"):
            h.prolongation(h.n_levels)

        b = array.array('d', [math.cos(0.1 * i) for i in range(n)])
        mg = h.multigrid(shift=1.0)
        self.assertEqual(mg.info["level_sizes"][0], n)
        x, report = mg.pcg(b)
        self.assertTrue(report["converged"])
        # The residual of (L + M) x = b, assembled from the COO triplets.
        mass = h.mass(0)
        ax = [mass[i] * x[i] for i in range(n)]
        for r, c, v in zip(lr, lc, lv):
            ax[r] += v * x[c]
        res = math.sqrt(sum((ax[i] - b[i]) ** 2 for i in range(n)))
        self.assertLess(res / math.sqrt(sum(t * t for t in b)), 1e-8)
        y, report2 = mg.solve(b)
        self.assertTrue(report2["converged"])
        self.assertLess(max(abs(x[i] - y[i]) for i in range(n)), 1e-6)

    def test_known_values_and_explicit_matrix(self):
        positions, triangles = icosphere(3)
        n = len(positions) // 3
        h = ddg_tools.build_hierarchy(positions, triangles, ratio=0.25, n_min=40,
                                      fixed=array.array('i', [0, 7]))
        lr, lc, lv = h.laplacian(0)
        mg = h.multigrid((lr, lc, lv), known=array.array('i', [0, 7]))
        b = array.array('d', [0.0] * n)
        x, report = mg.pcg(b, known_values=array.array('d', [0.0, 1.0]))
        self.assertTrue(report["converged"])
        self.assertAlmostEqual(x[0], 0.0, places=12)
        self.assertAlmostEqual(x[7], 1.0, places=12)
        self.assertTrue(all(-1e-9 <= t <= 1 + 1e-9 for t in x))  # discrete maximum principle
        with self.assertRaises(ddg_tools.DdgError):
            h.multigrid((lr, lc, lv), known=array.array('i', [10 ** 6]))

    def test_errors(self):
        positions, triangles = icosphere(2)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^hierarchy"):
            ddg_tools.build_hierarchy(positions, triangles, ratio=1.5)
        h = ddg_tools.build_hierarchy(positions, triangles, ratio=0.25, n_min=20)
        with self.assertRaises(ddg_tools.DdgError):
            h.multigrid(shift=-1.0)
        with self.assertRaises(TypeError):
            h.multigrid(array.array('d', [1.0]))
        for fresh in (False, True):
            g = ddg_tools.build_hierarchy(positions, triangles, ratio=0.25, n_min=20,
                                          fresh_levels=fresh)
            self.assertEqual(g.report["fresh_levels"], fresh)


class QemTest(unittest.TestCase):
    def test_sphere(self):
        positions, triangles = icosphere(3)
        n = len(positions) // 3
        q = ddg_tools.qem_simplify(positions, triangles, target=60)
        self.assertEqual(len(q.positions) // 3, 60)
        self.assertEqual(len(q.triangles) // 3, 2 * 60 - 4)
        self.assertEqual(len(q.vertex_map), n)
        self.assertEqual(len(q.projection_face), n)
        self.assertEqual(len(q.projection_bary), 3 * n)
        self.assertTrue(q.report["reached_target"])
        self.assertLess(q.report["area_error_percent"], 20.0)
        self.assertLess(q.report["max_projection_distance"], 0.2)

    def test_errors(self):
        positions, triangles = icosphere(1)
        with self.assertRaises(ddg_tools.DdgError):
            ddg_tools.qem_simplify(positions, triangles, target=10, boundary_weight=-1.0)
        with self.assertRaises(TypeError):
            ddg_tools.qem_simplify(positions, triangles)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
