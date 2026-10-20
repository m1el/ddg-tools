"""Tests of the developability functions of the built wheel of `ddg_tools`, without Blender.

Run through `scripts/py-test.sh` (see test_ddg_tools.py). Standard library only.
"""

import array
import math
import unittest

import ddg_tools


def grid(n, bump=0.0):
    """An n x n grid of quads in z = 0 facing +z, with a Gaussian bump of height `bump`."""
    positions = []
    for j in range(n + 1):
        for i in range(n + 1):
            x, y = i / n - 0.5, j / n - 0.5
            positions += [i / n, j / n, bump * math.exp(-8 * (x * x + y * y))]
    v = lambda i, j: j * (n + 1) + i
    corner_vertices, starts, triangles = [], [], []
    for j in range(n):
        for i in range(n):
            s = len(corner_vertices)
            starts.append(s)
            corner_vertices += [v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)]
            triangles += [s, s + 1, s + 2, s, s + 2, s + 3]
    return positions, corner_vertices, starts, triangles


def cylinder(m=12, h=4):
    positions = []
    for j in range(h + 1):
        for i in range(m):
            a = 2 * math.pi * i / m
            positions += [math.cos(a), math.sin(a), 0.3 * j]
    corner_vertices, starts, triangles = [], [], []
    for j in range(h):
        for i in range(m):
            s = len(corner_vertices)
            starts.append(s)
            i1 = (i + 1) % m
            corner_vertices += [j * m + i, j * m + i1, (j + 1) * m + i1, (j + 1) * m + i]
            triangles += [s, s + 1, s + 2, s, s + 2, s + 3]
    return positions, corner_vertices, starts, triangles


def cube():
    positions = []
    for i in range(8):
        positions += [float(i & 1), float((i >> 1) & 1), float((i >> 2) & 1)]
    quads = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    cv = [v for q in quads for v in q]
    starts = [4 * p for p in range(6)]
    tri = []
    for s in starts:
        tri += [s, s + 1, s + 2, s, s + 2, s + 3]
    return positions, cv, starts, tri


def arrays(mesh):
    p, cv, st, tr = mesh
    return (array.array('f', p), array.array('i', cv), array.array('i', st), array.array('i', tr))


def three_fins():
    return ([0, 0, 0, 1, 0, 0, 0.5, 1, 0, 0.5, -1, 0, 0.5, 0, 1], [0, 1, 2, 1, 0, 3, 0, 1, 4],
            [0, 3, 6], list(range(9)))


class DevelopabilityTest(unittest.TestCase):
    def test_flat_grid_has_no_energy_and_a_bump_has(self):
        flat = ddg_tools.developability(*arrays(grid(6)))
        self.assertIsInstance(flat.vertex_energy, array.array)
        self.assertEqual(len(flat.vertex_energy), 49)
        self.assertEqual(len(flat.counted), 49)
        self.assertLess(flat.energy, 1e-12)
        bumpy = ddg_tools.developability(*arrays(grid(8, 0.2)))
        self.assertGreater(bumpy.energy, 1e-4)
        self.assertAlmostEqual(bumpy.energy, sum(e for e, c in zip(bumpy.vertex_energy, bumpy.counted) if c))
        self.assertEqual(bumpy.report["n_counted"], sum(bumpy.counted))

    def test_every_energy_works_and_agrees_with_zero_on_flat(self):
        self.assertEqual(len(ddg_tools.ENERGIES), 6)
        for name in ddg_tools.ENERGIES:
            r = ddg_tools.developability(*arrays(grid(5, 0.2)), energy=name)
            self.assertTrue(all(math.isfinite(e) and e >= -1e-12 for e in r.vertex_energy), name)
            self.assertGreater(r.energy, 0, name)

    def test_options_change_the_domain(self):
        m = arrays(grid(6, 0.2))
        a = ddg_tools.developability(*m)
        b = ddg_tools.developability(*m, include_boundary=False)
        self.assertLess(sum(b.counted), sum(a.counted))
        self.assertLessEqual(b.energy, a.energy)
        tetra = ([1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1], [0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2],
                 [0, 3, 6, 9], list(range(12)))
        c = ddg_tools.developability(*arrays(tetra), omit_valence_three=False)
        d = ddg_tools.developability(*arrays(tetra), omit_valence_three=True)
        self.assertGreater(sum(c.counted), sum(d.counted))

    def test_deterministic_and_independent_of_threads(self):
        m = arrays(grid(10, 0.3))
        a = ddg_tools.developability(*m, threads=1)
        b = ddg_tools.developability(*m, threads=4)
        self.assertEqual(a.vertex_energy, b.vertex_energy)
        self.assertEqual(a.energy, b.energy)

    def test_errors(self):
        with self.assertRaises(ddg_tools.MeshError) as c:
            ddg_tools.developability(*arrays(three_fins()))
        self.assertIn("not manifold", str(c.exception))
        self.assertIn("(0, 1)", str(c.exception))
        with self.assertRaises(ddg_tools.InputError) as c:
            ddg_tools.developability(*arrays(cube()), energy="bogus")
        self.assertTrue(str(c.exception).startswith("energy"))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.developability(*arrays(cube()), threads=0)
        with self.assertRaises(TypeError):
            ddg_tools.developability(3, 4, 5, 6)
        p, cv, st, tr = grid(2)
        cv = list(cv)
        cv[4:8] = reversed(cv[4:8])
        with self.assertRaises(ddg_tools.MeshError) as c:
            ddg_tools.developability(*arrays((p, cv, st, tr)))
        self.assertIn("consistently oriented", str(c.exception))
        # A degenerate triangle names its polygon.
        p, cv, st, tr = grid(2)
        p = list(p)
        p[3:6] = p[0:3]
        with self.assertRaises(ddg_tools.DdgError):
            ddg_tools.developability(*arrays((p, cv, st, tr)))


class MakeDevelopableTest(unittest.TestCase):
    def test_unchanged_topology_keeps_polygons(self):
        mesh = arrays(grid(8, 0.2))
        r = ddg_tools.make_developable(*mesh, subdivisions=0, iterations=60, remesh=False)
        self.assertIsNone(r.faces)
        self.assertIsNone(r.face_polygon)
        self.assertEqual(len(r.positions), 81 * 3)
        self.assertEqual(list(r.vertex_source), list(range(81)))
        rep = r.report
        self.assertFalse(rep["topology_changed"])
        self.assertLess(rep["energy_final"], 0.5 * rep["energy_initial"])
        self.assertEqual(len(rep["rounds"]), 1)
        self.assertIsNotNone(rep["nonplanarity_after"])
        self.assertIsNotNone(rep["nonplanarity_polygon"])
        self.assertLess(rep["nonplanarity_after_relative"], 1.0)
        self.assertNotEqual(list(r.positions), [float(x) for x in mesh[0]])
        # The input is untouched.
        self.assertEqual(mesh[0][2], mesh[0][2])

    def test_return_faces_gives_the_input_triangulation(self):
        mesh = arrays(grid(4, 0.2))
        r = ddg_tools.make_developable(*mesh, subdivisions=0, iterations=5, remesh=False,
                                       return_faces=True)
        self.assertEqual(len(r.faces), 3 * 32)
        self.assertEqual(list(r.face_polygon), [p for p in range(16) for _ in range(2)])
        self.assertEqual(list(r.faces[:6]), [0, 1, 6, 0, 6, 5])

    def test_subdivision_branch(self):
        mesh = arrays(grid(4, 0.2))
        r = ddg_tools.make_developable(*mesh, subdivisions=1, iterations=20, remesh=False)
        self.assertEqual(len(r.faces), 3 * 4 * 32)
        self.assertEqual(len(r.face_polygon), 4 * 32)
        for k in range(0, len(r.face_polygon), 8):  # the 4 children of 2 triangles of a polygon
            self.assertEqual(len(set(r.face_polygon[k:k + 8])), 1)
        self.assertEqual(max(r.face_polygon), 15)
        nv = len(r.positions) // 3
        self.assertEqual(len(r.vertex_source), nv)
        self.assertEqual(list(r.vertex_source[:25]), list(range(25)))
        self.assertTrue(all(s == -1 for s in r.vertex_source[25:]))
        self.assertTrue(r.report["topology_changed"])
        self.assertEqual(len(r.report["rounds"]), 2)
        self.assertIsNone(r.report["nonplanarity_after"])
        self.assertTrue(all(0 <= v < nv for v in r.faces))

    def test_coarse_energy_sets_the_energy_of_round_zero(self):
        mesh = arrays(grid(4, 0.2))
        kw = dict(subdivisions=2, iterations=12, remesh=False)
        names = lambda r: [rd["energy"] for rd in r.report["rounds"]]
        default = ddg_tools.make_developable(*mesh, **kw)
        self.assertEqual(names(default), ["combinatorial", "covariance", "covariance"])
        same = ddg_tools.make_developable(*mesh, coarse_energy=None, **kw)
        self.assertEqual(names(same), ["covariance"] * 3)
        explicit = ddg_tools.make_developable(*mesh, coarse_energy="covariance", **kw)
        self.assertEqual(list(explicit.positions), list(same.positions))
        self.assertNotEqual(list(default.positions), list(same.positions))
        other = ddg_tools.make_developable(*mesh, energy="intrinsic", coarse_energy="covariance", **kw)
        self.assertEqual(names(other), ["covariance", "intrinsic", "intrinsic"])
        # Without a subdivision the one round uses `energy`.
        one = ddg_tools.make_developable(*mesh, subdivisions=0, iterations=5, remesh=False)
        self.assertEqual(names(one), ["covariance"])
        with self.assertRaises(ddg_tools.InputError) as c:
            ddg_tools.make_developable(*mesh, coarse_energy="bogus", **kw)
        self.assertIn("coarse_energy", str(c.exception))

    def test_fallback_operations_and_counts_in_the_report(self):
        p, cv, st, tr = grid(8, 0.3)
        p = [x * (1.0 if k % 3 else 6.0) for k, x in enumerate(p)]
        mesh = arrays((p, cv, st, tr))
        kw = dict(subdivisions=1, iterations=30, small_angle=0.5, min_face_angle=0.0)
        for fb in (True, False):
            r = ddg_tools.make_developable(*mesh, fallback_operations=fb, **kw)
            for rd in r.report["rounds"]:
                self.assertLessEqual(rd["fallback_flips"], rd["flips"])
                self.assertLessEqual(rd["fallback_collapses"], rd["collapses"])
                if not fb:
                    self.assertEqual((rd["fallback_flips"], rd["fallback_collapses"]), (0, 0))

    def test_face_polygon_is_exact_over_two_rounds(self):
        mesh = arrays(grid(4, 0.2))
        r = ddg_tools.make_developable(*mesh, subdivisions=2, iterations=10, remesh=False)
        self.assertEqual(len(r.face_polygon), 16 * 32)
        self.assertEqual(list(r.face_polygon), [p for p in range(16) for _ in range(32)])
        rep = r.report
        self.assertEqual(len(rep["rounds"]), 3)
        self.assertIsNone(rep["stopped"])
        self.assertEqual(rep["bad_faces"], 0)
        self.assertIsInstance(rep["inadmissible_trials"], int)
        self.assertGreater(rep["min_angle"], 0.1)
        for rd in rep["rounds"]:
            self.assertIn("inadmissible_trials", rd)
            self.assertIn("remesh_rejected", rd)
            self.assertNotIn("bad_triangles_left", rd)

    def test_face_polygon_survives_remeshing_exactly(self):
        # An aggressive small angle forces flips and collapses; every output triangle keeps the
        # polygon of its parent, which the library tracks (so the set of polygons is intact and
        # polygons far from the repairs are untouched).
        p, cv, st, tr = grid(8, 0.3)
        p = [x * (1.0 if k % 3 else 6.0) for k, x in enumerate(p)]
        r = ddg_tools.make_developable(*arrays((p, cv, st, tr)), subdivisions=1, iterations=30,
                                       small_angle=0.5, min_face_angle=0.0)
        self.assertEqual(len(r.face_polygon), len(r.faces) // 3)
        self.assertTrue(all(0 <= q < 64 for q in r.face_polygon))
        self.assertTrue(r.report["topology_changed"])
        flips = sum(rd["flips"] + rd["collapses"] for rd in r.report["rounds"])
        self.assertGreater(flips, 0)

    def test_min_face_angle(self):
        p, cv, st, tr = grid(2, 0.0)
        p = list(p)
        p[12], p[13] = 0.75, 1e-9  # the centre vertex next to the edge (0.5,0)-(1,0)
        mesh = arrays((p, cv, st, tr))
        with self.assertRaises(ddg_tools.MeshError) as c:
            ddg_tools.make_developable(*mesh, remesh=False)
        self.assertIn("nearly degenerate", str(c.exception))
        self.assertIn("min_face_angle", str(c.exception))
        r = ddg_tools.make_developable(*mesh, remesh=False, min_face_angle=0.0, subdivisions=0,
                                       iterations=2)
        self.assertEqual(len(r.positions), 27)
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.make_developable(*arrays(grid(3, 0.1)), small_angle=0.01, min_face_angle=0.02)

    def test_remeshing_keeps_polygons_valid(self):
        # A strongly sheared grid has small angles to repair.
        p, cv, st, tr = grid(8, 0.3)
        p = [x * (1.0 if k % 3 else 6.0) for k, x in enumerate(p)]
        r = ddg_tools.make_developable(*arrays((p, cv, st, tr)), subdivisions=0, iterations=40,
                                       small_angle=0.2)
        if r.faces is not None:
            self.assertEqual(len(r.face_polygon), len(r.faces) // 3)
            self.assertTrue(all(0 <= q < 64 for q in r.face_polygon))
            self.assertEqual(len(r.vertex_source), len(r.positions) // 3)

    def test_pinned_vertices_do_not_move(self):
        mesh = arrays(grid(6, 0.2))
        r = ddg_tools.make_developable(*mesh, subdivisions=0, iterations=30, remesh=False,
                                       pinned_vertices=array.array('i', [24, 25, 0]))
        for v in (24, 25, 0):
            for k in range(3):
                self.assertAlmostEqual(r.positions[3 * v + k], mesh[0][3 * v + k], places=6)

    def test_methods_and_line_searches(self):
        mesh = arrays(grid(6, 0.2))
        for kw in ({"method": "gradient_descent"}, {"method": "fixed_step", "step": 1e-3},
                   {"method": "lbfgs", "memory": 3, "line_search": "backtracking"},
                   {"energy": "combinatorial"}, {"energy": "intrinsic"}):
            r = ddg_tools.make_developable(*mesh, subdivisions=0, iterations=20, remesh=False, **kw)
            self.assertTrue(all(math.isfinite(x) for x in r.positions), kw)
            self.assertLessEqual(r.report["energy_final"], r.report["energy_initial"] * 1.0001, kw)

    def test_first_trial(self):
        # D37: "grow" starts each line search at 2x the previous step; a different path than "unit".
        mesh = arrays(grid(6, 0.2))
        kw = dict(subdivisions=0, iterations=30, remesh=False, energy="covariance_max")
        unit = ddg_tools.make_developable(*mesh, **kw)
        self.assertEqual(list(ddg_tools.make_developable(*mesh, first_trial="unit", **kw).positions),
                         list(unit.positions))
        grow = ddg_tools.make_developable(*mesh, first_trial="grow", **kw)
        self.assertTrue(all(math.isfinite(x) for x in grow.positions))
        self.assertLess(grow.report["energy_final"], grow.report["energy_initial"])
        self.assertNotEqual(list(grow.positions), list(unit.positions))
        capped = ddg_tools.make_developable(*mesh, first_trial="grow", first_trial_growth=1.0,
                                            first_trial_max=0.5, **kw)
        self.assertLess(capped.report["energy_final"], capped.report["energy_initial"])

    def test_deterministic_and_independent_of_threads(self):
        mesh = arrays(grid(8, 0.3))
        a = ddg_tools.make_developable(*mesh, subdivisions=1, iterations=25, threads=1)
        b = ddg_tools.make_developable(*mesh, subdivisions=1, iterations=25, threads=4)
        c = ddg_tools.make_developable(*mesh, subdivisions=1, iterations=25)
        for r in (b, c):
            self.assertEqual(a.positions, r.positions)
            self.assertEqual(a.faces, r.faces)
            self.assertEqual(a.face_polygon, r.face_polygon)
            self.assertEqual(a.vertex_source, r.vertex_source)

    def test_progress_and_cancellation(self):
        mesh = arrays(grid(8, 0.3))
        calls = []

        def stop(round_, iteration, energy, gradient_norm):
            calls.append((round_, iteration))
            self.assertTrue(math.isfinite(energy) and gradient_norm >= 0)
            return iteration < 4

        r = ddg_tools.make_developable(*mesh, subdivisions=1, iterations=100, progress=stop,
                                       progress_interval=0)
        self.assertTrue(r.report["cancelled"])
        self.assertEqual(len(r.report["rounds"]), 1)
        self.assertEqual(calls[-1], (0, 4))
        # None counts as "continue".
        n = []
        r = ddg_tools.make_developable(*mesh, subdivisions=0, iterations=8, remesh=False,
                                       progress=lambda *a: n.append(a), progress_interval=0)
        self.assertFalse(r.report["cancelled"])
        self.assertGreaterEqual(len(n), 3)

        def boom(*args):
            raise ValueError("boom")

        with self.assertRaises(ValueError):
            ddg_tools.make_developable(*mesh, subdivisions=0, iterations=8, progress=boom,
                                       progress_interval=0)

    def test_errors(self):
        mesh = arrays(grid(3, 0.1))
        for kw, exc in (({"energy": "x"}, ddg_tools.InputError), ({"method": "x"}, ddg_tools.InputError),
                        ({"line_search": "x"}, ddg_tools.InputError),
                        ({"first_trial": "x"}, ddg_tools.InputError),
                        ({"first_trial": "grow", "first_trial_growth": 0.5}, ddg_tools.InputError),
                        ({"first_trial": "grow", "first_trial_max": 0.0}, ddg_tools.InputError),
                        ({"iterations": 0}, ddg_tools.InputError),
                        ({"subdivisions": -1}, ddg_tools.InputError),
                        ({"max_displacement": -1.0}, ddg_tools.InputError),
                        ({"small_angle": 2.0}, ddg_tools.InputError),
                        ({"pinned_vertices": array.array('i', [99])}, ddg_tools.InputError)):
            with self.assertRaises(exc, msg=str(kw)):
                ddg_tools.make_developable(*mesh, **kw)
        with self.assertRaises(ddg_tools.MeshError):
            ddg_tools.make_developable(*arrays(three_fins()))


class SeamsAndUnwrapTest(unittest.TestCase):
    def test_cube_seams_are_polygon_edges_and_feed_the_unwrap(self):
        mesh = arrays(cube())
        s = ddg_tools.developable_seams(*mesh)
        pairs = list(zip(s.seam_edges[0::2], s.seam_edges[1::2]))
        self.assertGreaterEqual(len(pairs), 7)
        for a, b in pairs:
            self.assertEqual(bin(a ^ b).count("1"), 1, (a, b))  # cube edges, never diagonals
        self.assertEqual(len(s.terminals), 8)
        self.assertEqual(len(s.vertex_energy), 8)
        self.assertEqual(len(s.vertex_defect), 8)
        self.assertEqual(len(s.star_embeddedness), 8)
        self.assertTrue(all(abs(d - math.pi / 2) < 1e-9 for d in s.vertex_defect))
        u = ddg_tools.unwrap(*mesh, seam_edges=s.seam_edges)
        self.assertEqual(u.report["n_flipped_uv_triangles"], 0)

    def test_defect_rule_can_be_switched_off(self):
        self.assertAlmostEqual(ddg_tools.DEFAULT_DEFECT_TOLERANCE, 0.1)
        mesh = arrays(cube())
        a = ddg_tools.developable_seams(*mesh, defect_tolerance=None, epsilon=1e9)
        b = ddg_tools.developable_seams(*mesh, epsilon=1e9)
        self.assertEqual(len(a.terminals) , 2)  # the library adds two on a closed genus-0 piece
        self.assertEqual(len(b.terminals), 8)

    def test_user_seams_extra_terminals_and_errors(self):
        mesh = arrays(cube())
        s = ddg_tools.developable_seams(*mesh, seam_edges=array.array('i', [0, 1]),
                                        defect_tolerance=None, epsilon=1e9)
        self.assertIn((0, 1), list(zip(s.seam_edges[0::2], s.seam_edges[1::2])))
        s = ddg_tools.developable_seams(*mesh, extra_terminals=array.array('i', [0, 7]),
                                        defect_tolerance=None, epsilon=1e9)
        self.assertEqual(list(s.terminals), [0, 7])
        with self.assertRaises(ddg_tools.InputError):  # a diagonal of the z = 0 face
            ddg_tools.developable_seams(*mesh, seam_edges=array.array('i', [0, 3]))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.developable_seams(*mesh, extra_terminals=array.array('i', [8]))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.developable_seams(*mesh, extension="x")
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.developable_seams(*mesh, epsilon=-1.0)
        with self.assertRaises(ddg_tools.MeshError):
            ddg_tools.developable_seams(*arrays(three_fins()))

    def test_unwrap_of_a_cylinder_is_isometric_and_packed(self):
        mesh = arrays(cylinder())
        r = ddg_tools.developable_unwrap(*mesh)
        rep = r.report
        self.assertEqual(len(r.corner_uvs), 2 * 12 * 4 * 4)
        self.assertEqual(rep["flipped_triangles"], 0)
        self.assertTrue(rep["is_injective"], rep)
        self.assertEqual(rep["n_charts"], 1)
        self.assertLess(rep["mean_abs_log_scale"], 1e-6)
        self.assertEqual(len(rep["charts"]), 1)
        self.assertEqual(list(r.polygon_chart), [0] * 48)
        self.assertEqual(len(r.corner_log_scale), 192)
        self.assertGreater(len(r.seam_edges), 0)
        self.assertTrue(all(-1e-6 <= x <= 1 + 1e-6 for x in r.corner_uvs))
        r2 = ddg_tools.developable_unwrap(*mesh)
        self.assertEqual(r.corner_uvs, r2.corner_uvs)
        self.assertEqual(r.seam_edges, r2.seam_edges)

    def test_unwrap_options(self):
        mesh = arrays(cylinder())
        r = ddg_tools.developable_unwrap(*mesh, extension="conjugate", margin=0.0,
                                         alignment_weight=0.5, energy="combinatorial",
                                         overlap_samples_per_triangle=0.0)
        self.assertEqual(r.report["energy"], "combinatorial")
        self.assertEqual(r.report["pack_margin"], 0.0)
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.developable_unwrap(*mesh, margin=0.7)
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.developable_unwrap(*mesh, length_weight=0.0, alignment_weight=0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
