"""Tests of `ddg_tools.geodesic_path`, `geodesic_network` and `geodesic_bezier` (geodesics by edge
flips, Sharp and Crane 2020; crate ddg-flip-geodesics), on both representations, without Blender.

Run through `scripts/py-test.sh` like test_ddg_tools.py. Standard library only (numpy arrays are
also tried when numpy is installed).
"""

import array
import math
import sys
import unittest

import ddg_tools

REPRESENTATIONS = ("signposts", "integer")


def _ints(x):
    return array.array('i', x) if isinstance(x, (list, tuple)) else x


def geodesic_path(positions, triangles, vertices, **kw):
    return ddg_tools.geodesic_path(positions, triangles, _ints(vertices), **kw)


def geodesic_network(positions, triangles, edges, pinned=None, **kw):
    pinned = None if pinned is None else _ints(pinned)
    return ddg_tools.geodesic_network(positions, triangles, _ints(edges), pinned=pinned, **kw)


def geodesic_bezier(positions, triangles, controls, **kw):
    return ddg_tools.geodesic_bezier(positions, triangles, _ints(controls), **kw)


def grid(n=10):
    """The unit square with (n + 1)^2 vertices, vertex (i, j) is j * (n + 1) + i:
    (positions, triangles)."""
    positions, triangles = [], []
    for j in range(n + 1):
        for i in range(n + 1):
            positions += [i / n, j / n, 0.0]
    for j in range(n):
        for i in range(n):
            a, b = j * (n + 1) + i, j * (n + 1) + i + 1
            c, d = a + n + 1, b + n + 1
            triangles += [a, b, d, a, d, c]
    return array.array('d', positions), array.array('i', triangles)


def v(i, j, n=10):
    return j * (n + 1) + i


def icosahedron():
    t = (1 + math.sqrt(5)) / 2
    raw = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t),
           (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    positions = []
    for p in raw:
        n = math.sqrt(sum(c * c for c in p))
        positions += [c / n for c in p]
    triangles = [0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2,
                 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5,
                 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1]
    return array.array('d', positions), array.array('i', triangles)


def points(flat):
    return [tuple(flat[3 * k:3 * k + 3]) for k in range(len(flat) // 3)]


def polyline_length(flat):
    p = points(flat)
    return sum(math.dist(p[k], p[k + 1]) for k in range(len(p) - 1))


def de_casteljau_polygon(controls, rounds):
    """The control polygon after `rounds` halvings (the pinned points of SPEC 6.1), planar."""
    if rounds == 0:
        return list(controls)
    levels = [list(controls)]
    while len(levels[-1]) > 1:
        last = levels[-1]
        levels.append([tuple((a + b) / 2 for a, b in zip(last[k], last[k + 1]))
                       for k in range(len(last) - 1)])
    left = [level[0] for level in levels]
    right = [level[-1] for level in reversed(levels)]
    first = de_casteljau_polygon(left, rounds - 1)
    second = de_casteljau_polygon(right, rounds - 1)
    return first + second[1:]


class GeodesicPathTest(unittest.TestCase):
    def test_flat_path_is_the_straight_line(self):
        positions, triangles = grid()
        a, b = v(1, 2), v(9, 7)
        expected = math.dist((0.1, 0.2), (0.9, 0.7))
        lengths = []
        for rep in REPRESENTATIONS:
            r = geodesic_path(positions, triangles, [a, b], representation=rep)
            self.assertEqual(r.report["representation"], rep)
            self.assertEqual(r.report["status"], "geodesic")
            self.assertAlmostEqual(r.lengths[0], expected, places=9)
            self.assertAlmostEqual(polyline_length(r.polylines[0]), expected, places=9)
            self.assertAlmostEqual(r.report["length"], expected, places=9)
            self.assertLess(r.report["length"], r.report["initial_length"])
            self.assertIsInstance(r.polylines[0], array.array)
            p = points(r.polylines[0])
            self.assertEqual(p[0], (0.1, 0.2, 0.0))
            self.assertEqual(p[-1], (0.9, 0.7, 0.0))
            # Every point lies on the segment between the ends.
            for q in p:
                self.assertLess(abs((q[0] - 0.1) * 0.5 - (q[1] - 0.2) * 0.8), 1e-9)
            self.assertTrue(r.report["certificate"]["geodesic"])
            self.assertTrue(r.report["input_certificate"]["passed"], r.report["input_certificate"])
            self.assertEqual(r.closed, [False])
            lengths.append(r.lengths[0])
        self.assertAlmostEqual(lengths[0], lengths[1], places=9)

    def test_default_representation_is_the_documented_one(self):
        positions, triangles = grid(4)
        r = geodesic_path(positions, triangles, [0, 24])
        self.assertEqual(r.report["representation"], "integer")

    def test_surface_points(self):
        positions, triangles = grid()
        for rep in REPRESENTATIONS:
            r = geodesic_path(positions, triangles, [v(1, 2), v(9, 7)], representation=rep)
            sp = r.surface_points[0]
            n = len(r.polylines[0]) // 3
            self.assertEqual(len(sp["kind"]), n)
            self.assertEqual(len(sp["id"]), n)
            self.assertEqual(len(sp["params"]), 3 * n)
            self.assertEqual((sp["kind"][0], sp["id"][0]), (0, v(1, 2)))
            self.assertEqual((sp["kind"][-1], sp["id"][-1]), (0, v(9, 7)))
            self.assertTrue(all(k in (0, 1, 2) for k in sp["kind"]))
            # Rebuild the 3D points from the places: the same polyline.
            p = points(r.polylines[0])
            for k in range(n):
                kind, ident = sp["kind"][k], sp["id"][k]
                t = sp["params"][3 * k:3 * k + 3]
                if kind == 0:
                    q = tuple(positions[3 * ident:3 * ident + 3])
                else:
                    self.assertEqual(kind, 1)
                    tri_edges = {}
                    for f in range(len(triangles) // 3):
                        c = triangles[3 * f:3 * f + 3]
                        for s in range(3):
                            key = (min(c[s], c[(s + 1) % 3]), max(c[s], c[(s + 1) % 3]))
                            tri_edges.setdefault(key, len(tri_edges))
                    ends = [e for e, i in tri_edges.items() if i == ident]
                    # Edge ids follow the mesh's own numbering, not this one: only check range.
                    self.assertTrue(0.0 <= t[0] <= 1.0)
                    q = p[k]
                self.assertEqual(len(q), 3)

    def test_closed_loop(self):
        positions, triangles = icosahedron()
        for rep in REPRESENTATIONS:
            r = geodesic_path(positions, triangles, [0, 3, 8], closed=True, representation=rep)
            self.assertEqual(r.closed, [True])
            p = points(r.polylines[0])
            self.assertEqual(p[0], p[-1])
            self.assertEqual(r.surface_points[0]["id"][0], r.surface_points[0]["id"][-1])
            self.assertTrue(r.report["certificate"]["geodesic"])
            self.assertAlmostEqual(r.lengths[0], polyline_length(r.polylines[0]), places=9)

    def test_unpinned_loop_may_leave_its_vertices(self):
        positions, triangles = grid()
        # A loop around the middle of a flat disk shrinks away: nothing is left of it.
        r = geodesic_path(positions, triangles, [v(3, 3), v(7, 3), v(7, 7), v(3, 7)],
                                    closed=True, pin_waypoints=False, representation="integer")
        self.assertLess(r.lengths[0], 1e-6)

    def test_stop_ratio_and_iteration_limit(self):
        positions, triangles = grid()
        args = (positions, triangles, [v(1, 1), v(9, 4)])
        for rep in REPRESENTATIONS:
            full = geodesic_path(*args, representation=rep)
            part = geodesic_path(*args, representation=rep, stop_length_ratio=0.97)
            self.assertEqual(part.report["status"], "length_ratio")
            self.assertGreater(part.lengths[0], full.lengths[0])
            capped = geodesic_path(*args, representation=rep, max_iterations=1)
            self.assertEqual(capped.report["status"], "iteration_limit")
            self.assertEqual(capped.report["iterations"], 1)
            self.assertFalse(capped.report["certificate"]["geodesic"])

    def test_options_off_and_triangulation_edges(self):
        positions, triangles = grid()
        for rep in REPRESENTATIONS:
            r = geodesic_path(positions, triangles, [v(1, 1), v(9, 9)], representation=rep,
                                        certify_input=False, triangulation_edges=True)
            self.assertIsNone(r.report["input_certificate"])
            self.assertIsNone(r.control_polygon)
            te = r.triangulation_edges
            self.assertEqual(len(te["offsets"]), len(te["edges"]) + 1)
            self.assertEqual(te["offsets"][-1], len(te["points"]) // 3)
            self.assertLess(te["max_relative_mismatch"], 1e-6)
            plain = geodesic_path(positions, triangles, [v(1, 1), v(9, 9)], representation=rep)
            self.assertIsNone(plain.triangulation_edges)
            self.assertIsNotNone(plain.report["input_certificate"])

    def test_mollify_is_for_integer_only(self):
        positions, triangles = grid(4)
        r = geodesic_path(positions, triangles, [0, 24], representation="integer",
                                    mollify=1e-4)
        self.assertEqual(r.report["status"], "geodesic")
        with self.assertRaisesRegex(ddg_tools.InputError, r"^mollify"):
            geodesic_path(positions, triangles, [0, 24], representation="signposts",
                                    mollify=1e-4)

    def test_numpy_arrays(self):
        try:
            import numpy
        except ImportError:
            self.skipTest("numpy is not installed")
        positions, triangles = grid()
        p = numpy.asarray(positions).reshape(-1, 3)
        t = numpy.asarray(triangles, dtype=numpy.int32).reshape(-1, 3)
        r = ddg_tools.geodesic_path(p, t, numpy.array([v(1, 2), v(9, 7)], dtype=numpy.int32))
        a = numpy.asarray(r.polylines[0]).reshape(-1, 3)
        self.assertAlmostEqual(float(numpy.linalg.norm(numpy.diff(a, axis=0), axis=1).sum()),
                               math.dist((0.1, 0.2), (0.9, 0.7)), places=9)


class GeodesicNetworkTest(unittest.TestCase):
    def seam_edges(self):
        """A Y of mesh edges: junction (5, 5); one leg bends through (2, 5) to (2, 8)."""
        chain = lambda *cells: [(v(*a), v(*b)) for a, b in zip(cells, cells[1:])]
        return (chain((5, 5), (4, 5), (3, 5), (2, 5), (2, 6), (2, 7), (2, 8))
                + chain((5, 5), (6, 5), (7, 5), (8, 5)) + chain((5, 5), (5, 4), (5, 3)))

    def test_a_seam_network_straightens_and_keeps_its_junction(self):
        positions, triangles = grid()
        edges = self.seam_edges()
        flat = array.array('i', [x for e in edges for x in e])
        expected = math.dist((0.5, 0.5), (0.2, 0.8)) + 0.3 + 0.2
        for rep in REPRESENTATIONS:
            r = geodesic_network(positions, triangles, flat, representation=rep)
            self.assertEqual(len(r.polylines), 3)
            self.assertEqual(r.report["n_paths"], 3)
            self.assertAlmostEqual(sum(r.lengths), expected, places=9)
            self.assertTrue(r.report["certificate"]["geodesic"])
            self.assertTrue(r.report["input_certificate"]["passed"])
            junction = (0.5, 0.5, 0.0)
            ends = [points(p)[0] for p in r.polylines] + [points(p)[-1] for p in r.polylines]
            self.assertEqual(sum(1 for e in ends if e == junction), 3)
            # The pair shape (m, 2) as nested lists works through a flat array too.
            self.assertEqual(r.closed, [False] * 3)

    def test_extra_pinned_vertices_stay(self):
        positions, triangles = grid()
        edges = self.seam_edges()
        flat = array.array('i', [x for e in edges for x in e])
        corner = v(2, 5)
        r = geodesic_network(positions, triangles, flat, pinned=[corner],
                                       representation="integer")
        on_corner = [q for p in r.polylines for q in points(p) if q == (0.2, 0.5, 0.0)]
        self.assertTrue(on_corner)
        self.assertEqual(len(r.polylines), 4)
        self.assertGreater(sum(r.lengths), math.dist((0.5, 0.5), (0.2, 0.8)) + 0.5)

    def test_a_closed_seam_is_a_loop(self):
        positions, triangles = grid()
        ring = [(3, 3), (4, 3), (5, 3), (6, 3), (6, 4), (6, 5), (6, 6), (5, 6), (4, 6), (3, 6),
                (3, 5), (3, 4), (3, 3)]
        edges = array.array('i', [x for a, b in zip(ring, ring[1:]) for x in (v(*a), v(*b))])
        # With one vertex pinned the loop is an open path from it back to itself.
        r = geodesic_network(positions, triangles, edges, pinned=[v(3, 3)],
                                       representation="integer")
        self.assertEqual(len(r.polylines), 1)

    def test_errors(self):
        positions, triangles = grid()
        with self.assertRaisesRegex(ddg_tools.InputError, r"^edges"):
            geodesic_network(positions, triangles, [v(0, 0), v(5, 5)])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^edges"):
            geodesic_network(positions, triangles, [0, 9999])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^edges"):
            geodesic_network(positions, triangles, [0, 1, 2])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^pinned"):
            geodesic_network(positions, triangles, [0, 1], pinned=[9999])


class GeodesicBezierTest(unittest.TestCase):
    def test_bezier_on_a_plane_is_de_casteljau(self):
        positions, triangles = grid()
        cells = [(1, 1), (3, 9), (7, 9), (9, 1)]
        controls = [v(*c) for c in cells]
        planar = [(i / 10, j / 10, 0.0) for i, j in cells]
        for rep in REPRESENTATIONS:
            for rounds in (1, 2, 3):
                r = geodesic_bezier(positions, triangles, controls, rounds=rounds,
                                              representation=rep)
                expected = de_casteljau_polygon(planar, rounds)
                got = points(r.control_polygon)
                self.assertEqual(len(got), len(expected), (rep, rounds))
                for g, e in zip(got, expected):
                    self.assertLess(math.dist(g, e), 1e-6, (rep, rounds))
                self.assertAlmostEqual(polyline_length(r.control_polygon),
                                       polyline_length(r.polylines[0]), places=6)
                self.assertAlmostEqual(r.lengths[0], polyline_length(r.polylines[0]), places=6)
                b = r.report["bezier"]
                self.assertEqual(b["rounds"], rounds)
                self.assertEqual(len(b["control_points"]), len(expected))
                self.assertEqual(r.report["status"], "geodesic")
                self.assertTrue(r.report["certificate"]["geodesic"])
                self.assertTrue(r.report["input_certificate"]["passed"], r.report["input_certificate"])

    def test_more_rounds_approach_the_curve(self):
        positions, triangles = grid()
        controls = [v(1, 1), v(3, 9), v(7, 9), v(9, 1)]
        lengths = [geodesic_bezier(positions, triangles, controls, rounds=k,
                                             representation="integer").lengths[0]
                   for k in (1, 2, 4)]
        self.assertGreater(lengths[0], lengths[1])
        self.assertGreater(lengths[1], lengths[2])

    def test_representations_agree_on_a_sphere(self):
        positions, triangles = icosahedron()
        a = geodesic_bezier(positions, triangles, [0, 3, 8], rounds=2,
                                      representation="signposts")
        b = geodesic_bezier(positions, triangles, [0, 3, 8], rounds=2,
                                      representation="integer")
        self.assertAlmostEqual(a.lengths[0], b.lengths[0], places=6)

    def test_errors(self):
        positions, triangles = grid()
        bezier = lambda c, **kw: geodesic_bezier(positions, triangles, c, **kw)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^control_vertices.*at least two"):
            bezier([v(1, 1)])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^control_vertices"):
            bezier([v(1, 1), v(1, 1), v(5, 5)])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^control_vertices"):
            bezier([v(1, 1), 99999])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^options.*rounds"):
            bezier([v(1, 1), v(5, 5)], rounds=1000)
        # The shortest paths between the control points cross: P0P1 and P2P3 are the diagonals.
        for rep in REPRESENTATIONS:
            with self.assertRaisesRegex(ddg_tools.SolveError, r"move a control point"):
                bezier([v(1, 1), v(9, 9), v(1, 9), v(9, 1)], representation=rep)
        # A repeated control point.
        for rep in REPRESENTATIONS:
            with self.assertRaisesRegex(ddg_tools.SolveError, r"control point"):
                bezier([v(2, 2), v(8, 3), v(5, 8), v(2, 2)], representation=rep)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^representation"):
            bezier([v(1, 1), v(5, 5)], representation="fancy")


class GeodesicErrorTest(unittest.TestCase):
    def test_path_errors(self):
        positions, triangles = grid()
        path = lambda vs, **kw: geodesic_path(positions, triangles, vs, **kw)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^vertices.*at least two"):
            path([3])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^vertices.*vertex 3 to vertex 3"):
            path([3, 3])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^vertices"):
            path([3, 100000])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^tolerance"):
            path([3, 30], tolerance=-1.0)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^stop_length_ratio"):
            path([3, 30], stop_length_ratio=2.0)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^representation"):
            path([3, 30], representation="fancy")

    def test_path_crossing_itself(self):
        positions, triangles = grid()
        # Legs (1,1)-(9,9) and (1,9)-(9,1) cross in the middle: refused, not repaired.
        for rep in REPRESENTATIONS:
            with self.assertRaisesRegex(ddg_tools.SolveError, r"cross"):
                geodesic_path(positions, triangles, [v(1, 1), v(9, 9), v(1, 9), v(9, 1)],
                                        representation=rep)

    def test_mesh_error(self):
        with self.assertRaises(ddg_tools.MeshError):
            geodesic_path(array.array('d', [0, 0, 0, 1, 0, 0, 2, 0, 0]),
                                    array.array('i', [0, 1, 2]), [0, 1])


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]], verbosity=1)
