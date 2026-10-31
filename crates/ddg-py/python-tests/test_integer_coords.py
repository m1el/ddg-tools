"""Tests of `ddg_tools.IntrinsicTriangulation` (integer coordinates, Gillespie, Sharp, Crane
2021; crate ddg-integer-coords) and of the `representation=` option of the intrinsic bindings,
without Blender.

Run through `scripts/py-test.sh` like test_ddg_tools.py. Standard library only.
"""

import array
import math
import sys
import unittest

import ddg_tools


def icosahedron():
    """The icosahedron on the unit sphere: (positions, triangles)."""
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


def sheared_grid(n=8, shear=2.0):
    """A flat parallelogram of obtuse triangles (with boundary): (positions, triangles)."""
    positions, triangles = [], []
    for j in range(n + 1):
        for i in range(n + 1):
            positions += [(i + shear * j) / n, j / n, 0.0]
    for j in range(n):
        for i in range(n):
            a, b = j * (n + 1) + i, j * (n + 1) + i + 1
            c, d = a + n + 1, b + n + 1
            triangles += [a, b, d, a, d, c]
    return array.array('d', positions), array.array('i', triangles)


def prepared(mesh, representation, rules):
    t = ddg_tools.IntrinsicTriangulation(*mesh, representation=representation)
    t.flip_to_delaunay()
    t.refine(25.0, rules)
    return t


class IntrinsicTriangulationTest(unittest.TestCase):
    def test_a_fresh_triangulation_is_the_input(self):
        positions, triangles = icosahedron()
        t = ddg_tools.IntrinsicTriangulation(positions, triangles)
        self.assertEqual(t.representation, "integer")
        c = t.counts()
        self.assertEqual((c["n_vertices"], c["n_faces"], c["n_edges"]), (12, 20, 30))
        self.assertEqual(c["n_crossings"], 0)
        self.assertAlmostEqual(c["min_angle"], math.pi / 3, places=12)
        tris, lengths = t.triangles()
        self.assertEqual(list(tris), list(triangles))
        self.assertEqual(len(lengths), 30)
        self.assertEqual(t.certify()["n_faces"], 20)

    def test_flip_and_refine(self):
        mesh = sheared_grid()
        for rep in ("integer", "signposts"):
            for rules in ("chew", "ruppert"):
                t = ddg_tools.IntrinsicTriangulation(*mesh, representation=rep)
                self.assertEqual(t.flip_to_delaunay()["status"], "delaunay")
                report = t.refine(25.0, rules)
                self.assertEqual(report["status"], "done")
                self.assertGreater(report["n_inserted"], 0)
                # The parallelogram has narrow corners that refinement leaves; elsewhere 25 degrees.
                self.assertGreater(math.degrees(report["min_angle_valid"]), 24.9, (rep, rules))
                self.assertEqual(report["certified"], rep == "integer")
                if rep == "integer":
                    t.certify()

    def test_representations_agree(self):
        for mesh in (icosahedron(), sheared_grid()):
            for rules in ("chew", "ruppert"):
                a, b = prepared(mesh, "signposts", rules), prepared(mesh, "integer", rules)
                ca, cb = a.counts(), b.counts()
                for key in ("n_vertices", "n_edges", "n_faces"):
                    self.assertEqual(ca[key], cb[key])
                self.assertAlmostEqual(ca["min_angle"], cb["min_angle"], places=9)
                sa, sb = a.common_subdivision(), b.common_subdivision()
                self.assertEqual(list(sa["polygon_starts"]), list(sb["polygon_starts"]))
                self.assertEqual(list(sa["input_face"]), list(sb["input_face"]))
                self.assertEqual(list(sa["intrinsic_face"]), list(sb["intrinsic_face"]))
                worst = max(abs(x - y) for x, y in zip(sa["positions"], sb["positions"]))
                self.assertLess(worst, 1e-8)

    def test_common_subdivision(self):
        mesh = icosahedron()
        t = prepared(mesh, "integer", "chew")
        c, s = t.counts(), t.common_subdivision()
        nv = len(s["positions"]) // 3
        self.assertEqual(len(s["on_input"]["kind"]), nv)
        self.assertEqual(len(s["on_triangulation"]["id"]), nv)
        self.assertEqual(s["n_intrinsic_vertices"], c["n_vertices"])
        self.assertEqual(nv - c["n_vertices"], c["n_crossings"])
        n_poly = len(s["input_face"])
        self.assertEqual(len(s["polygon_starts"]), n_poly + 1)
        self.assertEqual(s["polygon_starts"][-1], len(s["corners"]))
        self.assertTrue(all(0 <= v < nv for v in s["corners"]))
        self.assertEqual(set(s["intrinsic_face"]), set(range(c["n_faces"])))
        self.assertEqual(set(s["input_face"]), set(range(c["n_input_faces"])))
        # Every polygon has at least three corners; the input vertices are first and sit at
        # their own positions.
        starts = s["polygon_starts"]
        self.assertTrue(all(starts[i + 1] - starts[i] >= 3 for i in range(n_poly)))
        for v in range(12):
            self.assertEqual(s["on_input"]["kind"][v], 0)
            self.assertEqual(s["on_input"]["id"][v], v)
        # Total area of the polygons: the area of the input (planar polygons on flat faces).
        pos = s["positions"]
        area = 0.0
        for i in range(n_poly):
            ids = s["corners"][starts[i]:starts[i + 1]]
            p = [pos[3 * k:3 * k + 3] for k in ids]
            for k in range(1, len(p) - 1):
                u = [p[k][a] - p[0][a] for a in range(3)]
                w = [p[k + 1][a] - p[0][a] for a in range(3)]
                cross = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                         u[0] * w[1] - u[1] * w[0]]
                area += 0.5 * math.sqrt(sum(x * x for x in cross))
        positions, triangles = mesh
        total = 0.0
        for f in range(20):
            a, b, d = (positions[3 * triangles[3 * f + k]:3 * triangles[3 * f + k] + 3]
                       for k in range(3))
            u = [b[i] - a[i] for i in range(3)]
            w = [d[i] - a[i] for i in range(3)]
            cross = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
            total += 0.5 * math.sqrt(sum(x * x for x in cross))
        self.assertAlmostEqual(area, total, places=9)

    def test_edge_polylines(self):
        t = prepared(icosahedron(), "integer", "chew")
        c = t.counts()
        p = t.edge_polylines()
        self.assertEqual(len(p["edges"]), c["n_edges"])
        self.assertEqual(len(p["offsets"]), c["n_edges"] + 1)
        self.assertEqual(p["offsets"][-1] * 3, len(p["points"]))
        self.assertEqual(sum(p["n_crossings"]), c["n_crossings"])
        self.assertLess(p["max_relative_mismatch"], 1e-6)
        places = t.vertex_places()
        self.assertEqual(len(places["positions"]), 3 * c["n_vertices"])
        # Each polyline runs from the place of one end to the place of the other.
        tris, _ = t.triangles()
        ends = {}
        for f in range(c["n_faces"]):
            for k in range(3):
                ends.setdefault(frozenset((tris[3 * f + k], tris[3 * f + (k + 1) % 3])), None)
        self.assertEqual(len(ends), c["n_edges"])
        some = t.edge_polylines(array.array('i', [4, 2]))
        self.assertEqual(list(some["edges"]), [4, 2])
        with self.assertRaisesRegex(ddg_tools.InputError, r"^edges"):
            t.edge_polylines(array.array('i', [10_000]))
        # Signposts give the same polylines.
        s = prepared(icosahedron(), "signposts", "chew").edge_polylines()
        self.assertEqual(list(s["offsets"]), list(p["offsets"]))
        self.assertLess(max(abs(x - y) for x, y in zip(s["points"], p["points"])), 1e-8)

    def test_transfer(self):
        t = prepared(icosahedron(), "integer", "chew")
        c = t.counts()
        down = t.to_input(array.array('d', [1.0] * c["n_vertices"]))
        self.assertEqual(len(down), 12)
        self.assertLess(max(abs(x - 1.0) for x in down), 1e-10)
        up = t.to_intrinsic(down)
        self.assertEqual(len(up), c["n_vertices"])
        self.assertLess(max(abs(x - 1.0) for x in up), 1e-10)
        with self.assertRaises(ddg_tools.SolveError):
            t.to_input(array.array('d', [1.0]))

    def test_mollify_and_errors(self):
        positions, triangles = sheared_grid(4)
        t = ddg_tools.IntrinsicTriangulation(positions, triangles, mollify=1e-4)
        t.certify()
        with self.assertRaisesRegex(ddg_tools.InputError, r"^mollify"):
            ddg_tools.IntrinsicTriangulation(positions, triangles, "signposts", 1e-4)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^representation"):
            ddg_tools.IntrinsicTriangulation(positions, triangles, "nope")
        with self.assertRaisesRegex(ddg_tools.InputError, r"^rules"):
            t.refine(25.0, "nope")
        with self.assertRaisesRegex(ddg_tools.InputError, r"^min_angle_degrees"):
            t.refine(45.0)
        s = ddg_tools.IntrinsicTriangulation(positions, triangles, "signposts")
        with self.assertRaisesRegex(ddg_tools.InputError, r"^representation"):
            s.certify()
        with self.assertRaises(ddg_tools.MeshError):
            ddg_tools.IntrinsicTriangulation(array.array('d', [0, 0, 0, 1, 0, 0, 2, 0, 0]),
                                             array.array('i', [0, 1, 2]))
        with self.assertRaises(TypeError):
            ddg_tools.IntrinsicTriangulation(positions, 3)


class RepresentationOptionTest(unittest.TestCase):
    def test_default_is_signposts_and_both_agree(self):
        positions, triangles = icosahedron()
        sources = array.array('i', [0])
        default = ddg_tools.geodesic_distance(positions, triangles, sources)
        a = ddg_tools.geodesic_distance(positions, triangles, sources, representation="signposts")
        b = ddg_tools.geodesic_distance(positions, triangles, sources, representation="integer")
        self.assertEqual(list(default.distances), list(a.distances))
        self.assertLess(max(abs(x - y) for x, y in zip(a.distances, b.distances)), 1e-10)
        self.assertEqual(a.report["n_triangulation_vertices"], b.report["n_triangulation_vertices"])

        ta = ddg_tools.steiner_tree(positions, triangles, representation="signposts")
        tb = ddg_tools.steiner_tree(positions, triangles, representation="integer")
        self.assertAlmostEqual(ta.report["length"], tb.report["length"], places=10)
        self.assertEqual(ta.report["n_steiner_points"], tb.report["n_steiner_points"])

        sa = ddg_tools.adaptive_solve(positions, triangles, "heat", source=3, time=0.05,
                                      max_vertices=400)
        sb = ddg_tools.adaptive_solve(positions, triangles, "heat", source=3, time=0.05,
                                      max_vertices=400, representation="integer")
        self.assertLess(max(abs(x - y) for x, y in zip(sa.values, sb.values)), 1e-8)

    def test_unknown_representation(self):
        positions, triangles = icosahedron()
        for call in (
            lambda: ddg_tools.geodesic_distance(positions, triangles, array.array('i', [0]),
                                                representation="mixed"),
            lambda: ddg_tools.steiner_tree(positions, triangles, representation="mixed"),
            lambda: ddg_tools.adaptive_solve(positions, triangles, "heat", representation="x"),
        ):
            with self.assertRaisesRegex(ddg_tools.InputError, r"^representation"):
                call()


if __name__ == "__main__":
    print("python", sys.version.split()[0], "| ddg_tools", ddg_tools.__version__,
          "from", ddg_tools.__file__)
    unittest.main(verbosity=2)
