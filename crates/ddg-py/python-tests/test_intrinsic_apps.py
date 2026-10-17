"""Tests of `ddg_tools.geodesic_distance` and `ddg_tools.steiner_tree` (applications of
intrinsic triangulations, crate ddg-intrinsic-apps), without Blender.

Run through `scripts/py-test.sh` like test_ddg_tools.py. Standard library only.
"""

import array
import math
import sys
import unittest

import ddg_tools


def octahedron():
    """The regular octahedron with vertices on the unit axes: (positions, triangles)."""
    positions = [1.0, 0, 0, -1.0, 0, 0, 0, 1.0, 0, 0, -1.0, 0, 0, 0, 1.0, 0, 0, -1.0]
    triangles = [0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5]
    return array.array('d', positions), array.array('i', triangles)


def sheared_grid(n=12, shear=3.0):
    """A flat parallelogram of obtuse triangles: (positions, triangles)."""
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


class GeodesicDistanceTest(unittest.TestCase):
    def test_octahedron(self):
        positions, triangles = octahedron()
        for triangulation in ("refined", "delaunay", "input"):
            result = ddg_tools.geodesic_distance(positions, triangles, array.array('i', [0]),
                                                 triangulation=triangulation)
            d = result.distances
            self.assertIsInstance(d, array.array)
            self.assertEqual(d.typecode, 'd')
            self.assertEqual(len(d), 6)
            self.assertAlmostEqual(d[0], 0.0, places=12)
            # Neighbours along an edge (length sqrt 2); the opposite vertex is the farthest.
            self.assertTrue(all(d[v] > 0.5 for v in range(1, 6)))
            self.assertEqual(max(range(6), key=lambda v: d[v]), 1)
            self.assertIn("seconds_factor", result.report)
            self.assertGreaterEqual(result.report["n_triangulation_vertices"], 6)

    def test_flat_grid_is_close_to_euclidean(self):
        positions, triangles = sheared_grid()
        source = 6 * 13 + 6
        result = ddg_tools.geodesic_distance(positions, triangles, array.array('i', [source]))
        p = positions
        errors = []
        for v in range(len(p) // 3):
            e = math.dist(p[3 * v:3 * v + 3], p[3 * source:3 * source + 3])
            if e > 0:
                errors.append(abs(result.distances[v] - e) / e)
        self.assertLess(sum(errors) / len(errors), 0.05)

    def test_errors(self):
        positions, triangles = octahedron()
        with self.assertRaisesRegex(ddg_tools.InputError, r"^sources"):
            ddg_tools.geodesic_distance(positions, triangles, array.array('i', [17]))
        with self.assertRaisesRegex(ddg_tools.InputError, r"^triangulation"):
            ddg_tools.geodesic_distance(positions, triangles, array.array('i', [0]),
                                        triangulation="fancy")
        with self.assertRaisesRegex(ddg_tools.InputError, r"^time_factor"):
            ddg_tools.geodesic_distance(positions, triangles, array.array('i', [0]),
                                        time_factor=0.0)
        with self.assertRaises(ddg_tools.MeshError):
            ddg_tools.geodesic_distance(array.array('d', [0, 0, 0, 1, 0, 0, 2, 0, 0]),
                                        array.array('i', [0, 1, 2]), array.array('i', [0]))
        flat, flat_triangles = sheared_grid(4)
        with self.assertRaisesRegex(ddg_tools.InputError, r"^sources"):
            ddg_tools.geodesic_distance(flat, flat_triangles, array.array('i', [0]),
                                        boundary="dirichlet")


class SteinerTreeTest(unittest.TestCase):
    def test_octahedron(self):
        positions, triangles = octahedron()
        mst = ddg_tools.steiner_tree(positions, triangles, steiner_points=False)
        steiner = ddg_tools.steiner_tree(positions, triangles)
        # The MST of the octahedron's edges: 5 edges of length sqrt 2.
        self.assertAlmostEqual(mst.report["length"], 5 * math.sqrt(2), places=12)
        self.assertLess(steiner.report["length"], mst.report["length"])
        self.assertGreater(steiner.report["n_steiner_points"], 0)
        for result in (mst, steiner):
            self.assertTrue(result.report["is_disk"])
            self.assertLess(result.report["max_interior_defect"], 1e-9)
            points, offsets = result.points, result.offsets
            self.assertEqual(offsets[-1] * 3, len(points))
            length = 0.0
            for k in range(len(offsets) - 1):
                for i in range(offsets[k], offsets[k + 1] - 1):
                    length += math.dist(points[3 * i:3 * i + 3], points[3 * i + 3:3 * i + 6])
            self.assertAlmostEqual(length, result.report["length"], places=9)

    def test_errors(self):
        with self.assertRaisesRegex(ddg_tools.InputError, r"^positions"):
            ddg_tools.steiner_tree(array.array('d', [0.0] * 4), array.array('i', [0, 1, 2]))
        with self.assertRaises(TypeError):
            ddg_tools.steiner_tree(octahedron()[0], 3)


if __name__ == "__main__":
    print("python", sys.version.split()[0], "| ddg_tools", ddg_tools.__version__,
          "from", ddg_tools.__file__)
    unittest.main(verbosity=2)
