"""Tests of `ddg_tools.map_to_sphere` on the built wheel, without Blender.

Run through `scripts/py-test.sh` (see test_ddg_tools.py). Standard library only.
"""

import array
import math
import unittest

import ddg_tools


def octahedron(level):
    """A regular octahedron with each face split into 4**level flat triangles: (positions,
    triangles) as flat lists. Closed, genus 0, outward oriented."""
    positions = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0],
                 [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]
    faces = [(0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4), (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5)]
    for _ in range(level):
        midpoint = {}

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in midpoint:
                midpoint[key] = len(positions)
                positions.append([(positions[a][i] + positions[b][i]) / 2 for i in range(3)])
            return midpoint[key]

        new = []
        for a, b, c in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
        faces = new
    return [x for p in positions for x in p], [v for f in faces for v in f]


def torus(n_major=8, n_minor=6, big=2.0, small=0.5):
    positions, triangles = [], []
    for i in range(n_major):
        for j in range(n_minor):
            u, v = 2 * math.pi * i / n_major, 2 * math.pi * j / n_minor
            r = big + small * math.cos(v)
            positions += [r * math.cos(u), r * math.sin(u), small * math.sin(v)]
    for i in range(n_major):
        for j in range(n_minor):
            a = i * n_minor + j
            b = ((i + 1) % n_major) * n_minor + j
            c = ((i + 1) % n_major) * n_minor + (j + 1) % n_minor
            d = i * n_minor + (j + 1) % n_minor
            triangles += [a, b, c, a, c, d]
    return positions, triangles


class MapToSphereTest(unittest.TestCase):
    def test_octahedron_maps_injectively_to_the_unit_sphere(self):
        positions, triangles = octahedron(3)
        result = ddg_tools.map_to_sphere(array.array('d', positions), array.array('i', triangles))
        self.assertIsInstance(result.positions, array.array)
        self.assertEqual(result.positions.typecode, 'd')
        self.assertEqual(len(result.positions), len(positions))
        for k in range(0, len(result.positions), 3):
            x, y, z = result.positions[k:k + 3]
            self.assertAlmostEqual(math.sqrt(x * x + y * y + z * z), 1.0, places=12)
        report = result.report
        self.assertTrue(report["is_injective"])
        self.assertEqual(report["flipped_triangles"], 0)
        self.assertEqual(report["degenerate_triangles"], 0)
        self.assertAlmostEqual(report["covering_degree"], 1.0, places=9)
        self.assertEqual(report["centering_status"], "converged")
        self.assertLess(report["center_of_mass"], 1e-9)
        self.assertGreaterEqual(report["quasi_conformal_mean"], 1.0)
        self.assertLess(report["quasi_conformal_mean"], 1.2)
        for key in ("pole", "quasi_conformal_max", "area_distortion", "log_scale_min",
                    "log_scale_max", "max_radius_error", "centering_iterations", "seconds",
                    "seconds_total"):
            self.assertIn(key, report)
        self.assertIn("SphereResult(", repr(result))

    def test_pole_and_no_centering(self):
        positions, triangles = octahedron(2)
        result = ddg_tools.map_to_sphere(positions_buffer(positions), array.array('i', triangles),
                                         pole=7, center=False)
        self.assertEqual(result.report["pole"], 7)
        self.assertIsNone(result.report["centering_status"])
        self.assertIsNone(result.report["center_of_mass"])
        self.assertEqual(list(result.positions[21:24]), [0.0, 0.0, -1.0])

    def test_errors(self):
        positions, triangles = torus()
        with self.assertRaises(ddg_tools.MeshError) as caught:
            ddg_tools.map_to_sphere(array.array('d', positions), array.array('i', triangles))
        self.assertIn("genus 1", str(caught.exception))
        positions, triangles = octahedron(1)
        with self.assertRaises(ddg_tools.MeshError):
            # Drop one triangle: a disk, not a sphere.
            ddg_tools.map_to_sphere(array.array('d', positions), array.array('i', triangles[3:]))
        with self.assertRaises(ddg_tools.InputError) as caught:
            ddg_tools.map_to_sphere(array.array('d', positions), array.array('i', triangles),
                                    pole=10_000)
        self.assertTrue(str(caught.exception).startswith("pole"), str(caught.exception))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.map_to_sphere(array.array('d', positions), array.array('i', triangles),
                                    pole=-1)
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.map_to_sphere(array.array('d', positions[:-1]), array.array('i', triangles))
        with self.assertRaises(TypeError):
            ddg_tools.map_to_sphere(3, array.array('i', triangles))


def positions_buffer(positions):
    return array.array('f', positions)


if __name__ == "__main__":
    unittest.main(verbosity=2)
