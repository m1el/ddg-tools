"""Tests of `ddg_tools.smoothest_direction_field` and `ddg_tools.adaptive_solve` (tangent fields
and adaptive refinement on intrinsic triangulations, crate ddg-intrinsic-apps), without Blender.

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


def sheared_grid(n=8, shear=1.5):
    """A flat parallelogram of obtuse triangles (a disk): (positions, triangles)."""
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


class DirectionFieldTest(unittest.TestCase):
    def test_index_sum_is_n_times_euler_characteristic(self):
        positions, triangles = octahedron()
        for n in (1, 2, 4):
            for triangulation in ("delaunay", "refined", "input"):
                for laplacian in ("fem", "fd"):
                    with self.subTest(n=n, triangulation=triangulation, laplacian=laplacian):
                        r = ddg_tools.smoothest_direction_field(
                            positions, triangles, n=n, triangulation=triangulation,
                            laplacian=laplacian)
                        self.assertEqual(r.report["euler_characteristic"], 2)
                        self.assertEqual(r.report["index_sum"], 2 * n)
                        self.assertEqual(sum(r.singularity_indices), 2 * n)
                        self.assertEqual(len(r.singularity_points),
                                         3 * len(r.singularity_indices))
                        self.assertGreaterEqual(r.report["n_triangulation_vertices"], 6)
                        self.assertIn(r.report["status"], ("converged", "iteration_limit"))

    def test_vectors_are_unit_and_typed(self):
        positions, triangles = octahedron()
        r = ddg_tools.smoothest_direction_field(positions, triangles, n=2)
        v = r.vectors
        self.assertIsInstance(v, array.array)
        self.assertEqual(v.typecode, 'd')
        self.assertEqual(len(v), 18)
        self.assertEqual(r.singularity_indices.typecode, 'i')
        for k in range(6):
            self.assertAlmostEqual(math.hypot(v[3 * k], v[3 * k + 1], v[3 * k + 2]), 1.0,
                                   places=9)
        self.assertGreater(r.report["energy"], 0.0)

    def test_flat_disk_has_a_constant_field(self):
        # A flat disk has a parallel field (energy 0): the smoothest field is constant.
        positions, triangles = sheared_grid()
        r = ddg_tools.smoothest_direction_field(positions, triangles, n=1, laplacian="fd")
        self.assertLess(abs(r.report["energy"]), 1e-8)
        v = r.vectors
        for k in range(len(v) // 3):
            self.assertAlmostEqual(v[3 * k] * v[0] + v[3 * k + 1] * v[1], 1.0, places=6)
            self.assertAlmostEqual(v[3 * k + 2], 0.0, places=12)

    def test_default_representation_is_integer(self):
        positions, triangles = octahedron()
        r = ddg_tools.smoothest_direction_field(positions, triangles)
        self.assertEqual(r.report["representation"], "integer")

    def test_both_representations_give_the_same_field(self):
        for mesh in (octahedron(), sheared_grid()):
            positions, triangles = mesh
            for triangulation in ("delaunay", "refined"):
                for n in (1, 4):
                    with self.subTest(triangulation=triangulation, n=n):
                        a, b = (ddg_tools.smoothest_direction_field(
                            positions, triangles, n=n, triangulation=triangulation,
                            representation=representation)
                            for representation in ("signposts", "integer"))
                        self.assertEqual(a.report["representation"], "signposts")
                        self.assertEqual(b.report["representation"], "integer")
                        self.assertEqual(list(a.singularity_indices),
                                         list(b.singularity_indices))
                        self.assertAlmostEqual(a.report["eigenvalue"], b.report["eigenvalue"],
                                               delta=1e-9 * max(1.0, abs(a.report["eigenvalue"])))
                        for x, y in zip(a.singularity_points, b.singularity_points):
                            self.assertAlmostEqual(x, y, places=9)

    def test_bad_options(self):
        positions, triangles = octahedron()
        for kwargs, name in (({"n": 0}, "n"), ({"triangulation": "nope"}, "triangulation"),
                             ({"laplacian": "spectral"}, "laplacian"),
                             ({"triangulation": "refined", "min_angle_degrees": 50.0},
                              "min_angle_degrees"),
                             ({"representation": "halfedges"}, "representation")):
            with self.subTest(**kwargs):
                with self.assertRaises(ddg_tools.InputError) as cm:
                    ddg_tools.smoothest_direction_field(positions, triangles, **kwargs)
                self.assertTrue(str(cm.exception).startswith(name), str(cm.exception))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.smoothest_direction_field(positions, array.array('i', [0, 1]))


class AdaptiveSolveTest(unittest.TestCase):
    def test_heat_kernel(self):
        positions, triangles = octahedron()
        r = ddg_tools.adaptive_solve(positions, triangles, "heat", source=4, time=0.05,
                                     max_vertices=800)
        self.assertEqual(r.n_components, 1)
        u = r.values
        self.assertIsInstance(u, array.array)
        self.assertEqual(len(u), 6)
        self.assertEqual(max(range(6), key=lambda v: u[v]), 4)
        self.assertEqual(min(range(6), key=lambda v: u[v]), 5)
        self.assertGreater(r.report["n_triangulation_vertices"], 6)
        self.assertIn(r.report["status"],
                      ("tolerance", "vertex_limit", "iteration_limit", "stalled"))
        self.assertGreater(r.report["iterations"], 0)

    def test_green_function_peaks_at_the_source(self):
        positions, triangles = octahedron()
        r = ddg_tools.adaptive_solve(positions, triangles, "green", source=0, max_vertices=600)
        u = r.values
        self.assertEqual(max(range(6), key=lambda v: u[v]), 0)
        self.assertEqual(min(range(6), key=lambda v: u[v]), 1)

    def test_disk_map_sends_the_boundary_to_the_circle(self):
        positions, triangles = sheared_grid()
        r = ddg_tools.adaptive_solve(positions, triangles, "disk", max_vertices=600)
        self.assertEqual(r.n_components, 2)
        uv = r.values
        n = 9
        self.assertEqual(len(uv), 2 * n * n)
        for j in range(n):
            for i in range(n):
                v = j * n + i
                radius = math.hypot(uv[2 * v], uv[2 * v + 1])
                if i in (0, n - 1) or j in (0, n - 1):
                    self.assertAlmostEqual(radius, 1.0, places=9)
                else:
                    self.assertLess(radius, 1.0)

    def test_errors(self):
        positions, triangles = octahedron()
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.adaptive_solve(positions, triangles, "heat", source=99)
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.adaptive_solve(positions, triangles, "wave")
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.adaptive_solve(positions, triangles, "green", bulk=0.0)
        with self.assertRaises(ddg_tools.MeshError):
            ddg_tools.adaptive_solve(positions, triangles, "disk")


if __name__ == "__main__":
    print("python", sys.version.split()[0], "| ddg_tools", ddg_tools.__version__,
          "from", ddg_tools.__file__)
    unittest.main(verbosity=2)
