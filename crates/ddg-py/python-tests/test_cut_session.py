"""Tests of `ddg_tools.CutSession` (the cut flow), run through `scripts/py-test.sh` next to
test_ddg_tools.py: the built wheel, without Blender, standard library only (numpy tests are
skipped when numpy is not installed).
"""

import array
import math
import sys
import unittest

import ddg_tools

try:
    import numpy
except ImportError:
    numpy = None


def uv_sphere(segments=24, rings=12, radius=1.0):
    """Blender's UV sphere (triangle fans at the poles, quads between), outward facing, fan
    triangulated: (positions, corner_vertices, polygon_starts, triangles) as lists."""
    positions = [0.0, 0.0, radius]
    for i in range(1, rings):
        theta = math.pi * i / rings
        for j in range(segments):
            phi = 2.0 * math.pi * j / segments
            positions += [radius * math.sin(theta) * math.cos(phi),
                          radius * math.sin(theta) * math.sin(phi),
                          radius * math.cos(theta)]
    positions += [0.0, 0.0, -radius]
    south = (len(positions) // 3) - 1

    def ring(i, j):
        return 1 + (i - 1) * segments + (j % segments)

    polygons = []
    for j in range(segments):
        polygons.append((0, ring(1, j), ring(1, j + 1)))
    for i in range(1, rings - 1):
        for j in range(segments):
            polygons.append((ring(i, j), ring(i + 1, j), ring(i + 1, j + 1), ring(i, j + 1)))
    for j in range(segments):
        polygons.append((south, ring(rings - 1, j + 1), ring(rings - 1, j)))
    corner_vertices = [v for polygon in polygons for v in polygon]
    polygon_starts, triangles, corner = [], [], 0
    for polygon in polygons:
        polygon_starts.append(corner)
        for k in range(1, len(polygon) - 1):
            triangles += [corner, corner + k, corner + k + 1]
        corner += len(polygon)
    return positions, corner_vertices, polygon_starts, triangles


def arrays(mesh):
    positions, corner_vertices, polygon_starts, triangles = mesh
    return (array.array('d', positions), array.array('i', corner_vertices),
            array.array('i', polygon_starts), array.array('i', triangles))


def stop_never(**extra):
    return dict(motion_tolerance=0.0, energy_tolerance=0.0, **extra)


class SessionTest(unittest.TestCase):
    def setUp(self):
        self.mesh = uv_sphere()
        self.buffers = arrays(self.mesh)
        self.n_vertices = len(self.mesh[0]) // 3
        self.n_polygons = len(self.mesh[2])

    def session(self, **kwargs):
        kwargs.setdefault("circle_vertex", 5)
        kwargs.setdefault("circle_radius", 0.9)
        return ddg_tools.CutSession(*self.buffers, **kwargs)

    def test_names_and_exceptions(self):
        for name in ("CutSession", "CutDescription", "InsertedMesh", "CutUnwrapResult",
                     "visibility_penalty", "length_weight_for_ripples", "NotAvailableError"):
            self.assertTrue(hasattr(ddg_tools, name), name)
        self.assertTrue(issubclass(ddg_tools.NotAvailableError, ddg_tools.DdgError))

    def test_sizes_and_state_shape(self):
        s = self.session()
        self.assertEqual(s.n_regions, 2)
        self.assertEqual(s.n_vertices, self.n_vertices)
        self.assertEqual(s.n_polygons, self.n_polygons)
        self.assertEqual(s.n_corners, len(self.mesh[1]))
        self.assertGreater(s.length_scale, 0.0)
        self.assertAlmostEqual(s.total_area * s.length_scale ** 2, 1.0, places=12)
        state = s.state
        self.assertEqual(state.shape, (2, s.n_corners))
        self.assertEqual(state.format, 'd')
        rows = state.tolist()
        self.assertEqual(rows[1], [-x for x in rows[0]])
        self.assertEqual(list(s.state_row(0)), rows[0])
        # Region 0 is the disk: negative at the centre vertex.
        centre = self.mesh[1].index(5)
        self.assertLess(rows[0][centre], 0.0)
        regions = s.corner_regions()
        self.assertEqual(regions[centre], 0)
        self.assertEqual(set(regions), {0, 1})

    @unittest.skipIf(numpy is None, "numpy is not installed")
    def test_state_is_a_2d_numpy_array(self):
        s = self.session()
        state = numpy.asarray(s.state)
        self.assertEqual(state.shape, (2, s.n_corners))
        s.state = state.astype(numpy.float32)
        self.assertLess(abs(numpy.asarray(s.state) - state).max(), 1e-6)

    def test_k_steps_twice_equal_2k_steps_once(self):
        once = self.session()
        once.run(length_weight=0.05, max_steps=8, **stop_never())
        twice = self.session()
        twice.run(length_weight=0.05, max_steps=4, **stop_never())
        twice.state = twice.state
        twice.run(length_weight=0.05, max_steps=4, **stop_never())
        a, b = once.state.tolist(), twice.state.tolist()
        worst = max(abs(x - y) for ra, rb in zip(a, b) for x, y in zip(ra, rb))
        self.assertLess(worst, 1e-12)
        # Through float32 storage (what a Blender attribute holds): equal to float32 accuracy.
        third = self.session()
        third.run(length_weight=0.05, max_steps=4, **stop_never())
        third.state = array.array('f', third.state_row(0).tolist() + third.state_row(1).tolist())
        third.run(length_weight=0.05, max_steps=4, **stop_never())
        c = third.state.tolist()
        worst = max(abs(x - y) for ra, rb in zip(a, c) for x, y in zip(ra, rb))
        self.assertLess(worst, 1e-3)

    def test_run_report(self):
        s = self.session()
        report = s.run(length_weight=0.05, max_steps=6, **stop_never())
        self.assertEqual(report["steps_taken"], 6)
        self.assertEqual(report["stop_reason"], "max_steps")
        self.assertLess(report["energy_after"], report["energy_before"])
        self.assertEqual(len(report["steps"]), 6)
        self.assertEqual(report["steps"][0]["energy_total"], report["energy_before"])
        for key in ("cut_length_before", "cut_length_after", "length_weight", "distortion_0",
                    "distortion_1", "area_0", "area_1", "time", "n_steps_total"):
            self.assertIn(key, report)
        self.assertAlmostEqual(report["area_0"] + report["area_1"], 1.0, places=9)
        self.assertEqual(s.n_steps, 6)
        self.assertGreater(s.time, 0.0)
        # cut_length is in the unit of the mesh: unit-area length divided by the scale.
        first = report["steps"][0]
        self.assertAlmostEqual(first["cut_length"] * s.length_scale, first["cut_length_unit_area"],
                               places=12)

    def test_step_report(self):
        s = self.session()
        report = s.step(length_weight=0.05)
        for key in ("time", "tau", "energy_total", "roughening_rate", "cut_motion",
                    "n_cut_points", "length_weight", "seconds_evaluate"):
            self.assertIn(key, report)
        self.assertTrue(math.isfinite(report["energy_total"]))
        self.assertEqual(report["length_weight"], 0.05)

    def test_stop_rule_stops_early(self):
        s = self.session()
        report = s.run(length_weight=0.05, max_steps=2000)
        self.assertLess(report["steps_taken"], 2000)
        self.assertIn(report["stop_reason"],
                      ("converged", "energy_stalled", "energy_rose", "cut_vanished"))

    def test_ripple_wavelength(self):
        s = self.session()
        omega = s.roughening_rate(energy="hencky")
        self.assertGreater(omega, 0.0)
        report = s.step(ripple_wavelength=0.25)
        self.assertAlmostEqual(report["length_weight"],
                               ddg_tools.length_weight_for_ripples(omega, 0.25), places=9)
        self.assertAlmostEqual(ddg_tools.length_weight_for_ripples(2.0, math.pi), 1.0)
        # Neither given: the default wavelength.
        again = self.session().step()
        self.assertAlmostEqual(again["length_weight"],
                               ddg_tools.length_weight_for_ripples(omega, 0.25), places=9)
        with self.assertRaisesRegex(ddg_tools.InputError, "exclude each other"):
            s.step(length_weight=1.0, ripple_wavelength=0.25)

    def test_every_control_runs(self):
        n = self.n_vertices
        # Mirror symmetry about the plane y = 0: ring i, segment j and -j form a class.
        segments, rings = 24, 12
        classes = array.array('i', [0] + [(i * segments + min(j, segments - j) % segments)
                                          for i in range(1, rings) for j in range(segments)] + [1])
        eta = array.array('d', [0.01 * (v % 3) for v in range(n)])
        cases = [
            dict(energy="dirichlet"),
            dict(energy="scaled_dirichlet"),
            dict(energy="scaled_dirichlet", feature_r_eps=0.5, feature_smoothing=0.05),
            dict(energy="scaled_dirichlet", scale_weights=array.array('d', [1.0] * n)),
            dict(energy="hencky", penalty=eta),
            dict(energy="hencky", penalty=eta, penalty_rule="paper_diagonal"),
            dict(distortion_weights=[1.0, 2.0], area_weights=[0.0, 0.5]),
            dict(symmetry=classes),
            dict(curvature_uniform=2.0),
            dict(curvature_per_triangle=array.array('d', [4.0] * (len(self.mesh[3]) // 3)),
                 triangulation="caller"),
            dict(stages="paper"),
            dict(flux_recovery="nodal", sliver_aspect_ratio=4.0, smoothing_passes=1),
            dict(flux_recovery="smoothed", flux_smoothing=0.7),
            dict(length_step="on_mesh", extension="harmonic", speed_scaling="unit",
                 redistance="keep_adjacent"),
            dict(redistance="off"),
            dict(step_rule="adaptive", time_scale=0.1, max_motion=0.2, tau_max=0.5),
            dict(step_rule="capped", tau_max=0.5, max_motion=0.2),
            dict(step_rule="fixed", tau=1e-4),
            dict(min_crossing_offset=0.01),
        ]
        for case in cases:
            with self.subTest(case=sorted(case)):
                # The symmetric start is the disk around the pole.
                kwargs = {"circle_vertex": 0} if "symmetry" in case else {}
                if "triangulation" in case:
                    kwargs["triangulation"] = case.pop("triangulation")
                s = self.session(**kwargs)
                case.setdefault("length_weight", 0.05)
                report = s.run(max_steps=2, **stop_never(), **case)
                self.assertEqual(report["steps_taken"], 2)
                self.assertTrue(math.isfinite(report["energy_after"]))

    def test_paper_stage_set_differs_from_default(self):
        a, b = self.session(), self.session()
        a.run(length_weight=0.05, max_steps=3, **stop_never())
        b.run(length_weight=0.05, max_steps=3, stages="paper", **stop_never())
        self.assertNotEqual(a.state.tolist(), b.state.tolist())

    def test_keyword_errors(self):
        s = self.session()
        with self.assertRaisesRegex(TypeError, "unexpected keyword argument 'bogus'"):
            s.step(bogus=1)
        with self.assertRaisesRegex(TypeError, "unexpected keyword argument 'max_steps'"):
            s.step(max_steps=3)
        with self.assertRaisesRegex(ddg_tools.InputError, "energy: \"bogus\" is not one of"):
            s.step(energy="bogus")
        with self.assertRaisesRegex(ddg_tools.InputError, "tau is required"):
            s.step(step_rule="fixed", length_weight=0.1)
        with self.assertRaisesRegex(ddg_tools.InputError, "scale_weights"):
            s.step(scale_weights=array.array('d', [1.0] * self.n_vertices), length_weight=0.1)
        with self.assertRaisesRegex(ddg_tools.InputError, "distortion_weights has 3 entries"):
            s.step(distortion_weights=[1.0, 1.0, 1.0], length_weight=0.1)
        with self.assertRaisesRegex(ddg_tools.InputError, "penalty has 3 entries"):
            s.step(penalty=array.array('d', [0.0] * 3), length_weight=0.1)
        with self.assertRaisesRegex(ddg_tools.InputError, "vertex"):
            s.step(penalty=array.array('d', [-1.0] * self.n_vertices), length_weight=0.1)
        with self.assertRaisesRegex(ddg_tools.InputError, "length weight"):
            s.step(length_weight=-1.0)
        with self.assertRaisesRegex(ddg_tools.InputError, "exclude each other"):
            s.step(curvature_uniform=1.0,
                   curvature_per_triangle=array.array('d', [0.0] * len(self.mesh[3])),
                   length_weight=0.1)

    def test_a_geodesic_circle_starts_two_regions_only(self):
        # The session takes more than two regions since the core has them; the circle start
        # does not (it has an inside and an outside).
        with self.assertRaisesRegex(ddg_tools.InputError, "two regions"):
            self.session(n_regions=3)
        with self.assertRaises(ddg_tools.InputError):
            self.session(n_regions=1)

    def test_starts(self):
        labels = array.array('i', [0 if p < self.n_polygons // 2 else 1
                                   for p in range(self.n_polygons)])
        s = self.session(circle_vertex=None, circle_radius=None, polygon_labels=labels)
        regions = s.corner_regions()
        self.assertEqual(set(regions), {0, 1})
        saved = s.state
        t = self.session(circle_vertex=None, circle_radius=None, corner_values=saved)
        self.assertEqual(t.state.tolist(), saved.tolist())
        s.restart(circle_vertex=100, circle_radius=0.5)
        self.assertEqual(s.n_steps, 0)
        s.restart(polygon_labels=labels)
        for bad in (dict(), dict(circle_vertex=1), dict(circle_vertex=1, circle_radius=1.0,
                                                          polygon_labels=labels)):
            with self.assertRaises(ddg_tools.InputError):
                s.restart(**bad)
        with self.assertRaisesRegex(ddg_tools.InputError, "radius"):
            s.restart(circle_vertex=1, circle_radius=-1.0)
        with self.assertRaisesRegex(ddg_tools.InputError, "out of range"):
            s.restart(circle_vertex=10 ** 6, circle_radius=1.0)
        with self.assertRaisesRegex(ddg_tools.InputError, "label"):
            s.restart(polygon_labels=array.array('i', [5] * self.n_polygons))

    def test_state_errors_name_the_problem(self):
        s = self.session()
        with self.assertRaisesRegex(ddg_tools.InputError, "state"):
            s.state = array.array('d', [0.0] * 7)
        flat = s.state_row(0).tolist() + s.state_row(1).tolist()
        flat[0] = float("nan")
        with self.assertRaisesRegex(ddg_tools.InputError, "is nan|NaN|not finite|nan"):
            s.state = array.array('d', flat)
        # Corners of one vertex that no seam separates must agree.
        flat = s.state_row(0).tolist() + s.state_row(1).tolist()
        corner = self.mesh[1].index(7, self.mesh[1].index(7) + 1)
        flat[corner] += 0.5
        with self.assertRaisesRegex(ddg_tools.InputError, "vertex 7"):
            s.state = array.array('d', flat)

    def test_preview(self):
        s = self.session()
        positions, starts, closed, regions = s.preview()
        self.assertEqual(positions.typecode, 'd')
        self.assertEqual(len(starts), len(closed) + 1)
        self.assertEqual(len(positions), 3 * starts[-1])
        self.assertEqual(list(regions), [0, 1])
        self.assertEqual(len(closed), 1)
        self.assertTrue(closed[0])
        with self.assertRaises(ddg_tools.InputError):
            s.preview(min_crossing_offset=0.7)

    def test_description_insert_and_unwrap(self):
        s = self.session()
        s.run(length_weight=0.05, max_steps=3, **stop_never())
        cut = s.describe_cut()
        self.assertGreater(cut.n_edge_points, 0)
        self.assertEqual(len(cut.edge_vertices), 2 * cut.n_edge_points)
        self.assertEqual(len(cut.edge_t), cut.n_edge_points)
        self.assertEqual(len(cut.edge_positions), 3 * cut.n_edge_points)
        self.assertEqual(len(cut.edge_side), cut.n_edge_points)
        self.assertTrue(all(0.0 < t < 1.0 for t in cut.edge_t))
        self.assertEqual(len(cut.path_start), cut.n_paths)
        self.assertEqual(len(cut.path_regions), 2 * cut.n_paths)
        self.assertEqual(len(cut.path_interior_starts), cut.n_paths + 1)
        n_interior = cut.path_interior_starts[-1]
        self.assertEqual(len(cut.interior_corners), 2 * n_interior)
        self.assertEqual(len(cut.interior_t), n_interior)
        self.assertEqual(len(cut.junction_positions), 0)
        polygons = list(cut.path_polygon)
        self.assertEqual(polygons, sorted(polygons))
        # A closed cut on a closed surface: every edge point is the end of two paths... or of
        # one path at each end.
        ends = list(cut.path_start) + list(cut.path_end)
        for point in range(cut.n_edge_points):
            self.assertEqual(ends.count(point), 2, point)

        inserted = s.insert(cut)
        n_new = len(inserted.positions) // 3
        self.assertEqual(n_new, self.n_vertices + cut.n_edge_points + n_interior)
        self.assertEqual(inserted.n_old_vertices, self.n_vertices)
        self.assertEqual(len(inserted.vertex_source), 3 * n_new)
        self.assertEqual(len(inserted.edge_point_vertex), cut.n_edge_points)
        self.assertEqual(len(inserted.polygon_parent), len(inserted.polygon_starts))
        self.assertGreater(len(inserted.polygon_starts), self.n_polygons)
        self.assertEqual(len(inserted.triangles),
                         3 * (len(inserted.corner_vertices) - 2 * len(inserted.polygon_starts)))
        seams = list(zip(inserted.new_seams[0::2], inserted.new_seams[1::2]))
        self.assertEqual(seams, sorted(set(seams)))
        self.assertTrue(all(a < b for a, b in seams))
        # Without an argument: the current cut, the same result.
        again = s.insert()
        self.assertEqual(again.new_seams.tolist(), inserted.new_seams.tolist())
        self.assertEqual(again.positions.tolist(), inserted.positions.tolist())
        # Every old polygon is untouched unless a path crosses it.
        crossed = set(polygons)
        parents = list(inserted.polygon_parent)
        for p in range(self.n_polygons):
            if p not in crossed:
                self.assertEqual(parents.count(p), 1, p)

        result = s.unwrap()
        self.assertEqual(len(result.uvs), 2 * len(result.inserted.corner_vertices))
        self.assertEqual(result.inserted.new_seams.tolist(), inserted.new_seams.tolist())
        self.assertEqual(len(result.polygon_charts), len(result.inserted.polygon_starts))
        self.assertEqual(max(result.polygon_charts) + 1, 2)
        self.assertEqual(len(result.cone_vertices), 0)
        self.assertTrue(all(math.isfinite(x) for x in result.uvs))
        self.assertIn("warnings", result.report)
        packed, info = ddg_tools.pack_charts(result.uvs, result.inserted.polygon_starts,
                                             result.polygon_charts)
        self.assertEqual(len(packed), len(result.uvs))
        self.assertEqual(info["n_charts"], 2)

    def test_empty_cut(self):
        s = self.session()
        s.state = array.array('d', [-1.0] * s.n_corners + [1.0] * s.n_corners)
        positions, starts, closed, regions = s.preview()
        self.assertEqual((len(positions), list(starts), list(closed), list(regions)), (0, [0], [], []))
        cut = s.describe_cut()
        self.assertEqual((cut.n_edge_points, cut.n_paths), (0, 0))
        self.assertEqual(len(s.insert(cut).positions), 3 * self.n_vertices)

    def test_seams_are_boundaries_of_the_flow(self):
        # A seam along one meridian: its two sides carry separate values.
        m = self.mesh
        segments, rings = 24, 12
        ring = lambda i, j: 1 + (i - 1) * segments + (j % segments)  # noqa: E731
        seam = []
        previous = 0
        for i in range(1, rings):
            seam += [previous, ring(i, 0)]
            previous = ring(i, 0)
        seam += [previous, len(m[0]) // 3 - 1]
        s = self.session(seam_edges=array.array('i', seam))
        self.assertEqual(len(s.seams_cut), len(seam))
        s.run(length_weight=0.05, max_steps=3, **stop_never())
        cut = s.describe_cut()
        # Points on a seam carry the side they are on.
        sided = [side for side in cut.edge_side if side >= 0]
        self.assertTrue(all(0 <= side < self.n_polygons for side in sided))
        inserted = s.insert(cut)
        self.assertGreater(len(inserted.new_seams), 0)

    def test_visibility_penalty(self):
        eta = ddg_tools.visibility_penalty(*self.buffers, (0.0, 0.0, 10.0))
        self.assertEqual(len(eta), self.n_vertices)
        self.assertEqual(eta.typecode, 'd')
        self.assertAlmostEqual(eta[0], 1.0, places=6)
        self.assertEqual(eta[-1], 0.0)
        self.assertTrue(all(0.0 <= x <= 1.0 for x in eta))
        none = ddg_tools.visibility_penalty(*self.buffers, (0.0, 0.0, 10.0),
                                            visibility=array.array('d', [0.0] * self.n_vertices))
        self.assertEqual(max(none), 0.0)
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.visibility_penalty(*self.buffers, (0.0, 0.0, 10.0),
                                         visibility=array.array('d', [1.0]))
        # As the penalty of a step.
        s = self.session()
        report = s.run(penalty=array.array('d', [0.05 * x for x in eta]), length_weight=0.05,
                       max_steps=2, **stop_never())
        self.assertEqual(report["steps_taken"], 2)

    def test_bad_meshes(self):
        positions = array.array('d', [0, 0, 0, 1, 0, 0, 0.5, 1, 0, 0.5, -1, 0, 0.5, 0, 1])
        fins = (positions, array.array('i', [0, 1, 2, 1, 0, 3, 0, 1, 4]),
                array.array('i', [0, 3, 6]), array.array('i', range(9)))
        with self.assertRaises(ddg_tools.MeshError):
            ddg_tools.CutSession(*fins, circle_vertex=0, circle_radius=0.1)
        bad = (self.buffers[0], self.buffers[1], self.buffers[2], array.array('i', [0, 1]))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.CutSession(*bad, circle_vertex=0, circle_radius=0.1)
        with self.assertRaises(TypeError):
            ddg_tools.CutSession(1, 2, 3, 4, circle_vertex=0, circle_radius=0.1)


class RegionsTest(unittest.TestCase):
    """More than two regions: the normal-cluster start, N-valued labels and corner values, the
    recommended stages, per-region weights, the network on polygons and its insertion."""

    def setUp(self):
        self.mesh = uv_sphere(32, 16)
        self.buffers = arrays(self.mesh)
        self.n_polygons = len(self.mesh[2])

    def clusters(self, n=4, **kwargs):
        return ddg_tools.CutSession(*self.buffers, n_regions=n, cluster_iterations=30, **kwargs)

    def test_normal_clusters_start(self):
        for n in (3, 4, 6):
            s = self.clusters(n)
            self.assertEqual(s.n_regions, n)
            self.assertEqual(s.state.shape, (n, s.n_corners))
            self.assertGreater(len(set(s.corner_regions())), 2)
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.CutSession(*self.buffers, circle_vertex=1, circle_radius=0.5,
                                 cluster_iterations=5)

    def test_labels_and_corner_values_with_n_regions(self):
        labels = array.array('i', [p % 3 for p in range(self.n_polygons)])
        s = ddg_tools.CutSession(*self.buffers, n_regions=3, polygon_labels=labels)
        self.assertEqual(s.state.shape, (3, s.n_corners))
        t = ddg_tools.CutSession(*self.buffers, n_regions=3, corner_values=s.state)
        self.assertEqual(t.state.tolist(), s.state.tolist())
        with self.assertRaisesRegex(ddg_tools.InputError, "label"):
            ddg_tools.CutSession(*self.buffers, n_regions=3,
                                 polygon_labels=array.array('i', [3] * self.n_polygons))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.CutSession(*self.buffers, n_regions=3,
                                 corner_values=array.array('d', [0.0] * 10))
        # NotAvailableError is kept as a class but not raised for this any more.
        self.assertTrue(issubclass(ddg_tools.NotAvailableError, ddg_tools.DdgError))

    def test_run_reports_every_region_and_uses_recommended_stages(self):
        s = self.clusters(4)
        report = s.run(length_weight=0.5, max_steps=4, motion_tolerance=0.0,
                       energy_tolerance=0.0)
        self.assertEqual(report["steps_taken"], 4)
        for k in range(4):
            self.assertIn("distortion_%d" % k, report)
            self.assertIn("area_%d" % k, report)
        self.assertAlmostEqual(sum(report["area_%d" % k] for k in range(4)), 1.0, places=9)
        # The recommended stages are the paper's with the flux recovery: not the default's.
        a, b = self.clusters(4), self.clusters(4)
        kw = dict(length_weight=0.5, max_steps=3, motion_tolerance=0.0, energy_tolerance=0.0)
        a.run(**kw)
        b.run(stages="default", **kw)
        self.assertNotEqual(a.state.tolist(), b.state.tolist())
        c = self.clusters(4)
        c.run(stages="recommended", **kw)
        self.assertEqual(a.state.tolist(), c.state.tolist())
        # Individual overrides still work.
        d = self.clusters(4)
        d.run(flux_recovery="nodal", **kw)

    def test_per_region_weights_have_length_n(self):
        s = self.clusters(4)
        kw = dict(length_weight=0.5, max_steps=2, motion_tolerance=0.0, energy_tolerance=0.0)
        s.run(distortion_weights=[1.0, 2.0, 0.5, 1.0], area_weights=[0.0, 0.1, 0.0, 0.0], **kw)
        with self.assertRaisesRegex(ddg_tools.InputError, "distortion_weights has 2 entries"):
            s.run(distortion_weights=[1.0, 1.0], **kw)

    def test_network_on_polygons_insert_and_unwrap(self):
        s = self.clusters(4)
        s.run(length_weight=0.5, max_steps=3, motion_tolerance=0.0, energy_tolerance=0.0)
        positions, starts, closed, regions = s.preview()
        self.assertEqual(len(regions), 2 * len(closed))
        self.assertTrue(all(0 <= r < 4 for r in regions))
        self.assertTrue(all(regions[2 * i] < regions[2 * i + 1] for i in range(len(closed))))
        cut = s.describe_cut()
        n_j = len(cut.junction_polygon)
        self.assertGreater(n_j, 0)
        self.assertEqual(len(cut.junction_positions), 3 * n_j)
        self.assertEqual(len(cut.junction_triangle), n_j)
        self.assertEqual(len(cut.junction_corners), 3 * n_j)
        self.assertEqual(len(cut.junction_weights), 3 * n_j)
        for j in range(n_j):
            w = cut.junction_weights[3 * j:3 * j + 3]
            self.assertAlmostEqual(sum(w), 1.0, places=9)
            corners = cut.junction_corners[3 * j:3 * j + 3]
            xs = [self.mesh[0][3 * self.mesh[1][c]:3 * self.mesh[1][c] + 3] for c in corners]
            for k in range(3):
                x = sum(w[i] * xs[i][k] for i in range(3))
                self.assertAlmostEqual(x, cut.junction_positions[3 * j + k], places=6)
        ends = list(cut.path_start) + list(cut.path_end)
        for j in range(n_j):
            self.assertEqual(ends.count(-(j + 1)), 3, j)
        inserted = s.insert(cut)
        self.assertEqual(len(inserted.junction_vertex), n_j)
        kinds = [inserted.vertex_source[3 * v] for v in range(len(inserted.positions) // 3)]
        self.assertEqual(kinds.count(3), n_j)
        for j, v in enumerate(inserted.junction_vertex):
            self.assertEqual(inserted.vertex_source[3 * v:3 * v + 3].tolist(), [3, j, 0])
        result = s.unwrap()
        self.assertTrue(all(math.isfinite(x) for x in result.uvs))
        self.assertEqual(len(result.uvs), 2 * len(result.inserted.corner_vertices))
        self.assertGreaterEqual(max(result.polygon_charts) + 1, 4)

    def test_two_region_description_has_no_junctions(self):
        s = ddg_tools.CutSession(*self.buffers, circle_vertex=5, circle_radius=0.9)
        s.run(length_weight=0.05, max_steps=2, motion_tolerance=0.0, energy_tolerance=0.0)
        cut = s.describe_cut()
        self.assertEqual(len(cut.junction_polygon), 0)
        self.assertEqual(len(s.insert(cut).junction_vertex), 0)



class IntrinsicTest(unittest.TestCase):
    """The triangulation the flow runs on: intrinsic (Delaunay flips, the default), refined, the
    caller's; the report, the state bytes, the stale-state error, the description's lines."""

    def setUp(self):
        self.mesh = uv_sphere(24, 12)
        self.buffers = arrays(self.mesh)

    def make(self, **kwargs):
        kwargs.setdefault("circle_vertex", 5)
        kwargs.setdefault("circle_radius", 0.9)
        return ddg_tools.CutSession(*self.buffers, **kwargs)

    @staticmethod
    def worst(a, b):
        return max(abs(x - y) for r, q in zip(a.state.tolist(), b.state.tolist())
                   for x, y in zip(r, q))

    def test_choices_and_report(self):
        default = self.make()
        self.assertEqual(default.triangulation, "intrinsic")
        report = default.intrinsic_report
        for key in ("n_added", "n_initial_flips", "min_angle_degrees", "max_angle_degrees",
                    "is_delaunay", "fingerprint", "seconds"):
            self.assertIn(key, report)
        self.assertTrue(report["is_delaunay"])
        self.assertEqual(report["n_added"], 0)
        self.assertTrue(default.state_is_per_corner)
        self.assertEqual(default.fingerprint, report["fingerprint"])
        caller = self.make(triangulation="caller")
        self.assertEqual((caller.triangulation, caller.intrinsic_report), ("caller", None))
        self.assertTrue(caller.state_is_per_corner)
        refined = self.make(triangulation="refined", refine_bound=30.0)
        self.assertEqual(refined.triangulation, "refined")
        self.assertGreater(refined.intrinsic_report["n_added"], 0)
        self.assertFalse(refined.state_is_per_corner)
        self.assertEqual(self.make(triangulation="refined").intrinsic_report["bound_degrees"],
                         25.0)
        for bad in (dict(triangulation="bogus"), dict(refine_bound=20.0),
                    dict(triangulation="refined", refine_bound=80.0)):
            with self.assertRaises(ddg_tools.InputError, msg=str(bad)):
                self.make(**bad)

    def test_optimal_delaunay_flag(self):
        off = self.make()
        self.assertIs(off.optimal_delaunay, False)
        self.assertEqual(off.intrinsic_report["n_optimized"], 0)
        self.assertNotIn("optimal_delaunay_status", off.intrinsic_report)
        self.assertIs(self.make(triangulation="caller").optimal_delaunay, False)
        kw = dict(length_weight=0.05, max_steps=3, motion_tolerance=0.0, energy_tolerance=0.0)
        for name in ("intrinsic", "refined"):
            with self.subTest(triangulation=name):
                s = self.make(triangulation=name, optimal_delaunay=True)
                self.assertIs(s.optimal_delaunay, True)
                self.assertEqual(s.triangulation, name)  # the names are unchanged
                report = s.intrinsic_report
                self.assertIn(report["optimal_delaunay_status"],
                              ("Converged", "IterationLimit", "InsertionLimit"))
                self.assertGreaterEqual(report["n_optimized"], 0)
                self.assertGreater(report["optimal_delaunay_iterations"], 0)
                self.assertTrue(math.isfinite(report["optimal_delaunay_last_motion"]))
                result = s.run(**kw)
                self.assertLess(result["energy_after"], result["energy_before"])
                bytes_ = s.snapshot()
                again = ddg_tools.CutSession(*self.buffers, state=bytes_)
                self.assertIs(again.optimal_delaunay, True)
                self.assertEqual(again.triangulation, name)
                self.assertEqual(again.snapshot(), bytes_)
                self.assertEqual(again.fingerprint, s.fingerprint)
                # A state made without the flag does not carry it.
                plain = self.make(triangulation=name)
                self.assertIs(ddg_tools.CutSession(
                    *self.buffers, state=plain.snapshot()).optimal_delaunay, False)
        with self.assertRaises(ddg_tools.InputError) as ctx:
            self.make(triangulation="caller", optimal_delaunay=True)
        self.assertIn("intrinsic triangulation", str(ctx.exception))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.CutSession(*self.buffers, state=off.snapshot(), optimal_delaunay=True)

    def test_the_old_expectations_hold_on_every_triangulation(self):
        kw = dict(length_weight=0.05, max_steps=6, motion_tolerance=0.0, energy_tolerance=0.0)
        for name in ("caller", "intrinsic", "refined"):
            with self.subTest(triangulation=name):
                s = self.make(triangulation=name)
                report = s.run(**kw)
                self.assertLess(report["energy_after"], report["energy_before"])
                self.assertEqual(s.state.shape, (2, s.n_corners))
                cut = s.describe_cut()
                self.assertGreater(cut.n_paths, 0)
                inserted = s.insert(cut)
                self.assertGreater(len(inserted.polygon_starts), len(self.mesh[2]))
                self.assertEqual(len(inserted.polygon_region), len(inserted.polygon_starts))
                positions, starts, closed, regions = s.preview()
                self.assertEqual(len(regions), 2 * len(closed))
                positions2, starts2, closed2, regions2 = cut.lines()
                self.assertEqual((len(closed2), sorted(regions2)), (len(closed), sorted(regions)))
                self.assertEqual(s.lines()[1].tolist(), starts2.tolist())
                result = s.unwrap()
                self.assertTrue(all(math.isfinite(x) for x in result.uvs))

    def test_k_steps_twice_equal_2k_steps_once_with_added_vertices(self):
        kw = dict(length_weight=0.05, motion_tolerance=0.0, energy_tolerance=0.0)
        once = self.make(triangulation="refined", refine_bound=30.0)
        once.run(max_steps=6, **kw)
        twice = self.make(triangulation="refined", refine_bound=30.0)
        twice.run(max_steps=3, **kw)
        bytes_ = twice.snapshot()
        self.assertIsInstance(bytes_, bytes)
        # Through the bytes into a new session on the same mesh, as the Blender operators do.
        again = ddg_tools.CutSession(*self.buffers, state=bytes_)
        self.assertEqual(again.triangulation, "refined")
        self.assertEqual(again.snapshot(), bytes_)
        again.run(max_steps=3, **kw)
        self.assertLess(self.worst(once, again), 1e-12)
        # restore() on a live session.
        twice.run(max_steps=3, **kw)
        other = self.make(triangulation="refined", refine_bound=30.0)
        other.restore(twice.snapshot())
        self.assertLess(self.worst(other, once), 1e-12)

    def test_stale_state(self):
        s = self.make(triangulation="refined", refine_bound=30.0)
        bytes_ = s.snapshot()
        # Another seam: the triangulation is another one.
        seam = array.array('i', [0, 1])
        with self.assertRaises(ddg_tools.StaleStateError):
            ddg_tools.CutSession(*self.buffers, state=bytes_, seam_edges=seam)
        other = ddg_tools.CutSession(*arrays(uv_sphere(20, 10)), circle_vertex=1,
                                     circle_radius=0.5)
        with self.assertRaises(ddg_tools.StaleStateError):
            other.restore(bytes_)
        self.assertTrue(issubclass(ddg_tools.StaleStateError, ddg_tools.DdgError))
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.CutSession(*self.buffers, state=b"not a state")
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.CutSession(*self.buffers, state=bytes_, circle_vertex=1, circle_radius=1.0)
        with self.assertRaisesRegex(ddg_tools.InputError, "n_regions"):
            ddg_tools.CutSession(*self.buffers, state=bytes_, n_regions=3)

    def test_merge_tolerance_and_the_description_report(self):
        kw = dict(length_weight=0.05, max_steps=6, motion_tolerance=0.0, energy_tolerance=0.0)
        merged = self.make()
        merged.run(**kw)
        report = merged.describe_cut().report
        for key in ("n_merged_nodes", "max_displacement", "max_displacement_relative",
                    "n_moved_points", "n_removed_excursions"):
            self.assertIn(key, report)
        self.assertLessEqual(report["max_displacement_relative"], 0.1 + 1e-12)
        raw = self.make(merge_tolerance=0.0)
        raw.run(**kw)
        r0 = raw.describe_cut().report
        self.assertEqual(r0["n_merged_nodes"], 0)
        self.assertGreaterEqual(report["n_merged_nodes"], 0)
        self.assertEqual(self.make(triangulation="caller").describe_cut().report["n_merged_nodes"], 0)
        # The tolerance can be changed when a session is rebuilt from the state bytes.
        again = ddg_tools.CutSession(*self.buffers, state=merged.snapshot(), merge_tolerance=0.0)
        self.assertEqual(again.describe_cut().report["n_merged_nodes"], 0)
        with self.assertRaises(ddg_tools.InputError):
            self.make(merge_tolerance=0.7)

    def test_patch_count_and_the_changed_polygons_of_an_insertion(self):
        for name in ("caller", "intrinsic", "refined"):
            with self.subTest(triangulation=name):
                s = self.make(triangulation=name)
                self.assertEqual(s.n_patches, 2)
                s.run(length_weight=0.05, max_steps=4, motion_tolerance=0.0,
                      energy_tolerance=0.0)
                self.assertEqual(s.n_patches, 2)
                ins = s.insert()
                n = len(ins.changed_polygons)
                self.assertGreater(n, 0)
                self.assertEqual(len(ins.changed_vertex_starts), n + 1)
                self.assertEqual(len(ins.changed_triangle_starts), n + 1)
                self.assertEqual(len(ins.changed_vertices), ins.changed_vertex_starts[-1])
                self.assertEqual(len(ins.changed_triangles), 3 * ins.changed_triangle_starts[-1])
                for k in range(n):
                    vertices = ins.changed_vertices[ins.changed_vertex_starts[k]:
                                                    ins.changed_vertex_starts[k + 1]]
                    triangles = ins.changed_triangles[3 * ins.changed_triangle_starts[k]:
                                                      3 * ins.changed_triangle_starts[k + 1]]
                    self.assertGreater(len(vertices), 3)
                    self.assertEqual(len(triangles), 3 * (len(vertices) - 2))
                    self.assertTrue(set(triangles) <= set(vertices))
        # Several regions: the clusters are more than two patches.
        many = ddg_tools.CutSession(*self.buffers, n_regions=4, cluster_iterations=30)
        self.assertGreaterEqual(many.n_patches, 4)

    def test_curvature_per_triangle_is_refused_on_an_intrinsic_triangulation(self):
        s = self.make()
        k = array.array('d', [0.0] * (len(self.mesh[3]) // 3))
        with self.assertRaises(ddg_tools.InputError):
            s.step(length_weight=0.1, curvature_per_triangle=k)
        self.make(triangulation="caller").step(length_weight=0.1, curvature_per_triangle=k)


if __name__ == "__main__":
    print("python", sys.version.split()[0], "| numpy",
          numpy.__version__ if numpy is not None else "absent",
          "| ddg_tools", ddg_tools.__version__, "from", ddg_tools.__file__)
    unittest.main(verbosity=2)
