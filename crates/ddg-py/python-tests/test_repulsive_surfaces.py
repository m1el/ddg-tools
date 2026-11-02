"""Tests of the repulsive-surfaces bindings of the built wheel of `ddg_tools`, without Blender.

Run through `scripts/py-test.sh` (see test_ddg_tools.py). Standard library only; the numpy tests
are skipped when numpy is not installed.
"""

import array
import ast
import math
import os
import tempfile
import threading
import time
import unittest

import ddg_tools

try:
    import numpy
except ImportError:  # pragma: no cover
    numpy = None

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------------------------
# Meshes
# ---------------------------------------------------------------------------------------------


def torus(n_major=24, n_minor=12, major=2.0, minor=0.7):
    """Outward oriented torus (genus 1)."""
    positions = array.array("d")
    for i in range(n_major):
        u = 2 * math.pi * i / n_major
        for j in range(n_minor):
            v = 2 * math.pi * j / n_minor
            r = major + minor * math.cos(v)
            positions += array.array("d", [r * math.cos(u), r * math.sin(u), minor * math.sin(v)])
    faces = array.array("i")
    for i in range(n_major):
        for j in range(n_minor):
            a = i * n_minor + j
            b = ((i + 1) % n_major) * n_minor + j
            c = ((i + 1) % n_major) * n_minor + (j + 1) % n_minor
            d = i * n_minor + (j + 1) % n_minor
            faces += array.array("i", [a, b, c, a, c, d])
    return positions, faces


def sphere(levels=2, radius=1.0, center=(0.0, 0.0, 0.0)):
    """Outward oriented icosphere-like sphere: an octahedron subdivided `levels` times."""
    verts = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    tris = [
        (0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
        (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5),
    ]
    for _ in range(levels):
        cache = {}

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                m = tuple((x + y) / 2 for x, y in zip(verts[a], verts[b]))
                n = math.sqrt(sum(x * x for x in m))
                verts.append(tuple(x / n for x in m))
                cache[key] = len(verts) - 1
            return cache[key]

        new = []
        for a, b, c in tris:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
        tris = new
    p = array.array("d")
    for v in verts:
        p += array.array("d", [center[k] + radius * v[k] for k in range(3)])
    f = array.array("i", [i for t in tris for i in t])
    return p, f


def volume(positions, faces):
    total = 0.0
    for k in range(0, len(faces), 3):
        a, b, c = (positions[3 * faces[k + i] : 3 * faces[k + i] + 3] for i in range(3))
        total += (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        ) / 6
    return total


def barycenter(positions):
    n = len(positions) // 3
    return [sum(positions[k::3]) / n for k in range(3)]


def write_obj(path, positions, faces):
    with open(path, "w", encoding="utf-8") as f:
        for k in range(0, len(positions), 3):
            f.write("v %.17g %.17g %.17g\n" % tuple(positions[k : k + 3]))
        for k in range(0, len(faces), 3):
            f.write("f %d %d %d\n" % tuple(i + 1 for i in faces[k : k + 3]))


SPHERE = sphere(2)
TORUS = torus(24, 12)
SMALL = {"threads": 2}


class Flow(unittest.TestCase):
    def test_torus_lowers_energy_keeps_barycenter_and_volume(self):
        p, f = TORUS
        e0, _, _ = ddg_tools.surface_energy(p, f)
        v0, b0 = volume(p, f), barycenter(p)
        q, g, report = ddg_tools.repulsive_surface_flow(
            p, f, steps=12, remesh=False,
            constraints=[{"type": "barycenter"}, {"type": "volume"}], **SMALL
        )
        e1, _, _ = ddg_tools.surface_energy(q, g)
        self.assertLess(e1, e0)
        self.assertAlmostEqual(report["energy_initial"], e0, delta=1e-6 * e0)
        self.assertEqual(report["iterations"], 12)
        self.assertEqual(len(report["steps"]), 12)
        self.assertAlmostEqual(volume(q, g), v0, delta=1e-6 * abs(v0))
        for a, b in zip(barycenter(q), b0):  # the constraint is the area-weighted barycenter
            self.assertAlmostEqual(a, b, delta=1e-3)
        for s in report["steps"]:
            self.assertLessEqual(s["energy_after"], s["energy_before"] * (1 + 1e-9))

    def test_sphere_area_growth_schedule(self):
        p, f = SPHERE
        area0 = 4 * math.pi
        q, g, _ = ddg_tools.repulsive_surface_flow(
            p, f, steps=6, remesh=False,
            constraints=[
                {"type": "barycenter"},
                {"type": "area", "factor": 1.2, "steps": 6},
            ],
            **SMALL,
        )
        self.assertEqual(len(q), len(p))
        area = 0.0
        for k in range(0, len(g), 3):
            a, b, c = (q[3 * g[k + i] : 3 * g[k + i] + 3] for i in range(3))
            u = [b[i] - a[i] for i in range(3)]
            w = [c[i] - a[i] for i in range(3)]
            cr = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
            area += 0.5 * math.sqrt(sum(x * x for x in cr))
        # The sphere mesh is an inscribed polyhedron: compare with its own start area.
        p0 = array.array("d", p)
        area_start = 0.0
        for k in range(0, len(f), 3):
            a, b, c = (p0[3 * f[k + i] : 3 * f[k + i] + 3] for i in range(3))
            u = [b[i] - a[i] for i in range(3)]
            w = [c[i] - a[i] for i in range(3)]
            cr = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
            area_start += 0.5 * math.sqrt(sum(x * x for x in cr))
        self.assertAlmostEqual(area, 1.2 * area_start, delta=1e-5 * area_start)
        self.assertLess(area_start, area0)

    def test_remeshing_changes_counts_and_the_snapshot_reflects_it(self):
        p, f = torus(40, 6, 2.0, 0.5)  # long thin quads: edges far from equal
        s = ddg_tools.RepulsiveSurfaceSession(p, f, **SMALL)
        snap0 = s.snapshot()
        self.assertEqual(list(snap0["faces"]), list(f))
        self.assertEqual(list(snap0["face_origin"]), list(range(len(f) // 3)))
        s.run(max_steps=4)
        snap = s.snapshot()
        nv, nf = len(snap["positions"]) // 3, len(snap["faces"]) // 3
        self.assertEqual((nv, nf), (s.n_vertices, s.n_faces))
        self.assertEqual(len(snap["face_origin"]), nf)
        self.assertNotEqual(nf, len(f) // 3)
        self.assertTrue(all(0 <= o < len(f) // 3 for o in snap["face_origin"]))
        self.assertTrue(all(0 <= v < nv for v in snap["faces"]))
        # Euler characteristic of the torus is kept by the remesher.
        edges = {tuple(sorted((snap["faces"][3 * k + i], snap["faces"][3 * k + (i + 1) % 3])))
                 for k in range(nf) for i in range(3)}
        self.assertEqual(nv - len(edges) + nf, 0)
        self.assertIsNotNone(snap["report"]["remesh"])
        self.assertGreater(
            snap["report"]["remesh"]["splits"] + snap["report"]["remesh"]["collapses"]
            + snap["report"]["remesh"]["flips"], 0)
        # Without remeshing the topology is untouched.
        q, g, _ = ddg_tools.repulsive_surface_flow(p, f, steps=3, remesh=False, **SMALL)
        self.assertEqual((len(q), list(g)), (len(p), list(f)))

    def test_pins_stay_through_remeshing(self):
        p, f = torus(40, 6, 2.0, 0.5)
        pins = [0, 7, 100]
        s = ddg_tools.RepulsiveSurfaceSession(p, f, pins=pins, **SMALL)
        before = {v: tuple(p[3 * v : 3 * v + 3]) for v in pins}
        s.run(max_steps=4)
        snap = s.snapshot()
        self.assertEqual(len(snap["pinned"]), 3)
        for old, new in zip(pins, snap["pinned"]):
            got = snap["positions"][3 * new : 3 * new + 3]
            for a, b in zip(before[old], got):
                self.assertAlmostEqual(a, b, delta=1e-6)
        self.assertNotEqual(s.n_faces, len(f) // 3)

    def test_moving_pins(self):
        p, f = SPHERE
        q, g, report = ddg_tools.repulsive_surface_flow(
            p, f, steps=5, remesh=False, pins=[0], pin_move=[0, 0, 0.1], pin_move_steps=5, **SMALL
        )
        self.assertAlmostEqual(q[2], p[2] + 0.1, delta=1e-6)
        self.assertEqual(list(report["pinned"]), [0])

    def test_obstacle_penalty_and_ccd(self):
        p, f = sphere(2, 1.0, (0.0, 0.0, 0.0))
        op, of = sphere(1, 0.5, (2.2, 0.0, 0.0))
        q, _, report = ddg_tools.repulsive_surface_flow(
            p, f, steps=3, remesh=False,
            penalties=[{"type": "mesh_obstacle", "vertices": op, "triangles": of, "weight": 1.0}],
            **SMALL,
        )
        self.assertEqual(len(report["steps"][0]["penalties_after"]), 1)
        self.assertGreater(report["steps"][0]["penalties_after"][0], 0.0)
        self.assertLess(max(q[0::3]), 2.0)  # still left of the obstacle


class Session(unittest.TestCase):
    def test_stepping_equals_the_one_call_flow_bitwise(self):
        p, f = TORUS
        for options in (
            {},
            {"constraints": [{"type": "barycenter"}, {"type": "volume", "factor": 1.1, "steps": 3}]},
            {"remesh": False, "metric": "h1"},
        ):
            q, g, report = ddg_tools.repulsive_surface_flow(p, f, steps=4, **SMALL, **options)
            s = ddg_tools.RepulsiveSurfaceSession(p, f, max_iterations=4, **SMALL, **options)
            for _ in range(4):
                s.step()
            self.assertEqual(s.positions.tobytes(), q.tobytes())
            self.assertEqual(s.faces.tobytes(), g.tobytes())
            snap = s.snapshot()
            self.assertEqual(snap["face_origin"].tobytes(), report["face_origin"].tobytes())
            self.assertEqual(snap["pinned"].tobytes(), report["pinned"].tobytes())
            # run() reaches the same state as the steps.
            s2 = ddg_tools.RepulsiveSurfaceSession(p, f, max_iterations=4, **SMALL, **options)
            end = s2.run()
            self.assertEqual(end["steps_taken"], 4)
            self.assertEqual(s2.positions.tobytes(), q.tobytes())

    def test_state_report_and_getters(self):
        p, f = SPHERE
        s = ddg_tools.RepulsiveSurfaceSession(p, f, **SMALL)
        self.assertEqual(s.iteration, 0)
        self.assertEqual(s.status, "running")
        self.assertIsNone(s.report)
        self.assertIsNone(s.error)
        r = s.step()
        for key in (
            "iteration", "status", "session_status", "n_vertices", "n_faces", "energy_before",
            "energy_after", "penalties_after", "gradient_norm", "tau", "halvings",
            "max_violation", "linear", "remesh", "seconds",
        ):
            self.assertIn(key, r)
        self.assertEqual(r["iteration"], 1)
        self.assertEqual(s.iteration, 1)
        self.assertEqual(s.energy, r["energy_after"])
        self.assertEqual(s.report["iteration"], 1)
        self.assertIn("RepulsiveSurfaceSession(", repr(s))
        snap = s.snapshot()
        self.assertEqual(
            set(snap),
            {"positions", "faces", "pinned", "face_origin", "iteration", "energy", "status",
             "error", "mean_edge_length0", "report"},
        )
        self.assertIsInstance(snap["positions"], array.array)
        self.assertEqual(snap["positions"].typecode, "d")
        self.assertEqual(snap["faces"].typecode, "i")

    def test_run_with_progress_and_limits(self):
        p, f = SPHERE
        s = ddg_tools.RepulsiveSurfaceSession(p, f, **SMALL)
        seen = []
        end = s.run(max_steps=5, progress=lambda r: seen.append(r["iteration"]) or len(seen) < 3)
        self.assertEqual(end["end"], "callback")
        self.assertEqual(seen, [1, 2, 3])
        end = s.run(max_steps=2)
        self.assertEqual((end["end"], end["steps_taken"], len(end["steps"])), ("max_steps", 2, 2))
        self.assertEqual(s.iteration, 5)
        s2 = ddg_tools.RepulsiveSurfaceSession(p, f, max_iterations=2, **SMALL)
        end = s2.run()
        self.assertEqual((end["end"], end["status"], end["steps_taken"]), ("final", "max_iterations", 2))
        self.assertEqual(s2.run()["steps_taken"], 0)
        with self.assertRaises(KeyError):  # an exception in progress() propagates
            ddg_tools.repulsive_surface_flow(
                p, f, steps=3, progress=lambda r: {}["boom"], **SMALL
            )

    def test_cancel_and_reset(self):
        p, f = SPHERE
        s = ddg_tools.RepulsiveSurfaceSession(p, f, **SMALL)
        s.step()
        s.cancel()
        self.assertEqual(s.status, "cancelled")
        r = s.step()
        self.assertEqual(r["status"], "cancelled")
        self.assertEqual(s.iteration, 1)
        self.assertEqual(s.run()["end"], "cancelled")
        s.reset()
        self.assertEqual((s.iteration, s.status), (0, "running"))
        self.assertEqual(s.positions.tobytes(), p.tobytes())
        s.step()
        self.assertEqual(s.iteration, 1)

    def test_cancel_from_another_thread_while_stepping(self):
        p, f = torus(60, 30)  # 3600 faces: a step takes a moment
        s = ddg_tools.RepulsiveSurfaceSession(p, f, max_iterations=100000)
        result = {}
        worker = threading.Thread(target=lambda: result.update(s.run()))
        worker.start()
        deadline = time.time() + 300
        ticks = 0
        while s.iteration < 2 and time.time() < deadline:
            ticks += 1  # the GIL is free: this thread polls without waiting for the step
            snap = s.snapshot()
            self.assertEqual(len(snap["faces"]) // 3, len(snap["face_origin"]))
            time.sleep(0.01)
        self.assertGreaterEqual(s.iteration, 2)
        self.assertGreater(ticks, 0)
        t0 = time.time()
        s.cancel()
        self.assertEqual(s.status, "cancelled")
        worker.join(120)
        self.assertFalse(worker.is_alive())
        self.assertLess(time.time() - t0, 60)
        self.assertEqual(result["end"], "cancelled")


class Defaults(unittest.TestCase):
    @staticmethod
    def stub_defaults(name, kind):
        path = os.path.join(HERE, "..", "ddg_tools.pyi")
        with open(path, encoding="utf-8") as fh:
            tree = ast.parse(fh.read())
        for node in tree.body:
            if kind == "class" and isinstance(node, ast.ClassDef) and node.name == name:
                fn = next(
                    i for i in node.body if isinstance(i, ast.FunctionDef) and i.name == "__init__"
                )
                break
            if kind == "func" and isinstance(node, ast.FunctionDef) and node.name == name:
                fn = node
                break
        else:
            raise AssertionError(f"{name} is not in the stub")
        a = fn.args
        return {
            arg.arg: ast.literal_eval(d)
            for arg, d in zip(a.kwonlyargs, a.kw_defaults)
            if d is not None
        }

    def test_defaults_are_the_papers_and_the_librarys(self):
        d = ddg_tools.repulsive_surface_defaults()
        self.assertEqual((d["alpha"], d["beta"]), (6.0, 12.0))
        self.assertEqual(d["metric"], "sobolev")
        self.assertTrue(d["low_order"])
        self.assertEqual(d["backend"], "auto")
        self.assertEqual(d["dense_max_faces"], 1000)
        self.assertEqual(d["energy_eval"], "barnes_hut")
        self.assertEqual((d["theta"], d["chi"]), (0.5, 0.5))
        self.assertTrue(d["remesh"])
        self.assertTrue(d["ccd"])
        self.assertEqual(d["ccd_gap"], 0.05)
        self.assertEqual(d["max_iterations"], 1000)
        self.assertEqual(d["method"], "projected_gradient")
        self.assertIsNone(d["constraints"])
        self.assertIsNone(d["penalties"])
        self.assertFalse(d["reference_compat"])
        self.assertEqual(d["remesh_ops"], ["split", "collapse", "flip", "smooth"])

    def test_stub_defaults_equal_the_library(self):
        d = ddg_tools.repulsive_surface_defaults()
        self.assertEqual(self.stub_defaults("RepulsiveSurfaceSession", "class"), d)
        self.assertEqual(self.stub_defaults("repulsive_surface_flow", "func"), d)
        self.assertEqual(self.stub_defaults("surface_energy", "func"), d)

    def test_passing_the_defaults_back_changes_nothing(self):
        p, f = SPHERE
        d = ddg_tools.repulsive_surface_defaults()
        a = ddg_tools.repulsive_surface_flow(p, f, steps=3, **SMALL)
        b = ddg_tools.repulsive_surface_flow(p, f, steps=3, **{**d, "threads": 2})
        self.assertEqual(a[0].tobytes(), b[0].tobytes())
        self.assertEqual(a[1].tobytes(), b[1].tobytes())


class Options(unittest.TestCase):
    """Every option class parses (a session is built, evaluating F once) and, where cheap, steps."""

    def build(self, mesh=SPHERE, step=False, **options):
        s = ddg_tools.RepulsiveSurfaceSession(*mesh, **{**SMALL, **options})
        if step:
            self.assertIn(s.step()["status"], ("accepted", "stuck"))
        return s

    def test_exponents_backend_metric_and_evaluation(self):
        self.build(step=True, alpha=4.0, beta=8.0)
        for metric in ("sobolev", "l2", "h1", "h2"):
            self.build(step=True, metric=metric)
        self.build(step=True, metric="l2", l2_kind="identity")
        self.build(step=True, low_order=False)
        for backend in ("auto", "dense", "hierarchical"):
            self.build(step=True, backend=backend)
        self.build(step=True, backend="hierarchical", chi=0.25, cluster_leaf_size=4,
                   gmres_restart=20, gmres_tolerance=1e-6, gmres_max_iterations=100,
                   gmres_side="left")
        self.build(step=True, energy_eval="exact")
        self.build(step=True, theta=0.25, leaf_size=4, tp_weight=2.0, dense_max_faces=10)
        self.build(step=True, energy="willmore", metric="h2")
        self.build(energy="coulomb", coulomb_alpha=2.0)
        self.build(step=True, reference_compat=True)
        self.build(step=True, threads=1)

    def test_constraints(self):
        self.build(step=True, constraints=[{"type": "barycenter", "target": [0, 0, 0]}])
        self.build(step=True, barycenter=False, pins=[0, 1, 2])
        self.build(step=True, constraints=[{"type": "pins", "vertices": [0, 1, 2]}])
        self.build(step=True, barycenter_mode="global", barycenter_shift=True)
        self.build(step=True, barycenter_mode="per_component")
        self.build(step=True, constraints=[{"type": "barycenter"}, {"type": "area", "target": 12.0}])
        self.build(step=True, constraints=[{"type": "barycenter"}, {"type": "volume"}],
                   remesh=False)
        self.build(step=True, constraints=[{"type": "barycenter"},
                                           {"type": "volume", "factor": 2.0, "offset": 0.1, "steps": 4}],
                   schedule_timing="after_projection")
        self.build(step=True, pins=[0], pin_move=[[0, 0, 0.1]], pin_move_steps=3)
        self.build(step=True, newton_steps=2)
        self.build(step=True, newton_max=3)

    def test_penalties_of_every_class(self):
        op, of = sphere(1, 0.3, (3, 0, 0))
        cloud = [3.0, 0.0, 0.0, 3.0, 1.0, 0.0]
        penalties = [
            {"type": "area_deviation", "weight": 0.1},
            {"type": "volume_deviation", "weight": 0.1, "target": 3.0},
            {"type": "total_area", "weight": 0.1, "power": 2.0},
            {"type": "total_volume", "weight": 0.1},
            {"type": "squared_error", "weight": 0.1},
            {"type": "boundary_length", "weight": 0.1, "target": 1.0},
            {"type": "boundary_curvature", "weight": 0.1, "dual_length": False},
            {"type": "willmore", "weight": 0.01},
            {"type": "mesh_obstacle", "vertices": op, "triangles": of, "weight": 0.1},
            {"type": "mesh_obstacle", "vertices": op, "triangles": of, "kind": "coulomb",
             "power": 2.0, "ccd": False},
            {"type": "point_cloud_obstacle", "points": cloud, "weight": 0.1},
        ]
        for surface in (
            {"type": "sphere", "center": [3, 0, 0], "radius": 0.5},
            {"type": "cylinder", "point": [0, 3, 0], "axis": [0, 0, 1], "radius": 0.5},
            {"type": "torus", "center": [0, 0, 3], "axis": [0, 0, 1], "major": 1.5, "minor": 0.4},
            {"type": "plane", "point": [0, 0, -3], "normal": [0, 0, 1]},
        ):
            penalties.append({"type": "implicit_obstacle", "surface": surface, "weight": 0.1,
                              "power": 2.0, "quadrature": "per_vertex"})
            penalties.append({"type": "implicit_attractor", "surface": surface, "weight": 0.1,
                              "vertices": [0, 1, 2]})
        s = self.build(penalties=penalties)
        # Every penalty has a value in the report of the first step.
        r = s.step()
        self.assertEqual(len(r["penalties_after"]), len(penalties))
        e, d, report = ddg_tools.surface_energy(*SPHERE, penalties=penalties)
        self.assertEqual([x["name"] for x in report["penalties"]][:3],
                         ["AreaDeviation", "VolumeDeviation", "TotalArea"])
        self.assertAlmostEqual(e, report["total"])
        self.assertAlmostEqual(report["total"], report["energy"] + sum(x["value"] for x in report["penalties"]))
        self.assertEqual(len(d), len(SPHERE[0]))

    def test_line_search_stop_remesh_and_limits(self):
        self.build(step=True, armijo="reference", armijo_c=1e-3)
        for start in ("grow_previous", "max_displacement", "reference"):
            self.build(step=True, start=start, start_lengths=0.5, start_factor=1.5)
        self.build(step=True, start="absolute", start_tau=0.01)
        self.build(step=True, tau_min_rel=1e-10, max_halvings=30)
        self.build(step=True, ccd=False)
        self.build(step=True, ccd=True, ccd_gap=0.1, ccd_max_iterations=50, ccd_close_factor=0.25)
        self.build(step=True, remesh=False)
        self.build(step=True, remesh_split=1.2, remesh_collapse=0.4, remesh_rounds=2,
                   remesh_rho=0.3, remesh_flip_cap=10, remesh_order="smooth_then_flip",
                   remesh_ops=["split", "smooth"])
        self.build(step=True, stop="energy_stall", stop_relative=1e-6, stop_window=5)
        self.build(step=True, stop="sobolev", stop_eps=1e-3)
        self.build(step=True, stop="off", max_iterations=7, time_limit=100.0, stuck_limit=3)

    def test_methods(self):
        self.build()  # projected gradient
        for method in ("aqp", "h1_lbfgs", "bcqn"):
            # The comparison optimizers parse; the library reports them as not available yet or
            # rejects the hard constraints: either way a ValueError that names them.
            with self.assertRaises(ValueError):
                self.build(method=method, aqp_kappa=50.0, lbfgs_memory=5, lbfgs_variant="reference")

    def test_array_and_numpy_inputs(self):
        p, f = SPHERE
        a = ddg_tools.repulsive_surface_flow(list(p), list(f), steps=1, **SMALL)
        b = ddg_tools.repulsive_surface_flow(p, f, steps=1, **SMALL)
        self.assertEqual(a[0].tobytes(), b[0].tobytes())
        nested = [list(p[k : k + 3]) for k in range(0, len(p), 3)]
        tri = [list(f[k : k + 3]) for k in range(0, len(f), 3)]
        c = ddg_tools.repulsive_surface_flow(nested, tri, steps=1, **SMALL)
        self.assertEqual(a[0].tobytes(), c[0].tobytes())

    @unittest.skipIf(numpy is None, "numpy is not installed")
    def test_numpy_inputs_and_option_arrays(self):
        p, f = SPHERE
        P = numpy.frombuffer(p, dtype=numpy.float64).reshape(-1, 3)
        F = numpy.frombuffer(f, dtype=numpy.int32).reshape(-1, 3)
        a = ddg_tools.repulsive_surface_flow(P, F, steps=2, **SMALL)
        b = ddg_tools.repulsive_surface_flow(p, f, steps=2, **SMALL)
        self.assertEqual(a[0].tobytes(), b[0].tobytes())
        op, of = sphere(1, 0.5, (3, 0, 0))
        OP = numpy.frombuffer(op, dtype=numpy.float64).reshape(-1, 3)
        OF = numpy.frombuffer(of, dtype=numpy.int32).reshape(-1, 3)
        pen = lambda v, t: [{"type": "mesh_obstacle", "vertices": v, "triangles": t}]
        a = ddg_tools.repulsive_surface_flow(P, F, steps=2, penalties=pen(OP, OF), **SMALL)
        b = ddg_tools.repulsive_surface_flow(p, f, steps=2, penalties=pen(op, of), **SMALL)
        self.assertEqual(a[0].tobytes(), b[0].tobytes())
        c = ddg_tools.repulsive_surface_flow(P, F, steps=2, pins=numpy.array([0, 1]), **SMALL)
        d = ddg_tools.repulsive_surface_flow(p, f, steps=2, pins=[0, 1], **SMALL)
        self.assertEqual(c[0].tobytes(), d[0].tobytes())


class Errors(unittest.TestCase):
    def assertInput(self, *args, message=None, **kwargs):
        with self.assertRaises(ddg_tools.RepulsiveSurfaceInputError) as cm:
            ddg_tools.RepulsiveSurfaceSession(*args, **{**SMALL, **kwargs})
        self.assertIsInstance(cm.exception, ValueError)
        self.assertIsInstance(cm.exception, ddg_tools.InputError)
        if message:
            self.assertIn(message, str(cm.exception))

    def test_unknown_option_is_named(self):
        self.assertInput(*SPHERE, message="alfa", alfa=3.0)
        with self.assertRaises(ValueError) as cm:
            ddg_tools.repulsive_surface_flow(*SPHERE, steps=1, remesh_splitt=1.0)
        self.assertIn("remesh_splitt", str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            ddg_tools.surface_energy(*SPHERE, bogus=1)
        self.assertIn("bogus", str(cm.exception))

    def test_bad_values(self):
        p = SPHERE
        self.assertInput(*p, message="alpha", alpha=0.5)
        self.assertInput(*p, message="beta", alpha=3.0, beta=3.0)
        self.assertInput(*p, message="metric", metric="l3")
        self.assertInput(*p, message="backend", backend="fast")
        self.assertInput(*p, message="split", remesh_split=-1.0)
        self.assertInput(*p, message="rho", remesh_rho=2.0)
        self.assertInput(*p, message="remesh_ops", remesh_ops=["split", "melt"])
        self.assertInput(*p, message="threads", threads=0)
        self.assertInput(*p, message="armijo", armijo_c=2.0)
        self.assertInput(*p, message="start_tau", start="absolute")
        self.assertInput(*p, message="stop_eps", stop="sobolev")
        self.assertInput(*p, message="coulomb_alpha", energy="coulomb")
        self.assertInput(*p, message="alpha", alpha="six")
        self.assertInput(*p, message="remesh", remesh=1.5)
        self.assertInput(*p, message="time_limit", time_limit=-1.0)
        self.assertInput(*p, message="max_iterations", max_iterations=-3)
        self.assertInput(*p, message="constraints", constraints=3)

    def test_bad_constraints_and_penalties(self):
        p = SPHERE
        self.assertInput(*p, message="constraints[0].type", constraints=[{"type": "knot"}])
        self.assertInput(*p, message="vertices", constraints=[{"type": "pins"}])
        self.assertInput(*p, message="unknown key 'colour'",
                         constraints=[{"type": "barycenter", "colour": 1}])
        self.assertInput(*p, message="steps", constraints=[{"type": "area", "factor": 2.0}])
        self.assertInput(*p, message="steps", pins=[0], pin_move=[0, 0, 1])
        self.assertInput(*p, message="pin", pin_move=[0, 0, 1], pin_move_steps=2)
        self.assertInput(*p, message="penalties[0].type", penalties=[{"type": "magic"}])
        self.assertInput(*p, message="unknown key 'wieght'",
                         penalties=[{"type": "willmore", "wieght": 1.0}])
        self.assertInput(*p, message="weight", penalties=[{"type": "willmore", "weight": -1.0}])
        self.assertInput(*p, message="surface",
                         penalties=[{"type": "implicit_obstacle"}])
        self.assertInput(*p, message="surface.type",
                         penalties=[{"type": "implicit_obstacle", "surface": {"type": "cube"}}])
        self.assertInput(*p, message="radius",
                         penalties=[{"type": "implicit_obstacle",
                                     "surface": {"type": "sphere", "radius": -1.0}}])
        self.assertInput(*p, message="triangles", penalties=[
            {"type": "mesh_obstacle", "vertices": [0, 0, 0, 1, 0, 0], "triangles": [0, 1]}])
        self.assertInput(*p, constraints=[{"type": "pins", "vertices": [10 ** 6]}])

    def test_bad_meshes(self):
        p, f = SPHERE
        self.assertInput(p, f[:-1])
        self.assertInput(p[:-1], f)
        self.assertInput(p, array.array("i", [0, 1, 99]))
        self.assertInput(p, array.array("i", [0, 0, 1]))
        bad = array.array("d", p)
        bad[4] = float("nan")
        self.assertInput(bad, f)
        # A non-manifold edge: three faces on one edge.
        self.assertInput(array.array("d", [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1]),
                         array.array("i", [0, 1, 2, 0, 1, 3, 0, 1, 4]))
        with self.assertRaises(ValueError):
            ddg_tools.surface_energy(p, f[:-1])
        with self.assertRaises(TypeError):
            ddg_tools.RepulsiveSurfaceSession(p)

    def test_session_methods_validate(self):
        s = ddg_tools.RepulsiveSurfaceSession(*SPHERE, **SMALL)
        with self.assertRaises(ValueError):
            s.run(max_steps=-1)
        with self.assertRaises(ValueError):
            ddg_tools.repulsive_surface_flow(*SPHERE, steps=-2)

    def test_volume_needs_a_closed_oriented_mesh(self):
        # A single triangle strip has a boundary: the volume constraint is refused.
        p = array.array("d", [0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0])
        f = array.array("i", [0, 1, 2, 1, 3, 2])
        with self.assertRaises(ValueError):
            ddg_tools.RepulsiveSurfaceSession(p, f, constraints=[{"type": "volume"}], **SMALL)


class Scenes(unittest.TestCase):
    def test_a_scene_file_runs(self):
        with tempfile.TemporaryDirectory() as d:
            p, f = SPHERE
            write_obj(os.path.join(d, "ball.obj"), p, f)
            op, of = sphere(1, 0.4, (2.5, 0.0, 0.0))
            write_obj(os.path.join(d, "wall.obj"), op, of)
            scene = os.path.join(d, "scene.txt")
            with open(scene, "w", encoding="utf-8") as fh:
                fh.write(
                    "repel_mesh ball.obj 6 12\n"
                    "constrain barycenter\n"
                    "constrain volume 1.2 5\n"
                    "obstacle wall.obj 1\n"
                    "implicit sphere repel 2 0.5 0.5 0 3 0\n"
                    "iteration_limit 3\n"
                )
            q, g, options = ddg_tools.repulsive_surface_scene(scene)
            self.assertEqual(q.tobytes(), p.tobytes())
            self.assertEqual(g.tobytes(), f.tobytes())
            self.assertEqual((options["alpha"], options["beta"]), (6.0, 12.0))
            self.assertEqual(options["max_iterations"], 3)
            kinds = [c["type"] for c in options["constraints"]]
            self.assertEqual(kinds, ["barycenter", "volume"])
            self.assertEqual(options["constraints"][1]["steps"], 5)
            self.assertEqual(len(options["penalties"]), 2)
            self.assertEqual(len(options["ccd_obstacles"]), 1)
            self.assertEqual(options["penalties"][0].__class__.__name__, "RepulsiveSurfaceHandle")
            # The options run the scene, and the pieces can be edited.
            out = ddg_tools.repulsive_surface_flow(q, g, **{**options, "threads": 2})
            self.assertEqual(out[2]["iterations"], 3)
            self.assertEqual(len(out[2]["steps"][0]["penalties_after"]), 2)
            out2 = ddg_tools.repulsive_surface_flow(q, g, **{**options, "threads": 2, "remesh": False,
                                                              "ccd": False, "ccd_obstacles": None})
            self.assertEqual(len(out2[0]), len(p))

    def test_scene_errors(self):
        with self.assertRaises(ValueError) as cm:
            ddg_tools.repulsive_surface_scene("/nonexistent/scene.txt")
        self.assertIn("scene.txt", str(cm.exception))
        with tempfile.TemporaryDirectory() as d:
            scene = os.path.join(d, "scene.txt")
            with open(scene, "w", encoding="utf-8") as fh:
                fh.write("repel_mesh ball.obj 6 12\nfrobnicate\n")
            with self.assertRaises(ValueError) as cm:
                ddg_tools.repulsive_surface_scene(scene)
            self.assertIn("frobnicate", str(cm.exception))

    def test_reference_scenes_when_present(self):
        root = os.path.normpath(os.path.join(HERE, "..", "..", "..", "data", "repulsive-surfaces", "scenes"))
        if not os.path.isdir(root):
            self.skipTest("data/repulsive-surfaces/scenes is not present (scripts/rs-fetch-scenes.sh)")
        loaded = []
        for name in sorted(os.listdir(root)):
            scene = os.path.join(root, name, "scene.txt")
            if not os.path.isfile(scene):
                continue
            try:
                q, g, options = ddg_tools.repulsive_surface_scene(scene)
            except ValueError:
                continue  # a scene whose mesh is not in the data set
            loaded.append((len(g), scene, q, g, options))
        if not loaded:
            self.skipTest("no scene could be loaded")
        loaded.sort(key=lambda t: t[0])
        n, scene, q, g, options = loaded[0]
        self.assertGreater(n // 3, 0)
        out = ddg_tools.repulsive_surface_flow(q, g, steps=1, **{**options, "threads": 4})
        self.assertEqual(out[2]["iterations"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
