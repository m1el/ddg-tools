"""Tests of the repulsive-curves bindings of the built wheel of `ddg_tools`, without Blender.

Run through `scripts/py-test.sh` (see test_ddg_tools.py). Standard library only; the numpy tests
are skipped when numpy is not installed.
"""

import ast
import array
import math
import os
import threading
import time
import unittest

import ddg_tools

try:
    import numpy
except ImportError:  # pragma: no cover
    numpy = None


def trefoil(n=64):
    positions, edges = [], []
    for i in range(n):
        t = 2 * math.pi * i / n
        positions += [
            math.sin(t) + 2 * math.sin(2 * t),
            math.cos(t) - 2 * math.cos(2 * t),
            -math.sin(3 * t),
        ]
        edges += [i, (i + 1) % n]
    return array.array("d", positions), array.array("i", edges)


def lengths(positions, edges):
    out = []
    for k in range(0, len(edges), 2):
        a, b = edges[k], edges[k + 1]
        out.append(math.dist(positions[3 * a : 3 * a + 3], positions[3 * b : 3 * b + 3]))
    return out


def max_z(positions):
    return max(positions[2::3])


PLANE = (
    array.array("d", [-8, -8, 1.2, 8, -8, 1.2, 8, 8, 1.2, -8, 8, 1.2]),
    array.array("i", [0, 1, 2, 0, 2, 3]),
)


class Flow(unittest.TestCase):
    def test_trefoil_lowers_energy_and_keeps_edge_lengths(self):
        p, e = trefoil()
        e0, _ = ddg_tools.curve_energy(p, e)
        q, f, report = ddg_tools.repulsive_flow(p, e, steps=15)
        e1, _ = ddg_tools.curve_energy(q, f)
        self.assertLess(e1, e0)
        self.assertAlmostEqual(report["energy_initial"], e0, delta=1e-9 * abs(e0))
        self.assertAlmostEqual(report["energy_final"], e1, delta=1e-9 * abs(e1))
        self.assertEqual(report["iterations"], 15)
        self.assertEqual(len(report["steps"]), 15)
        before, after = lengths(p, e), lengths(q, f)
        mean = sum(before) / len(before)
        for a, b in zip(before, after):
            self.assertAlmostEqual(a, b, delta=1e-3 * mean)
        # The energy never rises, step by step.
        for s in report["steps"]:
            self.assertLessEqual(s["energy_after"], s["energy_before"])
            self.assertEqual(s["status"], "accepted")

    def test_total_length_and_growth(self):
        p, e = trefoil(48)
        total0 = sum(lengths(p, e))
        q, f, _ = ddg_tools.repulsive_flow(p, e, steps=10, length="total")
        self.assertAlmostEqual(sum(lengths(q, f)), total0, delta=1e-3 * total0)
        # Growth by a factor: 10 steps of rate 0.02 make 20 percent.
        q, f, _ = ddg_tools.repulsive_flow(p, e, steps=10, grow_factor=1.5, grow_rate=0.02)
        self.assertAlmostEqual(sum(lengths(q, f)), 1.2 * total0, delta=2e-3 * total0)
        # Shrinking.
        q, f, _ = ddg_tools.repulsive_flow(p, e, steps=5, grow_factor=0.5, grow_rate=0.02)
        self.assertAlmostEqual(sum(lengths(q, f)), 0.9 * total0, delta=2e-3 * total0)

    def test_pins_stay(self):
        p, e = trefoil()
        pinned = [0, 21]
        q, _, _ = ddg_tools.repulsive_flow(p, e, steps=12, pinned=pinned)
        mean = sum(lengths(p, e)) / (len(e) // 2)
        for v in pinned:
            self.assertAlmostEqual(
                math.dist(p[3 * v : 3 * v + 3], q[3 * v : 3 * v + 3]), 0.0, delta=1e-3 * mean
            )
        free = math.dist(p[3 * 10 : 3 * 10 + 3], q[3 * 10 : 3 * 10 + 3])
        self.assertGreater(free, 1e-3 * mean)

    def test_fixed_tangent(self):
        n = 30
        p = array.array("d")
        for i in range(n):
            t = i / (n - 1)
            p += array.array("d", [3 * t, 0.4 * math.sin(5 * t), 0.4 * math.cos(5 * t)])
        e = array.array("i")
        for i in range(n - 1):
            e += array.array("i", [i, i + 1])

        def direction(q):
            d = [q[3 + k] - q[k] for k in range(3)]
            m = math.hypot(*d)
            return [x / m for x in d]

        q, _, _ = ddg_tools.repulsive_flow(p, e, steps=10, tangent_vertices=[0])
        self.assertGreater(sum(a * b for a, b in zip(direction(p), direction(q))), 1 - 1e-4)
        q2, _, _ = ddg_tools.repulsive_flow(p, e, steps=10)
        self.assertLess(sum(a * b for a, b in zip(direction(p), direction(q2))), 1 - 1e-4)

    def test_surface_constraint(self):
        n = 40
        p = array.array("d")
        for i in range(n):  # a wavy curve on the unit sphere
            t = 2 * math.pi * i / n
            lat = 0.5 * math.sin(2 * t)
            p += array.array(
                "d", [math.cos(lat) * math.cos(t), math.cos(lat) * math.sin(t), math.sin(lat)]
            )
        e = array.array("i", [x for i in range(n) for x in (i, (i + 1) % n)])
        q, _, report = ddg_tools.repulsive_flow(
            p, e, steps=8, surface={"type": "sphere", "radius": 1.0}, barycenter=False
        )
        for v in range(n):
            self.assertAlmostEqual(math.hypot(*q[3 * v : 3 * v + 3]), 1.0, delta=2e-3)
        self.assertGreater(len(report["steps"]), 0)

    def test_penalties_and_fields_run(self):
        p, e = trefoil(40)
        base = ddg_tools.repulsive_flow(p, e, steps=3)[2]
        self.assertEqual(base["steps"][0]["penalties_after"], [])
        for options, n_penalties in [
            ({"total_length_weight": 0.5, "length": "free"}, 1),
            ({"length_difference_weight": 0.5, "length": "free"}, 1),
            ({"field": {"type": "constant", "vector": (0, 0, 1)}, "field_weight": 0.5}, 1),
            ({"field": {"type": "circular"}}, 1),
            ({"obstacles": [PLANE], "obstacle_weight": 0.5}, 1),
        ]:
            report = ddg_tools.repulsive_flow(p, e, steps=3, **options)[2]
            self.assertEqual(len(report["steps"][0]["penalties_after"]), n_penalties, options)

    def test_planar(self):
        n = 24
        p = array.array("d")
        for i in range(n):
            t = 2 * math.pi * i / n
            p += array.array("d", [math.cos(t) * (1 + 0.3 * math.cos(3 * t)), math.sin(t), 0.5])
        e = array.array("i", [x for i in range(n) for x in (i, (i + 1) % n)])
        q, _, _ = ddg_tools.repulsive_flow(p, e, steps=5, planar_axis="z")
        for v in range(n):
            self.assertAlmostEqual(q[3 * v + 2], 0.5, delta=1e-12)

    def test_backends(self):
        p, e = trefoil(64)
        reference = ddg_tools.repulsive_flow(p, e, steps=2, backend="dense")[0]
        for backend in ["auto", "barnes_hut", "multigrid", "accelerated"]:
            q, _, report = ddg_tools.repulsive_flow(p, e, steps=2, backend=backend)
            d = max(abs(a - b) for a, b in zip(reference, q))
            self.assertLess(d, 5e-2, backend)
            self.assertLess(report["energy_final"], report["energy_initial"], backend)
        self.assertEqual(
            ddg_tools.repulsive_flow(p, e, steps=1, backend="accelerated")[2]["steps"][0][
                "solver_backend"
            ],
            "multigrid",
        )

    def test_obstacle_is_never_crossed(self):
        p, e = trefoil(48)
        options = dict(grow_factor=1.5, grow_rate=0.05)
        # Without the obstacle the growing knot passes the plane z = 1.2 ...
        free = ddg_tools.repulsive_flow(p, e, steps=12, **options)[0]
        self.assertGreater(max_z(free), 1.2)
        # ... with it, no step does.
        session = ddg_tools.RepulsiveSession(p, e, obstacles=[PLANE], **options)
        for _ in range(12):
            session.step()
            self.assertLess(max_z(session.positions), 1.2)
        self.assertGreater(max_z(session.positions), max_z(p))
        # Also in the one-shot flow, with the potential on.
        q = ddg_tools.repulsive_flow(
            p, e, steps=12, obstacles=[PLANE], obstacle_weight=0.1, **options
        )[0]
        self.assertLess(max_z(q), 1.2)


class Session(unittest.TestCase):
    def test_stepping_matches_flow_bitwise(self):
        p, e = trefoil(56)
        options = dict(pinned=[3], tangent_vertices=[], threads=3)
        q, f, _ = ddg_tools.repulsive_flow(p, e, steps=8, **options)
        s = ddg_tools.RepulsiveSession(p, e, max_steps=8, **options)
        for k in range(8):
            r = s.step()
            self.assertEqual(r["iteration"], k + 1)
        self.assertEqual(s.positions.tolist(), q.tolist())
        self.assertEqual(s.edges.tolist(), f.tolist())
        self.assertEqual(s.status, "max_steps")
        r = ddg_tools.RepulsiveSession(p, e, max_steps=8, **options)
        summary = r.run()
        self.assertEqual(summary["steps_taken"], 8)
        self.assertEqual(r.positions.tolist(), q.tolist())
        # Thread count does not change the curve.
        other = ddg_tools.repulsive_flow(p, e, steps=8, pinned=[3], threads=1)[0]
        self.assertEqual(other.tolist(), q.tolist())

    def test_state_and_report(self):
        p, e = trefoil(40)
        s = ddg_tools.RepulsiveSession(p, e)
        self.assertEqual((s.n_vertices, s.n_edges, s.iteration, s.status), (40, 40, 0, "running"))
        self.assertIsNone(s.report)
        e0 = s.energy
        r = s.step()
        for key in ["iteration", "status", "session_status", "energy_before", "energy_after",
                    "gradient_norm", "tau", "halvings", "max_violation", "linear", "retarget",
                    "subdivided", "seconds", "n_vertices", "penalties_after"]:
            self.assertIn(key, r)
        self.assertAlmostEqual(r["energy_before"], e0, delta=1e-9 * e0)
        self.assertEqual(s.energy, r["energy_after"])
        self.assertEqual(s.report["iteration"], 1)
        pos, edges, it, energy, status = s.snapshot()
        self.assertEqual((it, energy, status), (1, s.energy, "running"))
        self.assertEqual(len(pos), 120)
        self.assertEqual(len(edges), 80)
        self.assertIn("40 vertices", repr(s).replace("vertices=40", "40 vertices"))

    def test_run_with_progress(self):
        p, e = trefoil(40)
        s = ddg_tools.RepulsiveSession(p, e)
        seen = []

        def progress(report):
            seen.append(report["iteration"])
            return len(seen) < 3

        summary = s.run(progress=progress)
        self.assertEqual(summary["end"], "callback")
        self.assertEqual(summary["steps_taken"], 3)
        self.assertEqual(seen, [1, 2, 3])
        self.assertEqual(s.status, "running")
        summary = s.run(max_steps=2)
        self.assertEqual((summary["end"], summary["steps_taken"], s.iteration), ("max_steps", 2, 5))
        self.assertLess(summary["energy_after"], summary["energy_before"])
        # An exception in the callback propagates and stops the run.
        def boom(report):
            raise RuntimeError("boom")
        with self.assertRaises(RuntimeError):
            s.run(progress=boom)
        self.assertEqual(s.iteration, 6)
        # Same for the one-shot flow.
        with self.assertRaises(RuntimeError):
            ddg_tools.repulsive_flow(p, e, steps=5, progress=boom)
        calls = []
        ddg_tools.repulsive_flow(p, e, steps=5, progress=lambda r: calls.append(r["iteration"]))
        self.assertEqual(calls, [1, 2, 3, 4, 5])
        stopped = ddg_tools.repulsive_flow(p, e, steps=5, progress=lambda r: False)[2]
        self.assertEqual(stopped["status"], "cancelled")

    def test_converges(self):
        n = 24  # a circle with equal edges is at rest
        p = array.array("d")
        for i in range(n):
            t = 2 * math.pi * i / n
            p += array.array("d", [math.cos(t), math.sin(t), 0.0])
        e = array.array("i", [x for i in range(n) for x in (i, (i + 1) % n)])
        s = ddg_tools.RepulsiveSession(p, e, stop_threshold=1e-3)
        summary = s.run()
        self.assertEqual(summary["end"], "final")
        self.assertEqual(s.status, "converged")
        self.assertEqual(summary["steps"][-1]["status"], "converged")
        # A converged session does not step on its own.
        self.assertEqual(s.run()["steps_taken"], 0)

    def test_cancel_and_reset(self):
        p, e = trefoil(40)
        s = ddg_tools.RepulsiveSession(p, e)
        s.step()
        moved = s.positions.tolist()
        s.cancel()
        self.assertEqual(s.status, "cancelled")
        r = s.step()
        self.assertEqual(r["status"], "cancelled")
        self.assertEqual(s.positions.tolist(), moved)
        self.assertEqual(s.run()["end"], "cancelled")
        s.reset()
        self.assertEqual((s.status, s.iteration), ("running", 0))
        self.assertEqual(s.positions.tolist(), p.tolist())
        s.step()
        self.assertEqual(s.iteration, 1)

    def test_cancel_from_another_thread_while_stepping(self):
        p, e = trefoil(300)
        s = ddg_tools.RepulsiveSession(p, e, max_steps=100000)
        result = {}
        worker = threading.Thread(target=lambda: result.update(s.run()))
        worker.start()
        # The GIL is released: this thread keeps running and sees the progress.
        deadline = time.time() + 120
        ticks = 0
        while s.iteration < 2 and time.time() < deadline:
            ticks += 1
            time.sleep(0.01)
        self.assertGreaterEqual(s.iteration, 2)
        self.assertGreater(ticks, 0)
        s.cancel()
        worker.join(60)
        self.assertFalse(worker.is_alive())
        self.assertEqual(result["end"], "cancelled")
        self.assertEqual(s.status, "cancelled")

    def test_set_schedule(self):
        p, e = trefoil(40)
        total0 = sum(lengths(p, e))
        s = ddg_tools.RepulsiveSession(p, e)
        s.run(max_steps=2)
        s.set_schedule(1.2, 0.04)
        s.run(max_steps=5)
        self.assertAlmostEqual(s.total_length, 1.2 * total0, delta=2e-3 * total0)
        s.set_schedule(1.0)
        s.run(max_steps=2)
        self.assertAlmostEqual(s.total_length, 1.2 * total0, delta=2e-3 * total0)
        with self.assertRaises(ValueError):
            s.set_schedule(-1.0)
        free = ddg_tools.RepulsiveSession(p, e, length="free")
        with self.assertRaises(ValueError):
            free.set_schedule(1.2)

    def test_subdivision(self):
        n = 40
        p, e = trefoil(n)
        s = ddg_tools.RepulsiveSession(p, e, max_subdivisions=1, grow_factor=2.5, grow_rate=0.1)
        self.assertEqual(s.subdivisions, 0)
        self.assertEqual(s.vertex_source.tolist(), list(range(n)))
        steps = 0
        while s.n_vertices == n and steps < 60:
            r = s.step()
            steps += 1
            self.assertEqual(r["subdivided"], s.n_vertices != n)
        self.assertEqual(s.n_vertices, 2 * n)
        self.assertEqual(s.n_edges, 2 * n)
        self.assertEqual(s.subdivisions, 1)
        self.assertEqual(len(s.positions), 6 * n)
        self.assertEqual(len(s.edges), 4 * n)
        self.assertEqual(s.vertex_source.tolist(), list(range(n)) + [-1] * n)
        self.assertEqual(s.edge_source.tolist(), list(range(n)) + list(range(n)))
        edges = s.edges.tolist()
        for i in range(n):
            self.assertEqual(edges[2 * i : 2 * i + 2], [i, n + i])
            self.assertEqual(edges[2 * (n + i) : 2 * (n + i) + 2], [n + i, (i + 1) % n])
        # The midpoints are where the edge's halves meet: between their ends.
        pos = s.positions
        for i in range(0, n, 7):
            a, m, b = (pos[3 * v : 3 * v + 3] for v in (i, n + i, (i + 1) % n))
            self.assertLess(abs(math.dist(a, m) - math.dist(m, b)), 0.5 * math.dist(a, b))
        # The session stays usable, and no second subdivision happens (the limit is 1).
        s.run(max_steps=5)
        self.assertEqual(s.n_vertices, 2 * n)
        # The one-shot flow reports the same.
        q, f, report = ddg_tools.repulsive_flow(
            p, e, steps=steps + 3, max_subdivisions=1, grow_factor=2.5, grow_rate=0.1
        )
        self.assertEqual(len(q), 6 * n)
        self.assertEqual(len(f), 4 * n)
        self.assertEqual(report["n_subdivisions"], 1)
        self.assertEqual(report["vertex_source"], list(range(n)) + [-1] * n)
        # reset goes back to the input topology.
        s.reset()
        self.assertEqual((s.n_vertices, s.n_edges, s.subdivisions), (n, n, 0))

    def test_junction_and_open_curves(self):
        # A "Y": three arms meeting in vertex 0.
        pos, edges = [0.0, 0.0, 0.0], []
        for arm, d in enumerate([(1, 0, 0), (-0.5, 0.87, 0), (-0.5, -0.87, 0.2)]):
            prev = 0
            for k in range(1, 8):
                pos += [d[0] * k * 0.3, d[1] * k * 0.3, d[2] * k * 0.3 + 0.01 * k * k * (arm - 1)]
                edges += [prev, len(pos) // 3 - 1]
                prev = len(pos) // 3 - 1
        s = ddg_tools.RepulsiveSession(pos, edges, max_steps=3)
        self.assertEqual(s.n_vertices, 22)
        r = s.step()
        self.assertEqual(r["status"], "accepted")


class Contract(unittest.TestCase):
    @staticmethod
    def stub_defaults(name, kind):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ddg_tools.pyi")
        with open(path, encoding="utf-8") as f:
            tree = ast.parse(f.read())
        for node in tree.body:
            if kind == "class" and isinstance(node, ast.ClassDef) and node.name == name:
                for item in node.body:
                    if isinstance(item, ast.FunctionDef) and item.name == "__init__":
                        fn = item
                        break
                else:
                    raise AssertionError("no __init__")
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

    def test_stub_defaults_equal_the_library(self):
        defaults = ddg_tools.repulsive_defaults()
        self.assertEqual(set(defaults), set(self.stub_defaults("RepulsiveSession", "class")))
        self.assertEqual(self.stub_defaults("RepulsiveSession", "class"), defaults)
        flow = self.stub_defaults("repulsive_flow", "func")
        self.assertEqual(flow, defaults)
        # The numbers are the paper's and the library's (SPEC section 10).
        self.assertEqual((defaults["alpha"], defaults["beta"]), (3.0, 6.0))
        self.assertEqual(defaults["line_search"], "collision")
        self.assertEqual(defaults["length"], "edges")
        self.assertEqual(defaults["armijo"], 1e-4)
        self.assertEqual(defaults["grow_rate"], 0.01)
        self.assertEqual(defaults["max_steps"], 1000)
        self.assertEqual(defaults["max_subdivisions"], 0)
        self.assertEqual(defaults["projection_iterations"], 3)

    def test_energy_function_defaults(self):
        d = ddg_tools.repulsive_defaults()
        energy = self.stub_defaults("curve_energy", "func")
        for key in ["bh_theta_x", "bh_theta_t", "bh_leaf_size", "threads"]:
            self.assertEqual(energy[key], d[key], key)

    def test_explicit_defaults_change_nothing(self):
        p, e = trefoil(40)
        explicit = {k: v for k, v in ddg_tools.repulsive_defaults().items() if v is not None}
        a = ddg_tools.repulsive_flow(p, e, steps=3)[0]
        b = ddg_tools.repulsive_flow(p, e, steps=3, **explicit)[0]
        self.assertEqual(a.tolist(), b.tolist())
        s = ddg_tools.RepulsiveSession(p, e, **explicit)
        self.assertEqual(s.status, "running")

    def test_energy_matches_finite_differences(self):
        p, e = trefoil(32)
        energy, diff = ddg_tools.curve_energy(p, e)
        self.assertEqual(len(diff), len(p))
        h = 1e-6
        for k in [0, 5, 31, 50, 95]:
            up, dn = array.array("d", p), array.array("d", p)
            up[k] += h
            dn[k] -= h
            fd = (ddg_tools.curve_energy(up, e)[0] - ddg_tools.curve_energy(dn, e)[0]) / (2 * h)
            self.assertAlmostEqual(fd, diff[k], delta=1e-4 * max(1.0, abs(fd)))
        bh, bh_diff = ddg_tools.curve_energy(p, e, backend="barnes_hut")
        self.assertAlmostEqual(bh, energy, delta=0.05 * energy)
        # Scaling: E(c x) = c^(2 + alpha - beta) E(x).
        scaled = array.array("d", [2 * x for x in p])
        self.assertAlmostEqual(
            ddg_tools.curve_energy(scaled, e)[0], energy * 2.0 ** (2 + 3 - 6), delta=1e-9 * energy
        )

    def test_python_sequences_and_numpy_are_accepted(self):
        p, e = trefoil(24)
        a = ddg_tools.repulsive_flow(p, e, steps=2)[0]
        nested = [list(p[3 * i : 3 * i + 3]) for i in range(24)]
        pairs = [list(e[2 * i : 2 * i + 2]) for i in range(24)]
        self.assertEqual(ddg_tools.repulsive_flow(nested, pairs, steps=2)[0].tolist(), a.tolist())
        self.assertEqual(
            ddg_tools.repulsive_flow(list(p), list(e), steps=2)[0].tolist(), a.tolist()
        )
        if numpy is not None:
            q = numpy.array(nested).reshape(24, 3)
            f = numpy.array(pairs, dtype=numpy.int32)
            b = ddg_tools.repulsive_flow(q, f, steps=2)[0]
            self.assertEqual(b.tolist(), a.tolist())
            s = ddg_tools.RepulsiveSession(q, f)
            self.assertEqual(numpy.frombuffer(s.positions).reshape(-1, 3).shape, (24, 3))


class Multigrid(unittest.TestCase):
    """The multigrid solver choice, its options, the start rule and the fallback field."""

    SMALL = dict(mg_coarsest_vertices=16)

    def test_multigrid_options_and_their_defaults(self):
        d = ddg_tools.repulsive_defaults()
        self.assertEqual(d["backend"], "auto")
        self.assertEqual(d["auto_dense_max_vertices"], 2500)
        self.assertEqual(d["dense_fallback_max_vertices"], 5000)
        self.assertEqual(d["mg_coarsest_vertices"], 512)
        self.assertEqual(d["mg_near_field_factor"], 1.0)
        self.assertEqual(d["mg_smoother_factor"], 1.0)
        self.assertEqual(d["mg_max_fixed_fraction"], 0.4)
        self.assertIs(d["mg_energy_check"], True)
        self.assertEqual(d["mg_failure_residual"], 1e-2)
        self.assertEqual(d["line_search_start"], "grow_previous")
        self.assertEqual(d["line_search_start_factor"], 2.0)

    def test_multigrid_backend_runs_with_a_hierarchy(self):
        p, e = trefoil(64)
        for backend in ["multigrid", "accelerated"]:
            _, _, report = ddg_tools.repulsive_flow(
                p, e, steps=3, backend=backend, **self.SMALL
            )
            self.assertLess(report["energy_final"], report["energy_initial"])
            for step in report["steps"]:
                self.assertEqual(step["solver_backend"], "multigrid")
                self.assertEqual(step["linear"]["solver"], "multigrid")
                self.assertGreaterEqual(step["linear"]["levels"], 2)
                self.assertIsNone(step["solver_fallback"])
        self.assertEqual(report["steps"][0]["energy_backend"], "barnes-hut")

    def test_every_multigrid_option_is_accepted_and_validated(self):
        p, e = trefoil(64)
        options = dict(
            mg_coarsest_vertices=16,
            mg_max_levels=4,
            mg_smoothing_steps=6,
            mg_initial_smoothing_steps=12,
            mg_tolerance=1e-4,
            mg_max_cycles=30,
            mg_flexible_cg=False,
            mg_near_field_factor=0.5,
            mg_smoother_factor=0.0,
            mg_max_fixed_fraction=0.9,
            mg_energy_check=False,
            mg_failure_residual=0.1,
        )
        _, _, report = ddg_tools.repulsive_flow(p, e, steps=2, backend="multigrid", **options)
        self.assertEqual(report["steps"][0]["linear"]["solver"], "multigrid")
        for bad in [
            dict(mg_tolerance=0.0),
            dict(mg_failure_residual=1e-9),
            dict(mg_max_fixed_fraction=0.0),
            dict(mg_near_field_factor=-1.0),
            dict(mg_max_levels=0),
        ]:
            with self.assertRaises(ValueError, msg=str(bad)):
                ddg_tools.repulsive_flow(p, e, steps=1, backend="multigrid", **bad)

    def test_auto_follows_the_crossover_option(self):
        p, e = trefoil(64)
        small = ddg_tools.repulsive_flow(p, e, steps=1, auto_dense_max_vertices=10)[2]["steps"][0]
        self.assertEqual((small["energy_backend"], small["solver_backend"]), ("barnes-hut", "multigrid"))
        big = ddg_tools.repulsive_flow(p, e, steps=1)[2]["steps"][0]
        self.assertEqual((big["energy_backend"], big["solver_backend"]), ("exact", "direct"))

    def test_start_rules(self):
        p, e = trefoil(64)
        results = {}
        for rule in ["two_thirds", "grow_previous"]:
            q, _, report = ddg_tools.repulsive_flow(p, e, steps=6, line_search_start=rule)
            self.assertLess(report["energy_final"], report["energy_initial"], rule)
            results[rule] = q.tolist()
        self.assertNotEqual(results["two_thirds"], results["grow_previous"])
        # The default is the grow-previous rule with factor 2.
        a = ddg_tools.repulsive_flow(p, e, steps=6)[0].tolist()
        b = ddg_tools.repulsive_flow(
            p, e, steps=6, line_search_start="grow_previous", line_search_start_factor=2.0
        )[0].tolist()
        self.assertEqual(a, b)
        self.assertEqual(a, results["grow_previous"])
        with self.assertRaises(ValueError):
            ddg_tools.repulsive_flow(p, e, steps=1, line_search_start="sideways")
        with self.assertRaises(ValueError):
            ddg_tools.repulsive_flow(p, e, steps=1, line_search_start_factor=0.5)

    FAILING = dict(
        backend="multigrid",
        mg_coarsest_vertices=16,
        mg_tolerance=1e-14,
        mg_failure_residual=1e-13,
        mg_max_cycles=1,
    )

    def test_a_failed_multigrid_solve_falls_back_and_says_so(self):
        p, e = trefoil(64)
        dense = ddg_tools.RepulsiveSession(p, e, backend="dense")
        mg = ddg_tools.RepulsiveSession(p, e, **self.FAILING)
        for _ in range(2):
            rd = dense.step()
            rm = mg.step()
            fb = rm["solver_fallback"]
            self.assertIsNotNone(fb)
            self.assertEqual((fb["from"], fb["to"]), ("multigrid", "direct"))
            self.assertEqual(fb["cause"], "multigrid_not_converged")
            self.assertIn("did not converge", fb["message"])
            self.assertGreaterEqual(fb["dense_solves"], 1)
            self.assertIsNone(rd["solver_fallback"])
        # The fallback step is the dense step.
        self.assertEqual(mg.positions.tolist(), dense.positions.tolist())

    def test_the_fallback_can_be_refused(self):
        p, e = trefoil(64)
        with self.assertRaises(ddg_tools.SolveError) as cm:
            ddg_tools.RepulsiveSession(p, e, dense_fallback_max_vertices=0, **self.FAILING).step()
        self.assertIn("did not converge", str(cm.exception))

    def test_schedule_progress_is_in_the_step_report(self):
        p, e = trefoil(40)
        s = ddg_tools.RepulsiveSession(p, e, grow_factor=1.3, grow_rate=0.1)
        self.assertIsNone(ddg_tools.RepulsiveSession(p, e).step()["schedule"])
        last = None
        for _ in range(5):
            last = s.step()["schedule"]
        self.assertEqual(last["factor_target"], 1.3)
        self.assertAlmostEqual(last["factor_imposed"], 1.3, places=12)
        self.assertAlmostEqual(last["factor_measured"], 1.3, delta=1e-3)
        self.assertTrue(last["finished"])
        self.assertAlmostEqual(last["fraction"], 1.0)

    def test_min_edge_length_ratio_option(self):
        p, e = trefoil(40)
        self.assertEqual(ddg_tools.repulsive_defaults()["min_edge_length_ratio"], 1e-3)
        ddg_tools.repulsive_flow(p, e, steps=2, min_edge_length_ratio=None)
        for bad in (0.0, 1.0, -1.0):
            with self.assertRaises(ValueError):
                ddg_tools.repulsive_flow(p, e, steps=1, min_edge_length_ratio=bad)

    def test_dense_threads_is_an_option(self):
        p, e = trefoil(40)
        a = ddg_tools.repulsive_flow(p, e, steps=2, threads=2, dense_threads=1)[0].tolist()
        b = ddg_tools.repulsive_flow(p, e, steps=2, threads=3, dense_threads=2)[0].tolist()
        self.assertEqual(a, b)  # bitwise independent of the thread counts
        with self.assertRaises(ValueError):
            ddg_tools.repulsive_flow(p, e, steps=1, dense_threads=0)


class Errors(unittest.TestCase):
    def setUp(self):
        self.p, self.e = trefoil(30)

    def assertInput(self, *args, contains=None, **kwargs):
        for call in (ddg_tools.repulsive_flow, ddg_tools.RepulsiveSession):
            with self.assertRaises(ValueError) as cm:
                call(*args, **kwargs)
            self.assertIsInstance(cm.exception, ddg_tools.InputError)
            self.assertIsInstance(cm.exception, ddg_tools.RepulsiveInputError)
            self.assertIsInstance(cm.exception, ddg_tools.DdgError)
            if contains:
                self.assertIn(contains, str(cm.exception))

    def test_bad_curve(self):
        self.assertInput(self.p, [0, 99], contains="99")
        self.assertInput(self.p, [3, 3])
        self.assertInput(self.p, [0, 1, 1, 0])
        self.assertInput(self.p[:-1], self.e, contains="positions")
        self.assertInput(self.p, self.e[:-1], contains="edges")
        self.assertInput([0.0] * 6 + [float("nan")] + [0.0] * 2, [0, 1, 1, 2, 2, 0])
        coincident = array.array("d", self.p)
        coincident[3:6] = coincident[0:3]
        self.assertInput(coincident, self.e)
        self.assertInput(self.p, ["a", "b"])
        self.assertInput(self.p, [0, -1])

    def test_bad_options(self):
        a, e = self.p, self.e
        self.assertInput(a, e, alpha=0.5)
        self.assertInput(a, e, beta=float("nan"))
        self.assertInput(a, e, alpha="three", contains="alpha")
        self.assertInput(a, e, nonsense=1, contains="nonsense")
        self.assertInput(a, e, length="some")
        self.assertInput(a, e, line_search="magic", contains="line_search")
        self.assertInput(a, e, armijo=2.0)
        self.assertInput(a, e, shrink=0.0)
        self.assertInput(a, e, backend="gpu", contains="backend")
        self.assertInput(a, e, threads=0)
        self.assertInput(a, e, threads=-1)
        self.assertInput(a, e, max_steps=-1)
        self.assertInput(a, e, stop="never")
        self.assertInput(a, e, stop="stall", contains="stop_threshold")
        self.assertInput(a, e, stop_threshold=-1.0)
        self.assertInput(a, e, time_limit=-1.0)
        self.assertInput(a, e, planar_axis="w")
        self.assertInput(a, e, planar_axis="z", contains="plan")  # a trefoil is not planar
        self.assertInput(a, e, pinned=[99])
        self.assertInput(a, e, pinned=[-1])
        self.assertInput(a, e, tangent_vertices=[99])
        self.assertInput(a, e, grow_factor=-1.0)
        self.assertInput(a, e, grow_factor=2.0, length="free")
        self.assertInput(a, e, grow_rate=0.0, grow_factor=2.0)
        self.assertInput(a, e, total_length_weight=float("nan"))
        self.assertInput(a, e, obstacle_weight="heavy")
        self.assertInput(a, e, max_subdivisions=-1)
        self.assertInput(a, e, metric="l7")

    def test_bad_obstacles_surfaces_fields(self):
        a, e = self.p, self.e
        tri = array.array("d", [0, 0, 5, 1, 0, 5, 0, 1, 5])
        self.assertInput(a, e, obstacles=[(tri, array.array("i", [0, 1, 7]))])
        self.assertInput(a, e, obstacles=[(tri, array.array("i", []))])
        self.assertInput(a, e, obstacles=[tri])
        self.assertInput(a, e, obstacles=[(array.array("d", [0, 0, 5, 1, 0, 5, 2, 0, 5]), [0, 1, 2])])
        self.assertInput(a, e, surface={"type": "cube"}, contains="cube")
        self.assertInput(a, e, surface={"radius": 1.0})
        self.assertInput(a, e, surface={"type": "sphere", "radius": -1.0})
        self.assertInput(a, e, surface={"type": "sphere", "colour": "red"}, contains="colour")
        self.assertInput(a, e, surface={"type": "torus", "major": 2.0})
        self.assertInput(a, e, surface="sphere")
        self.assertInput(a, e, field={"type": "swirl"})
        self.assertInput(a, e, field={"type": "constant"})
        self.assertInput(a, e, field={"type": "constant", "vector": (0, 0, 0), "normalize": True})
        self.assertInput(a, e, field={"type": "constant", "vector": (1, 2)})
        self.assertInput(a, e, field_weight=-1.0, field={"type": "circular"})

    def test_infeasible_constraints_are_input_errors(self):
        # A surface the curve is not on: the library names the violation.
        self.assertInput(
            self.p, self.e, surface={"type": "sphere", "radius": 0.5}, barycenter=False
        )

    def test_session_methods_validate(self):
        s = ddg_tools.RepulsiveSession(self.p, self.e)
        with self.assertRaises(ValueError):
            s.run(max_steps=-1)
        with self.assertRaises(TypeError):
            ddg_tools.RepulsiveSession(self.p)
        with self.assertRaises(ValueError):
            ddg_tools.curve_energy(self.p, self.e, alpha=0.5)
        with self.assertRaises(ValueError):
            ddg_tools.curve_energy(self.p, self.e, backend="fast")
        with self.assertRaises(ValueError):
            ddg_tools.curve_energy(self.p, [0, 1])


if __name__ == "__main__":
    unittest.main(verbosity=2)
