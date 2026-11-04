"""Tests of the Points as Tori bindings of the built wheel of `ddg_tools`, without Blender.

Run through `scripts/py-test.sh` (see test_ddg_tools.py). Standard library only; the numpy tests
are skipped when numpy is not installed, the parity tests when `runs_tmp/pat-ref/` (the golden
outputs of the authors' program) is absent, the GPU test unless the wheel has CUDA, a device is
usable and no process is using the GPU.
"""

import array
import ast
import math
import os
import struct
import subprocess
import tempfile
import unittest

import ddg_tools
from ddg_tools import points_as_tori as pat

try:
    import numpy
except ImportError:  # pragma: no cover
    numpy = None

HERE = os.path.dirname(os.path.abspath(__file__))
GOLDEN = os.path.join(HERE, "..", "..", "..", "runs_tmp", "pat-ref")


# ---------------------------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------------------------


def fibonacci_sphere(n=512, radius=0.5, center=(0.0, 0.0, 0.0)):
    """n Fibonacci points on a sphere and their outward unit normals (flat arrays)."""
    golden = math.pi * (3.0 - math.sqrt(5.0))
    pts, nrm = array.array("d"), array.array("d")
    for i in range(n):
        z = 1.0 - 2.0 * (i + 0.5) / n
        r = math.sqrt(1.0 - z * z)
        a = golden * i
        d = (r * math.cos(a), r * math.sin(a), z)
        nrm.extend(d)
        pts.extend(radius * d[k] + center[k] for k in range(3))
    return pts, nrm


def circle(n=64, radius=1.0):
    pts, nrm = array.array("d"), array.array("d")
    for i in range(n):
        a = 2 * math.pi * i / n
        nrm.extend((math.cos(a), math.sin(a)))
        pts.extend((radius * math.cos(a), radius * math.sin(a)))
    return pts, nrm


def octahedron():
    v = array.array("d", [1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1])
    f = array.array("i", [0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5])
    return v, f


def read_npy(path):
    """float64 / int32 / float32 .npy arrays as (flat array.array, shape)."""
    with open(path, "rb") as fh:
        assert fh.read(6) == b"\x93NUMPY"
        major, _ = struct.unpack("BB", fh.read(2))
        hlen = struct.unpack("<H" if major == 1 else "<I", fh.read(2 if major == 1 else 4))[0]
        header = ast.literal_eval(fh.read(hlen).decode("latin1"))
        code = {"<f8": "d", "<f4": "f", "<i4": "i", "<i8": "q"}[header["descr"]]
        assert not header["fortran_order"]
        a = array.array(code)
        a.frombytes(fh.read())
        return a, header["shape"]


def max_abs_diff(a, b):
    assert len(a) == len(b)
    return max(abs(x - y) for x, y in zip(a, b))


def rows(a, w=3):
    return [tuple(a[i * w:(i + 1) * w]) for i in range(len(a) // w)]


def grid_queries(lo=-0.9, hi=0.9, n=4):
    out = array.array("d")
    for i in range(n):
        for j in range(n):
            for k in range(n):
                out.extend(lo + (hi - lo) * t / (n - 1) for t in (i, j, k))
    return out


SPHERE = fibonacci_sphere()
FIELD = ddg_tools.PointsAsTori(*SPHERE)
Q = grid_queries()


def exact_sphere(q, r=0.5):
    return [math.sqrt(x * x + y * y + z * z) - r for x, y, z in rows(q)]


# ---------------------------------------------------------------------------------------------


class Basics(unittest.TestCase):
    def test_submodule_is_the_same_objects(self):
        import importlib

        sub = importlib.import_module("ddg_tools.points_as_tori")
        self.assertIs(sub, pat)
        self.assertIs(pat.PointsAsTori, ddg_tools.PointsAsTori)
        self.assertIs(pat.reconstruct, ddg_tools.pat_reconstruct)
        self.assertIs(pat.CsgField, ddg_tools.PatCsgField)
        self.assertIs(pat.gpu_info, ddg_tools.gpu_info)

    def test_sphere_distance_gradient_laplacian(self):
        d = FIELD.signed_distance(Q)
        self.assertIsInstance(d, array.array)
        self.assertEqual(len(d), len(Q) // 3)
        self.assertLess(max_abs_diff(d, exact_sphere(Q)), 0.03)
        d2, g = FIELD.distance_and_gradient(Q)
        self.assertEqual(d.tobytes(), d2.tobytes())
        self.assertEqual(FIELD.gradient(Q).tobytes(), g.tobytes())
        self.assertEqual(FIELD.distance(Q).tobytes(), d.tobytes())
        for (x, y, z), (gx, gy, gz) in zip(rows(Q), rows(g)):
            r = math.sqrt(x * x + y * y + z * z)
            if r > 0.2:
                self.assertAlmostEqual(gx * x / r + gy * y / r + gz * z / r, 1.0, delta=0.1)
        lap = FIELD.laplacian(Q)
        self.assertEqual(len(lap), len(d))
        # Delta r = 2/r on a sphere: positive away from the origin
        for (x, y, z), l in zip(rows(Q), lap):
            r = math.sqrt(x * x + y * y + z * z)
            if r > 0.7:
                self.assertAlmostEqual(l, 2.0 / r, delta=1.0)

    def test_properties_and_reports(self):
        self.assertEqual(FIELD.n, 512)
        self.assertEqual(len(FIELD.points), 1536)
        self.assertEqual(len(FIELD.normals), 1536)
        br = FIELD.build_report
        self.assertEqual(br["n"], 512)
        self.assertEqual(br["predict_device"], "cpu")
        for key in ("knn_seconds", "forward_seconds", "saddle", "umbilic", "seconds"):
            self.assertIn(key, br)
        FIELD.signed_distance(Q)
        er = FIELD.eval_report
        self.assertEqual(er["queries"], len(Q) // 3)
        self.assertEqual(er["device"], "cpu")
        self.assertEqual(FIELD.cloud_report["n"], 512)
        self.assertEqual(FIELD.options["neighborhood"], "knn")
        self.assertIn("PointsAsTori", repr(FIELD))

    def test_nested_and_sequence_inputs(self):
        rows_ = [list(r) for r in rows(Q)[:5]]
        a = FIELD.signed_distance(rows_)
        b = FIELD.signed_distance(Q[:15])
        self.assertEqual(a.tobytes(), b.tobytes())
        f32 = array.array("f", Q[:15])
        c = FIELD.signed_distance(f32)
        self.assertLess(max_abs_diff(a, c), 1e-5)

    @unittest.skipIf(numpy is None, "numpy not installed")
    def test_numpy_in_and_out(self):
        p = numpy.asarray(SPHERE[0]).reshape(-1, 3)
        n = numpy.asarray(SPHERE[1]).reshape(-1, 3)
        f = ddg_tools.PointsAsTori(p, n)
        q = numpy.asarray(Q).reshape(-1, 3)
        d = numpy.asarray(f.signed_distance(q))
        self.assertEqual(d.shape, (len(q),))
        self.assertEqual(d.tobytes(), numpy.asarray(FIELD.signed_distance(Q)).tobytes())
        g = numpy.asarray(f.gradient(q)).reshape(-1, 3)
        self.assertEqual(g.shape, (len(q), 3))
        self.assertEqual(numpy.asarray(f.points).reshape(-1, 3).tolist(), p.tolist())

    def test_evaluation_options(self):
        self.assertEqual(FIELD.with_options(neighborhood="all").options["neighborhood"], "all")
        a = FIELD.with_options(neighborhood="all").signed_distance(Q)
        self.assertLess(max_abs_diff(a, FIELD.signed_distance(Q)), 0.05)
        c = FIELD.with_options(isovalue=0.1).signed_distance(Q)
        self.assertLess(max_abs_diff(c, [x - 0.1 for x in FIELD.signed_distance(Q)]), 1e-12)
        e = FIELD.with_options(gradient="exact").distance_and_gradient(Q)
        self.assertEqual(len(e[1]), len(Q))
        p = ddg_tools.PointsAsTori(*SPHERE, paper=True)
        o = p.options
        self.assertEqual((o["neighborhood"], o["lambda_rule"]), ("radius", "paper"))
        self.assertLess(max_abs_diff(p.signed_distance(Q), exact_sphere(Q)), 0.05)
        w = ddg_tools.PointsAsTori(*SPHERE, area_weights="barill")
        self.assertTrue(w.options["area_weights"])
        self.assertLess(max_abs_diff(w.signed_distance(Q), exact_sphere(Q)), 0.05)
        g = ddg_tools.PointsAsTori(*SPHERE, sign_rule="geometric", eigen_rule="robust", feature_order="reference", k=32)
        self.assertLess(max_abs_diff(g.signed_distance(Q), exact_sphere(Q)), 0.05)

    def test_normalize_policy(self):
        p, n = SPHERE
        n2 = array.array("d", (2.0 * x for x in n))
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.PointsAsTori(p, n2)
        f = ddg_tools.PointsAsTori(p, n2, normal_policy="normalize")
        self.assertEqual(f.cloud_report["normalized"], 512)
        self.assertLess(max_abs_diff(f.signed_distance(Q), FIELD.signed_distance(Q)), 1e-12)


class RoundTrips(unittest.TestCase):
    def test_save_and_load_tori(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "sphere.tori")
            FIELD.save_tori(path)
            g = ddg_tools.pat_load_tori(path, *SPHERE)
            self.assertIsNone(g.build_report)
            self.assertEqual(FIELD.sign_rule, "geometric")  # the default since 2026-11-04
            self.assertEqual(g.sign_rule, "geometric")  # the file records it
            self.assertEqual(g.signed_distance(Q).tobytes(), FIELD.signed_distance(Q).tobytes())
            self.assertEqual(g.gradient(Q).tobytes(), FIELD.gradient(Q).tobytes())
            self.assertEqual(g.tori["coefficients"].tobytes(), FIELD.tori["coefficients"].tobytes())
            # another cloud: refused
            other = fibonacci_sphere(512, 0.6)
            with self.assertRaises(pat.PointsAsToriError) as cm:
                ddg_tools.pat_load_tori(path, *other)
            self.assertIn("belongs to cloud", str(cm.exception))
            with open(path, "rb") as fh:
                data = fh.read()
            bad = os.path.join(d, "bad.tori")
            with open(bad, "wb") as fh:
                fh.write(data[:100])
            with self.assertRaises(pat.PointsAsToriError):
                ddg_tools.pat_load_tori(bad, *SPHERE)
            with open(bad, "wb") as fh:
                fh.write(b"NOTTORI\0" + data[8:])
            with self.assertRaises(pat.PointsAsToriError) as cm:
                ddg_tools.pat_load_tori(bad, *SPHERE)
            self.assertIn("magic", str(cm.exception))
            with self.assertRaises(pat.PointsAsToriError):
                ddg_tools.pat_load_tori(os.path.join(d, "missing.tori"), *SPHERE)

    def test_author_layout_export_and_import(self):
        t = FIELD.export_tori()
        self.assertEqual(sorted(t), ["axes", "centers", "major", "minor"])
        n = FIELD.n
        self.assertEqual((len(t["centers"]), len(t["axes"]), len(t["major"]), len(t["minor"])), (3 * n, 2 * n, n, n))
        full = FIELD.tori
        self.assertEqual(sorted(full), ["axes", "axis_vectors", "centers", "coefficients", "major", "minor", "sigma"])
        # the axis angles are (theta, phi) = (atan2(uy, ux), acos(uz))
        for i in (0, 7, 100):
            ux, uy, uz = full["axis_vectors"][3 * i:3 * i + 3]
            self.assertAlmostEqual(t["axes"][2 * i], math.atan2(uy, ux), places=12)
            self.assertAlmostEqual(t["axes"][2 * i + 1], math.acos(max(-1.0, min(1.0, uz))), places=12)
        g = ddg_tools.pat_from_tori_dict(*SPHERE, t)
        self.assertLess(max_abs_diff(g.signed_distance(Q), FIELD.signed_distance(Q)), 1e-8)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_from_tori_dict(*SPHERE, {"centers": t["centers"]})
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_from_tori_dict(*SPHERE, {**t, "major": t["major"][:-1]})


class CrateItems(unittest.TestCase):
    """WP15: the sign-rule marker, reconstruct's eval report, smooth = 0 and the progress hook."""

    def test_sign_rule_marker_survives_a_round_trip_and_old_files_load(self):
        geo = ddg_tools.PointsAsTori(*SPHERE, sign_rule="geometric")
        self.assertEqual(geo.sign_rule, "geometric")
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "g.tori")
            geo.save_tori(path)
            self.assertEqual(ddg_tools.pat_load_tori(path, *SPHERE).sign_rule, "geometric")
            with open(path, "rb") as fh:
                data = bytearray(fh.read())
            data[12:16] = b"\0\0\0\0"  # a file from before the marker
            old = os.path.join(d, "old.tori")
            with open(old, "wb") as fh:
                fh.write(data)
            self.assertIsNone(ddg_tools.pat_load_tori(old, *SPHERE).sign_rule)
            data[12] = 1  # a file marked Reference is read as Reference, not as the new default
            with open(old, "wb") as fh:
                fh.write(data)
            self.assertEqual(ddg_tools.pat_load_tori(old, *SPHERE).sign_rule, "reference")
            ref = ddg_tools.PointsAsTori(*SPHERE, sign_rule="reference")
            self.assertEqual(ref.sign_rule, "reference")
            data[12] = 9
            with open(old, "wb") as fh:
                fh.write(data)
            with self.assertRaises(pat.PointsAsToriError):
                ddg_tools.pat_load_tori(old, *SPHERE)
        imported = ddg_tools.pat_from_tori_dict(*SPHERE, FIELD.export_tori())
        self.assertIsNone(imported.sign_rule)

    def test_reconstruct_fills_the_eval_report(self):
        f = ddg_tools.PointsAsTori(*SPHERE)
        self.assertIsNone(f.eval_report)
        v, faces, rep = f.reconstruct(resolution=16, return_report=True)
        er = f.eval_report
        self.assertGreater(er["queries"], 0)
        self.assertEqual(er["queries"], rep["eval"]["queries"])
        self.assertGreaterEqual(er["terms"], 32 * (er["queries"] - 1))

    def test_smooth_zero_is_the_sharp_boolean(self):
        a, b = FIELD, ddg_tools.PointsAsTori(*fibonacci_sphere(512, 0.5, (0.3, 0.0, 0.0)))
        for op in ("union", "intersect", "subtract"):
            sharp = getattr(a, op)(b).signed_distance(Q)
            self.assertEqual(getattr(a, op)(b, smooth=0.0).signed_distance(Q).tobytes(), sharp.tobytes())

    def test_progress_counts_and_cancels(self):
        p = ddg_tools.PatProgress()
        self.assertEqual((p.stage, p.done, p.total, p.cancelled), ("idle", 0, 0, False))
        f = ddg_tools.PointsAsTori(*SPHERE, progress=p)
        self.assertEqual(f.signed_distance(Q).tobytes(), FIELD.signed_distance(Q).tobytes())
        self.assertEqual((p.stage, p.done, p.total), ("evaluate", len(Q) // 3, len(Q) // 3))
        self.assertEqual(p.fraction, 1.0)
        p.cancel()
        self.assertTrue(p.cancelled)
        with self.assertRaises(pat.PointsAsToriError) as cm:
            f.signed_distance(Q)
        self.assertIn("cancelled", str(cm.exception))
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.PointsAsTori(*SPHERE, progress=p)
        # with_options detaches the handle unless given again
        self.assertEqual(f.with_options().signed_distance(Q).tobytes(), FIELD.signed_distance(Q).tobytes())
        with self.assertRaises(pat.PointsAsToriError):
            f.with_options(progress=p).signed_distance(Q)
        p.reset()
        self.assertEqual(f.signed_distance(Q).tobytes(), FIELD.signed_distance(Q).tobytes())
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.PointsAsTori(*SPHERE, progress="not a progress")


class Errors(unittest.TestCase):
    def raises(self, fragment, fn, *args, **kwargs):
        with self.assertRaises(pat.PointsAsToriError) as cm:
            fn(*args, **kwargs)
        self.assertIn(fragment, str(cm.exception))
        self.assertIsInstance(cm.exception, ValueError)
        self.assertIsInstance(cm.exception, ddg_tools.DdgError)

    def test_cloud_errors(self):
        p, n = SPHERE
        self.raises("512 positions but 510 normals", ddg_tools.PointsAsTori, p, n[:-6])
        self.raises("at least 2", ddg_tools.PointsAsTori, p[:3], n[:3])
        bad = array.array("d", p)
        bad[7] = float("nan")
        self.raises("position 2 is not finite", ddg_tools.PointsAsTori, bad, n)
        z = array.array("d", n)
        z[3:6] = array.array("d", [0, 0, 0])
        self.raises("normal 1 has length zero", ddg_tools.PointsAsTori, p, z)
        self.raises("whole number of 3-vectors", ddg_tools.PointsAsTori, p[:-1], n[:-1])
        self.raises("expected an object with the buffer protocol", ddg_tools.PointsAsTori, "abc", n)

    def test_option_errors(self):
        p, n = SPHERE
        self.raises("unknown option", ddg_tools.PointsAsTori, p, n, bogus=1)
        self.raises("not one of", ddg_tools.PointsAsTori, p, n, sign_rule="nope")
        self.raises("expected a string", ddg_tools.PointsAsTori, p, n, device=3)
        self.raises("threads", ddg_tools.PointsAsTori, p, n, threads=0)
        self.raises("radius", ddg_tools.PointsAsTori, p, n, neighborhood="radius")
        self.raises("lambda_value", ddg_tools.PointsAsTori, p, n, lambda_rule="constant")
        self.raises("k = 0", ddg_tools.PointsAsTori, p, n, knn_k=0)
        self.raises("area weights", ddg_tools.PointsAsTori, p, n, area_weights=[1.0, 2.0])
        self.raises("weights: cannot read", ddg_tools.PointsAsTori, p, n, weights="/nonexistent/w.safetensors")
        self.raises("not one of", ddg_tools.pat_boolean, FIELD, FIELD, "xor")
        self.raises("expected a PointsAsTori", ddg_tools.pat_reconstruct, "nope")

    def test_query_errors(self):
        self.raises("query 1 is not finite", FIELD.signed_distance, [0, 0, 0, float("inf"), 0, 0])
        self.raises("whole number of 3-vectors", FIELD.signed_distance, [0.0, 1.0])

    def test_no_cuda_in_a_cpu_wheel(self):
        info = ddg_tools.gpu_info()
        if info["cuda_compiled"]:
            self.skipTest("this wheel has CUDA")
        self.raises("cuda", ddg_tools.PointsAsTori, *SPHERE, device="cuda")


class Reconstruction(unittest.TestCase):
    def test_tiny_reconstruction(self):
        v, f = FIELD.reconstruct(0.0, 24)
        self.assertIsInstance(v, array.array)
        self.assertEqual(len(v) % 3, 0)
        self.assertEqual(len(f) % 3, 0)
        self.assertGreater(len(f), 100)
        nv = len(v) // 3
        self.assertLess(max(f), nv)
        radii = [math.sqrt(x * x + y * y + z * z) for x, y, z in rows(v)]
        self.assertLess(max(abs(r - 0.5) for r in radii), 0.03)

    def test_report_normals_and_options(self):
        v, f, nrm, rep = FIELD.reconstruct(resolution=20, normals=True, return_report=True)
        self.assertEqual(len(nrm), len(v))
        self.assertEqual(rep["mesh"]["status"], "manifold")
        self.assertEqual(rep["mesh"]["components"], 1)
        self.assertEqual(rep["mesh"]["euler_characteristic"], 2)
        self.assertEqual(rep["mesh"]["boundary_edges"], 0)
        self.assertEqual(rep["vertices"], len(v) // 3)
        for x, y, z in rows(nrm)[:50]:
            self.assertAlmostEqual(x * x + y * y + z * z, 1.0, places=6)
        # outward normals: along the radius
        for (px, py, pz), (nx, ny, nz) in list(zip(rows(v), rows(nrm)))[:50]:
            self.assertGreater(px * nx + py * ny + pz * nz, 0.0)
        # explicit box and cells, dense band, an offset level
        w, g = FIELD.reconstruct(0.1, [16, 16, 16], [-1, -1, -1, 1, 1, 1], band="dense", check_mesh=False)
        self.assertLess(max(abs(math.sqrt(x * x + y * y + z * z) - 0.6) for x, y, z in rows(w)), 0.04)
        a = ddg_tools.pat_reconstruct(FIELD, 0.0, 20)
        self.assertEqual(a[0].tobytes(), v.tobytes())
        self.assertEqual(a[1].tobytes(), f.tobytes())
        for bad in ({"resolution": 0}, {"resolution": [1, 2]}, {"bbox": [0, 1, 2]}, {"band": "wide"}):
            with self.assertRaises(pat.PointsAsToriError):
                FIELD.reconstruct(**bad)

    def test_level_set_sampling(self):
        p, n, rep = FIELD.sample_level_set(0.1, 300, seed=3, return_report=True)
        self.assertEqual(len(p), 900)
        self.assertEqual(rep["produced"], 300)
        for (x, y, z), (nx, ny, nz) in zip(rows(p), rows(n)):
            self.assertAlmostEqual(math.sqrt(x * x + y * y + z * z), 0.6, delta=0.04)
            self.assertAlmostEqual(nx * nx + ny * ny + nz * nz, 1.0, places=6)
        q, _ = FIELD.sample_level_set(0.1, 300, seed=3)
        self.assertEqual(p.tobytes(), q.tobytes())
        r, _ = pat.sample_level_set(FIELD, 0.1, 300, seed=4)
        self.assertNotEqual(p.tobytes(), r.tobytes())


class Csg(unittest.TestCase):
    def setUp(self):
        self.other = ddg_tools.PointsAsTori(*fibonacci_sphere(512, 0.5, (0.4, 0.0, 0.0)))

    def test_boolean_values(self):
        a, b = FIELD, self.other
        da, db = a.signed_distance(Q), b.signed_distance(Q)
        for op, fn in (("union", min), ("intersect", max), ("subtract", lambda x, y: max(x, -y))):
            c = ddg_tools.pat_boolean(a, b, op)
            self.assertIsInstance(c, ddg_tools.PatCsgField)
            expect = [fn(x, y) for x, y in zip(da, db)]
            self.assertLess(max_abs_diff(c.signed_distance(Q), expect), 1e-12)
            m = getattr(a, op)(b)
            self.assertEqual(m.signed_distance(Q).tobytes(), c.signed_distance(Q).tobytes())
        s = a.union(b, smooth=0.2)
        self.assertTrue(all(x <= y + 1e-12 for x, y in zip(s.signed_distance(Q), a.union(b).signed_distance(Q))))
        d, g = a.union(b).distance_and_gradient(Q)
        self.assertEqual(len(g), len(Q))
        self.assertEqual(a.union(b).gradient(Q).tobytes(), g.tobytes())
        self.assertEqual(len(a.union(b).bbox), 6)

    def test_offset_and_reconstruction(self):
        o = FIELD.offset(0.1)
        self.assertLess(max_abs_diff(o.signed_distance(Q), [x - 0.1 for x in FIELD.signed_distance(Q)]), 1e-12)
        v, f = ddg_tools.pat_offset(FIELD, 0.1).reconstruct(0.0, 20)
        self.assertLess(max(abs(math.sqrt(x * x + y * y + z * z) - 0.6) for x, y, z in rows(v)), 0.04)
        u = ddg_tools.pat_boolean(FIELD, self.other, "union")
        v, f, rep = u.reconstruct(resolution=32, return_report=True)
        self.assertEqual(rep["mesh"]["components"], 1)
        self.assertEqual(rep["mesh"]["euler_characteristic"], 2)
        p, n = u.sample_level_set(0.0, 100)
        self.assertEqual(len(p), 300)
        self.assertLess(max_abs_diff(u.offset(0.05).offset(-0.05).signed_distance(Q), u.signed_distance(Q)), 1e-12)

    def test_csg_errors(self):
        with self.assertRaises(pat.PointsAsToriError):
            FIELD.union(FIELD, smooth=-1.0)
        with self.assertRaises(pat.PointsAsToriError):
            FIELD.offset(float("nan"))
        with self.assertRaises(pat.PointsAsToriError):
            FIELD.union("sphere")
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_boolean(FIELD, 1, "union")


class Applications(unittest.TestCase):
    def test_morphology(self):
        f, rep = ddg_tools.pat_morphology(FIELD, "dilate", 0.1, count=400, resolution=24)
        self.assertIsInstance(f, ddg_tools.PointsAsTori)
        self.assertEqual(len(rep["steps"]), 1)
        self.assertEqual(rep["steps"][0]["level"], 0.1)
        d = f.signed_distance([0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0])
        self.assertAlmostEqual(d[0], 0.0, delta=0.04)
        self.assertLess(d[1], -0.45)
        self.assertAlmostEqual(d[2], 0.2, delta=0.05)
        e, rep = pat.morphology(FIELD, "close", 0.05, count=400, resolution=24)
        self.assertEqual([s["level"] for s in rep["steps"]], [0.05, -0.05])
        self.assertAlmostEqual(e.signed_distance([0.5, 0.0, 0.0])[0], 0.0, delta=0.04)
        for bad in (("grow", 0.1), ("erode", -0.1), ("erode", float("nan"))):
            with self.assertRaises(pat.PointsAsToriError):
                ddg_tools.pat_morphology(FIELD, *bad)

    def test_orientation(self):
        p, n = SPHERE
        wrong = array.array("d", n)
        flipped = set(range(0, 512, 9))
        for i in flipped:
            for k in range(3):
                wrong[3 * i + k] = -wrong[3 * i + k]
        fixed, rep = ddg_tools.pat_orient_normals(p, wrong)
        self.assertEqual(len(fixed), len(n))
        self.assertEqual(rep["n"], 512)
        self.assertTrue(rep["converged"])
        bad = sum(1 for i in range(512) if sum(fixed[3 * i + k] * n[3 * i + k] for k in range(3)) < 0)
        self.assertLessEqual(bad, 5)
        self.assertGreater(rep["flipped_total"], 0)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_orient_normals(p, n, vote_k=0)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_orient_normals(p, n, epsilon=1.0, epsilon_scale=1.0)


class Baselines(unittest.TestCase):
    P, N = SPHERE
    INSIDE = [0.0, 0.0, 0.0, 0.1, 0.1, 0.0]
    OUTSIDE = [1.0, 0.0, 0.0, 0.0, 1.5, 0.0]
    QS = array.array("d", INSIDE + OUTSIDE)

    def test_point_areas(self):
        a = ddg_tools.pat_point_areas(self.P, self.N)
        self.assertEqual(len(a), 512)
        self.assertAlmostEqual(sum(a), 4 * math.pi * 0.25, delta=0.4)
        a2, rep = ddg_tools.pat_point_areas(self.P, self.N, k=12, return_report=True)
        self.assertEqual(rep["n"], 512)
        self.assertEqual(len(a2), 512)

    def test_winding_numbers(self):
        w, rep = ddg_tools.pat_gwn(self.P, self.N, self.QS, return_report=True)
        self.assertEqual(rep["queries"], 4)
        for x in w[:2]:
            self.assertAlmostEqual(x, 1.0, delta=0.15)
        for x in w[2:]:
            self.assertAlmostEqual(x, 0.0, delta=0.15)
        r = ddg_tools.pat_gwn(self.P, self.N, self.QS, regularization=0.1)
        self.assertAlmostEqual(r[0], 1.0, delta=0.2)
        f = ddg_tools.pat_fwn(self.P, self.N, self.QS)
        self.assertLess(max_abs_diff(f, w), 0.02)
        areas = ddg_tools.pat_point_areas(self.P, self.N)
        self.assertEqual(ddg_tools.pat_gwn(self.P, self.N, self.QS, areas=areas).tobytes(), w.tobytes())
        v, fc = octahedron()
        m = ddg_tools.pat_fwn_mesh(v, fc, [0, 0, 0, 3, 0, 0, 0.2, 0.1, 0.1])
        self.assertAlmostEqual(m[0], 1.0, places=2)
        self.assertAlmostEqual(m[1], 0.0, places=2)
        self.assertAlmostEqual(m[2], 1.0, places=2)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_gwn(self.P, self.N, self.QS, regularization=-1.0)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_fwn(self.P, self.N, self.QS, order=3)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_fwn_mesh(v, fc[:-1], [0, 0, 0])

    def test_kernel_baselines(self):
        exact = exact_sphere(self.QS)
        for fn in (ddg_tools.pat_shc, ddg_tools.pat_sspd):
            v = fn(self.P, self.N, self.QS)
            self.assertEqual(len(v), 4)
            self.assertTrue(all(x < 0 for x in v[:2]) and all(x > 0 for x in v[2:]), (fn, list(v)))
            self.assertLess(max_abs_diff(v, exact), 0.2)
            v2, g = fn(self.P, self.N, self.QS, gradient=True)
            self.assertEqual(v2.tobytes(), v.tobytes())
            self.assertEqual(len(g), 12)
            v3, g3, rep = fn(self.P, self.N, self.QS, gradient=True, return_report=True)
            self.assertEqual(rep["queries"], 4)
        lse = ddg_tools.pat_lse(self.P, self.N, self.QS)
        self.assertTrue(all(x > 0 for x in lse))
        self.assertAlmostEqual(lse[2], 0.5, delta=0.1)
        paper = ddg_tools.pat_shc(self.P, self.N, self.QS, paper=True)
        self.assertEqual(len(paper), 4)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_shc(self.P, self.N, self.QS, bogus=1)

    def test_mesh_truth(self):
        v, f = octahedron()
        t = ddg_tools.PatMeshTruth(v, f)
        self.assertTrue(t.mesh_report["closed"])
        self.assertAlmostEqual(t.mesh_report["signed_volume"], 4.0 / 3.0, places=9)
        q = [0, 0, 0, 2, 0, 0, 0.3, 0.1, 0.0]
        d = t.signed_distance(q)
        self.assertAlmostEqual(d[0], -1 / math.sqrt(3), places=12)
        self.assertAlmostEqual(d[1], 1.0, places=12)
        self.assertLess(d[2], 0.0)
        dd, w = t.evaluate(q)
        self.assertEqual(dd.tobytes(), d.tobytes())
        self.assertAlmostEqual(w[0], 1.0, places=6)
        self.assertAlmostEqual(w[1], 0.0, places=6)
        for kind in ("exact", "fast"):
            tt = ddg_tools.PatMeshTruth(v, f, sign=kind)
            self.assertLess(max_abs_diff(tt.signed_distance(q), d), 1e-9)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.PatMeshTruth(v, f[:-3])  # open: a face missing
        open_ok = ddg_tools.PatMeshTruth(v, f[:-3], require_closed=False)
        self.assertFalse(open_ok.mesh_report["closed"])
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.PatMeshTruth(v, array.array("i", [0, 1, 99]))

    def test_truth_against_pat_on_a_mesh_sphere(self):
        # the field against the exact SDF of the same sphere as a mesh (octahedron not round: use the analytic)
        err = max_abs_diff(FIELD.signed_distance(Q), exact_sphere(Q))
        self.assertLess(err, 0.03)


class Classical2D(unittest.TestCase):
    def test_circle(self):
        p, n = circle(256)
        f = ddg_tools.Pat2d(p, n)
        self.assertEqual(f.build_report["n"], 256)
        d = f.signed_distance([0, 0, 2, 0, 1, 0])
        self.assertAlmostEqual(d[0], -1.0, delta=0.05)
        self.assertAlmostEqual(d[1], 1.0, delta=0.05)
        self.assertAlmostEqual(d[2], 0.0, delta=0.05)
        d2, g = f.distance_and_gradient([2, 0])
        self.assertAlmostEqual(g[0], 1.0, delta=0.03)
        lines = f.contours([-1.5, -1.5, 1.5, 1.5], 60, 60)
        self.assertEqual(len(lines), 1)
        pts, closed = lines[0]
        self.assertTrue(closed)
        for i in range(len(pts) // 2):
            self.assertAlmostEqual(math.hypot(pts[2 * i], pts[2 * i + 1]), 1.0, delta=0.03)
        self.assertTrue(ddg_tools.Pat2d(p, n, paper=True).signed_distance([0, 0])[0] < -0.5)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.Pat2d(p, n[:-2])
        with self.assertRaises(pat.PointsAsToriError):
            f.contours([0, 0, 1], 5, 5)
        with self.assertRaises(pat.PointsAsToriError):
            f.contours([0, 0, 1, 1], 1, 5)


class Determinism(unittest.TestCase):
    def test_threads_and_repeats(self):
        a = ddg_tools.PointsAsTori(*SPHERE, threads=1)
        b = ddg_tools.PointsAsTori(*SPHERE, threads=4)
        c = ddg_tools.PointsAsTori(*SPHERE)
        for k in ("centers", "axes", "major", "minor", "coefficients"):
            self.assertEqual(a.tori[k].tobytes(), b.tori[k].tobytes(), k)
            self.assertEqual(a.tori[k].tobytes(), c.tori[k].tobytes(), k)
        for fn in ("signed_distance", "gradient", "laplacian"):
            self.assertEqual(getattr(a, fn)(Q).tobytes(), getattr(b, fn)(Q).tobytes(), fn)
        ma = a.reconstruct(resolution=20)
        mb = b.reconstruct(resolution=20, threads=3)
        self.assertEqual(ma[0].tobytes(), mb[0].tobytes())
        self.assertEqual(ma[1].tobytes(), mb[1].tobytes())

    def test_concurrent_evaluation_releases_the_gil(self):
        import threading

        big = grid_queries(-0.9, 0.9, 40)
        out = []
        ts = [threading.Thread(target=lambda: out.append(FIELD.signed_distance(big).tobytes())) for _ in range(3)]
        for t in ts:
            t.start()
        for t in ts:
            t.join()
        self.assertEqual(len(set(out)), 1)


class Rendering(unittest.TestCase):
    """WP11: sphere-traced images (Fig. 12)."""

    EYE, TARGET = (0.0, 0.4, 2.5), (0.0, 0.0, 0.0)

    def render(self, field=None, **kw):
        return ddg_tools.pat_render(field or FIELD, self.EYE, self.TARGET, fov=30.0, width=48, height=32, **kw)

    def test_image_shape_and_report(self):
        img, rep = self.render(return_report=True, threads=2)
        self.assertIsInstance(img, array.array)
        self.assertEqual(img.typecode, "B")
        self.assertEqual(len(img), 48 * 32 * 4)
        self.assertEqual(rep["rays"], 48 * 32)
        self.assertGreater(rep["hits"], 50)
        self.assertEqual(rep["step_limit"], 0)
        self.assertLessEqual(rep["steps_max"], 128)
        centre = (16 * 48 + 24) * 4
        self.assertEqual(img[centre + 3], 255)
        self.assertNotEqual(tuple(img[centre : centre + 3]), (255, 255, 255))
        if numpy is not None:
            self.assertEqual(numpy.asarray(img).reshape(32, 48, 4).shape, (32, 48, 4))

    def test_bitwise_at_any_thread_count(self):
        a = self.render(threads=1)
        b = self.render(threads=4)
        c = self.render(threads=3)
        self.assertEqual(a.tobytes(), b.tobytes())
        self.assertEqual(a.tobytes(), c.tobytes())

    def test_offset_boolean_and_cut_plane(self):
        big = FIELD.offset(0.1)
        _, r0 = self.render(return_report=True)
        _, r1 = self.render(big, return_report=True)
        self.assertGreater(r1["hits"], r0["hits"])
        other = ddg_tools.PointsAsTori(*fibonacci_sphere(256, 0.3, (0.4, 0.0, 0.0)))
        u = pat.boolean(FIELD, other, "union")
        _, ru = self.render(u, return_report=True)
        self.assertGreater(ru["hits"], 50)
        _, rc = self.render(return_report=True, cut_axis=2, cut_offset=-0.2, sky="white")
        self.assertGreater(rc["cut_plane_pixels"], 0)

    def test_bad_input_raises(self):
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_render(FIELD, self.TARGET, self.TARGET)
        with self.assertRaises(pat.PointsAsToriError):
            self.render(safety=2.0)
        with self.assertRaises(pat.PointsAsToriError):
            self.render(nonsense=1)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_render(FIELD, self.EYE, self.TARGET, width=0, height=4)
        with self.assertRaises(pat.PointsAsToriError):
            ddg_tools.pat_render(FIELD, (1, 2), self.TARGET)
        self.assertIs(pat.render, ddg_tools.pat_render)


class Defaults(unittest.TestCase):
    def test_stub_defaults_equal_the_library(self):
        path = os.path.join(HERE, "..", "ddg_tools.pyi")
        with open(path, encoding="utf-8") as fh:
            tree = ast.parse(fh.read())
        cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "PointsAsTori")
        fn = next(i for i in cls.body if isinstance(i, ast.FunctionDef) and i.name == "__init__")
        a = fn.args
        stub = {arg.arg: ast.literal_eval(d) for arg, d in zip(a.kwonlyargs, a.kw_defaults) if d is not None}
        self.assertEqual(stub, ddg_tools.pat_defaults())

    def test_defaults_are_the_decided_ones(self):
        d = ddg_tools.pat_defaults()
        self.assertEqual((d["neighborhood"], d["knn_k"]), ("knn", 32))
        self.assertEqual((d["lambda_rule"], d["lambda_c"]), ("per_query", 64.0))
        self.assertEqual(d["gpu_matmul"], "tf32")
        self.assertEqual((d["sign_rule"], d["eigen_rule"]), ("geometric", "reference"))
        self.assertEqual(d["gradient"], "frozen")
        self.assertEqual(d["precision"], "f64")
        self.assertEqual(d["device"], "cpu")
        self.assertEqual(ddg_tools.PointsAsTori(*SPHERE).options, {**FIELD.options})

    def test_every_pat_name_is_in_the_stub(self):
        path = os.path.join(HERE, "..", "ddg_tools.pyi")
        with open(path, encoding="utf-8") as fh:
            tree = ast.parse(fh.read())
        names = {n.name for n in tree.body if isinstance(n, (ast.ClassDef, ast.FunctionDef))}
        for name in dir(ddg_tools):
            if name.startswith("pat_") or name in ("PointsAsTori", "PatCsgField", "PatMeshTruth", "Pat2d", "PointsAsToriError", "PatProgress"):
                self.assertIn(name, names)


class Golden(unittest.TestCase):
    """Against the authors' program (`runs_tmp/pat-ref`, see its README): CPU network, 32 nearest
    neighbors per query. NOTES-WP1 measured end-to-end max 2.6e-7 on the grid."""

    @classmethod
    def setUpClass(cls):
        if not os.path.isdir(os.path.join(GOLDEN, "cpu_run1")):
            raise unittest.SkipTest("runs_tmp/pat-ref is absent")

    def case(self, name):
        inp = os.path.join(GOLDEN, "inputs", name)
        run = os.path.join(GOLDEN, "cpu_run1", name)
        p, _ = read_npy(os.path.join(inp, "points.npy"))
        n, _ = read_npy(os.path.join(inp, "normals.npy"))
        q, _ = read_npy(os.path.join(GOLDEN, "inputs", "queries_random.npy"))
        return p, n, q, run

    def check(self, name, tol_sdf, tol_grad):
        p, n, q, run = self.case(name)
        f = ddg_tools.PointsAsTori(p, n, feature_order="reference", sign_rule="reference", normal_policy="normalize")
        sdf, _ = read_npy(os.path.join(run, "query_sdf_accel.npy"))
        grad, _ = read_npy(os.path.join(run, "query_grad_accel.npy"))
        d, g = f.distance_and_gradient(q)
        self.assertLess(max_abs_diff(d, sdf), tol_sdf)
        self.assertLess(max_abs_diff(g, grad), tol_grad)
        # tori
        t = f.export_tori()
        for key, file in (("centers", "tori_centers"), ("major", "tori_major_radius"), ("minor", "tori_minor_radius")):
            ref, _ = read_npy(os.path.join(run, file + ".npy"))
            scale = max(1.0, max(abs(x) for x in ref if abs(x) < 1e6))
            diffs = sorted(abs(x - y) / max(1.0, abs(y)) for x, y in zip(t[key], ref))
            self.assertLess(diffs[int(0.99 * len(diffs))], 1e-3, key)
        # the same through the author-layout import
        h = ddg_tools.pat_from_tori_dict(p, n, {k: read_npy(os.path.join(run, f"tori_{k2}.npy"))[0] for k, k2 in (("centers", "centers"), ("axes", "axes_angles"), ("major", "major_radius"), ("minor", "minor_radius"))}, normal_policy="normalize")
        self.assertLess(max_abs_diff(h.signed_distance(q), sdf), 1e-9)
        self.assertLess(max_abs_diff(h.gradient(q), grad), 1e-7)

    def test_synth_sphere(self):
        self.check("synth_sphere", 1e-5, 1e-3)

    def test_tooth(self):
        self.check("tooth", 1e-5, 1e-3)

    def test_holey_bunny_grid_subset(self):
        # the saddle tori of the bunny make it the worst case; same tolerance as the Rust parity test
        p, n, q, run = self.case("holey-bunny")
        sdf, _ = read_npy(os.path.join(run, "query_sdf_accel.npy"))
        f = ddg_tools.PointsAsTori(p, n, feature_order="reference", sign_rule="reference", normal_policy="normalize")
        d = f.signed_distance(q)
        diffs = sorted(abs(x - y) for x, y in zip(d, sdf))
        self.assertLess(diffs[len(diffs) // 2], 1e-6)
        self.assertLess(diffs[-1], 1e-3)


class Gpu(unittest.TestCase):
    def test_device_cuda_matches_cpu(self):
        info = ddg_tools.gpu_info()
        if not (info["cuda_compiled"] and info["available"]):
            self.skipTest("no CUDA wheel or device")
        try:
            apps = subprocess.run(["nvidia-smi", "--query-compute-apps=pid", "--format=csv,noheader"], capture_output=True, text=True, timeout=20).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            apps = ""
        if apps:
            self.skipTest("a process is using the GPU")
        g = ddg_tools.PointsAsTori(*SPHERE, device="cuda", gpu_matmul="split3")
        self.assertTrue(g.build_report["predict_device"].startswith("cuda"))
        self.assertLess(max_abs_diff(g.signed_distance(Q), FIELD.signed_distance(Q)), 1e-4)
        t = ddg_tools.PointsAsTori(*SPHERE, device="auto")
        self.assertLess(max_abs_diff(t.signed_distance(Q), FIELD.signed_distance(Q)), 1e-2)

    def test_gpu_info_never_raises(self):
        info = ddg_tools.gpu_info()
        self.assertIn("cuda_compiled", info)


if __name__ == "__main__":
    unittest.main()
