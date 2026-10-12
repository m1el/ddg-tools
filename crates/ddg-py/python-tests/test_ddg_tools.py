"""Tests of the built wheel of `ddg_tools`, without Blender.

Run through `scripts/py-test.sh`, which installs the wheel into throw-away environments with
`uv run` (once without numpy, once with) and sets DDG_WHEEL to the wheel under test.
Standard library only; the numpy tests are skipped when numpy is not installed.
"""

import array
import hashlib
import json
import math
import os
import struct
import sys
import unittest
import zipfile

import ddg_tools

try:
    import numpy
except ImportError:
    numpy = None


def cube():
    """The unit cube as 6 outward-facing quads, fan-triangulated: (positions, corner_vertices,
    polygon_starts, triangles) as lists."""
    positions = []
    for i in range(8):
        positions += [float(i & 1), float((i >> 1) & 1), float((i >> 2) & 1)]
    quads = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    corner_vertices = [v for quad in quads for v in quad]
    polygon_starts = [4 * p for p in range(6)]
    triangles = []
    for start in polygon_starts:
        triangles += [start, start + 1, start + 2, start, start + 2, start + 3]
    return positions, corner_vertices, polygon_starts, triangles


def as_arrays(mesh, real='f', index='i'):
    positions, corner_vertices, polygon_starts, triangles = mesh
    return (
        array.array(real, positions),
        array.array(index, corner_vertices),
        array.array(index, polygon_starts),
        array.array(index, triangles),
    )


def three_fins():
    """Three triangles on the edge (0, 1): not a manifold."""
    positions = [0, 0, 0, 1, 0, 0, 0.5, 1, 0, 0.5, -1, 0, 0.5, 0, 1]
    return [float(x) for x in positions], [0, 1, 2, 1, 0, 3, 0, 1, 4], [0, 3, 6], list(range(9))


class ModuleTest(unittest.TestCase):
    def test_loaded_module_is_the_wheel_under_test(self):
        wheel = os.environ.get("DDG_WHEEL")
        if not wheel:
            self.skipTest("DDG_WHEEL is not set")
        package = os.path.dirname(ddg_tools.__file__)
        with zipfile.ZipFile(wheel) as archive:
            names = [n for n in archive.namelist() if n.endswith(".so") or n.endswith(".pyd")]
            self.assertEqual(len(names), 1, names)
            expected = hashlib.sha256(archive.read(names[0])).hexdigest()
        with open(os.path.join(package, os.path.basename(names[0])), "rb") as fh:
            self.assertEqual(hashlib.sha256(fh.read()).hexdigest(), expected)
        # One binary for every Python from 3.11 on.
        self.assertIn("-cp311-abi3-", os.path.basename(wheel))
        self.assertIn(".abi3.", names[0])

    def test_attributes(self):
        self.assertRegex(ddg_tools.__version__, r"^\d+\.\d+\.\d+")
        self.assertIsInstance(ddg_tools.PIPELINE, str)
        self.assertIsInstance(ddg_tools.IS_PLACEHOLDER, bool)
        for name in ("InputError", "MeshError", "SolveError"):
            error = getattr(ddg_tools, name)
            self.assertTrue(issubclass(error, ddg_tools.DdgError), name)
            self.assertEqual(error.__module__, "ddg_tools")
        self.assertTrue(issubclass(ddg_tools.DdgError, Exception))
        self.assertIn("lambda_=0.1", ddg_tools.unwrap.__text_signature__)
        self.assertIn("polygon corner", ddg_tools.unwrap.__doc__)

    def test_pipeline_is_announced(self):
        result = ddg_tools.unwrap(*as_arrays(cube()))
        self.assertEqual(result.report["pipeline"], ddg_tools.PIPELINE)
        self.assertIs(result.report["placeholder"], ddg_tools.IS_PLACEHOLDER)
        # The real pipeline is linked; the flag exists so that a stand-in could say so.
        self.assertIs(ddg_tools.IS_PLACEHOLDER, False)
        self.assertEqual(ddg_tools.PIPELINE, "ddg-uv")


class OutputTest(unittest.TestCase):
    def setUp(self):
        self.mesh = cube()
        self.result = ddg_tools.unwrap(*as_arrays(self.mesh))

    def test_result_arrays(self):
        r = self.result
        expected = {
            "uvs": ('f', 2 * 24),
            "polygon_charts": ('i', 6),
            "cut_edges": ('i', None),
            "cone_vertices": ('i', None),
            "cone_angles": ('d', None),
        }
        for name, (typecode, length) in expected.items():
            value = getattr(r, name)
            self.assertIs(type(value), array.array, name)
            self.assertEqual(value.typecode, typecode, name)
            if length is not None:
                self.assertEqual(len(value), length, name)
        self.assertEqual(len(r.cut_edges) % 2, 0)
        self.assertEqual(len(r.cone_vertices), len(r.cone_angles))
        self.assertTrue(all(math.isfinite(x) for x in r.uvs))
        self.assertIn("UnwrapResult(corners=24, polygons=6", repr(r))
        with self.assertRaises(AttributeError):
            r.uvs = None

    def test_report_is_a_dict_of_plain_values(self):
        report = self.result.report
        self.assertIs(type(report), dict)
        for key, value in report.items():
            self.assertIs(type(key), str)
            if key == "warnings":
                self.assertIs(type(value), list)
                self.assertTrue(all(type(w) is str for w in value))
            else:
                self.assertIn(type(value), (bool, int, float, str), key)
        json.dumps(report)
        for key, kind in {
            "pipeline": str, "placeholder": bool,
            "n_vertices": int, "n_polygons": int, "n_corners": int, "n_triangles": int,
            "n_seam_edges": int, "n_pinned": int, "n_prescribed": int, "n_forbidden": int,
            "lambda": float, "positive_cones_only": bool,
            "n_charts": int, "n_cones": int, "n_cut_edges": int,
            "n_degenerate_triangles": int, "n_collapsed_uv_triangles": int,
            "n_flipped_uv_triangles": int, "measured_triangles": int,
            "measured_angle_error_max": float, "measured_angle_error_mean": float,
            "measured_log_scale_std": float, "measured_log_scale_min": float,
            "measured_log_scale_max": float,
            "seconds_decode": float, "seconds_validate": float, "seconds_pipeline": float,
            "seconds_measure": float, "seconds_encode": float, "seconds_total": float,
            "warnings": list,
        }.items():
            self.assertIn(key, report)
            self.assertIs(type(report[key]), kind, key)
        self.assertEqual(
            (report["n_vertices"], report["n_polygons"], report["n_corners"],
             report["n_triangles"]),
            (8, 6, 24, 12))
        self.assertEqual(report["n_cones"], len(self.result.cone_vertices))
        self.assertEqual(report["n_cut_edges"], len(self.result.cut_edges) // 2)
        self.assertEqual(report["n_charts"], len(set(self.result.polygon_charts)))
        self.assertGreaterEqual(report["seconds_total"], report["seconds_pipeline"])

    def test_deterministic_and_inputs_untouched(self):
        arrays = as_arrays(self.mesh)
        before = [a.tobytes() for a in arrays]
        again = ddg_tools.unwrap(*arrays)
        self.assertEqual([a.tobytes() for a in arrays], before)
        for name in ("uvs", "polygon_charts", "cut_edges", "cone_vertices", "cone_angles"):
            self.assertEqual(getattr(again, name), getattr(self.result, name), name)


class BufferTypeTest(unittest.TestCase):
    """Every kind of buffer gives the same result as the reference call with array('f'/'i')."""

    def setUp(self):
        self.mesh = cube()
        self.reference = ddg_tools.unwrap(*as_arrays(self.mesh)).uvs

    def check(self, positions, corner_vertices, polygon_starts, triangles):
        result = ddg_tools.unwrap(positions, corner_vertices, polygon_starts, triangles)
        self.assertEqual(result.uvs, self.reference)

    def test_array_typecodes(self):
        for real in "fd":
            for index in "iIlLqQhH":
                with self.subTest(real=real, index=index):
                    self.check(*as_arrays(self.mesh, real, index))

    def test_bytes_bytearray_memoryview(self):
        arrays = as_arrays(self.mesh)
        self.check(*[a.tobytes() for a in arrays])
        self.check(*[bytearray(a.tobytes()) for a in arrays])
        self.check(*[memoryview(a) for a in arrays])
        self.check(*[memoryview(a.tobytes()) for a in arrays])
        self.check(*[memoryview(a.tobytes()).cast(a.typecode) for a in arrays])
        # Mixed, as a caller without numpy might do.
        self.check(arrays[0].tobytes(), arrays[1], memoryview(arrays[2]), bytearray(arrays[3]))

    def test_struct_packed_with_explicit_byte_order(self):
        positions, corner_vertices, polygon_starts, triangles = self.mesh
        native = "<" if sys.byteorder == "little" else ">"
        self.check(
            memoryview(struct.pack(native + "%df" % len(positions), *positions)).cast("f"),
            struct.pack(native + "%di" % len(corner_vertices), *corner_vertices),
            struct.pack("=%di" % len(polygon_starts), *polygon_starts),
            struct.pack("@%di" % len(triangles), *triangles),
        )

    def test_two_dimensional_memoryview(self):
        arrays = as_arrays(self.mesh)
        positions = memoryview(arrays[0]).cast("B").cast("f", (8, 3))
        triangles = memoryview(arrays[3]).cast("B").cast("i", (12, 3))
        self.check(positions, arrays[1], arrays[2], triangles)
        bad = memoryview(arrays[0]).cast("B").cast("f", (12, 2))
        with self.assertRaisesRegex(
                ddg_tools.InputError,
                r"positions: buffer has shape \(12, 2\); expected a flat buffer or one of "
                r"shape \(n, 3\)"):
            ddg_tools.unwrap(bad, arrays[1], arrays[2], arrays[3])

    def test_not_a_buffer_is_a_type_error(self):
        arrays = as_arrays(self.mesh)
        for position, value in enumerate([list(arrays[0]), None, 3, "abc"]):
            with self.subTest(value=type(value).__name__):
                with self.assertRaisesRegex(TypeError, "positions: expected an object that "
                                                       "supports the buffer protocol"):
                    ddg_tools.unwrap(value, arrays[1], arrays[2], arrays[3])
        with self.assertRaisesRegex(TypeError, "triangles: .* got list"):
            ddg_tools.unwrap(arrays[0], arrays[1], arrays[2], list(arrays[3]))
        with self.assertRaisesRegex(TypeError, "seam_edges: .* got tuple"):
            ddg_tools.unwrap(*arrays, seam_edges=(0, 1))

    def test_wrong_element_types(self):
        arrays = as_arrays(self.mesh)
        with self.assertRaisesRegex(ddg_tools.InputError, "positions: unsupported buffer element"):
            ddg_tools.unwrap(array.array('i', [0] * 24), arrays[1], arrays[2], arrays[3])
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    "corner_vertices: unsupported buffer element"):
            ddg_tools.unwrap(arrays[0], array.array('f', [0.0] * 24), arrays[2], arrays[3])
        with self.assertRaisesRegex(ddg_tools.InputError, "triangles: unsupported"):
            ddg_tools.unwrap(arrays[0], arrays[1], arrays[2], array.array('d', [0.0] * 36))
        with self.assertRaisesRegex(
                ddg_tools.InputError,
                "polygon_starts: a byte buffer is read as native-endian int32, but its length "
                "7 is not a multiple of 4"):
            ddg_tools.unwrap(arrays[0], arrays[1], b"\0" * 7, arrays[3])


class ShapeAndValueTest(unittest.TestCase):
    def setUp(self):
        self.mesh = cube()

    def unwrap_with(self, index, values, real='f', **kwargs):
        mesh = list(self.mesh)
        mesh[index] = values
        return ddg_tools.unwrap(*as_arrays(mesh, real), **kwargs)

    def assert_input_error(self, pattern, index, values, **kwargs):
        with self.assertRaisesRegex(ddg_tools.InputError, pattern) as caught:
            self.unwrap_with(index, values, **kwargs)
        self.assertIsInstance(caught.exception, ddg_tools.DdgError)

    def test_lengths(self):
        positions, corner_vertices, polygon_starts, triangles = self.mesh
        self.assert_input_error(
            r"positions: 23 items is not a multiple of 3 \(3 per vertex\)", 0, positions[:-1])
        self.assert_input_error(
            r"triangles: 35 items is not a multiple of 3 \(3 per triangle\)", 3, triangles[:-1])
        self.assert_input_error(
            r"triangles: got 11 triangles, expected 12 \(24 corners - 2 x 6 polygons\)",
            3, triangles[:-3])
        self.assert_input_error("polygon_starts is empty: the mesh has no faces", 2, [])
        self.assert_input_error(
            r"polygon 5 has 2 corners \(corners 20\.\.22\); a polygon needs at least 3",
            1, corner_vertices[:-2])

    def test_values(self):
        positions, corner_vertices, polygon_starts, triangles = self.mesh
        self.assert_input_error(
            r"corner_vertices\[3\] = 8 is out of range: there are 8 vertices",
            1, corner_vertices[:3] + [8] + corner_vertices[4:])
        self.assert_input_error(
            r"corner_vertices\[0\] = -1 is negative", 1, [-1] + corner_vertices[1:])
        self.assert_input_error(
            r"polygon_starts\[0\] = 1, but the first polygon must start at corner 0",
            2, [1] + polygon_starts[1:])
        self.assert_input_error(
            r"triangles\[2\] = 24 is out of range: there are 24 polygon corners",
            3, triangles[:2] + [24] + triangles[3:])
        self.assert_input_error(
            "triangles: triangle 0 = corners \\(0, 1, 4\\) spans polygons 0, 0, 1",
            3, [0, 1, 4] + triangles[3:])
        self.assert_input_error(
            "polygon 0 \\(4 corners\\) is not triangulated correctly: triangle 0 runs against",
            3, [0, 2, 1] + triangles[3:])
        for bad in (math.nan, math.inf):
            self.assert_input_error(
                r"positions\[4\] = (NaN|inf) is not a finite number",
                0, positions[:4] + [bad] + positions[5:])
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    r"corner_vertices\[0\] = 1099511627776 does not fit a "
                                    r"32-bit index"):
            ddg_tools.unwrap(
                array.array('f', positions), array.array('q', [2 ** 40] + corner_vertices[1:]),
                array.array('i', polygon_starts), array.array('i', triangles))

    def test_options(self):
        arrays = as_arrays(self.mesh)
        for bad in (-0.5, math.nan, math.inf):
            with self.assertRaisesRegex(ddg_tools.InputError, "lambda_ = .*: lambda must be"):
                ddg_tools.unwrap(*arrays, lambda_=bad)
        with self.assertRaises(TypeError):
            ddg_tools.unwrap(*arrays, lambda_="a lot")
        with self.assertRaises(TypeError):
            ddg_tools.unwrap(*arrays, array.array('i'))  # options are keyword-only
        with self.assertRaises(TypeError):
            ddg_tools.unwrap(*arrays, no_such_option=1)
        result = ddg_tools.unwrap(*arrays, lambda_=0.25, positive_cones_only=True)
        self.assertEqual(result.report["lambda"], 0.25)
        self.assertIs(result.report["positive_cones_only"], True)
        defaults = ddg_tools.unwrap(*arrays).report
        self.assertEqual(defaults["lambda"], 0.1)
        self.assertIs(defaults["positive_cones_only"], True)

    def test_seams(self):
        arrays = as_arrays(self.mesh)
        whole = ddg_tools.unwrap(*arrays)
        self.assertEqual(whole.report["n_charts"], 1)
        # The four edges of the top face (vertices 4, 5, 6, 7), in any orientation, repeated.
        seams = array.array('i', [4, 5, 7, 5, 7, 6, 6, 4, 5, 4])
        cut = ddg_tools.unwrap(*arrays, seam_edges=seams)
        self.assertEqual(cut.report["n_seam_edges"], 4)
        self.assertEqual(cut.report["n_charts"], 2)
        self.assertEqual(list(cut.polygon_charts), [0, 1, 0, 0, 0, 0])
        pairs = set(zip(cut.cut_edges[0::2], cut.cut_edges[1::2]))
        self.assertTrue({(4, 5), (5, 7), (6, 7), (4, 6)} <= pairs, pairs)
        self.assertEqual(ddg_tools.unwrap(*arrays, seam_edges=None).uvs, whole.uvs)
        self.assertEqual(ddg_tools.unwrap(*arrays, seam_edges=array.array('i')).uvs, whole.uvs)
        with self.assertRaisesRegex(ddg_tools.InputError, "seam_edges: edge 1 joins vertex 3"):
            ddg_tools.unwrap(*arrays, seam_edges=array.array('i', [0, 1, 3, 3]))
        with self.assertRaisesRegex(ddg_tools.InputError, "seam_edges: 3 items is not a multiple"):
            ddg_tools.unwrap(*arrays, seam_edges=array.array('i', [0, 1, 3]))
        with self.assertRaisesRegex(ddg_tools.InputError, r"seam_edges\[1\] = 8 is out of range"):
            ddg_tools.unwrap(*arrays, seam_edges=array.array('i', [0, 8]))

    def test_steering(self):
        arrays = as_arrays(self.mesh)
        result = ddg_tools.unwrap(
            *arrays,
            pinned_vertices=array.array('i', [6, 2]),
            prescribed_vertices=array.array('i', [4]),
            prescribed_angles=array.array('d', [1.25]),
            forbidden_vertices=array.array('i', [0, 1]),
        )
        report = result.report
        self.assertEqual(
            (report["n_pinned"], report["n_prescribed"], report["n_forbidden"]), (2, 1, 2))
        # Every prescribed cone is either imposed or reported as not applicable (it lies on a
        # boundary or a seam); the pipeline does not drop it silently.
        self.assertEqual(
            report["n_prescribed_applied"] + report["n_prescribed_not_applied"], 1)
        self.assertLessEqual(report["n_pinned_applied"], 2)
        self.assertLessEqual(report["n_forbidden_applied"], 2)
        with self.assertRaisesRegex(ddg_tools.InputError, "vertex 5 is both pinned and forbidden"):
            ddg_tools.unwrap(*arrays, pinned_vertices=array.array('i', [1, 5]),
                             forbidden_vertices=array.array('i', [5]))
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    "prescribed_vertices has 1 entries but prescribed_angles "
                                    "has 0"):
            ddg_tools.unwrap(*arrays, prescribed_vertices=array.array('i', [1]))
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    r"pinned_vertices\[0\] = 8 is out of range"):
            ddg_tools.unwrap(*arrays, pinned_vertices=array.array('i', [8]))
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    "prescribed_angles: unsupported buffer element"):
            ddg_tools.unwrap(*arrays, prescribed_vertices=array.array('i', [1]),
                             prescribed_angles=array.array('i', [1]))


class MeshErrorTest(unittest.TestCase):
    def test_non_manifold_edge(self):
        with self.assertRaises(ddg_tools.MeshError) as caught:
            ddg_tools.unwrap(*as_arrays(three_fins()))
        self.assertEqual(str(caught.exception), "edge (0, 1) is shared by 3 faces: 0, 1, 2")
        self.assertIsInstance(caught.exception, ddg_tools.DdgError)
        self.assertNotIsInstance(caught.exception, ddg_tools.InputError)

    def test_flipped_face(self):
        # Two triangles of a square; the second is wound the other way round.
        # The pipeline repairs the winding per connected piece and says so; it is not an error.
        positions = [0.0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0]
        result = ddg_tools.unwrap(array.array('f', positions),
                                  array.array('i', [0, 1, 2, 3, 2, 0]),
                                  array.array('i', [0, 3]), array.array('i', range(6)))
        self.assertEqual(result.report["n_flipped_faces"], 1)
        self.assertTrue(result.report["warnings"])
        self.assertEqual(len(result.uvs), 12)
        self.assertTrue(all(math.isfinite(x) for x in result.uvs))

    def test_bow_tie_vertex(self):
        positions = [-1.0, -1, 0, 0, 0, 0, -1, 1, 0, 1, -1, 0, 1, 1, 0]
        with self.assertRaisesRegex(ddg_tools.MeshError, "vertex 1 is non-manifold: its incident "
                                                         "faces form 2 separate fans"):
            ddg_tools.unwrap(array.array('f', positions), array.array('i', [0, 1, 2, 1, 3, 4]),
                             array.array('i', [0, 3]), array.array('i', range(6)))


class PackTest(unittest.TestCase):
    def test_single_chart(self):
        arrays = as_arrays(cube())
        result = ddg_tools.unwrap(*arrays)
        uvs, info = ddg_tools.pack_charts(result.uvs, arrays[2], result.polygon_charts, margin=0.125)
        self.assertIs(type(uvs), array.array)
        self.assertEqual((uvs.typecode, len(uvs)), ('f', 48))
        self.assertEqual((min(uvs), max(uvs)), (0.125, 0.875))
        self.assertIs(type(info), dict)
        self.assertEqual(info["n_charts"], 1)
        self.assertEqual(info["margin"], 0.125)
        # The cube unfolds without distortion into a net four unit squares long.
        self.assertEqual(info["scale"], 0.75 / 4)
        # The input is not modified.
        self.assertNotEqual(uvs, result.uvs)

    def test_charts_are_separated(self):
        arrays = as_arrays(cube())
        seams = array.array('i', [4, 5, 5, 7, 7, 6, 6, 4])
        result = ddg_tools.unwrap(*arrays, seam_edges=seams)
        margin = 0.02
        uvs, info = ddg_tools.pack_charts(result.uvs, arrays[2], result.polygon_charts,
                                          margin=margin)
        self.assertEqual(info["n_charts"], 2)
        boxes = {}
        for polygon, chart in enumerate(result.polygon_charts):
            box = boxes.setdefault(chart, [math.inf, math.inf, -math.inf, -math.inf])
            for corner in range(4 * polygon, 4 * polygon + 4):
                u, v = uvs[2 * corner], uvs[2 * corner + 1]
                box[0], box[1] = min(box[0], u), min(box[1], v)
                box[2], box[3] = max(box[2], u), max(box[3], v)
        a, b = boxes[0], boxes[1]
        tolerance = 1e-6
        for box in (a, b):
            self.assertGreaterEqual(min(box[:2]), margin - tolerance)
            self.assertLessEqual(max(box[2:]), 1 - margin + tolerance)
        gap = max(a[0] - b[2], b[0] - a[2], a[1] - b[3], b[1] - a[3])
        self.assertGreaterEqual(gap, 2 * margin - tolerance)

    def test_errors(self):
        arrays = as_arrays(cube())
        result = ddg_tools.unwrap(*arrays)
        with self.assertRaisesRegex(ddg_tools.InputError, "margin = 0.5: the margin must be"):
            ddg_tools.pack_charts(result.uvs, arrays[2], result.polygon_charts, margin=0.5)
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    "polygon_charts: got 5 rows, expected 6"):
            ddg_tools.pack_charts(result.uvs, arrays[2], result.polygon_charts[:5])
        with self.assertRaisesRegex(ddg_tools.InputError, "uvs: 47 items is not a multiple of 2"):
            ddg_tools.pack_charts(result.uvs[:47], arrays[2], result.polygon_charts)
        with self.assertRaisesRegex(ddg_tools.InputError, "polygon 5 has 2 corners"):
            ddg_tools.pack_charts(result.uvs[:44], arrays[2], result.polygon_charts)
        with self.assertRaisesRegex(TypeError, "uvs: expected an object that supports"):
            ddg_tools.pack_charts(list(result.uvs), arrays[2], result.polygon_charts)


class LargerMeshTest(unittest.TestCase):
    def test_grid(self):
        n = 200
        positions = array.array('f')
        for j in range(n + 1):
            for i in range(n + 1):
                positions.extend((i / n, j / n, 0.1 * math.sin(6.0 * i / n) * math.cos(5.0 * j / n)))
        corner_vertices = array.array('i')
        for j in range(n):
            for i in range(n):
                v = j * (n + 1) + i
                corner_vertices.extend((v, v + 1, v + n + 2, v + n + 1))
        polygon_starts = array.array('i', range(0, 4 * n * n, 4))
        triangles = array.array('i')
        for start in polygon_starts:
            triangles.extend((start, start + 1, start + 2, start, start + 2, start + 3))
        result = ddg_tools.unwrap(positions, corner_vertices, polygon_starts, triangles)
        self.assertEqual(len(result.uvs), 8 * n * n)
        self.assertEqual(result.report["n_charts"], 1)
        self.assertEqual(result.report["n_triangles"], 2 * n * n)
        self.assertEqual(result.report["measured_triangles"], 2 * n * n)
        self.assertEqual(result.report["n_flipped_uv_triangles"], 0)
        uvs, _ = ddg_tools.pack_charts(result.uvs, polygon_starts, result.polygon_charts)
        self.assertGreaterEqual(min(uvs), 0.0)
        self.assertLessEqual(max(uvs), 1.0)


@unittest.skipIf(numpy is None, "numpy is not installed")
class NumpyTest(unittest.TestCase):
    def setUp(self):
        positions, corner_vertices, polygon_starts, triangles = cube()
        self.lists = (positions, corner_vertices, polygon_starts, triangles)
        self.reference = ddg_tools.unwrap(*as_arrays(self.lists)).uvs

    def arrays(self, real=None, index=None):
        positions, corner_vertices, polygon_starts, triangles = self.lists
        return (
            numpy.array(positions, dtype=real or numpy.float32).reshape(-1, 3),
            numpy.array(corner_vertices, dtype=index or numpy.int32),
            numpy.array(polygon_starts, dtype=index or numpy.int32),
            numpy.array(triangles, dtype=index or numpy.int32).reshape(-1, 3),
        )

    def test_dtypes(self):
        for real in (numpy.float32, numpy.float64):
            for index in (numpy.int32, numpy.uint32, numpy.int64, numpy.uint64, numpy.intp,
                          numpy.int16, numpy.uint16):
                with self.subTest(real=real.__name__, index=index.__name__):
                    result = ddg_tools.unwrap(*self.arrays(real, index))
                    self.assertEqual(result.uvs, self.reference)

    def test_default_dtypes_and_flat_arrays(self):
        positions, corner_vertices, polygon_starts, triangles = self.lists
        # numpy's defaults: float64 and the platform integer.
        result = ddg_tools.unwrap(numpy.array(positions), numpy.array(corner_vertices),
                                  numpy.array(polygon_starts), numpy.array(triangles))
        self.assertEqual(result.uvs, self.reference)

    def test_results_convert_to_numpy_without_copying_semantics_surprises(self):
        result = ddg_tools.unwrap(*self.arrays())
        uvs = numpy.frombuffer(result.uvs, dtype=numpy.float32).reshape(-1, 2)
        self.assertEqual(uvs.shape, (24, 2))
        self.assertTrue(numpy.isfinite(uvs).all())
        self.assertEqual(numpy.asarray(result.polygon_charts).dtype, numpy.int32)
        self.assertEqual(numpy.asarray(result.cone_angles).dtype, numpy.float64)
        packed, _ = ddg_tools.pack_charts(uvs, self.arrays()[2], numpy.asarray(result.polygon_charts))
        self.assertEqual(len(packed), 48)

    def test_wrong_shapes(self):
        positions, corner_vertices, polygon_starts, triangles = self.arrays()
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    r"positions: buffer has shape \(12, 2\)"):
            ddg_tools.unwrap(positions.reshape(12, 2), corner_vertices, polygon_starts, triangles)
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    r"triangles: buffer has shape \(9, 4\)"):
            ddg_tools.unwrap(positions, corner_vertices, polygon_starts, triangles.reshape(9, 4))
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    r"corner_vertices: buffer has shape \(6, 4\); expected a "
                                    r"one-dimensional buffer"):
            ddg_tools.unwrap(positions, corner_vertices.reshape(6, 4), polygon_starts, triangles)
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    r"positions: buffer has shape \(2, 4, 3\)"):
            ddg_tools.unwrap(positions.reshape(2, 4, 3), corner_vertices, polygon_starts,
                             triangles)

    def test_non_contiguous_and_foreign_layouts(self):
        positions, corner_vertices, polygon_starts, triangles = self.arrays()
        doubled = numpy.repeat(positions, 2, axis=0)
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    "positions: the buffer is not C-contiguous"):
            ddg_tools.unwrap(doubled[::2], corner_vertices, polygon_starts, triangles)
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    "positions: the buffer is not C-contiguous"):
            ddg_tools.unwrap(numpy.asfortranarray(positions), corner_vertices, polygon_starts,
                             triangles)
        # The advice in the message works.
        result = ddg_tools.unwrap(numpy.ascontiguousarray(doubled[::2]), corner_vertices,
                                  polygon_starts, triangles)
        self.assertEqual(result.uvs, self.reference)
        swapped = positions.astype(positions.dtype.newbyteorder())
        with self.assertRaisesRegex(ddg_tools.InputError, "positions: unsupported buffer element"):
            ddg_tools.unwrap(swapped, corner_vertices, polygon_starts, triangles)
        with self.assertRaisesRegex(ddg_tools.InputError, "positions: unsupported buffer element"):
            ddg_tools.unwrap(positions.astype(numpy.float16), corner_vertices, polygon_starts,
                             triangles)
        with self.assertRaisesRegex(ddg_tools.InputError,
                                    "corner_vertices: unsupported buffer element"):
            ddg_tools.unwrap(positions, corner_vertices.astype(bool), polygon_starts, triangles)

    def test_uint8_indices_are_never_read_as_values(self):
        # One-byte buffers always mean raw bytes (so that `bytes` works). A uint8 index array
        # is therefore not accepted as indices; here its length gives it away.
        positions, corner_vertices, polygon_starts, triangles = self.arrays()
        with self.assertRaises(ddg_tools.InputError):
            ddg_tools.unwrap(positions, corner_vertices.astype(numpy.uint8)[:23],
                             polygon_starts, triangles)


if __name__ == "__main__":
    print("python", sys.version.split()[0], "| numpy",
          numpy.__version__ if numpy is not None else "absent",
          "| ddg_tools", ddg_tools.__version__, "from", ddg_tools.__file__)
    unittest.main(verbosity=2)
