# SPDX-License-Identifier: MIT
"""Cost of the Blender side of the variational cut against the library time. Runs inside Blender:

    blender -b --factory-startup --python blender/tests/bench_cut.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [--segments 256 --rings 128 --steps 10]

`scripts/blender-test.sh --bench` runs it. On a UV sphere, for Start, Flow (`--steps` steps)
and Apply: the wall time of the operator, the time spent in the native module (the session
build and the solve, which is what the library costs), and the rest (the add-on: reading the
mesh, writing the attributes, the preview object, bmesh).
"""

import argparse
import importlib
import sys
import time

import bpy

parser = argparse.ArgumentParser()
parser.add_argument("--module", required=True)
parser.add_argument("--segments", type=int, default=256)
parser.add_argument("--rings", type=int, default=128)
parser.add_argument("--steps", type=int, default=10)
ARGS = parser.parse_args(sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else [])
if ARGS.module not in bpy.context.preferences.addons.keys():
    assert bpy.ops.preferences.addon_enable(module=ARGS.module) == {'FINISHED'}
cut = importlib.import_module(ARGS.module + ".cut")
import ddg_tools

bpy.ops.mesh.primitive_uv_sphere_add(segments=ARGS.segments, ring_count=ARGS.rings)
obj = bpy.context.object
mesh = obj.data
print("BENCH BEGIN")
print("UV sphere: {} vertices, {} faces, {} corners".format(
    len(mesh.vertices), len(mesh.polygons), len(mesh.loops)))


def timed(label, function):
    start = time.perf_counter()
    result = function()
    return time.perf_counter() - start, result


# Library time of each call, measured by wrapping the native session methods.
native = {"seconds": 0.0}
original = ddg_tools.CutSession


class Timed:
    def __init__(self, *a, **k):
        start = time.perf_counter()
        self._s = original(*a, **k)
        native["seconds"] += time.perf_counter() - start

    def __getattr__(self, name):
        attr = getattr(self._s, name)
        if not callable(attr):
            return attr

        def call(*a, **k):
            start = time.perf_counter()
            try:
                return attr(*a, **k)
            finally:
                native["seconds"] += time.perf_counter() - start
        return call


ddg_tools.CutSession = Timed
results = {}
try:
    for choice in ("CALLER", "INTRINSIC", "REFINED"):
        for old in list(bpy.data.objects):
            bpy.data.objects.remove(old, do_unlink=True)
        for old in list(bpy.data.meshes):
            bpy.data.meshes.remove(old)
        bpy.ops.mesh.primitive_uv_sphere_add(segments=ARGS.segments, ring_count=ARGS.rings)
        obj = bpy.context.object
        mesh = obj.data
        rows = []
        native["seconds"] = 0.0
        t, _ = timed("start", lambda: bpy.ops.uv.ddg_cut_start(
            mode='CIRCLE', radius=0.5, vertex=0, triangulation=choice))
        rows.append(("Start", t, native["seconds"]))
        native["seconds"] = 0.0
        t, _ = timed("flow", lambda: bpy.ops.uv.ddg_cut_flow(
            steps=ARGS.steps, use_stop_rule=False, length_mode='EXPLICIT', length_weight=0.1))
        rows.append(("Flow ({} steps)".format(ARGS.steps), t, native["seconds"]))
        native["seconds"] = 0.0
        t, _ = timed("apply", lambda: bpy.ops.uv.ddg_cut_apply())
        rows.append(("Apply", t, native["seconds"]))
        results[choice] = rows
finally:
    ddg_tools.CutSession = original
print("  {:<12s}{:<18s}{:>12s}{:>14s}{:>14s}".format("triangulation", "operator", "total ms",
                                                      "library ms", "add-on ms"))
for choice, rows in results.items():
    for label, total, lib in rows:
        print("  {:<12s}{:<18s}{:12.1f}{:14.1f}{:14.1f}".format(
            choice.lower(), label, 1e3 * total, 1e3 * lib, 1e3 * (total - lib)))
print("BENCH END")
