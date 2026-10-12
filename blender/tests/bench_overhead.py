# SPDX-License-Identifier: MIT
"""Cost of the Blender side of the unwrap on meshes of about 300k triangles. Runs inside Blender:

    blender -b --factory-startup --python blender/tests/bench_overhead.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [--repeat 7]

`scripts/blender-test.sh --bench` runs it after the tests. Two meshes are measured: a grid of
quads (301 088 triangles after tessellation) and an icosphere (327 680 triangles). For each, the
steps of `core.unwrap_object` are timed separately, best and median of `--repeat` runs:

* extract: `foreach_get` of positions, face corners, face starts and loop triangles;
* seams / vertex groups: the optional inputs;
* native: `ddg_tools.unwrap`, with the binding's own share (decoding and validating the
  buffers, measuring distortion, building the result arrays) split from the pipeline's;
* pack: `ddg_tools.pack_charts`;
* write: `foreach_set` of the UVs and `mesh.update()`;
* mark seams: mapping cut edges back to edge indices;

and the whole operator from Object mode and from Edit mode (which adds leaving and re-entering
Edit mode).
"""

import argparse
import array
import importlib
import statistics
import sys
import time

import bpy


def timed(function, repeat):
    """(best, median) wall time of `function()` in seconds, and its last result."""
    times = []
    result = None
    for _ in range(repeat):
        start = time.perf_counter()
        result = function()
        times.append(time.perf_counter() - start)
    return min(times), statistics.median(times), result


def ms(seconds):
    return "{:9.2f}".format(1000.0 * seconds)


def row(label, best, median, note=""):
    print("  {:<44s}{} ms {} ms  {}".format(label, ms(best), ms(median), note).rstrip())


def bench(name, obj, core, ddg, repeat):
    mesh = obj.data
    print("{}: {} vertices, {} edges, {} faces, {} corners".format(
        name, len(mesh.vertices), len(mesh.edges), len(mesh.polygons), len(mesh.loops)))
    print("  {:<44s}{:>12s} {:>12s}".format("step", "best", "median"))

    start = time.perf_counter()
    buffers = core.mesh_buffers(mesh)
    first = time.perf_counter() - start
    n_triangles = len(buffers[3]) // 3
    row("extract, first call (Blender tessellates)", first, first,
        "{} triangles".format(n_triangles))
    best, median, buffers = timed(lambda: core.mesh_buffers(mesh), repeat)
    row("extract, tessellation cached", best, median)
    extract = best
    positions, corner_vertices, polygon_starts, triangles = buffers

    few = max(2, repeat // 3)
    signed = array.array('i', bytes(4 * len(mesh.loops)))
    best, median, _ = timed(lambda: mesh.loops.foreach_get("vertex_index", signed), few)
    row("  for comparison: loops.vertex_index -> 'i'", best, median, "(element API, signed)")
    signed = array.array('i', bytes(12 * n_triangles))
    best, median, _ = timed(lambda: mesh.loop_triangles.foreach_get("loops", signed), few)
    row("  for comparison: loop_triangles.loops -> 'i'", best, median, "(unsigned property)")
    unsigned = array.array('I', bytes(12 * n_triangles))
    best, median, _ = timed(lambda: mesh.loop_triangles.foreach_get("loops", unsigned), few)
    row("  for comparison: loop_triangles.loops -> 'I'", best, median)
    flags_view = memoryview(bytearray(len(mesh.edges))).cast('?')
    best, median, _ = timed(lambda: mesh.edges.foreach_get("use_seam", flags_view), few)
    row("  for comparison: edges.use_seam -> bool", best, median, "(element API)")

    best, median, seams = timed(lambda: core.seam_edges(mesh), repeat)
    row("seams: none marked", best, median)
    # Mark every 200th edge.
    flags = bytearray(len(mesh.edges))
    flags[::200] = b"\x01" * len(flags[::200])
    mesh.edges.foreach_set("use_seam", memoryview(flags).cast('?'))
    best, median, seams = timed(lambda: core.seam_edges(mesh), repeat)
    row("seams: {} marked".format(len(seams) // 2), best, median)
    seams_cost = best

    group = obj.vertex_groups.new(name="pin")
    group.add(list(range(0, len(mesh.vertices), 10)), 1.0, 'REPLACE')
    other = obj.vertex_groups.new(name="no cone")
    other.add(list(range(5, len(mesh.vertices), 10)), 1.0, 'REPLACE')
    best, median, pinned = timed(lambda: core.group_vertices(obj, "pin"), repeat)
    row("vertex group: {} of {} vertices".format(len(pinned), len(mesh.vertices)), best, median,
        "(Python loop over all vertices)")
    best, median, _ = timed(lambda: core.groups_vertices(obj, ["pin", "no cone"]), repeat)
    row("two vertex groups in one pass", best, median)
    group_cost = best

    def native():
        return ddg.unwrap(positions, corner_vertices, polygon_starts, triangles,
                          seam_edges=seams)

    best, median, result = timed(native, repeat)
    report = result.report
    row("native: ddg_tools.unwrap", best, median, "pipeline: " + report["pipeline"])
    binding = sum(report[k] for k in ("seconds_decode", "seconds_validate", "seconds_encode"))
    print("    of the last call: decode {} ms, validate {} ms, pipeline {} ms, measure {} ms, "
          "encode {} ms".format(*(ms(report[k]).strip() for k in (
              "seconds_decode", "seconds_validate", "seconds_pipeline", "seconds_measure",
              "seconds_encode"))))

    best, median, packed = timed(
        lambda: ddg.pack_charts(result.uvs, polygon_starts, result.polygon_charts,
                                margin=core.PACK_MARGIN), repeat)
    row("pack: ddg_tools.pack_charts", best, median,
        "{} charts".format(packed[1]["n_charts"]))
    pack = best
    uvs = packed[0]

    layer = mesh.uv_layers.active or mesh.uv_layers.new(name="UVMap", do_init=False)

    def write():
        mesh.uv_layers.active.uv.foreach_set("vector", uvs)
        mesh.update()

    best, median, _ = timed(write, repeat)
    row("write: foreach_set + mesh.update()", best, median)
    write_cost = best
    as_double = array.array('d', uvs)
    best, median, _ = timed(
        lambda: mesh.uv_layers.active.uv.foreach_set("vector", as_double), max(2, repeat // 3))
    row("  for comparison: foreach_set of array('d')", best, median,
        "(type mismatch, Blender converts item by item)")

    # Clearing the seams between runs is not part of what is measured.
    clear = memoryview(bytearray(len(mesh.edges))).cast('?')
    times = []
    for _ in range(repeat):
        mesh.attributes["uv_seam"].data.foreach_set("value", clear)
        start = time.perf_counter()
        marked = core.mark_seams(mesh, result.cut_edges)
        times.append(time.perf_counter() - start)
    assert marked == (len(result.cut_edges) // 2, 0), marked
    best, median = min(times), statistics.median(times)
    row("mark seams: {} cut edges".format(len(result.cut_edges) // 2), best, median,
        "(Python loop over all edges)")
    mark_cost = best
    mesh.edges.foreach_set("use_seam", memoryview(flags).cast('?'))

    print("  {:<44s}{} ms".format("Blender side, default options (extract+write)",
                                  ms(extract + write_cost)))
    print("  {:<44s}{} ms".format("Blender side, everything on",
                                  ms(extract + seams_cost + group_cost + write_cost
                                     + mark_cost)))
    print("  {:<44s}{} ms".format("binding (decode+validate+encode), last call", ms(binding)))
    print("  {:<44s}{} ms".format("packing", ms(pack)))

    bpy.context.view_layer.objects.active = obj
    best, median, _ = timed(lambda: bpy.ops.uv.ddg_conformal_unwrap(use_seams=False), repeat)
    row("operator, Object mode, default options", best, median)
    best, median, _ = timed(
        lambda: bpy.ops.uv.ddg_conformal_unwrap(
            pinned_group="pin", forbidden_group="no cone", mark_seams=True), max(2, repeat // 2))
    row("operator, Object mode, everything on", best, median)
    bpy.ops.object.mode_set(mode='EDIT')
    best, median, _ = timed(lambda: bpy.ops.uv.ddg_conformal_unwrap(use_seams=False),
                            max(2, repeat // 2))
    row("operator, Edit mode, default options", best, median)

    def toggle():
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')

    best, median, _ = timed(toggle, max(2, repeat // 2))
    row("  leaving and re-entering Edit mode alone", best, median, "(Blender, not the add-on)")
    bpy.ops.object.mode_set(mode='OBJECT')
    print()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", required=True)
    parser.add_argument("--repeat", type=int, default=7)
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    args = parser.parse_args(argv)

    assert bpy.ops.preferences.addon_enable(module=args.module) == {'FINISHED'}
    core = importlib.import_module(args.module + ".core")
    ddg = core.native()
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)

    print("BENCH BEGIN")
    try:
        import numpy
        numpy_version = numpy.__version__
    except ImportError:
        numpy_version = "absent"
    print("Blender {}, Python {}, numpy {}, ddg_tools {} ({})".format(
        bpy.app.version_string, sys.version.split()[0], numpy_version, ddg.__version__,
        ddg.PIPELINE))
    print("times are wall clock, best and median of {} runs\n".format(args.repeat))

    bpy.ops.mesh.primitive_grid_add(x_subdivisions=388, y_subdivisions=388, size=2.0)
    grid = bpy.context.view_layer.objects.active
    # Not flat, so that the stand-in's projection has something to measure.
    positions = array.array('f', bytes(12 * len(grid.data.vertices)))
    grid.data.vertices.foreach_get("co", positions)
    import math
    for i in range(0, len(positions), 3):
        positions[i + 2] = 0.2 * math.sin(3.0 * positions[i]) * math.cos(2.0 * positions[i + 1])
    grid.data.vertices.foreach_set("co", positions)
    grid.data.update()
    bench("quad grid 388 x 388", grid, core, ddg, args.repeat)
    bpy.data.objects.remove(grid, do_unlink=True)

    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=8)
    sphere = bpy.context.view_layer.objects.active
    bench("icosphere, 8 subdivisions", sphere, core, ddg, args.repeat)
    print("BENCH END")


main()
