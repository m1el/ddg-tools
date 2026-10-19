# SPDX-License-Identifier: MIT
"""End-to-end measurement of Make Developable followed by Developable Seams + UV on one mesh
(not a test: it prints one JSON line per run). Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/e2e_developable.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap --mesh PATH.obj | --primitive NAME \
        [--subdivisions N] [--iterations N] [--sweep] --out FILE.jsonl

The mesh is imported (an OBJ) or made (a primitive: cube, sphere, torus, monkey, cylinder, plane);
the flow runs once with the operator's defaults except `--subdivisions` / `--iterations` (use
`--subdivisions -1` for no flow at all); then the seams operator runs with the default threshold,
and with `--sweep` also with other thresholds, defect tolerances and lambdas (every setting
starts from the mesh after the flow, with the seams cleared). Each line has the flow report
(rounds, iterations, status, stopped, energies, triangles, time) and per setting the seam edges,
terminals, charts, cones, overlap, mirrored UV triangles and time. `scripts/e2e-developable.sh`
runs the cases in parallel.
"""

import argparse
import array
import importlib
import json
import sys
import time

import bmesh
import bpy

PRIMITIVES = {
    "cube": lambda: bpy.ops.mesh.primitive_cube_add(size=2.0),
    "sphere": lambda: bpy.ops.mesh.primitive_uv_sphere_add(segments=32, ring_count=16),
    "torus": lambda: bpy.ops.mesh.primitive_torus_add(major_segments=32, minor_segments=16),
    "monkey": lambda: bpy.ops.mesh.primitive_monkey_add(),
    "cylinder": lambda: bpy.ops.mesh.primitive_cylinder_add(vertices=32),
    "plane": lambda: bpy.ops.mesh.primitive_grid_add(x_subdivisions=16, y_subdivisions=16, size=2.0),
}


def make_object(args):
    if args.mesh:
        before = set(bpy.data.objects)
        bpy.ops.wm.obj_import(filepath=args.mesh)
        new = [o for o in bpy.data.objects if o not in before and o.type == 'MESH']
        assert len(new) == 1, new
        obj = new[0]
    else:
        PRIMITIVES[args.primitive]()
        obj = bpy.context.view_layer.objects.active
        if args.primitive == "cube":
            bm = bmesh.new()
            bm.from_mesh(obj.data)
            bmesh.ops.subdivide_edges(bm, edges=bm.edges[:], cuts=4, use_grid_fill=True)
            bm.to_mesh(obj.data)
            bm.free()
        if args.primitive == "plane":
            import random
            rng = random.Random(1)
            for v in obj.data.vertices:
                v.co.z = 0.15 * rng.uniform(-1.0, 1.0)
        obj.data.update()
    bpy.context.view_layer.objects.active = obj
    for o in bpy.context.view_layer.objects:
        o.select_set(o == obj)
    return obj


def mirrored(mesh, name):
    mesh.calc_loop_triangles()
    uv = array.array('f', bytes(8 * len(mesh.loops)))
    mesh.uv_layers[name].data.foreach_get("uv", uv)
    n = 0
    for t in mesh.loop_triangles:
        a, b, c = (uv[2 * i:2 * i + 2] for i in t.loops)
        if (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) < 0.0:
            n += 1
    return n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", default="ddg_uv_unwrap")
    parser.add_argument("--mesh", default="")
    parser.add_argument("--primitive", default="")
    parser.add_argument("--label", default="")
    parser.add_argument("--subdivisions", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=500)
    parser.add_argument("--sweep", action="store_true")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(sys.argv[sys.argv.index("--") + 1:])
    assert bpy.app.background and bpy.app.factory_startup
    assert bpy.ops.preferences.addon_enable(module=args.module) == {'FINISHED'}
    d = importlib.import_module(args.module + ".developable_operator")
    label = args.label or args.mesh.rsplit("/", 1)[-1] or args.primitive
    base = {"case": label, "subdivisions": args.subdivisions, "iterations": args.iterations}
    out = open(args.out, "a")

    def emit(record):
        out.write(json.dumps({**base, **record}) + "\n")
        out.flush()

    obj = make_object(args)
    mesh = obj.data
    base["polygons"] = len(mesh.polygons)
    base["vertices"] = len(mesh.vertices)
    if args.subdivisions >= 0:
        t0 = time.perf_counter()
        try:
            stats = d.make_developable_object(obj, iterations=args.iterations,
                                              subdivisions=args.subdivisions)
        except Exception as ex:  # a message for the user
            emit({"stage": "flow", "error": str(ex)})
            return
        r = stats["report"]
        emit({
            "stage": "flow",
            "triangles_in": r["n_input_triangles"], "triangles_out": r["n_triangles"],
            "energy_initial": r["energy_initial"], "energy_final": r["energy_final"],
            "iterations_done": r["iterations"],
            "rounds": [{"iterations": x["iterations"], "status": x["status"],
                        "flips": x["flips"], "collapses": x["collapses"]} for x in r["rounds"]],
            "stopped": r["stopped"], "bad_faces": r["bad_faces"],
            "inadmissible_trials": r["inadmissible_trials"], "min_angle": r["min_angle"],
            "seconds": time.perf_counter() - t0, "seconds_native": stats["seconds_native"],
            "warnings": r["warnings"],
        })
    settings = [{"name": "default"}]
    if args.sweep:
        settings += [
            {"name": "eps 1e-2", "epsilon": 1e-2},
            {"name": "eps 1e-4", "epsilon": 1e-4},
            {"name": "defect 0.3", "defect_tolerance": 0.3},
            {"name": "paper rule (no defect)", "use_defect_rule": False},
            {"name": "eps 1e-2, defect 0.3", "epsilon": 1e-2, "defect_tolerance": 0.3},
            {"name": "lambda 1", "lambda_": 1.0},
        ]
    for i, setting in enumerate(settings):
        for e in mesh.edges:
            e.use_seam = False
        kwargs = {k: v for k, v in setting.items() if k != "name"}
        uv_name = "E2E"  # Blender allows 8 UV maps; the operator reuses the layer
        t0 = time.perf_counter()
        try:
            s = d.developable_seams_uv(obj, uv_name=uv_name, **kwargs)
        except Exception as ex:
            emit({"stage": "seams", "setting": setting["name"], "error": str(ex)})
            continue
        seconds = time.perf_counter() - t0
        r = s["report"]
        emit({
            "stage": "seams", "setting": setting["name"],
            "seam_edges": s["n_seams"], "terminals": s["n_terminals"],
            "charts": r.get("n_charts"), "cones": s["n_cones"],
            "overlap": r.get("overlapped_fraction"),
            "mirrored": mirrored(mesh, uv_name),
            "mirrored_pipeline": r.get("n_flipped_uv_triangles"),
            "collapsed": r.get("n_collapsed_uv_triangles"),
            "angle_error_max_deg": (r.get("measured_angle_error_max") or 0) * 57.29577951308232,
            "seconds": seconds,
        })
    emit({"stage": "done"})


main()
