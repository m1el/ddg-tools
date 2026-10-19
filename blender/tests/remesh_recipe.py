# SPDX-License-Identifier: MIT
"""Blender's own way to a coarse, uniform input for Make Developable (blender/README.md,
"Starting mesh"): Remesh modifier (Voxel) to about PRE x the target, then Decimate (Collapse) to
the target, both applied. Writes the coarse mesh as OBJ and, with --flow, runs Make Developable
with its defaults on it and writes the result. Measured with crates/ddg-developability/examples/quality.py.

    blender -b --factory-startup --python blender/tests/remesh_recipe.py -- \
        --obj DENSE.obj --target 1000 --pre 4 --out coarse.obj [--flow RESULT.obj --module bl_ext.user_default.ddg_uv_unwrap]
"""
import argparse
import importlib
import math
import sys

import bpy
import bmesh


def counts(obj):
    return len(obj.data.vertices), len(obj.evaluated_get(bpy.context.evaluated_depsgraph_get()).data.polygons)


def triangles(obj):
    obj.data.calc_loop_triangles()
    return len(obj.data.loop_triangles)


def apply(obj, name):
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_apply(modifier=name)


def export(obj, path):
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.ops.wm.obj_export(filepath=path, export_selected_objects=True, export_materials=False,
                          export_uv=False, export_normals=False, export_triangulated_mesh=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--obj", required=True)
    p.add_argument("--target", type=int, default=1000)
    p.add_argument("--pre", type=float, default=4.0, help="remesh to this many times the target, then decimate")
    p.add_argument("--out", required=True)
    p.add_argument("--flow")
    p.add_argument("--module", default="bl_ext.user_default.ddg_uv_unwrap")
    p.add_argument("--subdivisions", type=int, default=2)
    a = p.parse_args(sys.argv[sys.argv.index("--") + 1:])
    bpy.ops.wm.obj_import(filepath=a.obj)
    obj = bpy.context.selected_objects[0]
    bpy.context.view_layer.objects.active = obj
    me = obj.data
    bm = bmesh.new()
    bm.from_mesh(me)
    area = sum(f.calc_area() for f in bm.faces)
    bm.free()
    # voxel size from the area: a voxel-remeshed quad has about the area of a voxel face.
    size = math.sqrt(area / (a.pre * a.target / 2.0))
    mod = obj.modifiers.new("Remesh", 'REMESH')
    mod.mode = 'VOXEL'
    mod.voxel_size = size
    mod.use_smooth_shade = False
    apply(obj, mod.name)
    n_remesh = triangles(obj)
    # a second pass of the guess: the count scales with 1/size^2
    ratio = min(1.0, a.target / max(n_remesh, 1))
    dec = obj.modifiers.new("Decimate", 'DECIMATE')
    dec.decimate_type = 'COLLAPSE'
    dec.ratio = ratio
    dec.use_collapse_triangulate = True
    apply(obj, dec.name)
    print("RECIPE voxel_size=%.5g remesh_tris=%d decimate_ratio=%.4f tris=%d" % (size, n_remesh, ratio, triangles(obj)))
    export(obj, a.out)
    if a.flow:
        assert bpy.ops.preferences.addon_enable(module=a.module) == {'FINISHED'}
        d = importlib.import_module(a.module + ".developable_operator")
        s = d.make_developable_object(obj, subdivisions=a.subdivisions)
        print("FLOW", d.summary(s))
        export(obj, a.flow)


main()
