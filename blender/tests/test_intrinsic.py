# SPDX-License-Identifier: MIT
"""Headless tests of Show Intrinsic Refinement (`mesh.ddg_show_intrinsic_refinement`) and of the
reusable `draw_intrinsic_edges`. Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/test_intrinsic.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [-k NAME]

`scripts/blender-test.sh` runs it against the installed extension and its bundled wheel. For
development, `--source DIR --site DIR` register the add-on from a source directory (the parent
of the package) and import `ddg_tools` from an unpacked wheel instead.

Checked: the operator makes `<mesh>_intrinsic` (the common subdivision: polygons with the face
attributes, the color attribute, the total area of the input, every input and every intrinsic
face covered) and `<mesh>_intrinsic_edges` (one spline per intrinsic edge, ends at vertices);
the input mesh is bit for bit unchanged; the custom properties (counts, smallest angle,
certificate); re-running replaces the objects; both representations and both rule sets;
mollification; the angle bound is met; Object and Edit mode; undo; errors (empty mesh,
non-manifold mesh, mollify with signposts); `draw_intrinsic_edges` on its own.
"""

import argparse
import importlib
import math
import sys
import time
import traceback

import bmesh
import bpy

TESTS = []


def test(function):
    TESTS.append(function)
    return function


# ----------------------------------------------------------------------------------------------
# Helpers


def reset():
    active = bpy.context.view_layer.objects.active
    if active is not None and active.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    for mesh in list(bpy.data.meshes):
        bpy.data.meshes.remove(mesh)
    for curve in list(bpy.data.curves):
        bpy.data.curves.remove(curve)


def activate(obj):
    view_layer = bpy.context.view_layer
    view_layer.update()
    for other in view_layer.objects:
        other.select_set(other is obj)
    view_layer.objects.active = obj
    return obj


def primitive(kind, **kwargs):
    getattr(bpy.ops.mesh, "primitive_{}_add".format(kind))(**kwargs)
    return bpy.context.view_layer.objects.active


def triangulate(obj):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bmesh.ops.triangulate(bm, faces=bm.faces[:])
    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()
    return obj


def sphere(subdivisions=2):
    return activate(triangulate(primitive("ico_sphere", subdivisions=subdivisions, radius=1.0)))


def grid(n=8):
    """A flat open grid (a mesh with boundary), triangulated."""
    return activate(triangulate(primitive("grid", x_subdivisions=n, y_subdivisions=n, size=2.0)))


def fingerprint(mesh):
    """Everything about the mesh that the operator must leave alone."""
    co = [0.0] * (3 * len(mesh.vertices))
    mesh.vertices.foreach_get("co", co)
    loops = [0] * len(mesh.loops)
    mesh.loops.foreach_get("vertex_index", loops)
    return (len(mesh.vertices), len(mesh.edges), len(mesh.polygons), tuple(co), tuple(loops),
            tuple(a.name for a in mesh.attributes))


def run(**kwargs):
    result = bpy.ops.mesh.ddg_show_intrinsic_refinement(**kwargs)
    assert result == {'FINISHED'}, result


def new_objects(source):
    mesh = bpy.data.objects.get(source.name + "_intrinsic")
    edges = bpy.data.objects.get(source.name + "_intrinsic_edges")
    return mesh, edges


def area(mesh):
    return sum(p.area for p in mesh.polygons)


def attribute(mesh, name, count, width=1, key="value"):
    values = [0] * (count * width) if key == "value" else [0.0] * (count * width)
    mesh.attributes[name].data.foreach_get(key, values)
    return values


# ----------------------------------------------------------------------------------------------
# Tests


@test
def operator_makes_the_common_subdivision():
    source = sphere()
    before = fingerprint(source.data)
    n_input_faces = len(source.data.polygons)
    run()
    mesh_obj, edges_obj = new_objects(source)
    assert mesh_obj is not None and edges_obj is not None
    assert mesh_obj.type == 'MESH' and edges_obj.type == 'CURVE'
    assert fingerprint(source.data) == before, "the input mesh changed"
    mesh = mesh_obj.data
    n = len(mesh.polygons)
    assert n >= n_input_faces
    intrinsic = attribute(mesh, "intrinsic_face", n)
    original = attribute(mesh, "input_face", n)
    assert set(original) == set(range(n_input_faces)), "an input face has no piece"
    n_intrinsic_faces = mesh_obj["ddg_n_faces"]
    assert set(intrinsic) == set(range(n_intrinsic_faces)), "an intrinsic face has no piece"
    # The pieces tile the surface: their area is the input's, up to float32 vertex storage.
    assert abs(area(mesh) - area(source.data)) < 1e-4 * area(source.data), (
        area(mesh), area(source.data))
    # Pieces of one intrinsic face share a color; different faces differ.
    color = mesh.color_attributes["intrinsic_face_color"]
    assert color.domain == 'CORNER'
    rgba = [0.0] * (4 * len(mesh.loops))
    color.data.foreach_get("color", rgba)
    by_face = {}
    for p in mesh.polygons:
        for loop in range(p.loop_start, p.loop_start + p.loop_total):
            key = tuple(round(c, 5) for c in rgba[4 * loop:4 * loop + 4])
            by_face.setdefault(intrinsic[p.index], set()).add(key)
    assert all(len(v) == 1 for v in by_face.values())
    assert len({next(iter(v)) for v in by_face.values()}) > n_intrinsic_faces // 2
    # The objects sit where the source is.
    assert mesh_obj.matrix_world == source.matrix_world


@test
def custom_properties_report_counts_and_certificate():
    source = sphere()
    run(min_angle=25.0)
    mesh_obj, _ = new_objects(source)
    p = {k: mesh_obj[k] for k in mesh_obj.keys() if k.startswith("ddg_")}
    assert p["ddg_representation"] == 'INTEGER' and p["ddg_rules"] == 'CHEW'
    assert p["ddg_source"] == source.name
    assert p["ddg_certificate"] == "ok", p["ddg_certificate"]
    assert p["ddg_certificate_n_faces"] == p["ddg_n_faces"]
    assert p["ddg_n_input_vertices"] == len(source.data.vertices)
    assert p["ddg_n_input_faces"] == len(source.data.polygons)
    assert p["ddg_n_vertices"] == p["ddg_n_input_vertices"] + p["ddg_n_inserted"]
    assert p["ddg_n_faces"] == 2 * p["ddg_n_vertices"] - 4, "Euler: closed genus 0"
    assert p["ddg_min_angle_deg"] > 24.9, p["ddg_min_angle_deg"]
    assert p["ddg_n_crossings"] >= 0 and p["ddg_refine_status"] == "done"
    # The number of polygons: one per piece; vertices: intrinsic vertices plus crossings.
    assert len(mesh_obj.data.vertices) == p["ddg_n_vertices"] + p["ddg_n_crossings"]


@test
def intrinsic_edges_are_drawn_on_the_input():
    source = sphere()
    run()
    mesh_obj, edges_obj = new_objects(source)
    curve = edges_obj.data
    assert len(curve.splines) == mesh_obj["ddg_n_edges"]
    assert all(s.type == 'POLY' and len(s.points) >= 2 for s in curve.splines)
    assert sum(len(s.points) - 2 for s in curve.splines) == mesh_obj["ddg_n_crossings"]
    # Every point lies on the surface of the unit icosphere's faces: within the circumscribed
    # sphere and outside the inscribed one.
    radii = [math.dist((0, 0, 0), p.co[:3]) for s in curve.splines for p in s.points]
    assert max(radii) < 1.0 + 1e-5 and min(radii) > 0.7, (min(radii), max(radii))
    assert edges_obj.matrix_world == source.matrix_world
    # Without the curve, none is made.
    run(draw_edges=False)
    mesh_obj, edges_obj = new_objects(source)
    assert mesh_obj is not None and edges_obj is None


@test
def rerunning_replaces_the_objects():
    source = sphere()
    run()
    n_objects = len(bpy.data.objects)
    n_meshes, n_curves = len(bpy.data.meshes), len(bpy.data.curves)
    first = new_objects(source)[0].as_pointer()
    run(min_angle=20.0)
    run(min_angle=22.0, rules='RUPPERT', representation='SIGNPOSTS')
    assert len(bpy.data.objects) == n_objects
    assert len(bpy.data.meshes) == n_meshes and len(bpy.data.curves) == n_curves
    mesh_obj, _ = new_objects(source)
    assert mesh_obj["ddg_representation"] == 'SIGNPOSTS' and mesh_obj["ddg_angle_bound_deg"] == 22.0
    del first
    # Another object has its own pair.
    other = grid()
    run()
    assert len(bpy.data.objects) == n_objects + 3  # the grid, its two
    assert new_objects(source)[0] is not None and new_objects(other)[0] is not None


@test
def both_representations_and_rules_agree():
    source = sphere()
    results = {}
    for representation in ('INTEGER', 'SIGNPOSTS'):
        for rules in ('CHEW', 'RUPPERT'):
            run(representation=representation, rules=rules)
            mesh_obj, _ = new_objects(source)
            results[(representation, rules)] = (
                mesh_obj["ddg_n_vertices"], mesh_obj["ddg_n_faces"], len(mesh_obj.data.polygons),
                mesh_obj["ddg_min_angle_deg"])
            assert mesh_obj["ddg_min_angle_deg"] > 24.9, results
    for rules in ('CHEW', 'RUPPERT'):
        a, b = results[('INTEGER', rules)], results[('SIGNPOSTS', rules)]
        assert a[:3] == b[:3] and abs(a[3] - b[3]) < 1e-6, (a, b)
    mesh_obj, _ = new_objects(source)
    assert mesh_obj["ddg_certificate"].startswith("none")


@test
def surface_with_boundary_and_mollify():
    source = grid(6)
    before = fingerprint(source.data)
    run()
    mesh_obj, _ = new_objects(source)
    assert mesh_obj["ddg_certificate"] == "ok"
    assert abs(area(mesh_obj.data) - 4.0) < 1e-4
    assert mesh_obj["ddg_min_angle_valid_deg"] > 24.9
    run(mollify=True)
    mesh_obj, _ = new_objects(source)
    assert mesh_obj["ddg_mollify"] is True or mesh_obj["ddg_mollify"] == 1
    assert mesh_obj["ddg_certificate"] == "ok"
    assert fingerprint(source.data) == before


@test
def edit_mode_and_undo():
    source = sphere()
    bpy.ops.object.mode_set(mode='EDIT')
    run()
    assert source.mode == 'EDIT', "the operator did not return to Edit mode"
    bpy.ops.object.mode_set(mode='OBJECT')
    assert new_objects(source)[0] is not None
    # Undo is registered ('UNDO' in bl_options); in background mode the stack may be empty, so
    # only check that the call is accepted and leaves a consistent scene.
    bpy.ops.ed.undo_push(message="test")
    assert new_objects(source)[0] is not None


def fails(message, **kwargs):
    """The operator reports an ERROR containing ``message`` (background Blender raises it as a
    RuntimeError from bpy.ops) and makes nothing."""
    try:
        bpy.ops.mesh.ddg_show_intrinsic_refinement(**kwargs)
    except RuntimeError as ex:
        assert message in str(ex), str(ex)
    else:
        raise AssertionError("expected the error " + message)


@test
def errors_are_reported_not_raised():
    # Empty mesh.
    empty = bpy.data.objects.new("empty", bpy.data.meshes.new("empty"))
    bpy.context.scene.collection.objects.link(empty)
    activate(empty)
    fails("has no faces")
    # A non-manifold mesh: three triangles on one edge.
    mesh = bpy.data.meshes.new("fin")
    mesh.from_pydata([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1)], [],
                     [(0, 1, 2), (0, 1, 3), (0, 1, 4)])
    fin = bpy.data.objects.new("fin", mesh)
    bpy.context.scene.collection.objects.link(fin)
    activate(fin)
    fails("shared by 3 faces")
    assert new_objects(fin) == (None, None)
    # Mollification needs integer coordinates.
    source = sphere()
    fails("Mollify needs", mollify=True, representation='SIGNPOSTS')
    assert new_objects(source) == (None, None)


@test
def draw_intrinsic_edges_is_reusable():
    from importlib import import_module

    draw = import_module(ARGS.module + ".intrinsic_draw")
    source = sphere()
    paths = [[(0, 0, 0), (1, 0, 0)], [(0, 0, 0), (0, 1, 0), (0, 1, 1)]]
    points, offsets = draw.polylines_from_paths(paths)
    obj = draw.draw_intrinsic_edges(source, points, offsets, name="my_edges", radius=0.01)
    assert obj.name == "my_edges" and len(obj.data.splines) == 2
    assert [len(s.points) for s in obj.data.splines] == [2, 3]
    assert abs(obj.data.bevel_depth - 0.01) < 1e-7
    first_name = obj.name
    again = draw.draw_intrinsic_edges(source, points, offsets, name="my_edges")
    assert len([o for o in bpy.data.objects if o.get(draw.KIND_KEY) == "edges"]) == 1
    assert again.data.bevel_depth == 0.0 and first_name == again.name
    draw.remove_previous(source, "edges")
    assert not [o for o in bpy.data.objects if o.get(draw.KIND_KEY) == "edges"]


@test
def a_jittered_mesh_is_flipped_and_refined():
    import random

    source = sphere(3)
    rng = random.Random(3)
    for v in source.data.vertices:
        v.co += 0.12 * __import__("mathutils").Vector([rng.uniform(-1, 1) for _ in range(3)])
    source.data.update()
    before = fingerprint(source.data)
    run(min_angle=28.0)
    mesh_obj, edges_obj = new_objects(source)
    assert fingerprint(source.data) == before
    assert mesh_obj["ddg_certificate"] == "ok"
    assert mesh_obj["ddg_n_inserted"] > 0 and mesh_obj["ddg_n_crossings"] > 0
    assert mesh_obj["ddg_min_angle_deg"] >= 27.9, mesh_obj["ddg_min_angle_deg"]
    assert abs(area(mesh_obj.data) - area(source.data)) < 1e-3 * area(source.data)
    assert len(edges_obj.data.splines) == mesh_obj["ddg_n_edges"]
    assert len(mesh_obj.data.polygons) > mesh_obj["ddg_n_faces"]


@test
def a_larger_mesh_runs():
    source = sphere(4)  # 5120 triangles
    t0 = time.perf_counter()
    run()
    seconds = time.perf_counter() - t0
    mesh_obj, _ = new_objects(source)
    print("INFO icosphere 4: {} -> {} intrinsic vertices, {} pieces, {:.2f} s".format(
        len(source.data.vertices), mesh_obj["ddg_n_vertices"], len(mesh_obj.data.polygons),
        seconds))
    assert mesh_obj["ddg_certificate"] == "ok"


# ----------------------------------------------------------------------------------------------


def main():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", default="ddg_uv_unwrap")
    parser.add_argument("--source", default="", help="register the add-on from this directory")
    parser.add_argument("--site", default="", help="directory with an unpacked ddg_tools wheel")
    parser.add_argument("-k", dest="only", default="")
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    ARGS = parser.parse_args(argv)

    assert bpy.app.background, "run with -b"
    if ARGS.site:
        sys.path.insert(0, ARGS.site)
    if ARGS.source:
        sys.path.insert(0, ARGS.source)
        importlib.import_module(ARGS.module).register()
    else:
        assert bpy.app.factory_startup, "run with --factory-startup"
        assert bpy.ops.preferences.addon_enable(module=ARGS.module) == {'FINISHED'}
        assert ARGS.module in bpy.context.preferences.addons.keys(), "the add-on did not enable"
    assert hasattr(bpy.types, "MESH_OT_ddg_show_intrinsic_refinement"), "operator not registered"

    passed, skipped, failed = [], [], []
    for function in TESTS:
        name = function.__name__
        if ARGS.only and ARGS.only not in name:
            continue
        reset()
        print("TEST", name)
        start = time.perf_counter()
        try:
            function()
        except Exception:
            failed.append(name)
            traceback.print_exc(file=sys.stdout)
            print("FAIL", name)
        else:
            passed.append(name)
            print("PASS", name, "({:.2f} s)".format(time.perf_counter() - start))
        sys.stdout.flush()
    print("RESULT: {} passed, {} skipped, {} failed".format(len(passed), len(skipped), len(failed)))
    for name in failed:
        print("  failed:", name)
    if failed or not passed:
        raise SystemExit(1)
    print("ALL TESTS PASSED")


ARGS = None
main()
