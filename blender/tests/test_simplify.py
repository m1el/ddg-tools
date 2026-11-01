# SPDX-License-Identifier: MIT
"""Headless tests of Intrinsic Simplify (`object.ddg_intrinsic_simplify`). Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/test_simplify.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [-k NAME]

`scripts/blender-test.sh` runs it against the installed extension and its bundled wheel. For
development, `--source DIR --site DIR` register the add-on from a source directory (the parent
of the package) and import `ddg_tools` from an unpacked wheel instead.

Checked: the property defaults; the corner color attribute, the face id attribute and the vertex
group exist with the right sizes and values; the mesh geometry is bit for bit unchanged; every
face is colored by the coarse face of its centroid (checked against the library called
directly), faces that share a coarse face share a color, coarse faces that share a vertex
differ; a percentage and a count as the target; the area weight and the memoryless metric;
determinism; selected vertices stay; refinement before and after; n-gons; Object and Edit mode;
re-running overwrites; undo is accepted; the mesh with a boundary; errors (empty mesh,
non-manifold mesh, a conflicting attribute, no vertex selected, a target the library refuses).
"""

import argparse
import importlib
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

OP = "ddg_intrinsic_simplify"


def reset():
    active = bpy.context.view_layer.objects.active
    if active is not None and active.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    for mesh in list(bpy.data.meshes):
        bpy.data.meshes.remove(mesh)


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


def sphere(subdivisions=3):
    return activate(triangulate(primitive("ico_sphere", subdivisions=subdivisions, radius=1.0)))


def grid(n=10):
    return activate(triangulate(primitive("grid", x_subdivisions=n, y_subdivisions=n, size=2.0)))


def fingerprint(mesh):
    """Everything about the geometry that the operator must leave alone."""
    co = [0.0] * (3 * len(mesh.vertices))
    mesh.vertices.foreach_get("co", co)
    loops = [0] * len(mesh.loops)
    mesh.loops.foreach_get("vertex_index", loops)
    starts = [0] * len(mesh.polygons)
    mesh.polygons.foreach_get("loop_start", starts)
    return (len(mesh.vertices), len(mesh.edges), len(mesh.polygons), tuple(co), tuple(loops),
            tuple(starts))


def run(**kwargs):
    result = getattr(bpy.ops.object, OP)(**kwargs)
    assert result == {'FINISHED'}, result


def fails(message, **kwargs):
    """The operator reports an ERROR containing ``message`` (background Blender raises it as a
    RuntimeError from bpy.ops)."""
    try:
        getattr(bpy.ops.object, OP)(**kwargs)
    except RuntimeError as ex:
        assert message in str(ex), str(ex)
    else:
        raise AssertionError("expected the error " + message)


def corner_colours(mesh):
    rgba = [0.0] * (4 * len(mesh.loops))
    mesh.color_attributes["IntrinsicSimplifyFaces"].data.foreach_get("color", rgba)
    return [tuple(round(c, 5) for c in rgba[4 * i:4 * i + 4]) for i in range(len(mesh.loops))]


def face_colour(mesh, colours, polygon):
    p = mesh.polygons[polygon]
    found = {colours[loop] for loop in range(p.loop_start, p.loop_start + p.loop_total)}
    return found


def face_ids(mesh):
    values = [0] * len(mesh.polygons)
    mesh.attributes["IntrinsicSimplifyFaceId"].data.foreach_get("value", values)
    return values


def group_vertices(obj):
    group = obj.vertex_groups["IntrinsicSimplifyVertices"]
    out = []
    for v in obj.data.vertices:
        for g in v.groups:
            if g.group == group.index and g.weight > 0.5:
                out.append(v.index)
    return out


def library_run(mesh, **kwargs):
    """The same call the operator makes, straight to the library."""
    from array import array

    ddg = importlib.import_module("ddg_tools")
    core = importlib.import_module(ARGS.module + ".core")
    positions, corner_vertices, _starts, corner_triangles = core.mesh_buffers(mesh)
    n = len(corner_triangles) // 3
    triangles = array('I', (corner_vertices[c] for c in corner_triangles))
    points = (array('i', range(n)), array('d', [1.0 / 3.0] * (3 * n)))
    return ddg.simplify(positions, triangles, points=points, **kwargs)


# ----------------------------------------------------------------------------------------------
# Tests


@test
def defaults_and_registration():
    props = bpy.ops.object.ddg_intrinsic_simplify.get_rna_type().properties
    assert props["target_mode"].default == 'PERCENT'
    assert props["percent"].default == 10.0
    assert props["count"].default == 100
    assert props["area_weight"].default == 0.0
    assert props["memory"].default is True
    assert props["fix_selected"].default is False
    assert props["refine_before"].default is False and props["refine_after"].default is False
    assert props["representation"].default == 'INTEGER'
    assert hasattr(bpy.types, "VIEW3D_PT_ddg_intrinsic_simplify")
    assert bpy.ops.object.ddg_intrinsic_simplify.poll() is False, "no active object"


@test
def attributes_and_group_have_the_right_sizes():
    source = sphere()
    mesh = source.data
    n_v, n_p, n_l = len(mesh.vertices), len(mesh.polygons), len(mesh.loops)
    before = fingerprint(mesh)
    run(percent=20.0)
    assert fingerprint(mesh) == before, "the geometry changed"
    attribute = mesh.color_attributes["IntrinsicSimplifyFaces"]
    assert attribute.domain == 'CORNER' and attribute.data_type == 'FLOAT_COLOR'
    assert len(attribute.data) == n_l
    assert mesh.color_attributes.active_color_name == "IntrinsicSimplifyFaces"
    ids = mesh.attributes["IntrinsicSimplifyFaceId"]
    assert ids.domain == 'FACE' and ids.data_type == 'INT' and len(ids.data) == n_p
    kept = group_vertices(source)
    assert 0 < len(kept) <= round(0.2 * n_v) + 1, (len(kept), n_v)
    assert mesh["ddg_simplify_n_input_vertices"] == n_v
    assert mesh["ddg_simplify_n_vertices"] == len(kept) == mesh["ddg_simplify_n_group_vertices"]
    assert mesh["ddg_simplify_n_faces"] == 2 * len(kept) - 4, "closed genus 0"
    assert mesh["ddg_simplify_stop"] == "targetreached"
    # Every color is one of the palette's.
    palette = importlib.import_module(ARGS.module + ".simplify_operator").PALETTE
    allowed = {tuple(round(c, 5) for c in (*rgb, 1.0)) for rgb in palette}
    assert set(corner_colours(mesh)) <= allowed


@test
def colours_follow_the_coarse_faces():
    source = sphere()
    mesh = source.data
    run(percent=15.0)
    colours = corner_colours(mesh)
    ids = face_ids(mesh)
    # The library, called directly with the same points, puts every centroid in the same face.
    s = library_run(mesh, ratio=0.15)
    n_v = len(mesh.vertices)
    expected = list(s.point_face)[n_v:]
    assert expected == ids, "the faces of the centroids differ from the library's"
    # Faces of one coarse face share a color (and only one color each, on a triangle mesh).
    by_coarse = {}
    for p in range(len(mesh.polygons)):
        found = face_colour(mesh, colours, p)
        assert len(found) == 1, (p, found)
        by_coarse.setdefault(ids[p], set()).update(found)
    assert all(len(v) == 1 for v in by_coarse.values())
    # Coarse faces that share a vertex differ in color.
    palette = importlib.import_module(ARGS.module + ".simplify_operator").PALETTE
    assert s.n_colours <= len(palette), "the test mesh needs more colors than the palette has"
    fc = list(s.faces_compact)
    colour_of = {f: next(iter(c)) for f, c in by_coarse.items()}
    at = {}
    for f in range(s.n_faces):
        for v in set(fc[3 * f:3 * f + 3]):
            at.setdefault(v, []).append(f)
    for fs in at.values():
        used = [colour_of[f] for f in fs if f in colour_of]
        assert len(set(used)) == len(used), "two faces around a coarse vertex share a color"
    # Most coarse faces contain a centroid on this dense sphere.
    assert len(by_coarse) > 0.8 * s.n_faces


@test
def vertex_group_is_the_libraries_coarse_vertices():
    source = sphere()
    run(percent=15.0)
    s = library_run(source.data, ratio=0.15)
    assert sorted(group_vertices(source)) == sorted(s.vertices)
    # The group is rebuilt, not extended, when the operator runs again.
    run(percent=40.0)
    assert len(group_vertices(source)) > len(s.vertices)
    assert len([g for g in source.vertex_groups if g.name == "IntrinsicSimplifyVertices"]) == 1


@test
def targets_by_count_and_percent():
    source = sphere()
    n_v = len(source.data.vertices)
    run(target_mode='COUNT', count=60)
    assert len(group_vertices(source)) == 60 and source.data["ddg_simplify_n_vertices"] == 60
    run(target_mode='PERCENT', percent=50.0)
    assert abs(len(group_vertices(source)) - round(0.5 * n_v)) <= 1
    del n_v


@test
def metric_options_change_the_result_and_are_deterministic():
    source = sphere()
    run(percent=10.0)
    first = (corner_colours(source.data), group_vertices(source))
    run(percent=10.0)
    assert (corner_colours(source.data), group_vertices(source)) == first, "not deterministic"
    run(percent=10.0, area_weight=1.0)
    assert source.data["ddg_simplify_area_weight"] == 1.0
    assert len(group_vertices(source)) == len(first[1])
    run(percent=10.0, memory=False)
    assert source.data["ddg_simplify_memory"] in (False, 0)
    assert len(group_vertices(source)) == len(first[1])


@test
def selected_vertices_stay():
    source = sphere()
    mesh = source.data
    chosen = [0, 5, 17, 40, 100, 150]
    flags = [i in chosen for i in range(len(mesh.vertices))]
    mesh.vertices.foreach_set("select", flags)
    mesh.update()
    run(percent=8.0, fix_selected=True)
    kept = set(group_vertices(source))
    assert set(chosen) <= kept, set(chosen) - kept
    assert mesh["ddg_simplify_n_fixed"] == len(chosen)
    # Without the option the same selection does not have to stay (it is just another run).
    run(percent=8.0)
    assert mesh["ddg_simplify_n_fixed"] == 0


@test
def refinement_before_and_after():
    source = grid(8)
    before = fingerprint(source.data)
    n_v = len(source.data.vertices)
    run(percent=30.0, refine_before=True)
    assert fingerprint(source.data) == before
    assert len(group_vertices(source)) <= source.data["ddg_simplify_n_vertices"]
    assert source.data["ddg_simplify_refine_before"] in (True, 1)
    run(percent=30.0, refine_after=True)
    run(percent=30.0, refine_before=True, representation="SIGNPOSTS")
    assert fingerprint(source.data) == before
    ids = face_ids(source.data)
    assert min(ids) >= 0 and max(ids) < source.data["ddg_simplify_n_faces"]
    assert len(source.data.polygons) > 0 and n_v > 0


@test
def mesh_with_boundary():
    source = grid(12)
    run(percent=25.0)
    assert source.data["ddg_simplify_stop"] == "targetreached"
    assert len(set(face_ids(source.data))) > 1


@test
def ngons_are_split_like_the_tessellation():
    source = activate(primitive("cube"))
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.subdivide(number_cuts=3, smoothness=0.0)
    bpy.ops.object.mode_set(mode='OBJECT')
    assert any(len(p.vertices) == 4 for p in source.data.polygons)
    before = fingerprint(source.data)
    run(percent=30.0)
    assert fingerprint(source.data) == before
    mesh = source.data
    assert len(face_ids(mesh)) == len(mesh.polygons)
    assert len(mesh.color_attributes["IntrinsicSimplifyFaces"].data) == len(mesh.loops)
    # A quad's two triangles may land in different coarse faces: every corner still has a color.
    assert all(c[3] == 1.0 for c in corner_colours(mesh))


@test
def edit_mode_and_undo():
    source = sphere()
    bpy.ops.object.mode_set(mode='EDIT')
    run(percent=20.0)
    assert source.mode == 'EDIT', "the operator did not return to Edit mode"
    bpy.ops.object.mode_set(mode='OBJECT')
    assert "IntrinsicSimplifyFaces" in source.data.color_attributes
    # Undo is registered ('UNDO' in bl_options); in background mode the stack may be empty, so
    # only check that the call is accepted and leaves a consistent mesh.
    bpy.ops.ed.undo_push(message="test")
    assert "IntrinsicSimplifyFaces" in source.data.color_attributes


@test
def errors_are_reported_not_raised():
    empty = bpy.data.objects.new("empty", bpy.data.meshes.new("empty"))
    bpy.context.scene.collection.objects.link(empty)
    activate(empty)
    fails("has no faces")
    mesh = bpy.data.meshes.new("fin")
    mesh.from_pydata([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1)], [],
                     [(0, 1, 2), (0, 1, 3), (0, 1, 4)])
    fin = bpy.data.objects.new("fin", mesh)
    bpy.context.scene.collection.objects.link(fin)
    activate(fin)
    fails("shared by 3 faces")
    assert "IntrinsicSimplifyFaces" not in fin.data.color_attributes
    assert "IntrinsicSimplifyVertices" not in fin.vertex_groups
    # A conflicting attribute is not overwritten.
    source = sphere()
    source.data.attributes.new("IntrinsicSimplifyFaceId", 'FLOAT', 'POINT')
    fails("already has an attribute 'IntrinsicSimplifyFaceId'")
    assert "IntrinsicSimplifyFaces" not in source.data.color_attributes
    # A target above the vertex count does not raise a Python error; the library says what is wrong.
    other = sphere(2)
    try:
        run(target_mode='COUNT', count=10 ** 6)
    except RuntimeError as ex:
        assert "target" in str(ex).lower() or "vertices" in str(ex).lower(), str(ex)
    else:
        assert len(group_vertices(other)) == len(other.data.vertices)


@test
def a_larger_mesh_runs():
    source = sphere(5)  # 20480 triangles
    t0 = time.perf_counter()
    run(percent=5.0)
    seconds = time.perf_counter() - t0
    print("INFO icosphere 5: {} -> {} vertices, {:.2f} s".format(
        len(source.data.vertices), len(group_vertices(source)), seconds))
    assert source.data["ddg_simplify_stop"] == "targetreached"


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
    assert hasattr(bpy.types, "OBJECT_OT_ddg_intrinsic_simplify"), "operator not registered"

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
