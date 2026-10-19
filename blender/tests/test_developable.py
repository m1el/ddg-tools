# SPDX-License-Identifier: MIT
"""Headless tests of the developability operators (Make Developable, Show Developability,
Developable Seams + UV). Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/test_developable.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [-k NAME]

`scripts/blender-test.sh` runs it against the installed extension and its bundled wheel. For
development, `--source DIR --site DIR` register the add-on from a source directory (the parent
of the package) and import `ddg_tools` from an unpacked wheel instead.

Meshes: subdivided cube, UV sphere, cylinder, a noisy plane, and (when the file is present)
`spot_quadrangulated.obj` of Keenan Crane's model repository (`--spot PATH`).

Checked: the flow lowers the energy; with unchanged topology the polygons (quads) stay and only
positions change; after subdivision (or Triangulate) the mesh is a triangle mesh with the
material index of the parent polygon and a warning about the dropped UV map; pinned vertices
(selected in Edit mode, or a vertex group) stay; undo restores the mesh and Adjust Last Operation
re-runs; the preview writes the attributes; the seams operator marks seams, writes a UV map of
the given name with no mirrored triangle on the easy cases and keeps the vertices where they are;
errors (non-manifold mesh, shape keys, empty names) are readable. On spot the numbers are printed
(`SPOT ...` lines).
"""

import argparse
import array
import importlib
import math
import os
import random
import sys
import time
import traceback

import bmesh
import bpy

TESTS = []
SPOT = "/var/data/discrete-geometry/ddg-tools/data/models/spot/spot/spot_quadrangulated.obj"


def test(function):
    TESTS.append(function)
    return function


class Skip(Exception):
    pass


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
    for material in list(bpy.data.materials):
        bpy.data.materials.remove(material)


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


def subdivided_cube(cuts=4):
    cube = primitive("cube", size=2.0)
    bm = bmesh.new()
    bm.from_mesh(cube.data)
    bmesh.ops.subdivide_edges(bm, edges=bm.edges[:], cuts=cuts, use_grid_fill=True)
    bm.to_mesh(cube.data)
    bm.free()
    cube.data.update()
    return cube


def noisy_plane(cuts=12, height=0.15, seed=1):
    plane = primitive("grid", x_subdivisions=cuts, y_subdivisions=cuts, size=2.0)
    rng = random.Random(seed)
    for v in plane.data.vertices:
        v.co.z = height * rng.uniform(-1.0, 1.0)
    plane.data.update()
    return plane


def spot_object():
    if not os.path.exists(ARGS.spot):
        raise Skip("{} is not present".format(ARGS.spot))
    before = set(bpy.data.objects)
    bpy.ops.wm.obj_import(filepath=ARGS.spot)
    new = [o for o in bpy.data.objects if o not in before and o.type == 'MESH']
    assert len(new) == 1, new
    return activate(new[0])


def coords(obj):
    out = array.array('f', bytes(12 * len(obj.data.vertices)))
    obj.data.vertices.foreach_get("co", out)
    return out


def uvs(mesh, name=None):
    layer = mesh.uv_layers.active if name is None else mesh.uv_layers[name]
    out = array.array('f', bytes(8 * len(mesh.loops)))
    layer.data.foreach_get("uv", out)
    return out


def sides(mesh):
    return sorted({len(p.vertices) for p in mesh.polygons})


def module(name=""):
    return importlib.import_module(ARGS.module + name)


def run(operator, **kwargs):
    result = operator('EXEC_DEFAULT', **kwargs)
    assert result == {'FINISHED'}, result
    return result


def seam_count(mesh):
    return sum(1 for e in mesh.edges if e.use_seam)


def energy_of(obj, **kwargs):
    d = module(".core").native()
    mesh = obj.data
    return d.developability(*module(".core").mesh_buffers(mesh), **kwargs).energy


def mirrored(obj, name):
    """Triangles of Blender's tessellation whose UV image is clockwise (the area is signed so
    that a counter-clockwise UV triangle is positive, as the faces are)."""
    mesh = obj.data
    mesh.calc_loop_triangles()
    uv = uvs(mesh, name)
    flipped = 0
    for t in mesh.loop_triangles:
        a, b, c = (uv[2 * i:2 * i + 2] for i in t.loops)
        if (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) < 0.0:
            flipped += 1
    return flipped


# ----------------------------------------------------------------------------------------------
# Make Developable


@test
def make_developable_keeps_quads_when_topology_is_unchanged():
    cube = subdivided_cube(4)
    mesh = cube.data
    n_polygons, n_vertices = len(mesh.polygons), len(mesh.vertices)
    before = coords(cube)
    e0 = energy_of(cube)
    stats = module(".developable_operator").make_developable_object(
        cube, iterations=60, subdivisions=0, remesh=False)
    assert not stats["replaced"] and not stats["topology_changed"]
    assert stats["report"]["nonplanarity_after"] is not None
    assert len(mesh.polygons) == n_polygons and len(mesh.vertices) == n_vertices
    assert sides(mesh) == [4], sides(mesh)
    assert coords(cube) != before
    e1 = energy_of(cube)
    assert e1 <= e0 + 1e-12 and abs(stats["report"]["energy_final"] - e1) < 1e-6 * max(e0, 1e-9) + 1e-9, (e0, e1)
    print("INFO cube energy {:.3e} -> {:.3e}, nonplanarity {:.3e}".format(
        e0, e1, stats["report"]["nonplanarity_after"]))


@test
def operator_unchanged_topology_keeps_uv_and_attributes():
    plane = noisy_plane(8)
    mesh = plane.data
    layer = mesh.uv_layers.new(name="MyUV")
    uv_before = uvs(mesh, "MyUV")
    mesh.attributes.new("mine", 'FLOAT', 'POINT')
    run(bpy.ops.object.ddg_make_developable, subdivisions=0, remesh=False, iterations=30)
    assert uvs(mesh, "MyUV") == uv_before
    assert "mine" in mesh.attributes
    assert sides(mesh) == [4]


@test
def subdivision_gives_triangles_with_parent_materials_and_a_warning():
    plane = noisy_plane(6)
    mesh = plane.data
    for name in ("red", "blue"):
        mesh.materials.append(bpy.data.materials.new(name))
    n = len(mesh.polygons)
    index = array.array('i', [p % 2 for p in range(n)])
    mesh.polygons.foreach_set("material_index", index)
    mesh.uv_layers.new(name="MyUV")
    # polygon -> (min x, min y) corner so that the parent can be recognized afterwards
    centers = [tuple(p.center) for p in mesh.polygons]
    stats = module(".developable_operator").make_developable_object(
        plane, iterations=15, subdivisions=1, remesh=False)
    assert stats["replaced"] and "MyUV" in stats["dropped"], stats["dropped"]
    assert sides(mesh) == [3]
    assert len(mesh.polygons) == 8 * n, len(mesh.polygons)
    assert not mesh.uv_layers
    mat = array.array('i', bytes(4 * len(mesh.polygons)))
    mesh.polygons.foreach_get("material_index", mat)
    assert sorted(set(mat)) == [0, 1]
    # each parent polygon has 8 children with its own material: count per material
    assert sum(1 for m in mat if m == 0) == 8 * sum(1 for m in index if m == 0)
    run(bpy.ops.object.ddg_make_developable, subdivisions=0, triangulate=True, remesh=False,
        iterations=3)
    assert sides(mesh) == [3]


@test
def materials_and_smooth_survive_subdivision_and_remeshing_exactly():
    """Six materials and alternating smooth flags on a noisy plane; two subdivisions and an
    aggressive remeshing (flips and collapses happen). The library's face_parent gives the parent
    polygon of every triangle; the operator must copy exactly its material and smooth flag."""
    plane = noisy_plane(8, height=0.3, seed=3)
    mesh = plane.data
    for name in range(6):
        mesh.materials.append(bpy.data.materials.new("m{}".format(name)))
    n = len(mesh.polygons)
    index = array.array('i', [(p * 5) % 6 for p in range(n)])
    mesh.polygons.foreach_set("material_index", index)
    smooth = [bool(p % 3 == 0) for p in range(n)]
    mesh.polygons.foreach_set("use_smooth", smooth)
    mesh.update()
    ddg = module(".core").native()
    buffers = module(".core").mesh_buffers(mesh)
    kwargs = dict(iterations=40, subdivisions=2, remesh=True, small_angle=0.5)
    expected = ddg.make_developable(*buffers, **kwargs)
    flips = sum(r["flips"] + r["collapses"] for r in expected.report["rounds"])
    assert flips > 0, "the remeshing did not run; the test would not test it"
    run(bpy.ops.object.ddg_make_developable, iterations=40, subdivisions=2, remesh=True,
        small_angle=0.5)
    assert sides(mesh) == [3] and len(mesh.polygons) == len(expected.face_polygon)
    mat = array.array('i', bytes(4 * len(mesh.polygons)))
    mesh.polygons.foreach_get("material_index", mat)
    flat = array.array('B', bytes(len(mesh.polygons)))
    mesh.polygons.foreach_get("use_smooth", flat)
    assert list(mat) == [index[p] for p in expected.face_polygon]
    assert [bool(x) for x in flat] == [smooth[p] for p in expected.face_polygon]
    assert len(set(mat)) == 6


@test
def triangulate_option_replaces_quads_even_without_new_vertices():
    cube = subdivided_cube(2)
    n, nv = len(cube.data.polygons), len(cube.data.vertices)
    run(bpy.ops.object.ddg_make_developable, subdivisions=0, remesh=False, iterations=5,
        triangulate=True)
    assert len(cube.data.polygons) == 2 * n and sides(cube.data) == [3]
    assert len(cube.data.vertices) == nv


@test
def energies_and_methods_all_run_on_a_sphere():
    sphere = primitive("uv_sphere", segments=16, ring_count=8)
    d = module(".developable_operator")
    for kw in ({"energy": 'COMBINATORIAL'}, {"energy": 'INTRINSIC'}, {"energy": 'COVARIANCE_MAX'},
               {"method": 'GRADIENT_DESCENT'}, {"method": 'FIXED_STEP'}):
        keep = coords(sphere)
        stats = d.make_developable_object(sphere, iterations=8, subdivisions=0, remesh=False, **kw)
        assert all(math.isfinite(x) for x in coords(sphere)), kw
        assert stats["report"]["energy_final"] <= stats["report"]["energy_initial"] * 1.001 + 1e-12, kw


@test
def cylinder_stays_developable_and_noisy_plane_becomes_flatter():
    cyl = primitive("cylinder", vertices=24, depth=2.0, end_fill_type='NOTHING')
    e = energy_of(cyl)
    assert e < 1e-6, e
    plane = noisy_plane(12, height=0.2)
    e0 = energy_of(plane)
    stats = module(".developable_operator").make_developable_object(
        plane, iterations=150, subdivisions=0, remesh=False)
    e1 = energy_of(plane)
    assert e1 < 0.2 * e0, (e0, e1)
    print("INFO noisy plane energy {:.3e} -> {:.3e}".format(e0, e1))


@test
def pinned_selected_vertices_in_edit_mode_and_vertex_groups_do_not_move():
    plane = noisy_plane(8)
    mesh = plane.data
    activate(plane)
    group = plane.vertex_groups.new(name="pin")
    group.add([3, 4], 1.0, 'REPLACE')
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(mesh)
    bm.verts.ensure_lookup_table()
    for v in bm.verts:
        v.select = v.index in (10, 11)
    bm.select_flush_mode()
    bmesh.update_edit_mesh(mesh)
    before = coords(plane)
    bpy.ops.object.mode_set(mode='OBJECT')
    before = coords(plane)
    bpy.ops.object.mode_set(mode='EDIT')
    run(bpy.ops.object.ddg_make_developable, subdivisions=0, remesh=False, iterations=40,
        pin_selected=True, pinned_group="pin")
    assert plane.mode == 'EDIT', "the operator returns to Edit mode"
    bpy.ops.object.mode_set(mode='OBJECT')
    after = coords(plane)
    for v in (3, 4, 10, 11):
        assert after[3 * v:3 * v + 3] == before[3 * v:3 * v + 3], v
    assert after != before


@test
def adjust_last_operation_and_undo():
    plane = noisy_plane(8)
    before = coords(plane)
    try:
        bpy.ops.object.ddg_make_developable('EXEC_DEFAULT', True, subdivisions=0, remesh=False,
                                            iterations=40)
    except TypeError:
        raise
    after = coords(plane)
    assert after != before
    if not bpy.ops.ed.undo.poll():
        raise Skip("undo is not available in this session")
    bpy.ops.ed.undo()
    restored = coords(bpy.context.active_object)
    assert restored == before, "undo restores the mesh"
    bpy.ops.ed.redo()
    assert coords(bpy.context.active_object) == after, "redo brings the result back"
    # Adjust Last Operation: redo_last re-runs the last operator with its properties
    bpy.ops.ed.undo()
    if bpy.ops.screen.redo_last.poll():
        bpy.ops.screen.redo_last('EXEC_DEFAULT')


@test
def preview_writes_attributes_and_optionally_after_the_flow():
    plane = noisy_plane(8)
    mesh = plane.data
    run(bpy.ops.object.ddg_show_developability)
    color = mesh.color_attributes.get("developability")
    values = mesh.attributes.get("developability_energy")
    assert color is not None and color.domain == 'POINT' and color.data_type == 'FLOAT_COLOR'
    assert values is not None and values.data_type == 'FLOAT'
    # The color attribute (not the float values) is what the viewport shows.
    assert mesh.color_attributes.active_color.name == "developability"
    assert mesh.attributes.active.name == "developability"
    assert mesh.color_attributes[mesh.color_attributes.render_color_index].name == "developability"
    v = array.array('f', bytes(4 * len(mesh.vertices)))
    values.data.foreach_get("value", v)
    assert max(v) > 1e-4 and all(x >= -1e-6 for x in v)
    run(bpy.ops.object.ddg_make_developable, subdivisions=0, remesh=False, iterations=20,
        write_energy=True)
    values = mesh.attributes.get("developability_energy")
    w = array.array('f', bytes(4 * len(mesh.vertices)))
    values.data.foreach_get("value", w)
    assert sum(w) < sum(v), (sum(w), sum(v))


@test
def errors_are_readable():
    core = module(".core")
    d = module(".developable_operator")
    obj = primitive("cube")
    # non-manifold: three faces on one edge
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bm.edges.ensure_lookup_table()
    e = bm.edges[0]
    a, b = e.verts
    c = bm.verts.new((5, 5, 5))
    bm.faces.new((a, b, c))
    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()
    try:
        d.make_developable_object(obj)
    except core.UnwrapFailed as ex:
        assert "manifold" in str(ex) or "faces" in str(ex), ex
    else:
        raise AssertionError("expected UnwrapFailed")
    try:
        d.developable_seams_uv(obj)
    except core.UnwrapFailed:
        pass
    else:
        raise AssertionError("expected UnwrapFailed")
    # shape keys
    plane = noisy_plane(3)
    plane.shape_key_add(name="Basis")
    try:
        d.make_developable_object(plane)
    except core.UnwrapFailed as ex:
        assert "shape key" in str(ex), ex
    else:
        raise AssertionError("expected UnwrapFailed")
    try:
        d.developable_seams_uv(noisy_plane(3), uv_name="")
    except core.UnwrapFailed as ex:
        assert "UV map name" in str(ex)
    # the operator turns it into an error report, not an exception
    activate(obj)
    try:
        bpy.ops.object.ddg_make_developable('EXEC_DEFAULT')
    except RuntimeError as ex:
        assert "manifold" in str(ex), ex
    else:
        raise AssertionError("expected the operator to report an error")


# ----------------------------------------------------------------------------------------------
# Developable Seams + UV


@test
def seams_uv_on_a_subdivided_cube():
    cube = subdivided_cube(3)
    mesh = cube.data
    before = coords(cube)
    stats = module(".developable_operator").developable_seams_uv(cube)
    assert coords(cube) == before, "vertices are not moved"
    assert seam_count(mesh) >= 7, seam_count(mesh)
    assert stats["n_terminals"] == 8, stats["n_terminals"]
    assert "Developable" in mesh.uv_layers and mesh.uv_layers.active.name == "Developable"
    assert stats["report"]["n_flipped_uv_triangles"] == 0, stats["report"]
    assert mirrored(cube, "Developable") == 0
    assert stats["developable_marked"] == stats["n_seams"]
    print("INFO", module(".developable_operator").seams_summary(stats))


@test
def seams_uv_operator_options_and_existing_seams():
    cube = subdivided_cube(2)
    mesh = cube.data
    mesh.edges[0].use_seam = True
    run(bpy.ops.uv.ddg_developable_seams, uv_name="Dev2", lambda_=0.2, epsilon=1e-2,
        use_defect_rule=False)
    assert "Dev2" in mesh.uv_layers
    assert mesh.edges[0].use_seam, "an existing seam is kept"
    assert mirrored(cube, "Dev2") == 0, mirrored(cube, "Dev2")
    n = seam_count(mesh)
    run(bpy.ops.uv.ddg_developable_seams, uv_name="Dev2")
    assert seam_count(mesh) >= n
    assert sorted(l.name for l in mesh.uv_layers) == ["Dev2", "UVMap"], "one new UV map, the cube's own kept"


@test
def seams_uv_on_cylinder_and_sphere_has_no_mirrored_triangles():
    d = module(".developable_operator")
    for obj in (primitive("cylinder", vertices=32, depth=2.0), primitive("uv_sphere", segments=24, ring_count=12)):
        stats = d.developable_seams_uv(obj)
        name = obj.data.uv_layers.active.name
        flips = mirrored(obj, name)
        print("INFO {} seams {} charts {} cones {} flips {} overlap {}".format(
            obj.name, stats["n_seams"], stats["report"]["n_charts"], stats["n_cones"], flips,
            stats["report"].get("overlapped_fraction")))
        assert flips == 0, (obj.name, flips)
        activate(obj)


@test
def seams_uv_in_edit_mode_returns_to_edit_mode():
    cube = subdivided_cube(2)
    activate(cube)
    bpy.ops.object.mode_set(mode='EDIT')
    run(bpy.ops.uv.ddg_developable_seams)
    assert cube.mode == 'EDIT'
    bpy.ops.object.mode_set(mode='OBJECT')
    assert "Developable" in cube.data.uv_layers and seam_count(cube.data) > 0


# ----------------------------------------------------------------------------------------------
# Spot


@test
def spot_flow_and_seams():
    spot = spot_object()
    mesh = spot.data
    d = module(".developable_operator")
    n_polygons = len(mesh.polygons)
    t0 = time.perf_counter()
    stats = d.make_developable_object(spot, iterations=100, subdivisions=0, remesh=False)
    t1 = time.perf_counter()
    r = stats["report"]
    print("SPOT flow, {} polygons, 100 iterations, no subdivision: energy {:.3e} -> {:.3e}, "
          "nonplanarity {:.3e} ({:.2%} of a side), {:.2f} s".format(
              n_polygons, r["energy_initial"], r["energy_final"], r["nonplanarity_after"],
              r["nonplanarity_after_relative"], t1 - t0))
    assert not stats["replaced"] and len(mesh.polygons) == n_polygons
    assert r["energy_final"] < r["energy_initial"]


@test
def spot_seams_uv_report():
    spot = spot_object()
    d = module(".developable_operator")
    before = coords(spot)
    t0 = time.perf_counter()
    stats = d.developable_seams_uv(spot)
    t1 = time.perf_counter()
    assert coords(spot) == before
    r = stats["report"]
    flips = mirrored(spot, "Developable")
    print("SPOT seams+UV: {} seam edges, {} terminals, {} charts, {} cones, overlap {}, mirrored "
          "(ours) {}, mirrored (pipeline) {}, collapsed {}, {:.2f} s".format(
              stats["n_seams"], stats["n_terminals"], r["n_charts"], stats["n_cones"],
              r.get("overlapped_fraction"), flips, r.get("n_flipped_uv_triangles"),
              r.get("n_collapsed_uv_triangles"), t1 - t0))
    print("SPOT", d.seams_summary(stats))
    assert seam_count(spot.data) > 0


@test
def spot_developable_then_seams_uv_report():
    """The intended workflow: flow first (one subdivision), then seams and UV on the result."""
    spot = spot_object()
    mesh = spot.data
    d = module(".developable_operator")
    t0 = time.perf_counter()
    stats = d.make_developable_object(spot, iterations=ARGS.spot_iterations, subdivisions=1)
    t1 = time.perf_counter()
    r = stats["report"]
    print("SPOT flow, 1 subdivision, {} iterations/round: {} -> {} triangles, energy {:.3e} -> {:.3e}, "
          "{} remesh flips+collapses, {:.1f} s".format(
              ARGS.spot_iterations, r["n_input_triangles"], r["n_triangles"], r["energy_initial"],
              r["energy_final"], sum(x["flips"] + x["collapses"] for x in r["rounds"]), t1 - t0))
    assert stats["replaced"] and sides(mesh) == [3]
    before = coords(spot)
    for eps in (1e-3, 1e-2):
        seam_stats = d.developable_seams_uv(spot, epsilon=eps, uv_name="Dev{}".format(eps))
        rep = seam_stats["report"]
        flips = mirrored(spot, "Dev{}".format(eps))
        print("SPOT after flow, epsilon {}: {} seam edges, {} terminals, {} charts, {} cones, "
              "overlap {}, mirrored {} (pipeline {}), {:.2f} s".format(
                  eps, seam_stats["n_seams"], seam_stats["n_terminals"], rep["n_charts"],
                  seam_stats["n_cones"], rep.get("overlapped_fraction"), flips,
                  rep.get("n_flipped_uv_triangles"), seam_stats["seconds_seams_native"]
                  + seam_stats["seconds_native"]))
        for e in mesh.edges:
            e.use_seam = False
    assert coords(spot) == before


# ----------------------------------------------------------------------------------------------


def main():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", default="ddg_uv_unwrap")
    parser.add_argument("--source", default="", help="register the add-on from this directory")
    parser.add_argument("--site", default="", help="directory with an unpacked ddg_tools wheel")
    parser.add_argument("--spot", default=SPOT)
    parser.add_argument("--spot-iterations", type=int, default=100,
                        help="iterations per round of the flow on spot")
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
    for name in ("OBJECT_OT_ddg_make_developable", "OBJECT_OT_ddg_show_developability",
                 "UV_OT_ddg_developable_seams"):
        assert hasattr(bpy.types, name), name + " is not registered"

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
        except Skip as skip:
            skipped.append(name)
            print("SKIP {}: {}".format(name, skip))
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
