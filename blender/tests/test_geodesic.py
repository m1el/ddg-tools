# SPDX-License-Identifier: MIT
"""Headless tests of Geodesic Path (`mesh.ddg_geodesic_path`) and Straighten Seams
(`mesh.ddg_straighten_seams`). Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/test_geodesic.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [-k NAME]

`scripts/blender-test.sh` runs it against the installed extension and its bundled wheel. For
development, `--source DIR --site DIR` register the add-on from a source directory (the parent
of the package) and import `ddg_tools` from an unpacked wheel instead.

Checked: every mode (Path, Smooth, Loop, Bezier) on a flat grid against the straight-line length
and de Casteljau, and on a sphere; the selection order from `bmesh.select_history`; the Curve
object (name, parent and transform, bevel depth, spline count, cyclic loops, custom properties);
Straighten Seams (junctions, extra pinned vertices, the stop ratio); Show triangulation (the
second curve object; it does not touch the objects of Show Intrinsic Refinement); replacement on
re-running; the mesh unchanged bit for bit; Edit mode kept; both representations; errors
(missing selection order, too few vertices, no seams, a crossing control polygon, mollify with
signposts, Object mode); the panel.
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


REPRESENTATIONS = ('SIGNPOSTS', 'INTEGER')


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


def grid(n=12):
    """A flat open grid in the plane z = 0, side 2, triangulated."""
    return activate(triangulate(primitive("grid", x_subdivisions=n, y_subdivisions=n, size=2.0)))


def fingerprint(mesh):
    """Everything about the mesh that the operators must leave alone."""
    co = [0.0] * (3 * len(mesh.vertices))
    mesh.vertices.foreach_get("co", co)
    loops = [0] * len(mesh.loops)
    mesh.loops.foreach_get("vertex_index", loops)
    seams = [False] * len(mesh.edges)
    mesh.edges.foreach_get("use_seam", seams)
    return (len(mesh.vertices), len(mesh.edges), len(mesh.polygons), tuple(co), tuple(loops),
            tuple(seams), tuple(a.name for a in mesh.attributes if not a.name.startswith(".")))


def positions(obj):
    return [tuple(v.co) for v in obj.data.vertices]


def vertex_at(obj, x, y, z=0.0):
    """The vertex nearest to (x, y, z) in object space."""
    return min(range(len(obj.data.vertices)),
               key=lambda i: math.dist(obj.data.vertices[i].co, (x, y, z)))


def select_in_order(obj, vertices):
    """Edit mode, vertex select mode, exactly these vertices selected, history in this order."""
    activate(obj)
    if obj.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
    # Deselect in Object mode first: a fresh primitive is "all selected" without selection
    # attributes, and Blender (headless) loses a selection made only in Edit mode.
    mesh = obj.data
    for collection in (mesh.polygons, mesh.edges, mesh.vertices):
        collection.foreach_set("select", [False] * len(collection))
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type='VERT')
    bm = bmesh.from_edit_mesh(obj.data)
    bm.verts.ensure_lookup_table()
    for v in bm.verts:
        v.select_set(False)
    bm.select_history.clear()
    for i in vertices:
        bm.verts[i].select_set(True)
        bm.select_history.add(bm.verts[i])
    bmesh.update_edit_mesh(obj.data)


def box_select(obj, vertices):
    """Selected without a history (as a box select leaves it)."""
    select_in_order(obj, vertices)
    bm = bmesh.from_edit_mesh(obj.data)
    bm.select_history.clear()
    bmesh.update_edit_mesh(obj.data)


def geodesic(obj, vertices, **kwargs):
    select_in_order(obj, vertices)
    result = bpy.ops.mesh.ddg_geodesic_path(**kwargs)
    assert result == {'FINISHED'}, result
    assert obj.mode == 'EDIT', "the operator did not return to Edit mode"
    return curve_of(obj)


def seams(obj, **kwargs):
    result = bpy.ops.mesh.ddg_straighten_seams(**kwargs)
    assert result == {'FINISHED'}, result
    return curve_of(obj)


def curve_of(obj):
    curve = bpy.data.objects.get(obj.name + "_geodesic")
    assert curve is not None, "no curve object"
    assert curve.type == 'CURVE'
    return curve


def spline_points(curve):
    return [[tuple(p.co)[:3] for p in s.points] for s in curve.data.splines]


def polyline_length(points):
    return sum(math.dist(points[k], points[k + 1]) for k in range(len(points) - 1))


def total_length(curve):
    return sum(polyline_length(p) for p in spline_points(curve))


def expect_error(function, *needles):
    try:
        function()
    except RuntimeError as ex:
        text = str(ex)
        for needle in needles:
            assert needle in text, "message {!r} lacks {!r}".format(text, needle)
        return text
    raise AssertionError("expected an error containing {}".format(needles))


def distance_to_polyline(point, points):
    best = math.inf
    for a, b in zip(points, points[1:]):
        ab = [y - x for x, y in zip(a, b)]
        ap = [y - x for x, y in zip(a, point)]
        L = sum(c * c for c in ab)
        t = 0.0 if L == 0 else max(0.0, min(1.0, sum(x * y for x, y in zip(ab, ap)) / L))
        q = [x + t * c for x, c in zip(a, ab)]
        best = min(best, math.dist(q, point))
    return best


def de_casteljau_polygon(controls, rounds):
    if rounds == 0:
        return list(controls)
    levels = [list(controls)]
    while len(levels[-1]) > 1:
        last = levels[-1]
        levels.append([tuple((a + b) / 2 for a, b in zip(last[k], last[k + 1]))
                       for k in range(len(last) - 1)])
    left = [level[0] for level in levels]
    right = [level[-1] for level in reversed(levels)]
    return de_casteljau_polygon(left, rounds - 1) + de_casteljau_polygon(right, rounds - 1)[1:]


def lattice_line(start, end, step=0.2):
    """Cells of the grid from start to end along an axis, one edge apart."""
    n = round(math.dist(start, end) / step)
    return [(round(start[0] + (end[0] - start[0]) * k / n, 6),
             round(start[1] + (end[1] - start[1]) * k / n, 6)) for k in range(n + 1)]


def mark_seam_path(obj, cells):
    """Marks the mesh edges between consecutive vertices (nearest to the (x, y) cells) as seams."""
    if obj.mode == 'EDIT':
        bpy.ops.object.mode_set(mode='OBJECT')
    ids = [vertex_at(obj, x, y) for x, y in cells]
    wanted = {frozenset(p) for p in zip(ids, ids[1:])}
    found = 0
    for edge in obj.data.edges:
        if frozenset(edge.vertices) in wanted:
            edge.use_seam = True
            found += 1
    assert found == len(wanted), "the cells are not joined by mesh edges"
    return ids


# ----------------------------------------------------------------------------------------------
# Tests


@test
def operators_and_panel_are_registered():
    assert hasattr(bpy.types, "MESH_OT_ddg_geodesic_path")
    assert hasattr(bpy.types, "MESH_OT_ddg_straighten_seams")
    assert hasattr(bpy.types, "VIEW3D_PT_ddg_geodesics")
    assert hasattr(bpy.ops.mesh, "ddg_geodesic_path") and hasattr(bpy.ops.mesh, "ddg_straighten_seams")


@test
def path_on_a_flat_grid_is_the_straight_line():
    source = grid()
    a, b = vertex_at(source, -0.8, -0.5), vertex_at(source, 0.7, 0.9)
    pa, pb = source.data.vertices[a].co, source.data.vertices[b].co
    expected = math.dist(pa, pb)
    for representation in REPRESENTATIONS:
        curve = geodesic(source, [a, b], mode='PATH', representation=representation)
        assert curve.name == source.name + "_geodesic"
        splines = spline_points(curve)
        assert len(splines) == 1
        assert abs(polyline_length(splines[0]) - expected) < 1e-5, (representation,)
        assert abs(curve["ddg_length"] - expected) < 1e-5
        assert curve["ddg_status"] == "geodesic"
        assert curve["ddg_certificate"] == "ok", curve["ddg_certificate"]
        assert curve["ddg_representation"] == representation
        assert curve["ddg_mode"] == 'PATH'
        assert curve["ddg_initial_length"] >= curve["ddg_length"]
        assert tuple(splines[0][0]) == tuple(pa) and tuple(splines[0][-1]) == tuple(pb)
        assert all(abs(p[2]) < 1e-6 for p in splines[0]), "left the surface"
        assert not curve.data.splines[0].use_cyclic_u


@test
def curve_object_is_parented_with_the_transform_and_a_bevel():
    source = grid()
    source.location = (1.0, 2.0, 3.0)
    source.rotation_euler = (0.3, 0.2, 0.9)
    source.scale = (2.0, 2.0, 2.0)
    activate(source)
    a, b = vertex_at(source, -0.8, -0.5), vertex_at(source, 0.7, 0.9)
    curve = geodesic(source, [a, b])
    bpy.context.view_layer.update()
    assert curve.parent is source
    for i in range(4):
        for j in range(4):
            assert abs(curve.matrix_world[i][j] - source.matrix_world[i][j]) < 1e-5
    assert curve.matrix_parent_inverse == curve.matrix_parent_inverse.__class__.Identity(4)
    diagonal = math.dist((-1, -1, 0), (1, 1, 0))
    assert abs(curve.data.bevel_depth - 0.002 * diagonal) < 1e-6, curve.data.bevel_depth
    assert curve.data.bevel_depth > 0 and curve.data.dimensions == '3D'
    # The curve lies on the surface in world space as well.
    world = source.matrix_world
    point = world @ source.data.vertices[a].co
    assert (curve.matrix_world @ __import__("mathutils").Vector(spline_points(curve)[0][0]) - point).length < 1e-4
    # Moving the mesh moves the curve.
    source.location = (5.0, 0.0, 0.0)
    bpy.context.view_layer.update()
    moved = curve.matrix_world.translation
    assert abs(moved.x - 5.0) < 1e-5


@test
def smooth_mode_only_seeds_the_path():
    source = grid()
    a = vertex_at(source, -0.8, -0.8)
    b = vertex_at(source, 0.8, -0.8)
    c = vertex_at(source, 0.8, 0.8)
    pa, pb, pc = (source.data.vertices[i].co for i in (a, b, c))
    for representation in REPRESENTATIONS:
        through = geodesic(source, [a, b, c], mode='PATH', representation=representation)
        assert abs(total_length(through) - (math.dist(pa, pb) + math.dist(pb, pc))) < 1e-5
        assert any(tuple(p) == tuple(pb) for p in spline_points(through)[0]), "left a pinned point"
        smooth = geodesic(source, [a, b, c], mode='SMOOTH', representation=representation)
        assert abs(total_length(smooth) - math.dist(pa, pc)) < 1e-5, representation
        assert smooth["ddg_mode"] == 'SMOOTH'


@test
def loop_mode_closes_the_curve():
    source = sphere(3)
    vertices = [vertex_at(source, *p) for p in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    lengths = {}
    for representation in REPRESENTATIONS:
        curve = geodesic(source, vertices, mode='LOOP', representation=representation)
        splines = curve.data.splines
        assert len(splines) == 1 and splines[0].use_cyclic_u
        points = spline_points(curve)[0]
        assert len(points) >= 3 and points[0] != points[-1], "a cyclic spline must not repeat its start"
        for v in vertices:
            assert any(tuple(p) == tuple(source.data.vertices[v].co) for p in points)
        assert curve["ddg_certificate"] == "ok", curve["ddg_certificate"]
        lengths[representation] = curve["ddg_length"]
        # A closed geodesic triangle of the octant: three quarter circles, about 3 pi / 2.
        assert 4.0 < curve["ddg_length"] < 5.0, curve["ddg_length"]
    assert abs(lengths['SIGNPOSTS'] - lengths['INTEGER']) < 1e-4, lengths
    unpinned = geodesic(source, vertices, mode='LOOP', pin_points=False)
    assert unpinned["ddg_length"] <= lengths['SIGNPOSTS'] + 1e-6


@test
def bezier_mode_matches_de_casteljau_on_a_plane():
    source = grid(16)
    cells = [(-0.75, -0.75), (-0.4, 0.75), (0.4, 0.75), (0.75, -0.75)]
    controls = [vertex_at(source, *c) for c in cells]
    planar = [tuple(source.data.vertices[i].co) for i in controls]
    for representation in REPRESENTATIONS:
        for rounds in (2, 3):
            curve = geodesic(source, controls, mode='BEZIER', rounds=rounds,
                             representation=representation)
            expected = de_casteljau_polygon(planar, rounds)
            points = spline_points(curve)[0]
            for e in expected:
                assert distance_to_polyline(e, points) < 2e-5, (representation, rounds, e)
            assert abs(polyline_length(points) - polyline_length(expected)) < 1e-4
            assert curve["ddg_rounds"] == rounds
            assert curve["ddg_n_control_points"] == 3 * 2 ** rounds + 1
            assert curve["ddg_mode"] == 'BEZIER'
            assert curve["ddg_certificate"] == "ok"
            assert curve["ddg_n_splits"] + curve["ddg_n_snapped"] >= 1


@test
def selection_order_decides():
    source = grid()
    a, b, c = (vertex_at(source, *p) for p in ((-0.8, -0.8), (0.8, -0.8), (0.4, 0.8)))
    forward = geodesic(source, [a, b, c], mode='PATH')
    first = spline_points(forward)[0]
    forward_length = forward["ddg_length"]
    assert tuple(first[0]) == tuple(source.data.vertices[a].co)
    assert tuple(first[-1]) == tuple(source.data.vertices[c].co)
    other = geodesic(source, [b, a, c], mode='PATH')
    assert tuple(spline_points(other)[0][0]) == tuple(source.data.vertices[b].co)
    assert abs(forward_length - other["ddg_length"]) > 0.1, "the order did not matter"
    # Two vertices without a history are fine (a straight line either way).
    box_select(source, [a, b])
    assert bpy.ops.mesh.ddg_geodesic_path() == {'FINISHED'}


@test
def straighten_seams_shortens_a_bent_seam():
    source = grid(10)
    cells = lattice_line((-0.8, -0.8), (0.8, -0.8)) + lattice_line((0.8, -0.8), (0.8, 0.8))[1:]
    ids = mark_seam_path(source, cells)
    corner = ids[8]
    before = fingerprint(source.data)
    p = [source.data.vertices[i].co for i in ids]
    for representation in REPRESENTATIONS:
        bpy.ops.object.mode_set(mode='OBJECT')
        curve = seams(source, representation=representation)
        assert len(curve.data.splines) == 1
        assert abs(total_length(curve) - math.dist(p[0], p[-1])) < 1e-5, representation
        assert curve["ddg_mode"] == 'SEAMS' and curve["ddg_n_paths"] == 1
        assert curve["ddg_status"] == "geodesic" and curve["ddg_certificate"] == "ok"
    assert fingerprint(source.data) == before, "the mesh changed"
    # The corner selected in Edit mode is pinned: the seam keeps passing through it.
    select_in_order(source, [corner])
    curve = seams(source)
    assert source.mode == 'EDIT'
    assert abs(total_length(curve) - (math.dist(p[0], p[8]) + math.dist(p[8], p[-1]))) < 1e-5
    assert len(curve.data.splines) == 2, "a pinned vertex splits the seam into two paths"
    # Stop at 0.9: stops early on a seam that has to shorten by more than that.
    bpy.ops.object.mode_set(mode='OBJECT')
    early = seams(source, stop_ratio=0.9)
    assert early["ddg_status"] == "length_ratio", early["ddg_status"]
    assert early["ddg_length"] > math.dist(p[0], p[-1]) + 1e-4
    assert early["ddg_length"] <= 0.9 * early["ddg_initial_length"] + 1e-6


@test
def straighten_seams_keeps_junctions():
    source = grid(10)
    ids = []
    for leg in (lattice_line((0.0, 0.0), (-0.8, 0.0)) + lattice_line((-0.8, 0.0), (-0.8, 0.8))[1:],
                lattice_line((0.0, 0.0), (0.8, 0.0)),
                lattice_line((0.0, 0.0), (0.0, -0.8))):
        ids.append(mark_seam_path(source, leg))
    junction = source.data.vertices[ids[0][0]].co
    for representation in REPRESENTATIONS:
        curve = seams(source, representation=representation)
        splines = spline_points(curve)
        assert len(splines) == 3
        ends = [s[0] for s in splines] + [s[-1] for s in splines]
        assert sum(1 for e in ends if tuple(e) == tuple(junction)) == 3, "the junction moved"
        end = source.data.vertices[ids[0][-1]].co
        expected = math.dist(junction, end) + 0.8 + 0.8
        assert abs(total_length(curve) - expected) < 1e-4, (total_length(curve), expected)


@test
def show_triangulation_adds_a_second_curve():
    source = grid()
    a, b = vertex_at(source, -0.8, -0.5), vertex_at(source, 0.7, 0.9)
    for representation in REPRESENTATIONS:
        curve = geodesic(source, [a, b], show_triangulation=True, representation=representation)
        edges = bpy.data.objects.get(source.name + "_geodesic_triangulation")
        assert edges is not None and edges.type == 'CURVE'
        assert len(edges.data.splines) == edges["ddg_n_edges"] >= len(source.data.edges)
        assert edges["ddg_intrinsic_of"] == source.name
        # The drawn edges lie on the surface of the input.
        assert all(abs(p[2]) < 1e-5 for s in spline_points(edges) for p in s)
        assert curve.parent is source
    # Running without it removes the second curve; the first stays.
    geodesic(source, [a, b], show_triangulation=False)
    assert bpy.data.objects.get(source.name + "_geodesic_triangulation") is None
    assert bpy.data.objects.get(source.name + "_geodesic") is not None


@test
def show_triangulation_leaves_the_intrinsic_refinement_alone():
    source = grid()
    bpy.ops.object.mode_set(mode='OBJECT')
    assert bpy.ops.mesh.ddg_show_intrinsic_refinement() == {'FINISHED'}
    assert bpy.data.objects.get(source.name + "_intrinsic_edges") is not None
    a, b = vertex_at(source, -0.8, -0.5), vertex_at(source, 0.7, 0.9)
    geodesic(source, [a, b], show_triangulation=True)
    assert bpy.data.objects.get(source.name + "_intrinsic_edges") is not None, "was replaced"
    assert bpy.data.objects.get(source.name + "_intrinsic") is not None
    assert bpy.data.objects.get(source.name + "_geodesic_triangulation") is not None
    bpy.ops.object.mode_set(mode='OBJECT')
    assert bpy.ops.mesh.ddg_show_intrinsic_refinement() == {'FINISHED'}
    assert bpy.data.objects.get(source.name + "_geodesic_triangulation") is not None, "was replaced"
    assert bpy.data.objects.get(source.name + "_geodesic") is not None


@test
def running_again_replaces_the_curve():
    source = grid(10)
    a, b, c = (vertex_at(source, *p) for p in ((-0.8, -0.8), (0.8, -0.8), (0.0, 0.8)))
    first = geodesic(source, [a, b])
    n_curves = len(bpy.data.curves)
    second = geodesic(source, [a, c], show_triangulation=True)
    third = geodesic(source, [b, c], show_triangulation=True)
    objects = [o for o in bpy.data.objects if o.get("ddg_intrinsic_kind") == "geodesic"]
    assert len(objects) == 1 and objects[0] is third
    assert len([o for o in bpy.data.objects if o.name.startswith(source.name + "_geodesic")]) == 2
    assert len(bpy.data.curves) == n_curves + 1, "curve data leaked"
    assert abs(third["ddg_length"] - math.dist(source.data.vertices[b].co,
                                               source.data.vertices[c].co)) < 1e-5
    # Straighten Seams replaces the same object.
    bpy.ops.object.mode_set(mode='OBJECT')
    mark_seam_path(source, lattice_line((-0.8, -0.8), (-0.4, -0.8)))
    fourth = seams(source)
    assert len([o for o in bpy.data.objects if o.get("ddg_intrinsic_kind") == "geodesic"]) == 1
    assert fourth["ddg_mode"] == 'SEAMS'


@test
def the_mesh_is_never_changed():
    source = sphere(2)
    ids = [vertex_at(source, *p) for p in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0))]
    select_in_order(source, ids)
    bpy.ops.object.mode_set(mode='OBJECT')
    before = fingerprint(source.data)
    for mode in ('PATH', 'SMOOTH', 'LOOP', 'BEZIER'):
        for representation in REPRESENTATIONS:
            select_in_order(source, ids)
            bpy.ops.mesh.ddg_geodesic_path(mode=mode, representation=representation,
                                           show_triangulation=True, rounds=2)
            bpy.ops.object.mode_set(mode='OBJECT')
            assert fingerprint(source.data) == before, (mode, representation)
    for e in source.data.edges:
        if e.index % 7 == 0:
            e.use_seam = True
    before = fingerprint(source.data)
    assert bpy.ops.mesh.ddg_straighten_seams(show_triangulation=True) == {'FINISHED'}
    assert fingerprint(source.data) == before


@test
def both_representations_agree_on_a_sphere():
    source = sphere(3)
    ids = [vertex_at(source, *p) for p in ((1, 0, 0), (0, 1, 0), (0, 0, -1))]
    for mode in ('PATH', 'SMOOTH', 'BEZIER'):
        lengths = [geodesic(source, ids, mode=mode, representation=r, rounds=3)["ddg_length"]
                   for r in REPRESENTATIONS]
        assert abs(lengths[0] - lengths[1]) < 1e-4, (mode, lengths)
    plain = geodesic(source, ids, mode='PATH', representation='INTEGER')["ddg_length"]
    mollified = geodesic(source, ids, mode='PATH', representation='INTEGER', mollify=True)
    assert mollified["ddg_mollify"]
    assert abs(mollified["ddg_length"] - plain) < 1e-3 * plain


@test
def errors_are_messages():
    source = grid()
    a, b, c = (vertex_at(source, *p) for p in ((-0.8, -0.8), (0.8, -0.8), (0.0, 0.8)))
    # Too few vertices.
    select_in_order(source, [a])
    expect_error(lambda: bpy.ops.mesh.ddg_geodesic_path(), "at least 2")
    select_in_order(source, [])
    expect_error(lambda: bpy.ops.mesh.ddg_geodesic_path(), "at least 2")
    assert source.mode == 'EDIT'
    # More than two vertices without a selection order.
    box_select(source, [a, b, c])
    expect_error(lambda: bpy.ops.mesh.ddg_geodesic_path(), "shift-click", "order")
    assert source.mode == 'EDIT'
    # A control polygon whose shortest paths cross.
    d = vertex_at(source, -0.8, 0.8)
    e = vertex_at(source, 0.8, 0.8)
    select_in_order(source, [a, e, d, b])
    for representation in REPRESENTATIONS:
        expect_error(lambda: bpy.ops.mesh.ddg_geodesic_path(
            mode='BEZIER', representation=representation), "move a control point")
        assert source.mode == 'EDIT'
    # Mollify needs integer coordinates.
    select_in_order(source, [a, b])
    expect_error(lambda: bpy.ops.mesh.ddg_geodesic_path(mollify=True, representation='SIGNPOSTS'),
                 "Mollify")
    # No seams.
    expect_error(lambda: bpy.ops.mesh.ddg_straighten_seams(), "no seams", "Mark Seam")
    assert source.mode == 'EDIT'
    # Object mode: the operator is not available.
    bpy.ops.object.mode_set(mode='OBJECT')
    expect_error(lambda: bpy.ops.mesh.ddg_geodesic_path(), "poll")
    # Nothing was made by any of the failures.
    assert bpy.data.objects.get(source.name + "_geodesic") is None


@test
def the_selection_survives_the_operators():
    source = grid()
    a, b, c = (vertex_at(source, *p) for p in ((-0.8, -0.8), (0.8, -0.8), (0.0, 0.8)))
    for work in (lambda: bpy.ops.mesh.ddg_geodesic_path(mode='PATH'),
                 lambda: expect_error(lambda: bpy.ops.mesh.ddg_geodesic_path(
                     mollify=True, representation='SIGNPOSTS'), "Mollify")):
        select_in_order(source, [c, a, b])
        work()
        bm = bmesh.from_edit_mesh(source.data)
        assert [v.index for v in bm.select_history] == [c, a, b], "the selection order was lost"
        assert sorted(v.index for v in bm.verts if v.select) == sorted([a, b, c])


@test
def a_failure_keeps_the_previous_curve():
    source = grid()
    a, b = vertex_at(source, -0.8, -0.8), vertex_at(source, 0.8, 0.8)
    curve = geodesic(source, [a, b])
    select_in_order(source, [a])
    expect_error(lambda: bpy.ops.mesh.ddg_geodesic_path(), "at least 2")
    assert bpy.data.objects.get(source.name + "_geodesic") is curve


@test
def a_closed_mesh_with_quads_and_ngons():
    # The operators triangulate as Blender does: a cube (quads) works and stays unchanged.
    source = activate(primitive("cube", size=2.0))
    ids = [vertex_at(source, *p) for p in ((-1, -1, -1), (1, 1, 1))]
    before = fingerprint(source.data)
    curve = geodesic(source, ids, mode='PATH')
    bpy.ops.object.mode_set(mode='OBJECT')
    assert fingerprint(source.data) == before
    assert curve["ddg_length"] > math.dist((-1, -1, -1), (1, 1, 1)) - 1e-6
    assert curve["ddg_certificate"] == "ok"


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
