# SPDX-License-Identifier: MIT
"""Headless tests of the Repulsive Curves operators. Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/test_repulsive.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [-k NAME]

`scripts/blender-test.sh` runs it against the installed extension and its bundled wheel. For
development, `--source DIR --site DIR` register the add-on from a source directory (the parent
of the package) and import `ddg_tools` from an unpacked wheel instead.

Checked: a trefoil edge mesh flows (energy down, edge lengths kept, undo and redo); a Bezier curve
round trip (poly write-back keeps the point count; Keep Bezier gives a Bezier spline; NURBS;
several splines; sampled Beziers); pins (vertex group, selected control points) stay; an
obstacle plane is never crossed, in world space with transformed objects, through the object and
the collection property; the live flow (``LiveFlow``: worker thread, snapshots, accept, Esc restore
to the exact original positions, also after live subdivision); live subdivision keeps vertex
groups; space and planar options; bad input (faces, empty mesh, coincident vertices, loose
vertices, a degenerate obstacle, wrong object type) gives a message and CANCELLED, not a crash.

NOT testable headless (see blender/README.md): the modal operator's event loop (timer events, Enter,
right click and Esc keys, the header text and the sidebar panel's Accept / Cancel buttons).
Timers do not run in `-b`; `ui_repulsive.py` drives them in a windowed Blender. What the modal
operator does on each event is `LiveFlow.poll` / `finish` / `cancel`, tested here with a real
worker thread.
"""

import argparse
import array
import importlib
import math
import sys
import time
import traceback

import bmesh
import bpy
from mathutils import Vector

TESTS = []


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
    for curve in list(bpy.data.curves):
        bpy.data.curves.remove(curve)
    for collection in list(bpy.data.collections):
        bpy.data.collections.remove(collection)
    settings = bpy.context.scene.ddg_repulsive
    for name in OPTIONS:
        settings.property_unset(name)


def module(name=""):
    return importlib.import_module(ARGS.module + name)


def ddg():
    return module(".core").native()


def link(obj):
    bpy.context.collection.objects.link(obj)
    view_layer = bpy.context.view_layer
    view_layer.update()
    for other in view_layer.objects:
        other.select_set(other is obj)
    view_layer.objects.active = obj
    return obj


def trefoil_points(n=60, scale=1.0):
    pts = []
    for i in range(n):
        t = 2.0 * math.pi * i / n
        pts.append((scale * (math.sin(t) + 2.0 * math.sin(2.0 * t)),
                    scale * (math.cos(t) - 2.0 * math.cos(2.0 * t)),
                    scale * -math.sin(3.0 * t)))
    return pts


def edge_mesh(points, edges, name="curve"):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(points, edges, [])
    mesh.update()
    return link(bpy.data.objects.new(name, mesh))


def ring_edges(n):
    return [(i, (i + 1) % n) for i in range(n)]


def trefoil(n=60, scale=1.0):
    return edge_mesh(trefoil_points(n, scale), ring_edges(n), "trefoil")


def coords(obj):
    out = array.array('f', bytes(12 * len(obj.data.vertices)))
    obj.data.vertices.foreach_get("co", out)
    return out


def pts_of(obj):
    c = coords(obj)
    return [tuple(c[i:i + 3]) for i in range(0, len(c), 3)]


def edge_array(obj):
    out = array.array('I', bytes(8 * len(obj.data.edges)))
    obj.data.edges.foreach_get("vertices", out)
    return [(out[i], out[i + 1]) for i in range(0, len(out), 2)]


def dist(a, b):
    return math.dist(a, b)


def edge_lengths(points, edges):
    return [dist(points[a], points[b]) for a, b in edges]


def ratio_error(a, b):
    """Largest |a_i / b_i - 1|."""
    return max(abs(x / y - 1.0) for x, y in zip(a, b))


def max_move(p, q, indices=None):
    indices = range(len(p)) if indices is None else indices
    return max(dist(p[i], q[i]) for i in indices)


def energy_of(points, edges, **kw):
    flat = array.array('d', [c for p in points for c in p])
    e, _ = ddg().curve_energy(flat, array.array('i', [v for pair in edges for v in pair]), **kw)
    return e


def world(obj, points):
    m = obj.matrix_world
    return [tuple(m @ Vector(p)) for p in points]


def run(**kwargs):
    result = bpy.ops.object.ddg_repulsive_flow('EXEC_DEFAULT', **kwargs)
    assert result == {'FINISHED'}, result


def refused(**kwargs):
    """The N-steps operator refuses: in a background Blender an operator that reports an error
    raises RuntimeError with the report. Returns the report text."""
    try:
        result = bpy.ops.object.ddg_repulsive_flow('EXEC_DEFAULT', **kwargs)
    except RuntimeError as ex:
        return str(ex)
    raise AssertionError("expected the operator to refuse, got {}".format(result))


def message_for(obj, **options):
    """The user message ``make_session`` raises for ``obj`` with the options (None: no error)."""
    op = module(".repulsive_operator")
    io = module(".repulsive_io")
    settings = bpy.context.scene.ddg_repulsive
    for k, v in options.items():
        setattr(settings, k, v)
    try:
        net, _session = op.make_session(obj, settings, bpy.context.evaluated_depsgraph_get())
        net.discard_backup()
    except io.CurveInputError as ex:
        return str(ex)
    return None


def make_bezier_curve(points, cyclic=True, kind='BEZIER'):
    curve = bpy.data.curves.new("curve", 'CURVE')
    curve.dimensions = '3D'
    spline = curve.splines.new(kind)
    n = len(points)
    if kind == 'BEZIER':
        spline.bezier_points.add(n - 1)
        for p, co in zip(spline.bezier_points, points):
            p.co = co
            p.handle_left_type = p.handle_right_type = 'AUTO'
    else:
        spline.points.add(n - 1)
        for p, co in zip(spline.points, points):
            p.co = (*co, 1.0)
    spline.use_cyclic_u = cyclic
    return link(bpy.data.objects.new("curve", curve))


def spline_points(spline):
    if spline.type == 'BEZIER':
        return [tuple(p.co) for p in spline.bezier_points]
    return [tuple(p.co)[:3] for p in spline.points]


def wait_for(flow, iterations, timeout=60.0):
    """Polls (as the timer does) until the session has completed ``iterations`` steps or ends."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        flow.poll()
        if flow.session.iteration >= iterations or not flow.running:
            return
        time.sleep(0.005)
    raise AssertionError("the flow did not reach step {} in {} s".format(iterations, timeout))


def session_points(session):
    flat = session.positions
    return [tuple(flat[i:i + 3]) for i in range(0, len(flat), 3)]


def f32(points):
    return [tuple(array.array('f', p)) for p in points]


# ----------------------------------------------------------------------------------------------
# Registration


@test
def operators_and_settings_are_registered():
    for name in ("OBJECT_OT_ddg_repulsive_flow", "OBJECT_OT_ddg_repulsive_live",
                 "OBJECT_OT_ddg_repulsive_live_finish", "VIEW3D_PT_ddg_repulsive"):
        assert hasattr(bpy.types, name), name + " is not registered"
    settings = bpy.context.scene.ddg_repulsive
    assert settings.alpha == 3.0 and settings.beta == 6.0 and settings.backend == 'auto'
    assert settings.subdivide is False, "live subdivision is off by default"
    assert settings.length == 'EDGES' and settings.space == 'WORLD'
    props = bpy.ops.object.ddg_repulsive_flow.get_rna_type().properties
    for name in OPTIONS:
        assert name in props, name


# ----------------------------------------------------------------------------------------------
# Trefoil: energy goes down, lengths are kept


@test
def trefoil_flows_energy_down_edge_lengths_kept():
    # squashed in z, so that the projection almost crosses itself: high energy to start with
    obj = edge_mesh([(x, y, 0.3 * z) for x, y, z in trefoil_points(40)], ring_edges(40))
    p0, e = pts_of(obj), edge_array(obj)
    e0 = energy_of(p0, e)
    l0 = edge_lengths(p0, e)
    run(steps=25)
    p1 = pts_of(obj)
    e1 = energy_of(p1, e)
    assert e1 < 0.8 * e0, (e0, e1)
    err = ratio_error(edge_lengths(p1, e), l0)
    assert err < 2e-3, err
    assert max_move(p0, p1) > 1e-3
    centroid = lambda p: [sum(q[k] for q in p) / len(p) for k in range(3)]
    assert dist(centroid(p0), centroid(p1)) < 1e-2, "the barycenter does not drift"
    print("INFO trefoil energy {:.4f} -> {:.4f}, edge length error {:.1e}".format(e0, e1, err))


@test
def to_convergence_with_zero_steps_and_status_report():
    obj = trefoil(24)
    p0, e = pts_of(obj), edge_array(obj)
    op = module(".repulsive_operator")
    settings = bpy.context.scene.ddg_repulsive
    settings.max_steps = 400
    stats = op.flow_object(obj, settings, bpy.context.evaluated_depsgraph_get(), steps=0)
    assert stats["status"] in {"converged", "max_steps", "stuck"}, stats["status"]
    assert stats["energy_after"] < stats["energy_before"]
    assert stats["steps"] > 5
    line = op.summary(stats)
    assert "energy" in line and "steps" in line, line
    assert energy_of(pts_of(obj), e) < energy_of(p0, e)


@test
def free_length_and_total_length_options():
    obj = trefoil(40)
    e = edge_array(obj)
    p0 = pts_of(obj)
    run(steps=10, length='TOTAL')
    p1 = pts_of(obj)
    assert abs(sum(edge_lengths(p1, e)) / sum(edge_lengths(p0, e)) - 1.0) < 2e-3
    reset()
    obj = trefoil(40)
    p0 = pts_of(obj)
    run(steps=10, length='EDGES', grow_factor=1.5, grow_rate=0.05)
    p1 = pts_of(obj)
    assert sum(edge_lengths(p1, e)) > 1.05 * sum(edge_lengths(p0, e))
    reset()
    trefoil(40)
    run(steps=5, length='FREE')


@test
def undo_redo_and_one_undo_step():
    trefoil(30)
    before = coords(bpy.context.active_object)
    run(steps=8)
    after = coords(bpy.context.active_object)
    assert after != before
    if not bpy.ops.ed.undo.poll():
        raise Skip("undo is not available in this session")
    bpy.ops.ed.undo()
    assert coords(bpy.context.active_object) == before, "one undo undoes all steps"
    bpy.ops.ed.redo()
    assert coords(bpy.context.active_object) == after
    bpy.ops.ed.undo()
    if bpy.ops.screen.redo_last.poll():
        bpy.ops.screen.redo_last('EXEC_DEFAULT')


@test
def invoke_takes_the_scene_settings():
    trefoil(30)
    settings = bpy.context.scene.ddg_repulsive
    settings.steps = 3
    settings.alpha = 2.5
    settings.beta = 5.0
    assert bpy.ops.object.ddg_repulsive_flow('INVOKE_DEFAULT') == {'FINISHED'}
    operators = bpy.context.window_manager.operators
    last = operators[-1] if operators else None
    if last is not None and last.bl_idname == "OBJECT_OT_ddg_repulsive_flow":
        assert last.properties.steps == 3 and abs(last.properties.alpha - 2.5) < 1e-6


@test
def native_backends_are_accepted():
    for backend in ('auto', 'dense', 'barnes_hut', 'multigrid', 'accelerated'):
        reset()
        trefoil(30)
        run(steps=2, backend=backend)


# ----------------------------------------------------------------------------------------------
# Curves


def bezier_trefoil(n=12):
    return make_bezier_curve(trefoil_points(n, 1.0), cyclic=True)


@test
def bezier_poly_writeback_keeps_the_point_count():
    obj = bezier_trefoil(12)
    before = spline_points(obj.data.splines[0])
    run(steps=6)
    spline = obj.data.splines[0]
    assert spline.type == 'POLY', spline.type
    assert len(spline.points) == 12 and spline.use_cyclic_u
    after = spline_points(spline)
    assert max_move(before, after) > 1e-4
    # exact: the polyline is the flow's output, with the lengths it kept
    e = ring_edges(12)
    assert ratio_error(edge_lengths(after, e), edge_lengths(before, e)) < 2e-3


@test
def keep_bezier_gives_a_bezier_spline_and_nurbs_stays_nurbs():
    obj = bezier_trefoil(12)
    before = spline_points(obj.data.splines[0])
    run(steps=6, keep_bezier=True)
    spline = obj.data.splines[0]
    assert spline.type == 'BEZIER' and len(spline.bezier_points) == 12 and spline.use_cyclic_u
    after = spline_points(spline)
    assert max_move(before, after) > 1e-4
    right = [tuple(p.handle_right) for p in spline.bezier_points]
    left = [tuple(p.handle_left) for p in spline.bezier_points]
    assert max(dist(r, p) for r, p in zip(right, after)) > 1e-3, "the handles are not degenerate"
    # the handles are aligned: left, point, right on one line (and the point between them)
    for l, p, r in zip(left, after, right):
        u = [r[k] - p[k] for k in range(3)]
        v = [p[k] - l[k] for k in range(3)]
        cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        assert max(abs(c) for c in cross) < 1e-3
    reset()
    obj = make_bezier_curve(trefoil_points(12), cyclic=True, kind='NURBS')
    run(steps=4, keep_bezier=True)
    assert obj.data.splines[0].type == 'NURBS' and len(obj.data.splines[0].points) == 12
    reset()
    obj = make_bezier_curve(trefoil_points(12), cyclic=True, kind='NURBS')
    run(steps=4)
    assert obj.data.splines[0].type == 'POLY' and len(obj.data.splines[0].points) == 12


@test
def poly_curve_roundtrip_and_several_splines():
    curve = bpy.data.curves.new("two", 'CURVE')
    curve.dimensions = '3D'
    for k, shift in enumerate((-3.0, 3.0)):
        spline = curve.splines.new('POLY')
        pts = trefoil_points(20, 0.5)
        spline.points.add(len(pts) - 1)
        for p, co in zip(spline.points, pts):
            p.co = (co[0] + shift, co[1], co[2], 1.0)
        spline.use_cyclic_u = True
        spline.material_index = k
    for name in ("a", "b"):
        curve.materials.append(bpy.data.materials.new(name))
    for k, spline in enumerate(curve.splines):
        spline.material_index = k
    link(bpy.data.objects.new("two", curve))
    run(steps=5)
    assert len(curve.splines) == 2
    for k, spline in enumerate(curve.splines):
        assert spline.type == 'POLY' and len(spline.points) == 20 and spline.use_cyclic_u
        assert spline.material_index == k
    # an open spline is an open polyline (n - 1 edges)
    reset()
    curve = bpy.data.curves.new("open", 'CURVE')
    curve.dimensions = '3D'
    spline = curve.splines.new('POLY')
    spline.points.add(29)
    for i, p in enumerate(spline.points):
        t = i / 29.0
        p.co = (3.0 * t, math.sin(6 * t), 0.3 * math.cos(9 * t), 1.0)
    link(bpy.data.objects.new("open", curve))
    run(steps=4)
    assert not curve.splines[0].use_cyclic_u and len(curve.splines[0].points) == 30


@test
def sampled_bezier_writes_more_points():
    obj = bezier_trefoil(8)
    run(steps=2, bezier_samples=3)
    spline = obj.data.splines[0]
    assert spline.type == 'POLY' and len(spline.points) == 8 * 3
    reset()
    obj = bezier_trefoil(8)
    run(steps=2, bezier_samples=3, keep_bezier=True)
    assert obj.data.splines[0].type == 'BEZIER' and len(obj.data.splines[0].bezier_points) == 24


@test
def selected_control_points_are_pinned():
    obj = bezier_trefoil(12)
    spline = obj.data.splines[0]
    for i, p in enumerate(spline.bezier_points):
        p.select_control_point = i in (0, 5)
    before = spline_points(spline)
    run(steps=10)
    after = spline_points(obj.data.splines[0])
    assert max_move(before, after, [0, 5]) < 1e-6
    assert max_move(before, after) > 1e-3
    # a NURBS curve with a kept tangent at the pin
    reset()
    obj = make_bezier_curve(trefoil_points(12), cyclic=True, kind='NURBS')
    for i, p in enumerate(obj.data.splines[0].points):
        p.select = i == 3
    before = spline_points(obj.data.splines[0])
    run(steps=4, pin_tangents=True)
    after = spline_points(obj.data.splines[0])
    assert max_move(before, after, [3]) < 1e-6


@test
def curve_in_edit_mode_is_handled():
    obj = bezier_trefoil(12)
    bpy.ops.object.mode_set(mode='EDIT')
    run(steps=2)
    assert obj.mode == 'EDIT', "the mode is restored"
    bpy.ops.object.mode_set(mode='OBJECT')


@test
def curve_live_subdivision_rebuilds_the_splines():
    obj = make_bezier_curve(trefoil_points(10, 1.0), cyclic=True, kind='POLY')
    run(steps=40, subdivide=True, max_subdivisions=1, length='EDGES', grow_factor=2.5,
        grow_rate=0.2)
    spline = obj.data.splines[0]
    assert len(spline.points) == 20 and spline.use_cyclic_u, len(spline.points)
    pts = spline_points(spline)
    assert all(math.isfinite(c) for p in pts for c in p)
    # consecutive points are consecutive along the curve: no long jumps
    steps = sorted(dist(pts[i], pts[(i + 1) % 20]) for i in range(20))
    assert steps[-1] < 3.0 * steps[10], (steps[-1], steps[10])


# ----------------------------------------------------------------------------------------------
# Pins on meshes


@test
def vertex_group_pins_stay_fixed():
    obj = trefoil(48)
    group = obj.vertex_groups.new(name="pins")
    group.add([0, 12, 30], 1.0, 'REPLACE')
    before = pts_of(obj)
    run(steps=20, pin_group="pins")
    after = pts_of(obj)
    assert max_move(before, after, [0, 12, 30]) < 1e-6
    assert max_move(before, after) > 1e-3, "the rest moved"
    # a weight of 0 is not a pin
    reset()
    obj = trefoil(48)
    group = obj.vertex_groups.new(name="pins")
    group.add([0], 0.0, 'REPLACE')
    group.add([12], 0.5, 'REPLACE')
    before = pts_of(obj)
    run(steps=10, pin_group="pins")
    after = pts_of(obj)
    assert max_move(before, after, [12]) < 1e-6
    assert max_move(before, after, [0]) > 1e-6


@test
def tangent_group_keeps_tangents():
    obj = trefoil(48)
    g = obj.vertex_groups.new(name="t")
    g.add([5], 1.0, 'REPLACE')
    before = pts_of(obj)
    run(steps=15, tangent_group="t")
    after = pts_of(obj)
    n = len(before)

    def tangent(p, i):
        d = [p[(i + 1) % n][k] - p[(i - 1) % n][k] for k in range(3)]
        length = math.sqrt(sum(c * c for c in d))
        return [c / length for c in d]

    dot = sum(a * b for a, b in zip(tangent(before, 5), tangent(after, 5)))
    assert dot > 0.999, dot


@test
def unknown_group_is_an_error():
    obj = trefoil(30)
    msg = message_for(obj, pin_group="nope")
    assert msg and "nope" in msg and "vertex group" in msg, msg
    refused(steps=2, pin_group="nope")


# ----------------------------------------------------------------------------------------------
# Obstacles


def plane_obstacle(z=0.0, size=40.0, location=None, name="plane"):
    bpy.ops.mesh.primitive_plane_add(size=size, location=(0, 0, z) if location is None else location)
    obj = bpy.context.view_layer.objects.active
    obj.name = name
    return obj


def trefoil_above(z0):
    obj = trefoil(60, 1.0)
    obj.location.z = z0 - min(p[2] for p in pts_of(obj))
    bpy.context.view_layer.update()
    return obj


@test
def an_obstacle_plane_is_never_crossed():
    plane = plane_obstacle(0.0)
    obj = trefoil_above(0.05)
    settings = bpy.context.scene.ddg_repulsive
    settings.obstacle_object = plane.name
    settings.grow_factor, settings.grow_rate = 1.5, 0.05  # pushes the curve into the plane
    op = module(".repulsive_operator")
    depsgraph = bpy.context.evaluated_depsgraph_get()
    net, session = op.make_session(obj, settings, depsgraph)
    lowest = []
    for _ in range(40):
        session.step()
        lowest.append(min(p[2] for p in session_points(session)))
    assert min(lowest) > 0.0, "crossed the plane: z = {}".format(min(lowest))
    assert session.iteration >= 30
    # Control: without the plane the same flow goes below it (so the obstacle acted)
    settings.obstacle_object = ""
    net2, session2 = op.make_session(obj, settings, depsgraph)
    for _ in range(40):
        session2.step()
    free = min(p[2] for p in session_points(session2))
    assert free < 0.0, "the control run does not reach the plane (z = {})".format(free)
    print("INFO lowest z with plane {:.4f}, without {:.4f}".format(min(lowest), free))


@test
def obstacle_through_the_operator_a_collection_and_world_transforms():
    # a tilted, moved plane in a collection; a rotated and scaled curve: all meet in world space
    plane = plane_obstacle(0.0, size=30.0, location=(0.5, 0.0, 0.0))
    plane.rotation_euler = (0.2, 0.0, 0.0)
    bpy.context.view_layer.update()
    coll = bpy.data.collections.new("obstacles")
    bpy.context.scene.collection.children.link(coll)
    coll.objects.link(plane)
    bpy.context.scene.collection.objects.unlink(plane)
    obj = trefoil(60, 1.0)
    obj.scale = (1.5, 1.5, 1.5)
    obj.rotation_euler = (0.3, 0.1, 0.0)
    bpy.context.view_layer.update()
    normal = plane.matrix_world.to_3x3() @ Vector((0, 0, 1))
    origin = plane.matrix_world.translation

    def sides():
        return [(Vector(p) - origin).dot(normal) for p in world(obj, pts_of(obj))]

    obj.location = obj.location + normal * (0.05 - min(sides()))
    bpy.context.view_layer.update()
    assert min(sides()) > 0.0
    e = edge_array(obj)
    lengths0 = edge_lengths(world(obj, pts_of(obj)), e)
    run(steps=30, obstacle_collection=coll.name, grow_factor=1.5, grow_rate=0.05)
    assert min(sides()) > -1e-5, "crossed the tilted plane in world space: {}".format(min(sides()))
    # the flow ran in world space: every world edge length grew by the same factor
    ratios = [a / b for a, b in zip(edge_lengths(world(obj, pts_of(obj)), e), lengths0)]
    assert max(ratios) / min(ratios) < 1.01 and 1.01 < min(ratios), (min(ratios), max(ratios))


@test
def a_degenerate_or_faceless_obstacle_gives_a_helpful_error():
    obj = trefoil(30)
    bad = bpy.data.meshes.new("degenerate")
    bad.from_pydata([(0, 0, 5), (1, 0, 5), (2, 0, 5), (0, 1, 5), (1, 1, 5)],
                    [], [(0, 1, 2), (0, 1, 3)])
    bad_obj = bpy.data.objects.new("bad_plane", bad)
    bpy.context.collection.objects.link(bad_obj)
    msg = message_for(obj, obstacle_object=bad_obj.name)
    assert msg and "bad_plane" in msg and "degenerate" in msg, msg
    assert "degenerate" in refused(obstacle_object=bad_obj.name, steps=2)
    empty = bpy.data.objects.new("wire", bpy.data.meshes.new("wire"))
    bpy.context.collection.objects.link(empty)
    msg = message_for(obj, obstacle_object=empty.name)
    assert msg and "wire" in msg and "no faces" in msg, msg


@test
def a_curve_that_starts_on_an_obstacle_does_not_crash():
    plane = plane_obstacle(0.0)
    obj = trefoil(30)  # passes through the plane z = 0
    before = coords(obj)
    msg = message_for(obj, obstacle_object=plane.name)
    # Either the library refuses at once (the message names the object) or the flow cannot make
    # a step ("stuck"); the object must not be damaged either way.
    try:
        result = bpy.ops.object.ddg_repulsive_flow(
            'EXEC_DEFAULT', obstacle_object=plane.name, steps=3)
        outcome = "ran: " + str(result)
    except RuntimeError as ex:
        outcome = "refused: " + str(ex)
        assert obj.name in str(ex) or "obstacle" in str(ex).lower(), ex
        assert coords(obj) == before
    print("INFO curve through obstacle: construction message {!r}; operator {}".format(msg, outcome))


# ----------------------------------------------------------------------------------------------
# Space and planar


@test
def planar_option_holds_one_axis():
    pts = []
    for i in range(50):
        t = 2.0 * math.pi * i / 50
        r = 1.0 + 0.3 * math.cos(3.0 * t)
        pts.append((r * math.cos(t), r * math.sin(t), 0.4))
    obj = edge_mesh(pts, ring_edges(50))
    run(steps=15, planar=True, planar_axis='Z')
    after = pts_of(obj)
    assert max(abs(p[2] - 0.4) for p in after) < 1e-5
    flat = lambda q: [(p[0], p[1], 0.0) for p in q]
    assert max_move(flat(pts), flat(after)) > 1e-3


@test
def world_and_local_space_differ_for_scaled_objects():
    results = {}
    e = ring_edges(40)
    for space in ('WORLD', 'LOCAL'):
        reset()
        obj = trefoil(40)
        obj.scale = (1.0, 1.0, 3.0)
        bpy.context.view_layer.update()
        run(steps=8, space=space)
        local = pts_of(obj)
        results[space] = (local, edge_lengths(local, e), edge_lengths(world(obj, local), e))
    # the constrained lengths are those of the chosen space
    p = trefoil_points(40)
    l_local = edge_lengths(p, e)
    l_world = edge_lengths([(x, y, 3.0 * z) for x, y, z in p], e)
    assert ratio_error(results['LOCAL'][1], l_local) < 2e-3
    assert ratio_error(results['WORLD'][2], l_world) < 2e-3
    assert max_move(results['LOCAL'][0], results['WORLD'][0]) > 1e-3


# ----------------------------------------------------------------------------------------------
# The live flow (everything the modal operator does, apart from events)


def live(obj, **options):
    op = module(".repulsive_operator")
    settings = bpy.context.scene.ddg_repulsive
    for k, v in options.items():
        setattr(settings, k, v)
    return op.LiveFlow.create(obj, settings, bpy.context.evaluated_depsgraph_get())


@test
def live_flow_runs_on_a_worker_thread_and_accepts():
    import threading

    obj = trefoil(60)
    p0, e = pts_of(obj), edge_array(obj)
    flow = live(obj)
    flow.start()
    assert flow.thread is not threading.current_thread() and flow.thread.is_alive()
    wait_for(flow, 15)
    stats = flow.finish()
    after = pts_of(obj)
    assert not flow.running
    assert stats["steps"] >= 5 and stats["energy_after"] < stats["energy_before"], stats
    assert energy_of(after, e) < energy_of(p0, e)
    # the object shows exactly the session's last completed step
    assert after == f32(session_points(flow.session))
    assert "Repulsive" in flow.status_line() and "energy" in flow.status_line()
    print("INFO live: {} steps, status {}".format(stats["steps"], stats["status"]))


@test
def live_flow_esc_restores_the_exact_original_positions():
    obj = trefoil(60)
    before = coords(obj)
    flow = live(obj)
    flow.start()
    wait_for(flow, 8)
    assert coords(obj) != before, "the object was updated while the flow ran"
    flow.cancel()
    assert coords(obj) == before, "bitwise the original positions"
    assert not flow.running


@test
def live_flow_esc_restores_a_transformed_object_and_a_curve():
    obj = trefoil(40)
    obj.location = (1.0, 2.0, 3.0)
    obj.scale = (2.0, 1.0, 0.5)
    obj.rotation_euler = (0.3, 0.2, 0.1)
    bpy.context.view_layer.update()
    before = coords(obj)
    flow = live(obj)
    flow.start()
    wait_for(flow, 5)
    flow.cancel()
    assert coords(obj) == before
    reset()
    obj = bezier_trefoil(12)
    sp = obj.data.splines[0]
    for i, p in enumerate(sp.bezier_points):
        p.handle_left_type = p.handle_right_type = 'ALIGNED'
        p.select_control_point = i == 2

    def state(spline):
        return [(tuple(p.co), tuple(p.handle_left), tuple(p.handle_right), p.select_control_point,
                 p.handle_left_type) for p in spline.bezier_points]

    snapshot = state(sp)
    flow = live(obj, keep_bezier=False)
    flow.start()
    wait_for(flow, 6)
    assert obj.data.splines[0].type == 'POLY', "the curve showed the poly result"
    flow.cancel()
    sp = obj.data.splines[0]
    assert sp.type == 'BEZIER' and len(sp.bezier_points) == 12 and sp.use_cyclic_u
    assert snapshot == state(sp)


@test
def live_flow_cancel_right_after_start_and_double_cancel_are_safe():
    obj = trefoil(60)
    before = coords(obj)
    flow = live(obj)
    flow.start()
    flow.cancel()
    assert coords(obj) == before
    flow.cancel()
    assert coords(obj) == before


@test
def live_flow_stops_on_convergence_and_waits():
    obj = trefoil(20)
    flow = live(obj, stop_tolerance=0.5, max_steps=500)
    flow.start()
    deadline = time.time() + 60
    while flow.running and time.time() < deadline:
        flow.poll()
        time.sleep(0.005)
    assert not flow.running, "the thread ends when the session is final"
    assert flow.session.status in {"converged", "stuck", "max_steps"}, flow.session.status
    n = flow.session.iteration
    time.sleep(0.05)
    assert flow.session.iteration == n, "no stepping after the end"
    stats = flow.finish()
    assert stats["steps"] == n


@test
def live_flow_pins_and_obstacles_hold_while_running():
    plane = plane_obstacle(0.0)
    obj = trefoil_above(0.05)
    group = obj.vertex_groups.new(name="pins")
    group.add([7], 1.0, 'REPLACE')
    before = pts_of(obj)
    flow = live(obj, obstacle_object=plane.name, pin_group="pins", grow_factor=1.5, grow_rate=0.05)
    flow.start()
    lowest = 1e9
    while flow.running and flow.session.iteration < 25:
        flow.poll()
        lowest = min(lowest, min(p[2] for p in world(obj, pts_of(obj))))
        time.sleep(0.005)
    flow.finish()
    after = pts_of(obj)
    assert max_move(before, after, [7]) < 1e-6
    final = min(p[2] for p in world(obj, after))
    # (the mesh stores float32: up to its rounding the curve does not reach the plane)
    assert lowest > -1e-6 and final > -1e-6, (lowest, final)
    assert min(p[2] for p in session_points(flow.session)) > 0.0, "the session is strictly above"


@test
def a_worker_error_is_kept_for_the_main_thread():
    obj = trefoil(30)
    flow = live(obj)

    def boom():
        raise RuntimeError("boom")

    flow.session = type("S", (), {"step": staticmethod(boom), "cancel": lambda self: None,
                                  "snapshot": lambda self: None})()
    flow.start()
    flow.thread.join()
    assert isinstance(flow.error, RuntimeError)


# ----------------------------------------------------------------------------------------------
# Live subdivision


def group_members(obj, name):
    index = obj.vertex_groups[name].index
    out = {}
    for v in obj.data.vertices:
        for g in v.groups:
            if g.group == index:
                out[v.index] = g.weight
    return out


def run_until_subdivided(flow, n):
    flow.start()
    deadline = time.time() + 60
    while flow.session.n_vertices == n and time.time() < deadline and flow.running:
        flow.poll()
        time.sleep(0.005)
    assert flow.session.n_vertices > n, "the session did not subdivide"
    flow.poll()


@test
def live_subdivision_keeps_vertex_groups():
    n = 12
    obj = trefoil(n, 1.0)
    a = obj.vertex_groups.new(name="a")
    a.add([0, 1, 2], 1.0, 'REPLACE')
    b = obj.vertex_groups.new(name="b")
    b.add([6], 0.25, 'REPLACE')
    flow = live(obj, subdivide=True, max_subdivisions=1, length='EDGES', grow_factor=3.0,
                grow_rate=0.3)
    run_until_subdivided(flow, n)
    flow.finish()
    mesh = obj.data
    assert len(mesh.vertices) == 2 * n and len(mesh.edges) == 2 * n, len(mesh.vertices)
    ga, gb = group_members(obj, "a"), group_members(obj, "b")
    # old vertices keep their indices and weights
    assert all(ga.get(v) == 1.0 for v in (0, 1, 2)), ga
    assert gb.get(6) == 0.25
    # the midpoints of edges (0,1) and (1,2) are vertices n + 0 and n + 1: both ends in the group
    assert ga.get(n + 0) == 1.0 and ga.get(n + 1) == 1.0, ga
    # a midpoint between a member and a non-member is not in the group
    assert (n + 2) not in ga and (n + 11) not in ga
    assert set(gb) == {6}
    # the edges are a closed loop through all vertices
    degree = [0] * (2 * n)
    for a_, b_ in edge_array(obj):
        degree[a_] += 1
        degree[b_] += 1
    assert degree == [2] * (2 * n)


@test
def live_subdivision_esc_restores_topology_groups_and_positions():
    n = 12
    obj = trefoil(n, 1.0)
    g = obj.vertex_groups.new(name="a")
    g.add([0, 3], 0.5, 'REPLACE')
    before = coords(obj)
    edges_before = edge_array(obj)
    flow = live(obj, subdivide=True, max_subdivisions=1, length='EDGES', grow_factor=3.0,
                grow_rate=0.3)
    run_until_subdivided(flow, n)
    assert len(obj.data.vertices) == 2 * n, "the object shows the subdivided curve"
    flow.cancel()
    assert len(obj.data.vertices) == n and coords(obj) == before
    assert edge_array(obj) == edges_before
    assert group_members(obj, "a") == {0: 0.5, 3: 0.5}
    assert not [m for m in bpy.data.meshes if m.name.startswith(".ddg_repulsive")], \
        "the backup is removed"


@test
def subdivision_is_off_by_default_so_vertices_stay():
    obj = trefoil(12)
    run(steps=30, grow_factor=3.0, grow_rate=0.3)
    assert len(obj.data.vertices) == 12


@test
def created_parents_undoes_rounds_in_order():
    io = module(".repulsive_io")
    # round 1: midpoints 4..7 of the ring 0-1-2-3 (edge I -> I and |E| + I)
    r1 = [(0, 4), (1, 5), (2, 6), (3, 7), (4, 1), (5, 2), (6, 3), (7, 0)]
    # round 2: midpoints 8..15
    r2 = []
    for i, (a, b) in enumerate(r1):
        r2 += [(a, 8 + i), (8 + i, b)]
    flat = [v for pair in r2 for v in pair]
    parents = io.created_parents(flat, 4, 16)
    assert parents[4] == (0, 1) and parents[7] == (0, 3), parents
    assert parents[8] == (0, 4) and parents[12] == (1, 4), parents


# ----------------------------------------------------------------------------------------------
# Meshes with extra stuff


@test
def mesh_with_faces_is_rejected_or_uses_its_loose_edges():
    bpy.ops.mesh.primitive_cube_add(size=2.0)
    cube = bpy.context.view_layer.objects.active
    msg = message_for(cube)
    assert msg and "6 faces" in msg and "Loose Edges Only" in msg, msg
    assert "faces" in refused(steps=2)
    # a cube with a loose loop next to it: with the option only the loop moves
    bm = bmesh.new()
    bm.from_mesh(cube.data)
    loop = [bm.verts.new((x + 6.0, y, z)) for x, y, z in trefoil_points(24)]
    for i in range(24):
        bm.edges.new((loop[i], loop[(i + 1) % 24]))
    bm.to_mesh(cube.data)
    bm.free()
    cube.data.update()
    before = pts_of(cube)
    run(steps=8, loose_edges_only=True)
    after = pts_of(cube)
    assert before[:8] == after[:8], "the cube did not move"
    assert max_move(before, after, range(8, 32)) > 1e-3
    assert message_for(cube, loose_edges_only=True) is None
    # subdivision cannot rebuild a partly used mesh
    msg = message_for(cube, loose_edges_only=True, subdivide=True)
    assert msg and "Subdivision" in msg, msg


@test
def loose_vertices_are_an_error_unless_ignored():
    pts = trefoil_points(30) + [(9.0, 9.0, 9.0)]
    obj = edge_mesh(pts, ring_edges(30))
    msg = message_for(obj)
    assert msg and "belong to no edge" in msg and "vertex 30" in msg, msg
    before = pts_of(obj)
    run(steps=5, ignore_loose=True)
    after = pts_of(obj)
    assert before[30] == after[30] and max_move(before, after, range(30)) > 1e-4


@test
def junctions_are_allowed():
    # a ring with a spoke: two degree-3 junctions
    pts = [(math.cos(2 * math.pi * i / 30), math.sin(2 * math.pi * i / 30), 0.0)
           for i in range(30)] + [(0.0, 0.0, 0.3)]
    obj = edge_mesh(pts, ring_edges(30) + [(0, 30), (15, 30)])
    before = pts_of(obj)
    run(steps=5)
    assert max_move(before, pts_of(obj)) > 1e-4


@test
def bad_input_gives_clear_messages_not_crashes():
    # empty mesh
    empty = link(bpy.data.objects.new("empty", bpy.data.meshes.new("empty")))
    msg = message_for(empty)
    assert msg and "empty" in msg, msg
    assert "empty" in refused()
    # vertices but no edges
    reset()
    obj = edge_mesh([(0, 0, 0), (1, 0, 0), (2, 1, 0)], [])
    msg = message_for(obj)
    assert msg and "empty" in msg, msg
    # coincident vertices
    reset()
    pts = trefoil_points(20)
    pts[10] = pts[3]
    obj = edge_mesh(pts, ring_edges(20))
    msg = message_for(obj)
    assert msg and "coincident" in msg and "3" in msg and "10" in msg and "Merge" in msg, msg
    assert "coincident" in refused(steps=2)
    # coincident across an edge (a zero-length edge)
    reset()
    pts = trefoil_points(20)
    pts[11] = pts[10]
    obj = edge_mesh(pts, ring_edges(20))
    msg = message_for(obj)
    assert msg and "coincident" in msg, msg
    # not finite
    reset()
    pts = trefoil_points(20)
    pts[4] = (float("nan"), 0.0, 0.0)
    obj = edge_mesh(pts, ring_edges(20))
    msg = message_for(obj)
    assert msg and "finite" in msg, msg
    # shape keys
    reset()
    obj = trefoil(20)
    obj.shape_key_add(name="Basis")
    msg = message_for(obj)
    assert msg and "shape keys" in msg, msg
    # a single edge: the library's message
    reset()
    obj = edge_mesh([(0, 0, 0), (1, 0, 0)], [(0, 1)])
    msg = message_for(obj)
    assert msg is not None
    print("INFO single edge:", msg)
    # wrong object types: the operators are not available
    reset()
    link(bpy.data.objects.new("just_an_empty", None))
    assert not bpy.ops.object.ddg_repulsive_flow.poll()
    # a curve without splines; a single-point spline
    reset()
    curve = bpy.data.curves.new("nothing", 'CURVE')
    obj = link(bpy.data.objects.new("nothing", curve))
    msg = message_for(obj)
    assert msg and "no splines" in msg, msg
    curve.splines.new('POLY')
    msg = message_for(obj)
    assert msg and "at least 2 points" in msg, msg


@test
def a_failed_run_leaves_the_object_untouched():
    pts = trefoil_points(20)
    pts[10] = pts[3]
    obj = edge_mesh(pts, ring_edges(20))
    before = coords(obj)
    refused(steps=3)
    assert coords(obj) == before


@test
def beta_not_above_alpha_is_a_readable_error():
    obj = trefoil(20)
    msg = message_for(obj, alpha=3.0, beta=2.5)
    assert msg is None or "beta" in msg.lower() or "alpha" in msg.lower(), msg
    print("INFO beta < alpha:", msg)


# ----------------------------------------------------------------------------------------------


OPTIONS = ()


def main():
    global ARGS, OPTIONS
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
    OPTIONS = module(".repulsive_operator").OPTION_NAMES

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
