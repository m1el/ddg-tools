# SPDX-License-Identifier: MIT
"""Headless tests of the Repulsive Surfaces operators. Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/test_repulsive_surface.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [-k NAME]

`scripts/blender-test.sh` runs it against the installed extension and its bundled wheel. For
development, `--source DIR --site DIR` register the add-on from a source directory (the parent
of the package) and import `ddg_tools` from an unpacked wheel instead.

Checked: N steps of the flow with remeshing off (energy down, positions only: polygons, UV maps and
vertex groups kept) and on (the mesh data is replaced; material index and smooth flag follow the
parent polygon; UVs and other groups dropped with a warning; the pin group kept on the pinned
vertices); `steps = 0` (to the step limit or convergence); pins that never move or vanish;
mesh obstacles from a collection; implicit obstacles and attractors from empties (shape by display
type, custom properties); the world and local space; collision detection switched off for an input
that intersects itself; undo and redo (one step); the live flow (``LiveFlow``: worker thread,
snapshots, accept, cancel restores the exact original, also after remeshing); bad input (an edge
shared by three faces, a mesh without faces, loose vertices, a wrong object type, a missing pin
group or collection) gives the library's located message and CANCELLED, not a crash.

NOT testable headless: the modal operator's event loop (timer events, Enter, right click, Esc, the
sidebar's Accept / Cancel buttons). `ui_repulsive_surface.py` drives them in a windowed Blender.
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
    for collection in list(bpy.data.collections):
        bpy.data.collections.remove(collection)
    settings = bpy.context.scene.ddg_repulsive_surface
    for name in OPTIONS:
        settings.property_unset(name)


def module(name=""):
    return importlib.import_module(ARGS.module + name)


def op():
    return module(".repulsive_surface_operator")


def io():
    return module(".repulsive_surface_io")


def link(obj):
    bpy.context.collection.objects.link(obj)
    view_layer = bpy.context.view_layer
    view_layer.update()
    for other in view_layer.objects:
        other.select_set(other is obj)
    view_layer.objects.active = obj
    return obj


def activate(obj):
    view_layer = bpy.context.view_layer
    for other in view_layer.objects:
        other.select_set(other is obj)
    view_layer.objects.active = obj


def mesh_object(bm, name="surface", location=(0.0, 0.0, 0.0)):
    mesh = bpy.data.meshes.new(name)
    bm.to_mesh(mesh)
    bm.free()
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    obj.location = location
    return link(obj)


def sphere(subdivisions=2, radius=1.0, location=(0.0, 0.0, 0.0), name="sphere"):
    bm = bmesh.new()
    bmesh.ops.create_icosphere(bm, subdivisions=subdivisions, radius=radius)
    return mesh_object(bm, name, location)


def quad_cube(cuts=3, name="cube", location=(0.0, 0.0, 0.0)):
    """A cube of quads (6 * (cuts + 1)^2 polygons) with a UV map that stores the x, y
    coordinates of its corners, and two materials (upper and lower half)."""
    bm = bmesh.new()
    bmesh.ops.create_cube(bm, size=2.0)
    bmesh.ops.subdivide_edges(bm, edges=bm.edges[:], cuts=cuts, use_grid_fill=True)
    uv = bm.loops.layers.uv.new("UVMap")
    for face in bm.faces:
        face.material_index = 0 if face.calc_center_median().z > 0 else 1
        face.smooth = True
        for loop in face.loops:
            loop[uv].uv = (loop.vert.co.x, loop.vert.co.y)
    return mesh_object(bm, name, location)


def coords(obj):
    out = array.array('f', bytes(12 * len(obj.data.vertices)))
    obj.data.vertices.foreach_get("co", out)
    return out


def points(obj):
    c = coords(obj)
    return [tuple(c[i:i + 3]) for i in range(0, len(c), 3)]


def polygons(obj):
    return [tuple(p.vertices) for p in obj.data.polygons]


def uv_values(obj):
    layer = obj.data.uv_layers.get("UVMap")
    if layer is None:
        return None
    out = array.array('f', bytes(8 * len(layer.data)))
    layer.data.foreach_get("uv", out)
    return list(out)


def materials(obj):
    return [p.material_index for p in obj.data.polygons]


def center_z(obj, p):
    cs = coords(obj)
    return sum(cs[3 * v + 2] for v in p.vertices) / len(p.vertices)


def dist(a, b):
    return math.dist(a, b)


def settings():
    return bpy.context.scene.ddg_repulsive_surface


def flow(obj, **options):
    """flow_object with the scene settings changed by ``options``; returns the statistics."""
    s = settings()
    for k, v in options.items():
        setattr(s, k, v)
    return op().flow_object(obj, s, bpy.context.evaluated_depsgraph_get())


def run(**kwargs):
    result = bpy.ops.object.ddg_repulsive_surface_flow('EXEC_DEFAULT', **kwargs)
    assert result == {'FINISHED'}, result


def refused(**kwargs):
    """The N-steps operator refuses: in a background Blender an operator that reports an error
    raises RuntimeError with the report. Returns the report text."""
    try:
        result = bpy.ops.object.ddg_repulsive_surface_flow('EXEC_DEFAULT', **kwargs)
    except RuntimeError as ex:
        return str(ex)
    raise AssertionError("expected the operator to refuse, got {}".format(result))


def message_for(obj, **options):
    """The user message ``make_session`` raises for ``obj`` (None: no error)."""
    s = settings()
    for k, v in options.items():
        setattr(s, k, v)
    try:
        surface, _session = op().make_session(obj, s, bpy.context.evaluated_depsgraph_get())
        surface.discard_backup()
    except io().SurfaceInputError as ex:
        return str(ex)
    return None


def new_collection(name, objects):
    collection = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(collection)
    for obj in objects:
        for c in list(obj.users_collection):
            c.objects.unlink(obj)
        collection.objects.link(obj)
    return collection


def plane_obstacle(z, size=6.0, cuts=8, name="plane"):
    bm = bmesh.new()
    bmesh.ops.create_grid(bm, x_segments=cuts, y_segments=cuts, size=size)
    for v in bm.verts:
        v.co.z = z
    obj = mesh_object(bm, name)
    return obj


def vertex_group(obj, name, indices, weight=1.0):
    group = obj.vertex_groups.new(name=name)
    group.add(list(indices), weight, 'REPLACE')
    return group


def group_members(obj, name):
    index = obj.vertex_groups[name].index
    return sorted(v.index for v in obj.data.vertices if any(g.group == index for g in v.groups))


def barycenter(pts):
    n = len(pts)
    return tuple(sum(p[k] for p in pts) / n for k in range(3))


# ----------------------------------------------------------------------------------------------
# Registration


@test
def registered():
    for name in ("OBJECT_OT_ddg_repulsive_surface_flow", "OBJECT_OT_ddg_repulsive_surface_live",
                 "OBJECT_OT_ddg_repulsive_surface_live_finish"):
        assert hasattr(bpy.types, name), name
    for name in ("VIEW3D_PT_ddg_repulsive_surface", "VIEW3D_PT_ddg_repulsive_surface_advanced"):
        assert hasattr(bpy.types, name), name
    props = bpy.ops.object.ddg_repulsive_surface_flow.get_rna_type().properties
    for name in ("steps", "alpha", "beta", "theta", "remesh", "ccd", "pin_group", "space",
                 "obstacle_collection", "implicit_collection", "metric", "w_willmore"):
        assert name in props, name
    assert settings().remesh is True and settings().ccd is True
    assert (settings().alpha, settings().beta) == (6.0, 12.0)


# ----------------------------------------------------------------------------------------------
# N steps


@test
def steps_without_remeshing_write_positions_only():
    obj = quad_cube()
    mesh = obj.data
    polys, uvs, mats = polygons(obj), uv_values(obj), materials(obj)
    vertex_group(obj, "side", range(10))
    before = points(obj)
    stats = flow(obj, steps=4, remesh=False)
    assert stats["steps"] == 4 and not stats["replaced"], stats
    assert stats["energy_after"] < stats["energy_before"], stats
    assert not stats["warnings"], stats["warnings"]
    assert obj.data is mesh
    assert polygons(obj) == polys, "the polygons (quads) are kept"
    assert uv_values(obj) == uvs, "the UV map is kept"
    assert materials(obj) == mats
    assert group_members(obj, "side") == list(range(10)), "vertex groups are kept"
    after = points(obj)
    assert max(dist(p, q) for p, q in zip(before, after)) > 1e-4, "the vertices moved"
    assert dist(barycenter(before), barycenter(after)) < 1e-3, "the barycenter stays"
    assert not [m for m in bpy.data.meshes if m.name.startswith(".ddg")], "no backup left"


@test
def operator_runs_n_steps():
    obj = sphere()
    before = coords(obj)
    run(steps=3, remesh=False)
    assert coords(obj) != before


@test
def steps_with_remeshing_replace_the_mesh_and_keep_materials():
    obj = quad_cube(cuts=3)
    n_before = len(obj.data.polygons)
    stats = flow(obj, steps=3, remesh=True)
    assert stats["replaced"], stats
    mesh = obj.data
    assert all(len(p.vertices) == 3 for p in mesh.polygons), "triangles"
    assert len(mesh.polygons) != n_before
    assert len(mesh.uv_layers) == 0, "UV maps are dropped"
    warnings = " ".join(stats["warnings"])
    assert "UVMap" in warnings and "replaced" in warnings, stats["warnings"]
    mats = materials(obj)
    assert set(mats) == {0, 1}
    agree = sum(1 for p in mesh.polygons if p.material_index == (0 if center_z(obj, p) > 0 else 1))
    assert agree >= 0.9 * len(mats), "material index follows the parent polygon ({}/{})".format(
        agree, len(mats))
    assert all(p.use_smooth for p in mesh.polygons), "the smooth flag follows the parent polygon"
    assert not [m for m in bpy.data.meshes if m.name.startswith(".ddg")], "no backup left"


@test
def remeshing_drops_other_vertex_groups_and_keeps_the_pin_group():
    obj = sphere(subdivisions=2)
    pins = [i for i, p in enumerate(points(obj)) if p[2] > 0.8]
    assert pins
    vertex_group(obj, "pins", pins, 0.75)
    vertex_group(obj, "other", range(5))
    original = [points(obj)[i] for i in pins]
    stats = flow(obj, steps=3, remesh=True, pin_group="pins")
    assert stats["replaced"]
    warnings = " ".join(stats["warnings"])
    assert "other" in warnings, stats["warnings"]
    members = group_members(obj, "pins")
    assert len(members) == len(pins), "pinned vertices never vanish"
    current = sorted(points(obj)[i] for i in members)
    assert all(dist(a, b) < 1e-5 for a, b in zip(current, sorted(original))), "pins did not move"
    index = obj.vertex_groups["pins"].index
    weights = {g.weight for v in obj.data.vertices for g in v.groups if g.group == index}
    assert all(abs(w - 0.75) < 1e-5 for w in weights), weights
    assert "other" not in obj.vertex_groups, "other groups are dropped"


@test
def pins_stay_without_remeshing():
    obj = sphere()
    pins = [0, 1, 2, 3]
    vertex_group(obj, "pins", pins)
    before = points(obj)
    flow(obj, steps=5, remesh=False, pin_group="pins")
    after = points(obj)
    assert all(dist(before[i], after[i]) < 1e-5 for i in pins)
    assert max(dist(before[i], after[i]) for i in range(len(before)) if i not in pins) > 1e-4


@test
def steps_zero_runs_to_the_limit_or_convergence():
    obj = sphere(subdivisions=1)
    stats = flow(obj, steps=0, max_steps=6, remesh=False)
    assert 1 <= stats["steps"] <= 6, stats
    assert stats["status"] in {"max_iterations", "converged", "stuck"}, stats["status"]


@test
def area_and_volume_constraints_and_penalties():
    obj = sphere()
    def area_volume(o):
        bm = bmesh.new()
        bm.from_mesh(o.data)
        a = sum(f.calc_area() for f in bm.faces)
        v = bm.calc_volume(signed=True)
        bm.free()
        return a, v
    a0, v0 = area_volume(obj)
    flow(obj, steps=4, remesh=False, area_constraint=True, volume_constraint=True,
         area_factor=1.0, volume_factor=1.0)
    a1, v1 = area_volume(obj)
    assert abs(a1 / a0 - 1) < 1e-3 and abs(v1 / v0 - 1) < 1e-3, (a0, a1, v0, v1)
    # a growth schedule: the volume moves towards twice the initial value
    for name in ("area_constraint", "area_factor", "volume_constraint", "volume_factor"):
        settings().property_unset(name)
    obj2 = sphere(name="grow")
    flow(obj2, steps=5, remesh=False, volume_constraint=True, volume_factor=2.0, volume_steps=5)
    assert area_volume(obj2)[1] > 1.3 * v0, area_volume(obj2)[1]
    # penalties of Sec. 7.3 run
    obj3 = sphere(name="pen")
    stats = flow(obj3, steps=2, remesh=False, w_area_deviation=1.0, w_willmore=0.1,
                 w_squared_error=0.1)
    assert stats["steps"] == 2


@test
def metric_and_exponents_options():
    for metric in ('L2', 'H1', 'H2'):
        obj = sphere(subdivisions=1, name=metric)
        flow(obj, steps=1, remesh=False, metric=metric)
    obj = sphere(subdivisions=1, name="exp")
    stats = flow(obj, steps=2, remesh=False, alpha=4.0, beta=8.0, theta=0.25)
    assert stats["steps"] == 2


# ----------------------------------------------------------------------------------------------
# Obstacles


def height(obj):
    zs = [p[2] for p in points(obj)]
    return max(zs) - min(zs)


def lowest(obj):
    return min(p[2] for p in points(obj))


@test
def mesh_obstacles_from_a_collection_are_not_crossed():
    free = sphere(name="free")
    flow(free, steps=4, remesh=False)
    reference = lowest(free)
    reset()
    obj = sphere()
    plane = plane_obstacle(-1.05)
    new_collection("Obstacles", [plane])
    activate(obj)
    flow(obj, steps=4, remesh=False, obstacle_collection="Obstacles", obstacle_weight=5.0)
    z = lowest(obj)
    assert z > -1.05, "the surface stays above the obstacle ({})".format(z)
    assert z > reference + 1e-4, "the obstacle pushes the surface up ({} vs {})".format(
        z, reference)
    # the obstacle is not changed
    assert all(abs(p[2] + 1.05) < 1e-5 for p in points(plane))


@test
def implicit_empties_repel_and_attract():
    base = sphere(name="base")
    flow(base, steps=4, remesh=False)
    reference = max(p[0] for p in points(base))
    reset()
    # a sphere empty (display type Sphere, size from the scale) next to the surface repels it
    obj = sphere()
    empty = bpy.data.objects.new("blocker", None)
    empty.empty_display_type = 'SPHERE'
    empty.empty_display_size = 0.5
    empty.location = (1.5, 0.0, 0.0)
    bpy.context.collection.objects.link(empty)
    new_collection("Empties", [empty])
    activate(obj)
    flow(obj, steps=4, remesh=False, implicit_collection="Empties", implicit_weight=2.0)
    xmax = max(p[0] for p in points(obj))
    assert xmax < reference - 1e-4, "the empty pushes the surface away ({} vs {})".format(
        xmax, reference)
    for p in points(obj):
        assert dist(p, (1.5, 0.0, 0.0)) > 0.5, "the surface stays outside the sphere"
    # the empty's scale gives the radius, a custom property switches it to an attractor
    mode = io().implicit_penalty(empty, _FakeSurface(), mode='REPEL', power=2.0, weight=1.0,
                                 minor_ratio=0.25)
    assert mode["type"] == "implicit_obstacle" and abs(mode["surface"]["radius"] - 0.5) < 1e-6
    empty.scale = (2.0, 2.0, 2.0)
    empty["ddg_mode"] = "attract"
    empty["ddg_power"] = 3.0
    empty["ddg_weight"] = 0.25
    bpy.context.view_layer.update()
    p = io().implicit_penalty(empty, _FakeSurface(), mode='REPEL', power=2.0, weight=1.0,
                              minor_ratio=0.25)
    assert p["type"] == "implicit_attractor" and p["power"] == 3.0 and p["weight"] == 0.25, p
    assert abs(p["surface"]["radius"] - 1.0) < 1e-6, p["surface"]


class _FakeSurface:
    world_to_space = (1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)


@test
def empty_display_types_choose_the_shape():
    shapes = {}
    for display, expected in (('SPHERE', "sphere"), ('CONE', "cylinder"), ('CIRCLE', "torus"),
                              ('PLAIN_AXES', "plane"), ('ARROWS', "plane"), ('CUBE', "plane")):
        empty = bpy.data.objects.new("e_" + display, None)
        empty.empty_display_type = display
        empty.empty_display_size = 2.0
        empty.location = (1.0, 2.0, 3.0)
        empty.rotation_euler = (math.pi / 2, 0.0, 0.0)
        bpy.context.collection.objects.link(empty)
        bpy.context.view_layer.update()
        p = io().implicit_penalty(empty, _FakeSurface(), mode='REPEL', power=2.0, weight=1.0,
                                  minor_ratio=0.25)
        shapes[display] = p["surface"]
        assert p["surface"]["type"] == expected, (display, p)
    assert shapes['SPHERE']["radius"] == 2.0 and shapes['SPHERE']["center"] == (1.0, 2.0, 3.0)
    axis = shapes['CONE']["axis"]  # local Z rotated by 90 degrees about X: -Y
    assert abs(axis[1] + 1.0) < 1e-6 and abs(axis[2]) < 1e-6, axis
    assert abs(shapes['CIRCLE']["major"] - 2.0) < 1e-6 and abs(shapes['CIRCLE']["minor"] - 0.5) < 1e-6
    assert shapes['PLAIN_AXES']["point"] == (1.0, 2.0, 3.0)


@test
def attractor_plane_flattens_the_surface():
    base = sphere(name="base")
    flow(base, steps=4, remesh=False)
    reference = height(base)
    reset()
    obj = sphere()
    empty = bpy.data.objects.new("floor", None)
    empty.empty_display_type = 'PLAIN_AXES'
    empty.location = (0.0, 0.0, -2.0)
    empty["ddg_mode"] = "attract"
    empty["ddg_weight"] = 1.0
    bpy.context.collection.objects.link(empty)
    new_collection("Empties", [empty])
    activate(obj)
    flow(obj, steps=4, remesh=False, implicit_collection="Empties")
    # all vertices are pulled to the plane while the barycenter stays: the surface flattens
    assert height(obj) < 0.7 * reference, (height(obj), reference)


# ----------------------------------------------------------------------------------------------
# Space, transform


@test
def object_transform_world_and_local_space():
    obj = sphere(location=(5.0, -3.0, 2.0))
    obj.scale = (1.5, 1.5, 1.5)
    bpy.context.view_layer.update()
    world0 = [tuple(obj.matrix_world @ Vector(p)) for p in points(obj)]
    flow(obj, steps=3, remesh=False, space='WORLD')
    world1 = [tuple(obj.matrix_world @ Vector(p)) for p in points(obj)]
    assert max(dist(a, b) for a, b in zip(world0, world1)) > 1e-4
    assert dist(barycenter(world0), barycenter(world1)) < 1e-3, "the world barycenter stays"
    reset()
    obj = sphere(location=(5.0, -3.0, 2.0))
    obj.scale = (1.5, 1.5, 1.5)
    bpy.context.view_layer.update()
    local0 = points(obj)
    flow(obj, steps=3, remesh=False, space='LOCAL')
    local1 = points(obj)
    assert max(dist(a, b) for a, b in zip(local0, local1)) > 1e-4
    assert dist(barycenter(local0), barycenter(local1)) < 1e-3, "the local barycenter stays"


@test
def self_intersecting_input_turns_ccd_off_with_a_warning():
    bm = bmesh.new()
    bmesh.ops.create_icosphere(bm, subdivisions=2, radius=1.0)
    other = bmesh.new()
    bmesh.ops.create_icosphere(other, subdivisions=2, radius=1.0)
    for v in other.verts:
        v.co.x += 0.5
    mesh = bpy.data.meshes.new("two")
    bm.to_mesh(mesh)
    bm.free()
    extra = bpy.data.meshes.new("extra")
    other.to_mesh(extra)
    other.free()
    mesh.vertices.add(0)
    # join the two spheres into one mesh
    joined = bmesh.new()
    joined.from_mesh(mesh)
    joined.from_mesh(extra)
    joined.to_mesh(mesh)
    joined.free()
    obj = link(bpy.data.objects.new("two", mesh))
    surface = io().read_object(obj, bpy.context.evaluated_depsgraph_get(), space='LOCAL')
    try:
        assert io().intersecting_pairs(surface.positions, surface.faces), "found the intersection"
    finally:
        surface.discard_backup()
    stats = flow(obj, steps=1, remesh=False, ccd=True, barycenter='ON')
    assert any("intersects itself" in w for w in stats["warnings"]), stats["warnings"]
    reset()
    clean = sphere()
    surface = io().read_object(clean, bpy.context.evaluated_depsgraph_get(), space='LOCAL')
    try:
        assert not io().intersecting_pairs(surface.positions, surface.faces), "no false positive"
    finally:
        surface.discard_backup()


# ----------------------------------------------------------------------------------------------
# Undo


@test
def undo_and_redo_are_one_step():
    if not bpy.ops.ed.undo.poll():
        raise Skip("undo is not available in this session")
    obj = sphere()
    bpy.ops.ed.undo_push(message="make sphere")
    before = coords(obj)
    run(steps=3, remesh=False)
    after = coords(bpy.context.active_object)
    assert after != before
    bpy.ops.ed.undo()
    assert coords(bpy.context.active_object) == before, "one undo undoes all steps"
    bpy.ops.ed.redo()
    assert coords(bpy.context.active_object) == after


@test
def undo_after_remeshing_restores_the_original_mesh():
    if not bpy.ops.ed.undo.poll():
        raise Skip("undo is not available in this session")
    obj = quad_cube()
    bpy.ops.ed.undo_push(message="make cube")
    n = len(obj.data.polygons)
    uvs = uv_values(obj)
    run(steps=2, remesh=True)
    assert len(bpy.context.active_object.data.polygons) != n
    bpy.ops.ed.undo()
    restored = bpy.context.active_object
    assert len(restored.data.polygons) == n and uv_values(restored) == uvs


# ----------------------------------------------------------------------------------------------
# The live flow


def live(obj, **options):
    s = settings()
    for k, v in options.items():
        setattr(s, k, v)
    return op().LiveFlow.create(obj, s, bpy.context.evaluated_depsgraph_get())


def wait_for(flow, iterations, timeout=60.0):
    t = time.time()
    while flow.session.iteration < iterations and flow.running and time.time() - t < timeout:
        time.sleep(0.02)


@test
def live_flow_accept_keeps_the_result():
    obj = sphere()
    before = coords(obj)
    flow_ = live(obj, remesh=False)
    flow_.start()
    wait_for(flow_, 3)
    assert flow_.poll(), "a new step was written"
    mid = coords(obj)
    assert mid != before
    stats = flow_.finish()
    assert not flow_.running
    assert stats["steps"] >= 3 and stats["energy_after"] < stats["energy_before"], stats
    assert coords(obj) != before
    assert not [m for m in bpy.data.meshes if m.name.startswith(".ddg")]
    assert flow_.error is None


@test
def live_flow_cancel_restores_the_exact_original():
    obj = sphere()
    before = coords(obj)
    flow_ = live(obj, remesh=False)
    flow_.start()
    wait_for(flow_, 3)
    flow_.poll()
    assert coords(obj) != before
    flow_.cancel()
    assert coords(obj) == before, "the exact original positions"
    assert not [m for m in bpy.data.meshes if m.name.startswith(".ddg")]


@test
def live_flow_cancel_after_remeshing_restores_polygons_and_uvs():
    obj = quad_cube()
    mesh = obj.data
    polys, uvs, mats, before = polygons(obj), uv_values(obj), materials(obj), coords(obj)
    vertex_group(obj, "side", range(7))
    flow_ = live(obj, remesh=True)
    flow_.start()
    wait_for(flow_, 3)
    flow_.poll()
    assert len(obj.data.polygons) != len(polys), "the mesh data was replaced while running"
    flow_.cancel()
    assert obj.data is mesh
    assert polygons(obj) == polys and uv_values(obj) == uvs and materials(obj) == mats
    assert coords(obj) == before
    assert group_members(obj, "side") == list(range(7)), "vertex groups come back"


@test
def live_flow_accept_with_pins_keeps_the_pin_group():
    obj = sphere()
    pins = [i for i, p in enumerate(points(obj)) if p[2] > 0.8]
    vertex_group(obj, "pins", pins)
    flow_ = live(obj, remesh=True, pin_group="pins")
    flow_.start()
    wait_for(flow_, 3)
    stats = flow_.finish()
    assert stats["replaced"]
    assert len(group_members(obj, "pins")) == len(pins)


# ----------------------------------------------------------------------------------------------
# Bad input


@test
def non_manifold_input_gives_the_librarys_located_message():
    mesh = bpy.data.meshes.new("fin")
    mesh.from_pydata([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)],
                     [], [(0, 1, 2), (0, 1, 3), (0, 1, 4)])
    mesh.update()
    link(bpy.data.objects.new("fin", mesh))
    text = refused(steps=1)
    assert "shared by 3 faces" in text or "3 faces" in text, text
    assert "fin" in text


@test
def loose_vertices_and_missing_things_are_refused():
    mesh = bpy.data.meshes.new("loose")
    mesh.from_pydata([(0, 0, 0), (1, 0, 0), (0, 1, 0), (5, 5, 5)], [], [(0, 1, 2)])
    mesh.update()
    obj = link(bpy.data.objects.new("loose", mesh))
    text = message_for(obj)
    assert text and "3" in text, text  # names the vertex
    reset()
    empty_mesh = bpy.data.meshes.new("edges")
    empty_mesh.from_pydata([(0, 0, 0), (1, 0, 0)], [(0, 1)], [])
    obj = link(bpy.data.objects.new("edges", empty_mesh))
    assert "no faces" in refused(steps=1)
    reset()
    obj = sphere()
    assert "no vertex group 'nope'" in message_for(obj, pin_group="nope")
    settings().property_unset("pin_group")
    assert "collection 'nowhere' does not exist" in message_for(obj, obstacle_collection="nowhere")
    settings().property_unset("obstacle_collection")
    cam = bpy.data.objects.new("empty", None)
    bpy.context.collection.objects.link(cam)
    bpy.context.view_layer.objects.active = cam
    assert not bpy.ops.object.ddg_repulsive_surface_flow.poll()
    assert not bpy.ops.object.ddg_repulsive_surface_live.poll()


@test
def bad_options_and_modifiers():
    obj = sphere()
    text = message_for(obj, alpha=6.0, beta=7.0)  # beta - alpha too small
    assert text is not None and "sphere" in text, text
    settings().property_unset("alpha")
    settings().property_unset("beta")
    # a modifier that changes the topology: the evaluated mesh is used and baked
    obj.modifiers.new("sub", 'SUBSURF').levels = 0
    obj.modifiers["sub"].render_levels = 0
    obj.modifiers.new("tri", 'TRIANGULATE')
    bpy.context.view_layer.update()
    stats = flow(obj, steps=1, remesh=False)
    assert any("modifiers" in w for w in stats["warnings"]) or not stats["replaced"], stats
    # shape keys are refused
    reset()
    obj = sphere()
    obj.shape_key_add(name="Basis")
    assert "shape keys" in message_for(obj)


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
    OPTIONS = module(".repulsive_surface_operator").OPTION_NAMES

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
