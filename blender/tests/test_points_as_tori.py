# SPDX-License-Identifier: MIT
"""Headless tests of the Points as Tori operator. Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/test_points_as_tori.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [-k NAME]

`scripts/blender-test.sh` runs it against the installed extension and its bundled wheel. For
development, `--source DIR --site DIR` register the add-on from a source directory (the parent
of the package) and import `ddg_tools` from an unpacked wheel instead.

Checked: Surface on a sphere cloud (radius, triangle mesh, closed, in the object's frame with
its transform: location, scale and rotation); Offset (radius moved) and a negative Level;
Boolean of two overlapping spheres (union, intersection, difference extents; partner by name and
by selection; partner errors); Orient Normals on a cloud with 30 % flipped normals (flips
written back, flips per iteration reported, undo restores); the cache (memory hit on the second
run, sidecar file next to a saved .blend and its hit after the memory is cleared, miss after a
change of one normal or one position, the sign rule is part of the key, Use Cache off, Clear Cache);
normals from a PLY import (nx, ny, nz), from the vertex normals of a mesh with faces; the Paper
evaluation rule and the Geometric sign rule; the resolution guard; bad input (no normals, a zero
normal, an empty mesh, a missing partner, CUDA in a CPU wheel) gives the message and CANCELLED;
undo and redo of Surface; the job on a worker thread and its cancel flag.

NOT testable headless: the modal operator's event loop (timer, Esc, the panel's Cancel button).
`ui_points_as_tori.py` covers the panel and properties in a windowed Blender.
"""

import argparse
import array
import importlib
import math
import os
import sys
import tempfile
import threading
import time
import traceback

import bmesh
import bpy

TESTS = []


def test(function):
    TESTS.append(function)
    return function


class Skip(Exception):
    pass


# ----------------------------------------------------------------------------------------------
# Helpers


def module(name=""):
    return importlib.import_module(ARGS.module + name)


def op():
    return module(".points_as_tori_operator")


def io():
    return module(".points_as_tori_io")


def reset():
    active = bpy.context.view_layer.objects.active
    if active is not None and active.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    for mesh in list(bpy.data.meshes):
        bpy.data.meshes.remove(mesh)
    settings = bpy.context.scene.ddg_points_as_tori
    for name in op().OPTION_NAMES:
        settings.property_unset(name)
    io().clear_cache(None)
    op().LAST.update(text="", error=False)


def activate(obj, *also):
    view_layer = bpy.context.view_layer
    view_layer.update()
    chosen = (obj,) + also
    for other in view_layer.objects:
        other.select_set(other in chosen)
    view_layer.objects.active = obj


def fibonacci(n, radius=0.5):
    golden = math.pi * (3.0 - math.sqrt(5.0))
    out = []
    for i in range(n):
        z = 1.0 - 2.0 * (i + 0.5) / n
        r = math.sqrt(1.0 - z * z)
        a = golden * i
        out.append((r * math.cos(a), r * math.sin(a), z))
    return [(radius * x, radius * y, radius * z) for x, y, z in out], out


def cloud(name="cloud", n=300, radius=0.5, location=(0.0, 0.0, 0.0), flip_every=0, attribute="normal"):
    """A faceless mesh of ``n`` points on a sphere with outward normals in the point attribute
    ``attribute`` (every ``flip_every``-th and a third of the others, deterministically, flipped
    when ``flip_every`` is not 0). Returns the object."""
    pts, nrm = fibonacci(n, radius)
    mesh = bpy.data.meshes.new(name)
    mesh.vertices.add(n)
    mesh.vertices.foreach_set("co", [c for p in pts for c in p])
    attr = mesh.attributes.new(attribute, 'FLOAT_VECTOR', 'POINT')
    values = []
    for i, d in enumerate(nrm):
        s = -1.0 if flip_every and (i * 7919) % 100 < 30 else 1.0
        values.extend(s * c for c in d)
    attr.data.foreach_set("vector", values)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    obj.location = location
    bpy.context.collection.objects.link(obj)
    bpy.context.view_layer.update()
    activate(obj)
    return obj


def settings():
    return bpy.context.scene.ddg_points_as_tori


def run(**kwargs):
    kwargs.setdefault("device", 'CPU')
    kwargs.setdefault("resolution", 40)
    before = set(bpy.data.objects.keys())
    result = bpy.ops.object.ddg_points_as_tori('EXEC_DEFAULT', **kwargs)
    assert result == {'FINISHED'}, result
    new = [o for o in bpy.data.objects if o.name not in before]
    assert len(new) <= 1, new
    return new[0] if new else None


def refused(**kwargs):
    """The operator refuses: in a background Blender an operator that reports an error raises
    RuntimeError with the report. Returns the report text."""
    kwargs.setdefault("device", 'CPU')
    kwargs.setdefault("resolution", 40)
    try:
        result = bpy.ops.object.ddg_points_as_tori('EXEC_DEFAULT', **kwargs)
    except RuntimeError as ex:
        return str(ex)
    raise AssertionError("expected the operator to refuse, got {}".format(result))


def last_stages():
    return op().LAST["stages"]


def coords(obj):
    out = array.array('f', bytes(12 * len(obj.data.vertices)))
    obj.data.vertices.foreach_get("co", out)
    return out


def rows(flat):
    return [tuple(flat[i:i + 3]) for i in range(0, len(flat), 3)]


def world_points(obj):
    m = obj.matrix_world
    return [tuple(m @ __import__("mathutils").Vector(p)) for p in rows(coords(obj))]


def radii(points, center=(0.0, 0.0, 0.0)):
    return [math.dist(p, center) for p in points]


def x_range(points):
    xs = [p[0] for p in points]
    return min(xs), max(xs)


def normals_of(obj, name="normal"):
    out = array.array('f', bytes(12 * len(obj.data.vertices)))
    obj.data.attributes[name].data.foreach_get("vector", out)
    return out


def assert_close(a, b, tol, what=""):
    assert abs(a - b) <= tol, "{}: {} vs {} (tolerance {})".format(what, a, b, tol)


# ----------------------------------------------------------------------------------------------
# Modes


@test
def surface_of_a_sphere_cloud():
    src = cloud()
    n_objects = len(bpy.data.objects)
    out = run()
    assert len(bpy.data.objects) == n_objects + 1
    assert out is not src and out.name == "cloud Surface", out.name
    assert all(len(p.vertices) == 3 for p in out.data.polygons)
    assert len(out.data.polygons) > 500
    r = radii(world_points(out))
    assert max(abs(x - 0.5) for x in r) < 0.03, (min(r), max(r))
    # closed manifold: every edge has two faces
    bm = bmesh.new()
    bm.from_mesh(out.data)
    assert all(len(e.link_faces) == 2 for e in bm.edges)
    assert len(bm.verts) - len(bm.edges) + len(bm.faces) == 2
    bm.free()
    stages = last_stages()
    assert [s["stage"] for s in stages] == ["build", "mesh"], stages
    assert stages[0]["cache"] == "miss" and stages[0]["device"] == "cpu"
    text = op().LAST["text"]
    assert "cache miss" in text and "triangles" in text and "CPU" in text, text


@test
def the_result_is_in_the_objects_frame():
    from mathutils import Euler
    src = cloud(location=(1.0, 2.0, 3.0))
    src.scale = (2.0, 2.0, 2.0)
    src.rotation_euler = Euler((0.3, 0.2, 0.1))
    bpy.context.view_layer.update()
    out = run()
    assert (out.matrix_world - src.matrix_world).median_scale < 1e-6 or \
        all(abs(a - b) < 1e-6 for ra, rb in zip(out.matrix_world, src.matrix_world)
            for a, b in zip(ra, rb))
    local = radii(rows(coords(out)))
    assert max(abs(x - 0.5) for x in local) < 0.03, "the local mesh is the unscaled sphere"
    world = radii(world_points(out), tuple(src.location))
    assert max(abs(x - 1.0) for x in world) < 0.06, "in the world the sphere has radius 2 * 0.5"


@test
def offset_moves_the_surface():
    cloud()
    out = run(mode='OFFSET', offset=0.1)
    assert out.name == "cloud Offset"
    assert max(abs(x - 0.6) for x in radii(world_points(out))) < 0.03
    out = run(mode='OFFSET', offset=-0.1)
    assert max(abs(x - 0.4) for x in radii(world_points(out))) < 0.03


@test
def level_sets():
    cloud()
    out = run(level=-0.1)
    assert max(abs(x - 0.4) for x in radii(world_points(out))) < 0.03
    out = run(level=0.2, margin=0.0)
    assert max(abs(x - 0.7) for x in radii(world_points(out))) < 0.03, "the box grows with the level"


@test
def boolean_of_two_spheres():
    a = cloud("a", location=(0.0, 0.0, 0.0))
    b = cloud("b", location=(0.5, 0.0, 0.0))
    activate(a, b)
    union = run(mode='BOOLEAN', boolean_op='UNION')
    assert union.name == "a Union b", union.name
    lo, hi = x_range(world_points(union))
    assert_close(lo, -0.5, 0.04, "union min x")
    assert_close(hi, 1.0, 0.04, "union max x")

    activate(a, b)
    inter = run(mode='BOOLEAN', boolean_op='INTERSECT')
    lo, hi = x_range(world_points(inter))
    assert_close(lo, 0.0, 0.04, "intersection min x")
    assert_close(hi, 0.5, 0.04, "intersection max x")

    activate(a, b)
    diff = run(mode='BOOLEAN', boolean_op='SUBTRACT', partner="b")
    assert diff.name == "a Difference b"
    lo, hi = x_range(world_points(diff))
    assert_close(lo, -0.5, 0.04, "difference min x")
    assert_close(hi, 0.25, 0.04, "difference max x (the plane through the circle of intersection)")
    # both tori sets are cached
    activate(a, b)
    run(mode='BOOLEAN', boolean_op='UNION')
    assert [s["cache"] for s in last_stages() if s["stage"] == "build"] == ["memory", "memory"]


@test
def boolean_needs_a_partner():
    a = cloud("a")
    message = refused(mode='BOOLEAN')
    assert "second point cloud" in message, message
    message = refused(mode='BOOLEAN', partner="nope")
    assert "'nope' does not exist" in message, message
    message = refused(mode='BOOLEAN', partner="a")
    assert "itself" in message, message
    plain = bpy.data.meshes.new("plain")
    plain.vertices.add(5)
    obj = bpy.data.objects.new("plain", plain)
    bpy.context.collection.objects.link(obj)
    activate(a)
    message = refused(mode='BOOLEAN', partner="plain")
    assert "'plain' has 5 vertices but no normals" in message, message


@test
def orient_normals_after_30_percent_flipped():
    obj = cloud(n=400, flip_every=1)
    truth = array.array('f', [c for d in fibonacci(400)[1] for c in d])
    before = normals_of(obj)
    wrong = sum(1 for i in range(400) if sum(before[3 * i + k] * truth[3 * i + k] for k in range(3)) < 0)
    assert 90 < wrong < 150, wrong
    n_objects = len(bpy.data.objects)
    run(mode='ORIENT')
    assert len(bpy.data.objects) == n_objects, "no new object"
    after = normals_of(bpy.data.objects["cloud"])
    bad = sum(1 for i in range(400) if sum(after[3 * i + k] * truth[3 * i + k] for k in range(3)) < 0)
    assert bad <= 8, "{} of {} normals still wrong".format(bad, wrong)
    # exact flips: every normal is the old one or its negation
    for i in range(400):
        for k in range(3):
            assert abs(abs(after[3 * i + k]) - abs(before[3 * i + k])) < 1e-6
    text = op().LAST["text"]
    assert "flips per iteration [" in text and "normals flipped in total" in text, text
    assert "attribute 'normal'" in text, text
    assert last_stages()[-1]["stage"] == "orient"


@test
def orient_normals_is_undoable():
    bpy.ops.ed.undo_push(message="start")
    obj = cloud(n=300, flip_every=1)
    bpy.ops.ed.undo_push(message="make cloud")
    if not bpy.ops.ed.undo.poll():
        raise Skip("undo is not available in this session")
    before = normals_of(obj)
    run(mode='ORIENT')
    bpy.ops.ed.undo_push(message="Points as Tori")  # as above
    after = normals_of(bpy.data.objects["cloud"])
    assert after != before
    bpy.ops.ed.undo()
    assert normals_of(bpy.data.objects["cloud"]) == before
    bpy.ops.ed.redo()
    assert normals_of(bpy.data.objects["cloud"]) == after


@test
def surface_is_undoable():
    bpy.ops.ed.undo_push(message="start")
    cloud()
    bpy.ops.ed.undo_push(message="make cloud")
    if not bpy.ops.ed.undo.poll():
        raise Skip("undo is not available in this session")
    run()
    # a script's operator call pushes no undo step in a background Blender (the viewport does:
    # ui_points_as_tori.py checks that); push the step the viewport would
    bpy.ops.ed.undo_push(message="Points as Tori")
    assert "cloud Surface" in bpy.data.objects
    bpy.ops.ed.undo()
    assert "cloud Surface" not in bpy.data.objects
    assert "cloud" in bpy.data.objects, "undo went back to the pushed state"
    bpy.ops.ed.redo()
    assert "cloud Surface" in bpy.data.objects


# ----------------------------------------------------------------------------------------------
# Cache


@test
def second_run_hits_the_cache():
    cloud()
    run()
    assert last_stages()[0]["cache"] == "miss"
    first = last_stages()[0]["seconds"]
    run(level=0.05)
    stage = last_stages()[0]
    assert stage["cache"] == "memory", stage
    assert stage["seconds"] < max(0.2 * first, 0.05), (stage["seconds"], first)
    assert "cache hit (memory)" in op().LAST["text"]


@test
def changing_the_cloud_invalidates_the_cache():
    obj = cloud()
    run()
    obj.data.vertices[0].co.x += 0.01
    run()
    assert last_stages()[0]["cache"] == "miss", "a moved point"
    run()
    assert last_stages()[0]["cache"] == "memory"
    values = list(obj.data.attributes["normal"].data[5].vector)
    obj.data.attributes["normal"].data[5].vector = [-v for v in values]
    run()
    assert last_stages()[0]["cache"] == "miss", "a flipped normal"
    activate(obj)
    run(sign_rule='GEOMETRIC')
    assert last_stages()[0]["cache"] == "miss", "the sign rule is part of the key"
    run(sign_rule='GEOMETRIC', evaluation='PAPER')
    assert last_stages()[0]["cache"] == "memory", "the evaluation rule is not"
    run(use_cache=False)
    assert last_stages()[0]["cache"] == "miss"


# ----------------------------------------------------------------------------------------------
# Input


@test
def normals_from_a_mesh_with_faces():
    bm = bmesh.new()
    bmesh.ops.create_icosphere(bm, subdivisions=3, radius=0.5)
    mesh = bpy.data.meshes.new("ico")
    bm.to_mesh(mesh)
    bm.free()
    obj = bpy.data.objects.new("ico", mesh)
    bpy.context.collection.objects.link(obj)
    activate(obj)
    out = run()
    assert max(abs(x - 0.5) for x in radii(world_points(out))) < 0.04
    c = io().read_cloud(obj)
    assert c.kind == 'vertex'


@test
def normals_from_a_ply_import():
    pts, nrm = fibonacci(300)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "sphere.ply")
        with open(path, "w") as f:
            f.write("ply\nformat ascii 1.0\nelement vertex 300\n")
            for p in ("x", "y", "z", "nx", "ny", "nz"):
                f.write("property float {}\n".format(p))
            f.write("end_header\n")
            for p, d in zip(pts, nrm):
                f.write("{:.7f} {:.7f} {:.7f} {:.7f} {:.7f} {:.7f}\n".format(*p, *d))
        bpy.ops.wm.ply_import(filepath=path)
    obj = bpy.context.active_object
    assert obj is not None and len(obj.data.vertices) == 300
    out = run()
    assert max(abs(x - 0.5) for x in radii(world_points(out))) < 0.03
    assert io().read_cloud(obj).kind == 'attribute'


@test
def normals_from_nx_ny_nz_attributes():
    obj = cloud()
    mesh = obj.data
    normal = array.array('f', normals_of(obj))
    mesh.attributes.remove(mesh.attributes["normal"])
    for k, name in enumerate(("nx", "ny", "nz")):
        attr = mesh.attributes.new(name, 'FLOAT', 'POINT')
        attr.data.foreach_set("value", normal[k::3])
    assert io().read_cloud(obj).kind == 'nxyz'
    out = run()
    assert max(abs(x - 0.5) for x in radii(world_points(out))) < 0.03
    # orientation writes them back
    for i in range(0, 300, 3):
        for name in ("nx", "ny", "nz"):
            mesh.attributes[name].data[i].value *= -1.0
    run(mode='ORIENT')
    assert "attributes nx, ny, nz" in op().LAST["text"]


@test
def a_named_normal_attribute():
    obj = cloud(attribute="my_normals")
    message = refused()
    assert "no normals" in message, message
    out = run(normal_attribute="my_normals")
    assert max(abs(x - 0.5) for x in radii(world_points(out))) < 0.03
    message = refused(normal_attribute="other")
    assert "no point attribute 'other'" in message, message


@test
def no_normals_is_an_error_with_a_message():
    mesh = bpy.data.meshes.new("bare")
    mesh.vertices.add(12345)
    obj = bpy.data.objects.new("bare", mesh)
    bpy.context.collection.objects.link(obj)
    activate(obj)
    message = refused()
    assert "object 'bare' has 12,345 vertices but no normals" in message, message
    assert "attribute named 'normal'" in message, message
    assert len(bpy.data.objects) == 1, "nothing was created"


@test
def zero_normals_and_empty_meshes_are_errors():
    obj = cloud()
    obj.data.attributes["normal"].data[3].vector = (0.0, 0.0, 0.0)
    message = refused()
    assert "1 of 300 normals" in message and "zero" in message, message
    mesh = bpy.data.meshes.new("empty")
    empty = bpy.data.objects.new("empty", mesh)
    bpy.context.collection.objects.link(empty)
    activate(empty)
    assert "has no vertices" in refused()


@test
def wrong_object_type_is_refused():
    bpy.ops.object.empty_add()
    assert not bpy.ops.object.ddg_points_as_tori.poll()


# ----------------------------------------------------------------------------------------------
# Options


@test
def paper_evaluation_and_geometric_sign_rule():
    cloud()
    out = run(evaluation='PAPER', sign_rule='GEOMETRIC')
    assert max(abs(x - 0.5) for x in radii(world_points(out))) < 0.05


@test
def resolution_guard():
    cloud()
    message = refused(resolution=1000)
    assert "lower the resolution" in message and "limit" in message, message
    assert len(bpy.data.objects) == 1


@test
def cuda_in_a_cpu_wheel_says_so():
    info = module(".repulsive_surface_operator").gpu_info()
    if info["cuda_compiled"]:
        raise Skip("this wheel has CUDA")
    cloud()
    message = refused(device='CUDA')
    assert "cuda" in message.lower(), message
    out = run(device='AUTO')  # falls back to the CPU
    assert last_stages()[0]["device"] == "cpu"
    assert out.name == "cloud Surface"


@test
def no_surface_at_the_level_is_an_error():
    cloud()
    message = refused(level=-5.0)
    assert "no surface at level -5" in message, message


@test
def smooth_boolean():
    a = cloud("a")
    b = cloud("b", location=(0.5, 0.0, 0.0))
    activate(a, b)
    out = run(mode='BOOLEAN', boolean_smooth=0.05)
    assert len(out.data.polygons) > 100


# ----------------------------------------------------------------------------------------------
# The worker thread


@test
def job_on_a_worker_thread():
    cloud()
    o = op()
    opts = o.snapshot_options(settings())
    opts.update(device='CPU', resolution=24)
    clouds, directory = o.prepare(bpy.context, opts)
    job = o.Job(clouds, opts, directory)
    thread = threading.Thread(target=job.run)
    thread.start()
    thread.join(120)
    assert job.done and job.error is None and job.progress == 1.0
    assert job.result["faces"] and len(job.stages) == 2


@test
def cancelled_job_reports_cancelled():
    cloud()
    o = op()
    opts = o.snapshot_options(settings())
    opts.update(device='CPU', resolution=24)
    clouds, directory = o.prepare(bpy.context, opts)
    job = o.Job(clouds, opts, directory)
    job.cancel.set()
    job.run()
    assert job.done and isinstance(job.error, o.Cancelled)


@test
def cancel_stops_the_native_build_and_caches_nothing():
    import time
    cloud(n=20000)
    o = op()
    opts = o.snapshot_options(settings())
    opts.update(device='CPU', resolution=24, use_cache=True)
    clouds, directory = o.prepare(bpy.context, opts)
    job = o.Job(clouds, opts, directory)
    if job.native_progress is None:
        return  # a wheel from before the cancel hook
    thread = threading.Thread(target=job.run)
    thread.start()
    deadline = time.time() + 60
    while job.native_progress.stage != "network" and time.time() < deadline and not job.done:
        time.sleep(0.005)
    job.cancel.set()
    thread.join(120)
    assert job.done and isinstance(job.error, o.Cancelled), job.error
    assert job.native_progress.cancelled
    assert io().memory_get((clouds[0].digest, "reference")) is None, "a cancelled build was cached"


# The last test: it saves the .blend, so the later tests would see a file path.


@test
def sidecar_file_next_to_the_blend():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "scene.blend")
        bpy.ops.wm.save_as_mainfile(filepath=path)
        assert io().sidecar_dir() == os.path.join(tmp, "scene.pat")
        cloud()
        run()
        assert last_stages()[0]["cache"] == "miss"
        files = os.listdir(os.path.join(tmp, "scene.pat"))
        assert len(files) == 1 and files[0].endswith("-reference.tori"), files
        io().clear_cache(None)  # forget the memory only
        run()
        assert last_stages()[0]["cache"] == "sidecar", last_stages()[0]
        assert "cache hit (sidecar file)" in op().LAST["text"]
        # the explicit operator removes the files
        assert bpy.ops.object.ddg_points_as_tori_clear_cache() == {'FINISHED'}
        assert not os.path.exists(os.path.join(tmp, "scene.pat"))
        run()
        assert last_stages()[0]["cache"] == "miss"
        # a corrupt sidecar file is ignored with a warning, and rebuilt
        side = os.path.join(tmp, "scene.pat")
        f = os.path.join(side, os.listdir(side)[0])
        with open(f, "wb") as fh:
            fh.write(b"garbage")
        io().clear_cache(None)
        run()
        assert last_stages()[0]["cache"] == "miss"
        assert any("ignored the sidecar" in w for w in last_stages()[0]["warnings"])


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
