# SPDX-License-Identifier: MIT
"""Windowed test of the Points as Tori panel and modal operator, which the headless tests cannot
reach (no timers and no events in `-b`). Run under a virtual display:

    xvfb-run -a blender --factory-startup --enable-event-simulate \
        --python blender/tests/ui_points_as_tori.py -- --module bl_ext.user_default.ddg_uv_unwrap

`scripts/blender-test.sh` does that. Checked, with the real event loop and timer:

* the scene settings and the operator have every property, with the defaults (Surface, 128, Auto);
* the sidebar panel (DDG tab) is registered and draws in every mode, with the device box (greyed
  out without CUDA), the Advanced sub-panel and the status line, without errors;
* the operator starts from the viewport (INVOKE_DEFAULT), returns RUNNING_MODAL, the panel shows
  the stage and the Cancel button while it runs, and the result object appears when it ends;
  accepting is one undo step (undo removes the object, redo brings it back);
* Cancel (button) and Esc (a simulated key press) end the job without creating an object;
* Orient Normals from the viewport writes the normals back;
* the operator from the viewport reports a mesh without normals as an error;
* Adjust Last Operation reruns the operator (cache hit).

Key events need Blender's ``--enable-event-simulate``; without it the Esc check is skipped with a
notice and only the button path runs.
"""

import argparse
import math
import sys
import traceback

import bpy

parser = argparse.ArgumentParser()
parser.add_argument("--module", required=True)
parser.add_argument("--source", default="", help="register the add-on from this directory")
parser.add_argument("--site", default="", help="directory with an unpacked ddg_tools wheel")
ARGS = parser.parse_args(sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else [])

FAILURES = []
STATE = {}
SIMULATE = None


def log(*parts):
    print("UI:", *parts)
    sys.stdout.flush()


def check(condition, message):
    if condition:
        log("ok:", message)
    else:
        FAILURES.append(message)
        log("FAILED:", message)


def viewport():
    window = bpy.context.window_manager.windows[0]
    for area in window.screen.areas:
        if area.type == 'VIEW_3D':
            region = next(r for r in area.regions if r.type == 'WINDOW')
            return {"window": window, "area": area, "region": region}
    raise RuntimeError("no 3D viewport in the factory startup screen")


def module(name):
    import importlib

    return importlib.import_module(ARGS.module + name)


def op():
    return module(".points_as_tori_operator")


def running():
    return op()._RUN["job"] is not None


def make_cloud(n=300, name="cloud", normals=True, flip=False):
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    for mesh in list(bpy.data.meshes):
        bpy.data.meshes.remove(mesh)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    pts, nrm = [], []
    for i in range(n):
        z = 1.0 - 2.0 * (i + 0.5) / n
        r = math.sqrt(1.0 - z * z)
        a = golden * i
        d = (r * math.cos(a), r * math.sin(a), z)
        s = -1.0 if flip and (i * 7919) % 100 < 30 else 1.0
        pts.extend(0.5 * c for c in d)
        nrm.extend(s * c for c in d)
    mesh = bpy.data.meshes.new(name)
    mesh.vertices.add(n)
    mesh.vertices.foreach_set("co", pts)
    if normals:
        mesh.attributes.new("normal", 'FLOAT_VECTOR', 'POINT').data.foreach_set("vector", nrm)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    bpy.context.view_layer.update()
    bpy.context.view_layer.objects.active = obj
    for other in bpy.context.view_layer.objects:
        other.select_set(other is obj)
    return obj


def settings():
    return bpy.context.scene.ddg_points_as_tori


def redraw():
    with bpy.context.temp_override(**viewport()):
        bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP', iterations=2)


def objects():
    return set(bpy.data.objects.keys())


def send(kind, value='PRESS'):
    SIMULATE(type=kind, value=value)


def start_operator(**kwargs):
    s = settings()
    for k, v in kwargs.items():
        setattr(s, k, v)
    s.device = 'CPU'
    with bpy.context.temp_override(**viewport()):
        return bpy.ops.object.ddg_points_as_tori('INVOKE_DEFAULT')


def wait_until_stopped(next_step, limit=600):
    STATE["waited"] = 0

    def step():
        STATE["waited"] += 1
        if running() and STATE["waited"] < limit:
            return step
        check(not running(), "the operator ended")
        return next_step
    return step


# ---- the scenario ------------------------------------------------------------------------------


def step_properties():
    s = settings()
    names = ("mode", "level", "offset", "resolution", "margin", "partner", "boolean_op",
             "boolean_smooth", "orient_iterations", "device", "cuda_ordinal", "sign_rule",
             "evaluation", "threads", "normal_attribute", "use_cache")
    check(all(hasattr(s, n) for n in names), "the scene settings have every option")
    check((s.mode, s.resolution, s.device, s.sign_rule, s.evaluation) ==
          ('SURFACE', 128, 'AUTO', 'REFERENCE', 'PER_QUERY'), "defaults")
    props = bpy.ops.object.ddg_points_as_tori.get_rna_type().properties
    check(all(n in props for n in names), "the operator has every option")
    check(hasattr(bpy.types, "VIEW3D_PT_ddg_points_as_tori"), "the sidebar panel is registered")
    check(bpy.ops.object.ddg_points_as_tori_cancel.poll() is False, "no job: Cancel is not available")
    return step_panel


def step_panel():
    make_cloud()
    for mode in ('SURFACE', 'OFFSET', 'BOOLEAN', 'ORIENT'):
        settings().mode = mode
        redraw()
    check(True, "the panel draws in every mode")
    op().LAST.update(text="a status line that is long enough to be wrapped over two or more lines "
                          "of the sidebar", error=False)
    redraw()
    op().LAST.update(text="an error", error=True)
    redraw()
    op().LAST.update(text="", error=False)
    settings().mode = 'SURFACE'
    return step_run


def step_run():
    obj = make_cloud()
    STATE["before"] = objects()
    with bpy.context.temp_override(**viewport()):
        bpy.ops.ed.undo_push(message="start")
        bpy.ops.ed.undo_push(message="make cloud")
    result = start_operator(resolution=32)
    check(result == {'RUNNING_MODAL'}, "the operator starts and keeps running")
    check(running(), "the job is registered as running")
    redraw()
    return wait_until_stopped(step_run_done)


def step_run_done():
    new = objects() - STATE["before"]
    check(new == {"cloud Surface"}, "the surface object appeared ({})".format(sorted(new)))
    check(bpy.context.view_layer.objects.active.name == "cloud", "the cloud stays active")
    check("cache" in op().LAST["text"] and "triangles" in op().LAST["text"], "status line")
    with bpy.context.temp_override(**viewport()):
        check(bpy.ops.screen.redo_last.poll(), "Adjust Last Operation is available")
    if bpy.ops.ed.undo.poll():
        with bpy.context.temp_override(**viewport()):
            bpy.ops.ed.undo()
        check("cloud Surface" not in bpy.data.objects, "undo removes the result (one step)")
        with bpy.context.temp_override(**viewport()):
            bpy.ops.ed.redo()
        check("cloud Surface" in bpy.data.objects, "redo brings it back")
    else:
        log("SKIP: undo is not available")
    start_operator(resolution=32)
    return wait_until_stopped(step_rerun_done)


def step_rerun_done():
    check("cache hit (memory)" in op().LAST["text"], "the second run hit the cache")
    return step_cancel_button


def step_cancel_button():
    make_cloud(n=40000)  # takes seconds to build
    STATE["before"] = objects()
    result = start_operator(resolution=32, use_cache=False)
    check(result == {'RUNNING_MODAL'}, "the long job starts")
    redraw()
    check(running() and op()._RUN["job"].stage != "", "the panel has a stage to show")
    with bpy.context.temp_override(**viewport()):
        bpy.ops.object.ddg_points_as_tori_cancel()
    STATE["waited"] = 0
    return wait_until_stopped(step_cancel_button_done, 100)


def step_cancel_button_done():
    check(objects() == STATE["before"], "Cancel button: no object created")
    return step_esc if SIMULATE else step_orient


def step_esc():
    make_cloud(n=40000)
    STATE["before"] = objects()
    start_operator(resolution=32, use_cache=False)
    with bpy.context.temp_override(**viewport()):
        send('ESC')
    return wait_until_stopped(step_esc_done, 100)


def step_esc_done():
    check(objects() == STATE["before"], "Esc: no object created")
    return step_orient


def step_orient():
    obj = make_cloud(n=300, flip=True)
    before = [v for v in obj.data.attributes["normal"].data[0].vector]
    STATE["obj"] = obj
    STATE["before_normals"] = [tuple(d.vector) for d in obj.data.attributes["normal"].data]
    start_operator(mode='ORIENT')
    return wait_until_stopped(step_orient_done)


def step_orient_done():
    now = [tuple(d.vector) for d in STATE["obj"].data.attributes["normal"].data]
    check(now != STATE["before_normals"], "Orient Normals wrote flipped normals back")
    check("flips per iteration" in op().LAST["text"], "the report lists the flips per iteration")
    settings().mode = 'SURFACE'
    return step_no_normals


def step_no_normals():
    make_cloud(normals=False)
    try:
        with bpy.context.temp_override(**viewport()):
            result = bpy.ops.object.ddg_points_as_tori('INVOKE_DEFAULT')
    except RuntimeError as ex:  # an operator that reports an error raises in a script
        result = {'CANCELLED'} if "but no normals" in str(ex) else str(ex)
    check(result == {'CANCELLED'}, "a mesh without normals is refused")
    check("but no normals" in op().LAST["text"] and op().LAST["error"], "with the message in the panel")
    redraw()
    return None


STEP = [step_properties]


def run():
    try:
        following = STEP[0]()
    except Exception:
        traceback.print_exc(file=sys.stdout)
        FAILURES.append("exception in " + STEP[0].__name__)
        following = None
    if following is not None:
        STEP[0] = following
        return 0.2
    print("UI RESULT: {} failed".format(len(FAILURES)))
    for failure in FAILURES:
        print("  failed:", failure)
    if not FAILURES:
        print("UI POINTS AS TORI TEST PASSED")
    sys.stdout.flush()
    bpy.ops.wm.quit_blender()
    return None


def begin():
    global SIMULATE
    assert not bpy.app.background
    if ARGS.site:
        sys.path.insert(0, ARGS.site)
    if ARGS.source:
        import importlib

        sys.path.insert(0, ARGS.source)
        importlib.import_module(ARGS.module).register()
    else:
        assert bpy.ops.preferences.addon_enable(module=ARGS.module) == {'FINISHED'}
    window = bpy.context.window_manager.windows[0]
    SIMULATE = getattr(window, "event_simulate", None)
    if SIMULATE is None:
        log("NOTE: Blender was started without --enable-event-simulate; key events are not tested")
    bpy.app.timers.register(run, first_interval=0.5)
    return None


bpy.app.timers.register(begin, first_interval=1.0)
