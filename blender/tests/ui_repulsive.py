# SPDX-License-Identifier: MIT
"""Windowed test of the live (modal) Repulsive Curves operator, which the headless tests cannot
reach (no timers and no events in `-b`). Run under a virtual display:

    xvfb-run -a blender --factory-startup --enable-event-simulate \
        --python blender/tests/ui_repulsive.py -- --module bl_ext.user_default.ddg_uv_unwrap

`scripts/blender-test.sh` does that. Checked, with the real event loop and timer:

* the operator starts from the viewport (INVOKE_DEFAULT), returns RUNNING_MODAL, and the object
  changes while it runs (worker thread + timer);
* the sidebar panel draws while it runs (with Accept / Cancel);
* Cancel (button) and Esc (a simulated key press) restore the exact original positions;
* Accept (button), Enter and right click (simulated key presses and mouse button) keep the result;
  accepting is one undo step: undo brings the original back, redo the result;
* the N-steps operator from the viewport: one undo step, redo, Adjust Last Operation;
* the operator ends cleanly each time (no running flow, no leftover backup mesh).

Key events need Blender's ``--enable-event-simulate``; without it those checks are skipped with a
notice and only the button paths run.
"""

import argparse
import array
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


def coords(obj):
    out = array.array('f', bytes(12 * len(obj.data.vertices)))
    obj.data.vertices.foreach_get("co", out)
    return out


def live_module():
    import importlib

    return importlib.import_module(ARGS.module + ".repulsive_operator")


def running():
    return live_module()._LIVE["flow"] is not None


def make_curve():
    import math

    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    n = 40
    pts = []
    for i in range(n):
        t = 2.0 * math.pi * i / n
        pts.append((math.sin(t) + 2.0 * math.sin(2.0 * t), math.cos(t) - 2.0 * math.cos(2.0 * t),
                    -0.3 * math.sin(3.0 * t)))
    mesh = bpy.data.meshes.new("trefoil")
    mesh.from_pydata(pts, [(i, (i + 1) % n) for i in range(n)], [])
    mesh.update()
    obj = bpy.data.objects.new("trefoil", mesh)
    bpy.context.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    for other in bpy.context.view_layer.objects:
        other.select_set(other is obj)
    return obj


def start():
    obj = STATE["obj"] = make_curve()
    with bpy.context.temp_override(**viewport()):
        bpy.ops.ed.undo_push(message="make curve")  # the state undo returns to
    STATE["before"] = coords(obj)
    with bpy.context.temp_override(**viewport()):
        result = bpy.ops.object.ddg_repulsive_live('INVOKE_DEFAULT', True)
    check(result == {'RUNNING_MODAL'}, "the live operator starts and keeps running")
    check(running(), "the live flow is registered as running")
    STATE["waited"] = 0


def wait_for_motion(next_step):
    def step():
        obj = STATE["obj"]
        STATE["waited"] += 1
        flow = live_module()._LIVE["flow"]
        moved = coords(obj) != STATE["before"] and flow is not None and flow.session.iteration >= 3
        if not moved and STATE["waited"] < 200:
            return step
        check(moved, "the object changes while the flow runs (step {})".format(
            flow.session.iteration if flow else None))
        return next_step
    return step


def draw_panel():
    # The sidebar panel draws with the Accept / Cancel buttons while the flow runs.
    with bpy.context.temp_override(**viewport()):
        bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP', iterations=2)
    check(running(), "the sidebar redraw did not stop the flow")


def finish_with_button(accept):
    with bpy.context.temp_override(**viewport()):
        bpy.ops.object.ddg_repulsive_live_finish('EXEC_DEFAULT', accept=accept)


def send(kind, value='PRESS'):
    SIMULATE(type=kind, value=value)


def wait_until_stopped(next_step):
    def step():
        STATE["waited"] += 1
        if running() and STATE["waited"] < 400:
            return step
        check(not running(), "the operator ended")
        return next_step
    return step


# ---- the scenario ------------------------------------------------------------------------------


def step_cancel_button():
    start()
    return wait_for_motion(step_cancel_button_2)


def step_cancel_button_2():
    draw_panel()
    STATE["waited"] = 0
    finish_with_button(False)
    return wait_until_stopped(step_cancel_button_3)


def step_cancel_button_3():
    check(coords(STATE["obj"]) == STATE["before"], "Cancel button: the exact original positions")
    check(not [m for m in bpy.data.meshes if m.name.startswith(".ddg")], "no backup mesh left")
    return step_esc if SIMULATE else step_accept_button


def step_esc():
    start()
    return wait_for_motion(step_esc_2)


def step_esc_2():
    STATE["waited"] = 0
    with bpy.context.temp_override(**viewport()):
        send('ESC')
    return wait_until_stopped(step_esc_3)


def step_esc_3():
    check(coords(STATE["obj"]) == STATE["before"], "Esc: the exact original positions")
    return step_accept_button


def step_accept_button():
    start()
    return wait_for_motion(step_accept_button_2)


def step_accept_button_2():
    STATE["waited"] = 0
    finish_with_button(True)
    return wait_until_stopped(step_accept_button_3)


def step_accept_button_3():
    obj = STATE["obj"]
    after = coords(obj)
    check(after != STATE["before"], "Accept button: the result is kept")
    STATE["after"] = after
    # one undo step
    if bpy.ops.ed.undo.poll():
        with bpy.context.temp_override(**viewport()):
            bpy.ops.ed.undo()
        check(coords(bpy.context.view_layer.objects.active) == STATE["before"],
              "undo restores the original curve (one step)")
        with bpy.context.temp_override(**viewport()):
            bpy.ops.ed.redo()
        check(coords(bpy.context.view_layer.objects.active) == after, "redo brings the result back")
    else:
        log("SKIP: undo is not available")
    return step_enter if SIMULATE else step_n_steps


def step_enter():
    start()
    return wait_for_motion(step_enter_2)


def step_enter_2():
    STATE["waited"] = 0
    with bpy.context.temp_override(**viewport()):
        send('RET')
    return wait_until_stopped(step_enter_3)


def step_enter_3():
    check(coords(STATE["obj"]) != STATE["before"], "Enter: the result is kept")
    return step_rmb


def step_rmb():
    start()
    return wait_for_motion(step_rmb_2)


def step_rmb_2():
    STATE["waited"] = 0
    with bpy.context.temp_override(**viewport()):
        send('RIGHTMOUSE')
    return wait_until_stopped(step_rmb_3)


def step_rmb_3():
    check(coords(STATE["obj"]) != STATE["before"], "Right click: the result is kept")
    return step_n_steps


def step_n_steps():
    """The N-steps operator from the viewport: one undo step, redo, Adjust Last Operation."""
    obj = STATE["obj"] = make_curve()
    before = coords(obj)
    with bpy.context.temp_override(**viewport()):
        bpy.ops.ed.undo_push(message="make curve")
        result = bpy.ops.object.ddg_repulsive_flow('INVOKE_DEFAULT', True)
    check(result == {'FINISHED'}, "the N-steps operator finished from the viewport")
    obj = bpy.context.view_layer.objects.active
    after = coords(obj)
    check(after != before, "the N-steps operator moved the curve")
    with bpy.context.temp_override(**viewport()):
        check(bpy.ops.screen.redo_last.poll(), "Adjust Last Operation is available")
        bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP', iterations=2)
        bpy.ops.ed.undo()
    check(coords(bpy.context.view_layer.objects.active) == before, "undo restores the curve")
    with bpy.context.temp_override(**viewport()):
        bpy.ops.ed.redo()
    check(coords(bpy.context.view_layer.objects.active) != before, "redo brings the result back")
    return None


STEP = [step_cancel_button]


def run():
    global SIMULATE
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
        print("UI REPULSIVE TEST PASSED")
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
