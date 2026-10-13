# SPDX-License-Identifier: MIT
"""Smoke test of the add-on in a windowed Blender (run under a virtual display):

    xvfb-run -a blender --factory-startup --python blender/tests/ui_smoke.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap

`scripts/blender-test-ui.sh` does that after installing the extension into a scratch user
directory. What the headless tests cannot reach is checked here:

* the operator runs from the 3D viewport the way a menu click runs it (INVOKE_DEFAULT, with an
  undo step) and is registered for "Adjust Last Operation";
* its redo panel is drawn (the operator's `draw` is called) without errors, in Object mode and
  in Edit mode;
* the Object menu and the Edit-mode UV menu draw the add-on's entry;
* undo restores the previous UVs, redo brings the new ones back.

The steps run from a timer so that the window exists and can redraw between them.
"""

import argparse
import array
import sys
import traceback

import bpy

parser = argparse.ArgumentParser()
parser.add_argument("--module", required=True)
ARGS = parser.parse_args(sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else [])

CALLS = {"draw": 0, "_draw_menu_entry": 0}
FAILURES = []


def profile(frame, event, _arg):
    # Counts calls of the add-on's UI functions without touching the add-on.
    if event == "call":
        code = frame.f_code
        if code.co_name in CALLS and "ddg_uv_unwrap" in code.co_filename:
            CALLS[code.co_name] += 1


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


def redraw():
    bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP', iterations=2)


def uvs(obj):
    mesh = obj.data
    values = array.array('f', bytes(8 * len(mesh.loops)))
    mesh.uv_layers.active.uv.foreach_get("vector", values)
    return values.tobytes()


def step_object_mode():
    assert not bpy.app.background
    assert bpy.ops.preferences.addon_enable(module=ARGS.module) == {'FINISHED'}
    sys.setprofile(profile)
    cube = bpy.context.view_layer.objects.active
    check(cube is not None and cube.type == 'MESH' and cube.mode == 'OBJECT',
          "the factory scene has an active mesh in Object mode")
    before = uvs(cube)
    with bpy.context.temp_override(**viewport()):
        # The second argument makes this an undo step, as a menu click would be.
        result = bpy.ops.uv.ddg_conformal_unwrap('INVOKE_DEFAULT', True, lambda_=0.3)
        check(result == {'FINISHED'}, "operator finished from the viewport in Object mode")
        last = bpy.context.window_manager.operators[-1] if bpy.context.window_manager.operators \
            else None
        check(last is not None and last.bl_idname == "UV_OT_ddg_conformal_unwrap",
              "operator is registered for Adjust Last Operation")
        check(last is not None and abs(last.properties.lambda_ - 0.3) < 1e-6,
              "its properties are kept for the redo panel")
        after = uvs(cube)
        check(after != before, "the UV map changed")
        redraw()
        check(bpy.ops.screen.redo_last.poll(), "Adjust Last Operation is available")
        bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        redraw()
        check(CALLS["draw"] > 0, "the operator's panel was drawn ({} calls)".format(CALLS["draw"]))

        check(bpy.ops.ed.undo.poll(), "undo is available")
        check(bpy.ops.ed.undo() == {'FINISHED'}, "undo ran")
        cube = bpy.context.view_layer.objects.active
        check(uvs(cube) == before, "undo restored the previous UVs")
        check(bpy.ops.ed.redo() == {'FINISHED'}, "redo ran")
        cube = bpy.context.view_layer.objects.active
        check(uvs(cube) == after, "redo brought the new UVs back")

        # Seams marked by the operator are undone with it.
        check(not any(e.use_seam for e in cube.data.edges), "the cube has no seams to begin with")
        result = bpy.ops.uv.ddg_conformal_unwrap('INVOKE_DEFAULT', True, mark_seams=True)
        cube = bpy.context.view_layer.objects.active
        marked = sum(e.use_seam for e in cube.data.edges)
        check(result == {'FINISHED'} and marked > 0, "cuts were marked as seams ({})".format(marked))
        bpy.ops.ed.undo()
        cube = bpy.context.view_layer.objects.active
        check(not any(e.use_seam for e in cube.data.edges), "undo removed the seams again")
        check(uvs(cube) == after, "and left the UVs of the step before")
    return step_edit_mode


def edit_uvs(obj):
    import bmesh
    edit_mesh = bmesh.from_edit_mesh(obj.data)
    layer = edit_mesh.loops.layers.uv.active
    values = array.array('f')
    for face in edit_mesh.faces:
        for loop in face.loops:
            values.extend(loop[layer].uv)
    return values.tobytes()


def step_edit_mode():
    draws = CALLS["draw"]
    with bpy.context.temp_override(**viewport()):
        # Entering Edit mode as an undo step of its own, as the Tab key does.
        bpy.ops.object.mode_set('EXEC_DEFAULT', True, mode='EDIT')
        cube = bpy.context.view_layer.objects.active
        # Give the UVs other values first (every face on the full square), as an undo step.
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.uv.reset('EXEC_DEFAULT', True)
        before = edit_uvs(cube)
        result = bpy.ops.uv.ddg_conformal_unwrap('INVOKE_DEFAULT', True, use_seams=False,
                                                 lambda_=0.7)
        check(result == {'FINISHED'} and cube.mode == 'EDIT',
              "operator finished from the viewport in Edit mode and stayed there")
        after = edit_uvs(cube)
        check(after != before, "the UV map changed in Edit mode")
        redraw()
        bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        redraw()
        check(CALLS["draw"] > draws, "the operator's panel was drawn in Edit mode")
        check(bpy.ops.ed.undo.poll(), "undo is available in Edit mode")
        check(bpy.ops.ed.undo() == {'FINISHED'}, "undo ran in Edit mode")
        cube = bpy.context.view_layer.objects.active
        check(cube.mode == 'EDIT', "undo stays in Edit mode")
        check(edit_uvs(cube) == before, "undo restored the previous UVs in Edit mode")
        check(bpy.ops.ed.redo() == {'FINISHED'}, "redo ran in Edit mode")
        cube = bpy.context.view_layer.objects.active
        check(cube.mode == 'EDIT' and edit_uvs(cube) == after,
              "redo brought the new UVs back in Edit mode")
    return step_menu_uv


def step_menu_uv():
    with bpy.context.temp_override(**viewport()):
        bpy.ops.wm.call_menu(name="VIEW3D_MT_uv_map")
        redraw()
    check(CALLS["_draw_menu_entry"] > 0, "the Edit-mode UV menu drew the add-on's entry")
    return step_menu_object


def step_menu_object():
    entries = CALLS["_draw_menu_entry"]
    with bpy.context.temp_override(**viewport()):
        bpy.ops.wm.call_menu(name="VIEW3D_MT_object")
        redraw()
    check(CALLS["_draw_menu_entry"] > entries, "the Object menu drew the add-on's entry")
    return step_cut


def step_cut():
    """Start, Flow (redo panel, undo and redo) and Apply of the variational cut from the viewport."""
    with bpy.context.temp_override(**viewport()):
        bpy.ops.object.mode_set('EXEC_DEFAULT', True, mode='OBJECT')
        bpy.ops.mesh.primitive_uv_sphere_add('EXEC_DEFAULT', True, segments=24, ring_count=12)
        sphere = bpy.context.view_layer.objects.active
        result = bpy.ops.uv.ddg_cut_start('INVOKE_DEFAULT', True, mode='CIRCLE', radius=0.9,
                                          vertex=5)
        check(result == {'FINISHED'}, "Cut: Start finished from the viewport")
        names = [a.name for a in sphere.data.attributes if a.name.startswith("ddg_cut_region_")]
        check(len(names) == 2, "Start made the region attributes")
        before = [list(sphere.data.attributes[n].data[3].value for n in [names[0]])]
        draws = CALLS["draw"]
        result = bpy.ops.uv.ddg_cut_flow('INVOKE_DEFAULT', True, steps=3, use_stop_rule=False,
                                         length_mode='EXPLICIT', length_weight=0.1)
        check(result == {'FINISHED'}, "Cut: Flow finished from the viewport")
        after = [list(sphere.data.attributes[n].data[3].value for n in [names[0]])]
        check(after != before, "Flow moved the cut")
        redraw()
        bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        redraw()
        check(CALLS["draw"] > draws, "the Flow panel was drawn")
        check(bpy.ops.ed.undo() == {'FINISHED'}, "undo of Flow ran")
        sphere = bpy.context.view_layer.objects.active
        again = [list(sphere.data.attributes[n].data[3].value for n in [names[0]])]
        check(again == before, "undo restored the cut")
        check(bpy.ops.ed.redo() == {'FINISHED'}, "redo of Flow ran")
        result = bpy.ops.uv.ddg_cut_apply('INVOKE_DEFAULT', True)
        sphere = bpy.context.view_layer.objects.active
        check(result == {'FINISHED'} and sum(e.use_seam for e in sphere.data.edges) > 0,
              "Cut: Apply marked seams")
        check(not [o for o in bpy.data.objects if o.name.endswith(" cut")], "preview removed")
        check(bpy.ops.ed.undo() == {'FINISHED'}, "undo of Apply ran")
        sphere = bpy.context.view_layer.objects.active
        check(not any(e.use_seam for e in sphere.data.edges), "undo removed the seams again")
    return None


STEP = [step_object_mode]


def run():
    try:
        following = STEP[0]()
    except Exception:
        traceback.print_exc(file=sys.stdout)
        FAILURES.append("exception in " + STEP[0].__name__)
        following = None
    if following is not None:
        STEP[0] = following
        return 0.3
    sys.setprofile(None)
    log("calls:", CALLS)
    print("UI RESULT: {} failed".format(len(FAILURES)))
    for failure in FAILURES:
        print("  failed:", failure)
    if not FAILURES:
        print("UI SMOKE TEST PASSED")
    sys.stdout.flush()
    bpy.ops.wm.quit_blender()
    return None


bpy.app.timers.register(run, first_interval=1.0)
