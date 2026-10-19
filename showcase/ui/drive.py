"""Drive a windowed Blender (under Xvfb) through the add-on and capture the screen.

    DISPLAY=:N blender --factory-startup blank.blend --python showcase/ui/drive.py -- <scene> <outdir>

Scenes: unwrap, cut, sphere, developable. Every step runs from a timer so the window redraws in between;
captures are taken from the X display with ImageMagick's `import` (Blender's own screenshot
operator returns a black image on llvmpipe). Captures: <outdir>/<scene>_<nn>_<tag>.png.
"""

import math
import subprocess
import sys
import traceback
from pathlib import Path

import bpy

ARGS = sys.argv[sys.argv.index("--") + 1:]
SCENE, OUT = ARGS[0], Path(ARGS[1])
OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path(__file__).resolve().parents[2]
SPOT = ROOT / "data/models/spot/spot/spot_quadrangulated.obj"
MOD = "bl_ext.user_default.ddg_uv_unwrap"
N = [0]


def log(*a):
    print("DRIVE:", *a, flush=True)


def win():
    return bpy.context.window_manager.windows[0]


def area(kind):
    w = win()
    for a in w.screen.areas:
        if a.type == kind:
            r = next(r for r in a.regions if r.type == 'WINDOW')
            return {"window": w, "screen": w.screen, "area": a, "region": r}
    raise RuntimeError("no area " + kind)


def redraw():
    bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP', iterations=3)


def capture(tag):
    redraw()
    path = OUT / "{}_{:02d}_{}.png".format(SCENE, N[0], tag)
    N[0] += 1
    subprocess.run(["import", "-window", "root", str(path)], check=True)
    log("captured", path)


def theme():
    t = bpy.context.preferences.themes[0].view_3d
    t.seam = (1.0, 0.12, 0.05)
    t.edge_width = 2
    t.wire_edit = (0.22, 0.22, 0.24)
    t.wire = (1.0, 0.72, 0.0)
    bpy.context.preferences.view.show_splash = False
    bpy.context.preferences.view.ui_line_width = 'THICK'


def set_workspace(name):
    win().workspace = bpy.data.workspaces[name]


def load_spot():
    for o in list(bpy.data.objects):
        bpy.data.objects.remove(o, do_unlink=True)
    bpy.ops.wm.obj_import(filepath=str(SPOT))
    obj = bpy.context.selected_objects[0]
    obj.name = "spot"
    bpy.context.view_layer.objects.active = obj
    me = obj.data
    while me.uv_layers:
        me.uv_layers.remove(me.uv_layers[0])
    me.materials.clear()
    for p in me.polygons:
        p.use_smooth = True
    return obj


def checker_material(obj):
    img = bpy.data.images.new("checker", 1024, 1024)
    img.generated_type = 'UV_GRID'
    mat = bpy.data.materials.new("checker")
    mat.use_nodes = True
    nt = mat.node_tree
    tex = nt.nodes.new("ShaderNodeTexImage")
    tex.image = img
    nt.links.new(tex.outputs["Color"], nt.nodes["Principled BSDF"].inputs["Base Color"])
    nt.nodes.active = tex
    obj.data.materials.append(mat)
    return img


def view3d_setup(texture=True, angle=(-35, 15), zoom_steps=0):
    for a in win().screen.areas:
        if a.type != 'VIEW_3D':
            continue
        sp = a.spaces.active
        sp.shading.type = 'SOLID'
        sp.shading.color_type = 'TEXTURE' if texture else 'MATERIAL'
        sp.shading.light = 'STUDIO'
        sp.overlay.show_floor = False
        sp.overlay.show_axis_x = sp.overlay.show_axis_y = False
        sp.show_gizmo_navigate = False
        r3d = sp.region_3d
        r3d.view_perspective = 'PERSP'
        from mathutils import Euler
        az, el = map(math.radians, angle)
        r3d.view_rotation = Euler((math.pi / 2 - el, 0, az)).to_quaternion()
        ctx = {"window": win(), "screen": win().screen, "area": a,
               "region": next(r for r in a.regions if r.type == 'WINDOW')}
        with bpy.context.temp_override(**ctx):
            bpy.ops.view3d.view_selected()
            for _ in range(zoom_steps):
                bpy.ops.view3d.zoom(delta=1)


def uv_editor_setup(img):
    for a in win().screen.areas:
        if a.type == 'IMAGE_EDITOR':
            a.spaces.active.image = img
            ctx = {"window": win(), "screen": win().screen, "area": a,
                   "region": next(r for r in a.regions if r.type == 'WINDOW')}
            with bpy.context.temp_override(**ctx):
                bpy.ops.image.view_all(fit_view=True)
                bpy.ops.image.view_zoom_out(location=(0.5, 0.5))


# ----------------------------------------------------------------------------- scenes
def edit_all_deselected():
    obj = bpy.data.objects["spot"]
    bpy.context.view_layer.objects.active = obj
    with bpy.context.temp_override(**area('VIEW_3D')):
        if obj.mode != 'EDIT':
            bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
    bpy.context.scene.tool_settings.use_uv_select_sync = True
    bpy.context.scene.tool_settings.mesh_select_mode = (False, True, False)
    for a in win().screen.areas:
        if a.type in ('VIEW_3D', 'IMAGE_EDITOR'):
            a.spaces.active.show_region_toolbar = False


def deselect():
    with bpy.context.temp_override(**area('VIEW_3D')):
        bpy.ops.mesh.select_all(action='DESELECT')


def mouse(x, y):
    win().event_simulate('MOUSEMOVE', 'NOTHING', x=x, y=y)


def key(k):
    win().event_simulate(k, 'PRESS')
    win().event_simulate(k, 'RELEASE')


def scene_unwrap():
    state = {}

    def s0():
        bpy.ops.preferences.addon_enable(module=MOD)
        theme()
        obj = load_spot()
        state["img"] = checker_material(obj)
        set_workspace("UV Editing")
        return 1.0

    def s1():
        edit_all_deselected()
        view3d_setup(texture=True, angle=(-50, 12))
        deselect()
        uv_editor_setup(state["img"])
        return 1.0

    def s2():
        mouse(1900, 20)
        with bpy.context.temp_override(**area('VIEW_3D')):
            res = bpy.ops.uv.ddg_conformal_unwrap('INVOKE_DEFAULT', True, mark_seams=True)
        log("unwrap", res)
        uv_editor_setup(state["img"])
        return 2.0

    def s3():
        capture("result")
        for a in win().screen.areas:
            if a.type == 'IMAGE_EDITOR':
                a.spaces.active.uv_editor.show_stretch = True
                a.spaces.active.uv_editor.display_stretch_type = 'ANGLE'
        return 1.0

    def s4():
        capture("stretch")
        for a in win().screen.areas:
            if a.type == 'IMAGE_EDITOR':
                a.spaces.active.uv_editor.show_stretch = False
        mouse(1180, 700)
        return 0.3

    def s5():
        with bpy.context.temp_override(**area('VIEW_3D')):
            bpy.ops.wm.call_menu(name="VIEW3D_MT_uv_map")
        return 0.5

    def s6():
        mouse(1240, 700 - 336)
        return 2.0

    def s7():
        capture("menu")
        key('ESC')
        return 0.5

    def s8():
        mouse(1180, 700)
        with bpy.context.temp_override(**area('VIEW_3D')):
            bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        return 1.5

    def s9():
        capture("options")
        key('ESC')
        return None

    return [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9]


def scene_cut():
    state = {}
    regions = int(ARGS[2]) if len(ARGS) > 2 else 2

    def s0():
        bpy.ops.preferences.addon_enable(module=MOD)
        theme()
        obj = load_spot()
        state["img"] = checker_material(obj)
        set_workspace("Layout")
        return 1.0

    def s1():
        view3d_setup(texture=False, angle=(-140 if regions == 2 else -50, 18))
        for a in win().screen.areas:
            if a.type == 'VIEW_3D':
                a.spaces.active.show_region_toolbar = False
                a.spaces.active.shading.color_type = 'SINGLE'
                a.spaces.active.shading.single_color = (0.85, 0.85, 0.87)
        mouse(1900, 20)
        with bpy.context.temp_override(**area('VIEW_3D')):
            res = bpy.ops.uv.ddg_cut_start('INVOKE_DEFAULT', True, mode='CLUSTERS', regions=regions)
        log("start", res)
        state["k"] = 0
        return 1.0

    def shot_flow():
        capture("flow{:02d}".format(state["k"]))
        return 0.1

    def flow():
        with bpy.context.temp_override(**area('VIEW_3D')):
            res = bpy.ops.uv.ddg_cut_flow("INVOKE_DEFAULT", True, steps=int(ARGS[4]) if len(ARGS) > 4 else 5, use_stop_rule=False)
        log("flow", state["k"], res, bpy.data.objects["spot"].data.get("ddg_cut_last_flow"))
        state["k"] += 1
        return 0.5

    def options():
        mouse(1180, 700)
        with bpy.context.temp_override(**area('VIEW_3D')):
            bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        return 1.5

    def shot_options():
        capture("options")
        key('ESC')
        return 0.5

    def apply():
        mouse(1900, 20)
        with bpy.context.temp_override(**area('VIEW_3D')):
            res = bpy.ops.uv.ddg_cut_apply('INVOKE_DEFAULT', True, unwrap_after=True)
        log("apply", res)
        set_workspace("UV Editing")
        return 1.0

    def show():
        edit_all_deselected()
        view3d_setup(texture=True, angle=(-140 if regions == 2 else -50, 18))
        deselect()
        uv_editor_setup(state["img"])
        return 1.0

    def shot_applied():
        capture("applied")
        return None

    steps = [s0, s1]
    for _ in range(int(ARGS[3]) if len(ARGS) > 3 else 12):
        steps += [shot_flow, flow]
    return steps + [shot_flow, options, shot_options, apply, show, shot_applied]


def scene_developable():
    """Make Developable (Layout workspace, with its Adjust Last Operation panel), Show
    Developability (the energy as vertex colors) and Developable Seams + UV (UV Editing
    workspace, seams in Edit mode)."""
    state = {}

    def s0():
        bpy.ops.preferences.addon_enable(module=MOD)
        theme()
        obj = load_spot()
        state["img"] = checker_material(obj)
        set_workspace("Layout")
        return 1.0

    def s1():
        view3d_setup(texture=False, angle=(-140, 16))
        for a in win().screen.areas:
            if a.type == 'VIEW_3D':
                a.spaces.active.show_region_toolbar = False
                a.spaces.active.shading.color_type = 'SINGLE'
                a.spaces.active.shading.single_color = (0.85, 0.85, 0.87)
        capture("input")
        mouse(1900, 20)
        return 0.5

    def s2():
        with bpy.context.temp_override(**area('VIEW_3D')):
            res = bpy.ops.object.ddg_make_developable('INVOKE_DEFAULT', True, subdivisions=0,
                                                      iterations=int(ARGS[2]) if len(ARGS) > 2 else 1500)
        log("make developable", res, bpy.data.objects["spot"].data.get("ddg_developable_last"))
        return 1.5

    def s3():
        capture("result")
        mouse(1180, 700)
        with bpy.context.temp_override(**area('VIEW_3D')):
            bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        return 1.5

    def s4():
        capture("options")
        key('ESC')
        return 0.5

    def s5():
        mouse(1900, 20)
        with bpy.context.temp_override(**area('VIEW_3D')):
            res = bpy.ops.object.ddg_show_developability('INVOKE_DEFAULT', True)
        log("show developability", res)
        me = bpy.data.objects["spot"].data
        log("color attrs", [(c.name, c.domain, c.data_type) for c in me.color_attributes],
            "active", me.color_attributes.active_color_name)
        c0 = me.color_attributes.get("developability")
        if c0 is not None:
            vals = [tuple(d.color) for d in c0.data[:2000]]
            log("color sample min/max", min(vals), max(vals))
            # The operator leaves the float attribute 'developability_energy' active, which the
            # viewport draws black; make the color attribute the active one.
            me.color_attributes.active_color = c0
        for a in win().screen.areas:
            if a.type == 'VIEW_3D':
                a.spaces.active.shading.color_type = 'VERTEX'
        return 1.5

    def s6():
        capture("energy")
        mouse(1180, 700)
        with bpy.context.temp_override(**area('VIEW_3D')):
            bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        return 1.5

    def s7():
        capture("energy_options")
        key('ESC')
        return 0.5

    def s8():
        set_workspace("UV Editing")
        return 1.0

    def s9():
        edit_all_deselected()
        view3d_setup(texture=True, angle=(-140, 16))
        deselect()
        for a in win().screen.areas:
            if a.type == 'VIEW_3D':
                a.spaces.active.shading.color_type = 'TEXTURE'
        uv_editor_setup(state["img"])
        return 1.0

    def s10():
        mouse(1900, 20)
        with bpy.context.temp_override(**area('VIEW_3D')):
            res = bpy.ops.uv.ddg_developable_seams('INVOKE_DEFAULT', True)
        log("seams", res)
        uv_editor_setup(state["img"])
        return 2.0

    def s11():
        capture("seams")
        mouse(1180, 700)
        with bpy.context.temp_override(**area('VIEW_3D')):
            bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        return 1.5

    def s12():
        capture("seams_options")
        key('ESC')
        return None

    return [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12]


def scene_sphere():
    state = {}

    def s0():
        bpy.ops.preferences.addon_enable(module=MOD)
        theme()
        obj = load_spot()
        set_workspace("Layout")
        # Spot's own texture through its OBJ UVs: re-import keeps them, so reload.
        for o in list(bpy.data.objects):
            bpy.data.objects.remove(o, do_unlink=True)
        bpy.ops.wm.obj_import(filepath=str(SPOT))
        obj = bpy.context.selected_objects[0]
        obj.name = "spot"
        bpy.context.view_layer.objects.active = obj
        for p in obj.data.polygons:
            p.use_smooth = True
        img = bpy.data.images.load(str(SPOT.parent / "spot_texture.png"))
        mat = obj.data.materials[0] if obj.data.materials else bpy.data.materials.new("spot")
        if not obj.data.materials:
            obj.data.materials.append(mat)
        mat.use_nodes = True
        nt = mat.node_tree
        tex = nt.nodes.new("ShaderNodeTexImage")
        tex.image = img
        bsdf = next(n for n in nt.nodes if n.type == 'BSDF_PRINCIPLED')
        nt.links.new(tex.outputs["Color"], bsdf.inputs["Base Color"])
        nt.nodes.active = tex
        return 1.0

    def s1():
        view3d_setup(texture=True, angle=(-140, 12))
        for a in win().screen.areas:
            if a.type == 'VIEW_3D':
                a.spaces.active.show_region_toolbar = False
        mouse(1900, 20)
        with bpy.context.temp_override(**area('VIEW_3D')):
            res = bpy.ops.object.ddg_map_to_sphere('INVOKE_DEFAULT', True, shape_key=True)
        log("sphere", res)
        state["k"] = 0
        return 1.0

    def options():
        mouse(1180, 700)
        with bpy.context.temp_override(**area('VIEW_3D')):
            bpy.ops.screen.redo_last('INVOKE_DEFAULT')
        return 1.5

    def shot_options():
        capture("options")
        key('ESC')
        if ALIGN:
            align_shape_key(bpy.data.objects["spot"], "Conformal Sphere")
        return 0.5

    def morph():
        obj = bpy.data.objects["spot"]
        kb = obj.data.shape_keys.key_blocks["Conformal Sphere"]
        t = min(state["k"], 30) / 30
        kb.value = t * t * (3 - 2 * t)
        state["k"] += 1
        return 0.05

    def shot():
        capture("morph{:02d}".format(state["k"]))
        return 0.05

    steps = [s0, s1, options, shot_options]
    for _ in range(40):
        steps += [morph, shot]
    return steps


def align_shape_key(obj, name):
    """Rotate a shape key about its centroid so it best matches the basis (polar
    decomposition of the cross-covariance; visualization only)."""
    from mathutils import Matrix, Vector
    me = obj.data
    kb = me.shape_keys.key_blocks[name]
    P = [v.co.copy() for v in me.vertices]
    Q = [d.co.copy() for d in kb.data]
    cp = sum(P, Vector()) / len(P)
    cq = sum(Q, Vector()) / len(Q)
    H = Matrix(((0.0,) * 3,) * 3)
    for p, q in zip(P, Q):
        a, b = q - cq, p - cp
        for i in range(3):
            for j in range(3):
                H[i][j] += b[i] * a[j]
    R = H.copy()
    for _ in range(100):
        R = (R + R.inverted().transposed()) * 0.5
    if R.determinant() < 0:
        log("align: reflection, skipped")
        return
    for d, q in zip(kb.data, Q):
        d.co = cp + R @ (q - cq)
    log("aligned shape key")


ALIGN = True




# ----------------------------------------------------------------------------- runner
def run(steps):
    it = iter(steps)

    def tick():
        try:
            f = next(it)
        except StopIteration:
            bpy.ops.wm.quit_blender()
            return None
        try:
            d = f()
        except Exception:
            traceback.print_exc()
            log("FAILED in", f.__name__)
            bpy.ops.wm.quit_blender()
            return None
        if d is None:
            bpy.ops.wm.quit_blender()
            return None
        return d

    bpy.app.timers.register(tick, first_interval=3.0)


run({"unwrap": scene_unwrap, "cut": scene_cut, "sphere": scene_sphere,
     "developable": scene_developable}[SCENE]())
