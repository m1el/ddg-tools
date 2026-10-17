"""Blender render library of the showcase: studio scene, meshes with colors, tubes, spheres,
camera fitting, frame rendering. Imported by the per-scene scripts, which Blender runs with
`blender -b --python <scene>.py -- <args>` (see showcase/render.sh)."""

import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHOW = HERE.parent
_site = SHOW / ".venv" / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}" / "site-packages"
if str(_site) not in sys.path:
    sys.path.insert(0, str(_site))

import bmesh  # noqa: E402
import bpy  # noqa: E402
import numpy as np  # noqa: E402
from mathutils import Matrix, Vector  # noqa: E402

WORK = SHOW / "work"
MEDIA = SHOW / "media"
FRAMES = WORK / "frames"

# ----------------------------------------------------------------------------- palette
PALETTE = {
    "surface": (0.80, 0.80, 0.82),
    "ink": (0.06, 0.07, 0.09),
    "red": (0.80, 0.12, 0.10),
    "blue": (0.10, 0.30, 0.78),
    "gold": (0.95, 0.62, 0.10),
    "teal": (0.05, 0.55, 0.55),
    "green": (0.20, 0.62, 0.25),
    "purple": (0.50, 0.25, 0.70),
}

# Colormaps as control points (sRGB), sampled linearly.
_CMAPS = {
    "viridis": ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89",
                "#35b779", "#6ece58", "#b5de2b", "#fde725"],
    "magma": ["#000004", "#1c1044", "#4f127b", "#812581", "#b5367a", "#e55964",
              "#fb8761", "#fec287", "#fcfdbf"],
    "rdbu": ["#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#e6e4e2", "#fddbc7",
             "#f4a582", "#d6604d", "#b2182b"],
    "ylgnbu": ["#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0",
               "#225ea8", "#253494", "#081d58"],
    "spectral": ["#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#fee08b",
                 "#fdae61", "#f46d43", "#d53e4f", "#9e0142"],
    "warm": ["#fff5eb", "#fdd0a2", "#fd8d3c", "#d94801", "#7f2704"],
}


def _hex(h):
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]


def srgb_to_linear(c):
    c = np.asarray(c, dtype=np.float64)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def cmap(values, name="viridis", lo=None, hi=None):
    """Linear-RGB colors (n, 3) of `values` through a colormap."""
    v = np.asarray(values, dtype=np.float64)
    lo = np.nanmin(v) if lo is None else lo
    hi = np.nanmax(v) if hi is None else hi
    t = np.clip((v - lo) / max(hi - lo, 1e-300), 0, 1)
    pts = np.array([_hex(h) for h in _CMAPS[name]])
    x = np.linspace(0, 1, len(pts))
    rgb = np.stack([np.interp(t, x, pts[:, k]) for k in range(3)], axis=-1)
    return srgb_to_linear(rgb)


def lin(rgb):
    """sRGB tuple -> linear tuple (Blender colors are linear)."""
    return tuple(float(x) for x in srgb_to_linear(rgb))


# ----------------------------------------------------------------------------- scene
def reset():
    bpy.ops.wm.read_factory_settings(use_empty=True)


def setup(width=1600, height=1000, samples=96, transparent=True, exposure=0.0):
    scene = bpy.context.scene
    scene.render.engine = "CYCLES"
    prefs = bpy.context.preferences.addons["cycles"].preferences
    try:
        prefs.compute_device_type = "CUDA"
        prefs.get_devices()
        for d in prefs.devices:
            d.use = d.type == "CUDA"
        scene.cycles.device = "GPU"
    except Exception as e:  # CPU fallback
        print("GPU unavailable:", e)
    scene.cycles.samples = samples
    scene.cycles.use_adaptive_sampling = True
    scene.cycles.adaptive_threshold = 0.02
    scene.cycles.use_denoising = True
    try:
        scene.cycles.denoiser = "OPENIMAGEDENOISE"
    except Exception:
        pass
    scene.cycles.max_bounces = 6
    scene.render.resolution_x = width
    scene.render.resolution_y = height
    scene.render.resolution_percentage = 100
    scene.render.film_transparent = transparent
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGBA" if transparent else "RGB"
    scene.view_settings.view_transform = "Standard"
    scene.view_settings.look = "None"
    scene.view_settings.exposure = exposure
    world = bpy.data.worlds.new("world")
    scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes["Background"]
    bg.inputs["Color"].default_value = (0.9, 0.92, 0.95, 1)
    bg.inputs["Strength"].default_value = 0.3
    return scene


def area_light(name, location, target, size, energy, color=(1, 1, 1)):
    data = bpy.data.lights.new(name, "AREA")
    data.size = size
    data.energy = energy
    data.color = color
    obj = bpy.data.objects.new(name, data)
    bpy.context.collection.objects.link(obj)
    obj.location = location
    look_at(obj, target)
    return obj


def look_at(obj, target):
    d = Vector(target) - Vector(obj.location)
    obj.rotation_euler = d.to_track_quat("-Z", "Y").to_euler()


def studio(center, radius, ground_z=None, key=1.0):
    """Three soft area lights around a subject of `radius` at `center`, and a shadow catcher
    under it (z = ground_z)."""
    c = Vector(center)
    r = radius
    area_light("key", c + Vector((-1.6, -2.2, 2.4)) * r, c, 2.5 * r, 90 * r * r * key)
    area_light("fill", c + Vector((2.4, -1.4, 0.8)) * r, c, 3.0 * r, 30 * r * r * key,
               color=(0.92, 0.95, 1.0))
    area_light("rim", c + Vector((0.6, 2.6, 1.8)) * r, c, 2.0 * r, 70 * r * r * key)
    if ground_z is not None:
        bpy.ops.mesh.primitive_plane_add(size=40 * r, location=(c.x, c.y, ground_z))
        g = bpy.context.active_object
        g.name = "ground"
        g.is_shadow_catcher = True
        g.data.materials.append(mat_solid((1, 1, 1), roughness=1.0))
        return g


def camera(location, target, lens=50, ortho_scale=None, name="camera"):
    data = bpy.data.cameras.new(name)
    data.lens = lens
    data.sensor_width = 36
    data.clip_start = 0.001
    data.clip_end = 1e5
    if ortho_scale is not None:
        data.type = "ORTHO"
        data.ortho_scale = ortho_scale
    obj = bpy.data.objects.new(name, data)
    bpy.context.collection.objects.link(obj)
    obj.location = location
    look_at(obj, target)
    bpy.context.scene.camera = obj
    return obj


def fit_camera(points, azimuth=-35, elevation=18, lens=50, margin=1.08, target=None,
               aspect=None):
    """Perspective camera that frames `points` (n, 3) seen from azimuth/elevation (degrees,
    azimuth 0 looks along +Y, i.e. the camera sits at -Y)."""
    pts = np.asarray(points, dtype=np.float64)
    lo, hi = pts.min(0), pts.max(0)
    ctr = (lo + hi) / 2 if target is None else np.asarray(target, dtype=np.float64)
    a, e = math.radians(azimuth), math.radians(elevation)
    d = np.array([math.sin(a) * math.cos(e), -math.cos(a) * math.cos(e), math.sin(e)])
    # Camera basis.
    fwd = -d
    right = np.cross(fwd, [0, 0, 1.0])
    right /= np.linalg.norm(right)
    up = np.cross(right, fwd)
    scene = bpy.context.scene
    if aspect is None:
        aspect = scene.render.resolution_x / scene.render.resolution_y
    half_w = 18.0 / lens  # tan of half horizontal fov (sensor 36 mm)
    half_h = half_w / aspect if aspect >= 1 else half_w
    if aspect < 1:
        half_w = half_h * aspect
    # Smallest distance such that every point projects inside the frame.
    rel = pts - ctr
    x, y, z = rel @ right, rel @ up, rel @ fwd
    dist = np.max(np.maximum(np.abs(x) / half_w, np.abs(y) / half_h) - z) * margin
    loc = ctr + d * dist
    cam = camera(tuple(loc), tuple(ctr), lens=lens)
    return cam, ctr, dist


# ----------------------------------------------------------------------------- materials
def _principled(mat):
    mat.use_nodes = True
    nt = mat.node_tree
    bsdf = nt.nodes["Principled BSDF"]
    return nt, bsdf


def mat_solid(color, roughness=0.45, name="solid", metallic=0.0, coat=0.0):
    mat = bpy.data.materials.new(name)
    nt, bsdf = _principled(mat)
    bsdf.inputs["Base Color"].default_value = (*color, 1)
    bsdf.inputs["Roughness"].default_value = roughness
    bsdf.inputs["Metallic"].default_value = metallic
    if coat:
        bsdf.inputs["Coat Weight"].default_value = coat
    return mat


def mat_attr(attr="Col", roughness=0.5, name="vcol", lines=None, lines_period=1.0,
             lines_width=0.08, lines_color=(0.02, 0.02, 0.03), lines_strength=0.85):
    """Base color from a color attribute; optional isolines of a float attribute `lines`
    every `lines_period`, `lines_width` the fraction of a period."""
    mat = bpy.data.materials.new(name)
    nt, bsdf = _principled(mat)
    bsdf.inputs["Roughness"].default_value = roughness
    a = nt.nodes.new("ShaderNodeAttribute")
    a.attribute_name = attr
    color_out = a.outputs["Color"]
    if lines is not None:
        f = nt.nodes.new("ShaderNodeAttribute")
        f.attribute_name = lines
        div = nt.nodes.new("ShaderNodeMath")
        div.operation = "DIVIDE"
        nt.links.new(f.outputs["Fac"], div.inputs[0])
        div.inputs[1].default_value = lines_period
        fr = nt.nodes.new("ShaderNodeMath")
        fr.operation = "FRACT"
        nt.links.new(div.outputs[0], fr.inputs[0])
        # distance to the nearest integer: 0.5 - |fract - 0.5|
        sub = nt.nodes.new("ShaderNodeMath")
        sub.operation = "SUBTRACT"
        nt.links.new(fr.outputs[0], sub.inputs[0])
        sub.inputs[1].default_value = 0.5
        ab = nt.nodes.new("ShaderNodeMath")
        ab.operation = "ABSOLUTE"
        nt.links.new(sub.outputs[0], ab.inputs[0])
        dist = nt.nodes.new("ShaderNodeMath")
        dist.operation = "SUBTRACT"
        dist.inputs[0].default_value = 0.5
        nt.links.new(ab.outputs[0], dist.inputs[1])
        ss = nt.nodes.new("ShaderNodeMapRange")
        ss.interpolation_type = "SMOOTHSTEP"
        nt.links.new(dist.outputs[0], ss.inputs["Value"])
        ss.inputs["From Min"].default_value = lines_width * 0.5
        ss.inputs["From Max"].default_value = lines_width * 0.5 + 0.02
        ss.inputs["To Min"].default_value = 1.0 - lines_strength
        ss.inputs["To Max"].default_value = 1.0
        mix = nt.nodes.new("ShaderNodeMix")
        mix.data_type = "RGBA"
        mix.blend_type = "MIX"
        nt.links.new(ss.outputs["Result"], mix.inputs["Factor"])
        mix.inputs["A"].default_value = (*lines_color, 1)
        nt.links.new(a.outputs["Color"], mix.inputs["B"])
        color_out = mix.outputs["Result"]
    nt.links.new(color_out, bsdf.inputs["Base Color"])
    return mat


def mat_checker(scale=16.0, c1=(0.95, 0.95, 0.95), c2=(0.12, 0.13, 0.16), uv="UVMap",
                roughness=0.5, name="checker", tint_attr=None):
    """A UV checkerboard; with `tint_attr` the light squares take that color attribute."""
    mat = bpy.data.materials.new(name)
    nt, bsdf = _principled(mat)
    bsdf.inputs["Roughness"].default_value = roughness
    uvn = nt.nodes.new("ShaderNodeUVMap")
    uvn.uv_map = uv
    ch = nt.nodes.new("ShaderNodeTexChecker")
    ch.inputs["Scale"].default_value = scale
    ch.inputs["Color1"].default_value = (*c1, 1)
    ch.inputs["Color2"].default_value = (*c2, 1)
    nt.links.new(uvn.outputs["UV"], ch.inputs["Vector"])
    out = ch.outputs["Color"]
    if tint_attr:
        a = nt.nodes.new("ShaderNodeAttribute")
        a.attribute_name = tint_attr
        nt.links.new(a.outputs["Color"], ch.inputs["Color1"])
    nt.links.new(out, bsdf.inputs["Base Color"])
    return mat


def mat_image(path, uv="UVMap", roughness=0.55, name="image"):
    mat = bpy.data.materials.new(name)
    nt, bsdf = _principled(mat)
    bsdf.inputs["Roughness"].default_value = roughness
    uvn = nt.nodes.new("ShaderNodeUVMap")
    uvn.uv_map = uv
    tex = nt.nodes.new("ShaderNodeTexImage")
    tex.image = bpy.data.images.load(str(path))
    nt.links.new(uvn.outputs["UV"], tex.inputs["Vector"])
    nt.links.new(tex.outputs["Color"], bsdf.inputs["Base Color"])
    return mat


def mat_emission(color, strength=1.0, name="emit"):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.remove(nt.nodes["Principled BSDF"])
    em = nt.nodes.new("ShaderNodeEmission")
    em.inputs["Color"].default_value = (*color, 1)
    em.inputs["Strength"].default_value = strength
    nt.links.new(em.outputs[0], nt.nodes["Material Output"].inputs["Surface"])
    return mat


# ----------------------------------------------------------------------------- geometry
def mesh_object(name, positions, poly_verts, poly_sizes=None, smooth=True, material=None):
    """A mesh from positions (n, 3) and polygons: either (m, k) array or flat list + sizes."""
    pos = np.asarray(positions, dtype=np.float32)
    pv = np.asarray(poly_verts, dtype=np.int32)
    if poly_sizes is None:
        poly_sizes = np.full(pv.shape[0], pv.shape[1], dtype=np.int32)
        pv = pv.ravel()
    poly_sizes = np.asarray(poly_sizes, dtype=np.int32)
    starts = np.concatenate([[0], np.cumsum(poly_sizes)[:-1]]).astype(np.int32)
    me = bpy.data.meshes.new(name)
    me.vertices.add(len(pos))
    me.vertices.foreach_set("co", pos.ravel())
    me.loops.add(len(pv))
    me.loops.foreach_set("vertex_index", pv)
    me.polygons.add(len(poly_sizes))
    me.polygons.foreach_set("loop_start", starts)
    me.update(calc_edges=True)
    me.validate(verbose=False)
    if smooth:
        me.shade_smooth()
    obj = bpy.data.objects.new(name, me)
    bpy.context.collection.objects.link(obj)
    if material is not None:
        obj.data.materials.append(material)
    return obj


def set_point_colors(obj, rgb, name="Col"):
    me = obj.data
    attr = me.color_attributes.new(name, "FLOAT_COLOR", "POINT")
    rgba = np.ones((len(me.vertices), 4), dtype=np.float32)
    rgba[:, :3] = rgb
    attr.data.foreach_set("color", rgba.ravel())
    return attr


def set_corner_colors(obj, rgb, name="Col"):
    me = obj.data
    attr = me.color_attributes.new(name, "FLOAT_COLOR", "CORNER")
    rgba = np.ones((len(me.loops), 4), dtype=np.float32)
    rgba[:, :3] = rgb
    attr.data.foreach_set("color", rgba.ravel())
    return attr


def set_face_colors(obj, rgb, name="Col"):
    me = obj.data
    sizes = np.zeros(len(me.polygons), dtype=np.int32)
    me.polygons.foreach_get("loop_total", sizes)
    return set_corner_colors(obj, np.repeat(np.asarray(rgb), sizes, axis=0), name)


def set_point_scalar(obj, values, name):
    me = obj.data
    attr = me.attributes.new(name, "FLOAT", "POINT")
    attr.data.foreach_set("value", np.asarray(values, dtype=np.float32))
    return attr


def set_uvs(obj, uv, name="UVMap"):
    """Per-corner UVs (n_loops, 2)."""
    me = obj.data
    layer = me.uv_layers.new(name=name)
    layer.data.foreach_set("uv", np.asarray(uv, dtype=np.float32).ravel())
    return layer


def tubes(name, polylines, radius, material, closed=None, resolution=4):
    """One curve object of round tubes along polylines (list of (k, 3) arrays)."""
    cu = bpy.data.curves.new(name, "CURVE")
    cu.dimensions = "3D"
    cu.bevel_depth = radius
    cu.bevel_resolution = resolution
    cu.use_fill_caps = True
    for i, pl in enumerate(polylines):
        pl = np.asarray(pl, dtype=np.float64)
        if len(pl) < 2:
            continue
        sp = cu.splines.new("POLY")
        sp.points.add(len(pl) - 1)
        co = np.ones((len(pl), 4))
        co[:, :3] = pl
        sp.points.foreach_set("co", co.ravel())
        if closed is not None and closed[i]:
            sp.use_cyclic_u = True
    obj = bpy.data.objects.new(name, cu)
    bpy.context.collection.objects.link(obj)
    obj.data.materials.append(material)
    return obj


_ICO = {}


def _ico(subdiv):
    if subdiv not in _ICO:
        bm = bmesh.new()
        bmesh.ops.create_icosphere(bm, subdivisions=subdiv, radius=1.0)
        v = np.array([x.co[:] for x in bm.verts])
        f = np.array([[x.index for x in face.verts] for face in bm.faces])
        bm.free()
        _ICO[subdiv] = (v, f)
    return _ICO[subdiv]


def spheres(name, centers, radii, colors=None, material=None, subdiv=3):
    """Many spheres as one mesh; `colors` (n, 3) linear RGB go to the attribute 'Col'."""
    c = np.asarray(centers, dtype=np.float64).reshape(-1, 3)
    r = np.broadcast_to(np.asarray(radii, dtype=np.float64), (len(c),))
    v, f = _ico(subdiv)
    nv = len(v)
    pos = (v[None, :, :] * r[:, None, None] + c[:, None, :]).reshape(-1, 3)
    faces = (f[None, :, :] + nv * np.arange(len(c))[:, None, None]).reshape(-1, 3)
    obj = mesh_object(name, pos, faces, smooth=True,
                      material=material if material is not None else mat_attr("Col", 0.35))
    if colors is not None:
        set_point_colors(obj, np.repeat(np.asarray(colors).reshape(-1, 3), nv, axis=0))
    return obj


def arrows(name, origins, directions, length, radius, material, head=0.35):
    """Cone-tipped arrows as one mesh."""
    o = np.asarray(origins, dtype=np.float64).reshape(-1, 3)
    d = np.asarray(directions, dtype=np.float64).reshape(-1, 3)
    n = np.linalg.norm(d, axis=1)
    keep = n > 1e-12
    o, d = o[keep], d[keep] / n[keep, None]
    k = 8
    ang = np.linspace(0, 2 * np.pi, k, endpoint=False)
    # local frame per arrow
    tmp = np.where(np.abs(d[:, 2:3]) < 0.9, [[0, 0, 1.0]], [[1.0, 0, 0]])
    u = np.cross(d, tmp)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    w = np.cross(d, u)
    L = length
    shaft = L * (1 - head)
    rings = []
    for t, rad in ((0.0, radius), (shaft, radius), (shaft, radius * 2.6), (L, 0.0)):
        ring = (o[:, None, :] + d[:, None, :] * t + rad * (np.cos(ang)[None, :, None] * u[:, None, :]
                + np.sin(ang)[None, :, None] * w[:, None, :]))
        rings.append(ring)
    pos = np.concatenate(rings, axis=1)  # (m, 4k, 3)
    faces = []
    for ri in range(3):
        for j in range(k):
            a, b = ri * k + j, ri * k + (j + 1) % k
            faces.append([a, b, b + k, a + k])
    faces = np.array(faces)
    m = len(o)
    allf = (faces[None] + (4 * k) * np.arange(m)[:, None, None]).reshape(-1, 4)
    return mesh_object(name, pos.reshape(-1, 3), allf, smooth=True, material=material)


# ----------------------------------------------------------------------------- transforms
def y_up_to_z_up(p):
    p = np.asarray(p, dtype=np.float64)
    return np.stack([p[..., 0], -p[..., 2], p[..., 1]], axis=-1)


def normalize(p, size=2.0, center=None):
    """Scale so the bounding box diagonal is `size`, bottom at z = 0, centered in xy.
    Returns (points, (scale, offset)) so other data can follow with `apply_norm`."""
    p = np.asarray(p, dtype=np.float64)
    lo, hi = p.min(0), p.max(0)
    s = size / np.linalg.norm(hi - lo)
    off = -np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, lo[2]])
    return (p + off) * s, (s, off)


def apply_norm(p, norm):
    s, off = norm
    return (np.asarray(p, dtype=np.float64) + off) * s


def rotz(p, deg):
    a = math.radians(deg)
    R = np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]])
    return np.asarray(p) @ R.T


# ----------------------------------------------------------------------------- output
def render_still(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    bpy.context.scene.render.filepath = str(path)
    bpy.ops.render.render(write_still=True)
    print("rendered", path)


def render_frames(name, n, update, start=0):
    """Render frames `start..n-1` of an animation: `update(i)` changes the scene first.
    Frames go to work/frames/<name>/0000.png."""
    out = FRAMES / name
    out.mkdir(parents=True, exist_ok=True)
    only = os.environ.get("SHOW_FRAMES")  # "a:b" to render a subrange
    lo, hi = start, n
    if only:
        a, b = only.split(":")
        lo, hi = int(a), int(b)
    for i in range(lo, hi):
        update(i)
        f = out / f"{i:04d}.png"
        if f.exists() and os.environ.get("SHOW_KEEP"):
            continue
        bpy.context.scene.render.filepath = str(f)
        bpy.ops.render.render(write_still=True)
        print("frame", i, "of", n, flush=True)
    return out


def args():
    argv = sys.argv
    return argv[argv.index("--") + 1:] if "--" in argv else []


def load(name):
    return dict(np.load(WORK / f"{name}.npz", allow_pickle=False))


def ease(t):
    t = min(max(t, 0.0), 1.0)
    return t * t * (3 - 2 * t)
