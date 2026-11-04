# SPDX-License-Identifier: MIT
"""Points as Tori: signed distance to an oriented point cloud, as a Blender operator.

The work is done by ``ddg_tools.points_as_tori`` (Rust, crates/ddg-points-as-tori and
crates/ddg-py), after

    Nicole Feng, Ioannis Gkioulekas, Keenan Crane, Points as Tori: Fast Pointwise Signed Distance
    for Point Clouds, ACM Transactions on Graphics 45(4), Article 53, 2026.

``object.ddg_points_as_tori`` works on the active mesh object, a point cloud (its vertices, with
normals; see ``points_as_tori_io.py`` for where the normals come from and for the cache), in four
modes:

* **Surface**: marching cubes on the level set ``phi = level`` of the blended field, on a grid of
  ``resolution`` cells along the longest side of the box; a new mesh object.
* **Offset**: the same at ``phi = offset`` (the surface moved by that distance, outwards when
  positive); a new mesh object.
* **Boolean**: union, intersection or difference of the active cloud with a second one (another
  selected mesh object, or the one named in the panel); a new mesh object.
* **Orient Normals**: the orientation iteration of the paper's Fig. 27; the flipped normals are
  written back into the object (undoable), with the flips per iteration in the report.

The tori are cached per cloud (hash of positions and normals; memory and a sidecar next to the
.blend). The modal ``invoke`` runs the work on a worker thread while a timer shows progress; Esc or
the panel's Cancel button stops it: the tori build and the marching-cubes evaluations check the
native ``PatProgress`` cancel flag between chunks (orientation and the other native calls run to
their end in the background and their result is dropped). ``execute`` (redo, scripts, headless Blender) runs the same job in the
calling thread. Conventions as for the Repulsive Surfaces operators: scene settings in
``scene.ddg_points_as_tori`` copied into the operator on invoke, sidebar panel in the DDG tab,
device box greyed out without CUDA, a status line.
"""

import os
import threading
from math import ceil
from time import perf_counter

import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
    StringProperty,
)

from . import core
from . import points_as_tori_io as io
from .repulsive_surface_operator import DEVICE_ITEMS, gpu_info

MODE_ITEMS = [
    ('SURFACE', "Surface", "Mesh of the level set phi = Level of the signed distance field"),
    ('OFFSET', "Offset", "Mesh of the surface moved by the Offset distance (phi = Offset)"),
    ('BOOLEAN', "Boolean", "Union, intersection or difference with a second point cloud"),
    ('ORIENT', "Orient Normals", "Flip the normals of the cloud until they are consistent "
     "(the paper's Fig. 27); written back into the object"),
]
BOOLEAN_ITEMS = [
    ('UNION', "Union", "Points inside either cloud"),
    ('INTERSECT', "Intersection", "Points inside both clouds"),
    ('SUBTRACT', "Difference", "The active cloud minus the partner"),
]
SIGN_RULE_ITEMS = [
    ('REFERENCE', "Reference", "The sign rule of the authors' program (the torus of a saddle "
     "point is the mirror image of the osculating torus); for parity with their results"),
    ('GEOMETRIC', "Geometric", "The default: the torus center from the signed curvature, the "
     "osculating torus at saddle points too. Differs from Reference only at saddle points"),
]
EVALUATION_ITEMS = [
    ('PER_QUERY', "Per Query", "Blend the 32 nearest tori with a lambda that follows the "
     "neighborhood (the library's default)"),
    ('PAPER', "Paper", "The paper's Sec. 4.2 rule: a fixed lambda = 1000 / (mean point spacing) "
     "over a radius neighborhood"),
]

# Largest grid the marching cubes accepts: 2^29 values (crates/ddg-points-as-tori, src/mc.rs).
MAX_GRID_POINTS = 1 << 29
OP_LABELS = {'UNION': "Union", 'INTERSECT': "Intersection", 'SUBTRACT': "Difference"}
OP_ARGUMENT = {'UNION': "union", 'INTERSECT': "intersect", 'SUBTRACT': "subtract"}


class PatError(Exception):
    """The job failed; the message is meant for the user."""


class Cancelled(Exception):
    pass


def option_properties():
    """The options shared by the operator and the scene settings (fresh property objects)."""
    return {
        "mode": EnumProperty(name="Mode", items=MODE_ITEMS, default='SURFACE',
                             description="What to compute from the point cloud"),
        "level": FloatProperty(
            name="Level", default=0.0, soft_min=-1.0, soft_max=1.0, unit='LENGTH',
            description="Level of the signed distance field to mesh (negative: inside)"),
        "offset": FloatProperty(
            name="Offset", default=0.05, soft_min=-1.0, soft_max=1.0, unit='LENGTH',
            description="Distance to move the surface (positive: outwards)"),
        "resolution": IntProperty(
            name="Resolution", default=128, min=4, max=1000, soft_max=512,
            description="Grid cells along the longest side of the bounding box. The dense grid "
                        "must stay under 2^29 values; larger requests are refused with a message"),
        "margin": FloatProperty(
            name="Margin", default=0.1, min=0.0, soft_max=0.5,
            description="Space around the cloud's bounding box, as a fraction of its diagonal "
                        "(an outward level or offset is added on top)"),
        "partner": StringProperty(
            name="Partner", default="",
            description="Second point cloud for Boolean (empty: the other selected mesh object)"),
        "boolean_op": EnumProperty(name="Operation", items=BOOLEAN_ITEMS, default='UNION'),
        "boolean_smooth": FloatProperty(
            name="Smooth", default=0.0, min=0.0, soft_max=0.5, unit='LENGTH',
            description="Blend radius of a smooth Boolean (0: the hard min / max)"),
        "orient_iterations": IntProperty(
            name="Iterations", default=10, min=1, soft_max=50,
            description="Maximum orientation iterations (the paper uses 10)"),
        "device": EnumProperty(
            name="Device", items=DEVICE_ITEMS, default='AUTO',
            description="Where the network and the field evaluation run. Auto uses an NVIDIA GPU "
                        "when the extension has CUDA support and one is usable. CUDA needs the "
                        "extension built with CUDA support"),
        "cuda_ordinal": IntProperty(name="GPU Index", default=0, min=0, soft_max=8,
                                    description="Which CUDA device to use (0: the first)"),
        "sign_rule": EnumProperty(name="Sign Rule", items=SIGN_RULE_ITEMS, default='GEOMETRIC',
                                  description="How the torus of a point is placed (changes the "
                                              "tori, so it is part of the cache key)"),
        "evaluation": EnumProperty(name="Evaluation", items=EVALUATION_ITEMS,
                                   default='PER_QUERY',
                                   description="How the tori are blended at a query point"),
        "threads": IntProperty(name="Threads", default=0, min=0, soft_max=64,
                               description="Native threads (0: the library's choice)"),
        "normal_attribute": StringProperty(
            name="Normals", default="",
            description="Point attribute (Vector) with the normals; empty: 'normal', nx/ny/nz, "
                        "custom normals, or the vertex normals of a mesh with faces"),
        "use_cache": BoolProperty(
            name="Use Cache", default=True,
            description="Reuse the tori of an unchanged cloud (memory and a sidecar directory "
                        "next to the .blend file)"),
    }


OPTION_NAMES = tuple(option_properties())


class DDG_PG_points_as_tori_settings(bpy.types.PropertyGroup):
    """The options of the sidebar panel; the operator starts from them."""

    __annotations__ = option_properties()


def copy_options(source, target):
    for name in OPTION_NAMES:
        setattr(target, name, getattr(source, name))


def snapshot_options(settings):
    """The options as a plain dict (the worker thread must not touch Blender properties)."""
    return {name: getattr(settings, name) for name in OPTION_NAMES}


# ----------------------------------------------------------------------------------------------
# The cache and the fields (no bpy below this line until the operators)


def effective_device(requested):
    """'cpu' or 'cuda' for the setting 'cpu' / 'cuda' / 'auto' (what the library will use)."""
    if requested == 'cuda':
        return 'cuda'
    if requested == 'auto':
        return 'cuda' if gpu_info()["available"] else 'cpu'
    return 'cpu'


def build_kwargs(opts):
    kwargs = {"device": opts["device"].lower(), "cuda_ordinal": opts["cuda_ordinal"],
              "sign_rule": opts["sign_rule"].lower()}
    if opts["threads"] > 0:
        kwargs["threads"] = opts["threads"]
    return kwargs


class CancelEvent(threading.Event):
    """A ``threading.Event`` that also cancels the native call in flight (``ddg_tools.PatProgress``),
    so Esc stops a build at its next chunk instead of letting it finish in the background."""

    def __init__(self):
        super().__init__()
        self.native = None

    def set(self):
        super().set()
        if self.native is not None:
            self.native.cancel()


def get_field(cloud, opts, directory, native_progress=None):
    """``(field, info)``: the blended field of ``cloud`` (from the cache or built) and a dict:
    ``cache`` ('memory', 'sidecar' or 'miss'), ``seconds``, ``device`` (where the tori were
    predicted: 'cpu' or 'cuda'), ``warnings``."""
    ddg = core.native()
    t0 = perf_counter()
    kwargs = build_kwargs(opts)
    key = (cloud.digest, kwargs["sign_rule"])
    use_cache = opts["use_cache"]
    warnings = []
    field, cache, device = None, "miss", None
    if use_cache:
        field = io.memory_get(key)
        if field is not None:
            cache = "memory"
    path = io.sidecar_path(directory, *key) if directory else None
    if field is None and use_cache and path and os.path.isfile(path):
        try:
            field = ddg.pat_load_tori(path, cloud.positions, cloud.normals, **kwargs)
            if field.sign_rule not in (None, kwargs["sign_rule"]):
                raise ValueError("the file was fitted with the {} sign rule".format(field.sign_rule))
            cache = "sidecar"
        except ValueError as ex:
            warnings.append("ignored the sidecar file {} ({})".format(os.path.basename(path), ex))
            field = None
    if field is None:
        try:
            field = ddg.PointsAsTori(cloud.positions, cloud.normals, progress=native_progress, **kwargs)
            if native_progress is not None:
                field = field.with_options(progress=None)  # the cache must not keep this job's handle
        except ValueError as ex:
            if native_progress is not None and native_progress.cancelled:
                raise Cancelled() from ex
            raise PatError("{}: {}".format(cloud.name, ex)) from ex
        report = field.build_report or {}
        device = report.get("predict_device")
        if use_cache and path:
            try:
                os.makedirs(directory, exist_ok=True)
                field.save_tori(path)
            except (OSError, ValueError) as ex:
                warnings.append("could not write the sidecar file {}: {}".format(path, ex))
    if use_cache:
        io.memory_put(key, field)
    # A cached field may have been made with another evaluation rule, device or thread count.
    eval_options = {"device": kwargs["device"], "cuda_ordinal": kwargs["cuda_ordinal"],
                    "paper": opts["evaluation"] == 'PAPER'}
    if native_progress is not None:
        eval_options["progress"] = native_progress
    if opts["threads"] > 0:
        eval_options["threads"] = opts["threads"]
    try:
        field = field.with_options(**eval_options)
    except ValueError as ex:
        raise PatError("{}: {}".format(cloud.name, ex)) from ex
    return field, {"cache": cache, "seconds": perf_counter() - t0, "device": device,
                   "warnings": warnings}


def box_for(clouds, opts, level):
    """Bounding box ``[x0, y0, z0, x1, y1, z1]`` of the clouds with the margin and an outward
    level added."""
    boxes = [c.bounds() for c in clouds]
    lo = [min(b[k] for b in boxes) for k in range(3)]
    hi = [max(b[3 + k] for b in boxes) for k in range(3)]
    diag = sum((h - l) ** 2 for l, h in zip(lo, hi)) ** 0.5
    pad = opts["margin"] * diag + max(level, 0.0) + (
        opts["boolean_smooth"] if opts["mode"] == 'BOOLEAN' else 0.0)
    if diag == 0.0 and pad == 0.0:
        pad = 1.0
    return [l - pad for l in lo] + [h + pad for h in hi]


def check_grid(box, resolution):
    """Raises ``PatError`` if the grid would exceed the library's dense-grid limit."""
    extent = [box[3 + k] - box[k] for k in range(3)]
    h = max(extent) / resolution
    cells = [max(1, ceil(e / h - 1e-9)) for e in extent]
    points = (cells[0] + 1) * (cells[1] + 1) * (cells[2] + 1)
    if points > MAX_GRID_POINTS:
        raise PatError("resolution {} makes a grid of {:,} points (limit {:,}, about {:.1f} GB "
                       "of values); lower the resolution".format(
                           resolution, points, MAX_GRID_POINTS, 8.0 * points / 1e9))
    return cells, points


class Job:
    """One run of the operator, everything but the Blender data.

    ``mode``, ``clouds`` (``io.Cloud`` list: the active cloud, then the partner for Boolean),
    ``opts`` (plain dict of the options), ``directory`` (sidecar directory or None). ``run()`` fills
    ``result`` (a dict) or raises; ``progress`` (0..1) and ``stage`` (text) can be read from
    another thread; ``cancel`` is an ``Event`` checked between stages."""

    def __init__(self, clouds, opts, directory):
        self.clouds = clouds
        self.opts = opts
        self.mode = opts["mode"]
        self.directory = directory
        self.cancel = CancelEvent()
        self.native_progress = None
        try:
            make = getattr(core.native(), "PatProgress", None)  # wheels before WP15 have none
        except Exception:  # a broken installation is reported by the job itself
            make = None
        if make is not None:
            self.native_progress = make()
            self.cancel.native = self.native_progress
        self.progress = 0.0
        self.stage = "starting"
        self.stages = []  # one dict per finished stage
        self.result = None
        self.error = None
        self.done = False
        self.started = perf_counter()

    def _check(self):
        if self.cancel.is_set():
            raise Cancelled()

    def _field(self, index, label, fraction):
        cloud = self.clouds[index]
        self.stage = "{}: tori of '{}' ({:,} points)".format(label, cloud.name, cloud.n)
        self._check()
        field, info = get_field(cloud, self.opts, self.directory, self.native_progress)
        info.update(stage="build", cloud=cloud.name, n=cloud.n)
        self.stages.append(info)
        self.progress = fraction
        self._check()
        return field

    def run(self):
        try:
            self._run()
        except Cancelled:
            self.error = Cancelled()
        except Exception as ex:  # reported by the main thread
            # a cancelled native call raises its own error: it is the user's Esc, not a failure
            self.error = Cancelled() if self.cancel.is_set() else ex
        finally:
            self.done = True

    def _run(self):
        if self.mode == 'ORIENT':
            self._orient()
        else:
            self._surface()

    def _orient(self):
        ddg = core.native()
        cloud = self.clouds[0]
        o = self.opts
        self.stage = "orienting the normals of '{}' ({:,} points)".format(cloud.name, cloud.n)
        t0 = perf_counter()
        kwargs = build_kwargs(o)
        kwargs["max_iterations"] = o["orient_iterations"]
        try:
            fixed, report = ddg.pat_orient_normals(cloud.positions, cloud.normals, **kwargs)
        except ValueError as ex:
            raise PatError("{}: {}".format(cloud.name, ex)) from ex
        self._check()
        flipped = [i for i in range(cloud.n)
                   if cloud.normals[3 * i] * fixed[3 * i] + cloud.normals[3 * i + 1] * fixed[3 * i + 1]
                   + cloud.normals[3 * i + 2] * fixed[3 * i + 2] < 0.0]
        self.stages.append({"stage": "orient", "cloud": cloud.name, "n": cloud.n,
                            "seconds": perf_counter() - t0, "device": effective_device(kwargs["device"])})
        self.progress = 1.0
        self.result = {"flipped": flipped, "report": report}

    def _surface(self):
        o = self.opts
        mode = self.mode
        total = 2 if mode != 'BOOLEAN' else 3
        field = self._field(0, "1/{}".format(total), 1.0 / total)
        level = o["level"]
        if mode == 'OFFSET':
            level = o["offset"]
        label = {'SURFACE': "Surface", 'OFFSET': "Offset"}.get(mode)
        if mode == 'BOOLEAN':
            other = self._field(1, "2/3", 1.0 / 3.0)
            smooth = o["boolean_smooth"]  # 0 is the hard min / max in the library
            try:
                field = getattr(field, OP_ARGUMENT[o["boolean_op"]])(other, smooth=smooth)
            except ValueError as ex:
                raise PatError("Boolean failed: {}".format(ex)) from ex
            label = OP_LABELS[o["boolean_op"]]
            level = o["level"]
        box = box_for(self.clouds, o, level)
        cells, points = check_grid(box, o["resolution"])
        self.stage = "marching cubes, {} x {} x {} cells at level {:.4g}".format(*cells, level)
        t0 = perf_counter()
        try:
            vertices, faces, report = field.reconstruct(
                level=level, resolution=o["resolution"], bbox=box, return_report=True)
        except ValueError as ex:
            raise PatError("marching cubes failed: {}".format(ex)) from ex
        self._check()
        self.stages.append({"stage": "mesh", "seconds": perf_counter() - t0,
                            "device": effective_device(o["device"].lower()),
                            "triangles": len(faces) // 3, "grid_points": points,
                            "components": (report.get("mesh") or {}).get("components")})
        self.progress = 1.0
        self.result = {"vertices": vertices, "faces": faces, "report": report, "label": label,
                       "level": level, "box": box}


def stage_text(stage):
    """One short phrase for a stage record."""
    if stage["stage"] == "build":
        hit = {"memory": "cache hit (memory)", "sidecar": "cache hit (sidecar file)",
               "miss": "cache miss, built"}[stage["cache"]]
        device = " on {}".format(stage["device"].upper()) if stage.get("device") else ""
        return "tori of '{}': {}{}, {:.2f} s".format(stage["cloud"], hit, device, stage["seconds"])
    if stage["stage"] == "orient":
        return "orientation on {}, {:.2f} s".format(stage["device"].upper(), stage["seconds"])
    return "mesh: {:,} triangles, evaluated on {}, {:.2f} s".format(
        stage["triangles"], stage["device"].upper(), stage["seconds"])


def summary(job, extra=""):
    """The status line: the stages, and the total time."""
    text = "; ".join(stage_text(s) for s in job.stages)
    if extra:
        text = extra + ": " + text
    return "{} (total {:.2f} s)".format(text, perf_counter() - job.started)


# What the panel shows after a run: {"text": ..., "stages": [...], "error": bool}.
LAST = {"text": "", "error": False}
# The running modal job, for the panel's Cancel button.
_RUN = {"job": None}


# ----------------------------------------------------------------------------------------------
# Blender side: preparing the job, applying the result


def prepare(context, opts, active=None):
    """Reads the clouds named by ``opts`` from the scene. Returns ``(clouds, directory)``.
    Raises ``io.CloudInputError``."""
    active = active or context.active_object
    clouds = [io.read_cloud(active, opts["normal_attribute"])]
    if opts["mode"] == 'BOOLEAN':
        name = opts["partner"]
        if name:
            partner = bpy.data.objects.get(name)
            if partner is None:
                raise io.CloudInputError("the Boolean partner '{}' does not exist".format(name))
        else:
            others = io.other_cloud_objects(context, active)
            if not others:
                raise io.CloudInputError(
                    "Boolean needs a second point cloud: select another mesh object or name it "
                    "as the Partner")
            partner = others[0]
        if partner is active:
            raise io.CloudInputError("the Boolean partner is the active object itself")
        clouds.append(io.read_cloud(partner, opts["normal_attribute"]))
    return clouds, io.sidecar_dir()


def apply_result(context, job):
    """Writes ``job.result`` into the scene. Returns the created object or None (orient)."""
    clouds, result = job.clouds, job.result
    cloud = clouds[0]
    if job.mode == 'ORIENT':
        if result["flipped"]:
            place = io.write_flips(cloud, result["flipped"])
            job.stages[-1]["written"] = place
        return None
    if not result["faces"]:
        raise PatError(
            "no surface at level {:.4g}: the field has no sign change around it inside the "
            "box (an empty Boolean, or a level outside the field's range); try another level, "
            "or a larger margin".format(result["level"]))
    label = result["label"]
    if job.mode == 'BOOLEAN':
        name = "{} {} {}".format(cloud.name, label, clouds[1].name)
    else:
        name = "{} {}".format(cloud.name, label)
    obj = io.new_surface_object(cloud, name, result["vertices"], result["faces"],
                                cloud.obj.matrix_world.copy(), list(cloud.obj.users_collection))
    # The selection stays on the cloud, so that a second run (Adjust Last Operation, another
    # level) starts from the cloud again.
    return obj


def orient_summary(job):
    report = job.result["report"]
    flips = report.get("flips_per_iteration") or []
    return "flips per iteration {}, {:,} normals flipped in total ({}converged)".format(
        flips, len(job.result["flipped"]), "" if report.get("converged") else "not ")


def finish_job(operator, context, job):
    """Applies a finished job and reports. Returns the operator result set."""
    try:
        if job.error is not None:
            raise job.error
        obj = apply_result(context, job)
    except Cancelled:
        LAST.update(text="cancelled", error=False)
        operator.report({'INFO'}, "Points as Tori cancelled")
        return {'CANCELLED'}
    except (PatError, io.CloudInputError) as ex:
        LAST.update(text=str(ex), error=True)
        operator.report({'ERROR'}, str(ex))
        return {'CANCELLED'}
    except Exception as ex:  # the library's own errors that were not caught as ValueError
        LAST.update(text="Points as Tori failed: {}".format(ex), error=True)
        operator.report({'ERROR'}, "Points as Tori failed: {}".format(ex))
        return {'CANCELLED'}
    for stage in job.stages:
        for warning in stage.get("warnings", ()):
            operator.report({'WARNING'}, warning)
    if job.mode == 'ORIENT':
        written = job.stages[-1].get("written")
        text = summary(job, "Oriented '{}'".format(job.clouds[0].name)) + "; " + orient_summary(job)
        if written:
            text += "; written to " + written
    else:
        text = summary(job, "'{}'".format(obj.name))
    LAST.update(text=text, error=False, stages=list(job.stages))
    operator.report({'INFO'}, text)
    return {'FINISHED'}


def run_sync(context, opts):
    """The job for ``opts`` run in this thread, without applying it; used by ``execute`` and the
    tests. Raises ``io.CloudInputError``."""
    clouds, directory = prepare(context, opts)
    job = Job(clouds, opts, directory)
    job.run()
    return job


# ----------------------------------------------------------------------------------------------
# Operators


def _poll_object(context):
    obj = context.active_object
    return obj is not None and obj.type == 'MESH' and obj.mode == 'OBJECT' \
        and _RUN["job"] is None


def draw_options(layout, settings, context, *, in_operator=False):
    """The options, for the redo panel and the sidebar."""
    layout.use_property_split = True
    layout.use_property_decorate = False
    layout.prop(settings, "mode")
    mode = settings.mode
    if mode == 'SURFACE':
        layout.prop(settings, "level")
    elif mode == 'OFFSET':
        layout.prop(settings, "offset")
    elif mode == 'BOOLEAN':
        layout.prop_search(settings, "partner", bpy.data, "objects", text="Partner")
        layout.prop(settings, "boolean_op")
        layout.prop(settings, "boolean_smooth")
        layout.prop(settings, "level")
    else:
        layout.prop(settings, "orient_iterations")
    if mode != 'ORIENT':
        layout.prop(settings, "resolution")
        layout.prop(settings, "margin")
    box = layout.box()
    box.label(text="Device")
    info = gpu_info()
    col = box.column()
    col.enabled = bool(info["available"])
    col.prop(settings, "device")
    sub = col.column()
    sub.active = settings.device != 'CPU'
    sub.prop(settings, "cuda_ordinal")
    if info["available"]:
        names = ", ".join("{}: {}".format(d["ordinal"], d["name"]) for d in info["devices"])
        box.label(text=names, icon='INFO')
    else:
        box.label(text="CUDA is not available; running on the CPU", icon='INFO')
        box.label(text=info["reason"] or "no usable GPU")
    header, body = layout.panel("ddg_points_as_tori_advanced_" + ("op" if in_operator else "ui"),
                                default_closed=True)
    header.label(text="Advanced")
    if body is not None:
        draw_advanced(body, settings)


def draw_advanced(layout, settings):
    layout.use_property_split = True
    layout.use_property_decorate = False
    layout.prop(settings, "sign_rule")
    layout.prop(settings, "evaluation")
    layout.prop(settings, "threads")
    layout.prop(settings, "normal_attribute")
    layout.prop(settings, "use_cache")


class OBJECT_OT_ddg_points_as_tori(bpy.types.Operator):
    bl_idname = "object.ddg_points_as_tori"
    bl_label = "Points as Tori"
    bl_description = (
        "Signed distance to the active point cloud (a mesh object: its vertices and normals) by "
        "blending fitted tori; method: Feng, Gkioulekas, Crane, Points as Tori (ACM TOG 2026). "
        "Modes: mesh of a level set (Surface), of an offset, a Boolean with a second cloud, or "
        "orient the normals. The tori are cached per cloud. Esc cancels"
    )
    bl_options = {'REGISTER', 'UNDO'}

    __annotations__ = option_properties()

    _timer = None
    _job = None
    _thread = None
    _area = None

    @classmethod
    def poll(cls, context):
        return _poll_object(context)

    def draw(self, context):
        draw_options(self.layout, self, context, in_operator=True)

    def _opts(self):
        return snapshot_options(self)

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        try:
            job = run_sync(context, self._opts())
        except io.CloudInputError as ex:
            LAST.update(text=str(ex), error=True)
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        return finish_job(self, context, job)

    def invoke(self, context, event):
        copy_options(context.scene.ddg_points_as_tori, self)
        try:
            clouds, directory = prepare(context, self._opts())
        except io.CloudInputError as ex:
            LAST.update(text=str(ex), error=True)
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        job = self._job = Job(clouds, self._opts(), directory)
        _RUN["job"] = job
        self._thread = threading.Thread(target=job.run, name="ddg-points-as-tori", daemon=True)
        self._thread.start()
        wm = context.window_manager
        wm.progress_begin(0, 100)
        self._timer = wm.event_timer_add(0.1, window=context.window)
        wm.modal_handler_add(self)
        self._area = context.area
        self._status(context)
        return {'RUNNING_MODAL'}

    def _status(self, context):
        job = self._job
        text = "Points as Tori: {} ({:.0f} s)   [Esc: cancel]".format(
            job.stage, perf_counter() - job.started)
        if self._area is not None:
            self._area.header_text_set(text)
        context.workspace.status_text_set(text)
        context.window_manager.progress_update(int(100 * job.progress))

    def _cleanup(self, context):
        wm = context.window_manager
        if self._timer is not None:
            wm.event_timer_remove(self._timer)
            self._timer = None
        wm.progress_end()
        if self._area is not None:
            self._area.header_text_set(None)
            self._area.tag_redraw()
        context.workspace.status_text_set(None)
        _RUN["job"] = None

    def modal(self, context, event):
        job = self._job
        if event.type == 'ESC' and event.value == 'PRESS' or job.cancel.is_set():
            job.cancel.set()
            self._cleanup(context)
            LAST.update(text="cancelled", error=False)
            self.report({'INFO'}, "Points as Tori cancelled (a native call without a cancel hook "
                                  "finishes in the background; its result is dropped)")
            return {'CANCELLED'}
        if event.type == 'TIMER':
            if job.done:
                self._cleanup(context)
                return finish_job(self, context, job)
            self._status(context)
            for area in context.window.screen.areas:
                if area.type == 'VIEW_3D':
                    area.tag_redraw()
            return {'RUNNING_MODAL'}
        return {'PASS_THROUGH'}

    def cancel(self, context):
        if self._job is not None and _RUN["job"] is self._job:
            self._job.cancel.set()
            self._cleanup(context)


class OBJECT_OT_ddg_points_as_tori_cancel(bpy.types.Operator):
    bl_idname = "object.ddg_points_as_tori_cancel"
    bl_label = "Cancel Points as Tori"
    bl_description = "Cancel the running Points as Tori job"
    bl_options = {'INTERNAL'}

    @classmethod
    def poll(cls, context):
        return _RUN["job"] is not None

    def execute(self, context):
        _RUN["job"].cancel.set()
        return {'FINISHED'}


class OBJECT_OT_ddg_points_as_tori_clear_cache(bpy.types.Operator):
    bl_idname = "object.ddg_points_as_tori_clear_cache"
    bl_label = "Clear Tori Cache"
    bl_description = ("Forget the cached tori in memory and delete the sidecar files next to the "
                      ".blend file")
    bl_options = {'REGISTER'}

    def execute(self, context):
        dropped, deleted = io.clear_cache(io.sidecar_dir())
        self.report({'INFO'}, "Cleared the tori cache: {} in memory, {} sidecar files".format(
            dropped, deleted))
        return {'FINISHED'}


class VIEW3D_PT_ddg_points_as_tori(bpy.types.Panel):
    bl_label = "Points as Tori"
    bl_idname = "VIEW3D_PT_ddg_points_as_tori"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "DDG"

    def draw(self, context):
        layout = self.layout
        job = _RUN["job"]
        if job is not None:
            layout.label(text=job.stage)
            layout.progress(factor=job.progress, text="{:.0f} %".format(100 * job.progress))
            layout.operator(OBJECT_OT_ddg_points_as_tori_cancel.bl_idname, text="Cancel",
                            icon='X')
            return
        settings = context.scene.ddg_points_as_tori
        label = "Orient Normals" if settings.mode == 'ORIENT' else "Compute Mesh"
        layout.operator(OBJECT_OT_ddg_points_as_tori.bl_idname, text=label, icon='PLAY')
        draw_options(layout, settings, context)
        if LAST["text"]:
            box = layout.box()
            col = box.column(align=True)
            for line in _wrap(LAST["text"], 44):
                col.label(text=line, icon='ERROR' if LAST["error"] else 'NONE')
        layout.operator(OBJECT_OT_ddg_points_as_tori_clear_cache.bl_idname, icon='TRASH')


def _wrap(text, width):
    lines, line = [], ""
    for word in text.split():
        if line and len(line) + 1 + len(word) > width:
            lines.append(line)
            line = word
        else:
            line = (line + " " + word).strip()
    if line:
        lines.append(line)
    return lines


CLASSES = (
    DDG_PG_points_as_tori_settings,
    OBJECT_OT_ddg_points_as_tori,
    OBJECT_OT_ddg_points_as_tori_cancel,
    OBJECT_OT_ddg_points_as_tori_clear_cache,
    VIEW3D_PT_ddg_points_as_tori,
)


def _draw_menu_entry(self, _context):
    self.layout.separator()
    self.layout.operator(OBJECT_OT_ddg_points_as_tori.bl_idname)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.Scene.ddg_points_as_tori = PointerProperty(type=DDG_PG_points_as_tori_settings)
    bpy.types.VIEW3D_MT_object.append(_draw_menu_entry)


def unregister():
    bpy.types.VIEW3D_MT_object.remove(_draw_menu_entry)
    del bpy.types.Scene.ddg_points_as_tori
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
