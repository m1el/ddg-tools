# SPDX-License-Identifier: MIT
"""Repulsive Surfaces operators: a live (modal) flow and an undoable N-steps flow.

The work is done by ``ddg_tools.RepulsiveSurfaceSession`` (Rust, crates/ddg-repulsive-surfaces and
crates/ddg-py), after

    Chris Yu, Caleb Brakensiek, Henrik Schumacher, Keenan Crane, Repulsive Surfaces,
    ACM Transactions on Graphics 40(6), 2021,

and the reference implementation by Nick Sharp (MIT) that the Rust port follows.

The flow decreases the tangent-point energy of the active triangle mesh (n-gons are split by
Blender's loop triangles) by fractional Sobolev descent, under barycenter, pin, area and volume
constraints, with penalties, mesh obstacles, implicit obstacles and attractors (empties),
continuous collision detection and remeshing. See ``repulsive_surface_io.py`` for what is read and
how the result is written back (positions only when the triangles are unchanged, else the mesh
data is replaced and attributes are dropped with a warning).

* ``object.ddg_repulsive_surface_flow``: N steps (or to convergence) in one go, one undo step,
  with a redo panel.
* ``object.ddg_repulsive_surface_live``: a worker thread steps the session while a timer (1/30 s)
  copies the last completed step into the object. Enter, right click, or the Accept button of the
  sidebar panel (``object.ddg_repulsive_surface_live_finish``) keeps the result as one undo step;
  Esc or the Cancel button restores the original mesh exactly. The options are those of the
  sidebar panel (``scene.ddg_repulsive_surface``).

``LiveFlow`` holds everything the modal operator does apart from events, header text and the
timer, so that it is testable in a background Blender (no timers there).
"""

import threading
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
from . import repulsive_surface_io as io

SPACE_ITEMS = [
    ('WORLD', "World", "Run in world space: the object's transform (scale!) counts as it would "
     "physically, obstacles are used as they stand"),
    ('LOCAL', "Local", "Run in the object's local space; obstacles are transformed into it"),
]
METRIC_ITEMS = [
    ('SOBOLEV', "Fractional Sobolev (Hs)", "The paper's metric: the tangent-point energy's own "
     "fractional Sobolev inner product (the default)"),
    ('L2', "L2", "Plain L2 gradient (for experiments: slow and stiff)"),
    ('H1', "H1", "H1 gradient (for experiments)"),
    ('H2', "H2", "H2 gradient (for experiments)"),
]
BARYCENTER_ITEMS = [
    ('AUTO', "Auto", "Fix the barycenter unless pins are used (the paper: pins make it "
     "unnecessary)"),
    ('ON', "On", "Always fix the barycenter"),
    ('OFF', "Off", "No barycenter constraint (needs pins, or the flow refuses: the surface "
     "would translate)"),
]
MODE_ITEMS = [
    ('REPEL', "Repel", "Empties repel the surface (an implicit obstacle)"),
    ('ATTRACT', "Attract", "Empties attract the surface (an implicit attractor)"),
]

DEVICE_ITEMS = [
    ('CPU', "CPU", "Run everything on the CPU (the reference; the default)"),
    ('CUDA', "CUDA", "Run the hierarchical products and the Barnes-Hut energy on an NVIDIA GPU; "
     "the flow fails with a message if no usable GPU is found"),
    ('AUTO', "Auto", "Use an NVIDIA GPU when one is usable, the CPU otherwise"),
]
DETERMINISM_ITEMS = [
    ('RUN_TO_RUN', "Run to Run", "Bitwise the same result on every run (the default); within "
     "round-off of the CPU result"),
    ('MATCH_CPU', "Match CPU", "Bit for bit the CPU result (integer exponents; slower)"),
]

_GPU_INFO = []


def gpu_info(refresh=False):
    """``ddg_tools.gpu_info()`` (cached: it opens the driver), or a dict saying the module is
    too old or cannot be loaded. Keys: ``cuda_compiled``, ``available``, ``devices``,
    ``reason``."""
    if refresh or not _GPU_INFO:
        try:
            info = core.native().gpu_info()
        except (core.UnwrapFailed, AttributeError) as ex:
            info = {"cuda_compiled": False, "available": False, "devices": [],
                    "reason": "the native module has no GPU support information ({})".format(ex)}
        _GPU_INFO[:] = [info]
    return _GPU_INFO[0]


def device_text(device, name, metric_device=None):
    """'CPU' or 'CUDA (name)' for a step report's ``device`` and ``device_name``; when the GPU is
    open but the metric was solved by the dense back end on the CPU (small meshes), says so."""
    if device == "cuda":
        text = "CUDA ({})".format(name) if name else "CUDA"
        if metric_device == "cpu":
            text += ", dense metric on CPU"
        return text
    return "CPU"


STATUS_LABELS = {
    "running": "running",
    "converged": "converged",
    "max_iterations": "step limit reached",
    "time_limit": "time limit reached",
    "stuck": "stuck (the line search found no step)",
    "cancelled": "stopped",
    "error": "error",
}

PENALTY_WEIGHTS = (
    # (property, penalty type)
    ("w_area_deviation", "area_deviation"),
    ("w_volume_deviation", "volume_deviation"),
    ("w_boundary_length", "boundary_length"),
    ("w_boundary_curvature", "boundary_curvature"),
    ("w_willmore", "willmore"),
    ("w_squared_error", "squared_error"),
)


def option_properties():
    """The options shared by the operators and the scene settings (fresh property objects)."""
    props = {
        "alpha": FloatProperty(name="Alpha", default=6.0, min=1.0, soft_max=16.0,
                               description="Exponent alpha of the tangent-point kernel"),
        "beta": FloatProperty(name="Beta", default=12.0, min=2.0, soft_max=24.0,
                              description="Exponent beta of the tangent-point kernel (the "
                                          "energy is finite for beta - alpha above 2 on "
                                          "surfaces; 6 and 12 are the paper's)"),
        "theta": FloatProperty(
            name="Accuracy (Theta)", default=0.5, min=0.0, max=2.0,
            description="Barnes-Hut opening angle of the energy (smaller: more accurate and "
                        "slower; 0: exact)"),
        "remesh": BoolProperty(
            name="Remesh", default=True,
            description="Split long edges, collapse short ones, flip to Delaunay and smooth "
                        "after every step. The mesh data is replaced by the result (UV maps, "
                        "colours, attributes and vertex groups other than the pin group are "
                        "dropped; material index and smooth flag follow the parent face). Off: "
                        "the triangles stay and positions are written into the mesh"),
        "ccd": BoolProperty(
            name="Collision Detection", default=True,
            description="Continuous collision detection in the line search: the surface never "
                        "passes through itself. Turned off automatically (with a warning) when "
                        "the input already intersects itself"),
        "metric": EnumProperty(name="Metric", items=METRIC_ITEMS, default='SOBOLEV',
                               description="The inner product of the descent"),
        "space": EnumProperty(name="Space", items=SPACE_ITEMS, default='WORLD',
                              description="The coordinates the flow runs in"),
        "pin_group": StringProperty(
            name="Pin Group", default="",
            description="Vertex group (weights above 0) of vertices that never move or vanish"),
        "barycenter": EnumProperty(name="Barycenter", items=BARYCENTER_ITEMS, default='AUTO',
                                   description="Keep the barycenter of the surface fixed"),
        "barycenter_shift": BoolProperty(
            name="Barycenter Shift", default=False,
            description="The paper's mean-motion shift of the barycenter (Sec. 7.1.1; "
                        "unit-inconsistent, for comparisons)"),
        "area_constraint": BoolProperty(
            name="Fix Area", default=False,
            description="Total area is a hard constraint (up to the growth factor)"),
        "area_factor": FloatProperty(
            name="Area Factor", default=1.0, min=0.01, soft_max=10.0,
            description="The area target moves to this factor times the initial area"),
        "area_steps": IntProperty(name="Area Steps", default=100, min=1, soft_max=1000,
                                  description="Steps over which the area target moves"),
        "volume_constraint": BoolProperty(
            name="Fix Volume", default=False,
            description="Enclosed volume is a hard constraint (needs a closed, consistently "
                        "oriented surface)"),
        "volume_factor": FloatProperty(
            name="Volume Factor", default=1.0, min=0.01, soft_max=10.0,
            description="The volume target moves to this factor times the initial volume"),
        "volume_steps": IntProperty(name="Volume Steps", default=100, min=1, soft_max=1000,
                                    description="Steps over which the volume target moves"),
        "obstacle_collection": StringProperty(
            name="Obstacles", default="",
            description="Name of a collection: its mesh objects are tangent-point obstacles"),
        "obstacle_weight": FloatProperty(
            name="Obstacle Weight", default=1.0, min=0.0, soft_max=10.0,
            description="Weight of the tangent-point obstacle energy"),
        "obstacle_ccd": BoolProperty(
            name="Obstacle Collision", default=True,
            description="The line search also keeps the surface out of the obstacles"),
        "implicit_collection": StringProperty(
            name="Empties", default="",
            description="Name of a collection: its empties are implicit obstacles or "
                        "attractors; the display type chooses the shape (Sphere: sphere, Cone: "
                        "cylinder, Circle: torus, others: plane), the size and scale its "
                        "dimensions. Custom properties ddg_mode, ddg_power, ddg_weight and "
                        "ddg_minor_ratio override the defaults below per empty"),
        "implicit_mode": EnumProperty(name="Mode", items=MODE_ITEMS, default='REPEL',
                                      description="Default for empties without ddg_mode"),
        "implicit_power": FloatProperty(
            name="Power", default=2.0, min=0.1, soft_max=8.0,
            description="Default exponent of the potential (d^-p for obstacles, |d|^p for "
                        "attractors)"),
        "implicit_weight": FloatProperty(
            name="Weight", default=1.0, min=0.0, soft_max=10.0,
            description="Default weight of the potential"),
        "implicit_minor_ratio": FloatProperty(
            name="Torus Tube Ratio", default=0.25, min=0.01, max=0.99,
            description="Default minor over major radius of the tori made from Circle empties"),
        "max_steps": IntProperty(name="Step Limit", default=1000, min=1, soft_max=100000,
                                 description="The flow stops after this many steps"),
        "steps": IntProperty(name="Steps", default=20, min=0, soft_max=500,
                             description="Steps to take; 0: until convergence (at most the "
                                         "step limit)"),
        "device": EnumProperty(
            name="Device", items=DEVICE_ITEMS, default='CPU',
            description="Where the heavy products of the flow run. CUDA needs the extension to be "
                        "built with CUDA support and an NVIDIA GPU with its driver"),
        "cuda_ordinal": IntProperty(name="GPU Index", default=0, min=0, soft_max=8,
                                    description="Which CUDA device to use (0: the first)"),
        "gpu_determinism": EnumProperty(
            name="Determinism", items=DETERMINISM_ITEMS, default='RUN_TO_RUN',
            description="How reproducible the GPU result is"),
        "threads": IntProperty(name="Threads", default=0, min=0, soft_max=64,
                               description="Native threads (0: the library's choice)"),
    }
    for name, label in (("w_area_deviation", "Area Deviation"),
                        ("w_volume_deviation", "Volume Deviation"),
                        ("w_boundary_length", "Boundary Length"),
                        ("w_boundary_curvature", "Boundary Curvature"),
                        ("w_willmore", "Willmore"),
                        ("w_squared_error", "Squared Error")):
        props[name] = FloatProperty(
            name=label, default=0.0, min=0.0, soft_max=100.0,
            description="Weight of the {} penalty of the paper's Sec. 7.3 (0: off)".format(
                label.lower()))
    return props


OPTION_NAMES = tuple(option_properties())


class DDG_PG_repulsive_surface_settings(bpy.types.PropertyGroup):
    """The options of the sidebar panel; the live operator and Repulsive Surface Flow start
    from them."""

    __annotations__ = option_properties()


def copy_options(source, target):
    for name in OPTION_NAMES:
        setattr(target, name, getattr(source, name))


# ----------------------------------------------------------------------------------------------
# Building the session


def session_options(settings, surface, penalties, steps=None, ccd=None):
    """Keyword arguments of ``RepulsiveSurfaceSession`` for ``settings`` and ``surface``."""
    steps = settings.steps if steps is None else steps
    kwargs = dict(
        alpha=settings.alpha, beta=settings.beta, theta=settings.theta,
        metric=settings.metric.lower(), remesh=bool(settings.remesh),
        ccd=settings.ccd if ccd is None else ccd,
        max_iterations=max(settings.max_steps, steps),
        barycenter_shift=bool(settings.barycenter_shift),
    )
    if settings.threads > 0:
        kwargs["threads"] = settings.threads
    if settings.device != 'CPU':
        kwargs["device"] = settings.device.lower()
        kwargs["cuda_ordinal"] = settings.cuda_ordinal
        kwargs["gpu_determinism"] = settings.gpu_determinism.lower()
    constraints = []
    if settings.barycenter == 'ON' or (settings.barycenter == 'AUTO' and not surface.pinned):
        constraints.append({"type": "barycenter"})
    if surface.pinned:
        constraints.append({"type": "pins", "vertices": list(surface.pinned)})
    if settings.area_constraint:
        constraints.append({"type": "area", "factor": settings.area_factor,
                            "steps": settings.area_steps})
    if settings.volume_constraint:
        constraints.append({"type": "volume", "factor": settings.volume_factor,
                            "steps": settings.volume_steps})
    kwargs["constraints"] = constraints
    if penalties:
        kwargs["penalties"] = penalties
    return kwargs


def penalty_list(settings, surface, depsgraph):
    """Soft penalties, mesh obstacles and implicit obstacles/attractors as the session's
    ``penalties`` dicts."""
    penalties = []
    for prop, kind in PENALTY_WEIGHTS:
        weight = getattr(settings, prop)
        if weight > 0.0:
            penalties.append({"type": kind, "weight": weight})
    objects = io.objects_of(settings.obstacle_collection, 'MESH')
    for pts, tris in io.read_mesh_obstacles(objects, surface, depsgraph):
        penalties.append({"type": "mesh_obstacle", "vertices": pts, "triangles": tris,
                          "kind": "tangent_point", "weight": settings.obstacle_weight,
                          "ccd": bool(settings.obstacle_ccd)})
    for empty in io.objects_of(settings.implicit_collection, 'EMPTY'):
        penalties.append(io.implicit_penalty(
            empty, surface, mode=settings.implicit_mode, power=settings.implicit_power,
            weight=settings.implicit_weight, minor_ratio=settings.implicit_minor_ratio))
    return penalties


def make_session(obj, settings, depsgraph, steps=None):
    """``(surface, session)`` for ``obj``. Raises ``io.SurfaceInputError`` (message for the user)
    for unusable input, including everything the library rejects."""
    ddg = core.native()
    surface = io.read_object(obj, depsgraph, space=settings.space, pin_group=settings.pin_group)
    try:
        penalties = penalty_list(settings, surface, depsgraph)
        ccd = settings.ccd
        if ccd and io.intersecting_pairs(surface.positions, surface.faces):
            ccd = False
            surface.warnings.append(
                "'{}' intersects itself: collision detection was turned off for this flow "
                "(it would keep the touching triangles from moving)".format(obj.name))
        options = session_options(settings, surface, penalties, steps, ccd)
        try:
            session = ddg.RepulsiveSurfaceSession(surface.positions, surface.faces, **options)
        except ValueError as ex:  # RepulsiveSurfaceInputError is a ValueError
            raise io.SurfaceInputError("{}: {}".format(obj.name, ex)) from ex
    except Exception:
        surface.discard_backup()
        raise
    return surface, session


def status_text(status):
    return STATUS_LABELS.get(status, status)


def summary(stats):
    """One line for the status bar."""
    return "{} steps, energy {:.6g} -> {:.6g}, {} ({:.2f} s native, {:.2f} s Blender), on {}".format(
        stats["steps"], stats["energy_before"], stats["energy_after"],
        status_text(stats["status"]), stats["seconds_native"], stats["seconds_blender"],
        stats.get("device", "CPU"))


def _write_snapshot(surface, snap):
    surface.write(snap["positions"], snap["faces"], snap["pinned"], snap["face_origin"])


# ----------------------------------------------------------------------------------------------
# N steps


def flow_object(obj, settings, depsgraph, *, steps=None, progress=None):
    """Runs the flow on ``obj`` for ``steps`` steps (default ``settings.steps``; 0: until
    convergence, at most ``settings.max_steps``) and writes the result. Returns a dict: ``steps``,
    ``energy_before``, ``energy_after``, ``status``, ``gradient_norm``, ``max_violation``,
    ``warnings``, ``replaced`` (the mesh data was replaced), ``n_vertices``, ``n_faces``,
    ``seconds_native``, ``seconds_blender``. Raises ``io.SurfaceInputError``."""
    steps = settings.steps if steps is None else steps
    t0 = perf_counter()
    surface, session = make_session(obj, settings, depsgraph, steps)
    try:
        energy0 = session.energy
        t1 = perf_counter()
        try:
            report = session.run(max_steps=steps if steps > 0 else None, progress=progress)
        except Exception as ex:
            raise io.SurfaceInputError("{}: the flow failed: {}".format(obj.name, ex)) from ex
        t2 = perf_counter()
        snap = session.snapshot()
        if snap["status"] == "error":
            raise io.SurfaceInputError("{}: the flow failed: {}".format(obj.name, snap["error"]))
        _write_snapshot(surface, snap)
        steps_list = report.get("steps") or []
        last = steps_list[-1] if steps_list else {}
        stats = {
            "steps": report["steps_taken"], "energy_before": energy0,
            "energy_after": snap["energy"], "status": snap["status"],
            "gradient_norm": last.get("gradient_norm"),
            "max_violation": last.get("max_violation"),
            "warnings": list(surface.warnings), "replaced": surface.replaced,
            "n_vertices": len(snap["positions"]) // 3, "n_faces": len(snap["faces"]) // 3,
            "seconds_native": t2 - t1, "seconds_blender": (t1 - t0) + (perf_counter() - t2),
            "description": surface.description,
            "device": device_text(last.get("device"), last.get("device_name"),
                                  last.get("metric_device")),
        }
    except Exception:
        # nothing was written before the failure, except possibly by write itself
        if surface.replaced:
            surface.restore()
        raise
    finally:
        surface.discard_backup()
    return stats


# ----------------------------------------------------------------------------------------------
# The live flow


class LiveFlow:
    """A session stepped on a worker thread and the object it is shown on.

    ``start()`` starts the thread; ``poll()`` (main thread, from a timer) writes the last
    completed step to the object if there is a new one; ``finish()`` stops the thread, writes the
    final state and returns the statistics; ``cancel()`` stops the thread and restores the
    original mesh. Only the main thread touches Blender data; the thread touches the session only.
    """

    def __init__(self, surface, session):
        self.surface = surface
        self.session = session
        self.energy0 = session.energy
        self.settings_device = "CPU"
        self.thread = None
        self.thread_done = False
        self.applied = -1
        self.last_report = None
        self.error = None
        self.blender_seconds = 0.0
        self.native_seconds = 0.0
        self.started = None
        self.stopped = threading.Event()

    @classmethod
    def create(cls, obj, settings, depsgraph):
        surface, session = make_session(obj, settings, depsgraph)
        flow = cls(surface, session)
        flow.settings_device = {'CPU': "CPU", 'CUDA': "CUDA", 'AUTO': "Auto"}[settings.device]
        return flow

    # -- thread

    def _work(self):
        try:
            while not self.stopped.is_set():
                report = self.session.step()
                self.last_report = report
                if report["session_status"] != "running":
                    break
        except BaseException as ex:  # reported by the main thread
            self.error = ex
        finally:
            self.thread_done = True

    def start(self):
        self.thread_done = False
        self.started = perf_counter()
        self.thread = threading.Thread(target=self._work, name="ddg-repulsive-surface",
                                       daemon=True)
        self.thread.start()

    @property
    def running(self):
        return self.thread is not None and not self.thread_done

    def stop(self):
        """Ends the worker thread (it finishes the line-search trial it is in)."""
        self.stopped.set()
        self.session.cancel()
        if self.thread is not None:
            self.thread.join()
            self.native_seconds = perf_counter() - self.started

    # -- main thread

    def poll(self):
        """Writes the last completed step to the object if it is new. Returns whether it did."""
        snap = self.session.snapshot()
        if snap["iteration"] == self.applied:
            return False
        t = perf_counter()
        _write_snapshot(self.surface, snap)
        self.blender_seconds += perf_counter() - t
        self.applied = snap["iteration"]
        return True

    def numbers(self):
        """``(iteration, energy, gradient norm or None, status)`` of the last completed step."""
        report = self.last_report
        gradient = report.get("gradient_norm") if report else None
        return self.session.iteration, self.session.energy, gradient, self.session.status

    def status_line(self):
        iteration, energy, gradient, status = self.numbers()
        text = "Repulsive surface flow: step {}   energy {:.6g}".format(iteration, energy)
        if gradient is not None:
            text += "   gradient {:.3g}".format(gradient)
        text += "   {}".format(status_text(status))
        if self.last_report:
            text += "   on {}".format(self.device_text())
        return text + "   [Enter / right click: accept, Esc: cancel]"

    def finish(self):
        """Stops, writes the final state and returns the statistics (as ``flow_object``)."""
        self.stop()
        if self.error is None:
            self.poll()
        self.surface.discard_backup()
        report = self.last_report or {}
        status = self.session.status
        if status == "cancelled" and self.last_report is not None:
            status = "stopped"
        snap = self.session.snapshot()
        return {
            "steps": self.session.iteration, "energy_before": self.energy0,
            "energy_after": self.session.energy, "status": status,
            "gradient_norm": report.get("gradient_norm"),
            "max_violation": report.get("max_violation"),
            "warnings": list(self.surface.warnings), "replaced": self.surface.replaced,
            "n_vertices": len(snap["positions"]) // 3, "n_faces": len(snap["faces"]) // 3,
            "seconds_native": self.native_seconds, "seconds_blender": self.blender_seconds,
            "description": self.surface.description,
            "device": self.device_text(),
        }

    def device_text(self):
        """'CPU' or 'CUDA (name)': the device of the last completed step (the session's
        setting before the first one)."""
        report = self.last_report
        if report:
            return device_text(report.get("device"), report.get("device_name"),
                               report.get("metric_device"))
        return self.settings_device

    def cancel(self):
        """Stops and puts the mesh back exactly as it was."""
        self.stop()
        self.surface.restore()
        self.surface.discard_backup()


# ----------------------------------------------------------------------------------------------
# Operators


def _poll_object(context):
    obj = context.active_object
    return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}


def draw_options(layout, settings, obj):
    """The options, for the redo panel and the sidebar."""
    layout.use_property_split = True
    layout.use_property_decorate = False
    col = layout.column(heading="Energy")
    col.prop(settings, "alpha")
    col.prop(settings, "beta")
    layout.prop(settings, "theta")
    layout.prop(settings, "remesh")
    layout.prop(settings, "ccd")
    layout.prop(settings, "space")
    box = layout.box()
    box.label(text="Pins")
    if obj is not None and obj.type == 'MESH':
        box.prop_search(settings, "pin_group", obj, "vertex_groups")
    else:
        box.prop(settings, "pin_group")
    box = layout.box()
    box.label(text="Constraints")
    box.prop(settings, "barycenter")
    row = box.row(heading="Area")
    row.prop(settings, "area_constraint", text="")
    sub = row.row(align=True)
    sub.active = settings.area_constraint
    sub.prop(settings, "area_factor", text="x")
    sub.prop(settings, "area_steps", text="steps")
    row = box.row(heading="Volume")
    row.prop(settings, "volume_constraint", text="")
    sub = row.row(align=True)
    sub.active = settings.volume_constraint
    sub.prop(settings, "volume_factor", text="x")
    sub.prop(settings, "volume_steps", text="steps")
    box = layout.box()
    box.label(text="Obstacles")
    box.prop_search(settings, "obstacle_collection", bpy.data, "collections")
    sub = box.column()
    sub.active = bool(settings.obstacle_collection)
    sub.prop(settings, "obstacle_weight")
    sub.prop(settings, "obstacle_ccd")
    box.prop_search(settings, "implicit_collection", bpy.data, "collections")
    sub = box.column()
    sub.active = bool(settings.implicit_collection)
    sub.prop(settings, "implicit_mode")
    sub.prop(settings, "implicit_power")
    sub.prop(settings, "implicit_weight")
    sub.prop(settings, "implicit_minor_ratio")
    layout.prop(settings, "max_steps")
    layout.prop(settings, "threads")


def draw_penalties(layout, settings):
    layout.use_property_split = True
    layout.use_property_decorate = False
    for prop, _kind in PENALTY_WEIGHTS:
        layout.prop(settings, prop)


def draw_advanced(layout, settings):
    layout.use_property_split = True
    layout.use_property_decorate = False
    layout.prop(settings, "metric")
    layout.prop(settings, "barycenter_shift")
    box = layout.box()
    box.label(text="GPU")
    info = gpu_info()
    col = box.column()
    col.enabled = bool(info["available"])
    col.prop(settings, "device")
    sub = col.column()
    sub.active = settings.device != 'CPU'
    sub.prop(settings, "cuda_ordinal")
    sub.prop(settings, "gpu_determinism")
    if info["available"]:
        names = ", ".join("{}: {}".format(d["ordinal"], d["name"]) for d in info["devices"])
        box.label(text=names, icon='INFO')
    else:
        box.label(text="CUDA is not available:", icon='ERROR')
        box.label(text=info["reason"] or "no usable GPU")
        if not info["cuda_compiled"]:
            box.label(text="(Build the extension with DDG_WHEEL_CUDA=1.)")


def _report_warnings(operator, warnings):
    for warning in warnings:
        operator.report({'WARNING'}, warning)


class _FlowProgress:
    def __init__(self, wm, total):
        self.wm = wm
        self.total = max(total, 1)
        wm.progress_begin(0, self.total)

    def __call__(self, report):
        self.wm.progress_update(min(report["iteration"], self.total))
        return True

    def end(self):
        self.wm.progress_end()


class OBJECT_OT_ddg_repulsive_surface_flow(bpy.types.Operator):
    bl_idname = "object.ddg_repulsive_surface_flow"
    bl_label = "Repulsive Surface Flow (Steps)"
    bl_description = (
        "Run N steps of the repulsive-surfaces flow (or until convergence) on the active mesh: "
        "the tangent-point energy decreases, the surface spreads out and untangles. Method: Yu, "
        "Brakensiek, Schumacher, Crane, Repulsive Surfaces (2021). With remeshing on, the mesh "
        "data is replaced (UVs and attributes are dropped). One undo step; Adjust Last "
        "Operation reruns from the original. The window waits until it is done"
    )
    bl_options = {'REGISTER', 'UNDO'}

    __annotations__ = option_properties()

    @classmethod
    def poll(cls, context):
        return _poll_object(context)

    def invoke(self, context, event):
        copy_options(context.scene.ddg_repulsive_surface, self)
        return self.execute(context)

    def draw(self, context):
        self.layout.prop(self, "steps")
        draw_options(self.layout, self, context.active_object)
        header, body = self.layout.panel("ddg_repulsive_surface_penalties", default_closed=True)
        header.label(text="Penalties")
        if body is not None:
            draw_penalties(body, self)
        header, body = self.layout.panel("ddg_repulsive_surface_advanced", default_closed=True)
        header.label(text="Advanced")
        if body is not None:
            draw_advanced(body, self)

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        was_edit_mode = obj.mode == 'EDIT'
        if was_edit_mode:
            bpy.ops.object.mode_set(mode='OBJECT')
        total = self.steps if self.steps > 0 else self.max_steps
        progress = _FlowProgress(context.window_manager, total)
        try:
            stats = flow_object(obj, self, context.evaluated_depsgraph_get(),
                                steps=self.steps, progress=progress)
        except io.SurfaceInputError as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        finally:
            progress.end()
            if was_edit_mode:
                bpy.ops.object.mode_set(mode='EDIT')
        _report_warnings(self, stats["warnings"])
        if stats["status"] == "stuck":
            self.report({'WARNING'}, "The line search found no step: the surface may be at a "
                                     "minimum, or touching an obstacle or itself")
        self.report({'INFO'}, summary(stats))
        return {'FINISHED'}


# The running live flow, for the Accept / Cancel buttons of the panel.
_LIVE = {"flow": None, "request": None}


class OBJECT_OT_ddg_repulsive_surface_live(bpy.types.Operator):
    bl_idname = "object.ddg_repulsive_surface_live"
    bl_label = "Repulsive Surface Flow (Live)"
    bl_description = (
        "Run the repulsive-surfaces flow live in the viewport until stopped: Enter, right click "
        "or the Accept button keep the result (one undo step), Esc or Cancel restore the "
        "original mesh. The options are those of the Repulsive Surfaces panel (sidebar, DDG tab)"
    )
    bl_options = {'REGISTER', 'UNDO'}

    _timer = None
    _flow = None
    _area = None

    @classmethod
    def poll(cls, context):
        return _poll_object(context) and _LIVE["flow"] is None

    def invoke(self, context, event):
        obj = context.active_object
        if obj.mode == 'EDIT':
            bpy.ops.object.mode_set(mode='OBJECT')
        try:
            flow = LiveFlow.create(obj, context.scene.ddg_repulsive_surface,
                                   context.evaluated_depsgraph_get())
        except io.SurfaceInputError as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        _report_warnings(self, flow.surface.warnings)
        self._flow = flow
        _LIVE["flow"], _LIVE["request"] = flow, None
        self._area = context.area
        wm = context.window_manager
        self._timer = wm.event_timer_add(1.0 / 30.0, window=context.window)
        wm.modal_handler_add(self)
        flow.start()
        self._header(context)
        return {'RUNNING_MODAL'}

    def _header(self, context):
        text = self._flow.status_line()
        if self._area is not None:
            self._area.header_text_set(text)
        context.workspace.status_text_set(text)

    def _cleanup(self, context):
        wm = context.window_manager
        if self._timer is not None:
            wm.event_timer_remove(self._timer)
            self._timer = None
        if self._area is not None:
            self._area.header_text_set(None)
            self._area.tag_redraw()
        context.workspace.status_text_set(None)
        _LIVE["flow"], _LIVE["request"] = None, None

    def _accept(self, context):
        flow = self._flow
        try:
            stats = flow.finish()
        finally:
            self._cleanup(context)
        _report_warnings(self, stats["warnings"])
        self.report({'INFO'}, summary(stats))
        return {'FINISHED'}

    def _cancel(self, context):
        try:
            self._flow.cancel()
        finally:
            self._cleanup(context)
        self.report({'INFO'}, "Repulsive surface flow cancelled: the original mesh is restored")
        return {'CANCELLED'}

    def modal(self, context, event):
        flow = self._flow
        request = _LIVE["request"]
        if event.type == 'ESC' and event.value == 'PRESS' or request == 'CANCEL':
            return self._cancel(context)
        if event.value == 'PRESS' and event.type in {'RET', 'NUMPAD_ENTER', 'RIGHTMOUSE'} \
                or request == 'ACCEPT':
            return self._accept(context)
        if event.type == 'TIMER':
            if flow.error is not None:
                message = "Repulsive surface flow failed: {}".format(flow.error)
                flow.cancel()
                self._cleanup(context)
                self.report({'ERROR'}, message + ". The original mesh is restored")
                return {'CANCELLED'}
            flow.poll()
            self._header(context)
            if self._area is not None:
                self._area.tag_redraw()
            return {'RUNNING_MODAL'}
        return {'PASS_THROUGH'}

    def cancel(self, context):
        # Blender ends the operator (window closed, file loaded): leave the mesh as it was.
        if self._flow is not None and _LIVE["flow"] is self._flow:
            try:
                self._flow.cancel()
            finally:
                self._cleanup(context)


class OBJECT_OT_ddg_repulsive_surface_live_finish(bpy.types.Operator):
    bl_idname = "object.ddg_repulsive_surface_live_finish"
    bl_label = "Finish Repulsive Surface Flow"
    bl_description = "Accept or cancel the running live repulsive surface flow"
    bl_options = {'INTERNAL'}

    accept: BoolProperty(name="Accept", default=True)

    @classmethod
    def poll(cls, context):
        return _LIVE["flow"] is not None

    def execute(self, context):
        _LIVE["request"] = 'ACCEPT' if self.accept else 'CANCEL'
        return {'FINISHED'}


class VIEW3D_PT_ddg_repulsive_surface(bpy.types.Panel):
    bl_label = "Repulsive Surfaces"
    bl_idname = "VIEW3D_PT_ddg_repulsive_surface"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "DDG"

    def draw(self, context):
        layout = self.layout
        flow = _LIVE["flow"]
        if flow is not None:
            layout.label(text=flow.status_line().split("   [")[0])
            row = layout.row()
            row.operator(OBJECT_OT_ddg_repulsive_surface_live_finish.bl_idname, text="Accept",
                         icon='CHECKMARK').accept = True
            row.operator(OBJECT_OT_ddg_repulsive_surface_live_finish.bl_idname, text="Cancel",
                         icon='X').accept = False
            return
        layout.operator(OBJECT_OT_ddg_repulsive_surface_live.bl_idname, icon='PLAY')
        row = layout.row(align=True)
        row.operator(OBJECT_OT_ddg_repulsive_surface_flow.bl_idname, icon='FF')
        row.prop(context.scene.ddg_repulsive_surface, "steps", text="Steps")
        draw_options(layout, context.scene.ddg_repulsive_surface, context.active_object)


class VIEW3D_PT_ddg_repulsive_surface_penalties(bpy.types.Panel):
    bl_label = "Penalties"
    bl_idname = "VIEW3D_PT_ddg_repulsive_surface_penalties"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "DDG"
    bl_parent_id = "VIEW3D_PT_ddg_repulsive_surface"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return _LIVE["flow"] is None

    def draw(self, context):
        draw_penalties(self.layout, context.scene.ddg_repulsive_surface)


class VIEW3D_PT_ddg_repulsive_surface_advanced(bpy.types.Panel):
    bl_label = "Advanced"
    bl_idname = "VIEW3D_PT_ddg_repulsive_surface_advanced"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "DDG"
    bl_parent_id = "VIEW3D_PT_ddg_repulsive_surface"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return _LIVE["flow"] is None

    def draw(self, context):
        draw_advanced(self.layout, context.scene.ddg_repulsive_surface)


CLASSES = (
    DDG_PG_repulsive_surface_settings,
    OBJECT_OT_ddg_repulsive_surface_flow,
    OBJECT_OT_ddg_repulsive_surface_live,
    OBJECT_OT_ddg_repulsive_surface_live_finish,
    VIEW3D_PT_ddg_repulsive_surface,
    VIEW3D_PT_ddg_repulsive_surface_penalties,
    VIEW3D_PT_ddg_repulsive_surface_advanced,
)


def _draw_menu_entry(self, _context):
    self.layout.separator()
    self.layout.operator(OBJECT_OT_ddg_repulsive_surface_live.bl_idname)
    self.layout.operator(OBJECT_OT_ddg_repulsive_surface_flow.bl_idname)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.Scene.ddg_repulsive_surface = PointerProperty(type=DDG_PG_repulsive_surface_settings)
    bpy.types.VIEW3D_MT_object.append(_draw_menu_entry)


def unregister():
    bpy.types.VIEW3D_MT_object.remove(_draw_menu_entry)
    del bpy.types.Scene.ddg_repulsive_surface
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
