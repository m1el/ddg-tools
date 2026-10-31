# SPDX-License-Identifier: MIT
"""Repulsive Curves operators: a live (modal) flow and an undoable N-steps flow.

The work is done by ``ddg_tools.RepulsiveSession`` (Rust, crates/ddg-repulsive-curves and
crates/ddg-py), after

    Chris Yu, Henrik Schumacher, Keenan Crane, Repulsive Curves,
    ACM Transactions on Graphics 40(2), 2021.

The flow decreases the tangent-point energy of a curve network (mesh edges or Curve splines;
``repulsive_io.py``) by fractional Sobolev descent, under length, pin, planarity and obstacle
constraints. The flow runs in world space (default) or in the object's local space ("Space");
obstacles are transformed into the same space.

* ``object.ddg_repulsive_flow``: N steps (or to convergence) in one go, one undo step, with a
  redo panel.
* ``object.ddg_repulsive_live``: a worker thread steps the session while a timer (1/30 s) copies
  the last completed step into the object. Enter, right click, or the Accept button of the
  sidebar panel (``object.ddg_repulsive_live_finish``) keeps the result as one undo step; Esc or
  the Cancel button restores the original object exactly. The options are those of the sidebar
  panel (``scene.ddg_repulsive``).

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
from . import repulsive_io as io

LENGTH_ITEMS = [
    ('EDGES', "Keep Edge Lengths", "Every edge keeps its length (up to the growth factor)"),
    ('TOTAL', "Keep Total Length", "The total length of the network is kept; edges may change"),
    ('FREE', "Free", "No length constraint: the curve may grow or shrink by itself"),
]
SPACE_ITEMS = [
    ('WORLD', "World", "Run in world space: the object's transform (scale!) counts as it would "
     "physically, obstacles are used as they stand"),
    ('LOCAL', "Local", "Run in the object's local space; obstacles are transformed into it"),
]
AXIS_ITEMS = [('X', "X", ""), ('Y', "Y", ""), ('Z', "Z", "")]
LINE_SEARCH_ITEMS = [
    ('collision', "Collision Guarded", "Backtracking line search with continuous collision "
     "detection: the curve never passes through itself or an obstacle"),
    ('armijo', "Armijo", "Plain backtracking on the energy; no collision check"),
]
BACKEND_ITEMS = [
    ('auto', "Auto", "Exact energy and dense solves up to about 2500 vertices, Barnes-Hut energy "
     "and multigrid solves above"),
    ('dense', "Dense", "Exact energy and dense solves"),
    ('barnes_hut', "Barnes-Hut", "Barnes-Hut energy and dense solves"),
    ('multigrid', "Multigrid", "Exact energy and multigrid solves (a step whose multigrid solve "
     "fails is repeated with the dense solver on curves of up to 5000 vertices)"),
    ('accelerated', "Accelerated", "Barnes-Hut energy and multigrid solves"),
]

METHOD_ITEMS = [
    ('PROJECTED_GRADIENT', "Projected Gradient", "The paper's flow: fractional Sobolev descent "
     "with the lengths (and pins) as hard constraints"),
    ('LBFGS', "L-BFGS", "L-BFGS on the paper's soft problem: no hard constraints, the edge "
     "lengths are a Hencky penalty. No pins or tangent constraints, no growth"),
    ('NCG', "Nonlinear CG", "Polak-Ribiere nonlinear conjugate gradients on the same soft "
     "problem as L-BFGS"),
]
SECOND_NEIGHBOUR_ITEMS = [
    ('FOLD_GATED', "Fold Gated", "Edges two apart are corrected near contact, and in folds "
     "(the default)"),
    ('ONE_SIDED', "One Sided", "Edges two apart are corrected but never below the paper's term"),
    ('TWO_SIDED', "Two Sided", "Edges two apart are corrected like every other pair (lowers "
     "the energy at corners)"),
]

STATUS_LABELS = {
    "running": "running",
    "converged": "converged",
    "max_steps": "step limit reached",
    "time_limit": "time limit reached",
    "stuck": "stuck (the line search found no step)",
    "cancelled": "stopped",
    "error": "error",
}


def option_properties():
    """The options shared by the operators and the scene settings (fresh property objects)."""
    return {
        "alpha": FloatProperty(name="Alpha", default=3.0, min=1.0, soft_max=10.0,
                               description="Exponent alpha of the tangent-point kernel"),
        "beta": FloatProperty(name="Beta", default=6.0, min=2.0, soft_max=12.0,
                              description="Exponent beta of the tangent-point kernel (beta > "
                                          "alpha; 3 and 6 are the paper's)"),
        "method": EnumProperty(name="Method", items=METHOD_ITEMS, default='PROJECTED_GRADIENT',
                               description="The descent method"),
        "hencky_weight": FloatProperty(
            name="Length Penalty", default=10.0, min=0.0, soft_max=100.0,
            description="Weight of the Hencky penalty that keeps the edge lengths of the soft "
                        "methods (the paper's evaluation uses 10)"),
        "length": EnumProperty(name="Length", items=LENGTH_ITEMS, default='EDGES',
                               description="What the flow keeps constant"),
        "grow_factor": FloatProperty(
            name="Grow/Shrink Factor", default=1.0, min=0.01, soft_max=4.0,
            description="The kept length moves to this factor times itself (1: off)"),
        "grow_rate": FloatProperty(
            name="Rate", default=0.01, min=0.0001, max=1.0,
            description="Fraction of the start length by which the target moves per step"),
        "subdivide": BoolProperty(
            name="Live Subdivision", default=False,
            description="Split every edge when the mean edge length has doubled. Off: the "
                        "vertices, pins and vertex groups stay as they are. On: old vertices "
                        "keep their indices, vertex groups are reapplied (a new vertex joins a "
                        "group if both neighbours it was inserted between are in it)"),
        "max_subdivisions": IntProperty(name="Subdivisions", default=2, min=1, soft_max=6,
                                        description="At most this many subdivisions"),
        "pin_group": StringProperty(
            name="Pin Group", default="",
            description="Mesh: vertex group (weights above 0) of vertices that do not move. "
                        "Curve objects: the selected control points are pinned"),
        "tangent_group": StringProperty(
            name="Tangent Group", default="",
            description="Mesh: vertex group of vertices whose tangent is kept"),
        "pin_tangents": BoolProperty(
            name="Keep Tangents at Pins", default=False,
            description="The pinned vertices also keep their tangent (do not pin neighbours "
                        "this way: the constraints would be dependent)"),
        "ignore_loose": BoolProperty(
            name="Ignore Loose Vertices", default=False,
            description="Mesh: leave out vertices that belong to no edge (otherwise they are "
                        "an error)"),
        "loose_edges_only": BoolProperty(
            name="Loose Edges Only", default=False,
            description="Mesh with faces: use only the edges that belong to no face (the rest "
                        "stays). Off: a mesh with faces is rejected"),
        "obstacle_object": StringProperty(
            name="Obstacle", default="",
            description="Name of a mesh object the curve may not cross (picked from the list)"),
        "obstacle_collection": StringProperty(
            name="Obstacles", default="",
            description="Name of a collection: its mesh objects are obstacles"),
        "obstacle_weight": FloatProperty(
            name="Obstacle Repulsion", default=0.0, min=0.0, soft_max=10.0,
            description="0: obstacles are only a collision check. Above 0: they also repel "
                        "the curve (a repulsive potential)"),
        "planar": BoolProperty(name="Planar", default=False,
                               description="Hold the coordinate of one axis fixed (flat drawings)"),
        "planar_axis": EnumProperty(name="Axis", items=AXIS_ITEMS, default='Z',
                                    description="The axis whose coordinate stays; of the world "
                                                "or the local axes, as Space says"),
        "space": EnumProperty(name="Space", items=SPACE_ITEMS, default='WORLD',
                              description="The coordinates the flow runs in"),
        "line_search": EnumProperty(name="Line Search", items=LINE_SEARCH_ITEMS,
                                    default='collision'),
        "backend": EnumProperty(name="Back End", items=BACKEND_ITEMS, default='auto',
                                description="How the energy and the solves are computed"),
        "threads": IntProperty(name="Threads", default=0, min=0, soft_max=64,
                               description="Native threads (0: the library's choice)"),
        "near_pair": BoolProperty(
            name="Near-Pair Quadrature", default=True,
            description="Accurate energy for edges that come close in their interiors (not in "
                        "the paper). Off: the paper's energy"),
        "second_neighbour_rule": EnumProperty(
            name="Second Neighbours", items=SECOND_NEIGHBOUR_ITEMS, default='FOLD_GATED',
            description="How edges with one edge between them are treated by the near-pair "
                        "quadrature"),
        "stuck_retry_chord": BoolProperty(
            name="Retry Stuck Steps with Chord", default=True,
            description="A stuck line search is repeated with the paper's chord projection"),
        "back_projection": BoolProperty(
            name="Back Projection", default=True,
            description="Project the trial points back onto the constraints (the paper). Off: "
                        "the paper's 'without projection' variant"),
        "newton_polish": BoolProperty(
            name="Newton Polish", default=True,
            description="Finish the projection with Newton steps once the violation is small "
                        "(not in the paper). Off: the paper's chord method throughout"),
        "newton_polish_threshold": FloatProperty(
            name="Polish Below", default=100.0, min=1.0, soft_max=1000.0,
            description="Newton steps start when the violation is at most this many times the "
                        "projection tolerance"),
        "projection_early_abort": BoolProperty(
            name="Early Abort", default=True,
            description="Give up a chord projection early when it cannot converge in time"),
        "second_order_correction": BoolProperty(
            name="Second-Order Correction", default=False,
            description="Correct the trial points to second order before projecting"),
        "mg_max_threads": IntProperty(
            name="Multigrid Max Threads", default=16, min=0, soft_max=64,
            description="Cap of the multigrid solver's threads (0: no cap; "
                        "16: the library's default)"),
        "stop_tolerance": FloatProperty(
            name="Stop Tolerance", default=1e-4, min=0.0, precision=6,
            description="The flow has converged when the gradient norm falls below this"),
        "max_steps": IntProperty(name="Step Limit", default=1000, min=1, soft_max=100000,
                                 description="The flow stops after this many steps"),
        "steps": IntProperty(name="Steps", default=50, min=0, soft_max=1000,
                             description="Steps to take; 0: until convergence (at most the "
                                         "step limit)"),
        "bezier_samples": IntProperty(
            name="Bezier Samples", default=1, min=1, soft_max=16,
            description="Points per Bezier segment. 1: the control points are the polyline "
                        "(the point count is kept). More: the curve is sampled, and writing "
                        "back creates that many points"),
        "keep_bezier": BoolProperty(
            name="Keep Bezier", default=False,
            description="Write Bezier splines back as Bezier splines, re-fitted through the "
                        "new points (approximate: Catmull-Rom tangents, aligned handles), and "
                        "NURBS splines as NURBS with moved control points. Off: poly splines "
                        "(exact)"),
    }


OPTION_NAMES = tuple(option_properties())


class DDG_PG_repulsive_settings(bpy.types.PropertyGroup):
    """The options of the sidebar panel; the live operator and Repulsive Flow start from them."""

    __annotations__ = option_properties()


def copy_options(source, target):
    for name in OPTION_NAMES:
        setattr(target, name, getattr(source, name))


# ----------------------------------------------------------------------------------------------
# Building the session


def session_options(settings, net, obstacles, steps=None):
    """Keyword arguments of ``RepulsiveSession`` for ``settings`` and the network ``net``."""
    steps = settings.steps if steps is None else steps
    kwargs = dict(
        alpha=settings.alpha, beta=settings.beta, length=settings.length.lower(),
        grow_factor=settings.grow_factor, grow_rate=settings.grow_rate,
        max_subdivisions=settings.max_subdivisions if settings.subdivide else 0,
        line_search=settings.line_search, backend=settings.backend,
        stop="gradient", stop_threshold=settings.stop_tolerance,
        max_steps=steps if steps > 0 else settings.max_steps,
    )
    if settings.threads > 0:
        kwargs["threads"] = settings.threads
    kwargs.update(
        method=settings.method.lower(), hencky_weight=settings.hencky_weight,
        near_pair=settings.near_pair,
        second_neighbour_rule=settings.second_neighbour_rule.lower(),
        stuck_retry_chord=settings.stuck_retry_chord,
        back_projection=settings.back_projection,
        newton_polish=settings.newton_polish_threshold if settings.newton_polish else None,
        projection_early_abort=settings.projection_early_abort,
        second_order_correction=settings.second_order_correction,
        mg_max_threads=settings.mg_max_threads if settings.mg_max_threads > 0 else None)
    if settings.planar:
        kwargs["planar_axis"] = settings.planar_axis.lower()
    if obstacles:
        kwargs["obstacles"] = obstacles
        kwargs["obstacle_weight"] = settings.obstacle_weight
    kwargs.update(io.session_arguments(net))
    return kwargs


def read_network(obj, settings):
    """The network of ``obj`` for ``settings``; raises ``io.CurveInputError``."""
    return io.read_object(
        obj, space=settings.space, pin_group=settings.pin_group,
        tangent_group=settings.tangent_group, pin_tangents=settings.pin_tangents,
        ignore_loose=settings.ignore_loose, loose_edges_only=settings.loose_edges_only,
        bezier_samples=settings.bezier_samples, keep_bezier=settings.keep_bezier,
        keep_groups=settings.subdivide, keep_backup=settings.subdivide)


def check_method(settings, net):
    """A soft method has no hard constraints: pins, tangents and growth are refused."""
    if settings.method == 'PROJECTED_GRADIENT':
        return
    label = dict((i[0], i[1]) for i in METHOD_ITEMS)[settings.method]
    if net.pinned or net.tangent:
        raise io.CurveInputError(
            "{}: pins and tangent constraints need Method: Projected Gradient ({} solves the "
            "soft problem without hard constraints)".format(net.obj.name, label))
    if settings.grow_factor != 1.0:
        raise io.CurveInputError(
            "Grow/Shrink Factor needs Method: Projected Gradient ({} has no hard length "
            "constraint to move)".format(label))


def make_session(obj, settings, depsgraph, steps=None):
    """``(network, session)`` for ``obj``. Raises ``io.CurveInputError`` (message for the user)
    for unusable input, including everything the library rejects."""
    ddg = core.native()
    net = read_network(obj, settings)
    try:
        if settings.subdivide and not net.whole_mesh and net.kind == 'MESH':
            raise io.CurveInputError(
                "Live Subdivision needs the whole mesh to be the curve network (no faces, no "
                "ignored loose vertices); turn it off or clean the mesh")
        check_method(settings, net)
        objects = io.obstacle_objects(settings.obstacle_object, settings.obstacle_collection)
        obstacles = io.read_obstacles(objects, net, depsgraph)
        options = session_options(settings, net, obstacles, steps)
        try:
            session = ddg.RepulsiveSession(net.positions, net.edges, **options)
        except ValueError as ex:  # RepulsiveInputError is a ValueError
            raise io.CurveInputError("{}: {}".format(net.obj.name, ex)) from ex
    except Exception:
        net.discard_backup()
        raise
    return net, session


def status_text(status):
    return STATUS_LABELS.get(status, status)


def summary(stats):
    """One line for the status bar."""
    return "{} steps, energy {:.6g} -> {:.6g}, {} ({:.2f} s native, {:.2f} s Blender)".format(
        stats["steps"], stats["energy_before"], stats["energy_after"],
        status_text(stats["status"]), stats["seconds_native"], stats["seconds_blender"])


# ----------------------------------------------------------------------------------------------
# N steps


def flow_object(obj, settings, depsgraph, *, steps=None, progress=None):
    """Runs the flow on ``obj`` for ``steps`` steps (default ``settings.steps``; 0: until
    convergence, at most ``settings.max_steps``) and writes the result. Returns a dict: ``steps``,
    ``energy_before``, ``energy_after``, ``status``, ``gradient_norm``, ``max_violation``,
    ``warnings``, ``seconds_native``, ``seconds_blender``, ``subdivisions``. Raises
    ``io.CurveInputError``."""
    t0 = perf_counter()
    net, session = make_session(obj, settings, depsgraph, steps)
    try:
        energy0 = session.energy
        t1 = perf_counter()
        report = session.run(progress=progress)
        t2 = perf_counter()
        net.write(session.positions, session.edges,
                  session.vertex_source if session.n_vertices != net.n_vertices else None)
        steps_list = report.get("steps") or []
        last = steps_list[-1] if steps_list else {}
        stats = {
            "steps": report["steps_taken"], "energy_before": energy0,
            "energy_after": session.energy, "status": session.status,
            "gradient_norm": last.get("gradient_norm"),
            "max_violation": last.get("max_violation"),
            "warnings": list(net.warnings), "subdivisions": session.subdivisions,
            "seconds_native": t2 - t1, "seconds_blender": (t1 - t0) + (perf_counter() - t2),
            "description": net.description,
        }
    finally:
        net.discard_backup()
    return stats


# ----------------------------------------------------------------------------------------------
# The live flow


class LiveFlow:
    """A session stepped on a worker thread and the object it is shown on.

    ``start()`` starts the thread; ``poll()`` (main thread, from a timer) writes the last
    completed step to the object if there is a new one; ``finish()`` stops the thread, writes the
    final state and returns the statistics; ``cancel()`` stops the thread and restores the
    original object. Only the main thread touches Blender data; the thread touches the session
    only.
    """

    def __init__(self, net, session):
        self.net = net
        self.session = session
        self.energy0 = session.energy
        self.thread = None
        self.applied = -1
        self.last_report = None
        self.error = None
        self.blender_seconds = 0.0
        self.native_seconds = 0.0
        self.started = None
        self.stopped = threading.Event()

    @classmethod
    def create(cls, obj, settings, depsgraph):
        net, session = make_session(obj, settings, depsgraph)
        return cls(net, session)

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
        self.thread = threading.Thread(target=self._work, name="ddg-repulsive", daemon=True)
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
        positions, edges, iteration, _energy, _status = self.session.snapshot()
        if iteration == self.applied:
            return False
        n = len(positions) // 3
        source = None
        if n != self.net.n_vertices or len(edges) // 2 != self.net.n_edges:
            source = self.session.vertex_source
            if len(source) != n:  # a step finished between the two reads; try again
                return False
        t = perf_counter()
        self.net.write(positions, edges, source)
        self.blender_seconds += perf_counter() - t
        self.applied = iteration
        return True

    def numbers(self):
        """``(iteration, energy, gradient norm or None, status)`` of the last completed step."""
        report = self.last_report
        gradient = report.get("gradient_norm") if report else None
        return self.session.iteration, self.session.energy, gradient, self.session.status

    def status_line(self):
        iteration, energy, gradient, status = self.numbers()
        text = "Repulsive flow: step {}   energy {:.6g}".format(iteration, energy)
        if gradient is not None:
            text += "   gradient {:.3g}".format(gradient)
        text += "   {}".format(status_text(status))
        return text + "   [Enter / right click: accept, Esc: cancel]"

    def finish(self):
        """Stops, writes the final state and returns the statistics (as ``flow_object``)."""
        self.stop()
        if self.error is None:
            self.poll()
        self.net.discard_backup()
        report = self.last_report or {}
        status = self.session.status
        if status == "cancelled" and self.last_report is not None:
            status = "stopped"
        return {
            "steps": self.session.iteration, "energy_before": self.energy0,
            "energy_after": self.session.energy, "status": status,
            "gradient_norm": report.get("gradient_norm"),
            "max_violation": report.get("max_violation"),
            "warnings": list(self.net.warnings), "subdivisions": self.session.subdivisions,
            "seconds_native": self.native_seconds, "seconds_blender": self.blender_seconds,
            "description": self.net.description,
        }

    def cancel(self):
        """Stops and puts the object back exactly as it was."""
        self.stop()
        self.net.restore()
        self.net.discard_backup()


# ----------------------------------------------------------------------------------------------
# Operators


def _poll_object(context):
    obj = context.active_object
    return obj is not None and obj.type in {'MESH', 'CURVE'} and obj.mode in {'OBJECT', 'EDIT'}


def draw_options(layout, settings, obj):
    """The options, for the redo panel and the sidebar."""
    layout.use_property_split = True
    layout.use_property_decorate = False
    is_mesh = obj is not None and obj.type == 'MESH'
    is_curve = obj is not None and obj.type == 'CURVE'
    col = layout.column(heading="Energy")
    col.prop(settings, "alpha")
    col.prop(settings, "beta")
    soft = settings.method != 'PROJECTED_GRADIENT'
    layout.prop(settings, "method")
    if soft:
        layout.prop(settings, "hencky_weight")
    sub = layout.column()
    sub.active = not soft
    sub.prop(settings, "length")
    sub = layout.column()
    sub.active = settings.length != 'FREE' and not soft
    sub.prop(settings, "grow_factor")
    sub.prop(settings, "grow_rate")
    layout.prop(settings, "space")
    box = layout.box()
    box.active = not soft
    box.label(text="Pins" if not soft else "Pins (Projected Gradient only)")
    if is_mesh:
        box.prop_search(settings, "pin_group", obj, "vertex_groups")
        box.prop_search(settings, "tangent_group", obj, "vertex_groups")
        box.prop(settings, "ignore_loose")
        box.prop(settings, "loose_edges_only")
    elif is_curve:
        box.label(text="Selected control points are pinned", icon='INFO')
        box.prop(settings, "bezier_samples")
        box.prop(settings, "keep_bezier")
    else:
        box.prop(settings, "pin_group")
    box.prop(settings, "pin_tangents")
    box = layout.box()
    box.label(text="Obstacles")
    box.prop_search(settings, "obstacle_object", bpy.data, "objects")
    box.prop_search(settings, "obstacle_collection", bpy.data, "collections")
    box.prop(settings, "obstacle_weight")
    row = layout.row(heading="Planar")
    row.prop(settings, "planar", text="")
    sub = row.row()
    sub.active = settings.planar
    sub.prop(settings, "planar_axis", text="")
    row = layout.row(heading="Live Subdivision")
    row.prop(settings, "subdivide", text="")
    sub = row.row()
    sub.active = settings.subdivide
    sub.prop(settings, "max_subdivisions", text="")
    layout.prop(settings, "stop_tolerance")
    layout.prop(settings, "max_steps")
    layout.prop(settings, "line_search")
    layout.prop(settings, "backend")
    layout.prop(settings, "threads")


def draw_advanced(layout, settings):
    """The advanced options (a sub-panel of the sidebar panel, a closed panel in the redo
    panel)."""
    layout.use_property_split = True
    layout.use_property_decorate = False
    layout.prop(settings, "near_pair")
    sub = layout.column()
    sub.active = settings.near_pair
    sub.prop(settings, "second_neighbour_rule")
    layout.prop(settings, "stuck_retry_chord")
    hard = settings.method == 'PROJECTED_GRADIENT'
    col = layout.column()
    col.active = hard
    col.prop(settings, "back_projection")
    row = col.row(heading="Newton Polish")
    row.prop(settings, "newton_polish", text="")
    sub = row.row()
    sub.active = settings.newton_polish
    sub.prop(settings, "newton_polish_threshold", text="")
    col.prop(settings, "projection_early_abort")
    col.prop(settings, "second_order_correction")
    layout.prop(settings, "mg_max_threads")


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


class OBJECT_OT_ddg_repulsive_flow(bpy.types.Operator):
    bl_idname = "object.ddg_repulsive_flow"
    bl_label = "Repulsive Flow (Steps)"
    bl_description = (
        "Run N steps of the repulsive-curves flow (or until convergence) on the active edge "
        "mesh or Curve object: the tangent-point energy decreases, the curve spreads out and "
        "untangles. Method: Yu, Schumacher, Crane, Repulsive Curves (2021). One undo step; "
        "Adjust Last Operation reruns from the original. The window waits until it is done"
    )
    bl_options = {'REGISTER', 'UNDO'}

    __annotations__ = option_properties()

    @classmethod
    def poll(cls, context):
        return _poll_object(context)

    def invoke(self, context, event):
        copy_options(context.scene.ddg_repulsive, self)
        return self.execute(context)

    def draw(self, context):
        self.layout.prop(self, "steps")
        draw_options(self.layout, self, context.active_object)
        header, body = self.layout.panel("ddg_repulsive_advanced", default_closed=True)
        header.label(text="Advanced")
        if body is not None:
            draw_advanced(body, self)

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type not in {'MESH', 'CURVE'}:
            self.report({'ERROR'}, "The active object is not an edge mesh or a Curve")
            return {'CANCELLED'}
        was_edit_mode = obj.mode == 'EDIT'
        if was_edit_mode:
            bpy.ops.object.mode_set(mode='OBJECT')
        total = self.steps if self.steps > 0 else self.max_steps
        progress = _FlowProgress(context.window_manager, total)
        try:
            stats = flow_object(obj, self, context.evaluated_depsgraph_get(),
                                steps=self.steps, progress=progress)
        except io.CurveInputError as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        finally:
            progress.end()
            if was_edit_mode:
                bpy.ops.object.mode_set(mode='EDIT')
        _report_warnings(self, stats["warnings"])
        if stats["status"] == "stuck":
            self.report({'WARNING'}, "The line search found no step: the curve may be at a "
                                     "minimum, or touching an obstacle or itself")
        self.report({'INFO'}, summary(stats))
        return {'FINISHED'}


# The running live flow, for the Accept / Cancel buttons of the panel.
_LIVE = {"flow": None, "request": None}


class OBJECT_OT_ddg_repulsive_live(bpy.types.Operator):
    bl_idname = "object.ddg_repulsive_live"
    bl_label = "Repulsive Flow (Live)"
    bl_description = (
        "Run the repulsive-curves flow live in the viewport until stopped: Enter, right click "
        "or the Accept button keep the result (one undo step), Esc or Cancel restore the "
        "original curve. The options are those of the Repulsive Curves panel (sidebar, DDG tab)"
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
            flow = LiveFlow.create(obj, context.scene.ddg_repulsive,
                                   context.evaluated_depsgraph_get())
        except io.CurveInputError as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        _report_warnings(self, flow.net.warnings)
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
        self.report({'INFO'}, summary(stats))
        return {'FINISHED'}

    def _cancel(self, context):
        try:
            self._flow.cancel()
        finally:
            self._cleanup(context)
        self.report({'INFO'}, "Repulsive flow cancelled: the original curve is restored")
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
                message = "Repulsive flow failed: {}".format(flow.error)
                flow.cancel()
                self._cleanup(context)
                self.report({'ERROR'}, message + ". The original curve is restored")
                return {'CANCELLED'}
            flow.poll()
            self._header(context)
            if self._area is not None:
                self._area.tag_redraw()
            return {'RUNNING_MODAL'}
        return {'PASS_THROUGH'}

    def cancel(self, context):
        # Blender ends the operator (window closed, file loaded): leave the object as it was.
        if self._flow is not None and _LIVE["flow"] is self._flow:
            try:
                self._flow.cancel()
            finally:
                self._cleanup(context)


class OBJECT_OT_ddg_repulsive_live_finish(bpy.types.Operator):
    bl_idname = "object.ddg_repulsive_live_finish"
    bl_label = "Finish Repulsive Flow"
    bl_description = "Accept or cancel the running live repulsive flow"
    bl_options = {'INTERNAL'}

    accept: BoolProperty(name="Accept", default=True)

    @classmethod
    def poll(cls, context):
        return _LIVE["flow"] is not None

    def execute(self, context):
        _LIVE["request"] = 'ACCEPT' if self.accept else 'CANCEL'
        return {'FINISHED'}


class VIEW3D_PT_ddg_repulsive(bpy.types.Panel):
    bl_label = "Repulsive Curves"
    bl_idname = "VIEW3D_PT_ddg_repulsive"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "DDG"

    def draw(self, context):
        layout = self.layout
        flow = _LIVE["flow"]
        if flow is not None:
            layout.label(text=flow.status_line().split("   [")[0])
            row = layout.row()
            row.operator(OBJECT_OT_ddg_repulsive_live_finish.bl_idname, text="Accept",
                         icon='CHECKMARK').accept = True
            row.operator(OBJECT_OT_ddg_repulsive_live_finish.bl_idname, text="Cancel",
                         icon='X').accept = False
            return
        layout.operator(OBJECT_OT_ddg_repulsive_live.bl_idname, icon='PLAY')
        row = layout.row(align=True)
        row.operator(OBJECT_OT_ddg_repulsive_flow.bl_idname, icon='FF')
        row.prop(context.scene.ddg_repulsive, "steps", text="Steps")
        draw_options(layout, context.scene.ddg_repulsive, context.active_object)


class VIEW3D_PT_ddg_repulsive_advanced(bpy.types.Panel):
    bl_label = "Advanced"
    bl_idname = "VIEW3D_PT_ddg_repulsive_advanced"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "DDG"
    bl_parent_id = "VIEW3D_PT_ddg_repulsive"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return _LIVE["flow"] is None

    def draw(self, context):
        draw_advanced(self.layout, context.scene.ddg_repulsive)


CLASSES = (
    DDG_PG_repulsive_settings,
    OBJECT_OT_ddg_repulsive_flow,
    OBJECT_OT_ddg_repulsive_live,
    OBJECT_OT_ddg_repulsive_live_finish,
    VIEW3D_PT_ddg_repulsive,
    VIEW3D_PT_ddg_repulsive_advanced,
)


def _draw_menu_entry(self, _context):
    self.layout.separator()
    self.layout.operator(OBJECT_OT_ddg_repulsive_live.bl_idname)
    self.layout.operator(OBJECT_OT_ddg_repulsive_flow.bl_idname)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.Scene.ddg_repulsive = PointerProperty(type=DDG_PG_repulsive_settings)
    bpy.types.VIEW3D_MT_object.append(_draw_menu_entry)


def unregister():
    bpy.types.VIEW3D_MT_object.remove(_draw_menu_entry)
    del bpy.types.Scene.ddg_repulsive
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
