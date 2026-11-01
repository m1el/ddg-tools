# SPDX-License-Identifier: MIT
"""Geodesic Path and Straighten Seams: geodesic paths, loops, Bezier curves and straightened
seams on the surface of a mesh, by flipping edges of an intrinsic triangulation (Sharp, Crane,
You Can Find Geodesic Paths in Triangle Meshes by Just Flipping Edges, ACM TOG 39(6), 2020).

* ``mesh.ddg_geodesic_path`` (Edit mode): the selected vertices, in the order they were selected
  (``bmesh.select_history``), are joined by a geodesic. Modes: *Path* (through the points),
  *Smooth* (the points only seed the path between the first and the last), *Loop* (closed through
  the points), *Bezier* (the points are control points).
* ``mesh.ddg_straighten_seams``: the edges marked as seams are shortened as a network (junctions
  stay; selected vertices are pinned too).

The result is a Curve object ``<mesh>_geodesic`` (poly splines, one per path), parented to the
mesh with its transform, with a round bevel of 0.2 % of the bounding-box diagonal and the report
in custom properties (``ddg_*``). Running again replaces it. With *Show triangulation* the edges
of the final intrinsic triangulation are drawn on the input mesh as a second Curve object
``<mesh>_geodesic_triangulation`` (``intrinsic_draw``). The mesh is never changed.

The default representation is integer coordinates, as for every other consumer, since they passed
their evaluation (Integer Coordinates SPEC Q3); ``REPRESENTATION_DEFAULT`` names it. Signposts stay
an option. Nothing falls back silently: when signposts fail to trace a path on the input mesh the
message says to choose Integer Coordinates.
"""

import math
from array import array
from time import perf_counter

import bmesh
import bpy
from bpy.props import BoolProperty, EnumProperty, FloatProperty, IntProperty
from mathutils import Matrix

from . import core, intrinsic_draw
from .intrinsic_operator import triangulated_buffers

PREFIX = "ddg_"
CURVE_SUFFIX = "_geodesic"
TRIANGULATION_SUFFIX = "_geodesic_triangulation"
KIND = "geodesic"
KIND_TRIANGULATION = "geodesic_triangulation"
BEVEL_FRACTION = 0.002
MOLLIFY = 1e-4

# The one place that names the default representation of these operators (see the module
# documentation).
REPRESENTATION_DEFAULT = 'INTEGER'
REPRESENTATION_ITEMS = (
    ('SIGNPOSTS', "Signposts", "Sharp, Soliman, Crane 2019: edge lengths and signpost angles"),
    ('INTEGER', "Integer Coordinates",
     "Gillespie, Sharp, Crane 2021: exact combinatorics; the fallback when signposts cannot "
     "trace a path on the mesh"),
)
MODE_ITEMS = (
    ('PATH', "Path", "A geodesic through the selected vertices, in selection order"),
    ('SMOOTH', "Smooth", "The selected vertices only seed the path: one geodesic between the "
     "first and the last"),
    ('LOOP', "Loop", "A closed geodesic through the selected vertices"),
    ('BEZIER', "Bezier", "The selected vertices are the control points of a geodesic Bezier "
     "curve (Morera et al. 2008)"),
)


class Result:
    """What ``geodesic_curve`` / ``straighten_seams`` made."""

    def __init__(self):
        self.curve_object = None
        self.triangulation_object = None
        self.report = {}
        self.lengths = []
        self.seconds_native = 0.0
        self.seconds_blender = 0.0


# ----------------------------------------------------------------------------------------------
# Selection (Edit mode)


def selected_vertices_in_order(obj):
    """The selected vertices of the mesh in Edit mode, in the order they were selected.

    Returns ``(vertices, complete)``; ``complete`` is False when the selection history does not
    list every selected vertex (box select, edge or face select mode), in which case the
    vertices are in index order. Raises ``core.UnwrapFailed`` outside Edit mode.
    """
    if obj is None or obj.type != 'MESH' or obj.mode != 'EDIT':
        raise core.UnwrapFailed("Geodesic Path works in Edit mode: select the vertices there")
    bm = bmesh.from_edit_mesh(obj.data)
    selected = [v.index for v in bm.verts if v.select]
    history = []
    for element in bm.select_history:
        if isinstance(element, bmesh.types.BMVert) and element.select \
                and element.index not in history:
            history.append(element.index)
    if sorted(history) == selected:
        return history, True
    return selected, False


def ordered_selection(obj, minimum=2):
    """The selection in selection order for the user's message rules (module documentation of
    the SPEC): at least ``minimum`` vertices; with more than two and no complete history, an
    error that asks to shift-click them in order."""
    vertices, complete = selected_vertices_in_order(obj)
    if len(vertices) < minimum:
        raise core.UnwrapFailed(
            "Select at least {} vertices; shift-click them one after the other in the order "
            "they should be visited".format(minimum))
    if not complete and len(vertices) > 2:
        raise core.UnwrapFailed(
            "The order of the {} selected vertices is not known: shift-click the vertices "
            "one after the other, in the order they should be visited (box and lasso "
            "selections do not keep the order)".format(len(vertices)))
    return vertices


# ----------------------------------------------------------------------------------------------
# The Curve object


def bounding_box_diagonal(mesh):
    n = len(mesh.vertices)
    if n == 0:
        return 0.0
    co = [0.0] * (3 * n)
    mesh.vertices.foreach_get("co", co)
    lo = [min(co[k::3]) for k in range(3)]
    hi = [max(co[k::3]) for k in range(3)]
    return math.dist(lo, hi)


def _add_splines(curve, polylines, closed):
    """One poly spline per polyline (``array('d')``, 3 per point); a closed path repeats its first
    point, the spline is cyclic instead. Returns how many splines were made."""
    made = 0
    for flat, is_closed in zip(polylines, closed):
        n = len(flat) // 3
        if is_closed and n > 1:
            n -= 1
        if n < 2:
            continue  # a contracted path has nothing to draw
        spline = curve.splines.new('POLY')
        spline.points.add(n - 1)
        coords = []
        for i in range(n):
            coords.extend((flat[3 * i], flat[3 * i + 1], flat[3 * i + 2], 1.0))
        spline.points.foreach_set("co", coords)
        spline.use_cyclic_u = bool(is_closed)
        made += 1
    return made


def make_curve_object(source, polylines, closed, *, name=None):
    """Creates ``<source>_geodesic`` (replacing the one made earlier for ``source``), parented to
    ``source`` with its transform, in the collections of ``source``."""
    intrinsic_draw.remove_previous(source, KIND)
    name = name or source.name + CURVE_SUFFIX
    curve = bpy.data.curves.new(name, 'CURVE')
    curve.dimensions = '3D'
    curve.resolution_u = 1
    curve.bevel_mode = 'ROUND'
    curve.bevel_depth = BEVEL_FRACTION * bounding_box_diagonal(source.data)
    curve.bevel_resolution = 2
    _add_splines(curve, polylines, closed)
    obj = bpy.data.objects.new(name, curve)
    intrinsic_draw._link_like(source, obj)
    obj.parent = source
    obj.matrix_parent_inverse = Matrix.Identity(4)
    obj[intrinsic_draw.SOURCE_KEY] = source.name
    obj[intrinsic_draw.KIND_KEY] = KIND
    return obj


def _certificate_text(certificate):
    if certificate is None:
        return "not computed"
    if "geodesic" in certificate:
        if certificate["geodesic"]:
            return "ok"
        return "FAILED: {} joint(s) not straight, smallest angle {:.4f} rad".format(
            len(certificate["failures"]), certificate["min_angle"])
    if certificate["passed"]:
        return "ok"
    return "FAILED: {} fault(s), {}".format(
        certificate["n_faults"], "; ".join(certificate["faults"][:3]))


def set_properties(target, result_report, lengths, *, mode, representation, mollify, source,
                   stop_ratio, rounds):
    r = result_report
    properties = {
        "source": source.name,
        "mode": mode,
        "representation": representation,
        "mollify": bool(mollify),
        "stop_ratio": float(stop_ratio),
        "length": float(r["length"]),
        "initial_length": float(r["initial_length"]),
        "n_paths": int(r["n_paths"]),
        "status": r["status"],
        "iterations": int(r["iterations"]),
        "flips": int(r["flips"]),
        "refused_flips": int(r["refused_flips"]),
        "stalled": int(r["stalled"]),
        "n_crossings": int(r["n_crossings"]),
        "n_triangulation_vertices": int(r["n_triangulation_vertices"]),
        "certificate": _certificate_text(r["certificate"]),
        "certificate_min_angle_deg": float(math.degrees(r["certificate"]["min_angle"])),
        "input_certificate": _certificate_text(r["input_certificate"]),
    }
    if r["bezier"] is not None:
        b = r["bezier"]
        properties.update({
            "rounds": int(rounds),
            "n_splits": int(b["n_splits"]),
            "n_snapped": int(b["n_snapped"]),
            "n_shortenings": int(b["n_shortenings"]),
            "n_not_geodesic": int(b["n_not_geodesic"]),
            "n_control_points": len(b["control_points"]),
        })
    for key, value in properties.items():
        target[PREFIX + key] = value


# ----------------------------------------------------------------------------------------------
# The work


def _common_options(representation, mollify, stop_ratio, certify_input, show_triangulation):
    if mollify and representation != 'INTEGER':
        raise core.UnwrapFailed("Mollify needs the integer coordinates representation")
    return {
        "representation": representation.lower(),
        "mollify": MOLLIFY if mollify else None,
        "stop_length_ratio": stop_ratio if stop_ratio and stop_ratio > 0.0 else None,
        "certify_input": bool(certify_input),
        "triangulation_edges": bool(show_triangulation),
    }


def _finish(source, out, *, mode, representation, mollify, stop_ratio, rounds,
            show_triangulation, start, native_seconds):
    result = Result()
    result.report = out.report
    result.lengths = list(out.lengths)
    result.seconds_native = native_seconds
    result.curve_object = make_curve_object(source, out.polylines, out.closed)
    set_properties(result.curve_object, out.report, out.lengths, mode=mode,
                   representation=representation, mollify=mollify, source=source,
                   stop_ratio=stop_ratio, rounds=rounds)
    intrinsic_draw.remove_previous(source, KIND_TRIANGULATION)
    if show_triangulation and out.triangulation_edges is not None:
        edges = out.triangulation_edges
        triangulation = intrinsic_draw.draw_intrinsic_edges(
            source, edges["points"], edges["offsets"], name=source.name + TRIANGULATION_SUFFIX,
            radius=0.0, replace=False, kind=KIND_TRIANGULATION)
        result.triangulation_object = triangulation
        triangulation["ddg_n_edges"] = len(edges["edges"])
        triangulation["ddg_representation"] = representation
    result.seconds_blender = perf_counter() - start - native_seconds
    return result


def geodesic_curve(obj, vertices, *, mode='PATH', representation=REPRESENTATION_DEFAULT,
                   mollify=False, stop_ratio=0.0, rounds=4, pin_points=True,
                   show_triangulation=False, certify_input=False):
    """Computes the geodesic through ``vertices`` (indices of the mesh, in order) and creates
    ``<obj>_geodesic``; returns a ``Result``. ``obj`` must not be in Edit mode. ``mode``:
    ``'PATH'``, ``'SMOOTH'``, ``'LOOP'`` or ``'BEZIER'``. ``stop_ratio``: 0 runs to a geodesic,
    otherwise stop at that fraction of the initial length. ``rounds``: Bezier subdivision rounds
    (1 to 8). ``pin_points``: for loops, keep the loop through its points. The mesh of ``obj`` is
    not changed; the objects of an earlier run are replaced. Raises ``core.UnwrapFailed`` with a
    message for the user."""
    start = perf_counter()
    ddg = core.native()
    mesh = obj.data
    if len(mesh.polygons) == 0:
        raise core.UnwrapFailed("Mesh '{}' has no faces".format(mesh.name))
    if mesh.is_editmode:
        raise core.UnwrapFailed("Internal error: the mesh is in Edit mode")
    options = _common_options(representation, mollify, stop_ratio, certify_input,
                              show_triangulation)
    if len(vertices) < 2:
        raise core.UnwrapFailed("Select at least two vertices")
    positions, triangles = triangulated_buffers(mesh)
    ids = array('i', vertices)
    native_start = perf_counter()
    try:
        if mode == 'BEZIER':
            out = ddg.geodesic_bezier(positions, triangles, ids, rounds=rounds, **options)
        elif mode in {'PATH', 'SMOOTH', 'LOOP'}:
            out = ddg.geodesic_path(
                positions, triangles, ids, closed=(mode == 'LOOP'),
                pin_waypoints=(mode == 'PATH') or (mode == 'LOOP' and pin_points), **options)
        else:
            raise core.UnwrapFailed("Unknown mode '{}'".format(mode))
    except ddg.DdgError as ex:
        raise core.UnwrapFailed(str(ex)) from ex
    native_seconds = perf_counter() - native_start
    return _finish(obj, out, mode=mode, representation=representation, mollify=mollify,
                   stop_ratio=stop_ratio, rounds=rounds, show_triangulation=show_triangulation,
                   start=start, native_seconds=native_seconds)


def straighten_seams(obj, pinned=(), *, representation=REPRESENTATION_DEFAULT, mollify=False,
                     stop_ratio=0.0, show_triangulation=False, certify_input=False):
    """Shortens the edges of ``obj`` that are marked as seams, as a network, and creates
    ``<obj>_geodesic``; returns a ``Result``. Vertices where other than two seams meet stay, and
    so do the vertices in ``pinned``. ``stop_ratio``: 0 runs to a geodesic (0.9 is the preset
    for segmentations). The mesh is not changed. Raises ``core.UnwrapFailed``."""
    start = perf_counter()
    ddg = core.native()
    mesh = obj.data
    if len(mesh.polygons) == 0:
        raise core.UnwrapFailed("Mesh '{}' has no faces".format(mesh.name))
    if mesh.is_editmode:
        raise core.UnwrapFailed("Internal error: the mesh is in Edit mode")
    options = _common_options(representation, mollify, stop_ratio, certify_input,
                              show_triangulation)
    seams = core.seam_edges(mesh)
    if len(seams) == 0:
        raise core.UnwrapFailed(
            "Mesh '{}' has no seams: mark edges as seams first (Edge > Mark Seam)".format(
                mesh.name))
    positions, triangles = triangulated_buffers(mesh)
    native_start = perf_counter()
    try:
        out = ddg.geodesic_network(
            positions, triangles, seams, pinned=array('i', pinned), **options)
    except ddg.DdgError as ex:
        raise core.UnwrapFailed(str(ex)) from ex
    native_seconds = perf_counter() - native_start
    return _finish(obj, out, mode='SEAMS', representation=representation, mollify=mollify,
                   stop_ratio=stop_ratio, rounds=0, show_triangulation=show_triangulation,
                   start=start, native_seconds=native_seconds)


def summary(result, name, what):
    r = result.report
    parts = ["{} path(s), length {:.5g}".format(r["n_paths"], r["length"])]
    if r["initial_length"] > 0.0:
        parts.append("{:.1f} % of the initial path".format(100.0 * r["length"] / r["initial_length"]))
    parts.append("{} flips".format(r["flips"]))
    parts.append("status {}".format(r["status"]))
    parts.append("{:.3f} s".format(result.seconds_native + result.seconds_blender))
    return "{} '{}': ".format(what, name) + " | ".join(parts)


# ----------------------------------------------------------------------------------------------
# Operators


def _warnings(operator, result):
    r = result.report
    for key, text in (("certificate", "the geodesic certificate"),
                      ("input_certificate", "the certificate on the input mesh")):
        value = result.curve_object.get(PREFIX + key, "")
        if value.startswith("FAILED"):
            operator.report({'WARNING'}, "{} does not hold: {}".format(text, value))
    if r["status"] not in {"geodesic", "length_ratio"}:
        operator.report({'WARNING'}, "The shortening stopped early ({}); the curve may not "
                        "be geodesic".format(r["status"]))


class _GeodesicProperties:
    representation: EnumProperty(
        name="Representation", items=REPRESENTATION_ITEMS, default=REPRESENTATION_DEFAULT)
    mollify: BoolProperty(
        name="Mollify", description="Intrinsic mollification of the input edge lengths first "
        "(relative margin 1e-4); integer coordinates only", default=False)
    stop_ratio: FloatProperty(
        name="Stop at Length Ratio", description="Stop once the length is below this fraction "
        "of the initial length (for example 0.9); 0 shortens to a geodesic",
        default=0.0, min=0.0, max=1.0)
    show_triangulation: BoolProperty(
        name="Show Triangulation", description="Also draw the edges of the final intrinsic "
        "triangulation on the input mesh, as a second curve object", default=False)
    certify_input: BoolProperty(
        name="Verify on Input Mesh", description="Also check the result on the input mesh "
        "alone (slower; the outcome is a custom property of the curve)", default=False)

    def _draw_common(self, layout):
        layout.prop(self, "representation")
        row = layout.row()
        row.enabled = self.representation == 'INTEGER'
        row.prop(self, "mollify")
        layout.prop(self, "stop_ratio")
        layout.prop(self, "show_triangulation")
        layout.prop(self, "certify_input")


def _run(operator, context, work):
    obj = context.active_object
    if obj is None or obj.type != 'MESH':
        operator.report({'ERROR'}, "The active object is not a mesh")
        return {'CANCELLED'}
    return work(obj)


class MESH_OT_ddg_geodesic_path(_GeodesicProperties, bpy.types.Operator):
    bl_idname = "mesh.ddg_geodesic_path"
    bl_label = "Geodesic Path"
    bl_description = (
        "Geodesic path, loop or Bezier curve through the selected vertices, in the order they "
        "were selected, as a new curve object on the surface. The mesh is not changed; running "
        "again replaces the curve. Sharp, Crane: You Can Find Geodesic Paths in Triangle Meshes "
        "by Just Flipping Edges (2020)"
    )
    bl_options = {'REGISTER', 'UNDO'}

    mode: EnumProperty(name="Mode", items=MODE_ITEMS, default='PATH')
    pin_points: BoolProperty(
        name="Pin Points", description="Loop: keep the loop through its points; otherwise they "
        "only seed it and it may slide off them", default=True)
    rounds: IntProperty(
        name="Rounds", description="Bezier: rounds of midpoint subdivision (the curve has "
        "2^rounds pieces)", default=4, min=1, max=8)

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode == 'EDIT'

    def draw(self, _context):
        layout = self.layout
        layout.prop(self, "mode")
        if self.mode == 'LOOP':
            layout.prop(self, "pin_points")
        if self.mode == 'BEZIER':
            layout.prop(self, "rounds")
        self._draw_common(layout)

    def execute(self, context):
        def work(obj):
            try:
                vertices = ordered_selection(obj)
            except core.UnwrapFailed as ex:
                self.report({'ERROR'}, str(ex))
                return {'CANCELLED'}
            bpy.ops.object.mode_set(mode='OBJECT')
            try:
                result = geodesic_curve(
                    obj, vertices, mode=self.mode, representation=self.representation,
                    mollify=self.mollify, stop_ratio=self.stop_ratio, rounds=self.rounds,
                    pin_points=self.pin_points, show_triangulation=self.show_triangulation,
                    certify_input=self.certify_input)
            except core.UnwrapFailed as ex:
                self.report({'ERROR'}, str(ex))
                return {'CANCELLED'}
            finally:
                bpy.ops.object.mode_set(mode='EDIT')
            _warnings(self, result)
            self.report({'INFO'}, summary(result, obj.name, "Geodesic " + self.mode.title()))
            return {'FINISHED'}

        return _run(self, context, work)


class MESH_OT_ddg_straighten_seams(_GeodesicProperties, bpy.types.Operator):
    bl_idname = "mesh.ddg_straighten_seams"
    bl_label = "Straighten Seams"
    bl_description = (
        "Shorten the edges marked as seams to geodesics on the surface, as a network that keeps "
        "its junctions (and the selected vertices), and show the result as a new curve object. "
        "The mesh is not changed; running again replaces the curve. Sharp, Crane: You Can Find "
        "Geodesic Paths in Triangle Meshes by Just Flipping Edges (2020)"
    )
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}

    def draw(self, _context):
        self._draw_common(self.layout)

    def execute(self, context):
        def work(obj):
            was_edit_mode = obj.mode == 'EDIT'
            pinned = []
            if was_edit_mode:
                bm = bmesh.from_edit_mesh(obj.data)
                pinned = [v.index for v in bm.verts if v.select]
                bpy.ops.object.mode_set(mode='OBJECT')
            try:
                result = straighten_seams(
                    obj, pinned, representation=self.representation, mollify=self.mollify,
                    stop_ratio=self.stop_ratio, show_triangulation=self.show_triangulation,
                    certify_input=self.certify_input)
            except core.UnwrapFailed as ex:
                self.report({'ERROR'}, str(ex))
                return {'CANCELLED'}
            finally:
                if was_edit_mode:
                    bpy.ops.object.mode_set(mode='EDIT')
            _warnings(self, result)
            self.report({'INFO'}, summary(result, obj.name, "Straightened seams of"))
            return {'FINISHED'}

        return _run(self, context, work)


class VIEW3D_PT_ddg_geodesics(bpy.types.Panel):
    bl_label = "Geodesics"
    bl_idname = "VIEW3D_PT_ddg_geodesics"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "DDG"

    def draw(self, context):
        layout = self.layout
        obj = context.active_object
        in_edit = obj is not None and obj.type == 'MESH' and obj.mode == 'EDIT'
        column = layout.column(align=True)
        column.enabled = in_edit
        column.label(text="Selected vertices, in selection order:")
        for mode, label, icon in (('PATH', "Path", 'CURVE_PATH'), ('SMOOTH', "Smooth", 'SMOOTHCURVE'),
                                  ('LOOP', "Loop", 'MESH_CIRCLE'), ('BEZIER', "Bezier", 'CURVE_BEZCURVE')):
            column.operator(MESH_OT_ddg_geodesic_path.bl_idname, text=label, icon=icon).mode = mode
        column = layout.column(align=True)
        column.label(text="Seams:")
        column.operator(MESH_OT_ddg_straighten_seams.bl_idname, text="Straighten Seams")
        column.operator(MESH_OT_ddg_straighten_seams.bl_idname,
                        text="Straighten Seams (stop at 0.9)").stop_ratio = 0.9


CLASSES = (MESH_OT_ddg_geodesic_path, MESH_OT_ddg_straighten_seams, VIEW3D_PT_ddg_geodesics)


def _draw_edit_entries(self, _context):
    self.layout.separator()
    self.layout.operator(MESH_OT_ddg_geodesic_path.bl_idname)
    self.layout.operator(MESH_OT_ddg_straighten_seams.bl_idname)


def _draw_object_entry(self, _context):
    self.layout.operator(MESH_OT_ddg_straighten_seams.bl_idname)


def _draw_uv_entry(self, _context):
    self.layout.operator(MESH_OT_ddg_straighten_seams.bl_idname)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.VIEW3D_MT_edit_mesh.append(_draw_edit_entries)
    bpy.types.VIEW3D_MT_object.append(_draw_object_entry)
    bpy.types.VIEW3D_MT_uv_map.append(_draw_uv_entry)


def unregister():
    bpy.types.VIEW3D_MT_uv_map.remove(_draw_uv_entry)
    bpy.types.VIEW3D_MT_object.remove(_draw_object_entry)
    bpy.types.VIEW3D_MT_edit_mesh.remove(_draw_edit_entries)
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
