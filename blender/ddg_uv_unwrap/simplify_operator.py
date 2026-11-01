# SPDX-License-Identifier: MIT
"""Intrinsic Simplify: a coarse-mesh view of the active mesh, after Liu, Gillespie, Chislett,
Sharp, Jacobson, Crane, Surface Simplification using Intrinsic Error Metrics (ACM TOG 42(4),
2023; Fig. 1 and Sec. 8.1).

The mesh is simplified by intrinsic vertex removal with the intrinsic error metric (library
``ddg_tools.simplify``), but its geometry is never changed: the coarse mesh is a coarse intrinsic
triangulation of the same surface. It is shown as

* a corner color attribute ``IntrinsicSimplifyFaces`` that colors every face of the mesh by the
  coarse face its centroid lands in (the coarse faces are colored greedily so that faces that
  share a vertex differ, from a fixed palette that is cycled; deterministic),
* a face integer attribute ``IntrinsicSimplifyFaceId`` with that coarse face's id, and
* a vertex group ``IntrinsicSimplifyVertices`` with the vertices that stay (weight 1).

One tracked point per face of Blender's tessellation (the centroid, in barycentric coordinates of
that triangle) is carried through the vertex removals; a corner takes the color of the triangle
of the tessellation that contains it (n-gons are split the way Blender splits them, so an n-gon
may show several colors). Counts and the report are custom properties of the mesh (``ddg_simplify_*``)
and a line in the info area. As with Show Developability the viewport color is not changed:
set it to Attribute to see the colors.
"""

from array import array
from time import perf_counter

import bpy
from bpy.props import BoolProperty, EnumProperty, FloatProperty, IntProperty

from . import core
from .developable_operator import selected_vertices

COLOR_ATTRIBUTE = "IntrinsicSimplifyFaces"
FACE_ID_ATTRIBUTE = "IntrinsicSimplifyFaceId"
VERTEX_GROUP = "IntrinsicSimplifyVertices"
PREFIX = "ddg_simplify_"

# A fixed qualitative palette (Tableau 20 without the grays and pale tints), cycled by color index.
PALETTE = (
    (0.122, 0.467, 0.706), (1.000, 0.498, 0.055), (0.173, 0.627, 0.173), (0.839, 0.153, 0.157),
    (0.580, 0.404, 0.741), (0.549, 0.337, 0.294), (0.890, 0.467, 0.761), (0.737, 0.741, 0.133),
    (0.090, 0.745, 0.812), (0.682, 0.780, 0.910), (1.000, 0.733, 0.471), (0.596, 0.875, 0.541),
)

TARGET_ITEMS = (
    ('PERCENT', "Percent", "Keep this percentage of the vertices"),
    ('COUNT', "Vertex Count", "Keep this many vertices"),
)
REPRESENTATION_ITEMS = (
    ('INTEGER', "Integer Coordinates", "Gillespie, Sharp, Crane 2021 (for the refinements)"),
    ('SIGNPOSTS', "Signposts", "Sharp, Soliman, Crane 2019 (for the refinements)"),
)

POINT_MAP_ITEMS = (
    ('NEGATIVE', "Negative (e^-u)",
     "The circumcircle-preserving map of Springborn et al. 2008 that the paper cites; much less "
     "distortion"),
    ('POSITIVE', "Positive (e^+u)", "Eq. 10 as printed in the paper, as the authors' code does"),
)


class Result:
    """What ``intrinsic_simplify`` made."""

    def __init__(self):
        self.report = {}
        self.n_input_vertices = 0
        self.n_vertices = 0
        self.n_faces = 0
        self.n_colours = 0
        self.n_group_vertices = 0
        self.n_fixed = 0
        self.coarse_face = None  # per face of the tessellation
        self.colours = None  # per face of the tessellation
        self.seconds_native = 0.0
        self.seconds = 0.0


def _check_attribute(mesh, name, domain, data_type):
    existing = mesh.attributes.get(name)
    if existing is not None and (existing.domain != domain or existing.data_type != data_type):
        raise core.UnwrapFailed(
            "Mesh '{}' already has an attribute '{}' that is not a {} {} attribute".format(
                mesh.name, name, domain.lower(), data_type.lower()))


def intrinsic_simplify(
        obj, *, target_mode='PERCENT', percent=10.0, count=100, area_weight=0.0, memory=True,
        fix_selected=False, refine_before=False, refine_after=False, representation='INTEGER',
        point_map='NEGATIVE'):
    """Simplifies the mesh of ``obj`` (library ``ddg_tools.simplify``) without changing it and
    writes the corner color attribute, the face id attribute and the vertex group. Returns a
    ``Result``. Raises ``core.UnwrapFailed`` with a message for the user."""
    start = perf_counter()
    ddg = core.native()
    mesh = obj.data
    if mesh.library is not None:
        raise core.UnwrapFailed(
            "Mesh '{}' is linked from another file and cannot be changed here".format(mesh.name))
    if len(mesh.polygons) == 0:
        raise core.UnwrapFailed("Mesh '{}' has no faces".format(mesh.name))
    if mesh.is_editmode:
        raise core.UnwrapFailed("Internal error: the mesh is in Edit mode")
    _check_attribute(mesh, COLOR_ATTRIBUTE, 'CORNER', 'FLOAT_COLOR')
    _check_attribute(mesh, FACE_ID_ATTRIBUTE, 'FACE', 'INT')
    positions, corner_vertices, _starts, corner_triangles = core.mesh_buffers(mesh)
    n_triangles = len(corner_triangles) // 3
    triangles = array('I', (corner_vertices[c] for c in corner_triangles))
    fixed = selected_vertices(mesh) if fix_selected else array('i')
    points = (array('i', range(n_triangles)), array('d', [1.0 / 3.0] * (3 * n_triangles)))
    options = dict(weights=float(area_weight), memory=memory, points=points,
                   refine_before=refine_before, refine_after=refine_after,
                   representation=representation.lower(), projective_sign=point_map.lower())
    if len(fixed):
        options["fixed"] = fixed
    if target_mode == 'PERCENT':
        options["ratio"] = percent / 100.0
    else:
        options["target"] = int(count)
    native_start = perf_counter()
    try:
        s = ddg.simplify(positions, triangles, **options)
    except ddg.DdgError as ex:
        raise core.UnwrapFailed(str(ex)) from ex
    result = Result()
    result.seconds_native = perf_counter() - native_start
    result.report = s.report
    result.n_input_vertices = s.n_input_vertices
    result.n_vertices = s.n_vertices
    result.n_faces = s.n_faces
    result.n_colours = s.n_colours
    result.n_fixed = len(fixed)
    n_vertices = len(mesh.vertices)
    coarse_face = list(s.point_face)[n_vertices:]
    point_colours = list(s.point_colours)[n_vertices:]
    if len(coarse_face) != n_triangles or min(coarse_face, default=0) < 0:
        raise core.UnwrapFailed("Internal error: a face centroid was not located on the coarse mesh")
    result.coarse_face = coarse_face
    result.colours = point_colours

    # Color of every corner: the color of the first triangle of the tessellation that has it.
    n_loops = len(mesh.loops)
    flat = array('f', bytes(4 * 4 * n_loops))
    done = bytearray(n_loops)
    polygon_of = array('i', bytes(4 * n_triangles))
    mesh.loop_triangles.foreach_get("polygon_index", polygon_of)
    face_ids = array('i', bytes(4 * len(mesh.polygons)))
    seen_polygon = bytearray(len(mesh.polygons))
    for t in range(n_triangles):
        r, g, b = PALETTE[point_colours[t] % len(PALETTE)]
        for loop in corner_triangles[3 * t:3 * t + 3]:
            if not done[loop]:
                done[loop] = 1
                flat[4 * loop:4 * loop + 4] = array('f', (r, g, b, 1.0))
        p = polygon_of[t]
        if not seen_polygon[p]:
            seen_polygon[p] = 1
            face_ids[p] = coarse_face[t]
    color = mesh.color_attributes.get(COLOR_ATTRIBUTE)
    if color is None:
        color = mesh.color_attributes.new(COLOR_ATTRIBUTE, 'FLOAT_COLOR', 'CORNER')
    color.data.foreach_set("color", flat)
    face_id = mesh.attributes.get(FACE_ID_ATTRIBUTE)
    if face_id is None:
        face_id = mesh.attributes.new(FACE_ID_ATTRIBUTE, 'INT', 'FACE')
    face_id.data.foreach_set("value", face_ids)
    # Creating the integer attribute invalidates the Python reference `color`, so the color
    # attribute is made active by name (as Show Developability does).
    mesh.color_attributes.active_color_name = COLOR_ATTRIBUTE
    mesh.color_attributes.render_color_index = mesh.color_attributes.active_color_index

    group = obj.vertex_groups.get(VERTEX_GROUP)
    if group is not None:
        obj.vertex_groups.remove(group)
    group = obj.vertex_groups.new(name=VERTEX_GROUP)
    kept = [v for v in s.vertices if v < n_vertices]
    result.n_group_vertices = len(kept)
    if kept:
        group.add(kept, 1.0, 'REPLACE')
    result.seconds = perf_counter() - start

    properties = {
        "n_input_vertices": n_vertices,
        "n_vertices": s.n_vertices,
        "n_group_vertices": len(kept),
        "n_faces": s.n_faces,
        "n_colours": s.n_colours,
        "n_fixed": len(fixed),
        "stop": result.report["stop"],
        "memory": bool(memory),
        "area_weight": float(area_weight),
        "refine_before": bool(refine_before),
        "refine_after": bool(refine_after),
        "seconds": result.seconds,
    }
    for key, value in properties.items():
        mesh[PREFIX + key] = value
    mesh.update()
    return result


def summary(result, name):
    r = result.report
    parts = [
        "{} -> {} vertices ({} faces)".format(
            result.n_input_vertices, result.n_vertices, result.n_faces),
        "{} interior, {} boundary, {} ear removals".format(
            r["n_interior"], r["n_boundary"], r["n_ear"]),
        "{} colors".format(result.n_colours),
    ]
    if result.n_fixed:
        parts.append("{} fixed".format(result.n_fixed))
    if result.n_vertices > result.n_group_vertices:
        parts.append("{} inserted vertices not in the group".format(
            result.n_vertices - result.n_group_vertices))
    parts.append("{:.3f} s".format(result.seconds))
    return "Intrinsic simplify of '{}': ".format(name) + " | ".join(parts)


class OBJECT_OT_ddg_intrinsic_simplify(bpy.types.Operator):
    bl_idname = "object.ddg_intrinsic_simplify"
    bl_label = "Intrinsic Simplify"
    bl_description = (
        "Simplify the active mesh by intrinsic vertex removal with the intrinsic error metric and "
        "show the coarse mesh on it: a corner color attribute colors every face by the coarse face "
        "its centroid lands in, and a vertex group holds the coarse vertices. The mesh geometry is "
        "not changed (set the viewport color to Attribute to see the colors). Liu, Gillespie, "
        "Chislett, Sharp, Jacobson, Crane: Surface Simplification using Intrinsic Error Metrics "
        "(2023)"
    )
    bl_options = {'REGISTER', 'UNDO'}

    target_mode: EnumProperty(name="Target", items=TARGET_ITEMS, default='PERCENT')
    percent: FloatProperty(
        name="Percent", description="Percentage of the vertices to keep", default=10.0,
        min=0.01, max=100.0, precision=2, subtype='PERCENTAGE')
    count: IntProperty(name="Vertices", description="Number of vertices to keep", default=100,
                       min=4)
    area_weight: FloatProperty(
        name="Area Weight", description="Curvature versus area in the error metric: 0 keeps "
        "curvature (the paper's metric), 1 keeps area (uniform coarse triangles)",
        default=0.0, min=0.0, max=1.0)
    memory: BoolProperty(
        name="Memory", description="Memory-based error metric (the paper's default); off: the "
        "memoryless metric", default=True)
    fix_selected: BoolProperty(
        name="Fix Selected Vertices", description="Keep the selected vertices (selection as "
        "stored: those selected in Edit mode)", default=False)
    refine_before: BoolProperty(
        name="Refine Before", description="Delaunay-refine the intrinsic triangulation before "
        "simplifying (Fig. 12); inserted vertices are removable and are not in the vertex group",
        default=False)
    refine_after: BoolProperty(
        name="Refine After", description="Delaunay-refine the coarse mesh at the end (Fig. 13); "
        "inserted vertices are not in the vertex group", default=False)
    representation: EnumProperty(
        name="Representation", items=REPRESENTATION_ITEMS, default='INTEGER',
        description="The intrinsic triangulation of the refinements")
    point_map: EnumProperty(
        name="Point Map", items=POINT_MAP_ITEMS, default='NEGATIVE',
        description="How the fine faces move through each flattening (paper Eq. 10)")

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}

    def draw(self, _context):
        layout = self.layout
        layout.prop(self, "target_mode", expand=True)
        layout.prop(self, "percent" if self.target_mode == 'PERCENT' else "count")
        layout.prop(self, "area_weight")
        layout.prop(self, "memory")
        layout.prop(self, "fix_selected")
        layout.prop(self, "refine_before")
        layout.prop(self, "refine_after")
        row = layout.row()
        row.enabled = self.refine_before or self.refine_after
        row.prop(self, "representation")
        layout.prop(self, "point_map")

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        was_edit_mode = obj.mode == 'EDIT'
        if was_edit_mode:
            bpy.ops.object.mode_set(mode='OBJECT')
        try:
            result = intrinsic_simplify(
                obj, target_mode=self.target_mode, percent=self.percent, count=self.count,
                area_weight=self.area_weight, memory=self.memory,
                fix_selected=self.fix_selected, refine_before=self.refine_before,
                refine_after=self.refine_after, representation=self.representation,
                point_map=self.point_map)
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        finally:
            if was_edit_mode:
                bpy.ops.object.mode_set(mode='EDIT')
        if self.fix_selected and result.n_fixed == 0:
            self.report({'WARNING'}, "Fix Selected Vertices: no vertex is selected")
        if result.report["stop"] != "targetreached":
            self.report({'WARNING'}, "Stopped at {} vertices: {}".format(
                result.n_vertices, result.report["stop"]))
        self.report({'INFO'}, summary(result, obj.name)
                    + "; color attribute '{}' (set the viewport color to Attribute)".format(
                        COLOR_ATTRIBUTE))
        return {'FINISHED'}


class VIEW3D_PT_ddg_intrinsic_simplify(bpy.types.Panel):
    bl_label = "Intrinsic Simplify"
    bl_idname = "VIEW3D_PT_ddg_intrinsic_simplify"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "DDG"

    def draw(self, context):
        obj = context.active_object
        column = self.layout.column(align=True)
        column.enabled = obj is not None and obj.type == 'MESH'
        column.label(text="Colors faces by their coarse face:")
        column.operator(OBJECT_OT_ddg_intrinsic_simplify.bl_idname, text="Simplify (10 %)",
                        icon='MOD_DECIM')


CLASSES = (OBJECT_OT_ddg_intrinsic_simplify, VIEW3D_PT_ddg_intrinsic_simplify)


def _draw_entry(self, _context):
    self.layout.operator(OBJECT_OT_ddg_intrinsic_simplify.bl_idname)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.VIEW3D_MT_object.append(_draw_entry)
    bpy.types.VIEW3D_MT_edit_mesh.append(_draw_entry)


def unregister():
    bpy.types.VIEW3D_MT_edit_mesh.remove(_draw_entry)
    bpy.types.VIEW3D_MT_object.remove(_draw_entry)
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
