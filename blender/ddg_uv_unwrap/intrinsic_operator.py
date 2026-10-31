# SPDX-License-Identifier: MIT
"""Show Intrinsic Refinement: the intrinsic Delaunay refinement of a mesh and its common
subdivision with the input, as new objects (Fig. 12 of Gillespie, Sharp, Crane, Integer
Coordinates for Intrinsic Geometry Processing, 2021).

The mesh is triangulated (as Blender's tessellation does), turned into an intrinsic
triangulation on integer coordinates (or on signposts, Sharp, Soliman, Crane 2019), flipped to
Delaunay and refined to an angle bound. The result is shown as

* a mesh object ``<mesh>_intrinsic`` whose faces are the pieces of the common subdivision (every
  piece lies in one face of the input and in one face of the intrinsic triangulation), with
  face attributes ``intrinsic_face`` and ``input_face`` and a color attribute
  ``intrinsic_face_color`` (a corner color attribute, one color per intrinsic face), and
* optionally a curve object ``<mesh>_intrinsic_edges`` with the intrinsic edges drawn on the
  input mesh (``intrinsic_draw.draw_intrinsic_edges``).

The input mesh is never changed. Running the operator again replaces both objects. Counts, the
smallest angle and the certificate are custom properties of the mesh object (``ddg_*``).
"""

from time import perf_counter

import bpy
from bpy.props import BoolProperty, EnumProperty, FloatProperty

from . import core, intrinsic_draw

REPRESENTATION_ITEMS = (
    ('INTEGER', "Integer Coordinates",
     "Gillespie, Sharp, Crane 2021: normal coordinates and roundabouts, exact combinatorics"),
    ('SIGNPOSTS', "Signposts", "Sharp, Soliman, Crane 2019: edge lengths and signpost angles"),
)
RULES_ITEMS = (
    ('CHEW', "Chew (paper)",
     "The integer paper's rules: split a boundary edge at its midpoint when a circumcenter path "
     "leaves through it, remove nearby inserted vertices"),
    ('RUPPERT', "Ruppert",
     "The signpost paper's rules: split encroached boundary edges, insert circumcenters"),
)

# Name of the property that carries the result's numbers.
PREFIX = "ddg_"
MESH_SUFFIX = "_intrinsic"


def _mollify_value(use):
    # The paper's relative margin (Gillespie, Sharp, Crane 2021, Sec. 7).
    return 1e-4 if use else None


def triangulated_buffers(mesh):
    """``(positions, triangles)`` of ``mesh`` as arrays for ``ddg_tools``: Blender's own
    tessellation, vertex indices per triangle."""
    from array import array

    positions, corner_vertices, _starts, corner_triangles = core.mesh_buffers(mesh)
    triangles = array('I', (corner_vertices[c] for c in corner_triangles))
    return positions, triangles


class Result:
    """What ``show_intrinsic_refinement`` made."""

    def __init__(self):
        self.mesh_object = None
        self.edges_object = None
        self.refine = {}
        self.counts = {}
        self.certificate = None
        self.certificate_error = None
        self.seconds_native = 0.0
        self.seconds_blender = 0.0


def _build_mesh_object(source, sub, n_intrinsic_faces):
    name = source.name + MESH_SUFFIX
    positions = sub["positions"]
    n = len(positions) // 3
    vertices = [(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]) for i in range(n)]
    starts, corners = sub["polygon_starts"], sub["corners"]
    faces = [tuple(corners[starts[i]:starts[i + 1]]) for i in range(len(starts) - 1)]
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    n_faces = len(faces)
    if len(mesh.polygons) != n_faces:
        raise core.UnwrapFailed("Internal error: Blender dropped faces of the common subdivision")

    intrinsic = list(sub["intrinsic_face"])
    original = list(sub["input_face"])
    attribute = mesh.attributes.new("intrinsic_face", 'INT', 'FACE')
    attribute.data.foreach_set("value", intrinsic)
    attribute = mesh.attributes.new("input_face", 'INT', 'FACE')
    attribute.data.foreach_set("value", original)
    palette = intrinsic_draw.face_colors(n_intrinsic_faces)
    # Color attributes live on points or corners: one color per corner of the piece.
    color = mesh.color_attributes.new("intrinsic_face_color", 'FLOAT_COLOR', 'CORNER')
    flat = []
    for i in range(n_faces):
        flat.extend(palette[intrinsic[i]] * (starts[i + 1] - starts[i]))
    color.data.foreach_set("color", flat)
    # Is the vertex a vertex of the intrinsic triangulation (1) or a crossing (0)?
    point = mesh.attributes.new("intrinsic_vertex", 'INT', 'POINT')
    n_intrinsic = sub["n_intrinsic_vertices"]
    point.data.foreach_set("value", [1 if i < n_intrinsic else 0 for i in range(n)])

    obj = bpy.data.objects.new(name, mesh)
    obj.matrix_world = source.matrix_world.copy()
    obj[intrinsic_draw.SOURCE_KEY] = source.name
    obj[intrinsic_draw.KIND_KEY] = "mesh"
    intrinsic_draw._link_like(source, obj)
    return obj


def show_intrinsic_refinement(
        obj, *, min_angle=25.0, rules='CHEW', representation='INTEGER', mollify=False,
        draw_edges=True, edge_radius=0.0, max_insertions=0):
    """Refines the triangulated mesh of ``obj`` and creates ``<name>_intrinsic`` (the common
    subdivision) and, with ``draw_edges``, ``<name>_intrinsic_edges``. Returns a ``Result``.
    The mesh of ``obj`` is not changed; the objects of an earlier run are replaced.

    ``min_angle``: angle bound in degrees (0 to 35). ``rules``: ``'CHEW'`` or ``'RUPPERT'``.
    ``representation``: ``'INTEGER'`` or ``'SIGNPOSTS'``. ``mollify``: intrinsic mollification
    (integer coordinates only). ``max_insertions``: 0 for the default limit (10 per vertex,
    at least 1000). Raises ``core.UnwrapFailed`` with a message for the user.
    """
    start = perf_counter()
    ddg = core.native()
    mesh = obj.data
    if len(mesh.polygons) == 0:
        raise core.UnwrapFailed("Mesh '{}' has no faces".format(mesh.name))
    if mesh.is_editmode:
        raise core.UnwrapFailed("Internal error: the mesh is in Edit mode")
    if mollify and representation != 'INTEGER':
        raise core.UnwrapFailed("Mollify needs the integer coordinates representation")
    positions, triangles = triangulated_buffers(mesh)
    result = Result()
    native_start = perf_counter()
    try:
        tri = ddg.IntrinsicTriangulation(
            positions, triangles, representation=representation.lower(),
            mollify=_mollify_value(mollify))
        tri.flip_to_delaunay()
        result.refine = tri.refine(
            min_angle_degrees=min_angle, rules=rules.lower(),
            max_insertions=max_insertions or None)
        result.counts = tri.counts()
        sub = tri.common_subdivision()
        polylines = tri.edge_polylines() if draw_edges else None
        if representation == 'INTEGER':
            try:
                result.certificate = tri.certify()
            except ddg.SolveError as ex:
                result.certificate_error = str(ex)
    except ddg.DdgError as ex:
        raise core.UnwrapFailed(str(ex)) from ex
    result.seconds_native = perf_counter() - native_start

    intrinsic_draw.remove_previous(obj, "mesh")
    intrinsic_draw.remove_previous(obj, "edges")
    result.mesh_object = _build_mesh_object(obj, sub, result.counts["n_faces"])
    _set_properties(result.mesh_object, obj, result, min_angle, rules, representation, mollify)
    if draw_edges:
        result.edges_object = intrinsic_draw.draw_intrinsic_edges(
            obj, polylines["points"], polylines["offsets"], radius=edge_radius, replace=False)
    result.seconds_blender = perf_counter() - start - result.seconds_native
    return result


def _set_properties(target, source, result, min_angle, rules, representation, mollify):
    c, r = result.counts, result.refine
    properties = {
        "source": source.name,
        "representation": representation,
        "rules": rules,
        "mollify": bool(mollify),
        "angle_bound_deg": float(min_angle),
        "min_angle_deg": float(c["min_angle"] * 180.0 / 3.141592653589793),
        "min_angle_valid_deg": float(r["min_angle_valid"] * 180.0 / 3.141592653589793),
        "n_below_bound": int(r["n_below_bound"]),
        "n_input_vertices": int(c["n_input_vertices"]),
        "n_input_faces": int(c["n_input_faces"]),
        "n_vertices": int(c["n_vertices"]),
        "n_edges": int(c["n_edges"]),
        "n_faces": int(c["n_faces"]),
        "n_inserted": int(r["n_inserted"]),
        "n_flips": int(r["n_flips"]),
        "n_crossings": int(c["n_crossings"]),
        "refine_status": r["status"],
    }
    if result.certificate is not None:
        properties["certificate"] = "ok"
        for key, value in result.certificate.items():
            properties["certificate_" + key] = int(value)
    elif result.certificate_error is not None:
        properties["certificate"] = "FAILED: " + result.certificate_error
    else:
        properties["certificate"] = "none (signposts)"
    for key, value in properties.items():
        target[PREFIX + key] = value


def summary(result, name):
    c, r = result.counts, result.refine
    parts = [
        "{} intrinsic vertices ({} inserted), {} faces".format(
            c["n_vertices"], r["n_inserted"], c["n_faces"]),
        "smallest angle {:.2f} deg".format(c["min_angle"] * 180.0 / 3.141592653589793),
    ]
    if c["n_crossings"]:
        parts.append("{} crossings".format(c["n_crossings"]))
    if result.certificate is not None:
        parts.append("certified")
    parts.append("{:.3f} s".format(result.seconds_native + result.seconds_blender))
    return "Intrinsic refinement of '{}': ".format(name) + " | ".join(parts)


class MESH_OT_ddg_show_intrinsic_refinement(bpy.types.Operator):
    bl_idname = "mesh.ddg_show_intrinsic_refinement"
    bl_label = "Show Intrinsic Refinement"
    bl_description = (
        "Refine an intrinsic Delaunay triangulation of the active mesh to an angle bound and show "
        "the common subdivision with the input as a new mesh object, colored by intrinsic face, "
        "and the intrinsic edges drawn on the input as a curve. The input mesh is not changed; "
        "running again replaces the new objects. Gillespie, Sharp, Crane: Integer Coordinates "
        "for Intrinsic Geometry Processing (2021)"
    )
    bl_options = {'REGISTER', 'UNDO'}

    min_angle: FloatProperty(
        name="Angle Bound", description="Refine until every corner is at least this angle "
        "(degrees); the provable range is up to 30", default=25.0, min=0.0, max=35.0)
    rules: EnumProperty(name="Rules", items=RULES_ITEMS, default='CHEW')
    representation: EnumProperty(
        name="Representation", items=REPRESENTATION_ITEMS, default='INTEGER')
    mollify: BoolProperty(
        name="Mollify", description="Intrinsic mollification of the input edge lengths first "
        "(relative margin 1e-4, as the paper); integer coordinates only", default=False)
    draw_edges: BoolProperty(
        name="Draw Intrinsic Edges", description="Also create a curve object with the "
        "intrinsic edges drawn on the input mesh", default=True)
    edge_radius: FloatProperty(
        name="Edge Radius", description="Bevel radius of the intrinsic edge curve; 0 draws "
        "thin lines", default=0.0, min=0.0, soft_max=0.05, unit='LENGTH')

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}

    def draw(self, _context):
        layout = self.layout
        layout.prop(self, "min_angle")
        layout.prop(self, "rules")
        layout.prop(self, "representation")
        row = layout.row()
        row.enabled = self.representation == 'INTEGER'
        row.prop(self, "mollify")
        layout.prop(self, "draw_edges")
        row = layout.row()
        row.enabled = self.draw_edges
        row.prop(self, "edge_radius")

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        was_edit_mode = obj.mode == 'EDIT'
        if was_edit_mode:
            bpy.ops.object.mode_set(mode='OBJECT')
        try:
            result = show_intrinsic_refinement(
                obj, min_angle=self.min_angle, rules=self.rules,
                representation=self.representation, mollify=self.mollify,
                draw_edges=self.draw_edges, edge_radius=self.edge_radius)
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        finally:
            if was_edit_mode:
                bpy.ops.object.mode_set(mode='EDIT')
        if result.certificate_error is not None:
            self.report({'WARNING'}, "The integer coordinates do not certify: {}".format(
                result.certificate_error))
        if result.refine["n_below_bound"]:
            self.report({'INFO'}, "{} corners stay below the bound (narrow cones, boundary "
                        "angles or the insertion limit)".format(result.refine["n_below_bound"]))
        self.report({'INFO'}, summary(result, obj.name))
        return {'FINISHED'}
