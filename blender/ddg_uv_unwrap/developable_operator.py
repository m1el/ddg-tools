# SPDX-License-Identifier: MIT
"""Developability operators: Make Developable, Show Developability, Developable Seams + UV.

The work is done by ``ddg_tools.make_developable``, ``ddg_tools.developability`` and
``ddg_tools.developable_seams`` (Rust, crates/ddg-developability and crates/ddg-py), after

    Oded Stein, Eitan Grinspun, Keenan Crane, Developability of Triangle Meshes,
    ACM Transactions on Graphics 37(4), 2018,

implemented clean-room from the text of the paper. The UV map of the third operator is made by
the same pipeline as Conformal Unwrap (``core.unwrap_object``), on a mesh whose seams were chosen
by the paper's cut.

The mesh is read through the same buffers as the other operators (``core.mesh_buffers``): the
triangulation is Blender's own (``mesh.loop_triangles``).

* Topology unchanged (no subdivision, no remeshing operation applied, "Triangulate" off): the new
  positions are written into the existing mesh; polygons, UVs and other data stay.
* Otherwise the mesh data is replaced by the triangle mesh the flow returns; the material index
  and the smooth flag of every triangle come from the polygon it descends from (the library
  tracks the parent of every triangle through subdivision, flips and collapses, so this is exact);
  UV maps and other attributes are dropped (the report warns).
* The report warns when the flow ended before a round because the energy is undefined there
  (``stopped``: the round and the reason; the mesh is the state before it) and when triangles with
  a corner angle below the small angle could not be repaired (``bad_faces``).
"""

from array import array
from math import log10, radians
from time import perf_counter

import bpy
from bpy.props import BoolProperty, EnumProperty, FloatProperty, IntProperty, StringProperty

from . import core

# (identifier, label, description) of the six energies; the library's name is the lower-case id.
ENERGY_ITEMS = [
    ('COVARIANCE', "Eigenvalue (Eλ)",
     "Smallest eigenvalue of the normal covariance matrix (Eq. 3-4): the default, smooth, and "
     "cheap"),
    ('COMBINATORIAL', "Combinatorial (E_P)",
     "Minimum over bipartitions of the summed normal variance of the parts (Eq. 1-2): sharpest "
     "creases, expensive on high valence"),
    ('COMBINATORIAL_MAX', "Combinatorial Max (π max)",
     "Minimum over bipartitions of the largest normal distance within a part (Eq. 5)"),
    ('COVARIANCE_MAX', "Eigenvalue Max (λ max)",
     "Smallest over directions of the largest squared normal component (Eq. 6)"),
    ('INTRINSIC', "Intrinsic (width)",
     "Eigenvalue energy of the log-mapped normals at the vertex normal (App. B.5)"),
    ('INTRINSIC_MAX', "Intrinsic Max",
     "Branching-free intrinsic width (App. B.5.1)"),
]

METHOD_ITEMS = [
    ('LBFGS', "L-BFGS", "Quasi-Newton minimization with a weak Wolfe line search (the paper's choice)"),
    ('GRADIENT_DESCENT', "Gradient Descent", "Steepest descent with a line search"),
    ('FIXED_STEP', "Fixed Step",
     "x <- x - step * (mean edge length)^2 * gradient, no line search: no decrease guarantee"),
]

# Names of the attributes the preview writes.
DEFAULT_ENERGY_ATTRIBUTE = "developability"
ENERGY_VALUES_SUFFIX = "_energy"

# Attribute names that Blender or this add-on keep for themselves and that survive a mesh
# replacement as they are (or are rebuilt from the faces).
_STRUCTURAL = {"position", "sharp_face", "sharp_edge", "material_index", "uv_seam"}


def energy_name(identifier):
    return identifier.lower()


def _check_writable(obj):
    mesh = obj.data
    if mesh.library is not None:
        raise core.UnwrapFailed(
            "Mesh '{}' is linked from another file and cannot be changed here".format(mesh.name))
    if len(mesh.polygons) == 0:
        raise core.UnwrapFailed("Mesh '{}' has no faces".format(mesh.name))
    return mesh


def selected_vertices(mesh):
    """The vertices selected in the mesh (as stored: after leaving Edit mode, those selected
    there), as ``array('i')``."""
    flags = bytearray(len(mesh.vertices))
    mesh.vertices.foreach_get("select", memoryview(flags).cast('?'))
    out = array('i')
    i = flags.find(1)
    while i >= 0:
        out.append(i)
        i = flags.find(1, i + 1)
    return out


def _write_positions(mesh, positions):
    """Writes ``positions`` (3 floats per vertex) as the vertex coordinates."""
    flat = array('f', positions)
    attribute = mesh.attributes.get("position")
    if attribute is not None and attribute.domain == 'POINT' and attribute.data_type == 'FLOAT_VECTOR' \
            and 3 * len(attribute.data) == len(flat):
        attribute.data.foreach_set("vector", flat)
    else:
        mesh.vertices.foreach_set("co", flat)
    mesh.update()


def _user_attributes(mesh):
    """Names of the data that a replacement of the mesh drops: UV maps, color attributes and other
    named attributes that are not structural."""
    names = [layer.name for layer in mesh.uv_layers]
    names += [a.name for a in mesh.color_attributes]
    skip = set(names) | _STRUCTURAL
    names += [a.name for a in mesh.attributes if not a.name.startswith(".") and a.name not in skip]
    return names


def _replace_mesh(mesh, positions, faces, face_polygon):
    """Replaces the geometry of ``mesh`` by the triangle mesh ``(positions, faces)``; the material
    index and smooth flag of triangle ``t`` are those of polygon ``face_polygon[t]``. Returns the
    names of the data that was dropped."""
    n_old = len(mesh.polygons)
    material = array('i', bytes(4 * n_old))
    mesh.polygons.foreach_get("material_index", material)
    smooth = bytearray(n_old)
    mesh.polygons.foreach_get("use_smooth", memoryview(smooth).cast('?'))
    dropped = _user_attributes(mesh)
    coordinates = [(positions[i], positions[i + 1], positions[i + 2])
                   for i in range(0, len(positions), 3)]
    triangles = [(faces[i], faces[i + 1], faces[i + 2]) for i in range(0, len(faces), 3)]
    mesh.clear_geometry()
    mesh.from_pydata(coordinates, [], triangles)
    mesh.polygons.foreach_set("material_index", array('i', (material[p] for p in face_polygon)))
    mesh.polygons.foreach_set(
        "use_smooth", memoryview(bytearray(smooth[p] for p in face_polygon)).cast('?'))
    mesh.update()
    return dropped


# ----------------------------------------------------------------------------------------------
# Make Developable


def make_developable_object(obj, *, energy='COVARIANCE', method='LBFGS', iterations=500,
                            subdivisions=2, remesh=True, small_angle=radians(1.0),
                            max_displacement=0.5, include_boundary=True, omit_valence_three=True,
                            pin_selected=False, pinned_group="", triangulate=False,
                            write_energy=False, energy_attribute=DEFAULT_ENERGY_ATTRIBUTE,
                            progress=None):
    """Runs the developability flow on the mesh of ``obj`` and writes the result.

    ``progress(round, iteration, energy, gradient_norm)`` is the library's progress function (it
    may return False to stop). Returns a dict: ``report`` (that of ``ddg_tools.make_developable``),
    ``topology_changed``, ``replaced`` (the mesh data was replaced), ``dropped`` (names of the
    data dropped by a replacement), ``n_pinned``, ``seconds_*``. Raises ``core.UnwrapFailed`` with
    a message for the user.
    """
    ddg = core.native()
    mesh = _check_writable(obj)
    if mesh.is_editmode:
        raise core.UnwrapFailed("Internal error: the mesh is in Edit mode")
    if mesh.shape_keys is not None:
        raise core.UnwrapFailed(
            "Mesh '{}' has shape keys; the flow would move only the basis. Apply or remove "
            "the shape keys first".format(mesh.name))
    start = perf_counter()
    positions, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
    pinned = set()
    if pin_selected:
        pinned.update(selected_vertices(mesh))
    if pinned_group:
        pinned.update(core.group_vertices(obj, pinned_group))
    pinned = array('i', sorted(pinned))
    extracted = perf_counter()
    try:
        result = ddg.make_developable(
            positions, corner_vertices, polygon_starts, triangles,
            energy=energy_name(energy),
            include_boundary=include_boundary,
            omit_valence_three=omit_valence_three,
            method=method.lower(),
            iterations=iterations,
            subdivisions=subdivisions,
            max_displacement=max_displacement,
            pinned_vertices=pinned,
            remesh=remesh,
            small_angle=small_angle,
            return_faces=triangulate,
            progress=progress,
        )
    except ddg.DdgError as ex:
        raise core.UnwrapFailed(str(ex)) from ex
    solved = perf_counter()
    report = result.report
    dropped = []
    replaced = result.faces is not None
    if replaced:
        dropped = _replace_mesh(mesh, result.positions, result.faces, result.face_polygon)
    else:
        _write_positions(mesh, result.positions)
    written = perf_counter()
    if write_energy:
        show_developability(obj, energy=energy, include_boundary=include_boundary,
                            omit_valence_three=omit_valence_three, name=energy_attribute)
    return {
        "report": report,
        "topology_changed": report["topology_changed"],
        "replaced": replaced,
        "dropped": dropped,
        "n_pinned": len(pinned),
        "seconds_extract": extracted - start,
        "seconds_native": solved - extracted,
        "seconds_write": written - solved,
    }


def summary(stats):
    """One line for the status bar from the dict ``make_developable_object`` returns."""
    r = stats["report"]
    parts = ["Made developable: energy {:.3g} → {:.3g}".format(
        r["energy_initial"], r["energy_final"])]
    parts.append("| {} round{}, {} iterations".format(
        len(r["rounds"]), "" if len(r["rounds"]) == 1 else "s", r["iterations"]))
    if r["n_vertices"] != r["n_input_vertices"] or r["n_triangles"] != r["n_input_triangles"]:
        parts.append("| {} → {} triangles".format(r["n_input_triangles"], r["n_triangles"]))
    if r["cancelled"]:
        parts.append("| STOPPED early")
    stopped = r.get("stopped")
    if stopped is not None:
        parts.append("| ENDED before round {}: {}".format(stopped["round"] + 1, stopped["reason"]))
    if r.get("bad_faces"):
        parts.append("| {} thin triangles left".format(r["bad_faces"]))
    if stats["n_pinned"]:
        parts.append("| {} pinned".format(stats["n_pinned"]))
    if r["nonplanarity_after"] is not None:
        parts.append("| largest polygon non-planarity {:.3g} ({:.1f} % of a side)".format(
            r["nonplanarity_after"], 100.0 * r["nonplanarity_after_relative"]))
    parts.append("| {:.2f} s (native {:.2f} s)".format(
        stats["seconds_extract"] + stats["seconds_native"] + stats["seconds_write"],
        stats["seconds_native"]))
    return " ".join(parts)


class _FlowProgress:
    """Reports the flow to Blender's progress cursor."""

    def __init__(self, wm, rounds, iterations):
        self.wm = wm
        self.iterations = iterations
        self.total = rounds * iterations
        wm.progress_begin(0, max(self.total, 1))

    def __call__(self, round_, iteration, energy, gradient_norm):
        self.wm.progress_update(min(self.total, round_ * self.iterations + iteration))
        return True

    def end(self):
        self.wm.progress_end()


class OBJECT_OT_ddg_make_developable(bpy.types.Operator):
    bl_idname = "object.ddg_make_developable"
    bl_label = "Make Developable"
    bl_description = (
        "Move the vertices of the active mesh so that the surface becomes developable (flat "
        "faces fold along straight lines): a gradient flow of the developability energy with "
        "subdivision and small-angle remeshing. Method: Stein, Grinspun, Crane, Developability "
        "of Triangle Meshes (2018). Runs in place and can be undone. The operator is "
        "synchronous: the Blender window waits until it is done"
    )
    bl_options = {'REGISTER', 'UNDO'}

    energy: EnumProperty(
        name="Energy", items=ENERGY_ITEMS, default='COVARIANCE',
        description="The developability energy that is minimized",
    )
    method: EnumProperty(
        name="Method", items=METHOD_ITEMS, default='LBFGS',
        description="The optimizer",
    )
    iterations: IntProperty(
        name="Iterations", default=500, min=1, soft_max=5000,
        description="Optimizer iterations per round",
    )
    subdivisions: IntProperty(
        name="Subdivision Rounds", default=2, min=0, soft_max=4,
        description=(
            "Number of 4-1 subdivisions (one extra round of optimization after each). Each one "
            "quadruples the number of triangles and turns the mesh into a triangle mesh"
        ),
    )
    remesh: BoolProperty(
        name="Remesh", default=True,
        description=(
            "Repair triangles with a corner angle below the small angle by edge flips and "
            "collapses while the flow runs. Changes the topology when it acts"
        ),
    )
    small_angle: FloatProperty(
        name="Small Angle", default=radians(1.0), min=radians(0.01), max=radians(59.0),
        subtype='ANGLE', unit='ROTATION',
        description="Corner angle below which a triangle is repaired by remeshing",
    )
    max_displacement: FloatProperty(
        name="Step Cap", default=0.5, min=0.001, soft_max=2.0,
        description="No vertex moves farther than this many mean edge lengths in one step",
    )
    include_boundary: BoolProperty(
        name="Include Boundary", default=True,
        description="Count the vertices on the boundary of the mesh in the energy",
    )
    omit_valence_three: BoolProperty(
        name="Omit Valence 3", default=True,
        description=(
            "Leave out interior vertices with exactly three faces (always developable in the "
            "sense of the energy)"
        ),
    )
    pin_selected: BoolProperty(
        name="Pin Selected", default=False,
        description="Vertices selected in Edit mode do not move",
    )
    pinned_group: StringProperty(
        name="Pinned Group", default="",
        description="Vertex group (weights above 0) of vertices that do not move",
    )
    triangulate: BoolProperty(
        name="Triangulate", default=False,
        description=(
            "Replace the mesh by the triangle mesh of the algorithm even if no vertex or edge "
            "was added. Off: when the topology is unchanged the polygons are kept (n-gons "
            "follow the flow only as far as they are planar)"
        ),
    )
    write_energy: BoolProperty(
        name="Write Energy Preview", default=False,
        description=(
            "Afterwards write the energy of every vertex as color attribute 'developability' "
            "(and the values as float attribute 'developability_energy')"
        ),
    )

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False
        layout.prop(self, "energy")
        layout.prop(self, "method")
        layout.prop(self, "iterations")
        layout.prop(self, "subdivisions")
        layout.prop(self, "remesh")
        sub = layout.row()
        sub.active = self.remesh
        sub.prop(self, "small_angle")
        layout.prop(self, "max_displacement")
        layout.prop(self, "include_boundary")
        layout.prop(self, "omit_valence_three")
        layout.prop(self, "pin_selected")
        obj = context.active_object
        if obj is not None and obj.type == 'MESH':
            layout.prop_search(self, "pinned_group", obj, "vertex_groups")
        else:
            layout.prop(self, "pinned_group")
        layout.prop(self, "triangulate")
        layout.prop(self, "write_energy")

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        was_edit_mode = obj.mode == 'EDIT'
        if was_edit_mode:
            bpy.ops.object.mode_set(mode='OBJECT')
        progress = _FlowProgress(context.window_manager, self.subdivisions + 1, self.iterations)
        try:
            stats = make_developable_object(
                obj, energy=self.energy, method=self.method, iterations=self.iterations,
                subdivisions=self.subdivisions, remesh=self.remesh, small_angle=self.small_angle,
                max_displacement=self.max_displacement, include_boundary=self.include_boundary,
                omit_valence_three=self.omit_valence_three, pin_selected=self.pin_selected,
                pinned_group=self.pinned_group, triangulate=self.triangulate,
                write_energy=self.write_energy, progress=progress)
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        except KeyboardInterrupt:
            self.report({'WARNING'}, "Interrupted")
            return {'CANCELLED'}
        finally:
            progress.end()
            if was_edit_mode:
                bpy.ops.object.mode_set(mode='EDIT')
        if stats["replaced"] and stats["dropped"]:
            self.report({'WARNING'}, (
                "The mesh was replaced by a triangle mesh; UV maps and attributes were dropped: "
                + ", ".join(stats["dropped"])))
        elif stats["replaced"]:
            self.report({'INFO'}, "The mesh was replaced by a triangle mesh")
        report = stats["report"]
        if report["nonplanarity_after_relative"] is not None \
                and report["nonplanarity_after_relative"] > 0.01:
            self.report({'WARNING'}, (
                "Polygons with more than three corners are not planar any more (largest "
                "deviation {:.1f} % of a side, polygon {}). The flow made the triangulation "
                "developable; use Triangulate to keep it exact").format(
                    100.0 * report["nonplanarity_after_relative"], report["nonplanarity_polygon"]))
        for warning in report["warnings"]:
            self.report({'WARNING'}, warning)
        self.report({'INFO'}, summary(stats))
        return {'FINISHED'}


# ----------------------------------------------------------------------------------------------
# Show Developability


def _ramp(energy):
    """Color of a vertex energy: white below 1e-6, through yellow to red at 1e-1 and above, on a
    logarithmic scale."""
    t = (log10(max(energy, 1e-12)) + 6.0) / 5.0
    t = min(1.0, max(0.0, t))
    if t < 0.5:
        u = t / 0.5
        return (1.0, 1.0, 1.0 - u, 1.0)
    u = (t - 0.5) / 0.5
    return (1.0, 1.0 - u, 0.0, 1.0)


def show_developability(obj, *, energy='COVARIANCE', include_boundary=True, omit_valence_three=True,
                        name=DEFAULT_ENERGY_ATTRIBUTE):
    """Writes the developability energy of every vertex of ``obj`` as the point color attribute
    ``name`` (white to red, logarithmic from 1e-6 to 1e-1) and as the float attribute
    ``name + "_energy"``. Returns the library's result (``vertex_energy``, ``energy``, ...)."""
    ddg = core.native()
    mesh = _check_writable(obj)
    if not name:
        raise core.UnwrapFailed("The attribute name is empty")
    if name.startswith(".") or name in _STRUCTURAL:
        raise core.UnwrapFailed("'{}' is a reserved attribute name".format(name))
    for attribute_name, domain, data_type in ((name, 'POINT', 'FLOAT_COLOR'),
                                              (name + ENERGY_VALUES_SUFFIX, 'POINT', 'FLOAT')):
        existing = mesh.attributes.get(attribute_name)
        if existing is not None and (existing.domain != domain or existing.data_type != data_type):
            raise core.UnwrapFailed(
                "Mesh '{}' already has an attribute '{}' that is not a {} {} attribute".format(
                    mesh.name, attribute_name, domain.lower(), data_type.lower()))
    positions, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
    try:
        result = ddg.developability(
            positions, corner_vertices, polygon_starts, triangles,
            energy=energy_name(energy), include_boundary=include_boundary,
            omit_valence_three=omit_valence_three)
    except ddg.DdgError as ex:
        raise core.UnwrapFailed(str(ex)) from ex
    values = array('f', result.vertex_energy)
    colors = array('f')
    for e in result.vertex_energy:
        colors.extend(_ramp(e))
    color_attribute = mesh.color_attributes.get(name)
    if color_attribute is None:
        color_attribute = mesh.color_attributes.new(name, 'FLOAT_COLOR', 'POINT')
    color_attribute.data.foreach_set("color", colors)
    value_attribute = mesh.attributes.get(name + ENERGY_VALUES_SUFFIX)
    if value_attribute is None:
        value_attribute = mesh.attributes.new(name + ENERGY_VALUES_SUFFIX, 'FLOAT', 'POINT')
    value_attribute.data.foreach_set("value", values)
    # Creating the float attribute invalidates the Python reference `color_attribute` (and makes
    # the float the active generic attribute, which Attribute shading then draws black), so the
    # color attribute is made active and the render color by name.
    mesh.color_attributes.active_color_name = name
    mesh.color_attributes.render_color_index = mesh.color_attributes.active_color_index
    mesh.attributes.active = mesh.attributes[name]
    mesh.update()
    return result


class OBJECT_OT_ddg_show_developability(bpy.types.Operator):
    bl_idname = "object.ddg_show_developability"
    bl_label = "Show Developability"
    bl_description = (
        "Write the developability energy of every vertex of the active mesh as a color attribute "
        "(white: developable, red: far from it) and a float attribute, to look at in the "
        "viewport. The mesh is not changed. Energies: Stein, Grinspun, Crane (2018)"
    )
    bl_options = {'REGISTER', 'UNDO'}

    energy: EnumProperty(name="Energy", items=ENERGY_ITEMS, default='COVARIANCE')
    include_boundary: BoolProperty(name="Include Boundary", default=True)
    omit_valence_three: BoolProperty(name="Omit Valence 3", default=True)
    attribute_name: StringProperty(
        name="Attribute", default=DEFAULT_ENERGY_ATTRIBUTE,
        description="Name of the color attribute; the float values are stored under this name "
                    "plus '_energy'",
    )

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}

    def execute(self, context):
        obj = context.active_object
        was_edit_mode = obj.mode == 'EDIT'
        if was_edit_mode:
            bpy.ops.object.mode_set(mode='OBJECT')
        try:
            result = show_developability(
                obj, energy=self.energy, include_boundary=self.include_boundary,
                omit_valence_three=self.omit_valence_three, name=self.attribute_name)
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        finally:
            if was_edit_mode:
                bpy.ops.object.mode_set(mode='EDIT')
        self.report({'INFO'}, (
            "Developability energy {:.4g} (largest vertex {:.3g}); color attribute '{}' "
            "(set the viewport color to Attribute)").format(
                result.energy, result.report["max_vertex_energy"], self.attribute_name))
        return {'FINISHED'}


# ----------------------------------------------------------------------------------------------
# Developable Seams + UV


def developable_seams_uv(obj, *, energy='COVARIANCE', epsilon=1e-3, use_defect_rule=True,
                         defect_tolerance=0.1, length_weight=1.0, alignment_weight=0.0,
                         use_existing_seams=True, uv_name="Developable", lambda_=0.1,
                         positive_cones_only=True, mark_all_cuts=False, split_overlaps=True):
    """Marks the developable seams of the mesh of ``obj`` and writes a UV map from the unwrap
    pipeline (``core.unwrap_object``) that cuts along the marked seams.

    Returns the dict of ``core.unwrap_object`` with ``seams`` (the library's report),
    ``n_seams`` (edges of the developable cut, existing seams included), ``developable_marked``
    (of them newly marked), ``seconds_seams_native``. Vertices are not moved. Raises
    ``core.UnwrapFailed`` with a message for the user.
    """
    ddg = core.native()
    mesh = _check_writable(obj)
    if mesh.is_editmode:
        raise core.UnwrapFailed("Internal error: the mesh is in Edit mode")
    if not uv_name:
        raise core.UnwrapFailed("The UV map name is empty")
    positions, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
    seams = core.seam_edges(mesh) if use_existing_seams else None
    start = perf_counter()
    try:
        cut = ddg.developable_seams(
            positions, corner_vertices, polygon_starts, triangles,
            seam_edges=seams,
            energy=energy_name(energy),
            epsilon=epsilon,
            defect_tolerance=defect_tolerance if use_defect_rule else None,
            length_weight=length_weight,
            alignment_weight=alignment_weight,
        )
    except ddg.DdgError as ex:
        raise core.UnwrapFailed(str(ex)) from ex
    seconds_seams_native = perf_counter() - start
    marked, not_found = core.mark_seams(mesh, cut.seam_edges)
    layer = mesh.uv_layers.get(uv_name)
    if layer is None:
        layer = mesh.uv_layers.new(name=uv_name, do_init=False)
        if layer is None:
            raise core.UnwrapFailed("Could not create the UV map '{}'".format(uv_name))
    mesh.uv_layers.active = layer
    stats = core.unwrap_object(
        obj, lambda_=lambda_, use_seams=True, positive_cones_only=positive_cones_only,
        mark_cuts_as_seams=mark_all_cuts, use_stored_triangulation=False,
        split_overlaps=split_overlaps)
    stats["seams"] = cut.report
    stats["n_seams"] = len(cut.seam_edges) // 2
    stats["developable_marked"] = marked
    stats["developable_not_found"] = not_found
    stats["n_terminals"] = len(cut.terminals)
    stats["seconds_seams_native"] = seconds_seams_native
    return stats


def seams_summary(stats):
    """One line for the status bar from the dict ``developable_seams_uv`` returns."""
    r = stats["report"]
    overlap = r.get("overlapped_fraction")
    parts = ["Developable seams + UV: {} seam edges ({} newly marked), {} terminals".format(
        stats["n_seams"], stats["developable_marked"], stats["n_terminals"])]
    parts.append("| {} chart{}, {} cone{}, {} cut edges in all".format(
        r.get("n_charts", 0), "" if r.get("n_charts") == 1 else "s",
        stats["n_cones"], "" if stats["n_cones"] == 1 else "s", stats["n_cut_edges"]))
    if r.get("measured_triangles"):
        from math import degrees, exp
        parts.append("| angle error max {:.2f}°, mean {:.2f}°, scale {:.2f}x to {:.2f}x".format(
            degrees(r["measured_angle_error_max"]), degrees(r["measured_angle_error_mean"]),
            exp(r["measured_log_scale_min"]), exp(r["measured_log_scale_max"])))
    if "n_splits" in r:
        parts.append("| {} split{} ({} undone)".format(
            r["n_splits"], "" if r["n_splits"] == 1 else "s", r["n_splits_rejected"]))
    parts.append("| overlap {}".format("n/a" if overlap is None else "{:.1f} %".format(100.0 * overlap)))
    parts.append("| {} mirrored, {} collapsed triangles".format(
        r.get("n_flipped_uv_triangles", 0), r.get("n_collapsed_uv_triangles", 0)))
    parts.append("| UV map '{}'".format(stats["uv_layer"]))
    return " ".join(parts)


class UV_OT_ddg_developable_seams(bpy.types.Operator):
    bl_idname = "uv.ddg_developable_seams"
    bl_label = "Developable Seams + UV"
    bl_description = (
        "Cut the active mesh where it is not developable and make a UV map: the seams run "
        "through the vertices of high developability energy (Stein, Grinspun, Crane, 2018, "
        "Sec. 4.4) along polygon edges and are marked as seams; the UV map is then made by "
        "the Conformal Unwrap pipeline, cutting along all marked seams. Vertices are not moved"
    )
    bl_options = {'REGISTER', 'UNDO'}

    energy: EnumProperty(
        name="Energy", items=ENERGY_ITEMS, default='COVARIANCE',
        description="The energy whose values choose the vertices the cut goes through",
    )
    epsilon: FloatProperty(
        name="Threshold", default=1e-3, min=0.0, soft_max=0.1, precision=5, step=1,
        description=(
            "A vertex whose energy exceeds this is cut through (the paper uses 0.01 to 0.001). "
            "Lower: more and longer seams"
        ),
    )
    use_defect_rule: BoolProperty(
        name="Cut Through Cone Vertices", default=True,
        description=(
            "Also cut through interior vertices whose angle defect exceeds the tolerance "
            "(corners of a cube, cone tips): the energy cannot see a defect on a star whose "
            "normals lie on a great circle. Off is the paper's rule"
        ),
    )
    defect_tolerance: FloatProperty(
        name="Defect Tolerance", default=0.1, min=0.0, max=3.0, subtype='ANGLE', unit='ROTATION',
        description="Angle defect (2π minus the sum of the corner angles) above which a "
                    "vertex is cut through",
    )
    length_weight: FloatProperty(
        name="Length Weight", default=1.0, min=0.0, soft_max=10.0,
        description="Weight of the edge length in the cost of a seam edge",
    )
    alignment_weight: FloatProperty(
        name="Alignment Weight", default=0.0, min=0.0, soft_max=10.0,
        description=(
            "Weight of the misalignment of an edge with the creases (rulings) at its ends in "
            "the cost of a seam edge. 0: length only"
        ),
    )
    use_existing_seams: BoolProperty(
        name="Keep Marked Seams", default=True,
        description="The edges already marked as seams are cut as well, and the cut is planned around them",
    )
    uv_name: StringProperty(
        name="UV Map", default="Developable",
        description="Name of the UV map that receives the result (created if it does not exist, "
                    "then made active)",
    )
    lambda_: FloatProperty(
        name="Lambda", default=0.1, min=0.0, soft_min=0.01, soft_max=1.0, step=1, precision=3,
        description="Trade-off between area distortion and extra cones in the unwrap. Larger "
                    "values give fewer cones and more distortion",
    )
    positive_cones_only: BoolProperty(
        name="Positive Cones Only", default=True,
        description="Allow only cones with a positive angle in the unwrap (see Conformal Unwrap)",
    )
    mark_all_cuts: BoolProperty(
        name="Mark All Cuts", default=False,
        description="Also mark the extra cuts the unwrap made (to reach its cones, and to split "
                    "overlapping charts) as seams",
    )
    split_overlaps: BoolProperty(
        name="Split Overlapping Charts", default=True,
        description=(
            "A chart whose flattening lies on top of itself is cut in two where the layers "
            "meet and flattened again, repeatedly. A cut is kept only if it lowers the overlap "
            "and mirrors no triangle. Off: the charts stay as the developable cut made them"
        ),
    )

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False
        layout.prop(self, "energy")
        layout.prop(self, "epsilon")
        layout.prop(self, "use_defect_rule")
        sub = layout.row()
        sub.active = self.use_defect_rule
        sub.prop(self, "defect_tolerance")
        layout.prop(self, "length_weight")
        layout.prop(self, "alignment_weight")
        layout.prop(self, "use_existing_seams")
        layout.prop(self, "uv_name")
        layout.prop(self, "lambda_")
        layout.prop(self, "positive_cones_only")
        layout.prop(self, "split_overlaps")
        layout.prop(self, "mark_all_cuts")

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        was_edit_mode = obj.mode == 'EDIT'
        if was_edit_mode:
            bpy.ops.object.mode_set(mode='OBJECT')
        try:
            stats = developable_seams_uv(
                obj, energy=self.energy, epsilon=self.epsilon,
                use_defect_rule=self.use_defect_rule, defect_tolerance=self.defect_tolerance,
                length_weight=self.length_weight, alignment_weight=self.alignment_weight,
                use_existing_seams=self.use_existing_seams, uv_name=self.uv_name,
                lambda_=self.lambda_, positive_cones_only=self.positive_cones_only,
                mark_all_cuts=self.mark_all_cuts, split_overlaps=self.split_overlaps)
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        finally:
            if was_edit_mode:
                bpy.ops.object.mode_set(mode='EDIT')
        for warning in stats["report"].get("warnings", ()):
            self.report({'WARNING'}, warning)
        if stats["developable_not_found"]:
            self.report({'WARNING'}, "{} seam edges are not edges of the mesh".format(
                stats["developable_not_found"]))
        self.report({'INFO'}, seams_summary(stats))
        return {'FINISHED'}
