# SPDX-License-Identifier: MIT
"""The cut flow in Blender: state attributes, preview, flow settings, insertion with bmesh.

"Variational Surface Cutting" (Sharp, Crane, ACM TOG 2018) is implemented clean-room from the
text of the paper in ``crates/ddg-variational-cuts``; ``ddg_tools.CutSession`` (``crates/ddg-py``)
runs it on the mesh buffers this add-on extracts. Everything here goes through Blender's data API
(``foreach_get`` / ``foreach_set``), the standard library, ``mathutils`` and ``bmesh``; there are
no operator calls and no numpy. The mesh must not be in Edit mode (the operators leave Edit mode
before calling in and re-enter it afterwards).

Data model.

* The state of the flow is one float attribute per region on the face-corner domain,
  ``ddg_cut_region_0``, ``ddg_cut_region_1``, ... (the two sides of a seam can differ). For two
  regions the second is the negative of the first.
* The preview is a child object of the mesh object that holds the cut as wire edges.
* Existing seams are fixed cuts: the flow cannot cross them, and they stay seams.
* Applying inserts the cut with bmesh: the vertices of the cut are put on the edges with
  ``bmesh.utils.edge_split`` and the polygons are split along the cut with
  ``bmesh.utils.face_split``, so that UV maps, vertex groups, attributes and materials are
  interpolated or inherited by Blender's own operations; only the points inside polygons, which
  Blender's operations do not interpolate, are interpolated here.
"""

import math
from array import array
from time import perf_counter

import bmesh
import bpy
from mathutils import Vector
from mathutils.kdtree import KDTree

from . import core

ATTRIBUTE_PREFIX = "ddg_cut_region_"
PREVIEW_KEY = "ddg_cut_preview"
LAST_FLOW_KEY = "ddg_cut_last_flow"

# The ripple wavelength of the length weight is a fraction of the square root of the surface
# area (the library works on the surface scaled to unit area).


class CutFailed(core.UnwrapFailed):
    """The cut operation could not be done; the message is meant for the user."""


# ----------------------------------------------------------------------------------------------
# The native module


def native():
    ddg = core.native()
    if not hasattr(ddg, "CutSession"):
        raise CutFailed(
            "The bundled 'ddg_tools' module has no CutSession (it is older than this add-on); "
            "reinstall the extension for this platform.")
    return ddg


def _translate(ex, hint=""):
    """A readable message from an exception of the native module."""
    message = str(ex)
    return CutFailed(message + hint)


# ----------------------------------------------------------------------------------------------
# State attributes


def state_names(mesh):
    """The names of the region attributes of ``mesh`` in order (a gap ends the list)."""
    names = []
    while mesh.attributes.get(ATTRIBUTE_PREFIX + str(len(names))) is not None:
        names.append(ATTRIBUTE_PREFIX + str(len(names)))
    return names


def has_state(mesh):
    return mesh.attributes.get(ATTRIBUTE_PREFIX + "0") is not None


def read_state(mesh):
    """The state of the mesh: ``(n_regions, array('d'))`` with ``n_regions * n_corners`` values,
    region after region."""
    names = state_names(mesh)
    if not names:
        raise CutFailed(
            "Mesh '{}' has no cut state (no attribute '{}0'). Run 'Cut: Start' first.".format(
                mesh.name, ATTRIBUTE_PREFIX))
    n_corners = len(mesh.loops)
    flat = array('d')
    for name in names:
        attribute = mesh.attributes[name]
        if attribute.domain != 'CORNER' or attribute.data_type != 'FLOAT':
            raise CutFailed(
                "Attribute '{}' must be a float attribute on the face-corner domain, but it is "
                "{} on {}. Remove it and run 'Cut: Start' again.".format(
                    name, attribute.data_type, attribute.domain))
        values = array('f', bytes(4 * n_corners))
        attribute.data.foreach_get("value", values)
        flat.extend(array('d', values))
    return len(names), flat


def write_state(mesh, session):
    """Writes the state of ``session`` to the region attributes, creating them as needed, and
    removes attributes of regions the session does not have."""
    for k in range(session.n_regions):
        name = ATTRIBUTE_PREFIX + str(k)
        attribute = mesh.attributes.get(name)
        if attribute is not None and (attribute.domain != 'CORNER'
                                      or attribute.data_type != 'FLOAT'):
            mesh.attributes.remove(attribute)
            attribute = None
        if attribute is None:
            attribute = mesh.attributes.new(name, 'FLOAT', 'CORNER')
        attribute.data.foreach_set("value", array('f', session.state_row(k)))
    k = session.n_regions
    while mesh.attributes.get(ATTRIBUTE_PREFIX + str(k)) is not None:
        mesh.attributes.remove(mesh.attributes[ATTRIBUTE_PREFIX + str(k)])
        k += 1


def remove_state(mesh):
    """Removes the region attributes (all of them, also behind a gap)."""
    for name in [a.name for a in mesh.attributes if a.name.startswith(ATTRIBUTE_PREFIX)]:
        mesh.attributes.remove(mesh.attributes[name])
    if LAST_FLOW_KEY in mesh:
        del mesh[LAST_FLOW_KEY]


# ----------------------------------------------------------------------------------------------
# Sessions


def check_mesh(obj):
    mesh = obj.data
    if obj.type != 'MESH':
        raise CutFailed("The active object is not a mesh")
    if mesh.is_editmode:
        raise CutFailed("Internal error: the mesh is in Edit mode")
    if mesh.library is not None:
        raise CutFailed(
            "Mesh '{}' is linked from another file ({}) and cannot be changed here".format(
                mesh.name, mesh.library.filepath))
    if len(mesh.polygons) == 0:
        raise CutFailed("Mesh '{}' has no faces".format(mesh.name))
    return mesh


def make_session(obj, start, n_regions=2):
    """A ``ddg_tools.CutSession`` on the mesh of ``obj`` (existing seams are fixed cuts).
    ``start`` is a dict of the start keywords of the constructor. Returns ``(session, seconds)``
    with the time spent extracting the buffers."""
    ddg = native()
    mesh = check_mesh(obj)
    begin = perf_counter()
    positions, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
    seams = core.seam_edges(mesh)
    extracted = perf_counter()
    try:
        session = ddg.CutSession(
            positions, corner_vertices, polygon_starts, triangles,
            n_regions=n_regions, seam_edges=seams if len(seams) else None, **start)
    except ddg.DdgError as ex:
        raise _translate(ex) from ex
    return session, extracted - begin


def selected_faces(mesh):
    """One byte (0 or 1) per face: whether the face is selected."""
    flags = bytearray(len(mesh.polygons))
    mesh.polygons.foreach_get("select", memoryview(flags).cast('?'))
    return flags


def start_labels_from_selection(mesh):
    """``array('i')`` with 0 for the selected faces and 1 for the others."""
    flags = selected_faces(mesh)
    n_selected = sum(flags)
    if n_selected == 0:
        raise CutFailed("No faces are selected. Select the faces that should be region 0, "
                        "or start from a geodesic circle.")
    if n_selected == len(flags):
        raise CutFailed("All faces are selected; select some faces, but not all of them.")
    return array('i', (0 if f else 1 for f in flags))


def start_cut(obj, *, mode, radius=0.0, vertex=-1, n_regions=2, active_vertex=None,
              selected_vertices=()):
    """Creates the state of a cut on ``obj`` and its preview.

    ``mode`` is ``'SELECTED'`` (the selected faces against the rest) or ``'CIRCLE'`` (the
    geodesic disk of ``radius`` object units, 0: a tenth of the bounding diagonal, around ``vertex``; with ``vertex < 0`` the
    ``active_vertex``, or the only one of ``selected_vertices``).

    Returns a dict with the session's sizes and timings.
    """
    mesh = check_mesh(obj)
    if mode == 'SELECTED':
        start = {"polygon_labels": start_labels_from_selection(mesh)}
        if n_regions != 2:
            raise CutFailed(
                "{} regions asked for, but the cut flow handles exactly 2 regions so far (N "
                "regions are not available yet).".format(n_regions))
    elif mode == 'CIRCLE':
        if vertex < 0:
            if active_vertex is not None:
                vertex = active_vertex
            elif len(selected_vertices) == 1:
                vertex = selected_vertices[0]
            else:
                raise CutFailed(
                    "There is no active vertex. Select the vertex the circle is centred on in "
                    "Edit mode (or set the Vertex property).")
        if vertex >= len(mesh.vertices):
            raise CutFailed("Vertex {} does not exist: the mesh has {} vertices".format(
                vertex, len(mesh.vertices)))
        if radius < 0.0:
            raise CutFailed("The radius must not be negative (it is {})".format(radius))
        if radius == 0.0:
            radius = 0.1 * _bounding_diagonal(mesh)
            if not radius > 0.0:
                raise CutFailed("The mesh has no extent, so no default radius can be chosen")
        start = {"circle_vertex": vertex, "circle_radius": radius}
    else:
        raise CutFailed("Unknown start mode {!r}".format(mode))
    session, seconds_extract = make_session(obj, start, n_regions)
    begin = perf_counter()
    write_state(mesh, session)
    if LAST_FLOW_KEY in mesh:
        del mesh[LAST_FLOW_KEY]
    mesh.update()
    written = perf_counter()
    update_preview(obj, session)
    previewed = perf_counter()
    return {
        "n_regions": session.n_regions,
        "n_corners": session.n_corners,
        "n_polygons": session.n_polygons,
        "vertex": start.get("circle_vertex"),
        "radius": start.get("circle_radius"),
        "warnings": list(session.warnings),
        "seconds_extract": seconds_extract,
        "seconds_write": written - begin,
        "seconds_preview": previewed - written,
    }


# ----------------------------------------------------------------------------------------------
# Preview


def preview_objects(obj):
    """The preview objects of ``obj`` (children marked as previews)."""
    return [child for child in obj.children if child.get(PREVIEW_KEY)]


def remove_preview(obj):
    """Removes the preview objects of ``obj`` and their meshes. Returns how many."""
    removed = 0
    for child in preview_objects(obj):
        mesh = child.data
        bpy.data.objects.remove(child, do_unlink=True)
        if mesh is not None and mesh.users == 0:
            bpy.data.meshes.remove(mesh)
        removed += 1
    return removed


def update_preview(obj, session):
    """Replaces the preview of ``obj`` by the current cut of ``session``: a wire object named
    after the mesh object, parented to it (so that it follows its transform). Returns the new
    object."""
    ddg = native()
    try:
        positions, starts, closed = session.preview()
    except ddg.DdgError as ex:
        raise _translate(ex) from ex
    remove_preview(obj)
    name = obj.name + " cut"
    mesh = bpy.data.meshes.new(name)
    n_points = starts[-1] if len(starts) else 0
    edges = array('i')
    for line in range(len(closed)):
        first, last = starts[line], starts[line + 1]
        for point in range(first, last - 1):
            edges.extend((point, point + 1))
        if closed[line] and last - first > 2:
            edges.extend((last - 1, first))
    mesh.vertices.add(n_points)
    if n_points:
        mesh.vertices.foreach_set("co", array('f', positions))
    mesh.edges.add(len(edges) // 2)
    if len(edges):
        mesh.edges.foreach_set("vertices", edges)
    mesh.update()
    mesh.validate()
    preview = bpy.data.objects.new(name, mesh)
    collections = obj.users_collection
    (collections[0] if collections else bpy.context.scene.collection).objects.link(preview)
    preview.parent = obj
    preview.display_type = 'WIRE'
    preview.show_in_front = True
    preview.hide_render = True
    preview.hide_select = True
    preview[PREVIEW_KEY] = True
    return preview


# ----------------------------------------------------------------------------------------------
# Penalty, symmetry, visibility


def group_weights(obj, name):
    """The weight of every vertex in the vertex group ``name`` (0 outside the group), as
    ``array('d')``."""
    group = obj.vertex_groups.get(name)
    if group is None:
        raise CutFailed("Vertex group '{}' does not exist on object '{}'".format(name, obj.name))
    weights = array('d', bytes(8 * len(obj.data.vertices)))
    index = group.index
    for vertex in obj.data.vertices:
        for element in vertex.groups:
            if element.group == index:
                weights[vertex.index] = element.weight
    return weights


def mirror_classes(mesh, axis, tolerance):
    """Symmetry classes for mirroring about the plane ``axis = 0`` of the object (``axis`` in
    ``'X'``, ``'Y'``, ``'Z'``): vertices that are mirror images of each other within
    ``tolerance`` (object units) share a class; a vertex on the plane is its own mirror image.

    Returns ``(classes, unmatched)``: ``array('I')`` with a class per vertex, and the list of
    vertices without a partner (each in a class of its own).
    """
    k = "XYZ".index(axis)
    n = len(mesh.vertices)
    positions = array('f', bytes(12 * n))
    mesh.vertices.foreach_get("co", positions)
    points = [Vector(positions[3 * i:3 * i + 3]) for i in range(n)]
    tree = KDTree(n)
    for i, p in enumerate(points):
        tree.insert(p, i)
    tree.balance()
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    unmatched = []
    for i, p in enumerate(points):
        mirrored = p.copy()
        mirrored[k] = -mirrored[k]
        _, partner, distance = tree.find(mirrored)
        if partner is not None and distance <= tolerance:
            a, b = find(i), find(partner)
            if a != b:
                parent[max(a, b)] = min(a, b)
        else:
            unmatched.append(i)
    return array('I', (find(i) for i in range(n))), unmatched


def vertex_visibility(obj, context, camera_position, tolerance_scale=1e-4):
    """Whether the camera sees each vertex of ``obj``: ``array('d')`` with 1 where a ray from
    ``camera_position`` (world space) reaches the vertex without hitting anything else in the
    scene, 0 where something is in the way. The rays are cast with the scene's ray cast against
    the evaluated geometry of every visible object (the mesh itself included)."""
    mesh = obj.data
    matrix = obj.matrix_world
    depsgraph = context.evaluated_depsgraph_get()
    scene = context.scene
    n = len(mesh.vertices)
    positions = array('f', bytes(12 * n))
    mesh.vertices.foreach_get("co", positions)
    visible = array('d', bytes(8 * n))
    origin = Vector(camera_position)
    for i in range(n):
        target = matrix @ Vector(positions[3 * i:3 * i + 3])
        direction = target - origin
        distance = direction.length
        if distance == 0.0:
            visible[i] = 1.0
            continue
        direction /= distance
        # Cast slightly beyond the vertex; anything that is hit clearly before it is in the way.
        hit, location, _normal, _index, _object, _matrix = scene.ray_cast(
            depsgraph, origin, direction, distance=distance * (1.0 + tolerance_scale) + 1e-6)
        if not hit or (location - origin).length >= distance * (1.0 - tolerance_scale) - 1e-6:
            visible[i] = 1.0
    return visible


# ----------------------------------------------------------------------------------------------
# The flow


ENERGIES = {'DIRICHLET': "dirichlet", 'SCALED_DIRICHLET': "scaled_dirichlet", 'HENCKY': "hencky"}


def flow_controls(obj, context, session, settings):
    """The keywords of ``session.run`` for the settings of the Flow operator, and notes.

    ``settings`` is a dict-like with the property names of the operator (see
    ``cut_operators.UV_OT_ddg_cut_flow``). Returns ``(controls, info)``; ``info`` holds the
    resolved length weight and messages for the user.
    """
    ddg = native()
    info = {"warnings": []}
    energy = ENERGIES[settings["energy"]]
    controls = {"energy": energy}
    if settings["stages"] == 'PAPER':
        controls["stages"] = "paper"
    if energy == "scaled_dirichlet":
        controls["feature_r_eps"] = settings["feature_r_eps"]
        controls["feature_smoothing"] = settings["feature_smoothing"]
    if settings["use_uniform_curvature"]:
        controls["curvature_uniform"] = settings["uniform_curvature"]
    controls["distortion_weights"] = tuple(settings["distortion_weights"])
    controls["area_weights"] = tuple(settings["area_weights"])

    step_rule = settings["step_rule"].lower()
    controls["step_rule"] = step_rule
    if step_rule == "adaptive":
        controls["time_scale"] = settings["time_scale"]
        controls["max_motion"] = settings["max_motion"]
        controls["tau_max"] = settings["tau_max"]
    elif step_rule == "capped":
        controls["max_motion"] = settings["max_motion"]
        controls["tau_max"] = settings["tau_max"]
    else:
        controls["tau"] = settings["tau"]
    for key, table in (
            ("flux_recovery", {'NODAL': "nodal", 'SMOOTHED': "smoothed"}),
            ("length_step", {'ON_CUT': "on_cut", 'ON_MESH': "on_mesh"}),
            ("extension", {'LOCAL': "local", 'HARMONIC': "harmonic"}),
            ("speed_scaling", {'MEAN_GRADIENT': "mean_gradient", 'UNIT': "unit"}),
            ("redistance", {'KEEP_ADJACENT': "keep_adjacent", 'OFF': "off"}),
            ("penalty_rule", {'GRADIENT': "gradient", 'PAPER_DIAGONAL': "paper_diagonal"})):
        value = settings[key]
        if value != 'DEFAULT':
            controls[key] = table[value]
    controls["min_crossing_offset"] = settings["min_crossing_offset"]

    mesh = obj.data
    if settings["symmetry"] != 'NONE':
        scale = _bounding_diagonal(mesh)
        tolerance = settings["symmetry_tolerance"] * scale
        classes, unmatched = mirror_classes(mesh, settings["symmetry"], tolerance)
        controls["symmetry"] = classes
        info["n_unmatched"] = len(unmatched)
        info["unmatched"] = unmatched
        if unmatched:
            info["warnings"].append(
                "{} of {} vertices have no mirror image about the {} axis within {:.3g} "
                "(first: {}); they are not constrained by the symmetry".format(
                    len(unmatched), len(mesh.vertices), settings["symmetry"], tolerance,
                    ", ".join(str(v) for v in unmatched[:8])
                    + (", ..." if len(unmatched) > 8 else "")))

    # The length weight: explicit, or from the roughening rate of the current state.
    if settings["length_mode"] == 'EXPLICIT':
        alpha = settings["length_weight"]
        info["roughening_rate"] = None
    else:
        wavelength = settings["ripple_wavelength"]
        if not wavelength > 0.0:
            raise CutFailed("The ripple wavelength must be larger than 0 (it is {})".format(
                wavelength))
        rate_controls = {k: v for k, v in controls.items()
                         if k not in ("symmetry", "min_crossing_offset")}
        try:
            omega = session.roughening_rate(**rate_controls)
        except ddg.DdgError as ex:
            raise _translate(ex) from ex
        alpha = ddg.length_weight_for_ripples(omega, wavelength)
        info["roughening_rate"] = omega
    controls["length_weight"] = alpha
    info["length_weight"] = alpha

    # The penalty eta: strength times the length weight times a weight per vertex in 0..1.
    eta = None
    if settings["penalty_group"]:
        weights = group_weights(obj, settings["penalty_group"])
        eta = array('d', (settings["penalty_strength"] * alpha * w for w in weights))
    if settings["hide_from_camera"]:
        camera = context.scene.camera
        if camera is None:
            raise CutFailed("The scene has no active camera. Add a camera and make it active, "
                            "or switch off 'Hide From Camera'.")
        begin = perf_counter()
        position = camera.matrix_world.translation
        visibility = vertex_visibility(obj, context, position)
        local = obj.matrix_world.inverted_safe() @ position
        positions, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
        try:
            facing = ddg.visibility_penalty(
                positions, corner_vertices, polygon_starts, triangles, tuple(local),
                visibility=visibility)
        except ddg.DdgError as ex:
            raise _translate(ex) from ex
        info["seconds_visibility"] = perf_counter() - begin
        info["n_hidden"] = sum(1 for v in visibility if v == 0.0)
        extra = array('d', (settings["hide_strength"] * alpha * f for f in facing))
        eta = extra if eta is None else array('d', (a + b for a, b in zip(eta, extra)))
    if eta is not None:
        controls["penalty"] = eta
    return controls, info


def _bounding_diagonal(mesh):
    n = len(mesh.vertices)
    positions = array('f', bytes(12 * n))
    mesh.vertices.foreach_get("co", positions)
    lows = [min(positions[k::3]) for k in range(3)]
    highs = [max(positions[k::3]) for k in range(3)]
    return math.sqrt(sum((h - l) ** 2 for l, h in zip(lows, highs)))


def run_flow(obj, context, settings):
    """Reads the state of ``obj``, runs the flow, writes the state back and updates the preview.
    Returns the report dict of ``CutSession.run`` plus notes (see keys ``info``, ``seconds_*``)."""
    ddg = native()
    mesh = check_mesh(obj)
    begin = perf_counter()
    n_regions, flat = read_state(mesh)
    session, seconds_extract = make_session(obj, {"corner_values": flat}, n_regions)
    state_read = perf_counter()
    controls, info = flow_controls(obj, context, session, settings)
    prepared = perf_counter()
    stop = {"max_steps": settings["steps"]}
    if settings["use_stop_rule"]:
        stop["motion_tolerance"] = settings["motion_tolerance"]
        stop["energy_tolerance"] = settings["energy_tolerance"]
        stop["energy_window"] = settings["energy_window"]
    else:
        stop["motion_tolerance"] = 0.0
        stop["energy_tolerance"] = 0.0
    try:
        report = session.run(**controls, **stop)
    except ddg.DdgError as ex:
        raise _translate(ex) from ex
    solved = perf_counter()
    write_state(mesh, session)
    mesh.update()
    written = perf_counter()
    update_preview(obj, session)
    previewed = perf_counter()
    summary = {k: v for k, v in report.items() if k != "steps"}
    summary["length_scale"] = session.length_scale
    if info.get("roughening_rate") is not None:
        summary["length_weight_roughening_rate"] = info["roughening_rate"]
    mesh[LAST_FLOW_KEY] = {k: v for k, v in summary.items()
                           if isinstance(v, (int, float, str)) and not isinstance(v, bool)}
    report["info"] = info
    report["warnings"] = list(session.warnings) + info["warnings"]
    report["seconds_extract"] = state_read - begin
    report["seconds_prepare"] = prepared - state_read
    report["seconds_native"] = solved - prepared
    report["seconds_write"] = written - solved
    report["seconds_preview"] = previewed - written
    return report


def flow_summary(report):
    """One line for Blender's status bar."""
    reason = {"max_steps": "step limit", "converged": "converged",
              "energy_stalled": "energy stalled", "cut_vanished": "the cut vanished"}[
        report["stop_reason"]]
    return ("Cut flow: {} step{} ({}) | energy {:.5g} -> {:.5g} | cut length {:.5g} "
            "(was {:.5g}) | length weight {:.4g}").format(
        report["steps_taken"], "" if report["steps_taken"] == 1 else "s", reason,
        report["energy_before"], report["energy_after"], report["cut_length_after"],
        report["cut_length_before"], report["length_weight"])


# ----------------------------------------------------------------------------------------------
# Inserting the cut with bmesh


def _loop_layers(bm):
    """The loop data layers whose values can be interpolated or copied: a list of
    ``(layer, interpolate)``; ``interpolate`` is False for layers that are copied from the
    nearer corner."""
    layers = []
    for collection, interpolate in ((bm.loops.layers.uv, True), (bm.loops.layers.float, True),
                                    (bm.loops.layers.float_vector, True),
                                    (bm.loops.layers.float_color, True),
                                    (bm.loops.layers.color, True), (bm.loops.layers.int, False),
                                    (bm.loops.layers.bool, False),
                                    (bm.loops.layers.string, False)):
        for layer in collection.values():
            layers.append((layer, interpolate))
    return layers


def _loop_value(loop, layer):
    value = loop[layer]
    if hasattr(value, "uv"):                 # BMLoopUV
        return Vector(value.uv)
    if isinstance(value, (bytes, bytearray)):
        return bytes(value)
    if isinstance(value, (int, float, bool)):
        return value
    return Vector(value)                     # float_vector, colors


def _set_loop_value(loop, layer, value):
    if hasattr(loop[layer], "uv"):
        loop[layer].uv = value
    else:
        loop[layer] = value


def _blend(a, b, t, interpolate):
    if not interpolate:
        return a if t < 0.5 else b
    return a * (1.0 - t) + b * t


def _vertex_layers(bm):
    layers = []
    for collection in (bm.verts.layers.float, bm.verts.layers.float_vector,
                       bm.verts.layers.float_color, bm.verts.layers.color,
                       bm.verts.layers.shape):
        for layer in collection.values():
            layers.append((layer, True))
    for collection in (bm.verts.layers.int, bm.verts.layers.bool, bm.verts.layers.string):
        for layer in collection.values():
            layers.append((layer, False))
    return layers


def _interpolate_vertex(bm, new, a, b, t, layers):
    """Vertex data of ``new``, a point at parameter ``t`` on the segment ``a``-``b``: vertex
    group weights and the vertex layers of ``a`` and ``b`` blended."""
    deform = bm.verts.layers.deform.active
    if deform is not None:
        wa, wb = dict(a[deform].items()), dict(b[deform].items())
        for group in set(wa) | set(wb):
            weight = (1.0 - t) * wa.get(group, 0.0) + t * wb.get(group, 0.0)
            if weight > 0.0:
                new[deform][group] = weight
    for layer, interpolate in layers:
        value = a[layer], b[layer]
        if isinstance(value[0], (int, float, bool, bytes)):
            new[layer] = _blend(value[0], value[1], t, interpolate)
        else:
            new[layer] = _blend(Vector(value[0]), Vector(value[1]), t, interpolate)


def insert_cut(mesh, cut, *, mark_seams=True):
    """Inserts the cut described by ``cut`` (a ``ddg_tools.CutDescription``) into ``mesh``, which
    must be the mesh the description was made on (same polygons and corners), and marks the
    edges of the cut as seams.

    The points on edges are made by splitting the edges (``bmesh.utils.edge_split``, which
    interpolates vertex and corner data), the polygons are split along the paths
    (``bmesh.utils.face_split``, which gives each new polygon the data of the old one); the points
    inside polygons get interpolated vertex-group weights, vertex data and corner data (UVs,
    attributes) from the corners of the diagonal they lie on. Where this differs from
    ``CutSession.insert`` (the Rust reference): the new polygons are not triangulated (only their
    shape is the same), and vertex numbers follow Blender's order of creation.

    Returns a dict with ``n_vertices``, ``n_polygons``, ``n_edge_points``, ``n_seams_marked``.
    """
    if mesh.shape_keys is not None:
        raise CutFailed(
            "Mesh '{}' has shape keys; applying a cut would change its vertices. Apply or remove "
            "the shape keys first.".format(mesh.name))
    bm = bmesh.new()
    try:
        bm.from_mesh(mesh)
        stats = _insert_cut(bm, cut, mark_seams)
        bm.to_mesh(mesh)
    finally:
        bm.free()
    mesh.update()
    return stats


def _insert_cut(bm, cut, mark_seams):
    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.faces.ensure_lookup_table()
    n_points = cut.n_edge_points
    n_paths = cut.n_paths
    face_starts = []
    corner = 0
    for face in bm.faces:
        face_starts.append(corner)
        corner += len(face.loops)

    loop_layers = _loop_layers(bm)
    vertex_layers = _vertex_layers(bm)

    # What the paths need to know about the old polygons, before anything is split.
    path_polygon = cut.path_polygon
    polygons_with_paths = sorted(set(path_polygon))
    old_corners = {}
    old_faces = {}
    for p in polygons_with_paths:
        face = old_faces[p] = bm.faces[p]
        old_corners[p] = [(loop.vert, [_loop_value(loop, layer) for layer, _ in loop_layers])
                          for loop in face.loops]

    # ----- points on edges -----
    edge_vertices = cut.edge_vertices
    edge_t = cut.edge_t
    by_edge = {}
    for i in range(n_points):
        key = (edge_vertices[2 * i], edge_vertices[2 * i + 1])
        by_edge.setdefault(key, []).append(i)
    point_vertex = [None] * n_points
    # The index tables are outdated as soon as the first edge is split.
    ends = {key: (bm.verts[key[0]], bm.verts[key[1]]) for key in by_edge}
    for (a, b), ids in by_edge.items():
        va, vb = ends[(a, b)]
        edge = bm.edges.get((va, vb))
        if edge is None:
            raise CutFailed("Internal error: the cut crosses ({}, {}), which is not an edge".format(
                a, b))
        ids = sorted(ids, key=lambda i: (edge_t[i], i))
        groups = []
        for i in ids:                         # the two sides of a seam can share a point
            if groups and edge_t[groups[-1][0]] == edge_t[i]:
                groups[-1].append(i)
            else:
                groups.append([i])
        current, from_vert, t_before = edge, va, 0.0
        for group in groups:
            t = edge_t[group[0]]
            fraction = (t - t_before) / (1.0 - t_before)
            new_edge, new_vert = bmesh.utils.edge_split(current, from_vert, fraction)
            for i in group:
                point_vertex[i] = new_vert
            current, from_vert, t_before = new_edge, new_vert, t

    # ----- paths -----
    path_start, path_end = cut.path_start, cut.path_end
    interior_starts = cut.path_interior_starts
    interior_corners, interior_t = cut.interior_corners, cut.interior_t
    pieces = {p: [old_faces[p]] for p in polygons_with_paths}
    cut_segments = []
    for k in range(n_paths):
        p = path_polygon[k]
        if path_start[k] < 0 or path_end[k] < 0:
            raise CutFailed("Paths that end at a junction of three or more regions cannot be "
                            "inserted yet")
        start, end = point_vertex[path_start[k]], point_vertex[path_end[k]]
        corners = old_corners[p]
        coords, params = [], []
        for j in range(interior_starts[k], interior_starts[k + 1]):
            c0 = interior_corners[2 * j] - face_starts[p]
            c1 = interior_corners[2 * j + 1] - face_starts[p]
            t = interior_t[j]
            coords.append(corners[c0][0].co * (1.0 - t) + corners[c1][0].co * t)
            params.append((c0, c1, t))
        face = next((f for f in pieces[p] if start in f.verts and end in f.verts), None)
        if face is None:
            raise CutFailed(
                "Internal error: the path {} of polygon {} does not lie in one piece of it".format(
                    k, p))
        before = set(face.verts)
        result = bmesh.utils.face_split(face, start, end, coords=coords)
        if not result or result[0] is None:
            raise CutFailed("Internal error: polygon {} could not be split along path {}".format(
                p, k))
        new_face = result[0]
        pieces[p].append(new_face)
        chain = [start]
        if coords:
            candidates = [v for v in new_face.verts if v not in before]
            for (c0, c1, t), co in zip(params, coords):
                vert = min(candidates, key=lambda v: (v.co - co).length_squared)
                candidates.remove(vert)
                chain.append(vert)
                _interpolate_vertex(bm, vert, corners[c0][0], corners[c1][0], t, vertex_layers)
                for loop in vert.link_loops:
                    for index, (layer, interpolate) in enumerate(loop_layers):
                        a, b = corners[c0][1][index], corners[c1][1][index]
                        _set_loop_value(loop, layer, _blend(a, b, t, interpolate))
        chain.append(end)
        cut_segments.extend(zip(chain, chain[1:]))

    n_marked = 0
    if mark_seams:
        for u, v in cut_segments:
            edge = bm.edges.get((u, v))
            if edge is None:
                raise CutFailed("Internal error: a segment of the cut is not an edge")
            if not edge.seam:
                edge.seam = True
                n_marked += 1
    return {
        "n_vertices": len(bm.verts),
        "n_polygons": len(bm.faces),
        "n_edge_points": n_points,
        "n_seams_marked": n_marked,
    }


def apply_cut(obj, *, mark_seams=True):
    """Applies the cut stored in the state attributes of ``obj``: inserts it into the mesh
    (``insert_cut``), marks it as seams, and removes the preview and the state attributes.
    Returns a dict with the insertion's counts and the time spent."""
    ddg = native()
    mesh = check_mesh(obj)
    begin = perf_counter()
    n_regions, flat = read_state(mesh)
    session, seconds_extract = make_session(obj, {"corner_values": flat}, n_regions)
    try:
        cut = session.describe_cut()
    except ddg.DdgError as ex:
        raise _translate(ex) from ex
    described = perf_counter()
    if cut.n_paths == 0:
        raise CutFailed("The cut is empty (one region covers the whole mesh), so there is "
                        "nothing to apply. Run 'Cut: Start' again.")
    n_polygons_before = len(mesh.polygons)
    n_vertices_before = len(mesh.vertices)
    seam_before = sum(core._seam_flags(mesh))
    stats = insert_cut(mesh, cut, mark_seams=mark_seams)
    remove_state(mesh)
    remove_preview(obj)
    done = perf_counter()
    stats.update({
        "n_polygons_before": n_polygons_before,
        "n_vertices_before": n_vertices_before,
        "n_paths": cut.n_paths,
        "seams_before": seam_before,
        "n_clamped_crossings": cut.n_clamped_crossings,
        "warnings": list(session.warnings),
        "seconds_extract": seconds_extract,
        "seconds_native": described - begin - seconds_extract,
        "seconds_insert": done - described,
    })
    return stats
