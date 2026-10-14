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

import base64
import bisect
import math
from array import array
from time import perf_counter

import bmesh
import bpy
from mathutils import Vector
from mathutils.kdtree import KDTree

from . import core

# The smallest distance of a point of the cut from a vertex, as a fraction of its edge, in the
# cut that is shown (preview) and inserted (Apply): the library's own default for the flow. It
# used to be 0.2 as a mitigation of flipped UV triangles; those come from Blender triangulating
# the split polygons again, and Apply now stores the library's triangulation instead (see
# `insert_cut`), after which no offset flips anything (blender/README.md).
DEFAULT_CUT_OFFSET = 1e-3

# The state of a cut. The per-corner attributes (below) are the whole state when the flow's
# triangulation has no vertex that is in no polygon. Otherwise the whole state is the bytes of
# `CutSession.snapshot()` in a mesh property (base64, as the property is a string), and the
# attributes are the view of the mesh's own vertices. The option says which triangulation the
# flow runs on: "intrinsic" (Delaunay flips), "refined:<degrees>" or "caller".
STATE_KEY = "ddg_cut_state"
OPTION_KEY = "ddg_cut_option"
DEFAULT_OPTION = "intrinsic"
DEFAULT_REFINE_BOUND = 25.0
# Bends of the cut that deviate from the straight line between their neighbours by less than
# this fraction of the longest edge of their triangle are removed before the cut is inserted
# (intrinsic triangulations; the library's default).
DEFAULT_MERGE_TOLERANCE = 0.1

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


def option_of(triangulation, refine_bound=DEFAULT_REFINE_BOUND):
    """The option string for the operator's enum value ('INTRINSIC', 'REFINED', 'CALLER')."""
    if triangulation == 'REFINED':
        return "refined:{:g}".format(refine_bound)
    return {'INTRINSIC': "intrinsic", 'CALLER': "caller"}[triangulation]


def option_keywords(option, offset=None, merge=None):
    """The keywords of ``CutSession`` for an option string; ``offset`` is the smallest distance
    of a point of the cut from a vertex of the caller's triangles (intrinsic sessions)."""
    kind, _, bound = option.partition(":")
    keywords = {"triangulation": kind}
    if kind == "refined":
        keywords["refine_bound"] = float(bound or DEFAULT_REFINE_BOUND)
    if kind != "caller" and offset is not None:
        keywords["min_offset"] = offset
    if kind != "caller" and merge is not None:
        keywords["merge_tolerance"] = merge
    return keywords


def write_state(mesh, session, option=None):
    """Writes the state of ``session`` to the region attributes, creating them as needed, and
    removes attributes of regions the session does not have. If the attributes are not the whole
    state (the triangulation added vertices) the state bytes go to the mesh property too."""
    if option is not None:
        mesh[OPTION_KEY] = option
    if session.state_is_per_corner:
        if STATE_KEY in mesh:
            del mesh[STATE_KEY]
    else:
        mesh[STATE_KEY] = base64.b64encode(session.snapshot()).decode("ascii")
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
    for key in (LAST_FLOW_KEY, STATE_KEY, OPTION_KEY):
        if key in mesh:
            del mesh[key]


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


def make_session(obj, start, n_regions=2, option=DEFAULT_OPTION, offset=None, merge=None):
    """A ``ddg_tools.CutSession`` on the mesh of ``obj`` (existing seams are fixed cuts).
    ``start`` is a dict of the start keywords of the constructor, ``option`` the triangulation
    (``option_keywords``). Returns ``(session, seconds)`` with the time spent extracting the
    buffers."""
    ddg = native()
    mesh = check_mesh(obj)
    begin = perf_counter()
    positions, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
    # Polygons that an earlier cut split keep the triangulation the library chose for them:
    # Blender's own can be refused (collinear vertices) and is worse shaped.
    triangles, _ = core.stored_triangulation(mesh, corner_vertices, polygon_starts, triangles)
    seams = core.seam_edges(mesh)
    extracted = perf_counter()
    keywords = dict(start)
    if "state" not in keywords:
        keywords["n_regions"] = n_regions
        keywords.update(option_keywords(option, offset, merge))
    elif not option.startswith("caller"):
        if offset is not None:
            keywords["min_offset"] = offset
        if merge is not None:
            keywords["merge_tolerance"] = merge
    try:
        session = ddg.CutSession(
            positions, corner_vertices, polygon_starts, triangles,
            seam_edges=seams if len(seams) else None, **keywords)
    except ddg.StaleStateError as ex:
        raise CutFailed(
            "The mesh or its seams changed since 'Cut: Start', so the stored cut no longer fits "
            "it ({}). Run 'Cut: Start' again.".format(ex)) from ex
    except ddg.DdgError as ex:
        raise _translate(ex) from ex
    return session, extracted - begin


def open_session(obj, offset=None, merge=None):
    """The session of the cut stored on the mesh of ``obj``: from the state bytes if the mesh has
    them, else from the per-corner attributes with the triangulation option stored at Start.
    Returns ``(session, seconds, option)``."""
    mesh = check_mesh(obj)
    option = mesh.get(OPTION_KEY, DEFAULT_OPTION)
    if STATE_KEY in mesh:
        try:
            state = base64.b64decode(mesh[STATE_KEY])
        except ValueError as ex:
            raise CutFailed("The stored cut state of mesh '{}' is damaged. Run 'Cut: Start' "
                            "again.".format(mesh.name)) from ex
        session, seconds = make_session(obj, {"state": state}, option=option, offset=offset,
                                         merge=merge)
        return session, seconds, option
    n_regions, flat = read_state(mesh)
    session, seconds = make_session(obj, {"corner_values": flat}, n_regions, option, offset, merge)
    return session, seconds, option


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


MAX_REGIONS = 16
"""The most regions the operators handle (the weights of the Flow operator are fixed-size
vectors)."""


def labels_from_values(values, what):
    """Region labels from one integer per face: the distinct values, in ascending order, become
    the regions 0, 1, ... Returns ``(array('i'), n_regions)``."""
    distinct = sorted(set(values))
    if len(distinct) < 2:
        raise CutFailed("{} has only one value, but a cut needs at least two regions".format(what))
    if len(distinct) > MAX_REGIONS:
        raise CutFailed("{} has {} distinct values, but at most {} regions are handled".format(
            what, len(distinct), MAX_REGIONS))
    index = {v: i for i, v in enumerate(distinct)}
    return array('i', (index[v] for v in values)), len(distinct)


def start_cut(obj, *, mode, radius=0.0, vertex=-1, n_regions=2, active_vertex=None,
              selected_vertices=(), cluster_iterations=50, attribute="",
              cut_offset=DEFAULT_CUT_OFFSET, merge_tolerance=DEFAULT_MERGE_TOLERANCE,
              triangulation='INTRINSIC',
              refine_bound=DEFAULT_REFINE_BOUND):
    """Creates the state of a cut on ``obj`` and its preview.

    ``mode`` is ``'SELECTED'`` (the selected faces against the rest, two regions), ``'CIRCLE'``
    (the geodesic disk of ``radius`` object units, 0: a tenth of the bounding diagonal, around
    ``vertex``; with ``vertex < 0`` the ``active_vertex``, or the only one of
    ``selected_vertices``; two regions), ``'CLUSTERS'`` (the vertex normals clustered into
    ``n_regions`` regions, at most ``cluster_iterations`` rounds), ``'MATERIALS'`` (the
    material slot of every face is its region) or ``'ATTRIBUTE'`` (the integer face attribute
    ``attribute``; its distinct values are the regions).

    Returns a dict with the session's sizes and timings.
    """
    mesh = check_mesh(obj)
    if mode == 'SELECTED':
        start = {"polygon_labels": start_labels_from_selection(mesh)}
        n_regions = 2
    elif mode == 'CLUSTERS':
        if not 2 <= n_regions <= MAX_REGIONS:
            raise CutFailed("The number of regions must be between 2 and {} (it is {})".format(
                MAX_REGIONS, n_regions))
        start = {"cluster_iterations": cluster_iterations}
    elif mode == 'MATERIALS':
        values = array('i', bytes(4 * len(mesh.polygons)))
        mesh.polygons.foreach_get("material_index", values)
        labels, n_regions = labels_from_values(values, "The material index of the faces")
        start = {"polygon_labels": labels}
    elif mode == 'ATTRIBUTE':
        attr = mesh.attributes.get(attribute) if attribute else None
        if attr is None:
            raise CutFailed("There is no attribute named '{}' on mesh '{}'".format(
                attribute, mesh.name))
        if attr.domain != 'FACE' or attr.data_type not in ('INT', 'INT8'):
            raise CutFailed(
                "Attribute '{}' must be an integer attribute on the face domain, but it is {} "
                "on {}".format(attribute, attr.data_type, attr.domain))
        values = array('i', bytes(4 * len(mesh.polygons)))
        attr.data.foreach_get("value", values)
        labels, n_regions = labels_from_values(values, "Attribute '{}'".format(attribute))
        start = {"polygon_labels": labels}
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
        n_regions = 2
    else:
        raise CutFailed("Unknown start mode {!r}".format(mode))
    option = option_of(triangulation, refine_bound)
    session, seconds_extract = make_session(obj, start, n_regions, option, cut_offset, merge_tolerance)
    begin = perf_counter()
    write_state(mesh, session, option)
    if LAST_FLOW_KEY in mesh:
        del mesh[LAST_FLOW_KEY]
    mesh.update()
    written = perf_counter()
    update_preview(obj, session, cut_offset)
    previewed = perf_counter()
    return {
        "n_regions": session.n_regions,
        "triangulation": option,
        "intrinsic": session.intrinsic_report,
        "state_is_per_corner": session.state_is_per_corner,
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


def update_preview(obj, session, offset=DEFAULT_CUT_OFFSET):
    """Replaces the preview of ``obj`` by the current cut of ``session``: a wire object named
    after the mesh object, parented to it (so that it follows its transform). Returns the new
    object."""
    ddg = native()
    try:
        positions, starts, closed, _regions = session.preview(min_crossing_offset=offset)
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
    controls["stages"] = settings["stages"].lower()
    if energy == "scaled_dirichlet":
        controls["feature_r_eps"] = settings["feature_r_eps"]
        controls["feature_smoothing"] = settings["feature_smoothing"]
    if settings["use_uniform_curvature"]:
        controls["curvature_uniform"] = settings["uniform_curvature"]
    n = session.n_regions
    controls["distortion_weights"] = tuple(settings["distortion_weights"])[:n]
    controls["area_weights"] = tuple(settings["area_weights"])[:n]

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
    session, seconds_extract, option = open_session(obj, settings["cut_offset"], settings["merge_tolerance"])
    n_regions = session.n_regions
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
    write_state(mesh, session, option)
    mesh.update()
    written = perf_counter()
    n_patches = count_patches(mesh, session.corner_regions())
    n_triple_points = 0
    merge = {}
    if n_regions > 2 or settings["show_advanced"]:
        try:
            described = session.describe_cut(min_crossing_offset=settings["cut_offset"])
        except ddg.DdgError as ex:
            raise _translate(ex) from ex
        n_triple_points = len(described.junction_polygon) if n_regions > 2 else 0
        merge = dict(described.report)
    update_preview(obj, session, settings["cut_offset"])
    previewed = perf_counter()
    summary = {k: v for k, v in report.items() if k != "steps"}
    summary["length_scale"] = session.length_scale
    summary["n_regions"] = n_regions
    summary["n_patches"] = n_patches
    summary["n_triple_points"] = n_triple_points
    report["n_patches"] = n_patches
    report["n_triple_points"] = n_triple_points
    report["n_regions"] = n_regions
    report["merge"] = merge
    intrinsic = session.intrinsic_report
    report["triangulation"] = option
    report["intrinsic"] = intrinsic
    summary["triangulation"] = option
    if intrinsic is not None:
        summary["vertices_added"] = intrinsic["n_added"]
        summary["min_angle_degrees"] = intrinsic["min_angle_degrees"]
        summary["max_angle_degrees"] = intrinsic["max_angle_degrees"]
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


def count_patches(mesh, regions):
    """The number of patches: connected pieces of the surface that are in one region. Two
    corners of a polygon that are next to each other and in one region are joined, and so are the
    two corners of one vertex on the two sides of an edge that is not a seam."""
    n = len(mesh.loops)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def join(a, b):
        if regions[a] == regions[b]:
            a, b = find(a), find(b)
            if a != b:
                parent[max(a, b)] = min(a, b)

    following = [0] * n
    by_edge = {}
    loops = mesh.loops
    for polygon in mesh.polygons:
        first = polygon.loop_start
        count = polygon.loop_total
        for k in range(count):
            c = first + k
            following[c] = first + (k + 1) % count
            join(c, following[c])
            by_edge.setdefault(loops[c].edge_index, []).append(c)
    seam = {e.index for e in mesh.edges if e.use_seam}
    for edge, cs in by_edge.items():
        if edge in seam or len(cs) != 2:
            continue
        a, b = cs
        if loops[a].vertex_index == loops[b].vertex_index:
            join(a, b)
            join(following[a], following[b])
        else:
            join(a, following[b])
            join(following[a], b)
    return len({find(c) for c in range(n)})


def triangulation_line(report):
    """One short line about the triangulation of a run (the advanced view)."""
    intrinsic = report.get("intrinsic")
    if intrinsic is None:
        return "triangulation: the mesh's own"
    line = "triangulation: {} | {} vertices added, angles {:.1f} to {:.1f} degrees".format(
        report["triangulation"], intrinsic["n_added"], intrinsic["min_angle_degrees"],
        intrinsic["max_angle_degrees"])
    merge = report.get("merge")
    if merge:
        line += " | {} bends merged, cut moved by at most {:.3g} of an edge".format(
            merge["n_merged_nodes"], merge["max_displacement_relative"])
    return line


def flow_summary(report, advanced=False):
    """One line for Blender's status bar."""
    reason = {"max_steps": "step limit", "converged": "converged",
              "energy_stalled": "energy stalled", "cut_vanished": "the cut vanished",
              "energy_rose": "energy rose, back at the lowest-energy state"}[
        report["stop_reason"]]
    network = (" | {} regions, {} patches, {} triple points".format(
        report["n_regions"], report["n_patches"], report["n_triple_points"])
        if report["n_regions"] > 2 else " | {} patches".format(report["n_patches"]))
    extra = " | " + triangulation_line(report) if advanced else ""
    return ("Cut flow: {} step{} ({}) | energy {:.5g} -> {:.5g} | cut length {:.5g} "
            "(was {:.5g}) | length weight {:.4g}" + network + extra).format(
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


def _interpolate_vertex(bm, new, sources, layers):
    """Vertex data of ``new``, a point that is the weighted mean of vertices: ``sources`` is a
    list of ``(vertex, weight)``, the weights add up to 1. Vertex group weights and the vertex
    layers are blended (layers that cannot be blended are taken from the heaviest source)."""
    deform = bm.verts.layers.deform.active
    if deform is not None:
        weights = {}
        for vert, w in sources:
            for group, value in vert[deform].items():
                weights[group] = weights.get(group, 0.0) + w * value
        for group, value in weights.items():
            if value > 0.0:
                new[deform][group] = value
    heaviest = max(sources, key=lambda source: source[1])[0]
    for layer, interpolate in layers:
        if not interpolate:
            new[layer] = heaviest[layer]
            continue
        values = [(vert[layer], w) for vert, w in sources]
        if isinstance(values[0][0], (int, float, bool)):
            new[layer] = sum(v * w for v, w in values)
        else:
            total = Vector(values[0][0]) * values[0][1]
            for v, w in values[1:]:
                total = total + Vector(v) * w
            new[layer] = total


def _blend_many(values, weights, interpolate):
    """The weighted mean of loop values (``_loop_value``), or the heaviest one."""
    if not interpolate:
        return values[max(range(len(values)), key=lambda i: weights[i])]
    total = values[0] * weights[0]
    for v, w in zip(values[1:], weights[1:]):
        total = total + v * w
    return total


def insert_cut(mesh, cut, *, mark_seams=True, reference=None):
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

    ``reference`` is the Rust insertion of the same cut (``CutSession.insert``). If it is given,
    the triangulation it chose for the polygons that were split is stored in the mesh property
    ``ddg_cut_triangulation`` (``core.TRIANGULATION_KEY``) for the unwrap: the surface the flow
    saw is made of the triangles of the original triangulation cut along the path, and this is
    the triangulation of the pieces; triangulating the new polygons again, as Blender does,
    gives slivers of the cut near vertices that the flattening can flip.

    Returns a dict with ``n_vertices``, ``n_polygons``, ``n_edge_points``, ``n_seams_marked``,
    ``n_triangulated`` (polygons with a stored triangulation).
    """
    if mesh.shape_keys is not None:
        raise CutFailed(
            "Mesh '{}' has shape keys; applying a cut would change its vertices. Apply or remove "
            "the shape keys first.".format(mesh.name))
    bm = bmesh.new()
    try:
        bm.from_mesh(mesh)
        old_entries = core.read_triangulation(mesh)
        stats, entries = _insert_cut(bm, cut, mark_seams, reference)
        if old_entries is not None:
            # Polygons the cut did not touch keep the triangulation they had.
            for key, triangles in old_entries[2].items():
                entries.setdefault(key, triangles)
        bm.to_mesh(mesh)
        if reference is not None or old_entries is not None:
            keys = {tuple(sorted(v.index for v in f.verts)) for f in bm.faces}
            core.write_triangulation(mesh, stats["n_vertices"], stats["n_polygons"],
                                     {k: t for k, t in entries.items() if k in keys})
            stats["n_triangulated"] = sum(1 for k in entries if k in keys)
    finally:
        bm.free()
    mesh.update()
    return stats


def _insert_cut(bm, cut, mark_seams, reference=None):
    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.faces.ensure_lookup_table()
    n_points = cut.n_edge_points
    n_paths = cut.n_paths
    face_starts = []
    old_sizes = []
    corner = 0
    for face in bm.faces:
        face_starts.append(corner)
        old_sizes.append(len(face.loops))
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
    edge_chain = {}
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
        chain_list = [va]
        for group in groups:
            t = edge_t[group[0]]
            fraction = (t - t_before) / (1.0 - t_before)
            _near, new_vert = bmesh.utils.edge_split(current, from_vert, fraction)
            for i in group:
                point_vertex[i] = new_vert
            chain_list.append(new_vert)
            # `current` now runs from the new vertex to the far end; `new_edge` is the near part.
            from_vert, t_before = new_vert, t
        edge_chain[(a, b)] = chain_list + [vb]

    # ----- paths: the network inside every polygon -----
    # A polygon's part of the cut is a forest whose leaves are edge points (vertices already) and
    # whose inner nodes are junctions (three paths each). It is embedded by splitting the polygon
    # along chains of paths that run from a node that is a vertex already to another one:
    # at first only edge points are, then also the junctions the chains passed through.
    path_start, path_end = cut.path_start, cut.path_end
    interior_starts = cut.path_interior_starts
    interior_corners, interior_t = cut.interior_corners, cut.interior_t
    junction_corners, junction_weights = cut.junction_corners, cut.junction_weights
    pieces = {p: [old_faces[p]] for p in polygons_with_paths}
    cut_segments = []
    junction_vertex = {}
    interior_vertex = {}
    paths_of_polygon = {}
    for k in range(n_paths):
        paths_of_polygon.setdefault(path_polygon[k], []).append(k)

    def node_of(code):
        return ("j", -code - 1) if code < 0 else ("e", code)

    def vertex_of(node):
        return junction_vertex.get(node[1]) if node[0] == "j" else point_vertex[node[1]]

    for p, path_ids in paths_of_polygon.items():
        corners = old_corners[p]
        pending = set(path_ids)
        ends = {k: (node_of(path_start[k]), node_of(path_end[k])) for k in path_ids}
        incident = {}
        for k, (u, v) in ends.items():
            incident.setdefault(u, []).append(k)
            incident.setdefault(v, []).append(k)

        def junction_position(j):
            c = [junction_corners[3 * j + i] - face_starts[p] for i in range(3)]
            w = [junction_weights[3 * j + i] for i in range(3)]
            return (corners[c[0]][0].co * w[0] + corners[c[1]][0].co * w[1]
                    + corners[c[2]][0].co * w[2]), c, w

        # A piece along a polygon edge: both ends are points of the same edge, nothing to split;
        # the sides of the chain of the edge between them are seams.
        for k in sorted(pending):
            if path_start[k] >= 0 and path_end[k] >= 0 and interior_starts[k] == interior_starts[k + 1]:
                a_, b_ = path_start[k], path_end[k]
                if tuple(edge_vertices[2 * a_:2 * a_ + 2]) == tuple(edge_vertices[2 * b_:2 * b_ + 2]):
                    chain_ = edge_chain[tuple(edge_vertices[2 * a_:2 * a_ + 2])]
                    i_, j_ = chain_.index(point_vertex[a_]), chain_.index(point_vertex[b_])
                    if i_ > j_:
                        i_, j_ = j_, i_
                    cut_segments.extend(zip(chain_[i_:j_], chain_[i_ + 1:j_ + 1]))
                    pending.discard(k)

        while pending:
            start_node = next(
                (n for n, ks in incident.items()
                 if vertex_of(n) is not None and any(k in pending for k in ks)), None)
            if start_node is None:
                raise CutFailed("Internal error: the paths of polygon {} do not reach an edge "
                                "point".format(p))
            # The shortest chain of paths from the start node, through junctions that are not
            # vertices yet, to another node that is a vertex.
            chain, used = _chain(start_node, incident, ends, pending, vertex_of)
            if chain is None:
                raise CutFailed("Internal error: the paths of polygon {} form a cycle that no "
                                "chain of splits can open".format(p))
            node = ends[chain[-1][0]][1] if chain[-1][1] else ends[chain[-1][0]][0]
            coords, labels = [], []        # labels: ("i", c0, c1, t) or ("j", j, corners, weights)
            walk_nodes = [start_node]
            for k, forward in chain:
                items = []
                for j in range(interior_starts[k], interior_starts[k + 1]):
                    c0 = interior_corners[2 * j] - face_starts[p]
                    c1 = interior_corners[2 * j + 1] - face_starts[p]
                    t = interior_t[j]
                    items.append((corners[c0][0].co * (1.0 - t) + corners[c1][0].co * t,
                                  ("i", c0, c1, t, k, j - interior_starts[k])))
                if not forward:
                    items.reverse()
                for co, label in items:
                    coords.append(co)
                    labels.append(label)
                next_node = ends[k][1] if forward else ends[k][0]
                walk_nodes.append(next_node)
                if vertex_of(next_node) is None:
                    position, c, w = junction_position(next_node[1])
                    coords.append(position)
                    labels.append(("j", next_node[1], c, w))
            first, last = vertex_of(start_node), vertex_of(walk_nodes[-1])
            candidates = [f for f in pieces[p] if first in f.verts and last in f.verts]
            if len(candidates) > 1 and coords:
                # Both ends are on the boundary of several pieces: the one the chain lies in.
                inside = [f for f in candidates
                          if bmesh.geometry.intersect_face_point(f, coords[0])]
                candidates = inside or candidates
            face = candidates[0] if candidates else None
            if face is None:
                raise CutFailed("Internal error: a chain of polygon {} does not lie in one "
                                "piece of it".format(p))
            before = set(face.verts)
            result = bmesh.utils.face_split(face, first, last, coords=coords)
            if not result or result[0] is None:
                raise CutFailed("Internal error: polygon {} could not be split along a chain of "
                                "paths".format(p))
            new_face = result[0]
            pieces[p].append(new_face)
            made = [v for v in new_face.verts if v not in before]
            order = _ordered(first, made) if coords else []
            if order is None or len(order) != len(coords):
                # Fall back to the nearest coordinates (points that coincide may swap).
                order = []
                for co in coords:
                    vert = min(made, key=lambda v: (v.co - co).length_squared)
                    made.remove(vert)
                    order.append(vert)
            vertices = [first]
            for vert, (co, label) in zip(order, zip(coords, labels)):
                vertices.append(vert)
                if label[0] == "i":
                    _, c0, c1, t, path_k, path_j = label
                    interior_vertex[(path_k, path_j)] = vert
                    sources = [(corners[c0][0], 1.0 - t), (corners[c1][0], t)]
                    sample = [corners[c0][1], corners[c1][1]]
                    weights = [1.0 - t, t]
                else:
                    _, j, c, w = label
                    junction_vertex[j] = vert
                    sources = [(corners[c[i]][0], w[i]) for i in range(3)]
                    sample = [corners[c[i]][1] for i in range(3)]
                    weights = list(w)
                _interpolate_vertex(bm, vert, sources, vertex_layers)
                for loop in vert.link_loops:
                    for index, (layer, interpolate) in enumerate(loop_layers):
                        _set_loop_value(loop, layer, _blend_many(
                            [values[index] for values in sample], weights, interpolate))
            vertices.append(last)
            cut_segments.extend(zip(vertices, vertices[1:]))
            pending -= used

    n_marked = 0
    if mark_seams:
        for u, v in cut_segments:
            edge = bm.edges.get((u, v))
            if edge is None:
                raise CutFailed("Internal error: a segment of the cut is not an edge")
            if not edge.seam:
                edge.seam = True
                n_marked += 1
    stats = {
        "n_vertices": len(bm.verts),
        "n_polygons": len(bm.faces),
        "n_edge_points": n_points,
        "n_seams_marked": n_marked,
        "n_triangulated": 0,
    }
    entries = {}
    if reference is not None:
        entries = _reference_triangulation(bm, reference, point_vertex, junction_vertex,
                                           interior_vertex, old_sizes)
    return stats, entries


def _chain(start, incident, ends, pending, vertex_of):
    """The shortest chain of pending paths from the node ``start`` (a vertex already) to another
    node that is a vertex, through nodes that are not (breadth first). Returns ``(chain, used)``
    with ``chain`` a list of ``(path, forward)``, or ``(None, None)``."""
    queue = [start]
    came = {start: None}
    while queue:
        node = queue.pop(0)
        for k in incident[node]:
            if k not in pending:
                continue
            u, v = ends[k]
            other = v if u == node else u
            if other in came or other == start:
                continue
            came[other] = (node, k, u == node)
            if vertex_of(other) is not None:
                chain = []
                at = other
                while came[at] is not None:
                    previous, path, forward = came[at]
                    chain.append((path, forward))
                    at = previous
                chain.reverse()
                return chain, {path for path, _ in chain}
            queue.append(other)
    return None, None


def _ordered(first, made):
    """The vertices ``made`` (new vertices of a chain), in order along it, starting from the
    neighbour of ``first``."""
    remaining = set(made)
    order = []
    current = first
    while remaining:
        following = [e.other_vert(current) for e in current.link_edges
                     if e.other_vert(current) in remaining]
        if len(following) != 1:
            return None
        current = following[0]
        order.append(current)
        remaining.discard(current)
    return order


def _reference_triangulation(bm, reference, point_vertex, junction_vertex, interior_vertex,
                             old_sizes):
    """The triangles of the polygons of the Rust insertion that were split or got extra vertices,
    as ``{ascending vertex tuple: flat vertex triples}`` in the numbering of ``bm``."""
    bm.verts.index_update()
    source = reference.vertex_source
    vertex_of = []
    for v in range(len(source) // 3):
        kind, a, b = source[3 * v:3 * v + 3]
        if kind == 0:
            vertex_of.append(a)
        elif kind == 1:
            vertex_of.append(point_vertex[a].index)
        elif kind == 2:
            vertex_of.append(interior_vertex[(a, b)].index)
        else:
            vertex_of.append(junction_vertex[a].index)
    starts = list(reference.polygon_starts)
    corners = reference.corner_vertices
    parents = reference.polygon_parent
    children = {}
    for parent in parents:
        children[parent] = children.get(parent, 0) + 1
    triangles = reference.triangles
    n_polygons = len(starts)
    by_polygon = {}
    for t in range(len(triangles) // 3):
        corners3 = triangles[3 * t:3 * t + 3]
        by_polygon.setdefault(bisect.bisect_right(starts, corners3[0]) - 1, []).extend(corners3)
    entries = {}
    for p in range(n_polygons):
        end = starts[p + 1] if p + 1 < n_polygons else len(corners)
        size = end - starts[p]
        if size <= 3:
            continue
        if children[parents[p]] == 1 and size == old_sizes[parents[p]]:
            continue
        vertices = [vertex_of[corners[c]] for c in range(starts[p], end)]
        flat = [vertex_of[corners[c]] for c in by_polygon[p]]
        entries[tuple(sorted(vertices))] = flat
    return entries


def apply_cut(obj, *, mark_seams=True, cut_offset=DEFAULT_CUT_OFFSET,
              merge_tolerance=DEFAULT_MERGE_TOLERANCE):
    """Applies the cut stored in the state attributes of ``obj``: inserts it into the mesh
    (``insert_cut``), marks it as seams, and removes the preview and the state attributes.
    Returns a dict with the insertion's counts and the time spent."""
    ddg = native()
    mesh = check_mesh(obj)
    begin = perf_counter()
    session, seconds_extract, option = open_session(obj, cut_offset, merge_tolerance)
    try:
        cut = session.describe_cut(min_crossing_offset=cut_offset)
    except ddg.DdgError as ex:
        raise _translate(ex) from ex
    described = perf_counter()
    if cut.n_paths == 0:
        raise CutFailed("The cut is empty (one region covers the whole mesh), so there is "
                        "nothing to apply. Run 'Cut: Start' again.")
    n_polygons_before = len(mesh.polygons)
    n_vertices_before = len(mesh.vertices)
    seam_before = sum(core._seam_flags(mesh))
    stats = insert_cut(mesh, cut, mark_seams=mark_seams, reference=session.insert(cut))
    remove_state(mesh)
    remove_preview(obj)
    done = perf_counter()
    stats.update({
        "n_polygons_before": n_polygons_before,
        "n_vertices_before": n_vertices_before,
        "n_paths": cut.n_paths,
        "seams_before": seam_before,
        "n_triangulated": stats.get("n_triangulated", 0),
        "n_clamped_crossings": cut.n_clamped_crossings,
        "merge": dict(cut.report),
        "warnings": list(session.warnings),
        "seconds_extract": seconds_extract,
        "seconds_native": described - begin - seconds_extract,
        "seconds_insert": done - described,
    })
    return stats
