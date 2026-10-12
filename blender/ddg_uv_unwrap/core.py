# SPDX-License-Identifier: MIT
"""Mesh to buffers, buffers to UVs.

Everything here goes through Blender's data API (``foreach_get`` / ``foreach_set``) and the
standard library; there are no operator calls and no numpy. The mesh must not be in Edit mode
(the operator leaves Edit mode before calling in and re-enters it afterwards).

The native module ``ddg_tools`` (crates/ddg-py of the ddg-tools repository) does the work. Its
array arguments are exactly the buffers extracted here.
"""

from array import array
from time import perf_counter

# Free space around every chart in the packed layout, as a fraction of the UV square.
PACK_MARGIN = 1.0 / 256.0

# Name given to the UV map when the mesh has none.
DEFAULT_UV_NAME = "UVMap"


class UnwrapFailed(Exception):
    """The unwrap could not be done; the message is meant for the user."""


def native():
    """The native module, imported on first use so that a broken installation is reported by
    the operator instead of preventing the add-on from loading."""
    try:
        import ddg_tools
    except ImportError as ex:
        raise UnwrapFailed(
            "The native module 'ddg_tools' could not be loaded ({}). It is bundled with the "
            "extension as a Python wheel; reinstall the extension for this platform.".format(ex)
        ) from ex
    return ddg_tools


def is_placeholder():
    """True when the bundled native module contains the stand-in pipeline, not the real unwrap.
    False when it cannot be loaded (the operator reports that when it runs)."""
    try:
        return bool(native().IS_PLACEHOLDER)
    except UnwrapFailed:
        return False


def _zeros(typecode, count):
    """A zero-filled ``array.array`` of ``count`` items of 4 bytes each ('f', 'i' or 'I')."""
    return array(typecode, bytes(4 * count))


# How mesh data is read, and why it matters: `foreach_get` copies memory only when the buffer's
# element type is exactly the property's (float32, signed or unsigned int32, bool) *and* Blender
# can reach the property as a plain array. Otherwise it converts item by item, which is 100 to
# 1000 times slower (measured in blender/README.md). So:
#
# * face starts, loop triangles and edge vertices are unsigned properties: `array('I')`;
# * positions, corner vertices and seam flags are stored as the attributes "position",
#   ".corner_vert" and "uv_seam", whose `data` is a plain array; the element collections
#   (`mesh.vertices`, `mesh.loops`, `mesh.edges`) reach the same values through functions.
#
# The attribute names are checked once per session against the documented element API on a
# two-triangle mesh (`_attribute_paths`); a name that does not check out is not used, and the
# element API, which is always correct, takes over.

_ATTRIBUTE_PATHS = None


def _attribute_paths():
    """Which attribute shortcuts are valid in this Blender: a dict with the keys "position",
    "corner_vert" and "seam", each an attribute name or None."""
    global _ATTRIBUTE_PATHS
    if _ATTRIBUTE_PATHS is not None:
        return _ATTRIBUTE_PATHS
    import bpy

    paths = {"position": None, "corner_vert": None, "seam": None}
    probe = bpy.data.meshes.new(".ddg_uv_unwrap_probe")
    try:
        probe.from_pydata(
            [(0.5, 0.25, 0.125), (1.5, 0.0, 2.0), (1.0, 3.0, -1.0), (-2.0, 1.0, 4.0)],
            [], [(0, 1, 2), (2, 1, 3)])
        probe.update()
        probe.edges[1].use_seam = True
        expected = _zeros('f', 12)
        probe.vertices.foreach_get("co", expected)
        got = _zeros('f', 12)
        if _attribute_get(probe, "position", 'POINT', 'FLOAT_VECTOR', "vector", got, 3) \
                and got == expected:
            paths["position"] = "position"
        expected = _zeros('I', 6)
        probe.loops.foreach_get("vertex_index", expected)
        got = _zeros('i', 6)
        if _attribute_get(probe, ".corner_vert", 'CORNER', 'INT', "value", got, 1) \
                and list(got) == list(expected):
            paths["corner_vert"] = ".corner_vert"
        expected = bytearray(len(probe.edges))
        probe.edges.foreach_get("use_seam", memoryview(expected).cast('?'))
        got = bytearray(len(probe.edges))
        if _attribute_get(probe, "uv_seam", 'EDGE', 'BOOLEAN', "value",
                          memoryview(got).cast('?'), 1) and got == expected and sum(got) == 1:
            paths["seam"] = "uv_seam"
    except Exception:
        # Any surprise here means: do not use shortcuts.
        paths = {"position": None, "corner_vert": None, "seam": None}
    finally:
        bpy.data.meshes.remove(probe)
    _ATTRIBUTE_PATHS = paths
    return paths


def _attribute_get(mesh, name, domain, data_type, field, values, width):
    """Fills ``values`` from the attribute ``name`` if the mesh has it with the expected domain,
    type and size. Returns whether it did."""
    attribute = mesh.attributes.get(name)
    if attribute is None or attribute.domain != domain or attribute.data_type != data_type:
        return False
    if width * len(attribute.data) != len(values):
        return False
    attribute.data.foreach_get(field, values)
    return True


def mesh_buffers(mesh):
    """The arrays that describe ``mesh`` to the native module.

    Returns ``(positions, corner_vertices, polygon_starts, triangles)``:

    * ``positions``: ``array('f')``, 3 per vertex (object space);
    * ``corner_vertices``: ``array('i')`` or ``array('I')``, the vertex of every face corner
      ("loop");
    * ``polygon_starts``: ``array('I')``, the first corner of every face;
    * ``triangles``: ``array('I')``, 3 corner indices per triangle of Blender's own
      tessellation (``mesh.loop_triangles``), so the native side triangulates n-gons exactly
      like Blender does.
    """
    paths = _attribute_paths()
    positions = _zeros('f', 3 * len(mesh.vertices))
    if not (paths["position"] and _attribute_get(
            mesh, paths["position"], 'POINT', 'FLOAT_VECTOR', "vector", positions, 3)):
        mesh.vertices.foreach_get("co", positions)
    corner_vertices = _zeros('i', len(mesh.loops))
    if not (paths["corner_vert"] and _attribute_get(
            mesh, paths["corner_vert"], 'CORNER', 'INT', "value", corner_vertices, 1)):
        corner_vertices = _zeros('I', len(mesh.loops))
        mesh.loops.foreach_get("vertex_index", corner_vertices)
    polygon_starts = _zeros('I', len(mesh.polygons))
    mesh.polygons.foreach_get("loop_start", polygon_starts)
    mesh.calc_loop_triangles()
    triangles = _zeros('I', 3 * len(mesh.loop_triangles))
    mesh.loop_triangles.foreach_get("loops", triangles)
    return positions, corner_vertices, polygon_starts, triangles


def _edge_vertices(mesh):
    """``array('I')`` with the 2 vertices of every edge."""
    edge_vertices = _zeros('I', 2 * len(mesh.edges))
    mesh.edges.foreach_get("vertices", edge_vertices)
    return edge_vertices


def _seam_flags(mesh):
    """One byte (0 or 1) per edge: whether the edge is marked as a seam."""
    flags = bytearray(len(mesh.edges))
    # `array` has no bool typecode; a bool view of a bytearray is what Blender copies into.
    view = memoryview(flags).cast('?')
    name = _attribute_paths()["seam"]
    if name is None:
        mesh.edges.foreach_get("use_seam", view)
    elif mesh.attributes.get(name) is not None:
        if not _attribute_get(mesh, name, 'EDGE', 'BOOLEAN', "value", view, 1):
            mesh.edges.foreach_get("use_seam", view)
    # else: the attribute only exists once an edge has been marked, so there are no seams.
    return flags


def _set_seam_flags(mesh, flags, first_new):
    """Writes one byte per edge to the seam flags. ``first_new`` is an edge that is being
    marked: marking it through the element API makes Blender create the seam attribute if the
    mesh has none yet."""
    view = memoryview(flags).cast('?')
    name = _attribute_paths()["seam"]
    if name is not None:
        mesh.edges[first_new].use_seam = True
        attribute = mesh.attributes.get(name)
        if attribute is not None and attribute.domain == 'EDGE' \
                and attribute.data_type == 'BOOLEAN' and len(attribute.data) == len(flags):
            attribute.data.foreach_set("value", view)
            return
    mesh.edges.foreach_set("use_seam", view)


def seam_edges(mesh):
    """The edges marked as seams, as ``array('I')`` with 2 vertex indices per edge."""
    flags = _seam_flags(mesh)
    out = array('I')
    i = flags.find(1)
    if i < 0:
        return out
    edge_vertices = _edge_vertices(mesh)
    while i >= 0:
        out.extend(edge_vertices[2 * i:2 * i + 2])
        i = flags.find(1, i + 1)
    return out


def mark_seams(mesh, cut_edges):
    """Marks the edges in ``cut_edges`` (2 vertex indices per edge) as seams. Existing seams are
    kept. Returns ``(newly_marked, not_found)``; ``not_found`` counts vertex pairs that are not
    an edge of the mesh."""
    wanted = set(zip(cut_edges[0::2], cut_edges[1::2]))
    wanted = {(min(a, b), max(a, b)) for a, b in wanted}
    if not wanted:
        return 0, 0
    both_ways = wanted | {(b, a) for a, b in wanted}
    edge_vertices = _edge_vertices(mesh)
    found = [edge for edge, pair in enumerate(zip(edge_vertices[0::2], edge_vertices[1::2]))
             if pair in both_ways]
    flags = _seam_flags(mesh)
    new = [edge for edge in found if not flags[edge]]
    if new:
        for edge in new:
            flags[edge] = 1
        _set_seam_flags(mesh, flags, new[0])
    return len(new), len(wanted) - len(found)


def groups_vertices(obj, group_names):
    """For each name in ``group_names``: the vertices of ``obj`` that are in that vertex group
    with a weight above zero, as ``array('i')`` in ascending order. One pass over the vertices
    for all groups together (there is no bulk access to vertex group weights)."""
    arrays = {}
    for name in group_names:
        group = obj.vertex_groups.get(name)
        if group is None:
            raise UnwrapFailed(
                "Vertex group '{}' does not exist on object '{}'".format(name, obj.name)
            )
        arrays.setdefault(group.index, array('i'))
    if arrays:
        lookup = arrays.get
        for vertex in obj.data.vertices:
            for element in vertex.groups:
                members = lookup(element.group)
                if members is not None and element.weight > 0.0:
                    members.append(vertex.index)
    return [arrays[obj.vertex_groups[name].index] for name in group_names]


def group_vertices(obj, group_name):
    """The vertices in one vertex group; see ``groups_vertices``."""
    return groups_vertices(obj, [group_name])[0]


def unwrap_object(obj, *, lambda_=0.5, use_seams=True, pinned_group="", forbidden_group="",
                  positive_cones_only=False, mark_cuts_as_seams=False):
    """Generates a UV map for the mesh of ``obj`` and writes it to the active UV map, creating
    one if the mesh has none. Vertices, edges and faces are not touched.

    Returns a dict:

    * ``report``: the report dict of ``ddg_tools.unwrap`` (counts, distortion, warnings);
    * ``pack``: the info dict of ``ddg_tools.pack_charts``;
    * ``uv_layer``: name of the UV map written, ``created_uv_layer``: whether it is new;
    * ``n_cones``, ``n_cut_edges``, ``seams_marked``, ``cut_edges_not_found``;
    * ``seconds_extract``, ``seconds_native``, ``seconds_pack``, ``seconds_write``,
      ``seconds_seams``: where the time went. ``seconds_native`` (``ddg_tools.unwrap``) and
      ``seconds_pack`` (``ddg_tools.pack_charts``) are spent in Rust; the other three are the
      cost of the Blender side: reading the mesh, writing the UVs, marking seams.

    Raises ``UnwrapFailed`` with a message for the user when the mesh cannot be unwrapped.
    """
    ddg = native()
    mesh = obj.data
    if mesh.is_editmode:
        raise UnwrapFailed("Internal error: the mesh is in Edit mode")

    start = perf_counter()
    positions, corner_vertices, polygon_starts, triangles = mesh_buffers(mesh)
    if len(polygon_starts) == 0:
        raise UnwrapFailed("Mesh '{}' has no faces".format(mesh.name))
    seams = seam_edges(mesh) if use_seams else None
    names = [name for name in (pinned_group, forbidden_group) if name]
    members = dict(zip(names, groups_vertices(obj, names)))
    pinned = members.get(pinned_group) if pinned_group else None
    forbidden = members.get(forbidden_group) if forbidden_group else None
    extracted = perf_counter()

    try:
        result = ddg.unwrap(
            positions, corner_vertices, polygon_starts, triangles,
            seam_edges=seams,
            pinned_vertices=pinned,
            forbidden_vertices=forbidden,
            lambda_=lambda_,
            positive_cones_only=positive_cones_only,
        )
        solved = perf_counter()
        uvs, pack_info = ddg.pack_charts(
            result.uvs, polygon_starts, result.polygon_charts, margin=PACK_MARGIN,
        )
    except ddg.DdgError as ex:
        raise UnwrapFailed(str(ex)) from ex
    packed = perf_counter()

    uv_layer = mesh.uv_layers.active
    created = uv_layer is None
    if created:
        uv_layer = mesh.uv_layers.new(name=DEFAULT_UV_NAME, do_init=False)
        if uv_layer is None:
            raise UnwrapFailed("Could not create a UV map on mesh '{}'".format(mesh.name))
    uv_layer.uv.foreach_set("vector", uvs)
    written = perf_counter()

    seams_marked = not_found = 0
    if mark_cuts_as_seams:
        seams_marked, not_found = mark_seams(mesh, result.cut_edges)
    mesh.update()
    done = perf_counter()

    return {
        "report": result.report,
        "pack": pack_info,
        "uv_layer": uv_layer.name,
        "created_uv_layer": created,
        "n_cones": len(result.cone_vertices),
        "n_cut_edges": len(result.cut_edges) // 2,
        "seams_marked": seams_marked,
        "cut_edges_not_found": not_found,
        "seconds_extract": extracted - start,
        "seconds_native": solved - extracted,
        "seconds_pack": packed - solved,
        "seconds_write": written - packed,
        "seconds_seams": done - written,
    }


def summary(stats):
    """One line for Blender's status bar from the dict ``unwrap_object`` returns."""
    from math import degrees, exp

    report = stats["report"]
    parts = []
    if report.get("placeholder"):
        parts.append("PLACEHOLDER BUILD (planar projection, not the conformal unwrap):")
    count = "{} chart{}, {} cone{}, {} cut edge{}".format(
        report.get("n_charts", 0), "" if report.get("n_charts") == 1 else "s",
        stats["n_cones"], "" if stats["n_cones"] == 1 else "s",
        stats["n_cut_edges"], "" if stats["n_cut_edges"] == 1 else "s",
    )
    parts.append(count)
    if report.get("measured_triangles"):
        parts.append("| angle error max {:.2f}°, mean {:.2f}°".format(
            degrees(report["measured_angle_error_max"]),
            degrees(report["measured_angle_error_mean"]),
        ))
        parts.append("| scale {:.2f}x to {:.2f}x of the mean".format(
            exp(report["measured_log_scale_min"]), exp(report["measured_log_scale_max"]),
        ))
    flipped = report.get("n_flipped_uv_triangles", 0)
    collapsed = report.get("n_collapsed_uv_triangles", 0)
    if flipped or collapsed:
        parts.append("| {} flipped, {} collapsed triangles".format(flipped, collapsed))
    blender_side = stats["seconds_extract"] + stats["seconds_write"] + stats["seconds_seams"]
    native_side = stats["seconds_native"] + stats["seconds_pack"]
    parts.append("| {:.3f} s (native {:.3f} s, Blender side {:.3f} s)".format(
        blender_side + native_side, native_side, blender_side,
    ))
    if stats["seams_marked"]:
        parts.append("| {} seams marked".format(stats["seams_marked"]))
    if stats["created_uv_layer"]:
        parts.append("| new UV map '{}'".format(stats["uv_layer"]))
    return " ".join(parts)
