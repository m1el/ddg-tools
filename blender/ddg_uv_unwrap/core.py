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
    """A zero-filled ``array.array`` of ``count`` items of 4 bytes each ('f' or 'i')."""
    return array(typecode, bytes(4 * count))


def mesh_buffers(mesh):
    """The arrays that describe ``mesh`` to the native module.

    Returns ``(positions, corner_vertices, polygon_starts, triangles)``:

    * ``positions``: ``array('f')``, 3 per vertex (object space);
    * ``corner_vertices``: ``array('i')``, the vertex of every face corner ("loop");
    * ``polygon_starts``: ``array('i')``, the first corner of every face;
    * ``triangles``: ``array('i')``, 3 corner indices per triangle of Blender's own
      tessellation (``mesh.loop_triangles``), so the native side triangulates n-gons exactly
      like Blender does.
    """
    positions = _zeros('f', 3 * len(mesh.vertices))
    mesh.vertices.foreach_get("co", positions)
    corner_vertices = _zeros('i', len(mesh.loops))
    mesh.loops.foreach_get("vertex_index", corner_vertices)
    polygon_starts = _zeros('i', len(mesh.polygons))
    mesh.polygons.foreach_get("loop_start", polygon_starts)
    mesh.calc_loop_triangles()
    triangles = _zeros('i', 3 * len(mesh.loop_triangles))
    mesh.loop_triangles.foreach_get("loops", triangles)
    return positions, corner_vertices, polygon_starts, triangles


def _edge_flags(mesh, attribute):
    """One byte (0 or 1) per edge for a boolean edge property."""
    flags = bytearray(len(mesh.edges))
    # A bool buffer is what `foreach_get` copies into directly; `array` has no bool typecode.
    mesh.edges.foreach_get(attribute, memoryview(flags).cast('?'))
    return flags


def seam_edges(mesh):
    """The edges marked as seams, as ``array('i')`` with 2 vertex indices per edge."""
    flags = _edge_flags(mesh, "use_seam")
    out = array('i')
    i = flags.find(1)
    if i < 0:
        return out
    edge_vertices = _zeros('i', 2 * len(mesh.edges))
    mesh.edges.foreach_get("vertices", edge_vertices)
    while i >= 0:
        out.extend(edge_vertices[2 * i:2 * i + 2])
        i = flags.find(1, i + 1)
    return out


def mark_seams(mesh, cut_edges):
    """Marks the edges in ``cut_edges`` (2 vertex indices per edge) as seams. Existing seams are
    kept. Returns ``(newly_marked, not_found)``; ``not_found`` counts vertex pairs that are not
    an edge of the mesh."""
    if len(cut_edges) == 0:
        return 0, 0
    n_edges = len(mesh.edges)
    edge_vertices = _zeros('i', 2 * n_edges)
    mesh.edges.foreach_get("vertices", edge_vertices)
    lookup = dict(zip(zip(edge_vertices[0::2], edge_vertices[1::2]), range(n_edges)))
    flags = _edge_flags(mesh, "use_seam")
    newly_marked = not_found = 0
    for k in range(0, len(cut_edges), 2):
        a, b = cut_edges[k], cut_edges[k + 1]
        edge = lookup.get((a, b))
        if edge is None:
            edge = lookup.get((b, a))
        if edge is None:
            not_found += 1
        elif not flags[edge]:
            flags[edge] = 1
            newly_marked += 1
    if newly_marked:
        mesh.edges.foreach_set("use_seam", memoryview(flags).cast('?'))
    return newly_marked, not_found


def group_vertices(obj, group_name):
    """The vertices of ``obj`` that are in the vertex group ``group_name`` with a weight above
    zero, as ``array('i')`` in ascending order."""
    group = obj.vertex_groups.get(group_name)
    if group is None:
        raise UnwrapFailed(
            "Vertex group '{}' does not exist on object '{}'".format(group_name, obj.name)
        )
    index = group.index
    out = array('i')
    for vertex in obj.data.vertices:
        for element in vertex.groups:
            if element.group == index:
                if element.weight > 0.0:
                    out.append(vertex.index)
                break
    return out


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
    pinned = group_vertices(obj, pinned_group) if pinned_group else None
    forbidden = group_vertices(obj, forbidden_group) if forbidden_group else None
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
