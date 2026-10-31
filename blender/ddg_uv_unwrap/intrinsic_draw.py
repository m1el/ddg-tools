# SPDX-License-Identifier: MIT
"""Drawing intrinsic triangulations on their input mesh in Blender. No operator here: the
functions are reused by Show Intrinsic Refinement (``intrinsic_operator``) and, later, by the
FlipOut operators' "Show triangulation" (``geodesic_operator``).

``draw_intrinsic_edges`` turns the polylines of ``ddg_tools.IntrinsicTriangulation.edge_polylines``
(or any ``points`` / ``offsets`` pair of the same layout, e.g. ``steiner_tree``'s) into a Curve
object. ``polylines_from_edge_paths`` is the same for paths given as lists of 3D points.
"""

import colorsys

import bpy

# Custom property that marks an object made by one of these functions and names the object it
# was made for; re-running replaces exactly these objects.
SOURCE_KEY = "ddg_intrinsic_of"
KIND_KEY = "ddg_intrinsic_kind"


def remove_previous(source, kind):
    """Removes the objects (and their data) that were made for ``source`` as ``kind``
    (``"mesh"``, ``"edges"`` or another ``kind``). Returns how many were removed."""
    removed = 0
    for obj in list(bpy.data.objects):
        if obj.get(SOURCE_KEY) == source.name and obj.get(KIND_KEY) == kind:
            data = obj.data
            bpy.data.objects.remove(obj, do_unlink=True)
            if data is not None and data.users == 0:
                if isinstance(data, bpy.types.Mesh):
                    bpy.data.meshes.remove(data)
                elif isinstance(data, bpy.types.Curve):
                    bpy.data.curves.remove(data)
            removed += 1
    return removed


def _link_like(source, obj):
    collections = list(source.users_collection) or [bpy.context.scene.collection]
    for collection in collections:
        collection.objects.link(obj)


def draw_intrinsic_edges(source, points, offsets, *, name=None, radius=0.0, replace=True,
                         kind="edges"):
    """Creates a Curve object with one poly spline per polyline and returns it.

    ``points``: flat sequence of x, y, z (object space of ``source``); ``offsets``: where each
    polyline starts in points (in points), one more entry at the end (the layout of
    ``IntrinsicTriangulation.edge_polylines``). The object gets ``source``'s transform, so the
    lines lie on ``source``'s surface. ``radius`` > 0 gives the curve a round bevel of that
    radius (otherwise it shows as a thin line). With ``replace`` the curve made earlier for
    ``source`` is removed first. ``source`` is never changed. ``kind`` is the value of the marker
    property (``"edges"`` for Show Intrinsic Refinement; another operator passes its own so that
    the two do not replace each other's objects).
    """
    name = name or source.name + "_intrinsic_edges"
    if replace:
        remove_previous(source, kind)
    curve = bpy.data.curves.new(name, 'CURVE')
    curve.dimensions = '3D'
    curve.bevel_depth = float(radius)
    curve.resolution_u = 1
    for k in range(len(offsets) - 1):
        a, b = offsets[k], offsets[k + 1]
        n = b - a
        if n < 2:
            continue
        spline = curve.splines.new('POLY')
        spline.points.add(n - 1)
        flat = []
        for i in range(a, b):
            flat.extend((points[3 * i], points[3 * i + 1], points[3 * i + 2], 1.0))
        spline.points.foreach_set("co", flat)
    obj = bpy.data.objects.new(name, curve)
    obj.matrix_world = source.matrix_world.copy()
    obj[SOURCE_KEY] = source.name
    obj[KIND_KEY] = kind
    obj.hide_select = False
    _link_like(source, obj)
    return obj


def polylines_from_paths(paths):
    """``(points, offsets)`` from a list of paths, each a list of (x, y, z)."""
    points, offsets = [], [0]
    for path in paths:
        for p in path:
            points.extend(p)
        offsets.append(offsets[-1] + len(path))
    return points, offsets


def face_colors(n, seed_offset=0.0):
    """``n`` distinct pleasant colors (RGBA tuples), by the golden-ratio walk around the hue
    circle."""
    golden = 0.6180339887498949
    colors = []
    for i in range(n):
        h = (seed_offset + i * golden) % 1.0
        s = 0.35 + 0.25 * ((i * 7) % 5) / 4.0
        v = 0.75 + 0.2 * ((i * 3) % 4) / 3.0
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        colors.append((r, g, b, 1.0))
    return colors
