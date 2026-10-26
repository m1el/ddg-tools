# SPDX-License-Identifier: MIT
"""Reading curve networks from Blender objects and writing the flow's result back.

Used by the Repulsive Curves operators (``repulsive_operator.py``). The algorithm is
Yu, Schumacher, Crane, Repulsive Curves, ACM TOG 40(2), 2021 (``ddg_tools.RepulsiveSession``).
Like the rest of the add-on this module uses only ``bpy``, ``mathutils`` and the standard library
(the system Blender of some distributions has no numpy).

Inputs
------
* **Mesh objects** with edges only: vertices and edges are the curve network (junctions allowed).
  A mesh with faces is rejected with a message that gives the count; "Loose Edges Only" uses only
  the edges that belong to no face (and says so). Vertices without an edge are rejected unless
  "Ignore Loose Vertices" is on.
* **Curve objects** (poly, Bezier, NURBS splines). Every spline is one component, a cyclic spline
  is a closed loop. The polyline is the control points (Bezier: optionally ``samples`` points per
  segment, evaluated on the curve; NURBS: the control polygon).

Space
-----
The flow runs in **world space** (the default; obstacles are transformed into it and the energy
sees the object's scale as it would physically) or in the **object's local space** (obstacles are
transformed into it). ``Network.matrix`` is the matrix from local coordinates to that space.

Write-back
----------
* Mesh: the new positions are written in local coordinates. When the topology changed (live
  subdivision) the mesh is rebuilt from vertices and edges and the vertex groups are reapplied
  through ``vertex_source``; a vertex created by subdivision gets the average weight of the two
  vertices it was inserted between when both are in the group.
* Curve: poly splines (exact, same point count) by default. "Keep Bezier": Bezier splines stay
  Bezier splines, re-fitted through the new polyline (Catmull-Rom tangents, aligned handles:
  approximate); NURBS splines stay NURBS with moved control points.
"""

from array import array
from math import isfinite, sqrt

import bpy
from mathutils import Matrix, Vector, geometry


class CurveInputError(Exception):
    """The object cannot be used; the message is meant for the user."""


# ----------------------------------------------------------------------------------------------
# Spaces and flat arrays (3 reals per point)

IDENTITY = (1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)


def space_matrix(obj, space):
    """The 3x4 matrix (12 floats, row major) from the local coordinates of ``obj`` to the space
    of the flow."""
    if space == 'WORLD':
        return matrix_rows(obj.matrix_world)
    return IDENTITY


def matrix_rows(m):
    return tuple(m[i][j] for i in range(3) for j in range(4))


def invert(m):
    """The inverse of a 3x4 affine matrix given as 12 floats."""
    full = Matrix((m[0:4], m[4:8], m[8:12], (0.0, 0.0, 0.0, 1.0))).inverted()
    return matrix_rows(full)


def transform(flat, m, typecode='d'):
    """The points of the flat sequence ``flat`` transformed by the 3x4 matrix ``m``; an
    ``array`` of ``typecode``."""
    m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11 = m
    xs, ys, zs = flat[0::3], flat[1::3], flat[2::3]
    out = array(typecode, bytes(array(typecode).itemsize * len(flat)))
    out[0::3] = array(typecode, [m0 * x + m1 * y + m2 * z + m3 for x, y, z in zip(xs, ys, zs)])
    out[1::3] = array(typecode, [m4 * x + m5 * y + m6 * z + m7 for x, y, z in zip(xs, ys, zs)])
    out[2::3] = array(typecode, [m8 * x + m9 * y + m10 * z + m11 for x, y, z in zip(xs, ys, zs)])
    return out


def point(flat, i):
    return (flat[3 * i], flat[3 * i + 1], flat[3 * i + 2])


def gather(flat, indices):
    """The points ``indices`` of the flat ``array`` ``flat``, as an array of the same type."""
    out = array(flat.typecode)
    for i in indices:
        out.extend(flat[3 * i:3 * i + 3])
    return out


def extent(flat):
    """The largest side of the bounding box of the points in ``flat``."""
    if not len(flat):
        return 0.0
    return max(max(flat[k::3]) - min(flat[k::3]) for k in range(3))


def distance(a, b):
    return sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


# ----------------------------------------------------------------------------------------------
# The network


class Network:
    """A curve network read from an object, with what is needed to write a result back.

    Attributes: ``obj``, ``kind`` ('MESH' or 'CURVE'), ``positions`` (flat ``array('d')`` in the
    space of the flow), ``edges`` (flat ``array('i')``), ``matrix`` / ``inverse`` (3x4, local ->
    flow space and back), ``pinned`` and ``tangent`` (sorted vertex indices), ``warnings``
    (strings for the user), ``description`` (what was read, for the report).
    """

    def __init__(self, obj, kind, space):
        self.obj = obj
        self.kind = kind
        self.space = space
        self.matrix = space_matrix(obj, space)
        self.inverse = invert(self.matrix)
        # world coordinates -> the space of the flow (for the obstacles)
        self.world_to_space = IDENTITY if space == 'WORLD' else invert(matrix_rows(obj.matrix_world))
        self.positions = array('d')
        self.edges = array('i')
        self.pinned = []
        self.tangent = []
        self.warnings = []
        self.description = ""
        self.n_vertices = 0  # of the object as currently written
        self.n_edges = 0
        self.n_vertices_original = 0
        self.n_edges_original = 0
        # mesh
        self.vertex_map = None  # network vertex -> mesh vertex
        self.whole_mesh = False  # the network is the whole mesh (may be rebuilt)
        self.groups = {}  # vertex group name -> {network vertex: weight}
        self.original_co = None  # float32 flat, all mesh vertices, before the flow
        self.backup = None  # a copy of the mesh, when the topology may change
        # curve
        self.splines = []
        self.keep_bezier = False
        self._built = False

    def label(self, v):
        """Index of network vertex ``v`` as the user sees it (mesh vertex or curve vertex)."""
        return str(int(self.vertex_map[v])) if self.vertex_map is not None else str(v)

    def write(self, positions, edges=None, vertex_source=None):
        """Writes the flow's ``positions`` (flat, in the space of the flow) to the object.
        ``edges`` (flat) and ``vertex_source`` are needed when the topology changed
        (subdivision); they are the session's."""
        if self.kind == 'MESH':
            self._write_mesh(positions, edges, vertex_source)
        else:
            self._write_curve(positions, edges)
        self.obj.update_tag()

    def restore(self):
        """Puts the object back exactly as it was read."""
        if self.kind == 'MESH':
            self._restore_mesh()
        else:
            self._restore_curve()
        self.obj.update_tag()

    def discard_backup(self):
        if self.backup is not None:
            bpy.data.meshes.remove(self.backup)
            self.backup = None

    # ------------------------------------------------------------------ mesh

    def _set_mesh_co(self, co):
        """Writes the flat float32 ``co`` as the vertex coordinates (through the position
        attribute when possible: ``foreach_set("co")`` may write into an array that an undo
        step still shares)."""
        mesh = self.obj.data
        attribute = mesh.attributes.get("position")
        if attribute is not None and attribute.domain == 'POINT' \
                and attribute.data_type == 'FLOAT_VECTOR' and 3 * len(attribute.data) == len(co):
            attribute.data.foreach_set("vector", co)
        else:
            mesh.vertices.foreach_set("co", co)
        mesh.update()

    def _write_mesh(self, positions, edges, vertex_source):
        mesh = self.obj.data
        n = len(positions) // 3
        local = transform(positions, self.inverse, 'f')
        changed = n != self.n_vertices or (edges is not None and len(edges) // 2 != self.n_edges)
        if not changed:
            if self.whole_mesh:
                self._set_mesh_co(local)
            else:
                co = array('f', self.original_co)
                for k, v in enumerate(self.vertex_map):
                    co[3 * v:3 * v + 3] = local[3 * k:3 * k + 3]
                self._set_mesh_co(co)
            return
        if not self.whole_mesh or edges is None or vertex_source is None:
            raise CurveInputError("Internal error: the topology changed on a mesh that is only "
                                  "partly used")
        mesh.clear_geometry()
        mesh.from_pydata([tuple(local[3 * i:3 * i + 3]) for i in range(n)],
                         [(edges[2 * i], edges[2 * i + 1]) for i in range(len(edges) // 2)], [])
        mesh.update()
        self.n_vertices, self.n_edges = n, len(edges) // 2
        self._reapply_groups(edges, vertex_source, n)

    def _reapply_groups(self, edges, vertex_source, n):
        if not self.groups:
            return
        parents = created_parents(edges, self.n_vertices_original, n)
        for name, original in self.groups.items():
            weights = {}
            for v in range(n):
                source = vertex_source[v]
                if source >= 0:
                    if source in original:
                        weights[v] = original[source]
                else:
                    pair = parents.get(v)
                    if pair is not None and pair[0] in weights and pair[1] in weights:
                        weights[v] = 0.5 * (weights[pair[0]] + weights[pair[1]])
            group = self.obj.vertex_groups.get(name)
            if group is None:
                group = self.obj.vertex_groups.new(name=name)
            by_weight = {}
            for v, weight in weights.items():
                by_weight.setdefault(weight, []).append(v)
            for weight, vertices in by_weight.items():
                group.add(vertices, weight, 'REPLACE')

    def _restore_mesh(self):
        mesh = self.obj.data
        if self.backup is not None and (len(mesh.vertices) != self.n_vertices_original
                                        or len(mesh.edges) != self.n_edges_original):
            import bmesh

            bm = bmesh.new()
            bm.from_mesh(self.backup)
            bm.to_mesh(mesh)
            bm.free()
            mesh.update()
        else:
            self._set_mesh_co(self.original_co)
        self.n_vertices, self.n_edges = self.n_vertices_original, self.n_edges_original

    # ------------------------------------------------------------------ curve

    def _target_type(self, spline):
        if self.keep_bezier and spline["type"] in ('BEZIER', 'NURBS'):
            return spline["type"]
        return 'POLY'

    def _write_curve(self, positions, edges):
        curve = self.obj.data
        n = len(positions) // 3
        if edges is not None and (n != self.n_vertices or len(edges) // 2 != self.n_edges):
            chains = self._chains_after_subdivision(edges)
        else:
            chains = [range(s["offset"], s["offset"] + s["count"]) for s in self.splines]
        local = transform(positions, self.inverse, 'd')
        if self._built and all(len(c) == s["built_count"] for c, s in zip(chains, self.splines)):
            for spline, chain, blender_spline in zip(self.splines, chains, curve.splines):
                _set_spline_points(blender_spline, spline, gather(local, chain))
            curve.update_tag()
            return
        curve.splines.clear()
        for spline, chain in zip(self.splines, chains):
            _build_spline(curve, spline, self._target_type(spline), gather(local, chain))
            spline["built_count"] = len(chain)
        self._built = True
        self.n_vertices = n
        if edges is not None:
            self.n_edges = len(edges) // 2

    def _chains_after_subdivision(self, edges):
        """Per spline, the vertices in order along it, from the edges of a subdivided network
        (old vertices keep their indices, so the first vertex of every spline is where it was)."""
        adjacency = {}
        for i in range(len(edges) // 2):
            a, b = edges[2 * i], edges[2 * i + 1]
            adjacency.setdefault(a, []).append(b)
            adjacency.setdefault(b, []).append(a)
        chains = []
        for spline in self.splines:
            start = spline["offset"]
            walk = [start]
            previous, current = None, start
            nxt = adjacency[start][0]
            while True:
                walk.append(nxt)
                previous, current = current, nxt
                options = [x for x in adjacency[current] if x != previous]
                if not options or options[0] == start:
                    break
                nxt = options[0]
            if spline["cyclic"]:
                # the walk may go either way round; the old vertices must come in order
                old = [v for v in walk if start <= v < start + spline["count"]]
                if len(old) > 1 and old[1] != start + 1:
                    walk = [walk[0]] + walk[1:][::-1]
            chains.append(walk)
        return chains

    def _restore_curve(self):
        curve = self.obj.data
        curve.splines.clear()
        for spline in self.splines:
            _rebuild_from_capture(curve, spline)
        self._built = False
        self.n_vertices, self.n_edges = self.n_vertices_original, self.n_edges_original


def created_parents(edges, n_old, n):
    """For every vertex ``v >= n_old`` created by edge midpoints, the two vertices it was
    inserted between: undo the rounds from the highest index down (the vertices of the last
    round are not adjacent to each other, and their neighbours have lower indices)."""
    adjacency = [set() for _ in range(n)]
    for i in range(len(edges) // 2):
        a, b = int(edges[2 * i]), int(edges[2 * i + 1])
        adjacency[a].add(b)
        adjacency[b].add(a)
    parents = {}
    for v in range(n - 1, n_old - 1, -1):
        neighbors = sorted(adjacency[v])
        if len(neighbors) != 2:
            continue
        a, b = neighbors
        parents[v] = (a, b)
        adjacency[a].discard(v)
        adjacency[b].discard(v)
        adjacency[a].add(b)
        adjacency[b].add(a)
    return parents


# ----------------------------------------------------------------------------------------------
# Splines


def _get(points, key, width, typecode='f'):
    buf = array(typecode, bytes(array(typecode).itemsize * width * len(points)))
    points.foreach_get(key, buf)
    return buf


def _capture_spline(spline):
    """Everything needed to rebuild ``spline`` exactly."""
    kind = spline.type
    info = {
        "type": kind,
        "cyclic": bool(spline.use_cyclic_u),
        "resolution": spline.resolution_u,
        "material": spline.material_index,
        "smooth": bool(spline.use_smooth),
        "order": spline.order_u,
        "endpoint": bool(spline.use_endpoint_u),
        "bezier": bool(spline.use_bezier_u),
    }
    if kind == 'BEZIER':
        pts = spline.bezier_points
        for key in ("co", "handle_left", "handle_right"):
            info[key] = _get(pts, key, 3)
        info["handle_left_type"] = [p.handle_left_type for p in pts]
        info["handle_right_type"] = [p.handle_right_type for p in pts]
        info["select"] = [bool(p.select_control_point) for p in pts]
        info["select_left"] = [bool(p.select_left_handle) for p in pts]
        info["select_right"] = [bool(p.select_right_handle) for p in pts]
    else:
        pts = spline.points
        info["co"] = _get(pts, "co", 4)
        info["select"] = [bool(p.select) for p in pts]
    for key in ("tilt", "radius", "weight_softbody"):
        info[key] = _get(pts, key, 1)
    info["n_control"] = len(pts)
    return info


def _new_spline(curve, info, kind, n):
    spline = curve.splines.new(kind)
    if kind == 'BEZIER':
        spline.bezier_points.add(n - 1)
    else:
        spline.points.add(n - 1)
    spline.use_cyclic_u = info["cyclic"]
    spline.resolution_u = info["resolution"]
    spline.material_index = info["material"]
    spline.use_smooth = info["smooth"]
    if kind == 'NURBS':
        spline.order_u = min(info["order"], max(n, 2))
        spline.use_endpoint_u = info["endpoint"]
        spline.use_bezier_u = info["bezier"]
    return spline


def _rebuild_from_capture(curve, info):
    kind = info["type"]
    spline = _new_spline(curve, info, kind, info["n_control"])
    if kind == 'BEZIER':
        pts = spline.bezier_points
        for i, p in enumerate(pts):
            p.handle_left_type = info["handle_left_type"][i]
            p.handle_right_type = info["handle_right_type"][i]
        for key in ("co", "handle_left", "handle_right"):
            pts.foreach_set(key, info[key])
        for i, p in enumerate(pts):
            p.select_control_point = info["select"][i]
            p.select_left_handle = info["select_left"][i]
            p.select_right_handle = info["select_right"][i]
    else:
        pts = spline.points
        pts.foreach_set("co", info["co"])
        for i, p in enumerate(pts):
            p.select = info["select"][i]
    for key in ("tilt", "radius", "weight_softbody"):
        pts.foreach_set(key, info[key])
    if kind == 'NURBS':
        spline.order_u = info["order"]


def _catmull_handles(points, cyclic):
    """Aligned Bezier handles through ``points`` (flat): tangent (next - previous) / 2, a third
    of it on each side (one-sided at the ends of an open spline). Returns the left and the right
    handles, flat ``array('f')``."""
    n = len(points) // 3
    left, right = array('f', bytes(4 * len(points))), array('f', bytes(4 * len(points)))
    for i in range(n):
        if n == 1:
            t = (0.0, 0.0, 0.0)
        elif not cyclic and i == 0:
            t = [points[3 + k] - points[k] for k in range(3)]
        elif not cyclic and i == n - 1:
            t = [points[3 * i + k] - points[3 * (i - 1) + k] for k in range(3)]
        else:
            j, h = (i + 1) % n, (i - 1) % n
            t = [0.5 * (points[3 * j + k] - points[3 * h + k]) for k in range(3)]
        for k in range(3):
            left[3 * i + k] = points[3 * i + k] - t[k] / 3.0
            right[3 * i + k] = points[3 * i + k] + t[k] / 3.0
    return left, right


def _build_spline(curve, info, kind, points):
    """Creates a spline of type ``kind`` through ``points`` (flat, local coordinates)."""
    n = len(points) // 3
    spline = _new_spline(curve, info, kind, n)
    if kind == 'BEZIER':
        pts = spline.bezier_points
        for p in pts:
            p.handle_left_type = 'ALIGNED'
            p.handle_right_type = 'ALIGNED'
        _set_spline_points(spline, info, points)
    else:
        pts = spline.points
        if info["type"] == 'NURBS' and info["n_control"] == n:
            w = info["co"][3::4]
        else:
            w = array('f', [1.0]) * n
        co = array('f', bytes(16 * n))
        for k in range(3):
            co[k::4] = array('f', points[k::3])
        co[3::4] = w
        pts.foreach_set("co", co)
    if info["n_control"] == n:
        for key in ("tilt", "radius", "weight_softbody"):
            pts.foreach_set(key, info[key])
    return spline


def _set_spline_points(spline, info, points):
    """Moves the points of an already built spline (same count) to ``points`` (flat)."""
    n = len(points) // 3
    if spline.type == 'BEZIER':
        left, right = _catmull_handles(points, info["cyclic"])
        pts = spline.bezier_points
        pts.foreach_set("co", array('f', points))
        pts.foreach_set("handle_left", left)
        pts.foreach_set("handle_right", right)
    else:
        pts = spline.points
        co = _get(pts, "co", 4)
        for k in range(3):
            co[k::4] = array('f', points[k::3])
        pts.foreach_set("co", co)
    return n


def _spline_polyline(spline, samples):
    """The points of one spline (flat, local coordinates) and the vertex of every control
    point."""
    cyclic = bool(spline.use_cyclic_u)
    if spline.type == 'BEZIER':
        pts = spline.bezier_points
        n = len(pts)
        co = _get(pts, "co", 3)
        if samples <= 1 or n < 2:
            return array('d', co), list(range(n))
        hl, hr = _get(pts, "handle_left", 3), _get(pts, "handle_right", 3)
        out, control = array('d'), []
        for i in range(n if cyclic else n - 1):
            j = (i + 1) % n
            control.append(len(out) // 3)
            curve_points = geometry.interpolate_bezier(
                Vector(co[3 * i:3 * i + 3]), Vector(hr[3 * i:3 * i + 3]),
                Vector(hl[3 * j:3 * j + 3]), Vector(co[3 * j:3 * j + 3]), samples + 1)
            for p in curve_points[:samples]:
                out.extend(p)
        if not cyclic:
            control.append(len(out) // 3)
            out.extend(co[-3:])
        return out, control
    pts = spline.points
    co = _get(pts, "co", 4)
    out = array('d', bytes(8 * 3 * len(pts)))
    for k in range(3):
        out[k::3] = array('d', co[k::4])
    return out, list(range(len(pts)))


# ----------------------------------------------------------------------------------------------
# Reading


def vertex_group_weights(obj, name):
    """``{vertex: weight}`` of a vertex group of a mesh object (members only)."""
    group = obj.vertex_groups.get(name)
    if group is None:
        raise CurveInputError("Object '{}' has no vertex group '{}'".format(obj.name, name))
    index = group.index
    out = {}
    for v in obj.data.vertices:
        for g in v.groups:
            if g.group == index:
                out[v.index] = g.weight
    return out


def read_object(obj, *, space='WORLD', pin_group="", tangent_group="", pin_tangents=False,
                ignore_loose=False, loose_edges_only=False, bezier_samples=1, keep_bezier=False,
                keep_groups=False, keep_backup=False):
    """Reads the curve network of ``obj`` (not in Edit mode). ``keep_groups``: record all vertex
    groups (for live subdivision); ``keep_backup``: copy the mesh so that Esc can restore a
    changed topology. Raises ``CurveInputError`` with a message for the user."""
    if obj is None:
        raise CurveInputError("There is no active object")
    if obj.type == 'MESH':
        net = _read_mesh(obj, space, pin_group, tangent_group, pin_tangents, ignore_loose,
                         loose_edges_only, keep_groups, keep_backup)
    elif obj.type == 'CURVE':
        net = _read_curve(obj, space, pin_tangents, bezier_samples, keep_bezier)
    else:
        raise CurveInputError(
            "Object '{}' is a {}; select an edge-only mesh or a Curve object (poly, Bezier or "
            "NURBS)".format(obj.name, obj.type.lower()))
    net.n_vertices_original = net.n_vertices
    net.n_edges_original = net.n_edges
    try:
        _check_network(net)
    except CurveInputError:
        net.discard_backup()
        raise
    return net


def _check_network(net):
    pos = net.positions
    n = len(pos) // 3
    for i in range(n):
        if not (isfinite(pos[3 * i]) and isfinite(pos[3 * i + 1]) and isfinite(pos[3 * i + 2])):
            raise CurveInputError("Vertex {} of '{}' has a position that is not finite".format(
                net.label(i), net.obj.name))
    if not len(net.edges):
        raise CurveInputError("Object '{}' has no edges to flow".format(net.obj.name))
    # Coincident vertices (and so zero-length edges): the flow cannot start; say where.
    tolerance = max(extent(pos), 1e-12) * 1e-9
    seen, clashes, first = {}, 0, None
    for i in range(n):
        key = (round(pos[3 * i] / tolerance), round(pos[3 * i + 1] / tolerance),
               round(pos[3 * i + 2] / tolerance))
        j = seen.setdefault(key, i)
        if j != i:
            clashes += 1
            if first is None:
                first = (j, i)
    if clashes:
        raise CurveInputError(
            "'{}' has coincident vertices ({} and {} are at the same place; {} vertices share a "
            "position with an earlier one). Merge them first (Mesh > Clean Up > Merge by "
            "Distance, or remove the duplicate spline points)".format(
                net.obj.name, net.label(first[0]), net.label(first[1]), clashes))


def _read_mesh(obj, space, pin_group, tangent_group, pin_tangents, ignore_loose, loose_edges_only,
               keep_groups, keep_backup):
    mesh = obj.data
    net = Network(obj, 'MESH', space)
    if mesh.library is not None:
        raise CurveInputError("Mesh '{}' is linked from another file".format(mesh.name))
    if mesh.shape_keys is not None:
        raise CurveInputError(
            "Mesh '{}' has shape keys; the flow would move only the basis. Apply or remove the "
            "shape keys first".format(mesh.name))
    n, m = len(mesh.vertices), len(mesh.edges)
    if n == 0 or m == 0:
        raise CurveInputError("Mesh '{}' is empty ({} vertices, {} edges): there is no curve to "
                              "flow".format(mesh.name, n, m))
    faces = len(mesh.polygons)
    if faces and not loose_edges_only:
        raise CurveInputError(
            "Mesh '{}' has {} faces. The flow needs an edge-only mesh (a curve network). Delete "
            "the faces (Edit mode: Delete > Only Faces), or turn on 'Loose Edges Only' to use "
            "just the edges that belong to no face".format(mesh.name, faces))
    co = array('f', bytes(12 * n))
    mesh.vertices.foreach_get("co", co)
    pairs = array('I', bytes(8 * m))
    mesh.edges.foreach_get("vertices", pairs)
    keep = [True] * m
    if faces:
        loop_edge = array('I', bytes(4 * len(mesh.loops)))
        mesh.loops.foreach_get("edge_index", loop_edge)
        for e in loop_edge:
            keep[e] = False
        if not any(keep):
            raise CurveInputError("Mesh '{}' has {} faces and no loose edges: nothing to flow"
                                  .format(mesh.name, faces))
        net.warnings.append("Using only the {} loose edges of '{}'; its {} faces and the {} "
                            "edges around them stay where they are".format(
                                sum(keep), mesh.name, faces, m - sum(keep)))
    edges = [(pairs[2 * i], pairs[2 * i + 1]) for i in range(m) if keep[i]]
    used = [False] * n
    for a, b in edges:
        used[a] = used[b] = True
    unused = [v for v in range(n) if not used[v]]
    if faces:
        vertex_map = [v for v in range(n) if used[v]]
    elif not unused:
        vertex_map = list(range(n))
    elif ignore_loose:
        vertex_map = [v for v in range(n) if used[v]]
        net.warnings.append("Ignoring {} vertices of '{}' that belong to no edge".format(
            len(unused), mesh.name))
    else:
        raise CurveInputError(
            "Mesh '{}' has {} vertices that belong to no edge (the first is vertex {}). Delete "
            "them (Mesh > Clean Up > Delete Loose), or turn on 'Ignore Loose Vertices'".format(
                mesh.name, len(unused), unused[0]))
    remap = {v: k for k, v in enumerate(vertex_map)}
    net.vertex_map = vertex_map
    net.whole_mesh = len(vertex_map) == n and not faces
    net.edges = array('i', [remap[v] for pair in edges for v in pair])
    net.positions = transform(gather(co, vertex_map), net.matrix, 'd')
    net.original_co = co
    net.n_vertices, net.n_edges = len(vertex_map), len(edges)
    net.description = "mesh '{}': {} vertices, {} edges".format(
        obj.name, net.n_vertices, net.n_edges)

    def pick(name):
        weights = vertex_group_weights(obj, name)
        return sorted(remap[v] for v, w in weights.items() if w > 0.0 and v in remap)

    net.pinned = pick(pin_group) if pin_group else []
    net.tangent = pick(tangent_group) if tangent_group else []
    if pin_tangents:
        net.tangent = sorted(set(net.tangent) | set(net.pinned))
    if keep_groups and net.whole_mesh:
        for group in obj.vertex_groups:
            net.groups[group.name] = vertex_group_weights(obj, group.name)
    if keep_backup and net.whole_mesh:
        net.backup = mesh.copy()
        net.backup.name = ".ddg_repulsive_backup"
    return net


def _read_curve(obj, space, pin_tangents, samples, keep_bezier):
    curve = obj.data
    net = Network(obj, 'CURVE', space)
    net.keep_bezier = keep_bezier
    if curve.library is not None:
        raise CurveInputError("Curve '{}' is linked from another file".format(curve.name))
    if len(curve.splines) == 0:
        raise CurveInputError(
            "Curve '{}' has no splines: there is nothing to flow".format(curve.name))
    points, edges, pinned = array('d'), [], []
    for k, spline in enumerate(curve.splines):
        info = _capture_spline(spline)
        polyline, control = _spline_polyline(spline, max(1, samples))
        n = len(polyline) // 3
        if n < 2:
            raise CurveInputError(
                "Spline {} of '{}' has {} point; a spline needs at least 2 points to be a curve"
                .format(k, curve.name, n))
        if info["cyclic"] and n < 3:
            info["cyclic"] = False
            net.warnings.append("Cyclic spline {} of '{}' has only {} points; treated as open"
                                .format(k, curve.name, n))
        offset = len(points) // 3
        info.update(offset=offset, count=n, built_count=n, control=control)
        points.extend(polyline)
        edges.extend(v for i in range(n - 1) for v in (offset + i, offset + i + 1))
        if info["cyclic"]:
            edges.extend((offset + n - 1, offset))
        pinned.extend(offset + control[i] for i, flag in enumerate(info["select"]) if flag)
        net.splines.append(info)
    net.positions = transform(points, net.matrix, 'd')
    net.edges = array('i', edges)
    net.pinned = sorted(pinned)
    if pin_tangents:
        net.tangent = list(net.pinned)
    net.n_vertices, net.n_edges = len(points) // 3, len(edges) // 2
    net.description = "curve '{}': {} splines, {} points".format(
        obj.name, len(net.splines), net.n_vertices)
    if curve.dimensions == '2D':
        net.warnings.append(
            "Curve '{}' is 2D: the z coordinate is dropped when it is displayed. Turn on Planar "
            "(axis Z, local space) to keep it flat".format(curve.name))
    return net


# ----------------------------------------------------------------------------------------------
# Obstacles


def read_obstacles(objects, net, depsgraph):
    """``[(positions, triangles)]`` (flat ``array('d')`` / ``array('i')``) of the mesh objects
    ``objects``, evaluated (modifiers applied) and transformed into the space of ``net``.
    Raises ``CurveInputError`` naming the object when a mesh is unusable."""
    out = []
    for ob in objects:
        if ob is net.obj:
            continue
        if ob.type != 'MESH':
            raise CurveInputError(
                "Obstacle '{}' is a {}, not a mesh".format(ob.name, ob.type.lower()))
        evaluated = ob.evaluated_get(depsgraph)
        mesh = evaluated.to_mesh()
        try:
            mesh.calc_loop_triangles()
            nv, nt = len(mesh.vertices), len(mesh.loop_triangles)
            if nt == 0:
                raise CurveInputError(
                    "Obstacle '{}' has no faces ({} vertices): an obstacle needs triangles or "
                    "polygons to block the curve".format(ob.name, nv))
            co = array('f', bytes(12 * nv))
            mesh.vertices.foreach_get("co", co)
            tri = array('I', bytes(12 * nt))
            mesh.loop_triangles.foreach_get("vertices", tri)
        finally:
            evaluated.to_mesh_clear()
        pts = transform(transform(co, matrix_rows(ob.matrix_world)), net.world_to_space, 'd')
        if not all(isfinite(x) for x in pts):
            raise CurveInputError("Obstacle '{}' has vertices that are not finite".format(ob.name))
        scale = max(extent(pts), 1e-12)
        degenerate = []
        for t in range(nt):
            a, b, c = (point(pts, tri[3 * t + k]) for k in range(3))
            u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
            v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
            cross = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                     u[0] * v[1] - u[1] * v[0])
            if 0.5 * sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2) <= 1e-12 * scale ** 2:
                degenerate.append(t)
        if degenerate:
            raise CurveInputError(
                "Obstacle '{}' has {} degenerate (zero-area) triangles, the first is triangle "
                "{}. Fix it first (Mesh > Clean Up > Degenerate Dissolve, or Merge by Distance)"
                .format(ob.name, len(degenerate), degenerate[0]))
        out.append((pts, array('i', tri)))
    return out


def obstacle_objects(object_name, collection_name):
    """The mesh objects named by the two obstacle properties (without duplicates)."""
    found = []
    if object_name:
        obj = bpy.data.objects.get(object_name)
        if obj is None:
            raise CurveInputError("The obstacle object '{}' does not exist".format(object_name))
        found.append(obj)
    if collection_name:
        collection = bpy.data.collections.get(collection_name)
        if collection is None:
            raise CurveInputError(
                "The obstacle collection '{}' does not exist".format(collection_name))
        found.extend(o for o in collection.all_objects if o.type == 'MESH')
    seen, out = set(), []
    for o in found:
        if o.name not in seen:
            seen.add(o.name)
            out.append(o)
    return out


def session_arguments(net):
    """The keyword arguments of ``RepulsiveSession`` that come from the network."""
    kwargs = {}
    if net.pinned:
        kwargs["pinned"] = array('i', net.pinned)
    if net.tangent:
        kwargs["tangent_vertices"] = array('i', net.tangent)
    return kwargs
