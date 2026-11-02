# SPDX-License-Identifier: MIT
"""Reading triangle surfaces from Blender objects and writing the repulsive-surface flow back.

Used by the Repulsive Surfaces operators (``repulsive_surface_operator.py``). The algorithm is
Yu, Brakensiek, Schumacher, Crane, Repulsive Surfaces, ACM TOG 40(6), 2021
(``ddg_tools.RepulsiveSurfaceSession``). Like the rest of the add-on this module uses only
``bpy``, ``mathutils`` and the standard library.

Input
-----
The active mesh object's *evaluated* mesh (modifiers applied), triangulated by Blender's own loop
triangles: n-gons are split exactly as Blender splits them, the library never changes the caller's
polygons. A vertex is a mesh vertex (no vertex is merged or dropped), so a vertex that belongs to
no face, a non-manifold edge and so on are rejected by the library with a located message.

Space
-----
The flow runs in world space (default: the object's scale counts as it would physically,
obstacles are used as they stand) or in the object's local space. ``SurfaceMesh.matrix`` is the
matrix from local coordinates to the flow's space.

Write-back
----------
* **Same triangles** (remeshing off, or no operation was accepted): the new positions are written
  into the existing mesh; polygons, UV maps, colours, vertex groups and everything else stay.
* **Topology changed** (remeshing): the mesh data is replaced by the result's triangles. The
  material index and the smooth flag of a triangle are those of the polygon that its ancestor
  triangle belonged to (``face_origin``). UV maps, colour attributes and other attributes are
  dropped, and so are vertex groups other than the pin group (which is rebuilt on the pinned
  vertices, which the library never moves or removes); the report says so (user decision: no
  attribute transfer, as for ``ddg-developability`` D24).
* If the evaluated mesh differs from the stored mesh (modifiers that add or remove vertices or
  faces) the result is always a replacement, and the report warns that the modifiers' result was
  baked into the mesh data.

Obstacles
---------
Mesh objects (``read_mesh_obstacles``) become ``mesh_obstacle`` penalties (tangent-point). Empties
(``read_empties``) become ``implicit_obstacle`` / ``implicit_attractor`` penalties, the shape
chosen by the empty's display type and its size by ``empty_display_size`` times its scale:

=================================  ===================================================
Sphere                             sphere, radius = size * mean scale
Cone                               infinite cylinder along local Z, radius = size * mean(x, y) scale
Circle                             torus around local Z: major radius = size * mean(x, y) scale,
                                   minor = major * ``ddg_minor_ratio``
Plain Axes, Arrows, Single Arrow,  plane through the origin, normal = local Z (the side the
Cube, Image                        normal points to is the free side)
=================================  ===================================================

Custom properties of the empty override the panel's defaults: ``ddg_mode`` (``"repel"`` or
``"attract"``), ``ddg_power``, ``ddg_weight``, ``ddg_minor_ratio``.
"""

from array import array
from math import isfinite, sqrt

import bpy
from mathutils import Matrix, Vector
from mathutils.bvhtree import BVHTree

from .repulsive_io import (
    IDENTITY,
    extent,
    invert,
    matrix_rows,
    point,
    transform,
)

BACKUP_NAME = ".ddg_repulsive_surface_backup"

# Mesh data that Blender or the add-on rebuilds from the faces (not "dropped by a replacement").
_STRUCTURAL = {"position", "sharp_face", "sharp_edge", "material_index", "uv_seam"}

EMPTY_SHAPES = {
    'SPHERE': "sphere",
    'CONE': "cylinder",
    'CIRCLE': "torus",
}


class SurfaceInputError(Exception):
    """The object cannot be used; the message is meant for the user."""


def space_matrix(obj, space):
    return matrix_rows(obj.matrix_world) if space == 'WORLD' else IDENTITY


def _user_attributes(mesh):
    names = [layer.name for layer in mesh.uv_layers]
    names += [a.name for a in mesh.color_attributes]
    skip = set(names) | _STRUCTURAL
    names += [a.name for a in mesh.attributes if not a.name.startswith(".") and a.name not in skip]
    return names


def group_weights(mesh, index):
    """``{vertex: weight}`` of the vertex group number ``index`` of ``mesh`` (members only)."""
    out = {}
    for v in mesh.vertices:
        for g in v.groups:
            if g.group == index:
                out[v.index] = g.weight
    return out


class SurfaceMesh:
    """A triangle surface read from an object, with what is needed to write a result back.

    Attributes: ``obj``, ``positions`` (flat ``array('d')`` in the flow's space), ``faces`` (flat
    ``array('i')``, 3 vertex ids per loop triangle), ``pinned`` (sorted vertex ids),
    ``pin_weights`` (their weights in the group), ``matrix`` / ``inverse`` (3x4, local -> flow
    space and back), ``world_to_space``, ``warnings`` (strings for the user), ``description``,
    ``baked`` (the evaluated mesh differs from the stored one: every write replaces the data),
    ``replaced`` (the object currently holds replaced data).
    """

    def __init__(self, obj, space):
        self.obj = obj
        self.space = space
        self.matrix = space_matrix(obj, space)
        self.inverse = invert(self.matrix)
        self.world_to_space = IDENTITY if space == 'WORLD' else invert(matrix_rows(obj.matrix_world))
        self.positions = array('d')
        self.faces = array('i')
        self.pinned = []
        self.pin_weights = []
        self.pin_group = ""
        self.warnings = []
        self.description = ""
        self.baked = False
        self.replaced = False
        self.n_vertices_original = 0
        self.n_polygons_original = 0
        self.original_co = None
        self.backup = None
        self.triangle_polygon = []  # polygon (of the evaluated mesh) of every input triangle
        self.material = []
        self.smooth = []
        self.dropped = []  # names of the data a replacement drops
        self.other_groups = []
        self.group_order = []  # vertex group names in their original order
        self.groups = {}  # name -> {vertex: weight} of the stored mesh, for a restore
        self.warned = False

    # ------------------------------------------------------------------ writing

    def write(self, positions, faces=None, pinned=None, face_origin=None):
        """Writes ``positions`` (flat, in the flow's space). ``faces``, ``pinned`` and
        ``face_origin`` come from the same ``snapshot()``; without ``faces`` the triangles are
        taken to be the input's."""
        local = transform(positions, self.inverse, 'f')
        n = len(positions) // 3
        same = (not self.baked and not self.replaced and n == self.n_vertices_original
                and (faces is None or faces == self.faces))
        if same:
            self._set_co(local)
        else:
            self._replace(local, faces, pinned, face_origin)
        self.obj.update_tag()

    def _set_co(self, co):
        """Vertex coordinates through the position attribute when possible (``foreach_set("co")``
        may write into an array that an undo step still shares)."""
        mesh = self.obj.data
        attribute = mesh.attributes.get("position")
        if attribute is not None and attribute.domain == 'POINT' \
                and attribute.data_type == 'FLOAT_VECTOR' and 3 * len(attribute.data) == len(co):
            attribute.data.foreach_set("vector", co)
        else:
            mesh.vertices.foreach_set("co", co)
        mesh.update()

    def _replace(self, local, faces, pinned, face_origin):
        if faces is None or face_origin is None:
            raise SurfaceInputError("Internal error: the topology changed without faces")
        mesh = self.obj.data
        n = len(local) // 3
        coordinates = [(local[i], local[i + 1], local[i + 2]) for i in range(0, len(local), 3)]
        triangles = [(faces[i], faces[i + 1], faces[i + 2]) for i in range(0, len(faces), 3)]
        mesh.clear_geometry()
        mesh.from_pydata(coordinates, [], triangles)
        count = len(triangles)
        source = [self.triangle_polygon[face_origin[t]] for t in range(count)]
        mesh.polygons.foreach_set("material_index", array('i', (self.material[p] for p in source)))
        mesh.polygons.foreach_set(
            "use_smooth", memoryview(bytearray(self.smooth[p] for p in source)).cast('?'))
        mesh.update()
        # Blender 5 stores the group names in the mesh, so clear_geometry may have removed them
        # (older versions keep empty groups): keep exactly the pin group, rebuilt on the pins.
        for group in list(self.obj.vertex_groups):
            if group.name != self.pin_group:
                self.obj.vertex_groups.remove(group)
        if self.pin_group and pinned is not None and len(pinned) == len(self.pin_weights):
            group = self.obj.vertex_groups.get(self.pin_group) \
                or self.obj.vertex_groups.new(name=self.pin_group)
            by_weight = {}
            for v, w in zip(pinned, self.pin_weights):
                by_weight.setdefault(w, []).append(int(v))
            for w, vertices in by_weight.items():
                group.add(vertices, w, 'REPLACE')
        self.replaced = True
        if not self.warned:
            self.warned = True
            if self.baked:
                self.warnings.append(
                    "'{}' has modifiers that change its vertices or faces: the evaluated mesh "
                    "was used and its result replaced the mesh data (the modifiers apply to it "
                    "again; remove or apply them)".format(self.obj.name))
            if self.dropped:
                self.warnings.append(
                    "The mesh data of '{}' was replaced by the remeshed triangles; dropped (no "
                    "attribute transfer): {}".format(self.obj.name, ", ".join(self.dropped)))
            else:
                self.warnings.append(
                    "The mesh data of '{}' was replaced by the remeshed triangles (material "
                    "index and smooth flag come from the parent polygon)".format(self.obj.name))
        return n

    def restore(self):
        """Puts the object back exactly as it was read."""
        mesh = self.obj.data
        if self.replaced or len(mesh.vertices) != self.n_vertices_original \
                or len(mesh.polygons) != self.n_polygons_original:
            if self.backup is None:
                raise SurfaceInputError("Internal error: no backup to restore from")
            import bmesh

            bm = bmesh.new()
            bm.from_mesh(self.backup)
            bm.to_mesh(mesh)
            bm.free()
            mesh.update()
            self._restore_groups()
            self.replaced = False
        else:
            self._set_co(self.original_co)
        self.obj.update_tag()

    def _restore_groups(self):
        obj = self.obj
        for group in list(obj.vertex_groups):
            obj.vertex_groups.remove(group)
        for name in self.group_order:
            group = obj.vertex_groups.new(name=name)
            by_weight = {}
            for v, w in self.groups.get(name, {}).items():
                by_weight.setdefault(w, []).append(v)
            for w, vertices in by_weight.items():
                group.add(vertices, w, 'REPLACE')

    def discard_backup(self):
        if self.backup is not None:
            bpy.data.meshes.remove(self.backup)
            self.backup = None


def read_object(obj, depsgraph, *, space='WORLD', pin_group=""):
    """The ``SurfaceMesh`` of the mesh object ``obj``. Raises ``SurfaceInputError``."""
    if obj is None or obj.type != 'MESH':
        raise SurfaceInputError("The active object is not a mesh")
    mesh = obj.data
    if mesh.library is not None:
        raise SurfaceInputError("Mesh '{}' is linked from another file".format(mesh.name))
    if mesh.shape_keys is not None:
        raise SurfaceInputError(
            "Mesh '{}' has shape keys; the flow would move only the basis. Apply or remove the "
            "shape keys first".format(mesh.name))
    if len(mesh.polygons) == 0:
        raise SurfaceInputError("Mesh '{}' has no faces ({} vertices): the flow needs a triangle "
                                "surface".format(mesh.name, len(mesh.vertices)))
    surface = SurfaceMesh(obj, space)
    surface.n_vertices_original = len(mesh.vertices)
    surface.n_polygons_original = len(mesh.polygons)
    co0 = array('f', bytes(12 * len(mesh.vertices)))
    mesh.vertices.foreach_get("co", co0)
    surface.original_co = co0
    pin_index = None
    if pin_group:
        group = obj.vertex_groups.get(pin_group)
        if group is None:
            raise SurfaceInputError("Object '{}' has no vertex group '{}'".format(
                obj.name, pin_group))
        pin_index = group.index
    evaluated = obj.evaluated_get(depsgraph)
    em = evaluated.to_mesh()
    try:
        em.calc_loop_triangles()
        nv, nt, npoly = len(em.vertices), len(em.loop_triangles), len(em.polygons)
        if nt == 0:
            raise SurfaceInputError("Mesh '{}' has no faces after its modifiers".format(obj.name))
        co = array('f', bytes(12 * nv))
        em.vertices.foreach_get("co", co)
        tri = array('I', bytes(12 * nt))
        em.loop_triangles.foreach_get("vertices", tri)
        tri_polygon = array('I', bytes(4 * nt))
        em.loop_triangles.foreach_get("polygon_index", tri_polygon)
        material = array('i', bytes(4 * npoly))
        em.polygons.foreach_get("material_index", material)
        smooth = bytearray(npoly)
        em.polygons.foreach_get("use_smooth", memoryview(smooth).cast('?'))
        weights = group_weights(em, pin_index) if pin_index is not None else {}
    finally:
        evaluated.to_mesh_clear()
    surface.baked = nv != surface.n_vertices_original or npoly != surface.n_polygons_original
    surface.positions = transform(co, surface.matrix, 'd')
    if not all(isfinite(x) for x in surface.positions):
        raise SurfaceInputError("Mesh '{}' has vertices that are not finite".format(obj.name))
    surface.faces = array('i', tri)
    surface.triangle_polygon = list(tri_polygon)
    surface.material = list(material)
    surface.smooth = list(smooth)
    surface.pin_group = pin_group
    pins = sorted(v for v, w in weights.items() if w > 0.0)
    surface.pinned = pins
    surface.pin_weights = [weights[v] for v in pins]
    surface.dropped = _user_attributes(mesh)
    surface.group_order = [g.name for g in obj.vertex_groups]
    surface.groups = {g.name: group_weights(mesh, g.index) for g in obj.vertex_groups}
    surface.other_groups = [g.name for g in obj.vertex_groups if g.name != pin_group]
    if surface.other_groups:
        surface.dropped += ["vertex group " + name for name in surface.other_groups]
    surface.description = "mesh '{}': {} vertices, {} triangles{}".format(
        obj.name, nv, nt, ", {} pinned".format(len(pins)) if pins else "")
    surface.backup = mesh.copy()
    surface.backup.name = BACKUP_NAME
    return surface


# ----------------------------------------------------------------------------------------------
# Self-intersection of the input (SPEC 16: CCD is off for inputs that already intersect)


def intersecting_pairs(positions, faces, limit=1):
    """Up to ``limit`` pairs ``(i, j)`` of triangles that overlap although they share no vertex
    (an approximation by the triangle BVH of Blender: pairs whose triangles intersect)."""
    n = len(positions) // 3
    verts = [point(positions, i) for i in range(n)]
    tris = [tuple(faces[i:i + 3]) for i in range(0, len(faces), 3)]
    tree = BVHTree.FromPolygons(verts, tris, epsilon=0.0)
    found = []
    for i, j in tree.overlap(tree):
        if i >= j:
            continue
        if set(tris[i]) & set(tris[j]):
            continue
        found.append((i, j))
        if len(found) >= limit:
            break
    return found


# ----------------------------------------------------------------------------------------------
# Obstacles


def objects_of(collection_name, kind):
    """The objects of type ``kind`` in the collection called ``collection_name`` (and its
    children), without duplicates; ``[]`` for no name."""
    if not collection_name:
        return []
    collection = bpy.data.collections.get(collection_name)
    if collection is None:
        raise SurfaceInputError("The collection '{}' does not exist".format(collection_name))
    return [o for o in collection.all_objects if o.type == kind]


def read_mesh_obstacles(objects, surface, depsgraph):
    """``[(positions, triangles)]`` (flat ``array('d')`` / ``array('i')``) of the mesh objects,
    evaluated and transformed into the flow's space."""
    out = []
    for ob in objects:
        if ob is surface.obj:
            continue
        evaluated = ob.evaluated_get(depsgraph)
        mesh = evaluated.to_mesh()
        try:
            mesh.calc_loop_triangles()
            nv, nt = len(mesh.vertices), len(mesh.loop_triangles)
            if nt == 0:
                raise SurfaceInputError(
                    "Obstacle '{}' has no faces ({} vertices)".format(ob.name, nv))
            co = array('f', bytes(12 * nv))
            mesh.vertices.foreach_get("co", co)
            tri = array('I', bytes(12 * nt))
            mesh.loop_triangles.foreach_get("vertices", tri)
        finally:
            evaluated.to_mesh_clear()
        pts = transform(transform(co, matrix_rows(ob.matrix_world)), surface.world_to_space, 'd')
        if not all(isfinite(x) for x in pts):
            raise SurfaceInputError("Obstacle '{}' has vertices that are not finite".format(ob.name))
        scale = max(extent(pts), 1e-12)
        for t in range(nt):
            a, b, c = (point(pts, tri[3 * t + k]) for k in range(3))
            u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
            v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
            cross = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                     u[0] * v[1] - u[1] * v[0])
            if 0.5 * sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2) <= 1e-12 * scale ** 2:
                raise SurfaceInputError(
                    "Obstacle '{}' has a degenerate (zero-area) triangle, the first is triangle "
                    "{}. Fix it first (Mesh > Clean Up > Degenerate Dissolve)".format(ob.name, t))
        out.append((pts, array('i', tri)))
    return out


def _custom(obj, key, default):
    value = obj.get(key, default)
    return default if value is None else value


def implicit_penalty(empty, surface, *, mode, power, weight, minor_ratio):
    """The penalty dict of one empty (see the module docstring); ``mode`` ('REPEL'/'ATTRACT'),
    ``power``, ``weight`` and ``minor_ratio`` are the defaults its custom properties override."""
    to_space = Matrix((surface.world_to_space[0:4], surface.world_to_space[4:8],
                       surface.world_to_space[8:12], (0.0, 0.0, 0.0, 1.0)))
    m = to_space @ empty.matrix_world
    center = m.translation
    axis = (m.to_3x3() @ Vector((0.0, 0.0, 1.0))).normalized()
    scale = m.to_scale()
    size = empty.empty_display_size
    radial = size * 0.5 * (abs(scale.x) + abs(scale.y))
    shape = EMPTY_SHAPES.get(empty.empty_display_type, "plane")
    if shape == "sphere":
        radius = size * (abs(scale.x) + abs(scale.y) + abs(scale.z)) / 3.0
        surf = {"type": "sphere", "center": tuple(center), "radius": radius}
    elif shape == "cylinder":
        surf = {"type": "cylinder", "point": tuple(center), "axis": tuple(axis), "radius": radial}
    elif shape == "torus":
        ratio = float(_custom(empty, "ddg_minor_ratio", minor_ratio))
        surf = {"type": "torus", "center": tuple(center), "axis": tuple(axis), "major": radial,
                "minor": radial * ratio}
    else:
        surf = {"type": "plane", "point": tuple(center), "normal": tuple(axis)}
    mode = str(_custom(empty, "ddg_mode", mode)).upper()
    if mode not in {"REPEL", "ATTRACT"}:
        raise SurfaceInputError("Empty '{}': custom property ddg_mode must be \"repel\" or "
                                "\"attract\", not {!r}".format(empty.name, mode.lower()))
    return {
        "type": "implicit_attractor" if mode == "ATTRACT" else "implicit_obstacle",
        "surface": surf,
        "power": float(_custom(empty, "ddg_power", power)),
        "weight": float(_custom(empty, "ddg_weight", weight)),
    }
