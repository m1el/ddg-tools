# SPDX-License-Identifier: MIT
"""Reading point clouds from Blender objects, the tori cache, and writing results back.

Used by the Points as Tori operator (``points_as_tori_operator.py``). The method is

    Feng, Gkioulekas, Crane, Signed Distance to Point Clouds by Blending Tori (Points as Tori),
    ACM Transactions on Graphics 45(4), SIGGRAPH 2026,

implemented in ``ddg_tools.points_as_tori`` (crates/ddg-points-as-tori, crates/ddg-py). Like the
rest of the add-on this module uses only ``bpy``, ``mathutils`` and the standard library.

Input
-----
A point cloud is a mesh object: its vertices count, edges and faces are ignored. The mesh data
itself is read (modifiers are not applied). The normals come from, in this order:

1. the point-domain FLOAT_VECTOR attribute named ``normal`` (what Blender's PLY importer makes of
   ``nx, ny, nz``, Blender 4.x and 5.x), or the name given in the panel;
2. three point-domain FLOAT attributes ``nx``, ``ny``, ``nz``;
3. the custom normals of the mesh (averaged per vertex);
4. the vertex normals of a mesh that has faces (computed by Blender from the faces).

Positions are taken to world space with the object's matrix, normals with the inverse transpose of
its linear part (and normalized). The result mesh is written in the object's frame.

The cache
---------
The tori are built once per cloud. The key is a BLAKE2 hash of the world-space positions and
normals (and the sign rule, which changes the fit). A hit is served from memory (this session,
the last few clouds) or from a sidecar file ``<name>.pat/<hash>-<sign rule>.tori`` next to the
.blend file (``save_tori``; ``pat_load_tori`` checks the cloud hash again). An unsaved file keeps
the tori in memory only. Changing positions or normals changes the hash, which is the whole
invalidation.
"""

import hashlib
import os
import threading
from array import array
from collections import OrderedDict
from math import isfinite, sqrt

import bpy
from mathutils import Matrix

from .repulsive_io import invert, matrix_rows, transform

SIDECAR_SUFFIX = ".pat"
MEMORY_ENTRIES = 4
NORMAL_NAME = "normal"


class CloudInputError(Exception):
    """The object cannot be used; the message is meant for the user."""


# ----------------------------------------------------------------------------------------------
# Reading


def _point_attribute(mesh, name, data_type):
    attr = mesh.attributes.get(name)
    if attr is None or attr.domain != 'POINT' or attr.data_type != data_type:
        return None
    return attr


def _vectors(attr, count):
    out = array('f', bytes(12 * count))
    attr.data.foreach_get("vector", out)
    return out


def _floats(attr, count):
    out = array('f', bytes(4 * count))
    attr.data.foreach_get("value", out)
    return out


def find_normals(mesh, name=""):
    """``(kind, label, local normals as flat array('f'))`` or ``None``. ``kind`` is 'attribute'
    (``label`` its name), 'nxyz', 'custom' or 'vertex'."""
    n = len(mesh.vertices)
    wanted = name or NORMAL_NAME
    attr = _point_attribute(mesh, wanted, 'FLOAT_VECTOR')
    if attr is not None:
        return 'attribute', wanted, _vectors(attr, n)
    if name:
        return None
    comps = [_point_attribute(mesh, c, 'FLOAT') for c in ("nx", "ny", "nz")]
    if all(c is not None for c in comps):
        xs, ys, zs = (_floats(c, n) for c in comps)
        out = array('f', bytes(12 * n))
        out[0::3], out[1::3], out[2::3] = xs, ys, zs
        return 'nxyz', "nx, ny, nz", out
    if mesh.has_custom_normals and len(mesh.polygons) > 0:
        out = array('f', bytes(12 * n))
        corner = array('f', bytes(12 * len(mesh.loops)))
        mesh.corner_normals.foreach_get("vector", corner)
        verts = array('i', bytes(4 * len(mesh.loops)))
        mesh.loops.foreach_get("vertex_index", verts)
        for i, v in enumerate(verts):
            out[3 * v] += corner[3 * i]
            out[3 * v + 1] += corner[3 * i + 1]
            out[3 * v + 2] += corner[3 * i + 2]
        return 'custom', "custom normals", out
    if len(mesh.polygons) > 0:
        out = array('f', bytes(12 * n))
        mesh.vertex_normals.foreach_get("vector", out)
        return 'vertex', "vertex normals", out
    return None


class Cloud:
    """An oriented point cloud read from an object.

    ``positions`` / ``normals``: flat ``array('d')`` in world space (normals unit length);
    ``local_normals``: the normals as stored (``array('f')``), for the write-back; ``kind`` /
    ``label``: where they came from; ``matrix`` / ``inverse``: 3x4 object <-> world;
    ``digest``: the cache key (hex).
    """

    def __init__(self, obj, positions, normals, local_normals, kind, label, matrix):
        self.obj = obj
        self.name = obj.name
        self.positions = positions
        self.normals = normals
        self.local_normals = local_normals
        self.kind = kind
        self.label = label
        self.matrix = matrix
        self.inverse = invert(matrix)
        self.n = len(positions) // 3
        self.digest = cloud_digest(positions, normals)

    def bounds(self):
        xs, ys, zs = self.positions[0::3], self.positions[1::3], self.positions[2::3]
        return (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs))


def cloud_digest(positions, normals):
    h = hashlib.blake2b(digest_size=16)
    h.update(len(positions).to_bytes(8, "little"))
    h.update(memoryview(positions).cast("B"))
    h.update(memoryview(normals).cast("B"))
    return h.hexdigest()


def read_cloud(obj, normal_attribute=""):
    """The point cloud of ``obj``. Raises ``CloudInputError`` with a message for the user."""
    if obj is None or obj.type != 'MESH':
        raise CloudInputError("{} is not a mesh object".format(
            "'{}'".format(obj.name) if obj is not None else "the active object"))
    mesh = obj.data
    n = len(mesh.vertices)
    if n == 0:
        raise CloudInputError("object '{}' has no vertices".format(obj.name))
    found = find_normals(mesh, normal_attribute)
    if found is None:
        if normal_attribute:
            raise CloudInputError(
                "object '{}' has {:,} vertices but no point attribute '{}' with normals".format(
                    obj.name, n, normal_attribute))
        raise CloudInputError(
            "object '{}' has {:,} vertices but no normals; add a point-domain Vector attribute "
            "named 'normal' (a PLY file with nx, ny, nz imports that way), or use a mesh with "
            "faces (its vertex normals are used)".format(obj.name, n))
    kind, label, local_normals = found
    local = array('f', bytes(12 * n))
    mesh.attributes["position"].data.foreach_get("vector", local)

    world = obj.matrix_world
    matrix = matrix_rows(world)
    positions = transform(local, matrix)
    linear = world.to_3x3()
    try:
        normal_matrix = linear.inverted().transposed()
    except ValueError as ex:
        raise CloudInputError("object '{}' has a singular transform (scale 0 on an axis)".format(
            obj.name)) from ex
    nm = Matrix.to_4x4(normal_matrix)
    normals = transform(local_normals, matrix_rows(nm))
    bad = 0
    for i in range(0, len(normals), 3):
        x, y, z = normals[i], normals[i + 1], normals[i + 2]
        length = sqrt(x * x + y * y + z * z)
        if not isfinite(length) or length < 1e-12:
            bad += 1
            continue
        normals[i], normals[i + 1], normals[i + 2] = x / length, y / length, z / length
    if bad:
        raise CloudInputError(
            "object '{}': {:,} of {:,} normals (from {}) are zero or not finite".format(
                obj.name, bad, n, label))
    if not all(isfinite(x) for x in positions[:3] + positions[-3:]):
        raise CloudInputError("object '{}' has non-finite positions".format(obj.name))
    return Cloud(obj, positions, normals, local_normals, kind, label, matrix)


def other_cloud_objects(context, active):
    """Selected mesh objects other than ``active``."""
    return [o for o in context.selected_objects if o is not active and o.type == 'MESH']


# ----------------------------------------------------------------------------------------------
# The cache

_LOCK = threading.Lock()
_MEMORY = OrderedDict()


def sidecar_dir():
    """The sidecar directory of the current .blend file, or None when it is not saved."""
    path = bpy.data.filepath
    if not path:
        return None
    base, _ext = os.path.splitext(path)
    return base + SIDECAR_SUFFIX


def sidecar_path(directory, digest, sign_rule):
    return os.path.join(directory, "{}-{}.tori".format(digest, sign_rule))


def memory_get(key):
    with _LOCK:
        field = _MEMORY.get(key)
        if field is not None:
            _MEMORY.move_to_end(key)
        return field


def memory_put(key, field):
    with _LOCK:
        _MEMORY[key] = field
        _MEMORY.move_to_end(key)
        while len(_MEMORY) > MEMORY_ENTRIES:
            _MEMORY.popitem(last=False)


def memory_entries():
    with _LOCK:
        return list(_MEMORY)


def clear_cache(directory=None):
    """Forgets the clouds in memory and deletes the ``.tori`` files of ``directory``. Returns
    ``(memory entries dropped, files deleted)``."""
    with _LOCK:
        dropped = len(_MEMORY)
        _MEMORY.clear()
    deleted = 0
    if directory and os.path.isdir(directory):
        for name in os.listdir(directory):
            if name.endswith(".tori"):
                os.remove(os.path.join(directory, name))
                deleted += 1
        try:
            os.rmdir(directory)
        except OSError:
            pass
    return dropped, deleted


# ----------------------------------------------------------------------------------------------
# Writing


def write_flips(cloud, flipped):
    """Negates the stored normals of the points in ``flipped`` (an iterable of vertex indices),
    in the place they were read from. A mesh whose normals were custom normals gets them back
    through ``normals_split_custom_set_from_vertices``; one that used the vertex normals computed
    from its faces gets a point attribute ``normal`` (which is read first from then on). Returns
    a short description of the place written."""
    mesh = cloud.obj.data
    flips = set(flipped)
    local = array('f', cloud.local_normals)
    for i in flips:
        local[3 * i] = -local[3 * i]
        local[3 * i + 1] = -local[3 * i + 1]
        local[3 * i + 2] = -local[3 * i + 2]
    if cloud.kind == 'attribute':
        mesh.attributes[cloud.label].data.foreach_set("vector", local)
        place = "attribute '{}'".format(cloud.label)
    elif cloud.kind == 'nxyz':
        for k, name in enumerate(("nx", "ny", "nz")):
            mesh.attributes[name].data.foreach_set("value", local[k::3])
        place = "attributes nx, ny, nz"
    elif cloud.kind == 'custom':
        # The custom normals are per corner; set them from the new per-vertex normals.
        vectors = [tuple(local[i:i + 3]) for i in range(0, len(local), 3)]
        mesh.normals_split_custom_set_from_vertices(vectors)
        place = "custom normals"
    else:
        attr = mesh.attributes.get(NORMAL_NAME)
        if attr is not None and (attr.domain != 'POINT' or attr.data_type != 'FLOAT_VECTOR'):
            raise CloudInputError(
                "cannot create the point attribute 'normal' on '{}': an attribute of that name "
                "exists with another type".format(cloud.name))
        if attr is None:
            attr = mesh.attributes.new(NORMAL_NAME, 'FLOAT_VECTOR', 'POINT')
        attr.data.foreach_set("vector", local)
        place = "a new point attribute 'normal' (it now takes precedence over the face normals)"
    mesh.update()
    return place


def new_surface_object(cloud, name, vertices, faces, matrix_world, collections):
    """A mesh object ``name`` holding the triangles (world-space ``vertices``, flat ``array('d')``;
    ``faces``, flat triangle indices) in the frame of ``cloud``'s object, linked next to it."""
    local = transform(vertices, cloud.inverse, 'f')
    n_faces = len(faces) // 3
    mesh = bpy.data.meshes.new(name)
    mesh.vertices.add(len(local) // 3)
    mesh.attributes["position"].data.foreach_set("vector", local)
    mesh.polygons.add(n_faces)
    mesh.loops.add(3 * n_faces)
    mesh.polygons.foreach_set("loop_start", array('i', range(0, 3 * n_faces, 3)))
    mesh.loops.foreach_set("vertex_index", array('i', faces))
    mesh.update(calc_edges=True)
    mesh.validate()
    mesh.shade_smooth()
    obj = bpy.data.objects.new(name, mesh)
    obj.matrix_world = matrix_world
    for collection in collections or [bpy.context.collection]:
        collection.objects.link(obj)
    return obj
