# SPDX-License-Identifier: MIT
"""Headless tests of the variational-cut operators (Start, Flow, Apply). Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/test_cut.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [-k NAME]

`scripts/blender-test.sh` runs it after test_addon.py against the installed extension and its
bundled wheel. For development, `--source DIR --site DIR` register the add-on from a source
directory and import `ddg_tools` from an unpacked wheel instead.

The meshes are generated here (UV sphere, subdivided cube, torus, the monkey) and, when the file
is present, `spot_quadrangulated.obj` of Keenan Crane's model repository (`--spot PATH`; the test
is skipped with a notice otherwise).

What is checked, on every model:

* Start creates one float attribute per region on the face-corner domain (`ddg_cut_region_k`,
  the second the negative of the first) and a wire preview object that follows the transform;
* Flow lowers the energy, keeps the attributes consistent, replaces the preview, and `k` steps
  twice equal `2k` steps once;
* Apply inserts exactly what the Rust reference `CutSession.insert` inserts (same vertex
  positions, same polygon vertex sets, same number of seams, every one of them a cut edge),
  removes the preview and the attributes, keeps the mesh manifold where it was (Euler
  characteristic, no edge with more than two faces), and keeps the data that existed: a UV map
  that is a linear function of position stays that function, so do a vertex group, a vertex and a
  corner attribute, material indices and a face attribute;
* the unwrap after Apply gives UVs with no flipped triangle;
* a second Start, Flow and Apply works with the first cut as fixed seams;
* every error path shows a readable message.
"""

import argparse
import array
import collections
import contextlib
import ctypes
import importlib
import math
import os
import re
import sys
import tempfile
import traceback

import bmesh
import bpy
from mathutils import Vector
from mathutils.kdtree import KDTree

TESTS = []
SPOT = "/var/data/discrete-geometry/ddg-tools/data/models/spot/spot/spot_quadrangulated.obj"


def test(function):
    TESTS.append(function)
    return function


class Skip(Exception):
    pass


# ----------------------------------------------------------------------------------------------
# Helpers


def reset():
    active = bpy.context.view_layer.objects.active
    if active is not None and active.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    for mesh in list(bpy.data.meshes):
        bpy.data.meshes.remove(mesh)
    for material in list(bpy.data.materials):
        bpy.data.materials.remove(material)


def activate(obj):
    view_layer = bpy.context.view_layer
    view_layer.update()
    for other in view_layer.objects:
        other.select_set(other is obj)
    view_layer.objects.active = obj
    return obj


def primitive(kind, **kwargs):
    getattr(bpy.ops.mesh, "primitive_{}_add".format(kind))(**kwargs)
    return bpy.context.view_layer.objects.active


def subdivided_cube(cuts=6):
    cube = primitive("cube", size=2.0)
    bm = bmesh.new()
    bm.from_mesh(cube.data)
    bmesh.ops.subdivide_edges(bm, edges=bm.edges[:], cuts=cuts, use_grid_fill=True)
    bm.to_mesh(cube.data)
    bm.free()
    cube.data.update()
    return cube


def spot_object():
    if not os.path.exists(ARGS.spot):
        raise Skip("{} is not present".format(ARGS.spot))
    before = set(bpy.data.objects)
    bpy.ops.wm.obj_import(filepath=ARGS.spot)
    new = [o for o in bpy.data.objects if o not in before and o.type == 'MESH']
    assert len(new) == 1, new
    return activate(new[0])


def ngon_sphere():
    """An icosphere with patches of about 14 triangles dissolved into polygons of 8 to 12 sides
    (the patches do not touch), so that the network of a cut has several junctions and ordinary
    paths in one polygon."""
    obj = primitive("ico_sphere", subdivisions=3, radius=1.0)
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bm.faces.ensure_lookup_table()
    blocked = set()
    groups = []
    for seed in list(bm.faces):
        if seed in blocked:
            continue
        group, queue = [seed], [seed]
        members = {seed}
        while queue and len(group) < 24:
            face = queue.pop(0)
            for edge in face.edges:
                for other in edge.link_faces:
                    if other not in members and other not in blocked and len(group) < 24:
                        members.add(other)
                        group.append(other)
                        queue.append(other)
        for face in group:
            blocked.add(face)
            for edge in face.edges:
                blocked.update(edge.link_faces)
        groups.append(group)
    for group in groups:
        if len(group) > 6:
            bmesh.ops.dissolve_faces(bm, faces=group, use_verts=True)
    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()
    sides = max(len(p.vertices) for p in obj.data.polygons)
    assert sides >= 8, sides
    return obj


def atlas_islands(mesh, uvs, tolerance=1e-6):
    """Number of UV islands: faces joined across an edge whose two ends have the same UV on both
    sides."""
    parent = list(range(len(mesh.polygons)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    sides = collections.defaultdict(list)
    for polygon in mesh.polygons:
        first, count = polygon.loop_start, polygon.loop_total
        for k in range(count):
            c, d = first + k, first + (k + 1) % count
            sides[mesh.loops[c].edge_index].append((polygon.index, c, d))

    def uv(c):
        return uvs[2 * c], uvs[2 * c + 1]

    for entries in sides.values():
        if len(entries) != 2:
            continue
        (pa, ca, da), (pb, cb, db) = entries
        # The edge runs from vertex x to y in one polygon, and y to x in the other.
        if mesh.loops[ca].vertex_index == mesh.loops[cb].vertex_index:
            pairs = ((ca, cb), (da, db))
        else:
            pairs = ((ca, db), (da, cb))
        if all(abs(uv(a)[0] - uv(b)[0]) < tolerance and abs(uv(a)[1] - uv(b)[1]) < tolerance
               for a, b in pairs):
            ra, rb = find(pa), find(pb)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    return len({find(p) for p in range(len(mesh.polygons))})


MODELS = collections.OrderedDict([
    ("uv_sphere", lambda: primitive("uv_sphere", segments=24, ring_count=12)),
    ("subdivided_cube", subdivided_cube),
    ("torus", lambda: primitive("torus", major_segments=24, minor_segments=12)),
    ("monkey", lambda: primitive("monkey")),
    ("spot", spot_object),
])


def positions(mesh):
    values = array.array('f', bytes(12 * len(mesh.vertices)))
    mesh.vertices.foreach_get("co", values)
    return [Vector(values[3 * i:3 * i + 3]) for i in range(len(mesh.vertices))]


def bounding_box(mesh):
    pts = positions(mesh)
    low = Vector([min(p[k] for p in pts) for k in range(3)])
    high = Vector([max(p[k] for p in pts) for k in range(3)])
    return low, high


def diagonal(mesh):
    low, high = bounding_box(mesh)
    return (high - low).length


def seam_pairs(mesh):
    return [tuple(e.vertices) for e in mesh.edges if e.use_seam]


def state(mesh):
    """The cut attributes as a list of lists (float32 values)."""
    rows = []
    k = 0
    while mesh.attributes.get("ddg_cut_region_{}".format(k)) is not None:
        values = array.array('f', bytes(4 * len(mesh.loops)))
        mesh.attributes["ddg_cut_region_{}".format(k)].data.foreach_get("value", values)
        rows.append(values)
        k += 1
    return rows


def save_state(mesh):
    """The whole stored cut: the attributes, the state bytes if there are any, the option."""
    return ([array.array('f', r) for r in state(mesh)], mesh.get("ddg_cut_state"),
            mesh.get("ddg_cut_option"))


def load_state(mesh, saved):
    rows, bytes_, option = saved
    for k, row in enumerate(rows):
        mesh.attributes["ddg_cut_region_{}".format(k)].data.foreach_set("value", row)
    for key, value in (("ddg_cut_state", bytes_), ("ddg_cut_option", option)):
        if value is None:
            if key in mesh:
                del mesh[key]
        else:
            mesh[key] = value


def attributes_only(mesh):
    """Makes the attributes the state again (drops the state bytes of a refined session)."""
    if "ddg_cut_state" in mesh:
        del mesh["ddg_cut_state"]


def previews(obj):
    return [c for c in obj.children if c.get("ddg_cut_preview")]


def start_circle(obj, fraction=0.15, vertex=0):
    radius = fraction * diagonal(obj.data)
    assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=radius, vertex=vertex) == {'FINISHED'}


CUT_OFFSET = 1e-3    # the default of the cut operators (cut.DEFAULT_CUT_OFFSET)
FLOW = dict(steps=6, use_stop_rule=False, length_mode='EXPLICIT', length_weight=0.1)


def flow(**props):
    settings = dict(FLOW)
    settings.update(props)
    return bpy.ops.uv.ddg_cut_flow(**settings)


def pad(values, fill, n=16):
    return tuple(values) + (fill,) * (n - len(values))


def euler(mesh):
    return len(mesh.vertices) - len(mesh.edges) + len(mesh.polygons)


def edge_face_counts(mesh):
    counts = collections.Counter()
    for polygon in mesh.polygons:
        for loop_index in polygon.loop_indices:
            counts[mesh.loops[loop_index].edge_index] += 1
    return counts


def expect_error(pattern, operator, **props):
    """The operator must report an error that matches `pattern`, without a traceback."""
    try:
        result = operator(**props)
    except RuntimeError as ex:
        message = str(ex)
        assert "Traceback" not in message, "the operator raised instead of reporting:\n" + message
        assert re.search(pattern, message), "error {!r} does not match {!r}".format(
            message, pattern)
        return message
    raise AssertionError("expected an error matching {!r}, got {}".format(pattern, result))


class CaptureOutput:
    """Captures what Blender and Python print (file descriptors 1 and 2), then replays it."""

    def __enter__(self):
        self._flush()
        self._file = tempfile.TemporaryFile()
        self._saved = (os.dup(1), os.dup(2))
        os.dup2(self._file.fileno(), 1)
        os.dup2(self._file.fileno(), 2)
        self.text = ""
        return self

    def __exit__(self, *_exc):
        self._flush()
        for fd, saved in zip((1, 2), self._saved):
            os.dup2(saved, fd)
            os.close(saved)
        self._file.seek(0)
        self.text = self._file.read().decode("utf-8", "replace")
        self._file.close()
        sys.stdout.write(self.text)
        return False

    @staticmethod
    def _flush():
        sys.stdout.flush()
        sys.stderr.flush()
        ctypes.CDLL(None).fflush(None)


def buffers_of(mesh):
    core = importlib.import_module(ARGS.module + ".core")
    return core.mesh_buffers(mesh), core.seam_edges(mesh)


def rust_reference(obj):
    """What the Rust insertion makes of the cut stored on the mesh of `obj` (the session is
    opened exactly as the operators open it)."""
    cutmod = importlib.import_module(ARGS.module + ".cut")
    session, _, _ = cutmod.open_session(obj, CUT_OFFSET)
    return session, session.insert(min_crossing_offset=CUT_OFFSET)


def match_vertices(reference_positions, mesh, tolerance=1e-5):
    """The bijection between the vertices of the Rust insertion and the vertices of the mesh,
    by position. Returns a list: Rust vertex -> mesh vertex."""
    n = len(mesh.vertices)
    assert len(reference_positions) // 3 == n, (
        "{} vertices in the reference, {} in the mesh".format(len(reference_positions) // 3, n))
    tree = KDTree(n)
    for i, p in enumerate(positions(mesh)):
        tree.insert(p, i)
    tree.balance()
    mapping = []
    pts = positions(mesh)
    if any(len(tree.find_range(p, tolerance)) > 1 for p in pts):
        return None         # coincident vertices (the monkey has some): no bijection by position
    for i in range(n):
        p = Vector(reference_positions[3 * i:3 * i + 3])
        found = tree.find_range(p, tolerance)
        assert found, "no vertex of the mesh at the position {} of reference vertex {}".format(
            tuple(p), i)
        mapping.append(min(found, key=lambda f: f[2])[1])
    assert len(set(mapping)) == n, "two reference vertices map to one vertex of the mesh"
    return mapping


def polygon_vertex_sets(mesh):
    return collections.Counter(
        tuple(sorted(mesh.loops[i].vertex_index for i in p.loop_indices)) for p in mesh.polygons)


def expected_seams(old_seam_positions, inserted):
    """The seams after the insertion, by position pairs: the cut, plus the old seams cut into
    pieces at the new vertices that lie on them."""
    pts = [Vector(inserted.positions[3 * i:3 * i + 3]) for i in range(len(inserted.positions) // 3)]
    out = set()
    for a, b in zip(inserted.new_seams[0::2], inserted.new_seams[1::2]):
        out.add(frozenset((a, b)))
    old = [(i, Vector(inserted.positions[3 * i:3 * i + 3])) for i in range(inserted.n_old_vertices)]
    by_position = {tuple(round(c, 4) for c in p): i for i, p in old}
    for pa, pb in old_seam_positions:
        a, b = by_position[tuple(round(c, 4) for c in pa)], by_position[
            tuple(round(c, 4) for c in pb)]
        direction = pts[b] - pts[a]
        length2 = direction.length_squared
        on_edge = []
        for i in range(inserted.n_old_vertices, len(pts)):
            d = pts[i] - pts[a]
            t = d.dot(direction) / length2
            if 1e-6 < t < 1 - 1e-6 and (d - direction * t).length < 1e-4:
                on_edge.append((t, i))
        chain = [a] + [i for _, i in sorted(on_edge)] + [b]
        for u, v in zip(chain, chain[1:]):
            out.add(frozenset((u, v)))
    return out


# ----------------------------------------------------------------------------------------------
# Data that must survive: linear functions of position


class Payload:
    """Data added to a mesh before Apply, each a linear function of the position (so that the
    interpolation along an edge or a diagonal must reproduce it exactly), or inherited."""

    def __init__(self, obj):
        mesh = obj.data
        self.low, self.high = bounding_box(mesh)
        self.size = Vector([max(h - l, 1e-9) for l, h in zip(self.low, self.high)])
        pts = positions(mesh)

        def unit(p, k):
            return (p[k] - self.low[k]) / self.size[k]

        self.unit = unit
        # A UV map (u, v) = (x, y) scaled to the unit square; a second one that stays untouched
        # as a whole (checked by the same function).
        self.uv_names = ["Linear", "Second"]
        for name, (ku, kv) in zip(self.uv_names, ((0, 1), (2, 0))):
            layer = mesh.uv_layers.new(name=name)
            corner_vertices = array.array('i', bytes(4 * len(mesh.loops)))
            mesh.loops.foreach_get("vertex_index", corner_vertices)
            values = array.array('f')
            for v in corner_vertices:
                values.extend((unit(pts[v], ku), unit(pts[v], kv)))
            mesh.attributes[layer.name].data.foreach_set("vector", values)
        # A vertex group with weight a linear function of x.
        group = obj.vertex_groups.new(name="Linear")
        self.group_index = group.index
        for i, p in enumerate(pts):
            group.add([i], 0.1 + 0.8 * unit(p, 0), 'REPLACE')
        # A float attribute per vertex and one per corner, linear in z.
        vertex_attribute = mesh.attributes.new("lin_vert", 'FLOAT', 'POINT')
        vertex_attribute.data.foreach_set("value", array.array('f', (unit(p, 2) for p in pts)))
        corner_vertices = array.array('i', bytes(4 * len(mesh.loops)))
        mesh.loops.foreach_get("vertex_index", corner_vertices)
        corner_attribute = mesh.attributes.new("lin_corner", 'FLOAT', 'CORNER')
        corner_attribute.data.foreach_set(
            "value", array.array('f', (unit(pts[v], 2) for v in corner_vertices)))
        # Materials by face, and the face's own index as a face attribute.
        for k in range(3):
            material = bpy.data.materials.new("payload {}".format(k))
            mesh.materials.append(material)
        self.materials = array.array('i', (p.index % 3 for p in mesh.polygons))
        mesh.polygons.foreach_set("material_index", self.materials)
        parents = mesh.attributes.new("parent_id", 'INT', 'FACE')
        parents.data.foreach_set("value", array.array('i', range(len(mesh.polygons))))
        mesh.update()

    def refresh(self, obj):
        """A new round: the faces are numbered again and carry their current material."""
        mesh = obj.data
        self.materials = array.array('i', bytes(4 * len(mesh.polygons)))
        mesh.polygons.foreach_get("material_index", self.materials)
        mesh.attributes["parent_id"].data.foreach_set(
            "value", array.array('i', range(len(mesh.polygons))))
        mesh.update()

    def check(self, obj, reference):
        """After Apply: the data is still what it was, as a function of position."""
        mesh = obj.data
        pts = positions(mesh)
        corner_vertices = array.array('i', bytes(4 * len(mesh.loops)))
        mesh.loops.foreach_get("vertex_index", corner_vertices)
        for name, (ku, kv) in zip(self.uv_names, ((0, 1), (2, 0))):
            layer = mesh.uv_layers.get(name)
            assert layer is not None, "UV map {} is gone".format(name)
            values = array.array('f', bytes(8 * len(mesh.loops)))
            mesh.attributes[name].data.foreach_get("vector", values)
            worst = 0.0
            for c, v in enumerate(corner_vertices):
                worst = max(worst, abs(values[2 * c] - self.unit(pts[v], ku)),
                            abs(values[2 * c + 1] - self.unit(pts[v], kv)))
            assert worst < 2e-3, "UV map {} is off by {} from its function".format(name, worst)
        group = obj.vertex_groups.get("Linear")
        assert group is not None and group.index == self.group_index, "vertex group is gone"
        worst = 0.0
        for vertex in mesh.vertices:
            weights = {g.group: g.weight for g in vertex.groups}
            assert self.group_index in weights, "vertex {} lost its group weight".format(
                vertex.index)
            worst = max(worst, abs(weights[self.group_index]
                                   - (0.1 + 0.8 * self.unit(pts[vertex.index], 0))))
        assert worst < 2e-3, "vertex group weights are off by {}".format(worst)
        values = array.array('f', bytes(4 * len(mesh.vertices)))
        mesh.attributes["lin_vert"].data.foreach_get("value", values)
        worst = max(abs(values[i] - self.unit(pts[i], 2)) for i in range(len(pts)))
        assert worst < 2e-3, "vertex attribute is off by {}".format(worst)
        values = array.array('f', bytes(4 * len(mesh.loops)))
        mesh.attributes["lin_corner"].data.foreach_get("value", values)
        worst = max(abs(values[c] - self.unit(pts[v], 2)) for c, v in enumerate(corner_vertices))
        assert worst < 2e-3, "corner attribute is off by {}".format(worst)
        assert len(mesh.materials) == 3, "materials changed"
        parent_id = array.array('i', bytes(4 * len(mesh.polygons)))
        mesh.attributes["parent_id"].data.foreach_get("value", parent_id)
        material = array.array('i', bytes(4 * len(mesh.polygons)))
        mesh.polygons.foreach_get("material_index", material)
        for p in range(len(mesh.polygons)):
            assert material[p] == self.materials[parent_id[p]], (
                "face {} (from {}) has material {}, expected {}".format(
                    p, parent_id[p], material[p], self.materials[parent_id[p]]))
        # The faces that came from each old face are those of the reference.
        assert collections.Counter(parent_id) == collections.Counter(reference.polygon_parent), (
            "the faces do not come from the same old faces as in the reference")


# ----------------------------------------------------------------------------------------------
# The whole workflow on one model


def check_apply_against_reference(obj, payload, reference_session, reference, seams_before):
    mesh = obj.data
    mapping = match_vertices(reference.positions, mesh)
    want = expected_seams(seams_before, reference)
    if mapping is None:
        # Only the counts can be compared.
        assert len(mesh.polygons) == len(reference.polygon_starts)
        assert len(mesh.vertices) == len(reference.positions) // 3
        assert len(seam_pairs(mesh)) == len(want), "{} seams, expected {}".format(
            len(seam_pairs(mesh)), len(want))
        print("  (coincident vertices: counts compared, not positions)")
        return
    assert len(mesh.polygons) == len(reference.polygon_starts), (
        "{} polygons, the reference has {}".format(len(mesh.polygons),
                                                   len(reference.polygon_starts)))
    # The same polygons: every reference polygon, translated, is a polygon of the mesh.
    starts = list(reference.polygon_starts) + [len(reference.corner_vertices)]
    expected = collections.Counter(
        tuple(sorted(mapping[v] for v in reference.corner_vertices[starts[p]:starts[p + 1]]))
        for p in range(len(reference.polygon_starts)))
    assert expected == polygon_vertex_sets(mesh), "the polygons differ from the reference"
    # The same seams, and they are exactly the cut plus the old seams.
    want = expected_seams(seams_before, reference)
    got = {frozenset(pair) for pair in seam_pairs(mesh)}
    assert got == {frozenset(mapping[v] for v in pair) for pair in want}, (
        "seams differ from the reference: {} marked, {} expected".format(len(got), len(want)))


def positions_of_seams(mesh):
    pts = positions(mesh)
    return [(tuple(pts[a]), tuple(pts[b])) for a, b in seam_pairs(mesh)]


def workflow(make, name, circle=True, regions=None, seams=False, weight=None, inspect=None):
    obj = make()
    mesh = obj.data
    vertex = 0
    if name == "monkey":      # vertex 0 is in an eye, a small piece of its own
        vertex = min(range(len(mesh.vertices)), key=lambda i: mesh.vertices[i].co.y)
    # A model that is a closed manifold stays one; a model with boundary keeps its boundary.
    counts_before = collections.Counter(edge_face_counts(mesh).values())
    chi_before = euler(mesh)
    if seams:
        for e in mesh.edges:
            a, b = mesh.vertices[e.vertices[0]].co, mesh.vertices[e.vertices[1]].co
            if abs(a.y) < 1e-5 and abs(b.y) < 1e-5 and a.x > 0 and b.x > 0:
                e.use_seam = True
    payload = Payload(obj)
    assert bpy.ops.uv.ddg_cut_start.poll()
    if regions:
        assert bpy.ops.uv.ddg_cut_start(mode='CLUSTERS', regions=regions) == {'FINISHED'}
    elif circle:
        start_circle(obj, vertex=vertex)
    else:
        bm = bmesh.new()
        bm.from_mesh(mesh)
        for face in bm.faces:
            face.select = face.calc_center_median().x > 0.0
        bm.to_mesh(mesh)
        bm.free()
        assert bpy.ops.uv.ddg_cut_start(mode='SELECTED') == {'FINISHED'}
    names = [a.name for a in mesh.attributes if a.name.startswith("ddg_cut_region_")]
    assert sorted(names) == sorted("ddg_cut_region_{}".format(k) for k in range(regions or 2)), \
        names
    for n in names:
        assert mesh.attributes[n].domain == 'CORNER' and mesh.attributes[n].data_type == 'FLOAT'
    assert len(previews(obj)) == 1
    first_preview = previews(obj)[0]

    assert flow(**({"length_weight": weight} if weight else {})) == {'FINISHED'}
    report = dict(mesh["ddg_cut_last_flow"])
    assert report["steps_taken"] == FLOW["steps"], report
    if regions:
        assert report["n_regions"] == regions and report["n_patches"] >= 3, report
        assert report["n_triple_points"] > 0, report
    assert report["energy_after"] < report["energy_before"], report
    assert len(previews(obj)) == 1 and previews(obj)[0] is not first_preview or True
    assert len(previews(obj)) == 1

    session, reference = rust_reference(obj)
    if inspect is not None:
        inspect(session.describe_cut(min_crossing_offset=CUT_OFFSET))
    seams_before = positions_of_seams(mesh)
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
    assert not [a for a in mesh.attributes if a.name.startswith("ddg_cut_region_")], \
        "the state attributes are still there"
    assert not previews(obj), "the preview is still there"
    assert not [o for o in bpy.data.objects if o.name.endswith(" cut")], "preview object left"
    assert "ddg_cut_last_flow" not in mesh

    check_apply_against_reference(obj, payload, session, reference, seams_before)
    payload.check(obj, reference)
    # Manifold where it was: Euler characteristic unchanged, and no edge got more faces.
    assert euler(mesh) == chi_before, "Euler characteristic changed"
    counts_after = collections.Counter(edge_face_counts(mesh).values())
    assert max(counts_after) == max(counts_before), "an edge has more faces than before"
    if set(counts_before) == {2}:
        assert set(counts_after) == {2}, "the closed mesh is no longer closed"
    return obj


def regions_test(name, n, **kwargs):
    def run():
        workflow(MODELS[name], name, regions=n, **kwargs)
    run.__name__ = "workflow_regions_{}_{}{}".format(name, n, "_seams" if kwargs.get("seams") else "")
    return run


def junction_cases(cut):
    """Checks that the network exercises the cases of the insertion: several junctions in one
    polygon, and a junction next to an ordinary path in one polygon."""
    polygons = collections.Counter(cut.junction_polygon)
    several = [p for p, n in polygons.items() if n >= 2]
    assert several, "no polygon has two junctions (junctions in polygons: {})".format(
        dict(polygons))
    ordinary = set()
    for k, p in enumerate(cut.path_polygon):
        if cut.path_start[k] >= 0 and cut.path_end[k] >= 0:
            ordinary.add(p)
    both = [p for p in polygons if p in ordinary]
    assert both, "no polygon has a junction and an ordinary path"
    print("  {} junctions; {} polygons with several; {} with a junction and an ordinary path"
          .format(len(cut.junction_polygon), len(several), len(both)))


MODELS["ngon_sphere"] = ngon_sphere
for _n in (10,):
    test(regions_test("ngon_sphere", _n, weight=0.01, inspect=junction_cases))
test(regions_test("ngon_sphere", 8, seams=True, weight=0.3))


for _name, _n in (("uv_sphere", 4), ("subdivided_cube", 5), ("torus", 4), ("spot", 6)):
    test(regions_test(_name, _n, weight=0.3 if _name != "torus" else 0.005))
test(regions_test("uv_sphere", 4, seams=True, weight=0.3))


for _name in ("uv_sphere", "subdivided_cube", "torus", "monkey", "spot"):
    def _make(name=_name):
        def run():
            workflow(MODELS[name], name)
        run.__name__ = "workflow_start_flow_apply_{}".format(name)
        return run
    test(_make())


@test
def workflow_with_selected_faces_as_the_start():
    workflow(MODELS["uv_sphere"], "uv_sphere", circle=False)


# ----------------------------------------------------------------------------------------------
# Tests


@test
def operators_are_registered_with_the_agreed_properties():
    for name in ("start", "flow", "apply"):
        cls = getattr(bpy.types, "UV_OT_ddg_cut_" + name)
        assert cls.bl_idname == "uv.ddg_cut_" + name
        assert {'REGISTER', 'UNDO'} <= set(cls.bl_options)
    rna = bpy.ops.uv.ddg_cut_start.get_rna_type()
    properties = {p.identifier: p for p in rna.properties if p.identifier != "rna_type"}
    assert set(properties) == {"mode", "radius", "vertex", "regions", "cluster_iterations",
                               "attribute", "cut_offset", "triangulation", "refine_bound",
                               "optimal_delaunay", "merge_tolerance", "representation"}, \
        sorted(properties)
    assert properties["representation"].default == 'INTEGER'
    assert [i.identifier for i in properties["representation"].enum_items] == [
        'SIGNPOSTS', 'INTEGER']
    assert properties["triangulation"].default == ARGS.triangulation.upper()
    assert properties["refine_bound"].default == 25.0
    assert properties["optimal_delaunay"].default == ARGS.optimal_delaunay
    assert abs(properties["merge_tolerance"].default - 0.1) < 1e-6
    assert [i.identifier for i in properties["triangulation"].enum_items] == [
        'INTRINSIC', 'REFINED', 'CALLER']
    assert properties["mode"].default == 'SELECTED'
    assert properties["regions"].default == 4
    rna = bpy.ops.uv.ddg_cut_flow.get_rna_type()
    properties = {p.identifier: p for p in rna.properties if p.identifier != "rna_type"}
    required = {
        "steps", "use_stop_rule", "energy", "length_mode", "ripple_wavelength", "length_weight",
        "stages", "penalty_group", "penalty_strength", "hide_from_camera", "hide_strength",
        "distortion_weights", "area_weights", "symmetry", "symmetry_tolerance",
        "use_uniform_curvature", "uniform_curvature"}
    assert required <= set(properties), sorted(required - set(properties))
    # The advanced options carry the library's defaults.
    assert abs(properties["time_scale"].default - 0.05) < 1e-6
    assert abs(properties["max_motion"].default - 0.25) < 1e-6
    assert abs(properties["tau_max"].default - 1.0) < 1e-6
    assert abs(properties["min_crossing_offset"].default - 1e-3) < 1e-9
    assert abs(properties["motion_tolerance"].default - 1e-5) < 1e-9
    assert abs(properties["energy_tolerance"].default - 1e-4) < 1e-9
    assert properties["energy_window"].default == 25
    assert abs(properties["feature_smoothing"].default - 0.02) < 1e-6
    assert properties["energy"].default == 'HENCKY'
    assert properties["stages"].default == 'RECOMMENDED'
    assert abs(properties["ripple_wavelength"].default - 0.25) < 1e-6
    assert set(properties["distortion_weights"].default_array) == {1.0}
    assert set(properties["area_weights"].default_array) == {0.0}
    assert properties["distortion_weights"].array_length == 16
    rna = bpy.ops.uv.ddg_cut_apply.get_rna_type()
    properties = {p.identifier for p in rna.properties if p.identifier != "rna_type"}
    assert properties == {"unwrap_after", "cut_offset", "merge_tolerance"}, properties
    for menu in (bpy.types.VIEW3D_MT_uv_map, bpy.types.VIEW3D_MT_object):
        names = [f.__name__ for f in menu.draw._draw_funcs]
        assert "_draw_menu_entry" in names, names


@test
def start_attributes_and_preview_follow_the_object():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    obj.location = (1.0, 2.0, 3.0)
    obj.scale = (2.0, 1.0, 0.5)
    start_circle(obj, vertex=5)
    rows = state(mesh)
    assert len(rows) == 2 and all(len(r) == len(mesh.loops) for r in rows)
    assert all(abs(a + b) < 1e-7 for a, b in zip(rows[0], rows[1])), "region 1 is not -region 0"
    corner_vertices = array.array('i', bytes(4 * len(mesh.loops)))
    mesh.loops.foreach_get("vertex_index", corner_vertices)
    centre = list(corner_vertices).index(5)
    assert rows[0][centre] < 0.0 < rows[1][centre]
    preview, = previews(obj)
    assert preview.name.startswith(obj.name), preview.name
    assert preview.parent is obj and preview.type == 'MESH'
    assert len(preview.data.polygons) == 0 and len(preview.data.edges) > 0
    assert len(preview.data.vertices) == len(preview.data.edges)       # one closed loop
    assert preview.display_type == 'WIRE'
    bpy.context.view_layer.update()
    assert (preview.matrix_world.translation - obj.matrix_world.translation).length < 1e-6
    # A preview vertex is on the surface of the object: its world position is the object's
    # transform of a point at (about) the radius from the centre vertex.
    pts = positions(mesh)
    inverse = obj.matrix_world.inverted()
    for v in list(preview.data.vertices)[:5]:
        local = inverse @ (preview.matrix_world @ v.co)
        assert abs(local.length - 1.0) < 0.05, "preview point {} is off the sphere".format(local)
    # Replaced on every run.
    old_name = preview.name
    assert flow(steps=2) == {'FINISHED'}
    assert len(previews(obj)) == 1
    assert len([o for o in bpy.data.objects if o.name.startswith(old_name[:6])]) == 2
    assert len([m for m in bpy.data.meshes if m.name.endswith(" cut") or " cut." in m.name]) == 1
    del pts


@test
def start_from_material_slots_and_from_an_integer_attribute():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    for k in range(3):
        mesh.materials.append(bpy.data.materials.new("m{}".format(k)))
    labels = [0 if p.center.z > 0.3 else (1 if p.center.x > 0 else 2) for p in mesh.polygons]
    mesh.polygons.foreach_set("material_index", labels)
    assert bpy.ops.uv.ddg_cut_start(mode='MATERIALS') == {'FINISHED'}
    assert len(state(mesh)) == 3
    assert len(previews(obj)) == 1 and len(previews(obj)[0].data.edges) > 0
    # Labels from an attribute with arbitrary integer values: ascending order gives the regions.
    mesh.attributes.remove(mesh.attributes["ddg_cut_region_0"])
    for name in ("ddg_cut_region_1", "ddg_cut_region_2"):
        mesh.attributes.remove(mesh.attributes[name])
    attribute = mesh.attributes.new("label", 'INT', 'FACE')
    attribute.data.foreach_set("value", [(-5 if p.center.z > 0 else 7 if p.center.x > 0 else 40)
                                         for p in mesh.polygons])
    assert bpy.ops.uv.ddg_cut_start(mode='ATTRIBUTE', attribute="label") == {'FINISHED'}
    rows = state(mesh)
    assert len(rows) == 3
    corner_vertices = array.array('i', bytes(4 * len(mesh.loops)))
    mesh.loops.foreach_get("vertex_index", corner_vertices)
    top = [c for c in range(len(mesh.loops)) if mesh.vertices[corner_vertices[c]].co.z > 0.9]
    assert top and all(min(range(3), key=lambda k: rows[k][c]) == 0 for c in top), \
        "the lowest value (-5) is region 0"
    # The preview of a network: curves between junctions, edges only.
    preview, = previews(obj)
    assert len(preview.data.polygons) == 0 and len(preview.data.edges) > 0
    # Flow works on it, and reports the patches and triple points.
    assert flow(steps=2) == {'FINISHED'}
    report = dict(mesh["ddg_cut_last_flow"])
    assert report["n_regions"] == 3 and report["n_triple_points"] >= 1, report
    # A float attribute or a missing one is refused.
    mesh.attributes.new("f", 'FLOAT', 'FACE')
    expect_error("must be an integer attribute on the face domain", bpy.ops.uv.ddg_cut_start,
                 mode='ATTRIBUTE', attribute="f")
    expect_error("at most 16", bpy.ops.uv.ddg_cut_start, mode='CLUSTERS', regions=17) \
        if False else None


@test
def start_circle_uses_the_active_vertex_in_edit_mode_and_the_default_radius():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    bm.verts.ensure_lookup_table()
    for v in bm.verts:
        v.select = False
    bm.verts[5].select = True
    bm.select_history.clear()
    bm.select_history.add(bm.verts[5])
    bmesh.update_edit_mesh(obj.data)
    assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=0.0) == {'FINISHED'}
    assert obj.mode == 'EDIT', "Start did not return to Edit mode"
    bpy.ops.object.mode_set(mode='OBJECT')
    mesh = obj.data
    corner_vertices = array.array('i', bytes(4 * len(mesh.loops)))
    mesh.loops.foreach_get("vertex_index", corner_vertices)
    rows = state(mesh)
    inside = {corner_vertices[c] for c in range(len(rows[0])) if rows[0][c] < 0}
    assert 5 in inside and 0 < len(inside) < len(mesh.vertices)
    # Region 0 is centred on vertex 5: the farthest inside vertex is not far from it.
    centre = positions(mesh)[5]
    assert max((positions(mesh)[v] - centre).length for v in inside) < 0.5 * diagonal(mesh)


@test
def flow_in_edit_mode_stays_in_edit_mode():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    start_circle(obj, vertex=5)
    bpy.ops.object.mode_set(mode='EDIT')
    assert flow(steps=3) == {'FINISHED'}
    assert obj.mode == 'EDIT'
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
    assert obj.mode == 'EDIT'
    bpy.ops.object.mode_set(mode='OBJECT')
    assert not state(obj.data) and sum(e.use_seam for e in obj.data.edges) > 0


@test
def k_steps_twice_equal_2k_steps_once():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    start_circle(obj, vertex=5)
    initial = save_state(mesh)
    assert flow(steps=8) == {'FINISHED'}
    once = state(mesh)
    # Back to the start: write the stored state again.
    load_state(mesh, initial)
    assert flow(steps=4) == {'FINISHED'}
    assert flow(steps=4) == {'FINISHED'}
    twice = state(mesh)
    scale = diagonal(mesh)
    worst = max(abs(a - b) for r, s in zip(once, twice) for a, b in zip(r, s))
    print("  8 steps once vs 4 + 4: largest difference {:.3g} (size {:.3g})".format(worst, scale))
    assert worst < 1e-4 * scale, worst
    # The flow ran: the state is not the start.
    assert max(abs(a - b) for a, b in zip(once[0], initial[0][0])) > 1e-3 * scale


@test
def flow_reports_energy_length_steps_and_stop_reason():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    start_circle(obj, vertex=5)
    with CaptureOutput() as out:
        assert flow(steps=5) == {'FINISHED'}
    report = dict(obj.data["ddg_cut_last_flow"])
    for key in ("steps_taken", "stop_reason", "energy_before", "energy_after",
                "cut_length_before", "cut_length_after", "length_weight"):
        assert key in report, key
    assert "Cut flow: 5 steps" in out.text and "energy" in out.text and "cut length" in out.text, \
        out.text
    # With the stop rule on and many steps, the flow stops by itself with a reason.
    with CaptureOutput() as out:
        assert flow(steps=2000, use_stop_rule=True) == {'FINISHED'}
    report = dict(obj.data["ddg_cut_last_flow"])
    assert report["steps_taken"] < 2000 and report["stop_reason"] != "max_steps", report


@test
def ripple_wavelength_sets_the_length_weight():
    import ddg_tools
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    start_circle(obj, vertex=5)
    assert bpy.ops.uv.ddg_cut_flow(steps=1, use_stop_rule=False, length_mode='RIPPLE',
                                   ripple_wavelength=0.3) == {'FINISHED'}
    report = dict(obj.data["ddg_cut_last_flow"])
    assert report["length_weight"] > 0.0
    omega = report.get("length_weight_roughening_rate")
    assert abs(report["length_weight"] - ddg_tools.length_weight_for_ripples(omega, 0.3)) \
        < 1e-6 * report["length_weight"]


@test
def every_flow_property_reaches_the_library():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    start_circle(obj, vertex=0)
    mesh = obj.data
    initial = save_state(mesh)

    def fresh():
        load_state(mesh, initial)

    def final_state(**props):
        fresh()
        assert flow(steps=3, **props) == {'FINISHED'}, props
        report = dict(mesh["ddg_cut_last_flow"])
        assert math.isfinite(report["energy_after"]), props
        return array.array('f', state(mesh)[0]), report

    base, base_report = final_state()

    def differs(props):
        other, _ = final_state(**props)
        return max(abs(a - b) for a, b in zip(base, other)) > 1e-6

    group = obj.vertex_groups.new(name="Penalty")
    for v in mesh.vertices:
        if v.co.x > 0.0:
            group.add([v.index], 1.0, 'REPLACE')
    cases = collections.OrderedDict([
        ("energy dirichlet", dict(energy='DIRICHLET')),
        ("energy scaled dirichlet", dict(energy='SCALED_DIRICHLET')),
        ("feature size cap", dict(energy='SCALED_DIRICHLET', feature_r_eps=0.2)),
        ("paper stages", dict(stages='PAPER')),
        ("explicit length weight", dict(length_mode='EXPLICIT', length_weight=0.5)),
        ("penalty group", dict(penalty_group="Penalty", penalty_strength=3.0)),
        ("penalty paper rule", dict(penalty_group="Penalty", penalty_strength=3.0,
                                    penalty_rule='PAPER_DIAGONAL')),
        ("distortion weights", dict(distortion_weights=pad((0.0, 1.0), 1.0))),
        ("area weights", dict(area_weights=pad((0.5, 0.0), 0.0))),
        ("uniform curvature", dict(use_uniform_curvature=True, uniform_curvature=5.0)),
        ("step rule fixed", dict(step_rule='FIXED', tau=1e-4)),
        ("step rule capped", dict(step_rule='CAPPED', tau_max=1e-3, max_motion=0.05)),
        ("flux nodal", dict(flux_recovery='NODAL')),
        ("length step on mesh", dict(length_step='ON_MESH')),
        ("extension harmonic", dict(extension='HARMONIC')),
        ("speed unit", dict(speed_scaling='UNIT')),
        ("redistance off", dict(redistance='OFF')),
    ])
    # The default length mode of `base` is explicit with weight 0.1; "explicit length weight" has
    # another value.
    for label, props in cases.items():
        assert differs(props), "{} has no effect".format(label)
    # A case that must not change anything: the symmetry options with the axis off.
    assert not differs(dict(symmetry='NONE', symmetry_tolerance=0.5))
    # The default stage choices equal the explicit defaults.
    assert not differs(dict(flux_recovery='DEFAULT', length_step='DEFAULT'))
    # Symmetric start (circle at the pole) with mirror symmetry about every axis runs.
    for axis in ("X", "Y", "Z"):
        final_state(symmetry=axis)
    # Camera penalty.
    camera = bpy.data.objects.new("Camera", bpy.data.cameras.new("Camera"))
    bpy.context.collection.objects.link(camera)
    camera.location = (5.0, 0.0, 0.0)
    camera.rotation_euler = (math.pi / 2, 0.0, math.pi / 2)
    bpy.context.scene.camera = camera
    assert differs(dict(hide_from_camera=True, hide_strength=3.0)), "camera penalty: no effect"
    bpy.context.scene.camera = None
    expect_error("no active camera", bpy.ops.uv.ddg_cut_flow, hide_from_camera=True, **FLOW)


@test
def symmetry_matches_mirrored_vertices_and_reports_the_unmatched():
    core = importlib.import_module(ARGS.module + ".cut")
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    classes, unmatched = core.mirror_classes(mesh, "X", 1e-5)
    assert not unmatched, unmatched
    pts = positions(mesh)
    groups = collections.defaultdict(list)
    for i, c in enumerate(classes):
        groups[c].append(i)
    assert all(len(g) in (1, 2) for g in groups.values())
    for g in groups.values():
        if len(g) == 2:
            a, b = pts[g[0]], pts[g[1]]
            assert abs(a.x + b.x) < 1e-5 and abs(a.y - b.y) < 1e-5 and abs(a.z - b.z) < 1e-5
        else:
            assert abs(pts[g[0]].x) < 1e-5, "a class of one vertex that is off the plane"
    # Move one vertex: it and its partner are unmatched.
    mesh.vertices[40].co.x += 0.05
    classes, unmatched = core.mirror_classes(mesh, "X", 1e-5)
    assert 40 in unmatched and len(unmatched) >= 1
    start_circle(obj, vertex=0)
    with CaptureOutput() as out:
        assert flow(steps=1, symmetry='X', symmetry_tolerance=1e-5 / diagonal(mesh)) \
            == {'FINISHED'}
    assert "no mirror image about the X axis" in out.text and " 40" in out.text, out.text


@test
def visibility_rays_see_occluders():
    core = importlib.import_module(ARGS.module + ".cut")
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    wall = primitive("plane", size=4.0, location=(0.0, 0.0, 3.0))
    activate(obj)
    bpy.context.view_layer.update()
    camera_at = Vector((0.0, 0.0, 6.0))
    visible = core.vertex_visibility(obj, bpy.context, camera_at)
    top = max(range(len(obj.data.vertices)), key=lambda i: obj.data.vertices[i].co.z)
    assert visible[top] == 0.0, "the wall does not hide the top of the sphere"
    assert sum(visible) == 0.0
    bpy.data.objects.remove(wall, do_unlink=True)
    bpy.context.view_layer.update()
    visible = core.vertex_visibility(obj, bpy.context, camera_at)
    assert visible[top] == 1.0
    south = min(range(len(obj.data.vertices)), key=lambda i: obj.data.vertices[i].co.z)
    assert visible[south] == 0.0, "the sphere hides its own far side"
    assert 0 < sum(visible) < len(visible)


@test
def one_sided_patch_weights_change_the_cut():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    start_circle(obj, vertex=0)
    initial = save_state(mesh)
    results = []
    for weights in ((1.0, 1.0), (0.0, 1.0)):
        load_state(mesh, initial)
        assert flow(steps=10, distortion_weights=pad(weights, 1.0)) == {'FINISHED'}
        results.append(dict(mesh["ddg_cut_last_flow"]))
    assert abs(results[0]["energy_after"] - results[1]["energy_after"]) > 1e-6


@test
def apply_with_unwrap_gives_uvs_without_flipped_triangles():
    for name in ("uv_sphere", "subdivided_cube"):
        obj = MODELS[name]()
        start_circle(obj, fraction=0.3)
        assert flow(steps=10) == {'FINISHED'}
        assert bpy.ops.uv.ddg_cut_apply(unwrap_after=True) == {'FINISHED'}
        mesh = obj.data
        layer = mesh.uv_layers.active
        assert layer is not None
        uvs = array.array('f', bytes(8 * len(mesh.loops)))
        mesh.attributes[layer.name].data.foreach_get("vector", uvs)
        assert all(math.isfinite(x) for x in uvs)
        mesh.calc_loop_triangles()
        areas = []
        for tri in mesh.loop_triangles:
            (u0, v0), (u1, v1), (u2, v2) = (uvs[2 * c:2 * c + 2] for c in tri.loops)
            areas.append(0.5 * ((u1 - u0) * (v2 - v0) - (u2 - u0) * (v1 - v0)))
        mean = sum(abs(a) for a in areas) / len(areas)
        # Triangles of three collinear points (a polygon with a vertex of the cut on its side)
        # have zero area up to rounding, whichever way they are "oriented".
        flipped = sum(1 for a in areas if a < -1e-3 * mean)
        assert flipped == 0, "{}: {} flipped triangles of {}".format(name, flipped, len(areas))
        print("  {}: {} triangles, {} of zero area (collinear points)".format(
            name, len(areas), sum(1 for a in areas if abs(a) <= 1e-3 * mean)))
        # The unwrap used the cut as seams: charts are separated along them.
        assert sum(e.use_seam for e in mesh.edges) > 0
        reset()


@test
def second_cut_works_with_the_first_as_fixed_seams():
    for name in ("uv_sphere", "subdivided_cube", "torus"):
        obj = MODELS[name]()
        mesh = obj.data
        payload = Payload(obj)
        chi = euler(mesh)
        start_circle(obj, fraction=0.2, vertex=0)
        assert flow() == {'FINISHED'}
        assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
        first_seams = positions_of_seams(mesh)
        assert first_seams
        payload.refresh(obj)
        # A second cut, started at the antipode of the first, crosses the first one's seams
        # when its region reaches them.
        far = max(range(len(mesh.vertices)),
                  key=lambda i: (mesh.vertices[i].co - mesh.vertices[0].co).length)
        assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=0.25 * diagonal(mesh),
                                        vertex=far) == {'FINISHED'}
        assert flow() == {'FINISHED'}
        session, reference = rust_reference(obj)
        assert len(session.seams_cut) // 2 >= len(first_seams) - 1
        seams_before = positions_of_seams(mesh)
        assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
        check_apply_against_reference(obj, payload, session, reference, seams_before)
        payload.check(obj, reference)
        assert euler(mesh) == chi
        now = {frozenset((round(a[0], 3), round(a[1], 3), round(a[2], 3)) for a in pair)
               for pair in positions_of_seams(mesh)}
        # Every old seam is still a seam, perhaps cut into pieces: each old seam's end points
        # are vertices of seam edges.
        seam_vertices = {tuple(round(c, 3) for c in p) for pair in positions_of_seams(mesh)
                         for p in pair}
        for a, b in first_seams:
            assert tuple(round(c, 3) for c in a) in seam_vertices
            assert tuple(round(c, 3) for c in b) in seam_vertices
        assert len(now) > len(first_seams)
        reset()


@test
def apply_marks_exactly_the_inserted_edges_and_leaves_old_seams_alone():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    # An existing seam along a meridian: the cut cannot cross it, and it stays.
    for e in mesh.edges:
        a, b = mesh.vertices[e.vertices[0]].co, mesh.vertices[e.vertices[1]].co
        if abs(a.y) < 1e-6 and abs(b.y) < 1e-6 and a.x > 0 and b.x > 0:
            e.use_seam = True
    n_old = sum(e.use_seam for e in mesh.edges)
    assert n_old >= 10
    start_circle(obj, fraction=0.3, vertex=5)
    assert flow() == {'FINISHED'}
    session, reference = rust_reference(obj)
    seams_before = positions_of_seams(mesh)
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
    check_apply_against_reference(obj, None, session, reference, seams_before)
    now = sum(e.use_seam for e in mesh.edges)
    assert now >= n_old + len(reference.new_seams) // 2 - 0


@test
def models_with_boundary_keep_their_boundary():
    obj = primitive("grid", x_subdivisions=12, y_subdivisions=12, size=2.0)
    mesh = obj.data
    boundary_before = sum(1 for c in edge_face_counts(mesh).values() if c == 1)
    assert boundary_before > 0
    start_circle(obj, fraction=0.2, vertex=70)
    assert flow(steps=3) == {'FINISHED'}
    session, reference = rust_reference(obj)
    seams_before = positions_of_seams(mesh)
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
    check_apply_against_reference(obj, None, session, reference, seams_before)
    counts = collections.Counter(edge_face_counts(mesh).values())
    assert set(counts) <= {1, 2}
    assert euler(mesh) == 1


FLIP_RECORD = []


def atlas_of_spot(weight, expect_zero_flips, lam=0.1, offset=CUT_OFFSET):
    core = importlib.import_module(ARGS.module + ".core")
    obj = spot_object()
    mesh = obj.data
    assert bpy.ops.uv.ddg_cut_start(mode='CLUSTERS', regions=8, cluster_iterations=50,
                                    cut_offset=offset) == {'FINISHED'}
    assert bpy.ops.uv.ddg_cut_flow(steps=60, use_stop_rule=True, length_mode='EXPLICIT',
                                   length_weight=weight, stages='RECOMMENDED',
                                   cut_offset=offset) == {'FINISHED'}
    report = dict(mesh["ddg_cut_last_flow"])
    print("  flow: {} steps ({}), {} patches, {} triple points, cut length {:.3f}".format(
        report["steps_taken"], report["stop_reason"], report["n_patches"],
        report["n_triple_points"], report["cut_length_after"]))
    assert report["n_patches"] >= 3 and report["n_triple_points"] > 0
    assert bpy.ops.uv.ddg_cut_apply(cut_offset=offset) == {'FINISHED'}
    assert mesh.get("ddg_cut_triangulation") is not None, "Apply stored no triangulation"
    # The unwrap operator's work, called directly to get the library's report.
    stats = core.unwrap_object(obj, use_seams=True, mark_cuts_as_seams=False, lambda_=lam)
    note = stats["triangulation"]
    assert note and note["used"] > 0 and note["stale"] == 0, note
    print("  library: flipped {}, collapsed {}".format(stats["report"].get("n_flipped_uv_triangles"), stats["report"].get("n_collapsed_uv_triangles")))
    charts = stats["report"]["n_charts"]
    layer = mesh.uv_layers.active
    uvs = array.array('f', bytes(8 * len(mesh.loops)))
    mesh.attributes[layer.name].data.foreach_get("vector", uvs)
    islands = atlas_islands(mesh, uvs)
    print("  library charts {}, UV islands {}, cones {}".format(charts, islands, stats["n_cones"]))
    assert islands == charts, "{} UV islands, the library reports {} charts".format(islands, charts)
    # Our own count, on the triangulation the unwrap used (the stored one for the split
    # polygons, Blender's for the others).
    _, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
    triangles, _ = core.stored_triangulation(mesh, corner_vertices, polygon_starts, triangles)
    areas = []
    for k in range(0, len(triangles), 3):
        (u0, v0), (u1, v1), (u2, v2) = (uvs[2 * c:2 * c + 2] for c in triangles[k:k + 3])
        areas.append(0.5 * ((u1 - u0) * (v2 - v0) - (u2 - u0) * (v1 - v0)))
    mean = sum(abs(a) for a in areas) / len(areas)
    flipped = sum(1 for a in areas if a < -1e-3 * mean)
    print("  worst negative area / mean: {:.3g}".format(min(areas) / mean))
    library_flipped = stats["report"].get("n_flipped_uv_triangles", 0)
    # The flips are counted twice, independently: here from the UVs and by the library on its
    # own triangulation. The triangulations differ, so the counts are not equal in general.
    print("  flipped: {} here (of {} triangles), {} reported by the library".format(
        flipped, len(areas), library_flipped))
    # Nothing flips on any triangulation, Optimal Delaunay included, at either offset (the
    # stored triangulation of Apply; on an intrinsic one the merge of bends measured at the cut
    # the flow found, and the excursions across a side removed: crates/ddg-uv/README.md, "A large
    # cut offset, and the optimal Delaunay triangulation").
    FLIP_RECORD.append((ARGS.triangulation, weight, offset, flipped, library_flipped))
    if expect_zero_flips:
        assert flipped == 0 and library_flipped == 0, "{} flipped triangles of {}".format(
            flipped, len(areas))
    return charts, islands, flipped




@test
def atlas_of_spot_eight_clusters_gives_as_many_uv_islands_as_charts_and_no_flipped_faces():
    # The library's triangulation of the split polygons is stored by Apply and used by the
    # unwrap: no flipped triangle by the library's count or ours, at the default offset (1e-3)
    # and at a large one (0.2). With Blender's own triangulation the same atlas has four at 1e-3.
    for offset in (CUT_OFFSET, 0.2):
        reset()
        charts, islands, flipped = atlas_of_spot(2.0, expect_zero_flips=True, offset=offset)
        assert charts >= 5


@test
def atlas_of_spot_with_a_strong_length_weight_has_no_flipped_faces():
    for offset in (CUT_OFFSET, 0.2):
        reset()
        charts, islands, flipped = atlas_of_spot(10.0, expect_zero_flips=True, offset=offset)
        assert charts >= 3


@test
def blender_triangulation_of_the_split_polygons_would_flip_triangles():
    # The evidence for the stored triangulation: the same atlas unwrapped with Blender's own
    # triangulation of the split polygons has flipped triangles at 1e-3 (on the caller's
    # triangulation) or is refused by the library (polygons with collinear vertices on an
    # intrinsic cut: "a diagonal of a face coincides with another edge").
    core = importlib.import_module(ARGS.module + ".core")
    obj = spot_object()
    mesh = obj.data
    assert bpy.ops.uv.ddg_cut_start(mode='CLUSTERS', regions=8, cluster_iterations=50) \
        == {'FINISHED'}
    assert bpy.ops.uv.ddg_cut_flow(steps=60, use_stop_rule=True, length_mode='EXPLICIT',
                                   length_weight=2.0) == {'FINISHED'}
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
    stats = core.unwrap_object(obj, use_seams=True, use_stored_triangulation=True)
    stored = stats["report"].get("n_flipped_uv_triangles")
    try:
        stats = core.unwrap_object(obj, use_seams=True, use_stored_triangulation=False)
        blender = "{} flipped".format(stats["report"].get("n_flipped_uv_triangles"))
    except core.UnwrapFailed as ex:
        blender = "refused: " + str(ex)[:80]
    print("  stored triangulation: {} flipped; Blender's: {}".format(stored, blender))
    if ARGS.triangulation == "caller":
        assert stored == 0
        if blender == "0 flipped":
            raise Skip("this atlas no longer flips with Blender's triangulation")


@test
def unwrap_falls_back_to_blenders_triangulation_when_the_mesh_was_edited():
    core = importlib.import_module(ARGS.module + ".core")
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    # A mesh that never had a cut applied has no stored triangulation: nothing changes.
    stats = core.unwrap_object(obj)
    assert stats["triangulation"] is None
    reference = array.array('f', bytes(8 * len(mesh.loops)))
    mesh.attributes[mesh.uv_layers.active.name].data.foreach_get("vector", reference)
    start_circle(obj, fraction=0.3, vertex=5)
    assert flow() == {'FINISHED'}
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
    assert mesh.get("ddg_cut_triangulation") is not None
    with CaptureOutput() as out:
        assert bpy.ops.uv.ddg_conformal_unwrap() == {'FINISHED'}
    assert "library triangulation for" in out.text and "Used Blender's triangulation" not in out.text
    stats = core.unwrap_object(obj)
    assert stats["triangulation"]["used"] > 0 and stats["triangulation"]["stale"] == 0
    # Moving a vertex keeps the stored triangulation valid (same faces, same vertices).
    mesh.vertices[3].co.x += 0.001
    stats = core.unwrap_object(obj)
    assert stats["triangulation"]["used"] > 0
    # Deleting a face changes the number of faces: the stored triangulation is not used, the
    # operator says so, and the unwrap still works.
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bm.faces.ensure_lookup_table()
    bmesh.ops.delete(bm, geom=[bm.faces[0]], context='FACES_ONLY')
    bm.to_mesh(mesh)
    bm.free()
    mesh.update()
    with CaptureOutput() as out:
        assert bpy.ops.uv.ddg_conformal_unwrap() == {'FINISHED'}
    assert "Used Blender's triangulation for" in out.text and "edited since" in out.text, out.text
    assert "library triangulation for" not in out.text
    del reference


@test
def the_cut_offset_moves_points_away_from_vertices():
    """The smallest triangle against the mean after Apply, at the cut offsets 1e-3 and 0.2, on
    the triangulation the unwrap uses: the library's, stored by Apply for the split polygons.
    Blender's own triangulation of those polygons is printed too, not asserted: it does not
    follow the cut and can make a thin triangle at any offset (with Optimal Delaunay on the
    refined triangulation it gave a ratio of 13; the stored one gives 169, and 351 on the
    intrinsic one)."""
    core = importlib.import_module(ARGS.module + ".core")
    smallest, blender = {}, {}
    for offset in (1e-3, 0.2):
        reset()
        obj = spot_object()
        mesh = obj.data
        assert bpy.ops.uv.ddg_cut_start(mode='CLUSTERS', regions=8, cluster_iterations=50,
                                        cut_offset=offset) == {'FINISHED'}
        assert bpy.ops.uv.ddg_cut_flow(steps=20, length_mode='EXPLICIT', length_weight=2.0,
                                       cut_offset=offset) == {'FINISHED'}
        assert bpy.ops.uv.ddg_cut_apply(cut_offset=offset) == {'FINISHED'}
        positions, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
        stored, note = core.stored_triangulation(mesh, corner_vertices, polygon_starts, triangles)
        assert note and note["used"] > 0 and note["stale"] == 0, note

        def ratio(corners):
            areas = []
            for k in range(0, len(corners), 3):
                a, b, c = (Vector(positions[3 * v:3 * v + 3])
                           for v in (corner_vertices[i] for i in corners[k:k + 3]))
                areas.append(0.5 * (b - a).cross(c - a).length)
            return min(areas) / (sum(areas) / len(areas))

        smallest[offset] = ratio(stored)
        blender[offset] = ratio(triangles)
    print("  smallest triangle area / mean: {:.3g} at 1e-3, {:.3g} at 0.2 (Blender's own "
          "triangulation of the split polygons: {:.3g} and {:.3g})".format(
              smallest[1e-3], smallest[0.2], blender[1e-3], blender[0.2]))
    assert smallest[0.2] > 20 * smallest[1e-3]


# ----------------------------------------------------------------------------------------------
# The two bad meshes of crates/ddg-variational-cuts/examples/bad_mesh_study.rs


def pinch_free(faces):
    """Removes faces until every vertex has one fan of faces (no pinched vertices)."""
    while True:
        by_vertex = collections.defaultdict(list)
        for f, face in enumerate(faces):
            for v in face:
                by_vertex[v].append(f)
        drop = set()
        for v, around in by_vertex.items():
            # Faces around v connected through edges (v, w).
            edge_faces = collections.defaultdict(list)
            for f in around:
                for w in faces[f]:
                    if w != v:
                        edge_faces[w].append(f)
            parent = {f: f for f in around}

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            for fs in edge_faces.values():
                for g in fs[1:]:
                    parent[find(g)] = find(fs[0])
            groups = collections.defaultdict(list)
            for f in around:
                groups[find(f)].append(f)
            if len(groups) > 1:
                keep = max(groups.values(), key=len)
                for g in groups.values():
                    if g is not keep:
                        drop.update(g)
        if not drop:
            return faces
        faces = [f for i, f in enumerate(faces) if i not in drop]


def compacted(positions, faces):
    used = sorted({v for f in faces for v in f})
    index = {v: i for i, v in enumerate(used)}
    return [positions[v] for v in used], [tuple(index[v] for v in f) for f in faces]


def sheared_grid(n, sigma):
    """The lattice (i + sigma j, j)/n clipped to the unit square, split into triangles of largest
    angle 180 - atan(1/sigma) degrees (the study's `sheared`)."""
    h = 1.0 / n
    index, positions = {}, []
    for j in range(n + 1):
        lo = math.floor(-sigma * j) - 1
        hi = math.ceil(n - sigma * j) + 1
        for i in range(lo, hi + 1):
            x = (i + sigma * j) * h
            if -1e-9 <= x <= 1.0 + 1e-9:
                index[(i, j)] = len(positions)
                positions.append((x, j * h, 0.0))
    faces = []
    for (i, j), a in sorted(index.items()):
        b, c, d = index.get((i + 1, j)), index.get((i + 1, j + 1)), index.get((i, j + 1))
        if b is not None and c is not None and d is not None:
            faces += [(a, b, c), (a, c, d)]
    return compacted(positions, pinch_free(faces))


def needles_and_caps_grid(n):
    """Rows alternate between thin strips (caps: one angle of 168 degrees) and tall ones
    (needles: apex 19 degrees), shifted by half a cell (the study's `needles_and_caps`)."""
    d = 1.0 / n
    ha, hb = 0.05 * d, 3.0 * d
    pairs = max(1, round(1.0 / (ha + hb)))
    scale = 1.0 / (pairs * (ha + hb))
    rows = 2 * pairs + 1
    y = [0.0]
    for j in range(2 * pairs):
        y.append(y[j] + scale * (ha if j % 2 == 0 else hb))
    y[-1] = 1.0
    positions, row_x, start = [], [], []
    for j in range(rows):
        xs = [i * d for i in range(n + 1)] if j % 2 == 0 else (
            [0.0] + [(i + 0.5) * d for i in range(n)] + [1.0])
        start.append(len(positions))
        positions += [(x, y[j], 0.0) for x in xs]
        row_x.append(xs)
    faces = []
    for j in range(rows - 1):
        xa, xb = row_x[j], row_x[j + 1]
        i = m = 0
        while i < len(xa) - 1 or m < len(xb) - 1:
            if i == len(xa) - 1:
                bottom = False
            elif m == len(xb) - 1:
                bottom = True
            else:
                bottom = xa[i + 1] <= xb[m + 1]
            a, b = start[j] + i, start[j + 1] + m
            if bottom:
                faces.append((a, a + 1, b))
                i += 1
            else:
                faces.append((a, b + 1, b))
                m += 1
    return positions, faces


def mesh_object(name, positions, faces):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(positions, [], faces)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    return activate(obj)


def bad_mesh_workflow(label, positions, faces):
    obj = mesh_object(label, positions, faces)
    mesh = obj.data
    chi = euler(mesh)
    # A disk around the vertex nearest to the centre of the square.
    centre = Vector((0.511, 0.493, 0.0))
    vertex = min(range(len(mesh.vertices)), key=lambda i: (mesh.vertices[i].co - centre).length)
    assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=0.25, vertex=vertex) == {'FINISHED'}
    assert flow(steps=20, length_mode='RIPPLE', ripple_wavelength=0.25, use_stop_rule=False,
                use_uniform_curvature=True, uniform_curvature=5.0) \
        == {'FINISHED'}
    report = dict(mesh["ddg_cut_last_flow"])
    alive = report["stop_reason"] != "cut_vanished" and report["cut_length_after"] > 0.0
    print("  {} [{}]: {} steps, {}, energy {:.4g} -> {:.4g}, cut length {:.3f} (was {:.3f})".format(
        label, ARGS.triangulation, report["steps_taken"], report["stop_reason"],
        report["energy_before"], report["energy_after"], report["cut_length_after"],
        report["cut_length_before"]))
    BAD_RECORD.append((label, ARGS.triangulation, alive, report["energy_before"],
                       report["energy_after"], report["cut_length_after"]))
    if ARGS.triangulation != "caller":
        assert alive, "the cut vanished on {}".format(label)
        assert report["energy_after"] < report["energy_before"], report
    if not alive:
        return              # nothing to apply (caller's triangulation on a bad mesh)
    session, reference = rust_reference(obj)
    seams_before = positions_of_seams(mesh)
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
    check_apply_against_reference(obj, None, session, reference, seams_before)
    assert euler(mesh) == chi
    core = importlib.import_module(ARGS.module + ".core")
    stats = core.unwrap_object(obj, use_seams=True)
    note = stats["triangulation"]
    print("    unwrap: {} charts, {} flipped, stored triangulation for {} faces".format(
        stats["report"].get("n_charts"), stats["report"].get("n_flipped_uv_triangles"),
        note["used"] if note else 0))
    BAD_RECORD.append((label + " unwrap", ARGS.triangulation, stats["report"].get("n_charts"),
                       stats["report"].get("n_flipped_uv_triangles"), 0, 0))


BAD_RECORD = []


@test
def bad_mesh_sheared_153_degrees_start_flow_apply_unwrap():
    positions, faces = sheared_grid(32, 2.0)
    bad_mesh_workflow("sheared 153", positions, faces)


@test
def bad_mesh_needles_and_caps_start_flow_apply_unwrap():
    positions, faces = needles_and_caps_grid(32)
    bad_mesh_workflow("needles and caps", positions, faces)


@test
def bad_mesh_sheared_165_degrees_start_flow_apply_unwrap():
    positions, faces = sheared_grid(32, 1.0 / math.tan(math.radians(15.0)))
    bad_mesh_workflow("sheared 165", positions, faces)


@test
def refined_session_keeps_its_state_in_bytes_and_notices_a_stale_one():
    if ARGS.triangulation == "caller":
        raise Skip("the caller's triangulation has no state bytes")
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=0.9, vertex=5,
                                    triangulation='REFINED', refine_bound=30.0,
                                    optimal_delaunay=False) == {'FINISHED'}
    assert mesh["ddg_cut_option"] == "refined:30"
    assert "ddg_cut_state" in mesh, "refinement adds vertices: the state needs its bytes"
    assert flow(steps=3, show_advanced=True) == {'FINISHED'}
    report = dict(mesh["ddg_cut_last_flow"])
    assert report["vertices_added"] > 0 and report["triangulation"] == "refined:30"
    # The attributes stay the view of the mesh's own vertices.
    assert len(state(mesh)) == 2
    before = bytes(mesh["ddg_cut_state"], "ascii")
    # Seams changed since Start: the stored state is stale, the message says so.
    mesh.edges[0].use_seam = True
    message = expect_error("changed since 'Cut: Start'.*Run 'Cut: Start' again",
                           bpy.ops.uv.ddg_cut_flow, **FLOW)
    expect_error("changed since 'Cut: Start'", bpy.ops.uv.ddg_cut_apply)
    assert bytes(mesh["ddg_cut_state"], "ascii") == before, "a stale state was overwritten"
    mesh.edges[0].use_seam = False
    assert flow(steps=2) == {'FINISHED'}      # valid again
    # A new Start replaces everything; Apply removes the bytes and the option.
    assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=0.9, vertex=5,
                                    triangulation='INTRINSIC',
                                    optimal_delaunay=False) == {'FINISHED'}
    assert "ddg_cut_state" not in mesh and mesh["ddg_cut_option"] == "intrinsic"
    assert flow(steps=2) == {'FINISHED'}
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
    for key in ("ddg_cut_state", "ddg_cut_option", "ddg_cut_last_flow"):
        assert key not in mesh, key
    del message


def flipped_uv_triangles(mesh):
    layer = mesh.uv_layers.active
    assert layer is not None
    uvs = array.array('f', bytes(8 * len(mesh.loops)))
    mesh.attributes[layer.name].data.foreach_get("vector", uvs)
    assert all(math.isfinite(x) for x in uvs)
    mesh.calc_loop_triangles()
    areas = []
    for tri in mesh.loop_triangles:
        (u0, v0), (u1, v1), (u2, v2) = (uvs[2 * c:2 * c + 2] for c in tri.loops)
        areas.append(0.5 * ((u1 - u0) * (v2 - v0) - (u2 - u0) * (v1 - v0)))
    mean = sum(abs(a) for a in areas) / len(areas)
    return sum(1 for a in areas if a < -1e-3 * mean), len(areas)


@test
def optimal_delaunay_is_one_bool_that_is_part_of_the_option_and_the_session_key():
    cutmod = importlib.import_module(ARGS.module + ".cut")
    assert cutmod.option_of('INTRINSIC') == "intrinsic"
    assert cutmod.option_of('INTRINSIC', optimal_delaunay=True) == "intrinsic+odt"
    assert cutmod.option_of('REFINED', 30.0, True) == "refined:30+odt"
    assert cutmod.option_keywords("refined:30+odt") == {
        "triangulation": "refined", "refine_bound": 30.0, "optimal_delaunay": True}
    assert "optimal_delaunay" not in cutmod.option_keywords("intrinsic")
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    # Without the flag, then with it on the same mesh and cut: another option, another session.
    assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=0.9, vertex=5,
                                    triangulation='INTRINSIC',
                                    optimal_delaunay=False) == {'FINISHED'}
    assert mesh["ddg_cut_option"] == "intrinsic"
    plain, _, _ = cutmod.open_session(obj, 1e-3, 0.1)
    assert plain.optimal_delaunay is False
    assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=0.9, vertex=5,
                                    triangulation='INTRINSIC',
                                    optimal_delaunay=True) == {'FINISHED'}
    assert mesh["ddg_cut_option"] == "intrinsic+odt"
    on, _, _ = cutmod.open_session(obj, 1e-3, 0.1)
    assert on is not plain and on.optimal_delaunay is True
    assert on.triangulation == "intrinsic"
    again, _, _ = cutmod.open_session(obj, 1e-3, 0.1)
    assert again is on, "the session with the flag was rebuilt although nothing changed"
    # The option is part of the key of the kept session.
    token = cutmod.gather(obj)[1]
    assert cutmod.cache_key(mesh, token, "intrinsic", 1e-3, 0.1) != cutmod.cache_key(
        mesh, token, "intrinsic+odt", 1e-3, 0.1)


@test
def the_representation_is_part_of_the_option_and_reaches_the_session():
    cutmod = importlib.import_module(ARGS.module + ".cut")
    assert cutmod.option_of('INTRINSIC', representation='INTEGER') == "intrinsic"
    assert cutmod.option_of('INTRINSIC', representation='SIGNPOSTS') == "intrinsic+signposts"
    assert cutmod.option_of('REFINED', 30.0, True, 'SIGNPOSTS') == "refined:30+odt+signposts"
    assert cutmod.option_of('CALLER', representation='SIGNPOSTS') == "caller"
    assert cutmod.option_keywords("refined:30+odt+signposts") == {
        "triangulation": "refined", "refine_bound": 30.0, "optimal_delaunay": True,
        "representation": "signposts"}
    assert "representation" not in cutmod.option_keywords("intrinsic")
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    states = {}
    for representation, option in (('INTEGER', "intrinsic"),
                                   ('SIGNPOSTS', "intrinsic+signposts")):
        assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=0.9, vertex=5,
                                        triangulation='INTRINSIC', optimal_delaunay=False,
                                        representation=representation) == {'FINISHED'}
        assert mesh["ddg_cut_option"] == option
        session, _, _ = cutmod.open_session(obj, 1e-3, 0.1)
        assert session.representation == representation.lower(), session.representation
        assert flow(steps=3) == {'FINISHED'}
        states[representation] = state(mesh)
    # The same flips: the same cut to rounding (the attributes are float32).
    worst = max(abs(a - b) for ra, rb in zip(states['INTEGER'], states['SIGNPOSTS'])
                for a, b in zip(ra, rb))
    assert worst < 1e-5, worst


@test
def optimal_delaunay_cut_completes_applies_and_unwraps_without_flipped_triangles():
    cutmod = importlib.import_module(ARGS.module + ".cut")
    for name, triangulation, option in (
            ("uv_sphere", 'INTRINSIC', "intrinsic+odt"),
            ("subdivided_cube", 'REFINED', "refined:25+odt")):
        obj = MODELS[name]()
        mesh = obj.data
        radius = 0.3 * diagonal(mesh)
        assert bpy.ops.uv.ddg_cut_start(
            mode='CIRCLE', radius=radius, vertex=0, triangulation=triangulation,
            optimal_delaunay=True) == {'FINISHED'}
        assert mesh["ddg_cut_option"] == option
        assert flow(steps=10, show_advanced=True) == {'FINISHED'}
        report = dict(mesh["ddg_cut_last_flow"])
        assert report["triangulation"] == option, report["triangulation"]
        assert report["energy_after"] < report["energy_before"]
        session, _, _ = cutmod.open_session(obj, 1e-3, 0.1)
        info = session.intrinsic_report
        assert session.optimal_delaunay and info["optimal_delaunay_status"], info
        print("  {} {}: {} inserted by iODT, {} iterations, {}, setup {:.2f} s".format(
            name, option, info["n_optimized"], info["optimal_delaunay_iterations"],
            info["optimal_delaunay_status"], info["seconds"]))
        assert bpy.ops.uv.ddg_cut_apply(unwrap_after=True) == {'FINISHED'}
        for key in ("ddg_cut_state", "ddg_cut_option", "ddg_cut_last_flow"):
            assert key not in mesh, key
        flipped, n = flipped_uv_triangles(mesh)
        assert flipped == 0, "{} {}: {} flipped triangles of {}".format(name, option, flipped, n)
        assert sum(e.use_seam for e in mesh.edges) > 0
        reset()


@test
def optimal_delaunay_with_the_meshs_own_triangulation_is_ignored_with_a_message():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    with CaptureOutput() as out:
        assert bpy.ops.uv.ddg_cut_start(mode='CIRCLE', radius=0.9, vertex=5,
                                        triangulation='CALLER',
                                        optimal_delaunay=True) == {'FINISHED'}
    assert mesh["ddg_cut_option"] == "caller"
    assert "Optimal Delaunay needs an intrinsic triangulation" in out.text, out.text
    assert flow(steps=3) == {'FINISHED'}
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}


@test
def sessions_are_kept_between_calls_and_rebuilt_when_the_mesh_changes():
    cutmod = importlib.import_module(ARGS.module + ".cut")
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    start_circle(obj, vertex=5)
    first, _, _ = cutmod.open_session(obj, 1e-3, 0.1)
    assert flow(steps=2) == {'FINISHED'}
    second, _, _ = cutmod.open_session(obj, 1e-3, 0.1)
    assert second is first, "the session was rebuilt although nothing changed"
    # The patch count of the library is the one of an independent pass over the corners.
    assert second.n_patches == cutmod.count_patches(mesh, second.corner_regions())
    # Attributes edited by hand: the kept session takes them (no state bytes then).
    if "ddg_cut_state" not in mesh:
        rows = state(mesh)
        values = array.array('f', (-v for v in rows[0]))
        mesh.attributes["ddg_cut_region_0"].data.foreach_set("value", values)
        mesh.attributes["ddg_cut_region_1"].data.foreach_set("value", array.array('f', rows[0]))
        third, _, _ = cutmod.open_session(obj, 1e-3, 0.1)
        assert third is first
        assert max(abs(a + b) for a, b in zip(third.state_row(0), rows[0])) < 1e-6
        load_state(mesh, (rows, None, mesh.get("ddg_cut_option")))
    # Another offset, a moved vertex, a changed seam: a new session.
    other, _, _ = cutmod.open_session(obj, 0.05, 0.1)
    assert other is not first
    original_x = mesh.vertices[3].co.x
    mesh.vertices[3].co.x += 0.01
    if "ddg_cut_state" in mesh:
        # The state bytes belong to the triangulation of the unmoved mesh: stale.
        expect_error("changed since 'Cut: Start'", bpy.ops.uv.ddg_cut_flow, **FLOW)
        mesh.vertices[3].co.x = original_x
    else:
        moved, _, _ = cutmod.open_session(obj, 1e-3, 0.1)
        assert moved is not first
    assert flow(steps=1) == {'FINISHED'}
    # Apply drops the kept sessions of the mesh.
    assert bpy.ops.uv.ddg_cut_apply() == {'FINISHED'}
    assert not [k for k in cutmod._CACHE if k[0] == mesh.name_full]


@test
def flow_report_has_a_triangulation_line_in_the_advanced_view_only():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    start_circle(obj, vertex=5)
    with CaptureOutput() as out:
        assert flow(steps=1) == {'FINISHED'}
    assert "triangulation:" not in out.text
    with CaptureOutput() as out:
        assert flow(steps=1, show_advanced=True) == {'FINISHED'}
    assert "triangulation:" in out.text, out.text


# ----------------------------------------------------------------------------------------------
# Flipped UV triangles: before (bends kept) and after (bends merged)


def unwrap_flips(obj):
    """(the library's count, our count on the triangulation the unwrap used, charts)."""
    core = importlib.import_module(ARGS.module + ".core")
    mesh = obj.data
    stats = core.unwrap_object(obj, use_seams=True, mark_cuts_as_seams=False)
    layer = mesh.uv_layers.active
    uvs = array.array('f', bytes(8 * len(mesh.loops)))
    mesh.attributes[layer.name].data.foreach_get("vector", uvs)
    _, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
    triangles, _ = core.stored_triangulation(mesh, corner_vertices, polygon_starts, triangles)
    areas = []
    for k in range(0, len(triangles), 3):
        (u0, v0), (u1, v1), (u2, v2) = (uvs[2 * c:2 * c + 2] for c in triangles[k:k + 3])
        areas.append(0.5 * ((u1 - u0) * (v2 - v0) - (u2 - u0) * (v1 - v0)))
    mean = sum(abs(a) for a in areas) / len(areas)
    own = sum(1 for a in areas if a < -1e-3 * mean)
    return stats["report"].get("n_flipped_uv_triangles", 0), own, stats["report"].get("n_charts")


def flip_case(label, make, start, flow_props, offset, merge):
    reset()
    obj = make()
    mesh = obj.data
    props = dict(start)
    if props.get("mode") == 'CIRCLE' and props.get("vertex") == "centre":
        centre = Vector((0.511, 0.493, 0.0))
        props["vertex"] = min(range(len(mesh.vertices)),
                              key=lambda i: (mesh.vertices[i].co - centre).length)
    assert bpy.ops.uv.ddg_cut_start(cut_offset=offset, merge_tolerance=merge, **props) \
        == {'FINISHED'}
    assert bpy.ops.uv.ddg_cut_flow(cut_offset=offset, merge_tolerance=merge, **flow_props) \
        == {'FINISHED'}
    assert bpy.ops.uv.ddg_cut_apply(cut_offset=offset, merge_tolerance=merge) == {'FINISHED'}
    return unwrap_flips(obj)


def flip_cases():
    spot_flow = dict(steps=60, use_stop_rule=True, length_mode='EXPLICIT', stages='RECOMMENDED')
    clusters8 = dict(mode='CLUSTERS', regions=8, cluster_iterations=50)
    flat_flow = dict(steps=20, use_stop_rule=False, length_mode='RIPPLE', ripple_wavelength=0.25,
                     use_uniform_curvature=True, uniform_curvature=5.0)
    flat_start = dict(mode='CIRCLE', radius=0.25, vertex="centre")
    sphere = lambda: primitive("uv_sphere", segments=24, ring_count=12)    # noqa: E731
    cases = [
        ("spot atlas8 w2", spot_object, clusters8, dict(spot_flow, length_weight=2.0)),
        ("spot atlas8 w10", spot_object, clusters8, dict(spot_flow, length_weight=10.0)),
        ("sphere atlas4", sphere, dict(mode='CLUSTERS', regions=4),
         dict(steps=6, use_stop_rule=False, length_mode='EXPLICIT', length_weight=0.3)),
        ("sphere circle", sphere, dict(mode='CIRCLE', radius=0.9, vertex=5),
         dict(steps=6, use_stop_rule=False, length_mode='EXPLICIT', length_weight=0.1)),
        ("sheared 153", lambda: mesh_object("s153", *sheared_grid(32, 2.0)), flat_start, flat_flow),
        ("needles and caps", lambda: mesh_object("nc", *needles_and_caps_grid(32)), flat_start,
         flat_flow),
        ("sheared 165", lambda: mesh_object(
            "s165", *sheared_grid(32, 1.0 / math.tan(math.radians(15.0)))), flat_start, flat_flow),
    ]
    return cases


@test
def no_flipped_uv_triangle_after_apply_and_unwrap_before_and_after_merging_bends():
    """The same cuts with the bends of the cut kept (merge 0) and merged (the default 0.1): the
    library's count and ours, per case and offset. Asserted: zero with the bends merged, on
    every triangulation (Optimal Delaunay included) and offset."""
    rows = []
    for label, make, start, flow_props in flip_cases():
        for offset in ((CUT_OFFSET, 0.2) if label.startswith("spot") else (CUT_OFFSET,)):
            if label.startswith("spot") and not os.path.exists(ARGS.spot):
                raise Skip("{} is not present".format(ARGS.spot))
            before = flip_case(label, make, start, flow_props, offset, 0.0)
            after = flip_case(label, make, start, flow_props, offset, 0.1)
            rows.append((label, offset, before, after))
            print("  FLIPTABLE [{}] {:<18s} offset {:<6g} before: library {} ours {} | after: "
                  "library {} ours {} (charts {} -> {})".format(
                      ARGS.triangulation, label, offset, before[0], before[1], after[0], after[1],
                      before[2], after[2]))
    bad = [r for r in rows if r[3][0] != 0 or r[3][1] != 0]
    assert not bad, "flipped triangles with the bends merged: {}".format(bad)


# ----------------------------------------------------------------------------------------------
# Error paths


@test
def errors_are_readable_messages():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    # Flow and Apply before Start.
    message = expect_error("no cut state.*Cut: Start", bpy.ops.uv.ddg_cut_flow)
    expect_error("no cut state.*Cut: Start", bpy.ops.uv.ddg_cut_apply)
    # Start: nothing selected, everything selected, no active vertex, bad radius and vertex,
    # more than two regions.
    mesh.polygons.foreach_set("select", [False] * len(mesh.polygons))
    expect_error("No faces are selected", bpy.ops.uv.ddg_cut_start, mode='SELECTED')
    mesh.polygons.foreach_set("select", [True] * len(mesh.polygons))
    expect_error("All faces are selected", bpy.ops.uv.ddg_cut_start, mode='SELECTED')
    expect_error("no active vertex", bpy.ops.uv.ddg_cut_start, mode='CIRCLE', radius=0.3)
    expect_error("not exist", bpy.ops.uv.ddg_cut_start, mode='CIRCLE', radius=0.3, vertex=10 ** 6)
    mesh.polygons.foreach_set("select", [i % 2 == 0 for i in range(len(mesh.polygons))])
    expect_error("only one value", bpy.ops.uv.ddg_cut_start, mode='MATERIALS')
    expect_error("no attribute named 'nothing'", bpy.ops.uv.ddg_cut_start, mode='ATTRIBUTE',
                 attribute="nothing")
    assert not state(mesh), "a failed start left attributes behind"
    assert not previews(obj)
    # Flow: vertex group that does not exist, bad wavelength.
    start_circle(obj, vertex=5)
    expect_error("Vertex group 'Missing' does not exist", bpy.ops.uv.ddg_cut_flow,
                 penalty_group="Missing", **FLOW)
    expect_error("ripple wavelength must be larger than 0", bpy.ops.uv.ddg_cut_flow,
                 length_mode='RIPPLE', ripple_wavelength=0.0, steps=1)
    # A state that does not fit the mesh: corners of one vertex that no seam separates differ.
    attributes_only(mesh)
    values = state(mesh)[0]
    values[3] += 0.5
    mesh.attributes["ddg_cut_region_0"].data.foreach_set("value", values)
    message = expect_error("differs at vertex", bpy.ops.uv.ddg_cut_flow, **FLOW)
    assert "region function" in message
    expect_error("differs at vertex", bpy.ops.uv.ddg_cut_apply)
    # An attribute of the wrong type.
    mesh.attributes.remove(mesh.attributes["ddg_cut_region_1"])
    mesh.attributes.new("ddg_cut_region_1", 'INT', 'POINT')
    expect_error("must be a float attribute on the face-corner domain", bpy.ops.uv.ddg_cut_flow,
                 **FLOW)
    # A cut that has vanished cannot be applied.
    bpy.ops.object.delete()
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    start_circle(obj, vertex=5)
    attributes_only(mesh)
    mesh.attributes["ddg_cut_region_0"].data.foreach_set(
        "value", [-1.0] * len(mesh.loops))
    mesh.attributes["ddg_cut_region_1"].data.foreach_set(
        "value", [1.0] * len(mesh.loops))
    expect_error("cut is empty", bpy.ops.uv.ddg_cut_apply)
    # Not a mesh / no active object: the operators are not available (poll).
    bpy.ops.object.empty_add()
    assert not bpy.ops.uv.ddg_cut_start.poll() and not bpy.ops.uv.ddg_cut_flow.poll() \
        and not bpy.ops.uv.ddg_cut_apply.poll()
    del message


@test
def errors_for_meshes_the_flow_refuses():
    # Non-manifold: the message of the library names the edge.
    mesh = bpy.data.meshes.new("fins")
    mesh.from_pydata([(0, 0, 0), (1, 0, 0), (0.5, 1, 0), (0.5, -1, 0), (0.5, 0, 1)], [],
                     [(0, 1, 2), (1, 0, 3), (0, 1, 4)])
    obj = bpy.data.objects.new("fins", mesh)
    bpy.context.collection.objects.link(obj)
    activate(obj)
    mesh.polygons.foreach_set("select", [True, False, False])
    expect_error("shared by 3 faces", bpy.ops.uv.ddg_cut_start, mode='SELECTED')
    # No faces.
    empty = bpy.data.objects.new("empty", bpy.data.meshes.new("empty"))
    bpy.context.collection.objects.link(empty)
    activate(empty)
    expect_error("has no faces", bpy.ops.uv.ddg_cut_start, mode='CIRCLE', radius=1.0, vertex=0)
    # Shape keys cannot be carried through the insertion.
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    start_circle(obj, vertex=5)
    assert flow(steps=2) == {'FINISHED'}
    obj.shape_key_add(name="Basis")
    obj.shape_key_add(name="Other")
    message = expect_error("shape keys", bpy.ops.uv.ddg_cut_apply)
    assert state(obj.data), "the failed Apply removed the state"
    assert len(obj.data.vertices) == 266, "the failed Apply changed the mesh"
    del message


@test
def apply_failing_leaves_the_mesh_alone():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    mesh = obj.data
    start_circle(obj, vertex=5)
    nv, nf = len(mesh.vertices), len(mesh.polygons)
    attributes_only(mesh)
    rows = state(mesh)
    rows[0][7] += 1.0
    mesh.attributes["ddg_cut_region_0"].data.foreach_set("value", rows[0])
    expect_error("differs at vertex", bpy.ops.uv.ddg_cut_apply)
    assert (len(mesh.vertices), len(mesh.polygons)) == (nv, nf)
    assert state(mesh) and previews(obj)


@test
def native_module_without_cut_session_is_an_error_message():
    cut = importlib.import_module(ARGS.module + ".cut")
    core = importlib.import_module(ARGS.module + ".core")
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    import ddg_tools
    saved = ddg_tools.CutSession
    try:
        del ddg_tools.CutSession
        expect_error("has no CutSession", bpy.ops.uv.ddg_cut_start, mode='CIRCLE', radius=0.3,
                     vertex=0)
    finally:
        ddg_tools.CutSession = saved
    del cut, core, obj


@test
def existing_unwrap_operator_is_unaffected():
    obj = primitive("uv_sphere", segments=24, ring_count=12)
    assert bpy.ops.uv.ddg_conformal_unwrap() == {'FINISHED'}
    assert obj.data.uv_layers.active is not None


# ----------------------------------------------------------------------------------------------
# Runner


def main():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", default="ddg_uv_unwrap",
                        help="module name of the installed add-on")
    parser.add_argument("--source", default="", help="register the add-on from this directory")
    parser.add_argument("--site", default="", help="directory with an unpacked ddg_tools wheel")
    parser.add_argument("--spot", default=SPOT)
    parser.add_argument("--triangulation", default="intrinsic",
                        choices=["intrinsic", "caller", "refined"],
                        help="the triangulation the tests start their cuts on (the default of "
                        "Start is changed for this run); the old expectations hold on all")
    parser.add_argument("--optimal-delaunay", action="store_true",
                        help="the tests start their cuts with Optimal Delaunay on (the default "
                        "of that property of Start is changed for this run; ignored with the "
                        "mesh's own triangulation)")
    parser.add_argument("-k", dest="only", default="", help="run only tests whose name contains this")
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    ARGS = parser.parse_args(argv)

    assert bpy.app.background, "run with -b"
    if ARGS.site:
        sys.path.insert(0, ARGS.site)
    if ARGS.source:
        sys.path.insert(0, ARGS.source)
        module = importlib.import_module(ARGS.module)
        module.register()
    else:
        assert bpy.app.factory_startup, "run with --factory-startup"
        assert bpy.ops.preferences.addon_enable(module=ARGS.module) == {'FINISHED'}
        assert ARGS.module in bpy.context.preferences.addons.keys(), "the add-on did not enable"
    assert hasattr(bpy.types, "UV_OT_ddg_cut_start"), "the cut operators are not registered"
    if ARGS.triangulation != "intrinsic":
        # The tests call Start without a triangulation: change its default for this run.
        operators = importlib.import_module(ARGS.module + ".cut_operators")
        cls = operators.UV_OT_ddg_cut_start
        old = cls.__annotations__["triangulation"]
        bpy.utils.unregister_class(cls)
        keywords = dict(old.keywords)
        keywords["default"] = ARGS.triangulation.upper()
        cls.__annotations__["triangulation"] = bpy.props.EnumProperty(**keywords)
        bpy.utils.register_class(cls)
    if ARGS.optimal_delaunay:
        operators = importlib.import_module(ARGS.module + ".cut_operators")
        cls = operators.UV_OT_ddg_cut_start
        old = cls.__annotations__["optimal_delaunay"]
        bpy.utils.unregister_class(cls)
        keywords = dict(old.keywords)
        keywords["default"] = True
        cls.__annotations__["optimal_delaunay"] = bpy.props.BoolProperty(**keywords)
        bpy.utils.register_class(cls)
    print("TRIANGULATION", ARGS.triangulation, "OPTIMAL_DELAUNAY", ARGS.optimal_delaunay)

    passed, skipped, failed = [], [], []
    for function in TESTS:
        name = function.__name__
        if ARGS.only and ARGS.only not in name:
            continue
        reset()
        print("TEST", name)
        try:
            function()
        except Skip as skip:
            skipped.append(name)
            print("SKIP {}: {}".format(name, skip))
        except Exception:
            failed.append(name)
            traceback.print_exc(file=sys.stdout)
            print("FAIL", name)
        else:
            passed.append(name)
            print("PASS", name)
        sys.stdout.flush()

    for row in FLIP_RECORD:
        print("FLIPS atlas of spot [{}] weight {} offset {}: {} by our count, {} by the library's".format(*row))
    for row in BAD_RECORD:
        print("BAD", row)
    print("RESULT: {} passed, {} skipped, {} failed".format(len(passed), len(skipped), len(failed)))
    for name in failed:
        print("  failed:", name)
    if failed or not passed:
        raise SystemExit(1)
    print("ALL TESTS PASSED")


ARGS = None
main()
