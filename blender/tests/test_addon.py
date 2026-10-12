# SPDX-License-Identifier: MIT
"""Headless tests of the installed extension. Runs inside Blender:

    blender -b --factory-startup --python-exit-code 1 --python blender/tests/test_addon.py -- \
        --module bl_ext.user_default.ddg_uv_unwrap [--zip dist/blender/ddg_uv_unwrap-X.Y.Z.zip]

`scripts/blender-test.sh` builds the wheel and the extension zip, installs the zip into a scratch
Blender user directory and then runs this file, so everything here goes through the installed
extension and its bundled wheel: the path a user takes.

Every unwrap goes through `run_checked`, which asserts that

* the operator finishes and leaves the object in the mode it was in;
* vertex, edge, face and corner counts, vertex positions, edge vertices, the corner order of
  every face and the seam flags are bit-for-bit what they were;
* the active UV map has one finite UV per face corner, inside the unit square;
* a second run gives bit-identical UVs.

The pipeline behind the operator is the real one (`ddg-uv`); `ddg_tools.IS_PLACEHOLDER` is
asserted to be False, so no test is conditional on a stand-in.
"""

import argparse
import array
import ctypes
import hashlib
import importlib
import math
import os
import re
import sys
import tempfile
import traceback
import zipfile

import bmesh
import bpy

OPERATOR = "uv.ddg_conformal_unwrap"
TESTS = []


def test(function):
    TESTS.append(function)
    return function


class Skip(Exception):
    pass


# ----------------------------------------------------------------------------------------------
# Helpers


def op(**props):
    return bpy.ops.uv.ddg_conformal_unwrap(**props)


def reset():
    """An empty scene; the add-on stays enabled."""
    active = bpy.context.view_layer.objects.active
    if active is not None and active.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    for mesh in list(bpy.data.meshes):
        bpy.data.meshes.remove(mesh)


def activate(obj):
    view_layer = bpy.context.view_layer
    view_layer.update()     # objects linked or removed since the last update
    for other in view_layer.objects:
        other.select_set(other is obj)
    view_layer.objects.active = obj
    return obj


def new_object(name, vertices, faces, edges=()):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(vertices, edges, faces)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    return activate(obj)


def primitive(kind, **kwargs):
    getattr(bpy.ops.mesh, "primitive_{}_add".format(kind))(**kwargs)
    return bpy.context.view_layer.objects.active


def remove_uv_layers(mesh):
    while mesh.uv_layers:
        mesh.uv_layers.remove(mesh.uv_layers[0])


def _get(collection, attribute, typecode, width):
    values = array.array(typecode, bytes(4 * width * len(collection)))
    collection.foreach_get(attribute, values)
    return values


def seam_flags(mesh):
    flags = bytearray(len(mesh.edges))
    mesh.edges.foreach_get("use_seam", memoryview(flags).cast('?'))
    return bytes(flags)


def snapshot(mesh):
    """Everything the operator must leave alone, as comparable values."""
    return {
        "counts (vertices, edges, faces, corners)": (
            len(mesh.vertices), len(mesh.edges), len(mesh.polygons), len(mesh.loops)),
        "vertex positions": _get(mesh.vertices, "co", 'f', 3).tobytes(),
        "edge vertices": _get(mesh.edges, "vertices", 'i', 2).tobytes(),
        "corner vertices": _get(mesh.loops, "vertex_index", 'i', 1).tobytes(),
        "corner edges": _get(mesh.loops, "edge_index", 'i', 1).tobytes(),
        "face corner starts": _get(mesh.polygons, "loop_start", 'i', 1).tobytes(),
        "face corner counts": _get(mesh.polygons, "loop_total", 'i', 1).tobytes(),
        "seams": seam_flags(mesh),
    }


def assert_same_mesh(before, after, ignore=()):
    for key, value in before.items():
        if key not in ignore:
            assert after[key] == value, "the operator changed the mesh: {}".format(key)


def uv_values(obj):
    """The active UV map of `obj`, 2 floats per face corner in corner order."""
    mesh = obj.data
    if obj.mode == 'EDIT':
        # In Edit mode the UV arrays of the datablock read as empty; the edit mesh has them.
        edit_mesh = bmesh.from_edit_mesh(mesh)
        layer = edit_mesh.loops.layers.uv.active
        assert layer is not None, "the edit mesh has no active UV map"
        values = array.array('f')
        for face in edit_mesh.faces:
            for loop in face.loops:
                values.extend(loop[layer].uv)
        return values
    layer = mesh.uv_layers.active
    assert layer is not None, "the mesh has no active UV map"
    assert len(layer.uv) == len(mesh.loops), "UV map size {} != {} corners".format(
        len(layer.uv), len(mesh.loops))
    return _get(layer.uv, "vector", 'f', 2)


def assert_valid_uvs(obj):
    uvs = uv_values(obj)
    n_corners = len(current_mesh(obj).loops)
    assert len(uvs) == 2 * n_corners, "{} UV values for {} corners".format(len(uvs), n_corners)
    assert all(math.isfinite(x) for x in uvs), "non-finite UV"
    low, high = min(uvs), max(uvs)
    assert low >= -1e-6 and high <= 1.0 + 1e-6, "UVs outside the unit square: {}..{}".format(
        low, high)
    return uvs


def current_mesh(obj):
    """The mesh datablock, brought up to date if the object is in Edit mode."""
    if obj.mode == 'EDIT':
        obj.update_from_editmode()
    return obj.data


def run_checked(obj, **props):
    """Runs the operator twice on `obj` and checks the invariants listed in the module
    docstring. Returns the UVs."""
    mode = obj.mode
    before = snapshot(current_mesh(obj))
    assert op(**props) == {'FINISHED'}
    assert obj.mode == mode, "mode changed from {} to {}".format(mode, obj.mode)
    assert_same_mesh(before, snapshot(current_mesh(obj)))
    first = assert_valid_uvs(obj).tobytes()

    assert op(**props) == {'FINISHED'}
    assert_same_mesh(before, snapshot(current_mesh(obj)))
    second = assert_valid_uvs(obj)
    assert second.tobytes() == first, "two runs gave different UVs"
    return second


def expect_error(pattern, **props):
    """The operator must report an error that matches `pattern`, without a traceback."""
    try:
        result = op(**props)
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


def face_components(mesh, skip_edges=frozenset()):
    """Connected components of faces across shared edges (not across `skip_edges`), computed
    here independently of the add-on. Returns one component id per face, numbered by lowest
    face."""
    parent = list(range(len(mesh.polygons)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    edge_face = {}
    for polygon in mesh.polygons:
        for loop_index in polygon.loop_indices:
            edge = mesh.loops[loop_index].edge_index
            if edge in skip_edges:
                continue
            other = edge_face.setdefault(edge, polygon.index)
            if other != polygon.index:
                a, b = find(other), find(polygon.index)
                parent[max(a, b)] = min(a, b)
    ids = {}
    return [ids.setdefault(find(p), len(ids)) for p in range(len(mesh.polygons))]


def chart_boxes(mesh, uvs, components):
    boxes = {}
    for polygon in mesh.polygons:
        box = boxes.setdefault(components[polygon.index],
                               [math.inf, math.inf, -math.inf, -math.inf])
        for corner in polygon.loop_indices:
            u, v = uvs[2 * corner], uvs[2 * corner + 1]
            box[0], box[1] = min(box[0], u), min(box[1], v)
            box[2], box[3] = max(box[2], u), max(box[3], v)
    return boxes


def assert_boxes_disjoint(boxes, gap):
    keys = sorted(boxes)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            p, q = boxes[a], boxes[b]
            distance = max(p[0] - q[2], q[0] - p[2], p[1] - q[3], q[1] - p[3])
            assert distance >= gap - 1e-6, "charts {} and {} are {} apart, expected {}".format(
                a, b, distance, gap)


# ----------------------------------------------------------------------------------------------
# Tests


@test
def packaged_module_is_the_bundled_wheel():
    import ddg_tools
    extensions = os.path.realpath(bpy.utils.user_resource('EXTENSIONS'))
    module_file = os.path.realpath(ddg_tools.__file__)
    assert module_file.startswith(extensions + os.sep), (
        "ddg_tools was imported from {}, not from the extension directory {}".format(
            module_file, extensions))
    addon = importlib.import_module(ARGS.module)
    assert os.path.realpath(addon.__file__).startswith(extensions + os.sep)
    print("  ddg_tools {} ({}), Python {}, Blender {}".format(
        ddg_tools.__version__, ddg_tools.PIPELINE, sys.version.split()[0],
        bpy.app.version_string))
    try:
        import numpy
        print("  numpy", numpy.__version__)
    except ImportError:
        print("  numpy absent")
    if ARGS.zip:
        with zipfile.ZipFile(ARGS.zip) as outer:
            wheels = [n for n in outer.namelist() if n.endswith(".whl")]
            assert wheels, "no wheel in " + ARGS.zip
            expected = set()
            for wheel in wheels:
                with zipfile.ZipFile(outer.open(wheel)) as inner:
                    for name in inner.namelist():
                        if name.endswith((".so", ".pyd")):
                            expected.add(hashlib.sha256(inner.read(name)).hexdigest())
        binaries = [f for f in os.listdir(os.path.dirname(module_file))
                    if f.endswith((".so", ".pyd"))]
        assert len(binaries) == 1, binaries
        with open(os.path.join(os.path.dirname(module_file), binaries[0]), "rb") as fh:
            assert hashlib.sha256(fh.read()).hexdigest() in expected, (
                "the loaded native module is not the one in " + ARGS.zip)


@test
def operator_is_registered_with_exactly_the_agreed_properties():
    cls = bpy.types.UV_OT_ddg_conformal_unwrap
    assert cls.bl_idname == OPERATOR
    assert {'REGISTER', 'UNDO'} <= set(cls.bl_options)
    rna = bpy.ops.uv.ddg_conformal_unwrap.get_rna_type()
    properties = {p.identifier: p for p in rna.properties if p.identifier != "rna_type"}
    assert set(properties) == {
        "lambda_", "use_seams", "pinned_group", "forbidden_group", "positive_cones_only",
        "mark_seams"}, sorted(properties)
    assert abs(properties["lambda_"].default - 0.1) < 1e-6
    assert properties["use_seams"].default is True
    assert properties["positive_cones_only"].default is True
    assert properties["mark_seams"].default is False
    assert properties["pinned_group"].default == "" and properties["forbidden_group"].default == ""
    for menu in (bpy.types.VIEW3D_MT_uv_map, bpy.types.VIEW3D_MT_object):
        modules = [f.__module__ for f in menu.draw._draw_funcs]
        assert ARGS.module in modules, "{} has no entry of the add-on".format(menu.__name__)


@test
def poll_wants_an_active_mesh_object_in_object_or_edit_mode():
    assert not bpy.ops.uv.ddg_conformal_unwrap.poll(), "poll passed without an active object"
    bpy.ops.object.empty_add()
    assert not bpy.ops.uv.ddg_conformal_unwrap.poll(), "poll passed for an empty"
    cube = primitive("cube")
    assert bpy.ops.uv.ddg_conformal_unwrap.poll()
    bpy.ops.object.mode_set(mode='EDIT')
    assert bpy.ops.uv.ddg_conformal_unwrap.poll()
    bpy.ops.object.mode_set(mode='SCULPT')
    assert not bpy.ops.uv.ddg_conformal_unwrap.poll(), "poll passed in Sculpt mode"
    bpy.ops.object.mode_set(mode='OBJECT')
    assert cube.mode == 'OBJECT'


@test
def cube_quads_existing_and_new_uv_map():
    cube = primitive("cube")
    mesh = cube.data
    assert [len(p.vertices) for p in mesh.polygons] == [4] * 6
    # Blender's cube comes with a UV map: it is overwritten, not duplicated.
    assert [layer.name for layer in mesh.uv_layers] == ["UVMap"]
    old = uv_values(cube).tobytes()
    new = run_checked(cube)
    assert [layer.name for layer in mesh.uv_layers] == ["UVMap"]
    assert new.tobytes() != old

    # Without any UV map, one is created.
    remove_uv_layers(mesh)
    assert len(mesh.uv_layers) == 0
    created = run_checked(cube)
    assert [layer.name for layer in mesh.uv_layers] == ["UVMap"]
    assert created.tobytes() == new.tobytes()

    # With two UV maps, only the active one is written.
    second = mesh.uv_layers.new(name="Second", do_init=False)
    second.uv.foreach_set("vector", array.array('f', [0.25] * (2 * len(mesh.loops))))
    mesh.uv_layers.active = mesh.uv_layers["Second"]
    first_before = _get(mesh.uv_layers["UVMap"].uv, "vector", 'f', 2).tobytes()
    run_checked(cube)
    assert mesh.uv_layers.active.name == "Second"
    assert _get(mesh.uv_layers["UVMap"].uv, "vector", 'f', 2).tobytes() == first_before
    assert _get(mesh.uv_layers["Second"].uv, "vector", 'f', 2).tobytes() == new.tobytes()


@test
def ngons_cylinder_caps_and_concave_faces():
    cylinder = primitive("cylinder", vertices=12, end_fill_type='NGON')
    sizes = sorted(len(p.vertices) for p in cylinder.data.polygons)
    assert sizes == [4] * 12 + [12, 12], sizes
    remove_uv_layers(cylinder.data)
    run_checked(cylinder)

    # An L-shaped (concave) hexagon next to a quad, and a non-planar pentagon on top.
    vertices = [(0, 0, 0), (2, 0, 0), (2, 1, 0), (1, 1, 0), (1, 2, 0), (0, 2, 0),
                (0, -1, 0), (2, -1, 0), (2.5, 0.5, 0.7), (3, 0, 0.2)]
    faces = [(0, 1, 2, 3, 4, 5), (6, 7, 1, 0), (1, 7, 9, 8, 2)]
    shape = new_object("ngons", vertices, faces)
    mesh = shape.data
    mesh.calc_loop_triangles()
    assert len(mesh.loop_triangles) == len(mesh.loops) - 2 * len(mesh.polygons) == 9
    run_checked(shape)


@test
def face_with_a_zero_area_triangle():
    # Vertices 1 and 2 are at the same place, so whichever way Blender splits the quad, one
    # triangle has no area. That must never crash: the pipeline mollifies the edge lengths (and
    # says so) and unwraps the mesh.
    import ddg_tools
    assert ddg_tools.IS_PLACEHOLDER is False
    vertices = [(0, 0, 0), (1, 0, 0), (1, 0, 0), (0, 1, 0), (0, 2, 0), (1, 2, 0)]
    shape = new_object("degenerate", vertices, [(0, 1, 2, 3), (3, 2, 5, 4)])
    with CaptureOutput() as output:
        run_checked(shape)
    assert "no area" in output.text and "mollification" in output.text, output.text
    core = importlib.import_module(ARGS.module + ".core")
    report = core.unwrap_object(shape)["report"]
    assert report["n_zero_area_triangles"] == 1 and report["mollified"] is True, report


@test
def monkey_three_components_do_not_overlap():
    monkey = primitive("monkey")
    mesh = monkey.data
    assert (len(mesh.vertices), len(mesh.polygons)) == (507, 500)
    uvs = run_checked(monkey)
    components = face_components(mesh)
    assert len(set(components)) == 3, "Suzanne has a head and two eyes"
    core = importlib.import_module(ARGS.module + ".core")
    stats = core.unwrap_object(monkey)
    assert stats["report"]["n_charts"] == 3, stats["report"]["n_charts"]
    assert_boxes_disjoint(chart_boxes(mesh, uvs, components), 2 * core.PACK_MARGIN)


@test
def torus():
    torus_object = primitive("torus")
    assert len(torus_object.data.polygons) == 576
    run_checked(torus_object)


@test
def uv_sphere_with_triangle_fans_at_the_poles():
    sphere = primitive("uv_sphere")
    sizes = {len(p.vertices) for p in sphere.data.polygons}
    assert sizes == {3, 4}, sizes
    run_checked(sphere)


@test
def buffers_agree_with_the_element_api():
    """`core` reads positions, corner vertices and seams through attributes where that is
    faster; whatever path it takes, the values must be the ones the element API gives."""
    core = importlib.import_module(ARGS.module + ".core")
    paths = core._attribute_paths()
    print("  attribute shortcuts in this Blender:", paths)
    if bpy.app.version[:2] == (5, 0):
        assert paths == {"position": "position", "corner_vert": ".corner_vert",
                         "seam": "uv_seam"}, paths
    assert not any(m.name.startswith(".ddg_uv_unwrap") for m in bpy.data.meshes), (
        "the probe mesh was left behind")
    for kind in ("monkey", "cylinder", "uv_sphere"):
        obj = primitive(kind)
        mesh = obj.data
        for index in range(0, len(mesh.edges), 7):
            mesh.edges[index].use_seam = True
        positions, corner_vertices, polygon_starts, triangles = core.mesh_buffers(mesh)
        assert positions.typecode == 'f'
        assert positions.tobytes() == _get(mesh.vertices, "co", 'f', 3).tobytes()
        assert list(corner_vertices) == [loop.vertex_index for loop in mesh.loops]
        assert list(polygon_starts) == [polygon.loop_start for polygon in mesh.polygons]
        mesh.calc_loop_triangles()
        assert list(triangles) == [c for t in mesh.loop_triangles for c in t.loops]
        expected = [v for edge in mesh.edges if edge.use_seam for v in edge.vertices]
        assert list(core.seam_edges(mesh)) == expected and expected
    # No seams at all: nothing is reported.
    assert len(core.seam_edges(primitive("cube").data)) == 0


@test
def marked_seams_are_cut_and_can_be_ignored():
    cube = primitive("cube")
    mesh = cube.data
    top = [e.index for e in mesh.edges
           if all(mesh.vertices[v].co.z > 0 for v in e.vertices)]
    assert len(top) == 4
    whole = run_checked(cube).tobytes()
    for index in top:
        mesh.edges[index].use_seam = True
    cut = run_checked(cube)                      # use_seams defaults to True
    assert cut.tobytes() != whole, "the seams had no effect"
    assert run_checked(cube, use_seams=False).tobytes() == whole

    core = importlib.import_module(ARGS.module + ".core")
    report = core.unwrap_object(cube)["report"]
    assert report["n_seam_edges"] == 4 and report["n_charts"] == 2, report
    assert core.unwrap_object(cube, use_seams=False)["report"]["n_seam_edges"] == 0
    # The two charts (lid, and the rest) are laid out apart from each other.
    components = face_components(mesh, skip_edges=frozenset(top))
    assert len(set(components)) == 2
    assert_boxes_disjoint(chart_boxes(mesh, cut, components), 2 * core.PACK_MARGIN)


@test
def cuts_can_be_marked_as_seams():
    import ddg_tools
    core = importlib.import_module(ARGS.module + ".core")
    sphere = primitive("uv_sphere", segments=16, ring_count=8)
    mesh = sphere.data
    mesh.edges[5].use_seam = True
    before = snapshot(mesh)
    # What the native module says it cut, asked independently of the operator.
    expected = ddg_tools.unwrap(*core.mesh_buffers(mesh), seam_edges=core.seam_edges(mesh))
    cut_pairs = {tuple(sorted(pair)) for pair in zip(expected.cut_edges[0::2],
                                                      expected.cut_edges[1::2])}
    assert cut_pairs, "a closed sphere cannot be flattened without a cut"

    assert op(mark_seams=False) == {'FINISHED'}
    assert seam_flags(mesh) == before["seams"], "seams changed although mark_seams is off"

    assert op(mark_seams=True) == {'FINISHED'}
    after = snapshot(mesh)
    assert_same_mesh(before, after, ignore=("seams",))
    assert_valid_uvs(sphere)
    marked = {tuple(sorted(e.vertices)) for e in mesh.edges if e.use_seam}
    assert marked == cut_pairs | {tuple(sorted(mesh.edges[5].vertices))}, (
        "seams after the run are not exactly the old seams plus the cut edges")
    assert len(marked) > 1

    # A mesh that never had a seam (Blender stores none of the seam data yet).
    torus = primitive("torus")
    mesh = torus.data
    assert not any(e.use_seam for e in mesh.edges)
    expected = ddg_tools.unwrap(*core.mesh_buffers(mesh))
    cut_pairs = {tuple(sorted(pair)) for pair in zip(expected.cut_edges[0::2],
                                                      expected.cut_edges[1::2])}
    assert cut_pairs
    before = snapshot(mesh)
    assert op(mark_seams=True) == {'FINISHED'}
    assert_same_mesh(before, snapshot(mesh), ignore=("seams",))
    assert {tuple(sorted(e.vertices)) for e in mesh.edges if e.use_seam} == cut_pairs
    # They are real seams for Blender too: Edit mode sees them.
    bpy.ops.object.mode_set(mode='EDIT')
    edit_mesh = bmesh.from_edit_mesh(mesh)
    assert {tuple(sorted(v.index for v in e.verts)) for e in edit_mesh.edges if e.seam} == cut_pairs
    bpy.ops.object.mode_set(mode='OBJECT')


@test
def edit_mode_uses_the_edited_mesh_and_stays_in_edit_mode():
    cube = primitive("cube")
    remove_uv_layers(cube.data)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.subdivide(number_cuts=2)       # the datablock is now out of date
    bpy.ops.mesh.select_all(action='DESELECT')
    edit_mesh = bmesh.from_edit_mesh(cube.data)
    edit_mesh.verts.ensure_lookup_table()
    edit_mesh.verts[3].select_set(True)
    edit_mesh.select_flush_mode()
    bmesh.update_edit_mesh(cube.data)
    counts = (len(edit_mesh.verts), len(edit_mesh.edges), len(edit_mesh.faces))
    assert counts == (56, 108, 54), counts

    uvs = run_checked(cube)
    assert cube.mode == 'EDIT'
    assert len(uvs) == 2 * 4 * 54

    edit_mesh = bmesh.from_edit_mesh(cube.data)
    assert (len(edit_mesh.verts), len(edit_mesh.edges), len(edit_mesh.faces)) == counts
    assert [v.index for v in edit_mesh.verts if v.select] == [3], "the selection changed"
    assert edit_mesh.loops.layers.uv.active is not None, "no UV layer in the edit mesh"

    # And the result is the same as unwrapping the same mesh from Object mode.
    bpy.ops.object.mode_set(mode='OBJECT')
    assert run_checked(cube).tobytes() == uvs.tobytes()


@test
def non_manifold_mesh_is_an_error_message_not_a_crash():
    fins = new_object(
        "fins",
        [(0, 0, 0), (1, 0, 0), (0.5, 1, 0), (0.5, -1, 0), (0.5, 0, 1)],
        [(0, 1, 2), (1, 0, 3), (0, 1, 4)])
    before = snapshot(fins.data)
    message = expect_error(r"edge \(0, 1\) is shared by 3 faces: 0, 1, 2")
    print("  reported:", message.strip())
    assert_same_mesh(before, snapshot(fins.data))
    assert len(fins.data.uv_layers) == 0, "a UV map was created although the unwrap failed"

    # The same from Edit mode: error, still in Edit mode, mesh untouched.
    bpy.ops.object.mode_set(mode='EDIT')
    expect_error(r"edge \(0, 1\) is shared by 3 faces: 0, 1, 2")
    assert fins.mode == 'EDIT'
    bpy.ops.object.mode_set(mode='OBJECT')
    assert_same_mesh(before, snapshot(fins.data))

    # A face glued into the middle of a cube: the typical accident in Blender.
    cube = primitive("cube")
    mesh = cube.data
    edit_mesh = bmesh.new()
    edit_mesh.from_mesh(mesh)
    edit_mesh.verts.ensure_lookup_table()
    diagonal = [v for v in edit_mesh.verts if abs(v.co.x - v.co.y) < 1e-6]
    assert len(diagonal) == 4
    bottom = sorted((v for v in diagonal if v.co.z < 0), key=lambda v: v.co.x)
    top = sorted((v for v in diagonal if v.co.z > 0), key=lambda v: v.co.x)
    edit_mesh.faces.new([bottom[0], bottom[1], top[1], top[0]])
    edit_mesh.to_mesh(mesh)
    edit_mesh.free()
    assert len(mesh.polygons) == 7
    activate(cube)
    # Either a cube edge with three faces, or a diagonal of the top or bottom face that
    # coincides with an edge of the new face; which one is found first depends on the indices.
    message = expect_error(r"edge \(\d+, \d+\) is shared by 3 (faces: \d+, \d+, \d+|"
                           r"triangles of the triangulation, in faces \d+, \d+)")
    print("  reported:", message.strip())

    # Two triangles that touch in one vertex only.
    bow = new_object(
        "bow",
        [(-1, -1, 0), (0, 0, 0), (-1, 1, 0), (1, -1, 0), (1, 1, 0)],
        [(0, 1, 2), (1, 3, 4)])
    message = expect_error(r"vertex 1 is non-manifold: its incident faces form 2 separate fans")
    print("  reported:", message.strip())
    assert len(bow.data.uv_layers) == 0


@test
def flipped_face_is_repaired():
    # Policy of the pipeline: inconsistent winding is repaired per connected piece (with a
    # warning), not refused. The mesh itself must stay as it is, flipped face included.
    cube = primitive("cube")
    edit_mesh = bmesh.new()
    edit_mesh.from_mesh(cube.data)
    edit_mesh.faces.ensure_lookup_table()
    bmesh.ops.reverse_faces(edit_mesh, faces=[edit_mesh.faces[2]])
    edit_mesh.to_mesh(cube.data)
    edit_mesh.free()
    before = snapshot(current_mesh(cube))
    assert op() == {'FINISHED'}
    assert_same_mesh(before, snapshot(current_mesh(cube)))
    assert_valid_uvs(cube)


@test
def mesh_without_faces_is_an_error_message():
    wire = new_object("wire", [(0, 0, 0), (1, 0, 0), (1, 1, 0)], [], edges=[(0, 1), (1, 2)])
    message = expect_error(r"has no faces")
    print("  reported:", message.strip())
    assert len(wire.data.uv_layers) == 0


@test
def linked_mesh_is_an_error_message():
    # A mesh linked from another .blend file is read-only.
    path = os.path.join(bpy.app.tempdir, "ddg_linked_source.blend")
    source = bpy.data.meshes.new("LinkedSource")
    source.from_pydata([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], [], [(0, 1, 2, 3)])
    bpy.data.libraries.write(path, {source})
    bpy.data.meshes.remove(source)
    with bpy.data.libraries.load(path, link=True) as (_data_from, data_to):
        data_to.meshes = ["LinkedSource"]
    linked = data_to.meshes[0]
    assert linked.library is not None
    obj = bpy.data.objects.new("uses linked mesh", linked)
    bpy.context.collection.objects.link(obj)
    activate(obj)
    message = expect_error(r"Mesh 'LinkedSource' is linked from another file")
    print("  reported:", message.strip())
    assert len(linked.uv_layers) == 0


@test
def missing_native_module_is_an_error_message():
    # What a user sees when the bundled wheel does not load on their system (for example a
    # wheel built for a newer C library): an error report, and the add-on stays usable.
    core = importlib.import_module(ARGS.module + ".core")
    primitive("cube")
    saved = sys.modules["ddg_tools"]
    sys.modules["ddg_tools"] = None          # makes `import ddg_tools` raise ImportError
    try:
        assert core.is_placeholder() is False
        message = expect_error(r"The native module 'ddg_tools' could not be loaded")
        print("  reported:", message.strip())
    finally:
        sys.modules["ddg_tools"] = saved
    assert op() == {'FINISHED'}


@test
def vertex_groups_reach_the_native_module():
    import ddg_tools
    core = importlib.import_module(ARGS.module + ".core")
    grid = primitive("grid", x_subdivisions=4, y_subdivisions=4)
    assert len(grid.data.vertices) == 25
    pin = grid.vertex_groups.new(name="pin")
    pin.add([1, 5, 12], 1.0, 'REPLACE')
    pin.add([6], 0.0, 'REPLACE')                 # weight 0 does not count
    nocone = grid.vertex_groups.new(name="no cone")
    nocone.add([0, 2, 24], 0.5, 'REPLACE')
    grid.vertex_groups.new(name="empty")
    assert list(core.group_vertices(grid, "pin")) == [1, 5, 12]
    assert list(core.group_vertices(grid, "no cone")) == [0, 2, 24]
    assert list(core.group_vertices(grid, "empty")) == []

    run_checked(grid, pinned_group="pin", forbidden_group="no cone")
    report = core.unwrap_object(grid, pinned_group="pin", forbidden_group="no cone")["report"]
    assert (report["n_pinned"], report["n_forbidden"]) == (3, 3), report
    assert core.unwrap_object(grid)["report"]["n_pinned"] == 0
    # Steering acts on interior vertices only: of the pinned 1, 5 and 12 only the centre of the
    # 5 x 5 grid is one, and all three forbidden vertices (0, 2, 24) are on the boundary.
    assert (report["n_pinned_applied"], report["n_forbidden_applied"]) == (1, 0), report

    # The same from Edit mode.
    bpy.ops.object.mode_set(mode='EDIT')
    run_checked(grid, pinned_group="pin", forbidden_group="no cone")
    bpy.ops.object.mode_set(mode='OBJECT')

    # A vertex in both groups is a contradiction the native side reports by index.
    nocone.add([12], 1.0, 'REPLACE')
    before = snapshot(grid.data)
    message = expect_error(r"vertex 12 is both pinned and forbidden",
                           pinned_group="pin", forbidden_group="no cone")
    print("  reported:", message.strip())
    assert_same_mesh(before, snapshot(grid.data))
    # A group that does not exist.
    message = expect_error(r"Vertex group 'nope' does not exist on object 'Grid'",
                           pinned_group="nope")
    print("  reported:", message.strip())


@test
def loose_vertices_and_wire_seams_do_not_stop_the_unwrap():
    """The pipeline ignores vertices that no face uses (and counts them), and a seam that is not
    an edge of a face (here a wire edge) is ignored with a warning."""
    core = importlib.import_module(ARGS.module + ".core")
    vertices = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (5, 5, 5), (6, 5, 5), (7, 7, 7)]
    shape = new_object("loose", vertices, [(0, 1, 2, 3)], edges=[(4, 5)])
    assert len(shape.data.edges) == 5
    wire = next(e for e in shape.data.edges if set(e.vertices) == {4, 5})
    wire.use_seam = True
    with CaptureOutput() as output:
        run_checked(shape)
    assert "1 seam edges are not edges of a face and were ignored" in output.text, output.text
    result = core.unwrap_object(shape)
    report = result["report"]
    assert report["n_unreferenced_vertices"] == 3, report
    assert report["n_seam_edges"] == 1 and report["n_seam_edges_ignored"] == 1, report
    assert report["n_seam_edges_cut"] == 0 and result["n_cut_edges"] == 0, report
    assert report["n_charts"] == 1 and result["n_cones"] == 0, report

    # With seams switched off the wire edge is not passed on at all.
    report = core.unwrap_object(shape, use_seams=False)["report"]
    assert report["n_seam_edges"] == 0 and report["n_seam_edges_ignored"] == 0, report
    assert report["n_unreferenced_vertices"] == 3, report


@test
def status_line_reports_charts_cones_distortion_and_timing():
    import ddg_tools
    primitive("monkey")
    with CaptureOutput() as output:
        assert op() == {'FINISHED'}
    lines = [line for line in output.text.splitlines() if "chart" in line]
    assert len(lines) == 1, output.text
    line = lines[0]
    pattern = (r"3 charts, \d+ cones?, \d+ cut edges? \| angle error max [\d.]+°, "
               r"mean [\d.]+° \| scale [\d.]+x to [\d.]+x of the mean .*"
               r"\| [\d.]+ s \(native [\d.]+ s, Blender side [\d.]+ s\)")
    assert re.search(pattern, line), line
    assert ddg_tools.IS_PLACEHOLDER is False
    assert line.startswith("Info: "), line


@test
def flat_grid_maps_every_corner_to_its_own_vertex():
    """Corner order and distortion, end to end: a flat sheet needs no cones and must come out
    undistorted, i.e. every corner's UV is the same similarity transform (rotation, uniform
    scale, translation; no mirroring) of the x and y of its vertex. A UV written to the wrong
    corner, or any distortion, shows as a residual."""
    core = importlib.import_module(ARGS.module + ".core")
    xs = [0.0, 0.3, 1.1, 2.0, 2.4]
    ys = [0.0, 0.5, 0.6, 1.5]
    vertices = [(x, y, 0.0) for y in ys for x in xs]
    faces = []
    for j in range(len(ys) - 1):
        for i in range(len(xs) - 1):
            v = j * len(xs) + i
            faces.append((v, v + 1, v + len(xs) + 1, v + len(xs)))
    # Mix in a triangle pair and rotate some faces so corner order is not regular.
    faces[0] = faces[0][1:] + faces[0][:1]
    faces[5] = faces[5][2:] + faces[5][:2]
    a, b, c, d = faces.pop()
    faces += [(a, b, c), (c, d, a)]
    sheet = new_object("sheet", vertices, faces)
    uvs = run_checked(sheet)
    check_flat_sheet(sheet.data, uvs)
    stats = core.unwrap_object(sheet)
    assert stats["n_cones"] == 0 and stats["n_cut_edges"] == 0 and stats["report"]["n_charts"] == 1


def check_flat_sheet(mesh, uvs, tolerance=1e-5):
    """Asserts that the UV of every corner is `a * (x + iy) + b` for one complex `a` (rotation
    and uniform scale) and one `b`: the least squares fit over all corners must leave a residual
    below `tolerance` (the UVs are float32 values in the unit square)."""
    z = [complex(mesh.vertices[loop.vertex_index].co.x, mesh.vertices[loop.vertex_index].co.y)
         for loop in mesh.loops]
    w = [complex(uvs[2 * loop.index], uvs[2 * loop.index + 1]) for loop in mesh.loops]
    z_mean, w_mean = sum(z) / len(z), sum(w) / len(w)
    spread = sum(abs(q - z_mean) ** 2 for q in z)
    a = sum((q - z_mean).conjugate() * (r - w_mean) for q, r in zip(z, w)) / spread
    assert abs(a) > 0.0, "the UVs collapsed to a point"
    b = w_mean - a * z_mean
    worst = max(range(len(z)), key=lambda i: abs(w[i] - (a * z[i] + b)))
    residual = abs(w[worst] - (a * z[worst] + b))
    assert residual < tolerance, (
        "corner {} (vertex {}): UV {} is {:.2e} away from the best similarity transform of "
        "the sheet (scale {:.4f})".format(
            worst, mesh.loops[worst].vertex_index, (w[worst].real, w[worst].imag), residual,
            abs(a)))


@test
def the_checks_themselves_catch_violations():
    """Negative controls: each check used above must fail when what it guards is broken."""
    import ddg_tools

    def must_fail(what, function):
        try:
            function()
        except AssertionError:
            return
        raise AssertionError("negative control passed: " + what)

    cube = primitive("cube")
    mesh = cube.data
    assert op() == {'FINISHED'}
    before = snapshot(mesh)

    # A moved vertex, a reordered face, an extra vertex, a new seam.
    mesh.vertices[0].co.x += 0.5
    must_fail("moved vertex", lambda: assert_same_mesh(before, snapshot(mesh)))
    mesh.vertices[0].co.x -= 0.5
    assert_same_mesh(before, snapshot(mesh))

    edit_mesh = bmesh.new()
    edit_mesh.from_mesh(mesh)
    edit_mesh.faces.ensure_lookup_table()
    face = edit_mesh.faces[0]
    corners = [v for v in face.verts]
    edit_mesh.faces.remove(face)
    edit_mesh.faces.new(corners[1:] + corners[:1])
    edit_mesh.to_mesh(mesh)
    edit_mesh.free()
    changed = snapshot(mesh)
    assert changed["counts (vertices, edges, faces, corners)"] == (8, 12, 6, 24)
    must_fail("rotated face corners", lambda: assert_same_mesh(before, changed))

    mesh.edges[0].use_seam = True
    must_fail("new seam", lambda: assert_same_mesh(changed, snapshot(mesh)))

    # UV checks: NaN, out of range, wrong count.
    layer = mesh.uv_layers.active
    good = _get(layer.uv, "vector", 'f', 2)
    bad = array.array('f', good)
    bad[7] = math.nan
    layer.uv.foreach_set("vector", bad)
    must_fail("NaN in the UVs", lambda: assert_valid_uvs(cube))
    bad[7] = 1.5
    layer.uv.foreach_set("vector", bad)
    must_fail("UV outside the unit square", lambda: assert_valid_uvs(cube))
    layer.uv.foreach_set("vector", good)
    assert_valid_uvs(cube)

    # Corner order and distortion: UVs shifted by one corner, or squeezed, on a flat sheet.
    sheet = new_object(
        "control", [(0, 0, 0), (2, 0, 0), (2, 1, 0), (0, 1, 0), (3, 0, 0), (3, 1, 0)],
        [(0, 1, 2, 3), (1, 4, 5, 2)])
    assert op() == {'FINISHED'}
    uvs = uv_values(sheet)
    check_flat_sheet(sheet.data, uvs)
    shifted = array.array('f', uvs[2:] + uvs[:2])
    must_fail("UVs shifted by one corner", lambda: check_flat_sheet(sheet.data, shifted))
    squeezed = array.array('f', uvs)
    for i in range(1, len(squeezed), 2):
        squeezed[i] *= 0.5
    must_fail("UVs squeezed in v", lambda: check_flat_sheet(sheet.data, squeezed))
    mirrored = array.array('f', uvs)
    for i in range(1, len(mirrored), 2):
        mirrored[i] = 1.0 - mirrored[i]
    must_fail("UVs mirrored", lambda: check_flat_sheet(sheet.data, mirrored))

    # Cut edges that are not edges of the mesh are counted, not marked.
    core = importlib.import_module(ARGS.module + ".core")
    quad = new_object("quad", [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], [(0, 1, 2, 3)])
    assert core.mark_seams(quad.data, array.array('i', [1, 0, 0, 2, 2, 3])) == (2, 1)
    assert [tuple(sorted(e.vertices)) for e in quad.data.edges if e.use_seam] == [(0, 1), (2, 3)]
    assert core.mark_seams(quad.data, array.array('i', [0, 1])) == (0, 0)

    # Chart separation.
    must_fail("overlapping boxes",
              lambda: assert_boxes_disjoint({0: [0, 0, 0.5, 0.5], 1: [0.4, 0.4, 1, 1]}, 0.01))


# ----------------------------------------------------------------------------------------------
# Runner


def main():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", required=True, help="module name of the installed add-on")
    parser.add_argument("--zip", default="", help="the extension zip that was installed")
    parser.add_argument("-k", dest="only", default="", help="run only tests whose name contains this")
    parser.add_argument("--bpy-module", action="store_true",
                        help="running in Blender-as-a-Python-module, which has no --factory-startup")
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    ARGS = parser.parse_args(argv)

    assert bpy.app.background, "run with -b"
    assert bpy.app.factory_startup or ARGS.bpy_module, "run with --factory-startup"
    # --factory-startup ignores the saved preferences, so the add-on is enabled here.
    assert ARGS.module not in bpy.context.preferences.addons.keys()
    assert bpy.ops.preferences.addon_enable(module=ARGS.module) == {'FINISHED'}
    assert ARGS.module in bpy.context.preferences.addons.keys(), "the add-on did not enable"
    assert hasattr(bpy.types, "UV_OT_ddg_conformal_unwrap"), "the operator is not registered"

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

    print("RESULT: {} passed, {} skipped, {} failed".format(len(passed), len(skipped), len(failed)))
    for name in failed:
        print("  failed:", name)
    if failed or not passed:
        raise SystemExit(1)
    print("ALL TESTS PASSED")


ARGS = None
main()
