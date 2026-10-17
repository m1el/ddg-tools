# SPDX-License-Identifier: MIT
"""Conformal Map to Sphere: store a conformal map of a closed genus-0 mesh to the unit sphere as
a point attribute (and, on request, as a shape key).

The map is ``ddg_tools.map_to_sphere`` (Boundary First Flattening, Sawhney and Crane 2017,
Sec. 7.5, with Moebius centering after Baden, Crane and Kazhdan 2018). The mesh itself is not
changed; only the attribute (and the shape key) are written.
"""

from array import array
from time import perf_counter

import bpy
from bpy.props import BoolProperty, StringProperty

from . import core

# Name of the shape key the operator writes when asked to.
SHAPE_KEY_NAME = "Conformal Sphere"


def map_object_to_sphere(obj, *, attribute_name="sphere", shape_key=False):
    """Maps the mesh of ``obj`` to the unit sphere and writes the result. Returns the report of
    ``ddg_tools.map_to_sphere`` with ``seconds_blender`` added. Raises ``core.UnwrapFailed`` with
    a message for the user."""
    start = perf_counter()
    mesh = obj.data
    if mesh.library is not None:
        raise core.UnwrapFailed("Mesh '{}' is linked from a library and cannot be written".format(
            mesh.name))
    if len(mesh.polygons) == 0:
        raise core.UnwrapFailed("Mesh '{}' has no faces".format(mesh.name))
    if not attribute_name:
        raise core.UnwrapFailed("The attribute name is empty")
    existing = mesh.attributes.get(attribute_name)
    if existing is not None and (existing.domain != 'POINT' or existing.data_type != 'FLOAT_VECTOR'):
        raise core.UnwrapFailed(
            "Mesh '{}' already has an attribute '{}' that is not a point vector".format(
                mesh.name, attribute_name))
    if attribute_name.startswith(".") or attribute_name in {"position", "sharp_face"}:
        raise core.UnwrapFailed("'{}' is a reserved attribute name".format(attribute_name))
    ddg = core.native()
    positions, corner_vertices, _starts, corner_triangles = core.mesh_buffers(mesh)
    # The native side wants vertex indices per triangle; Blender's triangles are corners.
    triangles = array('I', (corner_vertices[c] for c in corner_triangles))
    native_start = perf_counter()
    try:
        result = ddg.map_to_sphere(positions, triangles)
    except ddg.DdgError as ex:
        raise core.UnwrapFailed(str(ex)) from ex
    native_seconds = perf_counter() - native_start
    sphere = array('f', result.positions)
    attribute = mesh.attributes.get(attribute_name)
    if attribute is None:
        attribute = mesh.attributes.new(attribute_name, 'FLOAT_VECTOR', 'POINT')
    attribute.data.foreach_set("vector", sphere)
    if shape_key:
        if obj.data.shape_keys is None:
            obj.shape_key_add(name="Basis", from_mix=False)
        key = obj.data.shape_keys.key_blocks.get(SHAPE_KEY_NAME)
        if key is None:
            key = obj.shape_key_add(name=SHAPE_KEY_NAME, from_mix=False)
        # The unit sphere scaled to the mesh's area and centred on its vertex centroid, so that
        # the shape key morphs between shapes of comparable size.
        n = len(positions) // 3
        centre = [sum(positions[i::3]) / n for i in range(3)]
        area = sum(polygon.area for polygon in mesh.polygons)
        radius = (area / (4.0 * 3.141592653589793)) ** 0.5
        key.data.foreach_set("co", array('f', (
            centre[k % 3] + radius * sphere[k] for k in range(len(sphere)))))
    mesh.update()
    report = dict(result.report)
    report["seconds_native"] = native_seconds
    report["seconds_blender"] = perf_counter() - start - native_seconds
    return report


def summary(report):
    """One status line."""
    parts = [
        "pole vertex {}".format(report["pole"]),
        "angle distortion (sigma max/min) mean {:.3f}, max {:.3f}".format(
            report["quasi_conformal_mean"], report["quasi_conformal_max"]),
        "area distortion {:.3f}".format(report["area_distortion"]),
    ]
    bad = report["flipped_triangles"] + report["degenerate_triangles"]
    if bad:
        parts.append("{} flipped or collapsed triangles".format(bad))
    parts.append("{:.3f} s".format(report["seconds_native"] + report["seconds_blender"]))
    return "Mapped to sphere: " + " | ".join(parts)


class OBJECT_OT_ddg_map_to_sphere(bpy.types.Operator):
    bl_idname = "object.ddg_map_to_sphere"
    bl_label = "Conformal Map to Sphere"
    bl_description = (
        "Map the whole active mesh (closed, connected, genus 0) conformally to the unit sphere "
        "and store the result as a point attribute (optionally also as a shape key). The mesh "
        "itself is not changed. Method: Sawhney, Crane, Boundary First Flattening (2017)"
    )
    bl_options = {'REGISTER', 'UNDO'}

    attribute_name: StringProperty(
        name="Attribute",
        description="Name of the point attribute (float vector) that receives the sphere "
                    "coordinates",
        default="sphere",
    )
    shape_key: BoolProperty(
        name="Shape Key",
        description="Also write a shape key 'Conformal Sphere': the sphere, scaled to the "
                    "mesh's surface area and centred on its vertex centroid",
        default=False,
    )

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        was_edit_mode = obj.mode == 'EDIT'
        if was_edit_mode:
            bpy.ops.object.mode_set(mode='OBJECT')
        try:
            report = map_object_to_sphere(
                obj, attribute_name=self.attribute_name, shape_key=self.shape_key)
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        finally:
            if was_edit_mode:
                bpy.ops.object.mode_set(mode='EDIT')
        if not report["is_injective"]:
            self.report({'WARNING'}, "The map is not one-to-one: {} flipped or collapsed "
                        "triangles".format(report["flipped_triangles"]
                                           + report["degenerate_triangles"]))
        self.report({'INFO'}, summary(report))
        return {'FINISHED'}
