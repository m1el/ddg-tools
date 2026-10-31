# SPDX-License-Identifier: MIT
"""DDG Conformal UV Unwrap: an operator that generates a UV map for the active mesh, three
that cut a mesh open where it flattens best, and operators that make a mesh developable.

The work is done by the native module ``ddg_tools`` (Rust, bundled as a wheel). See
``blender/README.md`` in the ddg-tools repository.

Credits. UV unwrap: Soliman, Slepcev, Crane, Optimal Cone Singularities for Conformal Flattening
(2018); Sawhney, Crane, Boundary First Flattening (2017). Cut: Sharp, Crane, Variational Surface
Cutting (2018), implemented clean-room from the text of the paper. Map to sphere: Sawhney,
Crane, Boundary First Flattening (2017), Sec. 7.5; Baden, Crane, Kazhdan, Moebius Registration
(2018). Developability: Stein, Grinspun, Crane, Developability of Triangle Meshes (2018),
implemented clean-room from the text of the paper (Make Developable, Show Developability,
Developable Seams + UV). Repulsive curves: Yu, Schumacher, Crane, Repulsive Curves (2021): a live
flow and an N-steps flow of edge meshes and Curve objects. Show Intrinsic Refinement: Gillespie, Sharp, Crane, Integer
Coordinates for Intrinsic Geometry Processing (2021), Fig. 12; Sharp, Soliman, Crane (2019).
"""

import bpy

from . import core
from .cut_operators import UV_OT_ddg_cut_apply, UV_OT_ddg_cut_flow, UV_OT_ddg_cut_start
from .developable_operator import (
    OBJECT_OT_ddg_make_developable,
    OBJECT_OT_ddg_show_developability,
    UV_OT_ddg_developable_seams,
)
from .intrinsic_operator import MESH_OT_ddg_show_intrinsic_refinement
from .operator import UV_OT_ddg_conformal_unwrap
from . import repulsive_operator
from .sphere_operator import OBJECT_OT_ddg_map_to_sphere

_CUT_OPERATORS = (UV_OT_ddg_cut_start, UV_OT_ddg_cut_flow, UV_OT_ddg_cut_apply)


def _draw_menu_entry(self, _context):
    layout = self.layout
    layout.separator()
    if core.is_placeholder():
        layout.operator(
            UV_OT_ddg_conformal_unwrap.bl_idname,
            text=UV_OT_ddg_conformal_unwrap.bl_label + " [placeholder build]",
        )
    else:
        layout.operator(UV_OT_ddg_conformal_unwrap.bl_idname)
    for operator in _CUT_OPERATORS:
        layout.operator(operator.bl_idname)


def _draw_sphere_entry(self, _context):
    self.layout.operator(OBJECT_OT_ddg_map_to_sphere.bl_idname)


_DEVELOPABLE_OPERATORS = (
    OBJECT_OT_ddg_make_developable,
    OBJECT_OT_ddg_show_developability,
    UV_OT_ddg_developable_seams,
)


def _draw_developable_object_entries(self, _context):
    layout = self.layout
    layout.separator()
    for operator in _DEVELOPABLE_OPERATORS:
        layout.operator(operator.bl_idname)


def _draw_developable_edit_entries(self, _context):
    layout = self.layout
    layout.separator()
    layout.operator(OBJECT_OT_ddg_make_developable.bl_idname)
    layout.operator(OBJECT_OT_ddg_show_developability.bl_idname)


def _draw_developable_uv_entry(self, _context):
    self.layout.operator(UV_OT_ddg_developable_seams.bl_idname)


def _draw_intrinsic_entry(self, _context):
    self.layout.operator(MESH_OT_ddg_show_intrinsic_refinement.bl_idname)


# The UV menu of the 3D viewport in Edit mode, and the Object menu.
_MENUS = (bpy.types.VIEW3D_MT_uv_map, bpy.types.VIEW3D_MT_object)


def register():
    bpy.utils.register_class(UV_OT_ddg_conformal_unwrap)
    for operator in _CUT_OPERATORS:
        bpy.utils.register_class(operator)
    for menu in _MENUS:
        menu.append(_draw_menu_entry)
    bpy.utils.register_class(OBJECT_OT_ddg_map_to_sphere)
    bpy.types.VIEW3D_MT_object.append(_draw_sphere_entry)
    for operator in _DEVELOPABLE_OPERATORS:
        bpy.utils.register_class(operator)
    bpy.types.VIEW3D_MT_object.append(_draw_developable_object_entries)
    bpy.types.VIEW3D_MT_edit_mesh.append(_draw_developable_edit_entries)
    bpy.types.VIEW3D_MT_uv_map.append(_draw_developable_uv_entry)
    bpy.utils.register_class(MESH_OT_ddg_show_intrinsic_refinement)
    bpy.types.VIEW3D_MT_object.append(_draw_intrinsic_entry)
    bpy.types.VIEW3D_MT_edit_mesh.append(_draw_intrinsic_entry)
    repulsive_operator.register()


def unregister():
    repulsive_operator.unregister()
    bpy.types.VIEW3D_MT_edit_mesh.remove(_draw_intrinsic_entry)
    bpy.types.VIEW3D_MT_object.remove(_draw_intrinsic_entry)
    bpy.utils.unregister_class(MESH_OT_ddg_show_intrinsic_refinement)
    bpy.types.VIEW3D_MT_uv_map.remove(_draw_developable_uv_entry)
    bpy.types.VIEW3D_MT_edit_mesh.remove(_draw_developable_edit_entries)
    bpy.types.VIEW3D_MT_object.remove(_draw_developable_object_entries)
    for operator in reversed(_DEVELOPABLE_OPERATORS):
        bpy.utils.unregister_class(operator)
    bpy.types.VIEW3D_MT_object.remove(_draw_sphere_entry)
    bpy.utils.unregister_class(OBJECT_OT_ddg_map_to_sphere)
    for menu in _MENUS:
        menu.remove(_draw_menu_entry)
    for operator in reversed(_CUT_OPERATORS):
        bpy.utils.unregister_class(operator)
    bpy.utils.unregister_class(UV_OT_ddg_conformal_unwrap)
