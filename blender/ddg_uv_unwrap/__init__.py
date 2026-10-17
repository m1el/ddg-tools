# SPDX-License-Identifier: MIT
"""DDG Conformal UV Unwrap: an operator that generates a UV map for the active mesh, and three
that cut a mesh open where it flattens best.

The work is done by the native module ``ddg_tools`` (Rust, bundled as a wheel). See
``blender/README.md`` in the ddg-tools repository.

Credits. UV unwrap: Soliman, Slepcev, Crane, Optimal Cone Singularities for Conformal Flattening
(2018); Sawhney, Crane, Boundary First Flattening (2017). Cut: Sharp, Crane, Variational Surface
Cutting (2018), implemented clean-room from the text of the paper. Map to sphere: Sawhney,
Crane, Boundary First Flattening (2017), Sec. 7.5; Baden, Crane, Kazhdan, Moebius Registration
(2018).
"""

import bpy

from . import core
from .cut_operators import UV_OT_ddg_cut_apply, UV_OT_ddg_cut_flow, UV_OT_ddg_cut_start
from .operator import UV_OT_ddg_conformal_unwrap
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


def unregister():
    bpy.types.VIEW3D_MT_object.remove(_draw_sphere_entry)
    bpy.utils.unregister_class(OBJECT_OT_ddg_map_to_sphere)
    for menu in _MENUS:
        menu.remove(_draw_menu_entry)
    for operator in reversed(_CUT_OPERATORS):
        bpy.utils.unregister_class(operator)
    bpy.utils.unregister_class(UV_OT_ddg_conformal_unwrap)
