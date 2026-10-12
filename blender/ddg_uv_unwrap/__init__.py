# SPDX-License-Identifier: MIT
"""DDG Conformal UV Unwrap: one operator that generates a UV map for the active mesh.

The work is done by the native module ``ddg_tools`` (Rust, bundled as a wheel). See
``blender/README.md`` in the ddg-tools repository.
"""

import bpy

from . import core
from .operator import UV_OT_ddg_conformal_unwrap


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


# The UV menu of the 3D viewport in Edit mode, and the Object menu.
_MENUS = (bpy.types.VIEW3D_MT_uv_map, bpy.types.VIEW3D_MT_object)


def register():
    bpy.utils.register_class(UV_OT_ddg_conformal_unwrap)
    for menu in _MENUS:
        menu.append(_draw_menu_entry)


def unregister():
    for menu in _MENUS:
        menu.remove(_draw_menu_entry)
    bpy.utils.unregister_class(UV_OT_ddg_conformal_unwrap)
