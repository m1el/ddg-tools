# SPDX-License-Identifier: MIT
"""The one operator of this add-on: generate a UV map for the active mesh."""

import bpy
from bpy.props import BoolProperty, FloatProperty, StringProperty

from . import core


class UV_OT_ddg_conformal_unwrap(bpy.types.Operator):
    bl_idname = "uv.ddg_conformal_unwrap"
    bl_label = "Conformal Unwrap (Optimal Cones)"
    bl_description = (
        "Generate a UV map for the whole active mesh: place cone singularities that minimize "
        "area distortion, cut the surface open through them, and flatten it conformally. "
        "The mesh itself is not changed"
    )
    bl_options = {'REGISTER', 'UNDO'}

    # `lambda` is a Python keyword, hence the underscore.
    lambda_: FloatProperty(
        name="Lambda",
        description=(
            "Trade-off between area distortion and cones. "
            "Larger values give fewer cones and more distortion"
        ),
        default=0.5,
        min=0.0,
        soft_min=0.05,
        soft_max=2.0,
        step=5,
        precision=3,
    )
    use_seams: BoolProperty(
        name="Cut Marked Seams",
        description="Cut the surface along the edges marked as seams before placing cones",
        default=True,
    )
    pinned_group: StringProperty(
        name="Pinned Cones",
        description=(
            "Vertex group (weights above 0) of vertices where a cone may be placed at no cost"
        ),
        default="",
    )
    forbidden_group: StringProperty(
        name="No Cones",
        description="Vertex group (weights above 0) of vertices where no cone may be placed",
        default="",
    )
    positive_cones_only: BoolProperty(
        name="Positive Cones Only",
        description="Allow only cones with a positive angle (cone tips, no saddles)",
        default=False,
    )
    mark_seams: BoolProperty(
        name="Mark Cuts as Seams",
        description="Mark the edges the surface was cut along as seams",
        default=False,
    )

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False
        if core.is_placeholder():
            column = layout.column(align=True)
            column.label(text="Placeholder build", icon='ERROR')
            column.label(text="Planar projection, not the conformal unwrap")
        layout.prop(self, "lambda_")
        layout.prop(self, "positive_cones_only")
        layout.prop(self, "use_seams")
        obj = context.active_object
        if obj is not None and obj.type == 'MESH':
            layout.prop_search(self, "pinned_group", obj, "vertex_groups")
            layout.prop_search(self, "forbidden_group", obj, "vertex_groups")
        else:
            layout.prop(self, "pinned_group")
            layout.prop(self, "forbidden_group")
        layout.prop(self, "mark_seams")

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}

        # The mesh datablock is only current, and only writable, outside Edit mode.
        was_edit_mode = obj.mode == 'EDIT'
        if was_edit_mode:
            bpy.ops.object.mode_set(mode='OBJECT')
        try:
            stats = core.unwrap_object(
                obj,
                lambda_=self.lambda_,
                use_seams=self.use_seams,
                pinned_group=self.pinned_group,
                forbidden_group=self.forbidden_group,
                positive_cones_only=self.positive_cones_only,
                mark_cuts_as_seams=self.mark_seams,
            )
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        finally:
            if was_edit_mode:
                bpy.ops.object.mode_set(mode='EDIT')

        report = stats["report"]
        placeholder_warning = bool(report.get("placeholder"))
        for warning in report.get("warnings", ()):
            # The stand-in pipeline's own warning is already part of the summary line.
            if not (placeholder_warning and "stand-in" in warning):
                self.report({'WARNING'}, warning)
        if stats["cut_edges_not_found"]:
            self.report(
                {'WARNING'},
                "{} cut edges are not edges of the mesh and could not be marked as seams".format(
                    stats["cut_edges_not_found"]),
            )
        self.report({'WARNING'} if placeholder_warning else {'INFO'}, core.summary(stats))
        return {'FINISHED'}
