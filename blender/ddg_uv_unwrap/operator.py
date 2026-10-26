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
        "The mesh itself is not changed. "
        "Methods: Soliman, Slepcev, Crane, Optimal Cone Singularities for Conformal Flattening "
        "(2018); Sawhney, Crane, Boundary First Flattening (2017)"
    )
    bl_options = {'REGISTER', 'UNDO'}

    # `lambda` is a Python keyword, hence the underscore.
    lambda_: FloatProperty(
        name="Lambda",
        description=(
            "Trade-off between area distortion and cones. "
            "Larger values give fewer cones and more distortion"
        ),
        default=0.1,
        min=0.0,
        soft_min=0.01,
        soft_max=1.0,
        step=1,
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
        description=(
            "Allow only cones with a positive angle (cone tips, no saddles). "
            "Negative cones lower distortion but make the UV map overlap itself. "
            "A closed surface of genus 2 or more cannot be flattened with positive cones "
            "only: it gets cones of either sign, with a warning"
        ),
        default=True,
    )
    optimize_boundary: BoolProperty(
        name="Optimize Boundary",
        description=(
            "Let the solver stretch or shrink the border of the mesh and of every seam to "
            "lower the distortion. The two sides of a seam can then have different lengths "
            "(a texture is stretched differently on them), and on some meshes the result "
            "overlaps itself or contains mirrored triangles; the warnings say so. "
            "Off keeps every border edge at its length"
        ),
        default=False,
    )
    split_overlaps: BoolProperty(
        name="Split Overlapping Charts",
        description=(
            "A chart whose flattening lies on top of itself is cut in two where the layers "
            "meet and flattened again; a cut is kept only if it lowers the overlap and mirrors "
            "no triangle. Adds seams; off keeps the charts as the cones and seams made them"
        ),
        default=False,
    )
    cone_feedback: BoolProperty(
        name="Cone Feedback",
        description=(
            "Compare the cone solver's scale factor with the one the flattening realizes, "
            "ignore the distortion at vertices where they are far apart (the tips of thin "
            "spikes) and solve again. Costs a comparison where nothing is found, and a solve "
            "and a flattening per extra round where something is"
        ),
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
        layout.prop(self, "optimize_boundary")
        layout.prop(self, "split_overlaps")
        layout.prop(self, "cone_feedback")
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
                optimize_boundary=self.optimize_boundary,
                split_overlaps=self.split_overlaps,
                cone_feedback=self.cone_feedback,
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
            # A placeholder build's own warning is already part of the summary line. (The
            # shipped native module is never a placeholder; see core.is_placeholder.)
            if not (placeholder_warning and "stand-in" in warning):
                self.report({'WARNING'}, warning)
        if stats["cut_edges_not_found"]:
            self.report(
                {'WARNING'},
                "{} cut edges are not edges of the mesh and could not be marked as seams".format(
                    stats["cut_edges_not_found"]),
            )
        note = stats.get("triangulation")
        if note and (note["stale"] or note["reason"]):
            self.report({'WARNING'}, (
                "Used Blender's triangulation for {} of the {} faces whose triangulation was "
                "stored when the cut was applied{}").format(
                note["stale"], note["used"] + note["stale"],
                ": " + note["reason"] if note["reason"] else " (they were edited since)"))
        self.report({'WARNING'} if placeholder_warning else {'INFO'}, core.summary(stats))
        return {'FINISHED'}
