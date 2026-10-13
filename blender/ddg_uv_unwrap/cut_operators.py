# SPDX-License-Identifier: MIT
"""The three operators of the variational cut: Start, Flow, Apply.

Methods: Sharp, Crane, Variational Surface Cutting (2018), implemented clean-room from the text
of the paper (``crates/ddg-variational-cuts``). The work is in ``cut.py``.
"""

import bmesh
import bpy
from bpy.props import (BoolProperty, EnumProperty, FloatProperty, FloatVectorProperty,
                       IntProperty, StringProperty)

from . import cut, core

_CITATION = ("Method: Sharp, Crane, Variational Surface Cutting (2018); the implementation is "
             "clean-room, from the text of the paper")


def _poll(cls, context):
    obj = context.active_object
    return obj is not None and obj.type == 'MESH' and obj.mode in {'OBJECT', 'EDIT'}


class _ObjectMode:
    """Leaves Edit mode for the duration of the block (the mesh datablock is only current, and
    only writable, outside Edit mode) and returns to it afterwards."""

    def __init__(self, obj):
        self.obj = obj
        self.edit = obj.mode == 'EDIT'

    def __enter__(self):
        if self.edit:
            bpy.ops.object.mode_set(mode='OBJECT')
        return self

    def __exit__(self, *_exc):
        if self.edit:
            bpy.ops.object.mode_set(mode='EDIT')
        return False


def _edit_selection(obj):
    """``(active_vertex, selected_vertices)`` of the edit mesh, or ``(None, ())`` outside Edit
    mode."""
    if obj.mode != 'EDIT':
        return None, ()
    bm = bmesh.from_edit_mesh(obj.data)
    active = bm.select_history.active
    active_vertex = active.index if isinstance(active, bmesh.types.BMVert) else None
    return active_vertex, tuple(v.index for v in bm.verts if v.select)


class UV_OT_ddg_cut_start(bpy.types.Operator):
    bl_idname = "uv.ddg_cut_start"
    bl_label = "Cut: Start"
    bl_description = (
        "Start a variational cut of the active mesh: store the initial cut as one float "
        "attribute per region on the face corners, and show it as a wire object. "
        "Existing seams are fixed cuts. " + _CITATION)
    bl_options = {'REGISTER', 'UNDO'}

    mode: EnumProperty(
        name="Start From",
        description="Where the first cut comes from",
        items=[
            ('SELECTED', "Selected Faces", "The selected faces are region 0, the rest region 1"),
            ('CIRCLE', "Geodesic Circle",
             "A geodesic disk of the given radius around a vertex is region 0"),
        ],
        default='SELECTED',
    )
    radius: FloatProperty(
        name="Radius",
        description=(
            "Radius of the geodesic disk, in object units (the object's scale is not applied). "
            "Distances are measured along the surface, not across seams. 0: a tenth of the "
            "diagonal of the bounding box"),
        default=0.0, min=0.0, soft_max=100.0, subtype='DISTANCE', unit='LENGTH',
    )
    vertex: IntProperty(
        name="Vertex",
        description=(
            "Index of the vertex the circle is centred on. -1: the active vertex of the edit "
            "mesh (or the only selected vertex)"),
        default=-1, min=-1,
    )
    regions: IntProperty(
        name="Regions",
        description=(
            "Number of regions. The cut flow handles 2 so far; larger values are refused with "
            "a message"),
        default=2, min=2, soft_max=4,
    )

    @classmethod
    def poll(cls, context):
        return _poll(cls, context)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False
        layout.prop(self, "mode")
        if self.mode == 'CIRCLE':
            layout.prop(self, "radius")
            layout.prop(self, "vertex")
        layout.prop(self, "regions")

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        active_vertex, selected_vertices = _edit_selection(obj)
        try:
            with _ObjectMode(obj):
                stats = cut.start_cut(
                    obj, mode=self.mode, radius=self.radius, vertex=self.vertex,
                    n_regions=self.regions, active_vertex=active_vertex,
                    selected_vertices=selected_vertices)
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        for warning in stats["warnings"]:
            self.report({'WARNING'}, warning)
        if self.mode == 'SELECTED':
            what = "the selected faces"
        else:
            what = "a geodesic circle of radius {:.4g} around vertex {}".format(
                stats["radius"], stats["vertex"])
        self.report({'INFO'}, "Cut started from {}: {} regions on {} face corners. "
                    "Run 'Cut: Flow' to move the cut.".format(
                        what, stats["n_regions"], stats["n_corners"]))
        return {'FINISHED'}


_STAGE_ITEMS = {
    "flux_recovery": [('DEFAULT', "From Stage Set", "As the stage set says"),
                      ('NODAL', "Nodal", "The paper's: normal derivative per cut vertex"),
                      ('SMOOTHED', "Smoothed", "Smoothed along the cut before it is divided")],
    "length_step": [('DEFAULT', "From Stage Set", "As the stage set says"),
                    ('ON_CUT', "On the Cut", "The length term is integrated on the cut"),
                    ('ON_MESH', "On the Mesh", "The length term is integrated on the mesh")],
    "extension": [('DEFAULT', "From Stage Set", "As the stage set says"),
                  ('LOCAL', "Local", "Extended from the cut outward"),
                  ('HARMONIC', "Harmonic", "Harmonic extension, as in the paper")],
    "speed_scaling": [('DEFAULT', "From Stage Set", "As the stage set says"),
                      ('MEAN_GRADIENT', "Mean Gradient", "Scaled by the mean gradient norm"),
                      ('UNIT', "Unit", "The speed is the rate, as in the paper")],
    "redistance": [('DEFAULT', "From Stage Set", "As the stage set says"),
                   ('KEEP_ADJACENT', "Keep Adjacent", "Cut-adjacent vertices keep their values"),
                   ('OFF', "Off", "No redistancing")],
    "penalty_rule": [('DEFAULT', "Default", "The gradient rule"),
                     ('GRADIENT', "Gradient", "Full shape derivative of the length terms"),
                     ('PAPER_DIAGONAL', "Paper Diagonal", "The rule printed in the paper")],
}


class UV_OT_ddg_cut_flow(bpy.types.Operator):
    bl_idname = "uv.ddg_cut_flow"
    bl_label = "Cut: Flow"
    bl_description = (
        "Move the cut of the active mesh toward lower flattening distortion. Reads the cut from "
        "the mesh attributes, runs the flow, writes it back and updates the preview; running "
        "it again continues from the stored cut. " + _CITATION)
    bl_options = {'REGISTER', 'UNDO'}

    steps: IntProperty(
        name="Steps",
        description="Number of flow steps (at most, if the stop rule is on)",
        default=10, min=1, soft_max=1000,
    )
    use_stop_rule: BoolProperty(
        name="Stop Early",
        description=(
            "Stop before the given number of steps when the cut stops moving or the energy stops "
            "falling, or when the cut vanishes. Off: exactly the given number of steps"),
        default=True,
    )
    energy: EnumProperty(
        name="Energy",
        description="The distortion energy of each region",
        items=[
            ('DIRICHLET', "Dirichlet", "Integral of the squared gradient of the log scale factor"),
            ('SCALED_DIRICHLET', "Scaled Dirichlet",
             "Dirichlet energy weighted by the curvature-based feature size (paper Eq. 17)"),
            ('HENCKY', "Hencky", "Integral of the squared log scale factor (paper Eq. 19)"),
        ],
        default='HENCKY',
    )
    length_mode: EnumProperty(
        name="Length Weight",
        description="How the weight of the cut length is chosen",
        items=[
            ('RIPPLE', "From Ripple Wavelength",
             "Set from the roughening rate of the current cut, so that ripples longer than the "
             "wavelength grow and shorter ones decay"),
            ('EXPLICIT', "Explicit", "A fixed value"),
        ],
        default='RIPPLE',
    )
    ripple_wavelength: FloatProperty(
        name="Ripple Wavelength",
        description=(
            "Shortest ripple of the cut that the distortion may amplify, as a fraction of "
            "the square root of the surface area. The cut_flow example of the library uses 0.25. "
            "It must span several mesh edges, or the flow amplifies the noise of the mesh"),
        default=0.25, min=0.0, soft_max=1.0,
    )
    length_weight: FloatProperty(
        name="Weight",
        description=(
            "The weight of the cut length, on the surface scaled to unit area (the status "
            "line of 'Ripple Wavelength' runs reports the value it chose)"),
        default=1.0, min=0.0, soft_max=100.0,
    )
    stages: EnumProperty(
        name="Stages",
        description="The numerical stages of a step",
        items=[
            ('DEFAULT', "Default", "The library's defaults (the robust variants)"),
            ('PAPER', "Paper", "The stages as printed in the paper"),
        ],
        default='DEFAULT',
    )
    penalty_group: StringProperty(
        name="Penalty Group",
        description=(
            "Vertex group whose weights (0 to 1) are a penalty on cutting there (paper Sec. "
            "7.1): the penalty is Strength times the length weight times the weight"),
        default="",
    )
    penalty_strength: FloatProperty(
        name="Strength",
        description="Penalty of a vertex with weight 1, in units of the length weight",
        default=1.0, min=0.0, soft_max=10.0,
    )
    hide_from_camera: BoolProperty(
        name="Hide From Camera",
        description=(
            "Penalize cutting where the active camera sees the surface (paper Sec. 7.1): "
            "each vertex gets a penalty if a ray from the camera reaches it, times how much "
            "it faces the camera"),
        default=False,
    )
    hide_strength: FloatProperty(
        name="Strength",
        description="Penalty of a visible vertex facing the camera, in units of the length weight",
        default=1.0, min=0.0, soft_max=10.0,
    )
    distortion_weights: FloatVectorProperty(
        name="Distortion Weights",
        description=(
            "Weight of the distortion of region 0 and region 1. 0 for one region lets the "
            "cut ignore its distortion: a one-sided patch (paper Fig. 14)"),
        size=2, default=(1.0, 1.0), min=0.0, soft_max=10.0,
    )
    area_weights: FloatVectorProperty(
        name="Area Weights",
        description=(
            "Weight of the area of region 0 and region 1, on the surface scaled to unit area. "
            "Positive values make a region shrink, negative values make it grow"),
        size=2, default=(0.0, 0.0), soft_min=-10.0, soft_max=10.0,
    )
    symmetry: EnumProperty(
        name="Mirror Symmetry",
        description="Keep the cut symmetric about the plane through the object's origin",
        items=[
            ('NONE', "None", "No symmetry"),
            ('X', "X", "Mirror about the YZ plane"),
            ('Y', "Y", "Mirror about the XZ plane"),
            ('Z', "Z", "Mirror about the XY plane"),
        ],
        default='NONE',
    )
    symmetry_tolerance: FloatProperty(
        name="Tolerance",
        description=(
            "How far a vertex may be from the mirror image of its partner, as a fraction of "
            "the diagonal of the bounding box. Vertices without a partner are reported"),
        default=1e-4, min=0.0, soft_max=0.01, precision=5,
    )
    use_uniform_curvature: BoolProperty(
        name="Uniform Curvature",
        description="Replace the Gaussian curvature by the same curvature density everywhere "
                    "(paper Sec. 7.4)",
        default=False,
    )
    uniform_curvature: FloatProperty(
        name="Curvature",
        description="Curvature density per unit area of the surface scaled to unit area",
        default=0.0, soft_min=-100.0, soft_max=100.0,
    )

    # ----- advanced options: the library's defaults -----
    feature_r_eps: FloatProperty(
        name="Feature Size Cap", description="The cap on the feature size of the scaled "
        "Dirichlet energy, on the surface scaled to unit area", default=1.0, min=1e-6,
        soft_max=10.0)
    feature_smoothing: FloatProperty(
        name="Feature Smoothing", description="The smoothing constant of the feature size "
        "(the paper uses 0.02)", default=0.02, min=0.0, soft_max=1.0)
    step_rule: EnumProperty(
        name="Step Size Rule",
        description="How the time step is chosen",
        items=[('ADAPTIVE', "Adaptive", "Scaled by the roughening rate and the cut motion"),
               ('CAPPED', "Capped", "At most a maximum, and at most a fraction of an edge"),
               ('FIXED', "Fixed", "The same step every time")],
        default='ADAPTIVE',
    )
    time_scale: FloatProperty(name="Time Scale", description="Step in units of 1/roughening "
                              "rate (adaptive)", default=0.05, min=1e-6, soft_max=1.0)
    max_motion: FloatProperty(name="Max Motion", description="Largest motion of the cut in a "
                              "step, as a fraction of the shortest crossed edge", default=0.25,
                              min=1e-6, soft_max=1.0)
    tau_max: FloatProperty(name="Max Step", description="Largest time step", default=1.0,
                           min=1e-9, soft_max=10.0)
    tau: FloatProperty(name="Step", description="The time step (fixed rule)", default=1e-3,
                       min=1e-9, soft_max=1.0)
    min_crossing_offset: FloatProperty(
        name="Crossing Offset", description="Crossings closer than this fraction of an edge to "
        "a vertex are moved to that distance", default=1e-3, min=0.0, max=0.499, precision=5)
    motion_tolerance: FloatProperty(
        name="Motion Tolerance", description="Stop when the cut moves less than this in a step "
        "(on the surface scaled to unit area)", default=1e-5, min=0.0, precision=7)
    energy_tolerance: FloatProperty(
        name="Energy Tolerance", description="Stop when the energy fell by less than this "
        "fraction over the window", default=1e-4, min=0.0, precision=6)
    energy_window: IntProperty(name="Energy Window", description="Steps the energy is compared "
                               "over", default=25, min=1)
    flux_recovery: EnumProperty(name="Normal Derivative", items=_STAGE_ITEMS["flux_recovery"],
                                default='DEFAULT')
    length_step: EnumProperty(name="Length Step", items=_STAGE_ITEMS["length_step"],
                              default='DEFAULT')
    extension: EnumProperty(name="Extension", items=_STAGE_ITEMS["extension"], default='DEFAULT')
    speed_scaling: EnumProperty(name="Speed Scaling", items=_STAGE_ITEMS["speed_scaling"],
                                default='DEFAULT')
    redistance: EnumProperty(name="Redistancing", items=_STAGE_ITEMS["redistance"],
                             default='DEFAULT')
    penalty_rule: EnumProperty(name="Penalty Rule", items=_STAGE_ITEMS["penalty_rule"],
                               default='DEFAULT')

    show_advanced: BoolProperty(
        name="Advanced Options", description="Show the advanced options", default=False,
        options={'SKIP_SAVE'})

    @classmethod
    def poll(cls, context):
        return _poll(cls, context)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False
        layout.prop(self, "steps")
        layout.prop(self, "use_stop_rule")
        layout.prop(self, "energy")
        layout.prop(self, "length_mode")
        layout.prop(self, "ripple_wavelength" if self.length_mode == 'RIPPLE'
                    else "length_weight")
        layout.prop(self, "stages")
        obj = context.active_object
        box = layout.box()
        box.label(text="Penalty")
        if obj is not None and obj.type == 'MESH':
            box.prop_search(self, "penalty_group", obj, "vertex_groups")
        else:
            box.prop(self, "penalty_group")
        box.prop(self, "penalty_strength", text="Group Strength")
        box.prop(self, "hide_from_camera")
        if self.hide_from_camera:
            box.prop(self, "hide_strength", text="Camera Strength")
        layout.prop(self, "distortion_weights")
        layout.prop(self, "area_weights")
        layout.prop(self, "symmetry")
        if self.symmetry != 'NONE':
            layout.prop(self, "symmetry_tolerance")
        layout.prop(self, "use_uniform_curvature")
        if self.use_uniform_curvature:
            layout.prop(self, "uniform_curvature")
        layout.prop(self, "show_advanced")
        body = layout.box() if self.show_advanced else None
        if body is not None:
            body.prop(self, "feature_r_eps")
            body.prop(self, "feature_smoothing")
            body.prop(self, "step_rule")
            if self.step_rule == 'ADAPTIVE':
                body.prop(self, "time_scale")
            if self.step_rule != 'FIXED':
                body.prop(self, "max_motion")
                body.prop(self, "tau_max")
            else:
                body.prop(self, "tau")
            body.prop(self, "min_crossing_offset")
            body.prop(self, "motion_tolerance")
            body.prop(self, "energy_tolerance")
            body.prop(self, "energy_window")
            for name in ("flux_recovery", "length_step", "extension", "speed_scaling",
                         "redistance", "penalty_rule"):
                body.prop(self, name)

    def settings(self):
        keys = (
            "steps", "use_stop_rule", "energy", "length_mode", "ripple_wavelength",
            "length_weight", "stages", "penalty_group", "penalty_strength", "hide_from_camera",
            "hide_strength", "distortion_weights", "area_weights", "symmetry",
            "symmetry_tolerance", "use_uniform_curvature", "uniform_curvature", "feature_r_eps",
            "feature_smoothing", "step_rule", "time_scale", "max_motion", "tau_max", "tau",
            "min_crossing_offset", "motion_tolerance", "energy_tolerance", "energy_window",
            "flux_recovery", "length_step", "extension", "speed_scaling", "redistance",
            "penalty_rule")
        return {key: getattr(self, key) for key in keys}

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        try:
            with _ObjectMode(obj):
                report = cut.run_flow(obj, context, self.settings())
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        for warning in report["warnings"]:
            self.report({'WARNING'}, warning)
        self.report({'INFO'}, cut.flow_summary(report))
        return {'FINISHED'}


class UV_OT_ddg_cut_apply(bpy.types.Operator):
    bl_idname = "uv.ddg_cut_apply"
    bl_label = "Cut: Apply"
    bl_description = (
        "Insert the cut into the mesh of the active object and mark it as seams: crossed faces "
        "are split along the cut and other data (UV maps, vertex groups, materials, attributes) "
        "is interpolated or inherited. Removes the preview and the cut attributes. Existing "
        "seams stay. " + _CITATION)
    bl_options = {'REGISTER', 'UNDO'}

    unwrap_after: BoolProperty(
        name="Unwrap Afterwards",
        description="Run Conformal Unwrap (Optimal Cones) on the result, with its default "
                    "settings and the seams as they are now",
        default=False,
    )

    @classmethod
    def poll(cls, context):
        return _poll(cls, context)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False
        layout.prop(self, "unwrap_after")

    def execute(self, context):
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "The active object is not a mesh")
            return {'CANCELLED'}
        unwrap_error = None
        try:
            with _ObjectMode(obj):
                stats = cut.apply_cut(obj)
                if self.unwrap_after:
                    try:
                        bpy.ops.uv.ddg_conformal_unwrap(use_seams=True)
                    except RuntimeError as ex:
                        unwrap_error = str(ex).strip()
        except core.UnwrapFailed as ex:
            self.report({'ERROR'}, str(ex))
            return {'CANCELLED'}
        for warning in stats["warnings"]:
            self.report({'WARNING'}, warning)
        self.report({'INFO'}, (
            "Cut applied: {} vertices and {} faces added, {} edges marked as seams "
            "({} -> {} vertices, {} -> {} faces)").format(
            stats["n_vertices"] - stats["n_vertices_before"],
            stats["n_polygons"] - stats["n_polygons_before"], stats["n_seams_marked"],
            stats["n_vertices_before"], stats["n_vertices"], stats["n_polygons_before"],
            stats["n_polygons"]))
        if unwrap_error is not None:
            self.report({'ERROR'}, "The cut was applied, but the unwrap failed: " + unwrap_error)
        return {'FINISHED'}
