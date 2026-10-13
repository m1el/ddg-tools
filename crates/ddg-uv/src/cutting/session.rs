//! The cut session: the cut flow on a polygon mesh as Blender provides it.

use std::collections::HashMap;

use ddg_core::Geometry;
use ddg_core::distance::upwind_distance;
use ddg_core::topology::{NONE, TriMesh};
use ddg_variational_cuts::{
    CurvatureSource, Cut, CutProblem, CurveEnd, Energy, EvalOptions, Evaluation, Flow,
    FlowOptions, RegionFlow, RunReport, ScaleWeight, StepReport, StepSize, StopRule, init,
};

use super::describe::{
    CutOnPolygons, EdgePoint, InteriorPoint, JunctionPoint, PathEnd, PolygonPath,
};
use super::error::{CuttingError, internal};
use crate::input::{MeshInput, PolygonMesh};
use crate::options::{Degenerate, Orientation, UnwrapOptions};
use crate::prepare::prepare;

/// The default `δ` of the cut: crossings closer than this fraction of an edge to a vertex are
/// moved to that distance (the default of the flow, `EvalOptions::default`).
pub const DEFAULT_CROSSING_OFFSET: f64 = 1e-3;

/// How a session starts, or restarts ([`Session::restart`]). All at the level of the caller's
/// polygons.
#[derive(Debug, Clone)]
pub enum Start {
    /// The geodesic disk of the given radius (in the caller's length unit) around a vertex of
    /// the polygon mesh: region 0 is the disk. The distance is measured on the surface cut open
    /// along the seams, from every copy of the vertex: across a seam it does not reach.
    GeodesicCircle {
        /// The vertex of the polygon mesh.
        vertex: u32,
        /// The radius, positive and finite.
        radius: f64,
    },
    /// One label per polygon, each below the number of regions. The cut follows the interface
    /// between the labels as `ddg_variational_cuts::init::from_face_labels` does for the
    /// triangles of the polygons: a vertex takes the label that covers more of its angle.
    PolygonLabels(Vec<u32>),
    /// One region function per region and per polygon corner, as [`Session::state`] returns them.
    Corners(Vec<Vec<f64>>),
    /// The vertex normals of the flow mesh clustered into `n_regions` clusters
    /// (`ddg_variational_cuts::init::normal_clusters`, at most `max_iterations` rounds), the
    /// clusters as regions: the start of the atlas of the paper (Sec. 7). Works for any number
    /// of regions; the clusters are speckled and the flow is meant to clean them up. The
    /// clustering runs on the flow mesh, so a vertex that a seam splits may get a different
    /// cluster on each side.
    NormalClusters {
        /// The most rounds of the clustering.
        max_iterations: usize,
    },
}

/// Which flow a session runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FlowKind {
    /// The two-region flow for two regions, [`RegionFlow`] for more. The default of
    /// [`Session::new`].
    #[default]
    Auto,
    /// The two-region flow (`Flow`); an error for more than two regions.
    TwoRegion,
    /// [`RegionFlow`], also for two regions (where it takes the same steps as the two-region
    /// flow, to rounding, see the core crate).
    Regions,
}

/// A line of the cut for display: positions in space.
#[derive(Debug, Clone, PartialEq)]
pub struct PreviewPolyline {
    /// The positions, in order. A closed polyline does not repeat its first point.
    pub positions: Vec<[f64; 3]>,
    /// Whether the last position is joined to the first. An open polyline ends on the boundary
    /// of the surface or on a seam (the two sides of a seam give separate polylines).
    pub closed: bool,
    /// The two regions the line separates, `[lower id, higher id]` (`[0, 1]` for two regions).
    pub regions: [u32; 2],
}

/// The maps between the polygon mesh and the triangle mesh the flow runs on.
#[derive(Debug, Clone)]
struct Maps {
    n_regions: usize,
    n_vertices: usize,
    n_corners: usize,
    n_polygons: usize,
    /// Flow vertex to vertex of the polygon mesh.
    vertex_input: Vec<u32>,
    /// Flow face to triangle of the supplied triangulation.
    face_triangle: Vec<u32>,
    /// Triangle to polygon.
    triangle_polygon: Vec<u32>,
    /// Polygon corner to flow vertex.
    corner_vertex: Vec<u32>,
    /// First corner of every flow vertex.
    first_corner: Vec<u32>,
    /// Vertex of the polygon mesh to its flow vertices (compressed rows).
    copy_start: Vec<u32>,
    copies: Vec<u32>,
    /// Flow edge to the two polygon corners of its first face, as `(c, d)`.
    edge_corners: Vec<(u32, u32)>,
    /// Whether a flow edge is an edge of a polygon (and not a diagonal).
    edge_is_polygon: Vec<bool>,
    /// The seams that were cut (given seams and extensions of one-edge seams), `[min, max]`.
    used_seams: Vec<[u32; 2]>,
}

impl Maps {
    fn copies_of(&self, v: u32) -> &[u32] {
        &self.copies[self.copy_start[v as usize] as usize..self.copy_start[v as usize + 1] as usize]
    }

    /// The level-set function from the region functions of the corners.
    fn phi_from_regions<R: AsRef<[f64]>>(&self, regions: &[R]) -> Result<Vec<f64>, CuttingError> {
        self.check_regions(regions)?;
        // Two regions: the functions are φ and −φ; the level-set function is half the
        // difference.
        let (f0, f1) = (regions[0].as_ref(), regions[1].as_ref());
        Ok((0..self.vertex_input.len())
            .map(|fv| {
                let c = self.first_corner[fv] as usize;
                0.5 * (f0[c] - f1[c])
            })
            .collect())
    }

    /// The level-set functions of the flow mesh (one per region) from the region functions of
    /// the corners.
    fn phis_from_regions<R: AsRef<[f64]>>(
        &self,
        regions: &[R],
    ) -> Result<Vec<Vec<f64>>, CuttingError> {
        self.check_regions(regions)?;
        Ok(regions
            .iter()
            .map(|f| {
                let f = f.as_ref();
                (0..self.vertex_input.len())
                    .map(|fv| f[self.first_corner[fv] as usize])
                    .collect()
            })
            .collect())
    }

    /// The region functions of the corners from the level-set functions of the flow mesh.
    fn regions_from_phis(&self, phis: &[Vec<f64>]) -> Vec<Vec<f64>> {
        phis.iter()
            .map(|phi| self.corner_vertex.iter().map(|&v| phi[v as usize]).collect())
            .collect()
    }

    /// Checks the region functions: count, length, finite, equal across corners of a vertex.
    fn check_regions<R: AsRef<[f64]>>(&self, regions: &[R]) -> Result<(), CuttingError> {
        if regions.len() != self.n_regions {
            return Err(CuttingError::WrongLength {
                what: "region functions",
                per: "region",
                expected: self.n_regions,
                got: regions.len(),
            });
        }
        for (k, f) in regions.iter().enumerate() {
            let f = f.as_ref();
            if f.len() != self.n_corners {
                return Err(CuttingError::WrongLength {
                    what: "region function",
                    per: "polygon corner",
                    expected: self.n_corners,
                    got: f.len(),
                });
            }
            if let Some(c) = f.iter().position(|x| !x.is_finite()) {
                return Err(CuttingError::NonFinite {
                    what: "region function",
                    region: k,
                    position: c,
                    value: f[c],
                });
            }
            for (c, &fv) in self.corner_vertex.iter().enumerate() {
                let first = self.first_corner[fv as usize] as usize;
                if f[c] != f[first] {
                    return Err(CuttingError::CornersDisagree {
                        vertex: self.vertex_input[fv as usize],
                        region: k,
                        corners: [first as u32, c as u32],
                        values: [f[first], f[c]],
                    });
                }
            }
        }
        Ok(())
    }

    /// The region functions of the corners from the level-set function.
    fn regions_from_phi(&self, phi: &[f64]) -> Vec<Vec<f64>> {
        let at_corner: Vec<f64> = self
            .corner_vertex
            .iter()
            .map(|&v| phi[v as usize])
            .collect();
        let negated: Vec<f64> = at_corner.iter().map(|x| -x).collect();
        vec![at_corner, negated]
    }
}

/// A session of the cut flow on a polygon mesh.
///
/// It owns a copy of the polygon mesh with its triangulation, the triangle mesh the flow runs on
/// (the triangulation cut open along the seams, so the seams are boundary for the flow), the
/// maps between the two, and a [`Flow`]. See the [module](super) documentation.
///
/// **Policies fixed by the session.** The triangulated mesh must be consistently oriented and
/// without zero-area triangles: [`Orientation::Reject`] and [`Degenerate::Reject`] are used for
/// the validation, because the flow works on the true geometry and a repaired or mollified mesh
/// would no longer be the caller's. The seams are handled as `prepare` does: a seam that is not
/// an edge is ignored with a warning ([`Session::warnings`]), a seam of one edge between two
/// interior vertices is extended by one edge.
#[derive(Debug)]
pub struct Session {
    mesh: PolygonMesh,
    seams: Vec<[u32; 2]>,
    pub(super) pinned: Vec<u32>,
    pub(super) prescribed: Vec<ddg_core::io::Cone>,
    pub(super) forbidden: Vec<u32>,
    maps: Maps,
    flow_mesh: TriMesh,
    flow_geometry: Geometry,
    flow_positions: Vec<[f64; 3]>,
    /// The two-region flow; with `many` it is an unused one on the same mesh (see
    /// [`Session::flow`]).
    flow: Flow,
    many: Option<RegionFlow>,
    warnings: Vec<String>,
}

impl Session {
    /// Builds a session from a polygon mesh with its triangulation and seams, with
    /// `n_regions` regions, started as `start` says.
    ///
    /// # Errors
    /// [`CuttingError::TooFewRegions`], [`CuttingError::RegionsNotAvailable`] for more than two
    /// regions, [`CuttingError::Mesh`] (the mesh is refused; the message names polygons and
    /// vertices), and the errors of [`Session::restart`].
    pub fn new(
        input: &MeshInput<'_>,
        n_regions: usize,
        start: &Start,
    ) -> Result<Session, CuttingError> {
        Session::new_with(input, n_regions, start, FlowKind::Auto)
    }

    /// [`Session::new`] with the choice of the flow ([`FlowKind`]).
    ///
    /// # Errors
    /// As [`Session::new`]; [`CuttingError::InvalidOption`] for [`FlowKind::TwoRegion`] with
    /// more than two regions, and for a geodesic circle with more than two regions.
    pub fn new_with(
        input: &MeshInput<'_>,
        n_regions: usize,
        start: &Start,
        kind: FlowKind,
    ) -> Result<Session, CuttingError> {
        if n_regions < 2 {
            return Err(CuttingError::TooFewRegions { got: n_regions });
        }
        let many_flow = match kind {
            FlowKind::Auto => n_regions > 2,
            FlowKind::Regions => true,
            FlowKind::TwoRegion => {
                if n_regions > 2 {
                    return Err(CuttingError::InvalidOption {
                        name: "kind",
                        reason: format!("the two-region flow cannot run {n_regions} regions"),
                    });
                }
                false
            }
        };
        let options = UnwrapOptions {
            orientation: Orientation::Reject,
            degenerate: Degenerate::Reject,
            ..UnwrapOptions::default()
        };
        let prepared = prepare(input, &options)?;
        let n_vertices = input.positions.len();
        let n_corners = input.n_corners();
        let offsets = input.polygon_offsets;

        // The flow mesh: the pieces side by side.
        let mut faces: Vec<[u32; 3]> = Vec::new();
        let mut face_triangle: Vec<u32> = Vec::new();
        let mut vertex_input: Vec<u32> = Vec::new();
        for piece in &prepared.pieces {
            let shift = vertex_input.len() as u32;
            faces.extend(piece.mesh.faces().iter().map(|f| f.map(|v| v + shift)));
            face_triangle.extend_from_slice(&piece.face_triangle);
            vertex_input.extend_from_slice(&piece.vertex_input);
        }
        let flow_mesh = TriMesh::new(vertex_input.len(), faces)
            .map_err(|e| internal(format!("the flow mesh was refused: {e}")))?;
        let flow_positions: Vec<[f64; 3]> = vertex_input
            .iter()
            .map(|&v| input.positions[v as usize])
            .collect();
        let flow_geometry = Geometry::from_positions(&flow_mesh, &flow_positions)
            .map_err(|e| internal(format!("the geometry of the flow mesh was refused: {e}")))?;

        // Corner <-> flow vertex.
        let mut corner_polygon = vec![0u32; n_corners];
        for p in 0..input.n_polygons() {
            corner_polygon[offsets[p] as usize..offsets[p + 1] as usize].fill(p as u32);
        }
        let triangle_polygon: Vec<u32> = input
            .triangles
            .iter()
            .map(|t| corner_polygon[t[0] as usize])
            .collect();
        let mut corner_vertex = vec![NONE; n_corners];
        for (f, face) in flow_mesh.faces().iter().enumerate() {
            let corners = prepared.triangle_corners[face_triangle[f] as usize];
            for k in 0..3 {
                let slot = &mut corner_vertex[corners[k] as usize];
                if *slot != NONE && *slot != face[k] {
                    return Err(internal(format!(
                        "polygon corner {} belongs to two vertices of the flow mesh",
                        corners[k]
                    )));
                }
                *slot = face[k];
            }
        }
        if let Some(c) = corner_vertex.iter().position(|&v| v == NONE) {
            return Err(internal(format!("polygon corner {c} is in no triangle")));
        }
        let mut first_corner = vec![NONE; vertex_input.len()];
        for (c, &fv) in corner_vertex.iter().enumerate() {
            if first_corner[fv as usize] == NONE {
                first_corner[fv as usize] = c as u32;
            }
        }
        let mut copy_start = vec![0u32; n_vertices + 1];
        for &v in &vertex_input {
            copy_start[v as usize + 1] += 1;
        }
        for v in 0..n_vertices {
            copy_start[v + 1] += copy_start[v];
        }
        let mut copies = vec![0u32; vertex_input.len()];
        let mut fill = copy_start.clone();
        for (fv, &v) in vertex_input.iter().enumerate() {
            copies[fill[v as usize] as usize] = fv as u32;
            fill[v as usize] += 1;
        }

        // Flow edges: polygon corners of the first face, polygon edge or diagonal.
        let next_corner = |c: u32| {
            let p = corner_polygon[c as usize] as usize;
            if c + 1 == offsets[p + 1] {
                offsets[p]
            } else {
                c + 1
            }
        };
        let mut edge_corners = vec![(NONE, NONE); flow_mesh.n_edges()];
        let mut edge_is_polygon = vec![false; flow_mesh.n_edges()];
        for (f, &triangle) in face_triangle.iter().enumerate() {
            let corners = prepared.triangle_corners[triangle as usize];
            for k in 0..3 {
                let e = flow_mesh.edge(3 * f + k);
                if edge_corners[e].0 == NONE {
                    let (c, d) = (corners[k], corners[(k + 1) % 3]);
                    edge_corners[e] = (c, d);
                    edge_is_polygon[e] = next_corner(c) == d || next_corner(d) == c;
                }
            }
        }

        let maps = Maps {
            n_regions,
            n_vertices,
            n_corners,
            n_polygons: input.n_polygons(),
            vertex_input,
            face_triangle,
            triangle_polygon,
            corner_vertex,
            first_corner,
            copy_start,
            copies,
            edge_corners,
            edge_is_polygon,
            used_seams: prepared.seam_cut_edges.clone(),
        };
        let mesh = PolygonMesh {
            positions: input.positions.to_vec(),
            corner_vertices: input.corner_vertices.to_vec(),
            polygon_offsets: offsets.to_vec(),
            triangles: input.triangles.to_vec(),
        };
        let level = start_state(
            &maps,
            &flow_mesh,
            &flow_geometry,
            &flow_positions,
            start,
            many_flow,
        )?;
        let (flow, many) = if many_flow {
            let phis = level.into_many();
            let aux = Flow::new(&flow_mesh, &flow_geometry, &phis[0])?;
            (aux, Some(RegionFlow::new(&flow_mesh, &flow_geometry, &phis)?))
        } else {
            let phi = level.into_two();
            (Flow::new(&flow_mesh, &flow_geometry, &phi)?, None)
        };
        let mut seams: Vec<[u32; 2]> = input
            .seam_edges
            .iter()
            .map(|&[a, b]| [a.min(b), a.max(b)])
            .collect();
        seams.sort_unstable();
        seams.dedup();
        Ok(Session {
            mesh,
            seams,
            pinned: input.pinned_vertices.to_vec(),
            prescribed: input.prescribed_cones.to_vec(),
            forbidden: input.forbidden_vertices.to_vec(),
            maps,
            flow_mesh,
            flow_geometry,
            flow_positions,
            flow,
            many,
            warnings: prepared.warnings.clone(),
        })
    }

    /// Replaces the state by a new start (and resets the flow time and step count).
    ///
    /// # Errors
    /// [`CuttingError::VertexOutOfRange`] for a vertex that is not used by any polygon,
    /// [`CuttingError::InvalidOption`] for a radius that is not positive and finite,
    /// [`CuttingError::LabelOutOfRange`], and for corner functions
    /// [`CuttingError::WrongLength`], [`CuttingError::NonFinite`] and
    /// [`CuttingError::CornersDisagree`] (the vertex is named).
    pub fn restart(&mut self, start: &Start) -> Result<(), CuttingError> {
        let level = start_state(
            &self.maps,
            &self.flow_mesh,
            &self.flow_geometry,
            &self.flow_positions,
            start,
            self.many.is_some(),
        )?;
        match &mut self.many {
            Some(f) => f.set_level_sets(&level.into_many())?,
            None => self.flow.set_level_set(&level.into_two())?,
        }
        Ok(())
    }

    /// Whether the session runs the N-region flow ([`RegionFlow`]).
    pub fn uses_region_flow(&self) -> bool {
        self.many.is_some()
    }

    /// The N-region flow, if the session runs it.
    pub fn region_flow(&self) -> Option<&RegionFlow> {
        self.many.as_ref()
    }

    /// The N-region flow for the calls the session does not wrap.
    pub fn region_flow_mut(&mut self) -> Option<&mut RegionFlow> {
        self.many.as_mut()
    }

    /// The options the core crate recommends for this session's number of regions:
    /// `FlowOptions::new` for two regions (the default stages), `FlowOptions::for_regions` for
    /// more (the stages of the paper with the flux recovery of the core crate, SPEC D20: on a
    /// network the default stages come to rest with more patches and a higher energy from a
    /// speckled start). Options built with `FlowOptions::new` are accepted for any number of
    /// regions and are not overridden.
    pub fn recommended_options(&self, step_size: StepSize) -> FlowOptions<'static> {
        if self.maps.n_regions > 2 {
            FlowOptions::for_regions(step_size)
        } else {
            FlowOptions::new(step_size)
        }
    }

    /// A problem with default weights for each of this session's regions
    /// (`CutProblem::for_regions`).
    pub fn new_problem(&self, energy: Energy, length_weight: f64) -> CutProblem {
        CutProblem::for_regions(energy, length_weight, self.maps.n_regions)
    }

    /// The number of regions.
    pub fn n_regions(&self) -> usize {
        self.maps.n_regions
    }

    /// The polygon mesh the session was built from (a copy).
    pub fn mesh(&self) -> &PolygonMesh {
        &self.mesh
    }

    /// The seams as given (sorted, each once).
    pub fn seams(&self) -> &[[u32; 2]] {
        &self.seams
    }

    /// The seams that were cut: the given ones that are interior polygon edges, and the
    /// extensions of one-edge seams. Sorted `[min, max]`.
    pub fn seams_cut(&self) -> &[[u32; 2]] {
        &self.maps.used_seams
    }

    /// Warnings of the validation (ignored seams, extended seams).
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// The two-region flow. **Only the flow of a two-region session**: in a session that
    /// [`uses_region_flow`](Session::uses_region_flow) this is an unused flow on the same mesh
    /// whose level-set function is not the state; use [`Session::region_flow`] there.
    pub fn flow(&self) -> &Flow {
        &self.flow
    }

    /// The two-region flow, for the calls the session does not wrap (`redistance`,
    /// `roughening_rate`, ...); see [`Session::flow`] for when it is not the state. Its
    /// numbering is that of [`Session::flow_mesh`].
    pub fn flow_mut(&mut self) -> &mut Flow {
        &mut self.flow
    }

    /// The triangle mesh the flow runs on: the triangulation cut along the seams.
    pub fn flow_mesh(&self) -> &TriMesh {
        &self.flow_mesh
    }

    /// The vertex of the polygon mesh behind every vertex of the flow mesh.
    pub fn flow_vertex_to_vertex(&self) -> &[u32] {
        &self.maps.vertex_input
    }

    /// The flow vertex of every polygon corner.
    pub fn corner_to_flow_vertex(&self) -> &[u32] {
        &self.maps.corner_vertex
    }

    /// The triangle of the supplied triangulation behind every face of the flow mesh.
    pub fn flow_face_to_triangle(&self) -> &[u32] {
        &self.maps.face_triangle
    }

    /// The state: one function per region and per polygon corner, in the caller's length unit.
    /// For two regions they are `φ` and `−φ`; region 0 is where `φ ≤ 0`, in general a corner
    /// belongs to the region whose function is the smallest there (ties: the lowest).
    pub fn state(&self) -> Vec<Vec<f64>> {
        match &self.many {
            Some(f) => self.maps.regions_from_phis(&f.level_sets()),
            None => self.maps.regions_from_phi(&self.flow.level_set()),
        }
    }

    /// The region of every polygon corner.
    pub fn corner_regions(&self) -> Vec<u32> {
        if let Some(f) = &self.many {
            let regions = f.regions();
            return self
                .maps
                .corner_vertex
                .iter()
                .map(|&v| regions[v as usize])
                .collect();
        }
        let phi = self.flow.level_set();
        self.maps
            .corner_vertex
            .iter()
            .map(|&v| u32::from(phi[v as usize] > 0.0))
            .collect()
    }

    /// Replaces the state by region functions per polygon corner (resets time and step count).
    /// With two regions the level-set function is `φ = ½(φ₀ − φ₁)`; the flow takes it as it is.
    ///
    /// # Errors
    /// [`CuttingError::WrongLength`], [`CuttingError::NonFinite`], and
    /// [`CuttingError::CornersDisagree`], which names the vertex whose corners carry different
    /// values although no seam separates them.
    pub fn set_state<R: AsRef<[f64]>>(&mut self, regions: &[R]) -> Result<(), CuttingError> {
        if let Some(f) = &mut self.many {
            f.set_level_sets(&self.maps.phis_from_regions(regions)?)?;
            return Ok(());
        }
        let phi = self.maps.phi_from_regions(regions)?;
        self.flow.set_level_set(&phi)?;
        Ok(())
    }

    /// One step of the flow. `problem.penalty` is `η` per vertex of the polygon mesh,
    /// `problem.curvature` of the kind `PerFace` has one value per triangle of the supplied
    /// triangulation, `Energy::ScaledDirichlet(PerVertex)` one value per polygon-mesh vertex,
    /// and `options.symmetry` one class per polygon-mesh vertex (a vertex that a seam splits
    /// into several copies is left out of the symmetry). The report is the flow's.
    ///
    /// # Errors
    /// [`CuttingError::WrongLength`], [`CuttingError::InvalidOption`] for values the flow would
    /// refuse (named by polygon-mesh vertex), [`CuttingError::Flow`] otherwise.
    pub fn step(
        &mut self,
        problem: &CutProblem,
        options: &FlowOptions<'_>,
    ) -> Result<StepReport, CuttingError> {
        let problem = self.flow_problem(problem)?;
        let options = self.flow_options(options)?;
        match &mut self.many {
            Some(f) => Ok(f.step(&problem, &options)?),
            None => Ok(self.flow.step(&problem, &options)?),
        }
    }

    /// Steps until the stop rule says so. Arguments as for [`Session::step`].
    ///
    /// # Errors
    /// As [`Session::step`].
    pub fn run(
        &mut self,
        problem: &CutProblem,
        options: &FlowOptions<'_>,
        stop: &StopRule,
    ) -> Result<RunReport, CuttingError> {
        let problem = self.flow_problem(problem)?;
        let options = self.flow_options(options)?;
        match &mut self.many {
            Some(f) => Ok(f.run(&problem, &options, stop)?),
            None => Ok(self.flow.run(&problem, &options, stop)?),
        }
    }

    /// The roughening rate `ω` of the current state (the scale for the length weight, see
    /// `ddg_variational_cuts::length_weight_for_ripples`). `problem` as for [`Session::step`].
    ///
    /// # Errors
    /// As [`Session::step`].
    pub fn roughening_rate(
        &mut self,
        problem: &CutProblem,
        options: &EvalOptions<'_>,
    ) -> Result<f64, CuttingError> {
        let problem = self.flow_problem(problem)?;
        match &mut self.many {
            Some(f) => Ok(f.roughening_rate(&problem, options)?),
            None => Ok(self.flow.roughening_rate(&problem, options)?),
        }
    }

    /// Evaluates the current state without changing it. The evaluation is in the numbering of
    /// the flow mesh ([`Session::flow_vertex_to_vertex`], [`Session::flow_face_to_triangle`]).
    ///
    /// # Errors
    /// As [`Session::step`].
    pub fn evaluate(
        &mut self,
        problem: &CutProblem,
        options: &EvalOptions<'_>,
    ) -> Result<Evaluation, CuttingError> {
        let problem = self.flow_problem(problem)?;
        match &mut self.many {
            Some(f) => Ok(f.evaluate(&problem, options)?),
            None => Ok(self.flow.evaluate(&problem, options)?),
        }
    }

    /// The problem of the flow mesh from the problem of the polygon mesh.
    fn flow_problem(&self, problem: &CutProblem) -> Result<CutProblem, CuttingError> {
        let m = &self.maps;
        if problem.regions.len() != m.n_regions {
            return Err(CuttingError::WrongLength {
                what: "problem.regions",
                per: "region",
                expected: m.n_regions,
                got: problem.regions.len(),
            });
        }
        let per_vertex = |what: &'static str, values: &[f64], positive: bool| {
            if values.len() != m.n_vertices {
                return Err(CuttingError::WrongLength {
                    what,
                    per: "vertex of the polygon mesh",
                    expected: m.n_vertices,
                    got: values.len(),
                });
            }
            if let Some(v) = values
                .iter()
                .position(|x| !(x.is_finite() && if positive { *x > 0.0 } else { *x >= 0.0 }))
            {
                return Err(CuttingError::InvalidOption {
                    name: what,
                    reason: format!(
                        "vertex {v} has the value {}, but it must be finite and {}",
                        values[v],
                        if positive { "positive" } else { "not negative" }
                    ),
                });
            }
            Ok(m.vertex_input
                .iter()
                .map(|&v| values[v as usize])
                .collect::<Vec<f64>>())
        };
        let mut out = problem.clone();
        if let Some(eta) = &problem.penalty {
            out.penalty = Some(per_vertex("penalty", eta, false)?);
        }
        if let Energy::ScaledDirichlet(ScaleWeight::PerVertex(a)) = &problem.energy {
            out.energy = Energy::ScaledDirichlet(ScaleWeight::PerVertex(per_vertex(
                "scaled Dirichlet weight",
                a,
                true,
            )?));
        }
        if let CurvatureSource::PerFace(k) = &problem.curvature {
            if k.len() != self.mesh.triangles.len() {
                return Err(CuttingError::WrongLength {
                    what: "curvature override",
                    per: "triangle of the triangulation",
                    expected: self.mesh.triangles.len(),
                    got: k.len(),
                });
            }
            out.curvature =
                CurvatureSource::PerFace(m.face_triangle.iter().map(|&t| k[t as usize]).collect());
        }
        Ok(out)
    }

    /// The options of the flow mesh from the options of the polygon mesh.
    fn flow_options<'a>(&self, options: &FlowOptions<'a>) -> Result<FlowOptions<'a>, CuttingError> {
        let m = &self.maps;
        let mut out = options.clone();
        if let Some(classes) = &options.symmetry {
            if classes.len() != m.n_vertices {
                return Err(CuttingError::WrongLength {
                    what: "symmetry classes",
                    per: "vertex of the polygon mesh",
                    expected: m.n_vertices,
                    got: classes.len(),
                });
            }
            // Compress the labels, then give every copy of a split vertex a class of its own.
            let mut dense: HashMap<u32, u32> = HashMap::new();
            for &c in classes {
                let next = dense.len() as u32;
                dense.entry(c).or_insert(next);
            }
            let mut fresh = dense.len() as u32;
            let mapped: Vec<u32> = m
                .vertex_input
                .iter()
                .map(|&v| {
                    if m.copies_of(v).len() == 1 {
                        dense[&classes[v as usize]]
                    } else {
                        fresh += 1;
                        fresh - 1
                    }
                })
                .collect();
            out.symmetry = Some(mapped);
        }
        Ok(out)
    }

    /// The cut of the current state, as the flow sees it.
    fn current_cut(&self, min_crossing_offset: f64) -> Result<Cut, CuttingError> {
        Ok(match &self.many {
            Some(f) => Cut::from_regions(&self.flow_mesh, &f.level_sets(), min_crossing_offset)?,
            None => Cut::from_level_set(
                &self.flow_mesh,
                &self.flow.level_set(),
                min_crossing_offset,
            )?,
        })
    }

    /// The cut of the current state as polylines for display. Only the regions of the vertices
    /// matter; `min_crossing_offset` is `δ` of the flow (use [`DEFAULT_CROSSING_OFFSET`] to
    /// see the cut the flow sees).
    ///
    /// With more than two regions (or [`FlowKind::Regions`]) the lines are the curves of the
    /// network between its junctions (`Cut::curves`): each runs from a triple point, or from
    /// the boundary, to a triple point or the boundary, or is closed; a curve that ends at a
    /// triple point has the position of that triple point as its end. Every line carries the
    /// two regions it separates.
    ///
    /// # Errors
    /// [`CuttingError::Flow`] for a `δ` outside `[0, ½)`.
    pub fn preview(&self, min_crossing_offset: f64) -> Result<Vec<PreviewPolyline>, CuttingError> {
        let cut = self.current_cut(min_crossing_offset)?;
        if self.many.is_some() {
            return Ok(cut
                .curves(&self.flow_mesh, &self.flow_positions)?
                .into_iter()
                .map(|c| PreviewPolyline {
                    positions: c.positions,
                    closed: c.ends == [CurveEnd::Closed; 2],
                    regions: c.regions,
                })
                .collect());
        }
        Ok(cut
            .polylines(&self.flow_mesh, &self.flow_positions)?
            .into_iter()
            .map(|p| PreviewPolyline {
                positions: p.positions,
                closed: p.closed,
                regions: [0, 1],
            })
            .collect())
    }

    /// The cut points as edge points and interior points.
    fn classify_points(
        &self,
        cut: &Cut,
        at: &[[f64; 3]],
    ) -> (Vec<EdgePoint>, Vec<u32>, Vec<Option<InteriorPoint>>) {
        let m = &self.maps;
        let mut edge_points = Vec::new();
        let mut edge_point_of = vec![NONE; cut.points.len()];
        let mut interior_of: Vec<Option<InteriorPoint>> = vec![None; cut.points.len()];
        for (k, point) in cut.points.iter().enumerate() {
            let e = point.edge as usize;
            let [a, b] = self.flow_mesh.edge_vertices(e);
            if m.edge_is_polygon[e] {
                let (va, vb) = (m.vertex_input[a], m.vertex_input[b]);
                let t = if va < vb { point.t } else { 1.0 - point.t };
                let edge = [va.min(vb), va.max(vb)];
                let side = m.used_seams.binary_search(&edge).ok().map(|_| {
                    let face = self.flow_mesh.face(self.flow_mesh.edge_halfedge(e));
                    m.triangle_polygon[m.face_triangle[face] as usize]
                });
                edge_point_of[k] = edge_points.len() as u32;
                edge_points.push(EdgePoint {
                    edge,
                    t,
                    position: at[k],
                    side,
                });
            } else {
                let (c, d) = m.edge_corners[e];
                let (lo, hi) = (c.min(d), c.max(d));
                let t = if m.corner_vertex[lo as usize] as usize == a {
                    point.t
                } else {
                    1.0 - point.t
                };
                interior_of[k] = Some(InteriorPoint {
                    corners: [lo, hi],
                    t,
                    position: at[k],
                });
            }
        }
        (edge_points, edge_point_of, interior_of)
    }

    /// The cut of the current state expressed on the caller's polygons, for
    /// [`insert`](super::insert).
    ///
    /// With two regions on the two-region flow, every path runs between two edge points. On the
    /// N-region flow the cut is a network: wherever three regions meet in a triangle of the
    /// supplied triangulation there is a [`JunctionPoint`] (at the centroid of the three cut
    /// points on the sides of that triangle, the triple point of the core crate), and three
    /// paths end there, one through each side of the triangle. A path runs from an edge point
    /// or a junction to an edge point or a junction, through the diagonals it crosses; it is
    /// directed so that the region with the lower id is on its left (`regions` is
    /// `[lower, higher]`).
    ///
    /// # Errors
    /// [`CuttingError::Flow`] for a `δ` outside `[0, ½)`, [`CuttingError::Internal`] if the
    /// segments of the cut do not form paths (cannot happen for a valid flow mesh).
    pub fn describe_cut(&self, min_crossing_offset: f64) -> Result<CutOnPolygons, CuttingError> {
        let cut = self.current_cut(min_crossing_offset)?;
        if self.many.is_some() {
            self.describe_network(&cut)
        } else {
            self.describe_two(&cut)
        }
    }

    /// [`Session::describe_cut`] for the two-region flow.
    fn describe_two(&self, cut: &Cut) -> Result<CutOnPolygons, CuttingError> {
        let m = &self.maps;
        let at = cut.point_positions(&self.flow_mesh, &self.flow_positions)?;
        let (edge_points, edge_point_of, interior_of) = self.classify_points(cut, &at);
        let mut description = CutOnPolygons {
            n_clamped_crossings: cut.n_clamped,
            edge_points,
            ..CutOnPolygons::default()
        };
        // Segments by polygon; a path starts where a segment starts on a polygon edge.
        let polygon_of = |seg: usize| {
            m.triangle_polygon[m.face_triangle[cut.segments[seg].face as usize] as usize]
        };
        let mut order: Vec<usize> = (0..cut.segments.len()).collect();
        order.sort_by_key(|&s| (polygon_of(s), s));
        let mut by_first: HashMap<(u32, u32), usize> = HashMap::new();
        for &s in &order {
            by_first.insert((polygon_of(s), cut.segments[s].points[0]), s);
        }
        let mut used = vec![false; cut.segments.len()];
        for &s in &order {
            let polygon = polygon_of(s);
            let first = cut.segments[s].points[0];
            if edge_point_of[first as usize] == NONE {
                continue;
            }
            let mut interior = Vec::new();
            let mut current = s;
            let end = loop {
                used[current] = true;
                let last = cut.segments[current].points[1];
                if edge_point_of[last as usize] != NONE {
                    break PathEnd::EdgePoint(edge_point_of[last as usize]);
                }
                interior.push(interior_of[last as usize].clone().ok_or_else(|| {
                    internal("a cut point is neither on a polygon edge nor on a diagonal")
                })?);
                current = *by_first.get(&(polygon, last)).ok_or_else(|| {
                    internal(format!(
                        "the cut leaves polygon {polygon} through a diagonal point that has no \
                         continuation"
                    ))
                })?;
            };
            description.paths.push(PolygonPath {
                polygon,
                start: PathEnd::EdgePoint(edge_point_of[first as usize]),
                interior,
                end,
                regions: [0, 1],
            });
        }
        if let Some(s) = used.iter().position(|u| !u) {
            return Err(internal(format!(
                "a segment of the cut in polygon {} is on no path from a polygon edge",
                polygon_of(s)
            )));
        }
        description.paths.sort_by_key(|p| {
            let PathEnd::EdgePoint(i) = p.start else {
                unreachable!("two regions")
            };
            (p.polygon, i)
        });
        Ok(description)
    }

    /// [`Session::describe_cut`] for a network (the N-region flow).
    fn describe_network(&self, cut: &Cut) -> Result<CutOnPolygons, CuttingError> {
        let m = &self.maps;
        let n_points = cut.points.len();
        let at = cut.point_positions(&self.flow_mesh, &self.flow_positions)?;
        let (edge_points, edge_point_of, interior_of) = self.classify_points(cut, &at);
        let triple_at = cut.triple_positions(&self.flow_mesh, &self.flow_positions)?;
        let mut description = CutOnPolygons {
            n_clamped_crossings: cut.n_clamped,
            edge_points,
            ..CutOnPolygons::default()
        };
        for (mi, triple) in cut.triple_points.iter().enumerate() {
            let triangle = m.face_triangle[triple.face as usize];
            description.junctions.push(JunctionPoint {
                position: triple_at[mi],
                polygon: m.triangle_polygon[triangle as usize],
                triangle,
                corners: self.mesh.triangles[triangle as usize],
                barycentric: triple.barycentric,
            });
        }
        let polygon_of = |seg: usize| {
            m.triangle_polygon[m.face_triangle[cut.segments[seg].face as usize] as usize]
        };
        // The segments at every node, per polygon.
        let mut incident: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (s, segment) in cut.segments.iter().enumerate() {
            for node in segment.points {
                incident.entry((polygon_of(s), node)).or_default().push(s);
            }
        }
        // A path ends at an edge point or at a junction; diagonal points have two segments.
        let is_terminal = |node: u32| node as usize >= n_points || edge_point_of[node as usize] != NONE;
        let mut order: Vec<usize> = (0..cut.segments.len()).collect();
        order.sort_by_key(|&s| (polygon_of(s), s));
        let mut used = vec![false; cut.segments.len()];
        for &s in &order {
            if used[s] {
                continue;
            }
            let polygon = polygon_of(s);
            let [x, y] = cut.segments[s].points;
            let (start, forward) = if is_terminal(x) {
                (x, true)
            } else if is_terminal(y) {
                (y, false)
            } else {
                continue;
            };
            let mut nodes = vec![start];
            let (mut current, mut from) = (s, start);
            loop {
                used[current] = true;
                let [a, b] = cut.segments[current].points;
                let to = if a == from { b } else { a };
                nodes.push(to);
                if is_terminal(to) {
                    break;
                }
                let list = &incident[&(polygon, to)];
                if list.len() != 2 {
                    return Err(internal(format!(
                        "a diagonal point of the cut in polygon {polygon} has {} segments",
                        list.len()
                    )));
                }
                current = if list[0] == current { list[1] } else { list[0] };
                from = to;
            }
            // The regions on the left and on the right of the first segment, walked from
            // `start`.
            let [p0, p1] = cut.segments[s].points;
            let (left, right) = if p1 as usize >= n_points {
                let triple = &cut.triple_points[p1 as usize - n_points];
                let k = triple.points.iter().position(|&q| q == p0).ok_or_else(|| {
                    internal("a segment ends at a triple point that does not list its first end")
                })?;
                let face = self.flow_mesh.faces()[triple.face as usize];
                (
                    cut.vertex_region[face[k] as usize],
                    cut.vertex_region[face[(k + 1) % 3] as usize],
                )
            } else {
                let [lo, hi] = cut.segment_regions(&self.flow_mesh, s);
                (lo, hi)
            };
            let (mut left, mut right) = if forward { (left, right) } else { (right, left) };
            if left > right {
                nodes.reverse();
                std::mem::swap(&mut left, &mut right);
            }
            let end_of = |node: u32| {
                if node as usize >= n_points {
                    PathEnd::Junction(node - n_points as u32)
                } else {
                    PathEnd::EdgePoint(edge_point_of[node as usize])
                }
            };
            let interior = nodes[1..nodes.len() - 1]
                .iter()
                .map(|&node| {
                    interior_of[node as usize].clone().ok_or_else(|| {
                        internal("a cut point is neither on a polygon edge nor on a diagonal")
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            description.paths.push(PolygonPath {
                polygon,
                start: end_of(nodes[0]),
                interior,
                end: end_of(nodes[nodes.len() - 1]),
                regions: [left, right],
            });
        }
        if let Some(s) = used.iter().position(|u| !u) {
            return Err(internal(format!(
                "a segment of the cut in polygon {} is on no path from a polygon edge or a junction",
                polygon_of(s)
            )));
        }
        description.paths.sort_by_key(|p| {
            let key = |end: PathEnd| match end {
                PathEnd::EdgePoint(i) => (0, i),
                PathEnd::Junction(i) => (1, i),
            };
            (p.polygon, key(p.start), key(p.end))
        });
        Ok(description)
    }
}

/// A start: one level-set function (two-region flow) or one per region.
enum Level {
    Phi(Vec<f64>),
    Phis(Vec<Vec<f64>>),
}

impl Level {
    /// For the two-region flow: `φ`, or `½(φ₀ − φ₁)` of two functions.
    fn into_two(self) -> Vec<f64> {
        match self {
            Level::Phi(phi) => phi,
            Level::Phis(p) => p[0]
                .iter()
                .zip(&p[1])
                .map(|(a, b)| 0.5 * (a - b))
                .collect(),
        }
    }

    /// For the N-region flow: the functions, or `[φ, −φ]`.
    fn into_many(self) -> Vec<Vec<f64>> {
        match self {
            Level::Phi(phi) => {
                let negated = phi.iter().map(|x| -x).collect();
                vec![phi, negated]
            }
            Level::Phis(p) => p,
        }
    }
}

/// The state of the flow mesh for a start. `many` says whether the N-region flow will take it
/// (labels then give one function per region).
fn start_state(
    maps: &Maps,
    mesh: &TriMesh,
    geometry: &Geometry,
    positions: &[[f64; 3]],
    start: &Start,
    many: bool,
) -> Result<Level, CuttingError> {
    match start {
        Start::Corners(regions) => {
            if many {
                Ok(Level::Phis(maps.phis_from_regions(regions)?))
            } else {
                Ok(Level::Phi(maps.phi_from_regions(regions)?))
            }
        }
        Start::PolygonLabels(labels) => {
            let expected = maps.n_polygons;
            if labels.len() != expected {
                return Err(CuttingError::WrongLength {
                    what: "polygon labels",
                    per: "polygon",
                    expected,
                    got: labels.len(),
                });
            }
            if let Some(p) = labels.iter().position(|&l| l as usize >= maps.n_regions) {
                return Err(CuttingError::LabelOutOfRange {
                    polygon: p as u32,
                    label: labels[p],
                    n_regions: maps.n_regions,
                });
            }
            let face_labels: Vec<u32> = maps
                .face_triangle
                .iter()
                .map(|&t| labels[maps.triangle_polygon[t as usize] as usize])
                .collect();
            if many {
                Ok(Level::Phis(init::regions_from_face_labels(
                    mesh,
                    geometry,
                    &face_labels,
                    maps.n_regions,
                )?))
            } else {
                Ok(Level::Phi(init::from_face_labels(
                    mesh,
                    geometry,
                    &face_labels,
                )?))
            }
        }
        Start::NormalClusters { max_iterations } => {
            let labels = init::normal_clusters(mesh, positions, maps.n_regions, *max_iterations)?;
            let phis = init::regions_from_vertex_labels(mesh, geometry, &labels, maps.n_regions)?;
            Ok(Level::Phis(phis))
        }
        Start::GeodesicCircle { vertex, radius } => {
            if maps.n_regions != 2 {
                return Err(CuttingError::InvalidOption {
                    name: "start",
                    reason: format!(
                        "a geodesic circle starts two regions, but there are {}",
                        maps.n_regions
                    ),
                });
            }
            if *vertex as usize >= maps.n_vertices {
                return Err(CuttingError::VertexOutOfRange {
                    what: "center",
                    vertex: *vertex,
                    n_vertices: maps.n_vertices,
                });
            }
            if !(radius.is_finite() && *radius > 0.0) {
                return Err(CuttingError::InvalidOption {
                    name: "radius",
                    reason: format!("must be positive and finite, got {radius}"),
                });
            }
            let seeds: Vec<(u32, f64)> =
                maps.copies_of(*vertex).iter().map(|&v| (v, 0.0)).collect();
            if seeds.is_empty() {
                return Err(CuttingError::InvalidOption {
                    name: "center",
                    reason: format!("vertex {vertex} is not used by any polygon"),
                });
            }
            let distance = upwind_distance(mesh, geometry, &seeds, None)
                .map_err(|e| internal(format!("the distance from vertex {vertex} failed: {e}")))?;
            let mut phi: Vec<f64> = distance.values.iter().map(|d| d - radius).collect();
            // A vertex no front reached (another component) is outside.
            let far = phi
                .iter()
                .copied()
                .filter(|x| x.is_finite())
                .fold(*radius, f64::max);
            for x in &mut phi {
                if !x.is_finite() {
                    *x = far;
                }
            }
            Ok(Level::Phi(phi))
        }
    }
}
