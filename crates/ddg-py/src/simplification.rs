//! `simplify`, `build_hierarchy` and `qem_simplify`: surface simplification with intrinsic error
//! metrics (Liu, Gillespie, Chislett, Sharp, Jacobson, Crane, *Surface Simplification using
//! Intrinsic Error Metrics*, ACM TOG 42(4), 2023), from the crate `ddg-intrinsic-simplification`
//! (contract `SPEC.md` section 11).
//!
//! As in [`crate::intrinsic_apps`] the Python glue decodes buffers into plain vectors; the
//! functions here validate them (messages name the parameter), run the crate and return plain
//! vectors and the crate's own report types. No Python is involved, so `cargo test -p ddg-py`
//! covers all of it.
//!
//! * [`Simplification::new`]: the options of Alg. 1 (target, curvature/area weights, memory,
//!   fixed vertices, extra mass channel, extra tracked points, pre- and post-refinement,
//!   anisotropy) and the run; the result keeps the coarse mesh, its compacted output (faces,
//!   lengths, point locations, colouring), `P` and `P^vec`, and builds the applications of
//!   SPEC 10.1 (Poisson, geodesic distance, Voronoi) on first use.
//! * [`HierarchyHandle`] and [`MultigridHandle`]: the intrinsic hierarchy and the multigrid
//!   solver on it (SPEC 9).
//! * [`qem`]: the quadric error metric baseline (SPEC 10.4).

use ddg_core::sparse::ThreadPool;
use ddg_core::{Geometry, SparseMatrix, TriMesh};
use ddg_intrinsic_apps::exact::ExactOptions;
use ddg_intrinsic_apps::representation::Representation;
use ddg_intrinsic_simplification::anisotropy::{
    AnisotropyOptions, AnisotropyReport, PrincipalDirection, anisotropic_lengths,
    principal_direction_field,
};
use ddg_intrinsic_simplification::applications::{CoarseApplications, PoissonReport};
use ddg_intrinsic_simplification::baseline::{QemOptions, QemReport, qem_simplify};
use ddg_intrinsic_simplification::curvature::{CurvatureOptions, direction_angles};
use ddg_intrinsic_simplification::hierarchy::{
    Hierarchy, HierarchyOptions, HierarchyReport, build_hierarchy,
};
use ddg_intrinsic_simplification::multigrid::{
    Multigrid, MultigridInfo, MultigridOptions, SolveReport,
};
use ddg_intrinsic_simplification::refine::{RefinementOptions, RefinementReport};
use ddg_intrinsic_simplification::{
    ChannelSource, ChannelWeights, CoarseMesh, CoarseMeshOptions, CoarseOutput, CoarsePoint,
    ProjectiveSign,
    FineFrames, InputPoint, MetricOptions, SimplifyOptions, SimplifyReport, Target,
    VectorProlongationReport, simplify,
};

use crate::error::{InputError, UnwrapError};
use crate::intrinsic_apps::{Mesh, MeshInput, bad_option, build};

fn message(text: String) -> UnwrapError {
    InputError::Message(text).into()
}

fn iem(e: ddg_intrinsic_simplification::IemError) -> UnwrapError {
    use ddg_intrinsic_simplification::IemError as E;
    match e {
        E::Mesh(_) | E::Geometry(_) => UnwrapError::Mesh(e.to_string()),
        E::InvalidPoint { .. } | E::InvalidOptions { .. } | E::IndexOutOfRange { .. } => {
            message(e.to_string())
        }
        other => UnwrapError::Solve(other.to_string()),
    }
}

fn refine_error(e: ddg_intrinsic_simplification::refine::RefineError) -> UnwrapError {
    use ddg_intrinsic_simplification::refine::RefineError as E;
    match e {
        E::Coarse(inner) => iem(inner),
        other => UnwrapError::Solve(other.to_string()),
    }
}

fn app_error(e: ddg_intrinsic_simplification::applications::AppError) -> UnwrapError {
    use ddg_intrinsic_simplification::applications::AppError as E;
    match e {
        E::Invalid(text) => message(text),
        other => UnwrapError::Solve(other.to_string()),
    }
}

fn projective_sign(name: &str) -> Result<ProjectiveSign, UnwrapError> {
    match name {
        "negative" => Ok(ProjectiveSign::Negative),
        "positive" => Ok(ProjectiveSign::Positive),
        other => Err(bad_option(
            "projective_sign",
            format!("{other:?}"),
            "must be \"negative\" or \"positive\"",
        )),
    }
}

fn representation(name: &str) -> Result<Representation, UnwrapError> {
    match name {
        "integer" => Ok(Representation::Integer),
        "signposts" => Ok(Representation::Signposts),
        other => Err(bad_option(
            "representation",
            format!("{other:?}"),
            "must be \"integer\" or \"signposts\"",
        )),
    }
}

/// The principal direction (or the directions given) that anisotropic coarsening elongates
/// across (SPEC 8).
#[derive(Debug, Clone)]
pub enum AnisotropyDirection {
    /// The larger principal curvature's direction.
    Max,
    /// The smaller principal curvature's direction.
    Min,
    /// A 3D direction (a line) per vertex.
    Given(Vec<f64>),
}

/// The anisotropy of `simplify` (paper Sec. 8.6.2).
#[derive(Debug, Clone)]
pub struct AnisotropyArguments {
    /// The strength in `[0, 1]`.
    pub tau: f64,
    /// The direction field.
    pub direction: AnisotropyDirection,
    /// Lower `tau` to the largest valid value instead of failing.
    pub clamp: bool,
}

/// The arguments of [`Simplification::new`] (the Python keywords).
#[derive(Debug, Clone)]
pub struct SimplifyArguments {
    /// Stop at this many vertices (of the mesh being simplified: the refined one with
    /// `refine_before`).
    pub target: Option<usize>,
    /// Stop at `round(ratio * input vertices)`.
    pub ratio: Option<f64>,
    /// The curvature/area weight `w` in `[0, 1]`: `mixture(w)` (0: the paper's curvature
    /// metric, 1: the area).
    pub weights: f64,
    /// Memory-based metric (default) or the memoryless one.
    pub memory: bool,
    /// Fine vertices that stay.
    pub fixed: Vec<u32>,
    /// An extra per-vertex mass channel (one real per input vertex).
    pub masses: Option<Vec<f64>>,
    /// The weight of the `masses` channel.
    pub mass_weight: f64,
    /// Extra tracked points: faces of the input mesh ...
    pub point_faces: Vec<u32>,
    /// ... and 3 barycentric coordinates each.
    pub point_bary: Vec<f64>,
    /// Refine the input before simplifying (Fig. 12).
    pub refine_before: bool,
    /// Refine the coarse mesh afterwards (Fig. 13).
    pub refine_after: bool,
    /// Repair a non-simplicial coarse mesh by tracked flips before post-refinement.
    pub repair_simplicial: bool,
    /// Anisotropic coarsening.
    pub anisotropy: Option<AnisotropyArguments>,
    /// `"integer"` or `"signposts"`: the refinement's representation.
    pub representation: String,
    /// `"negative"` (default, `e^{-u}`) or `"positive"` (`e^{+u}`): the sign of Eq. 10, the map of the tracked points
    /// through a flattening (SPEC S25).
    pub projective_sign: String,
}

impl Default for SimplifyArguments {
    fn default() -> Self {
        SimplifyArguments {
            target: None,
            ratio: None,
            weights: 0.0,
            memory: true,
            fixed: Vec::new(),
            masses: None,
            mass_weight: 1.0,
            point_faces: Vec::new(),
            point_bary: Vec::new(),
            refine_before: false,
            refine_after: false,
            repair_simplicial: true,
            anisotropy: None,
            representation: "integer".to_owned(),
            projective_sign: "negative".to_owned(),
        }
    }
}

/// A sparse matrix as triplets.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Coo {
    /// Rows.
    pub nrows: usize,
    /// Columns.
    pub ncols: usize,
    /// Row index per entry.
    pub rows: Vec<u32>,
    /// Column index per entry.
    pub cols: Vec<u32>,
    /// Value per entry.
    pub values: Vec<f64>,
}

impl Coo {
    fn of(m: &SparseMatrix) -> Coo {
        let mut c = Coo {
            nrows: m.nrows(),
            ncols: m.ncols(),
            ..Coo::default()
        };
        for (r, col, v) in m.triplets() {
            c.rows.push(r);
            c.cols.push(col);
            c.values.push(v);
        }
        c
    }
}

/// `P^vec` as complex triplets.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ComplexCoo {
    /// Rows.
    pub nrows: usize,
    /// Columns.
    pub ncols: usize,
    /// Row per entry.
    pub rows: Vec<u32>,
    /// Column per entry.
    pub cols: Vec<u32>,
    /// Real parts.
    pub re: Vec<f64>,
    /// Imaginary parts.
    pub im: Vec<f64>,
    /// How the rows were found.
    pub report: VectorProlongationReport,
}

/// Where every tracked point is on the coarse mesh.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PointLocations {
    /// The coarse face of each point (a point that is a coarse vertex: its incident face of the
    /// lowest id), `-1` for none.
    pub face: Vec<i64>,
    /// Barycentric coordinates in that face's corners (3 per point; a vertex: 1 at its corner).
    pub bary: Vec<f64>,
    /// The compacted coarse vertex a point is, `-1` for a point inside a face.
    pub vertex: Vec<i64>,
}

/// What [`Simplification::poisson`] returns.
#[derive(Debug, Clone)]
pub struct PoissonOutput {
    /// The solution on the coarse vertices.
    pub coarse: Vec<f64>,
    /// The prolongated solution on the fine vertices.
    pub fine: Vec<f64>,
    /// The solver's report.
    pub report: PoissonReport,
}

/// What [`Simplification::voronoi`] returns.
#[derive(Debug, Clone)]
pub struct VoronoiOutput {
    /// Nearest source per fine vertex.
    pub fine_labels: Vec<u32>,
    /// Nearest source per coarse vertex.
    pub coarse_labels: Vec<u32>,
    /// Prolongated distance to it.
    pub fine_distance: Vec<f64>,
}

/// The result of a simplification: the coarse mesh in several forms, and what the applications
/// need.
pub struct Simplification {
    /// The compacted coarse mesh with the point locations.
    pub output: CoarseOutput,
    /// The report of Alg. 1.
    pub report: SimplifyReport,
    /// The pre-refinement's report.
    pub pre: Option<RefinementReport>,
    /// The post-refinement's report.
    pub post: Option<RefinementReport>,
    /// The anisotropy's report.
    pub anisotropy: Option<AnisotropyReport>,
    /// Vertices of the input mesh.
    pub n_input: usize,
    /// Seconds of the whole call (construction, refinements, simplification, compaction).
    pub seconds: f64,
    coarse: CoarseMesh,
    mesh: TriMesh,
    geometry: Geometry,
    plain_frames: bool,
    apps: Option<CoarseApplications>,
}

fn check_buffer_len(name: &'static str, got: usize, expected: usize) -> Result<(), UnwrapError> {
    if got == expected {
        Ok(())
    } else {
        Err(InputError::LengthMismatch {
            name,
            got,
            expected,
            reason: "the number of rows follows from the mesh or from the other array",
        }
        .into())
    }
}

impl Simplification {
    /// Builds the coarse mesh of the input (with the refinements and the anisotropy asked for),
    /// simplifies it and compacts the result.
    ///
    /// # Errors
    /// [`UnwrapError::Input`] for malformed arguments (naming the parameter), [`UnwrapError::Mesh`]
    /// for a mesh that is not a manifold triangle mesh, [`UnwrapError::Solve`] for a failure of a
    /// refinement or of the simplification (with the crate's message).
    pub fn new(input: &MeshInput, args: &SimplifyArguments) -> Result<Simplification, UnwrapError> {
        let t0 = std::time::Instant::now();
        let Mesh {
            mesh,
            geometry,
            positions,
        } = build(input)?;
        let n = positions.len();
        let rep = representation(&args.representation)?;
        if !(0.0..=1.0).contains(&args.weights) || args.weights.is_nan() {
            return Err(bad_option(
                "weights",
                args.weights.to_string(),
                "must be in [0, 1] (0: curvature, 1: area)",
            ));
        }
        if !(args.mass_weight.is_finite() && args.mass_weight >= 0.0) {
            return Err(bad_option(
                "mass_weight",
                args.mass_weight.to_string(),
                "must be finite and >= 0",
            ));
        }
        let target_count = match (args.target, args.ratio) {
            (Some(_), Some(_)) => {
                return Err(message(
                    "target and ratio: give one of them, not both".to_owned(),
                ));
            }
            (None, None) => {
                return Err(message(
                    "target or ratio: one of them is required".to_owned(),
                ));
            }
            (Some(t), None) => t,
            (None, Some(r)) => {
                if !(r.is_finite() && r > 0.0 && r <= 1.0) {
                    return Err(bad_option("ratio", r.to_string(), "must be in (0, 1]"));
                }
                ((r * n as f64).round() as usize).max(1)
            }
        };
        // Extra tracked points.
        check_buffer_len(
            "points (bary)",
            args.point_bary.len(),
            3 * args.point_faces.len(),
        )?;
        let points: Vec<InputPoint> = args
            .point_faces
            .iter()
            .enumerate()
            .map(|(k, &face)| InputPoint {
                face,
                bary: [
                    args.point_bary[3 * k],
                    args.point_bary[3 * k + 1],
                    args.point_bary[3 * k + 2],
                ],
            })
            .collect();
        // The metric: scaled by the anisotropy, else the input's.
        let (geometry_used, anisotropy_report) = match &args.anisotropy {
            None => (geometry.clone(), None),
            Some(a) => {
                let angles = match &a.direction {
                    AnisotropyDirection::Max | AnisotropyDirection::Min => {
                        let which = if matches!(a.direction, AnisotropyDirection::Max) {
                            PrincipalDirection::Max
                        } else {
                            PrincipalDirection::Min
                        };
                        principal_direction_field(
                            &positions,
                            mesh.faces(),
                            which,
                            &CurvatureOptions::default(),
                        )
                        .map_err(|e| UnwrapError::Solve(e.to_string()))?
                        .angles
                    }
                    AnisotropyDirection::Given(d) => {
                        check_buffer_len("anisotropy_direction", d.len(), 3 * n)?;
                        let dirs: Vec<[f64; 3]> =
                            d.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
                        direction_angles(&mesh, &positions, &dirs)
                            .map_err(|e| message(format!("anisotropy_direction: {e}")))?
                    }
                };
                let options = AnisotropyOptions {
                    tau: a.tau,
                    clamp_tau: a.clamp,
                    ..AnisotropyOptions::default()
                };
                let (lengths, report) =
                    anisotropic_lengths(&mesh, geometry.edge_lengths(), &angles, &options)
                        .map_err(|e| message(format!("anisotropy: {e}")))?;
                let g = Geometry::from_edge_lengths(&mesh, lengths)
                    .map_err(|e| UnwrapError::Mesh(format!("anisotropy: {e}")))?;
                (g, Some(report))
            }
        };
        let mesh_options = CoarseMeshOptions {
            projective_sign: projective_sign(&args.projective_sign)?,
            ..CoarseMeshOptions::default()
        };
        let refinement = RefinementOptions {
            representation: rep,
            make_simplicial_flips: args.repair_simplicial.then_some(1000),
            ..RefinementOptions::default()
        };
        let (mut coarse, pre) = if args.refine_before {
            let (c, r) =
                CoarseMesh::from_refined(&mesh, &geometry_used, &points, mesh_options, &refinement)
                    .map_err(refine_error)?;
            (c, Some(r))
        } else {
            (
                CoarseMesh::from_mesh(&mesh, &geometry_used, &points, mesh_options).map_err(iem)?,
                None,
            )
        };
        let mut channels = ChannelWeights::mixture(args.weights);
        if let Some(m) = &args.masses {
            check_buffer_len("masses", m.len(), n)?;
            channels = channels.with(ChannelSource::PerVertex(m.clone()), args.mass_weight);
        }
        let options = SimplifyOptions {
            target: Target::Vertices(target_count),
            channels,
            metric: MetricOptions {
                memory: args.memory,
                ..MetricOptions::default()
            },
            fixed: args.fixed.clone(),
            ..SimplifyOptions::default()
        };
        let report = simplify(&mut coarse, &options).map_err(iem)?;
        let (coarse, post) = if args.refine_after {
            let (c, r) = coarse.post_refine(&refinement).map_err(refine_error)?;
            (c, Some(r))
        } else {
            (coarse, None)
        };
        let output = coarse.compact().map_err(iem)?;
        let plain_frames = !args.refine_before && !args.refine_after && args.anisotropy.is_none();
        Ok(Simplification {
            output,
            report,
            pre,
            post,
            anisotropy: anisotropy_report,
            n_input: n,
            seconds: t0.elapsed().as_secs_f64(),
            coarse,
            mesh,
            geometry,
            plain_frames,
            apps: None,
        })
    }

    /// The faces as fine vertex ids (3 per face).
    pub fn faces_fine(&self) -> Vec<u32> {
        self.output.faces_fine().into_iter().flatten().collect()
    }

    /// The faces as compacted vertex ids (3 per face).
    pub fn faces_compact(&self) -> Vec<u32> {
        self.output.faces.iter().flatten().copied().collect()
    }

    /// The side lengths of every face (3 per face: slot `k` to slot `k + 1`).
    pub fn lengths(&self) -> Vec<f64> {
        self.output.face_lengths.iter().flatten().copied().collect()
    }

    /// Where each tracked point is (the fine vertices first, then the extra points).
    pub fn point_locations(&self) -> PointLocations {
        let mut out = PointLocations::default();
        for (p, loc) in self.output.points.iter().enumerate() {
            match *loc {
                CoarsePoint::Face { face, bary } => {
                    out.face.push(i64::from(face));
                    out.bary.extend_from_slice(&bary);
                    out.vertex.push(-1);
                }
                CoarsePoint::Vertex(v) => {
                    let face = self.output.face_of_point(p);
                    out.face.push(face.map_or(-1, i64::from));
                    let mut b = [0.0; 3];
                    if let Some(f) = face
                        && let Some(k) = self.output.faces[f as usize].iter().position(|&c| c == v)
                    {
                        b[k] = 1.0;
                    }
                    out.bary.extend_from_slice(&b);
                    out.vertex.push(i64::from(v));
                }
            }
        }
        out
    }

    /// The colour of every coarse face and the number of colours (SPEC 7.4).
    pub fn face_colours(&self) -> (Vec<u32>, u32) {
        let c = self.output.colour_faces();
        (c.colours, c.n_colours)
    }

    /// The colour of the coarse face of every tracked point, `-1` for none.
    pub fn point_colours(&self) -> Vec<i64> {
        self.output
            .point_colours()
            .into_iter()
            .map(|c| c.map_or(-1, i64::from))
            .collect()
    }

    /// `P`: fine vertices by coarse vertices (SPEC 7.2).
    pub fn prolongation(&self) -> Result<Coo, UnwrapError> {
        let p = self
            .output
            .scalar_prolongation()
            .map_err(|e| UnwrapError::Solve(e.to_string()))?;
        Ok(Coo::of(&p))
    }

    /// The same for every tracked point (the extra points, `points=`, included).
    pub fn point_prolongation(&self) -> Result<Coo, UnwrapError> {
        let p = self
            .output
            .point_prolongation()
            .map_err(|e| UnwrapError::Solve(e.to_string()))?;
        Ok(Coo::of(&p))
    }

    /// `P^vec` for `n_sym`-direction fields (SPEC 7.3), in the normalized frames of the input.
    ///
    /// # Errors
    /// [`UnwrapError::Input`] when the run used pre- or post-refinement or anisotropy (the frames
    /// of the coarse mesh are then not those of the input mesh).
    pub fn vector_prolongation(&self, n_sym: u32) -> Result<ComplexCoo, UnwrapError> {
        if !self.plain_frames {
            return Err(message(
                "vector prolongation: only without refine_before, refine_after and anisotropy \
                 (the coarse frames are then not those of the input mesh)"
                    .to_owned(),
            ));
        }
        if n_sym == 0 {
            return Err(bad_option("n", "0".to_owned(), "must be >= 1"));
        }
        let frames = FineFrames::from_mesh(&self.mesh, &self.geometry)
            .map_err(|e| UnwrapError::Solve(e.to_string()))?;
        let (m, report) = self
            .output
            .vector_prolongation(&frames, n_sym)
            .map_err(|e| UnwrapError::Solve(e.to_string()))?;
        let mut out = ComplexCoo {
            nrows: m.nrows,
            ncols: m.ncols,
            report,
            ..ComplexCoo::default()
        };
        for i in 0..m.nrows {
            let (cols, vals) = m.row(i);
            for (&c, v) in cols.iter().zip(vals) {
                out.rows.push(i as u32);
                out.cols.push(c);
                out.re.push(v.re);
                out.im.push(v.im);
            }
        }
        Ok(out)
    }

    fn apps(&mut self) -> Result<&CoarseApplications, UnwrapError> {
        if self.apps.is_none() {
            self.apps = Some(CoarseApplications::new(&self.coarse).map_err(app_error)?);
        }
        Ok(self.apps.as_ref().expect("just built"))
    }

    /// The coarse vertex of each fixed fine vertex.
    fn coarse_ids(&self, name: &'static str, fine: &[u32]) -> Result<Vec<u32>, UnwrapError> {
        fine.iter()
            .map(|&v| match self.output.coarse_of_fine.get(v as usize) {
                Some(Some(c)) => Ok(*c),
                Some(None) => Err(message(format!(
                    "{name}: vertex {v} was removed; keep it with simplify(fixed=...)"
                ))),
                None => Err(message(format!(
                    "{name}: vertex {v} is out of range ({} input vertices)",
                    self.n_input
                ))),
            })
            .collect()
    }

    /// Poisson's equation `L x = M f` on the coarse mesh, prolongated (Fig. 21). `f` has one
    /// value per input vertex (sampled at the coarse vertices); `fixed` are input vertices that
    /// are still coarse vertices, with `values`; without them the mean-zero problem is solved.
    pub fn poisson(
        &mut self,
        f: &[f64],
        fixed: &[u32],
        values: &[f64],
    ) -> Result<PoissonOutput, UnwrapError> {
        check_buffer_len("f", f.len(), self.n_input)?;
        check_buffer_len("values", values.len(), fixed.len())?;
        let coarse_fixed = self.coarse_ids("fixed", fixed)?;
        let apps = self.apps()?;
        let fbar = apps.sample(f).map_err(app_error)?;
        let (coarse, fine, report) = apps
            .poisson(&fbar, &coarse_fixed, values)
            .map_err(app_error)?;
        Ok(PoissonOutput {
            coarse,
            fine,
            report,
        })
    }

    /// Geodesic distance from tracked point `source` on the coarse mesh and prolongated:
    /// `"exact"` (Δ-complex windows) or `"heat"` (Crane et al. 2013).
    pub fn geodesic(
        &mut self,
        source: usize,
        method: &str,
        time_factor: f64,
    ) -> Result<(Vec<f64>, Vec<f64>), UnwrapError> {
        let apps = self.apps()?;
        match method {
            "exact" => apps
                .geodesic_exact(source, &ExactOptions::default())
                .map_err(app_error),
            "heat" => apps.geodesic_heat(source, time_factor).map_err(app_error),
            other => Err(bad_option(
                "method",
                format!("{other:?}"),
                "must be \"exact\" or \"heat\"",
            )),
        }
    }

    /// Geodesic Voronoi labels from tracked points `sources` (Fig. 25).
    pub fn voronoi(
        &mut self,
        sources: &[usize],
        threads: usize,
    ) -> Result<VoronoiOutput, UnwrapError> {
        let pool = ThreadPool::new(threads.clamp(1, 16))
            .map_err(|e| UnwrapError::Solve(format!("thread pool: {e}")))?;
        let apps = self.apps()?;
        let v = apps
            .voronoi(sources, &pool, &ExactOptions::default())
            .map_err(app_error)?;
        Ok(VoronoiOutput {
            fine_labels: v.fine_labels,
            coarse_labels: v.coarse_labels,
            fine_distance: v.fine_distance,
        })
    }
}

// ---------------------------------------------------------------------------------------------
// Hierarchy and multigrid
// ---------------------------------------------------------------------------------------------

/// The arguments of [`HierarchyHandle::new`].
#[derive(Debug, Clone)]
pub struct HierarchyArguments {
    /// Vertex count of a level over the one above, in `(0, 1)`.
    pub ratio: f64,
    /// Levels are added while the current one has more vertices than this.
    pub n_min: usize,
    /// One run per level with fresh masses.
    pub fresh_levels: bool,
    /// Curvature/area weight as in [`SimplifyArguments`].
    pub weights: f64,
    /// Memory-based metric.
    pub memory: bool,
    /// Vertices kept in every level.
    pub fixed: Vec<u32>,
}

impl Default for HierarchyArguments {
    fn default() -> Self {
        HierarchyArguments {
            ratio: 0.25,
            n_min: 500,
            fresh_levels: false,
            weights: 0.0,
            memory: true,
            fixed: Vec::new(),
        }
    }
}

/// An intrinsic hierarchy of a mesh (SPEC 9.1).
pub struct HierarchyHandle {
    hierarchy: Hierarchy,
}

/// The arguments of [`HierarchyHandle::multigrid`].
#[derive(Debug, Clone)]
pub struct MultigridArguments {
    /// The matrix as triplets (square, symmetric positive definite after elimination of
    /// `known`); `None`: `L + shift M` of the finest level.
    pub matrix: Option<Coo>,
    /// The shift of the default matrix.
    pub shift: f64,
    /// Vertices with known values.
    pub known: Vec<u32>,
    /// Smoothing sweeps before and after the coarse correction.
    pub pre_sweeps: usize,
    /// Backward sweeps after.
    pub post_sweeps: usize,
    /// Relative residual tolerance.
    pub tolerance: f64,
    /// Cap on V-cycles.
    pub max_cycles: usize,
}

impl Default for MultigridArguments {
    fn default() -> Self {
        let d = MultigridOptions::default();
        MultigridArguments {
            matrix: None,
            shift: 0.0,
            known: Vec::new(),
            pre_sweeps: d.pre_sweeps,
            post_sweeps: d.post_sweeps,
            tolerance: d.tolerance,
            max_cycles: d.max_cycles,
        }
    }
}

/// The multigrid solver on a hierarchy.
pub struct MultigridHandle {
    solver: Multigrid,
    n_known: usize,
}

fn multigrid_error(e: ddg_intrinsic_simplification::multigrid::MultigridError) -> UnwrapError {
    use ddg_intrinsic_simplification::multigrid::MultigridError as E;
    match e {
        E::SizeMismatch { .. }
        | E::BadKnown { .. }
        | E::NotSquare { .. }
        | E::NotSymmetric { .. }
        | E::BadOption { .. }
        | E::NonFinite { .. } => message(format!("multigrid: {e}")),
        other => UnwrapError::Solve(other.to_string()),
    }
}

impl HierarchyHandle {
    /// Builds the hierarchy of a mesh.
    pub fn new(
        input: &MeshInput,
        args: &HierarchyArguments,
    ) -> Result<HierarchyHandle, UnwrapError> {
        let Mesh { mesh, geometry, .. } = build(input)?;
        if !(0.0..=1.0).contains(&args.weights) || args.weights.is_nan() {
            return Err(bad_option(
                "weights",
                args.weights.to_string(),
                "must be in [0, 1] (0: curvature, 1: area)",
            ));
        }
        let coarse = CoarseMesh::from_mesh(&mesh, &geometry, &[], CoarseMeshOptions::default())
            .map_err(iem)?;
        let options = HierarchyOptions {
            ratio: args.ratio,
            n_min: args.n_min,
            fresh_levels: args.fresh_levels,
            simplify: SimplifyOptions {
                channels: ChannelWeights::mixture(args.weights),
                metric: MetricOptions {
                    memory: args.memory,
                    ..MetricOptions::default()
                },
                fixed: args.fixed.clone(),
                ..SimplifyOptions::default()
            },
        };
        let hierarchy = build_hierarchy(&coarse, &options).map_err(|e| {
            use ddg_intrinsic_simplification::hierarchy::HierarchyError as E;
            match e {
                E::BadOptions { reason } => message(format!("hierarchy: {reason}")),
                E::Simplify(inner) => iem(inner),
                other => UnwrapError::Solve(other.to_string()),
            }
        })?;
        Ok(HierarchyHandle { hierarchy })
    }

    /// Number of levels.
    pub fn n_levels(&self) -> usize {
        self.hierarchy.n_levels()
    }

    /// Vertices per level, finest first.
    pub fn level_sizes(&self) -> Vec<usize> {
        self.hierarchy
            .levels()
            .iter()
            .map(|l| l.n_vertices())
            .collect()
    }

    fn level(
        &self,
        k: usize,
    ) -> Result<&ddg_intrinsic_simplification::hierarchy::Level, UnwrapError> {
        self.hierarchy.levels().get(k).ok_or_else(|| {
            message(format!(
                "level = {k}: out of range ({} levels)",
                self.hierarchy.n_levels()
            ))
        })
    }

    /// The vertices of level `k` as input vertex ids.
    pub fn level_vertices(&self, k: usize) -> Result<Vec<u32>, UnwrapError> {
        Ok(self.level(k)?.vertices.clone())
    }

    /// The cotan Laplacian of level `k`.
    pub fn laplacian(&self, k: usize) -> Result<Coo, UnwrapError> {
        Ok(Coo::of(&self.level(k)?.laplacian))
    }

    /// The lumped masses of level `k`.
    pub fn mass(&self, k: usize) -> Result<Vec<f64>, UnwrapError> {
        Ok(self.level(k)?.mass.clone())
    }

    /// The faces of level `k` (compacted ids, 3 per face).
    pub fn faces(&self, k: usize) -> Result<Vec<u32>, UnwrapError> {
        Ok(self.level(k)?.faces.iter().flatten().copied().collect())
    }

    /// `P_k` (level `k` by level `k + 1`).
    pub fn prolongation(&self, k: usize) -> Result<Coo, UnwrapError> {
        self.hierarchy
            .prolongations()
            .get(k)
            .map(Coo::of)
            .ok_or_else(|| {
                message(format!(
                    "level = {k}: out of range ({} prolongations)",
                    self.hierarchy.prolongations().len()
                ))
            })
    }

    /// The report of the construction.
    pub fn report(&self) -> &HierarchyReport {
        self.hierarchy.report()
    }

    /// The multigrid solver for a system on the finest level.
    pub fn multigrid(&self, args: &MultigridArguments) -> Result<MultigridHandle, UnwrapError> {
        let n = self.hierarchy.levels()[0].n_vertices();
        let a = match &args.matrix {
            Some(c) => {
                if c.nrows != n || c.ncols != n {
                    return Err(message(format!(
                        "matrix: {} x {} for {n} vertices on the finest level",
                        c.nrows, c.ncols
                    )));
                }
                let triplets: Vec<(u32, u32, f64)> = (0..c.rows.len())
                    .map(|i| (c.rows[i], c.cols[i], c.values[i]))
                    .collect();
                if let Some(&(r, col, _)) = triplets
                    .iter()
                    .find(|t| t.0 as usize >= n || t.1 as usize >= n)
                {
                    return Err(message(format!(
                        "matrix: entry ({r}, {col}) is outside the {n} x {n} matrix"
                    )));
                }
                SparseMatrix::from_triplets(n, n, &triplets)
                    .map_err(|e| message(format!("matrix: {e}")))?
            }
            None => {
                if !(args.shift.is_finite() && args.shift >= 0.0) {
                    return Err(bad_option(
                        "shift",
                        args.shift.to_string(),
                        "must be finite and >= 0",
                    ));
                }
                let level = &self.hierarchy.levels()[0];
                let mut a = level.laplacian.clone();
                let d: Vec<f64> = level.mass.iter().map(|m| m * args.shift).collect();
                a.add_to_diagonal(&d)
                    .map_err(|e| UnwrapError::Solve(format!("shift: {e}")))?;
                a
            }
        };
        let options = MultigridOptions {
            pre_sweeps: args.pre_sweeps,
            post_sweeps: args.post_sweeps,
            tolerance: args.tolerance,
            max_cycles: args.max_cycles,
            ..MultigridOptions::default()
        };
        let solver = self
            .hierarchy
            .multigrid(&a, &args.known, options)
            .map_err(multigrid_error)?;
        Ok(MultigridHandle {
            solver,
            n_known: args.known.len(),
        })
    }
}

impl MultigridHandle {
    fn known_values(&self, v: Option<&[f64]>) -> Result<Vec<f64>, UnwrapError> {
        match v {
            Some(v) => Ok(v.to_vec()),
            None => Ok(vec![0.0; self.n_known]),
        }
    }

    /// V-cycles.
    pub fn solve(
        &self,
        b: &[f64],
        known_values: Option<&[f64]>,
        x0: Option<&[f64]>,
    ) -> Result<(Vec<f64>, SolveReport), UnwrapError> {
        let k = self.known_values(known_values)?;
        self.solver.solve(b, &k, x0).map_err(multigrid_error)
    }

    /// Conjugate gradients preconditioned by a V-cycle.
    pub fn pcg(
        &self,
        b: &[f64],
        known_values: Option<&[f64]>,
        x0: Option<&[f64]>,
        max_iterations: usize,
    ) -> Result<(Vec<f64>, SolveReport), UnwrapError> {
        let k = self.known_values(known_values)?;
        self.solver
            .pcg(b, &k, x0, max_iterations)
            .map_err(multigrid_error)
    }

    /// What the setup built.
    pub fn info(&self) -> &MultigridInfo {
        self.solver.info()
    }
}

// ---------------------------------------------------------------------------------------------
// QEM
// ---------------------------------------------------------------------------------------------

/// The arguments of [`qem`].
#[derive(Debug, Clone)]
pub struct QemArguments {
    /// Vertices to stop at.
    pub target: usize,
    /// Weight the planes by area.
    pub area_weighted: bool,
    /// Weight of the boundary planes.
    pub boundary_weight: f64,
}

/// The output of [`qem`].
#[derive(Debug, Clone)]
pub struct QemOutput {
    /// 3 reals per coarse vertex.
    pub positions: Vec<f64>,
    /// 3 vertex ids per face.
    pub triangles: Vec<u32>,
    /// Per fine vertex, the coarse vertex it collapsed into.
    pub vertex_map: Vec<u32>,
    /// Per fine vertex, the closest coarse face.
    pub projection_face: Vec<u32>,
    /// 3 barycentric coordinates per fine vertex.
    pub projection_bary: Vec<f64>,
    /// Distance per fine vertex.
    pub projection_distance: Vec<f64>,
    /// The report.
    pub report: QemReport,
}

/// Quadric error metric simplification (Garland, Heckbert 1997) with a closest-point map.
pub fn qem(input: &MeshInput, args: &QemArguments) -> Result<QemOutput, UnwrapError> {
    let Mesh {
        mesh, positions, ..
    } = build(input)?;
    let mut options = QemOptions::new(args.target);
    options.area_weighted = args.area_weighted;
    options.boundary_weight = args.boundary_weight;
    let r = qem_simplify(&positions, mesh.faces(), &options).map_err(|e| {
        use ddg_intrinsic_simplification::baseline::QemError as E;
        match e {
            E::InvalidMesh(m) => UnwrapError::Mesh(format!("triangles: {m}")),
            E::NonFinitePosition { .. } => message(e.to_string()),
            E::InvalidOptions(_) => message(e.to_string()),
        }
    })?;
    Ok(QemOutput {
        positions: r.positions.iter().flatten().copied().collect(),
        triangles: r.faces.iter().flatten().copied().collect(),
        vertex_map: r.vertex_map,
        projection_face: r.projection.iter().map(|p| p.face).collect(),
        projection_bary: r.projection.iter().flat_map(|p| p.bary).collect(),
        projection_distance: r.projection.iter().map(|p| p.distance).collect(),
        report: r.report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::generators::icosphere;

    fn sphere(level: u32) -> MeshInput {
        let m = icosphere(level).unwrap();
        MeshInput {
            positions: m.positions.iter().flatten().copied().collect(),
            triangles: m.faces.iter().flatten().copied().collect(),
        }
    }

    fn centroids(input: &MeshInput) -> (Vec<u32>, Vec<f64>) {
        let nf = input.triangles.len() / 3;
        ((0..nf as u32).collect(), vec![1.0 / 3.0; 3 * nf])
    }

    #[test]
    fn simplify_sphere_with_centroids() {
        let input = sphere(3);
        let n = input.positions.len() / 3;
        let (pf, pb) = centroids(&input);
        let args = SimplifyArguments {
            ratio: Some(0.25),
            point_faces: pf,
            point_bary: pb,
            ..SimplifyArguments::default()
        };
        let s = Simplification::new(&input, &args).unwrap();
        assert!(s.output.n_vertices() <= n / 4 + 1);
        let loc = s.point_locations();
        assert_eq!(loc.face.len(), n + input.triangles.len() / 3);
        assert!(loc.face.iter().all(|&f| f >= 0));
        let (colours, nc) = s.face_colours();
        assert_eq!(colours.len(), s.output.n_faces());
        assert!(nc >= 1);
        let p = s.prolongation().unwrap();
        let mut sums = vec![0.0; n];
        for (r, v) in p.rows.iter().zip(&p.values) {
            sums[*r as usize] += v;
        }
        assert!(sums.iter().all(|x| (x - 1.0).abs() < 1e-9));
        let v = s.vector_prolongation(1).unwrap();
        assert_eq!(v.report.n_default, 0);
    }

    #[test]
    fn target_and_ratio_errors_name_the_parameters() {
        let input = sphere(1);
        let e = Simplification::new(&input, &SimplifyArguments::default())
            .err()
            .unwrap();
        assert!(e.to_string().contains("target or ratio"), "{e}");
        let both = SimplifyArguments {
            target: Some(10),
            ratio: Some(0.5),
            ..SimplifyArguments::default()
        };
        assert!(
            Simplification::new(&input, &both)
                .err()
                .unwrap()
                .to_string()
                .contains("not both")
        );
        let bad = SimplifyArguments {
            ratio: Some(1.5),
            ..SimplifyArguments::default()
        };
        assert!(
            Simplification::new(&input, &bad)
                .err()
                .unwrap()
                .to_string()
                .starts_with("ratio")
        );
    }

    #[test]
    fn fixed_vertices_stay_and_poisson_runs() {
        let input = sphere(3);
        let n = input.positions.len() / 3;
        let args = SimplifyArguments {
            target: Some(n / 3),
            fixed: vec![0, 5, 7],
            ..SimplifyArguments::default()
        };
        let mut s = Simplification::new(&input, &args).unwrap();
        for v in [0usize, 5, 7] {
            assert!(s.output.coarse_of_fine[v].is_some());
        }
        let f: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
        let out = s.poisson(&f, &[], &[]).unwrap();
        assert_eq!(out.fine.len(), n);
        let (_, d) = s.geodesic(0, "heat", 1.0).unwrap();
        assert_eq!(d.len(), n);
        let v = s.voronoi(&[0, 5, 7], 2).unwrap();
        assert_eq!(v.fine_labels.len(), n);
    }

    #[test]
    fn hierarchy_and_multigrid_solve() {
        let input = sphere(3);
        let h = HierarchyHandle::new(
            &input,
            &HierarchyArguments {
                n_min: 20,
                ..HierarchyArguments::default()
            },
        )
        .unwrap();
        assert!(h.n_levels() >= 2);
        let n = h.level_sizes()[0];
        let mg = h
            .multigrid(&MultigridArguments {
                shift: 1.0,
                ..MultigridArguments::default()
            })
            .unwrap();
        let b: Vec<f64> = (0..n).map(|i| (i as f64).cos()).collect();
        let (x, report) = mg.pcg(&b, None, None, 100).unwrap();
        assert!(report.converged);
        assert_eq!(x.len(), n);
    }

    #[test]
    fn qem_runs() {
        let input = sphere(2);
        let o = qem(
            &input,
            &QemArguments {
                target: 20,
                area_weighted: false,
                boundary_weight: 1000.0,
            },
        )
        .unwrap();
        assert_eq!(o.positions.len() / 3, 20);
    }
}
