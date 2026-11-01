//! The intrinsic mesh hierarchy of the paper's Sec. 8.7 (SPEC §9.1, S24): a sequence of coarse
//! intrinsic triangulations of one input mesh, with the cotan Laplacian and lumped masses of
//! every level and the scalar prolongations between consecutive levels, ready for
//! [`Multigrid`].
//!
//! > Surface Simplification using Intrinsic Error Metrics (Liu, Gillespie, Chislett, Sharp,
//! > Jacobson, Crane 2023), Sec. 8.7, Fig. 29-31; and Liu, Zhang, Ben-Chen, Jacobson,
//! > *Surface Multigrid via Intrinsic Prolongation* (2021), Sec. 5 and App. D for the defaults
//! > (ratio 1/4, coarsest level of at most 500 vertices).
//!
//! # Construction
//!
//! [`build_hierarchy`] runs [`simplify`] **once** with
//! [`snapshots`](SimplifyOptions::snapshots) at `n₀ ρ, n₀ ρ², …` (each count the rounded `ρ ×`
//! the previous), so the masses and error vectors carry across the levels (S24). Level `k + 1`
//! is the [`Snapshot`] at that count; the prolongation `P_k ∈ ℝ^{n_k × n_{k+1}}` is read from
//! the locations of the level-`k` vertices in level `k + 1` (`Snapshot::previous`, tracked
//! through every flip and removal): a vertex that is also a vertex of level `k + 1` has one
//! entry `1`, any other vertex has the barycentric coordinates of its coarse face (summed on a
//! corner that repeats in a Δ-face), so every row sums to one and `P_k 1 = 1` (constants are
//! reproduced). Level `0` is the input as given (its Laplacian is that of the input
//! connectivity, before the initial Delaunay flipping of the run).
//!
//! The level counts follow the reference (`mg_precompute.cpp`; SPEC R33): levels are added
//! while the current one has more than `n_min` vertices, so the coarsest level is the first
//! one with at most `n_min` vertices (spot: 2930, 733, 183).
//!
//! With [`HierarchyOptions::fresh_levels`] the reference's variant is reproduced instead: one
//! `simplify` call per level, each starting from the previous level with **fresh** masses and
//! error vectors from the current curvature (a new call recomputes them; no memory crosses
//! levels, R33). The prolongations are read the same way (each call takes one snapshot at its
//! target, whose `previous` is the vertex set the call started with).
//!
//! # Level operators
//!
//! [`Level::laplacian`] is the cotan Laplacian of the level's lengths (SPEC §2.3, `CLAUDE.md`
//! conventions: symmetric positive semidefinite), [`Level::mass`] the lumped masses, both
//! computed from the snapshot's faces and lengths so that they are safe on Δ-complexes (a face
//! with a repeated vertex adds its cotangent to the diagonal of the repeated vertex only for the
//! non-degenerate edges; a self-edge contributes nothing, as in
//! [`CoarseMesh::cotan_laplacian`](crate::CoarseMesh::cotan_laplacian)).

use std::collections::HashMap;
use std::time::Instant;

use ddg_core::geometry::{corner_cotan, triangle_area};
use ddg_core::{SparseError, SparseMatrix};

use crate::PointLocation;
use crate::error::IemError;
use crate::mesh::CoarseMesh;
use crate::multigrid::{Multigrid, MultigridError, MultigridOptions};
use crate::simplify::{SimplifyOptions, SimplifyReport, Snapshot, Target, simplify};

/// Errors of [`build_hierarchy`].
#[derive(Debug, Clone, PartialEq)]
pub enum HierarchyError {
    /// An option is out of range, or the simplification refused its options or failed.
    Simplify(IemError),
    /// The coarsening ratio or the coarsest size is not valid.
    BadOptions {
        /// What is wrong.
        reason: String,
    },
    /// A level has a face that is not a triangle by its lengths, or a non-finite length.
    BadLevel {
        /// The level.
        level: usize,
        /// The face.
        face: usize,
        /// What is wrong.
        reason: String,
    },
    /// A snapshot is inconsistent (a located vertex that is not a vertex of the snapshot).
    Inconsistent {
        /// What is wrong.
        reason: String,
    },
    /// A failure of the sparse backend.
    Sparse(SparseError),
}

impl std::fmt::Display for HierarchyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HierarchyError::Simplify(e) => write!(f, "hierarchy: {e}"),
            HierarchyError::BadOptions { reason } => {
                write!(f, "invalid hierarchy options: {reason}")
            }
            HierarchyError::BadLevel {
                level,
                face,
                reason,
            } => write!(f, "hierarchy level {level}, face {face}: {reason}"),
            HierarchyError::Inconsistent { reason } => {
                write!(f, "inconsistent hierarchy snapshot: {reason}")
            }
            HierarchyError::Sparse(e) => write!(f, "hierarchy: sparse backend: {e}"),
        }
    }
}

impl std::error::Error for HierarchyError {}

impl From<IemError> for HierarchyError {
    fn from(e: IemError) -> Self {
        HierarchyError::Simplify(e)
    }
}

impl From<SparseError> for HierarchyError {
    fn from(e: SparseError) -> Self {
        HierarchyError::Sparse(e)
    }
}

/// Options of [`build_hierarchy`] (SPEC §9.1).
#[derive(Debug, Clone, PartialEq)]
pub struct HierarchyOptions {
    /// The vertex count of a level over the one above, in `(0, 1)` (default `0.25`, the paper's
    /// and the reference's).
    pub ratio: f64,
    /// Levels are added while the current one has more than this many vertices (default 500; the
    /// coarsest level is the first with at most `n_min`, as the reference).
    pub n_min: usize,
    /// One `simplify` call per level with fresh masses (the reference's variant, R33) instead
    /// of one run with snapshots (S24, default `false`).
    pub fresh_levels: bool,
    /// The options of the simplification (channels, metric, fixed vertices, ...). `target`,
    /// `snapshots` and `record_sequence` are set by [`build_hierarchy`].
    pub simplify: SimplifyOptions,
}

impl Default for HierarchyOptions {
    fn default() -> Self {
        HierarchyOptions {
            ratio: 0.25,
            n_min: 500,
            fresh_levels: false,
            simplify: SimplifyOptions::default(),
        }
    }
}

/// One level of the hierarchy.
#[derive(Debug, Clone)]
pub struct Level {
    /// The vertices as fine ids (index = row of this level's matrices).
    pub vertices: Vec<u32>,
    /// The faces (indices into `vertices`; a Δ-face may repeat one).
    pub faces: Vec<[u32; 3]>,
    /// Per face, `lengths[f][k]` is the length of the side from slot `k` to slot `k + 1`.
    pub lengths: Vec<[f64; 3]>,
    /// The cotan Laplacian (symmetric, positive semidefinite, zero row sums).
    pub laplacian: SparseMatrix,
    /// The lumped masses `⅓ Σ A_f` (unnormalized; they sum to the level's area).
    pub mass: Vec<f64>,
}

impl Level {
    /// Number of vertices.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }
}

/// What [`build_hierarchy`] did.
#[derive(Debug, Clone)]
pub struct HierarchyReport {
    /// The vertex counts asked for, finest first (the first is the input's).
    pub targets: Vec<usize>,
    /// The vertex counts of the levels built.
    pub counts: Vec<usize>,
    /// Whether the run was one run with snapshots (`false`) or one run per level (`true`).
    pub fresh_levels: bool,
    /// The report of each `simplify` call (one, or one per level).
    pub simplify: Vec<SimplifyReport>,
    /// Stored entries of each `P_k`.
    pub prolongation_nnz: Vec<usize>,
    /// Seconds in `simplify` (all calls).
    pub simplify_seconds: f64,
    /// Seconds building the level operators and prolongations.
    pub operators_seconds: f64,
}

/// The hierarchy: levels finest first and `P_k` (`n_k × n_{k+1}`) between them.
#[derive(Debug, Clone)]
pub struct Hierarchy {
    levels: Vec<Level>,
    prolongations: Vec<SparseMatrix>,
    report: HierarchyReport,
}

/// The vertex counts of the levels: `n₀`, then `round(ρ n)` while `n > n_min`.
fn level_targets(n0: usize, ratio: f64, n_min: usize) -> Vec<usize> {
    let mut out = vec![n0];
    let mut n = n0;
    while n > n_min {
        let next = ((n as f64) * ratio).round() as usize;
        if next >= n || next < 4 {
            break;
        }
        out.push(next);
        n = next;
    }
    out
}

/// The cotan Laplacian and lumped masses of a triangulation given by faces and side lengths
/// (`lengths[f][k]` from slot `k` to slot `k + 1`); self-edges contribute nothing.
///
/// [formula: cotan-laplacian] [formula: lumped-mass]
fn level_operators(
    n: usize,
    faces: &[[u32; 3]],
    lengths: &[[f64; 3]],
    level: usize,
) -> Result<(SparseMatrix, Vec<f64>), HierarchyError> {
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(6 * faces.len() + n);
    let mut diagonal = vec![0.0; n];
    let mut mass = vec![0.0; n];
    for (f, (face, l)) in faces.iter().zip(lengths).enumerate() {
        let area = triangle_area(l[0], l[1], l[2]).ok_or_else(|| HierarchyError::BadLevel {
            level,
            face: f,
            reason: format!("sides {l:?} do not form a triangle"),
        })?;
        for &v in face {
            mass[v as usize] += area / 3.0;
        }
        for k in 0..3 {
            let (i, j) = (face[k] as usize, face[(k + 1) % 3] as usize);
            if i == j {
                continue;
            }
            let cot = corner_cotan(l[k], l[(k + 1) % 3], l[(k + 2) % 3]).ok_or_else(|| {
                HierarchyError::BadLevel {
                    level,
                    face: f,
                    reason: format!("sides {l:?} do not form a triangle"),
                }
            })?;
            let w = 0.5 * cot;
            triplets.push((i as u32, j as u32, -w));
            triplets.push((j as u32, i as u32, -w));
            diagonal[i] += w;
            diagonal[j] += w;
        }
    }
    for (v, d) in diagonal.into_iter().enumerate() {
        triplets.push((v as u32, v as u32, d));
    }
    Ok((SparseMatrix::from_triplets(n, n, &triplets)?, mass))
}

/// `P` (`n_prev × n_this`) from the locations of the previous vertices in this snapshot.
fn prolongation_of(
    snapshot: &Snapshot,
    previous_rows: &HashMap<u32, usize>,
    n_prev: usize,
) -> Result<SparseMatrix, HierarchyError> {
    let column: HashMap<u32, u32> = snapshot
        .vertices
        .iter()
        .enumerate()
        .map(|(c, &v)| (v, c as u32))
        .collect();
    let col = |fine: u32| {
        column
            .get(&fine)
            .copied()
            .ok_or_else(|| HierarchyError::Inconsistent {
                reason: format!("fine vertex {fine} is located at a corner that is not a vertex"),
            })
    };
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(3 * n_prev);
    for (&v, loc) in snapshot.previous_vertices.iter().zip(&snapshot.previous) {
        let row = *previous_rows
            .get(&v)
            .ok_or_else(|| HierarchyError::Inconsistent {
                reason: format!("fine vertex {v} was not a vertex of the previous level"),
            })? as u32;
        match *loc {
            PointLocation::Vertex(w) => triplets.push((row, col(w)?, 1.0)),
            PointLocation::Face { corners, bary, .. } => {
                for k in 0..3 {
                    triplets.push((row, col(corners[k])?, bary[k]));
                }
            }
        }
    }
    let p = SparseMatrix::from_triplets(n_prev, snapshot.vertices.len(), &triplets)?;
    Ok(p)
}

/// Builds the intrinsic hierarchy of `mesh` (SPEC §9.1). `mesh` is not changed (the run is on a
/// copy). Its vertices at the start are level 0 (normally the whole input); fixed vertices of
/// the mesh and of `options.simplify.fixed` stay in every level.
///
/// # Errors
/// [`HierarchyError`]: a ratio outside `(0, 1)`, a failure or refusal of `simplify`, a level
/// whose lengths are not a triangulation, or a sparse backend failure.
pub fn build_hierarchy(
    mesh: &CoarseMesh,
    options: &HierarchyOptions,
) -> Result<Hierarchy, HierarchyError> {
    if !(options.ratio > 0.0 && options.ratio < 1.0) {
        return Err(HierarchyError::BadOptions {
            reason: format!("the ratio must be in (0, 1), got {}", options.ratio),
        });
    }
    let n0 = mesh.n_vertices();
    let targets = level_targets(n0, options.ratio, options.n_min);
    let mut work = mesh.clone();
    let mut sim = options.simplify.clone();
    sim.record_sequence = false;

    // Level 0: the input as given, before any flip.
    let t_ops = Instant::now();
    let mut operators_seconds = 0.0;
    let vertices0: Vec<u32> = work.vertices().to_vec();
    let faces0: Vec<[u32; 3]> = {
        let index: HashMap<u32, u32> = vertices0
            .iter()
            .enumerate()
            .map(|(c, &v)| (v, c as u32))
            .collect();
        work.faces().iter().map(|f| f.map(|v| index[&v])).collect()
    };
    let lengths0: Vec<[f64; 3]> = (0..work.n_faces())
        .map(|f| crate::mesh::face_lengths(work.complex(), work.edge_lengths(), f))
        .collect();
    let (l0, m0) = level_operators(n0, &faces0, &lengths0, 0)?;
    let mut levels = vec![Level {
        vertices: vertices0,
        faces: faces0,
        lengths: lengths0,
        laplacian: l0,
        mass: m0,
    }];
    operators_seconds += t_ops.elapsed().as_secs_f64();

    let mut prolongations: Vec<SparseMatrix> = Vec::new();
    let mut reports: Vec<SimplifyReport> = Vec::new();
    let mut simplify_seconds = 0.0;
    let mut snapshots: Vec<Snapshot> = Vec::new();
    if targets.len() > 1 {
        if options.fresh_levels {
            let mut n_now = n0;
            for _ in 1..targets.len() {
                let target = ((n_now as f64) * options.ratio).round() as usize;
                let mut s = sim.clone();
                s.target = Target::Vertices(target);
                s.snapshots = vec![target];
                let t = Instant::now();
                let report = simplify(&mut work, &s)?;
                simplify_seconds += t.elapsed().as_secs_f64();
                let got = report.snapshots.len();
                reports.push(report);
                let Some(snap) = reports.last().and_then(|r| r.snapshots.last()).cloned() else {
                    let _ = got;
                    break;
                };
                n_now = snap.vertices.len();
                snapshots.push(snap);
                if n_now <= options.n_min {
                    break;
                }
            }
        } else {
            let mut s = sim.clone();
            s.target = Target::Vertices(*targets.last().expect("non-empty"));
            s.snapshots = targets[1..].to_vec();
            let t = Instant::now();
            let report = simplify(&mut work, &s)?;
            simplify_seconds += t.elapsed().as_secs_f64();
            snapshots = report.snapshots.clone();
            reports.push(report);
        }
    }

    let t_ops = Instant::now();
    for (k, snap) in snapshots.iter().enumerate() {
        let prev_rows: HashMap<u32, usize> = levels[k]
            .vertices
            .iter()
            .enumerate()
            .map(|(r, &v)| (v, r))
            .collect();
        let p = prolongation_of(snap, &prev_rows, levels[k].n_vertices())?;
        let (lap, mass) = level_operators(snap.vertices.len(), &snap.faces, &snap.lengths, k + 1)?;
        prolongations.push(p);
        levels.push(Level {
            vertices: snap.vertices.clone(),
            faces: snap.faces.clone(),
            lengths: snap.lengths.clone(),
            laplacian: lap,
            mass,
        });
    }
    operators_seconds += t_ops.elapsed().as_secs_f64();
    let report = HierarchyReport {
        targets,
        counts: levels.iter().map(Level::n_vertices).collect(),
        fresh_levels: options.fresh_levels,
        prolongation_nnz: prolongations.iter().map(SparseMatrix::nnz).collect(),
        simplify: reports,
        simplify_seconds,
        operators_seconds,
    };
    Ok(Hierarchy {
        levels,
        prolongations,
        report,
    })
}

impl Hierarchy {
    /// The levels, finest (the input) first.
    pub fn levels(&self) -> &[Level] {
        &self.levels
    }

    /// Number of levels.
    pub fn n_levels(&self) -> usize {
        self.levels.len()
    }

    /// `P_k ∈ ℝ^{n_k × n_{k+1}}`, `k = 0 … n_levels − 2`.
    pub fn prolongations(&self) -> &[SparseMatrix] {
        &self.prolongations
    }

    /// What the construction did.
    pub fn report(&self) -> &HierarchyReport {
        &self.report
    }

    /// The Galerkin multigrid solver for the fine system `a` (`n_0 × n_0`, e.g. `L + εM` of
    /// level 0) with the fine vertices `known` fixed (SPEC §9.2).
    ///
    /// # Errors
    /// [`MultigridError`] as [`Multigrid::new`].
    pub fn multigrid(
        &self,
        a: &SparseMatrix,
        known: &[u32],
        options: MultigridOptions,
    ) -> Result<Multigrid, MultigridError> {
        Multigrid::new(a, &self.prolongations, known, options)
    }
}
