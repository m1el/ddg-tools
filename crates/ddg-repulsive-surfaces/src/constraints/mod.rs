//! Constraints: barycenter, pins, area, volume; values, Jacobians, the exact projection of the
//! positional ones, the translation and rank checks (`SPEC.md` §5.1, D8–D10).
//!
//! A [`ConstraintSet`] is built from the caller's [`Constraint`] list at the initial positions:
//! `None` targets become the current values, indices are checked, a volume needs consistently
//! oriented components, every component's translations must be fixed by a positional row
//! ([`RsError::TranslationNotFixed`]) and linearly dependent rows are dropped and reported
//! ([`rank`], D10). It then evaluates `Φ` ([`ConstraintSet::values`]) and the Jacobian
//! ([`ConstraintSet::jacobian`], `k × 3|V|`, column `3 v + axis`) at any positions.
//!
//! | constraint | rows | kind |
//! |---|---|---|
//! | `Barycenter` | 3 per component (axis fastest): `Σ a_v f_v / A_c − X₀`, Jacobian `a_v / A_c ⊗ I₃` with the areas frozen (D8) | positional |
//! | `Pins` | 3 per pin (axis fastest): `f_v − f_v⁰`, Jacobian `e_v ⊗ I₃` | positional |
//! | `Area` | 1: `Σ a_S − A₀`, Jacobian `Σ_{S∋v} ½ N_S × e_opp` | nonlinear (Schur) |
//! | `Volume` | 1: `⅙ Σ f_i·(f_j × f_k) − V₀`, Jacobian `⅙ Σ_{S∋v} f_j × f_k` (exact, D9) | nonlinear (Schur) |
//!
//! Rows are listed in constraint order; a barycenter with `component: None` gives one group of
//! three rows per component, in component order ([`BarycenterMode::PerComponent`], the paper)
//! or one group for the whole mesh ([`BarycenterMode::Global`], the reference).

pub mod rank;

pub use rank::{RANK_TOLERANCE, RankReport};

use ddg_core::SparseMatrix;

use crate::error::RsError;
use crate::geometry::v3::{cross, dot, scale};
use crate::geometry::{FaceGeometry, corners, face_area_gradient};
use crate::surface::Surface;

/// A hard constraint (`SPEC.md` §5.1, §9.1).
#[derive(Debug, Clone, PartialEq)]
pub enum Constraint {
    /// The area-weighted barycenter `Σ a_v f_v / Σ a_v` of a component (`None`: every component
    /// separately, or the whole mesh in [`BarycenterMode::Global`]). `target: None` = the
    /// initial value (per component); a target with `component: None` applies to every
    /// component.
    Barycenter {
        /// The component (numbered by lowest face index), or all.
        component: Option<u32>,
        /// The target, or the initial value.
        target: Option<[f64; 3]>,
    },
    /// Fixed vertex positions.
    Pins {
        /// The pinned vertices.
        vertices: Vec<u32>,
        /// One target per vertex, or the initial positions.
        targets: Option<Vec<[f64; 3]>>,
    },
    /// The total area of a component (`None`: the whole mesh).
    Area {
        /// The component, or all.
        component: Option<u32>,
        /// The target, or the initial value.
        target: Option<f64>,
    },
    /// The signed volume `⅙ Σ f_i·(f_j × f_k)` of a component (`None`: the whole mesh). Needs
    /// [`crate::Orientation::Consistent`] components. On an open component it is the signed
    /// volume of the cone over the origin.
    Volume {
        /// The component, or all.
        component: Option<u32>,
        /// The target, or the initial value.
        target: Option<f64>,
    },
}

/// How `Barycenter { component: None }` is resolved (`SPEC.md` §5.1, D8).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BarycenterMode {
    /// One barycenter per component (the paper, Sec. 7.1.1, Fig. 8).
    #[default]
    PerComponent,
    /// One barycenter of the whole mesh (the reference). On a mesh with several components the
    /// relative translations of the components stay in the kernel of `A′` (the reference
    /// regularizes its preconditioner with `1e-2 M` instead); the translation check accepts the
    /// global row for every component in this mode, and a dense solve then reports
    /// [`RsError::SingularSaddle`].
    Global,
}

/// What one row constrains.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowKind {
    /// One coordinate of a barycenter (`component: None` = the whole mesh, global mode).
    Barycenter {
        /// The component.
        component: Option<u32>,
        /// The axis (0, 1, 2).
        axis: u8,
    },
    /// One coordinate of a pinned vertex.
    Pin {
        /// The vertex.
        vertex: u32,
        /// The axis (0, 1, 2).
        axis: u8,
    },
    /// The area of a component (or the whole mesh).
    Area {
        /// The component, or all.
        component: Option<u32>,
    },
    /// The signed volume of a component (or the whole mesh).
    Volume {
        /// The component, or all.
        component: Option<u32>,
    },
}

impl RowKind {
    /// Positional rows live inside `A′`; the others go through the Schur complement (§5.1).
    pub fn is_positional(&self) -> bool {
        matches!(self, RowKind::Barycenter { .. } | RowKind::Pin { .. })
    }

    /// `"barycenter"`, `"pins"`, `"area"`, `"volume"`.
    pub fn name(&self) -> &'static str {
        match self {
            RowKind::Barycenter { .. } => "barycenter",
            RowKind::Pin { .. } => "pins",
            RowKind::Area { .. } => "area",
            RowKind::Volume { .. } => "volume",
        }
    }
}

/// One scalar constraint row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Row {
    /// Index of the constraint in the caller's list.
    pub constraint: usize,
    /// What the row constrains.
    pub kind: RowKind,
    /// The target value (`Φ = value − target`).
    pub target: f64,
}

/// The relative tolerance of the nonlinear rows (`SPEC.md` §6.2): `|Φ_r| ≤ 1e-8 |target|`, or
/// `1e-8` for a zero target.
pub const NONLINEAR_TOLERANCE: f64 = 1e-8;

/// The resolved constraints of one surface, their rows and the rows dropped by the rank check.
#[derive(Debug, Clone)]
pub struct ConstraintSet {
    n_vertices: usize,
    mode: BarycenterMode,
    faces: Vec<[u32; 3]>,
    face_component: Vec<u32>,
    vertex_component: Vec<u32>,
    n_components: usize,
    rows: Vec<Row>,
    active: Vec<bool>,
    rank: RankReport,
}

fn bad_target(name: &'static str, value: String) -> RsError {
    RsError::InvalidOption {
        name,
        value,
        requirement: "finite target values",
    }
}

impl ConstraintSet {
    /// Resolves `constraints` at `positions` (`SPEC.md` §5.1): `None` targets become the current
    /// values; indices and orientation are checked; every component's translations must be fixed
    /// by a barycenter or a pin; linearly dependent rows are dropped (pins first, then
    /// barycenters, then area/volume in list order are kept; [`rank`]).
    ///
    /// # Errors
    /// The errors of [`Surface::validate_positions`]; [`RsError::ConstraintOutOfRange`] for a
    /// component or vertex that does not exist; [`RsError::InvalidOption`] for non-finite
    /// targets or a pin target list of the wrong length; [`RsError::VolumeNeedsOrientation`];
    /// [`RsError::InfeasibleConstraint`] for a vertex pinned twice to different targets;
    /// [`RsError::TranslationNotFixed`].
    pub fn new(
        surface: &Surface,
        positions: &[[f64; 3]],
        constraints: &[Constraint],
        mode: BarycenterMode,
    ) -> Result<ConstraintSet, RsError> {
        surface.validate_positions(positions)?;
        let mut set = ConstraintSet {
            n_vertices: surface.n_vertices(),
            mode,
            faces: surface.faces().to_vec(),
            face_component: surface.face_components().to_vec(),
            vertex_component: surface.vertex_components().to_vec(),
            n_components: surface.n_components(),
            rows: Vec::new(),
            active: Vec::new(),
            rank: RankReport::default(),
        };
        let nc = set.n_components as u32;
        let check_comp = |i: usize, c: Option<u32>| -> Result<(), RsError> {
            match c {
                Some(c) if c >= nc => Err(RsError::ConstraintOutOfRange {
                    constraint: i,
                    index: c,
                }),
                _ => Ok(()),
            }
        };
        let geom = FaceGeometry::new(surface, positions)?;
        let vareas = crate::geometry::vertex_areas(surface, &geom.areas);
        for (i, c) in constraints.iter().enumerate() {
            match c {
                Constraint::Barycenter { component, target } => {
                    check_comp(i, *component)?;
                    if let Some(t) = target
                        && !t.iter().all(|x| x.is_finite())
                    {
                        return Err(bad_target("barycenter target", format!("{t:?}")));
                    }
                    let groups: Vec<Option<u32>> = match (component, mode) {
                        (Some(c), _) => vec![Some(*c)],
                        (None, BarycenterMode::PerComponent) => (0..nc).map(Some).collect(),
                        (None, BarycenterMode::Global) => vec![None],
                    };
                    for g in groups {
                        let cur = set.barycenter(&vareas, positions, g);
                        let t = target.unwrap_or(cur);
                        for axis in 0..3u8 {
                            set.rows.push(Row {
                                constraint: i,
                                kind: RowKind::Barycenter {
                                    component: g,
                                    axis,
                                },
                                target: t[axis as usize],
                            });
                        }
                    }
                }
                Constraint::Pins { vertices, targets } => {
                    if let Some(&v) = vertices.iter().find(|&&v| v as usize >= set.n_vertices) {
                        return Err(RsError::ConstraintOutOfRange {
                            constraint: i,
                            index: v,
                        });
                    }
                    if let Some(t) = targets {
                        if t.len() != vertices.len() {
                            return Err(RsError::InvalidOption {
                                name: "pin targets",
                                value: format!("{} targets for {} pins", t.len(), vertices.len()),
                                requirement: "one target per pinned vertex",
                            });
                        }
                        if let Some(p) = t.iter().find(|p| !p.iter().all(|x| x.is_finite())) {
                            return Err(bad_target("pin target", format!("{p:?}")));
                        }
                    }
                    for (k, &v) in vertices.iter().enumerate() {
                        let t = targets.as_ref().map_or(positions[v as usize], |t| t[k]);
                        for axis in 0..3u8 {
                            set.rows.push(Row {
                                constraint: i,
                                kind: RowKind::Pin { vertex: v, axis },
                                target: t[axis as usize],
                            });
                        }
                    }
                }
                Constraint::Area { component, target } => {
                    check_comp(i, *component)?;
                    if let Some(t) = target
                        && !t.is_finite()
                    {
                        return Err(bad_target("area target", t.to_string()));
                    }
                    let t = target.unwrap_or_else(|| set.area_of(&geom.areas, *component));
                    set.rows.push(Row {
                        constraint: i,
                        kind: RowKind::Area {
                            component: *component,
                        },
                        target: t,
                    });
                }
                Constraint::Volume { component, target } => {
                    check_comp(i, *component)?;
                    for (k, comp) in surface.components().iter().enumerate() {
                        if component.is_none_or(|c| c as usize == k)
                            && let Some(edge) = comp.inconsistent_edge()
                        {
                            return Err(RsError::VolumeNeedsOrientation {
                                component: k as u32,
                                edge,
                            });
                        }
                    }
                    if let Some(t) = target
                        && !t.is_finite()
                    {
                        return Err(bad_target("volume target", t.to_string()));
                    }
                    let t = target.unwrap_or_else(|| set.volume_of(positions, *component));
                    set.rows.push(Row {
                        constraint: i,
                        kind: RowKind::Volume {
                            component: *component,
                        },
                        target: t,
                    });
                }
            }
        }
        set.active = vec![true; set.rows.len()];
        set.analyze(surface, positions)?;
        Ok(set)
    }

    /// Re-runs the rank check (`SPEC.md` §5.1, D10) at `positions` (all rows active first),
    /// e.g. once per step; the report is [`ConstraintSet::rank_report`].
    ///
    /// # Errors
    /// [`RsError::InfeasibleConstraint`], [`RsError::TranslationNotFixed`],
    /// [`RsError::PositionCountMismatch`].
    pub fn analyze(&mut self, surface: &Surface, positions: &[[f64; 3]]) -> Result<(), RsError> {
        self.check_surface(surface, positions)?;
        self.active = vec![true; self.rows.len()];
        let report = rank::analyze(self, positions)?;
        for &r in &report.dropped_rows {
            self.active[r] = false;
        }
        self.rank = report;
        rank::check_translations(self)?;
        Ok(())
    }

    fn check_surface(&self, surface: &Surface, positions: &[[f64; 3]]) -> Result<(), RsError> {
        if surface.n_vertices() != self.n_vertices || surface.faces() != self.faces.as_slice() {
            return Err(RsError::InvalidOption {
                name: "surface",
                value: format!("{} vertices, {} faces", surface.n_vertices(), surface.n_faces()),
                requirement: "the surface the constraint set was built for",
            });
        }
        if positions.len() != self.n_vertices {
            return Err(RsError::PositionCountMismatch {
                expected: self.n_vertices,
                found: positions.len(),
            });
        }
        Ok(())
    }

    /// All rows (active and dropped), in constraint order.
    pub fn rows(&self) -> &[Row] {
        &self.rows
    }

    /// Number of rows (active and dropped).
    pub fn n_rows(&self) -> usize {
        self.rows.len()
    }

    /// `|V|`.
    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }

    /// The barycenter mode.
    pub fn mode(&self) -> BarycenterMode {
        self.mode
    }

    /// Whether row `r` is active (not dropped by the rank check).
    pub fn is_active(&self, r: usize) -> bool {
        self.active[r]
    }

    /// Indices of the active positional rows (inside `A′`), increasing.
    pub fn positional_rows(&self) -> Vec<usize> {
        (0..self.rows.len())
            .filter(|&r| self.active[r] && self.rows[r].kind.is_positional())
            .collect()
    }

    /// Indices of the active nonlinear rows (Schur), increasing.
    pub fn nonlinear_rows(&self) -> Vec<usize> {
        (0..self.rows.len())
            .filter(|&r| self.active[r] && !self.rows[r].kind.is_positional())
            .collect()
    }

    /// Indices of the rows dropped by the rank check (D10), increasing.
    pub fn dropped_rows(&self) -> &[usize] {
        &self.rank.dropped_rows
    }

    /// The last rank check.
    pub fn rank_report(&self) -> &RankReport {
        &self.rank
    }

    /// Sets the target of row `r` (schedules, D11; moving pins).
    ///
    /// # Errors
    /// [`RsError::ConstraintOutOfRange`] for a row that does not exist;
    /// [`RsError::InvalidOption`] for a non-finite value.
    pub fn set_target(&mut self, r: usize, value: f64) -> Result<(), RsError> {
        if r >= self.rows.len() {
            return Err(RsError::ConstraintOutOfRange {
                constraint: r,
                index: r as u32,
            });
        }
        if !value.is_finite() {
            return Err(bad_target("row target", value.to_string()));
        }
        self.rows[r].target = value;
        Ok(())
    }

    /// The tolerance of row `r` in the Newton projection: `1e-8 |target|`, or `1e-8` for a zero
    /// target (`SPEC.md` §6.2). Positional rows are projected exactly; the same rule is reported
    /// for them.
    pub fn tolerance(&self, r: usize) -> f64 {
        let t = self.rows[r].target.abs();
        if t > 0.0 {
            NONLINEAR_TOLERANCE * t
        } else {
            NONLINEAR_TOLERANCE
        }
    }

    fn in_component(&self, v: usize, c: Option<u32>) -> bool {
        c.is_none_or(|c| self.vertex_component[v] == c)
    }

    fn barycenter(&self, vareas: &[f64], positions: &[[f64; 3]], c: Option<u32>) -> [f64; 3] {
        let mut s = [0.0; 3];
        let mut tot = 0.0;
        for v in 0..self.n_vertices {
            if self.in_component(v, c) {
                for k in 0..3 {
                    s[k] += vareas[v] * positions[v][k];
                }
                tot += vareas[v];
            }
        }
        scale(s, 1.0 / tot)
    }

    fn area_of(&self, areas: &[f64], c: Option<u32>) -> f64 {
        areas
            .iter()
            .zip(&self.face_component)
            .filter(|(_, fc)| c.is_none_or(|c| **fc == c))
            .map(|(a, _)| a)
            .sum()
    }

    /// The signed volume `⅙ Σ f_i·(f_j × f_k)` of a component (or all).
    ///
    /// [formula: rs-volume-gradient]
    fn volume_of(&self, positions: &[[f64; 3]], c: Option<u32>) -> f64 {
        let mut s = 0.0;
        for (t, fc) in self.faces.iter().zip(&self.face_component) {
            if c.is_none_or(|c| *fc == c) {
                let p = corners(positions, t);
                s += dot(p[0], cross(p[1], p[2])) / 6.0;
            }
        }
        s
    }

    /// Per-face areas and vertex areas at `positions` (no validation: degenerate faces give
    /// zero area).
    fn areas(&self, positions: &[[f64; 3]]) -> (Vec<f64>, Vec<f64>) {
        let fa: Vec<f64> = self
            .faces
            .iter()
            .map(|t| crate::geometry::face_quantities(corners(positions, t)).0)
            .collect();
        let mut va = vec![0.0; self.n_vertices];
        for (t, a) in self.faces.iter().zip(&fa) {
            for &v in t {
                va[v as usize] += a / 3.0;
            }
        }
        (fa, va)
    }

    fn check_len(&self, positions: &[[f64; 3]]) -> Result<(), RsError> {
        if positions.len() != self.n_vertices {
            return Err(RsError::PositionCountMismatch {
                expected: self.n_vertices,
                found: positions.len(),
            });
        }
        Ok(())
    }

    /// `Φ` of every row (active and dropped) at `positions`.
    ///
    /// [formula: rs-barycenter]
    /// [formula: rs-volume-gradient]
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`].
    pub fn values(&self, positions: &[[f64; 3]]) -> Result<Vec<f64>, RsError> {
        self.check_len(positions)?;
        let (fa, va) = self.areas(positions);
        let mut bary_cache: Vec<(Option<u32>, [f64; 3])> = Vec::new();
        let mut out = Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            let value = match row.kind {
                RowKind::Barycenter { component, axis } => {
                    let b = match bary_cache.iter().find(|e| e.0 == component) {
                        Some(e) => e.1,
                        None => {
                            let b = self.barycenter(&va, positions, component);
                            bary_cache.push((component, b));
                            b
                        }
                    };
                    b[axis as usize]
                }
                RowKind::Pin { vertex, axis } => positions[vertex as usize][axis as usize],
                RowKind::Area { component } => self.area_of(&fa, component),
                RowKind::Volume { component } => self.volume_of(positions, component),
            };
            out.push(value - row.target);
        }
        Ok(out)
    }

    /// The Jacobian of the listed rows at `positions` (`rows.len() × 3|V|`, column
    /// `3 v + axis`; the barycenter's areas frozen at `positions`, D8).
    ///
    /// [formula: rs-barycenter]
    /// [formula: rs-area-gradient]
    /// [formula: rs-volume-gradient]
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`]; [`RsError::ConstraintOutOfRange`] for a row index
    /// that does not exist.
    pub fn jacobian(&self, positions: &[[f64; 3]], rows: &[usize]) -> Result<SparseMatrix, RsError> {
        self.check_len(positions)?;
        if let Some(&r) = rows.iter().find(|&&r| r >= self.rows.len()) {
            return Err(RsError::ConstraintOutOfRange {
                constraint: r,
                index: r as u32,
            });
        }
        let nv = self.n_vertices;
        let (_, va) = self.areas(positions);
        let mut trip: Vec<(u32, u32, f64)> = Vec::new();
        for (i, &r) in rows.iter().enumerate() {
            let i = i as u32;
            match self.rows[r].kind {
                RowKind::Barycenter { component, axis } => {
                    let tot: f64 = (0..nv)
                        .filter(|&v| self.in_component(v, component))
                        .map(|v| va[v])
                        .sum();
                    for v in 0..nv {
                        if self.in_component(v, component) {
                            trip.push((i, (3 * v) as u32 + axis as u32, va[v] / tot));
                        }
                    }
                }
                RowKind::Pin { vertex, axis } => {
                    trip.push((i, 3 * vertex + axis as u32, 1.0));
                }
                RowKind::Area { component } => {
                    let g = self.area_gradient(positions, component);
                    push_vertex_row(&mut trip, i, &g);
                }
                RowKind::Volume { component } => {
                    let g = self.volume_gradient(positions, component);
                    push_vertex_row(&mut trip, i, &g);
                }
            }
        }
        SparseMatrix::from_triplets(rows.len(), 3 * nv, &trip).map_err(|e| {
            RsError::InvalidOption {
                name: "constraint Jacobian",
                value: e.to_string(),
                requirement: "a representable sparse matrix",
            }
        })
    }

    /// `∇A` of a component (or all): `Σ_{S∋v} ½ N_S × e_opp(v)`.
    ///
    /// [formula: rs-area-gradient]
    pub fn area_gradient(&self, positions: &[[f64; 3]], component: Option<u32>) -> Vec<[f64; 3]> {
        let mut g = vec![[0.0; 3]; self.n_vertices];
        for (t, fc) in self.faces.iter().zip(&self.face_component) {
            if component.is_none_or(|c| *fc == c) {
                let d = face_area_gradient(corners(positions, t));
                for c in 0..3 {
                    let v = t[c] as usize;
                    for k in 0..3 {
                        g[v][k] += d[c][k];
                    }
                }
            }
        }
        g
    }

    /// `∇V` of a component (or all): `⅙ Σ_{S∋v} f_j × f_k` (`(v, j, k)` the face's cyclic
    /// order). On a closed consistently oriented fan it equals `⅓ Σ a_S N_S`; the reference's
    /// row `Σ a_S N_S` is 3 times this (D9).
    ///
    /// [formula: rs-volume-gradient]
    pub fn volume_gradient(&self, positions: &[[f64; 3]], component: Option<u32>) -> Vec<[f64; 3]> {
        let mut g = vec![[0.0; 3]; self.n_vertices];
        for (t, fc) in self.faces.iter().zip(&self.face_component) {
            if component.is_none_or(|c| *fc == c) {
                let p = corners(positions, t);
                for c in 0..3 {
                    let d = cross(p[(c + 1) % 3], p[(c + 2) % 3]);
                    let v = t[c] as usize;
                    for k in 0..3 {
                        g[v][k] += d[k] / 6.0;
                    }
                }
            }
        }
        g
    }

    /// The exact projection onto the active positional rows (`SPEC.md` §5.1, §6.2): every
    /// active barycenter group translates its component (or the whole mesh) by `−Φ` (a
    /// translation keeps the areas, so the barycenter then equals its target exactly), then every
    /// active pin is reset to its target (pins last, so they hold exactly; a component with both
    /// pins and a barycenter keeps the pins and its barycenter is off by the pins' share).
    ///
    /// [formula: rs-barycenter]
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`].
    pub fn project_positional(&self, positions: &mut [[f64; 3]]) -> Result<(), RsError> {
        self.check_len(positions)?;
        let mut r = 0;
        while r < self.rows.len() {
            if let RowKind::Barycenter { component, axis: 0 } = self.rows[r].kind
                && self.active[r]
            {
                let (_, va) = self.areas(positions);
                let b = self.barycenter(&va, positions, component);
                let d = [
                    b[0] - self.rows[r].target,
                    b[1] - self.rows[r + 1].target,
                    b[2] - self.rows[r + 2].target,
                ];
                for v in 0..self.n_vertices {
                    if self.in_component(v, component) {
                        for k in 0..3 {
                            positions[v][k] -= d[k];
                        }
                    }
                }
                r += 3;
                continue;
            }
            r += 1;
        }
        for (row, &act) in self.rows.iter().zip(&self.active) {
            if let RowKind::Pin { vertex, axis } = row.kind
                && act
            {
                positions[vertex as usize][axis as usize] = row.target;
            }
        }
        Ok(())
    }

    /// The vertex component numbers.
    pub(crate) fn vertex_component(&self) -> &[u32] {
        &self.vertex_component
    }

    /// Number of components.
    pub fn n_components(&self) -> usize {
        self.n_components
    }

    pub(crate) fn vertex_areas(&self, positions: &[[f64; 3]]) -> Vec<f64> {
        self.areas(positions).1
    }
}

fn push_vertex_row(trip: &mut Vec<(u32, u32, f64)>, i: u32, g: &[[f64; 3]]) {
    for (v, gv) in g.iter().enumerate() {
        for k in 0..3 {
            if gv[k] != 0.0 {
                trip.push((i, (3 * v + k) as u32, gv[k]));
            }
        }
    }
}
