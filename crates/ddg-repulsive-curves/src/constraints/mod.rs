//! Hard constraints `Φ(γ) = 0` of `SPEC.md` §5.1 and their rank checks (§5.2).
//!
//! A [`ConstraintSet`] is built once per session (or after a topology change) from the
//! caller's [`Constraint`] list: targets given as `None` take the values at the positions passed
//! to [`ConstraintSet::new`] ("keep"), indices and tangent sites are validated, and the known
//! redundancies of D9 are resolved (rows dropped and listed in
//! [`ConstraintSet::dropped_rows`]). It then evaluates `Φ` ([`ConstraintSet::values`]) and the
//! Jacobian `C = dΦ` ([`ConstraintSet::jacobian`]) at any positions, with rows in the order of
//! the list and columns in the interleaved layout of [`DofLayout`] (D29).
//!
//! Row units and tolerances (§5.1): length rows use `τ_len`, dimensionless (tangent) rows
//! `τ_len / ℓ̄₀`, surface rows `τ_len · |∇f(γ_v)|` ([`ConstraintSet::row_tolerances`]).
//!
//! The rank checks ([`ConstraintSet::analyze`], [`rank`]) give
//! [`RcError::TranslationNotFixed`] and [`RcError::RedundantConstraints`].

use std::fmt;
use std::sync::Arc;

use ddg_core::sparse::SparseMatrix;

use crate::error::RcError;
use crate::kernel::v3::{dot, norm, scale};
use crate::metric::DofLayout;
use crate::network::{Axis, CurveNetwork};
use crate::solve::{Tolerance, sparse_error};

pub mod jacobian;
pub mod rank;

pub use jacobian::tangent_frame;
pub use rank::{RankReport, check_component_translations, check_translations};

/// An implicit surface `f(x) = 0` for [`Constraint::Surface`] (`SPEC.md` §5.1, D8). `f` need
/// not be a signed distance. Built-in surfaces (sphere, torus, plane, union) implement this
/// trait in [`crate::surfaces`].
pub trait ImplicitSurface: Send + Sync {
    /// `f(x)`.
    fn value(&self, x: [f64; 3]) -> f64;
    /// `∇f(x)`.
    fn gradient(&self, x: [f64; 3]) -> [f64; 3];
}

/// Where a tangent constraint applies (D7).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TangentSite {
    /// The tangent `T_I` of edge `I`.
    Edge(u32),
    /// The tangent at a vertex of degree 1 (the edge oriented away from it) or 2 (the sum of
    /// the two unit edge vectors through it).
    Vertex(u32),
}

/// A hard constraint (`SPEC.md` §5.1, §9.1). `None` targets are the values at the positions
/// the [`ConstraintSet`] is built from.
#[derive(Clone)]
pub enum Constraint {
    /// `(Σ_I ℓ_I x_I)/(Σ_I ℓ_I) − x₀` (D6), 3 rows (2 in planar mode).
    Barycenter {
        /// `x₀`.
        target: Option<[f64; 3]>,
    },
    /// `L⁰ − Σ_I ℓ_I`, 1 row.
    TotalLength {
        /// `L⁰`.
        target: Option<f64>,
    },
    /// `ℓ⁰_I − ℓ_I`, one row per listed edge.
    EdgeLengths {
        /// The edges (`None`: all, in index order).
        edges: Option<Vec<u32>>,
        /// One target per listed edge (`None`: current lengths).
        targets: Option<Vec<f64>>,
    },
    /// `γ_v − x_v`, 3 rows (2 in planar mode).
    Point {
        /// The vertex.
        vertex: u32,
        /// `x_v`.
        target: Option<[f64; 3]>,
    },
    /// `f(γ_v)` for each listed vertex (D8), one row each.
    Surface {
        /// The vertices (`None`: all, in index order).
        vertices: Option<Vec<u32>>,
        /// The surface.
        surface: Arc<dyn ImplicitSurface>,
    },
    /// `(e₁·s, e₂·s)` with `e₁, e₂` completing `X/|X|` to an orthonormal frame (D7); 2 rows
    /// (1 in planar mode), plus the orientation requirement `s·X > 0`.
    Tangent {
        /// The site.
        site: TangentSite,
        /// `X` (`None`: the current `s`).
        direction: Option<[f64; 3]>,
    },
}

impl fmt::Debug for Constraint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Constraint::Barycenter { target } => f
                .debug_struct("Barycenter")
                .field("target", target)
                .finish(),
            Constraint::TotalLength { target } => f
                .debug_struct("TotalLength")
                .field("target", target)
                .finish(),
            Constraint::EdgeLengths { edges, targets } => f
                .debug_struct("EdgeLengths")
                .field("edges", edges)
                .field("targets", targets)
                .finish(),
            Constraint::Point { vertex, target } => f
                .debug_struct("Point")
                .field("vertex", vertex)
                .field("target", target)
                .finish(),
            Constraint::Surface { vertices, .. } => f
                .debug_struct("Surface")
                .field("vertices", vertices)
                .finish_non_exhaustive(),
            Constraint::Tangent { site, direction } => f
                .debug_struct("Tangent")
                .field("site", site)
                .field("direction", direction)
                .finish(),
        }
    }
}

/// The default constraint list of `SPEC.md` §10: `[Barycenter, EdgeLengths(all, current)]`.
pub fn default_constraints() -> Vec<Constraint> {
    vec![
        Constraint::Barycenter { target: None },
        Constraint::EdgeLengths {
            edges: None,
            targets: None,
        },
    ]
}

/// A constraint with every target and index resolved.
#[derive(Clone)]
pub enum ResolvedConstraint {
    /// See [`Constraint::Barycenter`].
    Barycenter {
        /// `x₀`.
        target: [f64; 3],
    },
    /// See [`Constraint::TotalLength`].
    TotalLength {
        /// `L⁰`.
        target: f64,
    },
    /// See [`Constraint::EdgeLengths`].
    EdgeLengths {
        /// The edges.
        edges: Vec<u32>,
        /// One target per edge.
        targets: Vec<f64>,
    },
    /// See [`Constraint::Point`].
    Point {
        /// The vertex.
        vertex: u32,
        /// `x_v`.
        target: [f64; 3],
    },
    /// See [`Constraint::Surface`].
    Surface {
        /// The vertices.
        vertices: Vec<u32>,
        /// The surface.
        surface: Arc<dyn ImplicitSurface>,
    },
    /// See [`Constraint::Tangent`].
    Tangent {
        /// The site.
        site: TangentSite,
        /// `X` as given (not normalized).
        direction: [f64; 3],
        /// The rows' unit vectors `e_k ⊥ X` (one in planar mode, then `frame[1]` is unused).
        frame: [[f64; 3]; 2],
        /// Number of rows (2, or 1 in planar mode).
        n_rows: usize,
    },
}

impl fmt::Debug for ResolvedConstraint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ResolvedConstraint::Barycenter { target } => f
                .debug_struct("Barycenter")
                .field("target", target)
                .finish(),
            ResolvedConstraint::TotalLength { target } => f
                .debug_struct("TotalLength")
                .field("target", target)
                .finish(),
            ResolvedConstraint::EdgeLengths { edges, targets } => f
                .debug_struct("EdgeLengths")
                .field("edges", edges)
                .field("targets", targets)
                .finish(),
            ResolvedConstraint::Point { vertex, target } => f
                .debug_struct("Point")
                .field("vertex", vertex)
                .field("target", target)
                .finish(),
            ResolvedConstraint::Surface { vertices, .. } => f
                .debug_struct("Surface")
                .field("vertices", vertices)
                .finish_non_exhaustive(),
            ResolvedConstraint::Tangent {
                site,
                direction,
                frame,
                n_rows,
            } => f
                .debug_struct("Tangent")
                .field("site", site)
                .field("direction", direction)
                .field("frame", frame)
                .field("n_rows", n_rows)
                .finish(),
        }
    }
}

/// The unit of a constraint row, which sets its projection tolerance (§5.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowUnit {
    /// A length: tolerance `τ_len`.
    Length,
    /// Dimensionless (tangent rows): tolerance `τ_len / ℓ̄₀`.
    Dimensionless,
    /// The unit of `f` (surface rows): tolerance `τ_len · |∇f(γ_v)|`.
    Surface,
}

/// What a constraint row is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowKind {
    /// Coordinate `coord` of the barycenter.
    Barycenter {
        /// 0, 1, 2 for x, y, z.
        coord: usize,
    },
    /// The total length.
    TotalLength,
    /// The length of one edge.
    EdgeLength {
        /// The edge.
        edge: u32,
        /// Position of the edge in the constraint's list.
        slot: usize,
    },
    /// Coordinate `coord` of a pinned vertex.
    Point {
        /// The vertex.
        vertex: u32,
        /// 0, 1, 2 for x, y, z.
        coord: usize,
    },
    /// `f` at one vertex.
    Surface {
        /// The vertex.
        vertex: u32,
        /// Position of the vertex in the constraint's list.
        slot: usize,
    },
    /// Row `k` of a tangent constraint.
    Tangent {
        /// 0 or 1.
        k: usize,
    },
}

impl RowKind {
    /// The row's unit.
    pub fn unit(&self) -> RowUnit {
        match self {
            RowKind::Tangent { .. } => RowUnit::Dimensionless,
            RowKind::Surface { .. } => RowUnit::Surface,
            _ => RowUnit::Length,
        }
    }
}

/// One constraint row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RowInfo {
    /// Index of the constraint in the list.
    pub constraint: usize,
    /// What the row is.
    pub kind: RowKind,
}

/// The resolved constraints of a session, their rows and the rows dropped by D9.
#[derive(Debug, Clone)]
pub struct ConstraintSet {
    n_vertices: usize,
    layout: DofLayout,
    resolved: Vec<ResolvedConstraint>,
    rows: Vec<RowInfo>,
    /// Indices into `rows` of the active rows, increasing.
    active: Vec<usize>,
    /// Indices into `rows` of the dropped rows, increasing.
    dropped: Vec<usize>,
}

fn finite3(x: [f64; 3]) -> bool {
    x.iter().all(|c| c.is_finite())
}

fn bad_option(name: &'static str, value: String, requirement: &'static str) -> RcError {
    RcError::InvalidOption {
        name,
        value,
        requirement,
    }
}

impl ConstraintSet {
    /// Resolves `constraints` at `positions` (`SPEC.md` §5.1): `None` targets become the current
    /// values, indices are checked, tangent frames are built (and `s·X > 0` is required), and
    /// the redundancies of D9 are dropped (an `EdgeLengths` row whose two endpoints are both
    /// `Point`-constrained; the `Barycenter` rows when every vertex is `Point`-constrained).
    /// Feasibility of the other targets is checked separately
    /// ([`ConstraintSet::check_feasible`]), the rank by [`ConstraintSet::analyze`].
    ///
    /// # Errors
    /// [`RcError::PositionCountMismatch`]; [`RcError::ConstraintOutOfRange`] for an edge or
    /// vertex that does not exist; [`RcError::InvalidTangentSite`] for a vertex site of degree
    /// other than 1 or 2; [`RcError::InvalidOption`] for non-finite targets, a target list of
    /// the wrong length, non-positive length targets, a zero tangent direction, or (planar
    /// mode) a tangent direction or point target off the plane; [`RcError::InfeasibleConstraint`]
    /// for a tangent site with `s·X ≤ 0` (violation `1 − cos∠(s, X)`).
    pub fn new(
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        constraints: &[Constraint],
        planar: Option<Axis>,
    ) -> Result<ConstraintSet, RcError> {
        let n = network.n_vertices();
        let m = network.n_edges();
        if positions.len() != n {
            return Err(RcError::PositionCountMismatch {
                expected: n,
                found: positions.len(),
            });
        }
        let layout = DofLayout::new(planar);
        let geom = network.edge_geometry(positions);
        let mut resolved = Vec::with_capacity(constraints.len());
        for (ci, c) in constraints.iter().enumerate() {
            let out_of_range = |index: u32| RcError::ConstraintOutOfRange {
                constraint: ci,
                index,
            };
            let r = match c {
                Constraint::Barycenter { target } => {
                    let target = match target {
                        Some(t) => {
                            if !finite3(*t) {
                                return Err(bad_option(
                                    "Barycenter.target",
                                    format!("{t:?}"),
                                    "finite coordinates",
                                ));
                            }
                            *t
                        }
                        None => jacobian::barycenter(network, &geom).0,
                    };
                    ResolvedConstraint::Barycenter { target }
                }
                Constraint::TotalLength { target } => {
                    let target = match target {
                        Some(t) => {
                            if !(t.is_finite() && *t > 0.0) {
                                return Err(bad_option(
                                    "TotalLength.target",
                                    t.to_string(),
                                    "a finite positive length",
                                ));
                            }
                            *t
                        }
                        None => geom.length.iter().sum(),
                    };
                    ResolvedConstraint::TotalLength { target }
                }
                Constraint::EdgeLengths { edges, targets } => {
                    let edges: Vec<u32> = match edges {
                        Some(e) => e.clone(),
                        None => (0..m as u32).collect(),
                    };
                    if let Some(&bad) = edges.iter().find(|&&e| e as usize >= m) {
                        return Err(out_of_range(bad));
                    }
                    let targets = match targets {
                        Some(t) => {
                            if t.len() != edges.len() {
                                return Err(bad_option(
                                    "EdgeLengths.targets",
                                    format!("{} targets for {} edges", t.len(), edges.len()),
                                    "one target per listed edge",
                                ));
                            }
                            if let Some(x) = t.iter().find(|x| !(x.is_finite() && **x > 0.0)) {
                                return Err(bad_option(
                                    "EdgeLengths.targets",
                                    x.to_string(),
                                    "finite positive lengths",
                                ));
                            }
                            t.clone()
                        }
                        None => edges.iter().map(|&e| geom.length[e as usize]).collect(),
                    };
                    ResolvedConstraint::EdgeLengths { edges, targets }
                }
                Constraint::Point { vertex, target } => {
                    if *vertex as usize >= n {
                        return Err(out_of_range(*vertex));
                    }
                    let target = match target {
                        Some(t) => {
                            if !finite3(*t) {
                                return Err(bad_option(
                                    "Point.target",
                                    format!("{t:?}"),
                                    "finite coordinates",
                                ));
                            }
                            if let Some(ax) = planar {
                                let k = ax.index();
                                if t[k] != positions[*vertex as usize][k] {
                                    return Err(bad_option(
                                        "Point.target",
                                        format!("{t:?}"),
                                        "in planar mode a point target must lie in the plane",
                                    ));
                                }
                            }
                            *t
                        }
                        None => positions[*vertex as usize],
                    };
                    ResolvedConstraint::Point {
                        vertex: *vertex,
                        target,
                    }
                }
                Constraint::Surface { vertices, surface } => {
                    let vertices: Vec<u32> = match vertices {
                        Some(v) => v.clone(),
                        None => (0..n as u32).collect(),
                    };
                    if let Some(&bad) = vertices.iter().find(|&&v| v as usize >= n) {
                        return Err(out_of_range(bad));
                    }
                    ResolvedConstraint::Surface {
                        vertices,
                        surface: surface.clone(),
                    }
                }
                Constraint::Tangent { site, direction } => {
                    match *site {
                        TangentSite::Edge(e) => {
                            if e as usize >= m {
                                return Err(out_of_range(e));
                            }
                        }
                        TangentSite::Vertex(v) => {
                            if v as usize >= n {
                                return Err(out_of_range(v));
                            }
                            let degree = network.degree(v);
                            if !(1..=2).contains(&degree) {
                                return Err(RcError::InvalidTangentSite { vertex: v, degree });
                            }
                        }
                    }
                    let s = jacobian::tangent_vector(network, positions, *site).0;
                    let x = match direction {
                        Some(x) => {
                            if !finite3(*x) || norm(*x) == 0.0 {
                                return Err(bad_option(
                                    "Tangent.direction",
                                    format!("{x:?}"),
                                    "a finite non-zero vector",
                                ));
                            }
                            *x
                        }
                        None => s,
                    };
                    let (frame, n_rows) = match planar {
                        None => (tangent_frame(x), 2),
                        Some(ax) => {
                            let k = ax.index();
                            if x[k].abs() > 1e-12 * norm(x) {
                                return Err(bad_option(
                                    "Tangent.direction",
                                    format!("{x:?}"),
                                    "in planar mode a tangent direction must lie in the plane",
                                ));
                            }
                            let mut axis = [0.0; 3];
                            axis[k] = 1.0;
                            let e = crate::kernel::v3::cross(axis, x);
                            ([scale(e, 1.0 / norm(e)), [0.0; 3]], 1)
                        }
                    };
                    let cosine = dot(s, x) / (norm(s) * norm(x));
                    if cosine.is_nan() || cosine <= 0.0 {
                        return Err(RcError::InfeasibleConstraint {
                            constraint: ci,
                            violation: 1.0 - if cosine.is_finite() { cosine } else { 0.0 },
                        });
                    }
                    ResolvedConstraint::Tangent {
                        site: *site,
                        direction: x,
                        frame,
                        n_rows,
                    }
                }
            };
            resolved.push(r);
        }
        let mut set = ConstraintSet {
            n_vertices: n,
            layout,
            resolved,
            rows: Vec::new(),
            active: Vec::new(),
            dropped: Vec::new(),
        };
        set.build_rows(network);
        Ok(set)
    }

    /// Row list and the D9 drops.
    fn build_rows(&mut self, network: &CurveNetwork) {
        let coords = self.layout.coords();
        let mut pinned = vec![false; self.n_vertices];
        for c in &self.resolved {
            if let ResolvedConstraint::Point { vertex, .. } = c {
                pinned[*vertex as usize] = true;
            }
        }
        let all_pinned = pinned.iter().all(|&p| p);
        let edges = network.edges();
        self.rows.clear();
        self.active.clear();
        self.dropped.clear();
        for (ci, c) in self.resolved.iter().enumerate() {
            let mut push = |kind: RowKind, drop: bool| {
                let r = self.rows.len();
                self.rows.push(RowInfo {
                    constraint: ci,
                    kind,
                });
                if drop {
                    self.dropped.push(r);
                } else {
                    self.active.push(r);
                }
            };
            match c {
                ResolvedConstraint::Barycenter { .. } => {
                    for &coord in coords {
                        push(RowKind::Barycenter { coord }, all_pinned);
                    }
                }
                ResolvedConstraint::TotalLength { .. } => push(RowKind::TotalLength, false),
                ResolvedConstraint::EdgeLengths { edges: list, .. } => {
                    for (slot, &e) in list.iter().enumerate() {
                        let [a, b] = edges[e as usize];
                        let drop = pinned[a as usize] && pinned[b as usize];
                        push(RowKind::EdgeLength { edge: e, slot }, drop);
                    }
                }
                ResolvedConstraint::Point { vertex, .. } => {
                    for &coord in coords {
                        push(
                            RowKind::Point {
                                vertex: *vertex,
                                coord,
                            },
                            false,
                        );
                    }
                }
                ResolvedConstraint::Surface { vertices, .. } => {
                    for (slot, &v) in vertices.iter().enumerate() {
                        push(RowKind::Surface { vertex: v, slot }, false);
                    }
                }
                ResolvedConstraint::Tangent { n_rows, .. } => {
                    for k in 0..*n_rows {
                        push(RowKind::Tangent { k }, false);
                    }
                }
            }
        }
    }

    /// The coordinate layout of the Jacobian's columns.
    pub fn layout(&self) -> DofLayout {
        self.layout
    }

    /// Number of vertices of the network the set was built for.
    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }

    /// The resolved constraints, in list order.
    pub fn constraints(&self) -> &[ResolvedConstraint] {
        &self.resolved
    }

    /// All rows (active and dropped), in constraint order.
    pub fn rows(&self) -> &[RowInfo] {
        &self.rows
    }

    /// Indices (into [`ConstraintSet::rows`]) of the active rows, increasing. Row `k` of
    /// [`ConstraintSet::values`] and [`ConstraintSet::jacobian`] is row `active_rows()[k]`.
    pub fn active_rows(&self) -> &[usize] {
        &self.active
    }

    /// Indices of the rows dropped by D9 (redundant and fixed by `Point` constraints).
    pub fn dropped_rows(&self) -> &[usize] {
        &self.dropped
    }

    /// Number of active rows `k`.
    pub fn n_active_rows(&self) -> usize {
        self.active.len()
    }

    fn check_vertices(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<(), RcError> {
        if network.n_vertices() != self.n_vertices || positions.len() != self.n_vertices {
            return Err(RcError::PositionCountMismatch {
                expected: self.n_vertices,
                found: positions.len().min(network.n_vertices()),
            });
        }
        Ok(())
    }

    /// `Φ(γ)` on all rows (dropped ones included) and, if asked, the Jacobian entries
    /// `(row, vertex, coordinate, value)` in 3-D coordinates.
    fn evaluate_all(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        with_jacobian: bool,
    ) -> (Vec<f64>, jacobian::JacobianEntries) {
        jacobian::evaluate(self, network, positions, with_jacobian)
    }

    /// `Φ(γ)` on the active rows.
    ///
    /// # Errors
    /// [`RcError::PositionCountMismatch`] if `network` or `positions` do not match the set.
    pub fn values(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<Vec<f64>, RcError> {
        self.check_vertices(network, positions)?;
        let (all, _) = self.evaluate_all(network, positions, false);
        Ok(self.active.iter().map(|&r| all[r]).collect())
    }

    /// The Jacobian `C = dΦ(γ)` on the active rows, `k × d|V|` (columns interleaved, held
    /// coordinate dropped in planar mode).
    ///
    /// # Errors
    /// [`RcError::PositionCountMismatch`]; resource errors of the sparse matrix.
    pub fn jacobian(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<SparseMatrix, RcError> {
        self.check_vertices(network, positions)?;
        let (_, entries) = self.evaluate_all(network, positions, true);
        let mut new_row = vec![usize::MAX; self.rows.len()];
        for (k, &r) in self.active.iter().enumerate() {
            new_row[r] = k;
        }
        let mut t = Vec::with_capacity(entries.len());
        for (r, v, c, val) in entries {
            let k = new_row[r];
            if k == usize::MAX {
                continue;
            }
            if let Some(slot) = self.layout.slot(c) {
                t.push((k as u32, self.layout.index(v as usize, slot) as u32, val));
            }
        }
        SparseMatrix::from_triplets(self.active.len(), self.layout.n_dofs(self.n_vertices), &t)
            .map_err(|e| sparse_error(e, "constraint Jacobian"))
    }

    /// The length tolerance `τ_len` of `tolerance` for the mean edge length `ℓ̄₀` at session
    /// start.
    pub fn length_tolerance(tolerance: Tolerance, mean_edge_length0: f64) -> f64 {
        match tolerance {
            Tolerance::Absolute(t) => t,
            Tolerance::MeanEdgeLength(t) => t * mean_edge_length0,
        }
    }

    /// Per active row, the tolerance of the projection (§5.1): `τ_len` for length rows,
    /// `τ_len / ℓ̄₀` for tangent rows, `τ_len · |∇f(γ_v)|` for surface rows (gradient over the
    /// active coordinates, at `positions`).
    ///
    /// # Errors
    /// [`RcError::PositionCountMismatch`].
    pub fn row_tolerances(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        tolerance: Tolerance,
        mean_edge_length0: f64,
    ) -> Result<Vec<f64>, RcError> {
        self.check_vertices(network, positions)?;
        let tau = Self::length_tolerance(tolerance, mean_edge_length0);
        Ok(self
            .active
            .iter()
            .map(|&r| {
                let row = self.rows[r];
                match row.kind {
                    RowKind::Surface { vertex, .. } => {
                        let ResolvedConstraint::Surface { surface, .. } =
                            &self.resolved[row.constraint]
                        else {
                            return tau;
                        };
                        let g = surface.gradient(positions[vertex as usize]);
                        let g2: f64 = self.layout.coords().iter().map(|&c| g[c] * g[c]).sum();
                        tau * g2.sqrt()
                    }
                    RowKind::Tangent { .. } => tau / mean_edge_length0,
                    _ => tau,
                }
            })
            .collect())
    }

    /// The first tangent constraint whose site has `s·X ≤ 0` at `positions`, with `s·X`.
    pub fn tangent_orientation_failure(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Option<(usize, f64)> {
        for (ci, c) in self.resolved.iter().enumerate() {
            if let ResolvedConstraint::Tangent {
                site, direction, ..
            } = c
            {
                let s = jacobian::tangent_vector(network, positions, *site).0;
                let d = dot(s, *direction);
                if d.is_nan() || d <= 0.0 {
                    return Some((ci, d));
                }
            }
        }
        None
    }

    /// Checks that every constraint holds at `positions` within the projection tolerance
    /// (all constraints must hold at the start, §5.1), that the rows dropped by D9 are
    /// consistent with the `Point` targets (an `EdgeLengths` target equals the distance of
    /// the two pinned targets within `τ_len`), and the tangent orientation.
    ///
    /// # Errors
    /// [`RcError::InfeasibleConstraint`] for the first violated row (its constraint index and
    /// the violation); [`RcError::PositionCountMismatch`].
    pub fn check_feasible(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        tolerance: Tolerance,
        mean_edge_length0: f64,
    ) -> Result<(), RcError> {
        self.check_vertices(network, positions)?;
        let tol = self.row_tolerances(network, positions, tolerance, mean_edge_length0)?;
        let (all, _) = self.evaluate_all(network, positions, false);
        for (k, &r) in self.active.iter().enumerate() {
            if all[r].is_nan() || all[r].abs() > tol[k] {
                return Err(RcError::InfeasibleConstraint {
                    constraint: self.rows[r].constraint,
                    violation: all[r],
                });
            }
        }
        let tau = Self::length_tolerance(tolerance, mean_edge_length0);
        let mut pin_target: Vec<Option<[f64; 3]>> = vec![None; self.n_vertices];
        for c in &self.resolved {
            if let ResolvedConstraint::Point { vertex, target } = c {
                pin_target[*vertex as usize] = Some(*target);
            }
        }
        for &r in &self.dropped {
            let row = self.rows[r];
            if let RowKind::EdgeLength { edge, slot } = row.kind {
                let ResolvedConstraint::EdgeLengths { targets, .. } =
                    &self.resolved[row.constraint]
                else {
                    continue;
                };
                let [a, b] = network.edges()[edge as usize];
                if let (Some(pa), Some(pb)) = (pin_target[a as usize], pin_target[b as usize]) {
                    let violation = targets[slot] - norm(crate::kernel::v3::sub(pb, pa));
                    if violation.is_nan() || violation.abs() > tau {
                        return Err(RcError::InfeasibleConstraint {
                            constraint: row.constraint,
                            violation,
                        });
                    }
                }
            }
        }
        if let Some((ci, d)) = self.tangent_orientation_failure(network, positions) {
            return Err(RcError::InfeasibleConstraint {
                constraint: ci,
                violation: d,
            });
        }
        Ok(())
    }

    fn constraint_mut(&mut self, constraint: usize) -> Result<&mut ResolvedConstraint, RcError> {
        let n = self.resolved.len();
        self.resolved
            .get_mut(constraint)
            .ok_or(RcError::ConstraintOutOfRange {
                constraint,
                index: n as u32,
            })
    }

    /// Sets the target `L⁰` of a `TotalLength` constraint (length schedules, §5.4).
    ///
    /// # Errors
    /// [`RcError::ConstraintOutOfRange`] for a bad index; [`RcError::InvalidOption`] if the
    /// constraint is not `TotalLength` or the target is not finite and positive.
    pub fn set_total_length_target(
        &mut self,
        constraint: usize,
        target: f64,
    ) -> Result<(), RcError> {
        if !(target.is_finite() && target > 0.0) {
            return Err(bad_option(
                "TotalLength.target",
                target.to_string(),
                "a finite positive length",
            ));
        }
        match self.constraint_mut(constraint)? {
            ResolvedConstraint::TotalLength { target: t } => {
                *t = target;
                Ok(())
            }
            _ => Err(bad_option(
                "constraint",
                constraint.to_string(),
                "the index of a TotalLength constraint",
            )),
        }
    }

    /// Sets the targets of an `EdgeLengths` constraint (one per listed edge).
    ///
    /// # Errors
    /// As [`ConstraintSet::set_total_length_target`], plus a wrong number of targets.
    pub fn set_edge_length_targets(
        &mut self,
        constraint: usize,
        targets: &[f64],
    ) -> Result<(), RcError> {
        if let Some(x) = targets.iter().find(|x| !(x.is_finite() && **x > 0.0)) {
            return Err(bad_option(
                "EdgeLengths.targets",
                x.to_string(),
                "finite positive lengths",
            ));
        }
        match self.constraint_mut(constraint)? {
            ResolvedConstraint::EdgeLengths { targets: t, .. } => {
                if t.len() != targets.len() {
                    return Err(bad_option(
                        "EdgeLengths.targets",
                        format!("{} targets for {} edges", targets.len(), t.len()),
                        "one target per listed edge",
                    ));
                }
                t.copy_from_slice(targets);
                Ok(())
            }
            _ => Err(bad_option(
                "constraint",
                constraint.to_string(),
                "the index of an EdgeLengths constraint",
            )),
        }
    }
}
