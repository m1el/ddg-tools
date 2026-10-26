//! The error type of the crate.
//!
//! Every variant names the element that is wrong (an edge, a vertex, a constraint row) and says
//! what is wrong with it. Bad input never panics. The enum is `#[non_exhaustive]`: later work
//! packages of `SPEC.md` §15 use the variants for their own checks (constraints, solves,
//! obstacles) and may be extended in the future.

use std::fmt;

/// All errors of the crate (`SPEC.md` §2.3, §9.1).
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum RcError {
    // ---- topology (CurveNetwork::new) ----
    /// An edge refers to a vertex index `>= n_vertices`.
    VertexOutOfRange {
        /// Edge index.
        edge: usize,
        /// The offending vertex index.
        vertex: u32,
        /// Number of vertices of the network.
        n_vertices: usize,
    },
    /// An edge joins a vertex with itself.
    SelfLoop {
        /// Edge index.
        edge: usize,
        /// The vertex.
        vertex: u32,
    },
    /// Two edges join the same unordered pair of vertices.
    DuplicateEdge {
        /// Smaller vertex index of the pair.
        a: u32,
        /// Larger vertex index of the pair.
        b: u32,
        /// The two edge indices.
        edges: [usize; 2],
    },
    /// A vertex belongs to no edge (no tangent, no mass).
    IsolatedVertex {
        /// The vertex.
        vertex: u32,
    },
    /// No two edges are non-adjacent (a single edge, a triangle, a star): the energy and the
    /// metric are identically zero.
    NoInteractingPairs {
        /// Number of edges of the network.
        n_edges: usize,
    },
    // ---- geometry (CurveNetwork::validate_positions) ----
    /// The position array has the wrong length.
    PositionCountMismatch {
        /// Number of vertices of the network.
        expected: usize,
        /// Length of the array that was passed.
        found: usize,
    },
    /// A coordinate is NaN or infinite.
    NonFinitePosition {
        /// The vertex.
        vertex: u32,
    },
    /// An edge has length zero (its endpoints are at the same point).
    ZeroLengthEdge {
        /// Edge index.
        edge: usize,
        /// First endpoint.
        a: u32,
        /// Second endpoint.
        b: u32,
    },
    /// Two distinct vertices are at the same point (the kernel is infinite there).
    CoincidentVertices {
        /// Smaller vertex index.
        a: u32,
        /// Larger vertex index.
        b: u32,
    },
    /// Two edges without a common vertex are at distance `<= 1e-12 * mean edge length`.
    EdgesIntersect {
        /// First edge index.
        a: usize,
        /// Second edge index.
        b: usize,
        /// Distance between the two segments.
        distance: f64,
    },
    /// Planar mode: a vertex does not have the common value of the held coordinate.
    NotPlanar {
        /// The vertex that differs from vertex 0.
        vertex: u32,
        /// The held axis (0, 1, 2 for X, Y, Z).
        axis: usize,
        /// The value of the coordinate at `vertex`.
        value: f64,
    },
    // ---- options ----
    /// `Exponents::new` rejected `(alpha, beta)`.
    InvalidExponents {
        /// The exponent `α` that was given.
        alpha: f64,
        /// The exponent `β` that was given.
        beta: f64,
        /// The requirement that is violated.
        requirement: &'static str,
    },
    /// An option has a value outside its range.
    InvalidOption {
        /// Name of the option.
        name: &'static str,
        /// The value, printed.
        value: String,
        /// The requirement that is violated.
        requirement: &'static str,
    },
    /// A feature that a later work package provides was asked for.
    Unsupported {
        /// What was asked for.
        what: &'static str,
    },
    // ---- constraints (WP2/WP3) ----
    /// The constraints do not remove the translations: `Ā` is singular on `ker C`.
    TranslationNotFixed,
    /// The constraint Jacobian is rank deficient.
    RedundantConstraints {
        /// Row indices (in constraint order) that are redundant.
        rows: Vec<usize>,
    },
    /// A constraint does not hold at the start of a session.
    InfeasibleConstraint {
        /// Index of the constraint in the list.
        constraint: usize,
        /// The violation, in the row's unit.
        violation: f64,
    },
    /// A tangent constraint at a vertex of degree other than 1 or 2.
    InvalidTangentSite {
        /// The vertex.
        vertex: u32,
        /// Its degree.
        degree: usize,
    },
    /// A constraint refers to an edge or vertex that does not exist.
    ConstraintOutOfRange {
        /// Index of the constraint in the list.
        constraint: usize,
        /// The offending index.
        index: u32,
    },
    /// A length schedule without a total-length or edge-length constraint.
    ScheduleWithoutLengthConstraint,
    /// L-BFGS or nonlinear CG asked for with hard constraints.
    MethodNeedsSoftProblem,
    // ---- solves (WP2/WP3/WP6) ----
    /// The metric is singular on the constraint's null space.
    SingularMetric {
        /// Number of connected components of the network.
        components: usize,
    },
    /// The saddle-point matrix has a zero pivot.
    SingularSaddle {
        /// Index of the pivot.
        pivot: usize,
    },
    /// The multigrid residual grew (or became non-finite): the cycles are not a contraction on
    /// this problem. The session falls back to the dense solver when the size permits
    /// (`SPEC.md` §17h).
    MultigridDiverged {
        /// Relative residuals per cycle (the first entry is the full-multigrid start).
        residuals: Vec<f64>,
    },
    /// The multigrid residual did not grow but stagnated above `failure_residual` when the
    /// cycle budget ended: the solve is too inaccurate to use, yet the cycles did not diverge.
    /// Handled like [`RcError::MultigridDiverged`] by the session (`SPEC.md` §17h).
    MultigridNotConverged {
        /// Relative residuals per cycle (the first entry is the full-multigrid start).
        residuals: Vec<f64>,
        /// The bound `MgOptions::failure_residual` the final residual stayed above.
        failure_residual: f64,
    },
    /// A step would shrink an edge below `FlowOptions::min_edge_length_ratio` times the mean edge
    /// length at the start of the session (`SPEC.md` §17j). The step is not taken: the curve is
    /// that of the last accepted step.
    EdgeCollapse {
        /// The edge.
        edge: usize,
        /// Its length after the rejected step.
        length: f64,
        /// Its length before.
        length_before: f64,
        /// The bound `ratio · ℓ̄₀`.
        threshold: f64,
        /// The step (1-based) that would have done it.
        step: usize,
    },
    /// The energy is NaN or infinite.
    NonFiniteEnergy {
        /// The step at which it happened.
        step: usize,
    },
    // ---- obstacles (WP7) ----
    /// An obstacle triangle refers to a vertex that does not exist.
    ObstacleTriangleOutOfRange {
        /// Triangle index.
        triangle: usize,
        /// The offending vertex index.
        vertex: u32,
        /// Number of obstacle vertices.
        n_vertices: usize,
    },
    /// An obstacle triangle has zero area.
    DegenerateObstacleTriangle {
        /// Triangle index.
        triangle: usize,
    },
    // ---- I/O ----
    /// A file could not be read or written.
    Io {
        /// The path.
        path: String,
        /// The operating system's message.
        message: String,
    },
    /// A syntax error in an input file.
    Parse {
        /// 1-based line number.
        line: usize,
        /// What is wrong.
        message: String,
    },
    /// The run was cancelled through the cancel flag.
    Cancelled,
}

impl fmt::Display for RcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use RcError::*;
        match self {
            VertexOutOfRange {
                edge,
                vertex,
                n_vertices,
            } => write!(
                f,
                "edge {edge} refers to vertex {vertex}, but the network has only {n_vertices} vertices"
            ),
            SelfLoop { edge, vertex } => {
                write!(
                    f,
                    "edge {edge} ({vertex}, {vertex}) is degenerate: it joins a vertex with itself"
                )
            }
            DuplicateEdge { a, b, edges } => write!(
                f,
                "edges {} and {} both join vertices {a} and {b}",
                edges[0], edges[1]
            ),
            IsolatedVertex { vertex } => {
                write!(
                    f,
                    "vertex {vertex} belongs to no edge (isolated vertices are not allowed)"
                )
            }
            NoInteractingPairs { n_edges } => write!(
                f,
                "no two of the {n_edges} edges are disjoint (a single edge, a triangle or a star): the energy and the metric are zero"
            ),
            PositionCountMismatch { expected, found } => {
                write!(f, "expected {expected} vertex positions, got {found}")
            }
            NonFinitePosition { vertex } => {
                write!(f, "vertex {vertex} has a non-finite coordinate")
            }
            ZeroLengthEdge { edge, a, b } => write!(
                f,
                "edge {edge} ({a}, {b}) has length zero: its endpoints coincide"
            ),
            CoincidentVertices { a, b } => {
                write!(f, "vertices {a} and {b} coincide (same position)")
            }
            EdgesIntersect { a, b, distance } => write!(
                f,
                "edges {a} and {b} intersect (distance {distance:e}) although they share no vertex"
            ),
            NotPlanar {
                vertex,
                axis,
                value,
            } => write!(
                f,
                "planar mode holds coordinate {} fixed, but vertex {vertex} has {value} there, unlike vertex 0",
                ["x", "y", "z"].get(*axis).copied().unwrap_or("?")
            ),
            InvalidExponents {
                alpha,
                beta,
                requirement,
            } => write!(
                f,
                "invalid exponents (alpha = {alpha}, beta = {beta}): {requirement}"
            ),
            InvalidOption {
                name,
                value,
                requirement,
            } => {
                write!(f, "option {name} = {value} is invalid: {requirement}")
            }
            Unsupported { what } => write!(f, "{what} is not available in this build"),
            TranslationNotFixed => write!(
                f,
                "the constraints do not fix the translations of the network (the metric is singular on their null space)"
            ),
            RedundantConstraints { rows } => {
                write!(
                    f,
                    "constraint rows {rows:?} are linearly dependent on the others"
                )
            }
            InfeasibleConstraint {
                constraint,
                violation,
            } => write!(
                f,
                "constraint {constraint} is violated by {violation:e} at the start"
            ),
            InvalidTangentSite { vertex, degree } => write!(
                f,
                "tangent constraint at vertex {vertex} of degree {degree}: only degrees 1 and 2 have a tangent"
            ),
            ConstraintOutOfRange { constraint, index } => write!(
                f,
                "constraint {constraint} refers to element {index}, which does not exist"
            ),
            ScheduleWithoutLengthConstraint => write!(
                f,
                "a length schedule needs a TotalLength or EdgeLengths constraint"
            ),
            MethodNeedsSoftProblem => write!(
                f,
                "L-BFGS and nonlinear CG need a problem without hard constraints"
            ),
            SingularMetric { components } => write!(
                f,
                "the metric is singular on the constraint null space ({components} components)"
            ),
            SingularSaddle { pivot } => {
                write!(
                    f,
                    "the saddle-point matrix has a zero pivot at index {pivot}"
                )
            }
            MultigridDiverged { residuals } => write!(
                f,
                "multigrid diverged: the relative residual grew or became non-finite \
                 (first {:e}, last {:e}); relative residuals {residuals:?}",
                residuals.first().copied().unwrap_or(f64::NAN),
                residuals.last().copied().unwrap_or(f64::NAN),
            ),
            MultigridNotConverged {
                residuals,
                failure_residual,
            } => write!(
                f,
                "multigrid did not converge: the relative residual stagnated at {:e} after {} \
                 residual evaluations, above the failure bound {failure_residual:e}; relative \
                 residuals {residuals:?}",
                residuals.last().copied().unwrap_or(f64::NAN),
                residuals.len(),
            ),
            EdgeCollapse {
                edge,
                length,
                length_before,
                threshold,
                step,
            } => write!(
                f,
                "step {step} would shrink edge {edge} from {length_before:e} to {length:e}, below the \
                 bound {threshold:e} (min_edge_length_ratio times the initial mean edge length): \
                 the edge is collapsing, typically an end edge with a fixed tangent and a free \
                 length (SPEC section 17j); fix the edge lengths (EdgeLengths constraint), or set \
                 min_edge_length_ratio to None to allow it"
            ),
            NonFiniteEnergy { step } => write!(f, "the energy is not finite at step {step}"),
            ObstacleTriangleOutOfRange {
                triangle,
                vertex,
                n_vertices,
            } => write!(
                f,
                "obstacle triangle {triangle} refers to vertex {vertex}, but the obstacle has only {n_vertices} vertices"
            ),
            DegenerateObstacleTriangle { triangle } => {
                write!(f, "obstacle triangle {triangle} has zero area")
            }
            Io { path, message } => write!(f, "{path}: {message}"),
            Parse { line, message } => write!(f, "line {line}: {message}"),
            Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for RcError {}
