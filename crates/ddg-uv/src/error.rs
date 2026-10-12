//! The error type of the pipeline.
//!
//! Every message is written in the **caller's** numbering: vertex indices are indices into the
//! caller's position array (before loose vertices were dropped), and "face" means a polygon of
//! the caller, never a triangle of the triangulation.

use std::fmt;

/// The three families of [`UvError`], matching the three exceptions of the Python binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorKind {
    /// An argument is malformed: wrong length, index out of range, a number that is not finite,
    /// contradictory steering, an option outside its range.
    Input,
    /// The arguments are well formed but the mesh cannot be unwrapped as it is: non-manifold,
    /// not orientable, degenerate under a policy that refuses, a closed surface under a policy
    /// that refuses, steering that leaves no feasible cone configuration.
    Mesh,
    /// The mesh was accepted but no valid map came out: the cone solver or the flattener
    /// failed, or the pipeline's own final checks did.
    Solve,
}

/// Why a mesh could not be unwrapped.
#[derive(Debug, Clone, PartialEq)]
pub enum UvError {
    // ----- input -----
    /// The mesh has no polygons.
    NoPolygons,
    /// More vertices or polygon corners than 32-bit ids can address.
    TooLarge {
        /// Which array.
        what: &'static str,
        /// Its number of entries.
        count: usize,
    },
    /// `polygon_offsets` is not a list of `n_polygons + 1` increasing corner offsets that
    /// starts at 0 and ends at the number of corners.
    PolygonOffsets {
        /// Position in `polygon_offsets` where the rule is broken.
        position: usize,
        /// The value found there.
        value: u32,
        /// What it should have been.
        requirement: String,
    },
    /// A polygon has fewer than three corners.
    PolygonTooSmall {
        /// The polygon.
        polygon: u32,
        /// Its number of corners.
        n_corners: usize,
    },
    /// An index is out of range.
    IndexOutOfRange {
        /// Which array.
        what: &'static str,
        /// Position in that array (for arrays of pairs or triples: the row).
        position: usize,
        /// The offending value.
        value: u32,
        /// Number of valid targets.
        count: usize,
        /// What the index refers to ("vertices", "polygon corners").
        of: &'static str,
    },
    /// A number is NaN or infinite.
    NonFinite {
        /// Which array.
        what: &'static str,
        /// Position in that array (for positions: the vertex).
        position: usize,
        /// The value.
        value: f64,
    },
    /// The number of triangles is not `corners − 2 · polygons`.
    TriangleCount {
        /// Number of triangles given.
        got: usize,
        /// Number required.
        expected: usize,
    },
    /// The triangles of a polygon do not triangulate it.
    BadTriangulation {
        /// The polygon.
        polygon: u32,
        /// What is wrong.
        problem: String,
    },
    /// A seam edge joins a vertex to itself.
    DegenerateSeamEdge {
        /// Row in `seam_edges`.
        position: usize,
        /// The vertex.
        vertex: u32,
    },
    /// A vertex is prescribed twice.
    DuplicatePrescribed {
        /// The vertex.
        vertex: u32,
    },
    /// A prescribed cone angle is `2π` or more: no such cone exists.
    PrescribedAngleTooLarge {
        /// The vertex.
        vertex: u32,
        /// The angle.
        angle: f64,
    },
    /// A vertex is in more than one of the pinned, prescribed and forbidden sets.
    ConflictingSteering {
        /// The vertex.
        vertex: u32,
        /// One set it is in.
        first: &'static str,
        /// Another set it is in.
        second: &'static str,
    },
    /// An option is outside its range.
    InvalidOption {
        /// The option.
        name: &'static str,
        /// Its value.
        value: String,
        /// What is required.
        requirement: &'static str,
    },

    // ----- mesh -----
    /// A polygon uses a vertex more than once.
    RepeatedVertexInPolygon {
        /// The polygon.
        polygon: u32,
        /// The vertex.
        vertex: u32,
    },
    /// Two polygons cover the same triangle.
    DuplicatePolygons {
        /// The lower-numbered polygon.
        first: u32,
        /// The higher-numbered polygon (equal to `first` if one polygon covers a triangle
        /// twice).
        second: u32,
        /// The three vertices of the triangle.
        vertices: [u32; 3],
    },
    /// An edge is shared by more than two polygons.
    NonManifoldEdge {
        /// The edge as `[min, max]`.
        edge: [u32; 2],
        /// The polygons at the edge, ascending.
        polygons: Vec<u32>,
    },
    /// A diagonal of the given triangulation coincides with another edge of the mesh, so the
    /// triangulated surface is not manifold there.
    DiagonalOnEdge {
        /// The edge as `[min, max]`.
        edge: [u32; 2],
        /// The polygons concerned, ascending.
        polygons: Vec<u32>,
        /// Number of triangles at the edge.
        n_triangles: usize,
    },
    /// The polygons around a vertex form more than one fan.
    NonManifoldVertex {
        /// The vertex.
        vertex: u32,
        /// Number of separate fans.
        fans: u32,
    },
    /// A connected component has no consistent orientation (a Möbius strip, a Klein bottle).
    NonOrientable {
        /// The edge at which the contradiction was found, as `[min, max]`. It can be a
        /// diagonal of the triangulation.
        edge: [u32; 2],
        /// The two polygons that cannot both be oriented consistently, ascending. They are
        /// the same polygon if the contradiction was found across one of its diagonals.
        polygons: [u32; 2],
    },
    /// Two polygons traverse their shared edge in the same direction, and the orientation
    /// policy is to refuse.
    InconsistentOrientation {
        /// The edge as `[min, max]`.
        edge: [u32; 2],
        /// The two polygons, ascending.
        polygons: [u32; 2],
    },
    /// The triangulation has triangles without area, and the policy is to refuse (or
    /// mollification did not make them usable).
    DegenerateTriangles {
        /// Number of degenerate triangles.
        count: usize,
        /// The polygon of the first one.
        polygon: u32,
        /// The vertices of the first one.
        vertices: [u32; 3],
        /// Its edge lengths `v0v1`, `v1v2`, `v2v0`.
        lengths: [f64; 3],
        /// Whether mollification was tried.
        after_mollification: bool,
    },
    /// A connected piece of the surface is closed (it has no boundary after the seams were
    /// cut) and the policy for closed surfaces is to refuse.
    ClosedSurface {
        /// Lowest-numbered polygon of the piece.
        polygon: u32,
        /// Number of polygons of the piece.
        n_polygons: usize,
        /// Its genus.
        genus: usize,
    },
    /// Only positive cones were asked for, on a closed piece whose cone angles must sum to a
    /// negative number.
    NegativeConesRequired {
        /// Lowest-numbered polygon of the piece.
        polygon: u32,
        /// Its genus.
        genus: usize,
    },
    /// The steering leaves no cone configuration that the topology of a piece allows.
    InfeasibleSteering {
        /// Lowest-numbered polygon of the piece.
        polygon: u32,
        /// What the cone solver said.
        message: String,
    },

    // ----- solve -----
    /// The cone solver refused the problem of a piece or failed.
    ConeSolver {
        /// Lowest-numbered polygon of the piece.
        polygon: u32,
        /// What it said.
        message: String,
    },
    /// The cone solver did not reach a certified optimum and `require_certified_cones` is set.
    ConesNotCertified {
        /// Lowest-numbered polygon of the piece.
        polygon: u32,
        /// The solver's status and message.
        message: String,
    },
    /// A cone angle of `2π` or more came out of the cone stage; no flattening has such a cone.
    UnrealizableCone {
        /// The vertex.
        vertex: u32,
        /// The angle.
        angle: f64,
    },
    /// The flattener refused a piece or failed.
    Flatten {
        /// Lowest-numbered polygon of the piece.
        polygon: u32,
        /// What it said.
        message: String,
    },
    /// A check of the pipeline on its own result failed. This is a bug in the pipeline.
    Internal {
        /// What was found.
        message: String,
    },
}

impl UvError {
    /// Which family the error belongs to.
    pub fn kind(&self) -> ErrorKind {
        use UvError::*;
        match self {
            NoPolygons
            | TooLarge { .. }
            | PolygonOffsets { .. }
            | PolygonTooSmall { .. }
            | IndexOutOfRange { .. }
            | NonFinite { .. }
            | TriangleCount { .. }
            | BadTriangulation { .. }
            | DegenerateSeamEdge { .. }
            | DuplicatePrescribed { .. }
            | PrescribedAngleTooLarge { .. }
            | ConflictingSteering { .. }
            | InvalidOption { .. } => ErrorKind::Input,
            RepeatedVertexInPolygon { .. }
            | DuplicatePolygons { .. }
            | NonManifoldEdge { .. }
            | DiagonalOnEdge { .. }
            | NonManifoldVertex { .. }
            | NonOrientable { .. }
            | InconsistentOrientation { .. }
            | DegenerateTriangles { .. }
            | ClosedSurface { .. }
            | NegativeConesRequired { .. }
            | InfeasibleSteering { .. } => ErrorKind::Mesh,
            ConeSolver { .. }
            | ConesNotCertified { .. }
            | UnrealizableCone { .. }
            | Flatten { .. }
            | Internal { .. } => ErrorKind::Solve,
        }
    }
}

/// Writes at most 16 items of a list, then the number of items left out.
fn write_list(f: &mut fmt::Formatter<'_>, items: &[u32]) -> fmt::Result {
    const MAX: usize = 16;
    for (i, item) in items.iter().take(MAX).enumerate() {
        if i > 0 {
            write!(f, ", ")?;
        }
        write!(f, "{item}")?;
    }
    if items.len() > MAX {
        write!(f, ", ... ({} more)", items.len() - MAX)?;
    }
    Ok(())
}

impl fmt::Display for UvError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use UvError::*;
        match self {
            NoPolygons => write!(f, "the mesh has no faces"),
            TooLarge { what, count } => write!(
                f,
                "{what} has {count} entries, more than 32-bit ids can address (limit {})",
                u32::MAX - 1
            ),
            PolygonOffsets {
                position,
                value,
                requirement,
            } => write!(f, "polygon_offsets[{position}] = {value}: {requirement}"),
            PolygonTooSmall { polygon, n_corners } => write!(
                f,
                "face {polygon} has {n_corners} corners; a face needs at least 3"
            ),
            IndexOutOfRange {
                what,
                position,
                value,
                count,
                of,
            } => {
                if *count == 0 {
                    write!(
                        f,
                        "{what}[{position}] = {value} is out of range: there are no {of}"
                    )
                } else {
                    write!(
                        f,
                        "{what}[{position}] = {value} is out of range: there are {count} {of} \
                         (valid indices 0..={})",
                        count - 1
                    )
                }
            }
            NonFinite {
                what,
                position,
                value,
            } => write!(f, "{what}[{position}] = {value} is not a finite number"),
            TriangleCount { got, expected } => write!(
                f,
                "triangles has {got} triangles, but the faces need {expected} \
                 (corners minus twice the number of faces)"
            ),
            BadTriangulation { polygon, problem } => write!(
                f,
                "the triangles of face {polygon} do not triangulate it: {problem}"
            ),
            DegenerateSeamEdge { position, vertex } => write!(
                f,
                "seam_edges[{position}] joins vertex {vertex} to itself"
            ),
            DuplicatePrescribed { vertex } => {
                write!(f, "vertex {vertex} has more than one prescribed cone")
            }
            PrescribedAngleTooLarge { vertex, angle } => write!(
                f,
                "the cone angle {angle} prescribed at vertex {vertex} is not below 2π \
                 ({}): no flattening has such a cone",
                std::f64::consts::TAU
            ),
            ConflictingSteering {
                vertex,
                first,
                second,
            } => write!(f, "vertex {vertex} is both {first} and {second}"),
            InvalidOption {
                name,
                value,
                requirement,
            } => write!(f, "option {name} = {value}: {requirement}"),

            RepeatedVertexInPolygon { polygon, vertex } => {
                write!(f, "face {polygon} uses vertex {vertex} more than once")
            }
            DuplicatePolygons {
                first,
                second,
                vertices: [a, b, c],
            } => {
                if first == second {
                    write!(
                        f,
                        "face {first} covers the triangle of vertices ({a}, {b}, {c}) twice"
                    )
                } else {
                    write!(
                        f,
                        "faces {first} and {second} are duplicates or overlap: both cover the \
                         triangle of vertices ({a}, {b}, {c})"
                    )
                }
            }
            NonManifoldEdge {
                edge: [a, b],
                polygons,
            } => {
                write!(f, "edge ({a}, {b}) is shared by {} faces: ", polygons.len())?;
                write_list(f, polygons)
            }
            DiagonalOnEdge {
                edge: [a, b],
                polygons,
                n_triangles,
            } => {
                write!(
                    f,
                    "edge ({a}, {b}) is shared by {n_triangles} triangles of the triangulation, \
                     in faces "
                )?;
                write_list(f, polygons)?;
                write!(f, ": a diagonal of a face coincides with another edge")
            }
            NonManifoldVertex { vertex, fans } => write!(
                f,
                "vertex {vertex} is non-manifold: its incident faces form {fans} separate fans"
            ),
            NonOrientable {
                edge: [a, b],
                polygons: [p, q],
            } => {
                if p == q {
                    // The search closed its loop inside one face, across a diagonal.
                    write!(
                        f,
                        "the surface is not orientable: the faces connected to face {p} have \
                         no consistent orientation (the contradiction closes at the diagonal \
                         ({a}, {b}) of that face)"
                    )
                } else {
                    write!(
                        f,
                        "the surface is not orientable: faces {p} and {q} cannot be oriented \
                         consistently across edge ({a}, {b})"
                    )
                }
            }
            InconsistentOrientation {
                edge: [a, b],
                polygons: [p, q],
            } => write!(
                f,
                "faces {p} and {q} traverse their shared edge ({a}, {b}) in the same \
                 direction: the mesh is not consistently oriented (one of them is flipped)"
            ),
            DegenerateTriangles {
                count,
                polygon,
                vertices: [a, b, c],
                lengths: [l0, l1, l2],
                after_mollification,
            } => {
                write!(
                    f,
                    "{count} {} of the triangulation {} no area",
                    if *count == 1 { "triangle" } else { "triangles" },
                    if *count == 1 { "has" } else { "have" },
                )?;
                if *after_mollification {
                    write!(f, " even after mollifying the edge lengths")?;
                }
                write!(
                    f,
                    "; the first is in face {polygon} (vertices {a}, {b}, {c}; edge lengths \
                     {l0}, {l1}, {l2})"
                )
            }
            ClosedSurface {
                polygon,
                n_polygons,
                genus,
            } => write!(
                f,
                "the piece of the surface that contains face {polygon} ({n_polygons} faces, \
                 genus {genus}) is closed: mark a seam on it"
            ),
            NegativeConesRequired { polygon, genus } => write!(
                f,
                "the piece of the surface that contains face {polygon} is closed with genus \
                 {genus}: its cone angles must sum to {} (2π times the Euler characteristic), \
                 which positive cones cannot do; allow negative cones or mark seams",
                std::f64::consts::TAU * (2.0 - 2.0 * *genus as f64)
            ),
            InfeasibleSteering { polygon, message } => write!(
                f,
                "the cone constraints on the piece of the surface that contains face \
                 {polygon} cannot be met: {message}"
            ),

            ConeSolver { polygon, message } => write!(
                f,
                "the cone solver failed on the piece of the surface that contains face \
                 {polygon}: {message}"
            ),
            ConesNotCertified { polygon, message } => write!(
                f,
                "the cones of the piece of the surface that contains face {polygon} are not a \
                 certified optimum: {message}"
            ),
            UnrealizableCone { vertex, angle } => write!(
                f,
                "the cone stage produced a cone angle of {angle} at vertex {vertex}, which is \
                 not below 2π: no flattening has such a cone"
            ),
            Flatten { polygon, message } => write!(
                f,
                "flattening failed on the piece of the surface that contains face {polygon}: \
                 {message}"
            ),
            Internal { message } => write!(f, "internal error: {message}"),
        }
    }
}

impl std::error::Error for UvError {}

/// Shorthand for an internal error.
pub(crate) fn internal(message: impl Into<String>) -> UvError {
    UvError::Internal {
        message: message.into(),
    }
}
