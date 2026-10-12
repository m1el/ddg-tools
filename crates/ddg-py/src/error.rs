//! Error types of the binding.
//!
//! [`InputError`] describes arguments that are malformed as arrays (wrong element type, wrong
//! length, an index out of range, contradictory steering). [`UnwrapError`] adds the two ways a
//! well-formed call can still fail: the mesh is not usable ([`UnwrapError::Mesh`]) or the pipeline
//! could not produce a map ([`UnwrapError::Solve`]). The Python module maps them to
//! `ddg_tools.InputError`, `ddg_tools.MeshError` and `ddg_tools.SolveError`, all subclasses of
//! `ddg_tools.DdgError`.

use std::fmt;

/// An argument is malformed. Every message starts with the name of the Python parameter.
#[derive(Debug, Clone, PartialEq)]
pub enum InputError {
    /// The buffer is not C-contiguous (for example a strided numpy slice).
    NotContiguous {
        /// Parameter name.
        name: &'static str,
    },
    /// The element type of a buffer is not one this parameter accepts.
    UnsupportedFormat {
        /// Parameter name.
        name: &'static str,
        /// The buffer's `struct`-module format string.
        format: String,
        /// The buffer's item size in bytes.
        item_size: usize,
        /// What the parameter accepts, as prose.
        expected: &'static str,
    },
    /// A buffer has more than two dimensions, or a second dimension of the wrong width.
    BadShape {
        /// Parameter name.
        name: &'static str,
        /// The buffer's shape.
        shape: Vec<usize>,
        /// Required number of items per row.
        width: usize,
    },
    /// A flat buffer's item count is not a multiple of the row width.
    LengthNotMultiple {
        /// Parameter name.
        name: &'static str,
        /// Number of items in the buffer.
        len: usize,
        /// Required number of items per row.
        width: usize,
        /// What one row is, for the message ("vertex", "triangle", ...).
        row: &'static str,
    },
    /// A byte buffer (`bytes`, `bytearray`, format `B`) whose length is not a multiple of the
    /// size of the element type it is reinterpreted as.
    RawByteLength {
        /// Parameter name.
        name: &'static str,
        /// Length of the buffer in bytes.
        len_bytes: usize,
        /// The element type the bytes stand for ("float32" or "int32").
        element: &'static str,
    },
    /// An index is negative.
    NegativeIndex {
        /// Parameter name.
        name: &'static str,
        /// Flat position in the buffer.
        position: usize,
        /// The value.
        value: i64,
    },
    /// An index does not fit the 32-bit ids used by the Rust side.
    IndexTooLarge {
        /// Parameter name.
        name: &'static str,
        /// Flat position in the buffer.
        position: usize,
        /// The value.
        value: u64,
    },
    /// An index is not below the number of elements it refers to.
    IndexOutOfRange {
        /// Parameter name.
        name: &'static str,
        /// Flat position in the buffer.
        position: usize,
        /// The value.
        value: u32,
        /// Number of elements of the kind the index refers to.
        count: usize,
        /// That kind, plural ("vertices", "polygon corners").
        of: &'static str,
    },
    /// A coordinate or an angle is NaN or infinite.
    NonFinite {
        /// Parameter name.
        name: &'static str,
        /// Flat position in the buffer.
        position: usize,
        /// The value.
        value: f64,
    },
    /// The mesh has more elements than 32-bit ids can address.
    TooLarge {
        /// Parameter name.
        name: &'static str,
        /// Number of rows in the buffer.
        count: usize,
    },
    /// `polygon_starts` is empty: there is nothing to unwrap.
    NoPolygons,
    /// `polygon_starts[0]` is not 0.
    FirstPolygonStart {
        /// The value found.
        start: u32,
    },
    /// A polygon has fewer than three corners (its start offset is too close to the next one,
    /// or to the end of `corner_vertices`), or the offsets decrease.
    PolygonTooSmall {
        /// The polygon.
        polygon: usize,
        /// Its start offset.
        start: u32,
        /// The next polygon's start offset, or the number of corners for the last polygon.
        end: u32,
    },
    /// The number of triangles is not `corners − 2 · polygons`.
    TriangleCount {
        /// Number of triangles given.
        got: usize,
        /// Number expected.
        expected: usize,
        /// Number of polygon corners.
        n_corners: usize,
        /// Number of polygons.
        n_polygons: usize,
    },
    /// The three corners of a triangle do not belong to one polygon.
    TriangleSpansPolygons {
        /// The triangle.
        triangle: usize,
        /// Its three corner indices.
        corners: [u32; 3],
        /// The polygons those corners belong to.
        polygons: [u32; 3],
    },
    /// The triangles of a polygon are not a triangulation of it.
    BadTriangulation {
        /// The polygon.
        polygon: u32,
        /// Number of corners of the polygon.
        n_corners: usize,
        /// What is wrong.
        problem: TriangulationProblem,
    },
    /// A seam edge joins a vertex to itself.
    DegenerateSeamEdge {
        /// Row in `seam_edges`.
        row: usize,
        /// The vertex.
        vertex: u32,
    },
    /// `prescribed_vertices` and `prescribed_angles` differ in length, or only one was given.
    PrescribedLengthMismatch {
        /// Length of `prescribed_vertices`.
        vertices: usize,
        /// Length of `prescribed_angles`.
        angles: usize,
    },
    /// A prescribed cone angle is not below 2π (a cone angle is 2π minus a positive angle sum).
    PrescribedAngleTooLarge {
        /// The vertex.
        vertex: u32,
        /// The angle.
        angle: f64,
    },
    /// A vertex is listed twice in `prescribed_vertices`.
    DuplicatePrescribed {
        /// The vertex.
        vertex: u32,
    },
    /// A vertex is in more than one of the pinned, prescribed and forbidden sets.
    ConflictingSteering {
        /// The vertex.
        vertex: u32,
        /// The first set it is in.
        first: &'static str,
        /// The second set it is in.
        second: &'static str,
    },
    /// An option is outside its domain.
    BadOption {
        /// Parameter name.
        name: &'static str,
        /// The value, formatted.
        value: String,
        /// The requirement, as prose.
        requirement: &'static str,
    },
    /// Two parameters that describe the same mesh disagree in length.
    LengthMismatch {
        /// Parameter name.
        name: &'static str,
        /// Its number of rows.
        got: usize,
        /// The number of rows required.
        expected: usize,
        /// Why that number is required.
        reason: &'static str,
    },
}

/// How the triangles of one polygon fail to triangulate it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriangulationProblem {
    /// The polygon has this many triangles instead of `corners − 2`.
    TriangleCount(usize),
    /// A triangle uses the same corner twice.
    RepeatedCorner {
        /// The triangle.
        triangle: usize,
    },
    /// A triangle runs against the winding of the polygon.
    ReversedTriangle {
        /// The triangle.
        triangle: usize,
    },
    /// A side of the polygon is covered by no triangle, or by more than one.
    SideCoverage {
        /// First corner of the side (the side runs to the next corner of the polygon).
        corner: u32,
        /// Number of triangles on it.
        count: usize,
    },
    /// A diagonal is not shared by exactly two triangles in opposite directions.
    UnpairedDiagonal {
        /// The two corners of the diagonal.
        corners: [u32; 2],
    },
}

impl fmt::Display for InputError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InputError::NotContiguous { name } => write!(
                f,
                "{name}: the buffer is not C-contiguous; pass a contiguous copy \
                 (for numpy: numpy.ascontiguousarray)"
            ),
            InputError::UnsupportedFormat {
                name,
                format,
                item_size,
                expected,
            } => write!(
                f,
                "{name}: unsupported buffer element type (format {format:?}, {item_size} \
                 bytes per item); expected {expected}"
            ),
            InputError::BadShape { name, shape, width } => {
                write!(f, "{name}: buffer has shape (")?;
                for (i, s) in shape.iter().enumerate() {
                    if i > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{s}")?;
                }
                if *width == 1 {
                    write!(f, "); expected a one-dimensional buffer")
                } else {
                    write!(f, "); expected a flat buffer or one of shape (n, {width})")
                }
            }
            InputError::LengthNotMultiple {
                name,
                len,
                width,
                row,
            } => write!(
                f,
                "{name}: {len} items is not a multiple of {width} ({width} per {row})"
            ),
            InputError::RawByteLength {
                name,
                len_bytes,
                element,
            } => write!(
                f,
                "{name}: a byte buffer is read as native-endian {element}, but its length \
                 {len_bytes} is not a multiple of 4"
            ),
            InputError::NegativeIndex {
                name,
                position,
                value,
            } => write!(f, "{name}[{position}] = {value} is negative"),
            InputError::IndexTooLarge {
                name,
                position,
                value,
            } => write!(
                f,
                "{name}[{position}] = {value} does not fit a 32-bit index (limit {})",
                u32::MAX - 1
            ),
            InputError::IndexOutOfRange {
                name,
                position,
                value,
                count,
                of,
            } => {
                if *count == 0 {
                    write!(
                        f,
                        "{name}[{position}] = {value} is out of range: there are no {of}"
                    )
                } else {
                    write!(
                        f,
                        "{name}[{position}] = {value} is out of range: there are {count} {of} \
                         (valid indices 0..={})",
                        count - 1
                    )
                }
            }
            InputError::NonFinite {
                name,
                position,
                value,
            } => write!(f, "{name}[{position}] = {value} is not a finite number"),
            InputError::TooLarge { name, count } => write!(
                f,
                "{name}: {count} rows are more than 32-bit ids can address (limit {})",
                u32::MAX - 1
            ),
            InputError::NoPolygons => write!(f, "polygon_starts is empty: the mesh has no faces"),
            InputError::FirstPolygonStart { start } => write!(
                f,
                "polygon_starts[0] = {start}, but the first polygon must start at corner 0"
            ),
            InputError::PolygonTooSmall {
                polygon,
                start,
                end,
            } => {
                if end < start {
                    write!(
                        f,
                        "polygon_starts[{polygon}] = {start} is followed by {end}: polygon start \
                         offsets must increase"
                    )
                } else {
                    write!(
                        f,
                        "polygon {polygon} has {} corners (corners {start}..{end}); a polygon \
                         needs at least 3",
                        end - start
                    )
                }
            }
            InputError::TriangleCount {
                got,
                expected,
                n_corners,
                n_polygons,
            } => write!(
                f,
                "triangles: got {got} triangles, expected {expected} \
                 ({n_corners} corners - 2 x {n_polygons} polygons)"
            ),
            InputError::TriangleSpansPolygons {
                triangle,
                corners: [a, b, c],
                polygons: [p, q, r],
            } => write!(
                f,
                "triangles: triangle {triangle} = corners ({a}, {b}, {c}) spans polygons \
                 {p}, {q}, {r}; the corners of a triangle must belong to one polygon"
            ),
            InputError::BadTriangulation {
                polygon,
                n_corners,
                problem,
            } => {
                write!(
                    f,
                    "triangles: polygon {polygon} ({n_corners} corners) is not triangulated \
                     correctly: "
                )?;
                match problem {
                    TriangulationProblem::TriangleCount(n) => {
                        write!(f, "it has {n} triangles instead of {}", n_corners - 2)
                    }
                    TriangulationProblem::RepeatedCorner { triangle } => {
                        write!(f, "triangle {triangle} uses a corner twice")
                    }
                    TriangulationProblem::ReversedTriangle { triangle } => write!(
                        f,
                        "triangle {triangle} runs against the winding of the polygon"
                    ),
                    TriangulationProblem::SideCoverage { corner, count } => write!(
                        f,
                        "the side starting at corner {corner} is covered by {count} triangles \
                         instead of 1"
                    ),
                    TriangulationProblem::UnpairedDiagonal { corners: [a, b] } => write!(
                        f,
                        "the diagonal between corners {a} and {b} is not shared by exactly two \
                         triangles in opposite directions"
                    ),
                }
            }
            InputError::DegenerateSeamEdge { row, vertex } => {
                write!(f, "seam_edges: edge {row} joins vertex {vertex} to itself")
            }
            InputError::PrescribedLengthMismatch { vertices, angles } => write!(
                f,
                "prescribed_vertices has {vertices} entries but prescribed_angles has {angles}; \
                 they must be given together and have equal length"
            ),
            InputError::PrescribedAngleTooLarge { vertex, angle } => write!(
                f,
                "prescribed_angles: cone angle {angle} at vertex {vertex} is not below 2\u{3c0}"
            ),
            InputError::DuplicatePrescribed { vertex } => write!(
                f,
                "prescribed_vertices: vertex {vertex} is listed more than once"
            ),
            InputError::ConflictingSteering {
                vertex,
                first,
                second,
            } => write!(
                f,
                "vertex {vertex} is both {first} and {second}; a vertex can be in only one of \
                 the pinned, prescribed and forbidden sets"
            ),
            InputError::BadOption {
                name,
                value,
                requirement,
            } => write!(f, "{name} = {value}: {requirement}"),
            InputError::LengthMismatch {
                name,
                got,
                expected,
                reason,
            } => write!(f, "{name}: got {got} rows, expected {expected} ({reason})"),
        }
    }
}

impl std::error::Error for InputError {}

/// Why an unwrap call failed.
#[derive(Debug, Clone, PartialEq)]
pub enum UnwrapError {
    /// An argument is malformed.
    Input(InputError),
    /// The arrays are well formed but do not describe a mesh the pipeline accepts
    /// (non-manifold, inconsistently oriented, ...). The message says what and where, with vertex
    /// and polygon indices of the caller's mesh.
    Mesh(String),
    /// The pipeline accepted the mesh but could not produce a valid map.
    Solve(String),
}

impl fmt::Display for UnwrapError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            UnwrapError::Input(e) => e.fmt(f),
            UnwrapError::Mesh(message) | UnwrapError::Solve(message) => f.write_str(message),
        }
    }
}

impl std::error::Error for UnwrapError {}

impl From<InputError> for UnwrapError {
    fn from(e: InputError) -> Self {
        UnwrapError::Input(e)
    }
}
