//! Intrinsic triangle geometry: corner angles, areas and cotangents computed from edge lengths.
//!
//! The algorithms of this workspace need only a positive length per edge (a piecewise Euclidean
//! metric), so [`Geometry`] is built from lengths; positions are just one way to obtain them
//! ([`Geometry::from_positions`]).
//!
//! The scalar formulas ([`triangle_area`], [`corner_angle`], [`corner_cotan`]) are the versions
//! recommended by W. Kahan, *Miscalculating Area and Angles of a Needle-like Triangle* (2014):
//! they stay accurate for needle and cap triangles, where the textbook forms lose all digits.
//!
//! # Degenerate triangles
//!
//! A face is rejected, with an error that names it, when its three lengths do not form a
//! triangle of strictly positive area: the strict triangle inequality fails in floating point
//! arithmetic, or the computed area or a computed cotangent is not a positive / finite number
//! (possible only for lengths near the limits of `f64`). There is no tolerance and nothing is
//! clamped: a thin but valid triangle is accepted and yields large, finite cotangents. Every
//! value stored in a [`Geometry`] is finite.

use std::fmt;

use crate::topology::TriMesh;

/// Why a geometry could not be built.
#[derive(Debug, Clone, PartialEq)]
pub enum GeometryError {
    /// The number of positions is not the number of mesh vertices.
    PositionCountMismatch {
        /// Number of vertices of the mesh.
        expected: usize,
        /// Number of positions passed.
        got: usize,
    },
    /// A position has a NaN or infinite coordinate.
    NonFinitePosition {
        /// The vertex.
        vertex: u32,
        /// Its position.
        position: [f64; 3],
    },
    /// The number of edge lengths is not the number of mesh edges.
    EdgeLengthCountMismatch {
        /// Number of edges of the mesh.
        expected: usize,
        /// Number of lengths passed.
        got: usize,
    },
    /// An edge length is negative, NaN or infinite.
    InvalidEdgeLength {
        /// The edge id.
        edge: u32,
        /// Its vertices `[min, max]`.
        vertices: [u32; 2],
        /// The offending length.
        length: f64,
        /// Total number of edges with an invalid length.
        count: usize,
    },
    /// One or more faces are degenerate; all of them are listed, in ascending order.
    DegenerateFaces {
        /// The degenerate faces.
        faces: Vec<DegenerateFace>,
    },
    /// A [`Geometry`] was used with a mesh it was not built for.
    MeshMismatch {
        /// `(n_edges, n_faces)` of the mesh.
        mesh: (usize, usize),
        /// `(n_edges, n_faces)` of the geometry.
        geometry: (usize, usize),
    },
}

/// A face whose edge lengths do not form a triangle of positive area.
#[derive(Debug, Clone, PartialEq)]
pub struct DegenerateFace {
    /// The face.
    pub face: u32,
    /// Its vertices.
    pub vertices: [u32; 3],
    /// Lengths of its edges `v0v1`, `v1v2`, `v2v0`.
    pub lengths: [f64; 3],
    /// What is wrong with it.
    pub reason: DegenerateReason,
}

/// How a triangle fails to have positive area.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegenerateReason {
    /// The longest side equals the sum of the other two: the corners are collinear or coincide.
    ZeroArea,
    /// The longest side exceeds the sum of the other two: no such triangle exists.
    TriangleInequality,
    /// The lengths are so small or so large that the area or a cotangent under- or overflows.
    OutOfRange,
}

impl fmt::Display for DegenerateFace {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let [a, b, c] = self.vertices;
        let [l0, l1, l2] = self.lengths;
        write!(f, "face {} (vertices {a}, {b}, {c}; edge lengths {l0}, {l1}, {l2}) ", self.face)?;
        match self.reason {
            DegenerateReason::ZeroArea => write!(f, "has zero area"),
            DegenerateReason::TriangleInequality => write!(f, "violates the triangle inequality"),
            DegenerateReason::OutOfRange => write!(f, "has lengths outside the range where its area is computable"),
        }
    }
}

impl fmt::Display for GeometryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GeometryError::PositionCountMismatch { expected, got } => {
                write!(f, "the mesh has {expected} vertices but {got} positions were given")
            }
            GeometryError::NonFinitePosition { vertex, position: [x, y, z] } => {
                write!(f, "vertex {vertex} has a non-finite position ({x}, {y}, {z})")
            }
            GeometryError::EdgeLengthCountMismatch { expected, got } => {
                write!(f, "the mesh has {expected} edges but {got} edge lengths were given")
            }
            GeometryError::InvalidEdgeLength { edge, vertices: [a, b], length, count } => {
                write!(f, "edge {edge} = ({a}, {b}) has invalid length {length}")?;
                if *count > 1 {
                    write!(f, " ({count} edges have a negative or non-finite length)")?;
                }
                Ok(())
            }
            GeometryError::DegenerateFaces { faces } => {
                const MAX: usize = 8;
                let n = faces.len();
                write!(f, "{n} degenerate {}: ", if n == 1 { "face" } else { "faces" })?;
                for (i, face) in faces.iter().take(MAX).enumerate() {
                    if i > 0 {
                        write!(f, "; ")?;
                    }
                    write!(f, "{face}")?;
                }
                if n > MAX {
                    write!(f, "; ... ({} more)", n - MAX)?;
                }
                Ok(())
            }
            GeometryError::MeshMismatch { mesh, geometry } => write!(
                f,
                "geometry does not belong to this mesh: the mesh has {} edges and {} faces, the \
                 geometry {} edges and {} faces",
                mesh.0, mesh.1, geometry.0, geometry.1
            ),
        }
    }
}

impl std::error::Error for GeometryError {}

/// Euclidean distance between two points, `ℓ = √(Δx² + Δy² + Δz²)`.
///
/// [formula: edge-length-euclidean]
#[inline]
pub fn edge_length(p: [f64; 3], q: [f64; 3]) -> f64 {
    let (dx, dy, dz) = (q[0] - p[0], q[1] - p[1], q[2] - p[2]);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Sorts three lengths in descending order and classifies the triangle they form.
/// `Ok` means all lengths are positive and finite and the strict triangle inequality holds.
fn sorted_valid(a: f64, b: f64, c: f64) -> Result<[f64; 3], DegenerateReason> {
    if !(a.is_finite() && b.is_finite() && c.is_finite()) {
        return Err(DegenerateReason::OutOfRange);
    }
    let mut s = [a, b, c];
    if s[0] < s[1] {
        s.swap(0, 1);
    }
    if s[1] < s[2] {
        s.swap(1, 2);
    }
    if s[0] < s[1] {
        s.swap(0, 1);
    }
    // x >= y >= z; the triangle inequality is z > x - y. The subtraction is exact whenever the
    // inequality is close to equality (Sterbenz), so this test has no rounding error of its own.
    // It also implies z > 0.
    let slack = s[2] - (s[0] - s[1]);
    if slack > 0.0 {
        Ok(s)
    } else if slack == 0.0 {
        Err(DegenerateReason::ZeroArea)
    } else {
        Err(DegenerateReason::TriangleInequality)
    }
}

/// Area of the triangle with side lengths `a`, `b`, `c` by Heron's formula
/// `16 A² = (a + b + c)(−a + b + c)(a − b + c)(a + b − c)`, evaluated in Kahan's stable
/// form: with the sides sorted `x ≥ y ≥ z`,
/// `A = ¼ √( (x + (y + z)) (z − (x − y)) (z + (x − y)) (x + (y − z)) )`.
///
/// Returns `None` unless the lengths are positive, finite, satisfy the strict triangle
/// inequality, and give a positive finite area.
///
/// [formula: triangle-area-heron]
pub fn triangle_area(a: f64, b: f64, c: f64) -> Option<f64> {
    let [x, y, z] = sorted_valid(a, b, c).ok()?;
    let area = 0.25 * ((x + (y + z)) * (z - (x - y)) * (z + (x - y)) * (x + (y - z))).sqrt();
    (area > 0.0 && area.is_finite()).then_some(area)
}

/// Interior angle of a triangle opposite the side of length `opposite`, the other two sides
/// having lengths `side1` and `side2`, from the half-angle formula
/// `tan²(θ/2) = ((a − b + c)(a + b − c)) / ((a + b + c)(−a + b + c))` with `a = opposite`,
/// evaluated in Kahan's stable form. The result lies in `(0, π)`.
///
/// Returns `None` unless the three lengths form a triangle of positive area
/// (see [`triangle_area`]).
///
/// [formula: corner-angle-from-lengths]
pub fn corner_angle(opposite: f64, side1: f64, side2: f64) -> Option<f64> {
    sorted_valid(opposite, side1, side2).ok()?;
    // Kahan's notation: the angle C is opposite c and lies between a >= b.
    let c = opposite;
    let (a, b) = if side1 >= side2 { (side1, side2) } else { (side2, side1) };
    let mu = if b >= c { c - (a - b) } else { b - (a - c) };
    let num = ((a - b) + c) * mu;
    let den = (a + (b + c)) * ((a - c) + b);
    let angle = 2.0 * num.sqrt().atan2(den.sqrt());
    (angle > 0.0 && angle.is_finite()).then_some(angle)
}

/// Cotangent of the interior angle opposite the side of length `opposite`:
/// `cot θ = (b² + c² − a²) / (4 A)` with `a = opposite`, `b`, `c` the other two sides and `A`
/// the triangle area ([`triangle_area`]). The numerator is evaluated as
/// `c² + (b − a)(b + a)` with `b ≥ c` to avoid cancellation.
///
/// Returns `None` unless the three lengths form a triangle of positive area and the result is
/// finite.
///
/// [formula: cotan-from-lengths]
pub fn corner_cotan(opposite: f64, side1: f64, side2: f64) -> Option<f64> {
    let area = triangle_area(opposite, side1, side2)?;
    let (b, c) = if side1 >= side2 { (side1, side2) } else { (side2, side1) };
    let cot = (c * c + (b - opposite) * (b + opposite)) / (4.0 * area);
    cot.is_finite().then_some(cot)
}

/// Edge lengths of a mesh and everything derived from them per face and per corner.
///
/// Corner data is stored per face as `[_; 3]`, where slot `k` belongs to the corner at vertex
/// `faces[f][k]`, i.e. to corner id `3 f + k` (see [`crate::topology`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Geometry {
    edge_lengths: Vec<f64>,
    corner_angles: Vec<[f64; 3]>,
    corner_cotans: Vec<[f64; 3]>,
    face_areas: Vec<f64>,
}

impl Geometry {
    /// Builds the geometry from vertex positions.
    ///
    /// # Errors
    /// [`GeometryError::PositionCountMismatch`], [`GeometryError::NonFinitePosition`],
    /// [`GeometryError::InvalidEdgeLength`] (only if a coordinate difference overflows), and
    /// [`GeometryError::DegenerateFaces`] naming every face without positive area — in
    /// particular faces with coincident or collinear vertices.
    pub fn from_positions(mesh: &TriMesh, positions: &[[f64; 3]]) -> Result<Geometry, GeometryError> {
        if positions.len() != mesh.n_vertices() {
            return Err(GeometryError::PositionCountMismatch { expected: mesh.n_vertices(), got: positions.len() });
        }
        if let Some(v) = positions.iter().position(|p| !p.iter().all(|x| x.is_finite())) {
            return Err(GeometryError::NonFinitePosition { vertex: v as u32, position: positions[v] });
        }
        let lengths = mesh
            .edges()
            .iter()
            .map(|&[a, b]| edge_length(positions[a as usize], positions[b as usize]))
            .collect();
        Geometry::from_edge_lengths(mesh, lengths)
    }

    /// Builds the geometry from one length per edge (`edge_lengths[e]` for edge id `e`, see
    /// [`TriMesh::edges`]). This is the intrinsic constructor: no embedding is needed.
    ///
    /// # Errors
    /// [`GeometryError::EdgeLengthCountMismatch`]; [`GeometryError::InvalidEdgeLength`] for a
    /// negative, NaN or infinite length; [`GeometryError::DegenerateFaces`] naming every face
    /// whose lengths do not form a triangle of positive area (a zero length makes its faces
    /// degenerate and is reported this way).
    pub fn from_edge_lengths(mesh: &TriMesh, edge_lengths: Vec<f64>) -> Result<Geometry, GeometryError> {
        if edge_lengths.len() != mesh.n_edges() {
            return Err(GeometryError::EdgeLengthCountMismatch { expected: mesh.n_edges(), got: edge_lengths.len() });
        }
        let invalid = |l: &f64| !(l.is_finite() && *l >= 0.0);
        if let Some(e) = edge_lengths.iter().position(invalid) {
            return Err(GeometryError::InvalidEdgeLength {
                edge: e as u32,
                vertices: mesh.edges()[e],
                length: edge_lengths[e],
                count: edge_lengths.iter().filter(|l| invalid(l)).count(),
            });
        }

        let n_faces = mesh.n_faces();
        let mut corner_angles = Vec::with_capacity(n_faces);
        let mut corner_cotans = Vec::with_capacity(n_faces);
        let mut face_areas = Vec::with_capacity(n_faces);
        let mut degenerate = Vec::new();
        for f in 0..n_faces {
            // l[k] is the length of halfedge 3f + k, from corner k to corner k + 1; the side
            // opposite corner k is halfedge 3f + (k + 1) % 3.
            let l = [0, 1, 2].map(|k| edge_lengths[mesh.edge(3 * f + k)]);
            match face_geometry(l) {
                Ok((area, angles, cotans)) => {
                    face_areas.push(area);
                    corner_angles.push(angles);
                    corner_cotans.push(cotans);
                }
                Err(reason) => {
                    degenerate.push(DegenerateFace { face: f as u32, vertices: mesh.faces()[f], lengths: l, reason })
                }
            }
        }
        if !degenerate.is_empty() {
            return Err(GeometryError::DegenerateFaces { faces: degenerate });
        }
        Ok(Geometry { edge_lengths, corner_angles, corner_cotans, face_areas })
    }

    /// Length of every edge, indexed by edge id.
    #[inline]
    pub fn edge_lengths(&self) -> &[f64] {
        &self.edge_lengths
    }

    /// Interior angles in radians: `corner_angles()[f][k]` is the angle of face `f` at vertex
    /// `faces[f][k]`. Each lies in `(0, π)`.
    #[inline]
    pub fn corner_angles(&self) -> &[[f64; 3]] {
        &self.corner_angles
    }

    /// Cotangents of the interior angles, laid out like [`Geometry::corner_angles`].
    /// `corner_cotans()[f][k]` is the weight that face `f` contributes to the edge opposite
    /// corner `k`, i.e. to halfedge `3 f + (k + 1) % 3`.
    #[inline]
    pub fn corner_cotans(&self) -> &[[f64; 3]] {
        &self.corner_cotans
    }

    /// Area of every face. Each is strictly positive.
    #[inline]
    pub fn face_areas(&self) -> &[f64] {
        &self.face_areas
    }

    /// Sum of all face areas.
    pub fn total_area(&self) -> f64 {
        self.face_areas.iter().sum()
    }

    /// Number of faces this geometry describes.
    #[inline]
    pub fn n_faces(&self) -> usize {
        self.face_areas.len()
    }

    /// Number of edges this geometry describes.
    #[inline]
    pub fn n_edges(&self) -> usize {
        self.edge_lengths.len()
    }

    /// Checks that this geometry has the edge and face counts of `mesh`.
    ///
    /// # Errors
    /// [`GeometryError::MeshMismatch`].
    pub fn check_mesh(&self, mesh: &TriMesh) -> Result<(), GeometryError> {
        if self.n_edges() == mesh.n_edges() && self.n_faces() == mesh.n_faces() {
            Ok(())
        } else {
            Err(GeometryError::MeshMismatch {
                mesh: (mesh.n_edges(), mesh.n_faces()),
                geometry: (self.n_edges(), self.n_faces()),
            })
        }
    }
}

/// Area, corner angles and corner cotangents of one face from the lengths of its halfedges.
fn face_geometry(l: [f64; 3]) -> Result<(f64, [f64; 3], [f64; 3]), DegenerateReason> {
    sorted_valid(l[0], l[1], l[2])?;
    let area = triangle_area(l[0], l[1], l[2]).ok_or(DegenerateReason::OutOfRange)?;
    let mut angles = [0.0; 3];
    let mut cotans = [0.0; 3];
    for k in 0..3 {
        // Corner k lies between halfedges k (leaving it) and k + 2 (arriving); opposite is k + 1.
        let (opposite, s1, s2) = (l[(k + 1) % 3], l[k], l[(k + 2) % 3]);
        angles[k] = corner_angle(opposite, s1, s2).ok_or(DegenerateReason::OutOfRange)?;
        cotans[k] = corner_cotan(opposite, s1, s2).ok_or(DegenerateReason::OutOfRange)?;
    }
    Ok((area, angles, cotans))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn right_triangle_3_4_5() {
        assert_eq!(triangle_area(3.0, 4.0, 5.0), Some(6.0));
        assert_eq!(triangle_area(5.0, 3.0, 4.0), Some(6.0));
        let right = corner_angle(5.0, 3.0, 4.0).unwrap();
        assert!((right - PI / 2.0).abs() < 1e-15);
        let a = corner_angle(3.0, 4.0, 5.0).unwrap();
        let b = corner_angle(4.0, 5.0, 3.0).unwrap();
        assert!((a - (3.0f64 / 4.0).atan()).abs() < 1e-15);
        assert!((a + b + right - PI).abs() < 1e-15);
        assert_eq!(corner_cotan(5.0, 3.0, 4.0), Some(0.0));
        assert!((corner_cotan(3.0, 4.0, 5.0).unwrap() - 4.0 / 3.0).abs() < 1e-15);
        assert!((corner_cotan(4.0, 3.0, 5.0).unwrap() - 0.75).abs() < 1e-15);
    }

    #[test]
    fn equilateral_triangle() {
        let s = 2.5;
        let area = triangle_area(s, s, s).unwrap();
        assert!((area - 3f64.sqrt() / 4.0 * s * s).abs() < 1e-14);
        assert!((corner_angle(s, s, s).unwrap() - PI / 3.0).abs() < 1e-15);
        assert!((corner_cotan(s, s, s).unwrap() - 1.0 / 3f64.sqrt()).abs() < 1e-15);
    }

    #[test]
    fn needle_triangle_keeps_its_digits() {
        // Kahan's example class: two long sides and a tiny one. Compare against the exact
        // isosceles formulas.
        let (long, short) = (1.0, 1e-9);
        let area = triangle_area(long, long, short).unwrap();
        let exact_area = 0.5 * short * (long * long - 0.25 * short * short).sqrt();
        assert!((area / exact_area - 1.0).abs() < 1e-15);
        let apex = corner_angle(short, long, long).unwrap();
        let exact_apex = 2.0 * (0.5 * short / long).asin();
        assert!((apex / exact_apex - 1.0).abs() < 1e-15);
        let base = corner_angle(long, long, short).unwrap();
        assert!((apex + 2.0 * base - PI).abs() < 1e-15);
        let cot_apex = corner_cotan(short, long, long).unwrap();
        assert!((cot_apex * exact_apex.tan() - 1.0).abs() < 1e-14);
        let cot_base = corner_cotan(long, long, short).unwrap();
        assert!((cot_base / (0.5 * exact_apex).tan() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn invalid_triangles_give_none() {
        for (a, b, c) in [
            (1.0, 2.0, 3.0),
            (1.0, 2.0, 3.5),
            (0.0, 1.0, 1.0),
            (0.0, 0.0, 0.0),
            (-1.0, 1.0, 1.0),
            (f64::NAN, 1.0, 1.0),
            (f64::INFINITY, 1.0, 1.0),
            (1e-200, 1e-200, 1e-200),
            (1e200, 1e200, 1e200),
        ] {
            assert_eq!(triangle_area(a, b, c), None, "{a} {b} {c}");
            assert_eq!(corner_cotan(a, b, c), None, "{a} {b} {c}");
        }
        assert_eq!(corner_angle(1.0, 2.0, 3.0), None);
        assert_eq!(corner_angle(3.0, 1.0, 2.0), None);
        assert_eq!(corner_angle(f64::NAN, 1.0, 2.0), None);
        assert_eq!(sorted_valid(1.0, 2.0, 3.0), Err(DegenerateReason::ZeroArea));
        assert_eq!(sorted_valid(1.0, 1.0, 0.0), Err(DegenerateReason::ZeroArea));
        assert_eq!(sorted_valid(0.0, 0.0, 0.0), Err(DegenerateReason::ZeroArea));
        assert_eq!(sorted_valid(1.0, 2.0, 0.0), Err(DegenerateReason::TriangleInequality));
        assert_eq!(sorted_valid(1.0, 2.0, 3.5), Err(DegenerateReason::TriangleInequality));
        assert_eq!(sorted_valid(1.0, f64::NAN, 3.5), Err(DegenerateReason::OutOfRange));
    }

    fn two_triangles() -> TriMesh {
        // A unit square split along the diagonal 0-2.
        TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap()
    }

    const SQUARE: [[f64; 3]; 4] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]];

    #[test]
    fn geometry_of_a_square() {
        let mesh = two_triangles();
        let g = Geometry::from_positions(&mesh, &SQUARE).unwrap();
        g.check_mesh(&mesh).unwrap();
        assert_eq!(mesh.edges(), &[[0, 1], [0, 2], [0, 3], [1, 2], [2, 3]]);
        let d = 2f64.sqrt();
        assert_eq!(g.edge_lengths(), &[1.0, d, 1.0, 1.0, 1.0]);
        assert_eq!((g.n_edges(), g.n_faces()), (5, 2));
        for f in 0..2 {
            assert!((g.face_areas()[f] - 0.5).abs() < 1e-15);
        }
        assert!((g.total_area() - 1.0).abs() < 1e-15);
        // Face 0 = (0, 1, 2): right angle at vertex 1.
        let a = g.corner_angles()[0];
        assert!((a[0] - PI / 4.0).abs() < 1e-15 && (a[1] - PI / 2.0).abs() < 1e-15 && (a[2] - PI / 4.0).abs() < 1e-15);
        let c = g.corner_cotans()[0];
        assert!((c[0] - 1.0).abs() < 1e-15 && c[1].abs() < 1e-15 && (c[2] - 1.0).abs() < 1e-15);
        // Face 1 = (0, 2, 3): right angle at vertex 3.
        let a = g.corner_angles()[1];
        assert!((a[2] - PI / 2.0).abs() < 1e-15);

        // The intrinsic constructor gives the same result.
        let g2 = Geometry::from_edge_lengths(&mesh, g.edge_lengths().to_vec()).unwrap();
        assert_eq!(g, g2);
    }

    #[test]
    fn error_counts_and_non_finite_input() {
        let mesh = two_triangles();
        let e = Geometry::from_positions(&mesh, &SQUARE[..3]).unwrap_err();
        assert_eq!(e.to_string(), "the mesh has 4 vertices but 3 positions were given");
        let mut p = SQUARE;
        p[2][1] = f64::NAN;
        let e = Geometry::from_positions(&mesh, &p).unwrap_err();
        assert_eq!(e.to_string(), "vertex 2 has a non-finite position (1, NaN, 0)");
        let e = Geometry::from_edge_lengths(&mesh, vec![1.0; 4]).unwrap_err();
        assert_eq!(e.to_string(), "the mesh has 5 edges but 4 edge lengths were given");
        let e = Geometry::from_edge_lengths(&mesh, vec![1.0, 1.0, -2.0, 1.0, f64::INFINITY]).unwrap_err();
        assert_eq!(
            e.to_string(),
            "edge 2 = (0, 3) has invalid length -2 (2 edges have a negative or non-finite length)"
        );
        let other = TriMesh::new(3, vec![[0, 1, 2]]).unwrap();
        let g = Geometry::from_positions(&mesh, &SQUARE).unwrap();
        assert_eq!(
            g.check_mesh(&other).unwrap_err().to_string(),
            "geometry does not belong to this mesh: the mesh has 3 edges and 1 faces, the geometry 5 edges and 2 faces"
        );
    }

    #[test]
    fn degenerate_faces_are_named() {
        let mesh = two_triangles();
        // Vertex 1 moved onto the diagonal: face 0 is collinear.
        let mut p = SQUARE;
        p[1] = [0.5, 0.5, 0.0];
        let e = Geometry::from_positions(&mesh, &p).unwrap_err();
        match &e {
            GeometryError::DegenerateFaces { faces } => {
                assert_eq!(faces.len(), 1);
                assert_eq!(faces[0].face, 0);
                assert_eq!(faces[0].vertices, [0, 1, 2]);
                assert_eq!(faces[0].reason, DegenerateReason::ZeroArea);
            }
            other => panic!("unexpected {other:?}"),
        }
        let h = 0.5f64.sqrt();
        assert_eq!(
            e.to_string(),
            format!(
                "1 degenerate face: face 0 (vertices 0, 1, 2; edge lengths {h}, {h}, {}) has zero area",
                2f64.sqrt()
            )
        );

        // Coincident vertices: both faces collapse.
        let mut p = SQUARE;
        p[2] = p[0];
        let e = Geometry::from_positions(&mesh, &p).unwrap_err();
        assert_eq!(
            e.to_string(),
            "2 degenerate faces: face 0 (vertices 0, 1, 2; edge lengths 1, 1, 0) has zero area; \
             face 1 (vertices 0, 2, 3; edge lengths 0, 1, 1) has zero area"
        );

        // Intrinsic lengths that violate the triangle inequality in face 1 only.
        let e = Geometry::from_edge_lengths(&mesh, vec![1.0, 1.5, 1.0, 1.0, 3.0]).unwrap_err();
        assert_eq!(
            e.to_string(),
            "1 degenerate face: face 1 (vertices 0, 2, 3; edge lengths 1.5, 3, 1) violates the triangle inequality"
        );
    }
}
