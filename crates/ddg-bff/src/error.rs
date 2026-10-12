//! The error type of this crate.

use std::fmt;

use ddg_core::geometry::GeometryError;
use ddg_core::sparse::SparseError;
use ddg_core::topology::MeshError;

/// Why a cut could not be planned or a mesh could not be flattened.
///
/// Every variant names the offending element (vertex, edge, face or chart) so that the message
/// can be shown to a user as is.
#[derive(Debug, Clone, PartialEq)]
pub enum BffError {
    /// A mesh built by this crate (a cut mesh) failed validation. This indicates a bug in the
    /// cutting code, not bad input; it is reported instead of panicking.
    Mesh(MeshError),
    /// The geometry does not belong to the mesh, or is otherwise unusable.
    Geometry(GeometryError),
    /// A sparse factorization or solve failed.
    Solve {
        /// Which linear system.
        system: &'static str,
        /// What the solver reported.
        error: SparseError,
    },
    /// A per-element input array has the wrong length.
    LengthMismatch {
        /// Which argument.
        what: &'static str,
        /// The required length.
        expected: usize,
        /// The length that was passed.
        got: usize,
    },
    /// An input value is NaN or infinite.
    NonFiniteInput {
        /// Which argument.
        what: &'static str,
        /// Index of the first offending entry.
        index: usize,
        /// The value.
        value: f64,
    },
    /// An edge id is not an edge of the mesh.
    EdgeOutOfRange {
        /// Which argument listed it.
        what: &'static str,
        /// Position in that list.
        position: usize,
        /// The offending id.
        edge: u32,
        /// Number of edges of the mesh.
        n_edges: usize,
    },
    /// A cut consists of a single edge between two interior vertices that touch no other cut
    /// edge. Such a slit cannot be represented by an indexed triangle mesh, because its two
    /// sides would have the same two end points; extend the cut by one more edge.
    IsolatedCutEdge {
        /// The edge.
        edge: u32,
        /// Its vertices `[min, max]`.
        vertices: [u32; 2],
    },
    /// A cone is placed on a vertex that does not exist.
    ConeVertexOutOfRange {
        /// Position in the cone list.
        position: usize,
        /// The offending vertex index.
        vertex: u32,
        /// Number of vertices of the mesh.
        n_vertices: usize,
    },
    /// Two cones are placed on the same vertex.
    DuplicateCone {
        /// The vertex.
        vertex: u32,
        /// Position of its first occurrence in the cone list.
        first: usize,
        /// Position of its second occurrence.
        second: usize,
    },
    /// A cone is placed on a boundary vertex. Cone angles are prescribed at interior vertices
    /// only; on the boundary the scale factor is fixed to zero and the boundary curvature
    /// follows from it.
    ConeOnBoundary {
        /// Position in the cone list.
        position: usize,
        /// The boundary vertex.
        vertex: u32,
    },
    /// A cone angle is not finite, or is `>= 2π` (the angle sum `2π − µ` around the cone would
    /// not be positive).
    InvalidConeAngle {
        /// Position in the cone list.
        position: usize,
        /// The vertex.
        vertex: u32,
        /// The angle.
        angle: f64,
    },
    /// The cone angles on a closed connected component do not sum to `2πχ`.
    GaussBonnetViolated {
        /// The component (numbered as in `TriMesh::face_components`).
        component: u32,
        /// Lowest-numbered face of the component.
        first_face: u32,
        /// Euler characteristic of the component.
        euler_characteristic: i64,
        /// Number of cones on the component.
        n_cones: usize,
        /// Sum of its cone angles.
        cone_sum: f64,
        /// `2πχ`.
        expected: f64,
        /// The tolerance that was exceeded.
        tolerance: f64,
    },
    /// A closed component has a discrepancy to distribute over its cones but no cone with a
    /// non-zero angle.
    NoConesToAdjust {
        /// The component.
        component: u32,
        /// `Σµ − 2πχ`.
        discrepancy: f64,
    },
    /// A closed connected component of genus 0 needs at least two cones: its cone angles sum to
    /// `4π` while each is below `2π`, and a cut needs two end points.
    SphereNeedsTwoCones {
        /// The component of the mesh after the required cuts were made.
        component: u32,
        /// Lowest-numbered face of the component.
        first_face: u32,
        /// Number of cones on it.
        n_cones: usize,
    },
    /// An edge is listed as required and is also marked as not cuttable.
    RequiredEdgeForbidden {
        /// The edge.
        edge: u32,
        /// Its vertices `[min, max]`.
        vertices: [u32; 2],
    },
    /// A cone lies at a vertex all of whose edges are marked as not cuttable, so no cut can
    /// reach it.
    ConeEnclosedByUncuttableEdges {
        /// The cone vertex.
        vertex: u32,
    },
    /// The edges that may not be cut hold together a set of faces that is not a topological
    /// disk, so no admissible cut makes the surface a disk.
    UncuttableRegionNotDisk {
        /// A face of the region (its lowest-numbered one).
        first_face: u32,
        /// Number of faces of the region.
        n_faces: usize,
        /// Its Euler characteristic once every cuttable edge is cut.
        euler_characteristic: i64,
        /// Its number of boundary loops once every cuttable edge is cut.
        boundary_loops: usize,
    },
    /// A chart of the cut surface is not a topological disk.
    ChartNotDisk {
        /// The chart (numbered by its lowest face).
        chart: u32,
        /// Lowest-numbered face of the chart.
        first_face: u32,
        /// Number of faces of the chart.
        n_faces: usize,
        /// Its Euler characteristic (a disk has 1).
        euler_characteristic: i64,
        /// Its number of boundary loops (a disk has 1).
        boundary_loops: usize,
        /// Its genus (a disk has 0).
        genus: usize,
    },
    /// After cutting, a cone is still an interior vertex of its chart.
    ConeNotOnCut {
        /// The cone vertex.
        vertex: u32,
        /// The chart that contains it.
        chart: u32,
    },
    /// The mesh passed to a disk flattening is not a topological disk.
    NotADisk {
        /// Number of connected components (a disk has 1).
        components: usize,
        /// Euler characteristic (a disk has 1).
        euler_characteristic: i64,
        /// Number of boundary loops (a disk has 1).
        boundary_loops: usize,
    },
    /// An anchor is not a corner of the mesh.
    AnchorOutOfRange {
        /// Position in the anchor list.
        position: usize,
        /// The offending corner id.
        corner: u32,
        /// Number of corners of the mesh (three per face).
        n_corners: usize,
    },
    /// An anchor corner lies at a vertex that is not on the boundary of its chart, so no
    /// boundary curve can start there.
    AnchorNotOnBoundary {
        /// The corner.
        corner: u32,
        /// Its vertex.
        vertex: u32,
    },
    /// Two anchors lie in the same chart.
    DuplicateAnchor {
        /// The chart.
        chart: u32,
        /// The first anchor corner in that chart.
        first: u32,
        /// The second one.
        second: u32,
    },
    /// The vertex at which the boundary of a disk should start is not a boundary vertex.
    StartNotOnBoundary {
        /// The vertex.
        vertex: u32,
    },
    /// Prescribed exterior angles of a disk do not sum to `2π`.
    ExteriorAnglesDoNotClose {
        /// Their sum.
        sum: f64,
        /// The tolerance on `|sum − 2π|` that was exceeded.
        tolerance: f64,
    },
    /// Scale factors passed in by the caller are not flat away from the listed cones: they
    /// imply curvature at a vertex that is interior to its chart.
    ScaleFactorsNotFlat {
        /// The vertex with the largest implied curvature among the chart-interior vertices.
        vertex: u32,
        /// The implied cone angle `Ω + (L u)` there.
        curvature: f64,
        /// The tolerance that was exceeded.
        tolerance: f64,
    },
    /// A computed quantity is not finite. Reported instead of returning NaN.
    NonFiniteResult {
        /// Which quantity.
        what: &'static str,
        /// The chart, if the quantity belongs to one.
        chart: Option<u32>,
    },
}

impl fmt::Display for BffError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BffError::Mesh(e) => write!(f, "internal error: the cut mesh is invalid: {e}"),
            BffError::Geometry(e) => write!(f, "{e}"),
            BffError::Solve { system, error } => write!(f, "solving {system} failed: {error}"),
            BffError::LengthMismatch {
                what,
                expected,
                got,
            } => write!(f, "{what} has {got} entries, expected {expected}"),
            BffError::NonFiniteInput { what, index, value } => {
                write!(f, "{what}[{index}] = {value} is not finite")
            }
            BffError::EdgeOutOfRange {
                what,
                position,
                edge,
                n_edges,
            } => write!(
                f,
                "{what}[{position}] is edge {edge}, but the mesh has {n_edges} edges"
            ),
            BffError::IsolatedCutEdge {
                edge,
                vertices: [a, b],
            } => write!(
                f,
                "cut edge {edge} = ({a}, {b}) is an isolated slit: neither end point is on the \
                 boundary or on another cut edge, so the two sides of the cut cannot be told \
                 apart; extend the cut by one more edge"
            ),
            BffError::ConeVertexOutOfRange {
                position,
                vertex,
                n_vertices,
            } => write!(
                f,
                "cone {position} is at vertex {vertex}, but the mesh has {n_vertices} vertices"
            ),
            BffError::DuplicateCone {
                vertex,
                first,
                second,
            } => write!(f, "cones {first} and {second} are both at vertex {vertex}"),
            BffError::ConeOnBoundary { position, vertex } => write!(
                f,
                "cone {position} is at boundary vertex {vertex}: cone angles can only be \
                 prescribed at interior vertices"
            ),
            BffError::InvalidConeAngle {
                position,
                vertex,
                angle,
            } => write!(
                f,
                "cone {position} at vertex {vertex} has angle {angle}: a cone angle must be \
                 finite and below 2π"
            ),
            BffError::GaussBonnetViolated {
                component,
                first_face,
                euler_characteristic,
                n_cones,
                cone_sum,
                expected,
                tolerance,
            } => write!(
                f,
                "the {n_cones} cone angles on closed component {component} (containing face \
                 {first_face}, Euler characteristic {euler_characteristic}) sum to {cone_sum}, \
                 but Gauss-Bonnet requires 2πχ = {expected}: discrepancy {:e} exceeds the \
                 tolerance {tolerance:e}",
                cone_sum - expected
            ),
            BffError::NoConesToAdjust {
                component,
                discrepancy,
            } => write!(
                f,
                "closed component {component} violates Gauss-Bonnet by {discrepancy:e}, and it \
                 has no cone with a non-zero angle over which the discrepancy could be distributed"
            ),
            BffError::SphereNeedsTwoCones {
                component,
                first_face,
                n_cones,
            } => write!(
                f,
                "component {component} (containing face {first_face}) is a closed surface of \
                 genus 0 with {n_cones} cone{}: at least two cones are needed to cut it into a disk",
                if *n_cones == 1 { "" } else { "s" }
            ),
            BffError::RequiredEdgeForbidden {
                edge,
                vertices: [a, b],
            } => write!(
                f,
                "edge {edge} = ({a}, {b}) is required to be cut and also marked as not cuttable"
            ),
            BffError::ConeEnclosedByUncuttableEdges { vertex } => write!(
                f,
                "the cone at vertex {vertex} cannot be reached by a cut: every edge at this \
                 vertex is marked as not cuttable"
            ),
            BffError::UncuttableRegionNotDisk {
                first_face,
                n_faces,
                euler_characteristic,
                boundary_loops,
            } => write!(
                f,
                "the edges that may not be cut hold together {n_faces} faces (among them face \
                 {first_face}) that do not form a disk: even with every other edge cut, the \
                 region has Euler characteristic {euler_characteristic} and {boundary_loops} \
                 boundary loops"
            ),
            BffError::ChartNotDisk {
                chart,
                first_face,
                n_faces,
                euler_characteristic,
                boundary_loops,
                genus,
            } => write!(
                f,
                "chart {chart} ({n_faces} faces, among them face {first_face}) is not a disk \
                 after cutting: Euler characteristic {euler_characteristic}, {boundary_loops} \
                 boundary loops, genus {genus}"
            ),
            BffError::ConeNotOnCut { vertex, chart } => write!(
                f,
                "the cone at vertex {vertex} is not on the cut: it is an interior vertex of \
                 chart {chart}"
            ),
            BffError::NotADisk {
                components,
                euler_characteristic,
                boundary_loops,
            } => write!(
                f,
                "the mesh is not a topological disk: {components} connected component{}, Euler \
                 characteristic {euler_characteristic}, {boundary_loops} boundary loop{} \
                 (a disk has 1, 1 and 1)",
                if *components == 1 { "" } else { "s" },
                if *boundary_loops == 1 { "" } else { "s" }
            ),
            BffError::AnchorOutOfRange {
                position,
                corner,
                n_corners,
            } => write!(
                f,
                "anchor {position} is corner {corner}, but the mesh has {n_corners} corners"
            ),
            BffError::AnchorNotOnBoundary { corner, vertex } => write!(
                f,
                "anchor corner {corner} is at vertex {vertex}, which is not on the boundary of \
                 its chart after cutting"
            ),
            BffError::DuplicateAnchor {
                chart,
                first,
                second,
            } => write!(
                f,
                "anchor corners {first} and {second} are both in chart {chart}; a chart can \
                 have one anchor"
            ),
            BffError::StartNotOnBoundary { vertex } => write!(
                f,
                "vertex {vertex} is not on the boundary, so the boundary cannot start there"
            ),
            BffError::ExteriorAnglesDoNotClose { sum, tolerance } => write!(
                f,
                "the exterior angles sum to {sum}, not to 2π (difference {:e}, tolerance \
                 {tolerance:e})",
                sum - 2.0 * std::f64::consts::PI
            ),
            BffError::ScaleFactorsNotFlat {
                vertex,
                curvature,
                tolerance,
            } => write!(
                f,
                "the scale factors imply a cone of angle {curvature:e} at vertex {vertex}, which \
                 is not on the cut (tolerance {tolerance:e}): list it as a cone or correct the \
                 scale factors"
            ),
            BffError::NonFiniteResult { what, chart } => match chart {
                Some(c) => write!(f, "{what} of chart {c} is not finite"),
                None => write!(f, "{what} is not finite"),
            },
        }
    }
}

impl std::error::Error for BffError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            BffError::Mesh(e) => Some(e),
            BffError::Geometry(e) => Some(e),
            BffError::Solve { error, .. } => Some(error),
            _ => None,
        }
    }
}

impl From<MeshError> for BffError {
    fn from(e: MeshError) -> Self {
        BffError::Mesh(e)
    }
}

impl From<GeometryError> for BffError {
    fn from(e: GeometryError) -> Self {
        BffError::Geometry(e)
    }
}

/// Checks that `values` has `expected` entries, all finite.
pub(crate) fn check_values(
    what: &'static str,
    values: &[f64],
    expected: usize,
) -> Result<(), BffError> {
    if values.len() != expected {
        return Err(BffError::LengthMismatch {
            what,
            expected,
            got: values.len(),
        });
    }
    if let Some(index) = values.iter().position(|x| !x.is_finite()) {
        return Err(BffError::NonFiniteInput {
            what,
            index,
            value: values[index],
        });
    }
    Ok(())
}
