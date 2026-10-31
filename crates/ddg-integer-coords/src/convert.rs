//! Reading integer coordinates off a signpost triangulation (`ddg_intrinsic`): the normal
//! coordinate of every edge from its traced path over the input mesh, the roundabouts from the
//! input faces in which the paths leave their ends.
//!
//! **Floating point decides this conversion.** Each edge is traced over the input mesh with
//! the signposts ([`IntrinsicTriangulation::trace_edges`]); the number of input edges the path
//! crosses is `n`, and the input face in which it leaves an input vertex gives the roundabout.
//! A trace that runs within rounding of an input vertex or edge can get either wrong. The
//! result is therefore certified ([`IntegerTriangulation::certify`]) before it is returned, and
//! a conversion that does not pass is an error, not a repaired state.
//!
//! Two readings are exact rather than metric, and are applied before certification:
//!
//! * a path with no crossing whose two ends lie on a common input edge (two of its end points,
//!   or a point of it) runs along that input edge: a straight segment inside one flat triangle
//!   between two points of the same side is that side (`n = −1`);
//! * a last crossing with the very input edge on which the far vertex lies is the arrival at
//!   that vertex, not a crossing (a straight segment inside a triangle cannot meet one of its
//!   sides twice).

use std::sync::Arc;

use ddg_core::TriMesh;
use ddg_intrinsic::{EdgeTraceOptions, IntrinsicTriangulation, SurfacePoint};

use crate::certify::CertifyReport;
use crate::coords::{IntegerTriangulation, LocalIndices, NO_ROUNDABOUT, Parts};
use crate::error::IntegerError;

/// What [`IntegerTriangulation::from_intrinsic`] did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConversionReport {
    /// Number of edges.
    pub n_edges: usize,
    /// Edges found to run along an input edge (`n = −1`), by the trace or because they are
    /// known to.
    pub n_along: usize,
    /// Of those, edges whose trace had no crossing and stayed in one input face between two
    /// points of one of its sides, read as running along that side.
    pub n_along_by_endpoints: usize,
    /// Last crossings dropped because they are the arrival at a far vertex on that input edge.
    pub n_arrival_crossings_dropped: usize,
    /// The largest mismatch of a trace at its far vertex, relative to the edge's length.
    pub max_relative_mismatch: f64,
    /// The certificate of the result.
    pub certificate: CertifyReport,
}

/// The input edges of input face `f` on which point `q` lies (an end point or a point of it).
fn sides_through(mesh: &TriMesh, f: usize, q: SurfacePoint) -> Vec<usize> {
    let sides = [0, 1, 2].map(|k| mesh.edge(3 * f + k));
    match q {
        SurfacePoint::Vertex(v) => sides
            .into_iter()
            .filter(|&e| mesh.edge_vertices(e).contains(&(v as usize)))
            .collect(),
        SurfacePoint::Edge { edge, .. } => {
            sides.into_iter().filter(|&e| e == edge as usize).collect()
        }
        SurfacePoint::Face { .. } => Vec::new(),
    }
}

/// The local index of the first input halfedge counter-clockwise after a direction that leaves
/// input vertex `a` strictly inside the corner of input face `f` at `a`.
fn roundabout_in_face(
    mesh: &TriMesh,
    local: &LocalIndices,
    a: usize,
    f: usize,
) -> Result<u32, String> {
    let k = (0..3)
        .find(|&k| mesh.faces()[f][k] as usize == a)
        .ok_or_else(|| format!("input face {f} has no corner at vertex {a}"))?;
    let i = local.of_halfedge(3 * f + k) + 1;
    let deg = local.degree(a);
    Ok(if local.is_boundary(a) { i } else { i % deg })
}

impl IntegerTriangulation {
    /// The integer coordinates of a signpost triangulation (SPEC WP1 item 5): same
    /// connectivity, lengths and vertex places; `n` of an edge from its traced path over the
    /// input mesh (the number of input edges it crosses, `−1` along an input edge), the
    /// roundabouts from where the paths leave their input-vertex ends. See the
    /// [module documentation](crate::convert) for what floating point decides.
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] when a trace fails its check
    /// ([`IntrinsicTriangulation::trace_edges`]), [`IntegerError::Conversion`] naming the edge
    /// when a trace cannot be read, and the error of [`certify`](Self::certify) when the
    /// integers read off the traces are not valid.
    pub fn from_intrinsic(
        signposts: &IntrinsicTriangulation,
        options: &EdgeTraceOptions,
    ) -> Result<(IntegerTriangulation, ConversionReport), IntegerError> {
        let input = Arc::new(signposts.input().clone());
        let mesh = input.mesh();
        let local = LocalIndices::new(mesh);
        let complex = signposts.complex().clone();
        let traces = signposts.trace_edges(options)?;
        let locations = signposts.vertex_locations().to_vec();
        let n_input = mesh.n_vertices();
        let mut normal = vec![0i64; complex.n_edges()];
        let mut roundabouts = vec![NO_ROUNDABOUT; complex.n_halfedges()];
        let mut report = ConversionReport {
            n_edges: complex.n_edges(),
            n_along: 0,
            n_along_by_endpoints: 0,
            n_arrival_crossings_dropped: 0,
            max_relative_mismatch: traces.report.max_relative_mismatch,
            certificate: CertifyReport {
                n_faces: 0,
                n_along: 0,
                n_crossings: 0,
                n_input_edges: 0,
                n_inserted_on_edges: 0,
                max_crossings_per_input_edge: 0,
            },
        };
        for (e, n_e) in normal.iter_mut().enumerate() {
            let fail = |reason: String| IntegerError::Conversion {
                edge: e as u32,
                reason,
            };
            let h = complex.edge_halfedge(e);
            let (a, b) = (complex.from_vertex(h), complex.to_vertex(h));
            let mut crossings = traces.crossings(e);
            // The arrival at a far vertex on an input edge is not a crossing.
            if let (SurfacePoint::Edge { edge, .. }, Some(last)) = (locations[b], crossings.last())
                && last.edge == edge
            {
                crossings = &crossings[..crossings.len() - 1];
                report.n_arrival_crossings_dropped += 1;
            }
            let mut along = traces.along_input_edge(e).map(|x| x as usize);
            let first_face = traces.first_face(e).map(|f| f as usize);
            if along.is_none()
                && crossings.is_empty()
                && let Some(f) = first_face
            {
                let at_a = sides_through(mesh, f, locations[a]);
                let at_b = sides_through(mesh, f, locations[b]);
                if let Some(&side) = at_a.iter().find(|s| at_b.contains(s)) {
                    along = Some(side);
                    report.n_along_by_endpoints += 1;
                }
            }
            // The input face in which the path arrives at b.
            let last_face = match crossings.last() {
                Some(c) => Some(
                    mesh.twin(c.halfedge as usize)
                        .map(|t| t / 3)
                        .ok_or_else(|| fail("a crossing with a boundary input edge".into()))?,
                ),
                None => first_face,
            };
            match along {
                Some(ie) => {
                    *n_e = -1;
                    report.n_along += 1;
                    for (v, g) in [(a, Some(h)), (b, complex.twin(h))] {
                        if let (true, Some(g)) = (v < n_input, g) {
                            roundabouts[g] = local.of_edge_at(mesh, ie, v).ok_or_else(|| {
                                fail(format!(
                                    "it runs along input edge {ie}, which does not end at its \
                                     end point {v}"
                                ))
                            })?;
                        }
                    }
                }
                None => {
                    *n_e = crossings.len() as i64;
                    let (Some(f0), Some(f1)) = (first_face, last_face) else {
                        return Err(fail("the trace has no first face".into()));
                    };
                    if a < n_input {
                        roundabouts[h] = roundabout_in_face(mesh, &local, a, f0).map_err(fail)?;
                    }
                    if let (true, Some(t)) = (b < n_input, complex.twin(h)) {
                        roundabouts[t] = roundabout_in_face(mesh, &local, b, f1).map_err(fail)?;
                    }
                }
            }
        }
        let result = IntegerTriangulation::from_parts(Parts {
            input,
            complex,
            lengths: signposts.edge_lengths().to_vec(),
            normal,
            roundabouts,
            locations,
        })?;
        report.certificate = result.certify()?;
        Ok((result, report))
    }
}
