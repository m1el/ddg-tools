//! Tracing the edges of the intrinsic triangulation over the input mesh (paper Sec. 3.4.2), with
//! the check that each trace arrives where its edge ends.
//!
//! An intrinsic edge is traced from the place of its tail on the input mesh, in the direction
//! of its signpost, for its length. While the triangles of the input mesh are unfolded along
//! the path (see [`crate::trace`]), the place of the edge's head is unfolded with every
//! triangle that contains it. The trace **succeeds** when that place comes to lie within the
//! tolerance of the end of the path; the path then ends there and the crossings up to that
//! triangle are the result. It **fails** when the path has been followed for the edge's whole
//! length (and the tolerance beyond it) without this happening, or leaves the surface, or runs
//! into another vertex. A failure is an error that names the edge; nothing is perturbed,
//! snapped or retried (SPEC N3).
//!
//! The paper's Appendix A describes this check as a topological one (the path must end with
//! a crossing into the star of the simplex that contains the target) followed by repairs:
//! crossings after that one are dropped, and a path that stops short of the star is connected
//! through a neighbouring triangle. Here the first is not a repair but the stopping rule
//! (the path ends when the target is met, so there are no later crossings), and the second
//! has its counterpart in following the path for the tolerance beyond the edge's length.
//! Both are bounded by a stated metric tolerance, which the paper's rule is not.

use crate::error::{EdgeTraceFailure, EdgeTraceFailureReason, IntrinsicError};
use crate::point::SurfacePoint;
use crate::trace::{Crossing, Target, TraceOptions, TraceStop, run};
use crate::triangulation::IntrinsicTriangulation;

/// Options of [`IntrinsicTriangulation::trace_edge`] and
/// [`IntrinsicTriangulation::trace_edges`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeTraceOptions {
    /// How far the unfolded place of the edge's far vertex may be from the end of the traced
    /// path, as a fraction of the edge's length. Finite and not negative.
    pub tolerance: f64,
    /// The limit on the number of input edges one trace may cross; see
    /// [`TraceOptions::max_crossings`].
    pub max_crossings: Option<usize>,
}

impl Default for EdgeTraceOptions {
    /// A tolerance of `1e-6` of the edge's length and the default crossing limit. On meshes
    /// without nearly degenerate triangles the mismatch is of the order of `1e-13`; the
    /// tolerance leaves room for the loss of accuracy in needle and cap triangles (the
    /// crate's README has the measured values).
    fn default() -> Self {
        EdgeTraceOptions {
            tolerance: 1e-6,
            max_crossings: None,
        }
    }
}

/// The path of one intrinsic edge over the input mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeTrace {
    /// The intrinsic edge.
    pub edge: u32,
    /// Its end points, in the direction of the trace (tail and head of its halfedge).
    pub vertices: [u32; 2],
    /// The input edges it crosses, in order from its tail. `halfedge` and `edge` are ids of
    /// the input mesh, `t` runs from the lower-numbered end of the input edge, and `distance`
    /// is measured along the intrinsic edge from its tail.
    pub crossings: Vec<Crossing>,
    /// The input edge the intrinsic edge lies along, if it does; it then has no crossings.
    pub along_input_edge: Option<u32>,
    /// The input face in which the path leaves the tail; `None` if it lies along an input
    /// edge.
    pub first_face: Option<u32>,
    /// The distance between the unfolded place of the far vertex and the end of the traced
    /// path (an absolute length). Zero for an edge that was not traced because it is known to
    /// lie along an input edge ([`IntrinsicTriangulation::edge_on_input`]).
    pub mismatch: f64,
    /// See [`Trace::clamp`](crate::Trace::clamp).
    pub clamp: f64,
}

/// What tracing all edges found.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeTraceReport {
    /// Number of edges of the triangulation.
    pub n_edges: usize,
    /// Number of edges that were not traced because they are known to lie along an input
    /// edge.
    pub n_on_input: usize,
    /// Number of edges that were traced.
    pub n_traced: usize,
    /// Total number of crossings with input edges.
    pub n_crossings: usize,
    /// The largest number of crossings of one edge.
    pub max_crossings: usize,
    /// The largest mismatch of a traced edge divided by its length.
    pub max_relative_mismatch: f64,
    /// The largest [`EdgeTrace::clamp`].
    pub max_clamp: f64,
}

/// The paths of all intrinsic edges over the input mesh. Returned by
/// [`IntrinsicTriangulation::trace_edges`].
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeTraces {
    offsets: Vec<u32>,
    crossings: Vec<Crossing>,
    along: Vec<Option<u32>>,
    first_faces: Vec<Option<u32>>,
    mismatches: Vec<f64>,
    /// Counts and the largest mismatch.
    pub report: EdgeTraceReport,
}

impl EdgeTraces {
    /// Number of edges.
    pub fn n_edges(&self) -> usize {
        self.along.len()
    }

    /// The crossings of intrinsic edge `e` with input edges, in order from the tail of
    /// [`edge_halfedge`](IntrinsicTriangulation::edge_halfedge); see [`EdgeTrace::crossings`].
    pub fn crossings(&self, e: usize) -> &[Crossing] {
        &self.crossings[self.offsets[e] as usize..self.offsets[e + 1] as usize]
    }

    /// The input edge that intrinsic edge `e` lies along, if it does.
    pub fn along_input_edge(&self, e: usize) -> Option<u32> {
        self.along[e]
    }

    /// The input face in which the path of edge `e` leaves its tail; `None` if the edge lies
    /// along an input edge.
    pub fn first_face(&self, e: usize) -> Option<u32> {
        self.first_faces[e]
    }

    /// The mismatch of edge `e` at its far vertex (an absolute length); see
    /// [`EdgeTrace::mismatch`].
    pub fn mismatch(&self, e: usize) -> f64 {
        self.mismatches[e]
    }
}

fn check_options(options: &EdgeTraceOptions) -> Result<(), IntrinsicError> {
    if options.tolerance.is_finite() && options.tolerance >= 0.0 {
        Ok(())
    } else {
        Err(IntrinsicError::InvalidOptions {
            reason: format!(
                "the tolerance of an edge trace must be finite and not negative, got {}",
                options.tolerance
            ),
        })
    }
}

impl IntrinsicTriangulation {
    /// The trace of edge `e`, which is in range, or how it failed.
    fn trace_edge_unchecked(
        &self,
        e: usize,
        options: &EdgeTraceOptions,
    ) -> Result<EdgeTrace, EdgeTraceFailure> {
        let h = self.complex.edge_halfedge[e] as usize;
        let (a, b) = (self.from_vertex(h), self.to_vertex(h));
        let vertices = [a as u32, b as u32];
        if let Some(input_edge) = self.edge_on_input(e) {
            return Ok(EdgeTrace {
                edge: e as u32,
                vertices,
                crossings: Vec::new(),
                along_input_edge: Some(input_edge as u32),
                first_face: None,
                mismatch: 0.0,
                clamp: 0.0,
            });
        }
        let length = self.edge_lengths[e];
        let tolerance = options.tolerance * length;
        let failure = |reason| EdgeTraceFailure {
            edge: e as u32,
            vertices,
            length,
            tolerance,
            reason,
        };
        let trace_options = TraceOptions {
            record_crossings: true,
            max_crossings: options.max_crossings,
        };
        let target = Target {
            point: self.locations[b],
            tolerance,
        };
        let outcome = run(
            &*self.input,
            self.locations[a],
            self.direction_on_input(a, self.signposts[h]),
            length,
            &trace_options,
            Some(target),
        )
        .map_err(|error| {
            failure(match error {
                IntrinsicError::CrossingLimit {
                    limit, travelled, ..
                } => EdgeTraceFailureReason::CrossingLimit { limit, travelled },
                other => EdgeTraceFailureReason::Query(Box::new(other)),
            })
        })?;
        let trace = outcome.trace;
        match outcome.mismatch {
            Some(mismatch) => Ok(EdgeTrace {
                edge: e as u32,
                vertices,
                crossings: trace.crossings,
                along_input_edge: trace.along_edge,
                first_face: trace.first_face,
                mismatch,
                clamp: trace.clamp,
            }),
            None => Err(failure(match (trace.stop, trace.end.as_vertex()) {
                (TraceStop::Boundary, _) => EdgeTraceFailureReason::HitBoundary {
                    at: trace.end,
                    travelled: trace.length,
                },
                // Arriving at the far vertex itself with length left over is a mismatch of
                // lengths, reported as such.
                (TraceStop::Vertex, Some(vertex))
                    if self.locations[b] != SurfacePoint::Vertex(vertex) =>
                {
                    EdgeTraceFailureReason::HitVertex {
                        vertex,
                        travelled: trace.length,
                    }
                }
                _ => EdgeTraceFailureReason::Missed {
                    end: trace.end,
                    mismatch: outcome.closest,
                    n_crossings: trace.n_crossings,
                },
            })),
        }
    }

    /// Traces edge `e` over the input mesh, from the tail of
    /// [`edge_halfedge`](Self::edge_halfedge) to its head, and checks that the trace arrives
    /// at the head (see the [module documentation](self) for the check).
    ///
    /// An edge that is known to lie along an input edge
    /// ([`edge_on_input`](Self::edge_on_input)) is not traced: it has no crossings.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`], [`IntrinsicError::InvalidOptions`], and
    /// [`IntrinsicError::EdgeTrace`] when the trace fails its check.
    pub fn trace_edge(
        &self,
        e: usize,
        options: &EdgeTraceOptions,
    ) -> Result<EdgeTrace, IntrinsicError> {
        self.check_edge(e)?;
        check_options(options)?;
        self.trace_edge_unchecked(e, options)
            .map_err(IntrinsicError::EdgeTrace)
    }

    /// Traces every edge over the input mesh ([`trace_edge`](Self::trace_edge)).
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidOptions`], and [`IntrinsicError::EdgeTraces`] listing every
    /// edge whose trace failed its check. No partial result is returned.
    pub fn trace_edges(&self, options: &EdgeTraceOptions) -> Result<EdgeTraces, IntrinsicError> {
        check_options(options)?;
        let n = self.n_edges();
        let mut offsets = Vec::with_capacity(n + 1);
        let mut crossings: Vec<Crossing> = Vec::new();
        let mut along = Vec::with_capacity(n);
        let mut first_faces = Vec::with_capacity(n);
        let mut mismatches = Vec::with_capacity(n);
        let mut failures = Vec::new();
        let mut report = EdgeTraceReport {
            n_edges: n,
            n_on_input: 0,
            n_traced: 0,
            n_crossings: 0,
            max_crossings: 0,
            max_relative_mismatch: 0.0,
            max_clamp: 0.0,
        };
        offsets.push(0);
        for e in 0..n {
            match self.trace_edge_unchecked(e, options) {
                Ok(trace) => {
                    if self.edge_on_input(e).is_some() {
                        report.n_on_input += 1;
                    } else {
                        report.n_traced += 1;
                        report.max_relative_mismatch = report
                            .max_relative_mismatch
                            .max(trace.mismatch / self.edge_lengths[e]);
                    }
                    report.n_crossings += trace.crossings.len();
                    report.max_crossings = report.max_crossings.max(trace.crossings.len());
                    report.max_clamp = report.max_clamp.max(trace.clamp);
                    crossings.extend_from_slice(&trace.crossings);
                    along.push(trace.along_input_edge);
                    first_faces.push(trace.first_face);
                    mismatches.push(trace.mismatch);
                }
                Err(failure) => {
                    failures.push(failure);
                    along.push(None);
                    first_faces.push(None);
                    mismatches.push(f64::INFINITY);
                }
            }
            offsets.push(crossings.len() as u32);
        }
        if !failures.is_empty() {
            return Err(IntrinsicError::EdgeTraces {
                n_edges: n,
                failures,
            });
        }
        Ok(EdgeTraces {
            offsets,
            crossings,
            along,
            first_faces,
            mismatches,
            report,
        })
    }
}

#[cfg(test)]
mod tests {
    //! The check must catch wrong data. A triangulation built through the public interface is
    //! never wrong, so these tests damage one from inside the crate.

    use ddg_core::generators::grid_disk;
    use ddg_core::{Geometry, TriMesh};

    use super::*;
    use crate::DelaunayOptions;

    /// A flat grid sheared by 1.5 cells per row, flipped to Delaunay, with its positions.
    fn sheared() -> (IntrinsicTriangulation, TriMesh) {
        let mut shape = grid_disk(8, 8).unwrap();
        for p in &mut shape.positions {
            p[0] += 1.5 * p[1];
        }
        let mesh = shape.topology().unwrap();
        let geometry = Geometry::from_positions(&mesh, &shape.positions).unwrap();
        let mut tri = IntrinsicTriangulation::new(&mesh, &geometry).unwrap();
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        (tri, mesh)
    }

    /// A flipped edge with at least two crossings whose end points are interior vertices.
    fn winding_edge(tri: &IntrinsicTriangulation, traces: &EdgeTraces) -> usize {
        (0..tri.n_edges())
            .find(|&e| {
                let [a, b] = tri.edge_vertices(e);
                traces.crossings(e).len() >= 2
                    && !tri.is_boundary_vertex(a)
                    && !tri.is_boundary_vertex(b)
            })
            .expect("the sheared grid has such an edge")
    }

    fn only_failure(tri: &IntrinsicTriangulation, e: usize) -> EdgeTraceFailure {
        let options = EdgeTraceOptions::default();
        let single = match tri.trace_edge(e, &options) {
            Err(IntrinsicError::EdgeTrace(failure)) => failure,
            other => panic!("expected a failed trace, got {other:?}"),
        };
        match tri.trace_edges(&options) {
            Err(IntrinsicError::EdgeTraces { n_edges, failures }) => {
                assert_eq!(n_edges, tri.n_edges());
                assert_eq!(failures, vec![single.clone()]);
            }
            other => panic!("expected failed traces, got {other:?}"),
        }
        assert_eq!(single.edge as usize, e);
        assert_eq!(single.vertices.map(|v| v as usize), tri.edge_vertices(e));
        assert_eq!(single.length, tri.edge_length(e));
        assert_eq!(single.tolerance, 1e-6 * tri.edge_length(e));
        single
    }

    #[test]
    fn a_wrong_signpost_is_a_failed_trace() {
        let (mut tri, _) = sheared();
        let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
        assert!(traces.report.max_relative_mismatch < 1e-13);
        let e = winding_edge(&tri, &traces);
        let h = tri.complex.edge_halfedge[e] as usize;
        // An interior vertex of a flat shape: normalized angles are angles. Turning the
        // signpost by 1e-4 moves the end of the path by 1e-4 of its length.
        tri.signposts[h] += 1e-4;
        let failure = only_failure(&tri, e);
        let length = tri.edge_length(e);
        match failure.reason {
            EdgeTraceFailureReason::Missed {
                mismatch: Some(d),
                n_crossings,
                end: SurfacePoint::Face { .. },
            } => {
                assert!((d / length - 1e-4).abs() < 1e-8, "{d} {length}");
                assert!(n_crossings >= 2);
            }
            ref other => panic!("{other:?}"),
        }
        let message = IntrinsicError::EdgeTrace(failure).to_string();
        let [a, b] = tri.edge_vertices(e);
        assert!(
            message.starts_with(&format!(
                "the trace of intrinsic edge {e} = ({a}, {b}) of length {length:?} over the input \
                 mesh ended at the point ("
            )),
            "{message}"
        );
        assert!(
            message.contains(&format!(
                "without reaching vertex {b}: the closest it came is"
            )),
            "{message}"
        );
        // With a tolerance that allows it the same trace is accepted and says how far off it
        // is.
        let loose = EdgeTraceOptions {
            tolerance: 1e-3,
            max_crossings: None,
        };
        let accepted = tri.trace_edge(e, &loose).unwrap();
        assert!((accepted.mismatch / length - 1e-4).abs() < 1e-8);
        assert_eq!(accepted.crossings.len(), traces.crossings(e).len());
    }

    #[test]
    fn a_wrong_length_is_a_failed_trace() {
        let (mut tri, _) = sheared();
        let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
        let e = winding_edge(&tri, &traces);
        let length = tri.edge_length(e);
        // Too short: the path ends before the vertex.
        tri.edge_lengths[e] = length * (1.0 - 1e-4);
        match only_failure(&tri, e).reason {
            EdgeTraceFailureReason::Missed {
                mismatch: Some(d), ..
            } => assert!((d / length - 1e-4).abs() < 1e-8, "{d}"),
            other => panic!("{other:?}"),
        }
        // Far too long: the path leaves the surface.
        tri.edge_lengths[e] = 100.0;
        match only_failure(&tri, e).reason {
            EdgeTraceFailureReason::HitBoundary {
                at: SurfacePoint::Edge { edge, .. },
                travelled,
            } => {
                assert!(tri.input().mesh().is_boundary_edge(edge as usize));
                assert!(travelled > length && travelled < 3.0);
            }
            other => panic!("{other:?}"),
        }
        let message = tri
            .trace_edge(e, &EdgeTraceOptions::default())
            .unwrap_err()
            .to_string();
        assert!(
            message.contains("left the surface at the point t = "),
            "{message}"
        );
    }

    #[test]
    fn a_path_into_another_vertex_is_a_failed_trace() {
        let (mut tri, mesh) = sheared();
        let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
        let e = winding_edge(&tri, &traces);
        let h = tri.complex.edge_halfedge[e] as usize;
        let [a, b] = tri.edge_vertices(e);
        // Point the signpost along an input edge at the tail that does not lead to the head,
        // and make the edge longer than that input edge.
        let along = mesh
            .outgoing_halfedges(a)
            .find(|&g| mesh.to_vertex(g) != b)
            .unwrap();
        let other = mesh.to_vertex(along);
        let input_length = tri.input().geometry().edge_lengths()[mesh.edge(along)];
        tri.signposts[h] = tri.input().signpost(along);
        tri.edge_lengths[e] = 1.5 * input_length;
        match only_failure(&tri, e).reason {
            EdgeTraceFailureReason::HitVertex { vertex, travelled } => {
                assert_eq!(vertex as usize, other);
                assert_eq!(travelled, input_length);
            }
            other => panic!("{other:?}"),
        }
        let message = tri
            .trace_edge(e, &EdgeTraceOptions::default())
            .unwrap_err()
            .to_string();
        assert!(
            message.ends_with(&format!(
                "ran into input vertex {other} after a length of {input_length:?}"
            )),
            "{message}"
        );
    }

    #[test]
    fn the_crossing_limit_is_a_failed_trace() {
        let (tri, _) = sheared();
        let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
        let e = winding_edge(&tri, &traces);
        let n = traces.crossings(e).len();
        let limited = EdgeTraceOptions {
            tolerance: 1e-6,
            max_crossings: Some(n - 1),
        };
        match tri.trace_edge(e, &limited) {
            Err(IntrinsicError::EdgeTrace(failure)) => match failure.reason {
                EdgeTraceFailureReason::CrossingLimit { limit, travelled } => {
                    assert_eq!(limit, n - 1);
                    assert_eq!(travelled, traces.crossings(e)[n - 1].distance);
                }
                other => panic!("{other:?}"),
            },
            other => panic!("{other:?}"),
        }
        let enough = EdgeTraceOptions {
            tolerance: 1e-6,
            max_crossings: Some(n),
        };
        assert_eq!(
            tri.trace_edge(e, &enough).unwrap().crossings,
            traces.crossings(e)
        );
    }
    /// SPEC N8: the path is followed for the tolerance beyond the edge's length. Here the far
    /// vertex sits just behind an input edge (a triangle with a height of 1e-4), and the edge
    /// is made a little too short, so that its path ends before that input edge.
    #[test]
    fn a_vertex_just_behind_an_input_edge_is_reached_within_the_tolerance() {
        // x below, u and w left and right, y barely above the segment u w.
        let positions = [
            [0.0, -1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1e-4, 0.0],
        ];
        let (x, u, w, y) = (0, 1, 2, 3);
        let mesh = TriMesh::new(4, vec![[0, 2, 1], [1, 2, 3]]).unwrap();
        let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
        let mut tri = IntrinsicTriangulation::new(&mesh, &geometry).unwrap();
        let e = mesh.find_edge(u, w).unwrap();
        tri.flip_edge(e).unwrap();
        if tri.edge_vertices(e) != [x, y] {
            // The other direction: two more flips turn the edge around.
            tri.flip_edge(e).unwrap();
            tri.flip_edge(e).unwrap();
        }
        assert_eq!(tri.edge_vertices(e), [x, y]);
        let length = tri.edge_length(e);
        // The triangle u w y has a height of 1e-4: its shape is known to 1e-12 only.
        assert!((length - 1.0001).abs() < 1e-11, "{length}");
        let options = EdgeTraceOptions {
            tolerance: 1e-3,
            max_crossings: None,
        };
        let exact = tri.trace_edge(e, &options).unwrap();
        assert_eq!(exact.crossings.len(), 1);
        assert!(exact.mismatch < 1e-11, "{}", exact.mismatch);
        assert!((exact.crossings[0].t - 0.5).abs() < 1e-11);
        assert!((exact.crossings[0].distance - 1.0).abs() < 1e-11);
        // Too short by 2e-4: the path ends 1e-4 before the input edge u w.
        tri.edge_lengths[e] = length - 2e-4;
        let short = tri.trace_edge(e, &options).unwrap();
        assert_eq!(short.crossings, exact.crossings);
        assert!((short.mismatch - 2e-4).abs() < 1e-11, "{}", short.mismatch);
        // With a tolerance below the error it is a failure: the path ends in the face below.
        let strict = EdgeTraceOptions {
            tolerance: 1e-5,
            max_crossings: None,
        };
        match tri.trace_edge(e, &strict) {
            Err(IntrinsicError::EdgeTrace(failure)) => match failure.reason {
                EdgeTraceFailureReason::Missed {
                    end: SurfacePoint::Face { face: 0, .. },
                    mismatch: None,
                    n_crossings: 0,
                } => {}
                other => panic!("{other:?}"),
            },
            other => panic!("{other:?}"),
        }
    }
}
