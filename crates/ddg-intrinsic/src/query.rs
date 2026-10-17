//! Point queries (paper Sec. 3.4.1, Alg. 11) and straight segments of the triangulation as
//! paths over the input mesh.
//!
//! A point of the triangulation is found on the input mesh by a tracing query over the input
//! mesh from a vertex of its face; a point of the input mesh is found on the triangulation by
//! a tracing query over the triangulation from an input vertex, which is a vertex of both with
//! the same frame.
//!
//! A straight segment inside a face of the triangulation is, on the input mesh, a path of
//! straight pieces that each lie in one input face. It is found by placing its start on the
//! input mesh and tracing from there; the trace is checked like the trace of an edge (SPEC
//! N8): it must arrive at the place of the segment's end, found independently, within a
//! tolerance. This is what carries a cut computed on the triangulation back to the input mesh
//! (SPEC §6).

use std::f64::consts::{PI, TAU};

use crate::edge_trace::EdgeTraceOptions;
use crate::error::IntrinsicError;
use crate::formulas::{signpost_update, triangle_layout};
use crate::insert::{EdgeSnap, LOCATION_SNAP};
use crate::point::{SurfacePoint, point_defect};
use crate::surface::Surface;
use crate::trace::{Target, TraceOptions, TraceStop, run};
use crate::triangulation::IntrinsicTriangulation;

/// A straight piece inside one face of the triangulation, between two points of a list.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Segment {
    /// The face of the triangulation the piece lies in.
    pub face: u32,
    /// Its two ends, as indices into the list of points; both are points of the closed face:
    /// a corner, a point of one of its edges, or a point inside it.
    pub nodes: [u32; 2],
}

/// A path over the input mesh: straight pieces that each lie in one input face.
#[derive(Debug, Clone, PartialEq)]
pub struct PathOnInput {
    /// The points of the path on the input mesh: its start, the points where it crosses input
    /// edges (each a [`SurfacePoint::Edge`]), its end.
    pub points: Vec<SurfacePoint>,
    /// For every piece (from `points[i]` to `points[i + 1]`) an input face that contains it.
    /// A piece that runs along an input edge lies in two faces; one of them is given.
    pub faces: Vec<u32>,
    /// The distance between the end of the traced path and the place of the segment's end
    /// (an absolute length).
    pub mismatch: f64,
}

/// The result of [`IntrinsicTriangulation::segments_on_input`].
#[derive(Debug, Clone, PartialEq)]
pub struct PathsOnInput {
    /// The place of every point of the list on the input mesh.
    pub nodes: Vec<SurfacePoint>,
    /// The path of every segment, from its first node to its second.
    pub paths: Vec<PathOnInput>,
    /// The largest mismatch of a path divided by the length of its segment.
    pub max_relative_mismatch: f64,
}

/// A point of the triangulation placed on the input mesh, with the rotation between the two
/// frames there.
#[derive(Debug, Clone, Copy)]
struct Located {
    /// The place on the input mesh.
    point: SurfacePoint,
    /// The direction, in the frame of `point` on the input mesh, of direction zero of the
    /// point's frame on the triangulation (along its edge, or along the first halfedge of its
    /// face). Not used for vertices.
    zero: f64,
}

fn wrap(angle: f64) -> f64 {
    let a = angle.rem_euclid(TAU);
    if a >= TAU { 0.0 } else { a }
}

impl IntrinsicTriangulation {
    fn check_point(&self, point: SurfacePoint) -> Result<(), IntrinsicError> {
        match point_defect(&point, self.n_vertices(), self.n_edges(), self.n_faces()) {
            Some(reason) => Err(IntrinsicError::InvalidPoint { point, reason }),
            None => Ok(()),
        }
    }

    fn degenerate(&self, f: usize) -> IntrinsicError {
        IntrinsicError::DegenerateFace {
            surface: Self::NAME,
            face: f as u32,
            lengths: self.face_lengths(f),
        }
    }

    /// The signpost direction at corner `k` of face `f` towards the point `p` of the face's
    /// layout `v`.
    fn direction_from_corner(&self, f: usize, k: usize, v: &[[f64; 2]; 3], p: [f64; 2]) -> f64 {
        let u = [v[(k + 1) % 3][0] - v[k][0], v[(k + 1) % 3][1] - v[k][1]];
        let w = [p[0] - v[k][0], p[1] - v[k][1]];
        let inside = (u[0] * w[1] - u[1] * w[0]).atan2(u[0] * w[0] + u[1] * w[1]);
        let vertex = self.faces[f][k] as usize;
        signpost_update(
            self.signposts[3 * f + k],
            inside.max(0.0),
            self.angle_sums[vertex],
            self.is_boundary_vertex(vertex),
        )
    }

    /// Follows a path over the input mesh from vertex `from` and returns its end with the
    /// direction back, or an error that says where it stopped.
    ///
    /// An end inside an input face within [`LOCATION_SNAP`] of the longest edge of that face,
    /// and within `snap_distance`, from an edge of the face is returned as the point of that
    /// input edge (SPEC N28): the end of a path to a point on an input edge comes out a
    /// rounding error inside one of the two faces of the edge, chosen by rounding, and a
    /// point given as a point of one face is met only by a path in that face. A point of the
    /// edge belongs to both.
    fn trace_from_vertex(
        &self,
        from: usize,
        phi: f64,
        length: f64,
        what: SurfacePoint,
        snap_distance: f64,
    ) -> Result<(SurfacePoint, f64), IntrinsicError> {
        let trace = run(
            &*self.input,
            self.locations[from],
            self.direction_on_input(from, phi),
            length,
            &TraceOptions::default(),
            None,
        )?
        .trace;
        // A point of a boundary edge: the path may reach the boundary a rounding error before
        // its full length. It has then arrived.
        let arrived = trace.stop == TraceStop::Completed
            || (trace.stop == TraceStop::Boundary && length - trace.length <= 1e-9 * length);
        if !arrived {
            return Err(IntrinsicError::InvalidPoint {
                point: what,
                reason: "the path to the point over the input mesh left the surface or ran into \
                         a vertex",
            });
        }
        Ok(
            match self.snap_onto_input_edge(
                trace.end,
                trace.back_direction,
                LOCATION_SNAP,
                snap_distance,
            ) {
                EdgeSnap::Onto { location, back, .. } => (location, back),
                // Exactly at a corner (two coordinates exactly zero) the point stays a point
                // of the face; a degenerate face cannot have been traced through.
                EdgeSnap::Far | EdgeSnap::Corner | EdgeSnap::Degenerate(_) => {
                    (trace.end, trace.back_direction)
                }
            },
        )
    }

    /// Places a point of an edge on the input mesh. Direction zero is along the edge.
    fn locate_on_edge(&self, e: usize, t: f64) -> Result<Located, IntrinsicError> {
        let h = self.edge_halfedge[e] as usize;
        let (a, b) = (self.from_vertex(h), self.to_vertex(h));
        if let Some(input_edge) = self.edge_on_input(e) {
            let mesh = self.input.mesh();
            let on = |v: usize| match self.locations[v] {
                SurfacePoint::Vertex(u) => {
                    let [lo, hi] = mesh.edge_vertices(input_edge);
                    (u as usize == lo)
                        .then_some(0.0)
                        .or((u as usize == hi).then_some(1.0))
                }
                SurfacePoint::Edge { edge, t } => (edge as usize == input_edge).then_some(t),
                SurfacePoint::Face { .. } => None,
            };
            if let (Some(ta), Some(tb)) = (on(a), on(b)) {
                return Ok(Located {
                    point: SurfacePoint::Edge {
                        edge: input_edge as u32,
                        t: (1.0 - t) * ta + t * tb,
                    },
                    zero: if tb > ta { 0.0 } else { PI },
                });
            }
        }
        let what = SurfacePoint::Edge { edge: e as u32, t };
        let (point, back) = self.trace_from_vertex(
            a,
            self.signposts[h],
            t * self.edge_lengths[e],
            what,
            LOCATION_SNAP * self.edge_lengths[e],
        )?;
        Ok(Located {
            point,
            zero: wrap(back + PI),
        })
    }

    /// Places a point inside a face on the input mesh. Direction zero is along the first
    /// halfedge of the face.
    fn locate_in_face(&self, f: usize, bary: [f64; 3]) -> Result<Located, IntrinsicError> {
        let v = triangle_layout(self.face_lengths(f)).ok_or_else(|| self.degenerate(f))?;
        let p = [
            bary[1] * v[1][0] + bary[2] * v[2][0],
            bary[1] * v[1][1] + bary[2] * v[2][1],
        ];
        // From the corner with the smallest weight: the one the point is farthest from in
        // relative terms, so that the direction is well defined.
        let k = (0..3)
            .min_by(|&a, &b| bary[a].total_cmp(&bary[b]))
            .unwrap_or(0);
        let to_corner = [v[k][0] - p[0], v[k][1] - p[1]];
        let length = to_corner[0].hypot(to_corner[1]);
        let phi = self.direction_from_corner(f, k, &v, p);
        let what = SurfacePoint::Face {
            face: f as u32,
            bary,
        };
        let l = self.face_lengths(f);
        let (point, back) = self.trace_from_vertex(
            self.faces[f][k] as usize,
            phi,
            length,
            what,
            LOCATION_SNAP * l[0].min(l[1]).min(l[2]),
        )?;
        // `back` points from the point to the corner; in the frame of the face that
        // direction has the angle of `to_corner`.
        Ok(Located {
            point,
            zero: wrap(back - to_corner[1].atan2(to_corner[0])),
        })
    }

    /// The place on the input mesh of a point of the triangulation (paper Alg. 11): the
    /// stored place for a vertex; for a point of an edge or a face, the end of a tracing
    /// query over the input mesh, or an interpolation along the input edge for an edge that
    /// lies along one. An end of a tracing query within [`LOCATION_SNAP`] of the longest edge
    /// of its input face, and within `LOCATION_SNAP` of the edge the point lies on (of the
    /// shortest side of the face the point lies in), from an edge of that input face is
    /// returned as the point of that input edge, which belongs to both faces of the edge
    /// (SPEC N29).
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`], and the errors of a tracing query.
    pub fn point_on_input(&self, point: SurfacePoint) -> Result<SurfacePoint, IntrinsicError> {
        self.check_point(point)?;
        Ok(match point {
            SurfacePoint::Vertex(v) => self.locations[v as usize],
            SurfacePoint::Edge { edge, t } if t == 0.0 || t == 1.0 => {
                let ends = self.edge_vertices(edge as usize);
                self.locations[ends[usize::from(t == 1.0)]]
            }
            SurfacePoint::Edge { edge, t } => self.locate_on_edge(edge as usize, t)?.point,
            SurfacePoint::Face { face, bary } => match bary.iter().position(|&b| b == 1.0) {
                Some(k) => self.locations[self.faces[face as usize][k] as usize],
                None => self.locate_in_face(face as usize, bary)?.point,
            },
        })
    }

    /// The place on the triangulation of a point of the input mesh (paper Alg. 11): the
    /// vertex itself for an input vertex; otherwise the end of a tracing query over the
    /// triangulation from an input vertex of the point's edge or face.
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`], and the errors of a tracing query.
    pub fn point_on_triangulation(
        &self,
        point: SurfacePoint,
    ) -> Result<SurfacePoint, IntrinsicError> {
        let mesh = self.input.mesh();
        if let Some(reason) =
            point_defect(&point, mesh.n_vertices(), mesh.n_edges(), mesh.n_faces())
        {
            return Err(IntrinsicError::InvalidPoint { point, reason });
        }
        let lengths = self.input.geometry().edge_lengths();
        let (vertex, phi, length) = match point {
            SurfacePoint::Vertex(v) => return Ok(SurfacePoint::Vertex(v)),
            SurfacePoint::Edge { edge, t } => {
                let [lo, hi] = mesh.edge_vertices(edge as usize);
                if t == 0.0 || t == 1.0 {
                    return Ok(SurfacePoint::Vertex(if t == 0.0 { lo } else { hi } as u32));
                }
                let h = mesh.edge_halfedge(edge as usize);
                let tail = mesh.from_vertex(h);
                let from_tail = if tail == lo { t } else { 1.0 - t };
                (
                    tail,
                    self.input.signpost(h),
                    from_tail * lengths[edge as usize],
                )
            }
            SurfacePoint::Face { face, bary } => {
                let f = face as usize;
                if let Some(k) = bary.iter().position(|&b| b == 1.0) {
                    return Ok(SurfacePoint::Vertex(mesh.faces()[f][k]));
                }
                let l = [0, 1, 2].map(|s| lengths[mesh.edge(3 * f + s)]);
                let v = triangle_layout(l).ok_or(IntrinsicError::DegenerateFace {
                    surface: "input mesh",
                    face,
                    lengths: l,
                })?;
                let p = [
                    bary[1] * v[1][0] + bary[2] * v[2][0],
                    bary[1] * v[1][1] + bary[2] * v[2][1],
                ];
                let k = (0..3)
                    .min_by(|&a, &b| bary[a].total_cmp(&bary[b]))
                    .unwrap_or(0);
                let u = [v[(k + 1) % 3][0] - v[k][0], v[(k + 1) % 3][1] - v[k][1]];
                let w = [p[0] - v[k][0], p[1] - v[k][1]];
                let inside = (u[0] * w[1] - u[1] * w[0]).atan2(u[0] * w[0] + u[1] * w[1]);
                let vertex = mesh.faces()[f][k] as usize;
                let phi = signpost_update(
                    self.input.signpost(3 * f + k),
                    inside.max(0.0),
                    self.input.vertex_angle_sum(vertex),
                    mesh.is_boundary_vertex(vertex),
                );
                (vertex, phi, w[0].hypot(w[1]))
            }
        };
        // A path along edges of the triangulation stops at every vertex on its way (an
        // inserted vertex on the input edge); it goes on in the same direction.
        let mut at = vertex;
        let mut direction = phi;
        let mut remaining = length;
        for _ in 0..=self.n_vertices() {
            let trace = run(
                self,
                SurfacePoint::Vertex(at as u32),
                direction,
                remaining,
                &TraceOptions::default(),
                None,
            )?
            .trace;
            match (trace.stop, trace.end) {
                (TraceStop::Completed, end) => return Ok(end),
                (TraceStop::Vertex, SurfacePoint::Vertex(v))
                    if self.angle_sums[v as usize]
                        == if self.is_boundary_vertex(v as usize) {
                            PI
                        } else {
                            TAU
                        } =>
                {
                    // Straight on through a flat vertex: the opposite of the way back.
                    remaining -= trace.length;
                    at = v as usize;
                    direction = if self.is_boundary_vertex(at) {
                        // Along the boundary: the way back is 0 or π, the way on the other.
                        PI - trace.back_direction
                    } else {
                        wrap(trace.back_direction + PI)
                    };
                    if remaining <= 0.0 || remaining.is_nan() {
                        return Ok(SurfacePoint::Vertex(v));
                    }
                }
                _ => break,
            }
        }
        Err(IntrinsicError::InvalidPoint {
            point,
            reason: "the path to the point over the triangulation left the surface or ran into \
                     a vertex",
        })
    }

    /// The position of a node of a segment in the layout `v` of face `f`, the slot of the
    /// side or corner it lies on, and its place on the input mesh.
    fn node_in_face(
        &self,
        f: usize,
        v: &[[f64; 2]; 3],
        node: SurfacePoint,
    ) -> Result<([f64; 2], NodeKind), IntrinsicError> {
        let outside = IntrinsicError::InvalidPoint {
            point: node,
            reason: "the point is not a point of the face of its segment",
        };
        match node {
            SurfacePoint::Vertex(u) => {
                let k = (0..3).find(|&k| self.faces[f][k] == u).ok_or(outside)?;
                Ok((v[k], NodeKind::Corner(k)))
            }
            SurfacePoint::Edge { edge, t } => {
                let k = (0..3)
                    .find(|&k| self.halfedge_edge[3 * f + k] == edge)
                    .ok_or(outside)?;
                let forward = self.edge_halfedge[edge as usize] as usize == 3 * f + k;
                let th = if forward { t } else { 1.0 - t };
                let (a, b) = (v[k], v[(k + 1) % 3]);
                Ok((
                    [a[0] + th * (b[0] - a[0]), a[1] + th * (b[1] - a[1])],
                    NodeKind::Side { slot: k, forward },
                ))
            }
            SurfacePoint::Face { face, bary } => {
                if face as usize != f {
                    return Err(outside);
                }
                Ok((
                    [
                        bary[1] * v[1][0] + bary[2] * v[2][0],
                        bary[1] * v[1][1] + bary[2] * v[2][1],
                    ],
                    NodeKind::Inside,
                ))
            }
        }
    }

    /// Carries straight segments of the triangulation to the input mesh (SPEC §4, §6).
    ///
    /// `nodes` are points of the triangulation; every segment joins two of them inside one
    /// face. Each node is placed on the input mesh once ([`point_on_input`](Self::point_on_input)),
    /// so that segments that share a node meet there exactly. Each segment is then traced over
    /// the input mesh from the place of its first node in the direction of its second, and
    /// must arrive at the place of the second within `options.tolerance` times its length.
    /// A segment along a side of its face that lies along an input edge is not traced: its
    /// path is the piece of that input edge between the places of its ends, whose length must
    /// be its length within the same tolerance (SPEC N29).
    ///
    /// A node must not be an end of an edge or a corner in disguise (`t` of 0 or 1, a
    /// barycentric coordinate of 1): use the vertex. In a face that uses a vertex or an edge
    /// twice, a node given as that vertex or as a point of that edge is taken at the first
    /// corner or side of the face that fits; a simplicial triangulation
    /// ([`export_simplicial`](Self::export_simplicial)) has no such faces, and
    /// [`segment_on_input`](Self::segment_on_input) takes barycentric coordinates.
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`] for a node that does not exist or does not lie in the
    /// face of a segment that uses it, [`IntrinsicError::IndexOutOfRange`],
    /// [`IntrinsicError::InvalidOptions`], and [`IntrinsicError::SegmentTrace`] for a segment
    /// whose trace does not arrive.
    pub fn segments_on_input(
        &self,
        nodes: &[SurfacePoint],
        segments: &[Segment],
        options: &EdgeTraceOptions,
    ) -> Result<PathsOnInput, IntrinsicError> {
        if !(options.tolerance.is_finite() && options.tolerance >= 0.0) {
            return Err(IntrinsicError::InvalidOptions {
                reason: format!(
                    "the tolerance of a segment trace must be finite and not negative, got {}",
                    options.tolerance
                ),
            });
        }
        let mut located = Vec::with_capacity(nodes.len());
        for &node in nodes {
            self.check_point(node)?;
            located.push(match node {
                SurfacePoint::Vertex(v) => Located {
                    point: self.locations[v as usize],
                    zero: 0.0,
                },
                SurfacePoint::Edge { edge, t } => {
                    if t == 0.0 || t == 1.0 {
                        return Err(IntrinsicError::InvalidPoint {
                            point: node,
                            reason: "the end of an edge must be given as the vertex",
                        });
                    }
                    self.locate_on_edge(edge as usize, t)?
                }
                SurfacePoint::Face { face, bary } => {
                    if bary.contains(&1.0) {
                        return Err(IntrinsicError::InvalidPoint {
                            point: node,
                            reason: "the corner of a face must be given as the vertex",
                        });
                    }
                    self.locate_in_face(face as usize, bary)?
                }
            });
        }
        let mut paths = Vec::with_capacity(segments.len());
        let mut max_relative_mismatch = 0.0f64;
        for segment in segments {
            let f = segment.face as usize;
            if f >= self.n_faces() {
                return Err(IntrinsicError::IndexOutOfRange {
                    what: "face",
                    index: f,
                    bound: self.n_faces(),
                });
            }
            for &n in &segment.nodes {
                if n as usize >= nodes.len() {
                    return Err(IntrinsicError::IndexOutOfRange {
                        what: "node",
                        index: n as usize,
                        bound: nodes.len(),
                    });
                }
            }
            let [from, to] = segment.nodes.map(|n| n as usize);
            let v = triangle_layout(self.face_lengths(f)).ok_or_else(|| self.degenerate(f))?;
            let (p, kind) = self.node_in_face(f, &v, nodes[from])?;
            let (q, kind_to) = self.node_in_face(f, &v, nodes[to])?;
            let path = self.segment_path(
                f,
                &v,
                (kind, p, located[from]),
                (kind_to, q, located[to]),
                options,
            )?;
            let d = [q[0] - p[0], q[1] - p[1]];
            max_relative_mismatch = max_relative_mismatch.max(path.mismatch / d[0].hypot(d[1]));
            paths.push(path);
        }
        Ok(PathsOnInput {
            nodes: located.iter().map(|l| l.point).collect(),
            paths,
            max_relative_mismatch,
        })
    }

    /// Traces the segment of face `f` (layout `v`) from the point `p` to the point `q` over
    /// the input mesh, from the place of `p` to the place of `q`.
    ///
    /// A segment along a side of the face that lies along an input edge
    /// ([`edge_on_input`](Self::edge_on_input)) is not traced: it is the piece of that input
    /// edge between the places of its ends (SPEC N28), as an edge along an input edge is not
    /// traced (SPEC N8). The check is the same: the length of that piece must be the length
    /// of the segment within the tolerance.
    fn segment_path(
        &self,
        f: usize,
        v: &[[f64; 2]; 3],
        (kind, p, from): (NodeKind, [f64; 2], Located),
        (kind_to, q, to): (NodeKind, [f64; 2], Located),
        options: &EdgeTraceOptions,
    ) -> Result<PathOnInput, IntrinsicError> {
        let mesh = self.input.mesh();
        let d = [q[0] - p[0], q[1] - p[1]];
        let length = d[0].hypot(d[1]);
        if length <= 0.0 || length.is_nan() {
            return Err(IntrinsicError::InvalidLength { length });
        }
        let tolerance = options.tolerance * length;
        if let Some(path) = self.path_along_input_edge(f, (kind, from), (kind_to, to), length) {
            if path.mismatch > tolerance {
                return Err(IntrinsicError::SegmentTrace {
                    face: f as u32,
                    length,
                    tolerance,
                    mismatch: Some(path.mismatch),
                });
            }
            return Ok(path);
        }
        // The direction of the segment in the frame of the face, then in the frame of the
        // place of its start on the input mesh.
        let in_face = d[1].atan2(d[0]);
        let direction = match kind {
            NodeKind::Corner(k) => {
                let vertex = self.faces[f][k] as usize;
                self.direction_on_input(vertex, self.direction_from_corner(f, k, v, q))
            }
            NodeKind::Side { slot, forward } => {
                let s = [
                    v[(slot + 1) % 3][0] - v[slot][0],
                    v[(slot + 1) % 3][1] - v[slot][1],
                ];
                let along = s[1].atan2(s[0]) + if forward { 0.0 } else { PI };
                wrap(from.zero + in_face - along)
            }
            NodeKind::Inside => wrap(from.zero + in_face),
        };
        let target = Target {
            point: to.point,
            tolerance,
        };
        let trace_options = TraceOptions {
            record_crossings: true,
            max_crossings: options.max_crossings,
        };
        let outcome = run(
            &*self.input,
            from.point,
            direction,
            length,
            &trace_options,
            Some(target),
        )?;
        let Some(mismatch) = outcome.mismatch else {
            return Err(IntrinsicError::SegmentTrace {
                face: f as u32,
                length,
                tolerance,
                mismatch: outcome.closest,
            });
        };
        let trace = outcome.trace;
        let mut points = Vec::with_capacity(trace.crossings.len() + 2);
        let mut faces = Vec::with_capacity(trace.crossings.len() + 1);
        points.push(from.point);
        let mut face = match (trace.first_face, trace.along_edge) {
            (Some(face), _) => face,
            (None, Some(edge)) => (mesh.edge_halfedge(edge as usize) / 3) as u32,
            (None, None) => 0,
        };
        for crossing in &trace.crossings {
            faces.push(face);
            points.push(SurfacePoint::Edge {
                edge: crossing.edge,
                t: crossing.t,
            });
            // The path goes on in the face on the other side of the crossed halfedge.
            face = mesh
                .twin(crossing.halfedge as usize)
                .map_or(face, |t| (t / 3) as u32);
        }
        faces.push(face);
        points.push(to.point);
        Ok(PathOnInput {
            points,
            faces,
            mismatch,
        })
    }

    /// The straight segment between two points of face `f` of the triangulation, given in
    /// barycentric coordinates, as a path over the input mesh: the ordered points where it
    /// crosses input edges, between the places of its two ends, and the input face of every
    /// piece. A point with a barycentric coordinate of exactly zero is a point of that side
    /// of the face, one with two zeros is that corner (in a face that uses a vertex twice the
    /// coordinates say which corner is meant, which a vertex id could not).
    ///
    /// The path is traced over the input mesh from the place of the first point and must
    /// arrive at the place of the second, found independently, within `options.tolerance`
    /// times the length of the segment.
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`], [`IntrinsicError::InvalidOptions`],
    /// [`IntrinsicError::InvalidLength`] if the two points coincide, and
    /// [`IntrinsicError::SegmentTrace`] if the trace does not arrive.
    pub fn segment_on_input(
        &self,
        f: usize,
        from: [f64; 3],
        to: [f64; 3],
        options: &EdgeTraceOptions,
    ) -> Result<PathOnInput, IntrinsicError> {
        if !(options.tolerance.is_finite() && options.tolerance >= 0.0) {
            return Err(IntrinsicError::InvalidOptions {
                reason: format!(
                    "the tolerance of a segment trace must be finite and not negative, got {}",
                    options.tolerance
                ),
            });
        }
        for bary in [from, to] {
            self.check_point(SurfacePoint::Face {
                face: f as u32,
                bary,
            })?;
        }
        let v = triangle_layout(self.face_lengths(f)).ok_or_else(|| self.degenerate(f))?;
        let end = |bary: [f64; 3]| -> Result<(NodeKind, [f64; 2], Located), IntrinsicError> {
            let p = [
                bary[1] * v[1][0] + bary[2] * v[2][0],
                bary[1] * v[1][1] + bary[2] * v[2][1],
            ];
            let zeros: Vec<usize> = (0..3).filter(|&k| bary[k] == 0.0).collect();
            Ok(match zeros[..] {
                [i, j] => {
                    let k = 3 - i - j;
                    let located = Located {
                        point: self.locations[self.faces[f][k] as usize],
                        zero: 0.0,
                    };
                    (NodeKind::Corner(k), v[k], located)
                }
                [k] => {
                    // On the side opposite corner k: halfedge k + 1.
                    let slot = (k + 1) % 3;
                    let h = 3 * f + slot;
                    let edge = self.halfedge_edge[h] as usize;
                    let th = bary[(k + 2) % 3] / (bary[(k + 1) % 3] + bary[(k + 2) % 3]);
                    let forward = self.edge_halfedge[edge] as usize == h;
                    let located = self.locate_on_edge(edge, if forward { th } else { 1.0 - th })?;
                    (NodeKind::Side { slot, forward }, p, located)
                }
                _ => (NodeKind::Inside, p, self.locate_in_face(f, bary)?),
            })
        };
        let (kind, p, from) = end(from)?;
        let (kind_to, q, to) = end(to)?;
        self.segment_path(f, &v, (kind, p, from), (kind_to, q, to), options)
    }
}

impl IntrinsicTriangulation {
    /// The path over the input mesh of a segment of face `f` of length `length` whose ends
    /// lie on the same side of the face, if that side lies along an input edge and the
    /// places of both ends are on that input edge (they are: a vertex of the side is placed
    /// on it, a point of the side is placed by interpolation along it); `None` otherwise.
    /// Its mismatch is the difference between the length of the piece of the input edge and
    /// `length`.
    fn path_along_input_edge(
        &self,
        f: usize,
        (kind, from): (NodeKind, Located),
        (kind_to, to): (NodeKind, Located),
        length: f64,
    ) -> Option<PathOnInput> {
        let side = kind
            .sides()
            .into_iter()
            .flatten()
            .find(|s| kind_to.sides().contains(&Some(*s)))?;
        let input_edge = self.edge_on_input(self.halfedge_edge[3 * f + side] as usize)?;
        let mesh = self.input.mesh();
        let parameter = |point: SurfacePoint| match point {
            SurfacePoint::Vertex(u) => {
                let [lo, hi] = mesh.edge_vertices(input_edge);
                if u as usize == lo {
                    Some(0.0)
                } else if u as usize == hi {
                    Some(1.0)
                } else {
                    None
                }
            }
            SurfacePoint::Edge { edge, t } => (edge as usize == input_edge).then_some(t),
            SurfacePoint::Face { .. } => None,
        };
        let (ta, tb) = (parameter(from.point)?, parameter(to.point)?);
        let piece = (tb - ta).abs() * self.input.geometry().edge_lengths()[input_edge];
        Some(PathOnInput {
            points: vec![from.point, to.point],
            faces: vec![(mesh.edge_halfedge(input_edge) / 3) as u32],
            mismatch: (piece - length).abs(),
        })
    }
}

/// Where a node of a segment lies in the face of the segment.
#[derive(Debug, Clone, Copy)]
enum NodeKind {
    /// At the corner in this slot.
    Corner(usize),
    /// On the side in this slot; `forward` if that halfedge gives the edge its direction.
    Side { slot: usize, forward: bool },
    /// Inside the face.
    Inside,
}

impl NodeKind {
    /// The sides of the face the node lies on: a corner on the two at it (the side leaving it
    /// and the side arriving at it), a point of a side on that one.
    fn sides(self) -> [Option<usize>; 2] {
        match self {
            NodeKind::Corner(k) => [Some(k), Some((k + 2) % 3)],
            NodeKind::Side { slot, .. } => [Some(slot), None],
            NodeKind::Inside => [None, None],
        }
    }
}
