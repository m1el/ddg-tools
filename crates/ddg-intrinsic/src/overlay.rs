//! The overlay: the common subdivision of the input mesh and the intrinsic triangulation
//! (paper Sec. 3.4.2).
//!
//! Its vertices are the vertices of the triangulation (which include all input vertices) and
//! the points where an edge of the triangulation crosses an input edge. Its faces are the
//! pieces into which the input edges cut the faces of the triangulation; each lies in one face
//! of the triangulation and in one input face, and is a convex polygon with three to six
//! corners there.
//!
//! # Construction
//!
//! The crossings come from the traces of the edges ([`EdgeTraces`]), which give, for every
//! edge, its crossings with input edges in order. Inside a face of the triangulation the
//! input edges are straight segments that do not cross each other, each from a side to a side
//! or from a corner to the opposite side (a corner is an input vertex, or a vertex inserted on
//! an input edge). Such a set of segments is determined by the numbers of crossings on the
//! three sides: with `n_k` crossings on side `k`, either `n_{k−1} + n_k − n_{k+1}` is an even
//! number `2 x_k ≥ 0` for every corner `k`, and `x_k` segments cut off corner `k`; or one
//! side has more crossings than the other two together, by `c`, and `c` segments run from the
//! opposite corner to it (the paper's cases 2 and 1). The segments are paired from the corners
//! inwards. The input edge of the two crossings a segment joins must be the same, and for a
//! segment from a corner that corner must lie on the input edge; if this fails, or the counts
//! fit neither case, the traces do not describe a subdivision and the overlay is refused
//! ([`IntrinsicError::InconsistentOverlay`]) instead of being built wrong.

use ddg_core::topology::NONE;

use crate::edge_trace::EdgeTraces;
use crate::error::IntrinsicError;
use crate::point::SurfacePoint;
use crate::triangulation::IntrinsicTriangulation;

/// A vertex of the overlay, as a point of both triangulations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OverlayVertex {
    /// The point on the input mesh: the place of a vertex of the triangulation, or a point of
    /// the input edge that is crossed.
    pub on_input: SurfacePoint,
    /// The point on the triangulation: a vertex, or a point of the edge that crosses.
    pub on_triangulation: SurfacePoint,
}

/// The common subdivision of the input mesh and the triangulation. Returned by
/// [`IntrinsicTriangulation::overlay`].
#[derive(Debug, Clone, PartialEq)]
pub struct Overlay {
    /// The vertices: first the vertices of the triangulation, with their ids, then one vertex
    /// per crossing of an edge of the triangulation with an input edge, in the order of the
    /// edges and along each edge from its tail.
    pub vertices: Vec<OverlayVertex>,
    /// For every face, the index of its first corner in [`corners`](Overlay::corners); one
    /// more entry at the end.
    pub face_starts: Vec<u32>,
    /// The corners of the faces, counter-clockwise, as vertex indices.
    pub corners: Vec<u32>,
    /// For every face, the input face it lies in.
    pub face_input: Vec<u32>,
    /// For every face, the face of the triangulation it lies in. Faces are in the order of
    /// the faces of the triangulation.
    pub face_triangulation: Vec<u32>,
}

impl Overlay {
    /// Number of faces.
    pub fn n_faces(&self) -> usize {
        self.face_input.len()
    }

    /// The corners of face `f`.
    pub fn face(&self, f: usize) -> &[u32] {
        &self.corners[self.face_starts[f] as usize..self.face_starts[f + 1] as usize]
    }

    /// The faces cut into triangles (every face is convex: a fan from its first corner), and
    /// for every triangle the face it is a part of.
    pub fn triangles(&self) -> (Vec<[u32; 3]>, Vec<u32>) {
        let mut triangles = Vec::with_capacity(self.corners.len());
        let mut parents = Vec::with_capacity(self.corners.len());
        for f in 0..self.n_faces() {
            let corners = self.face(f);
            for k in 1..corners.len() - 1 {
                triangles.push([corners[0], corners[k], corners[k + 1]]);
                parents.push(f as u32);
            }
        }
        (triangles, parents)
    }
}

/// The crossings on one side of a face, seen in the direction of that side.
struct Side<'a> {
    tri: &'a IntrinsicTriangulation,
    traces: &'a EdgeTraces,
    edge: usize,
    forward: bool,
    n: usize,
    /// Index of the edge's first crossing among all crossings.
    offset: usize,
}

impl Side<'_> {
    /// The index, in the direction of the edge, of the `i`-th crossing in the direction of
    /// the side.
    fn along_edge(&self, i: usize) -> usize {
        if self.forward { i } else { self.n - 1 - i }
    }

    /// The overlay vertex of the `i`-th crossing from the start of the side.
    fn vertex(&self, i: usize) -> u32 {
        (self.tri.n_vertices() + self.offset + self.along_edge(i)) as u32
    }

    /// The input edge of the `i`-th crossing from the start of the side.
    fn input_edge(&self, i: usize) -> u32 {
        self.traces.crossings(self.edge)[self.along_edge(i)].edge
    }

    /// The input face that contains the piece of the side between its crossings `i − 1` and
    /// `i` (from the start for `i = 0`, to the end for `i = n`), or, for a side that lies
    /// along an input edge, the input face on its left.
    fn input_face(&self, i: usize) -> Option<u32> {
        let mesh = self.tri.input().mesh();
        if let Some(input_edge) = self.traces.along_input_edge(self.edge) {
            // The halfedge of the input edge that runs the way the side does.
            let [a, b] = self.tri.edge_vertices(self.edge);
            let (tail, head) = if self.forward { (a, b) } else { (b, a) };
            let [lo, hi] = mesh.edge_vertices(input_edge as usize);
            let at = |v: usize| match self.tri.vertex_location(v) {
                SurfacePoint::Vertex(u) if u as usize == lo => Some(0.0),
                SurfacePoint::Vertex(u) if u as usize == hi => Some(1.0),
                SurfacePoint::Edge { edge, t } if edge == input_edge => Some(t),
                _ => None,
            };
            let ascending = at(tail)? < at(head)?;
            let (from, to) = if ascending { (lo, hi) } else { (hi, lo) };
            return mesh.find_halfedge(from, to).map(|h| (h / 3) as u32);
        }
        // Piece r of the edge, in its own direction, lies in the face entered by crossing
        // r − 1, or in the first face.
        let r = if self.forward { i } else { self.n - i };
        if r == 0 {
            return self.traces.first_face(self.edge);
        }
        let crossed = self.traces.crossings(self.edge)[r - 1].halfedge as usize;
        mesh.twin(crossed).map(|t| (t / 3) as u32)
    }
}

impl IntrinsicTriangulation {
    /// The overlay of the input mesh and the triangulation, from the traces of the
    /// triangulation's edges ([`trace_edges`](Self::trace_edges)), which must be those of
    /// this triangulation in its present state. See the [module documentation](self).
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidOptions`] if the traces are of another triangulation (a
    /// different number of edges), and [`IntrinsicError::InconsistentOverlay`] if the
    /// crossings around a face cannot be paired.
    pub fn overlay(&self, traces: &EdgeTraces) -> Result<Overlay, IntrinsicError> {
        if traces.n_edges() != self.n_edges() {
            return Err(IntrinsicError::InvalidOptions {
                reason: format!(
                    "the traces are of a triangulation with {} edges, this one has {}",
                    traces.n_edges(),
                    self.n_edges()
                ),
            });
        }
        let mesh = self.input.mesh();
        let nv = self.n_vertices();
        // Vertices.
        let mut offsets = Vec::with_capacity(self.n_edges());
        let mut vertices: Vec<OverlayVertex> = (0..nv)
            .map(|v| OverlayVertex {
                on_input: self.locations[v],
                on_triangulation: SurfacePoint::Vertex(v as u32),
            })
            .collect();
        for e in 0..self.n_edges() {
            offsets.push(vertices.len() - nv);
            let length = self.edge_lengths[e];
            for c in traces.crossings(e) {
                vertices.push(OverlayVertex {
                    on_input: SurfacePoint::Edge {
                        edge: c.edge,
                        t: c.t,
                    },
                    on_triangulation: SurfacePoint::Edge {
                        edge: e as u32,
                        t: (c.distance / length).clamp(0.0, 1.0),
                    },
                });
            }
        }
        // Whether a vertex of the triangulation lies on an input edge (as an end or inside).
        let on_input_edge = |v: usize, input_edge: u32| match self.locations[v] {
            SurfacePoint::Vertex(u) => mesh
                .edge_vertices(input_edge as usize)
                .contains(&(u as usize)),
            SurfacePoint::Edge { edge, .. } => edge == input_edge,
            SurfacePoint::Face { .. } => false,
        };

        let mut overlay = Overlay {
            vertices,
            face_starts: vec![0],
            corners: Vec::new(),
            face_input: Vec::new(),
            face_triangulation: Vec::new(),
        };
        for f in 0..self.n_faces() {
            let bad = |reason: String| IntrinsicError::InconsistentOverlay {
                face: f as u32,
                reason,
            };
            let sides: [Side<'_>; 3] = [0, 1, 2].map(|k| {
                let h = 3 * f + k;
                let edge = self.halfedge_edge[h] as usize;
                Side {
                    tri: self,
                    traces,
                    edge,
                    forward: self.edge_halfedge[edge] as usize == h,
                    n: traces.crossings(edge).len(),
                    offset: offsets[edge],
                }
            });
            let n = [sides[0].n, sides[1].n, sides[2].n];
            let corner = |k: usize| self.faces[f][k];
            // x[k]: segments that cut off corner k; fan: the corner that sends segments to
            // the opposite side, and how many.
            let mut x = [0usize; 3];
            let mut fan: Option<(usize, usize)> = None;
            for k in 0..3 {
                // The side opposite corner k is side k + 1.
                let (before, after, opposite) = (n[(k + 2) % 3], n[k], n[(k + 1) % 3]);
                if opposite > before + after {
                    fan = Some((k, opposite - before - after));
                }
            }
            match fan {
                Some((k, _)) => {
                    x[(k + 1) % 3] = n[k];
                    x[(k + 2) % 3] = n[(k + 2) % 3];
                }
                None => {
                    for k in 0..3 {
                        let twice = n[(k + 2) % 3] + n[k] - n[(k + 1) % 3];
                        if !twice.is_multiple_of(2) {
                            return Err(bad(format!(
                                "its sides have {}, {} and {} crossings, which cannot be joined \
                                 in pairs",
                                n[0], n[1], n[2]
                            )));
                        }
                        x[k] = twice / 2;
                    }
                }
            }
            // The two ends of the m-th segment at corner k, counted from the corner: on side
            // k from its start, and on side k − 1 from its end.
            let leaving = |k: usize, m: usize| (k, m);
            let arriving = |k: usize, m: usize| ((k + 2) % 3, n[(k + 2) % 3] - 1 - m);
            for (k, &segments) in x.iter().enumerate() {
                for m in 0..segments {
                    let (sa, ia) = leaving(k, m);
                    let (sb, ib) = arriving(k, m);
                    if sides[sa].input_edge(ia) != sides[sb].input_edge(ib) {
                        return Err(bad(format!(
                            "segment {m} at corner {k} would join crossings of the input edges \
                             {} and {}",
                            sides[sa].input_edge(ia),
                            sides[sb].input_edge(ib)
                        )));
                    }
                }
            }
            if let Some((k, c)) = fan {
                let side = &sides[(k + 1) % 3];
                for j in 0..c {
                    let input_edge = side.input_edge(x[(k + 1) % 3] + j);
                    if !on_input_edge(corner(k) as usize, input_edge) {
                        return Err(bad(format!(
                            "input edge {input_edge} crosses the side opposite vertex {} and \
                             neither of the other sides, but does not pass through that vertex",
                            corner(k)
                        )));
                    }
                }
            }

            let mut emit =
                |corners: &[u32], input_face: Option<u32>| -> Result<(), IntrinsicError> {
                    let Some(input_face) = input_face else {
                        return Err(bad("a piece has no input face".to_string()));
                    };
                    overlay.corners.extend_from_slice(corners);
                    overlay.face_starts.push(overlay.corners.len() as u32);
                    overlay.face_input.push(input_face);
                    overlay.face_triangulation.push(f as u32);
                    Ok(())
                };
            // The corner of the piece next to corner k that is left when its x[k] segments
            // are cut off: the ends of the last segment, or the corner itself.
            let cut_corner = |k: usize| -> Vec<u32> {
                if x[k] == 0 {
                    vec![corner(k)]
                } else {
                    let (sb, ib) = arriving(k, x[k] - 1);
                    let (sa, ia) = leaving(k, x[k] - 1);
                    vec![sides[sb].vertex(ib), sides[sa].vertex(ia)]
                }
            };
            // The pieces at the corners.
            for (k, &segments) in x.iter().enumerate() {
                for m in 0..segments {
                    let (sa, ia) = leaving(k, m);
                    let (sb, ib) = arriving(k, m);
                    let (a, b) = (sides[sa].vertex(ia), sides[sb].vertex(ib));
                    if m == 0 {
                        emit(&[corner(k), a, b], sides[k].input_face(0))?;
                    } else {
                        let (pa, pb) = (sides[sa].vertex(ia - 1), sides[sb].vertex(ib + 1));
                        emit(&[pa, a, b, pb], sides[k].input_face(m))?;
                    }
                }
            }
            match fan {
                None => {
                    // The piece in the middle.
                    let mut polygon = Vec::with_capacity(6);
                    for k in 0..3 {
                        polygon.extend(cut_corner(k));
                    }
                    // A piece of side k that bounds it: after the segments of corner k.
                    let face = (0..3).find_map(|k| sides[k].input_face(x[k]));
                    emit(&polygon, face)?;
                }
                Some((k, c)) => {
                    let (k1, k2) = ((k + 1) % 3, (k + 2) % 3);
                    let on_opposite = |j: usize| sides[k1].vertex(x[k1] + j);
                    // From corner k along side k to the first segment of the fan.
                    let mut first = vec![corner(k)];
                    first.extend(cut_corner(k1));
                    first.push(on_opposite(0));
                    emit(&first, sides[k].input_face(0))?;
                    for j in 1..c {
                        emit(
                            &[corner(k), on_opposite(j - 1), on_opposite(j)],
                            sides[k1].input_face(x[k1] + j),
                        )?;
                    }
                    // From the last segment of the fan back along side k − 1 to corner k.
                    let mut last = vec![corner(k), on_opposite(c - 1)];
                    last.extend(cut_corner(k2));
                    emit(&last, sides[k2].input_face(n[k2]))?;
                }
            }
        }
        debug_assert!(overlay.face_input.iter().all(|&g| g != NONE));
        Ok(overlay)
    }
}
