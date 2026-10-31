//! The common subdivision of the input mesh and the intrinsic triangulation, the interpolation
//! matrices `P₀` and `P₁`, and point queries from the input mesh to the triangulation
//! (SPEC §4.8, §4.9; paper Sec. 3.7, Sec. 4.3, Fig. 11, Supp. Alg. 9).
//!
//! # Construction (Supp. Alg. 9)
//!
//! The vertices of the subdivision are the vertices of `T¹` (with their ids) followed by one
//! vertex per crossing of an edge of `T¹` with an input edge, in the order of the edges and along
//! each edge from the tail of its stored halfedge: exactly the vertex order of the signpost
//! crate's [`Overlay`]. The crossings and their positions come from [`extract_all_edges`]
//! (Alg. 4) and the strip intersections of Alg. 2 for every input edge: a crossing
//! `(ij→, p)` of the curve of input halfedge `h` is the crossing number `p` (or `n − 1 − p` when
//! `ij→` is the twin of the stored halfedge) along its edge.
//!
//! The faces are made per face of `T¹` from the integers alone. With `c_k` (Eq. 1) the arcs that
//! cut corner `k` and `e_k` (Eq. 2) the arcs that leave corner `k` for the opposite side:
//!
//! * the `m`-th arc at corner `k` joins crossing `m` of side `k` to crossing `n_{k+2} − 1 − m` of
//!   side `k + 2` (counted from the starts of the sides); it cuts off the triangle or quadrilateral
//!   between itself and the previous arc (or the corner);
//! * case 1 (`e = 0`, Fig. 11 left): what is left is one polygon with up to six corners, one per
//!   end of the innermost arc at every corner (or the corner itself);
//! * case 2 (`e_κ > 0`, Fig. 11 right): the arcs from corner `κ` to side `κ + 1` fan the rest
//!   into `e_κ + 1` pieces, the first and the last bounded by the innermost arcs of corners
//!   `κ + 1` and `κ + 2`.
//!
//! The arcs are determined by the three side counts, so this is Alg. 9's per-face pairing without
//! first rotating the face so that `n_ij` is the largest (the pairing is the same; the oracle's
//! task `subdivide` states it). The input face of a piece is read off the crossings that bound it:
//! the geometric crossing knows on which side of its input edge the two parts of the `T¹` edge
//! lie. A face whose arcs fail the pairing (an odd numerator, a bound that does not add up, two
//! crossings joined by an arc that belong to different input edges, a fan arc whose input edge
//! misses the fan's corner) is an [`IntegerError::InvalidCoordinates`] naming the face, which
//! [`certify`](IntegerTriangulation::certify) catches earlier.
//!
//! # No cache
//!
//! The reference caches its subdivision and forgets to drop it after a face split (SPEC I19). Here
//! it is computed on request, and the result carries the [`fingerprint`](IntegerTriangulation::fingerprint)
//! of the state it was made from; [`IntegerTriangulation::interpolation_matrices`] refuses a
//! subdivision of another state.
//!
//! # Interpolation matrices (Sec. 4.3)
//!
//! `P₁` (`|V^S| × |V¹|`) writes every vertex of the subdivision as an affine combination of the
//! vertices of `T¹`: one entry for a vertex of `T¹`, two for a crossing (its end points on the
//! edge). `P₀` (`|V^S| × |V⁰|`) does the same with the vertices of `T⁰`: one entry for an input
//! vertex, two for a vertex on an input edge (an inserted vertex there, or a crossing), three for
//! an inserted vertex inside an input face. Rows sum to one.
//!
//! # Point queries `T⁰ → T¹` (SPEC §4.9, I10)
//!
//! [`PointLocator`] finds the pieces of an input face when asked, and keeps them. The reference
//! of SPEC I10 builds them from the three edges of the input face; that does not see the edges of
//! `T¹` that run inside the face between inserted vertices, so the pieces are found by a search
//! over the faces of `T¹` instead, from the first face along the first edge of the input face
//! (the face holding the start of `ExtractEdge`), across every side of a piece that is not
//! a part of an input edge. Each face of `T¹` met is cut as above, from
//! [`extract_intrinsic_edge`](IntegerTriangulation::extract_intrinsic_edge) of its three
//! sides. Cost: proportional to the crossings of the edges of the faces of `T¹` that meet the
//! input face.
//!
//! [`extract_all_edges`]: IntegerTriangulation::extract_all_edges

use std::collections::{HashMap, HashSet};

use ddg_core::SparseMatrix;
use ddg_core::topology::NONE;
use ddg_intrinsic::{EdgeTraceOptions, IntrinsicError, Overlay, OverlayVertex, SurfacePoint};

use crate::coords::IntegerTriangulation;
use crate::error::{IntegerError, Place, ValidityItem};
use crate::formulas::FaceCoordinates;
use crate::geometry::StripReport;
use crate::interface::to_intrinsic_error;
use crate::trace::PathPiece;

/// The common subdivision of a state, with the fingerprint of that state.
#[derive(Debug, Clone, PartialEq)]
pub struct CommonSubdivision {
    /// The subdivision, as the signpost crate's type: vertices (those of `T¹`, then the
    /// crossings), convex pieces with their input face and their face of `T¹`.
    pub overlay: Overlay,
    /// [`IntegerTriangulation::fingerprint`] of the state it was made from (SPEC I19).
    pub fingerprint: u64,
    /// What the construction did.
    pub report: SubdivisionReport,
}

/// What [`IntegerTriangulation::common_subdivision`] did.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SubdivisionReport {
    /// Vertices of the subdivision (`|V¹| + Σ n⁺`).
    pub n_vertices: usize,
    /// Pieces.
    pub n_faces: usize,
    /// Crossings, `Σ n⁺`.
    pub n_crossings: usize,
    /// Clamps and the length mismatch of the strip intersections of Alg. 2 (SPEC I5).
    pub strip: StripReport,
}

/// A crossing of an edge of `T¹` with an input edge, along the edge of `T¹` (in the direction of
/// its stored halfedge).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CrossInfo {
    /// The input edge, [`NONE`] while unset.
    input_edge: u32,
    /// The input face in which the edge lies before the crossing and after it.
    before: u32,
    after: u32,
    /// The parameter along the input edge, in its own direction.
    u: f64,
    /// The parameter along the edge of `T¹`, from its tail.
    v: f64,
}

const UNSET: CrossInfo = CrossInfo {
    input_edge: NONE,
    before: NONE,
    after: NONE,
    u: 0.0,
    v: 0.0,
};

/// The crossings on one side of a face, seen in the direction of that side.
#[derive(Clone, Copy)]
struct SideView<'a> {
    forward: bool,
    x: &'a [CrossInfo],
}

impl SideView<'_> {
    fn n(&self) -> usize {
        self.x.len()
    }
    /// The `i`-th crossing from the start of the side.
    fn crossing(&self, i: usize) -> &CrossInfo {
        &self.x[if self.forward { i } else { self.n() - 1 - i }]
    }
    /// The parameter of the `i`-th crossing along the side, from its start.
    fn along(&self, i: usize) -> f64 {
        let v = self.crossing(i).v;
        if self.forward { v } else { 1.0 - v }
    }
    /// The input face of the piece of the side between its crossings `i − 1` and `i` (`i = 0`:
    /// from the start; `i = n`: to the end). The side has a crossing.
    fn face_at(&self, i: usize) -> u32 {
        let j = if self.forward { i } else { self.n() - i };
        if j < self.n() {
            self.x[j].before
        } else {
            self.x[j - 1].after
        }
    }
}

/// A corner of a piece: a corner of the face, or the `i`-th crossing of a side from the side's
/// start.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Node {
    Corner(u8),
    Cross { side: u8, i: u32 },
}

/// A piece of a face of `T¹`: its corners counter-clockwise, and a piece of a side that bounds it
/// (side, index of the piece along it) which gives its input face, if any side has crossings.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Piece {
    nodes: [Node; 6],
    len: usize,
    on_side: Option<(u8, u32)>,
}

impl Piece {
    fn new(on_side: Option<(u8, u32)>) -> Piece {
        Piece {
            nodes: [Node::Corner(0); 6],
            len: 0,
            on_side,
        }
    }
    fn push(&mut self, node: Node) {
        self.nodes[self.len] = node;
        self.len += 1;
    }
    pub(crate) fn nodes(&self) -> &[Node] {
        &self.nodes[..self.len]
    }
}

/// Cuts a face with the side counts of `fc` into its pieces, calling `emit` for each: the pieces
/// at corner 0, 1, 2 from the corner outwards, then the middle piece (case 1) or the fan
/// (case 2). `fc` must satisfy the identities checked by `face_pieces`.
#[allow(clippy::needless_range_loop)]
fn cut_face(
    fc: &FaceCoordinates,
    mut emit: impl FnMut(&Piece) -> Result<(), IntegerError>,
) -> Result<(), IntegerError> {
    let (n, c, e) = (fc.n, fc.c, fc.e);
    let fan = (0..3).filter(|&k| e[k] > 0).collect::<Vec<_>>();
    // The last end of the arc `m − 1` at corner k, and the piece at corner k when its x[k]
    // arcs are cut off: the corner itself or the two ends of the innermost arc.
    let inner = |k: usize, piece: &mut Piece| {
        if c[k] == 0 {
            piece.push(Node::Corner(k as u8));
        } else {
            let k2 = (k + 2) % 3;
            piece.push(Node::Cross {
                side: k2 as u8,
                i: (n[k2] - c[k]) as u32,
            });
            piece.push(Node::Cross {
                side: k as u8,
                i: (c[k] - 1) as u32,
            });
        }
    };
    // The pieces at the corners.
    for k in 0..3 {
        let k2 = (k + 2) % 3;
        for m in 0..c[k] {
            let (a, b) = (
                Node::Cross {
                    side: k as u8,
                    i: m as u32,
                },
                Node::Cross {
                    side: k2 as u8,
                    i: (n[k2] - 1 - m) as u32,
                },
            );
            let mut piece = Piece::new(Some((k as u8, m as u32)));
            if m == 0 {
                piece.push(Node::Corner(k as u8));
                piece.push(a);
                piece.push(b);
            } else {
                piece.push(Node::Cross {
                    side: k as u8,
                    i: (m - 1) as u32,
                });
                piece.push(a);
                piece.push(b);
                piece.push(Node::Cross {
                    side: k2 as u8,
                    i: (n[k2] - m) as u32,
                });
            }
            emit(&piece)?;
        }
    }
    match fan.first() {
        None => {
            let on_side = (0..3).find(|&k| n[k] > 0).map(|k| (k as u8, c[k] as u32));
            let mut piece = Piece::new(on_side);
            for k in 0..3 {
                inner(k, &mut piece);
            }
            emit(&piece)?;
        }
        Some(&kappa) => {
            let (k1, k2) = ((kappa + 1) % 3, (kappa + 2) % 3);
            let on = |j: i64| Node::Cross {
                side: k1 as u8,
                i: (c[k1] + j) as u32,
            };
            // From the corner along the side to the first arc of the fan.
            let mut first = Piece::new(Some((k1 as u8, c[k1] as u32)));
            first.push(Node::Corner(kappa as u8));
            inner(k1, &mut first);
            first.push(on(0));
            emit(&first)?;
            for j in 1..e[kappa] {
                let mut piece = Piece::new(Some((k1 as u8, (c[k1] + j) as u32)));
                piece.push(Node::Corner(kappa as u8));
                piece.push(on(j - 1));
                piece.push(on(j));
                emit(&piece)?;
            }
            let mut last = Piece::new(Some((k1 as u8, (c[k1] + e[kappa]) as u32)));
            last.push(Node::Corner(kappa as u8));
            last.push(on(e[kappa] - 1));
            inner(k2, &mut last);
            emit(&last)?;
        }
    }
    Ok(())
}

fn bad_face(f: usize, detail: String) -> IntegerError {
    IntegerError::invalid(ValidityItem::Extraction, Place::Face(f as u32), detail)
}

impl IntegerTriangulation {
    /// A hash of the state: the faces, the length of every halfedge's edge, the normal
    /// coordinates, the roundabouts and the places of the vertices (FNV-1a over 64-bit words, as
    /// IT-SPEC N15). A [`CommonSubdivision`] carries it; any flip, split, insertion or removal
    /// changes it (up to hash collisions).
    pub fn fingerprint(&self) -> u64 {
        let mut hash = 0xcbf2_9ce4_8422_2325u64;
        let mut word = |w: u64| {
            for byte in w.to_le_bytes() {
                hash ^= u64::from(byte);
                hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        word(self.complex.n_vertices() as u64);
        word(self.complex.n_faces() as u64);
        for (f, face) in self.complex.faces().iter().enumerate() {
            for (k, &vertex) in face.iter().enumerate() {
                word(u64::from(vertex));
                word(self.lengths[self.complex.edge(3 * f + k)].to_bits());
            }
        }
        for &n in &self.normal {
            word(n as u64);
        }
        for &r in &self.roundabouts {
            word(u64::from(r));
        }
        for location in &self.locations {
            match *location {
                SurfacePoint::Vertex(v) => {
                    word(0);
                    word(u64::from(v));
                }
                SurfacePoint::Edge { edge, t } => {
                    word(1);
                    word(u64::from(edge));
                    word(t.to_bits());
                }
                SurfacePoint::Face { face, bary } => {
                    word(2);
                    word(u64::from(face));
                    for b in bary {
                        word(b.to_bits());
                    }
                }
            }
        }
        hash
    }

    /// Whether vertex `v` of `T¹` lies on input edge `input_edge` (as an end point or inside).
    fn vertex_on_input_edge(&self, v: usize, input_edge: u32) -> bool {
        match self.locations[v] {
            SurfacePoint::Vertex(u) => self
                .input
                .mesh()
                .edge_vertices(input_edge as usize)
                .contains(&(u as usize)),
            SurfacePoint::Edge { edge, .. } => edge == input_edge,
            SurfacePoint::Face { .. } => false,
        }
    }

    /// Cuts face `f` of `T¹` into its pieces (Supp. Alg. 9, per face; see the module
    /// documentation) and calls `emit` for each, in the order: the pieces at corner 0, 1, 2 from
    /// the corner outwards, then the middle piece (case 1) or the fan (case 2).
    ///
    /// The crossings of the sides are given in `sides`; their count must be `n⁺` of the side.
    fn face_pieces(
        &self,
        f: usize,
        sides: &[SideView<'_>; 3],
        emit: impl FnMut(&Piece) -> Result<(), IntegerError>,
    ) -> Result<(), IntegerError> {
        let counts = [0, 1, 2].map(|k| sides[k].n() as i64);
        let fc = FaceCoordinates::of(counts).map_err(|k| {
            bad_face(
                f,
                format!(
                    "the sides have {}, {} and {} crossings, and the numerator of Eq. 1 at \
                     corner {k} is odd",
                    counts[0], counts[1], counts[2]
                ),
            )
        })?;
        let (n, c, e) = (fc.n, fc.c, fc.e);
        let usz = |x: i64| x as usize;
        for k in 0..3 {
            if c[k] < 0 || e[k] < 0 || n[k] != c[k] + c[(k + 1) % 3] + e[(k + 2) % 3] {
                return Err(bad_face(
                    f,
                    format!(
                        "side {k} has {} crossings but the corners cut off {}, {} and \
                         {} arcs from the opposite corner",
                        n[k],
                        c[k],
                        c[(k + 1) % 3],
                        e[(k + 2) % 3]
                    ),
                ));
            }
        }
        let fan = (0..3).filter(|&k| e[k] > 0).collect::<Vec<_>>();
        if fan.len() > 1 {
            return Err(bad_face(
                f,
                format!("corners {fan:?} all emanate arcs, which would cross"),
            ));
        }
        let corner_vertex = |k: usize| self.complex.faces()[f][k] as usize;
        // The pairing of the arcs that cut corners.
        for k in 0..3 {
            let k2 = (k + 2) % 3;
            for m in 0..usz(c[k]) {
                let (a, b) = (
                    sides[k].crossing(m).input_edge,
                    sides[k2].crossing(usz(n[k2]) - 1 - m).input_edge,
                );
                if a != b {
                    return Err(bad_face(
                        f,
                        format!(
                            "arc {m} at corner {k} joins crossings of the input edges {a} and {b}"
                        ),
                    ));
                }
            }
        }
        if let Some(&kappa) = fan.first() {
            let k1 = (kappa + 1) % 3;
            for j in 0..e[kappa] {
                let ie = sides[k1].crossing(usz(c[k1] + j)).input_edge;
                if !self.vertex_on_input_edge(corner_vertex(kappa), ie) {
                    return Err(bad_face(
                        f,
                        format!(
                            "input edge {ie} crosses the side opposite vertex {} and no other \
                             side, but does not pass through that vertex",
                            corner_vertex(kappa)
                        ),
                    ));
                }
            }
        }
        cut_face(&fc, emit)
    }

    /// The input face that contains a side of a face of `T¹` with no crossings: the face on the
    /// left of the input edge it runs along, or the only face it crosses.
    fn whole_side_face(&self, f: usize, k: usize) -> Result<u32, IntegerError> {
        let h = 3 * f + k;
        let trace =
            self.extract_intrinsic_edge(self.complex.edge(h), &EdgeTraceOptions::default())?;
        let mesh = self.input.mesh();
        if let Some(ie) = trace.along_input_edge {
            let (tail, head) = (self.complex.from_vertex(h), self.complex.to_vertex(h));
            let ascending = self.parameter_on_input_edge(tail, ie as usize)?
                < self.parameter_on_input_edge(head, ie as usize)?;
            let [lo, hi] = mesh.edge_vertices(ie as usize);
            let (from, to) = if ascending { (lo, hi) } else { (hi, lo) };
            mesh.find_halfedge(from, to)
                .map(|g| (g / 3) as u32)
                .ok_or_else(|| {
                    bad_face(
                        f,
                        format!(
                            "side {k} runs along input edge {ie}, which has no face on its left"
                        ),
                    )
                })
        } else {
            trace
                .first_face
                .ok_or_else(|| bad_face(f, format!("side {k} lies in no input face")))
        }
    }

    /// The common subdivision of the input mesh and `T¹`, from the integers (SPEC §4.8,
    /// Supp. Alg. 9; see the [module documentation](self)).
    ///
    /// The state is assumed valid; [`certify`](Self::certify) first if it comes from outside.
    ///
    /// # Errors
    /// Those of [`extract_all_edges`](Self::extract_all_edges) and
    /// [`segment_geometry`](Self::segment_geometry); [`IntegerError::InvalidCoordinates`] naming
    /// the face whose arcs cannot be paired, or the edge with a crossing that no input edge
    /// passes through (or two).
    pub fn common_subdivision(&self) -> Result<CommonSubdivision, IntegerError> {
        let c = &self.complex;
        let mesh = self.input.mesh();
        let (nv, ne) = (c.n_vertices(), c.n_edges());
        let mut offsets = Vec::with_capacity(ne + 1);
        offsets.push(0usize);
        for e in 0..ne {
            let n = self.normal[e];
            if n < -1 {
                return Err(IntegerError::invalid(
                    ValidityItem::Ranges,
                    Place::Edge(e as u32),
                    format!("normal coordinate {n} is below −1"),
                ));
            }
            let total = offsets[e] as u128 + n.max(0) as u128;
            if total + nv as u128 >= u128::from(u32::MAX) {
                return Err(IntegerError::TooManyCrossings {
                    total,
                    limit: u64::from(u32::MAX),
                });
            }
            offsets.push(total as usize);
        }
        let n_crossings = offsets[ne];

        // The crossings, from the curves of the input edges (Alg. 4 and Alg. 2).
        let mut table = vec![UNSET; n_crossings];
        let mut strip = StripReport::default();
        let paths = self.extract_all_edges()?;
        for (ie, path) in paths.iter().enumerate() {
            let h_in = mesh.edge_halfedge(ie);
            let (right, left) = (
                mesh.twin(h_in).map_or(NONE, |t| (t / 3) as u32),
                (h_in / 3) as u32,
            );
            let vertex_u = path
                .vertices
                .iter()
                .map(|&v| self.parameter_on_input_edge(v as usize, ie))
                .collect::<Result<Vec<_>, _>>()?;
            for (i, piece) in path.pieces.iter().enumerate() {
                let PathPiece::Crossings(segment) = piece else {
                    continue;
                };
                let (found, report) =
                    self.segment_geometry(segment, ie, vertex_u[i], vertex_u[i + 1])?;
                strip.merge(&report);
                for z in found {
                    let h = z.crossing.halfedge as usize;
                    let e = c.edge(h);
                    let n = self.normal[e].max(0);
                    let forward = c.is_forward(h);
                    let p = if forward {
                        z.crossing.index
                    } else {
                        n - 1 - z.crossing.index
                    };
                    if p < 0 || p >= n {
                        return Err(IntegerError::invalid(
                            ValidityItem::Extraction,
                            z.crossing.place(),
                            format!(
                                "the crossing of input edge {ie} is not among the {n} of edge {e}"
                            ),
                        ));
                    }
                    let slot = &mut table[offsets[e] + p as usize];
                    if slot.input_edge != NONE {
                        return Err(IntegerError::invalid(
                            ValidityItem::Extraction,
                            z.crossing.place(),
                            format!(
                                "crossing {p} of edge {e} lies on input edges {} and {ie}",
                                slot.input_edge
                            ),
                        ));
                    }
                    // The halfedge of the triangulation points to the left of the input edge's
                    // direction (the curve leaves its face through it): its tail is on the right.
                    let (before, after) = if forward {
                        (right, left)
                    } else {
                        (left, right)
                    };
                    *slot = CrossInfo {
                        input_edge: ie as u32,
                        before,
                        after,
                        u: z.u,
                        v: if forward { z.v } else { 1.0 - z.v },
                    };
                }
            }
        }
        for e in 0..ne {
            let x = &table[offsets[e]..offsets[e + 1]];
            if let Some(p) = x.iter().position(|z| z.input_edge == NONE) {
                return Err(IntegerError::invalid(
                    ValidityItem::Extraction,
                    Place::Edge(e as u32),
                    format!(
                        "crossing {p} of the {} of the edge lies on no input edge",
                        x.len()
                    ),
                ));
            }
            if let Some(p) = x.windows(2).position(|w| w[0].after != w[1].before) {
                return Err(IntegerError::invalid(
                    ValidityItem::Extraction,
                    Place::Edge(e as u32),
                    format!(
                        "between its crossings {p} and {} the edge is in input faces {} and {}",
                        p + 1,
                        x[p].after,
                        x[p + 1].before
                    ),
                ));
            }
        }

        // Vertices.
        let mut vertices: Vec<OverlayVertex> = (0..nv)
            .map(|v| OverlayVertex {
                on_input: self.locations[v],
                on_triangulation: SurfacePoint::Vertex(v as u32),
            })
            .collect();
        vertices.reserve(n_crossings);
        for e in 0..ne {
            for z in &table[offsets[e]..offsets[e + 1]] {
                vertices.push(OverlayVertex {
                    on_input: SurfacePoint::Edge {
                        edge: z.input_edge,
                        t: z.u,
                    },
                    on_triangulation: SurfacePoint::Edge {
                        edge: e as u32,
                        t: z.v,
                    },
                });
            }
        }

        // Faces.
        let mut overlay = Overlay {
            vertices,
            face_starts: vec![0],
            corners: Vec::new(),
            face_input: Vec::new(),
            face_triangulation: Vec::new(),
        };
        for f in 0..c.n_faces() {
            let sides: [SideView<'_>; 3] = [0, 1, 2].map(|k| {
                let e = c.edge(3 * f + k);
                SideView {
                    forward: c.is_forward(3 * f + k),
                    x: &table[offsets[e]..offsets[e + 1]],
                }
            });
            let id = |node: Node| -> u32 {
                match node {
                    Node::Corner(k) => c.faces()[f][k as usize],
                    Node::Cross { side, i } => {
                        let side = side as usize;
                        let e = c.edge(3 * f + side);
                        let n = sides[side].n();
                        let along = if sides[side].forward {
                            i as usize
                        } else {
                            n - 1 - i as usize
                        };
                        (nv + offsets[e] + along) as u32
                    }
                }
            };
            self.face_pieces(f, &sides, |piece| {
                let input_face = match piece.on_side {
                    Some((side, i)) if sides[side as usize].n() > 0 => {
                        sides[side as usize].face_at(i as usize)
                    }
                    _ => self.whole_side_face(f, 0)?,
                };
                overlay
                    .corners
                    .extend(piece.nodes().iter().map(|&nd| id(nd)));
                overlay.face_starts.push(overlay.corners.len() as u32);
                overlay.face_input.push(input_face);
                overlay.face_triangulation.push(f as u32);
                Ok(())
            })?;
        }
        let report = SubdivisionReport {
            n_vertices: overlay.vertices.len(),
            n_faces: overlay.n_faces(),
            n_crossings,
            strip,
        };
        Ok(CommonSubdivision {
            overlay,
            fingerprint: self.fingerprint(),
            report,
        })
    }

    /// The interpolation matrices of a subdivision of this state (Sec. 4.3): `P₀`
    /// (`|V^S| × |V⁰|`) and `P₁` (`|V^S| × |V¹|`) with rows of one to three entries summing to
    /// one, such that for a function `x⁰` on the vertices of `T⁰` (linear on its faces), `P₀ x⁰`
    /// is its value at the vertices of the subdivision; and for `x¹` on the vertices of `T¹`,
    /// `P₁ x¹` is the function linear on the faces of `T¹` with those values. See the
    /// [module documentation](self) for the rows.
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] ([`IntrinsicError::InvalidOptions`]) if the subdivision was
    /// made from another state (its fingerprint differs: SPEC I19) or does not have the vertices
    /// of this one's.
    pub fn interpolation_matrices(
        &self,
        subdivision: &CommonSubdivision,
    ) -> Result<(SparseMatrix, SparseMatrix), IntegerError> {
        let stale =
            |reason: String| IntegerError::Intrinsic(IntrinsicError::InvalidOptions { reason });
        if subdivision.fingerprint != self.fingerprint() {
            return Err(stale(
                "the subdivision was made from another state (the fingerprints differ)".into(),
            ));
        }
        let vertices = &subdivision.overlay.vertices;
        if vertices.len() < self.complex.n_vertices() {
            return Err(stale(format!(
                "the subdivision has {} vertices, fewer than the {} of the triangulation",
                vertices.len(),
                self.complex.n_vertices()
            )));
        }
        let mesh = self.input.mesh();
        let (n0, n1) = (self.n_input_vertices(), self.complex.n_vertices());
        let mut t0 = Vec::with_capacity(2 * vertices.len());
        let mut t1 = Vec::with_capacity(2 * vertices.len());
        let out_of_range = |what: &'static str, index: usize, bound: usize| {
            IntegerError::IndexOutOfRange { what, index, bound }
        };
        for (s, vertex) in vertices.iter().enumerate() {
            let row = s as u32;
            // P0 from the place on the input mesh.
            match vertex.on_input {
                SurfacePoint::Vertex(v) => t0.push((row, v, 1.0)),
                SurfacePoint::Edge { edge, t } => {
                    if edge as usize >= mesh.n_edges() {
                        return Err(out_of_range("input edge", edge as usize, mesh.n_edges()));
                    }
                    let [a, b] = mesh.edge_vertices(edge as usize);
                    t0.push((row, a as u32, 1.0 - t));
                    t0.push((row, b as u32, t));
                }
                SurfacePoint::Face { face, bary } => {
                    if face as usize >= mesh.n_faces() {
                        return Err(out_of_range("input face", face as usize, mesh.n_faces()));
                    }
                    for (k, &v) in mesh.faces()[face as usize].iter().enumerate() {
                        t0.push((row, v, bary[k]));
                    }
                }
            }
            // P1 from the place on the triangulation.
            match vertex.on_triangulation {
                SurfacePoint::Vertex(v) => t1.push((row, v, 1.0)),
                SurfacePoint::Edge { edge, t } => {
                    if edge as usize >= self.complex.n_edges() {
                        return Err(out_of_range("edge", edge as usize, self.complex.n_edges()));
                    }
                    let [a, b] = self.complex.edge_vertices(edge as usize);
                    t1.push((row, a as u32, 1.0 - t));
                    t1.push((row, b as u32, t));
                }
                SurfacePoint::Face { face, bary } => {
                    if face as usize >= self.complex.n_faces() {
                        return Err(out_of_range("face", face as usize, self.complex.n_faces()));
                    }
                    for (k, &v) in self.complex.faces()[face as usize].iter().enumerate() {
                        t1.push((row, v, bary[k]));
                    }
                }
            }
        }
        let build = |rows, cols, triplets: &[(u32, u32, f64)]| {
            SparseMatrix::from_triplets(rows, cols, triplets)
                .map_err(|e| stale(format!("the interpolation matrix could not be built: {e}")))
        };
        Ok((
            build(vertices.len(), n0, &t0)?,
            build(vertices.len(), n1, &t1)?,
        ))
    }

    /// A [`PointLocator`] over this state.
    pub fn point_locator(&self) -> PointLocator<'_> {
        PointLocator {
            tri: self,
            sides: HashMap::new(),
            input_faces: HashMap::new(),
        }
    }

    /// The place on `T¹` of a point of the input mesh (SPEC §4.9 second half): the same
    /// signature and meaning as `ddg_intrinsic::IntrinsicTriangulation::point_on_triangulation`.
    /// Each call builds the pieces of one input face; ask a [`point_locator`](Self::point_locator)
    /// for many points.
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`] for a point that is not one of the input mesh; the errors
    /// of [`PointLocator::point_on_triangulation`] otherwise, as [`IntrinsicError::Representation`].
    pub fn point_on_triangulation(
        &self,
        point: SurfacePoint,
    ) -> Result<SurfacePoint, IntrinsicError> {
        self.point_locator().point_on_triangulation(point)
    }
}

// ----- point queries -----

/// A coordinate below this, in the barycentric coordinates of a result in a face of `T¹`, is
/// taken to be zero: the result is then a point of the edge or the vertex (the signposts snap
/// the same way, `1e-12` in IT-SPEC N29).
pub const LOCATE_SNAP: f64 = 1e-12;

/// How far outside every piece of the input face a query may fall (in barycentric coordinates of
/// the piece) before it is refused.
const LOCATE_SLACK: f64 = 1e-9;

/// A piece of an input face: its corners as barycentric coordinates in the input face and in
/// the face of `T¹` it lies in, counter-clockwise.
#[derive(Debug, Clone, PartialEq)]
pub struct LocatedPiece {
    /// The face of `T¹`.
    pub face: u32,
    /// The corners: (coordinates in the input face, coordinates in the face of `T¹`).
    pub corners: Vec<([f64; 3], [f64; 3])>,
}

/// Finds the place on `T¹` of points of the input mesh, building the pieces of each input face
/// when it is first asked for and keeping them (SPEC §4.9, I10). It borrows the state, so it
/// cannot outlive a change of it.
#[derive(Debug)]
pub struct PointLocator<'a> {
    tri: &'a IntegerTriangulation,
    /// The crossings of the edges of `T¹` met so far (`ExtractIntrinsicEdge`).
    sides: HashMap<u32, Vec<CrossInfo>>,
    input_faces: HashMap<u32, Vec<LocatedPiece>>,
}

impl PointLocator<'_> {
    /// The crossings of edge `e` of `T¹` along it, from `extract_intrinsic_edge`.
    fn edge_crossings(&mut self, e: usize) -> Result<&[CrossInfo], IntegerError> {
        if !self.sides.contains_key(&(e as u32)) {
            let tri = self.tri;
            let trace = tri.extract_intrinsic_edge(e, &EdgeTraceOptions::default())?;
            let mesh = tri.input.mesh();
            let length = tri.lengths[e];
            let mut x = Vec::with_capacity(trace.crossings.len());
            for z in &trace.crossings {
                let g = z.halfedge as usize;
                let Some(twin) = mesh.twin(g) else {
                    return Err(IntegerError::invalid(
                        ValidityItem::Extraction,
                        Place::Edge(e as u32),
                        format!(
                            "the edge leaves input face {} by the boundary halfedge {g}",
                            g / 3
                        ),
                    ));
                };
                x.push(CrossInfo {
                    input_edge: z.edge,
                    before: (g / 3) as u32,
                    after: (twin / 3) as u32,
                    u: z.t,
                    v: (z.distance / length).clamp(0.0, 1.0),
                });
            }
            self.sides.insert(e as u32, x);
        }
        Ok(&self.sides[&(e as u32)])
    }

    /// The pieces of input face `g` (found when first asked for; see the module documentation).
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; those of
    /// [`extract_intrinsic_edge`](IntegerTriangulation::extract_intrinsic_edge) and
    /// [`face_pieces`](IntegerTriangulation); [`IntegerError::InvalidCoordinates`] naming the
    /// face if the search leaves the pieces' faces inconsistent.
    pub fn pieces_of_input_face(&mut self, g: usize) -> Result<&[LocatedPiece], IntegerError> {
        let tri = self.tri;
        let mesh = tri.input.mesh();
        if g >= mesh.n_faces() {
            return Err(IntegerError::IndexOutOfRange {
                what: "input face",
                index: g,
                bound: mesh.n_faces(),
            });
        }
        if !self.input_faces.contains_key(&(g as u32)) {
            let pieces = self.build_input_face(g)?;
            self.input_faces.insert(g as u32, pieces);
        }
        Ok(&self.input_faces[&(g as u32)])
    }

    /// The face of `T¹` that holds the start of the curve of input halfedge `3g`, on the side of
    /// `g`.
    fn seed(&self, g: usize) -> Result<usize, IntegerError> {
        let tri = self.tri;
        let c = &tri.complex;
        let path = tri.extract_edge(3 * g)?;
        let fail = |detail: String| {
            IntegerError::invalid(
                ValidityItem::Extraction,
                Place::InputEdge(tri.input.mesh().edge(3 * g) as u32),
                detail,
            )
        };
        match path.pieces.first() {
            Some(PathPiece::Crossings(segment)) => segment
                .crossings
                .first()
                .map(|z| z.halfedge as usize / 3)
                .ok_or_else(|| fail("the first piece of the curve has no crossing".into())),
            Some(PathPiece::Along { edge }) => {
                let h = c.edge_halfedge(*edge as usize);
                let (from, to) = (path.vertices[0] as usize, path.vertices[1] as usize);
                let along = if c.from_vertex(h) == from && c.to_vertex(h) == to {
                    Some(h)
                } else {
                    c.twin(h)
                };
                along
                    .map(|h| h / 3)
                    .ok_or_else(|| fail("the edge along it has no face on its left".into()))
            }
            None => Err(fail("the curve has no piece".into())),
        }
    }

    fn build_input_face(&mut self, g: usize) -> Result<Vec<LocatedPiece>, IntegerError> {
        let tri = self.tri;
        let c = &tri.complex;
        let mesh = tri.input.mesh();
        let mut seen: HashSet<usize> = HashSet::new();
        let mut stack = vec![self.seed(g)?];
        seen.insert(stack[0]);
        let mut out: Vec<LocatedPiece> = Vec::new();
        while let Some(f) = stack.pop() {
            // The crossings of the sides (owned copies: the cache is borrowed mutably).
            let mut owned: Vec<Vec<CrossInfo>> = Vec::with_capacity(3);
            for k in 0..3 {
                owned.push(self.edge_crossings(c.edge(3 * f + k))?.to_vec());
            }
            let sides: [SideView<'_>; 3] = [0, 1, 2].map(|k| SideView {
                forward: c.is_forward(3 * f + k),
                x: &owned[k],
            });
            // The place of a node in the input face and in the face of T¹.
            let place = |node: Node| -> Result<([f64; 3], [f64; 3]), IntegerError> {
                match node {
                    Node::Corner(k) => {
                        let mut mu = [0.0; 3];
                        mu[k as usize] = 1.0;
                        let w = c.faces()[f][k as usize] as usize;
                        Ok((input_bary(tri, tri.locations[w], g, f)?, mu))
                    }
                    Node::Cross { side, i } => {
                        let s = side as usize;
                        let z = sides[s].crossing(i as usize);
                        let v = sides[s].along(i as usize);
                        let mut mu = [0.0; 3];
                        mu[s] = 1.0 - v;
                        mu[(s + 1) % 3] = v;
                        let on_input = SurfacePoint::Edge {
                            edge: z.input_edge,
                            t: z.u,
                        };
                        Ok((input_bary(tri, on_input, g, f)?, mu))
                    }
                }
            };
            let mut found: Vec<Piece> = Vec::new();
            let mut faces_of: Vec<u32> = Vec::new();
            tri.face_pieces(f, &sides, |piece| {
                let input_face = match piece.on_side {
                    Some((side, i)) if sides[side as usize].n() > 0 => {
                        sides[side as usize].face_at(i as usize)
                    }
                    _ => tri.whole_side_face(f, 0)?,
                };
                found.push(*piece);
                faces_of.push(input_face);
                Ok(())
            })?;
            for (piece, _) in found
                .iter()
                .zip(&faces_of)
                .filter(|&(_, &face)| face as usize == g)
            {
                let corners = piece
                    .nodes()
                    .iter()
                    .map(|&nd| place(nd))
                    .collect::<Result<Vec<_>, _>>()?;
                out.push(LocatedPiece {
                    face: f as u32,
                    corners,
                });
                // Across the sides of the face that the polygon has an edge along (and that
                // are not parts of input edges).
                let nodes = piece.nodes();
                for (a, b) in nodes.iter().zip(nodes.iter().cycle().skip(1)) {
                    let on_side = |nd: &Node, s: usize| match *nd {
                        Node::Corner(k) => k as usize == s || (k as usize + 2) % 3 == s,
                        Node::Cross { side, .. } => side as usize == s,
                    };
                    for s in 0..3 {
                        if !(on_side(a, s) && on_side(b, s)) {
                            continue;
                        }
                        // Two corners that are both on side s but lie on different sides'
                        // corners only: they are the ends of side s itself.
                        let h = 3 * f + s;
                        if tri.normal[c.edge(h)] == -1 {
                            continue;
                        }
                        if let Some(t) = c.twin(h)
                            && seen.insert(t / 3)
                        {
                            stack.push(t / 3);
                        }
                    }
                }
            }
        }
        if out.is_empty() {
            return Err(bad_face(
                self.seed(g)?,
                format!(
                    "no piece of the subdivision lies in input face {g}, whose first edge starts here"
                ),
            ));
        }
        let _ = mesh;
        Ok(out)
    }

    /// The place on `T¹` of a point of the input mesh; see
    /// [`IntegerTriangulation::point_on_triangulation`].
    ///
    /// A point on an input edge is located from the face on the left of the edge's stored
    /// halfedge. The result's coordinates below [`LOCATE_SNAP`] are taken to be zero.
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`]; [`IntrinsicError::Representation`] for the errors of
    /// the construction of the pieces, and for a point that lies in no piece of its face.
    pub fn point_on_triangulation(
        &mut self,
        point: SurfacePoint,
    ) -> Result<SurfacePoint, IntrinsicError> {
        let tri = self.tri;
        let mesh = tri.input.mesh();
        let invalid = |reason: &'static str| IntrinsicError::InvalidPoint { point, reason };
        let in_unit = |x: f64| x.is_finite() && (0.0..=1.0).contains(&x);
        let (g, lambda) = match point {
            SurfacePoint::Vertex(v) => {
                if v as usize >= mesh.n_vertices() {
                    return Err(invalid("the vertex does not exist"));
                }
                return Ok(SurfacePoint::Vertex(v));
            }
            SurfacePoint::Edge { edge, t } => {
                if edge as usize >= mesh.n_edges() {
                    return Err(invalid("the edge does not exist"));
                }
                if !in_unit(t) {
                    return Err(invalid("the parameter is not in [0, 1]"));
                }
                let [lo, hi] = mesh.edge_vertices(edge as usize);
                if t == 0.0 || t == 1.0 {
                    return Ok(SurfacePoint::Vertex(if t == 0.0 { lo } else { hi } as u32));
                }
                let g = mesh.edge_halfedge(edge as usize) / 3;
                let l = input_bary(tri, point, g, 0)
                    .map_err(|e| to_intrinsic_error(e, "point_on_triangulation"))?;
                (g, l)
            }
            SurfacePoint::Face { face, bary } => {
                if face as usize >= mesh.n_faces() {
                    return Err(invalid("the face does not exist"));
                }
                if !bary.iter().all(|&b| in_unit(b))
                    || (bary.iter().sum::<f64>() - 1.0).abs() > 1e-9
                {
                    return Err(invalid(
                        "the barycentric coordinates are not those of a point of a face",
                    ));
                }
                if let Some(k) = bary.iter().position(|&b| b == 1.0) {
                    return Ok(SurfacePoint::Vertex(mesh.faces()[face as usize][k]));
                }
                (face as usize, bary)
            }
        };
        let pieces = self
            .pieces_of_input_face(g)
            .map_err(|e| to_intrinsic_error(e, "point_on_triangulation"))?;
        // The fan triangle in which the point has the largest smallest coordinate.
        let mut best: Option<(f64, u32, [f64; 3])> = None;
        for piece in pieces {
            let p = &piece.corners;
            for k in 1..p.len() - 1 {
                let tri3 = [&p[0], &p[k], &p[k + 1]];
                let Some(s) = solve_bary(&[tri3[0].0, tri3[1].0, tri3[2].0], lambda) else {
                    continue;
                };
                let margin = s[0].min(s[1]).min(s[2]);
                if best.is_none_or(|(m, ..)| margin > m) {
                    let mut mu = [0.0; 3];
                    for (i, corner) in tri3.iter().enumerate() {
                        for (a, m) in mu.iter_mut().enumerate() {
                            *m += s[i] * corner.1[a];
                        }
                    }
                    best = Some((margin, piece.face, mu));
                }
            }
        }
        let Some((margin, face, mu)) = best.filter(|&(m, ..)| m >= -LOCATE_SLACK) else {
            return Err(IntrinsicError::Representation {
                representation: crate::interface::REPRESENTATION,
                operation: "point_on_triangulation",
                message: format!("the point lies in no piece of input face {g}"),
            });
        };
        let _ = margin;
        Ok(classify(tri, face as usize, mu))
    }
}

/// The barycentric coordinates in input face `g` of a point of the input mesh that lies in it
/// (closed); `f` names the face of `T¹` in errors.
fn input_bary(
    tri: &IntegerTriangulation,
    point: SurfacePoint,
    g: usize,
    f: usize,
) -> Result<[f64; 3], IntegerError> {
    let mesh = tri.input.mesh();
    let outside = || {
        bad_face(
            f,
            format!(
                "the point {point} is not in input face {g}, in which a piece of the face lies"
            ),
        )
    };
    let mut l = [0.0; 3];
    match point {
        SurfacePoint::Vertex(v) => {
            let k = mesh.faces()[g]
                .iter()
                .position(|&u| u == v)
                .ok_or_else(outside)?;
            l[k] = 1.0;
        }
        SurfacePoint::Edge { edge, t } => {
            let k = (0..3)
                .find(|&k| mesh.edge(3 * g + k) == edge as usize)
                .ok_or_else(outside)?;
            let forward = mesh.from_vertex(3 * g + k) == mesh.edge_vertices(edge as usize)[0];
            let th = if forward { t } else { 1.0 - t };
            l[k] = 1.0 - th;
            l[(k + 1) % 3] = th;
        }
        SurfacePoint::Face { face, bary } => {
            if face as usize != g {
                return Err(outside());
            }
            l = bary;
        }
    }
    Ok(l)
}

/// The coefficients `s` with `Σ s_i p_i = λ`, `Σ s_i = 1` (all points having coordinate sum
/// one), or `None` for a degenerate triangle.
fn solve_bary(p: &[[f64; 3]; 3], lambda: [f64; 3]) -> Option<[f64; 3]> {
    // Use the first two coordinates and the sum.
    let det3 = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let a = [
        [p[0][0], p[1][0], p[2][0]],
        [p[0][1], p[1][1], p[2][1]],
        [1.0, 1.0, 1.0],
    ];
    let rhs = [lambda[0], lambda[1], 1.0];
    let d = det3(a);
    if d.abs() < 1e-14 {
        return None;
    }
    let mut s = [0.0; 3];
    for (j, sj) in s.iter_mut().enumerate() {
        let mut m = a;
        for (i, row) in m.iter_mut().enumerate() {
            row[j] = rhs[i];
        }
        *sj = det3(m) / d;
    }
    Some(s)
}

/// A point of face `f` of `T¹` with barycentric coordinates `mu`, as a vertex, a point of an
/// edge or a point of the face, taking coordinates below [`LOCATE_SNAP`] to be zero.
fn classify(tri: &IntegerTriangulation, f: usize, mu: [f64; 3]) -> SurfacePoint {
    let c = &tri.complex;
    let zero = [0, 1, 2].map(|k| mu[k] <= LOCATE_SNAP);
    let n_zero = zero.iter().filter(|&&z| z).count();
    match n_zero {
        0 => SurfacePoint::Face {
            face: f as u32,
            bary: mu,
        },
        1 => {
            // The zero coordinate is that of the corner opposite side s.
            let z = zero.iter().position(|&z| z).unwrap_or(0);
            let s = (z + 1) % 3;
            let h = 3 * f + s;
            let v = mu[(s + 1) % 3] / (mu[s] + mu[(s + 1) % 3]);
            SurfacePoint::Edge {
                edge: c.edge(h) as u32,
                t: if c.is_forward(h) { v } else { 1.0 - v },
            }
        }
        _ => {
            let k = (0..3).max_by(|&a, &b| mu[a].total_cmp(&mu[b])).unwrap_or(0);
            SurfacePoint::Vertex(c.faces()[f][k])
        }
    }
}

/// Where the arc through a crossing ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArcEnd {
    /// At crossing `p` of side `k` (from the start of the side).
    Side {
        /// The side.
        k: u8,
        /// The crossing.
        p: u32,
    },
    /// At corner `k`.
    Corner(u8),
}

/// The arcs of a face with the side counts `n` (only `n⁺` enters): for every crossing
/// `(side k, p)`, in order of `k` then `p`, the other end of the arc through it. This is the
/// pairing of Supp. Alg. 9 as [`IntegerTriangulation::common_subdivision`] uses it (read off the
/// pieces it cuts), and what the oracle's task `subdivide` states.
///
/// # Errors
/// [`IntegerError::InvalidCoordinates`] (naming face 0) if the face is not valid (an odd
/// numerator of Eq. 1, or counts that do not add up).
pub fn face_arcs(n: [i64; 3]) -> Result<Vec<(u8, u32, ArcEnd)>, IntegerError> {
    let counts = n.map(|x| x.max(0));
    let fc = FaceCoordinates::of(counts)
        .map_err(|k| bad_face(0, format!("the numerator of Eq. 1 at corner {k} is odd")))?;
    for k in 0..3 {
        if fc.n[k] != fc.c[k] + fc.c[(k + 1) % 3] + fc.e[(k + 2) % 3] || fc.c[k] < 0 {
            return Err(bad_face(
                0,
                format!("the counts {counts:?} do not add up at side {k}"),
            ));
        }
    }
    let mut arcs: HashMap<(u8, u32), ArcEnd> = HashMap::new();
    cut_face(&fc, |piece| {
        let nodes = piece.nodes();
        for (a, b) in nodes.iter().zip(nodes.iter().cycle().skip(1)) {
            let on = |nd: &Node, s: u8| match *nd {
                Node::Corner(k) => k == s || (k + 2) % 3 == s,
                Node::Cross { side, .. } => side == s,
            };
            if (0..3u8).any(|s| on(a, s) && on(b, s)) {
                continue;
            }
            match (*a, *b) {
                (Node::Cross { side: s, i: p }, Node::Cross { side: t, i: q }) => {
                    arcs.insert((s, p), ArcEnd::Side { k: t, p: q });
                    arcs.insert((t, q), ArcEnd::Side { k: s, p });
                }
                (Node::Corner(k), Node::Cross { side: s, i: p })
                | (Node::Cross { side: s, i: p }, Node::Corner(k)) => {
                    arcs.insert((s, p), ArcEnd::Corner(k));
                }
                _ => {}
            }
        }
        Ok(())
    })?;
    let mut out: Vec<_> = arcs.into_iter().map(|((s, p), e)| (s, p, e)).collect();
    out.sort_by_key(|&(s, p, _)| (s, p));
    Ok(out)
}
