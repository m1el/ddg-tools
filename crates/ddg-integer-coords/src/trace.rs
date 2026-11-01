//! Reading the integers: the combinatorial part of Supp. Alg. 1–4 (SPEC §3.1).
//!
//! * [`trace_from`] is Alg. 1 (`TraceFrom`): from a combinatorial crossing, follow the curve
//!   face by face until it ends at a vertex (`ic-trace-step`).
//! * [`extract_curve`] is the combinatorial part of Alg. 2 (`ExtractCurve`): both directions,
//!   the crossing list from one end to the other.
//! * [`IntegerTriangulation::identify_traced_edge`] is Alg. 3 (`IdentifyTracedEdge`): the input
//!   edge a curve belongs to, walking back past inserted vertices that lie on it, with the
//!   halfedge `Next(Next(ij→)) = a→i` (SPEC §3.1; `ic-identify-edge`).
//! * [`IntegerTriangulation::extract_edge`] is the combinatorial part of Alg. 4 (`ExtractEdge`):
//!   the path of an input edge over the triangulation, continued through the inserted vertices
//!   on it.
//!
//! The geometry (unfolding the strip, the positions of the crossings, Alg. 2 lines 6–17) is work
//! package WP2; it consumes the types here ([`CurveSegment`], [`InputEdgePath`]) unchanged.
//!
//! **Crossings.** A combinatorial crossing `ζ = (h, p)` ([`Crossing`]) is the `p`-th crossing of
//! halfedge `h` counted from its tail, *with the curve leaving the face of `h` through `h`* (the
//! convention of Alg. 1 and of the oracle's `tracecurve` task). Its reversal is
//! `(twin(h), n − p − 1)` (`ic-reverse-crossing`). An *end* of a curve at a vertex is the first
//! crossing of the curve leaving the vertex, so [`trace_from`] of an end walks away from it.
//!
//! **Bad data never panics** (SPEC §3.1): a crossing index outside `[0, n⁺)`, a face whose
//! Eq. 1 is odd, a curve that leaves through a boundary edge, a crossing visited twice (SPEC
//! I18: a repeated `(halfedge, index)` pair in either direction, which no valid state has) are
//! [`IntegerError::InvalidCoordinates`] naming the face and the crossing.

use std::collections::HashSet;

use ddg_intrinsic::SurfacePoint;

use crate::coords::{IntegerTriangulation, face_coordinates_of};
use crate::error::{IntegerError, Place, ValidityItem};
use crate::formulas::{
    FaceCoordinates, Step, StepError, identify_index, n_plus, reverse_crossing, trace_step,
};

/// What tracing reads: a gluing of triangles (`h = 3f + k` from corner `k` to corner `k + 1`) with
/// one normal coordinate per halfedge, equal on twins. [`IntegerTriangulation`] is one; the tests
/// implement it for the hand-built Δ-complexes of the oracle's fixtures.
pub trait NormalGluing {
    /// Number of halfedges (three per face).
    fn n_halfedges(&self) -> usize;
    /// The twin of halfedge `h`, or `None` on the boundary.
    fn twin(&self, h: usize) -> Option<usize>;
    /// The vertex at the tail of halfedge `h` (its corner).
    fn tail(&self, h: usize) -> usize;
    /// The normal coordinate of the edge of halfedge `h`.
    fn normal(&self, h: usize) -> i64;
}

impl NormalGluing for IntegerTriangulation {
    #[inline]
    fn n_halfedges(&self) -> usize {
        self.complex.n_halfedges()
    }
    #[inline]
    fn twin(&self, h: usize) -> Option<usize> {
        self.complex.twin(h)
    }
    #[inline]
    fn tail(&self, h: usize) -> usize {
        self.complex.from_vertex(h)
    }
    #[inline]
    fn normal(&self, h: usize) -> i64 {
        self.halfedge_normal(h)
    }
}

/// A combinatorial crossing `(halfedge, index)`: the curve leaves the face of `halfedge` through
/// it at its `index`-th crossing counted from its tail (Sec. 2.5, Supp. Alg. 1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Crossing {
    /// The halfedge crossed.
    pub halfedge: u32,
    /// The crossing index, from the halfedge's tail, in `[0, n⁺)`.
    pub index: i64,
}

impl Crossing {
    /// The crossing `(halfedge, index)`.
    #[inline]
    pub fn new(halfedge: usize, index: i64) -> Crossing {
        Crossing {
            halfedge: halfedge as u32,
            index,
        }
    }

    /// The place, for error messages.
    #[inline]
    pub fn place(&self) -> Place {
        Place::Crossing {
            halfedge: self.halfedge,
            index: self.index,
            face: self.halfedge / 3,
        }
    }

    /// The same crossing seen from the other side, `(twin, n − p − 1)`; `None` on a boundary
    /// halfedge.
    ///
    /// [formula: ic-reverse-crossing]
    pub fn reversed<G: NormalGluing>(&self, g: &G) -> Option<Crossing> {
        let h = self.halfedge as usize;
        let t = g.twin(h)?;
        Some(Crossing::new(
            t,
            reverse_crossing(n_plus(g.normal(h)), self.index),
        ))
    }

    /// The corner (halfedge id) opposite the crossed halfedge in its face: the vertex an arc
    /// through this crossing comes from when the crossing is an end of a curve.
    #[inline]
    pub fn opposite_corner(&self) -> usize {
        let h = self.halfedge as usize;
        3 * (h / 3) + (h % 3 + 2) % 3
    }
}

/// Where a curve ended: the corner (halfedge id, `3f + k`) of the face in which it arrived at a
/// vertex, the vertex, and the arc by which it arrived.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Arrival {
    /// The corner: `face = corner / 3`, slot `corner % 3`.
    pub corner: u32,
    /// The vertex of that corner.
    pub vertex: u32,
    /// The arrival seen from the vertex: the crossing of the opposite side through which the arc
    /// leaves the face when followed away from the vertex.
    pub back: Crossing,
}

/// The result of [`trace_from`]: the crossings visited, starting with the given one, and where
/// the curve ended.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalfTrace {
    /// The crossings in order, the start first.
    pub crossings: Vec<Crossing>,
    /// The vertex where the curve ended.
    pub end: Arrival,
}

/// The decoded face of halfedge `h`'s face, or the error naming it.
fn face_of<G: NormalGluing>(g: &G, f: usize) -> Result<FaceCoordinates, IntegerError> {
    face_coordinates_of(
        [g.normal(3 * f), g.normal(3 * f + 1), g.normal(3 * f + 2)],
        f,
    )
}

fn invalid_at(z: Crossing, detail: String) -> IntegerError {
    IntegerError::invalid(ValidityItem::Extraction, z.place(), detail)
}

/// Supp. Alg. 1 without storing the crossings: calls `visit` on every crossing (the start
/// first; an error from `visit` stops the walk and is returned) and returns where the curve
/// ended and the number of crossings visited. `bound` limits that number (`None`: no limit; a
/// walk over finitely many crossings then ends at a vertex or repeats a crossing, which `visit`
/// must detect, as [`trace_from`] does). The order of the checks per crossing: its index is in
/// range, the bound, `visit`, the twin exists, the face entered is valid.
///
/// # Errors
/// As [`trace_from`], and whatever `visit` returns.
///
/// [formula: ic-trace-step]
pub fn trace_with<G: NormalGluing>(
    g: &G,
    start: Crossing,
    bound: Option<u64>,
    visit: &mut dyn FnMut(Crossing) -> Result<(), IntegerError>,
) -> Result<(Arrival, u64), IntegerError> {
    let nh = g.n_halfedges();
    let mut z = start;
    let mut steps = 0u64;
    loop {
        let h = z.halfedge as usize;
        if h >= nh {
            return Err(IntegerError::IndexOutOfRange {
                what: "halfedge",
                index: h,
                bound: nh,
            });
        }
        let n_ij = g.normal(h);
        if z.index < 0 || z.index >= n_plus(n_ij) {
            return Err(invalid_at(
                z,
                format!("crossing index {} is not in [0, {})", z.index, n_plus(n_ij)),
            ));
        }
        if bound.is_some_and(|b| steps >= b) {
            return Err(invalid_at(
                z,
                format!("no vertex reached after {steps} crossings"),
            ));
        }
        visit(z)?;
        steps += 1;
        let Some(t) = g.twin(h) else {
            return Err(invalid_at(
                z,
                "the curve leaves the surface through a boundary halfedge".to_string(),
            ));
        };
        let (f, s) = (t / 3, t % 3);
        let (h_ik, h_kj) = (3 * f + (s + 1) % 3, 3 * f + (s + 2) % 3);
        match trace_step(n_ij, g.normal(h_ik), g.normal(h_kj), z.index) {
            Ok(Step::ToIk(q)) => z = Crossing::new(h_ik, q),
            Ok(Step::ToKj(q)) => z = Crossing::new(h_kj, q),
            Ok(Step::Vertex) => {
                let corner = h_kj;
                let back = Crossing::new(t, reverse_crossing(n_plus(n_ij), z.index));
                let end = Arrival {
                    corner: corner as u32,
                    vertex: g.tail(corner) as u32,
                    back,
                };
                return Ok((end, steps));
            }
            Err(StepError::InvalidFace) => return Err(face_of(g, f).unwrap_err()),
            Err(StepError::OutOfRange) => unreachable!("the index was checked above"),
        }
    }
}

/// The key of a crossing independent of the direction in which it is crossed.
fn undirected<G: NormalGluing>(g: &G, z: Crossing) -> (u32, i64) {
    match z.reversed(g) {
        Some(r) if r.halfedge < z.halfedge => (r.halfedge, r.index),
        _ => (z.halfedge, z.index),
    }
}

/// Supp. Alg. 1, `TraceFrom(ζ)`: follows the curve through `start` into the face of the twin
/// and on, face by face (three cases per face, `ic-trace-step`), until it ends at a vertex.
/// `bound` limits the number of crossings visited (`None`: none; the walk then stops at the
/// first repeated crossing, SPEC I18).
///
/// # Errors
/// [`IntegerError::InvalidCoordinates`] naming the crossing or the face: an index out of range,
/// a face whose Eq. 1 is odd, a boundary halfedge, a crossing visited twice (the curve is
/// closed), the bound exceeded; [`IntegerError::IndexOutOfRange`] for a halfedge that does not
/// exist.
///
/// [formula: ic-trace-terminates]
pub fn trace_from<G: NormalGluing>(
    g: &G,
    start: Crossing,
    bound: Option<u64>,
) -> Result<HalfTrace, IntegerError> {
    let mut crossings = Vec::new();
    let mut seen = HashSet::new();
    let (end, _) = trace_with(g, start, bound, &mut |z| {
        crossings.push(z);
        if seen.insert(undirected(g, z)) {
            Ok(())
        } else {
            Err(invalid_at(
                z,
                "crossing visited twice: the curve is closed".to_string(),
            ))
        }
    })?;
    Ok(HalfTrace { crossings, end })
}

/// A curve segment between two vertices of the triangulation, crossing at least one edge: the
/// crossings in order (each oriented in the direction of travel) and the corners at both ends.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CurveSegment {
    /// The corner (halfedge id) at the start vertex, in the face of the first crossing.
    pub start: u32,
    /// The vertex at the start.
    pub start_vertex: u32,
    /// The crossings from the start to the end; never empty.
    pub crossings: Vec<Crossing>,
    /// Where the segment ends.
    pub end: Arrival,
}

impl CurveSegment {
    /// The end at the start vertex: the first crossing (`None` only for a segment built by
    /// hand without crossings).
    #[inline]
    pub fn first(&self) -> Option<Crossing> {
        self.crossings.first().copied()
    }
}

/// The combinatorial part of Supp. Alg. 2, `ExtractCurve(ζ)`: traces forwards from `start`
/// and backwards from its reversal and joins the two halves into the segment between two
/// vertices, in the direction of `start`.
///
/// # Errors
/// As [`trace_from`]; a crossing visited twice in the whole segment is an error.
pub fn extract_curve<G: NormalGluing>(
    g: &G,
    start: Crossing,
    bound: Option<u64>,
) -> Result<CurveSegment, IntegerError> {
    let mut seen = HashSet::new();
    let mut forward = Vec::new();
    let (end, _) = trace_with(g, start, bound, &mut |z| {
        forward.push(z);
        if seen.insert(undirected(g, z)) {
            Ok(())
        } else {
            Err(invalid_at(
                z,
                "crossing visited twice: the curve is closed".to_string(),
            ))
        }
    })?;
    let Some(back_start) = start.reversed(g) else {
        return Err(invalid_at(
            start,
            "the crossing lies on a boundary halfedge".to_string(),
        ));
    };
    let mut backward = Vec::new();
    let mut first = true;
    let (start_end, _) = trace_with(g, back_start, bound, &mut |z| {
        if std::mem::take(&mut first) {
            return Ok(());
        }
        backward.push(z);
        if seen.insert(undirected(g, z)) {
            Ok(())
        } else {
            Err(invalid_at(
                z,
                "crossing visited twice: the curve is closed".to_string(),
            ))
        }
    })?;
    let mut crossings: Vec<Crossing> = Vec::with_capacity(backward.len() + forward.len());
    for z in backward.iter().rev() {
        // Reversible: every crossing of `backward` was entered through its halfedge's twin.
        crossings.push(z.reversed(g).expect("a visited crossing has a twin"));
    }
    crossings.extend_from_slice(&forward);
    Ok(CurveSegment {
        start: start_end.corner,
        start_vertex: start_end.vertex,
        crossings,
        end,
    })
}

// ----- input edges -----

/// One end of an input-edge curve at a vertex of the triangulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CurveEnd {
    /// The curve runs along this edge of the triangulation (`n = −1`).
    Along {
        /// The edge.
        edge: u32,
    },
    /// The curve leaves the vertex inside the face of the crossing's halfedge and first crosses
    /// it there.
    Arc(Crossing),
}

/// A piece of an input edge between two consecutive vertices of the triangulation on it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PathPiece {
    /// Along an edge of the triangulation with `n = −1`.
    Along {
        /// The edge.
        edge: u32,
    },
    /// Across the faces of the triangulation.
    Crossings(CurveSegment),
}

/// The path of an input edge over the triangulation (combinatorial part of Supp. Alg. 4).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputEdgePath {
    /// The input edge.
    pub input_edge: u32,
    /// The vertices of the triangulation on it in order: its tail, the inserted vertices that
    /// lie on it, its head.
    pub vertices: Vec<u32>,
    /// The pieces between consecutive vertices (one fewer than vertices).
    pub pieces: Vec<PathPiece>,
}

/// The input edge a curve belongs to (Supp. Alg. 3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TracedInputEdge {
    /// The input edge.
    pub input_edge: u32,
    /// The input vertex the curve was walked back to.
    pub from: u32,
    /// The local index of the input edge at `from` (SPEC I3).
    pub local_index: u32,
}

impl IntegerTriangulation {
    /// The bound on the number of crossings of one walk (SPEC §3.1: `Σ n⁺ + |E¹|`).
    pub fn walk_bound(&self) -> u64 {
        let total = self.total_crossings() + self.complex.n_edges() as u128;
        u64::try_from(total).unwrap_or(u64::MAX)
    }

    pub(crate) fn check_crossing(&self, z: Crossing) -> Result<(), IntegerError> {
        let nh = self.complex.n_halfedges();
        if z.halfedge as usize >= nh {
            return Err(IntegerError::IndexOutOfRange {
                what: "halfedge",
                index: z.halfedge as usize,
                bound: nh,
            });
        }
        Ok(())
    }

    /// Supp. Alg. 1 on this triangulation ([`trace_from`]).
    ///
    /// # Errors
    /// As [`trace_from`].
    pub fn trace_from(&self, start: Crossing) -> Result<HalfTrace, IntegerError> {
        self.check_crossing(start)?;
        trace_from(self, start, None)
    }

    /// The combinatorial part of Supp. Alg. 2 on this triangulation ([`extract_curve`]).
    ///
    /// # Errors
    /// As [`trace_from`].
    pub fn extract_curve(&self, start: Crossing) -> Result<CurveSegment, IntegerError> {
        self.check_crossing(start)?;
        extract_curve(self, start, None)
    }

    /// The ends of the input edges at input vertex `a`, by local index (`deg₀(a)` of them,
    /// SPEC I3): each is an edge of the triangulation with `n = −1` along it or the arc that
    /// leaves `a` in the right corner. This is the window search of Supp. Alg. 4 lines 1–12 for
    /// all input edges at `a` at once, and it checks that the roundabouts and emanating counts
    /// around `a` agree: going counter-clockwise, the corner of each outgoing halfedge `h`
    /// holds `n⁻_h + e` input directions, and the next halfedge's roundabout is `r_h + n⁻_h + e`
    /// (modulo `deg₀(a)` at an interior vertex, where the counts add up to `deg₀(a)`; at a
    /// boundary vertex the first is the boundary halfedge with `r = 0` and the last corner ends
    /// at `deg₀(a) − 1`, the boundary edge that arrives).
    ///
    /// # Errors
    /// [`IntegerError::InvalidCoordinates`] (item [`ValidityItem::Extraction`]) naming the
    /// halfedge whose roundabout disagrees; item [`ValidityItem::FaceParity`] for a face around
    /// `a` whose Eq. 1 is odd; [`IntegerError::IndexOutOfRange`] if `a` is not an input vertex.
    pub fn input_directions(&self, a: usize) -> Result<Vec<CurveEnd>, IntegerError> {
        let n_input = self.n_input_vertices();
        if a >= n_input {
            return Err(IntegerError::IndexOutOfRange {
                what: "input vertex",
                index: a,
                bound: n_input,
            });
        }
        let deg = self.local.degree(a) as i64;
        let boundary = self.local.is_boundary(a);
        let bad = |h: usize, detail: String| {
            IntegerError::invalid(ValidityItem::Extraction, Place::Halfedge(h as u32), detail)
        };
        if boundary != self.complex.is_boundary_vertex(a) {
            return Err(IntegerError::invalid(
                ValidityItem::Ranges,
                Place::Vertex(a as u32),
                "the vertex is on the boundary of one triangulation and not of the other"
                    .to_string(),
            ));
        }
        let outs: Vec<usize> = self.complex.outgoing_halfedges(a).collect();
        // Per corner: (r, n⁻, e, the first emanating crossing).
        let mut wedges = Vec::with_capacity(outs.len());
        let mut total = 0i64;
        for &h in &outs {
            let Some(r) = self.roundabout(h) else {
                return Err(bad(
                    h,
                    "a halfedge leaving an input vertex has no roundabout".into(),
                ));
            };
            if r as i64 >= deg {
                return Err(bad(h, format!("roundabout {r} is not below deg₀ = {deg}")));
            }
            let coords = self.face_coordinates(h / 3)?;
            let slot = h % 3;
            let m = self.n_minus(h);
            let e = coords.e[slot];
            total = total.saturating_add(m).saturating_add(e);
            wedges.push((r as i64, m, e, coords.first_emanating(slot)));
        }
        let mut ends = vec![None; deg as usize];
        if boundary {
            let (r0, m0, ..) = wedges[0];
            if r0 != 0 || m0 != 1 {
                return Err(bad(
                    outs[0],
                    format!(
                        "the boundary halfedge leaving the vertex has roundabout {r0} and n⁻ = \
                         {m0}, expected 0 and 1"
                    ),
                ));
            }
        } else if total != deg {
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                Place::Vertex(a as u32),
                format!(
                    "the corners around the vertex hold {total} input directions (n⁻ plus \
                     emanating arcs), but deg₀ = {deg}"
                ),
            ));
        }
        for (t, &h) in outs.iter().enumerate() {
            let (r, m, e, first) = wedges[t];
            let next = if t + 1 < outs.len() {
                Some((outs[t + 1], wedges[t + 1].0))
            } else if boundary {
                None
            } else {
                Some((outs[0], wedges[0].0))
            };
            let expected = r + m + e;
            match next {
                Some((h_next, r_next)) => {
                    let ok = if boundary {
                        r_next == expected
                    } else {
                        r_next == expected.rem_euclid(deg)
                    };
                    if !ok {
                        return Err(bad(
                            h_next,
                            format!(
                                "roundabout {r_next}, but the previous halfedge {h} has \
                                 roundabout {r}, n⁻ = {m} and its corner emanates {e} arcs, so \
                                 {} was expected",
                                if boundary {
                                    expected
                                } else {
                                    expected.rem_euclid(deg)
                                }
                            ),
                        ));
                    }
                }
                None => {
                    if expected != deg - 1 {
                        return Err(bad(
                            h,
                            format!(
                                "the last corner at the boundary vertex ends at local index \
                                 {expected} (roundabout {r}, n⁻ = {m}, {e} emanating arcs), \
                                 expected deg₀ − 1 = {}",
                                deg - 1
                            ),
                        ));
                    }
                }
            }
            if boundary && expected > deg - 1 {
                return Err(bad(h, format!("local index {expected} beyond deg₀ − 1")));
            }
            let edge = self.complex.edge(h) as u32;
            let side = self.complex.next(h);
            let mut k = r;
            if m == 1 {
                ends[k.rem_euclid(deg) as usize] = Some(CurveEnd::Along { edge });
                k += 1;
            }
            for q in 0..e {
                let slot = &mut ends[(k + q).rem_euclid(deg) as usize];
                *slot = Some(CurveEnd::Arc(Crossing::new(side, first + q)));
            }
        }
        if boundary {
            // The boundary edge arriving at a: the edge of the halfedge before the last
            // outgoing one, which has no twin.
            let last = *outs.last().expect("a vertex has a corner");
            let arriving = self.complex.prev(last);
            if self.complex.twin(arriving).is_some() {
                return Err(bad(
                    arriving,
                    "expected the boundary halfedge arriving".into(),
                ));
            }
            ends[(deg - 1) as usize] = Some(CurveEnd::Along {
                edge: self.complex.edge(arriving) as u32,
            });
        }
        ends.into_iter()
            .enumerate()
            .map(|(k, end)| {
                end.ok_or_else(|| {
                    IntegerError::invalid(
                        ValidityItem::Extraction,
                        Place::Vertex(a as u32),
                        format!("no corner holds input direction {k}"),
                    )
                })
            })
            .collect()
    }

    /// The ends of input-edge curves at vertex `v`, counter-clockwise from its outgoing
    /// halfedge: edges with `n = −1` and the arcs that emanate from `v`. An inserted vertex has
    /// two (on an input edge, made by an edge split) or none (inside an input face).
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; item [`ValidityItem::FaceParity`] for a face around
    /// `v` whose Eq. 1 is odd; item [`ValidityItem::Extraction`] for more ends than a valid
    /// state has there (`deg₀` at an input vertex, 2 at an inserted one).
    pub fn curve_ends(&self, v: usize) -> Result<Vec<CurveEnd>, IntegerError> {
        let nv = self.complex.n_vertices();
        if v >= nv {
            return Err(IntegerError::IndexOutOfRange {
                what: "vertex",
                index: v,
                bound: nv,
            });
        }
        // A valid state has deg₀ ends at an input vertex and 0 or 2 at an inserted one; more
        // than that is reported before anything large is built.
        let limit = if self.is_input_vertex(v) {
            self.local.degree(v) as i64
        } else {
            2
        };
        let too_many = || {
            IntegerError::invalid(
                ValidityItem::Extraction,
                Place::Vertex(v as u32),
                format!("more than {limit} curve ends at the vertex"),
            )
        };
        let mut ends = Vec::new();
        let mut last = None;
        for h in self.complex.outgoing_halfedges(v) {
            if self.halfedge_normal(h) == -1 {
                ends.push(CurveEnd::Along {
                    edge: self.complex.edge(h) as u32,
                });
            }
            let coords = self.face_coordinates(h / 3)?;
            let slot = h % 3;
            let side = self.complex.next(h);
            if ends.len() as i64 + coords.e[slot] > limit {
                return Err(too_many());
            }
            for q in 0..coords.e[slot] {
                ends.push(CurveEnd::Arc(Crossing::new(
                    side,
                    coords.first_emanating(slot) + q,
                )));
            }
            last = Some(h);
        }
        if let Some(h) = last
            && self.complex.is_boundary_vertex(v)
        {
            let arriving = self.complex.prev(h);
            if self.halfedge_normal(arriving) == -1 {
                ends.push(CurveEnd::Along {
                    edge: self.complex.edge(arriving) as u32,
                });
            }
        }
        if ends.len() as i64 > limit {
            return Err(too_many());
        }
        Ok(ends)
    }

    /// Follows one end of a curve to the next vertex: the vertex, and the end there by which it
    /// arrived. `on_crossing` sees every crossing; `piece` receives the piece if asked for.
    pub(crate) fn follow(
        &self,
        from: usize,
        end: CurveEnd,
        bound: Option<u64>,
        on_crossing: &mut dyn FnMut(Crossing) -> Result<(), IntegerError>,
        piece: Option<&mut Option<PathPiece>>,
    ) -> Result<(usize, CurveEnd), IntegerError> {
        match end {
            CurveEnd::Along { edge } => {
                let [x, y] = self.complex.edge_vertices(edge as usize);
                let to = if x == from {
                    y
                } else if y == from {
                    x
                } else {
                    return Err(IntegerError::invalid(
                        ValidityItem::Extraction,
                        Place::Edge(edge),
                        format!("the edge does not end at vertex {from}"),
                    ));
                };
                if let Some(p) = piece {
                    *p = Some(PathPiece::Along { edge });
                }
                Ok((to, CurveEnd::Along { edge }))
            }
            CurveEnd::Arc(z) => {
                let mut collected = Vec::new();
                let collect = piece.is_some();
                let (arrival, _) = trace_with(self, z, bound, &mut |c| {
                    if collect {
                        collected.push(c);
                    }
                    on_crossing(c)
                })?;
                if let Some(p) = piece {
                    *p = Some(PathPiece::Crossings(CurveSegment {
                        start: z.opposite_corner() as u32,
                        start_vertex: from as u32,
                        crossings: collected,
                        end: arrival,
                    }));
                }
                Ok((arrival.vertex as usize, CurveEnd::Arc(arrival.back)))
            }
        }
    }

    /// At an inserted vertex `m` on the curve of input edge `input_edge`, entered by `arrived`:
    /// checks that `m` lies on that input edge and has exactly two ends, and returns the other.
    pub(crate) fn pass_inserted(
        &self,
        m: usize,
        input_edge: usize,
        arrived: CurveEnd,
    ) -> Result<CurveEnd, IntegerError> {
        let bad = |detail: String| {
            IntegerError::invalid(ValidityItem::Extraction, Place::Vertex(m as u32), detail)
        };
        match self.locations[m] {
            SurfacePoint::Edge { edge, .. } if edge as usize == input_edge => {}
            other => {
                return Err(bad(format!(
                    "the curve of input edge {input_edge} passes through this inserted vertex, \
                     which lies at {other}"
                )));
            }
        }
        let ends = self.curve_ends(m)?;
        if ends.len() != 2 {
            return Err(bad(format!(
                "an inserted vertex on an input edge must have 2 curve ends, it has {}",
                ends.len()
            )));
        }
        if ends[0] == arrived {
            Ok(ends[1])
        } else if ends[1] == arrived {
            Ok(ends[0])
        } else {
            Err(bad(format!(
                "the curve arrives by {arrived:?}, which is not one of its ends {ends:?}"
            )))
        }
    }

    /// The combinatorial part of Supp. Alg. 4, `ExtractEdge`: the path of input halfedge
    /// `input_halfedge` (an id of the input mesh) over the triangulation, from its tail to its
    /// head, through the inserted vertices that lie on it (SPEC §3.1, I17). Checks that the
    /// path arrives at the head by the end that the roundabouts there give the reverse
    /// direction (§2.4 item 4 for this edge).
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; [`IntegerError::InvalidCoordinates`] naming the
    /// vertex, halfedge, face or crossing where the integers are inconsistent.
    pub fn extract_edge(&self, input_halfedge: usize) -> Result<InputEdgePath, IntegerError> {
        let mesh = self.input.mesh();
        if input_halfedge >= mesh.n_halfedges() {
            return Err(IntegerError::IndexOutOfRange {
                what: "input halfedge",
                index: input_halfedge,
                bound: mesh.n_halfedges(),
            });
        }
        let (a, b) = (
            mesh.from_vertex(input_halfedge),
            mesh.to_vertex(input_halfedge),
        );
        let start = self.input_directions(a)?[self.local.of_halfedge(input_halfedge) as usize];
        let k_b = self
            .local
            .of_edge_at(mesh, mesh.edge(input_halfedge), b)
            .expect("the head of an input halfedge is an end point of its edge");
        let expected = self.input_directions(b)?[k_b as usize];
        self.extract_edge_between(input_halfedge, start, expected)
    }

    /// [`extract_edge`](Self::extract_edge) for the `edge_halfedge` of every input edge, in the
    /// order of the input edges, with `input_directions` computed once per input vertex instead
    /// of twice per edge.
    ///
    /// # Errors
    /// The first error of [`extract_edge`](Self::extract_edge) (edges in order), or of
    /// [`input_directions`](Self::input_directions).
    pub fn extract_all_edges(&self) -> Result<Vec<InputEdgePath>, IntegerError> {
        let mesh = self.input.mesh();
        let directions = (0..self.n_input_vertices())
            .map(|a| self.input_directions(a))
            .collect::<Result<Vec<_>, _>>()?;
        (0..mesh.n_edges())
            .map(|e| {
                let h = mesh.edge_halfedge(e);
                let (a, b) = (mesh.from_vertex(h), mesh.to_vertex(h));
                let k_b = self
                    .local
                    .of_edge_at(mesh, e, b)
                    .expect("the head of an input halfedge is an end point of its edge");
                self.extract_edge_between(
                    h,
                    directions[a][self.local.of_halfedge(h) as usize],
                    directions[b][k_b as usize],
                )
            })
            .collect()
    }

    /// The walk of [`extract_edge`](Self::extract_edge), given the end at the tail where the
    /// curve starts and the end at the head where it must arrive.
    fn extract_edge_between(
        &self,
        input_halfedge: usize,
        start: CurveEnd,
        expected: CurveEnd,
    ) -> Result<InputEdgePath, IntegerError> {
        let mesh = self.input.mesh();
        let input_edge = mesh.edge(input_halfedge);
        let (a, b) = (
            mesh.from_vertex(input_halfedge),
            mesh.to_vertex(input_halfedge),
        );
        let mut seen = HashSet::new();
        let mut on_crossing = |z: Crossing| {
            if seen.insert(undirected(self, z)) {
                Ok(())
            } else {
                Err(invalid_at(
                    z,
                    format!("crossing visited twice on input edge {input_edge}"),
                ))
            }
        };
        let mut vertices = vec![a as u32];
        let mut pieces = Vec::new();
        let (mut v, mut end) = (a, start);
        for _ in 0..=self.complex.n_vertices() {
            let mut piece = None;
            let (x, arrived) = self.follow(v, end, None, &mut on_crossing, Some(&mut piece))?;
            vertices.push(x as u32);
            pieces.push(piece.expect("follow fills the piece"));
            if self.is_input_vertex(x) {
                if x != b || arrived != expected {
                    return Err(IntegerError::invalid(
                        ValidityItem::Extraction,
                        Place::InputEdge(input_edge as u32),
                        format!(
                            "the curve from vertex {a} ends at vertex {x} by {arrived:?}, but the \
                             roundabouts at vertex {b} give {expected:?}"
                        ),
                    ));
                }
                return Ok(InputEdgePath {
                    input_edge: input_edge as u32,
                    vertices,
                    pieces,
                });
            }
            end = self.pass_inserted(x, input_edge, arrived)?;
            v = x;
        }
        Err(IntegerError::invalid(
            ValidityItem::Extraction,
            Place::InputEdge(input_edge as u32),
            "the curve passes more inserted vertices than there are".to_string(),
        ))
    }

    /// Supp. Alg. 3, `IdentifyTracedEdge`: the input edge that `segment` belongs to. If the
    /// segment starts at an inserted vertex, the curve is followed backwards through the other
    /// end there (an inserted vertex on an input edge has two) until an input vertex `a` is
    /// reached (Alg. 3 lines 1–5); then the local index at `a` is
    /// `r_ai + (p − c^{ja}_i) + n⁻_ai` for the first crossing `(ij→, p)` with `ai→ =
    /// Next(Next(ij→))` (line 9; SPEC §3.1), or the roundabout of the edge along which the curve
    /// leaves `a`.
    ///
    /// # Errors
    /// [`IntegerError::InvalidCoordinates`] naming the vertex or crossing where the walk fails.
    ///
    /// [formula: ic-identify-edge]
    pub fn identify_traced_edge(
        &self,
        segment: &CurveSegment,
    ) -> Result<TracedInputEdge, IntegerError> {
        let Some(&first) = segment.crossings.first() else {
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                Place::Vertex(segment.start_vertex),
                "a curve segment without crossings".to_string(),
            ));
        };
        self.check_crossing(first)?;
        let mut v = segment.start_vertex as usize;
        let nv = self.complex.n_vertices();
        if v >= nv || self.complex.from_vertex(first.opposite_corner()) != v {
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                first.place(),
                format!("the first crossing of the segment is not opposite its start vertex {v}"),
            ));
        }
        // The end by which the curve leaves v, looking forward.
        let mut leaving = CurveEnd::Arc(first);
        let mut seen = HashSet::new();
        let mut on_crossing = |z: Crossing| {
            if seen.insert(undirected(self, z)) {
                Ok(())
            } else {
                Err(invalid_at(z, "crossing visited twice".to_string()))
            }
        };
        for _ in 0..=self.complex.n_vertices() {
            if self.is_input_vertex(v) {
                return self.identify_at(v, leaving);
            }
            let ends = self.curve_ends(v)?;
            let back = match (ends.len(), ends.iter().position(|&x| x == leaving)) {
                (2, Some(i)) => ends[1 - i],
                _ => {
                    return Err(IntegerError::invalid(
                        ValidityItem::Extraction,
                        Place::Vertex(v as u32),
                        format!(
                            "walking back past an inserted vertex: it has {} curve ends {ends:?}, \
                             expected two, one of them {leaving:?}",
                            ends.len()
                        ),
                    ));
                }
            };
            let (x, arrived) = self.follow(v, back, None, &mut on_crossing, None)?;
            // Seen from x, the curve leaves x by the end it arrived by.
            v = x;
            leaving = arrived;
        }
        Err(IntegerError::invalid(
            ValidityItem::Extraction,
            Place::Vertex(segment.start_vertex),
            "walking back passes more inserted vertices than there are".to_string(),
        ))
    }

    /// The input direction at input vertex `a` of a curve leaving it by `leaving`.
    fn identify_at(&self, a: usize, leaving: CurveEnd) -> Result<TracedInputEdge, IntegerError> {
        let deg = self.local.degree(a);
        let modulus = self.local.modulus(a);
        let k = match leaving {
            CurveEnd::Along { edge } => {
                let h = self.complex.edge_halfedge(edge as usize);
                let h = if self.complex.from_vertex(h) == a {
                    Some(h)
                } else {
                    self.complex.twin(h)
                };
                match h {
                    Some(h) => self.roundabout(h).map(i64::from),
                    // A boundary edge arriving at a: the last direction.
                    None => Some(deg as i64 - 1),
                }
            }
            CurveEnd::Arc(z) => {
                let t = z.halfedge as usize;
                let ai = self.complex.next(self.complex.next(t));
                self.roundabout(ai).and_then(|r| {
                    identify_index(
                        r as i64,
                        self.halfedge_normal(ai),
                        self.halfedge_normal(t),
                        self.halfedge_normal(self.complex.next(t)),
                        z.index,
                        modulus,
                    )
                })
            }
        };
        match k {
            Some(k) if (0..deg as i64).contains(&k) => Ok(TracedInputEdge {
                input_edge: self.local.edge(a, k as u32) as u32,
                from: a as u32,
                local_index: k as u32,
            }),
            other => Err(IntegerError::invalid(
                ValidityItem::Extraction,
                Place::Vertex(a as u32),
                format!(
                    "the curve leaving by {leaving:?} has local index {other:?}, not in [0, {deg})"
                ),
            )),
        }
    }
}
