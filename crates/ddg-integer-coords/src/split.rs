//! The edge split (paper Sec. 3.4, Eq. 6; SPEC §4.5, I20).
//!
//! The refusals and the new lengths come from [`ddg_intrinsic::plan_edge_split`] (bitwise as
//! the signpost representation, IT-SPEC N10, FlipOut R4), the connectivity from
//! [`Complex::split_edge`](ddg_intrinsic::Complex::split_edge): the half from the old tail keeps
//! the edge's id, both halves keep its direction and its constrained flag, the faces keep their
//! ids (see [`EdgeSplit`](ddg_intrinsic::complex::EdgeSplit)).
//!
//! **Along an input edge** (`n = −1`, which every boundary edge is): both halves keep `n = −1`;
//! the edge from the new vertex `m` to the opposite corner `k` of each face gets Eq. 6,
//! `n_mk = max(n_ki, n_jk, 0)`; `q⁰_m = (1 − t) q⁰_i + t q⁰_j` on the input edge (the same
//! arithmetic as `ddg_intrinsic::IntrinsicTriangulation::split_edge`, so both representations
//! place the vertex bitwise at the same point when the signposts know the edge lies along that
//! input edge). The roundabout of the new halfedge `k → m` (for `k ∈ V⁰`) is that of `k → i`
//! plus `n⁻_ki`: no input edge leaves `k` into the face towards the side `ij`, which no input edge
//! crosses. The halfedges `i → m` and `j → m` leave in the directions of `i → j` and `j → i` and
//! keep their roundabouts.
//!
//! **Across input edges** (`n ≥ 0`, I20): the explicit split of the intrinsic edge (not the
//! paper's face split then flip, which fails for constrained edges). Floating point decides
//! two things, from the geometric crossings of the edge
//! ([`IntegerTriangulation::extract_intrinsic_edge`], WP2): how many crossings lie before the
//! split point (their `distance` along the edge below `t ℓ`), and the place `q⁰_m`, the point
//! of the input face between the two neighbouring crossings, interpolated linearly between
//! them (the edge is straight in that face). A point within rounding of a crossing (a
//! barycentric coordinate that is not strictly positive) is refused. The integers then follow
//! exactly ([`formulas::split_side`] on both sides): the halves get `before` and `n − before`,
//! the edges to the opposite corners count the arcs separated by the new chord, and the
//! roundabout of `k → m` is that of `k → i` plus `n⁻_ki` plus the arcs from `k` that reach
//! the edge before `m`. Whatever the rounding, the result is the exact integer data of the
//! vertex placed between those two crossings. [`IntegerTriangulation::split_edge_given`]
//! takes both decisions from the caller instead.

use ddg_intrinsic::formulas::triangle_layout;
use ddg_intrinsic::insert::LOCATION_SNAP;
use ddg_intrinsic::{
    EdgeSplitPlan, EdgeTraceOptions, InsertRefusal, Insertion, IntrinsicError, SurfacePoint,
    plan_edge_split,
};

use crate::coords::{IntegerTriangulation, NO_ROUNDABOUT};
use crate::error::{IntegerError, Place, ValidityItem};
use crate::formulas::{self, n_minus, split_edge_normal};

/// Where the split point is on the input mesh, from floating point (see
/// [`IntegerTriangulation::crossed_split_place`]).
enum CrossedPlace {
    /// Strictly inside an input face, between two crossings: `before` of them lie before it.
    Inside { before: i64, location: SurfacePoint },
    /// Within [`LOCATION_SNAP`] of crossing number `crossing` (0-based from the tail of the
    /// edge's halfedge), which the vertex is put on (SPEC I28).
    AtCrossing {
        crossing: i64,
        location: SurfacePoint,
        /// The barycentric coordinate of the point opposite the input edge, in its input face.
        coordinate: f64,
    },
}

/// The integers of the two faces of a split, computed before anything changes.
struct SplitIntegers {
    /// `n` of the two halves, from the old tail.
    halves: [i64; 2],
    /// `n` of the edges to the opposite corners `c` (the face of the edge's halfedge) and `d`
    /// (the other face; unused on the boundary).
    to_opposite: [i64; 2],
    /// The roundabouts of `c → m` and `d → m`.
    roundabouts: [u32; 2],
}

impl IntegerTriangulation {
    /// Splits edge `e` at the parameter `t` from the tail of its halfedge (`0 < t < 1`; paper
    /// Sec. 3.4): a vertex `m` is inserted there and joined to the one or two opposite
    /// vertices. Ids as [`Complex::split_edge`](ddg_intrinsic::Complex::split_edge);
    /// [`Insertion::halves`] gives both halves in order from the old tail. The lengths are those
    /// of [`plan_edge_split`] (`t ℓ`, `(1 − t) ℓ`, `it-point-distances`), bitwise as the signpost
    /// representation. A constrained edge gives two constrained halves.
    ///
    /// See the [module documentation](crate::split) for the two kinds of edges.
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] with the refusals of [`plan_edge_split`] (index, parameter,
    /// same face on both sides, too thin) and, for an edge that crosses input edges,
    /// [`InsertRefusal::NotLocated`] when the point is within rounding of a crossing (its
    /// place would not be strictly inside an input face); the errors of
    /// [`extract_intrinsic_edge`](Self::extract_intrinsic_edge);
    /// [`IntegerError::InvalidCoordinates`] when the integers of the edge do not say which input
    /// edge it runs along. Nothing is changed in every case.
    ///
    /// [formula: ic-split-edge-normal]
    pub fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntegerError> {
        let plan = plan_edge_split(&self.complex, &self.lengths, e, t)?;
        let n = self.normal[e];
        let (before, at_crossing, location, snapped) = if n >= 0 {
            match self.crossed_split_place(e, t, plan.halves[0])? {
                CrossedPlace::Inside { before, location } => (before, false, location, None),
                CrossedPlace::AtCrossing {
                    crossing,
                    location,
                    coordinate,
                } => (crossing, true, location, Some(coordinate)),
            }
        } else {
            (0, false, self.location_along_input_edge(e, t)?, None)
        };
        let integers = self.split_integers(e, before, at_crossing)?;
        let mut insertion = self.apply_split(e, plan, integers, location);
        insertion.snapped = snapped;
        Ok(insertion)
    }

    /// The edge split of [`split_edge`](Self::split_edge) for an edge with `n ≥ 0` (SPEC I20),
    /// with the two geometric decisions supplied by the caller: `before`, the number of the
    /// edge's crossings between the tail of its halfedge and the split point (`0 ≤ before ≤ n`;
    /// the point must not lie on a crossing), and `location`, the place `q⁰` of the new vertex
    /// on the input mesh (a point of an input face). The halves get `before` and `n − before`,
    /// the edges to the opposite corners and the roundabouts of their halfedges leaving input
    /// vertices come from [`formulas::split_side`] on each side. For an edge with `n = −1`,
    /// `before` must be `0` and the result equals [`split_edge`](Self::split_edge)'s except
    /// for `location`, which is taken as given.
    ///
    /// [`split_edge`](Self::split_edge) is this with both decided from the geometric crossings
    /// of the edge (`ExtractIntrinsicEdge`); this form lets tests and callers that know them
    /// otherwise (e.g. from a signpost trace) perform the split with their decisions.
    ///
    /// # Errors
    /// As [`split_edge`](Self::split_edge), and [`IntegerError::InvalidCoordinates`] for a
    /// `before` out of range or a location that is not a point of the input mesh of the right
    /// kind. Nothing is changed in every case.
    pub fn split_edge_given(
        &mut self,
        e: usize,
        t: f64,
        before: i64,
        location: SurfacePoint,
    ) -> Result<Insertion, IntegerError> {
        let plan = plan_edge_split(&self.complex, &self.lengths, e, t)?;
        let input = self.input.mesh();
        let ok = match location {
            SurfacePoint::Vertex(_) => false,
            SurfacePoint::Edge { edge, t } => {
                (edge as usize) < input.n_edges() && t > 0.0 && t < 1.0
            }
            SurfacePoint::Face { face, bary } => {
                (face as usize) < input.n_faces() && bary.iter().all(|&x| x > 0.0 && x < 1.0)
            }
        };
        if !ok {
            return Err(IntegerError::invalid(
                ValidityItem::Ranges,
                Place::Edge(e as u32),
                format!("{location:?} is not a place for a vertex inserted on an edge"),
            ));
        }
        let integers = self.split_integers(e, before, false)?;
        Ok(self.apply_split(e, plan, integers, location))
    }

    /// The two floating-point decisions of the split of edge `e` (`n ≥ 0`) at distance `at` from
    /// the tail of its halfedge: the number of crossings before the point, and the point's place
    /// on the input mesh (see the [module documentation](crate::split)).
    fn crossed_split_place(
        &self,
        e: usize,
        t: f64,
        at: f64,
    ) -> Result<CrossedPlace, IntegerError> {
        let trace = self.extract_intrinsic_edge(e, &EdgeTraceOptions::default())?;
        let mesh = self.input.mesh();
        let refused = |detail: String| {
            IntegerError::Intrinsic(IntrinsicError::InsertionRefused {
                at: SurfacePoint::Edge { edge: e as u32, t },
                reason: InsertRefusal::NotLocated { detail },
            })
        };
        let invalid = |detail: String| {
            IntegerError::invalid(ValidityItem::Extraction, Place::Edge(e as u32), detail)
        };
        let x = &trace.crossings;
        let before = x.iter().filter(|c| c.distance < at).count();
        // The input face between crossing `before − 1` and crossing `before`.
        let face = match (x.get(before), before.checked_sub(1).map(|k| &x[k])) {
            (Some(next), _) => next.halfedge as usize / 3,
            (None, Some(last)) => {
                mesh.twin(last.halfedge as usize)
                    .ok_or_else(|| invalid("its last crossing is a boundary input edge".into()))?
                    / 3
            }
            (None, None) => trace
                .first_face
                .ok_or_else(|| invalid("it has no crossing and no first input face".into()))?
                as usize,
        };
        let corners = mesh.faces()[face];
        let corner_of = |v: usize| corners.iter().position(|&c| c as usize == v);
        // Barycentric coordinates in `face` of a point on the input mesh.
        let bary = |q: SurfacePoint| -> Option<[f64; 3]> {
            let mut b = [0.0; 3];
            match q {
                SurfacePoint::Vertex(v) => b[corner_of(v as usize)?] = 1.0,
                SurfacePoint::Edge { edge, t } => {
                    let [p, r] = mesh.edge_vertices(edge as usize);
                    b[corner_of(p)?] = 1.0 - t;
                    b[corner_of(r)?] = t;
                }
                SurfacePoint::Face { face: g, bary } => {
                    if g as usize != face {
                        return None;
                    }
                    b = bary;
                }
            }
            Some(b)
        };
        let crossing_point = |k: usize| SurfacePoint::Edge {
            edge: x[k].edge,
            t: x[k].t,
        };
        let h = self.complex.edge_halfedge(e);
        let (a, b) = (self.complex.from_vertex(h), self.complex.to_vertex(h));
        let (p0, d0) = match before.checked_sub(1) {
            Some(k) => (crossing_point(k), x[k].distance),
            None => (self.locations[a], 0.0),
        };
        let (p1, d1) = match x.get(before) {
            Some(c) => (crossing_point(before), c.distance),
            None => (self.locations[b], self.lengths[e]),
        };
        let (Some(b0), Some(b1)) = (bary(p0), bary(p1)) else {
            return Err(invalid(format!(
                "the points {p0:?} and {p1:?} around the split point are not in input face {face}"
            )));
        };
        let s = ((at - d0) / (d1 - d0)).clamp(0.0, 1.0);
        let w = [0, 1, 2].map(|k| (1.0 - s) * b0[k] + s * b1[k]);
        // The rule of the signposts (`snap_onto_input_edge`, IT-SPEC N10): a point within
        // `LOCATION_SNAP` times the longest side of its input face from a side is on that side.
        let l = [0, 1, 2].map(|k| self.input.geometry().edge_lengths()[mesh.edge(3 * face + k)]);
        if let Some(v) = triangle_layout(l) {
            let double_area = v[1][0] * v[2][1];
            let distance = [0, 1, 2].map(|k| w[k] * double_area / l[(k + 1) % 3]);
            let k = (0..3)
                .min_by(|&a, &b| distance[a].total_cmp(&distance[b]))
                .unwrap_or(0);
            let longest = l[0].max(l[1]).max(l[2]);
            if distance[k] <= LOCATION_SNAP * longest {
                // On the side opposite corner `k`, which must be an input edge that `e` crosses
                // at one of the two crossings around the point.
                let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                let t_side = w[j] / (w[i] + w[j]);
                let input_edge = mesh.edge(3 * face + i);
                let crossing = [Some(before), before.checked_sub(1)]
                    .into_iter()
                    .flatten()
                    .find(|&c| x.get(c).is_some_and(|c| c.edge as usize == input_edge));
                return match crossing {
                    Some(c) if t_side > 0.0 && t_side < 1.0 => Ok(CrossedPlace::AtCrossing {
                        crossing: c as i64,
                        location: crossing_point(c),
                        coordinate: w[k],
                    }),
                    _ => Err(refused(format!(
                        "the point is within rounding of input edge {input_edge} (barycentric \
                         coordinates {w:?} in input face {face}), which edge {e} does not cross \
                         there, or of a corner of the face"
                    ))),
                };
            }
        } else if !w.iter().all(|&c| c > 0.0 && c < 1.0) {
            return Err(refused(format!(
                "the point is within rounding of an input edge (barycentric coordinates {w:?} in \
                 input face {face})"
            )));
        }
        Ok(CrossedPlace::Inside {
            before: before as i64,
            location: SurfacePoint::Face {
                face: face as u32,
                bary: w,
            },
        })
    }

    /// The integers of the split of edge `e` with `before` of its crossings on the half from its
    /// tail, from the current state.
    ///
    /// With `at_crossing` the new vertex lies on crossing number `before` (SPEC I28,
    /// [`formulas::split_side_at_crossing`]): the halves get `before` and `n − 1 − before`.
    fn split_integers(
        &self,
        e: usize,
        before: i64,
        at_crossing: bool,
    ) -> Result<SplitIntegers, IntegerError> {
        let c = &self.complex;
        let n = self.normal[e];
        let invalid = |detail: String| {
            IntegerError::invalid(ValidityItem::Ranges, Place::Edge(e as u32), detail)
        };
        if !(0..=n.max(0) - i64::from(at_crossing)).contains(&before) {
            return Err(invalid(format!(
                "{before} crossings before the split point of an edge with n = {n}"
            )));
        }
        let h0 = c.edge_halfedge(e);
        let h1 = c.twin(h0);
        // One side: halfedge g (from the side's tail), the number of crossings before the point
        // counted from g's tail.
        let side = |g: usize, before: i64| -> Result<(i64, u32), IntegerError> {
            let invalid_face = |f: usize, n_face: [i64; 3]| {
                IntegerError::invalid(
                    ValidityItem::FaceParity,
                    Place::Face(f as u32),
                    format!("face {f} with side counts {n_face:?} is not valid"),
                )
            };
            let f = g / 3;
            let n_face = self.face_normal(f);
            let s = g % 3;
            let (to_opposite, emanating_before) = if n == -1 {
                // Eq. 6, and no arc of the face reaches the side.
                (
                    split_edge_normal(n_face[(s + 2) % 3], n_face[(s + 1) % 3]),
                    0,
                )
            } else if at_crossing {
                let x = formulas::split_side_at_crossing(n_face, s, before)
                    .ok_or_else(|| invalid_face(f, n_face))?;
                (x.to_opposite, x.emanating_before)
            } else {
                let x = formulas::split_side(n_face, s, before)
                    .ok_or_else(|| invalid_face(f, n_face))?;
                (x.to_opposite, x.emanating_before)
            };
            // The halfedge from the opposite corner to the side's tail, just clockwise of the
            // new halfedge from that corner to the split point.
            let back = c.prev(g);
            let o = c.from_vertex(back);
            let r = if self.is_input_vertex(o) {
                let deg = self.local.degree(o);
                let r = self.roundabouts[back];
                let x = if r < deg {
                    formulas::roundabout_update(
                        r as i64,
                        emanating_before,
                        n_minus(self.normal[c.edge(back)]),
                        self.local.modulus(o),
                    )
                } else {
                    None
                };
                match x {
                    Some(x) if (0..deg as i64).contains(&x) => x as u32,
                    _ => {
                        return Err(IntegerError::invalid(
                            ValidityItem::Ranges,
                            Place::Halfedge(back as u32),
                            format!(
                                "the roundabout {r} of halfedge {back} at input vertex {o} \
                                 (deg₀ {deg}) gives no roundabout for the new edge"
                            ),
                        ));
                    }
                }
            } else {
                NO_ROUNDABOUT
            };
            Ok((to_opposite, r))
        };
        let (n_c, r_c) = side(h0, before)?;
        let (n_d, r_d) = match h1 {
            Some(h1) => side(h1, n.max(0) - before - i64::from(at_crossing))?,
            None => (0, NO_ROUNDABOUT),
        };
        let halves = if n == -1 {
            [-1, -1]
        } else {
            [before, n - before - i64::from(at_crossing)]
        };
        Ok(SplitIntegers {
            halves,
            to_opposite: [n_c, n_d],
            roundabouts: [r_c, r_d],
        })
    }

    /// The place on the input mesh of the point at parameter `t` of edge `e`, which runs along
    /// an input edge (`n = −1`): interpolated between the places of its ends on that input edge,
    /// as the signpost representation does.
    fn location_along_input_edge(&self, e: usize, t: f64) -> Result<SurfacePoint, IntegerError> {
        let c = &self.complex;
        let h0 = c.edge_halfedge(e);
        let (a, b) = (c.from_vertex(h0), c.to_vertex(h0));
        let invalid = |detail: String| {
            IntegerError::invalid(ValidityItem::AlongInputEdge, Place::Edge(e as u32), detail)
        };
        // The input edge: from the place of an inserted end, else from the roundabout at the
        // tail (both ends input vertices).
        let input_edge = match (self.locations[a], self.locations[b]) {
            (SurfacePoint::Edge { edge, .. }, _) | (_, SurfacePoint::Edge { edge, .. }) => {
                edge as usize
            }
            _ if self.is_input_vertex(a) => {
                let r = self.roundabouts[h0];
                if r >= self.local.degree(a) {
                    return Err(invalid(format!(
                        "the roundabout {r} of its halfedge {h0} is out of range"
                    )));
                }
                self.local.edge(a, r)
            }
            (qa, qb) => {
                return Err(invalid(format!(
                    "n = −1 but its ends are at {qa:?} and {qb:?}, not on an input edge"
                )));
            }
        };
        let (ta, tb) = match (
            self.split_parameter(a, input_edge),
            self.split_parameter(b, input_edge),
        ) {
            (Some(ta), Some(tb)) if ta != tb => (ta, tb),
            (ta, tb) => {
                return Err(invalid(format!(
                    "n = −1 but its ends are not two points of input edge {input_edge} \
                     (parameters {ta:?}, {tb:?})"
                )));
            }
        };
        let s = (1.0 - t) * ta + t * tb;
        if !(s > 0.0 && s < 1.0) {
            return Err(IntegerError::Intrinsic(IntrinsicError::InsertionRefused {
                at: SurfacePoint::Edge { edge: e as u32, t },
                reason: InsertRefusal::NotLocated {
                    detail: format!(
                        "the point is at parameter {s} of input edge {input_edge}, not strictly \
                         inside it"
                    ),
                },
            }));
        }
        Ok(SurfacePoint::Edge {
            edge: input_edge as u32,
            t: s,
        })
    }

    /// The parameter of vertex `v` along input edge `input_edge` (`0` at its first end point,
    /// `1` at its second), if `v` lies on it (as the signpost representation reads it).
    fn split_parameter(&self, v: usize, input_edge: usize) -> Option<f64> {
        match self.locations[v] {
            SurfacePoint::Vertex(u) => {
                let [a, b] = self.input.mesh().edge_vertices(input_edge);
                if u as usize == a {
                    Some(0.0)
                } else if u as usize == b {
                    Some(1.0)
                } else {
                    None
                }
            }
            SurfacePoint::Edge { edge, t } if edge as usize == input_edge => Some(t),
            _ => None,
        }
    }

    /// Changes the state: connectivity, lengths, integers, the new vertex's place.
    fn apply_split(
        &mut self,
        e: usize,
        plan: EdgeSplitPlan,
        integers: SplitIntegers,
        location: SurfacePoint,
    ) -> Insertion {
        let split = self
            .complex
            .split_edge(e)
            .expect("plan_edge_split refuses what Complex::split_edge refuses");
        split.move_halfedge_data(&mut self.roundabouts, NO_ROUNDABOUT);
        let ne = self.complex.n_edges();
        self.lengths.resize(ne, 0.0);
        self.normal.resize(ne, 0);
        let [_, e2] = split.halves.map(|x| x as usize);
        let [ec, ed] = split.to_opposite;
        self.lengths[e2] = plan.halves[1];
        self.lengths[ec as usize] = plan.to_opposite[0];
        self.lengths[e] = plan.halves[0];
        self.normal[e] = integers.halves[0];
        self.normal[e2] = integers.halves[1];
        self.normal[ec as usize] = integers.to_opposite[0];
        // The halfedges at the new vertex: those leaving it have no roundabout; those
        // arriving from the opposite corners get the computed ones.
        for h in split.outgoing() {
            self.roundabouts[h as usize] = NO_ROUNDABOUT;
        }
        let c_to_m = self.complex.edge_halfedge(ec as usize);
        self.roundabouts[c_to_m] = integers.roundabouts[0];
        if ed != ddg_core::topology::NONE {
            self.lengths[ed as usize] = plan.to_opposite[1];
            self.normal[ed as usize] = integers.to_opposite[1];
            let d_to_m = self.complex.edge_halfedge(ed as usize);
            self.roundabouts[d_to_m] = integers.roundabouts[1];
        }
        self.locations.push(location);
        Insertion {
            vertex: split.vertex,
            location,
            new_edges: split.new_edges(),
            new_faces: split
                .new_faces
                .into_iter()
                .filter(|&f| f != ddg_core::topology::NONE)
                .collect(),
            halves: Some(split.halves),
            snapped: None,
        }
    }
}
