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
//! **Across input edges** (`n ≥ 0`, I20): the explicit split, with the crossings of the edge
//! divided at the split point by their positions along it, needs the geometric crossings
//! (work package WP2) and is not available yet ([`IntegerError::Unsupported`]).
//! [`IntegerTriangulation::split_edge_given`] performs it when the caller supplies those
//! geometric decisions (the number of crossings before the point and its place on the input
//! mesh), with [`formulas::split_side`] on both sides.

use ddg_intrinsic::{
    EdgeSplitPlan, InsertRefusal, Insertion, IntrinsicError, SurfacePoint, plan_edge_split,
};

use crate::coords::{IntegerTriangulation, NO_ROUNDABOUT};
use crate::error::{IntegerError, Place, ValidityItem};
use crate::formulas::{self, n_minus, split_edge_normal};

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
    /// Only edges along an input edge (`n = −1`, every boundary edge among them) for now; see the
    /// [module documentation](crate::split).
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] with the refusals of [`plan_edge_split`] (index, parameter,
    /// same face on both sides, too thin); [`IntegerError::Unsupported`] for an edge with
    /// `n ≥ 0`; [`IntegerError::InvalidCoordinates`] when the integers of the edge do not say
    /// which input edge it runs along. Nothing is changed in every case.
    ///
    /// [formula: ic-split-edge-normal]
    pub fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntegerError> {
        let plan = plan_edge_split(&self.complex, &self.lengths, e, t)?;
        let n = self.normal[e];
        if n >= 0 {
            return Err(IntegerError::Unsupported {
                operation: "split_edge",
                reason: format!(
                    "edge {e} crosses {n} input edges; splitting it needs the positions of its \
                     crossings (SPEC I20, work package WP2), which are not available yet"
                ),
            });
        }
        let location = self.location_along_input_edge(e, t)?;
        let integers = self.split_integers(e, 0)?;
        Ok(self.apply_split(e, plan, integers, location))
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
    /// [`split_edge`](Self::split_edge) will call this once the geometric crossings (WP2)
    /// decide `before` and `location`; until then it lets tests and callers that know them
    /// (e.g. from a signpost trace) perform the split.
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
        let integers = self.split_integers(e, before)?;
        Ok(self.apply_split(e, plan, integers, location))
    }

    /// The integers of the split of edge `e` with `before` of its crossings on the half from its
    /// tail, from the current state.
    fn split_integers(&self, e: usize, before: i64) -> Result<SplitIntegers, IntegerError> {
        let c = &self.complex;
        let n = self.normal[e];
        let invalid = |detail: String| {
            IntegerError::invalid(ValidityItem::Ranges, Place::Edge(e as u32), detail)
        };
        if !(0..=n.max(0)).contains(&before) {
            return Err(invalid(format!(
                "{before} crossings before the split point of an edge with n = {n}"
            )));
        }
        let h0 = c.edge_halfedge(e);
        let h1 = c.twin(h0);
        // One side: halfedge g (from the side's tail), the number of crossings before the point
        // counted from g's tail.
        let side = |g: usize, before: i64| -> Result<(i64, u32), IntegerError> {
            let f = g / 3;
            let n_face = self.face_normal(f);
            let s = g % 3;
            let (to_opposite, emanating_before) = if n == -1 {
                // Eq. 6, and no arc of the face reaches the side.
                (
                    split_edge_normal(n_face[(s + 2) % 3], n_face[(s + 1) % 3]),
                    0,
                )
            } else {
                let x = formulas::split_side(n_face, s, before).ok_or_else(|| {
                    IntegerError::invalid(
                        ValidityItem::FaceParity,
                        Place::Face(f as u32),
                        format!("face {f} with side counts {n_face:?} is not valid"),
                    )
                })?;
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
            Some(h1) => side(h1, n.max(0) - before)?,
            None => (0, NO_ROUNDABOUT),
        };
        let halves = if n == -1 {
            [-1, -1]
        } else {
            [before, n - before]
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
