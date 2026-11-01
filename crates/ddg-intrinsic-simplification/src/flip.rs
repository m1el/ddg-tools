//! Intrinsic edge flips with point tracking (SPEC §3.2, S16) and Delaunay flipping (§3.3).

use std::collections::VecDeque;

use ddg_core::geometry::triangle_area;
use ddg_intrinsic::formulas::{barycentric, signpost_update, unfold_diamond};
use ddg_intrinsic::{DelaunayStatus, FlipLimit, FlipRefusal, plan_flip_with};

use crate::error::IemError;
use crate::formulas::wrap_angle;
use crate::mesh::{CoarseMesh, Location};

/// An edge flip, in fine ids: the end points before and after.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FlippedEdge {
    /// The end points before the flip (tail and head of the edge's halfedge).
    pub old: [u32; 2],
    /// The end points after it.
    pub new: [u32; 2],
}

/// What a Delaunay pass did (SPEC §3.3).
#[derive(Debug, Clone, PartialEq)]
pub struct DelaunayFlipReport {
    /// How it ended (as `ddg_intrinsic`'s Delaunay flipping).
    pub status: DelaunayStatus,
    /// The flips, in order.
    pub flips: Vec<FlippedEdge>,
    /// Edges taken from the queue and tested.
    pub n_tested: usize,
    /// Times an edge that failed the test was refused by the flip (skipped and counted).
    pub n_refused: usize,
    /// Edges that fail the test at the end among those that were refused.
    pub n_not_delaunay: usize,
    /// Points clamped onto a side during the flips (S16).
    pub n_clamped: usize,
}

/// Why a flip did not happen.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum FlipFailure {
    Refused(FlipRefusal),
    Internal(String),
}

/// A new location of a point computed by a flip: point, face, coordinates, clamped.
type Relocation = (u32, u32, [f64; 3], bool);

impl CoarseMesh {
    /// Flips edge `e` (internal id) with point tracking (SPEC §3.2): refused unless both
    /// faces form a strictly convex quadrilateral that passes the options' shape test
    /// (`ddg_intrinsic::plan_flip_with`). The new length is the other diagonal of the planar
    /// layout; the four outer halfedges keep their angular coordinates; the two new ones follow
    /// from their clockwise neighbours; a reference halfedge that is flipped passes to its
    /// clockwise neighbour; the points of the two faces are re-expressed through the layout.
    ///
    /// # Errors
    /// [`IemError::IndexOutOfRange`], [`IemError::RemovalPending`], [`IemError::FlipRefused`],
    /// and [`IemError::Internal`] if a point would lie outside both new faces (the mesh is
    /// unchanged then).
    pub fn flip_edge(&mut self, e: usize) -> Result<FlippedEdge, IemError> {
        if self.pending.is_some() {
            return Err(IemError::RemovalPending);
        }
        if e >= self.n_edges() {
            return Err(IemError::IndexOutOfRange {
                what: "edge",
                index: e,
                bound: self.n_edges(),
            });
        }
        let [a, b] = self.complex.edge_vertices(e);
        match self.flip_internal(e, true) {
            Ok((flipped, _)) => Ok(flipped),
            Err(FlipFailure::Refused(reason)) => Err(IemError::FlipRefused {
                edge: e as u32,
                vertices: [self.fine_of[a], self.fine_of[b]],
                reason,
            }),
            Err(FlipFailure::Internal(what)) => Err(IemError::Internal { what }),
        }
    }

    /// The flip of [`flip_edge`](Self::flip_edge), recorded on the undo log when a removal is
    /// pending; points are tracked only if `track`. Returns the flip in fine ids and the
    /// number of points clamped.
    pub(crate) fn flip_internal(
        &mut self,
        e: usize,
        track: bool,
    ) -> Result<(FlippedEdge, usize), FlipFailure> {
        let plan = plan_flip_with(&self.complex, &self.length, e, self.options.flip_test)
            .map_err(FlipFailure::Refused)?;
        let c = &self.complex;
        let h0 = c.edge_halfedge(e);
        let h1 = c
            .twin(h0)
            .ok_or(FlipFailure::Refused(FlipRefusal::BoundaryEdge))?;
        let (f0, f1) = (h0 / 3, h1 / 3);
        let (h0n, h1n) = (c.next(h0), c.next(h1));
        let (h0p, h1p) = (c.prev(h0), c.prev(h1));
        let (i, j) = (c.from_vertex(h0), c.from_vertex(h1));
        let (k, l) = (c.from_vertex(h0p), c.from_vertex(h1p));

        let mut relocations: Vec<Relocation> = Vec::new();
        let mut n_clamped = 0;
        if track && !(self.face_points[f0].is_empty() && self.face_points[f1].is_empty()) {
            relocations = self.flip_relocations(e, &plan.lengths, plan.new_length)?;
            n_clamped = relocations.iter().filter(|r| r.3).count();
        }

        // Connectivity.
        let undo = self.complex.flip_undo(e);
        self.undo.push_flip(undo);
        let (phi_il, phi_jk) = (self.phi[h1n], self.phi[h0n]);
        self.complex
            .flip_edge(e)
            .map_err(|err| FlipFailure::Internal(err.to_string()))?;
        // Lengths and angular coordinates: i → l moved to h0, j → k to h1; l → k (h0n) follows
        // l → j (h1p), k → l (h1n) follows k → i (h0p).
        self.set_length(e, plan.new_length);
        self.set_phi(h0, phi_il);
        self.set_phi(h1, phi_jk);
        let [at_l, at_k] = plan.new_corner_angles();
        let bl = self.complex.is_boundary_vertex(l);
        let bk = self.complex.is_boundary_vertex(k);
        let phi_lk = wrap_angle(
            signpost_update(self.phi[h1p], at_l, self.angle_sum[l], bl),
            bl,
        );
        let phi_kl = wrap_angle(
            signpost_update(self.phi[h0p], at_k, self.angle_sum[k], bk),
            bk,
        );
        self.set_phi(h0n, phi_lk);
        self.set_phi(h1n, phi_kl);
        // References: a flipped halfedge passes to its clockwise neighbour, which now sits in
        // the same slot (i → j's slot h0 holds i → l, j → i's slot h1 holds j → k); a moved one
        // follows its move.
        for v in [i, j] {
            let r = self.reference[v] as usize;
            if r == h1n {
                self.set_reference(v, h0 as u32);
            } else if r == h0n {
                self.set_reference(v, h1 as u32);
            }
        }
        // Points.
        for (p, face, bary, clamped) in relocations {
            let Location::Face { face: old, .. } = self.points[p as usize] else {
                continue;
            };
            if old == face {
                self.set_bary(p, bary);
            } else {
                self.move_point(p, face, bary);
            }
            if clamped {
                self.n_clamped += 1;
            }
        }
        let f = |v: usize| self.fine_of[v];
        Ok((
            FlippedEdge {
                old: [f(i), f(j)],
                new: [f(l), f(k)],
            },
            n_clamped,
        ))
    }

    /// The new locations of the points of the two faces of interior edge `e` after its flip
    /// (S16). The diamond is laid out with `i` at the origin, `j = (ℓ_ij, 0)`, `k` above and `l`
    /// below (`it-unfold-diamond`); the other diagonal's end points are expressed in each new
    /// face (`it-barycentric-coordinates`), and a point's coordinates in the new faces follow
    /// by linear substitution. Inside test: signed distances to the sides `≥ −slack · L` (`L`
    /// the mean of the five lengths); a point inside both goes to the face with the lower id,
    /// tiny negatives are clamped (counted), a point outside both is an internal error.
    fn flip_relocations(
        &self,
        e: usize,
        lengths: &[f64; 5],
        new_length: f64,
    ) -> Result<Vec<Relocation>, FlipFailure> {
        let c = &self.complex;
        let h0 = c.edge_halfedge(e);
        let h1 = c.twin(h0).unwrap_or(h0);
        let (f0, a) = (h0 / 3, h0 % 3);
        let (f1, b) = (h1 / 3, h1 % 3);
        let [l_ij, l_jk, l_ki, l_il, l_lj] = *lengths;
        let internal = |what: String| FlipFailure::Internal(what);
        let d = unfold_diamond(l_ij, l_jk, l_ki, l_il, l_lj)
            .ok_or_else(|| internal(format!("edge {e}: the diamond does not unfold")))?;
        let pi = [0.0, 0.0];
        let pj = [l_ij, 0.0];
        let pk = [d.x_k, d.h_k];
        let pl = [d.x_l, -d.h_l];
        // j in the new face (i, l, k); i in the new face (j, k, l).
        let bj = barycentric(&[pi, pl, pk], pj);
        let bi = barycentric(&[pj, pk, pl], pi);
        // Heights of the new faces over their sides, for the signed distances.
        let area0 = triangle_area(l_il, new_length, l_ki)
            .ok_or_else(|| internal(format!("edge {e}: new face (i, l, k) is degenerate")))?;
        let area1 = triangle_area(l_jk, new_length, l_lj)
            .ok_or_else(|| internal(format!("edge {e}: new face (j, k, l) is degenerate")))?;
        // Local order (i, l, k): opposite sides lk, ki, il. (j, k, l): kl, lj, jk.
        let height0 = [new_length, l_ki, l_il].map(|s| 2.0 * area0 / s);
        let height1 = [new_length, l_lj, l_jk].map(|s| 2.0 * area1 / s);
        let scale = (l_ij + l_jk + l_ki + l_il + l_lj) / 5.0;
        let slack = self.options.point_slack * scale;
        let min_distance =
            |x: [f64; 3], h: [f64; 3]| (x[0] * h[0]).min(x[1] * h[1]).min(x[2] * h[2]);
        let mut out = Vec::new();
        let rot = |s: usize, k: usize| (s + k) % 3;
        for (face, slot, from_f0) in [(f0, a, true), (f1, b, false)] {
            for &p in &self.face_points[face] {
                let Location::Face { bary, .. } = self.points[p as usize] else {
                    continue;
                };
                let (c0, c1) = if from_f0 {
                    // Old (i, j, k).
                    let (b_i, b_j, b_k) =
                        (bary[rot(slot, 0)], bary[rot(slot, 1)], bary[rot(slot, 2)]);
                    (
                        [b_i + b_j * bj[0], b_j * bj[1], b_k + b_j * bj[2]],
                        [b_j + b_i * bi[0], b_k + b_i * bi[1], b_i * bi[2]],
                    )
                } else {
                    // Old (j, i, l).
                    let (b_j, b_i, b_l) =
                        (bary[rot(slot, 0)], bary[rot(slot, 1)], bary[rot(slot, 2)]);
                    (
                        [b_i + b_j * bj[0], b_l + b_j * bj[1], b_j * bj[2]],
                        [b_j + b_i * bi[0], b_i * bi[1], b_l + b_i * bi[2]],
                    )
                };
                let m0 = min_distance(c0, height0);
                let m1 = min_distance(c1, height1);
                let ok0 = m0 >= -slack;
                let ok1 = m1 >= -slack;
                let use0 = match (ok0, ok1) {
                    (true, true) => f0 < f1,
                    (true, false) => true,
                    (false, true) => false,
                    (false, false) => {
                        return Err(internal(format!(
                            "point {p} lies outside both faces of the flip of edge {e} \
                             (signed distances {m0:e}, {m1:e}; slack {slack:e})"
                        )));
                    }
                };
                let (local, new_face, new_slot) = if use0 { (c0, f0, a) } else { (c1, f1, b) };
                let clamped = local.iter().any(|&x| x < 0.0);
                let local = if clamped {
                    let x = local.map(|x| x.max(0.0));
                    let s = x[0] + x[1] + x[2];
                    x.map(|x| x / s)
                } else {
                    local
                };
                let mut out_bary = [0.0; 3];
                for (k, &x) in local.iter().enumerate() {
                    out_bary[rot(new_slot, k)] = x;
                }
                out.push((p, new_face as u32, out_bary, clamped));
            }
        }
        Ok(out)
    }

    /// Flips every interior edge to Delaunay (SPEC §3.3; Bobenko–Springborn's greedy
    /// algorithm with the options' criterion and flip limit), tracking the points.
    ///
    /// # Errors
    /// [`IemError::RemovalPending`]; [`IemError::Internal`] from a flip (the flips before it
    /// stay).
    pub fn flip_to_delaunay(&mut self) -> Result<DelaunayFlipReport, IemError> {
        let seeds: Vec<u32> = (0..self.n_edges() as u32).collect();
        self.flip_to_delaunay_from(&seeds)
    }

    /// [`flip_to_delaunay`](Self::flip_to_delaunay) from the given edges (internal ids): the
    /// queue starts with them, and every flip pushes the four outer edges of its diamond.
    ///
    /// # Errors
    /// As [`flip_to_delaunay`](Self::flip_to_delaunay), and [`IemError::IndexOutOfRange`].
    pub fn flip_to_delaunay_from(&mut self, edges: &[u32]) -> Result<DelaunayFlipReport, IemError> {
        if self.pending.is_some() {
            return Err(IemError::RemovalPending);
        }
        if let Some(&e) = edges.iter().find(|&&e| e as usize >= self.n_edges()) {
            return Err(IemError::IndexOutOfRange {
                what: "edge",
                index: e as usize,
                bound: self.n_edges(),
            });
        }
        self.delaunay_pass(edges, true)
    }

    /// The Delaunay pass of [`flip_to_delaunay_from`](Self::flip_to_delaunay_from).
    pub(crate) fn delaunay_pass(
        &mut self,
        seeds: &[u32],
        track: bool,
    ) -> Result<DelaunayFlipReport, IemError> {
        let epoch = self.next_epoch();
        let mut queue: VecDeque<u32> = VecDeque::with_capacity(seeds.len());
        for &e in seeds {
            if self.edge_mark[e as usize] != epoch {
                self.edge_mark[e as usize] = epoch;
                queue.push_back(e);
            }
        }
        let limit = match self.options.max_delaunay_flips {
            FlipLimit::PerEdge(k) => k.saturating_mul(self.n_edges()),
            FlipLimit::Flips(n) => n,
            FlipLimit::Unlimited => usize::MAX,
        };
        let mut report = DelaunayFlipReport {
            status: DelaunayStatus::Delaunay,
            flips: Vec::new(),
            n_tested: 0,
            n_refused: 0,
            n_not_delaunay: 0,
            n_clamped: 0,
        };
        let mut refused: Vec<u32> = Vec::new();
        while let Some(e) = queue.pop_front() {
            let e = e as usize;
            self.edge_mark[e] = 0;
            if self.complex.is_boundary_edge(e) {
                continue;
            }
            report.n_tested += 1;
            let Some(excess) = self.delaunay_excess(e) else {
                continue;
            };
            if self.is_delaunay_with(e, excess) {
                continue;
            }
            if report.flips.len() >= limit {
                report.status = DelaunayStatus::FlipLimit;
                break;
            }
            match self.flip_internal(e, track) {
                Ok((flipped, clamped)) => {
                    report.flips.push(flipped);
                    report.n_clamped += clamped;
                    let c = &self.complex;
                    let h0 = c.edge_halfedge(e);
                    let h1 = c.twin(h0).unwrap_or(h0);
                    for h in [c.next(h0), c.prev(h0), c.next(h1), c.prev(h1)] {
                        let g = c.edge(h);
                        if self.edge_mark[g] != epoch {
                            self.edge_mark[g] = epoch;
                            queue.push_back(g as u32);
                        }
                    }
                }
                Err(FlipFailure::Refused(_)) => {
                    report.n_refused += 1;
                    refused.push(e as u32);
                }
                Err(FlipFailure::Internal(what)) => return Err(IemError::Internal { what }),
            }
        }
        if report.status == DelaunayStatus::Delaunay {
            refused.sort_unstable();
            refused.dedup();
            for &e in &refused {
                if let Some(excess) = self.delaunay_excess(e as usize)
                    && !self.is_delaunay_with(e as usize, excess)
                {
                    report.n_not_delaunay += 1;
                }
            }
            if report.n_not_delaunay > 0 {
                report.status = DelaunayStatus::Blocked;
            }
        }
        Ok(report)
    }
}
