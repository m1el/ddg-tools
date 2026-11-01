//! Vertex removal (paper Sec. 4.2; SPEC §3.5), tentative removal on the undo log (§3.6, S5).
//!
//! A removal runs in two phases. [`CoarseMesh::begin_removal`] flattens the vertex (Springborn's
//! flips if needed), flips it down to degree three (two on the boundary; an ear first gets the
//! flip of its opposite edge) and checks that the removal can complete; every write is
//! recorded. It returns a [`PreparedRemoval`] with what the error metric needs (the
//! neighbours, their curvature before and after, the edges and frames after flattening, the
//! flips). Then either [`CoarseMesh::revert_removal`] restores the mesh bitwise, or
//! [`CoarseMesh::finish_removal`] removes the vertex from the complex, re-expresses the points
//! of its star in the merged face, places the vertex there, and flips the edges among its
//! neighbours to Delaunay. A tentative removal ([`CoarseMesh::tentative_removal`]) is the
//! first phase without point tracking, then the revert: the same code computes the same
//! [`PreparedRemoval`] as the committed removal (S5).

use std::f64::consts::{PI, TAU};

use ddg_core::geometry::triangle_area;
use ddg_core::topology::NONE;
use ddg_intrinsic::plan_flip_with;

use crate::error::{IemError, RemovalFailure, RemovalRefusal};
use crate::flatten::{FlattenReport, StarFace};
use crate::flip::{DelaunayFlipReport, FlipFailure, FlippedEdge};
use crate::formulas::{boundary_merge, degree_three_barycentric};
use crate::mesh::{CoarseMesh, Location, ccw, corner_at, face_lengths};

/// What kind of vertex a removal removes (paper Sec. 4).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RemovalKind {
    /// An interior vertex: flattened, flipped to degree three, three faces merged.
    Interior,
    /// A boundary vertex on two or more faces: flattened, flipped to two faces, merged.
    Boundary,
    /// A boundary vertex on one face: its opposite edge flipped first, then as `Boundary`.
    Ear,
}

/// A neighbour `j` of the removed vertex `i` and what flattening `i` did to it: the data of
/// the transport (SPEC §4.2, §4.3), all after flattening.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NeighbourChange {
    /// The neighbour (fine id).
    pub vertex: u32,
    /// Whether it lies on the boundary (its curvature is then `κ`).
    pub boundary: bool,
    /// Its curvature `K_j` (or `κ_j`) before flattening.
    pub curvature_before: f64,
    /// Its curvature `K̃_j` (or `κ̃_j`) after flattening.
    pub curvature_after: f64,
    /// The length `ℓ_ij` of the shortest edge between `i` and `j` after flattening (ties: the
    /// lower edge id; S12).
    pub length: f64,
    /// The angular coordinate of that edge at `i`, `φ_ij` (`π` if it is the boundary edge
    /// arriving at `i`).
    pub phi_ij: f64,
    /// The angular coordinate of that edge at `j`, `φ_ji` (`π` if it is the boundary edge
    /// arriving at `j`).
    pub phi_ji: f64,
    /// `φ_ij` of the same edge before flattening (the frame the reference's tentative cost
    /// uses, R16, and its commit at `i`, R19; for `ReferenceCompat`).
    pub phi_ij_before: f64,
    /// `φ_ji` of the same edge before flattening (R16).
    pub phi_ji_before: f64,
}

impl NeighbourChange {
    /// The transport of a vector from `i`'s normalized frame to `j`'s along the edge (paper
    /// Sec. 3.2.2): `R_ij = e^{i(φ_ji + π − φ_ij)}`, as `(cos, sin)`.
    ///
    /// [formula: iem-transport-rotation]
    pub fn rotation(&self) -> [f64; 2] {
        let a = self.phi_ji + PI - self.phi_ij;
        [a.cos(), a.sin()]
    }

    /// The edge vector `e_ji` from `j` to `i` in `j`'s normalized frame, `ℓ_ij e^{i φ_ji}`.
    pub fn edge_vector(&self) -> [f64; 2] {
        [
            self.length * self.phi_ji.cos(),
            self.length * self.phi_ji.sin(),
        ]
    }
}

/// What the first phase of a removal did ([`CoarseMesh::begin_removal`]).
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedRemoval {
    /// The vertex (fine id).
    pub vertex: u32,
    /// What kind of vertex it is.
    pub kind: RemovalKind,
    /// Its curvature `K_i` (or `κ_i`) before flattening.
    pub curvature: f64,
    /// The flattening.
    pub flatten: FlattenReport,
    /// The flip of the opposite edge of an ear.
    pub ear_flip: Option<FlippedEdge>,
    /// Springborn's fallback flips (SPEC §3.4), in order.
    pub fallback_flips: Vec<FlippedEdge>,
    /// The distinct neighbours `N_i` at flattening (fine ids ascending) and their changes.
    pub neighbours: Vec<NeighbourChange>,
    /// The flips to degree three (two on the boundary), in order (S14).
    pub degree_flips: Vec<FlippedEdge>,
    /// Points clamped onto a side by these flips (S16; always 0 for a tentative removal).
    pub n_clamped: usize,
}

impl PreparedRemoval {
    /// Every flip of the first phase, in order: the ear flip, the fallback flips, the flips to
    /// degree three or two.
    pub fn flips(&self) -> impl Iterator<Item = &FlippedEdge> {
        self.ear_flip
            .iter()
            .chain(&self.fallback_flips)
            .chain(&self.degree_flips)
    }
}

/// What a completed removal did ([`CoarseMesh::finish_removal`]).
#[derive(Debug, Clone, PartialEq)]
pub struct RemovalReport {
    /// The first phase.
    pub prepared: PreparedRemoval,
    /// The merged face (internal id after the removal).
    pub face: u32,
    /// Its corners (fine ids, in its order).
    pub corners: [u32; 3],
    /// The coordinates of the removed vertex in it (`iem-degree-three-barycentric`, or on the
    /// boundary side by `iem-boundary-merge`).
    pub bary: [f64; 3],
    /// The Delaunay flips among the neighbours after the removal (S15), if enabled.
    pub delaunay: Option<DelaunayFlipReport>,
}

/// A removal between its two phases.
#[derive(Debug, Clone)]
pub(crate) struct Pending {
    vertex: usize,
    fine: u32,
    kind: RemovalKind,
    track: bool,
    prepared: PreparedRemoval,
    /// The clamp counter before the removal began (restored by a revert).
    n_clamped_before: u64,
}

/// The renumbering a removal from the complex reports.
struct Renumbering {
    halfedge_moves: Vec<(u32, u32)>,
    face_moves: Vec<(u32, u32)>,
    vertex_move: Option<(u32, u32)>,
    counts: [u32; 4],
}

fn permute<T: Clone + Default>(data: &mut Vec<T>, moves: &[(u32, u32)], len: u32) {
    let values: Vec<(usize, T)> = moves
        .iter()
        .map(|&(old, new)| (new as usize, std::mem::take(&mut data[old as usize])))
        .collect();
    for (new, value) in values {
        data[new] = value;
    }
    data.truncate(len as usize);
}

impl CoarseMesh {
    // ----- public API -----

    /// The first phase of removing fine vertex `v` with point tracking (see the
    /// [crate documentation](crate)). Afterwards a removal is pending: only
    /// [`finish_removal`](Self::finish_removal) or [`revert_removal`](Self::revert_removal)
    /// may change the mesh.
    ///
    /// # Errors
    /// A [`RemovalFailure`] (SPEC §3.5: the skip rules, an ear that cannot be flipped, a
    /// flattening that fails, a vertex stuck above degree three); the mesh is unchanged
    /// (bitwise).
    pub fn begin_removal(&mut self, v: u32) -> Result<PreparedRemoval, RemovalFailure> {
        self.begin(v, true)
    }

    /// Reverts the pending removal: every array is restored bitwise.
    ///
    /// # Errors
    /// [`IemError::NoRemovalPending`].
    pub fn revert_removal(&mut self) -> Result<(), IemError> {
        let Some(pending) = self.pending.take() else {
            return Err(IemError::NoRemovalPending);
        };
        self.n_clamped = pending.n_clamped_before;
        self.replay_undo();
        self.undo.end();
        Ok(())
    }

    /// The cost data of removing fine vertex `v`, without changing the mesh (SPEC §3.6): the
    /// first phase without point tracking, then the revert.
    ///
    /// # Errors
    /// As [`begin_removal`](Self::begin_removal).
    pub fn tentative_removal(&mut self, v: u32) -> Result<PreparedRemoval, RemovalFailure> {
        let prepared = self.begin(v, false)?;
        if let Some(pending) = self.pending.take() {
            self.n_clamped = pending.n_clamped_before;
        }
        self.replay_undo();
        self.undo.end();
        Ok(prepared)
    }

    /// Removes fine vertex `v`: [`begin_removal`](Self::begin_removal), then
    /// [`finish_removal`](Self::finish_removal).
    ///
    /// # Errors
    /// [`IemError::Removal`] if the vertex cannot be removed (the mesh is unchanged), and the
    /// errors of `finish_removal`.
    pub fn remove_vertex(&mut self, v: u32) -> Result<RemovalReport, IemError> {
        self.begin(v, true).map_err(IemError::Removal)?;
        self.finish_removal()
    }

    /// The second phase of the pending removal (see the [crate documentation](crate)).
    ///
    /// # Errors
    /// [`IemError::NoRemovalPending`]; [`IemError::Internal`] if the complex refuses the
    /// removal that the first phase checked (the removal is then reverted) or a Delaunay flip
    /// meets a point outside its faces (the removal stays; the flips before it too).
    pub fn finish_removal(&mut self) -> Result<RemovalReport, IemError> {
        let Some(pending) = self.pending.take() else {
            return Err(IemError::NoRemovalPending);
        };
        if !pending.track {
            self.n_clamped = pending.n_clamped_before;
            self.replay_undo();
            self.undo.end();
            return Err(IemError::Internal {
                what: "a removal begun without point tracking cannot be finished; it was \
                       reverted"
                    .into(),
            });
        }
        let result = match pending.kind {
            RemovalKind::Interior => self.remove_interior_now(pending.vertex, pending.fine),
            RemovalKind::Boundary | RemovalKind::Ear => {
                self.remove_boundary_now(pending.vertex, pending.fine)
            }
        };
        let (face, bary) = match result {
            Ok(x) => x,
            Err(what) => {
                self.n_clamped = pending.n_clamped_before;
                self.replay_undo();
                self.undo.end();
                return Err(IemError::Internal { what });
            }
        };
        self.undo.end();
        let corners = self.faces_fine(face as usize);
        let delaunay = if self.options.delaunay_after_removal {
            let fine: Vec<u32> = pending
                .prepared
                .neighbours
                .iter()
                .map(|n| n.vertex)
                .collect();
            let seeds = self.edges_among(&fine);
            Some(self.delaunay_pass(&seeds, true)?)
        } else {
            None
        };
        Ok(RemovalReport {
            prepared: pending.prepared,
            face,
            corners,
            bary,
            delaunay,
        })
    }

    /// Whether a removal is pending.
    pub fn removal_pending(&self) -> bool {
        self.pending.is_some()
    }

    // ----- phase one -----

    fn begin(&mut self, fine: u32, track: bool) -> Result<PreparedRemoval, RemovalFailure> {
        let failure = |reason| RemovalFailure {
            vertex: fine,
            reason,
        };
        if self.pending.is_some() {
            return Err(failure(RemovalRefusal::Internal {
                what: "another removal is pending".into(),
            }));
        }
        let Some(v) = self.complex_vertex(fine) else {
            return Err(failure(RemovalRefusal::NotAVertex));
        };
        if self.is_fixed(fine) {
            return Err(failure(RemovalRefusal::Fixed));
        }
        let boundary = self.complex.is_boundary_vertex(v);
        let fan: Vec<usize> = self.fan(v).collect();
        for &h in &fan {
            let f = self.complex.faces()[h / 3];
            if f[0] == f[1] && f[1] == f[2] {
                return Err(failure(RemovalRefusal::SelfFace {
                    face: (h / 3) as u32,
                }));
            }
        }
        if boundary && self.complex.to_vertex(self.reference[v] as usize) == v {
            return Err(failure(RemovalRefusal::BoundarySelfEdge));
        }
        if !boundary && fan.len() < 3 {
            return Err(failure(RemovalRefusal::InteriorDegree { faces: fan.len() }));
        }
        let kind = if !boundary {
            RemovalKind::Interior
        } else if fan.len() == 1 {
            RemovalKind::Ear
        } else {
            RemovalKind::Boundary
        };
        self.undo.begin();
        let n_clamped_before = self.n_clamped;
        match self.prepare(v, fine, kind, track) {
            Ok(prepared) => {
                self.pending = Some(Pending {
                    vertex: v,
                    fine,
                    kind,
                    track,
                    prepared: prepared.clone(),
                    n_clamped_before,
                });
                Ok(prepared)
            }
            Err(reason) => {
                self.n_clamped = n_clamped_before;
                self.replay_undo();
                self.undo.end();
                Err(failure(reason))
            }
        }
    }

    fn flip_for_removal(
        &mut self,
        e: usize,
        track: bool,
        clamped: &mut usize,
    ) -> Result<Option<FlippedEdge>, RemovalRefusal> {
        match self.flip_internal(e, track) {
            Ok((f, n)) => {
                *clamped += n;
                Ok(Some(f))
            }
            Err(FlipFailure::Refused(_)) => Ok(None),
            Err(FlipFailure::Internal(what)) => Err(RemovalRefusal::Internal { what }),
        }
    }

    fn prepare(
        &mut self,
        v: usize,
        fine: u32,
        kind: RemovalKind,
        track: bool,
    ) -> Result<PreparedRemoval, RemovalRefusal> {
        let boundary = kind != RemovalKind::Interior;
        let curvature = self.curvature(v);
        let mut n_clamped = 0;
        // An ear: flip the edge opposite it first (paper Sec. 4, "Special cases").
        let mut ear_flip = None;
        if kind == RemovalKind::Ear {
            let h = self.reference[v] as usize;
            let e = self.complex.edge(self.complex.next(h));
            match self.flip_internal(e, track) {
                Ok((f, n)) => {
                    ear_flip = Some(f);
                    n_clamped += n;
                }
                Err(FlipFailure::Refused(reason)) => {
                    return Err(RemovalRefusal::EarNotFlippable { reason });
                }
                Err(FlipFailure::Internal(what)) => return Err(RemovalRefusal::Internal { what }),
            }
        }
        // Flatten, with Springborn's flips while the flattened star is not valid.
        let target = if boundary { PI } else { TAU };
        let max_attempts = self.complex.vertex_degree(v) + 3;
        let mut report = FlattenReport::default();
        let mut fallback_flips = Vec::new();
        let (faces, solution) = loop {
            let faces = self.star_faces(v);
            report.attempts += 1;
            let solution = self.solve_flattening(&faces, target).map_err(|f| {
                RemovalRefusal::FlattenFailed {
                    residual: f.residual,
                    iterations: report.iterations + f.iterations,
                }
            })?;
            report.iterations += solution.iterations;
            report.bisections += solution.bisections;
            if solution.valid {
                break (faces, solution);
            }
            if report.attempts >= max_attempts {
                return Err(RemovalRefusal::NotFlattenable {
                    fallback_flips: fallback_flips.len(),
                });
            }
            let mut flipped = None;
            for e in self.fallback_candidates(&faces, solution.u) {
                if let Some(f) = self.flip_for_removal(e, track, &mut n_clamped)? {
                    flipped = Some(f);
                    break;
                }
            }
            match flipped {
                Some(f) => fallback_flips.push(f),
                None => {
                    return Err(RemovalRefusal::NotFlattenable {
                        fallback_flips: fallback_flips.len(),
                    });
                }
            }
        };
        report.u = solution.u;
        report.residual = solution.residual;

        // The neighbours and their curvature, before and after.
        let neighbours = self.star_neighbours(v, &faces);
        let before: Vec<f64> = neighbours.iter().map(|&w| self.curvature(w)).collect();
        // The angular coordinates of the spokes (both ends) before flattening.
        let mut phi_before: Vec<(u32, f64)> = Vec::with_capacity(16);
        {
            let c = &self.complex;
            let mut last = None;
            for h in self.fan(v) {
                last = Some(h);
                phi_before.push((h as u32, self.phi[h]));
                if let Some(t) = c.twin(h) {
                    phi_before.push((t as u32, self.phi[t]));
                }
            }
            if c.is_boundary_vertex(v)
                && let Some(h) = last
            {
                let g = c.prev(h);
                phi_before.push((g as u32, self.phi[g]));
            }
        }
        let lookup = |h: usize| {
            if h == NONE as usize {
                PI
            } else {
                phi_before
                    .iter()
                    .find(|x| x.0 as usize == h)
                    .map_or(f64::NAN, |x| x.1)
            }
        };
        self.commit_flattening(v, &faces, solution.u, track);
        let mut changes: Vec<NeighbourChange> = neighbours
            .iter()
            .zip(&before)
            .map(|(&w, &kb)| {
                let (length, phi_ij, phi_ji, h_ij, h_ji) = self.connecting_edge(v, w);
                NeighbourChange {
                    vertex: self.fine_of[w],
                    boundary: self.complex.is_boundary_vertex(w),
                    curvature_before: kb,
                    curvature_after: self.curvature(w),
                    length,
                    phi_ij,
                    phi_ji,
                    phi_ij_before: lookup(h_ij),
                    phi_ji_before: lookup(h_ji),
                }
            })
            .collect();
        changes.sort_by_key(|n| n.vertex);

        // Flip to degree three (two on the boundary), S14.
        let wanted = if boundary { 2 } else { 3 };
        let mut degree_flips = Vec::new();
        let cap = 4 * self.complex.vertex_degree(v) + 8;
        loop {
            let fan: Vec<usize> = self.fan(v).collect();
            let self_edge = fan.iter().any(|&h| self.complex.to_vertex(h) == v);
            if fan.len() == wanted && !self_edge {
                break;
            }
            if fan.len() <= wanted {
                return Err(if self_edge {
                    RemovalRefusal::SelfEdge
                } else {
                    RemovalRefusal::Stuck { faces: fan.len() }
                });
            }
            if degree_flips.len() >= cap {
                return Err(RemovalRefusal::Stuck { faces: fan.len() });
            }
            let Some(e) = self.degree_flip_choice(v, &fan) else {
                return Err(RemovalRefusal::Stuck { faces: fan.len() });
            };
            match self.flip_for_removal(e, track, &mut n_clamped)? {
                Some(f) => degree_flips.push(f),
                None => {
                    return Err(RemovalRefusal::Internal {
                        what: format!("edge {e} was flippable and then refused"),
                    });
                }
            }
        }

        // The removal must be able to complete.
        let c = &self.complex;
        let lengths = if boundary {
            if let Some(r) = c.boundary_removal_refusal(v) {
                return Err(match r {
                    ddg_intrinsic::RemovalRefusal::SelfEdge { .. } => RemovalRefusal::SelfEdge,
                    other => RemovalRefusal::Internal {
                        what: format!("the complex refuses the boundary removal: {other}"),
                    },
                });
            }
            let h_vc = c.outgoing_halfedge(v);
            let h_va = ccw(c, h_vc).unwrap_or(h_vc);
            let l = |h: usize| self.length[c.edge(h)];
            let (_, l_bc) = boundary_merge(l(c.prev(h_va)), l(h_vc));
            [l_bc, l(c.next(h_vc)), l(c.next(h_va))]
        } else {
            let o: Vec<usize> = c.outgoing_halfedges(v).collect();
            [0, 1, 2].map(|k| self.length[c.edge(c.next(o[k]))])
        };
        if triangle_area(lengths[0], lengths[1], lengths[2]).is_none() {
            return Err(RemovalRefusal::Degenerate { lengths });
        }
        Ok(PreparedRemoval {
            vertex: fine,
            kind,
            curvature,
            flatten: report,
            ear_flip,
            fallback_flips,
            neighbours: changes,
            degree_flips,
            n_clamped,
        })
    }

    /// The distinct vertices other than `v` of the star `faces`, in order of appearance.
    fn star_neighbours(&mut self, v: usize, faces: &[StarFace]) -> Vec<usize> {
        let epoch = self.next_epoch();
        let mut out = Vec::new();
        for sf in faces {
            for &x in &self.complex.faces()[sf.face] {
                let x = x as usize;
                if x != v && self.vertex_mark[x] != epoch {
                    self.vertex_mark[x] = epoch;
                    out.push(x);
                }
            }
        }
        out
    }

    /// The shortest edge between `v` and its neighbour `w` (ties: lower edge id): its length,
    /// its angular coordinate at `v` and at `w` (`π` for the virtual direction of an arriving
    /// boundary edge), and the halfedges these coordinates belong to (`NONE` for the virtual
    /// direction).
    fn connecting_edge(&self, v: usize, w: usize) -> (f64, f64, f64, usize, usize) {
        let c = &self.complex;
        let none = NONE as usize;
        let mut best: Option<(f64, usize, usize, usize)> = None;
        let mut consider = |e: usize, h_v: usize, h_w: usize| {
            let l = self.length[e];
            let better = match best {
                None => true,
                Some((bl, be, _, _)) => l < bl || (l == bl && e < be),
            };
            if better {
                best = Some((l, e, h_v, h_w));
            }
        };
        let mut last = None;
        for h in self.fan(v) {
            last = Some(h);
            if c.to_vertex(h) != w {
                continue;
            }
            consider(c.edge(h), h, c.twin(h).unwrap_or(none));
        }
        if c.is_boundary_vertex(v)
            && let Some(h) = last
        {
            // The boundary edge arriving at v.
            let g = c.prev(h);
            if c.from_vertex(g) == w {
                consider(c.edge(g), none, g);
            }
        }
        let phi = |h: usize| if h == none { PI } else { self.phi[h] };
        match best {
            Some((l, _, a, b)) => (l, phi(a), phi(b), a, b),
            None => (f64::NAN, f64::NAN, f64::NAN, none, none),
        }
    }

    /// The spoke of `v` to flip next on the way to degree three or two (S14): among the
    /// flippable spokes, self-edges first, then the largest opposite-angle sum, ties by the
    /// lower edge id. `None` if no spoke can be flipped.
    fn degree_flip_choice(&self, v: usize, fan: &[usize]) -> Option<usize> {
        let c = &self.complex;
        // Scores first (cheap), then the flip test in the order of preference: the first
        // flippable spoke is the best flippable one.
        let mut candidates: Vec<(bool, f64, usize)> = Vec::with_capacity(fan.len());
        for &h in fan {
            let e = c.edge(h);
            let Some(t) = c.twin(h) else { continue };
            let is_self = c.to_vertex(h) == v;
            if is_self && c.edge_halfedge(e) != h {
                continue; // each self-edge once
            }
            let score = corner_at(c, &self.length, c.prev(h)).unwrap_or(0.0)
                + corner_at(c, &self.length, c.prev(t)).unwrap_or(0.0);
            candidates.push((is_self, score, e));
        }
        candidates.sort_by(|a, b| b.0.cmp(&a.0).then(b.1.total_cmp(&a.1)).then(a.2.cmp(&b.2)));
        candidates
            .into_iter()
            .find(|&(_, _, e)| plan_flip_with(c, &self.length, e, self.options.flip_test).is_ok())
            .map(|b| b.2)
    }

    // ----- phase two -----

    /// Removes interior vertex `v` of degree three from the complex; returns the merged face
    /// and the vertex's coordinates in it.
    fn remove_interior_now(&mut self, v: usize, fine: u32) -> Result<(u32, [f64; 3]), String> {
        let c = &self.complex;
        let o: Vec<usize> = c.outgoing_halfedges(v).collect();
        if o.len() != 3 {
            return Err(format!("vertex {fine} has {} faces, not three", o.len()));
        }
        let f = [o[0] / 3, o[1] / 3, o[2] / 3];
        let areas: Vec<f64> = f
            .iter()
            .map(|&g| {
                let l = face_lengths(c, &self.length, g);
                triangle_area(l[0], l[1], l[2]).unwrap_or(0.0)
            })
            .collect();
        // Corner m of the merged face is x_m; the sub-triangle opposite it is f_{m+1}.
        let bary = degree_three_barycentric([areas[1], areas[2], areas[0]])
            .ok_or_else(|| format!("vertex {fine}: a face of its star has no area"))?;
        let keep = f[0].min(f[1]).min(f[2]);
        // The points of the star, by linear substitution (S16).
        let mut merged: Vec<(u32, [f64; 3])> = Vec::new();
        for k in 0..3 {
            let (sv, sx, sy) = (o[k] % 3, c.next(o[k]) % 3, c.prev(o[k]) % 3);
            for &p in &self.face_points[f[k]] {
                if let Location::Face { bary: b, .. } = self.points[p as usize] {
                    let mut nb = bary.map(|x| x * b[sv]);
                    nb[k] += b[sx];
                    nb[(k + 1) % 3] += b[sy];
                    merged.push((p, nb));
                }
            }
        }
        // References that die pass to their clockwise neighbour: x_{k+1} → v (prev(o_k))
        // becomes x_{k+1} → x_{k+2} (next(o_{k+1})).
        let mut fixes = Vec::new();
        for k in 0..3 {
            let dying = c.prev(o[k]);
            let x = c.from_vertex(dying);
            if self.reference[x] as usize == dying {
                fixes.push((x, c.next(o[(k + 1) % 3]) as u32));
            }
        }
        let removal = self
            .complex
            .remove_vertex(v)
            .map_err(|e| format!("the complex refuses the removal of vertex {fine}: {e}"))?;
        for (x, h) in fixes {
            self.reference[x] = h;
        }
        for &g in &f {
            self.face_points[g].clear();
        }
        removal.apply_edge_data(&mut self.length);
        removal.apply_halfedge_data(&mut self.phi);
        let renumbering = Renumbering {
            halfedge_moves: removal.halfedge_moves.clone(),
            face_moves: removal.face_moves.clone(),
            vertex_move: removal.vertex_move,
            counts: removal.counts,
        };
        let face = self.apply_renumbering(&renumbering, fine, keep, merged, bary);
        Ok((face, bary))
    }

    /// Removes boundary vertex `v` on two faces from the complex; returns the merged face and
    /// the vertex's coordinates in it (on the merged boundary side).
    fn remove_boundary_now(&mut self, v: usize, fine: u32) -> Result<(u32, [f64; 3]), String> {
        let c = &self.complex;
        let h_vc = c.outgoing_halfedge(v);
        let h_va = ccw(c, h_vc).ok_or_else(|| format!("vertex {fine} has one face"))?;
        let (fa, fb) = (h_vc / 3, h_va / 3);
        let l = |h: usize| self.length[c.edge(h)];
        let (t, l_bc) = boundary_merge(l(c.prev(h_va)), l(h_vc));
        let keep = fa.min(fb);
        // Slots of b, c, a in the merged face.
        let (sb, sc, sa) = if keep == fa {
            (h_vc % 3, c.next(h_vc) % 3, c.prev(h_vc) % 3)
        } else {
            (c.prev(h_va) % 3, h_va % 3, c.next(h_va) % 3)
        };
        let mut at_v = [0.0; 3];
        at_v[sb] += 1.0 - t;
        at_v[sc] += t;
        let mut merged: Vec<(u32, [f64; 3])> = Vec::new();
        for (g, h, s_next, s_prev) in [(fa, h_vc, sc, sa), (fb, h_va, sa, sb)] {
            let (sv, sx, sy) = (h % 3, c.next(h) % 3, c.prev(h) % 3);
            for &p in &self.face_points[g] {
                if let Location::Face { bary: b, .. } = self.points[p as usize] {
                    let mut nb = at_v.map(|x| x * b[sv]);
                    nb[s_next] += b[sx];
                    nb[s_prev] += b[sy];
                    merged.push((p, nb));
                }
            }
        }
        // a → v dies: a's reference passes to a → b.
        let dying = c.prev(h_vc);
        let a = c.from_vertex(dying);
        let fix = (self.reference[a] as usize == dying).then_some((a, c.next(h_va) as u32));
        let merge = self
            .complex
            .remove_boundary_vertex(v)
            .map_err(|e| format!("the complex refuses the removal of vertex {fine}: {e}"))?;
        if let Some((x, h)) = fix {
            self.reference[x] = h;
        }
        self.face_points[fa].clear();
        self.face_points[fb].clear();
        self.length[merge.halves[0] as usize] = l_bc;
        merge.apply_edge_data(&mut self.length);
        merge.apply_halfedge_data(&mut self.phi);
        let renumbering = Renumbering {
            halfedge_moves: merge.halfedge_moves.clone(),
            face_moves: merge.face_moves.clone(),
            vertex_move: merge.vertex_move,
            counts: merge.counts,
        };
        let face = self.apply_renumbering(&renumbering, fine, keep, merged, at_v);
        Ok((face, at_v))
    }

    /// Brings the per-vertex, per-face and point data along a removal from the complex (the
    /// lengths and angular coordinates are done by the caller), fills the merged face `keep`
    /// (old id) with `merged`, places the removed vertex, and refreshes the merged face's
    /// corners. Returns the merged face's new id.
    fn apply_renumbering(
        &mut self,
        r: &Renumbering,
        fine: u32,
        keep: usize,
        merged: Vec<(u32, [f64; 3])>,
        vertex_bary: [f64; 3],
    ) -> u32 {
        // Vertices.
        let moves: Vec<(u32, u32)> = r.vertex_move.into_iter().collect();
        let n_vertices = r.counts[3];
        permute(&mut self.reference, &moves, n_vertices);
        permute(&mut self.angle_sum, &moves, n_vertices);
        permute(&mut self.fine_of, &moves, n_vertices);
        self.complex_of[fine as usize] = NONE;
        if let Some((_, new)) = r.vertex_move {
            let moved_fine = self.fine_of[new as usize];
            self.complex_of[moved_fine as usize] = new;
        }
        // Reference halfedges that moved (read all before writing).
        let updates: Vec<(usize, u32)> = r
            .halfedge_moves
            .iter()
            .filter_map(|&(old, new)| {
                let t = self.complex.from_vertex(new as usize);
                (self.reference[t] == old).then_some((t, new))
            })
            .collect();
        for (t, new) in updates {
            self.reference[t] = new;
        }
        // Faces.
        let mut keep_points: Vec<u32> = merged.iter().map(|m| m.0).collect();
        for (p, bary) in merged {
            self.points[p as usize] = Location::Face {
                face: keep as u32,
                bary,
            };
        }
        if (fine as usize) < self.n_input as usize {
            self.points[fine as usize] = Location::Face {
                face: keep as u32,
                bary: vertex_bary,
            };
            keep_points.push(fine);
        }
        self.face_points[keep] = keep_points;
        permute(&mut self.face_points, &r.face_moves, r.counts[2]);
        let new_face = r
            .face_moves
            .iter()
            .find(|m| m.0 as usize == keep)
            .map_or(keep as u32, |m| m.1);
        let mut moved_faces: Vec<u32> = r.face_moves.iter().map(|m| m.1).collect();
        moved_faces.push(new_face);
        moved_faces.sort_unstable();
        moved_faces.dedup();
        for f in moved_faces {
            for slot in 0..self.face_points[f as usize].len() {
                let p = self.face_points[f as usize][slot] as usize;
                self.point_slot[p] = slot as u32;
                if let Location::Face { bary, .. } = self.points[p] {
                    self.points[p] = Location::Face { face: f, bary };
                }
            }
        }
        // The corners of the merged face: angle sums and angular coordinates.
        let corners = self.complex.faces()[new_face as usize];
        let mut seen: Vec<usize> = Vec::with_capacity(3);
        for x in corners {
            let x = x as usize;
            if !seen.contains(&x) {
                seen.push(x);
                self.refresh_vertex(x);
            }
        }
        new_face
    }

    /// The interior edges whose two end points are among the fine vertices `fine` (S15).
    fn edges_among(&mut self, fine: &[u32]) -> Vec<u32> {
        let epoch = self.next_epoch();
        let vertices: Vec<usize> = fine
            .iter()
            .filter_map(|&f| self.complex_vertex(f))
            .collect();
        for &w in &vertices {
            self.vertex_mark[w] = epoch;
        }
        let mut seeds = Vec::new();
        for &w in &vertices {
            for h in self.fan(w) {
                let x = self.complex.to_vertex(h);
                if self.vertex_mark[x] == epoch && self.complex.twin(h).is_some() {
                    seeds.push(self.complex.edge(h) as u32);
                }
            }
        }
        seeds.sort_unstable();
        seeds.dedup();
        seeds
    }
}
