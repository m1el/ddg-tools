//! Successive self-parameterization for the extrinsic baseline (SPEC §10.4, S23): the map from
//! the fine vertices to the QEM-simplified mesh that Liu, Zhang, Ben-Chen, Jacobson,
//! *Surface Multigrid via Intrinsic Prolongation*, ACM TOG 40(4), 2021, Sec. 4 and App. C build
//! while an extrinsic decimation runs (the paper's Fig. 18 and 21 baseline; the method of
//! Lee et al. 1998 "MAPS" and Cohen et al. 2003 "successive mappings" with their joint
//! flattening). Read from `papers/liu2021_surface_multigrid.pdf`, Sec. 4.1–4.4 and App. C.
//!
//! # What is built
//!
//! Every tracked point (here: the fine vertices) is kept as a face and barycentric coordinates
//! on the *current* mesh. For each collapse of the edge `(u, v)` into `k` (the mesh
//! `M^l → M^{l+1}`), the patch `N^l(u, v)` (the faces around `u` and `v`) before and the patch
//! `N^{l+1}(k)` (the faces around `k`) after have the same boundary polygon (the link of the
//! edge). Both are flattened into one shared planar domain `U` by minimizing
//! `E(V^l, F^l, U) + E(V^{l+1}, F^{l+1}, U)` (Liu et al. Eq. 3–5), the boundary vertices being
//! one set of unknowns (the "joint variable"); `E` is the **least-squares conformal** energy
//! (Lévy et al. 2002; the variant Liu et al. found best, Fig. 14), with the usual two pinned
//! boundary vertices to remove the similarity (the energy is invariant under it). A point is
//! lifted from its face of `N^l` to `U` and dropped into the face of `N^{l+1}` that contains
//! it (Fig. 8). The maps of all collapses compose, since the points are re-expressed on every
//! collapse. Bijectivity needs both flattenings to be embeddings of the same polygon: it is
//! checked (all signed areas positive).
//!
//! # Gaps against the paper
//!
//! * **Closed meshes only.** The boundary edges of Sec. 4.3 (three colinearity cases, snapping
//!   an interior vertex to the boundary curve) are not implemented; [`crate::baseline::qem_simplify`]
//!   returns an error for a mesh with boundary when the successive map is requested.
//! * **Validity fallback.** LSCM with a free boundary can fold. Liu et al. refer to their
//!   earlier paper (Liu et al. 2020, App. B, not available here) for the validity test; ours:
//!   if the LSCM layout has a non-positive triangle, the polygon is put on a circle (chord
//!   length parameter) and the interior vertices solved by uniform Tutte weights, which is an
//!   embedding for any triangulated disc (Tutte 1963; Floater 2003). If that fails too (a
//!   degenerate patch) the collapse is **rejected** and counted ([`SuccessiveStats`]). The
//!   share of fallbacks is reported.
//! * **Energy.** Only LSCM (no ARAP); the 3D positions of the patch are the plane triangles of
//!   `E` (the coarse patch uses the collapse position `p`).
//! * **Points on the shared boundary** land on the boundary of an after-face; they go to the
//!   after-face with the largest minimal barycentric coordinate. A point whose best face has a
//!   coordinate below `−1e-9` (rounding of a valid embedding) is clamped and counted.

use std::collections::BTreeSet;

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// Counts of a successive-map run.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SuccessiveStats {
    /// Collapses whose patches were flattened by the joint LSCM.
    pub n_lscm: usize,
    /// Collapses that needed the Tutte fallback (LSCM folded).
    pub n_tutte: usize,
    /// Collapses rejected because neither layout was an embedding (or the patch was degenerate).
    pub n_rejected: usize,
    /// Points that fell outside every after-face by more than rounding and were clamped.
    pub n_clamped: usize,
    /// The largest number of points moved by one collapse.
    pub max_points_moved: usize,
    /// The total number of point relocations.
    pub n_relocations: usize,
}

/// A read-only view of the QEM working mesh.
pub(crate) struct View<'a> {
    pub pos: &'a [V3],
    pub faces: &'a [[u32; 3]],
    pub vfaces: &'a [BTreeSet<u32>],
}

#[derive(Debug)]
pub(crate) enum PlanError {
    /// A vertex of the patch is on the boundary (unsupported).
    Degenerate,
}

/// The effect of one collapse on the tracked points, computed before the mesh is changed.
pub(crate) struct Plan {
    /// `(point, face, bary)` in the faces of the collapsed mesh (same face ids).
    moves: Vec<(u32, u32, [f64; 3])>,
    /// The faces of the patch before the collapse (their point lists are rebuilt).
    patch: Vec<u32>,
    tutte: bool,
    clamped: usize,
}

/// The tracked points through successive collapses.
pub(crate) struct SuccessiveMap {
    /// Per point: its face on the current mesh and barycentric coordinates.
    loc: Vec<(u32, [f64; 3])>,
    /// Per face: the points in it.
    face_points: Vec<Vec<u32>>,
    pub stats: SuccessiveStats,
}

/// Solves `a x = b` (dense, `n × n` row-major) by Gaussian elimination with partial pivoting.
fn solve_dense(a: &mut [f64], b: &mut [f64], n: usize) -> bool {
    for c in 0..n {
        let mut piv = c;
        for r in c + 1..n {
            if a[r * n + c].abs() > a[piv * n + c].abs() {
                piv = r;
            }
        }
        let d = a[piv * n + c];
        if !(d.abs() > 1e-300) || !d.is_finite() {
            return false;
        }
        if piv != c {
            for k in 0..n {
                a.swap(c * n + k, piv * n + k);
            }
            b.swap(c, piv);
        }
        for r in c + 1..n {
            let f = a[r * n + c] / d;
            if f != 0.0 {
                for k in c..n {
                    a[r * n + k] -= f * a[c * n + k];
                }
                b[r] -= f * b[c];
            }
        }
    }
    for c in (0..n).rev() {
        let mut s = b[c];
        for k in c + 1..n {
            s -= a[c * n + k] * b[k];
        }
        b[c] = s / a[c * n + c];
    }
    b.iter().all(|x| x.is_finite())
}

/// The plane coordinates of a 3D triangle: first corner at the origin, second on the x axis,
/// third with `y ≥ 0` (counter-clockwise for the vertex order).
fn local_2d(p: [V3; 3]) -> Option<[[f64; 2]; 3]> {
    let e1 = sub(p[1], p[0]);
    let e2 = sub(p[2], p[0]);
    let l1 = norm(e1);
    if !(l1 > 0.0) {
        return None;
    }
    let x = dot(e1, e2) / l1;
    let y = norm(cross(e1, e2)) / l1;
    (y > 0.0 && x.is_finite() && y.is_finite()).then_some([[0.0, 0.0], [l1, 0.0], [x, y]])
}

struct PatchFace {
    idx: [usize; 3],
    xy: [[f64; 2]; 3],
}

fn signed_area(q: [[f64; 2]; 3]) -> f64 {
    0.5 * ((q[1][0] - q[0][0]) * (q[2][1] - q[0][1]) - (q[2][0] - q[0][0]) * (q[1][1] - q[0][1]))
}

impl SuccessiveMap {
    /// The fine vertices as points: vertex `v` in the first face that has it.
    pub fn new(n_vertices: usize, faces: &[[u32; 3]]) -> SuccessiveMap {
        let mut loc = vec![(u32::MAX, [0.0; 3]); n_vertices];
        let mut face_points = vec![Vec::new(); faces.len()];
        for (fi, f) in faces.iter().enumerate() {
            for (k, &v) in f.iter().enumerate() {
                if loc[v as usize].0 == u32::MAX {
                    let mut b = [0.0; 3];
                    b[k] = 1.0;
                    loc[v as usize] = (fi as u32, b);
                    face_points[fi].push(v);
                }
            }
        }
        SuccessiveMap {
            loc,
            face_points,
            stats: SuccessiveStats::default(),
        }
    }

    /// The location of every point.
    pub fn locations(&self) -> &[(u32, [f64; 3])] {
        &self.loc
    }

    /// Plans the collapse of `(u, v)` (v removed, `u` moves to `p`; `ef` the faces on the edge).
    /// `Err` rejects the collapse.
    pub fn plan(&self, w: &View, u: u32, v: u32, p: V3, ef: &[u32]) -> Result<Plan, PlanError> {
        let mut patch: Vec<u32> = w.vfaces[u as usize]
            .union(&w.vfaces[v as usize])
            .copied()
            .collect();
        patch.sort_unstable();
        // Ring vertices: indices 0..R; then u, v (before) and k (after).
        let mut ring: Vec<u32> = Vec::new();
        for &f in &patch {
            for &x in &w.faces[f as usize] {
                if x != u && x != v && !ring.contains(&x) {
                    ring.push(x);
                }
            }
        }
        let r = ring.len();
        let (iu, iv, ik) = (r, r + 1, r + 2);
        let n = r + 3;
        let idx_before = |x: u32| -> usize {
            if x == u {
                iu
            } else if x == v {
                iv
            } else {
                ring.iter().position(|&y| y == x).expect("ring vertex")
            }
        };
        let idx_after = |x: u32| -> usize {
            if x == u || x == v {
                ik
            } else {
                ring.iter().position(|&y| y == x).expect("ring vertex")
            }
        };
        let pos_before = |x: u32| w.pos[x as usize];
        let pos_after = |x: u32| {
            if x == u || x == v {
                p
            } else {
                w.pos[x as usize]
            }
        };
        let mut before = Vec::new();
        let mut after = Vec::new();
        for &f in &patch {
            let t = w.faces[f as usize];
            let xy = local_2d(t.map(pos_before)).ok_or(PlanError::Degenerate)?;
            before.push(PatchFace {
                idx: t.map(idx_before),
                xy,
            });
            if !ef.contains(&f) {
                let xy = local_2d(t.map(pos_after)).ok_or(PlanError::Degenerate)?;
                after.push((
                    f,
                    PatchFace {
                        idx: t.map(idx_after),
                        xy,
                    },
                ));
            }
        }
        let after_faces: Vec<&PatchFace> = after.iter().map(|(_, f)| f).collect();

        // --- joint LSCM ---
        let mut uv: Option<Vec<[f64; 2]>> = None;
        let mut tutte = false;
        if let Some(sol) = self.joint_lscm(w, &ring, n, &before, &after_faces)
            && valid(&sol, &before, &after_faces)
        {
            uv = Some(sol);
        }
        if uv.is_none() {
            tutte = true;
            if let Some(sol) = tutte_layout(w, u, v, &ring, n, &before, &after_faces)
                && valid(&sol, &before, &after_faces)
            {
                uv = Some(sol);
            }
        }
        let uv = uv.ok_or(PlanError::Degenerate)?;

        // --- relocate the points ---
        let mut moves = Vec::new();
        let mut clamped = 0;
        for (&f, bf) in patch.iter().zip(&before) {
            for &pt in &self.face_points[f as usize] {
                let b = self.loc[pt as usize].1;
                let mut q = [0.0; 2];
                for a in 0..3 {
                    q[0] += b[a] * uv[bf.idx[a]][0];
                    q[1] += b[a] * uv[bf.idx[a]][1];
                }
                let mut best: Option<(f64, u32, [f64; 3])> = None;
                for (g, af) in &after {
                    let t = af.idx.map(|i| uv[i]);
                    let area = signed_area(t);
                    let bb = [
                        signed_area([q, t[1], t[2]]) / area,
                        signed_area([t[0], q, t[2]]) / area,
                        signed_area([t[0], t[1], q]) / area,
                    ];
                    let m = bb[0].min(bb[1]).min(bb[2]);
                    if best.is_none_or(|(bm, _, _)| m > bm) {
                        best = Some((m, *g, bb));
                    }
                }
                let (m, g, mut bb) = best.ok_or(PlanError::Degenerate)?;
                if m < 0.0 {
                    if m < -1e-9 {
                        clamped += 1;
                    }
                    for x in bb.iter_mut() {
                        *x = x.max(0.0);
                    }
                    let s: f64 = bb.iter().sum();
                    for x in bb.iter_mut() {
                        *x /= s;
                    }
                }
                moves.push((pt, g, bb));
            }
        }
        Ok(Plan {
            moves,
            patch,
            tutte,
            clamped,
        })
    }

    /// Applies a plan (call after the collapse is committed on the mesh; face ids unchanged).
    pub fn commit(&mut self, plan: Plan) {
        for &f in &plan.patch {
            self.face_points[f as usize].clear();
        }
        self.stats.n_relocations += plan.moves.len();
        self.stats.max_points_moved = self.stats.max_points_moved.max(plan.moves.len());
        for (pt, g, b) in plan.moves {
            self.loc[pt as usize] = (g, b);
            self.face_points[g as usize].push(pt);
        }
        if plan.tutte {
            self.stats.n_tutte += 1;
        } else {
            self.stats.n_lscm += 1;
        }
        self.stats.n_clamped += plan.clamped;
    }

    fn joint_lscm(
        &self,
        w: &View,
        ring: &[u32],
        n: usize,
        before: &[PatchFace],
        after: &[&PatchFace],
    ) -> Option<Vec<[f64; 2]>> {
        let r = ring.len();
        if r < 2 {
            return None;
        }
        // Pins: the two ring vertices farthest apart in space.
        let (mut a, mut b, mut dmax) = (0, 1, -1.0);
        for i in 0..r {
            for j in i + 1..r {
                let d = norm(sub(w.pos[ring[i] as usize], w.pos[ring[j] as usize]));
                if d > dmax {
                    (a, b, dmax) = (i, j, d);
                }
            }
        }
        if !(dmax > 0.0) {
            return None;
        }
        let m = 2 * n;
        let mut h = vec![0.0; m * m];
        for f in before.iter().chain(after.iter().copied()) {
            let q = f.xy;
            let area = signed_area(q);
            if !(area > 0.0) {
                return None;
            }
            let g = [
                [
                    (q[1][1] - q[2][1]) / (2.0 * area),
                    (q[2][0] - q[1][0]) / (2.0 * area),
                ],
                [
                    (q[2][1] - q[0][1]) / (2.0 * area),
                    (q[0][0] - q[2][0]) / (2.0 * area),
                ],
                [
                    (q[0][1] - q[1][1]) / (2.0 * area),
                    (q[1][0] - q[0][0]) / (2.0 * area),
                ],
            ];
            let sw = area.sqrt();
            // Two residual rows over the six unknowns (U_a, V_a).
            for row in 0..2 {
                let mut cols = [0usize; 6];
                let mut co = [0.0; 6];
                for k in 0..3 {
                    let (gx, gy) = (g[k][0], g[k][1]);
                    cols[2 * k] = 2 * f.idx[k];
                    cols[2 * k + 1] = 2 * f.idx[k] + 1;
                    if row == 0 {
                        co[2 * k] = sw * gx;
                        co[2 * k + 1] = -sw * gy;
                    } else {
                        co[2 * k] = sw * gy;
                        co[2 * k + 1] = sw * gx;
                    }
                }
                for s in 0..6 {
                    for t in 0..6 {
                        h[cols[s] * m + cols[t]] += co[s] * co[t];
                    }
                }
            }
        }
        let pinned = [2 * a, 2 * a + 1, 2 * b, 2 * b + 1];
        let pin_val = [0.0, 0.0, dmax, 0.0];
        let free: Vec<usize> = (0..m).filter(|i| !pinned.contains(i)).collect();
        let nf = free.len();
        let trace: f64 = (0..m).map(|i| h[i * m + i]).sum();
        let mut hf = vec![0.0; nf * nf];
        let mut rhs = vec![0.0; nf];
        for (ri, &i) in free.iter().enumerate() {
            for (ci, &j) in free.iter().enumerate() {
                hf[ri * nf + ci] = h[i * m + j];
            }
            hf[ri * nf + ri] += 1e-12 * trace / m as f64;
            for (pk, &pj) in pinned.iter().enumerate() {
                rhs[ri] -= h[i * m + pj] * pin_val[pk];
            }
        }
        if !solve_dense(&mut hf, &mut rhs, nf) {
            return None;
        }
        let mut x = vec![0.0; m];
        for (pk, &pj) in pinned.iter().enumerate() {
            x[pj] = pin_val[pk];
        }
        for (ri, &i) in free.iter().enumerate() {
            x[i] = rhs[ri];
        }
        Some((0..n).map(|i| [x[2 * i], x[2 * i + 1]]).collect())
    }
}

/// All faces of both patches positively oriented with a non-negligible area.
fn valid(uv: &[[f64; 2]], before: &[PatchFace], after: &[&PatchFace]) -> bool {
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    for p in uv {
        for k in 0..2 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let diag2 = (hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2);
    if !(diag2 > 0.0 && diag2.is_finite()) {
        return false;
    }
    before
        .iter()
        .chain(after.iter().copied())
        .all(|f| signed_area(f.idx.map(|i| uv[i])) > 1e-12 * diag2)
}

/// The Tutte layout: the ring on the unit circle by chord length, the interior vertices
/// (`u`, `v` before, `k` after) at the mean of their neighbours.
fn tutte_layout(
    w: &View,
    u: u32,
    v: u32,
    ring: &[u32],
    n: usize,
    before: &[PatchFace],
    after: &[&PatchFace],
) -> Option<Vec<[f64; 2]>> {
    let r = ring.len();
    let (iu, iv, ik) = (r, r + 1, r + 2);
    // The cyclic order of the ring: the side of each face of the patch opposite its centre.
    let mut next = vec![usize::MAX; r];
    for f in before {
        let c: Vec<usize> = (0..3).filter(|&k| f.idx[k] >= r).collect();
        if c.len() == 1 {
            let k = c[0];
            next[f.idx[(k + 1) % 3]] = f.idx[(k + 2) % 3];
        }
    }
    let mut order = Vec::with_capacity(r);
    let mut cur = 0usize;
    for _ in 0..r {
        order.push(cur);
        cur = next[cur];
        if cur == usize::MAX {
            return None;
        }
    }
    if cur != 0 {
        return None;
    }
    let _ = (u, v);
    let mut total = 0.0;
    let mut cum = Vec::with_capacity(r);
    for k in 0..r {
        cum.push(total);
        total += norm(sub(
            w.pos[ring[order[k]] as usize],
            w.pos[ring[order[(k + 1) % r]] as usize],
        ));
    }
    if !(total > 0.0) {
        return None;
    }
    let mut uv = vec![[0.0; 2]; n];
    for k in 0..r {
        let t = std::f64::consts::TAU * cum[k] / total;
        uv[order[k]] = [t.cos(), t.sin()];
    }
    // Interior: iu, iv, ik with neighbours from the faces.
    let mut nb: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    for f in before.iter().chain(after.iter().copied()) {
        for k in 0..3 {
            let (a, b) = (f.idx[k], f.idx[(k + 1) % 3]);
            if a != b {
                nb[a].insert(b);
                nb[b].insert(a);
            }
        }
    }
    let inter = [iu, iv, ik];
    // iu and iv are in the before patch, ik in the after patch; all neighbour sets are
    // restricted to their own patch because the index sets differ (ik never meets iu or iv).
    let mut a = vec![0.0; 9];
    let mut bx = vec![0.0; 3];
    let mut by = vec![0.0; 3];
    for (row, &i) in inter.iter().enumerate() {
        let d = nb[i].len() as f64;
        if d == 0.0 {
            return None;
        }
        a[row * 3 + row] = d;
        for &j in &nb[i] {
            match inter.iter().position(|&x| x == j) {
                Some(col) => a[row * 3 + col] -= 1.0,
                None => {
                    bx[row] += uv[j][0];
                    by[row] += uv[j][1];
                }
            }
        }
    }
    let mut a2 = a.clone();
    if !solve_dense(&mut a, &mut bx, 3) || !solve_dense(&mut a2, &mut by, 3) {
        return None;
    }
    for (row, &i) in inter.iter().enumerate() {
        uv[i] = [bx[row], by[row]];
    }
    Some(uv)
}
