//! Dynamic remeshing: split, collapse, Delaunay flips, tangential smoothing (`SPEC.md` §10;
//! D16, D17; paper Sec. 6; ported in part from the reference's `remeshing.cpp`, MIT,
//! "Copyright (c) 2019 Nick Sharp", see `LICENSES/repulsive-surfaces-MIT.txt`).
//!
//! # One pass ([`remesh`])
//!
//! With `L₀` given by the caller (the mean edge length of the run's initial mesh, fixed for
//! the run):
//!
//! 1. **Edge lengths.** Every edge of the input longer than `split · L₀` (default `3L₀/2`) is
//!    split at its midpoint (in lexicographic order of the input edges; the halves are not
//!    split again in this pass). Then edges shorter than `collapse · L₀` (default `L₀/2`) are
//!    collapsed, shortest first (ties by vertex ids), each new short edge at a survivor queued
//!    again: to the pinned endpoint, else to the single boundary endpoint, else to the
//!    midpoint; refused by the link condition, both ends pinned, a repeated face, a
//!    non-manifold result, a foldover (a surviving face turning its normal by more than
//!    `π − 2` or becoming degenerate) and the optional guard (the CCD hook, `None` until WP8).
//! 2. `rounds` times (default 5): **Delaunay flips** — a FIFO queue of all edges in
//!    lexicographic order; an interior edge whose opposite angles sum to more than `π`
//!    (`rs-delaunay-test`, `f64`) is flipped unless `c d` is already an edge or a new face
//!    turns its normal by more than `π − 2` (refusals counted); the four edges of its quad are
//!    queued again; at most `flip_cap_per_vertex · |V|` flips per round (hitting the cap is
//!    reported) — and **tangential smoothing** — for every vertex that is neither pinned nor
//!    on the boundary, `u = ρ Σ_S a_S (c_S − f_v)/Σ_S a_S` with `c_S` the circumcenter
//!    (`rs-circumcenter`; the barycenter for faces with a boundary edge, as the reference),
//!    projected onto the plane orthogonal to the angle-weighted vertex normal (face normals
//!    oriented consistently around the fan, so non-orientable surfaces work). All `u` are
//!    computed from the positions at the start of the smoothing (Jacobi); they are then
//!    applied in vertex order, each guarded: if a face of the star would turn its normal by
//!    more than `π − 2` against its normal at the start of the smoothing (or become
//!    degenerate), the move is halved, up to 3 times, then refused (counted). Order
//!    [`RemeshOrder::FlipThenSmooth`] (the paper, default) or
//!    [`RemeshOrder::SmoothThenFlip`] (the reference).
//! 3. The result is revalidated as a [`Surface`] with positions, and returned
//!    with the [`VertexMap`] (new vertex → old vertices and weights), the face parents and the
//!    remapped pins ([`Remeshed`]); per-vertex data follow with [`VertexMap::transport`].
//!
//! The pass may renumber freely; only the vertex map and the face parents are contractual.
//! Single operations ([`EditMesh::split`], [`EditMesh::collapse`], [`EditMesh::flip`]) follow
//! the numbering of the oracle's `remesh-ops` (`SPEC.md` §10).

mod geom;
mod mesh;
mod transport;

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet, VecDeque};
use std::time::Instant;

pub use geom::{
    COS_FOLDOVER_LIMIT, FOLDOVER_LIMIT, angle_between, barycenter, circumcenter, corner_angle,
    folds, smoothing_displacement, tri_normal, violates_delaunay,
};
pub use mesh::{
    CollapseChecks, CollapseGuard, CollapseQuery, EditMesh, Fan, FlipChecks, FlipOutcome,
    RefusalCounts, RefusalKind, Refused, Remeshed,
};
pub use transport::VertexMap;

use crate::error::RsError;
use crate::geometry::v3::{add, dot, norm, scale, sub};
use crate::surface::Surface;

/// The order of the two parts of a round (`SPEC.md` §10 step 2, D16).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RemeshOrder {
    /// Delaunay flips, then smoothing (the paper; default).
    #[default]
    FlipThenSmooth,
    /// Smoothing, then flips (the reference, `ref: dynamic_remesher.cpp:92–96`).
    SmoothThenFlip,
}

/// Which parts of the pass run (all by default).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnabledOps {
    /// Split long edges.
    pub split: bool,
    /// Collapse short edges.
    pub collapse: bool,
    /// Delaunay flips.
    pub flip: bool,
    /// Tangential smoothing.
    pub smooth: bool,
}

impl Default for EnabledOps {
    fn default() -> Self {
        EnabledOps {
            split: true,
            collapse: true,
            flip: true,
            smooth: true,
        }
    }
}

/// Options of a remeshing pass (`SPEC.md` §10). `L₀` is passed to [`remesh`] separately.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RemeshOptions {
    /// Split edges longer than `split · L₀` (default 1.5).
    pub split: f64,
    /// Collapse edges shorter than `collapse · L₀` (default 0.5).
    pub collapse: f64,
    /// Number of flip/smooth rounds `N` (default 5).
    pub rounds: usize,
    /// Smoothing step `ρ` (default 0.5).
    pub rho: f64,
    /// At most `flip_cap_per_vertex · |V|` flips per round (default 100).
    pub flip_cap_per_vertex: usize,
    /// Flip then smooth (default) or the reference's smooth then flip.
    pub order: RemeshOrder,
    /// Which parts run.
    pub enabled_ops: EnabledOps,
}

impl Default for RemeshOptions {
    fn default() -> Self {
        RemeshOptions {
            split: 1.5,
            collapse: 0.5,
            rounds: 5,
            rho: 0.5,
            flip_cap_per_vertex: 100,
            order: RemeshOrder::FlipThenSmooth,
            enabled_ops: EnabledOps::default(),
        }
    }
}

impl RemeshOptions {
    /// Checks the ranges.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] naming the option: `split` finite and positive, `collapse`
    /// finite and non-negative, `split > 2 · collapse` (the halves of a split edge must not be
    /// collapsible, or a pass could undo its own splits), `rho ∈ (0, 1]`,
    /// `flip_cap_per_vertex ≥ 1`.
    pub fn validate(&self) -> Result<(), RsError> {
        let bad = |name, value: f64, requirement| {
            Err(RsError::InvalidOption {
                name,
                value: value.to_string(),
                requirement,
            })
        };
        if !(self.split.is_finite() && self.split > 0.0) {
            return bad("split", self.split, "must be finite and positive");
        }
        if !(self.collapse.is_finite() && self.collapse >= 0.0) {
            return bad("collapse", self.collapse, "must be finite and non-negative");
        }
        if self.split <= 2.0 * self.collapse {
            return bad(
                "split",
                self.split,
                "must exceed 2 x collapse (the halves of a split edge must not be collapsible)",
            );
        }
        if !(self.rho > 0.0 && self.rho <= 1.0) {
            return bad("rho", self.rho, "must lie in (0, 1]");
        }
        if self.flip_cap_per_vertex == 0 {
            return bad("flip_cap_per_vertex", 0.0, "must be at least 1");
        }
        Ok(())
    }
}

/// What a pass did (`SPEC.md` §9.1 `StepReport::remesh`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RemeshReport {
    /// Edges split.
    pub splits: u64,
    /// Edges collapsed.
    pub collapses: u64,
    /// Edges flipped (all rounds).
    pub flips: u64,
    /// Long edges whose split was refused, by reason (only a degenerate half).
    pub refused_splits: RefusalCounts,
    /// Short edges whose collapse was refused, by reason.
    pub refused_collapses: RefusalCounts,
    /// Non-Delaunay edges whose flip was refused, by reason (all rounds).
    pub refused_flips: RefusalCounts,
    /// Smoothing moves applied (all rounds), including halved ones.
    pub smoothing_moves: u64,
    /// Smoothing moves applied after at least one halving.
    pub smoothing_halved: u64,
    /// Smoothing moves refused after three halvings.
    pub smoothing_refused: u64,
    /// Vertices skipped by the smoothing because no vertex normal or no center was defined.
    pub smoothing_skipped: u64,
    /// Number of rounds whose flips stopped at the cap `flip_cap_per_vertex · |V|`.
    pub flip_cap_hits: u32,
    /// Vertices after the pass.
    pub n_vertices: usize,
    /// Faces after the pass.
    pub n_faces: usize,
    /// Wall time of the pass.
    pub seconds: f64,
}

/// One remeshing pass (`SPEC.md` §10; see the module documentation).
///
/// `pinned` vertices are never removed or moved (a collapse with one pinned end goes to the
/// pinned vertex; both pinned: refused); [`Remeshed::pinned`] gives their new ids in the same
/// order. `guard` is the optional extra collapse check (the CCD gap test of `SPEC.md` §10,
/// WP8; pass `None`).
///
/// # Errors
/// [`RsError::InvalidOption`] for bad options, a non-positive or non-finite `l0` or a pin
/// that is not a vertex; the errors of [`Surface::validate_positions`] for the input; the
/// errors of [`Surface::new`] / [`Surface::validate_positions`] if the result failed its
/// revalidation (an internal error, never silently repaired).
pub fn remesh(
    surface: &Surface,
    positions: &[[f64; 3]],
    pinned: &[u32],
    l0: f64,
    options: &RemeshOptions,
    mut guard: Option<&mut dyn CollapseGuard>,
) -> Result<(Remeshed, RemeshReport), RsError> {
    let t0 = Instant::now();
    options.validate()?;
    if !(l0.is_finite() && l0 > 0.0) {
        return Err(RsError::InvalidOption {
            name: "l0",
            value: l0.to_string(),
            requirement: "the target edge length must be finite and positive",
        });
    }
    let mut mesh = EditMesh::new(surface, positions, pinned)?;
    let mut rep = RemeshReport::default();

    if options.enabled_ops.split {
        let thr = options.split * l0;
        for &[a, b] in surface.edges() {
            if edge_len(&mesh, a, b) > thr {
                match mesh.split(a, b) {
                    Ok(_) => rep.splits += 1,
                    Err(r) => rep.refused_splits.add(r.kind),
                }
            }
        }
    }
    if options.enabled_ops.collapse {
        collapse_short(&mut mesh, options.collapse * l0, &mut guard, &mut rep);
    }
    for _ in 0..options.rounds {
        match options.order {
            RemeshOrder::FlipThenSmooth => {
                if options.enabled_ops.flip {
                    delaunay_flips(&mut mesh, options.flip_cap_per_vertex, &mut rep);
                }
                if options.enabled_ops.smooth {
                    smooth(&mut mesh, options.rho, &mut rep);
                }
            }
            RemeshOrder::SmoothThenFlip => {
                if options.enabled_ops.smooth {
                    smooth(&mut mesh, options.rho, &mut rep);
                }
                if options.enabled_ops.flip {
                    delaunay_flips(&mut mesh, options.flip_cap_per_vertex, &mut rep);
                }
            }
        }
    }
    let out = mesh.finish()?;
    rep.n_vertices = out.surface.n_vertices();
    rep.n_faces = out.surface.n_faces();
    rep.seconds = t0.elapsed().as_secs_f64();
    Ok((out, rep))
}

fn edge_len(mesh: &EditMesh, a: u32, b: u32) -> f64 {
    let p = mesh.positions();
    norm(sub(p[a as usize], p[b as usize]))
}

type Short = Reverse<(u64, u32, u32)>;

fn push_short(heap: &mut BinaryHeap<Short>, mesh: &EditMesh, a: u32, b: u32, thr: f64) {
    let l = edge_len(mesh, a, b);
    if l < thr {
        // lengths are non-negative, so their bit patterns order like the values
        heap.push(Reverse((l.to_bits(), a.min(b), a.max(b))));
    }
}

fn collapse_short<'g>(
    mesh: &mut EditMesh,
    thr: f64,
    guard: &mut Option<&mut (dyn CollapseGuard + 'g)>,
    rep: &mut RemeshReport,
) {
    let mut heap = BinaryHeap::new();
    for [a, b] in mesh.edges() {
        push_short(&mut heap, mesh, a, b, thr);
    }
    let checks = CollapseChecks { foldover: true };
    while let Some(Reverse((bits, a, b))) = heap.pop() {
        if !mesh.is_vertex_alive(a) || !mesh.is_vertex_alive(b) {
            continue;
        }
        if mesh.edge_faces(a, b).is_empty() {
            continue;
        }
        let l = edge_len(mesh, a, b);
        if l.to_bits() != bits || l >= thr {
            continue; // stale entry: the edge moved since it was queued
        }
        let g: Option<&mut (dyn CollapseGuard + 'g)> = match guard {
            Some(g) => Some(&mut **g),
            None => None,
        };
        match mesh.collapse(a, b, checks, g) {
            Ok(keep) => {
                rep.collapses += 1;
                for w in mesh.neighbors(keep) {
                    push_short(&mut heap, mesh, keep, w, thr);
                }
            }
            Err(r) => rep.refused_collapses.add(r.kind),
        }
    }
}

fn edge_key(a: u32, b: u32) -> u64 {
    ((a.min(b) as u64) << 32) | a.max(b) as u64
}

fn delaunay_flips(mesh: &mut EditMesh, cap_per_vertex: usize, rep: &mut RemeshReport) {
    let cap = (cap_per_vertex as u64).saturating_mul(mesh.n_vertices() as u64);
    let mut queue: VecDeque<[u32; 2]> = mesh.edges().into();
    let mut in_queue: HashSet<u64> = queue.iter().map(|&[a, b]| edge_key(a, b)).collect();
    let mut flips = 0u64;
    let checks = FlipChecks {
        delaunay: false,
        foldover: true,
    };
    while let Some([a, b]) = queue.pop_front() {
        in_queue.remove(&edge_key(a, b));
        let inc = mesh.edge_faces(a, b);
        if inc.len() != 2 {
            continue;
        }
        let quad = |f: u32| {
            let t = mesh.face(f);
            t.into_iter().find(|&v| v != a && v != b).unwrap_or(a)
        };
        let (c, d) = (quad(inc[0]), quad(inc[1]));
        let p = mesh.positions();
        let (xa, xb, xc, xd) = (p[a as usize], p[b as usize], p[c as usize], p[d as usize]);
        if !violates_delaunay(xa, xb, xc, xd) {
            continue;
        }
        if flips >= cap {
            rep.flip_cap_hits += 1;
            break;
        }
        match mesh.flip(a, b, checks) {
            Ok(_) => {
                flips += 1;
                for [x, y] in [[a, c], [c, b], [b, d], [d, a]] {
                    if in_queue.insert(edge_key(x, y)) {
                        queue.push_back([x.min(y), x.max(y)]);
                    }
                }
            }
            Err(r) => rep.refused_flips.add(r.kind),
        }
    }
    rep.flips += flips;
}

fn smooth(mesh: &mut EditMesh, rho: f64, rep: &mut RemeshReport) {
    let nf = mesh.n_face_slots();
    let nv = mesh.n_vertex_slots();
    // boundary vertices and faces with a boundary edge (an edge with one face)
    let mut boundary = vec![false; nv];
    let mut face_on_boundary = vec![false; nf];
    for ([a, b], n) in mesh.edges_with_counts() {
        if n == 1 {
            boundary[a as usize] = true;
            boundary[b as usize] = true;
            for f in mesh.edge_faces(a, b) {
                face_on_boundary[f as usize] = true;
            }
        }
    }
    let mut normal0 = vec![[0.0; 3]; nf];
    let mut center = vec![[0.0; 3]; nf];
    let mut area = vec![0.0; nf];
    for f in 0..nf as u32 {
        if !mesh.is_face_alive(f) {
            continue;
        }
        let p = mesh.corners(f);
        let n = tri_normal(p);
        normal0[f as usize] = n;
        area[f as usize] = 0.5 * norm(n);
        center[f as usize] = if face_on_boundary[f as usize] {
            barycenter(p)
        } else {
            circumcenter(p)
        };
    }
    let mut moves: Vec<(u32, [f64; 3])> = Vec::new();
    for v in 0..nv as u32 {
        if !mesh.is_vertex_alive(v) || mesh.is_pinned(v) || boundary[v as usize] {
            continue;
        }
        let x = mesh.positions()[v as usize];
        let Some(fan) = mesh.fan(v) else {
            rep.smoothing_skipped += 1;
            continue;
        };
        // angle-weighted normal with the faces oriented consistently around the fan
        let mut nrm = [0.0; 3];
        for (&f, &fl) in fan.faces.iter().zip(&fan.flip) {
            let t = mesh.face(f);
            let k = t.iter().position(|&w| w == v).unwrap_or(0);
            let p = mesh.positions();
            let ang = corner_angle(x, p[t[(k + 1) % 3] as usize], p[t[(k + 2) % 3] as usize]);
            let n = normal0[f as usize];
            let len = norm(n);
            let s = if fl { -ang / len } else { ang / len };
            nrm = add(nrm, scale(n, s));
        }
        let len = norm(nrm);
        let disp = smoothing_displacement(
            x,
            rho,
            fan.faces
                .iter()
                .map(|&f| (area[f as usize], center[f as usize])),
        );
        let (Some(u), true) = (disp, len > 0.0 && len.is_finite()) else {
            rep.smoothing_skipped += 1;
            continue;
        };
        let n = scale(nrm, 1.0 / len);
        let u = sub(u, scale(n, dot(n, u)));
        if u.iter().all(|c| c.is_finite()) {
            moves.push((v, u));
        } else {
            rep.smoothing_skipped += 1;
        }
    }
    for (v, u) in moves {
        let x = mesh.positions()[v as usize];
        let mut step = u;
        let mut applied = false;
        for halvings in 0..=3 {
            let y = add(x, step);
            let ok = mesh.vertex_faces(v).iter().all(|&f| {
                let t = mesh.face(f);
                let mut p = mesh.corners(f);
                for k in 0..3 {
                    if t[k] == v {
                        p[k] = y;
                    }
                }
                !folds(normal0[f as usize], tri_normal(p))
            });
            if ok {
                mesh.set_position(v, y);
                rep.smoothing_moves += 1;
                if halvings > 0 {
                    rep.smoothing_halved += 1;
                }
                applied = true;
                break;
            }
            step = scale(step, 0.5);
        }
        if !applied {
            rep.smoothing_refused += 1;
        }
    }
}
