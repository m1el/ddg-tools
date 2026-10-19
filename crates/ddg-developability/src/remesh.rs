//! Small-angle remeshing (SPEC §5.3, D21): flips of cap triangles and collapses of needle
//! triangles, with the legality rules that keep the mesh a valid manifold without degenerate or
//! folded triangles.
//!
//! A triangle is *bad* if one of its corner angles is below `small_angle`. With
//! `small_angle < π/3` it has one or two such angles [formula: dev-small-angle-classification]:
//!
//! * two small angles (a cap, the third angle is `> π − 2·small_angle`): flip the edge opposite
//!   the large angle [formula: dev-flip-manifold];
//! * one small angle (a needle): collapse the opposite (shortest) edge
//!   [formula: dev-link-condition].
//!
//! Additions to the SPEC text, all of which only refuse more operations: a flip is refused when a new
//! face's normal turns by 90 degrees or more against either old face of the quad (the quad is not
//! convex; this is the flip's "normal turn" rule), a pinned vertex is never moved and survives a
//! collapse (so its id is kept) and the report contains the new pinned ids.

use std::collections::HashSet;

use ddg_core::topology::{TriMesh, compact_vertices};

use crate::error::DevError;
use crate::optimize::DEFAULT_MIN_FACE_ANGLE;

/// Options of [`remesh`] (SPEC §5.3).
#[derive(Debug, Clone, PartialEq)]
pub struct RemeshOptions {
    /// A triangle with a corner angle below this (radians, `0 < · < π/3`) is bad. Default 1 degree.
    pub small_angle: f64,
    /// Whether boundary edges may be collapsed (to their midpoint). Default `true`.
    pub collapse_boundary_edges: bool,
    /// Maximum number of passes (`>= 1`). Default 10.
    pub max_passes: usize,
    /// For a flow: run the check every this many accepted iterations (`>= 1`). Default 1.
    pub check_interval: usize,
    /// Admissibility floor (D31): a collapse is refused if a face of the merged star would get a
    /// smallest corner angle below this (radians) and below the angle that face has now. `0`
    /// switches the rule off. Default [`DEFAULT_MIN_FACE_ANGLE`].
    pub min_face_angle: f64,
    /// Rebuild and validate the whole mesh after every operation (O(|F|) each; for tests).
    /// Default `false`.
    pub validate_each_operation: bool,
}

impl Default for RemeshOptions {
    fn default() -> Self {
        RemeshOptions {
            small_angle: std::f64::consts::PI / 180.0,
            collapse_boundary_edges: true,
            max_passes: 10,
            check_interval: 1,
            min_face_angle: DEFAULT_MIN_FACE_ANGLE,
            validate_each_operation: false,
        }
    }
}

impl RemeshOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`DevError::InvalidOption`] naming the first violation.
    pub fn validate(&self) -> Result<(), DevError> {
        let a = self.small_angle;
        if !(a.is_finite() && a > 0.0 && a < std::f64::consts::PI / 3.0) {
            return Err(DevError::InvalidOption { name: "small_angle", value: a, requirement: "in (0, pi/3)" });
        }
        if self.max_passes == 0 {
            return Err(DevError::InvalidOption { name: "max_passes", value: 0.0, requirement: ">= 1" });
        }
        if self.check_interval == 0 {
            return Err(DevError::InvalidOption { name: "check_interval", value: 0.0, requirement: ">= 1" });
        }
        let m = self.min_face_angle;
        if !(m.is_finite() && (0.0..a).contains(&m)) {
            return Err(DevError::InvalidOption { name: "min_face_angle", value: m, requirement: "in [0, small_angle)" });
        }
        Ok(())
    }
}

/// What a remeshing did (SPEC §5.3).
#[derive(Debug, Clone, PartialEq)]
pub struct RemeshReport {
    /// Edge flips performed.
    pub flips: usize,
    /// Edge collapses performed.
    pub collapses: usize,
    /// Flips refused by a legality rule (counted per attempt, over all passes).
    pub refused_flips: usize,
    /// Collapses refused by a legality rule (counted per attempt, over all passes).
    pub refused_collapses: usize,
    /// Bad triangles left at the end.
    pub unresolved_bad_triangles: usize,
    /// Number of passes run.
    pub passes: usize,
    /// Smallest corner angle (radians) before.
    pub min_angle_before: f64,
    /// Smallest corner angle (radians) after.
    pub min_angle_after: f64,
    /// Set by the flow driver: the remeshing was computed but not applied because its result
    /// would not be a valid start for the optimizer (the energy undefined or a face below the
    /// optimizer's angle floor); the reason. The counts above describe the rejected remeshing.
    pub rejected: Option<&'static str>,
    /// Energy before; filled in by the flow driver (this function does not know the energy).
    pub energy_before: Option<f64>,
    /// Energy after; filled in by the flow driver.
    pub energy_after: Option<f64>,
}

/// Result of [`remesh`].
#[derive(Debug, Clone, PartialEq)]
pub struct Remeshed {
    /// The new mesh (always valid; vertices compacted).
    pub mesh: TriMesh,
    /// Positions of the new mesh.
    pub positions: Vec<[f64; 3]>,
    /// For each new vertex, the id it had in the input (survivors keep their identity).
    pub vertex_to_old: Vec<u32>,
    /// For each new face, the input face it descends from. A collapse removes two faces and the
    /// others keep their identity; a flip rewrites the two faces of the quad in place: the new
    /// face `(a, q, p)` inherits the flagged cap triangle `f` and `(q, b, p)` inherits its
    /// neighbour `g` (deterministic: the faces keep their slots, and slots keep their order).
    pub face_to_old: Vec<u32>,
    /// The input's pinned vertices in the new numbering (ascending, deduplicated; two pinned
    /// vertices are never merged).
    pub pinned: Vec<u32>,
    /// Counts and angles.
    pub report: RemeshReport,
}

impl Remeshed {
    /// Whether any operation was applied.
    pub fn changed(&self) -> bool {
        self.report.flips + self.report.collapses > 0
    }
}

type P3 = [f64; 3];

fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: P3, b: P3) -> P3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: P3) -> f64 {
    dot(a, a).sqrt()
}

/// Corner angles of the triangle `(a, b, c)` at `a`, `b`, `c` (atan2 form, accurate when small).
fn angles(a: P3, b: P3, c: P3) -> [f64; 3] {
    let at = |p: P3, q: P3, r: P3| {
        let (u, v) = (sub(q, p), sub(r, p));
        norm(cross(u, v)).atan2(dot(u, v))
    };
    [at(a, b, c), at(b, c, a), at(c, a, b)]
}

/// Normal `(b - a) x (c - a)` and whether the triangle is non-degenerate (relative tolerance).
fn normal_and_ok(a: P3, b: P3, c: P3) -> (P3, bool) {
    let n = cross(sub(b, a), sub(c, a));
    let l2 = [norm(sub(b, a)), norm(sub(c, b)), norm(sub(a, c))].into_iter().fold(0.0, f64::max).powi(2);
    let nn = norm(n);
    (n, nn.is_finite() && nn > 1e-12 * l2 && nn > 0.0)
}

struct State {
    pos: Vec<P3>,
    faces: Vec<[u32; 3]>,
    alive: Vec<bool>,
    vfaces: Vec<Vec<u32>>,
    pinned: Vec<bool>,
    merged_into: Vec<u32>,
}

impl State {
    fn new(pos: Vec<P3>, faces: Vec<[u32; 3]>, pinned: Vec<bool>) -> State {
        let nv = pos.len();
        let mut vfaces = vec![Vec::new(); nv];
        for (f, tri) in faces.iter().enumerate() {
            for &v in tri {
                vfaces[v as usize].push(f as u32);
            }
        }
        State { alive: vec![true; faces.len()], pos, faces, vfaces, pinned, merged_into: (0..nv as u32).collect() }
    }
    fn tri_pos(&self, f: usize) -> [P3; 3] {
        self.faces[f].map(|v| self.pos[v as usize])
    }
    fn face_angles(&self, f: usize) -> [f64; 3] {
        let [a, b, c] = self.tri_pos(f);
        angles(a, b, c)
    }
    fn edge_faces(&self, a: u32, b: u32) -> Vec<u32> {
        self.vfaces[a as usize].iter().copied().filter(|&f| self.faces[f as usize].contains(&b)).collect()
    }
    fn neighbors(&self, v: u32) -> HashSet<u32> {
        let mut s = HashSet::new();
        for &f in &self.vfaces[v as usize] {
            for &w in &self.faces[f as usize] {
                if w != v {
                    s.insert(w);
                }
            }
        }
        s
    }
    fn is_boundary_vertex(&self, v: u32) -> bool {
        self.neighbors(v).into_iter().any(|w| self.edge_faces(v, w).len() == 1)
    }
    fn set_face(&mut self, f: usize, tri: [u32; 3]) {
        for v in self.faces[f] {
            self.vfaces[v as usize].retain(|&g| g as usize != f);
        }
        self.faces[f] = tri;
        for v in tri {
            self.vfaces[v as usize].push(f as u32);
        }
    }
    fn kill_face(&mut self, f: usize) {
        for v in self.faces[f] {
            self.vfaces[v as usize].retain(|&g| g as usize != f);
        }
        self.alive[f] = false;
    }
    /// Is the component of `v` larger than 4 vertices?
    fn component_exceeds_four(&self, v: u32) -> bool {
        let mut seen = HashSet::from([v]);
        let mut stack = vec![v];
        while let Some(u) = stack.pop() {
            for w in self.neighbors(u) {
                if seen.insert(w) {
                    if seen.len() > 4 {
                        return true;
                    }
                    stack.push(w);
                }
            }
        }
        false
    }
    fn live_faces(&self) -> Vec<[u32; 3]> {
        self.faces.iter().zip(&self.alive).filter(|(_, a)| **a).map(|(f, _)| *f).collect()
    }
    fn validate(&self) -> Result<(), DevError> {
        let faces = self.live_faces();
        let c = compact_vertices(self.pos.len(), &faces)?;
        TriMesh::new(c.new_to_old.len(), c.faces)?;
        Ok(())
    }
}

/// Why an operation was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Refused {
    Rule(&'static str),
}

/// Flips the edge `a -> b` of face `f` (opposite vertex `p`) if legal.
///
/// [formula: dev-flip-manifold]
fn try_flip(st: &mut State, f: usize, a: u32, b: u32, p: u32) -> Result<(), Refused> {
    let ef = st.edge_faces(a, b);
    if ef.len() != 2 {
        return Err(Refused::Rule("flip: boundary edge"));
    }
    let g = if ef[0] as usize == f { ef[1] } else { ef[0] } as usize;
    // the other face is (b, a, q) up to rotation
    let q = *st.faces[g].iter().find(|&&v| v != a && v != b).ok_or(Refused::Rule("flip: inconsistent faces"))?;
    if p == q || st.neighbors(p).contains(&q) {
        return Err(Refused::Rule("flip: opposite vertices adjacent"));
    }
    // (implied by the adjacency rule for interior vertices; kept as a safety net)
    for v in [a, b] {
        let k = st.vfaces[v as usize].len();
        let min = if st.is_boundary_vertex(v) { 1 } else { 3 };
        if k - 1 < min {
            return Err(Refused::Rule("flip: valence"));
        }
    }
    let (pa, pb, pp, pq) = (st.pos[a as usize], st.pos[b as usize], st.pos[p as usize], st.pos[q as usize]);
    let old_min = st.face_angles(f).into_iter().chain(st.face_angles(g)).fold(f64::INFINITY, f64::min);
    // new faces (a, q, p) and (q, b, p)
    let (n1, ok1) = normal_and_ok(pa, pq, pp);
    let (n2, ok2) = normal_and_ok(pq, pb, pp);
    if !(ok1 && ok2) {
        return Err(Refused::Rule("flip: degenerate new face"));
    }
    // normal turn against the old faces (checked only against non-degenerate old faces)
    let old_normals = [normal_and_ok(pa, pb, pp), normal_and_ok(pb, pa, pq)];
    for n in [n1, n2] {
        if old_normals.iter().any(|&(o, ok)| ok && dot(n, o) <= 0.0) {
            return Err(Refused::Rule("flip: normal turn"));
        }
    }
    let new_min = angles(pa, pq, pp).into_iter().chain(angles(pq, pb, pp)).fold(f64::INFINITY, f64::min);
    if new_min.partial_cmp(&old_min) != Some(std::cmp::Ordering::Greater) {
        return Err(Refused::Rule("flip: smallest angle does not increase"));
    }
    st.set_face(f, [a, q, p]);
    st.set_face(g, [q, b, p]);
    Ok(())
}

/// Collapses the edge `(u, w)` if legal.
///
/// [formula: dev-link-condition]
fn try_collapse(st: &mut State, u: u32, w: u32, allow_boundary: bool, floor: f64) -> Result<(), Refused> {
    let ef = st.edge_faces(u, w);
    if ef.is_empty() || ef.len() > 2 {
        return Err(Refused::Rule("collapse: not an edge"));
    }
    let (bu, bw) = (st.is_boundary_vertex(u), st.is_boundary_vertex(w));
    let boundary_edge = ef.len() == 1;
    if boundary_edge && !allow_boundary {
        return Err(Refused::Rule("collapse: boundary edges disabled"));
    }
    if !boundary_edge && bu && bw {
        return Err(Refused::Rule("collapse: interior edge between boundary vertices"));
    }
    if st.pinned[u as usize] && st.pinned[w as usize] {
        return Err(Refused::Rule("collapse: both endpoints pinned"));
    }
    // link condition: common neighbours are exactly the opposite vertices of the edge's faces
    let opposite: HashSet<u32> = ef
        .iter()
        .filter_map(|&f| st.faces[f as usize].iter().copied().find(|&v| v != u && v != w))
        .collect();
    let (nu, nw) = (st.neighbors(u), st.neighbors(w));
    let common: HashSet<u32> = nu.intersection(&nw).copied().collect();
    if common != opposite {
        return Err(Refused::Rule("collapse: link condition"));
    }
    if !st.component_exceeds_four(u) {
        return Err(Refused::Rule("collapse: component has at most 4 vertices"));
    }
    // survivor, position
    let (s, r) = if st.pinned[w as usize] || (!st.pinned[u as usize] && w < u) { (w, u) } else { (u, w) };
    let (ps, pr) = (st.pos[s as usize], st.pos[r as usize]);
    let x = if st.pinned[s as usize] {
        ps
    } else if bu != bw {
        if bu { st.pos[u as usize] } else { st.pos[w as usize] }
    } else {
        [0, 1, 2].map(|k| 0.5 * (ps[k] + pr[k]))
    };
    let _ = pr;
    // surviving faces of the merged star
    let mut seen = HashSet::new();
    for &f in st.vfaces[s as usize].iter().chain(&st.vfaces[r as usize]) {
        if !seen.insert(f) {
            continue;
        }
        let tri = st.faces[f as usize];
        if tri.contains(&s) && tri.contains(&r) {
            continue; // removed
        }
        let old = st.tri_pos(f as usize);
        let new_tri = tri.map(|v| if v == r { s } else { v });
        let new_pos = new_tri.map(|v| if v == s { x } else { st.pos[v as usize] });
        let (nn, ok) = normal_and_ok(new_pos[0], new_pos[1], new_pos[2]);
        if !ok {
            return Err(Refused::Rule("collapse: degenerate face in the merged star"));
        }
        if floor > 0.0 {
            let new_min = angles(new_pos[0], new_pos[1], new_pos[2]).into_iter().fold(f64::INFINITY, f64::min);
            let old_min = angles(old[0], old[1], old[2]).into_iter().fold(f64::INFINITY, f64::min);
            if new_min < floor && new_min < old_min {
                return Err(Refused::Rule("collapse: face below the angle floor"));
            }
        }
        let (on, old_ok) = normal_and_ok(old[0], old[1], old[2]);
        if old_ok && dot(nn, on) <= 0.0 {
            return Err(Refused::Rule("collapse: normal turn of the merged star"));
        }
    }
    // valence rules
    let removed = ef.len();
    for &v in &opposite {
        let k = st.vfaces[v as usize].len() - 1;
        let min = if st.is_boundary_vertex(v) { 1 } else { 3 };
        if k < min {
            return Err(Refused::Rule("collapse: valence"));
        }
    }
    let merged_faces = st.vfaces[s as usize].len() + st.vfaces[r as usize].len() - 2 * removed;
    let merged_boundary = bu || bw;
    if merged_faces < if merged_boundary { 1 } else { 3 } {
        return Err(Refused::Rule("collapse: valence of the merged vertex"));
    }
    // apply
    for f in ef {
        st.kill_face(f as usize);
    }
    let moved: Vec<u32> = st.vfaces[r as usize].clone();
    for f in moved {
        let tri = st.faces[f as usize].map(|v| if v == r { s } else { v });
        st.set_face(f as usize, tri);
    }
    st.pos[s as usize] = x;
    st.pinned[s as usize] |= st.pinned[r as usize];
    st.merged_into[r as usize] = s;
    Ok(())
}

fn min_angle(st: &State) -> f64 {
    (0..st.faces.len())
        .filter(|&f| st.alive[f])
        .flat_map(|f| st.face_angles(f))
        .fold(f64::INFINITY, f64::min)
}

/// Small-angle remeshing of a triangle mesh (SPEC §5.3).
///
/// A pass collects the bad triangles, sorts them by `(smallest angle, face id)` and handles each
/// that is still present and still bad; passes repeat until a pass changes nothing or
/// `max_passes`. Vertices are then compacted. Nothing is moved except the surviving vertex of a
/// collapse. Pinned vertices never move; the pinned vertex survives a collapse (both pinned: the
/// collapse is refused).
///
/// # Errors
/// Invalid options, wrong array lengths, non-finite positions, out-of-range pinned ids, and
/// [`DevError::Internal`] if the result fails validation (a bug).
pub fn remesh(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    opts: &RemeshOptions,
    pinned: &[u32],
) -> Result<Remeshed, DevError> {
    opts.validate()?;
    let nv = mesh.n_vertices();
    if positions.len() != nv {
        return Err(DevError::LengthMismatch { what: "positions", expected: nv, got: positions.len() });
    }
    if let Some(v) = positions.iter().position(|p| p.iter().any(|c| !c.is_finite())) {
        return Err(DevError::NonFinitePosition { vertex: v as u32 });
    }
    let mut is_pinned = vec![false; nv];
    for &p in pinned {
        if p as usize >= nv {
            return Err(DevError::VertexOutOfRange { what: "pinned", vertex: p, n_vertices: nv });
        }
        is_pinned[p as usize] = true;
    }
    let mut st = State::new(positions.to_vec(), mesh.faces().to_vec(), is_pinned);
    let mut rep = RemeshReport {
        flips: 0,
        collapses: 0,
        refused_flips: 0,
        refused_collapses: 0,
        unresolved_bad_triangles: 0,
        passes: 0,
        min_angle_before: min_angle(&st),
        min_angle_after: 0.0,
        rejected: None,
        energy_before: None,
        energy_after: None,
    };
    let small = opts.small_angle;

    for _ in 0..opts.max_passes {
        rep.passes += 1;
        let mut bad: Vec<(f64, usize)> = (0..st.faces.len())
            .filter(|&f| st.alive[f])
            .filter_map(|f| {
                let m = st.face_angles(f).into_iter().fold(f64::INFINITY, f64::min);
                (m < small).then_some((m, f))
            })
            .collect();
        bad.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        let mut changed = false;
        for (_, f) in bad {
            if !st.alive[f] {
                continue;
            }
            let ang = st.face_angles(f);
            let n_small = ang.iter().filter(|&&a| a < small).count();
            if n_small == 0 {
                continue;
            }
            let tri = st.faces[f];
            if n_small >= 2 {
                // cap: flip the edge opposite the largest angle
                let k = (0..3).max_by(|&i, &j| ang[i].total_cmp(&ang[j])).unwrap_or(0);
                let (p, a, b) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
                match try_flip(&mut st, f, a, b, p) {
                    Ok(()) => {
                        rep.flips += 1;
                        changed = true;
                    }
                    Err(Refused::Rule(_)) => rep.refused_flips += 1,
                }
            } else {
                let k = (0..3).find(|&i| ang[i] < small).unwrap_or(0);
                let (u, w) = (tri[(k + 1) % 3], tri[(k + 2) % 3]);
                match try_collapse(&mut st, u, w, opts.collapse_boundary_edges, opts.min_face_angle) {
                    Ok(()) => {
                        rep.collapses += 1;
                        changed = true;
                    }
                    Err(Refused::Rule(_)) => rep.refused_collapses += 1,
                }
            }
            if changed && opts.validate_each_operation {
                st.validate().map_err(|e| DevError::Internal(format!("remesh produced an invalid mesh: {e}")))?;
            }
        }
        if !changed {
            break;
        }
    }

    rep.unresolved_bad_triangles = (0..st.faces.len())
        .filter(|&f| st.alive[f] && st.face_angles(f).into_iter().fold(f64::INFINITY, f64::min) < small)
        .count();
    rep.min_angle_after = min_angle(&st);

    let faces = st.live_faces();
    let face_to_old: Vec<u32> = (0..st.faces.len() as u32).filter(|&f| st.alive[f as usize]).collect();
    let c = compact_vertices(nv, &faces)?;
    let new_mesh = TriMesh::new(c.new_to_old.len(), c.faces)
        .map_err(|e| DevError::Internal(format!("remesh produced an invalid mesh: {e}")))?;
    let new_pos: Vec<P3> = c.new_to_old.iter().map(|&o| st.pos[o as usize]).collect();
    let resolve = |mut v: u32| {
        while st.merged_into[v as usize] != v {
            v = st.merged_into[v as usize];
        }
        v
    };
    let mut new_pinned: Vec<u32> = pinned
        .iter()
        .map(|&p| c.old_to_new[resolve(p) as usize])
        .collect();
    new_pinned.sort_unstable();
    new_pinned.dedup();
    Ok(Remeshed { mesh: new_mesh, positions: new_pos, vertex_to_old: c.new_to_old, face_to_old, pinned: new_pinned, report: rep })
}

#[cfg(test)]
mod tests {
    //! One test per legality rule (SPEC §10.6), calling the primitive operations directly.
    use super::*;
    use ddg_core::generators::{SplitMix64, grid_disk, perturb};
    use std::collections::BTreeSet;

    fn state(pos: Vec<P3>, faces: Vec<[u32; 3]>) -> State {
        let n = pos.len();
        State::new(pos, faces, vec![false; n])
    }

    fn from_gen(g: &ddg_core::generators::GeneratedMesh) -> State {
        state(g.positions.clone(), g.faces.clone())
    }

    /// Flip of the directed edge `a -> b` (the face that contains it supplies `p`).
    fn flip(st: &mut State, a: u32, b: u32) -> Result<(), Refused> {
        let directed = |st: &State, a: u32, b: u32| {
            st.vfaces[a as usize].iter().map(|&f| f as usize).find(|&f| {
                let t = st.faces[f];
                (0..3).any(|k| t[k] == a && t[(k + 1) % 3] == b)
            })
        };
        let (a, b, f) = match directed(st, a, b) {
            Some(f) => (a, b, f),
            None => (b, a, directed(st, b, a).expect("directed edge")),
        };
        let p = *st.faces[f].iter().find(|&&v| v != a && v != b).unwrap();
        try_flip(st, f, a, b, p)
    }

    fn reason(r: Result<(), Refused>) -> &'static str {
        match r {
            Ok(()) => "ok",
            Err(Refused::Rule(s)) => s,
        }
    }

    #[test]
    fn cap_flip_is_legal_and_valid() {
        let mut st = state(
            vec![[0., 0., 0.], [10., 0., 0.], [5., 0.05, 0.], [5., -3., 0.]],
            vec![[0, 1, 2], [1, 0, 3]],
        );
        let before = min_angle(&st);
        assert_eq!(reason(flip(&mut st, 0, 1)), "ok");
        st.validate().unwrap();
        assert!(min_angle(&st) > before);
        assert_eq!(st.live_faces(), vec![[0, 3, 2], [3, 1, 2]]);
    }

    #[test]
    fn flip_refusals() {
        // boundary edge: the cap's long edge (0, 1) is on the boundary
        let mut st = state(
            vec![[0., 0., 0.], [10., 0., 0.], [5., 0.05, 0.], [7., 3., 0.], [3., 3., 0.]],
            vec![[0, 1, 2], [2, 1, 3], [0, 2, 4]],
        );
        assert_eq!(reason(flip(&mut st, 0, 1)), "flip: boundary edge");
        // an interior vertex with three faces: the opposite vertices of its edges are adjacent
        let s3 = 0.866_025_403_784_438_6;
        let mut st = state(
            vec![[0., 0., 0.], [1., 0., 0.], [-0.5, s3, 0.], [-0.5, -s3, 0.]],
            vec![[0, 1, 2], [0, 2, 3], [0, 3, 1]],
        );
        assert_eq!(reason(flip(&mut st, 0, 1)), "flip: opposite vertices adjacent");
        // degenerate new face: a, q, p collinear
        let mut st = state(
            vec![[0., 0., 0.], [2., 0., 0.], [1., 1., 0.], [-1., -1., 0.]],
            vec![[0, 1, 2], [1, 0, 3]],
        );
        assert_eq!(reason(flip(&mut st, 0, 1)), "flip: degenerate new face");
        // normal turn: reflex quad
        let mut st = state(
            vec![[-3., 0., 0.], [0.5, 0., 0.], [1., 1., 0.], [1., -1., 0.]],
            vec![[0, 1, 2], [1, 0, 3]],
        );
        assert_eq!(reason(flip(&mut st, 0, 1)), "flip: normal turn");
        // a well shaped quad: the flip would make the smallest angle smaller
        let mut st = state(
            vec![[0., 0., 0.], [1., 0., 0.], [0.5, 0.9, 0.], [0.5, -0.9, 0.]],
            vec![[0, 1, 2], [1, 0, 3]],
        );
        assert_eq!(reason(flip(&mut st, 0, 1)), "flip: smallest angle does not increase");
    }

    fn collapse(st: &mut State, u: u32, w: u32, allow: bool) -> &'static str {
        reason(try_collapse(st, u, w, allow, 0.0))
    }

    #[test]
    fn collapse_refusals_and_success() {
        let g = grid_disk(4, 4).unwrap();
        let m = TriMesh::new(g.positions.len(), g.faces.clone()).unwrap();
        // interior edge between interior vertices: legal, valid result
        let (u, w) = (6u32, 7u32); // (1,1) and (1,2) of a 5x5 vertex grid
        assert!(m.find_edge(u as usize, w as usize).is_some());
        let mut st = from_gen(&g);
        assert_eq!(collapse(&mut st, u, w, true), "ok");
        st.validate().unwrap();
        assert_eq!(st.live_faces().len(), g.faces.len() - 2);
        assert_eq!(st.pos[6], [0.5 * (0.25 + 0.25), 0.5 * (0.25 + 0.5), 0.0]);
        // boundary edges disabled
        let mut st = from_gen(&g);
        assert_eq!(collapse(&mut st, 0, 1, false), "collapse: boundary edges disabled");
        // boundary edge allowed: survivor stays put when one end is interior? both boundary: midpoint
        let mut st = from_gen(&g);
        assert_eq!(collapse(&mut st, 0, 1, true), "ok");
        st.validate().unwrap();
        // both pinned
        let mut st = from_gen(&g);
        st.pinned[6] = true;
        st.pinned[7] = true;
        assert_eq!(collapse(&mut st, 6, 7, true), "collapse: both endpoints pinned");
        // one pinned: the pinned vertex stays, at its position, and survives even with the higher id
        let mut st = from_gen(&g);
        st.pinned[7] = true;
        let p7 = st.pos[7];
        assert_eq!(collapse(&mut st, 6, 7, true), "ok");
        assert_eq!(st.pos[7], p7);
        assert_eq!(st.merged_into[6], 7);
        assert!(st.pinned[7]);
        // pinch: the diagonal of a quad joins two boundary vertices
        let g1 = grid_disk(1, 1).unwrap();
        let mut st = from_gen(&g1);
        let m1 = TriMesh::new(4, g1.faces.clone()).unwrap();
        let diag = m1.edges().iter().find(|e| m1.find_edge(e[0] as usize, e[1] as usize).is_some_and(|ed| !m1.is_boundary_edge(ed))).unwrap();
        assert_eq!(collapse(&mut st, diag[0], diag[1], true), "collapse: interior edge between boundary vertices");
        // link condition: a stacked tetrahedron (vertex 4 inside face (0,2,1)); edge (0,1) has the
        // extra common neighbour 2
        let mut st = state(
            vec![[0., 0., 0.], [1., 0., 0.], [0.5, 1., 0.], [0.5, 0.4, 1.], [0.5, 0.4, -0.2]],
            vec![[0, 2, 4], [2, 1, 4], [1, 0, 4], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
        );
        st.validate().unwrap();
        assert_eq!(collapse(&mut st, 0, 1, true), "collapse: link condition");
        // a tetrahedron passes the link condition but is too small
        let mut st = state(
            vec![[0., 0., 0.], [1., 0., 0.], [0.5, 1., 0.], [0.5, 0.4, 1.]],
            vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
        );
        assert_eq!(collapse(&mut st, 0, 1, true), "collapse: component has at most 4 vertices");
    }

    #[test]
    fn collapse_geometric_refusals_are_hit() {
        // sweep all edges of jittered planar and bumpy grids and collect the reasons
        let mut reasons = BTreeSet::new();
        let mut rng = SplitMix64::new(5);
        for seed in 0..40u64 {
            let mut g = grid_disk(6, 6).unwrap();
            perturb(&mut g.positions, [0.09 + 0.002 * rng.next_f64(), 0.09, 0.02], seed);
            let m = TriMesh::new(g.positions.len(), g.faces.clone()).unwrap();
            for e in m.edges() {
                let mut st = from_gen(&g);
                reasons.insert(collapse(&mut st, e[0], e[1], true));
                let mut st = from_gen(&g);
                reasons.insert(reason(flip(&mut st, e[0], e[1])));
                st.validate().ok();
            }
        }
        eprintln!("{reasons:?}");
        assert!(reasons.contains("collapse: normal turn of the merged star"), "{reasons:?}");
        assert!(reasons.contains("flip: normal turn"), "{reasons:?}");
    }
}
