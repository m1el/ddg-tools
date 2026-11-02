//! The mutable mesh of the remesher and its single operations: split, collapse, flip
//! (`SPEC.md` §10, "single-operation conventions").
//!
//! # Representation
//!
//! Faces are vertex triples in slots that are never reused (a collapse marks two slots dead, a
//! split appends); every vertex keeps the list of its live faces. A half-edge is a face slot
//! and a side `k` (from corner `k` to corner `k + 1`); its twin is the side of the other face
//! on the same edge, found through the vertex–face lists. No orientation is assumed: the two
//! half-edges of an interior edge may run the same way (non-orientable components, `D1`), and
//! the operations that care (the flip, the vertex normal of the smoothing) reorder locally.
//! Dead vertices and faces are dropped by [`EditMesh::finish`], which keeps the relative order
//! of the survivors, so a single operation numbers exactly as the oracle's `remesh-ops`.
//!
//! # Legality
//!
//! Every operation either succeeds or returns a [`Refused`] naming the reason and the
//! elements; nothing panics and nothing is changed by a refused operation. After every
//! collapse and flip the stars of the touched vertices are checked (one fan, at most two faces
//! per edge, no repeated face, no isolated vertex) and the operation undone if they fail
//! ([`RefusalKind::NonManifoldResult`]).

use std::fmt;

use crate::error::RsError;
use crate::geometry::v3::{add, scale};
use crate::surface::{NONE, Surface};

use super::geom::{folds, tri_normal, violates_delaunay};
use super::transport::VertexMap;

/// Why an operation was refused (`SPEC.md` §10). The oracle's error kinds map as follows:
/// `NoSuchEdge`, `BothPinned`, `LinkCondition` (also an interior edge between two boundary
/// vertices), `Foldover`, `NotDelaunay`, `BoundaryEdge`, `FlipNonManifold` = [`ExistingEdge`],
/// `NonManifoldResult` = [`RepeatedFace`] or [`NonManifoldResult`]. The oracle's
/// `InconsistentOrientation` is not a refusal here: the flip reorders locally
/// ([`FlipOutcome::reordered`]).
///
/// [`ExistingEdge`]: RefusalKind::ExistingEdge
/// [`RepeatedFace`]: RefusalKind::RepeatedFace
/// [`NonManifoldResult`]: RefusalKind::NonManifoldResult
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum RefusalKind {
    /// The two vertices share no edge (or a vertex does not exist).
    NoSuchEdge,
    /// Collapse: both endpoints are pinned.
    BothPinned,
    /// Collapse: the endpoints have other common neighbours than the edge's opposite
    /// vertices, or the edge is interior but both endpoints are on the boundary.
    LinkCondition,
    /// Collapse, flip, smoothing: a face would turn its normal by more than `π − 2` or become
    /// degenerate.
    Foldover,
    /// Collapse: two faces would have the same three vertices.
    RepeatedFace,
    /// Collapse, flip: the result would not be a manifold (a star is not one fan, an edge has
    /// three faces, a vertex loses all its faces).
    NonManifoldResult,
    /// Flip: the edge has only one face.
    BoundaryEdge,
    /// Flip with the Delaunay check: the opposite angles do not exceed `π`.
    NotDelaunay,
    /// Flip: the two opposite vertices already share an edge.
    ExistingEdge,
    /// Collapse: the caller's guard (the CCD check of `SPEC.md` §10, WP8) said no.
    Guard,
    /// Split, collapse, flip: a face would have zero (or non-finite) area. Always checked, so that
    /// the result revalidates; with the foldover check on, this is reported as
    /// [`Foldover`](RefusalKind::Foldover) (as the oracle does).
    DegenerateFace,
}

impl RefusalKind {
    /// All kinds, in the order of [`RefusalCounts`].
    pub const ALL: [RefusalKind; 11] = [
        RefusalKind::NoSuchEdge,
        RefusalKind::BothPinned,
        RefusalKind::LinkCondition,
        RefusalKind::Foldover,
        RefusalKind::RepeatedFace,
        RefusalKind::NonManifoldResult,
        RefusalKind::BoundaryEdge,
        RefusalKind::NotDelaunay,
        RefusalKind::ExistingEdge,
        RefusalKind::Guard,
        RefusalKind::DegenerateFace,
    ];

    /// The kind's position in [`RefusalKind::ALL`].
    pub fn index(self) -> usize {
        RefusalKind::ALL
            .iter()
            .position(|&k| k == self)
            .unwrap_or(0)
    }

    /// The name, as the oracle prints its error kinds where they coincide.
    pub fn name(self) -> &'static str {
        match self {
            RefusalKind::NoSuchEdge => "NoSuchEdge",
            RefusalKind::BothPinned => "BothPinned",
            RefusalKind::LinkCondition => "LinkCondition",
            RefusalKind::Foldover => "Foldover",
            RefusalKind::RepeatedFace => "RepeatedFace",
            RefusalKind::NonManifoldResult => "NonManifoldResult",
            RefusalKind::BoundaryEdge => "BoundaryEdge",
            RefusalKind::NotDelaunay => "NotDelaunay",
            RefusalKind::ExistingEdge => "ExistingEdge",
            RefusalKind::Guard => "Guard",
            RefusalKind::DegenerateFace => "DegenerateFace",
        }
    }
}

/// Counts of refused operations by [`RefusalKind`] (indexed by [`RefusalKind::index`]).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RefusalCounts(pub [u64; 11]);

impl RefusalCounts {
    /// Adds one refusal.
    pub fn add(&mut self, kind: RefusalKind) {
        self.0[kind.index()] += 1;
    }
    /// The count of one kind.
    pub fn get(&self, kind: RefusalKind) -> u64 {
        self.0[kind.index()]
    }
    /// The total over all kinds.
    pub fn total(&self) -> u64 {
        self.0.iter().sum()
    }
}

/// A refused operation: the reason and a message naming the elements (current ids of the
/// [`EditMesh`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Refused {
    /// The reason.
    pub kind: RefusalKind,
    /// What is wrong, naming vertices and faces.
    pub message: String,
}

impl fmt::Display for Refused {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.kind.name(), self.message)
    }
}

fn refuse<T>(kind: RefusalKind, message: String) -> Result<T, Refused> {
    Err(Refused { kind, message })
}

/// Optional checks of a single collapse (the link condition, both-pinned, repeated faces and
/// manifoldness are always checked).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CollapseChecks {
    /// Refuse if a surviving face turns its normal by more than `π − 2` or becomes degenerate.
    pub foldover: bool,
}

/// Optional checks of a single flip (interior edge, no existing edge `c d` and manifoldness are
/// always checked).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct FlipChecks {
    /// Refuse unless the opposite angles sum to more than `π` (`rs-delaunay-test`).
    pub delaunay: bool,
    /// Refuse if a new face turns its normal by more than `π − 2` (against the face it
    /// replaces) or is degenerate.
    pub foldover: bool,
}

/// The result of a successful flip.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FlipOutcome {
    /// The new edge `c d`.
    pub new_edge: [u32; 2],
    /// The two faces traversed the edge in the same direction (non-orientable surface), and
    /// the second was reordered before the flip (`SPEC.md` §10; the oracle refuses such a
    /// flip with `InconsistentOrientation`).
    pub reordered: bool,
}

/// What a collapse guard is asked (the CCD hook of `SPEC.md` §10, filled by WP8): may the
/// edge `edge` collapse into `survivor` at `position`, given the faces of the survivor's new
/// star? Positions are those before the collapse (`mesh.positions()`).
pub struct CollapseQuery<'a> {
    /// The mesh before the collapse.
    pub mesh: &'a EditMesh,
    /// The edge, as passed to [`EditMesh::collapse`].
    pub edge: [u32; 2],
    /// The surviving vertex.
    pub survivor: u32,
    /// Its new position.
    pub position: [f64; 3],
    /// The faces of the survivor's star after the collapse (current vertex ids, the removed
    /// vertex already renamed to the survivor).
    pub new_star: &'a [[u32; 3]],
}

/// A check a collapse must pass in addition to the built-in ones (`SPEC.md` §10: the CCD gap
/// test, only when CCD is on). `None` everywhere until WP8.
pub trait CollapseGuard {
    /// `false` refuses the collapse ([`RefusalKind::Guard`]).
    fn allow(&mut self, query: &CollapseQuery<'_>) -> bool;
}

/// The output of a remeshing: a revalidated surface, positions, and the maps back to the
/// input (`SPEC.md` §10).
#[derive(Debug, Clone)]
pub struct Remeshed {
    /// The new surface (validated by [`Surface::new`]).
    pub surface: Surface,
    /// One position per new vertex (validated by [`Surface::validate_positions`]).
    pub positions: Vec<[f64; 3]>,
    /// New vertex → old vertices and weights (for per-vertex data, [`VertexMap::transport`]).
    pub vertex_map: VertexMap,
    /// New face → the old face it was cut from (or is).
    pub face_parents: Vec<u32>,
    /// The new index of every input pin, in the caller's order (pins are never removed).
    pub pinned: Vec<u32>,
}

/// The mutable mesh of the remesher (see the module documentation).
#[derive(Debug, Clone)]
pub struct EditMesh {
    pos: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    face_alive: Vec<bool>,
    parent: Vec<u32>,
    vfaces: Vec<Vec<u32>>,
    v_alive: Vec<bool>,
    pinned: Vec<bool>,
    /// The input pin a vertex carries (old vertex id), or `NONE`.
    pin_origin: Vec<u32>,
    sources: Vec<Vec<(u32, f64)>>,
    pins_in: Vec<u32>,
    n_vertices_in: usize,
}

/// `(p, q, c)`: the side of `f` on the edge `{a, b}` in the face's cyclic order, and the
/// opposite vertex.
fn orient(f: [u32; 3], a: u32, b: u32) -> Option<(u32, u32, u32)> {
    for k in 0..3 {
        let (p, q, c) = (f[k], f[(k + 1) % 3], f[(k + 2) % 3]);
        if (p == a && q == b) || (p == b && q == a) {
            return Some((p, q, c));
        }
    }
    None
}

/// Does `f` run from `a` to `b` along one of its sides?
fn traverses(f: [u32; 3], a: u32, b: u32) -> bool {
    (0..3).any(|k| f[k] == a && f[(k + 1) % 3] == b)
}

fn sorted3(f: [u32; 3]) -> [u32; 3] {
    let mut s = f;
    s.sort_unstable();
    s
}

/// The faces around a vertex in rotational order, with the relative orientation of each face
/// (`flip[i]`: face `i` is traversed against the first face's orientation).
#[derive(Debug, Clone, PartialEq)]
pub struct Fan {
    /// Faces in rotational order (a boundary fan starts and ends at a face with a boundary
    /// edge at the vertex).
    pub faces: Vec<u32>,
    /// Per face: its vertex order must be reversed to agree with the first face's.
    pub flip: Vec<bool>,
    /// The fan is a half-disk.
    pub boundary: bool,
}

impl EditMesh {
    /// Builds the mutable mesh from a validated surface, positions and the pinned vertices.
    ///
    /// # Errors
    /// The errors of [`Surface::validate_positions`]; [`RsError::InvalidOption`] (`pinned`) for
    /// a pin that is not a vertex.
    pub fn new(surface: &Surface, positions: &[[f64; 3]], pinned: &[u32]) -> Result<Self, RsError> {
        surface.validate_positions(positions)?;
        let nv = surface.n_vertices();
        let mut pin = vec![false; nv];
        let mut pin_origin = vec![NONE; nv];
        for &p in pinned {
            if p as usize >= nv {
                return Err(RsError::InvalidOption {
                    name: "pinned",
                    value: format!("{p} (the surface has {nv} vertices)"),
                    requirement: "every pinned vertex must be a vertex of the surface",
                });
            }
            pin[p as usize] = true;
            pin_origin[p as usize] = p;
        }
        let faces = surface.faces().to_vec();
        let mut vfaces = vec![Vec::new(); nv];
        for (f, t) in faces.iter().enumerate() {
            for &v in t {
                vfaces[v as usize].push(f as u32);
            }
        }
        Ok(EditMesh {
            pos: positions.to_vec(),
            face_alive: vec![true; faces.len()],
            parent: (0..faces.len() as u32).collect(),
            faces,
            vfaces,
            v_alive: vec![true; nv],
            pinned: pin,
            pin_origin,
            sources: (0..nv as u32).map(|v| vec![(v, 1.0)]).collect(),
            pins_in: pinned.to_vec(),
            n_vertices_in: nv,
        })
    }

    // ---------------------------------------------------------------- queries

    /// Positions of all vertex slots (dead vertices keep their last position).
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.pos
    }
    /// Number of vertex slots (live and dead).
    pub fn n_vertex_slots(&self) -> usize {
        self.pos.len()
    }
    /// Number of face slots (live and dead).
    pub fn n_face_slots(&self) -> usize {
        self.faces.len()
    }
    /// Number of live vertices.
    pub fn n_vertices(&self) -> usize {
        self.v_alive.iter().filter(|&&a| a).count()
    }
    /// Number of live faces.
    pub fn n_faces(&self) -> usize {
        self.face_alive.iter().filter(|&&a| a).count()
    }
    /// The vertices of face slot `f`.
    pub fn face(&self, f: u32) -> [u32; 3] {
        self.faces[f as usize]
    }
    /// Is face slot `f` live?
    pub fn is_face_alive(&self, f: u32) -> bool {
        self.face_alive[f as usize]
    }
    /// Is vertex slot `v` live?
    pub fn is_vertex_alive(&self, v: u32) -> bool {
        self.v_alive[v as usize]
    }
    /// Is `v` pinned?
    pub fn is_pinned(&self, v: u32) -> bool {
        self.pinned[v as usize]
    }
    /// The live faces of `v` (unordered).
    pub fn vertex_faces(&self, v: u32) -> &[u32] {
        &self.vfaces[v as usize]
    }
    /// The corner positions of face slot `f`.
    pub fn corners(&self, f: u32) -> [[f64; 3]; 3] {
        let t = self.faces[f as usize];
        [
            self.pos[t[0] as usize],
            self.pos[t[1] as usize],
            self.pos[t[2] as usize],
        ]
    }
    pub(crate) fn set_position(&mut self, v: u32, x: [f64; 3]) {
        self.pos[v as usize] = x;
    }

    fn valid_vertex(&self, v: u32) -> bool {
        (v as usize) < self.pos.len() && self.v_alive[v as usize]
    }

    /// The live faces on the edge `{a, b}`, ascending.
    pub fn edge_faces(&self, a: u32, b: u32) -> Vec<u32> {
        if !self.valid_vertex(a) || !self.valid_vertex(b) || a == b {
            return Vec::new();
        }
        let mut r: Vec<u32> = self.vfaces[a as usize]
            .iter()
            .copied()
            .filter(|&f| self.faces[f as usize].contains(&b))
            .collect();
        r.sort_unstable();
        r
    }

    /// All live edges as `[min, max]`, ascending.
    pub fn edges(&self) -> Vec<[u32; 2]> {
        let mut e = Vec::with_capacity(self.faces.len() * 3 / 2 + 1);
        for (f, t) in self.faces.iter().enumerate() {
            if self.face_alive[f] {
                for k in 0..3 {
                    let (a, b) = (t[k], t[(k + 1) % 3]);
                    e.push([a.min(b), a.max(b)]);
                }
            }
        }
        e.sort_unstable();
        e.dedup();
        e
    }

    /// The neighbours of `v`, ascending.
    pub fn neighbors(&self, v: u32) -> Vec<u32> {
        let mut r = Vec::with_capacity(8);
        for &f in &self.vfaces[v as usize] {
            for &w in &self.faces[f as usize] {
                if w != v {
                    r.push(w);
                }
            }
        }
        r.sort_unstable();
        r.dedup();
        r
    }

    /// Is `v` on the boundary (some edge at `v` has one face)?
    pub fn is_boundary_vertex(&self, v: u32) -> bool {
        let mut nb: Vec<u32> = Vec::with_capacity(12);
        for &f in &self.vfaces[v as usize] {
            for &w in &self.faces[f as usize] {
                if w != v {
                    nb.push(w);
                }
            }
        }
        nb.sort_unstable();
        let mut i = 0;
        while i < nb.len() {
            let mut j = i;
            while j < nb.len() && nb[j] == nb[i] {
                j += 1;
            }
            if j - i == 1 {
                return true;
            }
            i = j;
        }
        false
    }

    /// The ordered fan of `v` with the relative orientations of its faces, or `None` if the
    /// faces around `v` are not one fan (or an edge at `v` has more than two faces).
    pub fn fan(&self, v: u32) -> Option<Fan> {
        let fs = &self.vfaces[v as usize];
        if fs.is_empty() {
            return None;
        }
        // wedges: the two other vertices of each face
        let wedge = |f: u32| -> [u32; 2] {
            let t = self.faces[f as usize];
            let k = t.iter().position(|&x| x == v).unwrap_or(0);
            [t[(k + 1) % 3], t[(k + 2) % 3]]
        };
        // neighbour -> faces
        let mut nbr: Vec<(u32, u32)> = Vec::with_capacity(2 * fs.len());
        for &f in fs {
            for w in wedge(f) {
                nbr.push((w, f));
            }
        }
        nbr.sort_unstable();
        let faces_of = |w: u32| -> &[(u32, u32)] {
            let lo = nbr.partition_point(|&(x, _)| x < w);
            let hi = nbr.partition_point(|&(x, _)| x <= w);
            &nbr[lo..hi]
        };
        let mut ends: Vec<(u32, u32)> = Vec::new(); // (face, neighbour) with a boundary edge
        let mut i = 0;
        while i < nbr.len() {
            let mut j = i;
            while j < nbr.len() && nbr[j].0 == nbr[i].0 {
                j += 1;
            }
            match j - i {
                1 => ends.push((nbr[i].1, nbr[i].0)),
                2 => {}
                _ => return None,
            }
            i = j;
        }
        let boundary = match ends.len() {
            0 => false,
            2 => true,
            _ => return None,
        };
        let (start, mut entry) = if boundary {
            *ends.iter().min().unwrap_or(&ends[0])
        } else {
            let f0 = *fs.iter().min().unwrap_or(&fs[0]);
            (f0, wedge(f0)[0])
        };
        let mut out = Fan {
            faces: vec![start],
            flip: vec![false],
            boundary,
        };
        let mut cur = start;
        let mut cur_flip = false;
        loop {
            let w = wedge(cur);
            let exit = if w[0] == entry { w[1] } else { w[0] };
            let next = faces_of(exit).iter().map(|&(_, f)| f).find(|&f| f != cur);
            let Some(next) = next else { break };
            if next == start {
                break;
            }
            if out.faces.len() > fs.len() {
                return None;
            }
            let same = traverses(self.faces[cur as usize], v, exit)
                == traverses(self.faces[next as usize], v, exit);
            cur_flip ^= same;
            out.faces.push(next);
            out.flip.push(cur_flip);
            cur = next;
            entry = exit;
        }
        if out.faces.len() != fs.len() {
            return None;
        }
        Some(out)
    }

    /// The local manifold check of a live vertex: one fan, at most two faces per edge, no two
    /// faces with the same vertices.
    fn vertex_ok(&self, v: u32) -> bool {
        if !self.v_alive[v as usize] {
            return true;
        }
        if self.fan(v).is_none() {
            return false;
        }
        let fs = &self.vfaces[v as usize];
        for i in 0..fs.len() {
            let si = sorted3(self.faces[fs[i] as usize]);
            for j in i + 1..fs.len() {
                if si == sorted3(self.faces[fs[j] as usize]) {
                    return false;
                }
            }
        }
        true
    }

    fn remove_vface(&mut self, v: u32, f: u32) {
        let l = &mut self.vfaces[v as usize];
        if let Some(i) = l.iter().position(|&x| x == f) {
            l.swap_remove(i);
        }
    }

    // ---------------------------------------------------------------- split

    /// Splits the edge `{a, b}` at its midpoint (`SPEC.md` §10): the midpoint is appended as a
    /// new vertex `m`; every face on the edge, in increasing index, with side `p → q` and
    /// opposite vertex `c`, becomes `(p, m, c)` and `(m, q, c)` is appended (its parent is the
    /// cut face's parent). The midpoint is not pinned; its vertex-map entry is `½ a + ½ b`.
    /// Refused only if the edge does not exist or a half would be degenerate (a nearly
    /// degenerate face, where rounding puts the midpoint on the line). Returns `m`.
    pub fn split(&mut self, a: u32, b: u32) -> Result<u32, Refused> {
        let inc = self.edge_faces(a, b);
        if inc.is_empty() {
            return refuse(
                RefusalKind::NoSuchEdge,
                format!("({a}, {b}) is not an edge"),
            );
        }
        let m = self.pos.len() as u32;
        let mid = scale(add(self.pos[a as usize], self.pos[b as usize]), 0.5);
        for &f in &inc {
            if let Some((p, q, c)) = orient(self.faces[f as usize], a, b) {
                let (xp, xq, xc) = (
                    self.pos[p as usize],
                    self.pos[q as usize],
                    self.pos[c as usize],
                );
                if degenerate([xp, mid, xc]) || degenerate([mid, xq, xc]) {
                    return refuse(
                        RefusalKind::DegenerateFace,
                        format!("splitting ({a}, {b}) cuts face {f} into a degenerate face"),
                    );
                }
            }
        }
        self.pos.push(mid);
        self.v_alive.push(true);
        self.pinned.push(false);
        self.pin_origin.push(NONE);
        let src = combine(
            &self.sources[a as usize],
            0.5,
            &self.sources[b as usize],
            0.5,
        );
        self.sources.push(src);
        self.vfaces.push(Vec::with_capacity(6));
        for f in inc {
            let Some((p, q, c)) = orient(self.faces[f as usize], a, b) else {
                continue;
            };
            let nf = self.faces.len() as u32;
            self.faces[f as usize] = [p, m, c];
            self.faces.push([m, q, c]);
            self.face_alive.push(true);
            self.parent.push(self.parent[f as usize]);
            self.remove_vface(q, f);
            self.vfaces[q as usize].push(nf);
            self.vfaces[c as usize].push(nf);
            self.vfaces[m as usize].push(f);
            self.vfaces[m as usize].push(nf);
        }
        Ok(m)
    }

    // ---------------------------------------------------------------- collapse

    /// The position a collapse of `{a, b}` moves the survivor to: the pinned endpoint, else
    /// the single boundary endpoint, else the midpoint (`SPEC.md` §10), with the weights of
    /// `a` and `b` in the vertex map.
    fn collapse_target(&self, a: u32, b: u32) -> ([f64; 3], f64, f64) {
        let (pa, pb) = (self.pinned[a as usize], self.pinned[b as usize]);
        let (ba, bb) = (self.is_boundary_vertex(a), self.is_boundary_vertex(b));
        if pa {
            (self.pos[a as usize], 1.0, 0.0)
        } else if pb {
            (self.pos[b as usize], 0.0, 1.0)
        } else if ba && !bb {
            (self.pos[a as usize], 1.0, 0.0)
        } else if bb && !ba {
            (self.pos[b as usize], 0.0, 1.0)
        } else {
            (
                scale(add(self.pos[a as usize], self.pos[b as usize]), 0.5),
                0.5,
                0.5,
            )
        }
    }

    /// Collapses the edge `{a, b}` (`SPEC.md` §10): `min(a, b)` survives at the position of
    /// [the collapse rule](Self::collapse_target), `max(a, b)` dies, the faces on the edge are
    /// deleted. Checked in this order: the edge exists, not both pinned, the link condition
    /// (the common neighbours are exactly the opposite vertices; an interior edge between two
    /// boundary vertices fails), the foldover test if asked, no repeated face, the guard if
    /// given, and the stars of the touched vertices afterwards. Returns the survivor.
    pub fn collapse<'g>(
        &mut self,
        a: u32,
        b: u32,
        checks: CollapseChecks,
        guard: Option<&mut (dyn CollapseGuard + 'g)>,
    ) -> Result<u32, Refused> {
        let inc = self.edge_faces(a, b);
        if inc.is_empty() {
            return refuse(
                RefusalKind::NoSuchEdge,
                format!("({a}, {b}) is not an edge"),
            );
        }
        if self.pinned[a as usize] && self.pinned[b as usize] {
            return refuse(
                RefusalKind::BothPinned,
                format!("both ends of ({a}, {b}) are pinned"),
            );
        }
        let na = self.neighbors(a);
        let nb = self.neighbors(b);
        let common = na.iter().filter(|w| nb.binary_search(w).is_ok()).count();
        if common != inc.len() {
            return refuse(
                RefusalKind::LinkCondition,
                format!(
                    "the ends of ({a}, {b}) have {common} common neighbours but the edge has {} \
                     faces",
                    inc.len()
                ),
            );
        }
        if inc.len() == 2 && self.is_boundary_vertex(a) && self.is_boundary_vertex(b) {
            return refuse(
                RefusalKind::LinkCondition,
                format!("({a}, {b}) is an interior edge between two boundary vertices"),
            );
        }
        let (keep, gone) = (a.min(b), a.max(b));
        let (x, wa, wb) = self.collapse_target(a, b);

        // the faces touching exactly one end, ascending
        let mut star: Vec<u32> = self.vfaces[a as usize]
            .iter()
            .chain(self.vfaces[b as usize].iter())
            .copied()
            .filter(|f| !inc.contains(f))
            .collect();
        star.sort_unstable();
        star.dedup();
        if checks.foldover {
            for &f in &star {
                let t = self.faces[f as usize];
                let before = tri_normal(self.corners(f));
                let mut p = self.corners(f);
                for k in 0..3 {
                    if t[k] == a || t[k] == b {
                        p[k] = x;
                    }
                }
                if folds(before, tri_normal(p)) {
                    return refuse(
                        RefusalKind::Foldover,
                        format!(
                            "collapsing ({a}, {b}) turns the normal of face {f} by more than \
                             pi - 2 or makes it degenerate"
                        ),
                    );
                }
            }
        }
        for &f in &star {
            let t = self.faces[f as usize];
            let mut p = self.corners(f);
            for k in 0..3 {
                if t[k] == a || t[k] == b {
                    p[k] = x;
                }
            }
            if degenerate(p) {
                return refuse(
                    RefusalKind::DegenerateFace,
                    format!("collapsing ({a}, {b}) makes face {f} degenerate"),
                );
            }
        }
        let new_star: Vec<[u32; 3]> = star
            .iter()
            .map(|&f| self.faces[f as usize].map(|v| if v == gone { keep } else { v }))
            .collect();
        for i in 0..new_star.len() {
            let si = sorted3(new_star[i]);
            for j in i + 1..new_star.len() {
                if si == sorted3(new_star[j]) {
                    return refuse(
                        RefusalKind::RepeatedFace,
                        format!(
                            "faces {} and {} would have the same three vertices",
                            star[i], star[j]
                        ),
                    );
                }
            }
        }
        if let Some(g) = guard {
            let q = CollapseQuery {
                mesh: self,
                edge: [a, b],
                survivor: keep,
                position: x,
                new_star: &new_star,
            };
            if !g.allow(&q) {
                return refuse(
                    RefusalKind::Guard,
                    format!("the collapse guard refused ({a}, {b})"),
                );
            }
        }

        // ---- apply, with a snapshot for the undo
        let mut touched_v: Vec<u32> = vec![a, b];
        for &f in &inc {
            touched_v.extend_from_slice(&self.faces[f as usize]);
        }
        touched_v.sort_unstable();
        touched_v.dedup();
        let snap_v: Vec<(u32, Vec<u32>)> = touched_v
            .iter()
            .map(|&v| (v, self.vfaces[v as usize].clone()))
            .collect();
        let snap_f: Vec<(u32, [u32; 3])> =
            star.iter().map(|&f| (f, self.faces[f as usize])).collect();
        let old_keep = (
            self.pos[keep as usize],
            self.pinned[keep as usize],
            self.pin_origin[keep as usize],
            self.sources[keep as usize].clone(),
        );

        for &f in &inc {
            self.face_alive[f as usize] = false;
            for v in self.faces[f as usize] {
                self.remove_vface(v, f);
            }
        }
        let moved = std::mem::take(&mut self.vfaces[gone as usize]);
        for &f in &moved {
            for v in self.faces[f as usize].iter_mut() {
                if *v == gone {
                    *v = keep;
                }
            }
            self.vfaces[keep as usize].push(f);
        }
        self.v_alive[gone as usize] = false;
        self.pos[keep as usize] = x;
        let pinned = self.pinned[a as usize] || self.pinned[b as usize];
        let origin = if self.pinned[a as usize] {
            self.pin_origin[a as usize]
        } else if self.pinned[b as usize] {
            self.pin_origin[b as usize]
        } else {
            NONE
        };
        self.pinned[keep as usize] = pinned;
        self.pin_origin[keep as usize] = origin;
        self.sources[keep as usize] =
            combine(&self.sources[a as usize], wa, &self.sources[b as usize], wb);

        let mut check = self.neighbors(keep);
        check.push(keep);
        check.extend(touched_v.iter().copied().filter(|&v| v != gone));
        check.sort_unstable();
        check.dedup();
        if let Some(&bad) = check.iter().find(|&&v| !self.vertex_ok(v)) {
            // undo
            for (v, l) in snap_v {
                self.vfaces[v as usize] = l;
            }
            for (f, t) in snap_f {
                self.faces[f as usize] = t;
            }
            for &f in &inc {
                self.face_alive[f as usize] = true;
            }
            self.v_alive[gone as usize] = true;
            let (p, pin, o, s) = old_keep;
            self.pos[keep as usize] = p;
            self.pinned[keep as usize] = pin;
            self.pin_origin[keep as usize] = o;
            self.sources[keep as usize] = s;
            return refuse(
                RefusalKind::NonManifoldResult,
                format!("collapsing ({a}, {b}) leaves vertex {bad} with a non-manifold star"),
            );
        }
        Ok(keep)
    }

    // ---------------------------------------------------------------- flip

    /// Flips the interior edge `{a, b}` (`SPEC.md` §10): with faces `F1 < F2`, `F1` running
    /// `p → q` with opposite vertex `c` and `d` opposite in `F2`, `F1 := (p, d, c)` and
    /// `F2 := (d, q, c)`. If `F2` also runs `p → q` (non-orientable surface), it is reordered
    /// first ([`FlipOutcome::reordered`]). Refused on a boundary edge, if `c d` is already an
    /// edge, by the optional Delaunay and foldover checks, and if the result is not a
    /// manifold.
    pub fn flip(&mut self, a: u32, b: u32, checks: FlipChecks) -> Result<FlipOutcome, Refused> {
        let inc = self.edge_faces(a, b);
        if inc.is_empty() {
            return refuse(
                RefusalKind::NoSuchEdge,
                format!("({a}, {b}) is not an edge"),
            );
        }
        if inc.len() != 2 {
            return refuse(
                RefusalKind::BoundaryEdge,
                format!(
                    "({a}, {b}) has {} face(s); only interior edges flip",
                    inc.len()
                ),
            );
        }
        let (f1, f2) = (inc[0], inc[1]);
        let (t1, t2) = (self.faces[f1 as usize], self.faces[f2 as usize]);
        let (Some((p, q, c)), Some((p2, _, d))) = (orient(t1, a, b), orient(t2, a, b)) else {
            return refuse(
                RefusalKind::NoSuchEdge,
                format!("({a}, {b}) is not an edge"),
            );
        };
        let reordered = p2 == p;
        if !self.edge_faces(c, d).is_empty() || c == d {
            return refuse(
                RefusalKind::ExistingEdge,
                format!("vertices {c} and {d} already share an edge"),
            );
        }
        let (xp, xq, xc, xd) = (
            self.pos[p as usize],
            self.pos[q as usize],
            self.pos[c as usize],
            self.pos[d as usize],
        );
        if checks.delaunay && !violates_delaunay(xp, xq, xc, xd) {
            return refuse(
                RefusalKind::NotDelaunay,
                format!("the opposite angles of ({a}, {b}) do not exceed pi"),
            );
        }
        if checks.foldover
            && flip_folds(
                self.corners(f1),
                self.corners(f2),
                reordered,
                xp,
                xq,
                xc,
                xd,
            )
        {
            return refuse(
                RefusalKind::Foldover,
                format!("flipping ({a}, {b}) turns a normal by more than pi - 2"),
            );
        }
        if degenerate([xp, xd, xc]) || degenerate([xd, xq, xc]) {
            return refuse(
                RefusalKind::DegenerateFace,
                format!("flipping ({a}, {b}) makes a degenerate face"),
            );
        }
        let snap: Vec<(u32, Vec<u32>)> = [p, q, c, d]
            .iter()
            .map(|&v| (v, self.vfaces[v as usize].clone()))
            .collect();
        self.faces[f1 as usize] = [p, d, c];
        self.faces[f2 as usize] = [d, q, c];
        self.remove_vface(p, f2);
        self.remove_vface(q, f1);
        self.vfaces[c as usize].push(f2);
        self.vfaces[d as usize].push(f1);
        if let Some(&bad) = [p, q, c, d].iter().find(|&&v| !self.vertex_ok(v)) {
            self.faces[f1 as usize] = t1;
            self.faces[f2 as usize] = t2;
            for (v, l) in snap {
                self.vfaces[v as usize] = l;
            }
            return refuse(
                RefusalKind::NonManifoldResult,
                format!("flipping ({a}, {b}) leaves vertex {bad} with a non-manifold star"),
            );
        }
        Ok(FlipOutcome {
            new_edge: [c, d],
            reordered,
        })
    }

    // ---------------------------------------------------------------- output

    /// Drops dead vertices and faces (keeping the relative order of the survivors),
    /// revalidates the result as a [`Surface`] with its positions, and returns the maps back
    /// to the input.
    ///
    /// # Errors
    /// The errors of [`Surface::new`] and [`Surface::validate_positions`] if the result is not
    /// a valid surface (an internal error of the remesher: the operations keep every star a
    /// manifold fan and refuse degenerate faces, so this is not expected).
    pub fn finish(self) -> Result<Remeshed, RsError> {
        let mut new_id = vec![NONE; self.pos.len()];
        let mut positions = Vec::new();
        let mut offsets = vec![0u32];
        let mut src = Vec::new();
        let mut wts = Vec::new();
        for v in 0..self.pos.len() {
            if self.v_alive[v] {
                new_id[v] = positions.len() as u32;
                positions.push(self.pos[v]);
                for &(o, w) in &self.sources[v] {
                    src.push(o);
                    wts.push(w);
                }
                offsets.push(src.len() as u32);
            }
        }
        let mut faces = Vec::new();
        let mut parents = Vec::new();
        for f in 0..self.faces.len() {
            if self.face_alive[f] {
                faces.push(self.faces[f].map(|v| new_id[v as usize]));
                parents.push(self.parent[f]);
            }
        }
        let surface = Surface::new(positions.len(), faces)?;
        surface.validate_positions(&positions)?;
        let mut pin_new = vec![NONE; self.n_vertices_in];
        for v in 0..self.pos.len() {
            if self.v_alive[v] && self.pin_origin[v] != NONE {
                pin_new[self.pin_origin[v] as usize] = new_id[v];
            }
        }
        let pinned = self.pins_in.iter().map(|&p| pin_new[p as usize]).collect();
        Ok(Remeshed {
            surface,
            positions,
            vertex_map: VertexMap::from_parts(offsets, src, wts, self.n_vertices_in),
            face_parents: parents,
            pinned,
        })
    }
}

/// `wa · sa + wb · sb` of two sparse combinations of old vertices (first-appearance order:
/// `sa`'s entries, then `sb`'s new ones; zero weights are kept so that the map lists every old
/// vertex an operation merged).
fn combine(sa: &[(u32, f64)], wa: f64, sb: &[(u32, f64)], wb: f64) -> Vec<(u32, f64)> {
    let mut r: Vec<(u32, f64)> = sa.iter().map(|&(o, w)| (o, w * wa)).collect();
    for &(o, w) in sb {
        if let Some(e) = r.iter_mut().find(|e| e.0 == o) {
            e.1 += w * wb;
        } else {
            r.push((o, w * wb));
        }
    }
    r
}

/// A triangle with zero or non-finite area or normal (as [`Surface::validate_positions`]).
fn degenerate(p: [[f64; 3]; 3]) -> bool {
    let (area, normal, _) = crate::geometry::face_quantities(p);
    !(area > 0.0 && area.is_finite() && normal.iter().all(|x| x.is_finite()))
}

/// The foldover test of a flip: new `F1 = (p, d, c)` against old `F1`, new `F2 = (d, q, c)`
/// against old `F2` (reversed first if it was reordered).
fn flip_folds(
    old1: [[f64; 3]; 3],
    old2: [[f64; 3]; 3],
    reordered: bool,
    xp: [f64; 3],
    xq: [f64; 3],
    xc: [f64; 3],
    xd: [f64; 3],
) -> bool {
    let n1 = tri_normal(old1);
    let mut n2 = tri_normal(old2);
    if reordered {
        n2 = scale(n2, -1.0);
    }
    folds(n1, tri_normal([xp, xd, xc])) || folds(n2, tri_normal([xd, xq, xc]))
}
