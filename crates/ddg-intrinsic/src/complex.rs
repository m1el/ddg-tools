//! The connectivity of an intrinsic triangulation, shared by its representations (SPEC N50).
//!
//! [`Complex`] is a halfedge Δ-complex: the faces as vertex triples, the twin of every
//! halfedge, the edge of every halfedge, one halfedge per edge and per vertex, and a mask of
//! *constrained* edges. It knows nothing of lengths, signposts or integers. Its topological
//! operations ([`Complex::flip_edge`], [`Complex::split_face`], [`Complex::split_edge`],
//! [`Complex::remove_vertex`], [`Complex::merge_edge`], [`Complex::remove_boundary_vertex`])
//! change only these arrays and return what they did to the ids (which halfedges moved where,
//! which edges, faces and vertices are new, renumbered or gone), so that a representation can
//! bring its own per-halfedge, per-edge, per-face and per-vertex data along.
//! [`IntrinsicTriangulation`](crate::IntrinsicTriangulation) (signposts) embeds one; the integer
//! coordinates of Gillespie, Sharp and Crane (2021) embed the same type.
//!
//! # Conventions (SPEC N9)
//!
//! Those of [`ddg_core::topology`], so that a complex made from a `TriMesh` has its ids:
//!
//! * **Halfedge** `h = 3 f + k` runs from `faces[f][k]` to `faces[f][(k + 1) % 3]`; `next` and
//!   `prev` are arithmetic. A corner and the halfedge leaving it share one id. A boundary
//!   edge has one halfedge only; there are no halfedges outside the surface.
//! * **Edges** have ids of their own (two edges may join the same two vertices).
//!   [`Complex::edge_halfedge`] gives every edge a direction.
//! * Every vertex has an outgoing halfedge ([`Complex::outgoing_halfedge`]); at a boundary
//!   vertex it is the boundary halfedge leaving it, and the outgoing halfedges of a vertex form
//!   one fan counter-clockwise from it.
//!
//! What may occur (paper Sec. 2.1, Fig. 6): two edges between the same two vertices, an edge
//! from a vertex to itself, a face that uses a vertex two or three times or is glued to itself
//! along an edge, vertices of degree one. It is always an oriented manifold.
//!
//! # What the operations keep
//!
//! * Vertex, edge and face ids are stable under flips and splits (new ones are appended);
//!   halfedge ids are not (a flip moves two halfedges between its two faces). Every result
//!   lists the halfedges whose data moved (`(from, to)` pairs): a halfedge *moved* when it
//!   still leaves the same vertex in the same direction but has another id now.
//! * **Directions.** An operation never changes the direction of an edge it does not flip or
//!   split ([`Complex::edge_halfedge`] may change id, never tail and head). A flipped edge
//!   turns counter-clockwise; the two halves of a split edge keep its direction.
//! * Removal of a vertex ([`Complex::remove_vertex`], [`Complex::remove_boundary_vertex`]) and
//!   the merge of an edge's two halves ([`Complex::merge_edge`], the inverse of a split) are the
//!   only operations that renumber: the last vertex, the last edges and the last faces take the
//!   ids that became free.

use ddg_core::TriMesh;
use ddg_core::topology::NONE;

use crate::error::{FlipRefusal, InsertRefusal, IntrinsicError, RemovalRefusal};
use crate::point::SurfacePoint;

/// The halfedge after `h` in its face.
#[inline]
pub(crate) fn next(h: usize) -> usize {
    if h % 3 == 2 { h - 2 } else { h + 1 }
}

/// The halfedge before `h` in its face.
#[inline]
pub(crate) fn prev(h: usize) -> usize {
    if h.is_multiple_of(3) { h + 2 } else { h - 1 }
}

/// The connectivity of an intrinsic triangulation: a halfedge Δ-complex with a mask of
/// constrained edges. See the [module documentation](self).
#[derive(Debug, Clone, PartialEq)]
pub struct Complex {
    pub(crate) faces: Vec<[u32; 3]>,
    /// Per halfedge: the opposite halfedge, or `NONE`.
    pub(crate) twin: Vec<u32>,
    /// Per halfedge: its edge.
    pub(crate) halfedge_edge: Vec<u32>,
    /// Per edge: the halfedge that gives it its direction (the only one on the boundary).
    pub(crate) edge_halfedge: Vec<u32>,
    /// Per vertex: an outgoing halfedge; the one without twin at a boundary vertex.
    pub(crate) vertex_halfedge: Vec<u32>,
    /// Per edge: whether flips leave it alone.
    pub(crate) constrained: Vec<bool>,
}

/// The halfedges leaving a vertex, counter-clockwise. See [`Complex::outgoing_halfedges`].
#[derive(Debug, Clone)]
pub struct OutgoingHalfedges<'a> {
    pub(crate) complex: &'a Complex,
    pub(crate) start: u32,
    pub(crate) current: u32,
}

impl Iterator for OutgoingHalfedges<'_> {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        if self.current == NONE {
            return None;
        }
        let h = self.current as usize;
        let t = self.complex.twin[prev(h)];
        self.current = if t == self.start { NONE } else { t };
        Some(h)
    }
}

/// What [`Complex::flip_edge`] did. Before: the edge's halfedge `i → j` in face `f0 = (i, j, k)`
/// and its twin `j → i` in `f1 = (j, i, l)`. After: `f0 = (i, l, k)`, `f1 = (j, k, l)` and the
/// edge runs from `l` to `k`. The faces keep their ids, and so does the edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComplexFlip {
    /// The edge.
    pub edge: u32,
    /// The two faces: that of the edge's halfedge before the flip, then the other.
    pub faces: [u32; 2],
    /// The end points before the flip: `[i, j]`.
    pub old_vertices: [u32; 2],
    /// The end points after the flip: `[l, k]`, tail and head of the edge's new halfedge.
    pub new_vertices: [u32; 2],
    /// The two halfedges of the edge after the flip: `l → k` (the edge's halfedge), `k → l`.
    /// Their data is new; what was stored in these slots before belongs to no halfedge.
    pub new_halfedges: [u32; 2],
    /// The halfedges `i → l` and `j → k` changed faces: `(from, to)`.
    pub moved: [(u32, u32); 2],
}

impl ComplexFlip {
    /// Moves the per-halfedge data of the two halfedges that changed faces. The two slots of
    /// [`new_halfedges`](Self::new_halfedges) are left as they are, for the caller to set.
    #[inline]
    pub fn move_halfedge_data<T: Copy>(&self, data: &mut [T]) {
        for (from, to) in self.moved {
            data[to as usize] = data[from as usize];
        }
    }
}

/// What a flip of the complex overwrites ([`Complex::flip_undo`]), so that
/// [`Complex::undo_flip`] restores the complex exactly, ids and all: the two faces, the twins of
/// their six halfedges and of the two outer halfedges whose twins change, the edges of the four
/// halfedges that change edge, the halfedges of three edges and of two vertices. A removal that
/// fails after some flips undoes them with these in reverse order instead of restoring a copy of
/// the whole complex (IC-SPEC §4.6; a copy costs `O(n)` per removal).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComplexFlipUndo {
    faces: [(u32, [u32; 3]); 2],
    twin: [(u32, u32); 6],
    halfedge_edge: [(u32, u32); 4],
    edge_halfedge: [(u32, u32); 3],
    vertex_halfedge: [(u32, u32); 2],
}

/// What [`Complex::split_face`] did. Before: face `f = (a, b, c)`. After: `f = (a, b, p)`,
/// `f1 = (b, c, p)`, `f2 = (c, a, p)`, with the new vertex `p`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceSplit {
    /// The new vertex `p` (the last vertex id).
    pub vertex: u32,
    /// The face that was split; it keeps its id and becomes `(a, b, p)`.
    pub face: u32,
    /// The two new faces `(b, c, p)` and `(c, a, p)`.
    pub new_faces: [u32; 2],
    /// The three new edges, from corner `k` of the old face (`a`, `b`, `c`) to the new vertex;
    /// each is directed towards `p`.
    pub new_edges: [u32; 3],
    /// The halfedges leaving `p`, towards `a`, `b`, `c` (counter-clockwise).
    pub outgoing: [u32; 3],
    /// The halfedges `b → c` and `c → a` moved to the new faces: `(from, to)`.
    pub moved: [(u32, u32); 2],
    /// The number of halfedges after the split.
    pub n_halfedges: u32,
}

impl FaceSplit {
    /// Brings per-halfedge data along: grows `data` to the new number of halfedges, filling
    /// with `fill`, and moves the data of the two halfedges that moved. The six halfedges at
    /// the new vertex are left with `fill` or stale values, for the caller to set.
    pub fn move_halfedge_data<T: Copy>(&self, data: &mut Vec<T>, fill: T) {
        data.resize(self.n_halfedges as usize, fill);
        for (from, to) in self.moved {
            data[to as usize] = data[from as usize];
        }
    }
}

/// What [`Complex::split_edge`] did. Before: edge `e` with halfedge `a → b` in face
/// `f0 = (a, b, c)` and, if it is interior, the twin `b → a` in `f1 = (b, a, d)`. After: the
/// vertex `p` lies on it; `f0 = (a, p, c)`, the new face `f2 = (p, b, c)`, and on the other
/// side `f1 = (p, a, d)` and the new face `f3 = (b, p, d)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeSplit {
    /// The new vertex `p` (the last vertex id).
    pub vertex: u32,
    /// The two halves in order from the old edge's tail, both directed like it: `[e, e2]`
    /// with `e = a → p` (the old id) and `e2 = p → b` (new).
    pub halves: [u32; 2],
    /// The new edges `c → p` and `d → p` (directed towards `p`); the second is [`NONE`] on
    /// the boundary.
    pub to_opposite: [u32; 2],
    /// The faces of the old edge, `f0` and `f1` (`NONE` on the boundary); they keep their ids
    /// and lie on the half `e`.
    pub faces: [u32; 2],
    /// The new faces `f2` and `f3` (`NONE` on the boundary); they lie on the half `e2`.
    pub new_faces: [u32; 2],
    /// Halfedges whose data moved, `(from, to)`: `b → c` to `f2`; and, for an interior edge,
    /// `d → b` to `f3` and the old twin `b → a` (now `p → a`) to `b → p` in `f3` (it leaves
    /// `b` in the same direction).
    moved: [(u32, u32); 3],
    n_moved: u8,
    /// The number of halfedges after the split.
    pub n_halfedges: u32,
    p_to_c: u32,
    p_to_a: u32,
}

impl EdgeSplit {
    /// The halfedges whose data moved, `(from, to)`.
    pub fn moved(&self) -> &[(u32, u32)] {
        &self.moved[..self.n_moved as usize]
    }

    /// The new edges in the order of [`Insertion::new_edges`](crate::Insertion::new_edges):
    /// the half `e2`, then the edges to the opposite vertices.
    pub fn new_edges(&self) -> Vec<u32> {
        let mut out = vec![self.halves[1], self.to_opposite[0]];
        if self.to_opposite[1] != NONE {
            out.push(self.to_opposite[1]);
        }
        out
    }

    /// The halfedges leaving `p` counter-clockwise: towards `b`, `c`, `a`, `d` (two fewer on
    /// the boundary, where the first is the boundary halfedge `p → b`).
    pub fn outgoing(&self) -> Vec<u32> {
        // p → b is halfedge 0 of f2; p → c sits in f0 where b was; p → a is the old twin's
        // slot; p → d is halfedge 1 of f3.
        let mut out = vec![3 * self.new_faces[0], self.p_to_c];
        if self.faces[1] != NONE {
            out.push(self.p_to_a);
            out.push(3 * self.new_faces[1] + 1);
        }
        out
    }

    /// Brings per-halfedge data along: grows `data` to the new number of halfedges, filling
    /// with `fill`, and moves the data of the halfedges that moved. The halfedges at the new
    /// vertex (both directions of every edge at it, except `a → p`, which keeps the slot and
    /// direction of `a → b`) are left for the caller to set.
    pub fn move_halfedge_data<T: Copy>(&self, data: &mut Vec<T>, fill: T) {
        data.resize(self.n_halfedges as usize, fill);
        for &(from, to) in self.moved() {
            data[to as usize] = data[from as usize];
        }
    }
}

/// What [`Complex::remove_vertex`] did: the vertex of degree three and its three edges are
/// gone, its three faces are one. Removing frees one vertex id, three edge ids and two face
/// ids; the elements with the highest ids take them (`*_moves`, `(old, new)`), so ids stay
/// contiguous.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VertexRemoval {
    /// The removed vertex (its id before the removal; after it, the id belongs to the vertex
    /// in [`vertex_move`](Self::vertex_move), if any).
    pub vertex: u32,
    /// The face that the three faces became (after renumbering).
    pub face: u32,
    /// Its corners (vertex ids after renumbering), in the order of its halfedges.
    pub corners: [u32; 3],
    /// The removed edges (ids before the removal).
    pub removed_edges: [u32; 3],
    /// The removed faces (ids before the removal); the third face of the star is
    /// [`face`](Self::face).
    pub removed_faces: [u32; 2],
    /// Every surviving halfedge whose id changed, `(old, new)`.
    pub halfedge_moves: Vec<(u32, u32)>,
    /// Every surviving edge whose id changed, `(old, new)`.
    pub edge_moves: Vec<(u32, u32)>,
    /// Every surviving face whose id changed, `(old, new)`.
    pub face_moves: Vec<(u32, u32)>,
    /// The vertex that took the freed id, `(old, new)`.
    pub vertex_move: Option<(u32, u32)>,
    /// The numbers of halfedges, edges, faces and vertices after the removal.
    pub counts: [u32; 4],
}

/// Applies `(old, new)` moves to `data` as one permutation, then truncates it.
fn permute<T: Copy>(data: &mut Vec<T>, moves: &[(u32, u32)], len: u32) {
    let values: Vec<(usize, T)> = moves
        .iter()
        .map(|&(old, new)| (new as usize, data[old as usize]))
        .collect();
    for (new, value) in values {
        data[new] = value;
    }
    data.truncate(len as usize);
}

impl VertexRemoval {
    /// Brings per-halfedge data along (`data` indexed by the old halfedge ids). The three
    /// halfedges of the merged face are the outer halfedges of the star, unchanged in
    /// direction, so their data is kept as it was.
    pub fn apply_halfedge_data<T: Copy>(&self, data: &mut Vec<T>) {
        permute(data, &self.halfedge_moves, self.counts[0]);
    }

    /// Brings per-edge data along (`data` indexed by the old edge ids).
    pub fn apply_edge_data<T: Copy>(&self, data: &mut Vec<T>) {
        permute(data, &self.edge_moves, self.counts[1]);
    }

    /// Brings per-face data along (`data` indexed by the old face ids). The merged face keeps
    /// the data of the face whose id it has.
    pub fn apply_face_data<T: Copy>(&self, data: &mut Vec<T>) {
        permute(data, &self.face_moves, self.counts[2]);
    }

    /// Brings per-vertex data along (`data` indexed by the old vertex ids).
    pub fn apply_vertex_data<T: Copy>(&self, data: &mut Vec<T>) {
        let moves: Vec<(u32, u32)> = self.vertex_move.into_iter().collect();
        permute(data, &moves, self.counts[3]);
    }
}

/// What [`Complex::merge_edge`] did: the vertex between the two halves of an edge is gone, the
/// halves are one edge again, and on each side the two faces at the vertex are one face. The
/// inverse of [`EdgeSplit`]. Merging frees one vertex id, three edge ids and two face ids (two
/// edge ids and one face id on the boundary); the elements with the highest ids take them
/// (`*_moves`, `(old, new)`), as in [`VertexRemoval`], so ids stay contiguous.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgeMerge {
    /// The removed vertex (its id before the merge; after it, the id belongs to the vertex in
    /// [`vertex_move`](Self::vertex_move), if any).
    pub vertex: u32,
    /// The merged edge (its id after renumbering).
    pub edge: u32,
    /// Tail and head of the merged edge's halfedge (vertex ids after renumbering).
    pub ends: [u32; 2],
    /// The two halves before the merge: the one whose id and direction the merged edge kept
    /// (the smaller id), then the removed one.
    pub halves: [u32; 2],
    /// The merged faces, one on each side of the merged edge (ids after renumbering): that of
    /// the merged edge's halfedge ([`Complex::edge_halfedge`]) first, the second [`NONE`] on
    /// the boundary.
    pub faces: [u32; 2],
    /// The removed edges (ids before the merge): the other half, then the edges from the
    /// vertex to the opposite corners, ascending (the last [`NONE`] on the boundary).
    pub removed_edges: [u32; 3],
    /// The removed faces (ids before the merge), ascending; the second [`NONE`] on the boundary.
    pub removed_faces: [u32; 2],
    /// Every surviving halfedge whose id changed, `(old, new)`.
    pub halfedge_moves: Vec<(u32, u32)>,
    /// Every surviving edge whose id changed, `(old, new)`.
    pub edge_moves: Vec<(u32, u32)>,
    /// Every surviving face whose id changed, `(old, new)`.
    pub face_moves: Vec<(u32, u32)>,
    /// The vertex that took the freed id, `(old, new)`.
    pub vertex_move: Option<(u32, u32)>,
    /// The numbers of halfedges, edges, faces and vertices after the merge.
    pub counts: [u32; 4],
}

impl EdgeMerge {
    /// Brings per-halfedge data along (`data` indexed by the old halfedge ids). Every surviving
    /// halfedge leaves its tail in the same direction as before (the halfedge from the far end
    /// of a half towards the vertex is the merged edge's halfedge from that end), so its data
    /// is kept as it was.
    pub fn apply_halfedge_data<T: Copy>(&self, data: &mut Vec<T>) {
        permute(data, &self.halfedge_moves, self.counts[0]);
    }

    /// Brings per-edge data along (`data` indexed by the old edge ids). The merged edge gets
    /// what is stored for the kept half ([`halves`](Self::halves)`[0]`): set the merged
    /// values there before calling this.
    pub fn apply_edge_data<T: Copy>(&self, data: &mut Vec<T>) {
        permute(data, &self.edge_moves, self.counts[1]);
    }

    /// Brings per-face data along (`data` indexed by the old face ids). Each merged face keeps
    /// the data of the face whose id it has.
    pub fn apply_face_data<T: Copy>(&self, data: &mut Vec<T>) {
        permute(data, &self.face_moves, self.counts[2]);
    }

    /// Brings per-vertex data along (`data` indexed by the old vertex ids).
    pub fn apply_vertex_data<T: Copy>(&self, data: &mut Vec<T>) {
        let moves: Vec<(u32, u32)> = self.vertex_move.into_iter().collect();
        permute(data, &moves, self.counts[3]);
    }
}

/// Records that the element now at `from` moves to `to`, composing with earlier moves.
fn record_move(moves: &mut Vec<(u32, u32)>, from: u32, to: u32) {
    if let Some(m) = moves.iter_mut().find(|m| m.1 == from) {
        m.1 = to;
    } else {
        moves.push((from, to));
    }
}

impl Complex {
    /// The complex of a `TriMesh`: same vertices, edges, faces and halfedges, with the same
    /// ids; no edge constrained.
    ///
    /// [formula: it-complex-from-mesh]
    pub fn from_tri_mesh(mesh: &TriMesh) -> Complex {
        let n_halfedges = mesh.n_halfedges();
        Complex {
            faces: mesh.faces().to_vec(),
            twin: (0..n_halfedges)
                .map(|h| mesh.twin(h).map_or(NONE, |t| t as u32))
                .collect(),
            halfedge_edge: (0..n_halfedges).map(|h| mesh.edge(h) as u32).collect(),
            edge_halfedge: (0..mesh.n_edges())
                .map(|e| mesh.edge_halfedge(e) as u32)
                .collect(),
            vertex_halfedge: (0..mesh.n_vertices())
                .map(|v| mesh.outgoing_halfedge(v) as u32)
                .collect(),
            constrained: vec![false; mesh.n_edges()],
        }
    }

    // ----- counts -----

    /// Number of vertices.
    #[inline]
    pub fn n_vertices(&self) -> usize {
        self.vertex_halfedge.len()
    }

    /// Number of edges.
    #[inline]
    pub fn n_edges(&self) -> usize {
        self.edge_halfedge.len()
    }

    /// Number of faces.
    #[inline]
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Number of halfedges, which is also the number of corners: `3 * n_faces`.
    #[inline]
    pub fn n_halfedges(&self) -> usize {
        self.twin.len()
    }

    /// Euler characteristic `χ = V − E + F`. No operation changes it.
    ///
    /// [formula: it-complex-euler]
    pub fn euler_characteristic(&self) -> i64 {
        self.n_vertices() as i64 - self.n_edges() as i64 + self.n_faces() as i64
    }

    // ----- halfedges -----

    /// The faces as vertex triples. A face may list a vertex more than once.
    #[inline]
    pub fn faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    /// The next halfedge in the same face.
    #[inline]
    pub fn next(&self, h: usize) -> usize {
        next(h)
    }

    /// The previous halfedge in the same face.
    #[inline]
    pub fn prev(&self, h: usize) -> usize {
        prev(h)
    }

    /// The oppositely directed halfedge of the same edge, or `None` on the boundary. It can
    /// lie in the same face.
    #[inline]
    pub fn twin(&self, h: usize) -> Option<usize> {
        let t = self.twin[h];
        (t != NONE).then_some(t as usize)
    }

    /// The face of halfedge (or corner) `h`: `h / 3`.
    #[inline]
    pub fn face(&self, h: usize) -> usize {
        h / 3
    }

    /// The vertex halfedge `h` starts at; equivalently the vertex of corner `h`.
    #[inline]
    pub fn from_vertex(&self, h: usize) -> usize {
        self.faces[h / 3][h % 3] as usize
    }

    /// The vertex halfedge `h` points to.
    #[inline]
    pub fn to_vertex(&self, h: usize) -> usize {
        self.faces[h / 3][(h + 1) % 3] as usize
    }

    /// The edge halfedge `h` belongs to.
    #[inline]
    pub fn edge(&self, h: usize) -> usize {
        self.halfedge_edge[h] as usize
    }

    /// Whether halfedge `h` runs in the direction of its edge (it is
    /// [`edge_halfedge`](Self::edge_halfedge) of its edge).
    #[inline]
    pub fn is_forward(&self, h: usize) -> bool {
        self.edge_halfedge[self.halfedge_edge[h] as usize] == h as u32
    }

    // ----- edges -----

    /// The halfedge that gives edge `e` its direction. Its twin, if any, is the other halfedge
    /// of the edge.
    #[inline]
    pub fn edge_halfedge(&self, e: usize) -> usize {
        self.edge_halfedge[e] as usize
    }

    /// The end points of edge `e`: tail and head of [`edge_halfedge`](Self::edge_halfedge).
    /// They are equal for a self-edge.
    #[inline]
    pub fn edge_vertices(&self, e: usize) -> [usize; 2] {
        let h = self.edge_halfedge[e] as usize;
        [self.from_vertex(h), self.to_vertex(h)]
    }

    /// Whether edge `e` has only one incident triangle side.
    #[inline]
    pub fn is_boundary_edge(&self, e: usize) -> bool {
        self.twin[self.edge_halfedge[e] as usize] == NONE
    }

    /// Whether edge `e` is constrained: flips and Delaunay flipping leave it alone, as they
    /// leave boundary edges (a boundary edge need not be marked).
    #[inline]
    pub fn is_constrained(&self, e: usize) -> bool {
        self.constrained[e]
    }

    /// The constrained-edge mask, indexed by edge id.
    #[inline]
    pub fn constrained_edges(&self) -> &[bool] {
        &self.constrained
    }

    /// Number of constrained edges.
    pub fn n_constrained(&self) -> usize {
        self.constrained.iter().filter(|&&c| c).count()
    }

    /// Marks edge `e` as constrained or not. The mask follows the edge through operations:
    /// both halves of a split constrained edge are constrained, new edges are not.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`].
    pub fn set_constrained(&mut self, e: usize, constrained: bool) -> Result<(), IntrinsicError> {
        self.check_index("edge", e, self.n_edges())?;
        self.constrained[e] = constrained;
        Ok(())
    }

    // ----- vertices -----

    /// An outgoing halfedge of vertex `v`: the first one of
    /// [`outgoing_halfedges`](Self::outgoing_halfedges). At a boundary vertex it is the
    /// boundary halfedge leaving `v`.
    #[inline]
    pub fn outgoing_halfedge(&self, v: usize) -> usize {
        self.vertex_halfedge[v] as usize
    }

    /// Whether vertex `v` lies on the boundary.
    #[inline]
    pub fn is_boundary_vertex(&self, v: usize) -> bool {
        self.twin[self.vertex_halfedge[v] as usize] == NONE
    }

    /// The halfedges leaving `v`, one per incident corner, counter-clockwise. At a boundary
    /// vertex the iteration starts with the boundary halfedge leaving `v`; the last boundary
    /// neighbour is then reached only by an incoming halfedge. A self-edge at `v` appears with
    /// each of its halfedges.
    pub fn outgoing_halfedges(&self, v: usize) -> OutgoingHalfedges<'_> {
        OutgoingHalfedges {
            complex: self,
            start: self.vertex_halfedge[v],
            current: self.vertex_halfedge[v],
        }
    }

    /// The degree of `v`: the number of edge ends at `v` (a self-edge counts twice).
    pub fn vertex_degree(&self, v: usize) -> usize {
        self.outgoing_halfedges(v).count() + usize::from(self.is_boundary_vertex(v))
    }

    fn check_index(
        &self,
        what: &'static str,
        index: usize,
        bound: usize,
    ) -> Result<(), IntrinsicError> {
        if index < bound {
            Ok(())
        } else {
            Err(IntrinsicError::IndexOutOfRange { what, index, bound })
        }
    }

    // ----- consistency -----

    /// Checks the structure: array lengths; twins are mutual, run in opposite directions and
    /// share their edge; every edge has one or two halfedges and owns its
    /// [`edge_halfedge`](Self::edge_halfedge); the corners of every vertex form one fan that
    /// [`outgoing_halfedge`](Self::outgoing_halfedge) starts, at a boundary vertex with the
    /// boundary halfedge.
    ///
    /// # Errors
    /// [`IntrinsicError::Inconsistent`] with the first violation found.
    pub fn check(&self) -> Result<(), IntrinsicError> {
        let bad = |what: String| Err(IntrinsicError::Inconsistent { what });
        let (nv, ne, nf, nh) = (
            self.n_vertices(),
            self.n_edges(),
            self.n_faces(),
            self.n_halfedges(),
        );
        if nh != 3 * nf || self.halfedge_edge.len() != nh || self.constrained.len() != ne {
            return bad("array lengths do not match the element counts".to_string());
        }
        for (f, face) in self.faces.iter().enumerate() {
            if let Some(&v) = face.iter().find(|&&v| v as usize >= nv) {
                return bad(format!("face {f} references vertex {v} of {nv}"));
            }
        }
        let mut halfedges_of_edge = vec![0u8; ne];
        for h in 0..nh {
            let e = self.halfedge_edge[h];
            if e as usize >= ne {
                return bad(format!("halfedge {h} references edge {e} of {ne}"));
            }
            halfedges_of_edge[e as usize] = halfedges_of_edge[e as usize].saturating_add(1);
            let t = self.twin[h];
            if t == NONE {
                continue;
            }
            let t = t as usize;
            if t >= nh || t == h || self.twin[t] != h as u32 {
                return bad(format!("the twin of halfedge {h} is not mutual"));
            }
            if self.halfedge_edge[t] != e {
                return bad(format!(
                    "halfedge {h} and its twin {t} have different edges"
                ));
            }
            if self.from_vertex(t) != self.to_vertex(h) || self.to_vertex(t) != self.from_vertex(h)
            {
                return bad(format!(
                    "halfedge {h} and its twin {t} do not run in opposite directions"
                ));
            }
        }
        for (e, &n_halfedges) in halfedges_of_edge.iter().enumerate() {
            let h = self.edge_halfedge[e] as usize;
            if h >= nh || self.halfedge_edge[h] != e as u32 {
                return bad(format!("edge {e} does not own its halfedge {h}"));
            }
            let expected = if self.twin[h] == NONE { 1 } else { 2 };
            if n_halfedges != expected {
                return bad(format!(
                    "edge {e} has {n_halfedges} halfedges, expected {expected}"
                ));
            }
        }
        let mut corners_at = vec![0u32; nv];
        for face in &self.faces {
            for &v in face {
                corners_at[v as usize] += 1;
            }
        }
        for (v, &n_corners) in corners_at.iter().enumerate() {
            let h0 = self.vertex_halfedge[v] as usize;
            if h0 >= nh || self.from_vertex(h0) != v {
                return bad(format!("vertex {v} does not own its halfedge {h0}"));
            }
            let boundary = self.twin[h0] == NONE;
            let mut count = 0u32;
            let mut h = h0;
            loop {
                count += 1;
                if count > n_corners {
                    return bad(format!("the corners of vertex {v} do not form a fan"));
                }
                if self.from_vertex(h) != v {
                    return bad(format!(
                        "halfedge {h} in the fan of vertex {v} leaves another vertex"
                    ));
                }
                if h != h0 && self.twin[h] == NONE {
                    return bad(format!(
                        "boundary halfedge {h} is not the first halfedge of vertex {v}"
                    ));
                }
                let t = self.twin[prev(h)];
                if t == NONE {
                    if !boundary {
                        return bad(format!("the fan of interior vertex {v} is open"));
                    }
                    break;
                }
                h = t as usize;
                if h == h0 {
                    if boundary {
                        return bad(format!("the fan of boundary vertex {v} is closed"));
                    }
                    break;
                }
            }
            if count != n_corners {
                return bad(format!(
                    "vertex {v} has {n_corners} corners but its fan has {count}"
                ));
            }
        }
        Ok(())
    }

    // ----- flips -----

    /// The topological reasons why edge `e` (in range) cannot be flipped, in this order: a
    /// boundary edge, a constrained edge, an end point of degree one (both sides of the edge
    /// are the same face). Geometric reasons are the representation's.
    pub fn flip_refusal(&self, e: usize) -> Option<FlipRefusal> {
        let h0 = self.edge_halfedge[e] as usize;
        let h1 = self.twin[h0];
        if h1 == NONE {
            return Some(FlipRefusal::BoundaryEdge);
        }
        if self.constrained[e] {
            return Some(FlipRefusal::Constrained);
        }
        let h1 = h1 as usize;
        if h0 / 3 == h1 / 3 {
            // The face is glued to itself along the edge. Its third halfedge is a loop at one
            // end point; the other end point has this one corner only.
            let loop_halfedge = 3 * (h0 / 3) + (3 - h0 % 3 - h1 % 3);
            let at_loop = self.from_vertex(loop_halfedge);
            let [a, b] = [self.from_vertex(h0), self.to_vertex(h0)];
            return Some(FlipRefusal::DegreeOneEndpoint {
                vertex: (if a == at_loop { b } else { a }) as u32,
            });
        }
        None
    }

    /// Flips edge `e` (paper Fig. 6): its two faces `ijk`, `jil` become `ilk`, `jkl`. See
    /// [`ComplexFlip`] for the ids. Only the topology is checked ([`flip_refusal`]); whether
    /// the flip is geometrically possible is for the representation to decide first.
    ///
    /// [`flip_refusal`]: Self::flip_refusal
    ///
    /// [formula: it-complex-flip]
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`] and [`IntrinsicError::FlipRefused`]; the complex is
    /// unchanged then.
    pub fn flip_edge(&mut self, e: usize) -> Result<ComplexFlip, IntrinsicError> {
        self.check_index("edge", e, self.n_edges())?;
        if let Some(reason) = self.flip_refusal(e) {
            let [a, b] = self.edge_vertices(e);
            return Err(IntrinsicError::FlipRefused {
                edge: e as u32,
                vertices: [a as u32, b as u32],
                reason,
            });
        }
        Ok(self.flip_unchecked(e))
    }

    /// [`flip_edge`](Self::flip_edge) for an edge that [`flip_refusal`](Self::flip_refusal)
    /// accepts.
    pub(crate) fn flip_unchecked(&mut self, e: usize) -> ComplexFlip {
        let h0 = self.edge_halfedge[e] as usize;
        let h1 = self.twin[h0] as usize;
        let (f0, a) = (h0 / 3, h0 % 3);
        let (f1, b) = (h1 / 3, h1 % 3);
        // Before:  f0 = (i, j, k) with halfedges h0 = i→j, h0n = j→k, h0p = k→i,
        //          f1 = (j, i, l) with halfedges h1 = j→i, h1n = i→l, h1p = l→j.
        // After:   f0 = (i, l, k) with halfedges h0 = i→l, h0n = l→k, h0p = k→i,
        //          f1 = (j, k, l) with halfedges h1 = j→k, h1n = k→l, h1p = l→j.
        // The halfedges k→i and l→j stay where they are; i→l and j→k change faces.
        let (h0n, h1n) = (3 * f0 + (a + 1) % 3, 3 * f1 + (b + 1) % 3);
        let (i, j) = (self.faces[f0][a], self.faces[f1][b]);
        let (k, l) = (self.faces[f0][(a + 2) % 3], self.faces[f1][(b + 2) % 3]);

        let moved = |h: u32| {
            if h == h0n as u32 {
                h1 as u32
            } else if h == h1n as u32 {
                h0 as u32
            } else {
                h
            }
        };
        let twin_il = moved(self.twin[h1n]);
        let twin_jk = moved(self.twin[h0n]);
        let (edge_il, edge_jk) = (self.halfedge_edge[h1n], self.halfedge_edge[h0n]);

        self.faces[f0][(a + 1) % 3] = l;
        self.faces[f1][(b + 1) % 3] = k;

        self.twin[h0] = twin_il;
        if twin_il != NONE {
            self.twin[twin_il as usize] = h0 as u32;
        }
        self.twin[h1] = twin_jk;
        if twin_jk != NONE {
            self.twin[twin_jk as usize] = h1 as u32;
        }
        self.twin[h0n] = h1n as u32;
        self.twin[h1n] = h0n as u32;

        self.halfedge_edge[h0] = edge_il;
        self.halfedge_edge[h1] = edge_jk;
        self.halfedge_edge[h0n] = e as u32;
        self.halfedge_edge[h1n] = e as u32;
        if self.edge_halfedge[edge_il as usize] == h1n as u32 {
            self.edge_halfedge[edge_il as usize] = h0 as u32;
        }
        if self.edge_halfedge[edge_jk as usize] == h0n as u32 {
            self.edge_halfedge[edge_jk as usize] = h1 as u32;
        }
        self.edge_halfedge[e] = h0n as u32;

        if self.vertex_halfedge[i as usize] == h1n as u32 {
            self.vertex_halfedge[i as usize] = h0 as u32;
        }
        if self.vertex_halfedge[j as usize] == h0n as u32 {
            self.vertex_halfedge[j as usize] = h1 as u32;
        }
        ComplexFlip {
            edge: e as u32,
            faces: [f0 as u32, f1 as u32],
            old_vertices: [i, j],
            new_vertices: [l, k],
            new_halfedges: [h0n as u32, h1n as u32],
            moved: [(h1n as u32, h0 as u32), (h0n as u32, h1 as u32)],
        }
    }

    /// What [`flip_edge`](Self::flip_edge) of edge `e` would overwrite, for
    /// [`undo_flip`](Self::undo_flip). Call it right before the flip, with `e` an edge that
    /// [`flip_refusal`](Self::flip_refusal) accepts (interior).
    pub fn flip_undo(&self, e: usize) -> ComplexFlipUndo {
        let h0 = self.edge_halfedge[e] as usize;
        let h1 = self.twin[h0] as usize;
        let (f0, a) = (h0 / 3, h0 % 3);
        let (f1, b) = (h1 / 3, h1 % 3);
        let (h0n, h1n) = (3 * f0 + (a + 1) % 3, 3 * f1 + (b + 1) % 3);
        // The outer halfedges whose twins the flip sets (as in `flip_unchecked`).
        let moved = |h: u32| {
            if h == h0n as u32 {
                h1 as u32
            } else if h == h1n as u32 {
                h0 as u32
            } else {
                h
            }
        };
        let outer = |t: u32| {
            let t = moved(t);
            if t == NONE { h0 as u32 } else { t }
        };
        let (twin_il, twin_jk) = (outer(self.twin[h1n]), outer(self.twin[h0n]));
        let twin = [
            h0 as u32, h1 as u32, h0n as u32, h1n as u32, twin_il, twin_jk,
        ]
        .map(|h| (h, self.twin[h as usize]));
        let halfedge_edge = [h0, h1, h0n, h1n].map(|h| (h as u32, self.halfedge_edge[h]));
        let edge_halfedge = [self.halfedge_edge[h1n], self.halfedge_edge[h0n], e as u32]
            .map(|x| (x, self.edge_halfedge[x as usize]));
        let vertex_halfedge =
            [self.faces[f0][a], self.faces[f1][b]].map(|v| (v, self.vertex_halfedge[v as usize]));
        ComplexFlipUndo {
            faces: [f0, f1].map(|f| (f as u32, self.faces[f])),
            twin,
            halfedge_edge,
            edge_halfedge,
            vertex_halfedge,
        }
    }

    /// Restores what a flip overwrote ([`flip_undo`](Self::flip_undo), taken right before it).
    /// Flips undone in the reverse order of the flips restore the complex exactly.
    pub fn undo_flip(&mut self, undo: &ComplexFlipUndo) {
        for &(f, corners) in &undo.faces {
            self.faces[f as usize] = corners;
        }
        for &(h, t) in &undo.twin {
            self.twin[h as usize] = t;
        }
        for &(h, e) in &undo.halfedge_edge {
            self.halfedge_edge[h as usize] = e;
        }
        for &(e, h) in &undo.edge_halfedge {
            self.edge_halfedge[e as usize] = h;
        }
        for &(v, h) in &undo.vertex_halfedge {
            self.vertex_halfedge[v as usize] = h;
        }
    }

    // ----- splits -----

    /// Appends a face; its halfedges are set by the caller.
    fn push_face(&mut self, corners: [u32; 3]) -> usize {
        self.faces.push(corners);
        self.twin.extend([NONE; 3]);
        self.halfedge_edge.extend([NONE; 3]);
        self.faces.len() - 1
    }

    /// Appends an edge whose direction is that of `halfedge`.
    fn push_edge(&mut self, halfedge: usize, constrained: bool) -> u32 {
        self.edge_halfedge.push(halfedge as u32);
        self.constrained.push(constrained);
        (self.edge_halfedge.len() - 1) as u32
    }

    /// Moves the twin and the edge of halfedge `from` to the slot `to`, whose face already
    /// lists the tail, and repairs what pointed at `from`.
    fn move_halfedge(&mut self, from: usize, to: usize) {
        let t = self.twin[from];
        self.twin[to] = t;
        if t != NONE {
            self.twin[t as usize] = to as u32;
        }
        let e = self.halfedge_edge[from];
        self.halfedge_edge[to] = e;
        if self.edge_halfedge[e as usize] == from as u32 {
            self.edge_halfedge[e as usize] = to as u32;
        }
        let tail = self.faces[to / 3][to % 3] as usize;
        if self.vertex_halfedge[tail] == from as u32 {
            self.vertex_halfedge[tail] = to as u32;
        }
    }

    /// Makes `h` and `t` the two halfedges of `edge` (`t` may be [`NONE`]).
    fn pair(&mut self, h: usize, t: u32, edge: u32) {
        self.twin[h] = t;
        self.halfedge_edge[h] = edge;
        if t != NONE {
            self.twin[t as usize] = h as u32;
            self.halfedge_edge[t as usize] = edge;
        }
    }

    /// Splits face `f` at a new vertex joined to its three corners (paper Alg. 7, Fig. 7).
    /// See [`FaceSplit`] for the ids: the face keeps its id, two faces, three edges and one
    /// vertex are appended.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`].
    pub fn split_face(&mut self, f: usize) -> Result<FaceSplit, IntrinsicError> {
        self.check_index("face", f, self.n_faces())?;
        Ok(self.split_face_unchecked(f))
    }

    pub(crate) fn split_face_unchecked(&mut self, f: usize) -> FaceSplit {
        // Before: f = (a, b, c). After: f = (a, b, p), f1 = (b, c, p), f2 = (c, a, p).
        let [a, b, c] = self.faces[f];
        let p = self.n_vertices() as u32;
        let f1 = self.push_face([b, c, p]);
        let f2 = self.push_face([c, a, p]);
        self.move_halfedge(3 * f + 1, 3 * f1);
        self.move_halfedge(3 * f + 2, 3 * f2);
        self.faces[f] = [a, b, p];
        // b → p in f and p → b in f1; c → p in f1 and p → c in f2; a → p in f2 and p → a in f.
        let ea = self.push_edge(3 * f2 + 1, false);
        self.pair(3 * f2 + 1, (3 * f + 2) as u32, ea);
        let eb = self.push_edge(3 * f + 1, false);
        self.pair(3 * f + 1, (3 * f1 + 2) as u32, eb);
        let ec = self.push_edge(3 * f1 + 1, false);
        self.pair(3 * f1 + 1, (3 * f2 + 2) as u32, ec);
        self.vertex_halfedge.push((3 * f + 2) as u32);
        FaceSplit {
            vertex: p,
            face: f as u32,
            new_faces: [f1 as u32, f2 as u32],
            new_edges: [ea, eb, ec],
            outgoing: [(3 * f + 2) as u32, (3 * f1 + 2) as u32, (3 * f2 + 2) as u32],
            moved: [
                ((3 * f + 1) as u32, (3 * f1) as u32),
                ((3 * f + 2) as u32, (3 * f2) as u32),
            ],
            n_halfedges: self.n_halfedges() as u32,
        }
    }

    /// Splits edge `e` at a new vertex joined to the one or two opposite vertices (paper
    /// Alg. 8). See [`EdgeSplit`] for the ids: the half from the tail of the edge's halfedge
    /// keeps the id `e`, and the faces of the edge keep theirs. Both halves of a constrained
    /// edge are constrained.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`], and [`IntrinsicError::InsertionRefused`] with
    /// [`InsertRefusal::SameFaceOnBothSides`] for an edge with the same face on both sides.
    pub fn split_edge(&mut self, e: usize) -> Result<EdgeSplit, IntrinsicError> {
        self.check_index("edge", e, self.n_edges())?;
        if self.split_edge_refusal(e).is_some() {
            return Err(IntrinsicError::InsertionRefused {
                at: SurfacePoint::Edge {
                    edge: e as u32,
                    t: 0.5,
                },
                reason: InsertRefusal::SameFaceOnBothSides,
            });
        }
        Ok(self.split_edge_unchecked(e))
    }

    /// Why edge `e` (in range) cannot be split topologically.
    pub(crate) fn split_edge_refusal(&self, e: usize) -> Option<InsertRefusal> {
        let h0 = self.edge_halfedge[e] as usize;
        let h1 = self.twin[h0];
        (h1 != NONE && h1 as usize / 3 == h0 / 3).then_some(InsertRefusal::SameFaceOnBothSides)
    }

    pub(crate) fn split_edge_unchecked(&mut self, e: usize) -> EdgeSplit {
        let h0 = self.edge_halfedge[e] as usize;
        let (f0, k0) = (h0 / 3, h0 % 3);
        let h1 = self.twin[h0];
        let b = self.faces[f0][(k0 + 1) % 3];
        let constrained = self.constrained[e];
        // Before: f0 = (a, b, c) with a → b in slot k0. After: f0 = (a, p, c) and the new
        // face f2 = (p, b, c).
        let p = self.n_vertices() as u32;
        let c = self.faces[f0][(k0 + 2) % 3];
        let f2 = self.push_face([p, b, c]);
        let mut moved = [(NONE, NONE); 3];
        let mut n_moved = 0u8;
        let b_to_c = 3 * f0 + (k0 + 1) % 3;
        self.move_halfedge(b_to_c, 3 * f2 + 1);
        moved[0] = (b_to_c as u32, (3 * f2 + 1) as u32);
        n_moved += 1;
        self.faces[f0][(k0 + 1) % 3] = p;
        let p_to_c = b_to_c;
        // The half p → b, and its twin b → p in the new face f3 = (b, p, d).
        let e2 = self.push_edge(3 * f2, constrained);
        let (mut f3, mut ed, mut p_to_a) = (NONE, NONE, NONE);
        if h1 == NONE {
            self.pair(3 * f2, NONE, e2);
        } else {
            let h1 = h1 as usize;
            let (f1, k1) = (h1 / 3, h1 % 3);
            let d = self.faces[f1][(k1 + 2) % 3];
            // Before: f1 = (b, a, d) with b → a in slot k1. After: f1 = (p, a, d).
            let g3 = self.push_face([b, p, d]);
            let d_to_b = 3 * f1 + (k1 + 2) % 3;
            self.move_halfedge(d_to_b, 3 * g3 + 2);
            moved[1] = (d_to_b as u32, (3 * g3 + 2) as u32);
            moved[2] = (h1 as u32, (3 * g3) as u32);
            n_moved += 2;
            self.faces[f1][k1] = p;
            self.pair(3 * f2, (3 * g3) as u32, e2);
            if self.vertex_halfedge[b as usize] == h1 as u32 {
                self.vertex_halfedge[b as usize] = (3 * g3) as u32;
            }
            // d → p in f1, p → d in f3.
            let d_to_p = d_to_b;
            ed = self.push_edge(d_to_p, false);
            self.pair(d_to_p, (3 * g3 + 1) as u32, ed);
            f3 = g3 as u32;
            p_to_a = h1 as u32;
        }
        // c → p in f2, p → c in f0.
        let ec = self.push_edge(3 * f2 + 2, false);
        self.pair(3 * f2 + 2, p_to_c as u32, ec);
        self.vertex_halfedge.push((3 * f2) as u32);
        EdgeSplit {
            vertex: p,
            halves: [e as u32, e2],
            to_opposite: [ec, ed],
            faces: [f0 as u32, if h1 == NONE { NONE } else { h1 / 3 }],
            new_faces: [f2 as u32, f3],
            moved,
            n_moved,
            n_halfedges: self.n_halfedges() as u32,
            p_to_c: p_to_c as u32,
            p_to_a,
        }
    }

    // ----- removal -----

    /// Why vertex `v` (in range) cannot be removed by [`remove_vertex`](Self::remove_vertex),
    /// and otherwise its three outgoing halfedges counter-clockwise from
    /// [`outgoing_halfedge`](Self::outgoing_halfedge).
    pub(crate) fn removal_star(&self, v: usize) -> Result<[usize; 3], RemovalRefusal> {
        let h0 = self.vertex_halfedge[v] as usize;
        if self.twin[h0] == NONE {
            return Err(RemovalRefusal::BoundaryVertex);
        }
        let mut out = [0usize; 3];
        let mut n = 0usize;
        for h in self.outgoing_halfedges(v) {
            if n < 3 {
                out[n] = h;
            }
            n += 1;
            if n > 3 {
                break;
            }
        }
        if n != 3 {
            return Err(RemovalRefusal::Degree {
                degree: self.vertex_degree(v),
            });
        }
        for &h in &out {
            if self.to_vertex(h) == v {
                return Err(RemovalRefusal::SelfEdge {
                    edge: self.halfedge_edge[h],
                });
            }
        }
        for &h in &out {
            let e = self.halfedge_edge[h];
            if self.constrained[e as usize] {
                return Err(RemovalRefusal::ConstrainedEdge { edge: e });
            }
        }
        Ok(out)
    }

    /// Removes vertex `v`, an interior vertex of degree three (paper Sec. 3.5 of Gillespie,
    /// Sharp and Crane 2021; Alg. 8 flips a vertex down to degree three first, which is the
    /// representation's business): its three edges are deleted and its three faces become one
    /// face whose halfedges are the outer halfedges of the star, unchanged in direction. The
    /// merged face takes the smallest of the three face ids; the freed ids are filled by the
    /// last vertex, the last edges and the last faces (see [`VertexRemoval`]).
    ///
    /// Not supported: a boundary vertex, a vertex of another degree, a vertex with an edge to
    /// itself, a vertex with a constrained edge.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`], and [`IntrinsicError::RemovalRefused`] with the
    /// reason; the complex is unchanged then.
    pub fn remove_vertex(&mut self, v: usize) -> Result<VertexRemoval, IntrinsicError> {
        self.check_index("vertex", v, self.n_vertices())?;
        let o = self
            .removal_star(v)
            .map_err(|reason| IntrinsicError::RemovalRefused {
                vertex: v as u32,
                reason,
            })?;
        Ok(self.remove_vertex_unchecked(v, o))
    }

    pub(crate) fn remove_vertex_unchecked(&mut self, v: usize, o: [usize; 3]) -> VertexRemoval {
        // Face f_i = (v, x_i, x_{i+1}) holds o_i = v → x_i and the outer halfedge
        // q_i = x_i → x_{i+1}.
        let f = o.map(|h| h / 3);
        let q = o.map(next);
        let x = q.map(|h| self.faces[h / 3][h % 3]);
        let spokes = o.map(|h| self.halfedge_edge[h]);
        let keep = f[0].min(f[1]).min(f[2]);
        let nq = [0, 1, 2].map(|i| 3 * keep + i);
        let q_twin = q.map(|h| self.twin[h]);
        let q_edge = q.map(|h| self.halfedge_edge[h]);
        let q_edge_halfedge = q_edge.map(|e| self.edge_halfedge[e as usize]);
        let in_star = |h: u32| h != NONE && f.contains(&(h as usize / 3));
        let renamed = |h: u32| match q.iter().position(|&g| g as u32 == h) {
            Some(i) => nq[i] as u32,
            None => h,
        };
        // Outgoing halfedges of the corners that are about to disappear: x → v becomes the
        // outer halfedge leaving x.
        let corner_halfedge = x.map(|xi| {
            let h = self.vertex_halfedge[xi as usize];
            if !in_star(h) {
                return h;
            }
            match q.iter().position(|&g| g as u32 == h) {
                Some(i) => nq[i] as u32,
                // h = x → v, the halfedge before o_{i} in face f_{i-1}: its tail is x_i.
                None => {
                    let i = (0..3).find(|&i| x[i] == xi).unwrap_or(0);
                    nq[i] as u32
                }
            }
        });

        // Merge the star into face `keep`.
        self.faces[keep] = x;
        for i in 0..3 {
            self.twin[nq[i]] = renamed(q_twin[i]);
            self.halfedge_edge[nq[i]] = q_edge[i];
            if q_twin[i] != NONE && !in_star(q_twin[i]) {
                self.twin[q_twin[i] as usize] = nq[i] as u32;
            }
            self.edge_halfedge[q_edge[i] as usize] = renamed(q_edge_halfedge[i]);
            self.vertex_halfedge[x[i] as usize] = corner_halfedge[i];
        }
        let mut halfedge_moves = Vec::new();
        for i in 0..3 {
            if q[i] != nq[i] {
                halfedge_moves.push((q[i] as u32, nq[i] as u32));
            }
        }
        // Halfedges of the removed faces that are not outer ones are dead; a move into the
        // merged face from a dead slot was recorded above.
        let dead_faces: Vec<usize> = f.iter().copied().filter(|&g| g != keep).collect();
        let (dead_faces, face_moves, edge_moves, vertex_move) =
            self.compact(v, dead_faces, spokes.to_vec(), &mut halfedge_moves);
        let last = self.n_vertices();
        let rename_vertex = |u: u32| if u as usize == last { v as u32 } else { u };
        let mut removed_faces = [0u32; 2];
        for (slot, &d) in removed_faces.iter_mut().zip(&dead_faces) {
            *slot = d as u32;
        }
        removed_faces.sort_unstable();
        let mut removed_edges = spokes;
        removed_edges.sort_unstable();
        VertexRemoval {
            vertex: v as u32,
            face: keep as u32,
            corners: x.map(rename_vertex),
            removed_edges,
            removed_faces,
            halfedge_moves,
            edge_moves,
            face_moves,
            vertex_move,
            counts: [
                self.n_halfedges() as u32,
                self.n_edges() as u32,
                self.n_faces() as u32,
                self.n_vertices() as u32,
            ],
        }
    }

    /// Deletes the faces `dead_faces`, the edges `dead_edges` and the vertex `v`, which nothing
    /// alive refers to any more, and fills the freed ids with the last faces, edges and vertex
    /// (in decreasing order of the freed ids, so ids stay contiguous). The halfedge moves of the
    /// faces that moved are composed into `halfedge_moves`. Returns the dead faces in decreasing
    /// order, the face moves, the edge moves and the vertex move (all `(old, new)`).
    #[allow(clippy::type_complexity)]
    fn compact(
        &mut self,
        v: usize,
        mut dead_faces: Vec<usize>,
        mut dead_edges: Vec<u32>,
        halfedge_moves: &mut Vec<(u32, u32)>,
    ) -> (
        Vec<usize>,
        Vec<(u32, u32)>,
        Vec<(u32, u32)>,
        Option<(u32, u32)>,
    ) {
        dead_faces.sort_unstable_by(|a, b| b.cmp(a));
        let mut face_moves = Vec::new();
        for &d in &dead_faces {
            let last = self.n_faces() - 1;
            if d != last {
                self.move_face(last, d, halfedge_moves);
                record_move(&mut face_moves, last as u32, d as u32);
            }
            self.faces.truncate(last);
            self.twin.truncate(3 * last);
            self.halfedge_edge.truncate(3 * last);
        }
        dead_edges.sort_unstable_by(|a, b| b.cmp(a));
        let mut edge_moves = Vec::new();
        for &d in &dead_edges {
            let last = self.n_edges() - 1;
            if d as usize != last {
                let h = self.edge_halfedge[last] as usize;
                self.edge_halfedge[d as usize] = h as u32;
                self.constrained[d as usize] = self.constrained[last];
                self.halfedge_edge[h] = d;
                if self.twin[h] != NONE {
                    self.halfedge_edge[self.twin[h] as usize] = d;
                }
                record_move(&mut edge_moves, last as u32, d);
            }
            self.edge_halfedge.truncate(last);
            self.constrained.truncate(last);
        }
        let last = self.n_vertices() - 1;
        let mut vertex_move = None;
        if v != last {
            let fan: Vec<usize> = self.outgoing_halfedges(last).collect();
            for h in fan {
                self.faces[h / 3][h % 3] = v as u32;
            }
            self.vertex_halfedge[v] = self.vertex_halfedge[last];
            vertex_move = Some((last as u32, v as u32));
        }
        self.vertex_halfedge.truncate(last);
        (dead_faces, face_moves, edge_moves, vertex_move)
    }

    // ----- merge (the inverse of an edge split) -----

    /// The star of vertex `v` (in range) as [`merge_edge`](Self::merge_edge) of the half `e`
    /// needs it, or why it is not one. Interior: the four outgoing halfedges `[u0, u1, u2, u3]`
    /// counter-clockwise with `e` the edge of `u0`; the halves are the edges of `u0` and `u2`.
    /// Boundary: `[u0, u1]`, `u0` the boundary halfedge leaving `v`; the halves are the edges of
    /// `u0` and of the halfedge before `u1` (the arriving boundary halfedge), `e` one of them.
    pub(crate) fn merge_star(&self, v: usize, e: usize) -> Result<Vec<usize>, RemovalRefusal> {
        let boundary = self.is_boundary_vertex(v);
        let out: Vec<usize> = self.outgoing_halfedges(v).take(5).collect();
        let wanted = if boundary { 2 } else { 4 };
        if out.len() != wanted {
            return Err(RemovalRefusal::MergeDegree {
                degree: self.vertex_degree(v),
                boundary,
            });
        }
        for &h in &out {
            if self.to_vertex(h) == v {
                return Err(RemovalRefusal::SelfEdge {
                    edge: self.halfedge_edge[h],
                });
            }
        }
        let star = if boundary {
            let halves = [self.halfedge_edge[out[0]], self.halfedge_edge[prev(out[1])]];
            if !halves.contains(&(e as u32)) {
                return Err(RemovalRefusal::NotAHalf { edge: e as u32 });
            }
            out
        } else {
            let Some(i) = out
                .iter()
                .position(|&h| self.halfedge_edge[h] as usize == e)
            else {
                return Err(RemovalRefusal::NotAHalf { edge: e as u32 });
            };
            (0..4).map(|k| out[(i + k) % 4]).collect()
        };
        // The spokes that disappear must not be constrained; the halves must agree.
        let n_spokes = if boundary { 1 } else { 2 };
        for s in 0..n_spokes {
            let spoke = self.halfedge_edge[star[2 * s + 1]];
            if self.constrained[spoke as usize] {
                return Err(RemovalRefusal::ConstrainedEdge { edge: spoke });
            }
        }
        let [h0, h1] = self.merge_halves(&star);
        if self.constrained[h0] != self.constrained[h1] {
            let edge = if self.constrained[h0] { h0 } else { h1 };
            return Err(RemovalRefusal::ConstrainedEdge { edge: edge as u32 });
        }
        Ok(star)
    }

    /// The two halves of a star of [`merge_star`](Self::merge_star): the edge of `u0` and the
    /// edge of the halfedge before `u1`.
    pub(crate) fn merge_halves(&self, star: &[usize]) -> [usize; 2] {
        [
            self.halfedge_edge[star[0]] as usize,
            self.halfedge_edge[prev(star[1])] as usize,
        ]
    }

    /// Merges the two halves of an edge at vertex `v` back into one edge and deletes `v`: the
    /// inverse of [`split_edge`](Self::split_edge) (SPEC N53). `v` is an interior vertex of
    /// degree four, `e` one of its edges and the other half the edge opposite it (two steps
    /// around `v`); or a boundary vertex of degree three, `e` one of its two boundary edges and
    /// the other half the other one. On each side of the merged edge the two faces at `v`
    /// become one triangle, and the edge from `v` to the opposite corner is deleted. Whether the
    /// halves are straight is for the representation to decide first
    /// ([`plan_edge_merge`](crate::merge::plan_edge_merge)).
    ///
    /// **Ids.** The merged edge takes the smaller of the two halves' ids and the direction of
    /// that half; on each side the merged face takes the smaller of the two face ids, with the
    /// corner `v` replaced by the far end of the other half (its other two corners keep their
    /// slots). Then, as in [`remove_vertex`](Self::remove_vertex), the freed ids (one vertex,
    /// three edges and two faces; two edges and one face on the boundary) are filled by the last
    /// vertex, the last edges and the last faces, and [`EdgeMerge`] lists every move. The
    /// halfedges that survive leave their tails in the same directions: the outer halfedges of
    /// the four (two) faces, and the halfedge from the far end of each half towards `v`, which
    /// becomes the merged edge's halfedge from that end; [`EdgeMerge::apply_halfedge_data`]
    /// carries their data. Merging the vertex of a [`split_edge`](Self::split_edge) right after
    /// the split gives back the complex exactly, ids included. Two constrained halves give a
    /// constrained edge.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`], and [`IntrinsicError::RemovalRefused`] with
    /// [`RemovalRefusal::MergeDegree`] (not degree four, or three on the boundary),
    /// [`RemovalRefusal::SelfEdge`], [`RemovalRefusal::NotAHalf`] (`e` is not an edge at `v`, or
    /// not a boundary edge at a boundary vertex), [`RemovalRefusal::ConstrainedEdge`] (an edge
    /// that would disappear is constrained, or only one half is). The complex is unchanged then.
    pub fn merge_edge(&mut self, v: usize, e: usize) -> Result<EdgeMerge, IntrinsicError> {
        self.check_index("vertex", v, self.n_vertices())?;
        self.check_index("edge", e, self.n_edges())?;
        let star = self
            .merge_star(v, e)
            .map_err(|reason| IntrinsicError::RemovalRefused {
                vertex: v as u32,
                reason,
            })?;
        Ok(self.merge_edge_unchecked(v, &star))
    }

    /// Removes boundary vertex `v` that lies on exactly two faces (SPEC N54; Liu et al. 2023,
    /// Sec. 4.2, inset right): with `b` the boundary neighbour before `v` (the boundary edge
    /// `b → v` arrives at `v`), `c` the one after it (the boundary edge `v → c` leaves it) and
    /// `a` the far end of the one interior edge at `v`, the faces `vab` and `vca` become the one
    /// face `abc`, the spoke `va` is deleted and the boundary edges `bv`, `vc` become the one
    /// boundary edge `bc`. The outer edges `ab`, `ca` are kept. The length of `bc` is the
    /// caller's business (`ℓ_bv + ℓ_vc` if the corner angles at `v` sum to `π`).
    ///
    /// This is [`merge_edge`](Self::merge_edge) of the boundary edge leaving `v`, with the
    /// edge chosen here: the merged edge keeps the smaller of the two boundary edges' ids, the
    /// merged face the smaller face id, and the freed ids (one vertex, two edges, one face) are
    /// filled by the last ones as in [`remove_vertex`](Self::remove_vertex); the returned
    /// [`EdgeMerge`] lists every move ([`faces`](EdgeMerge::faces)`[0]` is the merged face,
    /// [`edge`](EdgeMerge::edge) the boundary edge `bc`). Every surviving halfedge leaves its
    /// tail in the same direction as before; the merged boundary halfedge `b → c` takes the
    /// data of `b → v`. Removing the vertex of a boundary [`split_edge`](Self::split_edge)
    /// right after the split gives back the complex exactly, ids included.
    ///
    /// Not supported: an interior vertex, a vertex on another number of faces, a vertex with an
    /// edge to itself (which includes the only vertex of a boundary loop), a vertex whose spoke
    /// is constrained or one of whose two boundary edges is constrained and the other not.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`], and [`IntrinsicError::RemovalRefused`] with
    /// [`RemovalRefusal::InteriorVertex`], [`RemovalRefusal::BoundaryFaces`],
    /// [`RemovalRefusal::SelfEdge`] or [`RemovalRefusal::ConstrainedEdge`]; the complex is
    /// unchanged then.
    pub fn remove_boundary_vertex(&mut self, v: usize) -> Result<EdgeMerge, IntrinsicError> {
        self.check_index("vertex", v, self.n_vertices())?;
        let star =
            self.boundary_removal_star(v)
                .map_err(|reason| IntrinsicError::RemovalRefused {
                    vertex: v as u32,
                    reason,
                })?;
        Ok(self.merge_edge_unchecked(v, &star))
    }

    /// Why vertex `v` (in range) cannot be removed by
    /// [`remove_boundary_vertex`](Self::remove_boundary_vertex), if it cannot.
    pub fn boundary_removal_refusal(&self, v: usize) -> Option<RemovalRefusal> {
        self.boundary_removal_star(v).err()
    }

    fn boundary_removal_star(&self, v: usize) -> Result<Vec<usize>, RemovalRefusal> {
        let h0 = self.vertex_halfedge[v] as usize;
        if self.twin[h0] != NONE {
            return Err(RemovalRefusal::InteriorVertex);
        }
        let n_faces = self.outgoing_halfedges(v).take(3).count();
        if n_faces != 2 {
            return Err(RemovalRefusal::BoundaryFaces {
                faces: self.outgoing_halfedges(v).count(),
            });
        }
        self.merge_star(v, self.halfedge_edge[h0] as usize)
    }

    /// [`merge_edge`](Self::merge_edge) for a star that [`merge_star`](Self::merge_star)
    /// returned.
    pub(crate) fn merge_edge_unchecked(&mut self, v: usize, star: &[usize]) -> EdgeMerge {
        let interior = star.len() == 4;
        // Side s has the faces of a = u_{2s}, (v, X0, X1), and of b = u_{2s+1}, (v, X1, X2)
        // (X0..X2 named per side). They become (X0, X1, X2) with the halfedges next(a),
        // next(b) and prev(b); the last, X2 → v, becomes the merged edge's halfedge X2 → X0.
        let n_sides = if interior { 2 } else { 1 };
        let halves = self.merge_halves(star).map(|x| x as u32);
        let keep_edge = halves[0].min(halves[1]);
        let dead_half = halves[0].max(halves[1]);
        let mut dead_edges = vec![dead_half];
        for s in 0..n_sides {
            dead_edges.push(self.halfedge_edge[star[2 * s + 1]]);
        }
        // (old slot, new slot) of every surviving halfedge, three per side, the merged
        // halfedge last; the corner of the kept face that changes from v; the dead halfedge
        // X1 → v with the survivor X1 → X2 that replaces it as X1's outgoing halfedge.
        let mut survivors: Vec<(usize, usize)> = Vec::with_capacity(6);
        let mut kept_faces = [NONE; 2];
        let mut dead_faces = Vec::with_capacity(2);
        let mut corner_change = Vec::with_capacity(2);
        let mut apex_replacement = Vec::with_capacity(2);
        for s in 0..n_sides {
            let (a, b) = (star[2 * s], star[2 * s + 1]);
            let (fa, fb) = (a / 3, b / 3);
            let (s1, s2, s3) = (next(a), next(b), prev(b));
            let x0 = self.faces[s1 / 3][s1 % 3];
            let x2 = self.faces[s3 / 3][s3 % 3];
            if fa < fb {
                // Face fa = (v, X0, X1) keeps next(a) in its slot; v becomes X2.
                survivors.extend([(s1, s1), (s2, prev(a)), (s3, a)]);
                kept_faces[s] = fa as u32;
                dead_faces.push(fb);
                corner_change.push((a, x2));
            } else {
                // Face fb = (v, X1, X2) keeps next(b) and prev(b) in their slots; v becomes X0.
                survivors.extend([(s1, b), (s2, s2), (s3, s3)]);
                kept_faces[s] = fb as u32;
                dead_faces.push(fa);
                corner_change.push((b, x0));
            }
            apex_replacement.push((prev(a), s2));
        }
        let merged_slots: Vec<usize> = (0..n_sides).map(|s| survivors[3 * s + 2].1).collect();
        let star_faces: Vec<usize> = star.iter().map(|&h| h / 3).collect();
        let in_star = |h: u32| h != NONE && star_faces.contains(&(h as usize / 3));
        let renamed = |h: u32| -> u32 {
            match survivors.iter().find(|m| m.0 as u32 == h) {
                Some(&(_, n)) => n as u32,
                None => h,
            }
        };
        // Everything is read before anything is written.
        let old: Vec<(u32, u32, u32)> = survivors
            .iter()
            .map(|&(o, _)| {
                let e = self.halfedge_edge[o];
                (self.twin[o], e, self.edge_halfedge[e as usize])
            })
            .collect();
        // The merged edge keeps the direction of the kept half. The half of u0 (towards X0 of
        // side 0) directed away from v runs like side 0's merged halfedge X2 → X0, directed
        // towards v like side 1's; the half before u1 (towards X2 of side 0) the other way
        // round. On the boundary both halves run like side 0's.
        let keep_halfedge = self.edge_halfedge[keep_edge as usize] as usize;
        let away_from_v = self.faces[keep_halfedge / 3][keep_halfedge % 3] as usize == v;
        let side = if keep_edge == halves[0] {
            usize::from(!away_from_v)
        } else {
            usize::from(away_from_v)
        };
        let merged_direction = merged_slots[side.min(n_sides - 1)];
        let vertex_fix: Vec<(usize, u32)> = {
            let mut corners: Vec<usize> = survivors
                .iter()
                .map(|&(o, _)| self.faces[o / 3][o % 3] as usize)
                .collect();
            corners.sort_unstable();
            corners.dedup();
            corners
                .into_iter()
                .map(|x| {
                    let h = self.vertex_halfedge[x];
                    let h = match apex_replacement.iter().find(|m| m.0 as u32 == h) {
                        Some(&(_, r)) => r as u32,
                        None => h,
                    };
                    (x, renamed(h))
                })
                .collect()
        };

        for &(slot, corner) in &corner_change {
            self.faces[slot / 3][slot % 3] = corner;
        }
        for (&(_, n), &(t, e, eh)) in survivors.iter().zip(&old) {
            if merged_slots.contains(&n) {
                continue;
            }
            let t_new = if t == NONE { NONE } else { renamed(t) };
            self.twin[n] = t_new;
            self.halfedge_edge[n] = e;
            if t != NONE && !in_star(t) {
                self.twin[t as usize] = n as u32;
            }
            self.edge_halfedge[e as usize] = renamed(eh);
        }
        if interior {
            let [m0, m1] = [merged_slots[0], merged_slots[1]];
            self.twin[m0] = m1 as u32;
            self.twin[m1] = m0 as u32;
            self.halfedge_edge[m1] = keep_edge;
        } else {
            self.twin[merged_slots[0]] = NONE;
        }
        self.halfedge_edge[merged_slots[0]] = keep_edge;
        self.edge_halfedge[keep_edge as usize] = merged_direction as u32;
        for &(x, h) in &vertex_fix {
            self.vertex_halfedge[x] = h;
        }
        let mut halfedge_moves: Vec<(u32, u32)> = survivors
            .iter()
            .filter(|m| m.0 != m.1)
            .map(|&(o, n)| (o as u32, n as u32))
            .collect();
        let mut removed_edges = [NONE; 3];
        removed_edges[..dead_edges.len()].copy_from_slice(&dead_edges);
        removed_edges[1..].sort_unstable();
        let mut removed_faces = [NONE; 2];
        for (slot, &d) in removed_faces.iter_mut().zip(&dead_faces) {
            *slot = d as u32;
        }
        removed_faces.sort_unstable();
        let (_, face_moves, edge_moves, vertex_move) =
            self.compact(v, dead_faces, dead_edges, &mut halfedge_moves);
        let moved = |moves: &[(u32, u32)], x: u32| {
            if x == NONE {
                NONE
            } else {
                moves.iter().find(|m| m.0 == x).map_or(x, |m| m.1)
            }
        };
        let edge = moved(&edge_moves, keep_edge);
        let h = self.edge_halfedge[edge as usize] as usize;
        // Independent of which half was given: the face of the merged edge's halfedge first,
        // the moves by old id.
        let mut faces = kept_faces.map(|f| moved(&face_moves, f));
        if faces[0] as usize != h / 3 {
            faces.swap(0, 1);
        }
        halfedge_moves.sort_unstable();
        EdgeMerge {
            vertex: v as u32,
            edge,
            ends: [self.from_vertex(h) as u32, self.to_vertex(h) as u32],
            halves: [keep_edge, dead_half],
            faces,
            removed_edges,
            removed_faces,
            halfedge_moves,
            edge_moves,
            face_moves,
            vertex_move,
            counts: [
                self.n_halfedges() as u32,
                self.n_edges() as u32,
                self.n_faces() as u32,
                self.n_vertices() as u32,
            ],
        }
    }

    /// Moves face `from` with its halfedges to the free slot `to` and repairs what pointed at
    /// them; the halfedge moves are recorded.
    fn move_face(&mut self, from: usize, to: usize, moves: &mut Vec<(u32, u32)>) {
        self.faces[to] = self.faces[from];
        let renamed = |h: u32| {
            if h != NONE && h as usize / 3 == from {
                (3 * to + h as usize % 3) as u32
            } else {
                h
            }
        };
        for k in 0..3 {
            let (old, new) = (3 * from + k, 3 * to + k);
            let t = self.twin[old];
            self.twin[new] = renamed(t);
            if t != NONE && t as usize / 3 != from {
                self.twin[t as usize] = new as u32;
            }
            let e = self.halfedge_edge[old];
            self.halfedge_edge[new] = e;
            if self.edge_halfedge[e as usize] == old as u32 {
                self.edge_halfedge[e as usize] = new as u32;
            }
            let tail = self.faces[to][k] as usize;
            if self.vertex_halfedge[tail] == old as u32 {
                self.vertex_halfedge[tail] = new as u32;
            }
            record_move(moves, old as u32, new as u32);
        }
    }
}
