//! What the tracing kernel needs from a triangulation. Implemented by the input mesh and by the
//! intrinsic triangulation, so that both are traced over by the same code.

use ddg_core::topology::NONE;

/// A triangulated surface with halfedges `h = 3 f + k`, edge lengths and signposts.
///
/// The halfedge `3 f + k` runs from corner `k` of face `f` to corner `k + 1`; `next` and `prev`
/// are arithmetic. A triangle may use a vertex or an edge more than once.
pub(crate) trait Surface {
    /// Name used in error messages.
    const NAME: &'static str;
    fn n_vertices(&self) -> usize;
    fn n_edges(&self) -> usize;
    fn n_faces(&self) -> usize;
    /// The opposite halfedge, or [`NONE`] on the boundary.
    fn twin(&self, h: usize) -> u32;
    /// The vertex a halfedge starts at.
    fn tail(&self, h: usize) -> u32;
    /// The edge of a halfedge.
    fn edge(&self, h: usize) -> u32;
    /// The length of an edge.
    fn length(&self, e: usize) -> f64;
    /// An outgoing halfedge of a vertex; the boundary one at a boundary vertex.
    fn outgoing(&self, v: usize) -> u32;
    /// The normalized direction of a halfedge at its tail.
    fn signpost(&self, h: usize) -> f64;
    /// The sum of the corner angles at a vertex.
    fn angle_sum(&self, v: usize) -> f64;
    /// Whether a halfedge runs in the direction of its edge (the direction edge parameters and
    /// edge frames refer to).
    fn is_forward(&self, h: usize) -> bool;
    /// A halfedge of an edge (the only one for a boundary edge).
    fn edge_halfedge(&self, e: usize) -> u32;
}

#[inline]
pub(crate) fn next(h: usize) -> usize {
    if h % 3 == 2 { h - 2 } else { h + 1 }
}

#[inline]
pub(crate) fn prev(h: usize) -> usize {
    if h.is_multiple_of(3) { h + 2 } else { h - 1 }
}

/// The lengths of the three halfedges of a face, by slot.
#[inline]
pub(crate) fn face_lengths<S: Surface>(s: &S, f: usize) -> [f64; 3] {
    [
        s.length(s.edge(3 * f) as usize),
        s.length(s.edge(3 * f + 1) as usize),
        s.length(s.edge(3 * f + 2) as usize),
    ]
}

#[inline]
pub(crate) fn is_boundary_vertex<S: Surface>(s: &S, v: usize) -> bool {
    s.twin(s.outgoing(v) as usize) == NONE
}

/// The parameter of a point of a halfedge in the direction of its edge.
#[inline]
pub(crate) fn edge_parameter<S: Surface>(s: &S, h: usize, t_from_tail: f64) -> f64 {
    if s.is_forward(h) {
        t_from_tail
    } else {
        1.0 - t_from_tail
    }
}
