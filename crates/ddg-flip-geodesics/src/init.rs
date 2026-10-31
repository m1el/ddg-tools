//! Initial edge paths (SPEC §2.2): along a given vertex sequence, and the shortest path in the
//! edge graph of the current triangulation (Dijkstra, F15).

use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use ddg_intrinsic::{Complex, Triangulation};

use crate::error::FlipError;
use crate::network::{EdgePath, Segment};

/// One way to leave a vertex along an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Step {
    segment: Segment,
    head: usize,
    /// The tie key (F15): the id of the halfedge that runs like the step, or
    /// `n_halfedges + edge` for a boundary edge walked against its only halfedge.
    key: usize,
}

/// The steps out of vertex `v`: one per outgoing halfedge, counter-clockwise, and at a
/// boundary vertex the incoming boundary edge walked backwards.
fn steps(complex: &Complex, v: usize) -> impl Iterator<Item = Step> + '_ {
    let out = complex.outgoing_halfedges(v).map(move |h| Step {
        segment: Segment::new(complex.edge(h), complex.is_forward(h)),
        head: complex.to_vertex(h),
        key: h,
    });
    let back = complex.is_boundary_vertex(v).then(|| {
        let last = complex
            .outgoing_halfedges(v)
            .last()
            .expect("a vertex has an outgoing halfedge");
        let g = complex.prev(last);
        Step {
            segment: Segment::new(complex.edge(g), false),
            head: complex.from_vertex(g),
            key: complex.n_halfedges() + complex.edge(g),
        }
    });
    out.chain(back)
}

/// A distance ordered by `f64::total_cmp`.
#[derive(Debug, Clone, Copy)]
struct Dist(f64);

impl PartialEq for Dist {
    fn eq(&self, other: &Self) -> bool {
        self.0.total_cmp(&other.0) == Ordering::Equal
    }
}
impl Eq for Dist {}
impl PartialOrd for Dist {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Dist {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

fn check_vertex(complex: &Complex, v: usize) -> Result<(), FlipError> {
    if v < complex.n_vertices() {
        Ok(())
    } else {
        Err(FlipError::IndexOutOfRange {
            what: "vertex",
            index: v,
            bound: complex.n_vertices(),
        })
    }
}

/// The shortest edge path from `from` to `to` in the edge graph of the current triangulation,
/// weighted by the edge lengths (paper Sec. 5.1; `gc: mga.cpp:14-84`).
///
/// Dijkstra's algorithm over steps along edges: queue entries are `(distance, tie key)`,
/// smallest first, distances compared by `f64::total_cmp`; a vertex is settled by the first
/// entry popped for it. The tie key is the id of the halfedge the step runs along (the
/// reference ties by halfedge, `gc: mga.cpp:29-30`); a boundary edge walked against its only
/// halfedge has the key `n_halfedges + edge`, after every halfedge (F15).
///
/// # Errors
/// [`FlipError::IndexOutOfRange`]; [`FlipError::SameVertex`] if `from == to`;
/// [`FlipError::NotConnected`] if `to` is in another component (the reference returns an
/// empty path, `gc: fg.cpp:275-278`).
pub fn dijkstra_path<T: Triangulation>(
    tri: &T,
    from: usize,
    to: usize,
) -> Result<EdgePath, FlipError> {
    let complex = tri.complex();
    check_vertex(complex, from)?;
    check_vertex(complex, to)?;
    if from == to {
        return Err(FlipError::SameVertex { vertex: from });
    }
    let n = complex.n_vertices();
    // Per vertex: the step that settled it.
    let mut incoming: Vec<Option<Segment>> = vec![None; n];
    let mut settled = vec![false; n];
    settled[from] = true;
    let mut heap: BinaryHeap<Reverse<(Dist, usize, Segment, usize)>> = BinaryHeap::new();
    let push = |heap: &mut BinaryHeap<_>, settled: &[bool], v: usize, d: f64| {
        for s in steps(complex, v) {
            if !settled[s.head] {
                let l = tri.edge_length(s.segment.edge as usize);
                heap.push(Reverse((Dist(d + l), s.key, s.segment, s.head)));
            }
        }
    };
    push(&mut heap, &settled, from, 0.0);
    while let Some(Reverse((Dist(d), _, segment, v))) = heap.pop() {
        if settled[v] {
            continue;
        }
        settled[v] = true;
        incoming[v] = Some(segment);
        if v == to {
            let mut segments = Vec::new();
            let mut w = to;
            while w != from {
                let s = incoming[w].expect("every settled vertex but the start has a step");
                segments.push(s);
                w = s.tail(complex);
            }
            segments.reverse();
            return Ok(EdgePath::open(segments));
        }
        push(&mut heap, &settled, v, d);
    }
    Err(FlipError::NotConnected { from, to })
}

/// The edge path through the given vertices, consecutive ones joined by an edge. Where
/// several edges join two vertices (a Δ-complex may have double edges and self-edges), the
/// shortest is taken, ties by the key of [`dijkstra_path`]. A closed path joins the last
/// vertex back to the first; do not repeat the first vertex at the end. A closed path of one
/// vertex is a loop along a self-edge.
///
/// # Errors
/// [`FlipError::TooFewVertices`], [`FlipError::IndexOutOfRange`],
/// [`FlipError::NotAdjacent`].
pub fn path_from_vertices<T: Triangulation>(
    tri: &T,
    vertices: &[usize],
    closed: bool,
) -> Result<EdgePath, FlipError> {
    let complex = tri.complex();
    let n = vertices.len();
    if n == 0 || (!closed && n < 2) {
        return Err(FlipError::TooFewVertices { n, closed });
    }
    for &v in vertices {
        check_vertex(complex, v)?;
    }
    let n_steps = if closed { n } else { n - 1 };
    let mut segments = Vec::with_capacity(n_steps);
    for i in 0..n_steps {
        let (u, w) = (vertices[i], vertices[(i + 1) % n]);
        let best = steps(complex, u)
            .filter(|s| s.head == w)
            .min_by(|a, b| {
                let la = tri.edge_length(a.segment.edge as usize);
                let lb = tri.edge_length(b.segment.edge as usize);
                la.total_cmp(&lb).then(a.key.cmp(&b.key))
            })
            .ok_or(FlipError::NotAdjacent {
                index: i,
                from: u,
                to: w,
            })?;
        segments.push(best.segment);
    }
    Ok(EdgePath { segments, closed })
}
