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

/// Initial paths and the vertices to pin for them (SPEC §2.2), as returned by
/// [`waypoint_path`] and [`network_from_edges`]; pass both to
/// [`GeodesicNetwork::new`](crate::GeodesicNetwork::new).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InitialPaths {
    /// The paths, in the order they were found.
    pub paths: Vec<EdgePath>,
    /// The vertices to pin, ascending, without repeats. The ends of open paths are pinned by
    /// the network anyway; they are listed too, so that the list says where the network may
    /// not bend.
    pub pinned: Vec<usize>,
}

/// The edge path through the waypoints, each consecutive pair joined by the shortest edge path
/// ([`dijkstra_path`]), concatenated (`gc: fg.cpp:283-312`). A closed path joins the last
/// waypoint back to the first; do not repeat the first waypoint at the end.
///
/// With `pin_waypoints` every waypoint is pinned, so the path keeps passing through it (the
/// reference's `markInterior`, `gc: fg.cpp:293-296`); without, only the ends of an open path
/// are pinned and the waypoints only seed the path (FlipOut may move it off them).
///
/// The concatenation may revisit a vertex or run along an edge twice (the legs are
/// independent shortest paths); whether that is a valid network is decided by
/// [`GeodesicNetwork::new`](crate::GeodesicNetwork::new), which refuses crossings.
///
/// # Errors
/// [`FlipError::TooFewVertices`] for fewer than two waypoints (a closed path through one
/// vertex has no shortest path to itself); [`FlipError::WaypointLeg`] naming the leg and its
/// two vertices when a leg is out of range, joins a vertex to itself (consecutive equal
/// waypoints), or the two are not connected.
pub fn waypoint_path<T: Triangulation>(
    tri: &T,
    points: &[usize],
    closed: bool,
    pin_waypoints: bool,
) -> Result<InitialPaths, FlipError> {
    let n = points.len();
    if n < 2 {
        return Err(FlipError::TooFewVertices { n, closed });
    }
    let n_legs = if closed { n } else { n - 1 };
    let mut segments = Vec::new();
    for leg in 0..n_legs {
        let (from, to) = (points[leg], points[(leg + 1) % n]);
        let part = dijkstra_path(tri, from, to).map_err(|cause| FlipError::WaypointLeg {
            leg,
            from,
            to,
            cause: Box::new(cause),
        })?;
        segments.extend(part.segments);
    }
    let mut pinned = if pin_waypoints {
        points.to_vec()
    } else {
        Vec::new()
    };
    pinned.sort_unstable();
    pinned.dedup();
    Ok(InitialPaths {
        paths: vec![EdgePath { segments, closed }],
        pinned,
    })
}

/// The paths and loops formed by a set of edges (`gc: fg.cpp:315-396`).
///
/// A vertex is an **end** if it is in `extra_pinned` or if the number of masked edge ends at it
/// is not `0` and not `2` (a self-edge counts twice at its vertex). From every end, in the order
/// of the halfedge ids (then the boundary edges walked against their halfedge, by edge id),
/// every masked edge not yet walked starts a path that follows the masked edges through
/// vertices of degree two until it reaches an end (possibly the one it started from: an open
/// path from `x` to `x`). Masked edges left over form closed loops, in the same order. At a
/// degree-two vertex the path continues along the other masked edge, the first in the
/// counter-clockwise order of the outgoing halfedges.
///
/// The result pins every end that lies on a masked edge, so junctions keep their topology (Fig. 6,
/// 19, 22); a path ends nowhere else. Every masked edge appears in exactly one segment of exactly
/// one path.
///
/// # Errors
/// [`FlipError::MaskSize`] if the mask does not have one entry per edge;
/// [`FlipError::IndexOutOfRange`] for an `extra_pinned` vertex out of range;
/// [`FlipError::Internal`] if a walk does not end (a bug).
pub fn network_from_edges<T: Triangulation>(
    tri: &T,
    edge_mask: &[bool],
    extra_pinned: &[usize],
) -> Result<InitialPaths, FlipError> {
    let complex = tri.complex();
    let (n_edges, n_vertices) = (complex.n_edges(), complex.n_vertices());
    if edge_mask.len() != n_edges {
        return Err(FlipError::MaskSize {
            len: edge_mask.len(),
            n_edges,
        });
    }
    for &v in extra_pinned {
        check_vertex(complex, v)?;
    }
    let mut count = vec![0u32; n_vertices];
    for e in (0..n_edges).filter(|&e| edge_mask[e]) {
        let [a, b] = complex.edge_vertices(e);
        count[a] += 1;
        count[b] += 1;
    }
    let mut is_end: Vec<bool> = count.iter().map(|&c| c != 0 && c != 2).collect();
    for &v in extra_pinned {
        is_end[v] = true;
    }
    // Directed edges in the reference's order: halfedges by id, then boundary edges the other
    // way round.
    let directed = (0..complex.n_halfedges())
        .map(|h| Segment::new(complex.edge(h), complex.is_forward(h)))
        .chain(
            (0..n_edges)
                .filter(|&e| complex.is_boundary_edge(e))
                .map(|e| Segment::new(e, false)),
        )
        .filter(|s| edge_mask[s.edge as usize]);
    let directed: Vec<Segment> = directed.collect();
    let mut walked = vec![false; n_edges];
    let mut paths = Vec::new();
    let next_step = |current: Segment| -> Option<Segment> {
        let at = current.head(complex);
        steps(complex, at)
            .map(|s| s.segment)
            .find(|&s| s != current.reversed() && edge_mask[s.edge as usize])
    };
    let guard = |what: &str, start: Segment| FlipError::Internal {
        what: format!(
            "the walk from edge {} did not end ({what})",
            start.edge
        ),
    };
    for &start in &directed {
        if walked[start.edge as usize] || !is_end[start.tail(complex)] {
            continue;
        }
        let mut segments = Vec::new();
        let mut current = start;
        loop {
            walked[current.edge as usize] = true;
            segments.push(current);
            if is_end[current.head(complex)] {
                break;
            }
            current = next_step(current).ok_or_else(|| guard("no continuation", start))?;
            if segments.len() > n_edges {
                return Err(guard("too long", start));
            }
        }
        paths.push(EdgePath::open(segments));
    }
    for &start in &directed {
        if walked[start.edge as usize] {
            continue;
        }
        let mut segments = Vec::new();
        let mut current = start;
        loop {
            walked[current.edge as usize] = true;
            segments.push(current);
            current = next_step(current).ok_or_else(|| guard("no continuation", start))?;
            if current == start {
                break;
            }
            if segments.len() > n_edges {
                return Err(guard("too long", start));
            }
        }
        paths.push(EdgePath::closed(segments));
    }
    let pinned = (0..n_vertices)
        .filter(|&v| is_end[v] && count[v] != 0)
        .collect();
    Ok(InitialPaths { paths, pinned })
}
