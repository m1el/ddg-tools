//! Edge networks: paths and loops as sequences of segments on the edges of an intrinsic
//! triangulation, the layering of segments that share an edge, pinned vertices, and their
//! validity (SPEC §2.1).
//!
//! A [`Network`] is pure combinatorics: it stores edge ids and refers to the connectivity
//! ([`Complex`]) only to know where an edge starts and ends. [`GeodesicNetwork`] owns a
//! triangulation `T` together with a network on it.
//!
//! # Conventions
//!
//! * A **segment** is a directed edge `(edge, forward)`, `forward` relative to
//!   [`Complex::edge_halfedge`] (F8). It is stored by edge id, which flips keep; an edge that
//!   is not flipped keeps its direction (SPEC N50), and edges that carry segments are never
//!   flipped (F7).
//! * Every segment of the network has an id ([`SegmentId`]) that never changes and is never
//!   reused. The segments of the paths given to [`Network::new`] get the ids `0, 1, …` in path
//!   order.
//! * The **layering** of edge `e` lists the segments on `e` ordered across the edge, from the
//!   face of `edge_halfedge(e)` to the face of its twin (the reference's deque,
//!   `gc: fg.h:133`). The segment *outermost on the side of halfedge `h`* is the first entry if
//!   `h` is `edge_halfedge(e)`, the last otherwise (`gc: fg.cpp:1319-1351`). On a boundary edge
//!   the "other side" is outside the surface; the order is kept all the same.
//! * A **joint** is a pair of consecutive segments `(ab, bc)`; it is named by its two segment
//!   ids and lies at the tail `b` of the outgoing one. A closed path also has the joint
//!   last → first.

use std::collections::VecDeque;

use ddg_intrinsic::{Complex, Insertion, Triangulation};

use crate::certify::{CertifyOptions, GeodesicCertificate};
use crate::error::FlipError;
use crate::wedge::{JointGeometry, VertexStar};

/// The "no segment" marker of the linked lists.
pub(crate) const NONE: u32 = u32::MAX;

/// A directed edge of the triangulation: `forward` means in the direction of
/// [`Complex::edge_halfedge`] of `edge` (F8).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Segment {
    /// The edge.
    pub edge: u32,
    /// Whether the segment runs from the tail to the head of the edge's halfedge.
    pub forward: bool,
}

impl Segment {
    /// The segment on edge `e`, forward or backward.
    pub fn new(edge: usize, forward: bool) -> Segment {
        Segment {
            edge: edge as u32,
            forward,
        }
    }

    /// The same edge in the other direction.
    pub fn reversed(self) -> Segment {
        Segment {
            edge: self.edge,
            forward: !self.forward,
        }
    }

    /// The vertex the segment starts at.
    pub fn tail(self, complex: &Complex) -> usize {
        let [a, b] = complex.edge_vertices(self.edge as usize);
        if self.forward { a } else { b }
    }

    /// The vertex the segment ends at.
    pub fn head(self, complex: &Complex) -> usize {
        let [a, b] = complex.edge_vertices(self.edge as usize);
        if self.forward { b } else { a }
    }

    /// The halfedge that runs like the segment: `edge_halfedge` if forward, its twin
    /// otherwise; `None` for a backward segment on a boundary edge.
    pub fn halfedge(self, complex: &Complex) -> Option<usize> {
        let h = complex.edge_halfedge(self.edge as usize);
        if self.forward {
            Some(h)
        } else {
            complex.twin(h)
        }
    }
}

/// A path or loop given as its segments, in order (F9: `closed` is explicit).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct EdgePath {
    /// The segments; consecutive ones share a vertex.
    pub segments: Vec<Segment>,
    /// Whether the last segment is followed by the first (a loop). An open path may start
    /// and end at the same vertex (a loop through a pinned point).
    pub closed: bool,
}

impl EdgePath {
    /// An open path.
    pub fn open(segments: Vec<Segment>) -> EdgePath {
        EdgePath {
            segments,
            closed: false,
        }
    }

    /// A closed path.
    pub fn closed(segments: Vec<Segment>) -> EdgePath {
        EdgePath {
            segments,
            closed: true,
        }
    }

    /// The vertices: the tail of every segment, then the head of the last one for an open
    /// path (`n + 1` vertices for `n` segments; `n` for a closed path).
    pub fn vertices(&self, complex: &Complex) -> Vec<usize> {
        let mut v: Vec<usize> = self.segments.iter().map(|s| s.tail(complex)).collect();
        if !self.closed
            && let Some(last) = self.segments.last()
        {
            v.push(last.head(complex));
        }
        v
    }

    /// The length: the sum of the edge lengths in path order.
    pub fn length<T: Triangulation>(&self, tri: &T) -> f64 {
        self.segments
            .iter()
            .map(|s| tri.edge_length(s.edge as usize))
            .fold(0.0, |s, l| s + l)
    }

    /// The same path traversed backwards.
    pub fn reversed(&self) -> EdgePath {
        EdgePath {
            segments: self.segments.iter().rev().map(|s| s.reversed()).collect(),
            closed: self.closed,
        }
    }
}

/// The id of a segment of a network. Ids are never reused.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SegmentId(pub u32);

/// A segment named by its path and its position in the path, for explicit layerings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SegmentRef {
    /// The path.
    pub path: usize,
    /// The position of the segment in the path.
    pub index: usize,
}

/// An explicit layering: per edge id, the segments on the edge ordered from the face of
/// `edge_halfedge(e)` to the face of its twin. Every segment of every path appears exactly
/// once, in the list of its own edge.
pub type Layering = Vec<Vec<SegmentRef>>;

/// What the network knows about one segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SegmentInfo {
    /// The path it belongs to.
    pub path: usize,
    /// The directed edge.
    pub segment: Segment,
    /// The segment before it in its path (`None` at the start of an open path).
    pub prev: Option<SegmentId>,
    /// The segment after it in its path (`None` at the end of an open path).
    pub next: Option<SegmentId>,
}

impl SegmentInfo {
    /// Whether this is the first or last segment of an open path (the reference's
    /// `isEndpoint`, `gc: fg.cpp:1700-1705`).
    pub fn is_path_end(&self) -> bool {
        self.prev.is_none() || self.next.is_none()
    }
}

/// A joint `γ_abc = (ab, bc)`: two consecutive segments of a path. It lies at the tail of
/// `outgoing`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Joint {
    /// The path.
    pub path: usize,
    /// The segment `ab` that arrives at `b`.
    pub incoming: SegmentId,
    /// The segment `bc` that leaves `b`.
    pub outgoing: SegmentId,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Node {
    pub(crate) path: u32,
    pub(crate) segment: Segment,
    pub(crate) prev: u32,
    pub(crate) next: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PathData {
    pub(crate) closed: bool,
    /// The first segment (any segment of a loop); `NONE` once the path has contracted.
    pub(crate) first: u32,
    pub(crate) n_segments: u32,
}

/// Paths and loops on the edges of a triangulation, with the layering of every edge and the
/// pinned vertices (SPEC §2.1). See the [module documentation](self).
#[derive(Debug, Clone, PartialEq)]
pub struct Network {
    /// By segment id; `None` for a removed segment.
    pub(crate) nodes: Vec<Option<Node>>,
    pub(crate) paths: Vec<PathData>,
    /// Per edge: segment ids from the face of `edge_halfedge` to the face of its twin.
    pub(crate) layers: Vec<VecDeque<u32>>,
    pub(crate) pinned: Vec<bool>,
}

impl Network {
    /// A network of the given paths on the triangulation with connectivity `complex`.
    ///
    /// `pinned` lists vertices where no joint is ever shortened (F10); the two ends of every
    /// open path are pinned as well. Without a `layering` each segment, in path order, is put
    /// outermost on the side of its own halfedge (`gc: fg.cpp:29-32`, F11): forward segments
    /// at the front of their edge's list, backward ones at the back. Then the network is
    /// checked (paths connected, loops closed, the layering complete, no crossing: see
    /// [`check`](Self::check)). Inputs that cross are refused, never repaired.
    ///
    /// # Errors
    /// [`FlipError::IndexOutOfRange`] for an edge, vertex or layering entry out of range,
    /// [`FlipError::EmptyPath`], [`FlipError::PathNotConnected`],
    /// [`FlipError::LoopNotClosed`], the `Layering*` errors for a layering that does not list
    /// every segment once on its own edge, [`FlipError::Crossing`].
    pub fn new(
        complex: &Complex,
        paths: &[EdgePath],
        pinned: &[usize],
        layering: Option<&Layering>,
    ) -> Result<Network, FlipError> {
        let n_edges = complex.n_edges();
        let n_vertices = complex.n_vertices();
        for (p, path) in paths.iter().enumerate() {
            if path.segments.is_empty() {
                return Err(FlipError::EmptyPath { path: p });
            }
            for s in &path.segments {
                if s.edge as usize >= n_edges {
                    return Err(FlipError::IndexOutOfRange {
                        what: "edge",
                        index: s.edge as usize,
                        bound: n_edges,
                    });
                }
            }
            check_connected(complex, p, path)?;
        }
        for &v in pinned {
            if v >= n_vertices {
                return Err(FlipError::IndexOutOfRange {
                    what: "vertex",
                    index: v,
                    bound: n_vertices,
                });
            }
        }

        // Segment ids in path order.
        let mut nodes = Vec::new();
        let mut path_data = Vec::with_capacity(paths.len());
        let mut first_id = Vec::with_capacity(paths.len());
        for (p, path) in paths.iter().enumerate() {
            let base = nodes.len() as u32;
            let n = path.segments.len() as u32;
            first_id.push(base);
            for (i, &segment) in path.segments.iter().enumerate() {
                let i = i as u32;
                let prev = if i > 0 {
                    base + i - 1
                } else if path.closed {
                    base + n - 1
                } else {
                    NONE
                };
                let next = if i + 1 < n {
                    base + i + 1
                } else if path.closed {
                    base
                } else {
                    NONE
                };
                nodes.push(Some(Node {
                    path: p as u32,
                    segment,
                    prev,
                    next,
                }));
            }
            path_data.push(PathData {
                closed: path.closed,
                first: base,
                n_segments: n,
            });
        }

        let mut layers: Vec<VecDeque<u32>> = vec![VecDeque::new(); n_edges];
        match layering {
            None => {
                for (id, node) in nodes.iter().enumerate() {
                    let node = node.expect("all nodes are live after construction");
                    let layer = &mut layers[node.segment.edge as usize];
                    if node.segment.forward {
                        layer.push_front(id as u32);
                    } else {
                        layer.push_back(id as u32);
                    }
                }
            }
            Some(layering) => {
                if layering.len() != n_edges {
                    return Err(FlipError::LayeringSize {
                        len: layering.len(),
                        n_edges,
                    });
                }
                let mut seen = vec![false; nodes.len()];
                for (e, list) in layering.iter().enumerate() {
                    for r in list {
                        if r.path >= paths.len() {
                            return Err(FlipError::IndexOutOfRange {
                                what: "path",
                                index: r.path,
                                bound: paths.len(),
                            });
                        }
                        let len = paths[r.path].segments.len();
                        if r.index >= len {
                            return Err(FlipError::IndexOutOfRange {
                                what: "segment of the path",
                                index: r.index,
                                bound: len,
                            });
                        }
                        let segment_edge = paths[r.path].segments[r.index].edge as usize;
                        if segment_edge != e {
                            return Err(FlipError::LayeringWrongEdge {
                                edge: e,
                                path: r.path,
                                index: r.index,
                                segment_edge,
                            });
                        }
                        let id = first_id[r.path] + r.index as u32;
                        if std::mem::replace(&mut seen[id as usize], true) {
                            return Err(FlipError::LayeringDuplicate {
                                edge: e,
                                path: r.path,
                                index: r.index,
                            });
                        }
                        layers[e].push_back(id);
                    }
                }
                for (p, path) in paths.iter().enumerate() {
                    for (i, s) in path.segments.iter().enumerate() {
                        if !seen[(first_id[p] + i as u32) as usize] {
                            return Err(FlipError::LayeringMissing {
                                edge: s.edge as usize,
                                path: p,
                                index: i,
                            });
                        }
                    }
                }
            }
        }

        let mut mask = vec![false; n_vertices];
        for &v in pinned {
            mask[v] = true;
        }
        for path in paths.iter().filter(|p| !p.closed) {
            mask[path.segments[0].tail(complex)] = true;
            mask[path.segments[path.segments.len() - 1].head(complex)] = true;
        }

        let network = Network {
            nodes,
            paths: path_data,
            layers,
            pinned: mask,
        };
        network.check(complex)?;
        Ok(network)
    }

    // ----- paths -----

    /// Number of paths (including paths that have contracted to nothing, which keep their
    /// index).
    pub fn n_paths(&self) -> usize {
        self.paths.len()
    }

    /// Whether path `p` is closed.
    ///
    /// # Panics
    /// If `p` is out of range.
    pub fn is_closed(&self, p: usize) -> bool {
        self.paths[p].closed
    }

    /// Number of segments of path `p` (`0` once it has contracted).
    ///
    /// # Panics
    /// If `p` is out of range.
    pub fn n_path_segments(&self, p: usize) -> usize {
        self.paths[p].n_segments as usize
    }

    /// Total number of segments.
    pub fn n_segments(&self) -> usize {
        self.paths.iter().map(|p| p.n_segments as usize).sum()
    }

    /// The segment ids of path `p` in order, from its first segment.
    ///
    /// # Panics
    /// If `p` is out of range.
    pub fn path_segment_ids(&self, p: usize) -> Vec<SegmentId> {
        let data = self.paths[p];
        let mut out = Vec::with_capacity(data.n_segments as usize);
        let mut id = data.first;
        while id != NONE && out.len() < data.n_segments as usize {
            out.push(SegmentId(id));
            id = self.node(id).next;
            if id == data.first {
                break;
            }
        }
        out
    }

    /// Path `p` as an [`EdgePath`].
    ///
    /// # Panics
    /// If `p` is out of range.
    pub fn path(&self, p: usize) -> EdgePath {
        EdgePath {
            segments: self
                .path_segment_ids(p)
                .into_iter()
                .map(|id| self.node(id.0).segment)
                .collect(),
            closed: self.paths[p].closed,
        }
    }

    /// All paths as [`EdgePath`]s, in path order.
    pub fn paths(&self) -> Vec<EdgePath> {
        (0..self.paths.len()).map(|p| self.path(p)).collect()
    }

    /// The joints of path `p` in order: for an open path at its interior vertices, for a
    /// closed path at every vertex, starting with the joint last → first.
    ///
    /// # Panics
    /// If `p` is out of range.
    pub fn path_joints(&self, p: usize) -> Vec<Joint> {
        self.path_segment_ids(p)
            .into_iter()
            .filter_map(|id| {
                let prev = self.node(id.0).prev;
                (prev != NONE).then_some(Joint {
                    path: p,
                    incoming: SegmentId(prev),
                    outgoing: id,
                })
            })
            .collect()
    }

    /// All joints, path by path.
    pub fn joints(&self) -> Vec<Joint> {
        (0..self.paths.len())
            .flat_map(|p| self.path_joints(p))
            .collect()
    }

    /// The vertex of a joint: the tail of its outgoing segment.
    ///
    /// # Panics
    /// If the joint's segments are not live segments of this network.
    pub fn joint_vertex(&self, complex: &Complex, joint: Joint) -> usize {
        self.node(joint.outgoing.0).segment.tail(complex)
    }

    // ----- segments -----

    /// Segment `id`, or `None` if there is no such segment (any more).
    pub fn segment(&self, id: SegmentId) -> Option<SegmentInfo> {
        let node = (*self.nodes.get(id.0 as usize)?)?;
        let opt = |x: u32| (x != NONE).then_some(SegmentId(x));
        Some(SegmentInfo {
            path: node.path as usize,
            segment: node.segment,
            prev: opt(node.prev),
            next: opt(node.next),
        })
    }

    /// One more than the largest segment id ever given out.
    pub fn segment_id_bound(&self) -> usize {
        self.nodes.len()
    }

    /// Whether segment `id` is the first or last segment of an open path.
    ///
    /// # Panics
    /// If `id` is not a live segment.
    pub fn is_path_end(&self, id: SegmentId) -> bool {
        let n = self.node(id.0);
        n.prev == NONE || n.next == NONE
    }

    pub(crate) fn node(&self, id: u32) -> Node {
        self.nodes[id as usize].expect("a live segment id")
    }

    // ----- layering -----

    /// The segments on edge `e`, ordered from the face of `edge_halfedge(e)` to the face of
    /// its twin.
    ///
    /// # Panics
    /// If `e` is out of range.
    pub fn layer(&self, e: usize) -> impl ExactSizeIterator<Item = SegmentId> + '_ {
        self.layers[e].iter().map(|&id| SegmentId(id))
    }

    /// Whether edge `e` carries a segment. Such an edge is never flipped (F7).
    ///
    /// # Panics
    /// If `e` is out of range.
    pub fn carries_segment(&self, e: usize) -> bool {
        !self.layers[e].is_empty()
    }

    /// The segment outermost on the side of halfedge `h` of its edge (`gc: fg.cpp:1319-1332`):
    /// the first of the edge's layer if `h` is its `edge_halfedge`, the last otherwise.
    /// `None` if the edge carries no segment.
    pub fn outermost(&self, complex: &Complex, h: usize) -> Option<SegmentId> {
        let e = complex.edge(h);
        let layer = &self.layers[e];
        let id = if complex.edge_halfedge(e) == h {
            layer.front()
        } else {
            layer.back()
        };
        id.map(|&id| SegmentId(id))
    }

    /// The current layering in the form [`Network::new`] takes, with segment positions
    /// counted from each path's first segment. `Network::new(complex, &self.paths(), …,
    /// Some(&self.layering()))` rebuilds this network up to segment ids.
    pub fn layering(&self) -> Layering {
        let mut position = vec![usize::MAX; self.nodes.len()];
        for p in 0..self.paths.len() {
            for (i, id) in self.path_segment_ids(p).into_iter().enumerate() {
                position[id.0 as usize] = i;
            }
        }
        self.layers
            .iter()
            .map(|layer| {
                layer
                    .iter()
                    .map(|&id| SegmentRef {
                        path: self.node(id).path as usize,
                        index: position[id as usize],
                    })
                    .collect()
            })
            .collect()
    }

    // ----- pinning -----

    /// Whether vertex `v` is pinned: no joint there is ever shortened (F10).
    ///
    /// # Panics
    /// If `v` is out of range.
    pub fn is_pinned(&self, v: usize) -> bool {
        self.pinned[v]
    }

    /// The pinned vertices, ascending.
    pub fn pinned_vertices(&self) -> Vec<usize> {
        (0..self.pinned.len()).filter(|&v| self.pinned[v]).collect()
    }

    /// Pins vertex `v`. Vertices added to the triangulation after the network was made can be
    /// pinned too (the mask grows).
    pub fn pin(&mut self, v: usize) {
        if v >= self.pinned.len() {
            self.pinned.resize(v + 1, false);
        }
        self.pinned[v] = true;
    }

    /// Unpins vertex `v`.
    ///
    /// # Errors
    /// [`FlipError::UnpinnedEnd`] if `v` is an end of an open path (ends stay pinned);
    /// nothing is changed then.
    pub fn unpin(&mut self, complex: &Complex, v: usize) -> Result<(), FlipError> {
        for p in 0..self.paths.len() {
            if let Some([a, b]) = self.path_ends(complex, p)
                && (a == v || b == v)
            {
                return Err(FlipError::UnpinnedEnd { path: p, vertex: v });
            }
        }
        if v < self.pinned.len() {
            self.pinned[v] = false;
        }
        Ok(())
    }

    /// The two end vertices of open path `p`; `None` for a closed or contracted path.
    ///
    /// # Panics
    /// If `p` is out of range.
    pub fn path_ends(&self, complex: &Complex, p: usize) -> Option<[usize; 2]> {
        let data = self.paths[p];
        if data.closed || data.first == NONE {
            return None;
        }
        let first = self.node(data.first);
        let mut last = data.first;
        while self.node(last).next != NONE {
            last = self.node(last).next;
        }
        Some([
            first.segment.tail(complex),
            self.node(last).segment.head(complex),
        ])
    }

    // ----- validity -----

    /// Checks the network against the connectivity `complex` (SPEC §2.1): every live segment
    /// lies on an existing edge; the linked lists of the paths are consistent and connected,
    /// loops are closed, open paths have two ends and both are pinned; every segment appears
    /// exactly once in the layering of its own edge; and nothing crosses: at every vertex, the
    /// segment ends listed in the cyclic order of the edges around it, each edge's segments
    /// in layer order, pair up into joints none of which interleave (path ends pair nothing).
    ///
    /// # Errors
    /// [`FlipError::Crossing`], [`FlipError::UnpinnedEnd`], [`FlipError::PathNotConnected`],
    /// [`FlipError::LoopNotClosed`], or [`FlipError::Internal`] for broken bookkeeping.
    pub fn check(&self, complex: &Complex) -> Result<(), FlipError> {
        let internal = |what: String| Err(FlipError::Internal { what });
        let n_edges = complex.n_edges();
        if self.layers.len() != n_edges {
            return internal(format!(
                "the network has layers for {} edges, the triangulation has {n_edges}",
                self.layers.len()
            ));
        }
        if self.pinned.len() < complex.n_vertices() {
            return internal(format!(
                "the pinned mask covers {} vertices, the triangulation has {}",
                self.pinned.len(),
                complex.n_vertices()
            ));
        }
        // Paths: linked lists, membership, counts.
        let mut owner = vec![NONE; self.nodes.len()];
        for (p, data) in self.paths.iter().enumerate() {
            if data.first == NONE {
                if data.n_segments != 0 {
                    return internal(format!("path {p} has no first segment but a count"));
                }
                continue;
            }
            let ids = self.walk_checked(p)?;
            for &id in &ids {
                owner[id as usize] = p as u32;
            }
            let path = EdgePath {
                segments: ids.iter().map(|&id| self.node(id).segment).collect(),
                closed: data.closed,
            };
            for s in &path.segments {
                if s.edge as usize >= n_edges {
                    return internal(format!(
                        "path {p} has a segment on edge {}, out of range ({n_edges} edges)",
                        s.edge
                    ));
                }
            }
            check_connected(complex, p, &path)?;
            if !data.closed {
                for v in [
                    path.segments[0].tail(complex),
                    path.segments[path.segments.len() - 1].head(complex),
                ] {
                    if !self.pinned[v] {
                        return Err(FlipError::UnpinnedEnd { path: p, vertex: v });
                    }
                }
            }
        }
        for (id, node) in self.nodes.iter().enumerate() {
            if let Some(node) = node
                && owner[id] != node.path
            {
                return internal(format!(
                    "segment {id} claims path {} but is not on its list",
                    node.path
                ));
            }
        }
        // Layering: every live segment once, on its own edge.
        let mut seen = vec![false; self.nodes.len()];
        for (e, layer) in self.layers.iter().enumerate() {
            for &id in layer {
                let Some(Some(node)) = self.nodes.get(id as usize) else {
                    return internal(format!("edge {e} lists segment {id}, which does not exist"));
                };
                if node.segment.edge as usize != e {
                    return internal(format!(
                        "edge {e} lists segment {id}, which lies on edge {}",
                        node.segment.edge
                    ));
                }
                if std::mem::replace(&mut seen[id as usize], true) {
                    return internal(format!("edge {e} lists segment {id} twice"));
                }
            }
        }
        for (id, node) in self.nodes.iter().enumerate() {
            if let Some(node) = node
                && !seen[id]
            {
                return internal(format!(
                    "segment {id} (path {}) is missing from the layering of edge {}",
                    node.path, node.segment.edge
                ));
            }
        }
        // Non-crossing, at every vertex with a joint.
        let mut vertices: Vec<usize> = self
            .nodes
            .iter()
            .flatten()
            .filter(|n| n.prev != NONE)
            .map(|n| n.segment.tail(complex))
            .collect();
        vertices.sort_unstable();
        vertices.dedup();
        for v in vertices {
            VertexStar::new(complex, self, v).check_non_crossing(self)?;
        }
        Ok(())
    }

    /// The segment ids of path `p` in order, checking the links.
    fn walk_checked(&self, p: usize) -> Result<Vec<u32>, FlipError> {
        let data = self.paths[p];
        let internal = |what: String| Err(FlipError::Internal { what });
        let n = data.n_segments as usize;
        let mut ids = Vec::with_capacity(n);
        let mut id = data.first;
        let mut prev = NONE;
        if !data.closed {
            match self.nodes.get(id as usize).copied().flatten() {
                Some(node) if node.prev == NONE => {}
                _ => {
                    return internal(format!(
                        "path {p} is open but its first segment {id} has a predecessor or does not exist"
                    ));
                }
            }
        }
        loop {
            let Some(Some(node)) = self.nodes.get(id as usize) else {
                return internal(format!(
                    "path {p} links to segment {id}, which does not exist"
                ));
            };
            if node.path as usize != p {
                return internal(format!(
                    "segment {id} is on path {p} but claims path {}",
                    node.path
                ));
            }
            if prev != NONE && node.prev != prev {
                return internal(format!(
                    "segment {id} of path {p} follows {prev}, but its predecessor is {}",
                    node.prev
                ));
            }
            ids.push(id);
            if ids.len() > n {
                return internal(format!(
                    "path {p} has more than its {n} segments (a cycle in an open path?)"
                ));
            }
            prev = id;
            id = node.next;
            if id == NONE {
                if data.closed {
                    return internal(format!(
                        "path {p} is closed but segment {prev} has no successor"
                    ));
                }
                break;
            }
            if id == data.first {
                if !data.closed {
                    return internal(format!("path {p} is open but its segments form a cycle"));
                }
                if self.node(id).prev != prev {
                    return internal(format!(
                        "loop {p}: the first segment {id} does not follow the last {prev}"
                    ));
                }
                break;
            }
        }
        if ids.len() != n {
            return internal(format!(
                "path {p} has {} segments, its count says {n}",
                ids.len()
            ));
        }
        Ok(ids)
    }
}

/// What [`Network::replace_run`] did.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Replaced {
    /// The ids of the new segments, in path order.
    pub(crate) new_ids: Vec<u32>,
    /// The segment that now follows the new ones (when nothing was added: the one that
    /// followed the removed run), if any: the joint at its tail changed.
    pub(crate) after: Option<u32>,
    /// The segments that were next to a removed one in the layer of its edge: the wedges at
    /// their ends may have become flexible.
    pub(crate) layer_neighbours: Vec<u32>,
    /// Whether the path has no segment left (it contracted, F12).
    pub(crate) contracted: bool,
}

impl Network {
    /// Replaces a run of one or two consecutive segments of one path by `new` (in path
    /// order), writing what it did into `out` (whose buffers are reused), each new segment
    /// pushed outermost on the side of the face of the halfedge given with it (SPEC §4.1 step
    /// 2: the wedge side; §4.3: the side of the face crossed).
    ///
    /// The removed segments are taken out of their layers wherever they are (Alg. 1 removes
    /// segments outermost on the wedge side, the contraction of §4.1 step 3 two adjacent
    /// ones). When the run is the whole of a closed path the new segments form the loop by
    /// themselves. When nothing is added the neighbours of the run are joined; a path left
    /// without segments has contracted (`first = NONE`, F12).
    ///
    /// # Errors
    /// [`FlipError::Internal`] if `run` is not one or two consecutive live segments of one
    /// path, or a removed segment is missing from its layer.
    pub(crate) fn replace_run(
        &mut self,
        complex: &Complex,
        run: &[u32],
        new: &[(Segment, usize)],
        out: &mut Replaced,
    ) -> Result<(), FlipError> {
        out.new_ids.clear();
        out.layer_neighbours.clear();
        let internal = |what: String| Err(FlipError::Internal { what });
        let live = |id: u32| self.nodes.get(id as usize).copied().flatten();
        let Some(first) = run.first().and_then(|&id| live(id)) else {
            return internal(format!("replace_run: bad run {run:?}"));
        };
        let p = first.path as usize;
        if run.len() > 2
            || (run.len() == 2
                && (first.next != run[1] || live(run[1]).is_none_or(|n| n.prev != run[0])))
        {
            return internal(format!(
                "replace_run: {run:?} are not consecutive segments of path {p}"
            ));
        }
        let last_id = *run.last().expect("non-empty");
        let last = live(last_id).expect("checked");
        let data = self.paths[p];
        let whole = data.closed && data.n_segments as usize == run.len();
        let (before, after) = if whole {
            (NONE, NONE)
        } else {
            (first.prev, last.next)
        };

        // Out of the layers.
        let layer_neighbours = &mut out.layer_neighbours;
        for &id in run {
            let e = live(id).expect("checked").segment.edge as usize;
            let layer = &mut self.layers[e];
            let Some(pos) = layer.iter().position(|&x| x == id) else {
                return internal(format!(
                    "segment {id} of path {p} is not in the layer of its edge {e}"
                ));
            };
            layer.remove(pos);
            for q in [pos.checked_sub(1), Some(pos)].into_iter().flatten() {
                if let Some(&x) = layer.get(q)
                    && !run.contains(&x)
                {
                    layer_neighbours.push(x);
                }
            }
        }
        for &id in run {
            self.nodes[id as usize] = None;
        }

        // The new segments.
        let new_ids = &mut out.new_ids;
        for &(segment, h) in new {
            let id = self.nodes.len() as u32;
            self.nodes.push(Some(Node {
                path: p as u32,
                segment,
                prev: NONE,
                next: NONE,
            }));
            let e = segment.edge as usize;
            if complex.edge_halfedge(e) == h {
                self.layers[e].push_front(id);
            } else {
                self.layers[e].push_back(id);
            }
            new_ids.push(id);
        }
        for w in new_ids.windows(2) {
            Network::link_nodes(&mut self.nodes, w[0], w[1]);
        }
        let new_ids = &out.new_ids;
        if whole {
            if let (Some(&a), Some(&b)) = (new_ids.first(), new_ids.last()) {
                self.link(b, a);
            }
        } else {
            match (new_ids.first(), new_ids.last()) {
                (Some(&a), Some(&b)) => {
                    if before != NONE {
                        self.link(before, a);
                    } else {
                        self.nodes[a as usize].as_mut().expect("new").prev = NONE;
                    }
                    if after != NONE {
                        self.link(b, after);
                    }
                }
                _ => {
                    if before != NONE && after != NONE {
                        self.link(before, after);
                    } else if before != NONE {
                        self.nodes[before as usize].as_mut().expect("live").next = NONE;
                    } else if after != NONE {
                        self.nodes[after as usize].as_mut().expect("live").prev = NONE;
                    }
                }
            }
        }

        // The path's bookkeeping.
        let n = data.n_segments as usize - run.len() + new.len();
        let data = &mut self.paths[p];
        data.n_segments = n as u32;
        if n == 0 {
            data.first = NONE;
        } else if run.contains(&data.first) {
            data.first = if !data.closed && before != NONE {
                // An open path whose first segment was removed has no `before`.
                return internal(format!("path {p}: its first segment has a predecessor"));
            } else {
                new_ids.first().copied().unwrap_or(after)
            };
        }
        out.after = (after != NONE).then_some(after);
        out.contracted = n == 0;
        Ok(())
    }

    /// After the triangulation split edge `e` into the halves `[e, e2]` (in order from the old
    /// tail, both directed like the old edge: [`Insertion::halves`]): every segment on `e`
    /// becomes two segments, one on each half, and keeps its place across the edge (SPEC §6.1,
    /// F19). The segment keeps its id for the half it runs along first (`e` if it is forward,
    /// `e2` if backward); the other half gets a new id, linked after it in its path. Both
    /// halves list their segments in the old order of `e`: an unchanged direction means an
    /// unchanged side of `edge_halfedge`, so the layering of the halves is that of `e`.
    ///
    /// The layers and the pinned mask grow to the triangulation's new size first. Returns the
    /// pairs `(first, second)` of ids in path order, in the old layer order of `e`.
    ///
    /// # Errors
    /// [`FlipError::Internal`] if `e2` already carries segments or a segment on `e` is not
    /// live.
    pub(crate) fn split_segments(
        &mut self,
        complex: &Complex,
        e: usize,
        e2: usize,
    ) -> Result<Vec<(SegmentId, SegmentId)>, FlipError> {
        self.layers.resize(complex.n_edges(), VecDeque::new());
        if self.pinned.len() < complex.n_vertices() {
            self.pinned.resize(complex.n_vertices(), false);
        }
        if !self.layers[e2].is_empty() {
            return Err(FlipError::Internal {
                what: format!("split of edge {e}: the new half {e2} already carries segments"),
            });
        }
        let old = std::mem::take(&mut self.layers[e]);
        let mut on_e = VecDeque::with_capacity(old.len());
        let mut on_e2 = VecDeque::with_capacity(old.len());
        let mut pairs = Vec::with_capacity(old.len());
        for &id in &old {
            let Some(node) = self.nodes.get(id as usize).copied().flatten() else {
                return Err(FlipError::Internal {
                    what: format!("split of edge {e}: its layer lists segment {id}, not live"),
                });
            };
            let forward = node.segment.forward;
            let new = self.nodes.len() as u32;
            // In path order: forward a → m (e) then m → b (e2); backward b → m (e2) then
            // m → a (e).
            let (first_edge, second_edge) = if forward { (e, e2) } else { (e2, e) };
            self.nodes.push(Some(Node {
                path: node.path,
                segment: Segment::new(second_edge, forward),
                prev: id,
                next: node.next,
            }));
            let first = self.nodes[id as usize].as_mut().expect("live");
            first.segment = Segment::new(first_edge, forward);
            first.next = new;
            // `node.next` is read now, so a run of segments of one path on `e` (out and back)
            // links through the halves made before.
            let after = node.next;
            if after != NONE {
                self.nodes[after as usize].as_mut().expect("live").prev = new;
            }
            if forward {
                on_e.push_back(id);
                on_e2.push_back(new);
            } else {
                on_e.push_back(new);
                on_e2.push_back(id);
            }
            self.paths[node.path as usize].n_segments += 1;
            pairs.push((SegmentId(id), SegmentId(new)));
        }
        self.layers[e] = on_e;
        self.layers[e2] = on_e2;
        Ok(pairs)
    }

    /// Links `a` → `b`.
    fn link(&mut self, a: u32, b: u32) {
        Network::link_nodes(&mut self.nodes, a, b);
    }

    fn link_nodes(nodes: &mut [Option<Node>], a: u32, b: u32) {
        nodes[a as usize].as_mut().expect("live").next = b;
        nodes[b as usize].as_mut().expect("live").prev = a;
    }
}

/// Checks that consecutive segments of `path` share a vertex, and that a loop is closed.
fn check_connected(complex: &Complex, p: usize, path: &EdgePath) -> Result<(), FlipError> {
    let s = &path.segments;
    for i in 0..s.len().saturating_sub(1) {
        let (end, start) = (s[i].head(complex), s[i + 1].tail(complex));
        if end != start {
            return Err(FlipError::PathNotConnected {
                path: p,
                index: i,
                end,
                start,
            });
        }
    }
    if path.closed {
        let (first, last) = (s[0].tail(complex), s[s.len() - 1].head(complex));
        if first != last {
            return Err(FlipError::LoopNotClosed {
                path: p,
                first,
                last,
            });
        }
    }
    Ok(())
}

/// A triangulation `T` together with an edge network on it (SPEC §2.3): what FlipOut
/// shortens.
#[derive(Debug, Clone)]
pub struct GeodesicNetwork<T: Triangulation> {
    pub(crate) tri: T,
    pub(crate) network: Network,
    /// The running sum of the segment lengths (F14).
    pub(crate) length: f64,
}

impl<T: Triangulation> GeodesicNetwork<T> {
    /// A network of `paths` on `tri` (see [`Network::new`] for `pinned`, `layering` and the
    /// checks). The triangulation is taken over; to keep it when the input is refused, check
    /// the input first with [`Network::new`] on `tri.complex()` and use
    /// [`from_parts`](Self::from_parts).
    ///
    /// # Errors
    /// As [`Network::new`].
    pub fn new(
        tri: T,
        paths: &[EdgePath],
        pinned: &[usize],
        layering: Option<&Layering>,
    ) -> Result<GeodesicNetwork<T>, FlipError> {
        let network = Network::new(tri.complex(), paths, pinned, layering)?;
        Ok(GeodesicNetwork::assemble(tri, network))
    }

    /// A triangulation and a network made for it.
    ///
    /// # Errors
    /// As [`Network::check`] on `tri.complex()`.
    pub fn from_parts(tri: T, network: Network) -> Result<GeodesicNetwork<T>, FlipError> {
        network.check(tri.complex())?;
        Ok(GeodesicNetwork::assemble(tri, network))
    }

    fn assemble(tri: T, network: Network) -> GeodesicNetwork<T> {
        let mut out = GeodesicNetwork {
            tri,
            network,
            length: 0.0,
        };
        out.length = out.recompute_length();
        out
    }

    /// The triangulation.
    pub fn triangulation(&self) -> &T {
        &self.tri
    }

    /// The triangulation, every path an edge path of it (SPEC Fig. 1 bottom).
    pub fn into_triangulation(self) -> T {
        self.tri
    }

    /// The triangulation and the network.
    pub fn into_parts(self) -> (T, Network) {
        (self.tri, self.network)
    }

    /// The network.
    pub fn network(&self) -> &Network {
        &self.network
    }

    /// The paths, in order.
    pub fn paths(&self) -> Vec<EdgePath> {
        self.network.paths()
    }

    /// The total length of all paths, kept as a running sum (F14). Initially the sum over
    /// paths in order of the segment lengths in path order.
    pub fn length(&self) -> f64 {
        self.length
    }

    /// The total length summed afresh, in the order of [`length`](Self::length)'s initial
    /// value.
    pub fn recompute_length(&self) -> f64 {
        (0..self.network.n_paths())
            .map(|p| self.path_length(p))
            .fold(0.0, |s, l| s + l)
    }

    /// The length of path `p`.
    ///
    /// # Panics
    /// If `p` is out of range.
    pub fn path_length(&self, p: usize) -> f64 {
        self.network
            .path_segment_ids(p)
            .into_iter()
            .map(|id| {
                self.tri
                    .edge_length(self.network.node(id.0).segment.edge as usize)
            })
            .fold(0.0, |s, l| s + l)
    }

    /// Checks the network against the current triangulation ([`Network::check`]).
    ///
    /// # Errors
    /// As [`Network::check`].
    pub fn check_invariants(&self) -> Result<(), FlipError> {
        self.network.check(self.tri.complex())
    }

    /// The wedges of a joint ([`JointGeometry`]).
    ///
    /// # Errors
    /// [`FlipError::IndexOutOfRange`] if the joint's segments are not consecutive live
    /// segments of the network.
    pub fn joint_geometry(&self, joint: Joint) -> Result<JointGeometry, FlipError> {
        JointGeometry::new(&self.tri, &self.network, joint)
    }

    /// The geodesic certificate of the current state (SPEC §7 layer 1a).
    ///
    /// # Errors
    /// As [`crate::certify::certify`].
    pub fn certify(&self, options: &CertifyOptions) -> Result<GeodesicCertificate, FlipError> {
        crate::certify::certify(&self.tri, &self.network, options)
    }

    /// Splits edge `e` of the triangulation at the parameter `t` from the tail of its
    /// halfedge ([`Triangulation::split_edge`], paper Alg. 8) and carries the network along
    /// (SPEC §6.1, F19): every segment on `e`, however many share it, becomes two segments,
    /// one on each half ([`Insertion::halves`]), in the same place across the edge as before.
    /// The new vertex is not pinned; the joint there is straight (both sides `π` up to
    /// rounding). The length is summed afresh (the halves' lengths `t ℓ` and `(1 − t) ℓ` add up
    /// to `ℓ` only up to rounding).
    ///
    /// # Errors
    /// [`FlipError::IndexOutOfRange`] for an edge out of range; the refusals of the
    /// triangulation's split as [`FlipError::Intrinsic`] (a parameter not strictly between 0
    /// and 1, an edge with the same face on both sides, a split too thin to be described by
    /// its lengths, ...). Nothing is changed then. [`FlipError::Internal`] if the split did
    /// not give two halves directed like the old edge (a bug of the triangulation).
    pub fn split_edge(&mut self, e: usize, t: f64) -> Result<NetworkSplit, FlipError> {
        let n_edges = self.tri.complex().n_edges();
        if e >= n_edges {
            return Err(FlipError::IndexOutOfRange {
                what: "edge",
                index: e,
                bound: n_edges,
            });
        }
        let [a, b] = self.tri.complex().edge_vertices(e);
        let insertion = self.tri.split_edge(e, t)?;
        let m = insertion.vertex as usize;
        let complex = self.tri.complex();
        let halves = match insertion.halves {
            Some([h0, h1])
                if h0 as usize == e
                    && complex.edge_vertices(h0 as usize) == [a, m]
                    && complex.edge_vertices(h1 as usize) == [m, b] =>
            {
                [h0, h1]
            }
            other => {
                return Err(FlipError::Internal {
                    what: format!(
                        "the split of edge {e} ({a} → {b}) at vertex {m} gave the halves \
                         {other:?}, not [{e}, _] directed like the old edge"
                    ),
                });
            }
        };
        let segments = self
            .network
            .split_segments(complex, e, halves[1] as usize)?;
        self.length = self.recompute_length();
        Ok(NetworkSplit {
            insertion,
            segments,
        })
    }
}

/// What [`GeodesicNetwork::split_edge`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct NetworkSplit {
    /// What the triangulation did: the new vertex, its place on the input mesh, the new edges
    /// and faces, the two halves.
    pub insertion: Insertion,
    /// Per segment that lay on the split edge, in its old layer order: the ids of its two
    /// parts in path order (the first keeps the old id).
    pub segments: Vec<(SegmentId, SegmentId)>,
}
