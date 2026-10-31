//! The geometry of a joint (SPEC §3, paper Sec. 3.2.3–3.2.5): the two wedges, their angles,
//! the classification straight / left / right, and whether a wedge is flexible or blocked by
//! a path end.
//!
//! # The star of a vertex
//!
//! [`VertexStar`] lists, around a vertex `b`, its **edge ends** counter-clockwise (one per
//! outgoing halfedge, plus at a boundary vertex the incoming boundary halfedge last; a
//! self-edge at `b` has two ends there), and in each edge end the **segment ends** of the
//! segments on that edge, in the order met when turning counter-clockwise (the layering
//! read from the face on the clockwise side to the face on the counter-clockwise side).
//! Between edge end `i` and edge end `i + 1` lies the corner of outgoing halfedge `i`; at a
//! boundary vertex the boundary gap lies between the last edge end and the first.
//!
//! A joint `γ_abc` pairs the end of `ab` at `b` (position `P`) with the start of `bc` (`Q`).
//! Its **right wedge** is everything swept counter-clockwise from `P` to `Q`, its **left
//! wedge** everything swept counter-clockwise from `Q` to `P` (the paper's `W_abc`, "the wedge
//! to the left of `γ_abc`", Sec. 3.2.3; the reference's `leftAngle`, `gc: fg.cpp:489-522`).
//! A wedge contains the corners it sweeps over. This handles the Δ-complex cases without
//! special code: a self-edge at `b` inside a wedge contributes two edge ends (it is incident
//! on the wedge twice, App. A), and when `ab` and `bc` lie on the same edge end (out and back
//! along one edge, e.g. at a vertex of degree one) the side between the two segments in the
//! layer order sweeps no corner (`α = 0`, `deg W = 0`, SPEC §4.1 step 3) and the other side
//! sweeps all of them.

use std::f64::consts::PI;

use ddg_intrinsic::{Complex, Triangulation};

use crate::error::FlipError;
use crate::network::{Joint, NONE, Network, SegmentId};

/// A side of a joint. Ordered `Left < Right` (F4).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Side {
    /// The wedge to the left of the path (`W_abc`).
    Left,
    /// The wedge to the right of the path.
    Right,
}

impl Side {
    /// The other side.
    pub fn opposite(self) -> Side {
        match self {
            Side::Left => Side::Right,
            Side::Right => Side::Left,
        }
    }
}

/// The classification of a joint (SPEC §3; `gc: fg.cpp:524-587`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Turn {
    /// The smaller wedge angle exceeds `π − tolerance`: locally shortest (Eq. 2).
    Straight,
    /// The left wedge is the smaller one and below `π − tolerance`.
    Left,
    /// The right wedge is the smaller one (or both are equal) and below `π − tolerance`.
    Right,
}

/// An end of an edge at a vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeEnd {
    /// The edge.
    pub edge: u32,
    /// Whether the vertex is the tail of the edge's `edge_halfedge` (else its head).
    pub at_tail: bool,
    /// The halfedge of this end: outgoing from the vertex, or for the last end of a boundary
    /// vertex the incoming boundary halfedge.
    pub halfedge: u32,
    /// The corner counter-clockwise after this end (the outgoing halfedge, which names its
    /// corner); `None` for the last end of a boundary vertex, after which lies the boundary
    /// gap.
    pub corner: Option<u32>,
}

/// An end of a segment at a vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SegmentEnd {
    /// The segment.
    pub segment: SegmentId,
    /// Whether the segment starts here (else it ends here).
    pub start: bool,
    /// The index of its edge end in [`VertexStar::ends`].
    pub edge_end: u32,
}

/// The edge ends and segment ends around a vertex, counter-clockwise (see the
/// [module documentation](self)).
#[derive(Debug, Clone, PartialEq)]
pub struct VertexStar {
    /// The vertex.
    pub vertex: usize,
    /// Whether it lies on the boundary.
    pub boundary: bool,
    /// The edge ends, counter-clockwise from the vertex's `outgoing_halfedge`.
    pub ends: Vec<EdgeEnd>,
    /// The segment ends, counter-clockwise: edge end by edge end, each in layer order as
    /// met turning counter-clockwise.
    pub entries: Vec<SegmentEnd>,
    /// Per edge end, the index of its first entry; one more element at the end.
    pub first_entry: Vec<u32>,
}

/// A sweep counter-clockwise around a vertex from one segment end to another.
#[derive(Debug, Clone, PartialEq)]
pub struct Wedge {
    /// Its angle: the sum of its corner angles, summed in `corners` order (Eq. 1); `+∞` if it
    /// crosses the boundary gap (Sec. 3.2.3).
    pub angle: f64,
    /// The corners swept, counter-clockwise from the start of the sweep, as the outgoing
    /// halfedges that name them. A corner can appear only once; for a boundary wedge the
    /// corners on both sides of the gap are listed.
    pub corners: Vec<u32>,
    /// Whether the sweep crosses the boundary gap.
    pub boundary: bool,
    /// Where the sweep starts (index into [`VertexStar::entries`]).
    pub from: u32,
    /// Where it ends.
    pub to: u32,
}

impl Wedge {
    /// The degree `deg W`: the number of corners (triangles) in the wedge.
    pub fn degree(&self) -> usize {
        self.corners.len()
    }
}

impl VertexStar {
    /// The star of vertex `v` with the segment ends of `network`.
    ///
    /// # Panics
    /// If `v` is out of range or the network's layers do not fit `complex`.
    pub fn new(complex: &Complex, network: &Network, v: usize) -> VertexStar {
        let mut ends: Vec<EdgeEnd> = complex
            .outgoing_halfedges(v)
            .map(|h| EdgeEnd {
                edge: complex.edge(h) as u32,
                at_tail: complex.is_forward(h),
                halfedge: h as u32,
                corner: Some(h as u32),
            })
            .collect();
        let boundary = complex.is_boundary_vertex(v);
        if boundary {
            let last = ends
                .last()
                .expect("a vertex has an outgoing halfedge")
                .halfedge as usize;
            let g = complex.prev(last);
            ends.push(EdgeEnd {
                edge: complex.edge(g) as u32,
                at_tail: false,
                halfedge: g as u32,
                corner: None,
            });
        }
        let mut entries = Vec::new();
        let mut first_entry = Vec::with_capacity(ends.len() + 1);
        for (i, end) in ends.iter().enumerate() {
            first_entry.push(entries.len() as u32);
            let layer = &network.layers[end.edge as usize];
            // Index 0 of a layer is next to the face of edge_halfedge. At the tail of
            // edge_halfedge that face is counter-clockwise after the edge, at its head
            // clockwise before it.
            let mut push = |id: u32| {
                let forward = network.node(id).segment.forward;
                entries.push(SegmentEnd {
                    segment: SegmentId(id),
                    start: if end.at_tail { forward } else { !forward },
                    edge_end: i as u32,
                });
            };
            if end.at_tail {
                layer.iter().rev().for_each(|&id| push(id));
            } else {
                layer.iter().for_each(|&id| push(id));
            }
        }
        first_entry.push(entries.len() as u32);
        VertexStar {
            vertex: v,
            boundary,
            ends,
            entries,
            first_entry,
        }
    }

    /// The position in [`entries`](Self::entries) of the start (`start = true`) or the end of
    /// segment `id`, if it is at this vertex.
    pub fn position(&self, id: SegmentId, start: bool) -> Option<usize> {
        self.entries
            .iter()
            .position(|e| e.segment == id && e.start == start)
    }

    /// The sweep counter-clockwise from entry `from` to entry `to` (distinct).
    ///
    /// # Panics
    /// If an index is out of range.
    pub fn sweep<T: Triangulation>(&self, tri: &T, from: usize, to: usize) -> Wedge {
        let (i, j) = (
            self.entries[from].edge_end as usize,
            self.entries[to].edge_end as usize,
        );
        let n = self.ends.len();
        let mut corners = Vec::new();
        let mut boundary = false;
        if !(i == j && to > from) {
            let mut k = i;
            loop {
                match self.ends[k].corner {
                    Some(h) => corners.push(h),
                    None => boundary = true,
                }
                k = (k + 1) % n;
                if k == j {
                    break;
                }
            }
        }
        let angle = if boundary {
            f64::INFINITY
        } else {
            wedge_angle(tri, &corners)
        };
        Wedge {
            angle,
            corners,
            boundary,
            from: from as u32,
            to: to as u32,
        }
    }

    /// The entries strictly inside a wedge, counter-clockwise.
    pub fn inside(&self, wedge: &Wedge) -> impl Iterator<Item = &SegmentEnd> + '_ {
        let m = self.entries.len();
        let (from, to) = (wedge.from as usize, wedge.to as usize);
        let count = (to + m - from) % m;
        (1..count).map(move |k| &self.entries[(from + k) % m])
    }

    /// Whether no segment end lies strictly inside the wedge: the joint's own segments are
    /// outermost on the wedge side of their edges and no edge inside the wedge carries a
    /// segment. With `from`, `to` the joint's own ends, this is the paper's *flexible*
    /// (Sec. 3.2.3; `gc: fg.cpp:589-658`).
    pub fn is_clear(&self, wedge: &Wedge) -> bool {
        self.inside(wedge).next().is_none()
    }

    /// Whether the wedge is blocked by a path end (`gc: fg.cpp:662-730`): a segment on an
    /// edge end strictly inside the wedge, or the segment outermost on the wedge side of the
    /// two bounding edge ends (if not the wedge's own), is the first or last segment of an
    /// open path. For a wedge of degree zero (both ends on one edge end) the segments between
    /// the two are checked.
    pub fn is_blocked_by_end(&self, network: &Network, wedge: &Wedge) -> bool {
        let is_end = |e: &SegmentEnd| network.is_path_end(e.segment);
        let (from, to) = (wedge.from as usize, wedge.to as usize);
        let (i, j) = (
            self.entries[from].edge_end as usize,
            self.entries[to].edge_end as usize,
        );
        if i == j && to > from {
            return self.inside(wedge).any(is_end);
        }
        let n = self.ends.len();
        // Edge ends strictly inside.
        let mut k = (i + 1) % n;
        while k != j {
            let range = self.first_entry[k] as usize..self.first_entry[k + 1] as usize;
            if self.entries[range].iter().any(is_end) {
                return true;
            }
            k = (k + 1) % n;
        }
        // The bounding edge ends: the wedge lies counter-clockwise after `from` and before
        // `to`, so their outermost segments there are the last entry of end i and the first
        // of end j.
        let last_i = self.first_entry[i + 1] as usize - 1;
        let first_j = self.first_entry[j] as usize;
        (last_i != from && is_end(&self.entries[last_i]))
            || (first_j != to && is_end(&self.entries[first_j]))
    }

    /// Checks that the joints at this vertex do not cross (SPEC §2.1): every joint pairs the
    /// end of its incoming segment with the start of its outgoing one, and no two pairs
    /// interleave in the cyclic order of [`entries`](Self::entries). Path ends pair nothing.
    ///
    /// # Errors
    /// [`FlipError::Crossing`] naming the vertex and the paths of two interleaving joints.
    pub fn check_non_crossing(&self, network: &Network) -> Result<(), FlipError> {
        let m = self.entries.len();
        // Partner keys, then positions by binary search.
        let mut index: Vec<(u32, bool, u32)> = self
            .entries
            .iter()
            .enumerate()
            .map(|(k, e)| (e.segment.0, e.start, k as u32))
            .collect();
        index.sort_unstable();
        let find = |id: u32, start: bool| -> Option<usize> {
            index
                .binary_search_by(|&(s, st, _)| (s, st).cmp(&(id, start)))
                .ok()
                .map(|k| index[k].2 as usize)
        };
        let mut partner = vec![usize::MAX; m];
        for (k, e) in self.entries.iter().enumerate() {
            let node = network.node(e.segment.0);
            let other = if e.start {
                (node.prev != NONE).then(|| find(node.prev, false))
            } else {
                (node.next != NONE).then(|| find(node.next, true))
            };
            match other {
                None => {}
                Some(Some(q)) => partner[k] = q,
                Some(None) => {
                    return Err(FlipError::Internal {
                        what: format!(
                            "segment {} has a joint at vertex {}, but its neighbour's end is \
                             not there",
                            e.segment.0, self.vertex
                        ),
                    });
                }
            }
        }
        let mut stack: Vec<usize> = Vec::new();
        for (k, &q) in partner.iter().enumerate() {
            if q == usize::MAX {
                continue;
            }
            if q > k {
                stack.push(k);
            } else if stack.last() == Some(&q) {
                stack.pop();
            } else {
                let path = |k: usize| network.node(self.entries[k].segment.0).path as usize;
                let top = *stack.last().expect("an open chord precedes a closing one");
                return Err(FlipError::Crossing {
                    vertex: self.vertex,
                    paths: [path(top), path(k)],
                });
            }
        }
        Ok(())
    }

    /// The sum of the corner angles at this vertex (`Θ_b` up to rounding), in edge-end order.
    pub fn angle_sum<T: Triangulation>(&self, tri: &T) -> f64 {
        let corners: Vec<u32> = self.ends.iter().filter_map(|e| e.corner).collect();
        wedge_angle(tri, &corners)
    }
}

/// The angle of a wedge, Eq. 1: the sum of its corner angles `θ` (from the edge lengths,
/// `ddg_core::geometry::corner_angle` through [`Triangulation::corner_angle`]), in the order
/// given (F2).
pub fn wedge_angle<T: Triangulation>(tri: &T, corners: &[u32]) -> f64 {
    corners.iter().map(|&h| tri.corner_angle(h as usize)).sum()
}

/// The two wedges of a joint, with the star of its vertex.
#[derive(Debug, Clone, PartialEq)]
pub struct JointGeometry {
    /// The joint.
    pub joint: Joint,
    /// The star of its vertex `b`.
    pub star: VertexStar,
    /// The left wedge `W_abc`: counter-clockwise from the start of `bc` to the end of `ab`.
    pub left: Wedge,
    /// The right wedge: counter-clockwise from the end of `ab` to the start of `bc`.
    pub right: Wedge,
}

impl JointGeometry {
    /// The wedges of `joint` in `network` on `tri`.
    ///
    /// # Errors
    /// [`FlipError::IndexOutOfRange`] if the joint is not two consecutive live segments.
    pub fn new<T: Triangulation>(
        tri: &T,
        network: &Network,
        joint: Joint,
    ) -> Result<JointGeometry, FlipError> {
        let bad = || FlipError::IndexOutOfRange {
            what: "joint (live consecutive segments)",
            index: joint.outgoing.0 as usize,
            bound: network.segment_id_bound(),
        };
        let outgoing = network.segment(joint.outgoing).ok_or_else(bad)?;
        if outgoing.prev != Some(joint.incoming) {
            return Err(bad());
        }
        let complex = tri.complex();
        let v = outgoing.segment.tail(complex);
        let star = VertexStar::new(complex, network, v);
        let p = star.position(joint.incoming, false).ok_or_else(bad)?;
        let q = star.position(joint.outgoing, true).ok_or_else(bad)?;
        let right = star.sweep(tri, p, q);
        let left = star.sweep(tri, q, p);
        Ok(JointGeometry {
            joint,
            star,
            left,
            right,
        })
    }

    /// The vertex `b`.
    pub fn vertex(&self) -> usize {
        self.star.vertex
    }

    /// The wedge on `side`.
    pub fn wedge(&self, side: Side) -> &Wedge {
        match side {
            Side::Left => &self.left,
            Side::Right => &self.right,
        }
    }

    /// The angle on `side`.
    pub fn angle(&self, side: Side) -> f64 {
        self.wedge(side).angle
    }

    /// The side of the smaller angle: left if strictly smaller, else right (as the reference,
    /// `gc: fg.cpp:530-543`).
    pub fn min_side(&self) -> Side {
        if self.left.angle < self.right.angle {
            Side::Left
        } else {
            Side::Right
        }
    }

    /// `α^min_b`.
    pub fn min_angle(&self) -> f64 {
        self.angle(self.min_side())
    }

    /// `α^max_b`.
    pub fn max_angle(&self) -> f64 {
        self.angle(self.min_side().opposite())
    }

    /// The classification (SPEC §3): [`Turn::Straight`] if `α^min > π − tolerance`
    /// (`gc: fg.cpp:534-543`, default tolerance `1e-5`, F3), otherwise the side of the
    /// smaller wedge.
    pub fn turn(&self, tolerance: f64) -> Turn {
        if self.min_angle() > PI - tolerance {
            Turn::Straight
        } else {
            match self.min_side() {
                Side::Left => Turn::Left,
                Side::Right => Turn::Right,
            }
        }
    }

    /// The sides to queue (SPEC §3, `gc: fg.cpp:998-1007`): none if straight; the smaller side,
    /// and the larger one too if it is also below `π − tolerance` (`Θ_b < 2π`).
    pub fn sides_to_shorten(&self, tolerance: f64) -> Vec<(Side, f64)> {
        let min = self.min_side();
        [min, min.opposite()]
            .into_iter()
            .map(|s| (s, self.angle(s)))
            .take_while(|&(_, a)| a <= PI - tolerance)
            .collect()
    }

    /// Whether the wedge on `side` is flexible: the joint's segments are outermost on that
    /// side of their edges and no edge inside it carries a segment (Sec. 3.2.3).
    pub fn is_flexible(&self, side: Side) -> bool {
        self.star.is_clear(self.wedge(side))
    }

    /// Whether the wedge on `side` is blocked by a path end ([`VertexStar::is_blocked_by_end`]):
    /// the joint is then locally shortest in its isotopy class (Sec. 4.2, case (b)).
    pub fn is_blocked_by_end(&self, network: &Network, side: Side) -> bool {
        self.star.is_blocked_by_end(network, self.wedge(side))
    }

    /// The corners of the wedge on `side` ordered from `ba` to `bc` (the order of the outer
    /// arc `a = n₀, …, n_k = c` of Alg. 1): clockwise for the left wedge, counter-clockwise
    /// for the right one.
    pub fn corners_from_incoming(&self, side: Side) -> Vec<u32> {
        match side {
            Side::Right => self.right.corners.clone(),
            Side::Left => self.left.corners.iter().rev().copied().collect(),
        }
    }
}
