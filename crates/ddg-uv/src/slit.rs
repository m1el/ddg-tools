//! The automatic slit of a closed piece: a simple path of polygon edges.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use ddg_core::operators::angle_defects;
use ddg_core::topology::NONE;

use crate::error::{UvError, internal};
use crate::options::SlitPlacement;
use crate::prepare::Piece;

#[derive(PartialEq)]
struct Entry {
    distance: f64,
    vertex: u32,
}

impl Eq for Entry {}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        // A min-heap on (distance, vertex).
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.vertex.cmp(&self.vertex))
    }
}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Shortest paths from `source` along the polygon edges of the piece, avoiding blocked
/// vertices: `(distance, predecessor)` per vertex, `∞` and `NONE` where unreachable.
fn dijkstra(piece: &Piece, blocked: &[bool], source: usize) -> (Vec<f64>, Vec<u32>) {
    let mesh = &piece.mesh;
    let lengths = piece.geometry.edge_lengths();
    let mut distance = vec![f64::INFINITY; mesh.n_vertices()];
    let mut parent = vec![NONE; mesh.n_vertices()];
    let mut heap = BinaryHeap::new();
    distance[source] = 0.0;
    heap.push(Entry {
        distance: 0.0,
        vertex: source as u32,
    });
    while let Some(Entry { distance: d, vertex }) = heap.pop() {
        let v = vertex as usize;
        if d > distance[v] {
            continue;
        }
        // A closed piece: every edge at `v` is the edge of an outgoing halfedge.
        for h in mesh.outgoing_halfedges(v) {
            let e = mesh.edge(h);
            let w = mesh.to_vertex(h);
            if !piece.is_polygon_edge[e] || blocked[w] {
                continue;
            }
            let candidate = d + lengths[e];
            if candidate < distance[w] || (candidate == distance[w] && vertex < parent[w]) {
                let improved = candidate < distance[w];
                distance[w] = candidate;
                parent[w] = vertex;
                if improved {
                    heap.push(Entry {
                        distance: candidate,
                        vertex: w as u32,
                    });
                }
            }
        }
    }
    (distance, parent)
}

/// The reachable vertex farthest from the source (ties: lowest index).
fn farthest(distance: &[f64]) -> usize {
    let mut best = 0usize;
    let mut best_distance = f64::NEG_INFINITY;
    for (v, &d) in distance.iter().enumerate() {
        if d.is_finite() && d > best_distance {
            best = v;
            best_distance = d;
        }
    }
    best
}

/// The path from the source of `parent` to `target`, as vertices in order.
fn path_to(parent: &[u32], target: usize) -> Vec<usize> {
    let mut path = vec![target];
    let mut v = target;
    while parent[v] != NONE {
        v = parent[v] as usize;
        path.push(v);
    }
    path.reverse();
    path
}

/// A slit: its edges (ids of the piece's mesh, ascending) and its length.
pub(crate) struct Slit {
    pub edges: Vec<u32>,
    pub length: f64,
}

/// Chooses the slit of a closed piece. The slit is a simple path of at least two polygon
/// edges that avoids the vertices in `avoid` if it can.
pub(crate) fn choose_slit(
    piece: &Piece,
    avoid: &[bool],
    placement: SlitPlacement,
    relative_length: f64,
) -> Result<Slit, UvError> {
    match slit_avoiding(piece, avoid, placement, relative_length) {
        Some(slit) => Ok(slit),
        None => {
            let nothing = vec![false; piece.mesh.n_vertices()];
            slit_avoiding(piece, &nothing, placement, relative_length).ok_or_else(|| {
                internal(format!(
                    "no slit of two edges was found on the closed piece that contains face {}",
                    piece.first_polygon
                ))
            })
        }
    }
}

fn slit_avoiding(
    piece: &Piece,
    blocked: &[bool],
    placement: SlitPlacement,
    relative_length: f64,
) -> Option<Slit> {
    let mesh = &piece.mesh;
    let lengths = piece.geometry.edge_lengths();
    let defects = angle_defects(mesh, &piece.geometry).ok()?;
    // The vertex of largest absolute angle defect (ties: lowest index).
    let mut start: Option<usize> = None;
    for v in (0..mesh.n_vertices()).filter(|&v| !blocked[v]) {
        if start.is_none_or(|s| defects[v].abs() > defects[s].abs()) {
            start = Some(v);
        }
    }
    let start = start?;
    let (distance, parent) = dijkstra(piece, blocked, start);
    let a = farthest(&distance);
    let mut path = match placement {
        SlitPlacement::AtLargestCurvature => path_to(&parent, a),
        SlitPlacement::DiameterMiddle => {
            let (distance, parent) = dijkstra(piece, blocked, a);
            path_to(&parent, farthest(&distance))
        }
    };
    // A path of fewer than two edges (tiny meshes): extend it by the shortest polygon edge at
    // its end that leads to a new vertex.
    while path.len() < 3 {
        let end = *path.last()?;
        let mut best: Option<(f64, usize)> = None;
        for h in mesh.outgoing_halfedges(end) {
            let (e, w) = (mesh.edge(h), mesh.to_vertex(h));
            if piece.is_polygon_edge[e] && !blocked[w] && !path.contains(&w) {
                let key = (lengths[e], w);
                if best.is_none_or(|b| key < b) {
                    best = Some(key);
                }
            }
        }
        path.push(best?.1);
    }
    let edge_of = |i: usize| mesh.find_edge(path[i], path[i + 1]);
    let n = path.len() - 1;
    let mut edge_lengths = Vec::with_capacity(n);
    for i in 0..n {
        edge_lengths.push(lengths[edge_of(i)?]);
    }
    let total: f64 = edge_lengths.iter().sum();
    let target = relative_length * total;

    // The range of path edges `first..last` that is kept.
    let (first, last) = match placement {
        SlitPlacement::AtLargestCurvature => {
            let mut last = 0;
            let mut covered = 0.0;
            while last < n && (last < 2 || covered < target) {
                covered += edge_lengths[last];
                last += 1;
            }
            (0, last)
        }
        SlitPlacement::DiameterMiddle => {
            // The edges whose midpoints lie within half the target length of the middle of
            // the path; at least the two edges around the vertex nearest to the middle.
            let mut position = 0.0;
            let mut middle_vertex = 1usize;
            let mut best_gap = f64::INFINITY;
            let (mut first, mut last) = (n, 0usize);
            for (i, &length) in edge_lengths.iter().enumerate() {
                let midpoint = position + 0.5 * length;
                if (midpoint - 0.5 * total).abs() <= 0.5 * target {
                    first = first.min(i);
                    last = last.max(i + 1);
                }
                position += length;
                // Vertex i + 1 lies at `position`; interior vertices of the path only.
                if i + 1 < n && (position - 0.5 * total).abs() < best_gap {
                    best_gap = (position - 0.5 * total).abs();
                    middle_vertex = i + 1;
                }
            }
            if first >= last {
                (first, last) = (middle_vertex - 1, middle_vertex + 1);
            }
            (first.min(middle_vertex - 1), last.max(middle_vertex + 1))
        }
    };
    let mut edges = Vec::with_capacity(last - first);
    let mut length = 0.0;
    for i in first..last {
        edges.push(edge_of(i)? as u32);
        length += edge_lengths[i];
    }
    edges.sort_unstable();
    Some(Slit { edges, length })
}
