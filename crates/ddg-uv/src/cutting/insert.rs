//! Inserting a cut into a polygon mesh: the crossed polygons are split along the paths.
//!
//! This is the operation an editor performs with its own mesh structure (in Blender: with
//! bmesh), written once in Rust so that the tests and the benchmarks can run it, and so that the
//! result can be compared with what the editor does.
//!
//! # What is built
//!
//! * Every [`EdgePoint`] becomes a vertex on its polygon edge, so both polygons of the edge get
//!   it as a corner. On a seam edge the two sides have their own points; each point becomes a
//!   vertex on the edge and **both** polygons get both vertices, in order of the parameter
//!   (two points with exactly the same parameter on the two sides are one vertex).
//!   A polygon that only receives such vertices is not split but has more corners; they are
//!   collinear with its edge.
//! * Every [`InteriorPoint`] becomes a vertex inside the cut.
//! * A polygon crossed by paths is split along them into one new polygon per part. A path cuts
//!   each triangle of the polygon's triangulation in a straight segment; the new polygons are
//!   unions of the resulting pieces of triangles, so their triangulation is made of those
//!   pieces (a piece with more than three corners is split by its better diagonal).
//! * Old vertices keep their numbers; new vertices follow (edge points in the order of
//!   [`CutOnPolygons::edge_points`], then interior points in the order of the paths). Polygons
//!   are numbered in the order of their parents; the parts of one polygon follow each other,
//!   starting with the part that has the first corner of the parent.
//!
//! Nothing is normalized: a description that does not fit the mesh is an error that names the
//! polygon (or the edge point).

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use ddg_core::topology::NONE;

use super::describe::{CutOnPolygons, PathEnd};
use super::error::{CuttingError, internal, invalid};
use crate::input::PolygonMesh;

/// Where a vertex of the new mesh comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VertexSource {
    /// A vertex of the old mesh (same index).
    Old(u32),
    /// An edge point (an index into [`CutOnPolygons::edge_points`]); for a point shared by the
    /// two sides of a seam, the first of the two.
    EdgePoint(u32),
    /// A junction (an index into [`CutOnPolygons::junctions`]): a vertex inside a triangle of
    /// its polygon where three paths meet.
    Junction(u32),
    /// An interior point: the `index`-th interior point of path `path`.
    Interior {
        /// The path (an index into [`CutOnPolygons::paths`]).
        path: u32,
        /// The position in [`PolygonPath::interior`](super::PolygonPath::interior).
        index: u32,
    },
}

/// The result of [`insert`].
#[derive(Debug, Clone, PartialEq)]
pub struct Inserted {
    /// The new polygon mesh with a triangulation of its polygons.
    pub mesh: PolygonMesh,
    /// The edges of the cut as vertex pairs `[min, max]` of the new mesh, ascending, each once:
    /// exactly the pieces between consecutive points of the paths.
    pub new_seams: Vec<[u32; 2]>,
    /// For every new polygon, the polygon of the old mesh it is a part of.
    pub polygon_parent: Vec<u32>,
    /// For every new vertex, where it comes from.
    pub vertex_source: Vec<VertexSource>,
    /// For every edge point, its vertex.
    pub edge_point_vertex: Vec<u32>,
    /// For every junction, its vertex.
    pub junction_vertex: Vec<u32>,
    /// For every path, the vertices of its points in order, from the start to the end.
    pub path_vertices: Vec<Vec<u32>>,
    /// Number of vertices of the old mesh.
    pub n_old_vertices: usize,
    /// For every old edge `[min, max]` that received points: its vertices in the new mesh in
    /// order from `min` to `max`, ends included.
    pub edge_chains: BTreeMap<[u32; 2], Vec<u32>>,
}

impl Inserted {
    /// The seams of the new mesh: the given seams of the old mesh (an old edge that received
    /// points is replaced by the pieces of its chain) and the cut. Sorted, each once.
    pub fn seams_after(&self, old_seams: &[[u32; 2]]) -> Vec<[u32; 2]> {
        let mut out: BTreeSet<[u32; 2]> = self.new_seams.iter().copied().collect();
        for &[a, b] in old_seams {
            let key = [a.min(b), a.max(b)];
            match self.edge_chains.get(&key) {
                Some(chain) => {
                    for w in chain.windows(2) {
                        out.insert([w[0].min(w[1]), w[0].max(w[1])]);
                    }
                }
                None => {
                    out.insert(key);
                }
            }
        }
        out.into_iter().collect()
    }
}

/// A vertex of a ring: the new vertex id and the set of lines (sides of the triangle, and the
/// cut chord) it lies on. Three consecutive ring vertices are collinear iff they share a line.
#[derive(Debug, Clone, Copy)]
struct RingVertex {
    id: u32,
    lines: u8,
}

const CHORD: u8 = 8;

/// The better triangulation of a convex ring of vertices on the boundary of a triangle: ears
/// are cut off one by one, never one with three collinear corners; of the ears the one with the
/// best shape (for a quadrilateral: of the two diagonals the one whose worse triangle is better).
fn triangulate_ring(
    ring: &[RingVertex],
    positions: &[[f64; 3]],
    polygon: u32,
) -> Result<Vec<[u32; 3]>, CuttingError> {
    let quality = |a: u32, b: u32, c: u32| -> f64 {
        let (pa, pb, pc) = (
            positions[a as usize],
            positions[b as usize],
            positions[c as usize],
        );
        let sub = |p: [f64; 3], q: [f64; 3]| [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
        let (u, v) = (sub(pb, pa), sub(pc, pa));
        let cross = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let area2 = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        let dot = |x: [f64; 3]| x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        let w = sub(pc, pb);
        let sum = dot(u) + dot(v) + dot(w);
        if sum > 0.0 {
            2.0 * 3f64.sqrt() * area2 / sum
        } else {
            0.0
        }
    };
    let mut rest: Vec<RingVertex> = ring.to_vec();
    let mut out = Vec::with_capacity(ring.len().saturating_sub(2));
    while rest.len() > 3 {
        let m = rest.len();
        let mut best: Option<(f64, usize)> = None;
        for i in 0..m {
            let (a, b, c) = (rest[(i + m - 1) % m], rest[i], rest[(i + 1) % m]);
            if a.lines & b.lines & c.lines != 0 {
                continue;
            }
            let mut score = quality(a.id, b.id, c.id);
            if m == 4 {
                let d = rest[(i + 2) % m];
                if a.lines & c.lines & d.lines != 0 {
                    continue;
                }
                score = score.min(quality(a.id, c.id, d.id));
            }
            if best.is_none_or(|(s, _)| score > s) {
                best = Some((score, i));
            }
        }
        let Some((_, i)) = best else {
            return Err(internal(format!(
                "polygon {polygon}: no ear without three collinear corners in a piece of {m} \
                 corners"
            )));
        };
        let (a, b, c) = (rest[(i + m - 1) % m], rest[i], rest[(i + 1) % m]);
        out.push([a.id, b.id, c.id]);
        rest.remove(i);
    }
    if rest[0].lines & rest[1].lines & rest[2].lines != 0 {
        return Err(internal(format!(
            "polygon {polygon}: a piece of a triangle is left with three collinear corners"
        )));
    }
    out.push([rest[0].id, rest[1].id, rest[2].id]);
    Ok(out)
}

struct Dsu(Vec<u32>);

impl Dsu {
    fn find(&mut self, mut x: u32) -> u32 {
        while self.0[x as usize] != x {
            self.0[x as usize] = self.0[self.0[x as usize] as usize];
            x = self.0[x as usize];
        }
        x
    }
    fn union(&mut self, a: u32, b: u32) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            self.0[a.max(b) as usize] = a.min(b);
        }
    }
}

/// A point of a path, seen from its polygon.
#[derive(Debug, Clone, Copy)]
struct PathPointRef {
    vertex: u32,
    /// The triangulation edge it lies on, as two polygon corners `(lower, higher)`; `None` for
    /// a junction.
    edge: Option<(u32, u32)>,
    /// For a junction, its triangle (an index into the polygon's triangles).
    junction_triangle: Option<usize>,
}

/// A triangle with a junction inside: the vertex `q[k]` of the path on its side `k` (from
/// `corners[k]` to `corners[k + 1]`) and the junction `j`.
#[derive(Debug, Clone, Copy)]
struct TripleTriangle {
    corners: [u32; 3],
    q: [u32; 3],
    j: u32,
}

/// A triangle of the polygon that a path cuts: corners rotated so that the apex is first, the
/// vertex on side `(a, b)` and the vertex on side `(c, a)`.
#[derive(Debug, Clone, Copy)]
struct CutTriangle {
    corners: [u32; 3],
    p: u32,
    q: u32,
}

/// One new polygon before it is written out.
struct Part {
    vertices: Vec<u32>,
    triangles: Vec<[u32; 3]>,
    key: f64,
}

/// Inserts `cut` into `mesh`. See the module documentation of `insert` (the top of this file) for what is built.
///
/// # Errors
/// [`CuttingError::Mesh`] if `mesh` is not a valid polygon mesh with a triangulation,
/// [`CuttingError::InvalidDescription`] if the description does not fit it (an edge point that
/// is not on a polygon edge, a path whose points are not in consecutive triangles of its
/// polygon, ...), [`CuttingError::Unsupported`] for a path that ends at a junction.
pub fn insert(mesh: &PolygonMesh, cut: &CutOnPolygons) -> Result<Inserted, CuttingError> {
    mesh.input().validate()?;
    let n_vertices = mesh.positions.len();
    let n_polygons = mesh.n_polygons();
    let n_corners = mesh.corner_vertices.len();
    let offsets = &mesh.polygon_offsets;
    let mut corner_polygon = vec![0u32; n_corners];
    for p in 0..n_polygons {
        corner_polygon[mesh.polygon_corners(p)].fill(p as u32);
    }
    {
        let mut stamp = vec![NONE; n_vertices];
        for p in 0..n_polygons {
            for &v in mesh.polygon(p) {
                if stamp[v as usize] == p as u32 {
                    return Err(CuttingError::Mesh(
                        crate::UvError::RepeatedVertexInPolygon {
                            polygon: p as u32,
                            vertex: v,
                        },
                    ));
                }
                stamp[v as usize] = p as u32;
            }
        }
    }
    // The polygons of every edge.
    let mut edge_polygons: HashMap<[u32; 2], Vec<u32>> = HashMap::new();
    for p in 0..n_polygons {
        let poly = mesh.polygon(p);
        for (i, &a) in poly.iter().enumerate() {
            let b = poly[(i + 1) % poly.len()];
            edge_polygons
                .entry([a.min(b), a.max(b)])
                .or_default()
                .push(p as u32);
        }
    }

    // ----- edge points -----
    for (i, point) in cut.edge_points.iter().enumerate() {
        let [a, b] = point.edge;
        let polygons = edge_polygons.get(&point.edge).ok_or_else(|| {
            invalid(
                None,
                format!("edge point {i} is on ({a}, {b}), which is not an edge of any polygon (or its vertices are not in ascending order)"),
            )
        })?;
        if !(point.t > 0.0 && point.t < 1.0) {
            return Err(invalid(
                None,
                format!(
                    "edge point {i} on ({a}, {b}) has the parameter {}; it must be strictly between 0 and 1",
                    point.t
                ),
            ));
        }
        if let Some(q) = point.side
            && !polygons.contains(&q)
        {
            return Err(invalid(
                Some(q),
                format!(
                    "edge point {i} on ({a}, {b}) is said to be on the side of polygon {q}, which does not have that edge"
                ),
            ));
        }
    }
    let mut by_edge: BTreeMap<[u32; 2], Vec<u32>> = BTreeMap::new();
    for (i, point) in cut.edge_points.iter().enumerate() {
        by_edge.entry(point.edge).or_default().push(i as u32);
    }
    let mut edge_point_vertex = vec![NONE; cut.edge_points.len()];
    let mut vertex_source: Vec<VertexSource> =
        (0..n_vertices as u32).map(VertexSource::Old).collect();
    let mut chain_of_edge: HashMap<[u32; 2], Vec<u32>> = HashMap::new(); // inner vertices, ascending t
    {
        // Merge points of the two sides with the same parameter; assign ids in point order.
        let mut partner = vec![NONE; cut.edge_points.len()];
        for (edge, list) in &mut by_edge {
            list.sort_by(|&x, &y| {
                let (tx, ty) = (cut.edge_points[x as usize].t, cut.edge_points[y as usize].t);
                tx.partial_cmp(&ty)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then(x.cmp(&y))
            });
            for w in list.windows(2) {
                let (x, y) = (
                    &cut.edge_points[w[0] as usize],
                    &cut.edge_points[w[1] as usize],
                );
                if x.t == y.t {
                    match (x.side, y.side) {
                        (Some(s), Some(u)) if s != u => {
                            partner[w[0] as usize] = w[1];
                            partner[w[1] as usize] = w[0];
                        }
                        _ => {
                            return Err(invalid(
                                None,
                                format!(
                                    "edge points {} and {} are both at t = {} on ({}, {}) \
                                     and not on the two sides of a seam",
                                    w[0], w[1], x.t, edge[0], edge[1]
                                ),
                            ));
                        }
                    }
                }
            }
        }
        for i in 0..cut.edge_points.len() {
            let pt = partner[i];
            if pt != NONE && (pt as usize) < i {
                edge_point_vertex[i] = edge_point_vertex[pt as usize];
            } else {
                edge_point_vertex[i] = vertex_source.len() as u32;
                vertex_source.push(VertexSource::EdgePoint(i as u32));
            }
        }
        for (edge, list) in &by_edge {
            let mut ids: Vec<u32> = list
                .iter()
                .map(|&i| edge_point_vertex[i as usize])
                .collect();
            ids.dedup();
            chain_of_edge.insert(*edge, ids);
        }
    }
    let mut positions = mesh.positions.clone();
    for (i, point) in cut.edge_points.iter().enumerate() {
        if vertex_source[edge_point_vertex[i] as usize] == VertexSource::EdgePoint(i as u32) {
            let [xa, xb] = point.edge.map(|v| mesh.positions[v as usize]);
            positions.push([0, 1, 2].map(|k| (1.0 - point.t) * xa[k] + point.t * xb[k]));
        }
    }

    // ----- junctions -----
    let mut junction_vertex: Vec<u32> = Vec::with_capacity(cut.junctions.len());
    for (j, junction) in cut.junctions.iter().enumerate() {
        let pol = junction.polygon;
        let bad = |problem: String| invalid(Some(pol), format!("junction {j}: {problem}"));
        if pol as usize >= n_polygons {
            return Err(invalid(
                None,
                format!("junction {j} is in polygon {pol}, but there are {n_polygons} polygons"),
            ));
        }
        let Some(triangle) = mesh.triangles.get(junction.triangle as usize) else {
            return Err(bad(format!(
                "the triangle {} does not exist ({} triangles)",
                junction.triangle,
                mesh.triangles.len()
            )));
        };
        if *triangle != junction.corners || corner_polygon[triangle[0] as usize] != pol {
            return Err(bad(format!(
                "triangle {} has the corners {:?} in polygon {}, not {:?} in polygon {pol}",
                junction.triangle, triangle, corner_polygon[triangle[0] as usize], junction.corners
            )));
        }
        let w = junction.barycentric;
        if !(w.iter().all(|x| x.is_finite() && *x > 0.0)
            && (w.iter().sum::<f64>() - 1.0).abs() <= 1e-9)
        {
            return Err(bad(format!(
                "the barycentric coordinates {w:?} must be positive and sum to 1 (to 1e-9)"
            )));
        }
        junction_vertex.push(vertex_source.len() as u32);
        vertex_source.push(VertexSource::Junction(j as u32));
        let x = junction
            .corners
            .map(|c| mesh.positions[mesh.corner_vertices[c as usize] as usize]);
        positions.push([0, 1, 2].map(|d| w[0] * x[0][d] + w[1] * x[1][d] + w[2] * x[2][d]));
    }

    // ----- paths -----
    let mut path_vertices: Vec<Vec<u32>> = Vec::with_capacity(cut.paths.len());
    let mut paths_of_polygon: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (k, path) in cut.paths.iter().enumerate() {
        let pol = path.polygon;
        if pol as usize >= n_polygons {
            return Err(invalid(
                None,
                format!("path {k} is in polygon {pol}, but there are {n_polygons} polygons"),
            ));
        }
        let mut ids = Vec::with_capacity(path.interior.len() + 2);
        for (name, end) in [("start", path.start), ("end", path.end)] {
            match end {
                PathEnd::EdgePoint(i) => {
                    let point = cut.edge_points.get(i as usize).ok_or_else(|| {
                        invalid(
                            Some(pol),
                            format!("path {k} {name}s at edge point {i}, which does not exist"),
                        )
                    })?;
                    let ok = edge_polygons
                        .get(&point.edge)
                        .is_some_and(|ps| ps.contains(&pol))
                        && point.side.is_none_or(|s| s == pol);
                    if !ok {
                        return Err(invalid(
                            Some(pol),
                            format!(
                                "path {k} {name}s at edge point {i} on ({}, {}), which is not \
                                 on an edge of this polygon (or on the side of another polygon)",
                                point.edge[0], point.edge[1]
                            ),
                        ));
                    }
                    if name == "start" {
                        ids.push(edge_point_vertex[i as usize]);
                    }
                }
                PathEnd::Junction(i) => {
                    let junction = cut.junctions.get(i as usize).ok_or_else(|| {
                        invalid(
                            Some(pol),
                            format!("path {k} {name}s at junction {i}, which does not exist"),
                        )
                    })?;
                    if junction.polygon != pol {
                        return Err(invalid(
                            Some(pol),
                            format!(
                                "path {k} {name}s at junction {i}, which is in polygon {}",
                                junction.polygon
                            ),
                        ));
                    }
                    if name == "start" {
                        ids.push(junction_vertex[i as usize]);
                    }
                }
            }
        }
        // interior ids are assigned below, the end id after them
        path_vertices.push(ids);
        paths_of_polygon.entry(pol).or_default().push(k);
    }
    for (k, path) in cut.paths.iter().enumerate() {
        for (j, point) in path.interior.iter().enumerate() {
            let [c0, c1] = point.corners;
            if !(c0 < c1
                && (c1 as usize) < n_corners
                && corner_polygon[c0 as usize] == path.polygon
                && corner_polygon[c1 as usize] == path.polygon)
            {
                return Err(invalid(
                    Some(path.polygon),
                    format!(
                        "path {k}, interior point {j}: the corners ({c0}, {c1}) must be two \
                         ascending corners of polygon {}",
                        path.polygon
                    ),
                ));
            }
            if !(point.t > 0.0 && point.t < 1.0) {
                return Err(invalid(
                    Some(path.polygon),
                    format!(
                        "path {k}, interior point {j} has the parameter {}; it must be strictly between 0 and 1",
                        point.t
                    ),
                ));
            }
            path_vertices[k].push(vertex_source.len() as u32);
            vertex_source.push(VertexSource::Interior {
                path: k as u32,
                index: j as u32,
            });
            let [xa, xb] =
                [c0, c1].map(|c| mesh.positions[mesh.corner_vertices[c as usize] as usize]);
            positions.push([0, 1, 2].map(|d| (1.0 - point.t) * xa[d] + point.t * xb[d]));
        }
        match path.end {
            PathEnd::EdgePoint(i) => path_vertices[k].push(edge_point_vertex[i as usize]),
            PathEnd::Junction(i) => path_vertices[k].push(junction_vertex[i as usize]),
        }
    }

    // ----- which polygons change -----
    let mut touched = vec![false; n_polygons];
    for &p in paths_of_polygon.keys() {
        touched[p as usize] = true;
    }
    for edge in chain_of_edge.keys() {
        for &p in &edge_polygons[edge] {
            touched[p as usize] = true;
        }
    }
    let mut triangles_of_polygon: Vec<Vec<[u32; 3]>> = vec![Vec::new(); n_polygons];
    for t in &mesh.triangles {
        let p = corner_polygon[t[0] as usize] as usize;
        if touched[p] {
            triangles_of_polygon[p].push(*t);
        }
    }

    // ----- polygon by polygon -----
    let mut out = PolygonMesh {
        positions: Vec::new(),
        corner_vertices: Vec::with_capacity(n_corners + 2 * cut.edge_points.len()),
        polygon_offsets: vec![0],
        triangles: Vec::with_capacity(mesh.triangles.len() + 2 * cut.paths.len()),
    };
    let mut polygon_parent: Vec<u32> = Vec::with_capacity(n_polygons);
    let ctx = PolygonContext {
        mesh,
        cut,
        chain_of_edge: &chain_of_edge,
        positions: &positions,
        path_vertices: &path_vertices,
    };
    // For every old polygon, the first corner of its first new polygon.
    let mut new_start = vec![0u32; n_polygons];
    for p in 0..n_polygons {
        new_start[p] = out.corner_vertices.len() as u32;
        if !touched[p] {
            out.corner_vertices.extend_from_slice(mesh.polygon(p));
            out.polygon_offsets.push(out.corner_vertices.len() as u32);
            polygon_parent.push(p as u32);
            continue;
        }
        let parts = ctx.process(
            p as u32,
            &triangles_of_polygon[p],
            paths_of_polygon
                .get(&(p as u32))
                .map_or(&[][..], |v| &v[..]),
        )?;
        for part in parts {
            let start = out.corner_vertices.len() as u32;
            let local: HashMap<u32, u32> = part
                .vertices
                .iter()
                .enumerate()
                .map(|(i, &v)| (v, start + i as u32))
                .collect();
            out.corner_vertices.extend_from_slice(&part.vertices);
            out.polygon_offsets.push(out.corner_vertices.len() as u32);
            polygon_parent.push(p as u32);
            for t in &part.triangles {
                let mut mapped = [0u32; 3];
                for (k, v) in t.iter().enumerate() {
                    mapped[k] = *local.get(v).ok_or_else(|| {
                        internal(format!(
                            "polygon {p}: a triangle uses vertex {v}, which is not on the \
                             boundary of its part"
                        ))
                    })?;
                }
                out.triangles.push(mapped);
            }
        }
    }
    // Untouched polygons keep their triangulation (corner indices shifted).
    for t in &mesh.triangles {
        let p = corner_polygon[t[0] as usize] as usize;
        if !touched[p] {
            let shift = |c: u32| new_start[p] + (c - offsets[p]);
            out.triangles.push([shift(t[0]), shift(t[1]), shift(t[2])]);
        }
    }
    out.positions = positions;

    // ----- seams, chains -----
    let mut seams: BTreeSet<[u32; 2]> = BTreeSet::new();
    for ids in &path_vertices {
        for w in ids.windows(2) {
            seams.insert([w[0].min(w[1]), w[0].max(w[1])]);
        }
    }
    let mut edge_chains = BTreeMap::new();
    for (edge, inner) in &chain_of_edge {
        let mut chain = vec![edge[0]];
        chain.extend_from_slice(inner);
        chain.push(edge[1]);
        edge_chains.insert(*edge, chain);
    }
    let inserted = Inserted {
        mesh: out,
        new_seams: seams.into_iter().collect(),
        polygon_parent,
        vertex_source,
        edge_point_vertex,
        junction_vertex,
        path_vertices,
        n_old_vertices: n_vertices,
        edge_chains,
    };
    inserted
        .mesh
        .input()
        .validate()
        .map_err(|e| internal(format!("the new mesh is not a valid polygon mesh: {e}")))?;
    Ok(inserted)
}

struct PolygonContext<'a> {
    mesh: &'a PolygonMesh,
    cut: &'a CutOnPolygons,
    chain_of_edge: &'a HashMap<[u32; 2], Vec<u32>>,
    positions: &'a [[f64; 3]],
    path_vertices: &'a [Vec<u32>],
}

impl PolygonContext<'_> {
    /// Splits polygon `p`; returns its parts, ordered by the first boundary position.
    fn process(
        &self,
        p: u32,
        triangles: &[[u32; 3]],
        paths: &[usize],
    ) -> Result<Vec<Part>, CuttingError> {
        let mesh = self.mesh;
        let range = mesh.polygon_corners(p as usize);
        let (lo, n) = (range.start as u32, range.len() as u32);
        let verts: Vec<u32> = mesh.polygon(p as usize).to_vec();
        let next = |c: u32| lo + (c - lo + 1) % n;
        let is_polygon_edge = |c: u32, d: u32| next(c) == d || next(d) == c;
        // Vertices strictly between two consecutive corners, in the direction c -> d.
        let extras = |c: u32, d: u32| -> Vec<u32> {
            let (a, b) = (verts[(c - lo) as usize], verts[(d - lo) as usize]);
            match self.chain_of_edge.get(&[a.min(b), a.max(b)]) {
                None => Vec::new(),
                Some(inner) => {
                    if a < b {
                        inner.clone()
                    } else {
                        inner.iter().rev().copied().collect()
                    }
                }
            }
        };
        // Triangulation edges to triangles.
        let mut tri_of_edge: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (ti, t) in triangles.iter().enumerate() {
            for k in 0..3 {
                let (c, d) = (t[k], t[(k + 1) % 3]);
                tri_of_edge
                    .entry((c.min(d), c.max(d)))
                    .or_default()
                    .push(ti);
            }
        }
        // The polygon edge (as corners) of an edge point.
        let edge_corners = |edge: [u32; 2]| -> Option<(u32, u32)> {
            (0..n).find_map(|i| {
                let (c, d) = (lo + i, next(lo + i));
                let (a, b) = (verts[i as usize], verts[(d - lo) as usize]);
                (a.min(b) == edge[0] && a.max(b) == edge[1]).then_some((c.min(d), c.max(d)))
            })
        };

        // Segments: which triangle each path crosses between two consecutive points.
        let mut cut_tri: Vec<Option<CutTriangle>> = vec![None; triangles.len()];
        let mut triple_tri: Vec<Option<TripleTriangle>> = vec![None; triangles.len()];
        for &k in paths {
            let path = &self.cut.paths[k];
            let ids = &self.path_vertices[k];
            let mut points: Vec<PathPointRef> = Vec::with_capacity(ids.len());
            for (j, &vertex) in ids.iter().enumerate() {
                let mut junction_triangle = None;
                let edge = if j == 0 || j + 1 == ids.len() {
                    let end = if j == 0 { path.start } else { path.end };
                    match end {
                        PathEnd::EdgePoint(i) => Some(
                            edge_corners(self.cut.edge_points[i as usize].edge).ok_or_else(
                                || {
                                    invalid(
                                        Some(p),
                                        format!(
                                            "path {k}: edge point {i} is not on an edge of the polygon"
                                        ),
                                    )
                                },
                            )?,
                        ),
                        PathEnd::Junction(i) => {
                            let corners = self.cut.junctions[i as usize].corners;
                            junction_triangle = triangles.iter().position(|t| *t == corners);
                            if junction_triangle.is_none() {
                                return Err(invalid(
                                    Some(p),
                                    format!("path {k}: the triangle of junction {i} is not a triangle of the polygon"),
                                ));
                            }
                            None
                        }
                    }
                } else {
                    let [c0, c1] = path.interior[j - 1].corners;
                    if !tri_of_edge.contains_key(&(c0, c1)) || is_polygon_edge(c0, c1) {
                        return Err(invalid(
                            Some(p),
                            format!(
                                "path {k}, interior point {}: ({c0}, {c1}) is not a diagonal of \
                                 the triangulation of the polygon",
                                j - 1
                            ),
                        ));
                    }
                    Some((c0, c1))
                };
                points.push(PathPointRef {
                    vertex,
                    edge,
                    junction_triangle,
                });
            }
            for w in points.windows(2) {
                let (x, y) = (w[0], w[1]);
                // A segment that ends at a junction: one of the three arms of its triangle.
                if x.junction_triangle.is_some() || y.junction_triangle.is_some() {
                    let (jp, other) = if x.junction_triangle.is_some() {
                        (x, y)
                    } else {
                        (y, x)
                    };
                    let (Some(ti), Some(edge)) = (jp.junction_triangle, other.edge) else {
                        return Err(invalid(
                            Some(p),
                            format!("path {k}: a segment joins two junctions"),
                        ));
                    };
                    let t = triangles[ti];
                    let side = (0..3).find(|&s| {
                        let (c, d) = (t[s], t[(s + 1) % 3]);
                        edge == (c.min(d), c.max(d))
                    });
                    let Some(side) = side else {
                        return Err(invalid(
                            Some(p),
                            format!(
                                "path {k}: the point next to the junction is on the edge {edge:?}, which is not a side of the junction's triangle {t:?}"
                            ),
                        ));
                    };
                    let entry = triple_tri[ti].get_or_insert(TripleTriangle {
                        corners: t,
                        q: [NONE; 3],
                        j: jp.vertex,
                    });
                    if entry.j != jp.vertex || entry.q[side] != NONE {
                        return Err(invalid(
                            Some(p),
                            format!(
                                "path {k}: the junction's triangle {t:?} has two arms through its side {side}"
                            ),
                        ));
                    }
                    entry.q[side] = other.vertex;
                    continue;
                }
                let (xe, ye) = (
                    x.edge.expect("not a junction"),
                    y.edge.expect("not a junction"),
                );
                if xe == ye {
                    return Err(invalid(
                        Some(p),
                        format!("path {k}: two consecutive points lie on the same edge {xe:?}"),
                    ));
                }
                let tx = &tri_of_edge[&xe];
                let ty = &tri_of_edge[&ye];
                let common: Vec<usize> = tx.iter().copied().filter(|t| ty.contains(t)).collect();
                if common.len() != 1 {
                    return Err(invalid(
                        Some(p),
                        format!(
                            "path {k}: consecutive points on the edges {xe:?} and {ye:?} do not \
                             share a triangle of the triangulation"
                        ),
                    ));
                }
                let ti = common[0];
                if cut_tri[ti].is_some() || triple_tri[ti].is_some() {
                    return Err(invalid(
                        Some(p),
                        format!(
                            "path {k}: triangle {:?} is crossed by more than one segment",
                            triangles[ti]
                        ),
                    ));
                }
                let t = triangles[ti];
                let apex = (0..3)
                    .find(|&i| {
                        let c = t[i];
                        [xe.0, xe.1].contains(&c) && [ye.0, ye.1].contains(&c)
                    })
                    .ok_or_else(|| internal("two edges of a triangle without a common corner"))?;
                let arr = [t[apex], t[(apex + 1) % 3], t[(apex + 2) % 3]];
                let side_ab = (arr[0].min(arr[1]), arr[0].max(arr[1]));
                let (px, qy) = if xe == side_ab { (x, y) } else { (y, x) };
                cut_tri[ti] = Some(CutTriangle {
                    corners: arr,
                    p: px.vertex,
                    q: qy.vertex,
                });
            }
        }
        for (ti, entry) in triple_tri.iter().enumerate() {
            if let Some(e) = entry
                && (e.q.contains(&NONE) || cut_tri[ti].is_some())
            {
                return Err(invalid(
                    Some(p),
                    format!(
                        "the junction in triangle {:?} needs three arms through its three sides and no other segment in the triangle",
                        triangles[ti]
                    ),
                ));
            }
        }

        // Nodes: the pieces of triangles.
        let mut base = vec![0u32; triangles.len()];
        let mut n_nodes = 0u32;
        for (ti, c) in cut_tri.iter().enumerate() {
            base[ti] = n_nodes;
            n_nodes += if c.is_some() {
                2
            } else if triple_tri[ti].is_some() {
                3
            } else {
                1
            };
        }
        let mut rings: Vec<Vec<RingVertex>> = vec![Vec::new(); n_nodes as usize];
        let corner_vertex = |c: u32| verts[(c - lo) as usize];
        let ring_extras = |ring: &mut Vec<RingVertex>, list: &[u32], line: u8| {
            for &id in list {
                ring.push(RingVertex { id, lines: line });
            }
        };
        for (ti, t) in triangles.iter().enumerate() {
            if let Some(tt) = triple_tri[ti] {
                let sides: Vec<Vec<u32>> = (0..3)
                    .map(|k| {
                        let (c, d) = (tt.corners[k], tt.corners[(k + 1) % 3]);
                        if is_polygon_edge(c, d) {
                            extras(c, d)
                        } else {
                            Vec::new()
                        }
                    })
                    .collect();
                // Position of q[k] in its side's list (MAX: the side is a diagonal).
                let mut at = [usize::MAX; 3];
                for k in 0..3 {
                    if !sides[k].is_empty() {
                        at[k] = sides[k].iter().position(|&x| x == tt.q[k]).ok_or_else(|| {
                            invalid(
                                Some(p),
                                format!(
                                    "the cut point of vertex {} is not on polygon edge {k} of triangle {t:?}",
                                    tt.q[k]
                                ),
                            )
                        })?;
                    }
                }
                for k in 0..3 {
                    let prev = (k + 2) % 3;
                    let ring = &mut rings[base[ti] as usize + k];
                    ring.push(RingVertex {
                        id: corner_vertex(tt.corners[k]),
                        lines: (1 << k) | (1 << prev),
                    });
                    if at[k] != usize::MAX {
                        ring_extras(ring, &sides[k][..at[k]], 1 << k);
                    }
                    ring.push(RingVertex {
                        id: tt.q[k],
                        lines: 1 << k,
                    });
                    ring.push(RingVertex { id: tt.j, lines: 0 });
                    ring.push(RingVertex {
                        id: tt.q[prev],
                        lines: 1 << prev,
                    });
                    if at[prev] != usize::MAX {
                        ring_extras(ring, &sides[prev][at[prev] + 1..], 1 << prev);
                    }
                }
                continue;
            }
            match cut_tri[ti] {
                None => {
                    let node = base[ti] as usize;
                    let ring = &mut rings[node];
                    for k in 0..3 {
                        ring.push(RingVertex {
                            id: corner_vertex(t[k]),
                            lines: (1 << k) | (1 << ((k + 2) % 3)),
                        });
                        let (c, d) = (t[k], t[(k + 1) % 3]);
                        if is_polygon_edge(c, d) {
                            ring_extras(ring, &extras(c, d), 1 << k);
                        }
                    }
                }
                Some(ct) => {
                    let [a, b, c] = ct.corners;
                    let (l0, l1, l2) = (
                        if is_polygon_edge(a, b) {
                            extras(a, b)
                        } else {
                            Vec::new()
                        },
                        if is_polygon_edge(b, c) {
                            extras(b, c)
                        } else {
                            Vec::new()
                        },
                        if is_polygon_edge(c, a) {
                            extras(c, a)
                        } else {
                            Vec::new()
                        },
                    );
                    let split = |list: &[u32],
                                 id: u32,
                                 side: &str|
                     -> Result<usize, CuttingError> {
                        if list.is_empty() {
                            return Ok(usize::MAX);
                        }
                        list.iter().position(|&x| x == id).ok_or_else(|| {
                            invalid(
                                Some(p),
                                format!("the cut point of vertex {id} is not on the {side} polygon edge of triangle {t:?}"),
                            )
                        })
                    };
                    let ip = split(&l0, ct.p, "first")?;
                    let iq = split(&l2, ct.q, "last")?;
                    let (before_p, after_p): (&[u32], &[u32]) = if ip == usize::MAX {
                        (&[], &[])
                    } else {
                        (&l0[..ip], &l0[ip + 1..])
                    };
                    let (before_q, after_q): (&[u32], &[u32]) = if iq == usize::MAX {
                        (&[], &[])
                    } else {
                        (&l2[..iq], &l2[iq + 1..])
                    };
                    // The part of the apex a: a, side (a,b) up to p, chord to q, side (c,a) from q.
                    let ring = &mut rings[base[ti] as usize];
                    ring.push(RingVertex {
                        id: corner_vertex(a),
                        lines: 1 | 4,
                    });
                    ring_extras(ring, before_p, 1);
                    ring.push(RingVertex {
                        id: ct.p,
                        lines: 1 | CHORD,
                    });
                    ring.push(RingVertex {
                        id: ct.q,
                        lines: 4 | CHORD,
                    });
                    ring_extras(ring, after_q, 4);
                    // The quadrilateral: p, side to b, side (b,c), side (c,a) up to q.
                    let ring = &mut rings[base[ti] as usize + 1];
                    ring.push(RingVertex {
                        id: ct.p,
                        lines: 1 | CHORD,
                    });
                    ring_extras(ring, after_p, 1);
                    ring.push(RingVertex {
                        id: corner_vertex(b),
                        lines: 1 | 2,
                    });
                    ring_extras(ring, &l1, 2);
                    ring.push(RingVertex {
                        id: corner_vertex(c),
                        lines: 2 | 4,
                    });
                    ring_extras(ring, before_q, 4);
                    ring.push(RingVertex {
                        id: ct.q,
                        lines: 4 | CHORD,
                    });
                }
            }
        }

        // Which node holds the part of the triangle next to a corner of an edge.
        let node_near = |ti: usize, edge: (u32, u32), near: u32| -> u32 {
            if let Some(tt) = triple_tri[ti] {
                for k in 0..3 {
                    let (c, d) = (tt.corners[k], tt.corners[(k + 1) % 3]);
                    if edge == (c.min(d), c.max(d)) {
                        return base[ti]
                            + if near == c {
                                k as u32
                            } else {
                                ((k + 1) % 3) as u32
                            };
                    }
                }
            }
            match cut_tri[ti] {
                None => base[ti],
                Some(ct) => {
                    let [a, b, c] = ct.corners;
                    let has = |x: u32, y: u32| edge == (x.min(y), x.max(y));
                    if (has(a, b) || has(c, a)) && near == a {
                        base[ti]
                    } else {
                        base[ti] + 1
                    }
                }
            }
        };
        // The vertex a segment puts on an edge of a triangle, if it crosses that edge.
        let crossing = |ti: usize, edge: (u32, u32)| -> Option<u32> {
            if let Some(tt) = triple_tri[ti] {
                return (0..3).find_map(|k| {
                    let (c, d) = (tt.corners[k], tt.corners[(k + 1) % 3]);
                    (edge == (c.min(d), c.max(d))).then_some(tt.q[k])
                });
            }
            let ct = cut_tri[ti]?;
            let [a, b, c] = ct.corners;
            let has = |x: u32, y: u32| edge == (x.min(y), x.max(y));
            if has(a, b) {
                Some(ct.p)
            } else if has(c, a) {
                Some(ct.q)
            } else {
                None
            }
        };
        let mut dsu = Dsu((0..n_nodes).collect());
        for (&edge, tris) in &tri_of_edge {
            if tris.len() == 2 && !is_polygon_edge(edge.0, edge.1) {
                let (t1, t2) = (tris[0], tris[1]);
                if crossing(t1, edge) != crossing(t2, edge) {
                    return Err(invalid(
                        Some(p),
                        format!(
                            "the diagonal {edge:?} is crossed at vertex {:?} in one triangle and \
                             at vertex {:?} in the other",
                            crossing(t1, edge),
                            crossing(t2, edge)
                        ),
                    ));
                }
                for near in [edge.0, edge.1] {
                    dsu.union(node_near(t1, edge, near), node_near(t2, edge, near));
                }
            } else if tris.len() == 1 && !is_polygon_edge(edge.0, edge.1) {
                return Err(invalid(
                    Some(p),
                    format!("the diagonal {edge:?} is in one triangle only"),
                ));
            }
        }

        // Triangulate the nodes.
        let mut node_triangles: Vec<Vec<[u32; 3]>> = Vec::with_capacity(rings.len());
        for ring in &rings {
            node_triangles.push(triangulate_ring(ring, self.positions, p)?);
        }

        // Keys along the boundary of the polygon.
        let mut key: HashMap<u32, f64> = HashMap::new();
        for i in 0..n {
            let (c, d) = (lo + i, next(lo + i));
            key.insert(verts[i as usize], f64::from(i));
            let list = extras(c, d);
            let m = list.len() as f64;
            for (j, id) in list.iter().enumerate() {
                key.entry(*id)
                    .or_insert(f64::from(i) + (j as f64 + 1.0) / (m + 1.0));
            }
        }

        // The parts.
        let mut groups: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
        for node in 0..n_nodes {
            groups
                .entry(dsu.find(node))
                .or_default()
                .push(node as usize);
        }
        let mut parts = Vec::with_capacity(groups.len());
        for nodes in groups.values() {
            let mut directed: HashSet<(u32, u32)> = HashSet::new();
            for &node in nodes {
                let ring = &rings[node];
                for i in 0..ring.len() {
                    let (u, v) = (ring[i].id, ring[(i + 1) % ring.len()].id);
                    if !directed.remove(&(v, u)) && !directed.insert((u, v)) {
                        return Err(internal(format!(
                            "polygon {p}: the side {u} -> {v} occurs twice"
                        )));
                    }
                }
            }
            let mut next_of: HashMap<u32, u32> = HashMap::with_capacity(directed.len());
            for &(u, v) in &directed {
                if next_of.insert(u, v).is_some() {
                    return Err(invalid(
                        Some(p),
                        format!("a part of the polygon touches itself at vertex {u}"),
                    ));
                }
            }
            let start = *next_of
                .keys()
                .min_by(|&&x, &&y| {
                    let (kx, ky) = (
                        key.get(&x).copied().unwrap_or(f64::INFINITY),
                        key.get(&y).copied().unwrap_or(f64::INFINITY),
                    );
                    kx.partial_cmp(&ky)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then(x.cmp(&y))
                })
                .ok_or_else(|| internal(format!("polygon {p}: a part without a boundary")))?;
            let mut vertices = vec![start];
            let mut v = next_of[&start];
            while v != start {
                vertices.push(v);
                v = *next_of.get(&v).ok_or_else(|| {
                    internal(format!(
                        "polygon {p}: the boundary of a part is open at vertex {v}"
                    ))
                })?;
                if vertices.len() > next_of.len() {
                    return Err(internal(format!(
                        "polygon {p}: the boundary of a part does not close"
                    )));
                }
            }
            if vertices.len() != next_of.len() {
                return Err(invalid(
                    Some(p),
                    "a part of the polygon has more than one boundary loop (the paths do not cut the polygon into discs)",
                ));
            }
            let triangles: Vec<[u32; 3]> = nodes
                .iter()
                .flat_map(|&n| node_triangles[n].iter().copied())
                .collect();
            let part_key = vertices
                .iter()
                .map(|v| key.get(v).copied().unwrap_or(f64::INFINITY))
                .fold(f64::INFINITY, f64::min);
            parts.push(Part {
                vertices,
                triangles,
                key: part_key,
            });
        }
        parts.sort_by(|a, b| {
            a.key
                .partial_cmp(&b.key)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        Ok(parts)
    }
}
