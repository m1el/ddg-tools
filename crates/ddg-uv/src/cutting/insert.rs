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
//! * Every [`JunctionPoint`](super::JunctionPoint) becomes a vertex inside its triangle (see
//!   below).
//! * Old vertices keep their numbers; new vertices follow (edge points in the order of
//!   [`CutOnPolygons::edge_points`], then junctions in the order of
//!   [`CutOnPolygons::junctions`], then interior points in the order of the paths). Polygons
//!   are numbered in the order of their parents; the parts of one polygon follow each other,
//!   starting with the part that has the first corner of the parent.
//!
//! # The construction, exactly (what an editor must mirror)
//!
//! Work polygon by polygon; a polygon that no path, junction or extra edge vertex touches is
//! copied as it is (corners and triangulation).
//!
//! 1. **Vertices on polygon edges.** For every polygon edge `{a, b}` collect the vertices of
//!    all edge points on it, sorted by the parameter `t` from the lower vertex (the points of
//!    the two sides of a seam together; equal parameters from the two sides are one vertex).
//!    Both polygons of the edge get all of them as corners between `a` and `b`.
//! 2. **Pieces of triangles.** Take the polygon's triangles. A triangle `T = (c0, c1, c2)`
//!    (polygon corners; side `k` runs from `c_k` to `c_{k+1}`) is
//!    * *untouched* if no path segment crosses it: one piece, `T` itself;
//!    * *cut by a segment* from the vertex `p` on a side to the vertex `q` on another side (the
//!      two sides share a corner `a`; walk `T` as `(a, b, c)` with `p` on `(a, b)` and `q` on
//!      `(c, a)`): two pieces, the corner piece `(a, p, q)` and the quadrilateral
//!      `(p, b, c, q)`;
//!    * *a junction triangle*: the triangle that contains a junction `J` (its position is
//!      `Σ barycentric[k] · x[corners[k]]`). Exactly three paths end at `J`, each as the
//!      segment from the vertex `q_k` on side `k` (from `c_k` to `c_{k+1}`) to `J`. Three
//!      pieces, **piece `k` is the quadrilateral `(c_k, q_k, J, q_{k−1})`** (indices modulo 3),
//!      the part of corner `c_k`. `J` is a vertex of those three pieces only.
//!
//!    The vertex of a path on a side is, on a polygon edge, an edge-point vertex and, on a
//!    diagonal of the triangulation, the interior point of the path there. Every piece lists,
//!    in order along its sides, the vertices from step 1 that lie on its polygon-edge sides
//!    (for a side crossed at `p`: those before `p` go to the piece of the first corner, those
//!    after `p` to the other piece).
//! 3. **Parts.** Two pieces on the two sides of an uncrossed diagonal belong to the same part;
//!    on a diagonal crossed at a vertex `p`, the pieces next to each end of the diagonal are
//!    joined (the piece of the first triangle that contains the diagonal's end `c` with the
//!    piece of the second triangle that contains `c`). A part is a connected set of pieces.
//!    Its polygon is the boundary: the sides of its pieces that are not shared with another
//!    piece of the part, walked in the winding of the polygon; it starts at the vertex with the
//!    lowest position along the boundary of the parent (the parent's first corner, then the
//!    edge vertices after it, ...), and a part without such a vertex starts at the lowest id.
//!    Parts are ordered by that start.
//! 4. **Triangulation of a part.** Each piece is triangulated by ear clipping (never an ear
//!    with three collinear vertices): of the ears the one whose triangle has the best shape
//!    (`4√3 · area / Σ squared lengths`) is clipped first; for a quadrilateral the diagonal
//!    whose worse triangle is better. The triangles of all pieces of a part are the
//!    triangulation of its polygon.
//!
//! The seams are the sides between consecutive points of each path (including the sides that
//! end at a junction); `Inserted::seams_after` adds the old seams.
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
    /// For every new polygon, the region it lies in where the description says so: a part
    /// next to a piece of the cut gets the region on its side of the piece (the region with
    /// the lower id is on the left of a path). `None` for a polygon the cut does not touch
    /// (the description does not know its region) and for the parts of a polygon that only
    /// received vertices.
    pub polygon_region: Vec<Option<u32>>,
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

/// One new polygon before it is written out.
struct Part {
    vertices: Vec<u32>,
    triangles: Vec<[u32; 3]>,
    key: f64,
    /// The region the part lies in, from the sides of the pieces of the cut on its boundary.
    region: Option<u32>,
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

    // The parameter of every vertex on a polygon edge or on a diagonal, from the lower end.
    let mut vertex_t: HashMap<u32, f64> = HashMap::new();
    for (i, point) in cut.edge_points.iter().enumerate() {
        vertex_t.insert(edge_point_vertex[i], point.t);
    }
    for (k, path) in cut.paths.iter().enumerate() {
        for (j, point) in path.interior.iter().enumerate() {
            vertex_t.insert(path_vertices[k][j + 1], point.t);
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
    let mut polygon_region: Vec<Option<u32>> = Vec::with_capacity(n_polygons);
    let ctx = PolygonContext {
        mesh,
        cut,
        chain_of_edge: &chain_of_edge,
        positions: &positions,
        path_vertices: &path_vertices,
        vertex_t: &vertex_t,
    };
    // For every old polygon, the first corner of its first new polygon.
    let mut new_start = vec![0u32; n_polygons];
    for p in 0..n_polygons {
        new_start[p] = out.corner_vertices.len() as u32;
        if !touched[p] {
            out.corner_vertices.extend_from_slice(mesh.polygon(p));
            out.polygon_offsets.push(out.corner_vertices.len() as u32);
            polygon_parent.push(p as u32);
            polygon_region.push(None);
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
            polygon_region.push(part.region);
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
    let mut edge_chains = BTreeMap::new();
    for (edge, inner) in &chain_of_edge {
        let mut chain = vec![edge[0]];
        chain.extend_from_slice(inner);
        chain.push(edge[1]);
        edge_chains.insert(*edge, chain);
    }
    let mut seams: BTreeSet<[u32; 2]> = BTreeSet::new();
    for ids in &path_vertices {
        for w in ids.windows(2) {
            let (u, v) = (w[0], w[1]);
            // A piece along a polygon edge covers every side of the chain between its ends.
            let along = edge_chains.values().find_map(|chain| {
                let (i, j) = (
                    chain.iter().position(|&x| x == u)?,
                    chain.iter().position(|&x| x == v)?,
                );
                Some(chain[i.min(j)..=i.max(j)].to_vec())
            });
            match along {
                Some(run) if run.len() > 2 => {
                    for w in run.windows(2) {
                        seams.insert([w[0].min(w[1]), w[0].max(w[1])]);
                    }
                }
                _ => {
                    seams.insert([u.min(v), u.max(v)]);
                }
            }
        }
    }
    let inserted = Inserted {
        mesh: out,
        new_seams: seams.into_iter().collect(),
        polygon_parent,
        polygon_region,
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
    /// The parameter of every vertex on a polygon edge or on a diagonal, from the lower end.
    vertex_t: &'a HashMap<u32, f64>,
}

/// Where a point of a path lies, seen from its polygon.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Loc {
    /// On a side of triangles: a polygon edge or a diagonal, as two polygon corners
    /// `(lower, higher)`.
    Side((u32, u32)),
    /// Inside the triangle with this index (into the polygon's triangles).
    Inner(usize),
}

/// A straight piece of the cut inside one triangle.
#[derive(Debug, Clone, Copy)]
struct PieceRec {
    a: u32,
    b: u32,
    tri: usize,
    /// The regions on the left and on the right of `a -> b`.
    left: u32,
    right: u32,
}

/// What a side of a cell is.
#[derive(Debug, Clone, Copy, PartialEq)]
enum SideKind {
    /// A part of a side of the triangle; `diag` for a diagonal.
    Boundary { diag: bool },
    /// A piece of the cut, walked in its direction (`forward`) or against it.
    Piece { rec: usize, forward: bool },
}

/// A face of the subdivision of one triangle by the pieces of the cut.
struct Cell {
    ids: Vec<u32>,
    xy: Vec<[f64; 2]>,
    sides: Vec<SideKind>,
    region: Option<u32>,
}

const EAR_EPS: f64 = 1e-12;

fn cross2(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Triangulates a simple polygon in the plane of a triangle (counter-clockwise `xy`, in
/// barycentric coordinates, which preserve straightness and orientation) by ear clipping: an ear
/// is a convex corner whose triangle holds no other vertex; of the ears the one with the best
/// shape in space (`4√3 · area / Σ squared lengths`) is clipped first, and for a
/// quadrilateral the one whose worse triangle is better. Collinear corners are never ears.
fn triangulate_cell(
    ids: &[u32],
    xy: &[[f64; 2]],
    positions: &[[f64; 3]],
    polygon: u32,
) -> Result<Vec<[u32; 3]>, CuttingError> {
    let quality = |a: usize, b: usize, c: usize| -> f64 {
        let p = [a, b, c].map(|i| positions[ids[i] as usize]);
        let sub = |x: [f64; 3], y: [f64; 3]| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
        let (u, v, w) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[2], p[1]));
        let cr = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let area2 = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
        let sq = |x: [f64; 3]| x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        let sum = sq(u) + sq(v) + sq(w);
        if sum > 0.0 {
            2.0 * 3f64.sqrt() * area2 / sum
        } else {
            0.0
        }
    };
    let n = ids.len();
    let mut rest: Vec<usize> = (0..n).collect();
    let mut out = Vec::with_capacity(n.saturating_sub(2));
    while rest.len() > 3 {
        let m = rest.len();
        let mut best: Option<(f64, usize)> = None;
        for i in 0..m {
            let (a, b, c) = (rest[(i + m - 1) % m], rest[i], rest[(i + 1) % m]);
            if cross2(xy[a], xy[b], xy[c]) <= EAR_EPS {
                continue;
            }
            let blocked = rest.iter().any(|&w| {
                w != a
                    && w != b
                    && w != c
                    && cross2(xy[a], xy[b], xy[w]) >= -EAR_EPS
                    && cross2(xy[b], xy[c], xy[w]) >= -EAR_EPS
                    && cross2(xy[c], xy[a], xy[w]) >= -EAR_EPS
            });
            if blocked {
                continue;
            }
            let mut score = quality(a, b, c);
            if m == 4 {
                let d = rest[(i + 2) % m];
                if cross2(xy[a], xy[c], xy[d]) <= EAR_EPS {
                    continue;
                }
                score = score.min(quality(a, c, d));
            }
            if best.is_none_or(|(s, _)| score > s) {
                best = Some((score, i));
            }
        }
        let Some((_, i)) = best else {
            return Err(internal(format!(
                "polygon {polygon}: no ear in a cell of {m} corners"
            )));
        };
        let (a, b, c) = (rest[(i + m - 1) % m], rest[i], rest[(i + 1) % m]);
        out.push([ids[a], ids[b], ids[c]]);
        rest.remove(i);
    }
    if cross2(xy[rest[0]], xy[rest[1]], xy[rest[2]]) <= EAR_EPS {
        return Err(internal(format!(
            "polygon {polygon}: a cell is left with three collinear corners"
        )));
    }
    out.push([ids[rest[0]], ids[rest[1]], ids[rest[2]]]);
    Ok(out)
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

        // The points of the paths on diagonals, per diagonal, by their parameter from the
        // lower corner.
        let mut diag_points: HashMap<(u32, u32), Vec<(f64, u32)>> = HashMap::new();
        for &k in paths {
            for (j, point) in self.cut.paths[k].interior.iter().enumerate() {
                let [c0, c1] = point.corners;
                if !tri_of_edge.contains_key(&(c0, c1)) || is_polygon_edge(c0, c1) {
                    return Err(invalid(
                        Some(p),
                        format!(
                            "path {k}, interior point {j}: ({c0}, {c1}) is not a diagonal of the \
                             triangulation of the polygon"
                        ),
                    ));
                }
                let id = self.path_vertices[k][j + 1];
                let list = diag_points.entry((c0, c1)).or_default();
                if !list.iter().any(|&(_, v)| v == id) {
                    list.push((point.t, id));
                }
            }
        }
        for (key, list) in &mut diag_points {
            list.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            if list.windows(2).any(|w| w[0].0 == w[1].0) {
                return Err(invalid(
                    Some(p),
                    format!("two points of the cut at the same place of the diagonal {key:?}"),
                ));
            }
        }
        // The points on a side of a triangle from corner c to corner d with their parameter
        // along it.
        let side_points = |c: u32, d: u32| -> Vec<(u32, f64)> {
            if is_polygon_edge(c, d) {
                let (a, b) = (verts[(c - lo) as usize], verts[(d - lo) as usize]);
                extras(c, d)
                    .into_iter()
                    .map(|id| {
                        let t = self.vertex_t[&id];
                        (id, if a < b { t } else { 1.0 - t })
                    })
                    .collect()
            } else {
                let mut v: Vec<(u32, f64)> = diag_points
                    .get(&(c.min(d), c.max(d)))
                    .map(|l| {
                        l.iter()
                            .map(|&(t, id)| (id, if c < d { t } else { 1.0 - t }))
                            .collect()
                    })
                    .unwrap_or_default();
                v.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
                v
            }
        };

        // Pieces of the cut: between consecutive points of the paths, inside one triangle.
        let mut inner_of_triangle: Vec<Vec<(u32, [f64; 3])>> = vec![Vec::new(); triangles.len()];
        let mut pieces: Vec<PieceRec> = Vec::new();
        for &k in paths {
            let path = &self.cut.paths[k];
            let ids = &self.path_vertices[k];
            let mut locs: Vec<Loc> = Vec::with_capacity(ids.len());
            for (j, &vertex) in ids.iter().enumerate() {
                let loc = if j == 0 || j + 1 == ids.len() {
                    let end = if j == 0 { path.start } else { path.end };
                    match end {
                        PathEnd::EdgePoint(i) => Loc::Side(
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
                            let junction = &self.cut.junctions[i as usize];
                            let ti = triangles
                                .iter()
                                .position(|t| *t == junction.corners)
                                .ok_or_else(|| {
                                    invalid(
                                        Some(p),
                                        format!(
                                            "path {k}: the triangle of junction {i} is not a triangle of the polygon"
                                        ),
                                    )
                                })?;
                            if !inner_of_triangle[ti].iter().any(|&(v, _)| v == vertex) {
                                inner_of_triangle[ti].push((vertex, junction.barycentric));
                            }
                            Loc::Inner(ti)
                        }
                    }
                } else {
                    let [c0, c1] = path.interior[j - 1].corners;
                    Loc::Side((c0, c1))
                };
                locs.push(loc);
            }
            for j in 0..ids.len() - 1 {
                let (x, y) = (ids[j], ids[j + 1]);
                let ti = match (locs[j], locs[j + 1]) {
                    (Loc::Inner(a), Loc::Inner(b)) => {
                        if a != b {
                            return Err(invalid(
                                Some(p),
                                format!("path {k}: a piece joins points inside two triangles"),
                            ));
                        }
                        a
                    }
                    (Loc::Inner(t), Loc::Side(e)) | (Loc::Side(e), Loc::Inner(t)) => {
                        if !tri_of_edge[&e].contains(&t) {
                            return Err(invalid(
                                Some(p),
                                format!(
                                    "path {k}: a piece joins a point inside triangle {:?} to the edge {e:?}, which is not a side of it",
                                    triangles[t]
                                ),
                            ));
                        }
                        t
                    }
                    (Loc::Side(e1), Loc::Side(e2)) => {
                        if e1 == e2 {
                            if !is_polygon_edge(e1.0, e1.1) {
                                return Err(CuttingError::Unsupported {
                                    what: format!(
                                        "path {k} of polygon {p} runs along the diagonal {e1:?}"
                                    ),
                                });
                            }
                            // Along a polygon edge: the boundary of the polygon is already
                            // cut there; only the seam is new.
                            continue;
                        }
                        let common: Vec<usize> = tri_of_edge[&e1]
                            .iter()
                            .copied()
                            .filter(|t| tri_of_edge[&e2].contains(t))
                            .collect();
                        if common.len() != 1 {
                            return Err(invalid(
                                Some(p),
                                format!(
                                    "path {k}: consecutive points on the edges {e1:?} and {e2:?} do not share a triangle of the triangulation"
                                ),
                            ));
                        }
                        common[0]
                    }
                };
                pieces.push(PieceRec {
                    a: x,
                    b: y,
                    tri: ti,
                    left: path.regions[0],
                    right: path.regions[1],
                });
            }
        }

        // The subdivision of every triangle by its pieces: cells.
        let mut cells: Vec<Cell> = Vec::new();
        for (ti, t) in triangles.iter().enumerate() {
            let corner_xy = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
            // The boundary of the triangle, counter-clockwise, with every point on it.
            let mut ring: Vec<(u32, [f64; 2], bool, bool)> = Vec::new(); // id, xy, diag side after
            for k in 0..3 {
                let (c, d) = (t[k], t[(k + 1) % 3]);
                let diag = !is_polygon_edge(c, d);
                ring.push((verts[(c - lo) as usize], corner_xy[k], diag, true));
                for (id, s) in side_points(c, d) {
                    let xy = [
                        (1.0 - s) * corner_xy[k][0] + s * corner_xy[(k + 1) % 3][0],
                        (1.0 - s) * corner_xy[k][1] + s * corner_xy[(k + 1) % 3][1],
                    ];
                    ring.push((id, xy, diag, false));
                }
            }
            let mut index: HashMap<u32, usize> = HashMap::new();
            let mut vid: Vec<u32> = Vec::new();
            let mut vxy: Vec<[f64; 2]> = Vec::new();
            for &(id, xy, _, _) in &ring {
                index.insert(id, vid.len());
                vid.push(id);
                vxy.push(xy);
            }
            for &(id, bary) in &inner_of_triangle[ti] {
                index.insert(id, vid.len());
                vid.push(id);
                vxy.push([bary[1], bary[2]]);
            }
            // Edges: the boundary (both directions, the inner one first), then the pieces.
            let mut from: Vec<usize> = Vec::new();
            let mut to: Vec<usize> = Vec::new();
            let mut kind: Vec<SideKind> = Vec::new();
            for i in 0..ring.len() {
                let (a, b) = (i, (i + 1) % ring.len());
                let diag = ring[i].2;
                from.push(a);
                to.push(b);
                kind.push(SideKind::Boundary { diag });
                from.push(b);
                to.push(a);
                kind.push(SideKind::Boundary { diag });
            }
            let mut mine: Vec<usize> = Vec::new();
            for (r, piece) in pieces.iter().enumerate() {
                if piece.tri == ti {
                    mine.push(r);
                    let (a, b) = (index[&piece.a], index[&piece.b]);
                    from.push(a);
                    to.push(b);
                    kind.push(SideKind::Piece {
                        rec: r,
                        forward: true,
                    });
                    from.push(b);
                    to.push(a);
                    kind.push(SideKind::Piece {
                        rec: r,
                        forward: false,
                    });
                }
            }
            // Pieces must not cross.
            for (i, &r1) in mine.iter().enumerate() {
                for &r2 in &mine[i + 1..] {
                    let (p1, p2) = (index[&pieces[r1].a], index[&pieces[r1].b]);
                    let (q1, q2) = (index[&pieces[r2].a], index[&pieces[r2].b]);
                    if [p1, p2].contains(&q1) || [p1, p2].contains(&q2) {
                        continue;
                    }
                    let d1 = cross2(vxy[p1], vxy[p2], vxy[q1]);
                    let d2 = cross2(vxy[p1], vxy[p2], vxy[q2]);
                    let d3 = cross2(vxy[q1], vxy[q2], vxy[p1]);
                    let d4 = cross2(vxy[q1], vxy[q2], vxy[p2]);
                    if d1 * d2 < -EAR_EPS * EAR_EPS && d3 * d4 < -EAR_EPS * EAR_EPS {
                        return Err(invalid(
                            Some(p),
                            format!("two pieces of the cut cross inside triangle {t:?}"),
                        ));
                    }
                }
            }
            // Outgoing half-edges of every vertex, by angle.
            let mut outgoing: Vec<Vec<usize>> = vec![Vec::new(); vid.len()];
            for h in 0..from.len() {
                outgoing[from[h]].push(h);
            }
            for list in &mut outgoing {
                list.sort_by(|&a, &b| {
                    let angle = |h: usize| {
                        (vxy[to[h]][1] - vxy[from[h]][1]).atan2(vxy[to[h]][0] - vxy[from[h]][0])
                    };
                    angle(a)
                        .partial_cmp(&angle(b))
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            }
            for (v, list) in outgoing.iter().enumerate() {
                if list.len() == 1 {
                    return Err(invalid(
                        Some(p),
                        format!(
                            "the cut ends inside triangle {t:?} at vertex {} (a dead end)",
                            vid[v]
                        ),
                    ));
                }
            }
            // Everything must be connected to the boundary.
            {
                let mut seen = vec![false; vid.len()];
                let mut stack = vec![0usize];
                seen[0] = true;
                while let Some(v) = stack.pop() {
                    for &h in &outgoing[v] {
                        if !seen[to[h]] {
                            seen[to[h]] = true;
                            stack.push(to[h]);
                        }
                    }
                }
                if let Some(v) = seen.iter().position(|s| !s) {
                    return Err(CuttingError::Unsupported {
                        what: format!(
                            "a closed piece of the cut lies inside the triangle {t:?} of polygon {p} (it would make a hole in a polygon); vertex {}",
                            vid[v]
                        ),
                    });
                }
            }
            // Faces: follow half-edges, turning at each vertex to the next one clockwise from
            // the way back.
            let mut visited = vec![false; from.len()];
            for h0 in 0..from.len() {
                if visited[h0] {
                    continue;
                }
                let mut cycle = Vec::new();
                let mut h = h0;
                while !visited[h] {
                    visited[h] = true;
                    cycle.push(h);
                    let twin = h ^ 1;
                    let list = &outgoing[to[h]];
                    let pos = list.iter().position(|&x| x == twin).expect("twin");
                    h = list[(pos + list.len() - 1) % list.len()];
                }
                if h != h0 {
                    return Err(internal(format!("polygon {p}: a face walk did not close")));
                }
                let area2: f64 = cycle
                    .iter()
                    .map(|&h| {
                        let (a, b) = (vxy[from[h]], vxy[to[h]]);
                        a[0] * b[1] - a[1] * b[0]
                    })
                    .sum();
                if area2 <= 0.0 {
                    continue; // the outer face
                }
                let mut region = None;
                for &h in &cycle {
                    if let SideKind::Piece { rec, forward } = kind[h] {
                        let r = if forward {
                            pieces[rec].left
                        } else {
                            pieces[rec].right
                        };
                        if *region.get_or_insert(r) != r {
                            return Err(invalid(
                                Some(p),
                                format!(
                                    "a part of triangle {t:?} lies on the left of a piece of region {r} and on the left of another of region {}",
                                    region.unwrap()
                                ),
                            ));
                        }
                    }
                }
                cells.push(Cell {
                    ids: cycle.iter().map(|&h| vid[from[h]]).collect(),
                    xy: cycle.iter().map(|&h| vxy[from[h]]).collect(),
                    sides: cycle.iter().map(|&h| kind[h]).collect(),
                    region,
                });
            }
        }

        // Parts: cells joined across the diagonals.
        let mut dsu = Dsu((0..cells.len() as u32).collect());
        let mut on_diagonal: HashMap<(u32, u32), usize> = HashMap::new();
        for (ci, cell) in cells.iter().enumerate() {
            for i in 0..cell.ids.len() {
                if cell.sides[i] == (SideKind::Boundary { diag: true }) {
                    on_diagonal.insert((cell.ids[i], cell.ids[(i + 1) % cell.ids.len()]), ci);
                }
            }
        }
        for (&(u, v), &ci) in &on_diagonal {
            if let Some(&cj) = on_diagonal.get(&(v, u)) {
                dsu.union(ci as u32, cj as u32);
            }
        }
        let mut node_triangles: Vec<Vec<[u32; 3]>> = Vec::with_capacity(cells.len());
        for cell in &cells {
            node_triangles.push(triangulate_cell(&cell.ids, &cell.xy, self.positions, p)?);
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

        let mut groups: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
        for ci in 0..cells.len() {
            groups.entry(dsu.find(ci as u32)).or_default().push(ci);
        }
        let mut parts = Vec::with_capacity(groups.len());
        for members in groups.values() {
            let mut directed: HashSet<(u32, u32)> = HashSet::new();
            let mut region: Option<u32> = None;
            for &ci in members {
                let cell = &cells[ci];
                if let Some(r) = cell.region {
                    if *region.get_or_insert(r) != r {
                        return Err(invalid(
                            Some(p),
                            format!(
                                "a part of the polygon lies in the regions {r} and {}",
                                region.unwrap()
                            ),
                        ));
                    }
                }
                for i in 0..cell.ids.len() {
                    let (u, v) = (cell.ids[i], cell.ids[(i + 1) % cell.ids.len()]);
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
            let triangles: Vec<[u32; 3]> = members
                .iter()
                .flat_map(|&ci| node_triangles[ci].iter().copied())
                .collect();
            let part_key = vertices
                .iter()
                .map(|v| key.get(v).copied().unwrap_or(f64::INFINITY))
                .fold(f64::INFINITY, f64::min);
            parts.push(Part {
                vertices,
                triangles,
                key: part_key,
                region,
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
