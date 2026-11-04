//! A regular 2D grid of field values and marching squares.

use super::field::Pat2dField;
use crate::PatError;
use std::collections::BTreeMap;

/// Values of a function on the nodes of a regular grid, `values[j * nx + i]` at
/// `origin + (i hx, j hy)`.
#[derive(Debug, Clone, PartialEq)]
pub struct Grid2d {
    /// Position of node `(0, 0)`.
    pub origin: [f64; 2],
    /// Spacing along x and y.
    pub spacing: [f64; 2],
    /// Nodes along x (≥ 2).
    pub nx: usize,
    /// Nodes along y (≥ 2).
    pub ny: usize,
    /// Row-major values.
    pub values: Vec<f64>,
}

impl Grid2d {
    /// Position of node `(i, j)`.
    pub fn node(&self, i: usize, j: usize) -> [f64; 2] {
        [self.origin[0] + self.spacing[0] * i as f64, self.origin[1] + self.spacing[1] * j as f64]
    }

    /// All node positions in storage order.
    pub fn nodes(&self) -> Vec<[f64; 2]> {
        (0..self.ny).flat_map(|j| (0..self.nx).map(move |i| (i, j))).map(|(i, j)| self.node(i, j)).collect()
    }
}

/// Evaluate the field at the `nx × ny` nodes spanning `bbox` (both ends included).
pub fn evaluate_grid(field: &Pat2dField, bbox: [[f64; 2]; 2], nx: usize, ny: usize) -> Result<Grid2d, PatError> {
    if nx < 2 || ny < 2 || bbox.iter().flatten().any(|c| !c.is_finite()) || bbox[1][0] <= bbox[0][0] || bbox[1][1] <= bbox[0][1] {
        return Err(PatError::InvalidOption { what: format!("evaluate_grid: {nx} x {ny} nodes on the box {bbox:?}") });
    }
    if nx.saturating_mul(ny) > 1 << 28 {
        return Err(PatError::InvalidOption { what: format!("evaluate_grid: {nx} x {ny} nodes is too many") });
    }
    let mut g = Grid2d {
        origin: bbox[0],
        spacing: [(bbox[1][0] - bbox[0][0]) / (nx - 1) as f64, (bbox[1][1] - bbox[0][1]) / (ny - 1) as f64],
        nx,
        ny,
        values: vec![0.0; nx * ny],
    };
    let nodes = g.nodes();
    field.evaluate(&nodes, Some(&mut g.values), None)?;
    Ok(g)
}

/// A segment of a contour, oriented so that values below the level are on its left.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment2 {
    /// Start.
    pub a: [f64; 2],
    /// End.
    pub b: [f64; 2],
    /// Vertex ids (grid edges) of the endpoints, for chaining.
    pub ids: [u64; 2],
}

/// A chain of contour vertices; `closed` repeats no point.
#[derive(Debug, Clone, PartialEq)]
pub struct Polyline {
    /// Vertices.
    pub points: Vec<[f64; 2]>,
    /// Whether the last vertex joins the first.
    pub closed: bool,
}

/// Marching squares at `level`. A cell with four crossings is resolved by the bilinear value at
/// its center (the asymptotic decider of Nielson and Hamann for a bilinear interpolant). Vertices
/// are linear interpolations along grid edges, computed from the same two values from both cells,
/// so shared vertices agree bitwise. Cells with a NaN corner are skipped.
pub fn marching_squares(grid: &Grid2d, level: f64) -> Vec<Segment2> {
    let (nx, ny) = (grid.nx, grid.ny);
    let v = |i: usize, j: usize| grid.values[j * nx + i];
    // vertex on the grid edge from (i0,j0) to (i1,j1), with (i0,j0) <= (i1,j1) in storage order
    let vert = |i0: usize, j0: usize, horizontal: bool| -> ([f64; 2], u64) {
        let (i1, j1) = if horizontal { (i0 + 1, j0) } else { (i0, j0 + 1) };
        let (va, vb) = (v(i0, j0), v(i1, j1));
        let t = (level - va) / (vb - va);
        let (pa, pb) = (grid.node(i0, j0), grid.node(i1, j1));
        ([pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1])], 2 * (j0 * nx + i0) as u64 + u64::from(!horizontal))
    };
    let mut out = Vec::new();
    for j in 0..ny - 1 {
        for i in 0..nx - 1 {
            let c = [v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)];
            if c.iter().any(|x| x.is_nan()) {
                continue;
            }
            let ins = [c[0] < level, c[1] < level, c[2] < level, c[3] < level];
            let count = ins.iter().filter(|&&b| b).count();
            if count == 0 || count == 4 {
                continue;
            }
            // edges: e0 = c0-c1 (horizontal at (i,j)), e1 = c1-c2 (vertical at (i+1,j)),
            // e2 = c3-c2 (horizontal at (i,j+1)), e3 = c0-c3 (vertical at (i,j))
            let edge = |e: usize| match e {
                0 => vert(i, j, true),
                1 => vert(i + 1, j, false),
                2 => vert(i, j + 1, true),
                _ => vert(i, j, false),
            };
            // corner k touches edges (e_{k}, e_{k+3 mod 4}) in this labelling: c0:e0,e3 c1:e0,e1 c2:e1,e2 c3:e2,e3
            let corner_edges = [(0, 3), (0, 1), (1, 2), (2, 3)];
            let mut pairs: Vec<(usize, usize)> = Vec::new();
            let alternating = count == 2 && ins[0] == ins[2];
            if alternating {
                let center = 0.25 * (c[0] + c[1] + c[2] + c[3]);
                let center_inside = center < level;
                // inside connected through the center: cut off the outside corners, else the inside ones
                for k in 0..4 {
                    if ins[k] != center_inside {
                        pairs.push(corner_edges[k]);
                    }
                }
            } else {
                // one cut-off corner (count 1 or 3) or a straight cut (count 2, adjacent)
                if count == 1 || count == 3 {
                    let k = (0..4).find(|&k| ins[k] == (count == 1)).unwrap_or(0);
                    pairs.push(corner_edges[k]);
                } else {
                    let crossing: Vec<usize> = (0..4).filter(|&e| ins[corner_of(e).0] != ins[corner_of(e).1]).collect();
                    if crossing.len() == 2 {
                        pairs.push((crossing[0], crossing[1]));
                    }
                }
            }
            let corner_pos = [grid.node(i, j), grid.node(i + 1, j), grid.node(i + 1, j + 1), grid.node(i, j + 1)];
            for (ea, eb) in pairs {
                let ((pa, ia), (pb, ib)) = (edge(ea), edge(eb));
                let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
                // nearest inside corner decides the orientation
                let mut best: Option<(f64, [f64; 2])> = None;
                for k in 0..4 {
                    if ins[k] {
                        let d = (corner_pos[k][0] - mid[0]).powi(2) + (corner_pos[k][1] - mid[1]).powi(2);
                        if best.is_none_or(|b| d < b.0) {
                            best = Some((d, corner_pos[k]));
                        }
                    }
                }
                let (a, b, ids) = match best {
                    Some((_, q)) => {
                        let cross = (pb[0] - pa[0]) * (q[1] - pa[1]) - (pb[1] - pa[1]) * (q[0] - pa[0]);
                        if cross > 0.0 { (pa, pb, [ia, ib]) } else { (pb, pa, [ib, ia]) }
                    }
                    None => (pa, pb, [ia, ib]),
                };
                out.push(Segment2 { a, b, ids });
            }
        }
    }
    out
}

/// The two corners of edge `e` of a cell (see [`marching_squares`]).
fn corner_of(e: usize) -> (usize, usize) {
    match e {
        0 => (0, 1),
        1 => (1, 2),
        2 => (3, 2),
        _ => (0, 3),
    }
}

/// Chain oriented segments into polylines (closed loops, or open chains ending at the grid
/// border). Deterministic: chains start at the smallest vertex id.
pub fn chain_segments(segments: &[Segment2]) -> Vec<Polyline> {
    let mut next: BTreeMap<u64, usize> = BTreeMap::new();
    let mut has_in: BTreeMap<u64, usize> = BTreeMap::new();
    for (s, seg) in segments.iter().enumerate() {
        next.insert(seg.ids[0], s);
        has_in.insert(seg.ids[1], s);
    }
    let mut used = vec![false; segments.len()];
    let mut out = Vec::new();
    let walk = |start: usize, used: &mut Vec<bool>| -> Polyline {
        let mut pts = vec![segments[start].a];
        let mut s = start;
        let mut closed = false;
        loop {
            used[s] = true;
            let seg = &segments[s];
            match next.get(&seg.ids[1]) {
                Some(&t) if !used[t] => {
                    pts.push(seg.b);
                    s = t;
                }
                Some(_) => {
                    closed = true;
                    break;
                }
                None => {
                    pts.push(seg.b);
                    break;
                }
            }
        }
        Polyline { points: pts, closed }
    };
    // open chains first: start at vertices without an incoming segment
    let mut starts: Vec<(u64, usize)> = next.iter().filter(|(id, _)| !has_in.contains_key(id)).map(|(&id, &s)| (id, s)).collect();
    starts.sort();
    for (_, s) in starts {
        if !used[s] {
            out.push(walk(s, &mut used));
        }
    }
    for s in 0..segments.len() {
        if !used[s] {
            out.push(walk(s, &mut used));
        }
    }
    out
}
