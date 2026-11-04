//! Marching cubes with the asymptotic decider and a narrow band (SPEC.md §8.6, D18). Owned by WP6.
//!
//! Our own implementation, on a regular grid over a box:
//!
//! 1. **Narrow band.** The grid is cut into blocks of `B³` cells (default 8, "masking voxels not
//!    within 8 voxels of the contour", paper Sec. 5.3). `φ − iso` is evaluated at the block
//!    corners; a block is skipped when some corner has `|φ − iso| > L · diag`, where `diag` is the
//!    block's diagonal and `L` the Lipschitz bound (default **2.0**: the blend is not exactly
//!    1-Lipschitz, `‖∇φ‖` is 1 up to a modest excess; the report states the bound used and
//!    measures the largest slope `|Δφ|/h` seen along evaluated grid edges). All grid values of
//!    the remaining blocks are then evaluated in large batches. **Growth:** a kept block whose
//!    face towards a skipped neighbour has a sign change keeps that neighbour too (repeated), so
//!    a bound that is too small can only lose isolated pieces of the surface that lie wholly in
//!    skipped blocks, never leave cracks. (The fitted blend has `‖∇φ‖` well above 1 in places,
//!    measured slopes of 5–50 per grid edge near jumps of the evaluation neighborhood on the
//!    reference clouds, so the bound is a heuristic; see `NOTES-WP6.md`.)
//! 2. **Cells.** A cell is *inside* at a corner where `φ − iso < 0` (zero counts as outside).
//!    Instead of a 256-entry table the crossing polygons are built from the faces: every face
//!    with two crossings gives one segment, a face with four crossings (the ambiguous case) is
//!    resolved by the **asymptotic decider** of Nielson and Hamann (1991): the inside diagonal is
//!    connected iff `f_P1 f_P2 > f_Q1 f_Q2` (the saddle value of the bilinear interpolant has the
//!    sign of the inside corners). The decision uses only the four face values in a form that is
//!    invariant under rotation and reflection of the face, so the two cells sharing a face agree
//!    bitwise, and the surface is crack-free and manifold for generic data. The segments close to
//!    loops on the cube boundary; loops are triangulated by the best fan whose diagonals stay inside the cell (a centroid star otherwise).
//! 3. **Vertices** live on grid edges and are shared (keyed by the edge). Linear interpolation
//!    `t = f₀/(f₀ − f₁)` `[formula: pat-mc-interpolation]`; optionally a few safeguarded Newton
//!    steps along the edge ([`McOptions::refine_steps`]).
//!
//! Everything is bitwise independent of the number of threads: every grid value is a pure
//! function of its position, the vertices are numbered by sorted edge key, and triangles follow
//! the block/cell order.

use crate::csg::ImplicitField;
use crate::pool::make_pool;
use crate::PatError;
use ddg_core::TriMesh;
use rayon::prelude::*;
use std::time::Instant;

/// Grid resolution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Resolution {
    /// `n` cells along the longest side of the box; cubic cells, the grid starts at the lower
    /// corner and may extend beyond the upper ones by less than a cell.
    Longest(usize),
    /// Cells along each axis; the grid fills the box exactly.
    Cells([usize; 3]),
}

/// Which parts of the grid are evaluated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Band {
    /// Every grid point (reference for tests).
    Dense,
    /// Blocks of `block³` cells whose corner values exceed `bound ×` the block diagonal are skipped.
    Lipschitz {
        /// The Lipschitz bound `L` (finite, `> 0`).
        bound: f64,
        /// Cells per block side (`≥ 1`).
        block: usize,
    },
}

/// Options of [`marching_cubes`].
#[derive(Debug, Clone, PartialEq)]
pub struct McOptions {
    /// The box (`None`: the field's [`ImplicitField::bbox`]).
    pub bbox: Option<[[f64; 3]; 2]>,
    /// Grid resolution.
    pub resolution: Resolution,
    /// The level extracted: `φ = isovalue` (an offset when `≠ 0`).
    pub isovalue: f64,
    /// Narrow band.
    pub band: Band,
    /// Safeguarded Newton steps along each edge after the linear interpolation (0: none).
    pub refine_steps: usize,
    /// Also return unit normals `∇φ/‖∇φ‖` at the vertices (zero vector where `∇φ = 0`).
    pub normals: bool,
    /// Run `TriMesh::new` on the output and fill [`McReport::mesh`].
    pub check_mesh: bool,
    /// Threads (`None`: all cores).
    pub threads: Option<usize>,
}

impl Default for McOptions {
    fn default() -> Self {
        McOptions {
            bbox: None,
            resolution: Resolution::Longest(64),
            isovalue: 0.0,
            band: Band::Lipschitz { bound: 2.0, block: 8 },
            refine_steps: 0,
            normals: false,
            check_mesh: true,
            threads: None,
        }
    }
}

/// The extracted surface, oriented so that normals point towards increasing `φ` (outward).
#[derive(Debug, Clone, PartialEq)]
pub struct McMesh {
    /// Vertex positions, ordered by grid edge (z-major).
    pub positions: Vec<[f64; 3]>,
    /// Triangles.
    pub faces: Vec<[u32; 3]>,
    /// Unit normals from `∇φ` if requested.
    pub normals: Option<Vec<[f64; 3]>>,
}

impl McMesh {
    /// The mesh as a validated [`TriMesh`] (fails if it is not an oriented manifold).
    pub fn to_trimesh(&self) -> Result<TriMesh, ddg_core::topology::MeshError> {
        TriMesh::new(self.positions.len(), self.faces.clone())
    }

    /// Enclosed signed volume (divergence theorem; meaningful for closed meshes).
    pub fn signed_volume(&self) -> f64 {
        let mut v = 0.0;
        for f in &self.faces {
            let (a, b, c) = (self.positions[f[0] as usize], self.positions[f[1] as usize], self.positions[f[2] as usize]);
            v += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
        }
        v / 6.0
    }

    /// Triangle areas.
    pub fn areas(&self) -> Vec<f64> {
        self.faces
            .iter()
            .map(|f| {
                let (a, b, c) = (self.positions[f[0] as usize], self.positions[f[1] as usize], self.positions[f[2] as usize]);
                let (u, w) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
                let n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
                0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
            })
            .collect()
    }
}

/// Topology of the output as seen by `ddg_core::TriMesh::new`.
#[derive(Debug, Clone, PartialEq)]
pub enum MeshCheck {
    /// Not checked.
    Skipped,
    /// An oriented manifold mesh (edge- and vertex-manifold, consistent orientation).
    Manifold {
        /// Edges used by one triangle (nonzero only where the surface leaves the box).
        boundary_edges: usize,
        /// `V − E + F`.
        euler_characteristic: i64,
        /// Connected components.
        components: usize,
    },
    /// Not a valid mesh; the message names the edge or vertex.
    Invalid(String),
}

/// What a run did.
#[derive(Debug, Clone, PartialEq)]
pub struct McReport {
    /// Cells per axis.
    pub cells: [usize; 3],
    /// Spacing per axis.
    pub spacing: [f64; 3],
    /// The grid origin and the box requested.
    pub origin: [f64; 3],
    /// Blocks in total.
    pub blocks: usize,
    /// Blocks added by the growth rule (surface crossing into a skipped neighbour) and the rounds it took.
    pub blocks_grown: usize,
    /// Evaluation rounds of the narrow band (1: no growth).
    pub band_rounds: usize,
    /// Blocks kept by the narrow band.
    pub blocks_active: usize,
    /// Cells in the kept blocks (the cells visited).
    pub cells_visited: usize,
    /// Cells with a sign change.
    pub cells_crossing: usize,
    /// Faces (counted per cell) with four crossings that the decider resolved.
    pub ambiguous_cell_faces: usize,
    /// Grid values evaluated (coarse corners + kept blocks).
    pub grid_values: usize,
    /// Grid points of the whole grid.
    pub grid_points: usize,
    /// Extra field evaluations by edge refinement and normals.
    pub extra_evaluations: usize,
    /// Vertices.
    pub vertices: usize,
    /// Triangles.
    pub triangles: usize,
    /// The Lipschitz bound used (`None`: dense).
    pub lipschitz_bound: Option<f64>,
    /// Largest `|Δφ|/spacing` along grid edges inside the kept blocks.
    pub max_edge_slope: f64,
    /// `max_edge_slope > lipschitz_bound`: the bound was too small for this field; blocks may have been skipped wrongly.
    pub bound_violated: bool,
    /// Normals that were zero because `∇φ = 0`.
    pub zero_normals: usize,
    /// Topology check.
    pub mesh: MeshCheck,
    /// Seconds: coarse pass and band, grid evaluation, cells, vertices and refinement, total.
    pub seconds: [f64; 5],
    /// What the fitted fields behind the surface evaluated during this run (the difference of
    /// [`ImplicitField::eval_totals`] before and after; batch evaluations only). `None` for
    /// analytic fields.
    pub eval: Option<crate::EvalReport>,
}

// ---------------------------------------------------------------------------------------------
// Cube tables

const fn corner_offset(c: usize) -> [usize; 3] {
    [c & 1, (c >> 1) & 1, (c >> 2) & 1]
}

/// The 12 edges as `(lower corner, axis)`, ordered by axis, then by corner.
const fn edge_list() -> [(usize, usize); 12] {
    let mut out = [(0, 0); 12];
    let mut n = 0;
    let mut a = 0;
    while a < 3 {
        let mut c = 0;
        while c < 8 {
            if c & (1 << a) == 0 {
                out[n] = (c, a);
                n += 1;
            }
            c += 1;
        }
        a += 1;
    }
    out
}
const EDGES: [(usize, usize); 12] = edge_list();

/// `EDGE_ID[lower corner][axis]` (255 where the corner has the axis bit set).
const fn edge_ids() -> [[u8; 3]; 8] {
    let mut out = [[255u8; 3]; 8];
    let e = edge_list();
    let mut i = 0;
    while i < 12 {
        out[e[i].0][e[i].1] = i as u8;
        i += 1;
    }
    out
}
const EDGE_ID: [[u8; 3]; 8] = edge_ids();

/// The six faces, corners counter-clockwise seen from outside the cube.
const fn face_list() -> [[usize; 4]; 6] {
    let mut out = [[0; 4]; 6];
    let mut n = 0;
    let mut a = 0;
    while a < 3 {
        let mut s = 0;
        while s < 2 {
            let u = (a + 1) % 3;
            let v = (a + 2) % 3;
            let q00 = s << a;
            let q10 = (s << a) | (1 << u);
            let q11 = (s << a) | (1 << u) | (1 << v);
            let q01 = (s << a) | (1 << v);
            out[n] = if s == 1 { [q00, q10, q11, q01] } else { [q00, q01, q11, q10] };
            n += 1;
            s += 1;
        }
        a += 1;
    }
    out
}
const FACES: [[usize; 4]; 6] = face_list();

const NONE: u8 = 255;

fn edge_between(c0: usize, c1: usize) -> usize {
    let axis = (c0 ^ c1).trailing_zeros() as usize;
    let lower = c0 & c1;
    EDGE_ID[lower][axis] as usize
}

/// Link the crossings of one cell: `next[e]` is the crossing after edge `e` on its loop, every
/// segment going from the edge where a face boundary (counter-clockwise) enters the inside to
/// the edge where it leaves it. Returns the number of ambiguous faces.
fn build_segments(v: &[f64; 8], inside: &[bool; 8], next: &mut [u8; 12]) -> usize {
    let mut ambiguous = 0;
    for q in &FACES {
        let ins = [inside[q[0]], inside[q[1]], inside[q[2]], inside[q[3]]];
        let ee = [edge_between(q[0], q[1]), edge_between(q[1], q[2]), edge_between(q[2], q[3]), edge_between(q[3], q[0])];
        let cross = [ins[0] != ins[1], ins[1] != ins[2], ins[2] != ins[3], ins[3] != ins[0]];
        let count = cross.iter().filter(|c| **c).count();
        if count == 2 {
            let (mut enter, mut exit) = (NONE, NONE);
            for m in 0..4 {
                if cross[m] {
                    if ins[(m + 1) % 4] {
                        enter = ee[m] as u8;
                    } else {
                        exit = ee[m] as u8;
                    }
                }
            }
            next[enter as usize] = exit;
        } else if count == 4 {
            ambiguous += 1;
            // The inside pair P is {q0, q2} or {q1, q3}. The decider (Nielson–Hamann): the inside
            // diagonal is connected iff f_P1 f_P2 > f_Q1 f_Q2.
            let (p, o) = if ins[0] { ([0, 2], [1, 3]) } else { ([1, 3], [0, 2]) };
            let connected = v[q[p[0]]] * v[q[p[1]]] > v[q[o[0]]] * v[q[o[1]]];
            let e = ee;
            let set = |next: &mut [u8; 12], from: usize, to: usize| next[e[from]] = e[to] as u8;
            if ins[0] {
                // e0, e2 exit; e1, e3 enter
                if connected {
                    set(next, 1, 0);
                    set(next, 3, 2);
                } else {
                    set(next, 3, 0);
                    set(next, 1, 2);
                }
            } else {
                // e0, e2 enter; e1, e3 exit
                if connected {
                    set(next, 0, 3);
                    set(next, 2, 1);
                } else {
                    set(next, 0, 1);
                    set(next, 2, 3);
                }
            }
        }
    }
    ambiguous
}

// ---------------------------------------------------------------------------------------------
// Grid

struct Grid {
    n: [usize; 3],
    lo: [f64; 3],
    h: [f64; 3],
    /// Strides of grid points per axis.
    stride: [usize; 3],
}

impl Grid {
    fn npoints(&self) -> usize {
        (self.n[0] + 1) * (self.n[1] + 1) * (self.n[2] + 1)
    }
    fn gp(&self, i: usize, j: usize, k: usize) -> usize {
        i + j * self.stride[1] + k * self.stride[2]
    }
    fn decode(&self, gp: usize) -> [usize; 3] {
        [gp % self.stride[1], (gp / self.stride[1]) % (self.n[1] + 1), gp / self.stride[2]]
    }
    fn coord(&self, ijk: [usize; 3]) -> [f64; 3] {
        [self.lo[0] + ijk[0] as f64 * self.h[0], self.lo[1] + ijk[1] as f64 * self.h[1], self.lo[2] + ijk[2] as f64 * self.h[2]]
    }
    /// Position on the edge `key` at parameter `t`.
    fn vertex(&self, key: u64, t: f64) -> [f64; 3] {
        let (gp, axis) = ((key / 3) as usize, (key % 3) as usize);
        let ijk = self.decode(gp);
        let mut p = self.coord(ijk);
        p[axis] = self.lo[axis] + (ijk[axis] as f64 + t) * self.h[axis];
        p
    }
    /// `[formula: pat-mc-interpolation]` the root of the linear interpolant on the edge `key`.
    fn linear_t(&self, vals: &[f64], key: u64) -> f64 {
        let (gp, axis) = ((key / 3) as usize, (key % 3) as usize);
        let (fa, fb) = (vals[gp], vals[gp + self.stride[axis]]);
        fa / (fa - fb)
    }
}

fn invalid(what: String) -> PatError {
    PatError::InvalidOption { what }
}

fn make_grid(bbox: [[f64; 3]; 2], res: Resolution) -> Result<Grid, PatError> {
    if !bbox.iter().flatten().all(|c| c.is_finite()) {
        return Err(invalid(format!("marching cubes box {bbox:?} is not finite")));
    }
    let ext = [0, 1, 2].map(|a| bbox[1][a] - bbox[0][a]);
    if ext.iter().any(|e| *e <= 0.0) {
        return Err(invalid(format!("marching cubes box {bbox:?} has an empty axis")));
    }
    let (n, h) = match res {
        Resolution::Longest(n) => {
            if n == 0 {
                return Err(invalid("Resolution::Longest(0)".into()));
            }
            let h = ext.iter().cloned().fold(0.0, f64::max) / n as f64;
            let cnt = ext.map(|e| ((e / h - 1e-9).ceil() as usize).max(1));
            (cnt, [h; 3])
        }
        Resolution::Cells(c) => {
            if c.contains(&0) {
                return Err(invalid(format!("Resolution::Cells({c:?}) has an empty axis")));
            }
            (c, [0, 1, 2].map(|a| ext[a] / c[a] as f64))
        }
    };
    let total = (n[0] as u128 + 1) * (n[1] as u128 + 1) * (n[2] as u128 + 1);
    if total > (1u128 << 29) {
        return Err(invalid(format!("grid {n:?} has {total} points (limit 2^29 = 536,870,912 values, 4.3 GB); lower the resolution")));
    }
    Ok(Grid { n, lo: bbox[0], h, stride: [1, n[0] + 1, (n[0] + 1) * (n[1] + 1)] })
}

/// Evaluate the grid points `gps` (all distinct) in batches and store them in `vals`.
fn evaluate_points(field: &dyn ImplicitField, grid: &Grid, gps: &[usize], vals: &mut [f64], pool: &rayon::ThreadPool) -> Result<(), PatError> {
    const BATCH: usize = 1 << 18;
    let mut xs: Vec<[f64; 3]> = Vec::new();
    let mut out: Vec<f64> = Vec::new();
    for chunk in gps.chunks(BATCH) {
        xs.clear();
        xs.extend(chunk.iter().map(|g| grid.coord(grid.decode(*g))));
        out.clear();
        out.resize(chunk.len(), 0.0);
        pool.install(|| field.distances(&xs, &mut out))?;
        for (g, v) in chunk.iter().zip(&out) {
            if !v.is_finite() {
                return Err(PatError::NonFinite { index: *g, what: "field value" });
            }
            vals[*g] = *v;
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------------------------

#[derive(Default)]
struct BlockOut {
    tris: Vec<[u64; 3]>,
    centroids: Vec<(u64, [f64; 3])>,
    crossing: usize,
    ambiguous: usize,
    slope: f64,
}

/// Triangulate one loop of crossings: the fan whose triangles have the best worst quality among
/// those whose diagonals do not join two vertices of a common cube face (the neighbour across
/// that face could draw the same diagonal and make the edge non-manifold); a star around the
/// centroid (a vertex private to the cell) when no such fan exists.
fn emit_loop(grid: &Grid, vals: &[f64], keys: &[u64], lp: &[u8], centroid_key: u64, out: &mut BlockOut) {
    let m = lp.len();
    if m < 3 {
        return;
    }
    if m == 3 {
        out.tris.push([keys[lp[0] as usize], keys[lp[1] as usize], keys[lp[2] as usize]]);
        return;
    }
    let pos: Vec<[f64; 3]> = lp
        .iter()
        .map(|e| {
            let k = keys[*e as usize];
            grid.vertex(k, grid.linear_t(vals, k))
        })
        .collect();
    let quality = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
        let (u, w) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        let n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        let area2 = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let d = |x: [f64; 3], y: [f64; 3]| (x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2);
        let s = d(a, b) + d(b, c) + d(c, a);
        if s > 0.0 { 2.0 * 3f64.sqrt() * area2 / s } else { 0.0 }
    };
    let face_masks: [u16; 6] = std::array::from_fn(|f| {
        let q = FACES[f];
        (0..4).fold(0u16, |acc, m| acc | (1 << edge_between(q[m], q[(m + 1) % 4])))
    });
    let shares_face = |a: u8, b: u8| {
        let t = (1u16 << a) | (1u16 << b);
        face_masks.iter().any(|fm| fm & t == t)
    };
    // a diagonal (i, j) is allowed if the loop neighbours or no common face
    let diag_ok = |i: usize, j: usize| {
        let adj = (i + 1) % m == j || (j + 1) % m == i;
        adj || !shares_face(lp[i], lp[j])
    };
    let mut best: Option<(usize, f64)> = None;
    for a in 0..m {
        let mut worst = f64::INFINITY;
        let mut ok = true;
        for i in 1..m - 1 {
            let (b, c) = ((a + i) % m, (a + i + 1) % m);
            ok &= diag_ok(a, b) && diag_ok(a, c);
            worst = worst.min(quality(pos[a], pos[b], pos[c]));
        }
        if ok && best.is_none_or(|(_, q)| worst > q) {
            best = Some((a, worst));
        }
    }
    if let Some((a, _)) = best {
        for i in 1..m - 1 {
            out.tris.push([keys[lp[a] as usize], keys[lp[(a + i) % m] as usize], keys[lp[(a + i + 1) % m] as usize]]);
        }
    } else {
        let mut c = [0.0; 3];
        for p in &pos {
            for i in 0..3 {
                c[i] += p[i] / m as f64;
            }
        }
        out.centroids.push((centroid_key, c));
        for i in 0..m {
            out.tris.push([centroid_key, keys[lp[i] as usize], keys[lp[(i + 1) % m] as usize]]);
        }
    }
}

fn process_block(grid: &Grid, vals: &[f64], edge_limit: u64, lo: [usize; 3], hi: [usize; 3]) -> BlockOut {
    let mut out = BlockOut::default();
    let offs: [usize; 8] = std::array::from_fn(|c| {
        let o = corner_offset(c);
        o[0] * grid.stride[0] + o[1] * grid.stride[1] + o[2] * grid.stride[2]
    });
    // slopes along the grid edges of the block
    let mut slope: f64 = 0.0;
    for k in lo[2]..=hi[2] {
        for j in lo[1]..=hi[1] {
            for i in lo[0]..=hi[0] {
                let g = grid.gp(i, j, k);
                let v = vals[g];
                for a in 0..3 {
                    let ijk = [i, j, k];
                    if ijk[a] < hi[a] {
                        slope = slope.max((vals[g + grid.stride[a]] - v).abs() / grid.h[a]);
                    }
                }
            }
        }
    }
    out.slope = slope;
    for k in lo[2]..hi[2] {
        for j in lo[1]..hi[1] {
            for i in lo[0]..hi[0] {
                let base = grid.gp(i, j, k);
                let v: [f64; 8] = std::array::from_fn(|c| vals[base + offs[c]]);
                let inside: [bool; 8] = std::array::from_fn(|c| v[c] < 0.0);
                let n_in = inside.iter().filter(|b| **b).count();
                if n_in == 0 || n_in == 8 {
                    continue;
                }
                out.crossing += 1;
                let mut next = [NONE; 12];
                out.ambiguous += build_segments(&v, &inside, &mut next);
                let keys: [u64; 12] = std::array::from_fn(|e| {
                    let (c0, a) = EDGES[e];
                    (base + offs[c0]) as u64 * 3 + a as u64
                });
                let mut seen = [false; 12];
                for start in 0..12 {
                    if next[start] == NONE || seen[start] {
                        continue;
                    }
                    let mut lp: Vec<u8> = Vec::with_capacity(12);
                    let mut e = start as u8;
                    loop {
                        seen[e as usize] = true;
                        lp.push(e);
                        e = next[e as usize];
                        if e == NONE || e as usize == start || seen[e as usize] {
                            break;
                        }
                    }
                    emit_loop(grid, vals, &keys, &lp, edge_limit + base as u64 * 12 + start as u64, &mut out);
                }
            }
        }
    }
    out
}

/// Is the box `[lo, hi]` inside the union of the balls (sufficient test: some single ball covers it,
/// or all eight octants are covered, down to `depth` levels)?
fn covered(lo: [f64; 3], hi: [f64; 3], balls: &[([f64; 3], f64)], depth: usize) -> bool {
    for (c, r) in balls {
        let far2: f64 = (0..3).map(|a| (c[a] - lo[a]).abs().max((c[a] - hi[a]).abs()).powi(2)).sum();
        if far2 <= r * r {
            return true;
        }
    }
    if depth == 0 {
        return false;
    }
    let mid = [0, 1, 2].map(|a| 0.5 * (lo[a] + hi[a]));
    (0..8).all(|o| {
        let l = [0, 1, 2].map(|a| if (o >> a) & 1 == 1 { mid[a] } else { lo[a] });
        let h = [0, 1, 2].map(|a| if (o >> a) & 1 == 1 { hi[a] } else { mid[a] });
        covered(l, h, balls, depth - 1)
    })
}

/// Extract `{φ = isovalue}` from `field`.
///
/// Returns the welded triangle mesh (outward oriented) and a report. The mesh is closed except
/// where the surface leaves the box; for generic data it is an oriented manifold (checked by
/// `TriMesh::new` when [`McOptions::check_mesh`]).
pub fn marching_cubes(field: &dyn ImplicitField, opts: &McOptions) -> Result<(McMesh, McReport), PatError> {
    let t_start = Instant::now();
    let totals_before = field.eval_totals();
    if !opts.isovalue.is_finite() {
        return Err(invalid(format!("marching cubes isovalue = {}", opts.isovalue)));
    }
    if opts.refine_steps > 64 {
        return Err(invalid(format!("refine_steps = {} (at most 64)", opts.refine_steps)));
    }
    let (bound, block) = match opts.band {
        Band::Dense => (None, usize::MAX),
        Band::Lipschitz { bound, block } => {
            if !(bound.is_finite() && bound > 0.0) || block == 0 {
                return Err(invalid(format!("Band::Lipschitz bound = {bound}, block = {block}")));
            }
            (Some(bound), block)
        }
    };
    let pool = make_pool(opts.threads)?;
    let bbox = opts.bbox.unwrap_or_else(|| field.bbox());
    let grid = make_grid(bbox, opts.resolution)?;
    let np = grid.npoints();
    let iso = opts.isovalue;

    // values are stored as φ − iso; NaN: not requested, +inf: queued
    let mut vals = vec![f64::NAN; np];
    // The field is queried through a shifted view.
    struct Shifted<'a>(&'a dyn ImplicitField, f64);
    impl ImplicitField for Shifted<'_> {
        fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
            Ok(self.0.distance(x)? - self.1)
        }
        fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
            let (d, g) = self.0.distance_gradient(x)?;
            Ok((d - self.1, g))
        }
        fn bbox(&self) -> [[f64; 3]; 2] {
            self.0.bbox()
        }
        fn eval_totals(&self) -> Option<crate::EvalReport> {
            self.0.eval_totals()
        }
        fn distances(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<(), PatError> {
            self.0.distances(xs, out)?;
            out.iter_mut().for_each(|v| *v -= self.1);
            Ok(())
        }
        fn distance_gradients(&self, xs: &[[f64; 3]], d: &mut [f64], g: &mut [[f64; 3]]) -> Result<(), PatError> {
            self.0.distance_gradients(xs, d, g)?;
            d.iter_mut().for_each(|v| *v -= self.1);
            Ok(())
        }
    }
    let shifted = Shifted(field, iso);
    let f: &dyn ImplicitField = &shifted;

    // blocks
    let bs = block.min(grid.n.iter().cloned().max().unwrap_or(1)).max(1);
    let nb = [0, 1, 2].map(|a| grid.n[a].div_ceil(bs));
    let n_blocks = nb[0] * nb[1] * nb[2];
    let coarse_idx = |a: usize, b: usize| (b * bs).min(grid.n[a]);

    // coarse corners
    let mut coarse: Vec<usize> = Vec::new();
    for bz in 0..=nb[2] {
        for by in 0..=nb[1] {
            for bx in 0..=nb[0] {
                let g = grid.gp(coarse_idx(0, bx), coarse_idx(1, by), coarse_idx(2, bz));
                if vals[g].is_nan() {
                    vals[g] = f64::INFINITY;
                    coarse.push(g);
                }
            }
        }
    }
    let t_coarse = Instant::now();
    evaluate_points(f, &grid, &coarse, &mut vals, &pool)?;
    let mut grid_values = coarse.len();

    // band
    let mut active: Vec<[usize; 3]> = Vec::new();
    for bz in 0..nb[2] {
        for by in 0..nb[1] {
            for bx in 0..nb[0] {
                active.push([bx, by, bz]);
            }
        }
    }
    if let Some(l) = bound {
        let keep: Vec<bool> = pool.install(|| {
            active
                .par_iter()
                .map(|b| {
                    // The surface cannot lie within |v_c|/L of a corner c; the block is skipped
                    // when these balls cover it (checked on the block and on two levels of octants).
                    let lo = [0, 1, 2].map(|a| grid.lo[a] + coarse_idx(a, b[a]) as f64 * grid.h[a]);
                    let hi = [0, 1, 2].map(|a| grid.lo[a] + coarse_idx(a, b[a] + 1) as f64 * grid.h[a]);
                    let balls: [([f64; 3], f64); 8] = std::array::from_fn(|c| {
                        let o = corner_offset(c);
                        let g = grid.gp(coarse_idx(0, b[0] + o[0]), coarse_idx(1, b[1] + o[1]), coarse_idx(2, b[2] + o[2]));
                        let centre = [0, 1, 2].map(|a| if o[a] == 1 { hi[a] } else { lo[a] });
                        (centre, vals[g].abs() / l)
                    });
                    !covered(lo, hi, &balls, 2)
                })
                .collect()
        });
        let mut i = 0;
        active.retain(|_| {
            i += 1;
            keep[i - 1]
        });
    }
    let block_range = |b: [usize; 3]| ([0, 1, 2].map(|a| coarse_idx(a, b[a])), [0, 1, 2].map(|a| coarse_idx(a, b[a] + 1)));
    let bindex = |b: [usize; 3]| b[0] + nb[0] * (b[1] + nb[1] * b[2]);
    let mut is_active = vec![false; n_blocks];
    for b in &active {
        is_active[bindex(*b)] = true;
    }
    // Fine values of the kept blocks, then growth: a kept block whose face towards a skipped
    // neighbour has a sign change continues the surface into that neighbour, so the neighbour is
    // kept as well (repeat). The Lipschitz bound therefore only decides which *isolated* pieces of
    // the surface can be missed; the surface that is found has no cracks except at the box.
    let t_band = Instant::now();
    let mut newly: Vec<[usize; 3]> = active.clone();
    let mut band_rounds = 0usize;
    let mut blocks_grown = 0usize;
    while !newly.is_empty() {
        band_rounds += 1;
        let mut todo: Vec<usize> = Vec::new();
        for b in &newly {
            let (lo, hi) = block_range(*b);
            for k in lo[2]..=hi[2] {
                for j in lo[1]..=hi[1] {
                    for i in lo[0]..=hi[0] {
                        let g = grid.gp(i, j, k);
                        if vals[g].is_nan() {
                            vals[g] = f64::INFINITY;
                            todo.push(g);
                        }
                    }
                }
            }
        }
        evaluate_points(f, &grid, &todo, &mut vals, &pool)?;
        grid_values += todo.len();
        let mut grow: Vec<[usize; 3]> = Vec::new();
        if bound.is_some() {
            for b in &newly {
                let (lo, hi) = block_range(*b);
                for axis in 0..3 {
                    for side in 0..2 {
                        let nbk = if side == 0 { b[axis].checked_sub(1) } else { Some(b[axis] + 1).filter(|v| *v < nb[axis]) };
                        let Some(nbk) = nbk else { continue };
                        let mut nbb = *b;
                        nbb[axis] = nbk;
                        if is_active[bindex(nbb)] {
                            continue;
                        }
                        // sign change on the face?
                        let fixed = if side == 0 { lo[axis] } else { hi[axis] };
                        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
                        let mut any_neg = false;
                        let mut any_pos = false;
                        for iv in lo[v]..=hi[v] {
                            for iu in lo[u]..=hi[u] {
                                let mut ijk = [0; 3];
                                ijk[axis] = fixed;
                                ijk[u] = iu;
                                ijk[v] = iv;
                                if vals[grid.gp(ijk[0], ijk[1], ijk[2])] < 0.0 {
                                    any_neg = true;
                                } else {
                                    any_pos = true;
                                }
                            }
                        }
                        if any_neg && any_pos {
                            is_active[bindex(nbb)] = true;
                            grow.push(nbb);
                        }
                    }
                }
            }
        }
        blocks_grown += grow.len();
        newly = grow;
    }
    let ranges: Vec<([usize; 3], [usize; 3])> = (0..n_blocks)
        .filter(|bi| is_active[*bi])
        .map(|bi| block_range([bi % nb[0], (bi / nb[0]) % nb[1], bi / (nb[0] * nb[1])]))
        .collect();
    let cells_visited: usize = ranges.iter().map(|(lo, hi)| (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])).sum();
    drop(coarse);
    let t_eval = Instant::now();

    // cells
    let edge_limit = np as u64 * 3;
    let outs: Vec<BlockOut> = pool.install(|| ranges.par_iter().map(|(lo, hi)| process_block(&grid, &vals, edge_limit, *lo, *hi)).collect());
    let t_cells = Instant::now();
    let (mut crossing, mut ambiguous, mut max_slope) = (0, 0, 0.0f64);
    let mut ntri = 0;
    for o in &outs {
        crossing += o.crossing;
        ambiguous += o.ambiguous;
        max_slope = max_slope.max(o.slope);
        ntri += o.tris.len();
    }
    let mut tri_keys: Vec<[u64; 3]> = Vec::with_capacity(ntri);
    let mut centroids: Vec<(u64, [f64; 3])> = Vec::new();
    for o in outs {
        tri_keys.extend(o.tris);
        centroids.extend(o.centroids);
    }
    centroids.sort_by_key(|c| c.0);
    let mut keys: Vec<u64> = tri_keys.iter().flatten().copied().collect();
    pool.install(|| keys.par_sort_unstable());
    keys.dedup();
    if keys.len() >= u32::MAX as usize {
        return Err(invalid(format!("{} vertices do not fit in u32 indices", keys.len())));
    }
    let faces: Vec<[u32; 3]> = pool.install(|| {
        tri_keys
            .par_iter()
            .map(|t| t.map(|k| keys.binary_search(&k).unwrap_or(0) as u32))
            .collect()
    });
    drop(tri_keys);

    // vertices
    let n_edge_keys = keys.partition_point(|k| *k < edge_limit);
    let mut ts: Vec<f64> = pool.install(|| keys[..n_edge_keys].par_iter().map(|k| grid.linear_t(&vals, *k)).collect());
    let mut extra = 0usize;
    if opts.refine_steps > 0 {
        extra += refine_edges(f, &grid, &keys[..n_edge_keys], &vals, &mut ts, opts.refine_steps, &pool)?;
    }
    let mut positions: Vec<[f64; 3]> = pool.install(|| keys[..n_edge_keys].par_iter().zip(ts.par_iter()).map(|(k, t)| grid.vertex(*k, *t)).collect());
    for k in &keys[n_edge_keys..] {
        let c = centroids.binary_search_by_key(k, |c| c.0).map_or([0.0; 3], |i| centroids[i].1);
        positions.push(c);
    }
    let mut zero_normals = 0;
    let normals = if opts.normals {
        let mut d = vec![0.0; positions.len()];
        let mut g = vec![[0.0; 3]; positions.len()];
        const BATCH: usize = 1 << 18;
        for ((xs, dd), gg) in positions.chunks(BATCH).zip(d.chunks_mut(BATCH)).zip(g.chunks_mut(BATCH)) {
            pool.install(|| f.distance_gradients(xs, dd, gg))?;
        }
        extra += positions.len();
        let nrm: Vec<[f64; 3]> = g
            .iter()
            .map(|v| {
                let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                if l > 0.0 && l.is_finite() {
                    [v[0] / l, v[1] / l, v[2] / l]
                } else {
                    zero_normals += 1;
                    [0.0; 3]
                }
            })
            .collect();
        Some(nrm)
    } else {
        None
    };
    let t_end = Instant::now();

    let mesh = McMesh { positions, faces, normals };
    let check = if opts.check_mesh && !mesh.faces.is_empty() {
        match mesh.to_trimesh() {
            Ok(m) => MeshCheck::Manifold {
                boundary_edges: m.n_boundary_edges(),
                euler_characteristic: m.euler_characteristic(),
                components: m.n_components(),
            },
            Err(e) => MeshCheck::Invalid(e.to_string()),
        }
    } else {
        MeshCheck::Skipped
    };
    let report = McReport {
        cells: grid.n,
        spacing: grid.h,
        origin: grid.lo,
        blocks: n_blocks,
        blocks_grown,
        band_rounds,
        blocks_active: ranges.len(),
        cells_visited,
        cells_crossing: crossing,
        ambiguous_cell_faces: ambiguous,
        grid_values,
        grid_points: np,
        extra_evaluations: extra,
        vertices: mesh.positions.len(),
        triangles: mesh.faces.len(),
        lipschitz_bound: bound,
        max_edge_slope: max_slope,
        bound_violated: bound.is_some_and(|l| max_slope > l),
        zero_normals,
        mesh: check,
        seconds: [
            (t_coarse - t_start).as_secs_f64() + (t_band - t_coarse).as_secs_f64() * 0.0,
            (t_eval - t_coarse).as_secs_f64(),
            (t_cells - t_eval).as_secs_f64(),
            (t_end - t_cells).as_secs_f64(),
            (t_end - t_start).as_secs_f64(),
        ],
        eval: field.eval_totals().zip(totals_before).map(|(after, before)| after.since(&before)),
    };
    Ok((mesh, report))
}

/// Safeguarded Newton steps along every edge: keep the bracket `[lo, hi]` of `t`, take the Newton
/// step of `t ↦ φ(p(t))` and fall back to bisection when it leaves the bracket. Returns the
/// number of evaluations.
fn refine_edges(
    f: &dyn ImplicitField,
    grid: &Grid,
    keys: &[u64],
    vals: &[f64],
    ts: &mut [f64],
    steps: usize,
    pool: &rayon::ThreadPool,
) -> Result<usize, PatError> {
    const BATCH: usize = 1 << 18;
    struct St {
        lo: f64,
        hi: f64,
        neg_lo: bool,
    }
    let mut st: Vec<St> = keys
        .iter()
        .map(|k| {
            let gp = (*k / 3) as usize;
            St { lo: 0.0, hi: 1.0, neg_lo: vals[gp] < 0.0 }
        })
        .collect();
    let mut evals = 0;
    let mut xs: Vec<[f64; 3]> = Vec::new();
    let mut d: Vec<f64> = Vec::new();
    let mut g: Vec<[f64; 3]> = Vec::new();
    for _ in 0..steps {
        for start in (0..keys.len()).step_by(BATCH) {
            let end = (start + BATCH).min(keys.len());
            xs.clear();
            xs.extend((start..end).map(|i| grid.vertex(keys[i], ts[i])));
            d.clear();
            d.resize(end - start, 0.0);
            g.clear();
            g.resize(end - start, [0.0; 3]);
            pool.install(|| f.distance_gradients(&xs, &mut d, &mut g))?;
            evals += end - start;
            for i in start..end {
                let (fv, gv) = (d[i - start], g[i - start]);
                if !fv.is_finite() {
                    return Err(PatError::NonFinite { index: i, what: "field value (edge refinement)" });
                }
                let axis = (keys[i] % 3) as usize;
                let s = &mut st[i];
                if fv == 0.0 {
                    s.lo = ts[i];
                    s.hi = ts[i];
                    continue;
                }
                if (fv < 0.0) == s.neg_lo {
                    s.lo = ts[i];
                } else {
                    s.hi = ts[i];
                }
                let slope = gv[axis] * grid.h[axis];
                let newton = ts[i] - fv / slope;
                ts[i] = if newton.is_finite() && newton > s.lo && newton < s.hi { newton } else { 0.5 * (s.lo + s.hi) };
            }
        }
    }
    Ok(evals)
}
