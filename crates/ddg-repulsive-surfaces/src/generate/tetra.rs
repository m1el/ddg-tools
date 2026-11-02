//! Marching tetrahedra on a regular grid (Kuhn decomposition: six tetrahedra per cube around the
//! main diagonal, which is conforming across cubes). The surface is a manifold by construction
//! (vertices lie on grid edges, one per sign change; the result is validated by [`Mesh::validated`]).
//!
//! Cleanup (no sliver triangles): an edge vertex is kept at least [`CLAMP`]`·h` away from the
//! grid vertex it would hug, then a few rounds of Laplacian smoothing, each followed by
//! Newton steps back onto the zero set of the field. Triangles stay oriented outward (normal
//! from the negative to the positive side of the field).

use std::collections::HashMap;

use rayon::prelude::*;

use super::Mesh;
use super::field::Field;
use crate::error::RsError;
use crate::geometry::v3::{add, cross, dot, norm, scale, sub};

/// Minimal distance, in cell widths along a grid edge, of a surface vertex from a grid vertex.
pub const CLAMP: f64 = 0.15;
const BLOCK: usize = 4;

/// Options of the smoothing pass.
#[derive(Debug, Clone, Copy)]
pub struct Cleanup {
    /// Rounds of smoothing (0: none).
    pub rounds: usize,
    /// Step towards the neighbour average.
    pub lambda: f64,
}

impl Default for Cleanup {
    fn default() -> Self {
        Cleanup {
            rounds: 6,
            lambda: 0.5,
        }
    }
}

/// Meshes `{f = 0}` inside the box `[lo, hi]` with cell width `h`. The field must be positive on
/// the boundary of the (rounded-up) grid.
///
/// # Errors
/// [`RsError::InvalidOption`] for `h ≤ 0`, an empty box, a grid over 600 cells per axis or a
/// field that is negative on the grid boundary; mesh errors if the result is not a valid surface.
pub fn marching_tetrahedra(
    f: &Field,
    lo: [f64; 3],
    hi: [f64; 3],
    h: f64,
    cleanup: Cleanup,
) -> Result<Mesh, RsError> {
    if !(h.is_finite() && h > 0.0) {
        return Err(super::bad("h", h, "a positive cell width"));
    }
    let mut n = [0usize; 3];
    for a in 0..3 {
        let c = ((hi[a] - lo[a]) / h).ceil();
        if !(2.0..=600.0).contains(&c) {
            return Err(super::bad("grid", c, "between 2 and 600 cells per axis"));
        }
        n[a] = c as usize;
    }
    let nv = [n[0] + 1, n[1] + 1, n[2] + 1];
    let vid = |i: usize, j: usize, k: usize| (i * nv[1] + j) * nv[2] + k;
    let pos = |i: usize, j: usize, k: usize| {
        [
            lo[0] + i as f64 * h,
            lo[1] + j as f64 * h,
            lo[2] + k as f64 * h,
        ]
    };
    let total = nv[0] * nv[1] * nv[2];
    // Values: far blocks get the value at their centre (right sign, no crossing inside), near
    // blocks are evaluated exactly.
    let mut val = vec![f64::NAN; total];
    let nb = [
        n[0].div_ceil(BLOCK),
        n[1].div_ceil(BLOCK),
        n[2].div_ceil(BLOCK),
    ];
    let blocks: Vec<(usize, usize, usize)> = (0..nb[0])
        .flat_map(|a| (0..nb[1]).flat_map(move |b| (0..nb[2]).map(move |c| (a, b, c))))
        .collect();
    let reach = 0.5 * BLOCK as f64 * h * 3f64.sqrt() * 1.3;
    let centers: Vec<f64> = blocks
        .par_iter()
        .map(|&(a, b, c)| {
            let p = pos(a * BLOCK + BLOCK / 2, b * BLOCK + BLOCK / 2, c * BLOCK + BLOCK / 2);
            f(p)
        })
        .collect();
    let mut near = vec![false; total];
    for (bi, &(a, b, c)) in blocks.iter().enumerate() {
        let near_block = centers[bi].abs() <= reach || !centers[bi].is_finite();
        for i in a * BLOCK..((a + 1) * BLOCK + 1).min(nv[0]) {
            for j in b * BLOCK..((b + 1) * BLOCK + 1).min(nv[1]) {
                for k in c * BLOCK..((c + 1) * BLOCK + 1).min(nv[2]) {
                    let v = vid(i, j, k);
                    if near_block {
                        near[v] = true;
                    } else if !near[v] {
                        val[v] = centers[bi];
                    }
                }
            }
        }
    }
    let todo: Vec<usize> = (0..total).filter(|&v| near[v]).collect();
    let computed: Vec<f64> = todo
        .par_iter()
        .map(|&v| {
            let k = v % nv[2];
            let j = (v / nv[2]) % nv[1];
            let i = v / (nv[2] * nv[1]);
            f(pos(i, j, k))
        })
        .collect();
    for (&v, &x) in todo.iter().zip(&computed) {
        val[v] = x;
    }
    // boundary check
    for i in 0..nv[0] {
        for j in 0..nv[1] {
            for k in 0..nv[2] {
                let edge = i == 0 || j == 0 || k == 0 || i == n[0] || j == n[1] || k == n[2];
                if edge && val[vid(i, j, k)] < 0.0 {
                    return Err(super::bad(
                        "field",
                        format!("negative at grid vertex ({i}, {j}, {k})"),
                        "positive on the boundary of the grid (enlarge the box)",
                    ));
                }
            }
        }
    }

    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut faces: Vec<[u32; 3]> = Vec::new();
    let mut edge_vertex: HashMap<(u32, u32), u32> = HashMap::new();
    let perms: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                let mut any_in = false;
                let mut any_out = false;
                for c in 0..8 {
                    let v = val[vid(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1))];
                    if v < 0.0 {
                        any_in = true;
                    } else {
                        any_out = true;
                    }
                }
                if !(any_in && any_out) {
                    continue;
                }
                for perm in perms {
                    let mut idx = [[i, j, k]; 4];
                    let mut cur = [i, j, k];
                    for (s, &ax) in perm.iter().enumerate() {
                        cur[ax] += 1;
                        idx[s + 1] = cur;
                    }
                    let ids: [u32; 4] = idx.map(|c| vid(c[0], c[1], c[2]) as u32);
                    let ps: [[f64; 3]; 4] = idx.map(|c| pos(c[0], c[1], c[2]));
                    let fv: [f64; 4] = ids.map(|v| val[v as usize]);
                    let inside: Vec<usize> = (0..4).filter(|&s| fv[s] < 0.0).collect();
                    let outside: Vec<usize> = (0..4).filter(|&s| fv[s] >= 0.0).collect();
                    if inside.is_empty() || outside.is_empty() {
                        continue;
                    }
                    let mut ev = |a: usize, b: usize| -> u32 {
                        let key = if ids[a] < ids[b] {
                            (ids[a], ids[b])
                        } else {
                            (ids[b], ids[a])
                        };
                        *edge_vertex.entry(key).or_insert_with(|| {
                            let (ai, bi) = if fv[a] < 0.0 { (a, b) } else { (b, a) };
                            let t = (fv[ai] / (fv[ai] - fv[bi])).clamp(CLAMP, 1.0 - CLAMP);
                            let p = add(ps[ai], scale(sub(ps[bi], ps[ai]), t));
                            positions.push(p);
                            (positions.len() - 1) as u32
                        })
                    };
                    let centroid = |set: &[usize]| {
                        let mut c = [0.0; 3];
                        for &s in set {
                            c = add(c, ps[s]);
                        }
                        scale(c, 1.0 / set.len() as f64)
                    };
                    let dir = sub(centroid(&outside), centroid(&inside));
                    let mut tris: Vec<[u32; 3]> = Vec::new();
                    match (inside.len(), outside.len()) {
                        (1, 3) => {
                            let x = inside[0];
                            tris.push([ev(x, outside[0]), ev(x, outside[1]), ev(x, outside[2])]);
                        }
                        (3, 1) => {
                            let x = outside[0];
                            tris.push([ev(x, inside[0]), ev(x, inside[1]), ev(x, inside[2])]);
                        }
                        _ => {
                            let (i0, i1, o0, o1) = (inside[0], inside[1], outside[0], outside[1]);
                            let e00 = ev(i0, o0);
                            let e01 = ev(i0, o1);
                            let e11 = ev(i1, o1);
                            let e10 = ev(i1, o0);
                            tris.push([e00, e01, e11]);
                            tris.push([e00, e11, e10]);
                        }
                    }
                    // orientation from the first triangle (the quad's two halves agree)
                    let t0 = tris[0];
                    let nrm = cross(
                        sub(positions[t0[1] as usize], positions[t0[0] as usize]),
                        sub(positions[t0[2] as usize], positions[t0[0] as usize]),
                    );
                    let flip = dot(nrm, dir) < 0.0;
                    for t in tris {
                        faces.push(if flip { [t[0], t[2], t[1]] } else { t });
                    }
                }
            }
        }
    }
    if faces.is_empty() {
        return Err(super::bad("field", "no sign change", "a surface inside the box"));
    }
    let mut mesh = Mesh { positions, faces };
    smooth_and_project(&mut mesh, f, h, cleanup);
    mesh.validated()
}

/// Laplacian smoothing, each round followed by two Newton steps onto `{f = 0}`.
pub fn smooth_and_project(mesh: &mut Mesh, f: &Field, h: f64, cleanup: Cleanup) {
    let n = mesh.positions.len();
    let mut nbrs: Vec<Vec<u32>> = vec![Vec::new(); n];
    for t in &mesh.faces {
        for c in 0..3 {
            let (a, b) = (t[c], t[(c + 1) % 3]);
            nbrs[a as usize].push(b);
            nbrs[b as usize].push(a);
        }
    }
    for l in &mut nbrs {
        l.sort_unstable();
        l.dedup();
    }
    let eps = 0.25 * h;
    let project = |p: [f64; 3]| -> [f64; 3] {
        let mut p = p;
        for _ in 0..2 {
            let v = f(p);
            let g = [
                f(add(p, [eps, 0.0, 0.0])) - f(sub(p, [eps, 0.0, 0.0])),
                f(add(p, [0.0, eps, 0.0])) - f(sub(p, [0.0, eps, 0.0])),
                f(add(p, [0.0, 0.0, eps])) - f(sub(p, [0.0, 0.0, eps])),
            ]
            .map(|x| x / (2.0 * eps));
            let g2 = dot(g, g);
            if g2 > 1e-12 && v.is_finite() {
                // never move more than a cell
                let step = scale(g, v / g2);
                let m = norm(step);
                p = sub(p, if m > h { scale(step, h / m) } else { step });
            }
        }
        p
    };
    for _ in 0..cleanup.rounds {
        let old = mesh.positions.clone();
        let new: Vec<[f64; 3]> = (0..n)
            .into_par_iter()
            .map(|v| {
                let mut avg = [0.0; 3];
                for &w in &nbrs[v] {
                    avg = add(avg, old[w as usize]);
                }
                avg = scale(avg, 1.0 / nbrs[v].len() as f64);
                let p = add(old[v], scale(sub(avg, old[v]), cleanup.lambda));
                project(p)
            })
            .collect();
        mesh.positions = new;
    }
}

/// Meshes `{f = 0}` with about `target` faces by adjusting the cell width (up to five tries;
/// the closest result is returned). The grid is anchored at `lo`.
///
/// # Errors
/// As [`marching_tetrahedra`]; `target < 100` is [`RsError::InvalidOption`].
pub fn mesh_with_faces(
    f: &Field,
    lo: [f64; 3],
    hi: [f64; 3],
    target: usize,
) -> Result<Mesh, RsError> {
    if target < 100 {
        return Err(super::bad("faces", target, "at least 100"));
    }
    let ext = (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max);
    let mut h = ext / 60.0;
    let mut best: Option<Mesh> = None;
    for _ in 0..5 {
        let m = marching_tetrahedra(f, lo, hi, h, Cleanup::default())?;
        let nf = m.faces.len() as f64;
        let better = best.as_ref().is_none_or(|b| {
            (nf - target as f64).abs() < (b.faces.len() as f64 - target as f64).abs()
        });
        let done = ((nf - target as f64) / target as f64).abs() < 0.01;
        if better {
            best = Some(m);
        }
        if done {
            break;
        }
        h *= (nf / target as f64).sqrt().clamp(0.5, 2.0);
    }
    best.ok_or_else(|| super::bad("faces", target, "reachable"))
}
