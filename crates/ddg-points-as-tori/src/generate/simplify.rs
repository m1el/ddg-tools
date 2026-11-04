//! Quadric edge-collapse simplification of closed triangle meshes (Garland and Heckbert, "Surface
//! simplification using quadric error metrics", SIGGRAPH 1997), standing in for the reference's
//! `pyfqmr` (Fast-Quadric-Mesh-Simplification) step that brings blobs down to 2048 faces.
//!
//! A collapse is refused unless it keeps the mesh a closed manifold (link condition) and does not
//! flip or collapse any incident face, so a closed oriented manifold stays one. Placement: the
//! minimizer of the summed quadric when well conditioned, else the best of the two endpoints and
//! the midpoint. Deterministic (ties broken by vertex indices).

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

type Q = [f64; 10]; // symmetric 4x4: a b c d / e f g / h i / j

fn plane_q(n: [f64; 3], d: f64, w: f64) -> Q {
    [n[0] * n[0], n[0] * n[1], n[0] * n[2], n[0] * d, n[1] * n[1], n[1] * n[2], n[1] * d, n[2] * n[2], n[2] * d, d * d].map(|x| x * w)
}

fn add(a: &Q, b: &Q) -> Q {
    let mut r = [0.0; 10];
    for i in 0..10 {
        r[i] = a[i] + b[i];
    }
    r
}

fn eval(q: &Q, p: [f64; 3]) -> f64 {
    let (x, y, z) = (p[0], p[1], p[2]);
    q[0] * x * x + 2.0 * q[1] * x * y + 2.0 * q[2] * x * z + 2.0 * q[3] * x + q[4] * y * y + 2.0 * q[5] * y * z + 2.0 * q[6] * y + q[7] * z * z + 2.0 * q[8] * z + q[9]
}

fn optimum(q: &Q) -> Option<[f64; 3]> {
    let a = [[q[0], q[1], q[2]], [q[1], q[4], q[5]], [q[2], q[5], q[7]]];
    let b = [-q[3], -q[6], -q[8]];
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    let scale = (a[0][0].abs() + a[1][1].abs() + a[2][2].abs()).max(1e-300);
    if det.abs() < 1e-9 * scale * scale * scale {
        return None;
    }
    let inv = |m: [[f64; 3]; 3], rhs: [f64; 3]| -> [f64; 3] {
        let col = |c: usize| {
            let mut mm = m;
            for r in 0..3 {
                mm[r][c] = rhs[r];
            }
            mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1]) - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0]) + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0])
        };
        [col(0) / det, col(1) / det, col(2) / det]
    };
    let x = inv(a, b);
    x.iter().all(|c| c.is_finite()).then_some(x)
}

#[derive(PartialEq)]
struct Cand {
    cost: f64,
    u: u32,
    v: u32,
    ver: (u32, u32),
    pos: [f64; 3],
}
impl Eq for Cand {}
impl Ord for Cand {
    fn cmp(&self, o: &Self) -> Ordering {
        // min-heap on (cost, u, v)
        o.cost.total_cmp(&self.cost).then(o.u.cmp(&self.u)).then(o.v.cmp(&self.v))
    }
}
impl PartialOrd for Cand {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

fn normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let (u, w) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
    [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
}

/// Simplify a closed manifold mesh to at most `target_faces` faces (fewer collapses happen if no
/// legal collapse remains). Returns compacted vertices and faces.
pub fn simplify_closed(vertices: &[[f64; 3]], faces: &[[u32; 3]], target_faces: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let nv = vertices.len();
    let mut pos = vertices.to_vec();
    let mut tri: Vec<[u32; 3]> = faces.to_vec();
    let mut alive = vec![true; tri.len()];
    let mut nalive = tri.len();
    let mut vfaces: Vec<Vec<u32>> = vec![Vec::new(); nv];
    let mut q: Vec<Q> = vec![[0.0; 10]; nv];
    for (f, t) in tri.iter().enumerate() {
        for &v in t {
            vfaces[v as usize].push(f as u32);
        }
        let (a, b, c) = (pos[t[0] as usize], pos[t[1] as usize], pos[t[2] as usize]);
        let n = normal(a, b, c);
        let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if l > 0.0 {
            let nn = [n[0] / l, n[1] / l, n[2] / l];
            let d = -(nn[0] * a[0] + nn[1] * a[1] + nn[2] * a[2]);
            let pq = plane_q(nn, d, 0.5 * l);
            for &v in t {
                q[v as usize] = add(&q[v as usize], &pq);
            }
        }
    }
    let mut version = vec![0u32; nv];
    let mut vdead = vec![false; nv];
    let mut heap = BinaryHeap::new();
    let edge_cost = |q: &Vec<Q>, pos: &Vec<[f64; 3]>, u: u32, v: u32, version: &Vec<u32>| -> Cand {
        let qs = add(&q[u as usize], &q[v as usize]);
        let (pu, pv) = (pos[u as usize], pos[v as usize]);
        let mid = [0, 1, 2].map(|k| 0.5 * (pu[k] + pv[k]));
        let mut best = (eval(&qs, mid), mid);
        for c in [pu, pv] {
            let e = eval(&qs, c);
            if e < best.0 {
                best = (e, c);
            }
        }
        if let Some(o) = optimum(&qs) {
            // accept the optimum only when it stays near the edge
            let el = ((pu[0] - pv[0]).powi(2) + (pu[1] - pv[1]).powi(2) + (pu[2] - pv[2]).powi(2)).sqrt();
            let dm = ((o[0] - mid[0]).powi(2) + (o[1] - mid[1]).powi(2) + (o[2] - mid[2]).powi(2)).sqrt();
            let e = eval(&qs, o);
            if dm < 2.0 * el && e <= best.0 {
                best = (e, o);
            }
        }
        Cand { cost: best.0.max(0.0), u, v, ver: (version[u as usize], version[v as usize]), pos: best.1 }
    };
    for t in &tri {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if a < b {
                heap.push(edge_cost(&q, &pos, a, b, &version));
            }
        }
    }
    // each undirected edge is pushed from the face where it runs low->high; the other face
    // pushes it when it runs high->low, so push those too (duplicates are harmless)
    for t in &tri {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if a > b {
                heap.push(edge_cost(&q, &pos, b, a, &version));
            }
        }
    }
    let neighbors = |vfaces: &Vec<Vec<u32>>, tri: &Vec<[u32; 3]>, alive: &Vec<bool>, v: u32| -> HashSet<u32> {
        let mut s = HashSet::new();
        for &f in &vfaces[v as usize] {
            if alive[f as usize] {
                for &w in &tri[f as usize] {
                    if w != v {
                        s.insert(w);
                    }
                }
            }
        }
        s
    };
    while nalive > target_faces {
        let Some(c) = heap.pop() else { break };
        let (u, v) = (c.u, c.v);
        if vdead[u as usize] || vdead[v as usize] || c.ver != (version[u as usize], version[v as usize]) {
            continue;
        }
        // faces containing both
        let shared: Vec<u32> = vfaces[u as usize].iter().copied().filter(|&f| alive[f as usize] && tri[f as usize].contains(&v)).collect();
        if shared.len() != 2 {
            continue;
        }
        // link condition
        let (nu, nv_) = (neighbors(&vfaces, &tri, &alive, u), neighbors(&vfaces, &tri, &alive, v));
        let common: HashSet<u32> = nu.intersection(&nv_).copied().collect();
        let opp: HashSet<u32> = shared.iter().flat_map(|&f| tri[f as usize].iter().copied().filter(|&w| w != u && w != v)).collect();
        if common != opp || nu.len() + nv_.len() < 8 {
            continue;
        }
        // flip / degeneracy check on all other incident faces
        let mut ok = true;
        for &w in [u, v].iter() {
            for &f in &vfaces[w as usize] {
                if !alive[f as usize] || shared.contains(&f) {
                    continue;
                }
                let t = tri[f as usize];
                let old = normal(pos[t[0] as usize], pos[t[1] as usize], pos[t[2] as usize]);
                let np = |x: u32| if x == u || x == v { c.pos } else { pos[x as usize] };
                let new = normal(np(t[0]), np(t[1]), np(t[2]));
                let (lo, ln) = ((old[0] * old[0] + old[1] * old[1] + old[2] * old[2]).sqrt(), (new[0] * new[0] + new[1] * new[1] + new[2] * new[2]).sqrt());
                if ln <= 1e-12 * lo.max(1e-300) || (old[0] * new[0] + old[1] * new[1] + old[2] * new[2]) <= 0.2 * lo * ln {
                    ok = false;
                    break;
                }
            }
            if !ok {
                break;
            }
        }
        if !ok {
            continue;
        }
        // collapse v into u
        for &f in &shared {
            alive[f as usize] = false;
            nalive -= 1;
        }
        let vf = std::mem::take(&mut vfaces[v as usize]);
        for f in vf {
            if alive[f as usize] {
                for x in tri[f as usize].iter_mut() {
                    if *x == v {
                        *x = u;
                    }
                }
                vfaces[u as usize].push(f);
            }
        }
        vfaces[u as usize].retain(|&f| alive[f as usize]);
        pos[u as usize] = c.pos;
        q[u as usize] = add(&q[u as usize], &q[v as usize]);
        vdead[v as usize] = true;
        version[u as usize] += 1;
        for w in neighbors(&vfaces, &tri, &alive, u) {
            let (a, b) = if u < w { (u, w) } else { (w, u) };
            heap.push(edge_cost(&q, &pos, a, b, &version));
        }
    }
    // compact
    let mut remap = vec![u32::MAX; nv];
    let mut out_v = Vec::new();
    let mut out_f = Vec::new();
    for (f, t) in tri.iter().enumerate() {
        if !alive[f] {
            continue;
        }
        let mut nt = [0u32; 3];
        for k in 0..3 {
            let v = t[k] as usize;
            if remap[v] == u32::MAX {
                remap[v] = out_v.len() as u32;
                out_v.push(pos[v]);
            }
            nt[k] = remap[v];
        }
        out_f.push(nt);
    }
    (out_v, out_f)
}
