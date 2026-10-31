//! Fig. 14 of the paper (Sec. 4.3): the accuracy of a Poisson solution on low-quality meshes of
//! the unit square, on the original mesh, on its intrinsic Delaunay triangulation (IDT) copied
//! back to the original mesh, and on the IDT transferred by the `L²` projection (Eq. 13-14).
//!
//! **Setup** (the caption of Fig. 14; the generator and the seeds are ours). `Δf = g`,
//! `g = sin(πx) sin(πy)` on `[0,1]²`, `f = 0` on the boundary; exact solution
//! `f = −g / (2π²)`. With the workspace's positive semidefinite cotan Laplacian `L` (`Δ ≈ −M⁻¹L`)
//! the discrete problem is `L f = −M g` with the Galerkin mass matrix `M`, `g` sampled at the
//! vertices, Dirichlet rows removed. Meshes: the square cut by a diagonal, then random edge
//! splits until `2^(k/2)` vertices for `k = 12 … 22` (the figure's x axis: `2⁶ … 2¹¹` in half
//! octaves; the paper marks about `2^7.5`, k = 15); 100 trials per size, trial `t` of size `k` uses
//! the seed `1000 k + t` of `SplitMix64`. The generator is ours (the paper gives none): a uniform
//! random point of the square picks its face, the edge of that face nearest to the point is split at
//! the projection of the point onto it (clamped to `[0.05, 0.95]`; boundary edges are split too).
//! `cargo run … fig14 100 1` instead splits the longest of 1 random edges (a random edge, at a
//! uniform position; `… 100 4` the longest of 4). The IDT has the same vertices; its `L` and `M` are assembled
//! here from the edge lengths of the integer triangulation (no `TriMesh`, no `ddg-intrinsic-apps`).
//!
//! **What is measured**, all as `‖f − f̃‖_{L²(unit square)}` (absolute, `‖f‖ = 0.0253`) by a
//! degree-4 quadrature on the pieces of the common subdivision `S` (which refines the original
//! mesh, so every function below is piecewise linear on the pieces):
//!
//! 1. `original`: `f̃` = the P1 solution on the original mesh `T⁰` (`P₀ f₀`);
//! 2. `IDT - copy values`: the IDT solution's nodal values used as nodal values of `T⁰`
//!    (same vertices, same ids), `P₀ f_copy`;
//! 3. `IDT - L²`: `P₀ f̂` with `f̂ = Transfer::to_input(f_IDT)` (Eq. 14);
//! 4. extra, `IDT itself`: the IDT's own piecewise linear function `P₁ f_IDT` (not in the paper).
//!
//! The paper plots the average over the 100 trials; this prints the mean and the median.
//!
//! ```text
//! cargo run --release -p ddg-integer-coords --example fig14 [trials [generator]]
//! ```

use std::f64::consts::PI;

use ddg_core::generators::SplitMix64;
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::{Geometry, SparseMatrix, TriMesh};
use ddg_integer_coords::{IntegerTriangulation, Transfer};
use ddg_intrinsic::{DelaunayOptions, Triangulation};

/// The unit square cut by a diagonal, then random edge splits until there are `n` vertices.
fn random_square(n: usize, seed: u64, tournament: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    use std::collections::HashMap;
    let mut rng = SplitMix64::new(seed);
    let mut pos = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let mut faces: Vec<[u32; 3]> = vec![[0, 1, 2], [0, 2, 3]];
    // Directed edge a -> b to the face that has it counter-clockwise.
    let mut he: HashMap<(u32, u32), usize> = HashMap::new();
    for (f, t) in faces.iter().enumerate() {
        for k in 0..3 {
            he.insert((t[k], t[(k + 1) % 3]), f);
        }
    }
    let mut edges: Vec<(u32, u32)> = vec![(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)];
    while pos.len() < n {
        let len2 = |e: usize, pos: &[[f64; 3]]| {
            let (p, q) = (pos[edges[e].0 as usize], pos[edges[e].1 as usize]);
            (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)
        };
        let mut e: usize;
        let t: f64;
        if tournament == 0 {
            // A uniform random point of the square: the edge of its face nearest to it is split
            // at the projection of the point onto it.
            let p = [rng.next_f64(), rng.next_f64()];
            let mut best = (f64::INFINITY, 0usize, 0.5);
            for tri in &faces {
                let q: [[f64; 2]; 3] = tri.map(|v| [pos[v as usize][0], pos[v as usize][1]]);
                let side = |a: [f64; 2], b: [f64; 2]| {
                    (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
                };
                if side(q[0], q[1]) < 0.0 || side(q[1], q[2]) < 0.0 || side(q[2], q[0]) < 0.0 {
                    continue;
                }
                for k in 0..3 {
                    let (a, b) = (q[k], q[(k + 1) % 3]);
                    let d = [b[0] - a[0], b[1] - a[1]];
                    let s = (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1])
                        / (d[0] * d[0] + d[1] * d[1]))
                        .clamp(0.05, 0.95);
                    let dist = (a[0] + s * d[0] - p[0]).hypot(a[1] + s * d[1] - p[1]);
                    if dist < best.0 {
                        let (u, v) = (tri[k], tri[(k + 1) % 3]);
                        let idx = edges
                            .iter()
                            .position(|&x| x == (u, v) || x == (v, u))
                            .unwrap();
                        // Parameter along the stored direction of the edge.
                        let along = if edges[idx] == (u, v) { s } else { 1.0 - s };
                        best = (dist, idx, along);
                    }
                }
                break;
            }
            (e, t) = (best.1, best.2);
        } else {
            e = (rng.next_u64() % edges.len() as u64) as usize;
            for _ in 1..tournament {
                let c = (rng.next_u64() % edges.len() as u64) as usize;
                if len2(c, &pos) > len2(e, &pos) {
                    e = c;
                }
            }
            t = rng.uniform(0.1, 0.9);
        }
        let (a, b) = edges[e];
        let m = pos.len() as u32;
        let (pa, pb) = (pos[a as usize], pos[b as usize]);
        pos.push(std::array::from_fn(|c| (1.0 - t) * pa[c] + t * pb[c]));
        edges[e] = (a, m);
        edges.push((m, b));
        // The face on the left of p -> q is cut into (p, m, c) and (m, q, c).
        for (p, q) in [(a, b), (b, a)] {
            let Some(f) = he.remove(&(p, q)) else {
                continue;
            };
            let old = faces[f];
            let c = old.iter().copied().find(|&v| v != p && v != q).unwrap();
            for k in 0..3 {
                he.remove(&(old[k], old[(k + 1) % 3]));
            }
            faces[f] = [p, m, c];
            faces.push([m, q, c]);
            let g = faces.len() - 1;
            for (id, t) in [(f, faces[f]), (g, faces[g])] {
                for k in 0..3 {
                    he.insert((t[k], t[(k + 1) % 3]), id);
                }
            }
            edges.push((m, c));
        }
    }
    (pos, faces)
}

/// Cotan Laplacian (positive semidefinite) and Galerkin mass matrix from the edge lengths of the
/// faces of a triangulation given as `faces` and the lengths `l[f][k]` of the side from corner
/// `k` to corner `k + 1`.
fn assemble(
    n: usize,
    faces: &[[u32; 3]],
    lengths: impl Fn(usize) -> [f64; 3],
) -> (SparseMatrix, SparseMatrix) {
    let mut lt = Vec::new();
    let mut mt = Vec::new();
    for (f, t) in faces.iter().enumerate() {
        let l = lengths(f);
        let s = 0.5 * (l[0] + l[1] + l[2]);
        let area = (s * (s - l[0]) * (s - l[1]) * (s - l[2])).max(0.0).sqrt();
        for c in 0..3 {
            // Corner c is opposite the side l[(c + 1) % 3].
            let (opp, s1, s2) = (l[(c + 1) % 3], l[c], l[(c + 2) % 3]);
            let cot = (s1 * s1 + s2 * s2 - opp * opp) / (4.0 * area);
            let (a, b) = (t[(c + 1) % 3], t[(c + 2) % 3]);
            let w = 0.5 * cot;
            lt.extend([(a, a, w), (b, b, w), (a, b, -w), (b, a, -w)]);
            for d in 0..3 {
                mt.push((t[c], t[d], area / 12.0 * if c == d { 2.0 } else { 1.0 }));
            }
        }
    }
    (
        SparseMatrix::from_triplets(n, n, &lt).expect("L"),
        SparseMatrix::from_triplets(n, n, &mt).expect("M"),
    )
}

/// `L f = -M g` with `f = 0` at the boundary vertices.
fn solve(l: &SparseMatrix, m: &SparseMatrix, g: &[f64], boundary: &[bool]) -> Vec<f64> {
    let inner: Vec<u32> = (0..g.len() as u32)
        .filter(|&v| !boundary[v as usize])
        .collect();
    let a = l.submatrix(&inner, &inner).expect("submatrix");
    let mut mg = vec![0.0; g.len()];
    m.mul_vec(g, &mut mg).unwrap();
    let mut rhs: Vec<f64> = inner.iter().map(|&v| -mg[v as usize]).collect();
    let fact = Symbolic::analyze(&a, Method::Llt)
        .and_then(|s| s.factor(&a, Parallelism::Sequential))
        .expect("Cholesky");
    fact.solve_in_place(&mut rhs, Parallelism::Sequential)
        .expect("solve");
    let mut f = vec![0.0; g.len()];
    for (&v, x) in inner.iter().zip(rhs) {
        f[v as usize] = x;
    }
    f
}

struct Trial {
    errors: [f64; 4],
    n_s: usize,
    transfer_residual: f64,
}

fn trial(n: usize, seed: u64, tournament: usize) -> Trial {
    let (pos, faces) = random_square(n, seed, tournament);
    let mesh = TriMesh::new(pos.len(), faces).expect("mesh");
    let geometry = Geometry::from_positions(&mesh, &pos).expect("geometry");
    let boundary: Vec<bool> = (0..n).map(|v| mesh.is_boundary_vertex(v)).collect();
    let g_of = |p: [f64; 3]| (PI * p[0]).sin() * (PI * p[1]).sin();
    let exact = |p: [f64; 3]| -g_of(p) / (2.0 * PI * PI);
    let g: Vec<f64> = pos.iter().map(|&p| g_of(p)).collect();

    // (1) The original mesh.
    let (l0, m0) = assemble(n, mesh.faces(), |f| {
        let p = mesh.faces()[f];
        let d = |a: u32, b: u32| {
            let (u, v) = (pos[a as usize], pos[b as usize]);
            ((u[0] - v[0]).powi(2) + (u[1] - v[1]).powi(2)).sqrt()
        };
        [d(p[0], p[1]), d(p[1], p[2]), d(p[2], p[0])]
    });
    let f_orig = solve(&l0, &m0, &g, &boundary);

    // (2), (3): the IDT.
    let mut tri = IntegerTriangulation::new(&mesh, &geometry).expect("init");
    tri.flip_to_delaunay(&DelaunayOptions::default())
        .expect("delaunay");
    let (l1, m1) = assemble(n, tri.complex().faces(), |f| tri.face_lengths(f));
    let f_idt = solve(&l1, &m1, &g, &boundary);
    let tr = Transfer::from_state(&tri).expect("transfer");
    let (f_hat, rep) = tr.to_input_report(&f_idt).expect("to_input");
    let f_copy = tr.copy_to_input(&f_idt).expect("copy");
    let spos = tr.subdivision_positions(&pos).expect("positions");
    let ns = spos.len();
    let err = |p: &SparseMatrix, x: &[f64]| {
        let mut y = vec![0.0; ns];
        p.mul_vec(x, &mut y).unwrap();
        tr.l2_error_to_function(&spos, &y, exact).expect("error")
    };
    Trial {
        errors: [
            err(tr.p0(), &f_orig),
            err(tr.p0(), &f_copy),
            err(tr.p0(), &f_hat),
            err(tr.p1(), &f_idt),
        ],
        n_s: ns,
        transfer_residual: rep.residual,
    }
}

fn par_map<T: Sync, R: Send>(items: &[T], threads: usize, f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let next = std::sync::atomic::AtomicUsize::new(0);
    let out: Vec<std::sync::Mutex<Option<R>>> =
        items.iter().map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if i >= items.len() {
                        break;
                    }
                    *out[i].lock().unwrap() = Some(f(&items[i]));
                }
            });
        }
    });
    out.into_iter()
        .map(|m| m.into_inner().unwrap().unwrap())
        .collect()
}

fn main() {
    let trials: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(100);
    let tournament: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    println!(
        "Poisson on the unit square, {trials} trials per size; absolute L2 errors of the solution on the original mesh"
    );
    println!(
        "{:>5} | {:>10} {:>10} {:>10} {:>10} | {:>8} {:>8} | {:>8} {:>8} | {:>8} {:>9}",
        "|V|",
        "original",
        "IDT copy",
        "IDT L2",
        "IDT own",
        "orig/cpy",
        "orig/L2",
        "med o/c",
        "med o/L2",
        "mean |S|",
        "max resid"
    );
    for k in 12..=22u32 {
        let n = (2f64.powf(f64::from(k) / 2.0)).round() as usize;
        let specs: Vec<u64> = (0..trials as u64)
            .map(|t| 1000 * u64::from(k) + t)
            .collect();
        let results = par_map(&specs, 12, |&seed| trial(n, seed, tournament));
        let mean =
            |i: usize| results.iter().map(|r| r.errors[i]).sum::<f64>() / results.len() as f64;
        let median = |i: usize| {
            let mut v: Vec<f64> = results.iter().map(|r| r.errors[i]).collect();
            v.sort_by(f64::total_cmp);
            v[v.len() / 2]
        };
        let ratios = |a: usize, b: usize| {
            let mut v: Vec<f64> = results.iter().map(|r| r.errors[a] / r.errors[b]).collect();
            v.sort_by(f64::total_cmp);
            v[v.len() / 2]
        };
        let resid = results
            .iter()
            .map(|r| r.transfer_residual)
            .fold(0.0, f64::max);
        let ns = results.iter().map(|r| r.n_s as f64).sum::<f64>() / results.len() as f64;
        println!(
            "{:>5} | {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e} | {:>8.2} {:>8.2} | {:>8.2} {:>8.2} | {:>8.0} {:>9.1e}",
            n,
            mean(0),
            mean(1),
            mean(2),
            mean(3),
            mean(0) / mean(1),
            mean(0) / mean(2),
            ratios(0, 1),
            ratios(0, 2),
            ns,
            resid
        );
        println!(
            "      | medians    {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e}",
            median(0),
            median(1),
            median(2),
            median(3)
        );
    }
}
