//! WP16a: consistency study of the paper's Sec. 8.1 / Fig. 10 (SPEC section 14.1).
//!
//! Modes (see `NOTES-WP16A.md`):
//! * `check`                      embeddedness / Jacobian sanity of the test surfaces.
//! * `smooth SURF A B XRULE DRULE THREADS`   the smooth energy by cubature.
//! * `discrete SURF FACES THREADS [THETAS]`  discrete energies, optionally Barnes-Hut errors.
//!
//! `SURF` is `pert` (perturbed torus) or `torus` (plain torus of revolution).
//! Everything prints `KEY k=v ...` lines on stdout.

use std::f64::consts::{PI, TAU};
use std::time::Instant;

use ddg_repulsive_surfaces::energy::barnes_hut::BhTree;
use ddg_repulsive_surfaces::energy::exact::{energy_value, make_pool};
use ddg_repulsive_surfaces::generate::Rng;
use ddg_repulsive_surfaces::generate::check::self_intersections;
use ddg_repulsive_surfaces::generate::perturbed::PerturbedTorus;
use ddg_repulsive_surfaces::geometry::{FaceGeometry, mean_edge_length};
use ddg_repulsive_surfaces::remesh::{RemeshOptions, remesh};
use ddg_repulsive_surfaces::{BhOptions, Exponents, Surface};
use rayon::prelude::*;

type V3 = [f64; 3];

const MAJOR: f64 = 1.0;
const MINOR: f64 = 1.0 / 3.0;
const DEGREE: u32 = 2;
const AMPLITUDE: f64 = 0.2;
const SEED: u64 = 16;

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// The surface with analytic first derivatives: `f = f0 + Phi(f0)`; the terms of `Phi` are
/// rebuilt exactly as `PerturbedTorus::new` does (checked against `PerturbedTorus::point`).
struct Surf {
    terms: Vec<(V3, V3, V3, f64)>,
}

impl Surf {
    fn new(amplitude: f64) -> Surf {
        let mut rng = Rng::new(SEED);
        let d = DEGREE as i32;
        let mut terms = Vec::new();
        for a in -d..=d {
            for b in -d..=d {
                for c in -d..=d {
                    let first = [a, b, c].into_iter().find(|&x| x != 0);
                    if first.is_none_or(|x| x < 0) {
                        continue;
                    }
                    let k = [a as f64, b as f64, c as f64];
                    let w = amplitude / (1.0 + dot(k, k));
                    let av = [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)];
                    let bv = [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)];
                    terms.push((k, av, bv, w));
                }
            }
        }
        Surf { terms }
    }

    /// `(f, f_u, f_v)`.
    fn eval(&self, u: f64, v: f64) -> (V3, V3, V3) {
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let ring = MAJOR + MINOR * cv;
        let p = [ring * cu, ring * su, MINOR * sv];
        let pu = [-ring * su, ring * cu, 0.0];
        let pv = [-MINOR * sv * cu, -MINOR * sv * su, MINOR * cv];
        let mut phi = [0.0; 3];
        // dphi[i][j] = d Phi_i / d x_j
        let mut dphi = [[0.0; 3]; 3];
        for (k, a, b, w) in &self.terms {
            let (s, c) = dot(*k, p).sin_cos();
            for i in 0..3 {
                phi[i] += w * (a[i] * c + b[i] * s);
                let g = w * (-a[i] * s + b[i] * c);
                for j in 0..3 {
                    dphi[i][j] += g * k[j];
                }
            }
        }
        let f = [p[0] + phi[0], p[1] + phi[1], p[2] + phi[2]];
        let push = |q: V3| -> V3 {
            [
                q[0] + dphi[0][0] * q[0] + dphi[0][1] * q[1] + dphi[0][2] * q[2],
                q[1] + dphi[1][0] * q[0] + dphi[1][1] * q[1] + dphi[1][2] * q[2],
                q[2] + dphi[2][0] * q[0] + dphi[2][1] * q[1] + dphi[2][2] * q[2],
            ]
        };
        (f, push(pu), push(pv))
    }
}

fn surf_for(name: &str) -> Surf {
    match name {
        "pert" => Surf::new(AMPLITUDE),
        "torus" => Surf::new(0.0),
        _ => panic!("surface must be pert or torus"),
    }
}

// ---------------------------------------------------------------- quadrature

/// Gauss-Legendre nodes and weights on `[-1, 1]`.
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n {
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, z);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            // p1 = P_n(z), p0 = P_{n-1}(z) (n >= 2)
            dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            let dz = p1 / dp;
            z -= dz;
            if dz.abs() < 1e-16 {
                break;
            }
        }
        x[i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
    (x, w)
}

/// GL nodes/weights on `[a, b]`.
fn gl(a: f64, b: f64, n: usize) -> Vec<(f64, f64)> {
    let (x, w) = gauss_legendre(n);
    x.iter().zip(&w).map(|(x, w)| (0.5 * (a + b) + 0.5 * (b - a) * x, 0.5 * (b - a) * w)).collect()
}

/// `P` equal panels of `n` GL nodes on `[a, b]`.
fn gl_panels(a: f64, b: f64, panels: usize, n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for p in 0..panels {
        let lo = a + (b - a) * p as f64 / panels as f64;
        let hi = a + (b - a) * (p + 1) as f64 / panels as f64;
        out.extend(gl(lo, hi, n));
    }
    out
}

type Pts = Vec<(f64, f64, f64)>;

/// Rules for the base point `x = (u, v)`.
fn x_rule(spec: &str) -> Pts {
    let parts: Vec<&str> = spec.split(':').collect();
    match parts[0] {
        // M x M trapezoid on the periodic torus (spectrally accurate for smooth periodic data)
        "trap" => {
            let m: usize = parts[1].parse().unwrap();
            let h = TAU / m as f64;
            (0..m).flat_map(|i| (0..m).map(move |j| (i as f64 * h, j as f64 * h, h * h))).collect()
        }
        // u-symmetry of the torus of revolution: u = 0 only, weight 2 pi
        "traps" => {
            let m: usize = parts[1].parse().unwrap();
            let h = TAU / m as f64;
            (0..m).map(|j| (0.0, j as f64 * h, TAU * h)).collect()
        }
        // tensor Gauss-Legendre, P panels of n nodes per direction
        "gl" => {
            let p: usize = parts[1].parse().unwrap();
            let n: usize = parts[2].parse().unwrap();
            let g = gl_panels(0.0, TAU, p, n);
            g.iter().flat_map(|&(u, wu)| g.iter().map(move |&(v, wv)| (u, v, wu * wv))).collect()
        }
        _ => panic!("bad x rule {spec}"),
    }
}

/// Rules for the offset `delta = y - x` over `[-pi, pi]^2` (the diagonal is the origin).
fn d_rule(spec: &str) -> Pts {
    let parts: Vec<&str> = spec.split(':').collect();
    let signs = [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)];
    let mut out = Vec::new();
    match parts[0] {
        // Duffy: each quadrant is two triangles (s, t) = (r, r w) and (r w, r); the Jacobian r
        // cancels the direction dependence, the transformed integrand is smooth.
        "duffy" => {
            let pr: usize = parts[1].parse().unwrap();
            let nr: usize = parts[2].parse().unwrap();
            let nw: usize = parts[3].parse().unwrap();
            let rs = gl_panels(0.0, PI, pr, nr);
            let ws = gl_panels(0.0, 1.0, 2, nw);
            for &(ss, st) in &signs {
                for &(r, wr) in &rs {
                    for &(w, ww) in &ws {
                        out.push((ss * r, st * r * w, wr * ww * r));
                        out.push((ss * r * w, st * r, wr * ww * r));
                    }
                }
            }
        }
        // geometric subdivision of the cells touching the diagonal: per quadrant and level l,
        // three squares of the L-shape around the corner cell [0, L/2]^2, plus the final
        // corner cell; plain tensor GL (n x n) on every cell, no change of variables.
        "cells" => {
            let depth: usize = parts[1].parse().unwrap();
            let n: usize = parts[2].parse().unwrap();
            let mut cells: Vec<[f64; 4]> = Vec::new();
            for l in 0..depth {
                let big = PI / (1u64 << l) as f64;
                let h = 0.5 * big;
                cells.push([h, big, 0.0, h]);
                cells.push([0.0, h, h, big]);
                cells.push([h, big, h, big]);
            }
            let last = PI / (1u64 << depth) as f64;
            cells.push([0.0, last, 0.0, last]);
            for &(ss, st) in &signs {
                for c in &cells {
                    let gs = gl(c[0], c[1], n);
                    let gt = gl(c[2], c[3], n);
                    for &(s, ws) in &gs {
                        for &(t, wt) in &gt {
                            out.push((ss * s, st * t, ws * wt));
                        }
                    }
                }
            }
        }
        _ => panic!("bad delta rule {spec}"),
    }
    out
}

/// `sum_x w_x J(x) sum_delta w_d J(y) |N(x).(f(x)-f(y))|^alpha |f(x)-f(y)|^-beta`.
fn smooth_energy(pool: &rayon::ThreadPool, surf: &Surf, alpha: f64, beta: f64, xs: &Pts, ds: &Pts) -> f64 {
    let (ha, hb) = (0.5 * alpha, 0.5 * beta);
    let pw = |x: f64, e: f64| if e.fract() == 0.0 { x.powi(e as i32) } else { x.powf(e) };
    let rows: Vec<f64> = pool.install(|| {
        xs.par_iter()
            .map(|&(u, v, wx)| {
                let (fx, fu, fv) = surf.eval(u, v);
                let c = cross(fu, fv);
                let jx = dot(c, c).sqrt();
                let n = [c[0] / jx, c[1] / jx, c[2] / jx];
                let mut acc = 0.0;
                for &(s, t, wd) in ds {
                    let (fy, gu, gv) = surf.eval(u + s, v + t);
                    let cy = cross(gu, gv);
                    let jy = dot(cy, cy).sqrt();
                    let d = sub(fx, fy);
                    let nd = dot(n, d);
                    let k = pw(nd * nd, ha) / pw(dot(d, d), hb);
                    acc += wd * jy * k;
                }
                wx * jx * acc
            })
            .collect()
    });
    rows.iter().sum()
}

// ------------------------------------------------------------- discretization

/// Closest point on the surface to `p` by damped Gauss-Newton in parameter space.
fn project(surf: &Surf, p: V3, mut u: f64, mut v: f64) -> (f64, f64) {
    let cost = |u: f64, v: f64| {
        let r = sub(surf.eval(u, v).0, p);
        dot(r, r)
    };
    let mut c0 = cost(u, v);
    let mut lambda = 1e-9;
    for _ in 0..100 {
        let (f, fu, fv) = surf.eval(u, v);
        let r = sub(f, p);
        let (g0, g1) = (dot(fu, r), dot(fv, r));
        if g0.hypot(g1) < 1e-14 {
            break;
        }
        let (a, b, d) = (dot(fu, fu), dot(fu, fv), dot(fv, fv));
        let mut accepted = false;
        for _ in 0..30 {
            let (a2, d2) = (a * (1.0 + lambda) + 1e-18, d * (1.0 + lambda) + 1e-18);
            let det = a2 * d2 - b * b;
            let (du, dv) = (-(d2 * g0 - b * g1) / det, -(a2 * g1 - b * g0) / det);
            let c1 = cost(u + du, v + dv);
            if c1 <= c0 {
                u += du;
                v += dv;
                c0 = c1;
                lambda = (lambda * 0.2).max(1e-12);
                accepted = true;
                break;
            }
            lambda *= 4.0;
        }
        if !accepted {
            break;
        }
    }
    (u, v)
}

struct Disc {
    surface: Surface,
    positions: Vec<V3>,
    max_move: f64,
}

/// Squeezed grid, then `passes` x (remesh, transport parameters, project back).
fn build_discrete(pool: &rayon::ThreadPool, name: &str, faces: usize, passes: usize) -> Disc {
    let surf = surf_for(name);
    let pt = PerturbedTorus::new(MAJOR, MINOR, DEGREE, if name == "pert" { AMPLITUDE } else { 0.0 }, SEED).unwrap();
    // the surf rebuild must agree with the library's
    for &(u, v) in &[(0.3, 1.1), (4.0, 5.5), (2.2, 0.1)] {
        let a = surf.eval(u, v).0;
        let b = pt.point(u, v);
        assert!(sub(a, b).iter().all(|x| x.abs() < 1e-13), "surface mismatch");
    }
    let (su, sv) = (0.4, 0.4);
    let nv = ((faces as f64 / 6.0).sqrt().round() as usize).max(4);
    let nu = 3 * nv;
    let mesh = pt.sample(nu, nv, su, sv).unwrap();
    let warp = |x: f64, s: f64| TAU * (x - s * (TAU * x).sin() / TAU);
    let mut uv: Vec<[f64; 4]> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let u = warp(i as f64 / nu as f64, su);
            let v = warp(j as f64 / nv as f64, sv);
            uv.push([u.cos(), u.sin(), v.cos(), v.sin()]);
        }
    }
    let mut surface = mesh.surface().unwrap();
    let mut positions = mesh.positions.clone();
    let area = FaceGeometry::new(&surface, &positions).unwrap().total_area();
    // equilateral triangles of edge l0 have area sqrt(3)/4 l0^2
    let l0 = (area / (faces as f64 * 3f64.sqrt() / 4.0)).sqrt();
    let opts = RemeshOptions::default();
    let mut max_move: f64 = 0.0;
    for _ in 0..passes {
        let (out, _rep) = remesh(&surface, &positions, &[], l0, &opts, None).unwrap();
        let comp = |k: usize| -> Vec<f64> {
            let col: Vec<f64> = uv.iter().map(|q| q[k]).collect();
            out.vertex_map.transport(&col).unwrap()
        };
        let cols = [comp(0), comp(1), comp(2), comp(3)];
        let n = out.positions.len();
        let res: Vec<(V3, [f64; 4], f64)> = pool.install(|| {
            (0..n)
                .into_par_iter()
                .map(|i| {
                    let u0 = cols[1][i].atan2(cols[0][i]);
                    let v0 = cols[3][i].atan2(cols[2][i]);
                    let (u, v) = project(&surf, out.positions[i], u0, v0);
                    let f = surf.eval(u, v).0;
                    let d = sub(f, out.positions[i]);
                    (f, [u.cos(), u.sin(), v.cos(), v.sin()], dot(d, d).sqrt())
                })
                .collect()
        });
        positions = res.iter().map(|r| r.0).collect();
        uv = res.iter().map(|r| r.1).collect();
        max_move = max_move.max(res.iter().map(|r| r.2).fold(0.0, f64::max));
        surface = out.surface;
    }
    Disc { surface, positions, max_move }
}

fn edge_stats(s: &Surface, p: &[V3]) -> (f64, f64, f64) {
    let (mut mx, mut mn, mut sum) = (0.0f64, f64::INFINITY, 0.0);
    for e in s.edges() {
        let d = sub(p[e[0] as usize], p[e[1] as usize]);
        let l = dot(d, d).sqrt();
        mx = mx.max(l);
        mn = mn.min(l);
        sum += l;
    }
    (mx, mn, sum / s.n_edges() as f64)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("check") => {
            let surf = surf_for("pert");
            let (mut jmin, mut jmax) = (f64::INFINITY, 0.0f64);
            for i in 0..400 {
                for j in 0..400 {
                    let (_, fu, fv) = surf.eval(TAU * i as f64 / 400.0, TAU * j as f64 / 400.0);
                    let c = cross(fu, fv);
                    let j2 = dot(c, c).sqrt();
                    jmin = jmin.min(j2);
                    jmax = jmax.max(j2);
                }
            }
            let pt = PerturbedTorus::new(MAJOR, MINOR, DEGREE, AMPLITUDE, SEED).unwrap();
            let mesh = pt.sample(300, 120, 0.0, 0.0).unwrap();
            let si = self_intersections(&mesh, 5);
            println!(
                "CHECK jacobian_min={jmin:.4} jacobian_max={jmax:.4} self_intersections={} area={:.6}",
                si.len(),
                mesh.area()
            );
            // finite-difference check of the analytic derivatives
            let (u, v, h) = (1.3, 2.9, 1e-5);
            let (_, fu, fv) = surf.eval(u, v);
            let fdu = sub(surf.eval(u + h, v).0, surf.eval(u - h, v).0).map(|x| x / (2.0 * h));
            let fdv = sub(surf.eval(u, v + h).0, surf.eval(u, v - h).0).map(|x| x / (2.0 * h));
            let e = |a: V3, b: V3| dot(sub(a, b), sub(a, b)).sqrt();
            println!("CHECK derivative_fd_error u={:.2e} v={:.2e}", e(fu, fdu), e(fv, fdv));
        }
        Some("smooth") => {
            let (name, a, b) = (&args[1], args[2].parse::<f64>().unwrap(), args[3].parse::<f64>().unwrap());
            let (xr, dr) = (&args[4], &args[5]);
            let pool = make_pool(Some(args[6].parse().unwrap())).unwrap();
            let surf = surf_for(name);
            let (xs, ds) = (x_rule(xr), d_rule(dr));
            let t0 = Instant::now();
            let val = smooth_energy(&pool, &surf, a, b, &xs, &ds);
            println!(
                "SMOOTH surf={name} alpha={a} beta={b} x={xr} d={dr} nx={} nd={} value={val:.15e} secs={:.1}",
                xs.len(),
                ds.len(),
                t0.elapsed().as_secs_f64()
            );
        }
        Some("discrete") => {
            let (name, faces) = (&args[1], args[2].parse::<usize>().unwrap());
            let pool = make_pool(Some(args[3].parse().unwrap())).unwrap();
            let thetas: Vec<f64> =
                args.get(4).map(|s| s.split(',').map(|x| x.parse().unwrap()).collect()).unwrap_or_default();
            let t0 = Instant::now();
            let d = build_discrete(&pool, name, faces, 4);
            let (hmax, hmin, hmean) = edge_stats(&d.surface, &d.positions);
            let geom = FaceGeometry::new(&d.surface, &d.positions).unwrap();
            let l0 = mean_edge_length(&d.surface, &d.positions);
            let mut line = format!(
                "DISC surf={name} target_faces={faces} faces={} vertices={} h={hmax:.8e} edge_min={hmin:.4e} edge_mean={hmean:.4e} l0={l0:.4e} area={:.8} proj_max_move={:.2e} mesh_secs={:.1}",
                d.surface.n_faces(),
                d.surface.n_vertices(),
                geom.total_area(),
                d.max_move,
                t0.elapsed().as_secs_f64()
            );
            for (a, b) in [(6.0, 12.0), (4.0, 8.0)] {
                let e = Exponents::new(a, b).unwrap();
                let t1 = Instant::now();
                let val = energy_value(&pool, &geom, &e);
                line += &format!(" E_{a}_{b}={val:.15e} secs_{a}_{b}={:.1}", t1.elapsed().as_secs_f64());
            }
            println!("{line}");
            if !thetas.is_empty() {
                for (a, b) in [(6.0, 12.0), (4.0, 8.0)] {
                    let e = Exponents::new(a, b).unwrap();
                    let exact = energy_value(&pool, &geom, &e);
                    for &theta in &thetas {
                        let t1 = Instant::now();
                        let val = if theta == 0.0 {
                            exact
                        } else {
                            let tree = BhTree::build(&pool, &d.surface, &d.positions, BhOptions { theta, leaf_size: 8 })
                                .unwrap();
                            tree.energy(&pool, &e)
                        };
                        println!(
                            "THETA surf={name} faces={} alpha={a} beta={b} theta={theta} bh={val:.15e} exact={exact:.15e} rel_err={:.6e} secs={:.1}",
                            d.surface.n_faces(),
                            (val - exact).abs() / exact,
                            t1.elapsed().as_secs_f64()
                        );
                    }
                }
            }
        }
        _ => eprintln!("usage: consistency check | smooth SURF A B XRULE DRULE THREADS | discrete SURF FACES THREADS [THETAS]"),
    }
}
