//! Flow fidelity against the paper's minimizers (Sec. 9.1, Fig. 1, 14, 15, 20); `NOTES-FIDELITY.md`.
//!
//! ```text
//! fidelity gen DIR
//!     Start meshes as OBJ: the symmetric torus of revolution (64 x 32, R = 1.0555, r = 0.4222,
//!     the `apps control_torus` start), a coarser one (48 x 24), genus 2 and 3 (`apps genus02/03`),
//!     linked and unlinked handcuffs (`apps handcuffs`).
//! fidelity flow OBJ [--steps N] [--threads T] [--exact | --theta TH] [--dense | --hier]
//!                   [--no-remesh] [--no-ccd] [--plane W] [--perturb EPS] [--every K]
//!                   [--converge REL --window W] [--swept] [--grid NM NN] [--exact-every K]
//!                   [--csv PATH] [--obj PATH] [--frames DIR --frame-every K]
//!                   [--volume-schedule N] [--volume-penalty W] [--boundary [--blen-scale S]] [--tau-start T]
//!     (--tau-start: a fixed first trial step, `StartRule::Absolute`, instead of doubling the last one)
//! fidelity symmetrize OBJ NM NN [--out PATH] [--threads T]
//!     A structured torus (vertex (i, j) = i NN + j, no remeshing) averaged over the rotations
//!     i -> i + 1 about the axis of least moment (rho, z and the angle offset per ring j), rescaled
//!     to the input's area; prints the exact energy along X_sym + t (X - X_sym).
//!     Area + barycenter constraints (as the applications), StopRule::Off. Prints every K steps:
//!     F, the exact energy (every --exact-every steps), area, inertia eigenvalues and the
//!     revolution indicators. --converge stops when the mean of F over the last W steps differs
//!     from the mean over the W steps before by less than REL (relative). --swept (needs
//!     --no-remesh) runs the library CCD with a tiny gap (1e-4 l0, close factor 0) on the linear
//!     motion from each step's start to its end: a contact is a (near-)crossing.
//! fidelity measure OBJ...
//!     Exact energy, area, volume, inertia indicators, revolution indicators, self-intersection.
//! fidelity scan --area A --faces F [--threads T] Q...
//!     Exact energy of tori of revolution with area A, R/r = Q, about F faces (n_major/n_minor ~ Q).
//! ```
#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::manual_is_multiple_of,
    clippy::collapsible_if
)]

#[path = "apps/apps.rs"]
mod apps;

use std::io::Write as _;
use std::sync::Arc;
use std::time::Instant;

use ddg_repulsive_surfaces::ccd::{self, CcdOptions};
use ddg_repulsive_surfaces::generate::{Mesh, primitives, shapes};
use ddg_repulsive_surfaces::implicit::Plane;
use ddg_repulsive_surfaces::step::{CSV_HEADER, csv_line};
use ddg_repulsive_surfaces::{
    Backend, BhOptions, Constraint, EnergyEval, Exponents, FlowOptions, HierarchicalOptions,
    ImplicitSurface, Penalty, Quadrature, Schedule, Session, StartRule, StopRule, Surface, energy, read_obj, write_obj,
};

type P = [f64; 3];

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: P, b: P) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: P, b: P) -> P {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn norm(a: P) -> f64 {
    dot(a, a).sqrt()
}

/// Jacobi eigen-decomposition of a symmetric 3x3 matrix: ascending eigenvalues, eigenvectors as
/// columns `v[k][i]` = component k of vector i.
fn eigsym(mut a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..100 {
        let (mut p, mut q, mut m) = (0, 1, 0.0);
        for i in 0..3 {
            for j in i + 1..3 {
                if a[i][j].abs() > m {
                    m = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if m < 1e-15 * (a[0][0].abs() + a[1][1].abs() + a[2][2].abs()).max(1e-300) {
            break;
        }
        let th = 0.5 * (2.0 * a[p][q]).atan2(a[q][q] - a[p][p]);
        let (c, s) = (th.cos(), th.sin());
        for k in 0..3 {
            let (akp, akq) = (a[k][p], a[k][q]);
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
        }
        for k in 0..3 {
            let (apk, aqk) = (a[p][k], a[q][k]);
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
        }
        for k in 0..3 {
            let (vkp, vkq) = (v[k][p], v[k][q]);
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
        }
    }
    let mut idx = [0usize, 1, 2];
    idx.sort_by(|&i, &j| a[i][i].partial_cmp(&a[j][j]).unwrap());
    let e = idx.map(|i| a[i][i]);
    let mut w = [[0.0; 3]; 3];
    for k in 0..3 {
        for (c, &i) in idx.iter().enumerate() {
            w[k][c] = v[k][i];
        }
    }
    (e, w)
}

/// Indicators of a closed surface.
struct Ind {
    area: f64,
    volume: f64,
    eig: [f64; 3],
    /// `(e2 - e1) / e2` for the two largest moments (0 for a torus of revolution about the
    /// least-moment axis).
    aniso_top: f64,
    /// `(e1 - e0) / e1` (0 for a body of revolution about the largest-moment axis, prolate).
    aniso_bot: f64,
    /// Angular Fourier amplitudes `|Σ a ρ e^{ikφ}| / Σ a ρ`, k = 1..4, about the least-moment
    /// axis through the area centroid.
    rho_modes: [f64; 4],
    /// Same for the height `z'` along the axis: `|Σ a z' e^{ikφ}| / Σ a |z'|`.
    z_modes: [f64; 4],
}

fn indicators(faces: &[[u32; 3]], x: &[P]) -> Ind {
    let (area, volume) = apps::area_volume(faces, x);
    let mut c = [0.0; 3];
    let mut cen = Vec::with_capacity(faces.len());
    for f in faces {
        let [p, q, r] = f.map(|i| x[i as usize]);
        let a = 0.5 * norm(cross(sub(q, p), sub(r, p)));
        let m = [0, 1, 2].map(|k| (p[k] + q[k] + r[k]) / 3.0);
        for k in 0..3 {
            c[k] += a * m[k];
        }
        cen.push((a, m));
    }
    for k in 0..3 {
        c[k] /= area;
    }
    let mut t = [[0.0; 3]; 3];
    for (a, m) in &cen {
        let d = sub(*m, c);
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] += a * d[i] * d[j] / area;
            }
        }
    }
    let (eig, v) = eigsym(t);
    let axis = [v[0][0], v[1][0], v[2][0]];
    let e1 = [v[0][1], v[1][1], v[2][1]];
    let e2 = [v[0][2], v[1][2], v[2][2]];
    let (mut rm, mut zm) = ([[0.0f64; 2]; 4], [[0.0f64; 2]; 4]);
    let (mut sr, mut sz) = (0.0, 0.0);
    for (a, m) in &cen {
        let d = sub(*m, c);
        let (u, w, z) = (dot(d, e1), dot(d, e2), dot(d, axis));
        let rho = u.hypot(w);
        let phi = w.atan2(u);
        for k in 0..4 {
            let kk = (k + 1) as f64;
            rm[k][0] += a * rho * (kk * phi).cos();
            rm[k][1] += a * rho * (kk * phi).sin();
            zm[k][0] += a * z * (kk * phi).cos();
            zm[k][1] += a * z * (kk * phi).sin();
        }
        sr += a * rho;
        sz += a * z.abs();
    }
    Ind {
        area,
        volume,
        eig,
        aniso_top: (eig[2] - eig[1]) / eig[2],
        aniso_bot: (eig[1] - eig[0]) / eig[1],
        rho_modes: rm.map(|m| m[0].hypot(m[1]) / sr),
        z_modes: zm.map(|m| m[0].hypot(m[1]) / sz.max(1e-300)),
    }
}

impl Ind {
    fn line(&self) -> String {
        format!(
            "area={:.6} volume={:.5} eig=[{:.5},{:.5},{:.5}] aniso_top={:.3e} aniso_bot={:.3e} rho_modes=[{:.1e},{:.1e},{:.1e},{:.1e}] z_modes=[{:.1e},{:.1e},{:.1e},{:.1e}]",
            self.area,
            self.volume,
            self.eig[0],
            self.eig[1],
            self.eig[2],
            self.aniso_top,
            self.aniso_bot,
            self.rho_modes[0],
            self.rho_modes[1],
            self.rho_modes[2],
            self.rho_modes[3],
            self.z_modes[0],
            self.z_modes[1],
            self.z_modes[2],
            self.z_modes[3]
        )
    }
}

/// Revolution defect of a structured torus (vertex `(i, j)` = `i nn + j`): the largest over the
/// minor index `j` of the spread (max - min) over `i` of `ρ` and `z`, about the z axis through the
/// vertex mean, relative to the minor radius estimate.
fn grid_defect(x: &[P], nm: usize, nn: usize) -> f64 {
    if x.len() != nm * nn {
        return f64::NAN;
    }
    let mut c = [0.0; 3];
    for p in x {
        for k in 0..3 {
            c[k] += p[k] / x.len() as f64;
        }
    }
    let mut worst: f64 = 0.0;
    let (mut zmin, mut zmax) = (f64::MAX, f64::MIN);
    for j in 0..nn {
        let (mut rlo, mut rhi, mut zlo, mut zhi) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for i in 0..nm {
            let p = sub(x[i * nn + j], c);
            let r = p[0].hypot(p[1]);
            rlo = rlo.min(r);
            rhi = rhi.max(r);
            zlo = zlo.min(p[2]);
            zhi = zhi.max(p[2]);
            zmin = zmin.min(p[2]);
            zmax = zmax.max(p[2]);
        }
        worst = worst.max(rhi - rlo).max(zhi - zlo);
    }
    worst / (0.5 * (zmax - zmin))
}

fn exact_energy(s: &Surface, x: &[P], threads: usize) -> f64 {
    energy(s, x, Exponents::default(), &EnergyEval::Exact, &[], Some(threads)).map(|r| r.energy).unwrap_or(f64::NAN)
}

fn plane(w: f64) -> Penalty {
    let s: Arc<dyn ImplicitSurface> = Arc::new(Plane::new([0.0; 3], [0.0, 0.0, 1.0]).unwrap());
    Penalty::ImplicitAttractor { weight: w, surface: s, power: 2.0, quadrature: Quadrature::DualArea, vertices: None }
}

fn unlinked_handcuffs(radius: f64, faces: usize) -> Mesh {
    let n = 160;
    let ring = |cx: f64| -> Vec<P> {
        (0..n)
            .map(|i| {
                let t = std::f64::consts::TAU * i as f64 / n as f64;
                [cx + 1.5 * t.cos(), 1.5 * t.sin(), 0.0]
            })
            .collect()
    };
    let bar: Vec<P> = vec![[-0.5, 0.0, 0.0], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
    shapes::graph_tube(&[(ring(-2.0), true), (ring(2.0), true), (bar, false)], radius, faces).unwrap()
}

fn cmd_gen(dir: &str) {
    std::fs::create_dir_all(dir).unwrap();
    let save = |name: &str, m: &Mesh| {
        let s = m.surface().unwrap();
        write_obj(format!("{dir}/{name}.obj"), &s, &m.positions).unwrap();
        println!("{name}: {} faces, chi {}", s.n_faces(), s.euler_characteristic());
    };
    save("torus64", &primitives::torus(64, 32, 1.0555, 0.4222).unwrap());
    save("torus48", &primitives::torus(48, 24, 1.0555, 0.4222).unwrap());
    save("torus32", &primitives::torus(32, 16, 1.0555, 0.4222).unwrap());
    for g in [1u32, 2, 3] {
        save(&format!("genus{g:02}"), &shapes::genus_surface(g, 3000 + 600 * g as usize).unwrap());
    }
    let linked = shapes::handcuffs(1, 9000).unwrap();
    let unl = unlinked_handcuffs(linked.radius, 9000);
    let unl = unl.clone().transformed((linked.mesh.area() / unl.area()).sqrt(), [0.0; 3]);
    save("handcuffs_linked", &linked.mesh);
    save("handcuffs_unlinked", &unl);
    save("punctured_torus", &primitives::punctured_torus(64, 32, 2.0, 0.8, 0.55).unwrap());
    save("punctured_torus_big", &primitives::punctured_torus(64, 32, 2.0, 0.8, 1.3).unwrap());
}

fn cmd_eval(paths: &[String], threads: usize) {
    for p in paths {
        let (s, x) = read_obj(p).unwrap();
        let e = exact_energy(&s, &x, threads);
        let ind = indicators(s.faces(), &x);
        let pairs = ccd::intersecting_pairs(&s, &x).unwrap();
        let inter = format!("{} pairs={} {:?}", !pairs.is_empty(), pairs.len(), &pairs[..pairs.len().min(4)]);
        println!(
            "MEASURE {p} faces={} genus={} exact_energy={e:.6} {} intersects={inter}",
            s.n_faces(),
            1 - s.euler_characteristic() / 2,
            ind.line()
        );
    }
}

fn cmd_scan(area: f64, faces: usize, threads: usize, qs: &[f64]) {
    for &q in qs {
        // A = 4 pi^2 R r, R = q r
        let r = (area / (4.0 * std::f64::consts::PI.powi(2) * q)).sqrt();
        let big = q * r;
        let nn = ((faces as f64 / (2.0 * q)).sqrt().round() as usize).max(3);
        let nm = ((q * nn as f64).round() as usize).max(3);
        let m = primitives::torus(nm, nn, big, r).unwrap();
        let s = m.surface().unwrap();
        // rescale to the exact area target (the discrete area is below the smooth one)
        let (a, _) = apps::area_volume(s.faces(), &m.positions);
        let k = (area / a).sqrt();
        let x: Vec<P> = m.positions.iter().map(|p| p.map(|c| c * k)).collect();
        let e = exact_energy(&s, &x, threads);
        println!("SCAN q={q:.3} R={:.4} r={:.4} grid={nm}x{nn} faces={} exact_energy={e:.6}", big * k, r * k, s.n_faces());
    }
}

fn cmd_symmetrize(path: &str, nm: usize, nn: usize, out: Option<String>, threads: usize) {
    let (s, x) = read_obj(path).unwrap();
    assert_eq!(x.len(), nm * nn, "not a {nm} x {nn} grid");
    let ind = indicators(s.faces(), &x);
    // frame: area centroid and principal axes; axis = least moment
    let mut c = [0.0; 3];
    let mut t = [[0.0; 3]; 3];
    let (mut at, mut cen) = (0.0, Vec::new());
    for f in s.faces() {
        let [p, q, r] = f.map(|i| x[i as usize]);
        let a = 0.5 * norm(cross(sub(q, p), sub(r, p)));
        let m = [0, 1, 2].map(|k| (p[k] + q[k] + r[k]) / 3.0);
        for k in 0..3 {
            c[k] += a * m[k];
        }
        at += a;
        cen.push((a, m));
    }
    for k in 0..3 {
        c[k] /= at;
    }
    for (a, m) in &cen {
        let d = sub(*m, c);
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] += a * d[i] * d[j] / at;
            }
        }
    }
    let (_, v) = eigsym(t);
    let ax = [v[0][0], v[1][0], v[2][0]];
    let e1 = [v[0][1], v[1][1], v[2][1]];
    let e2 = cross(ax, e1);
    let loc: Vec<P> = x
        .iter()
        .map(|p| {
            let d = sub(*p, c);
            [dot(d, e1), dot(d, e2), dot(d, ax)]
        })
        .collect();
    // per ring j: mean rho, z and angle offset relative to 2 pi i / nm
    let tau = std::f64::consts::TAU;
    let mut sym = vec![[0.0; 3]; x.len()];
    for j in 0..nn {
        let (mut rho, mut z, mut off) = (0.0, 0.0, [0.0, 0.0]);
        for i in 0..nm {
            let p = loc[i * nn + j];
            rho += p[0].hypot(p[1]) / nm as f64;
            z += p[2] / nm as f64;
            let ang = p[1].atan2(p[0]) - tau * i as f64 / nm as f64;
            off[0] += ang.cos();
            off[1] += ang.sin();
        }
        let o = off[1].atan2(off[0]);
        for i in 0..nm {
            let ang = o + tau * i as f64 / nm as f64;
            let l = [rho * ang.cos(), rho * ang.sin(), z];
            sym[i * nn + j] = [0, 1, 2].map(|k| c[k] + l[0] * e1[k] + l[1] * e2[k] + l[2] * ax[k]);
        }
    }
    let (a_sym, _) = apps::area_volume(s.faces(), &sym);
    let k = (ind.area / a_sym).sqrt();
    let sym: Vec<P> = sym.iter().map(|p| [0, 1, 2].map(|q| c[q] + k * (p[q] - c[q]))).collect();
    println!("SYM input {} | symmetrized {}", ind.line(), indicators(s.faces(), &sym).line());
    for tt in [-0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5] {
        let y: Vec<P> = sym.iter().zip(&x).map(|(p, q)| [0, 1, 2].map(|k| p[k] + tt * (q[k] - p[k]))).collect();
        let (ay, _) = apps::area_volume(s.faces(), &y);
        let sc = (ind.area / ay).sqrt();
        let y: Vec<P> = y.iter().map(|p| [0, 1, 2].map(|q| c[q] + sc * (p[q] - c[q]))).collect();
        println!("SEGMENT t={tt:+.2} exact_energy={:.8} (area rescaled by {sc:.6})", exact_energy(&s, &y, threads));
    }
    if let Some(p) = out {
        write_obj(p, &s, &sym).unwrap();
    }
}

struct FlowArgs {
    steps: usize,
    threads: usize,
    exact: bool,
    theta: f64,
    backend: String,
    remesh: bool,
    ccd: bool,
    plane: Option<f64>,
    perturb: f64,
    every: usize,
    exact_every: usize,
    converge: Option<f64>,
    window: usize,
    swept: bool,
    grid: Option<(usize, usize)>,
    csv: Option<String>,
    obj: Option<String>,
    frames: Option<String>,
    frame_every: usize,
    /// Fig. 20: a hard volume constraint moved to 0 over this many steps.
    volume_schedule: Option<usize>,
    /// Fig. 20: `w V^2` (the `apps eversion` setup uses 1e4).
    volume_penalty: Option<f64>,
    /// Fig. 20: boundary length (weight 30, initial length) and curvature (weight 1) penalties.
    boundary: bool,
    /// Fig. 20: the boundary-length target as a multiple of the initial length.
    blen_scale: f64,
    /// A fixed start `τ₀` (`StartRule::Absolute`) instead of the doubling default.
    tau_start: Option<f64>,
}

fn cmd_flow(path: &str, a: &FlowArgs) {
    let (s, mut x) = read_obj(path).unwrap();
    if a.perturb > 0.0 {
        let (a0, _) = apps::area_volume(s.faces(), &x);
        // deterministic pseudo-random displacement, a.perturb times the mean edge length
        let l: f64 = s.edges().iter().map(|e| norm(sub(x[e[0] as usize], x[e[1] as usize]))).sum::<f64>()
            / s.n_edges() as f64;
        let mut st: u64 = 0x9E3779B97F4A7C15;
        let mut rnd = || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            (st >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
        };
        for p in x.iter_mut() {
            for k in 0..3 {
                p[k] += a.perturb * l * rnd();
            }
        }
        // keep the area (the energy scales like area^-1 at (6, 12)): rescale about the mean
        let (a1, _) = apps::area_volume(s.faces(), &x);
        let k = (a0 / a1).sqrt();
        let c = x.iter().fold([0.0; 3], |c, p| [0, 1, 2].map(|i| c[i] + p[i] / x.len() as f64));
        for p in x.iter_mut() {
            *p = [0, 1, 2].map(|i| c[i] + k * (p[i] - c[i]));
        }
    }
    let mut o = FlowOptions {
        max_iterations: a.steps,
        stop: StopRule::Off,
        threads: Some(a.threads),
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        ..FlowOptions::default()
    };
    if a.exact {
        o.energy_eval = EnergyEval::Exact;
    } else {
        o.energy_eval = EnergyEval::BarnesHut(BhOptions { theta: a.theta, ..BhOptions::default() });
    }
    match a.backend.as_str() {
        "dense" => o.backend = Backend::Dense,
        "hier" => o.backend = Backend::Hierarchical(HierarchicalOptions::default()),
        _ => {}
    }
    if !a.remesh {
        o.remesh = None;
    }
    if !a.ccd {
        o.line_search.ccd = None;
    }
    if let Some(t) = a.tau_start {
        o.line_search.start = StartRule::Absolute { tau: t };
    }
    if let Some(n) = a.volume_schedule {
        o.constraints.push(Constraint::Volume { component: None, target: None });
        o.schedules.push(Schedule { constraint: o.constraints.len() - 1, factor: 0.0, offset: 0.0, steps: n });
    }
    if let Some(w) = a.volume_penalty {
        o.penalties.push(Penalty::TotalVolume { weight: w, power: 2.0 });
    }
    if a.boundary {
        let blen: f64 = s
            .edges()
            .iter()
            .enumerate()
            .filter(|(e, _)| s.is_boundary_edge(*e as u32))
            .map(|(_, e)| norm(sub(x[e[0] as usize], x[e[1] as usize])))
            .sum();
        o.penalties.push(Penalty::BoundaryLength { weight: 30.0, target: a.blen_scale * blen });
        o.penalties.push(Penalty::BoundaryCurvature { weight: 1.0, dual_length: true });
    }
    if let Some(w) = a.plane {
        o.penalties.push(plane(w));
    }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(a.threads).build().unwrap();
    let mut csv = a.csv.as_ref().map(|p| {
        let mut f = std::io::BufWriter::new(std::fs::File::create(p).unwrap());
        writeln!(f, "{CSV_HEADER},exact_energy,swept_contact_t").unwrap();
        f
    });
    if let Some(d) = &a.frames {
        std::fs::create_dir_all(d).unwrap();
    }
    let t0 = Instant::now();
    let e_ex0 = exact_energy(&s, &x, a.threads);
    let ind0 = indicators(s.faces(), &x);
    let gd0 = a.grid.map_or(f64::NAN, |(m, n)| grid_defect(&x, m, n));
    println!("STEP 0 F=nan exact={e_ex0:.6} grid_defect={gd0:.3e} {}", ind0.line());
    let mut sess = Session::new(s, x, o).unwrap();
    let l0 = sess.mean_edge_length0();
    let mut hist: Vec<f64> = Vec::new();
    let mut swept_hits = 0usize;
    let mut status = "MaxIterations".to_string();
    let mut capped = 0usize;
    let mut exempt = 0usize;
    while sess.iteration() < a.steps {
        let p0 = sess.positions().to_vec();
        let nf0 = sess.surface().n_faces();
        let r = match sess.step() {
            Ok(r) => r,
            Err(e) => {
                status = format!("Error({e})");
                break;
            }
        };
        let it = sess.iteration();
        if let Some(c) = &r.ccd {
            if let Some(ct) = &c.contact {
                capped += 1;
                // the limiting pair, its distance at the step start in l0, boundary flags
                let (s, x) = (sess.surface(), &p0);
                let bv = |v: u32| s.is_boundary_vertex(v);
                let desc = match ct.primitive {
                    ccd::Primitive::VertexTriangle { vertex, triangle } if (vertex as usize) < x.len() && (triangle as usize) < s.n_faces() => {
                        let t = s.faces()[triangle as usize];
                        let d = ddg_repulsive_curves::ccd::point_triangle_distance(x[vertex as usize], x[t[0] as usize], x[t[1] as usize], x[t[2] as usize]);
                        format!("VT v{vertex}(bd {}) f{triangle}(bd {}) d0={:.3e}l0", bv(vertex), t.iter().any(|&v| bv(v)), d / l0)
                    }
                    ccd::Primitive::EdgeEdge { a, b } if (a as usize) < s.n_edges() && (b as usize) < s.n_edges() => {
                        let (ea, eb) = (s.edges()[a as usize], s.edges()[b as usize]);
                        let d = ddg_repulsive_curves::ccd::segment_segment_distance(x[ea[0] as usize], x[ea[1] as usize], x[eb[0] as usize], x[eb[1] as usize]);
                        format!("EE {ea:?}(bd {}) {eb:?}(bd {}) d0={:.3e}l0", s.is_boundary_edge(a), s.is_boundary_edge(b), d / l0)
                    }
                    other => format!("{other:?}"),
                };
                if capped <= 200 {
                    println!("CAP step={} tau_max={:.3e} {desc}", sess.iteration(), c.tau_max);
                }
            }
            exempt += c.exempt_pairs as usize;
        }
        let mut swept_t = f64::NAN;
        if a.swept && sess.surface().n_faces() == nf0 && sess.positions().len() == p0.len() {
            let v: Vec<P> = p0.iter().zip(sess.positions()).map(|(p, q)| sub(*q, *p)).collect();
            let opts = CcdOptions { gap: 1e-4, reference_length: Some(l0), close_factor: 0.0, ..CcdOptions::default() };
            if let Ok(rep) = ccd::first_contact(&pool, sess.surface(), &p0, &v, 1.0, &[], &opts) {
                if let Some(cn) = rep.contact {
                    swept_hits += 1;
                    swept_t = rep.tau_max;
                    println!("SWEPT step={it} t={:.4} contact={cn:?} max_move={:.3}l0", rep.tau_max,
                        v.iter().map(|d| norm(*d)).fold(0.0, f64::max) / l0);
                }
            }
        }
        hist.push(r.energy_after);
        let want_exact = a.exact_every > 0 && it % a.exact_every == 0;
        let e_ex = if want_exact { exact_energy(sess.surface(), sess.positions(), a.threads) } else { f64::NAN };
        if let Some(f) = csv.as_mut() {
            writeln!(f, "{},{e_ex:.10e},{swept_t:.6e}", csv_line(&r)).unwrap();
            if it % 20 == 0 {
                f.flush().unwrap();
            }
        }
        if let Some(d) = &a.frames {
            if a.frame_every > 0 && it % a.frame_every == 0 {
                write_obj(format!("{d}/frame_{it:05}.obj"), sess.surface(), sess.positions()).unwrap();
            }
        }
        if a.every > 0 && it % a.every == 0 {
            let ind = indicators(sess.surface().faces(), sess.positions());
            let gd = a.grid.map_or(f64::NAN, |(m, n)| grid_defect(sess.positions(), m, n));
            println!(
                "STEP {it} F={:.8} exact={e_ex:.6} tau={:.3e} faces={} capped={capped} exempt={exempt} swept_hits={swept_hits} grid_defect={gd:.3e} {} t={:.0}s",
                r.energy_after,
                r.tau,
                sess.surface().n_faces(),
                ind.line(),
                t0.elapsed().as_secs_f64()
            );
        }
        if sess.status().is_terminal() {
            status = format!("{:?}", sess.status());
            break;
        }
        if let Some(rel) = a.converge {
            let w = a.window;
            if hist.len() >= 2 * w && hist.len() % 10 == 0 {
                let n = hist.len();
                let m1: f64 = hist[n - w..].iter().sum::<f64>() / w as f64;
                let m0: f64 = hist[n - 2 * w..n - w].iter().sum::<f64>() / w as f64;
                if (m0 - m1).abs() <= rel * m1.abs() {
                    status = format!("Converged(mean change {:.2e} over {w})", (m0 - m1) / m1);
                    break;
                }
            }
        }
    }
    if let Some(f) = csv.as_mut() {
        f.flush().unwrap();
    }
    let e_ex = exact_energy(sess.surface(), sess.positions(), a.threads);
    let ind = indicators(sess.surface().faces(), sess.positions());
    let gd = a.grid.map_or(f64::NAN, |(m, n)| grid_defect(sess.positions(), m, n));
    let inter = ccd::intersects(sess.surface(), sess.positions()).unwrap();
    if let Some(p) = &a.obj {
        write_obj(p, sess.surface(), sess.positions()).unwrap();
    }
    println!(
        "FINAL steps={} status={status} F={:.8} exact={e_ex:.6} faces={} capped={capped} exempt={exempt} swept_hits={swept_hits} grid_defect={gd:.3e} intersects={inter} {} seconds={:.0}",
        sess.iteration(),
        hist.last().copied().unwrap_or(f64::NAN),
        sess.surface().n_faces(),
        ind.line(),
        t0.elapsed().as_secs_f64()
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| args.iter().position(|a| a == k).map(|i| args[i + 1].clone());
    let has = |k: &str| args.iter().any(|a| a == k);
    let threads: usize = get("--threads").map_or(4, |v| v.parse().unwrap());
    match args.first().map(String::as_str) {
        Some("gen") => cmd_gen(&args[1]),
        Some("symmetrize") => {
            let nm: usize = args[2].parse().unwrap();
            let nn: usize = args[3].parse().unwrap();
            cmd_symmetrize(&args[1], nm, nn, get("--out"), threads)
        }
        Some("measure") => {
            let paths: Vec<String> = args[1..].iter().take_while(|a| !a.starts_with("--")).cloned().collect();
            cmd_eval(&paths, threads)
        }
        Some("scan") => {
            let area: f64 = get("--area").unwrap().parse().unwrap();
            let faces: usize = get("--faces").unwrap().parse().unwrap();
            let mut qs = Vec::new();
            let mut i = 1;
            while i < args.len() {
                if args[i].starts_with("--") {
                    i += 2;
                    continue;
                }
                qs.push(args[i].parse().unwrap());
                i += 1;
            }
            cmd_scan(area, faces, threads, &qs)
        }
        Some("flow") => {
            let grid = args.iter().position(|a| a == "--grid").map(|i| (args[i + 1].parse().unwrap(), args[i + 2].parse().unwrap()));
            let fa = FlowArgs {
                steps: get("--steps").map_or(800, |v| v.parse().unwrap()),
                threads,
                exact: has("--exact"),
                theta: get("--theta").map_or(0.5, |v| v.parse().unwrap()),
                backend: if has("--dense") { "dense".into() } else if has("--hier") { "hier".into() } else { "auto".into() },
                remesh: !has("--no-remesh"),
                ccd: !has("--no-ccd"),
                plane: get("--plane").map(|v| v.parse().unwrap()),
                perturb: get("--perturb").map_or(0.0, |v| v.parse().unwrap()),
                every: get("--every").map_or(50, |v| v.parse().unwrap()),
                exact_every: get("--exact-every").map_or(0, |v| v.parse().unwrap()),
                converge: get("--converge").map(|v| v.parse().unwrap()),
                window: get("--window").map_or(100, |v| v.parse().unwrap()),
                swept: has("--swept"),
                grid,
                csv: get("--csv"),
                obj: get("--obj"),
                frames: get("--frames"),
                frame_every: get("--frame-every").map_or(0, |v| v.parse().unwrap()),
                volume_schedule: get("--volume-schedule").map(|v| v.parse().unwrap()),
                volume_penalty: get("--volume-penalty").map(|v| v.parse().unwrap()),
                boundary: has("--boundary"),
                blen_scale: get("--blen-scale").map_or(1.0, |v| v.parse().unwrap()),
                tau_start: get("--tau-start").map(|v| v.parse().unwrap()),
            };
            cmd_flow(&args[1], &fa)
        }
        _ => eprintln!("usage: see the header of examples/fidelity.rs"),
    }
}
