//! Genus 11 (Fig. 15 of "Repulsive Surfaces", Yu, Brakensiek, Schumacher, Crane 2021): does the flow reach the
//! dodecahedral minimizer? `NOTES-DODECA.md`.
//!
//! ```text
//! dodeca gen ROW_OBJ OUT_OBJ [--faces F] [--radius R]
//!     A tube (radius R times the edge length, default 0.3) around the 30 edges of a regular dodecahedron
//!     (genus 30 - 20 + 1 = 11), about F faces (default 10136), scaled to the area of ROW_OBJ.
//! dodeca flow OBJ --out DIR [--threads T] [--max-steps S] [--window W] [--rel R] [--frame-every K]
//!                 [--perturb EPS --seed N] [--compat] [--no-ccd] [--armijo-ref] [--start-ref]
//!                 [--precond-ref] [--gmres-left] [--gmres-tol X] [--theta X] [--smooth-first]
//!                 [--bary global|percomp] [--newton1] [--no-remesh] [--no-area]
//!     Area + barycenter constraints (`--no-area`: barycenter only), defaults unless stated; the stop rule of `ishii genus` (mean F over the
//!     last W = 200 steps within R = 1e-5 relative of the W before; cap S = 8000). Writes DIR/steps.csv,
//!     DIR/frame_NNNNN.obj (every K = 100), DIR/final.obj; prints PROGRESS every 25 steps and SYM (moments,
//!     Q2 Q4 Q6) every 100.
//!     --perturb: every vertex moved by EPS mean edge lengths in a random direction (xorshift seeded by N),
//!     then rescaled to the input area.
//! dodeca ladder ROW_OBJ OUT_OBJ [--radius R] [--faces F]
//!     A row of 11 rings (tube of radius R = 0.22 around the ladder graph with 12 rungs, unit spacing),
//!     scaled to the area of ROW_OBJ.
//! dodeca chain ROW_OBJ OUT_OBJ [--hole H] [--outer D] [--waist W] [--rings N] [--faces F]
//!     The paper video's genus-11 start: N = 11 round tori (outer diameter D = 0.86, hole H = 0.3, in units of
//!     the ring spacing) in a row, fused by a smooth union with a neck W = 0.3 wide; about F = 12000 faces,
//!     scaled to the area of ROW_OBJ.
//! dodeca perturb IN_OBJ OUT_OBJ --perturb EPS --seed N
//!     The perturbed start of `flow --perturb EPS --seed N`, as OBJ.
//! dodeca measure OBJ... [--area A] [--threads T]
//!     Exact all-pairs energy, area, the energy at area A (E * area / A: E ~ area^-1 at (6, 12)), moments,
//!     Q2 Q4 Q6, genus, self-intersection.
//! ```
//! Symmetry indicators: eigenvalues of the area-weighted second-moment tensor about the area centroid, and
//! `Q_l = Σ_{f,g} a_f a_g P_l(u_f · u_g) / A²` (u the unit direction of a face centroid from the area centroid;
//! the rotation-invariant power of degree l of the area distribution on the sphere of directions). A shape
//! with icosahedral symmetry has `Q2 = Q4 = 0` (no invariant harmonics of degree 2 and 4) and `Q6 > 0`; a
//! cube has `Q2 = 0`, `Q4 > 0`; an elongated shape `Q2 > 0`.
#![allow(clippy::needless_range_loop)]

#[path = "apps/apps.rs"]
mod apps;

use std::io::Write as _;
use std::time::Instant;

use ddg_repulsive_surfaces::generate::shapes;
use ddg_repulsive_surfaces::metric::FractionalForm;
use ddg_repulsive_surfaces::remesh::RemeshOrder;
use ddg_repulsive_surfaces::solve::precond::PreconditionerForm;
use ddg_repulsive_surfaces::step::{CSV_HEADER, csv_line};
use ddg_repulsive_surfaces::{
    ArmijoRule, Backend, BarycenterMode, BhOptions, Constraint, EnergyEval, Exponents, FlowOptions, GmresSide,
    NewtonSteps, Session, StartRule, StopRule, Surface, ccd, energy, read_obj, write_obj,
};
use rayon::prelude::*;

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

/// (Q2, Q4, Q6) of the area distribution over the directions from the area centroid.
fn harmonic_power(faces: &[[u32; 3]], x: &[P]) -> [f64; 3] {
    let mut a = Vec::with_capacity(faces.len());
    let mut c = Vec::with_capacity(faces.len());
    for f in faces {
        let (p, q, r) = (x[f[0] as usize], x[f[1] as usize], x[f[2] as usize]);
        a.push(0.5 * norm(cross(sub(q, p), sub(r, p))));
        c.push([0, 1, 2].map(|k| (p[k] + q[k] + r[k]) / 3.0));
    }
    let at: f64 = a.iter().sum();
    let g = (0..faces.len()).fold([0.0; 3], |g, i| [0, 1, 2].map(|k| g[k] + a[i] * c[i][k] / at));
    let u: Vec<P> = c
        .iter()
        .map(|p| {
            let d = sub(*p, g);
            let n = norm(d).max(1e-300);
            d.map(|v| v / n)
        })
        .collect();
    let s = (0..faces.len())
        .into_par_iter()
        .map(|i| {
            let mut acc = [0.0; 3];
            for j in 0..faces.len() {
                let t = dot(u[i], u[j]);
                let t2 = t * t;
                let p2 = 0.5 * (3.0 * t2 - 1.0);
                let p4 = (35.0 * t2 * t2 - 30.0 * t2 + 3.0) / 8.0;
                let p6 = (231.0 * t2 * t2 * t2 - 315.0 * t2 * t2 + 105.0 * t2 - 5.0) / 16.0;
                let w = a[i] * a[j];
                acc[0] += w * p2;
                acc[1] += w * p4;
                acc[2] += w * p6;
            }
            acc
        })
        .reduce(|| [0.0; 3], |p, q| [p[0] + q[0], p[1] + q[1], p[2] + q[2]]);
    s.map(|v| v / (at * at))
}

fn sym_line(faces: &[[u32; 3]], x: &[P]) -> String {
    let (_, _, desc) = apps::inertia_indicators(faces, x);
    let q = harmonic_power(faces, x);
    format!("{desc} Q2={:.3e} Q4={:.3e} Q6={:.3e}", q[0], q[1], q[2])
}

fn exact_energy(s: &Surface, x: &[P], threads: usize) -> f64 {
    energy(s, x, Exponents::default(), &EnergyEval::Exact, &[], Some(threads)).map(|r| r.energy).unwrap_or(f64::NAN)
}

fn dodecahedron() -> (Vec<P>, Vec<[usize; 2]>) {
    let phi = (1.0 + 5f64.sqrt()) / 2.0;
    let ip = 1.0 / phi;
    let mut v: Vec<P> = Vec::new();
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                v.push([sx, sy, sz]);
            }
        }
    }
    for s1 in [-1.0, 1.0] {
        for s2 in [-1.0, 1.0] {
            v.push([0.0, s1 * ip, s2 * phi]);
            v.push([s1 * ip, s2 * phi, 0.0]);
            v.push([s1 * phi, 0.0, s2 * ip]);
        }
    }
    let el = 2.0 * ip;
    let mut e = Vec::new();
    for i in 0..v.len() {
        for j in i + 1..v.len() {
            if (norm(sub(v[i], v[j])) - el).abs() < 1e-9 {
                e.push([i, j]);
            }
        }
    }
    assert_eq!((v.len(), e.len()), (20, 30));
    (v, e)
}

fn cmd_gen(row: &str, out: &str, faces: usize, radius: f64) {
    let (s0, x0) = read_obj(row).unwrap();
    let (a0, _) = apps::area_volume(s0.faces(), &x0);
    let (v, e) = dodecahedron();
    let el = norm(sub(v[e[0][0]], v[e[0][1]]));
    let curves: Vec<(Vec<P>, bool)> = e
        .iter()
        .map(|&[i, j]| ((0..=8).map(|k| [0, 1, 2].map(|c| v[i][c] + (v[j][c] - v[i][c]) * k as f64 / 8.0)).collect(), false))
        .collect();
    let m = shapes::graph_tube(&curves, radius * el, faces).unwrap();
    let k = (a0 / m.area()).sqrt();
    let m = m.transformed(k, [0.0; 3]);
    let s = m.surface().unwrap();
    let inter = ccd::intersects(&s, &m.positions).unwrap();
    write_obj(out, &s, &m.positions).unwrap();
    println!(
        "GEN {out} faces={} genus={} area={:.6} (row {a0:.6}) scale={k:.4} intersects={inter} {}",
        s.n_faces(),
        1 - s.euler_characteristic() / 2,
        m.area(),
        sym_line(s.faces(), &m.positions)
    );
}

/// The paper video's genus-11 start (2:22.5): a chain of `n` round tori (axis z) along the x axis at unit
/// spacing, outer diameter `outer`, hole (throat) diameter `hole`, so major radius (outer + hole) / 4 and tube
/// radius (outer − hole) / 4; neighbours are fused by a smooth union whose blend width is chosen so the neck
/// between two rings is `waist` wide (in y, at the mid-plane x = 1/2). Meshed by marching tetrahedra with
/// about `faces` triangles and scaled to the area of ROW_OBJ.
fn cmd_chain(row: &str, out: &str, hole: f64, outer: f64, waist: f64, n: usize, faces: usize) {
    use ddg_repulsive_surfaces::generate::field::{Field, smooth_min};
    use ddg_repulsive_surfaces::generate::tetra::mesh_with_faces;
    let (s0, x0) = read_obj(row).unwrap();
    let (a0, _) = apps::area_volume(s0.faces(), &x0);
    let (big, r) = ((outer + hole) / 4.0, (outer - hole) / 4.0);
    // neck: at (1/2, waist/2, 0) the blended field is 0: d - k/4 = 0 (both tori at the same distance d)
    let d = ((0.25 + waist * waist / 4.0).sqrt() - big) - r;
    let k = (4.0 * d).max(0.0);
    let cx: Vec<f64> = (0..n).map(|i| i as f64 - (n as f64 - 1.0) / 2.0).collect();
    let cxf = cx.clone();
    let f: Field = std::sync::Arc::new(move |p: P| {
        let mut v = f64::INFINITY;
        for &c in &cxf {
            let rho = ((p[0] - c).powi(2) + p[1] * p[1]).sqrt();
            let t = ((rho - big).powi(2) + p[2] * p[2]).sqrt() - r;
            v = if v.is_finite() { smooth_min(v, t, k) } else { t };
        }
        v
    });
    let m = 0.25;
    let hx = cx[n - 1] + outer / 2.0 + m;
    let hy = outer / 2.0 + m;
    let hz = r + k / 4.0 + m;
    let mesh = mesh_with_faces(&f, [-hx, -hy, -hz], [hx, hy, hz], faces).unwrap();
    let sc = (a0 / mesh.area()).sqrt();
    let mesh = mesh.transformed(sc, [0.0; 3]);
    let s = mesh.surface().unwrap();
    let inter = ccd::intersects(&s, &mesh.positions).unwrap();
    write_obj(out, &s, &mesh.positions).unwrap();
    println!(
        "CHAIN {out} rings={n} hole={hole} outer={outer} waist={waist} R={big:.4} r={r:.4} blend={k:.4} faces={} \
         components={} genus={} area={:.6} (row {a0:.6}) scale={sc:.4} intersects={inter} {}",
        s.n_faces(),
        s.n_components(),
        1 - s.euler_characteristic() / 2,
        mesh.area(),
        sym_line(s.faces(), &mesh.positions)
    );
}

fn cmd_measure(paths: &[String], area: Option<f64>, threads: usize) {
    for p in paths {
        let (s, x) = read_obj(p).unwrap();
        let e = exact_energy(&s, &x, threads);
        let (a, _) = apps::area_volume(s.faces(), &x);
        let en = area.map_or(f64::NAN, |a0| e * a / a0);
        let inter = ccd::intersects(&s, &x).unwrap();
        println!(
            "MEASURE {p} faces={} genus={} area={a:.6} exact_energy={e:.6} energy_at_area={en:.6} intersects={inter} {}",
            s.n_faces(),
            1 - s.euler_characteristic() / 2,
            sym_line(s.faces(), &x)
        );
    }
}

struct FlowArgs {
    out: String,
    threads: usize,
    max_steps: usize,
    window: usize,
    rel: f64,
    frame_every: usize,
    perturb: f64,
    seed: u64,
}

fn perturb(s: &Surface, x: &mut [P], eps: f64, seed: u64) {
    let (a0, _) = apps::area_volume(s.faces(), x);
    let l: f64 = s.edges().iter().map(|e| norm(sub(x[e[0] as usize], x[e[1] as usize]))).sum::<f64>() / s.n_edges() as f64;
    let mut st: u64 = 0x9E3779B97F4A7C15 ^ seed.wrapping_mul(0xD1B54A32D192ED03);
    let mut rnd = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    };
    for _ in 0..10 {
        rnd();
    }
    for p in x.iter_mut() {
        let d = loop {
            let d = [rnd(), rnd(), rnd()];
            let n = norm(d);
            if n > 1e-3 && n <= 1.0 {
                break d.map(|v| v / n);
            }
        };
        for k in 0..3 {
            p[k] += eps * l * d[k];
        }
    }
    let (a1, _) = apps::area_volume(s.faces(), x);
    let k = (a0 / a1).sqrt();
    let c = x.iter().fold([0.0; 3], |c, p| [0, 1, 2].map(|i| c[i] + p[i] / x.len() as f64));
    for p in x.iter_mut() {
        *p = [0, 1, 2].map(|i| c[i] + k * (p[i] - c[i]));
    }
}

fn cmd_flow(path: &str, a: &FlowArgs, args: &[String]) {
    let has = |k: &str| args.iter().any(|x| x == k);
    let get = |k: &str| args.iter().position(|x| x == k).map(|i| args[i + 1].clone());
    let (s, mut x) = read_obj(path).unwrap();
    if a.perturb > 0.0 {
        perturb(&s, &mut x, a.perturb, a.seed);
    }
    let mut o = FlowOptions {
        max_iterations: a.max_steps,
        stop: StopRule::Off,
        threads: Some(a.threads),
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        ..FlowOptions::default()
    };
    if has("--no-area") {
        o.constraints.retain(|c| !matches!(c, Constraint::Area { .. }));
    }
    o.reference_compat = has("--compat");
    if has("--no-ccd") {
        o.line_search.ccd = None;
    }
    if has("--armijo-ref") {
        o.line_search.armijo = ArmijoRule::Reference { c: 1e-2 };
    }
    if has("--start-ref") {
        o.line_search.start = StartRule::Reference;
    }
    if let Some(t) = get("--theta") {
        o.energy_eval = EnergyEval::BarnesHut(BhOptions { theta: t.parse().unwrap(), ..BhOptions::default() });
    }
    if has("--smooth-first") {
        o.remesh.as_mut().unwrap().order = RemeshOrder::SmoothThenFlip;
    }
    match get("--bary").as_deref() {
        Some("global") => o.barycenter_mode = Some(BarycenterMode::Global),
        Some("percomp") => o.barycenter_mode = Some(BarycenterMode::PerComponent),
        Some(m) => panic!("unknown barycenter mode {m}"),
        None => {}
    }
    if has("--no-remesh") {
        o.remesh = None;
    }
    if has("--newton1") {
        o.newton_steps = NewtonSteps::Fixed(1);
    }
    if let Backend::Auto { hierarchical: h, .. } = &mut o.backend {
        if has("--precond-ref") {
            h.preconditioner.fractional = FractionalForm::FacePairs;
            h.preconditioner.form = PreconditionerForm::Congruence;
            h.preconditioner.eps_multi_rel = 1e-2;
        }
        if has("--gmres-left") {
            h.gmres.side = GmresSide::Left;
        }
        if let Some(t) = get("--gmres-tol") {
            h.gmres.tolerance = t.parse().unwrap();
        }
    }
    std::fs::create_dir_all(&a.out).unwrap();
    let mut csv = std::io::BufWriter::new(std::fs::File::create(format!("{}/steps.csv", a.out)).unwrap());
    writeln!(csv, "{CSV_HEADER}").unwrap();
    let t0 = Instant::now();
    println!("START {path} faces={} {}", s.n_faces(), sym_line(s.faces(), &x));
    write_obj(format!("{}/frame_00000.obj", a.out), &s, &x).unwrap();
    let mut sess = Session::new(s, x, o).unwrap();
    let mut hist: Vec<f64> = Vec::new();
    let mut status = "MaxSteps".to_string();
    while sess.iteration() < a.max_steps {
        let r = match sess.step() {
            Ok(r) => r,
            Err(e) => {
                status = format!("Error({e})");
                break;
            }
        };
        writeln!(csv, "{}", csv_line(&r)).unwrap();
        hist.push(r.energy_after);
        let it = sess.iteration();
        if it % 25 == 0 {
            csv.flush().unwrap();
            println!("PROGRESS step={it} faces={} energy={:.6e} t={:.0}s", sess.surface().n_faces(), r.energy_after, t0.elapsed().as_secs_f64());
        }
        if it % 100 == 0 {
            println!("SYM step={it} energy={:.6e} {}", r.energy_after, sym_line(sess.surface().faces(), sess.positions()));
        }
        if a.frame_every > 0 && it % a.frame_every == 0 {
            write_obj(format!("{}/frame_{it:05}.obj", a.out), sess.surface(), sess.positions()).unwrap();
        }
        if sess.status().is_terminal() {
            status = format!("{:?}", sess.status());
            break;
        }
        let n = hist.len();
        let w = a.window;
        if n >= 2 * w {
            let m = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
            let (last, prev) = (m(&hist[n - w..]), m(&hist[n - 2 * w..n - w]));
            if (prev - last).abs() <= a.rel * last.abs() {
                status = "Converged".into();
                break;
            }
        }
    }
    csv.flush().unwrap();
    write_obj(format!("{}/final.obj", a.out), sess.surface(), sess.positions()).unwrap();
    let e_ex = exact_energy(sess.surface(), sess.positions(), a.threads);
    let inter = ccd::intersects(sess.surface(), sess.positions()).unwrap();
    println!(
        "RESULT steps={} status={status} F={:.6} exact={e_ex:.6} faces={} genus={} intersects={inter} seconds={:.0} {}",
        sess.iteration(),
        hist.last().copied().unwrap_or(f64::NAN),
        sess.surface().n_faces(),
        1 - sess.surface().euler_characteristic() / 2,
        t0.elapsed().as_secs_f64(),
        sym_line(sess.surface().faces(), sess.positions())
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| args.iter().position(|a| a == k).map(|i| args[i + 1].clone());
    let threads: usize = get("--threads").map_or(4, |v| v.parse().unwrap());
    match args.first().map(String::as_str) {
        Some("gen") => cmd_gen(
            &args[1],
            &args[2],
            get("--faces").map_or(10136, |v| v.parse().unwrap()),
            get("--radius").map_or(0.3, |v| v.parse().unwrap()),
        ),
        Some("ladder") => {
            // a chain of 11 rings: the ladder graph with 12 rungs (24 vertices, 34 edges, genus 11), unit spacing
            let r = get("--radius").map_or(0.22, |v| v.parse().unwrap());
            let (s0, x0) = read_obj(&args[1]).unwrap();
            let (a0, _) = apps::area_volume(s0.faces(), &x0);
            let seg = |a: P, b: P| ((0..=8).map(|k| [0, 1, 2].map(|c| a[c] + (b[c] - a[c]) * k as f64 / 8.0)).collect(), false);
            let mut curves: Vec<(Vec<P>, bool)> = Vec::new();
            for i in 0..12 {
                let x = i as f64 - 5.5;
                curves.push(seg([x, -0.5, 0.0], [x, 0.5, 0.0]));
                if i < 11 {
                    curves.push(seg([x, -0.5, 0.0], [x + 1.0, -0.5, 0.0]));
                    curves.push(seg([x, 0.5, 0.0], [x + 1.0, 0.5, 0.0]));
                }
            }
            let m = shapes::graph_tube(&curves, r, get("--faces").map_or(10136, |v| v.parse().unwrap())).unwrap();
            let m = m.clone().transformed((a0 / m.area()).sqrt(), [0.0; 3]);
            let s = m.surface().unwrap();
            write_obj(&args[2], &s, &m.positions).unwrap();
            println!("LADDER {} faces={} genus={} area={:.6} intersects={}", args[2], s.n_faces(), 1 - s.euler_characteristic() / 2, m.area(), ccd::intersects(&s, &m.positions).unwrap());
        }
        Some("chain") => cmd_chain(
            &args[1],
            &args[2],
            get("--hole").map_or(0.3, |v| v.parse().unwrap()),
            get("--outer").map_or(0.86, |v| v.parse().unwrap()),
            get("--waist").map_or(0.3, |v| v.parse().unwrap()),
            get("--rings").map_or(11, |v| v.parse().unwrap()),
            get("--faces").map_or(12000, |v| v.parse().unwrap()),
        ),
        Some("perturb") => {
            // the start of `flow --perturb EPS --seed N`, written as OBJ (for the reference program)
            let (s, mut x) = read_obj(&args[1]).unwrap();
            perturb(&s, &mut x, get("--perturb").unwrap().parse().unwrap(), get("--seed").unwrap().parse().unwrap());
            write_obj(&args[2], &s, &x).unwrap();
        }
        Some("measure") => {
            let paths: Vec<String> = args[1..].iter().take_while(|a| !a.starts_with("--")).cloned().collect();
            cmd_measure(&paths, get("--area").map(|v| v.parse().unwrap()), threads)
        }
        Some("flow") => {
            let fa = FlowArgs {
                out: get("--out").expect("--out DIR"),
                threads,
                max_steps: get("--max-steps").map_or(8000, |v| v.parse().unwrap()),
                window: get("--window").map_or(200, |v| v.parse().unwrap()),
                rel: get("--rel").map_or(1e-5, |v| v.parse().unwrap()),
                frame_every: get("--frame-every").map_or(100, |v| v.parse().unwrap()),
                perturb: get("--perturb").map_or(0.0, |v| v.parse().unwrap()),
                seed: get("--seed").map_or(0, |v| v.parse().unwrap()),
            };
            cmd_flow(&args[1], &fa, &args)
        }
        _ => eprintln!("usage: see the header of examples/dodeca.rs"),
    }
}
