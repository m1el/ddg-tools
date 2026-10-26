//! Multigrid on curve networks with high-degree junctions and pinned endpoints (`SPEC.md` §8.4,
//! §17g).

#![allow(clippy::needless_range_loop)]

use std::sync::Arc;
use std::time::Instant;

use ddg_core::sparse::SparseMatrix;
use ddg_repulsive_curves::constraints::{Constraint, ConstraintSet};
use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::metric::GramQuadrature;
use ddg_repulsive_curves::multigrid::vcycle::Hierarchy;
use ddg_repulsive_curves::multigrid::{HierarchyInput, MgOptions, MultigridSaddle};
use ddg_repulsive_curves::scene::{Scene, read_scene};
use ddg_repulsive_curves::solve::{dense_sobolev_solver, projected_gradient};
use ddg_repulsive_curves::spatial::bct::BctOptions;
use ddg_repulsive_curves::{CurveNetwork, Exponents, RcError, read_curve_obj};

type P = [f64; 3];

fn data(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/repulsive-curves")
        .join(rel)
}

fn pool(n: usize) -> Arc<rayon::ThreadPool> {
    Arc::new(make_pool(Some(n)).unwrap())
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn rel_err(a: &[f64], b: &[f64]) -> f64 {
    let d: Vec<f64> = a.iter().zip(b).map(|(p, q)| p - q).collect();
    norm(&d) / norm(b)
}

fn mul(c: &SparseMatrix, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; c.nrows()];
    c.mul_vec(x, &mut y).unwrap();
    y
}

struct Problem {
    net: CurveNetwork,
    pos: Vec<P>,
    ex: Exponents,
    set: ConstraintSet,
    c: SparseMatrix,
    diff: Vec<P>,
}

fn problem(net: CurveNetwork, pos: Vec<P>, ex: Exponents, cons: &[Constraint]) -> Problem {
    let set = ConstraintSet::new(&net, &pos, cons, None).unwrap();
    set.analyze(&net, &pos).unwrap();
    let c = set.jacobian(&net, &pos).unwrap();
    let p = pool(16);
    let (_, diff) = energy_and_differential(&p, &net, &pos, ex);
    Problem {
        net,
        pos,
        ex,
        set,
        c,
        diff,
    }
}

fn mg_solver(pr: &Problem, threads: usize, o: MgOptions) -> Result<MultigridSaddle, String> {
    let p = pool(threads);
    let input = HierarchyInput {
        network: &pr.net,
        positions: &pr.pos,
        exponents: pr.ex,
        quadrature: GramQuadrature::Midpoint,
        layout: pr.set.layout(),
        jacobian: &pr.c,
        constraints: Some(&pr.set),
    };
    MultigridSaddle::new(&p, &p, &input, o).map_err(|e| e.to_string())
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }
}

fn probe_smoother() -> f64 {
    std::env::var("RC_PROBE_SMOOTHER")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1.0)
}

fn probe_bct() -> BctOptions {
    let theta: f64 = std::env::var("RC_PROBE_THETA")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.25);
    BctOptions {
        theta_x: theta,
        theta_t: theta,
        ..BctOptions::default()
    }
}

fn capillary(high: bool) -> Option<Scene> {
    let dir = data("scenes/CapillaryNetwork");
    if !dir.join("scene.txt").exists() {
        eprintln!("CapillaryNetwork scene not present; skipping");
        return None;
    }
    let scene = read_scene(dir.join("scene.txt")).unwrap();
    if high {
        return Some(scene);
    }
    let (net, pos) = read_curve_obj(dir.join("two-networks.obj")).unwrap();
    Some(Scene::new(scene.name.clone(), net, pos, scene.spec.clone()).unwrap())
}

fn describe(pr: &Problem) {
    let n = pr.net.n_vertices();
    let mut hist = [0usize; 8];
    for v in 0..n as u32 {
        hist[pr.net.incident_edges(v).len().min(7)] += 1;
    }
    eprintln!(
        "{} vertices, {} edges, {} components, degree histogram {:?}, {} rows ({} dropped)",
        n,
        pr.net.n_edges(),
        pr.net.components().len(),
        hist,
        pr.c.nrows(),
        pr.set.dropped_rows().len()
    );
}

fn describe_levels(mg: &MultigridSaddle) {
    for (l, lv) in mg.hierarchy().levels().iter().enumerate() {
        let n = lv.network().n_vertices();
        let mut hist = [0usize; 8];
        for v in 0..n as u32 {
            hist[lv.network().incident_edges(v).len().min(7)] += 1;
        }
        eprintln!(
            "  level {l}: {n} vertices, {} edges, rows {} (dropped {}), ker dim {}, degrees {:?}",
            lv.network().n_edges(),
            lv.jacobian().nrows(),
            lv.dropped_rows().len(),
            3 * n as i64 - lv.jacobian().nrows() as i64,
            hist
        );
    }
}

/// Rayleigh quotients of the Galerkin coarse operator against the rediscretized one:
/// `(Πx)ᵀ H_ℓ (Πx) / xᵀ H_{ℓ+1} x` for random `x` in `ker C_{ℓ+1}`.
fn level_ratios(mg: &MultigridSaddle) {
    let h = mg.hierarchy();
    for l in 0..h.n_levels() - 1 {
        let mut rng = Rng(7 + l as u64);
        let mut ratios = Vec::new();
        for _ in 0..6 {
            let x: Vec<f64> = (0..h.len(l + 1)).map(|_| rng.next()).collect();
            let x = h.levels()[l + 1].projector().project(&x).unwrap();
            let hx = h.apply(l + 1, &x).unwrap();
            let px = h.prolong(l, &x).unwrap();
            let hpx = h.apply(l, &px).unwrap();
            let num: f64 = px.iter().zip(&hpx).map(|(a, b)| a * b).sum();
            let den: f64 = x.iter().zip(&hx).map(|(a, b)| a * b).sum();
            ratios.push((num / den, norm(&px) / norm(&x)));
        }
        eprintln!(
            "  level {}: {} corrected pairs",
            l + 1,
            h.levels()[l + 1].near_correction().pairs().len()
        );
        eprintln!(
            "  level {l}->{}: Galerkin/rediscretized, |Πx|/|x|: {ratios:.3?}",
            l + 1
        );
    }
}

fn smooth_ratios(mg: &MultigridSaddle, pr: &Problem) {
    let layout = pr.set.layout();

    // Smooth test vectors: approximate coarse solutions for the restricted right side.
    let h = mg.hierarchy();
    let a = layout.gather(&pr.diff);
    let mut b = h.levels()[0].projector().project(&a).unwrap();
    for l in 0..h.n_levels() - 1 {
        let bc = h.restrict(l, &b).unwrap();
        let x = ddg_repulsive_curves::multigrid::vcycle::cg_smooth(
            h,
            l + 1,
            &bc,
            vec![0.0; bc.len()],
            300,
        )
        .unwrap();
        let e_c: f64 = x
            .iter()
            .zip(&h.apply(l + 1, &x).unwrap())
            .map(|(p, q)| p * q)
            .sum();
        let px = h.prolong(l, &x).unwrap();
        let e_f: f64 = px
            .iter()
            .zip(&h.apply(l, &px).unwrap())
            .map(|(p, q)| p * q)
            .sum();
        let jx = h.levels()[l]
            .coarsening()
            .unwrap()
            .prolong(&x, layout.dim())
            .unwrap();
        let e_j: f64 = jx
            .iter()
            .zip(&h.apply_metric(l, &jx).unwrap())
            .map(|(p, q)| p * q)
            .sum();
        let pjx = h.levels()[l].projector().project(&jx).unwrap();
        eprintln!(
            "  smooth {l}->{}: Galerkin/redisc {:.3}, unprojected J/redisc {:.3}, |Jx - PJx|/|Jx| {:.2e}",
            l + 1,
            e_f / e_c,
            e_j / e_c,
            rel_err(&pjx, &jx)
        );
        b = bc;
    }
}

#[test]
#[ignore]
fn probe_capillary() {
    let which = std::env::var("RC_PROBE_SIZE").unwrap_or_else(|_| "low".into());
    let Some(scene) = capillary(which == "high") else {
        return;
    };
    let t = Instant::now();
    let pr = problem(
        scene.network.clone(),
        scene.positions.clone(),
        scene.spec.exponents,
        &scene.constraints(),
    );
    describe(&pr);
    let layout = pr.set.layout();
    let theta: f64 = std::env::var("RC_PROBE_THETA")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.25);
    let env = |k: &str, d: usize| -> usize {
        std::env::var(k)
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(d)
    };
    let o = MgOptions {
        bct: BctOptions {
            theta_x: theta,
            theta_t: theta,
            ..BctOptions::default()
        },
        max_levels: env("RC_PROBE_LEVELS", 20),
        coarsest_vertices: env("RC_PROBE_COARSEST", 512),
        smoothing_steps: env("RC_PROBE_NU", 12),
        max_cycles: env("RC_PROBE_CYCLES", 20),
        near_field_factor: std::env::var("RC_PROBE_KAPPA")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(1.0),
        smoother_factor: probe_smoother(),
        energy_check: std::env::var("RC_PROBE_ENERGY").map_or(true, |v| v == "1"),
        max_fixed_fraction: std::env::var("RC_PROBE_FIXED")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(MgOptions::default().max_fixed_fraction),
        ..MgOptions::default()
    };
    let mut mg = match mg_solver(&pr, 32, o) {
        Ok(m) => m,
        Err(e) => panic!("build: {e}"),
    };
    eprintln!("hierarchy built in {:.1} s", t.elapsed().as_secs_f64());
    if std::env::var("RC_PROBE_VERBOSE").is_ok() {
        describe_levels(&mg);
        level_ratios(&mg);
        smooth_ratios(&mg, &pr);
    }
    if std::env::var("RC_PROBE_TIME").is_ok() {
        let h = mg.hierarchy();
        for l in 0..h.n_levels() {
            let x: Vec<f64> = (0..h.len(l)).map(|i| (i as f64 * 0.37).sin()).collect();
            let t = Instant::now();
            for _ in 0..5 {
                std::hint::black_box(h.apply_metric(l, &x).unwrap());
            }
            let tm = t.elapsed().as_secs_f64() / 5.0;
            let t = Instant::now();
            for _ in 0..5 {
                std::hint::black_box(h.levels()[l].projector().project(&x).unwrap());
            }
            let tp = t.elapsed().as_secs_f64() / 5.0;
            let t = Instant::now();
            for _ in 0..5 {
                std::hint::black_box(h.precondition(l, &x).unwrap());
            }
            let tq = t.elapsed().as_secs_f64() / 5.0;
            eprintln!(
                "  level {l}: metric {:.1} ms, project {:.1} ms, precondition {:.1} ms",
                tm * 1e3,
                tp * 1e3,
                tq * 1e3
            );
        }
    }
    let t = Instant::now();
    let gm = projected_gradient(&mut mg, layout, &pr.diff);
    eprintln!("mg solve {:.2} s", t.elapsed().as_secs_f64());
    match &gm {
        Ok(g) => eprintln!(
            "mg: {} cycles, residual {:.2e}, cycle residuals {:?}",
            g.linear.iterations,
            g.linear.relative_residual,
            mg.last_reports()
                .iter()
                .map(|r| r.residuals.clone())
                .collect::<Vec<_>>()
        ),
        Err(e) => eprintln!("mg: {e}"),
    }
    if let Ok(path) = std::env::var("RC_PROBE_DENSE") {
        let xd: Vec<f64> = match std::fs::read_to_string(&path) {
            Ok(text) => text.lines().map(|l| l.parse().unwrap()).collect(),
            Err(_) => {
                let t = Instant::now();
                let mut dense = dense_sobolev_solver(
                    pool(32),
                    &pr.net,
                    &pr.pos,
                    pr.ex,
                    GramQuadrature::Midpoint,
                    &pr.set,
                )
                .unwrap();
                let gd = projected_gradient(&mut dense, layout, &pr.diff).unwrap();
                eprintln!(
                    "dense: {:.1} s, residual {:.2e}",
                    t.elapsed().as_secs_f64(),
                    gd.linear.relative_residual
                );
                let xd = layout.gather(&gd.direction);
                let text: Vec<String> = xd.iter().map(|v| format!("{v:e}")).collect();
                let _ = std::fs::write(&path, text.join("\n"));
                xd
            }
        };
        if let Ok(g) = gm {
            let xm = layout.gather(&g.direction);
            eprintln!(
                "|g_mg - g_d|/|g_d| = {:.2e}, |C g_mg| = {:.2e}",
                rel_err(&xm, &xd),
                norm(&mul(&pr.c, &xm))
            );
        }
    }
}

/// The hierarchy with Galerkin coarse operators `Π_ℓᵀ H_ℓ Π_ℓ` (matrix-free) and a CG
/// coarsest solve, for comparison with the rediscretized one.
struct Galerkin<'a> {
    h: &'a ddg_repulsive_curves::multigrid::CurveHierarchy,
}

impl Hierarchy for Galerkin<'_> {
    fn n_levels(&self) -> usize {
        self.h.n_levels()
    }
    fn len(&self, level: usize) -> usize {
        self.h.len(level)
    }
    fn apply(&self, level: usize, x: &[f64]) -> Result<Vec<f64>, ddg_repulsive_curves::RcError> {
        if level == 0 {
            return self.h.apply(0, x);
        }
        let y = self.h.prolong(level - 1, x)?;
        let z = self.apply(level - 1, &y)?;
        self.h.restrict(level - 1, &z)
    }
    fn prolong(
        &self,
        level: usize,
        coarse: &[f64],
    ) -> Result<Vec<f64>, ddg_repulsive_curves::RcError> {
        self.h.prolong(level, coarse)
    }
    fn restrict(
        &self,
        level: usize,
        fine: &[f64],
    ) -> Result<Vec<f64>, ddg_repulsive_curves::RcError> {
        self.h.restrict(level, fine)
    }
    fn coarse_solve(&self, b: &[f64]) -> Result<Vec<f64>, ddg_repulsive_curves::RcError> {
        let l = self.n_levels() - 1;
        let bn = norm(b);
        let mut x = vec![0.0; b.len()];
        let mut r = b.to_vec();
        let mut p = r.clone();
        let mut rr: f64 = r.iter().map(|v| v * v).sum();
        for _ in 0..3000 {
            if rr.sqrt() <= 1e-10 * bn {
                break;
            }
            let q = self.apply(l, &p)?;
            let pq: f64 = p.iter().zip(&q).map(|(a, b)| a * b).sum();
            let alpha = rr / pq;
            for i in 0..x.len() {
                x[i] += alpha * p[i];
                r[i] -= alpha * q[i];
            }
            let rn: f64 = r.iter().map(|v| v * v).sum();
            for i in 0..p.len() {
                p[i] = r[i] + rn / rr * p[i];
            }
            rr = rn;
        }
        Ok(x)
    }
}

#[test]
#[ignore]
fn probe_galerkin() {
    let which = std::env::var("RC_PROBE_SIZE").unwrap_or_else(|_| "low".into());
    let pr = if which == "hub" {
        let (net, pos) = hub_network(40, 6, 0.05);
        let cons = pinned_edge_lengths(&net, &pos);
        problem(net, pos, Exponents::default(), &cons)
    } else {
        let Some(scene) = capillary(which == "high") else {
            return;
        };
        problem(
            scene.network.clone(),
            scene.positions.clone(),
            scene.spec.exponents,
            &scene.constraints(),
        )
    };
    let env = |k: &str, d: usize| -> usize {
        std::env::var(k)
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(d)
    };
    let o = MgOptions {
        max_levels: env("RC_PROBE_LEVELS", 3),
        coarsest_vertices: 1,
        bct: probe_bct(),
        smoother_factor: probe_smoother(),
        near_field_factor: std::env::var("RC_PROBE_KAPPA")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(1.0),
        ..MgOptions::default()
    };
    let mg = mg_solver(&pr, 32, o).unwrap();
    describe_levels(&mg);
    let layout = pr.set.layout();
    let a = layout.gather(&pr.diff);
    let h = mg.hierarchy();
    let b = h.levels()[0].projector().project(&a).unwrap();
    let cyc = ddg_repulsive_curves::multigrid::CycleOptions {
        smoothing_steps: 12,
        initial_smoothing_steps: 24,
        tolerance: 1e-6,
        max_cycles: 20,
        flexible_cg: true,
    };
    let t = Instant::now();
    match ddg_repulsive_curves::multigrid::vcycle::solve(h, &b, &cyc) {
        Ok((_, r)) => eprintln!(
            "rediscretized: {:.1} s {:?}",
            t.elapsed().as_secs_f64(),
            r.residuals
        ),
        Err(e) => eprintln!("rediscretized: {e}"),
    }
    let t = Instant::now();
    let g = Galerkin { h };
    match ddg_repulsive_curves::multigrid::vcycle::solve(&g, &b, &cyc) {
        Ok((_, r)) => eprintln!(
            "galerkin: {:.1} s {:?}",
            t.elapsed().as_secs_f64(),
            r.residuals
        ),
        Err(e) => eprintln!("galerkin: {e}"),
    }
    // Plain CG on the finest level for reference.
    let t = Instant::now();
    let mut x = vec![0.0; b.len()];
    let mut r = b.clone();
    let mut p = r.clone();
    let mut rr: f64 = r.iter().map(|v| v * v).sum();
    let bn = norm(&b);
    let mut hist = Vec::new();
    for k in 0..2000 {
        let q = h.apply(0, &p).unwrap();
        let pq: f64 = p.iter().zip(&q).map(|(a, b)| a * b).sum();
        let alpha = rr / pq;
        for i in 0..x.len() {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }
        let rn: f64 = r.iter().map(|v| v * v).sum();
        for i in 0..p.len() {
            p[i] = r[i] + rn / rr * p[i];
        }
        rr = rn;
        if k % 100 == 99 {
            hist.push(rr.sqrt() / bn);
        }
        if rr.sqrt() < 1e-6 * bn {
            eprintln!("cg converged in {k}");
            break;
        }
    }
    eprintln!(
        "plain CG: {:.1} s, every 100: {hist:?}",
        t.elapsed().as_secs_f64()
    );
}

// ---------------------------------------------------------------- synthetic network

/// Builder of a network made of straight rods between given points.
struct Builder {
    pos: Vec<P>,
    edges: Vec<[u32; 2]>,
}

impl Builder {
    fn vertex(&mut self, p: P) -> u32 {
        self.pos.push(p);
        self.pos.len() as u32 - 1
    }
    /// A rod from vertex `a` to `b` (or to a new end point) with `n` segments, bent by
    /// `bulge · sin(π t)`.
    fn rod(&mut self, a: u32, b: Option<u32>, end: P, n: usize, bulge: P) -> u32 {
        let p = self.pos[a as usize];
        let q = b.map_or(end, |b| self.pos[b as usize]);
        let mut prev = a;
        for k in 1..n {
            let t = k as f64 / n as f64;
            let w = (std::f64::consts::PI * t).sin();
            let v = self.vertex([
                p[0] + t * (q[0] - p[0]) + w * bulge[0],
                p[1] + t * (q[1] - p[1]) + w * bulge[1],
                p[2] + t * (q[2] - p[2]) + w * bulge[2],
            ]);
            self.edges.push([prev, v]);
            prev = v;
        }
        let last = b.unwrap_or_else(|| self.vertex(q));
        self.edges.push([prev, last]);
        last
    }
}

/// Three hubs of degree 5, 4 and 3 joined by rods, free arms with pinned ends, and `crossings`
/// open rods that pass over the A–B rod at height `gap · h` exactly above edge midpoints (the
/// near-contacts of the capillary scene). `h = 1/m` is the edge length of the A–B rod.
fn hub_network(m: usize, crossings: usize, gap: f64) -> (CurveNetwork, Vec<P>) {
    hub_network_at(m, crossings, gap, 0.5, 0.0)
}

/// [`hub_network`] with the crossings at `frac` along the A–B edge and the crossing rods'
/// edge midpoints shifted by `shift · h` from the crossing.
fn hub_network_at(
    m: usize,
    crossings: usize,
    gap: f64,
    frac: f64,
    shift: f64,
) -> (CurveNetwork, Vec<P>) {
    let h = 1.0 / m as f64;
    let mut b = Builder {
        pos: Vec::new(),
        edges: Vec::new(),
    };
    let a = b.vertex([0.0, 0.0, 0.0]);
    let bb = b.vertex([1.0, 0.0, 0.0]);
    let c = b.vertex([0.5, 0.8, 0.1]);
    // A–B is straight in the plane z = 0 (the crossings sit over its midpoints) and bends
    // in y only.
    b.rod(a, Some(bb), [0.0; 3], m, [0.0, -0.05, 0.0]);
    b.rod(a, Some(c), [0.0; 3], m, [-0.1, 0.05, 0.05]);
    b.rod(bb, Some(c), [0.0; 3], m, [0.1, 0.05, -0.05]);
    let arm = |b: &mut Builder, from: u32, deg: f64, len: f64, lift: f64| {
        let p = b.pos[from as usize];
        let t = deg.to_radians();
        let end = [p[0] + len * t.cos(), p[1] + len * t.sin(), p[2] + lift];
        let bulge = [0.1 * t.sin(), -0.1 * t.cos(), 0.05];
        b.rod(from, None, end, (len * m as f64).round() as usize, bulge);
    };
    for (deg, lift) in [(150.0, 0.1), (200.0, -0.1), (250.0, 0.05)] {
        arm(&mut b, a, deg, 0.6, lift);
    }
    for (deg, lift) in [(-60.0, 0.1), (20.0, -0.05)] {
        arm(&mut b, bb, deg, 0.6, lift);
    }
    arm(&mut b, c, 95.0, 0.5, 0.1);
    // Crossing rods over midpoints of the A–B rod (vertices 3.. are its interior vertices),
    // alternately above and below, bent in x so that they are not taut.
    let ab: Vec<P> = std::iter::once(b.pos[a as usize])
        .chain((0..m - 1).map(|k| b.pos[3 + k]))
        .chain(std::iter::once(b.pos[bb as usize]))
        .collect();
    for k in 0..crossings {
        let i = (k + 1) * m / (crossings + 1);
        let (p, q) = (ab[i], ab[i + 1]);
        let mid = [
            p[0] + frac * (q[0] - p[0]),
            p[1] + frac * (q[1] - p[1]) + shift * h,
            0.0,
        ];
        let z = if k % 2 == 0 { gap * h } else { -gap * h };
        // 2·half + 1 segments of length ≈ h along y, symmetric about the midpoint.
        let half = 9usize;
        let n = 2 * half + 1;
        let y0 = mid[1] - (half as f64 + 0.5) * h;
        let y1 = mid[1] + (half as f64 + 0.5) * h;
        let s0 = b.vertex([mid[0] - 0.1, y0, z]);
        // bulge 0.1 · sin(π t) brings x back to mid[0] at t = ½ (and by symmetry the edge
        // around t = ½ has its midpoint above `mid`).
        b.rod(s0, None, [mid[0] - 0.1, y1, z], n, [0.1, 0.0, 0.0]);
    }
    let n = b.pos.len();
    (CurveNetwork::new(n, b.edges).unwrap(), b.pos)
}

/// Edge lengths and pinned endpoints (the capillary scene's constraints).
fn pinned_edge_lengths(net: &CurveNetwork, pos: &[P]) -> Vec<Constraint> {
    let mut c = vec![Constraint::EdgeLengths {
        edges: None,
        targets: None,
    }];
    for v in 0..net.n_vertices() as u32 {
        if net.incident_edges(v).len() == 1 {
            c.push(Constraint::Point {
                vertex: v,
                target: Some(pos[v as usize]),
            });
        }
    }
    c
}

#[test]
#[ignore]
fn probe_hub_network() {
    let m: usize = std::env::var("RC_PROBE_M")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(40);
    let k: usize = std::env::var("RC_PROBE_X")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    let gap: f64 = std::env::var("RC_PROBE_GAP")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.05);
    let frac: f64 = std::env::var("RC_PROBE_FRAC")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.5);
    let shift: f64 = std::env::var("RC_PROBE_SHIFT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.0);
    let (net, pos) = hub_network_at(m, k, gap, frac, shift);
    let cons = pinned_edge_lengths(&net, &pos);
    let pr = problem(net, pos, Exponents::default(), &cons);
    describe(&pr);
    let layout = pr.set.layout();
    let mut dense = dense_sobolev_solver(
        pool(16),
        &pr.net,
        &pr.pos,
        pr.ex,
        GramQuadrature::Midpoint,
        &pr.set,
    )
    .unwrap();
    let gd = projected_gradient(&mut dense, layout, &pr.diff).unwrap();
    let xd = layout.gather(&gd.direction);
    for kappa in [0.0, 1.0] {
        let o = MgOptions {
            coarsest_vertices: 40,
            near_field_factor: kappa,
            bct: probe_bct(),
            smoother_factor: probe_smoother(),
            energy_check: std::env::var("RC_PROBE_ENERGY").map_or(true, |v| v == "1"),
            ..MgOptions::default()
        };
        let mut mg = mg_solver(&pr, 16, o).unwrap();
        if std::env::var("RC_PROBE_VERBOSE").is_ok() {
            describe_levels(&mg);
            level_ratios(&mg);
            smooth_ratios(&mg, &pr);
        }
        match projected_gradient(&mut mg, layout, &pr.diff) {
            Ok(g) => {
                let xm = layout.gather(&g.direction);
                eprintln!(
                    "kappa {kappa}: error {:.2e}, levels {}, {} cycles, final residual {:.1e}",
                    rel_err(&xm, &xd),
                    g.linear.levels,
                    g.linear.iterations,
                    mg.last_reports()[0].final_residual(),
                );
            }
            Err(e) => eprintln!("kappa {kappa}: {e}"),
        }
    }
}

/// Closed curves (Freedman–He–Wang) with the network options on and off: cycles, time, error
/// against dense. Env `RC_SIZES` (e.g. `600,2400,2400:1`), `RC_CONFIGS` (`k,s,nu` triples
/// separated by `;`).
#[test]
#[ignore]
fn probe_closed_curves() {
    use ddg_repulsive_curves::constraints::default_constraints;
    use ddg_repulsive_curves::subdivide::subdivide;
    let var = |k: &str, d: &str| std::env::var(k).unwrap_or_else(|_| d.to_string());
    let base = data("FreedmanHeWang");
    if !base.exists() {
        return;
    }
    let ex = Exponents::new(2.0, 4.5).unwrap();
    for size in var("RC_SIZES", "600,2400").split(',') {
        let mut it = size.split(':');
        let b: usize = it.next().unwrap().parse().unwrap();
        let times: usize = it.next().map_or(0, |k| k.parse().unwrap());
        let (mut net, mut pos) =
            read_curve_obj(base.join(format!("freedman-he-wang-{b}.obj"))).unwrap();
        for _ in 0..times {
            let (n, p, _) = subdivide(&net, &pos).unwrap();
            net = n;
            pos = p;
        }
        let pr = problem(net, pos, ex, &default_constraints());
        let layout = pr.set.layout();
        let n = pr.net.n_vertices();
        let xd = (n <= 2400).then(|| {
            let mut dense = dense_sobolev_solver(
                pool(16),
                &pr.net,
                &pr.pos,
                pr.ex,
                GramQuadrature::Midpoint,
                &pr.set,
            )
            .unwrap();
            layout.gather(
                &projected_gradient(&mut dense, layout, &pr.diff)
                    .unwrap()
                    .direction,
            )
        });
        for cfg in var("RC_CONFIGS", "0,0,12;1,0,12;1,1,12;1,1,24").split(';') {
            let v: Vec<f64> = cfg.split(',').map(|x| x.parse().unwrap()).collect();
            let o = MgOptions {
                near_field_factor: v[0],
                smoother_factor: v[1],
                smoothing_steps: v[2] as usize,
                initial_smoothing_steps: 2 * v[2] as usize,
                ..MgOptions::default()
            };
            let t = Instant::now();
            let mut mg = mg_solver(&pr, 16, o).unwrap();
            let tb = t.elapsed().as_secs_f64();
            let t = Instant::now();
            let g = projected_gradient(&mut mg, layout, &pr.diff).unwrap();
            let ts = t.elapsed().as_secs_f64();
            let err = xd
                .as_ref()
                .map_or(f64::NAN, |xd| rel_err(&layout.gather(&g.direction), xd));
            eprintln!(
                "freedman {size} ({n}): kappa {} smoother {} nu {}: levels {}, {} cycles, build {tb:.2} s, solve {ts:.2} s, error {err:.2e}",
                v[0], v[1], v[2], g.linear.levels, g.linear.iterations
            );
        }
    }
}

// ---------------------------------------------------------------- tests

/// The pair terms of `near` (Eq. 31 per unordered pair) sum to the dense metric
/// (`rc-kernel-matrix-form`), for both quadratures.
#[test]
fn pair_blocks_sum_to_the_dense_metric() {
    use ddg_repulsive_curves::metric::{GramGeometry, sobolev_gram};
    use ddg_repulsive_curves::multigrid::near::pair_block;
    let (net, pos) = hub_network(8, 2, 0.3);
    let ex = Exponents::default();
    let n = net.n_vertices();
    let edges = net.edges();
    for q in [GramQuadrature::Midpoint, GramQuadrature::Trapezoid] {
        let g = sobolev_gram(&pool(4), &net, &pos, ex, q);
        let geom = GramGeometry::new(&net, &pos, ex);
        let mut a = vec![0.0; n * n];
        for i in 0..edges.len() {
            for j in i + 1..edges.len() {
                let b = pair_block(&net, &geom, i, j, q);
                let v = [edges[i][0], edges[i][1], edges[j][0], edges[j][1]];
                for r in 0..4 {
                    for c in 0..4 {
                        a[v[r] as usize * n + v[c] as usize] += b[r][c];
                    }
                }
            }
        }
        let scale = g.a.max_abs();
        let err = (0..n * n)
            .map(|k| (a[k] - g.a.get(k / n, k % n)).abs())
            .fold(0.0, f64::max);
        assert!(err <= 1e-12 * scale, "{q:?}: {err:e} vs {scale:e}");
    }
}

/// With every coarse pair corrected (a huge near-field factor), the coarse metric is exactly
/// the Galerkin product `Jᵀ A J` of the fine metric (`rc-prolongation`).
#[test]
fn full_correction_gives_the_galerkin_coarse_metric() {
    use ddg_repulsive_curves::metric::sobolev_gram;
    let (net, pos) = hub_network(10, 3, 0.2);
    let cons = pinned_edge_lengths(&net, &pos);
    let pr = problem(net, pos, Exponents::default(), &cons);
    let o = MgOptions {
        coarsest_vertices: 1,
        max_levels: 2,
        max_fixed_fraction: 1.0,
        near_field_factor: 1e6,
        smoother_factor: 0.0,
        ..MgOptions::default()
    };
    let mg = mg_solver(&pr, 4, o).unwrap();
    let h = mg.hierarchy();
    assert_eq!(h.n_levels(), 2);
    let cg = h.levels()[0].coarsening().unwrap();
    let (nf, nc) = (cg.n_fine(), cg.n_coarse());
    let af = sobolev_gram(&pool(4), &pr.net, &pr.pos, pr.ex, GramQuadrature::Midpoint);
    let j = cg.prolongation_matrix().unwrap();
    // Jᵀ A J column by column, against the coarse level's metric (x slot of each vertex).
    let mut worst: f64 = 0.0;
    let mut scale: f64 = 0.0;
    for c in 0..nc {
        let mut e = vec![0.0; nc];
        e[c] = 1.0;
        let je = mul(&j, &e);
        let aje: Vec<f64> = (0..nf)
            .map(|u| (0..nf).map(|v| af.a.get(u, v) * je[v]).sum())
            .collect();
        let gal = mul(&j.transpose(), &aje);
        let mut x = vec![0.0; 3 * nc];
        x[3 * c] = 1.0;
        let y = h.apply_metric(1, &x).unwrap();
        for u in 0..nc {
            worst = worst.max((y[3 * u] - gal[u]).abs());
            scale = scale.max(gal[u].abs());
        }
    }
    assert!(worst <= 1e-10 * scale, "{worst:e} vs {scale:e}");
}

/// The network of the capillary scene in small: junctions of degree 3, 4, 5, pinned endpoints,
/// all edge lengths fixed, and branches crossing at a tenth of an edge length. The multigrid
/// gradient agrees with the dense one with a small cycle count.
#[test]
fn multigrid_matches_dense_on_a_network_with_junctions_crossings_and_pins() {
    for (frac, shift) in [(0.1, 0.2), (0.5, 0.0)] {
        let (net, pos) = hub_network_at(40, 12, 0.1, frac, shift);
        let mut hist = [0usize; 6];
        for v in 0..net.n_vertices() as u32 {
            hist[net.incident_edges(v).len().min(5)] += 1;
        }
        assert!(hist[3] > 0 && hist[4] > 0 && hist[5] > 0, "{hist:?}");
        let cons = pinned_edge_lengths(&net, &pos);
        let pr = problem(net, pos, Exponents::default(), &cons);
        let layout = pr.set.layout();
        let mut dense = dense_sobolev_solver(
            pool(8),
            &pr.net,
            &pr.pos,
            pr.ex,
            GramQuadrature::Midpoint,
            &pr.set,
        )
        .unwrap();
        let xd = layout.gather(
            &projected_gradient(&mut dense, layout, &pr.diff)
                .unwrap()
                .direction,
        );
        // Exact products (θ = 0), so that the error is the solver's. The default tolerance
        // bounds the residual; the Euclidean error of the direction is a few times larger
        // (the near-contacts make the metric ill conditioned); a tight tolerance shows that
        // the iteration converges to the dense solution.
        let runs: &[(f64, f64, usize)] = if frac == 0.5 {
            &[(1e-3, 2e-2, 15)]
        } else {
            &[(1e-3, 2e-2, 15), (1e-8, 1e-6, 40)]
        };
        for &(tol, bound, max_cycles) in runs {
            let o = MgOptions {
                coarsest_vertices: 40,
                tolerance: tol,
                failure_residual: 10.0 * tol,
                max_cycles,
                bct: BctOptions {
                    theta_x: 0.0,
                    theta_t: 0.0,
                    ..BctOptions::default()
                },
                ..MgOptions::default()
            };
            let mut mg = mg_solver(&pr, 8, o).unwrap();
            let g = projected_gradient(&mut mg, layout, &pr.diff).unwrap();
            let xm = layout.gather(&g.direction);
            let err = rel_err(&xm, &xd);
            let rep = &mg.last_reports()[0];
            eprintln!(
                "crossings at {frac}/{shift}, tolerance {tol:e}: {} levels, {} cycles, residuals {:?}, error {err:.2e}",
                g.linear.levels, g.linear.iterations, rep.residuals
            );
            assert!(g.linear.levels >= 3);
            assert!(rep.converged);
            assert!(g.linear.iterations <= max_cycles, "{}", g.linear.iterations);
            assert!(err <= bound, "{err:e}");
            assert!(norm(&mul(&pr.c, &xm)) <= 1e-10 * norm(&xm));
        }
    }
}

/// A solve that cannot reach its tolerance within the budget and ends above
/// `failure_residual` is an error that names it, not a silent result.
#[test]
fn unconverged_solve_above_the_failure_bound_is_an_error() {
    let (net, pos) = hub_network(20, 4, 0.2);
    let cons = pinned_edge_lengths(&net, &pos);
    let pr = problem(net, pos, Exponents::default(), &cons);
    let o = MgOptions {
        coarsest_vertices: 20,
        tolerance: 1e-14,
        failure_residual: 1e-13,
        max_cycles: 1,
        ..MgOptions::default()
    };
    let mut mg = mg_solver(&pr, 2, o).unwrap();
    let r = projected_gradient(&mut mg, pr.set.layout(), &pr.diff);
    assert!(
        matches!(r, Err(RcError::MultigridDiverged { .. })),
        "{:?}",
        r.map(|g| g.linear)
    );
    // Bad options are errors.
    for o in [
        MgOptions {
            near_field_factor: -1.0,
            ..MgOptions::default()
        },
        MgOptions {
            smoother_factor: f64::NAN,
            ..MgOptions::default()
        },
        MgOptions {
            max_fixed_fraction: 0.0,
            ..MgOptions::default()
        },
        MgOptions {
            failure_residual: 1e-9,
            ..MgOptions::default()
        },
    ] {
        assert!(matches!(o.validate(), Err(RcError::InvalidOption { .. })));
    }
}

/// The solve is bitwise identical at 1, 3 and 16 threads on a network with near-contacts (D25).
#[test]
fn network_solve_is_bitwise_identical_across_threads() {
    let (net, pos) = hub_network_at(30, 6, 0.1, 0.1, 0.2);
    let cons = pinned_edge_lengths(&net, &pos);
    let pr = problem(net, pos, Exponents::default(), &cons);
    let o = MgOptions {
        coarsest_vertices: 30,
        ..MgOptions::default()
    };
    let mut first: Option<Vec<f64>> = None;
    for t in [1, 3, 16] {
        let mut mg = mg_solver(&pr, t, o).unwrap();
        let g = projected_gradient(&mut mg, pr.set.layout(), &pr.diff).unwrap();
        let x = pr.set.layout().gather(&g.direction);
        match &first {
            None => first = Some(x),
            Some(f) => assert!(f.iter().zip(&x).all(|(a, b)| a.to_bits() == b.to_bits())),
        }
    }
}

/// The authors' CapillaryNetwork scene (both curves; data from
/// `showcase/repulsive-curves/fetch_scenes.sh`): the multigrid gradient against the dense one,
/// then `RC_FLOW_STEPS` steps (default 20 on the 2146-vertex curve, 3 on the 4296-vertex one)
/// of the exact-energy flow with each solver. Release only; slow (the dense solve of the large
/// curve takes about a minute per step).
/// `cargo test --release --test multigrid_networks -- --ignored --nocapture capillary_scene`.
#[test]
#[ignore]
fn capillary_scene_multigrid_matches_dense() {
    use ddg_repulsive_curves::scene::BackendChoice;
    for high in [false, true] {
        let Some(scene) = capillary(high) else {
            return;
        };
        let pr = problem(
            scene.network.clone(),
            scene.positions.clone(),
            scene.spec.exponents,
            &scene.constraints(),
        );
        let layout = pr.set.layout();
        let mut mg = mg_solver(&pr, 32, MgOptions::default()).unwrap();
        let t = Instant::now();
        let gm = projected_gradient(&mut mg, layout, &pr.diff).unwrap();
        let tm = t.elapsed().as_secs_f64();
        let mut dense = dense_sobolev_solver(
            pool(32),
            &pr.net,
            &pr.pos,
            pr.ex,
            GramQuadrature::Midpoint,
            &pr.set,
        )
        .unwrap();
        let gd = projected_gradient(&mut dense, layout, &pr.diff).unwrap();
        let err = rel_err(&layout.gather(&gm.direction), &layout.gather(&gd.direction));
        let rep = &mg.last_reports()[0];
        eprintln!(
            "{} vertices: levels {:?}, {} cycles, residuals {:?}, solve {tm:.1} s, error {err:.2e}",
            pr.net.n_vertices(),
            mg.hierarchy()
                .levels()
                .iter()
                .map(|l| l.network().n_vertices())
                .collect::<Vec<_>>(),
            gm.linear.iterations,
            rep.residuals
        );
        assert!(rep.converged);
        assert!(gm.linear.iterations <= 12);
        assert!(err <= 1e-2, "{err:e}");

        let steps: usize = std::env::var("RC_FLOW_STEPS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(if high { 3 } else { 20 });
        let mut sd = scene.session(BackendChoice::Exact, Some(32)).unwrap();
        let mut sm = scene
            .session(BackendChoice::ExactMultigrid, Some(32))
            .unwrap();
        let l0 = pr.net.mean_edge_length(&pr.pos);
        for k in 1..=steps {
            let rd = sd.step().unwrap();
            let rm = sm.step().unwrap();
            let de = (rm.energy_after - rd.energy_after).abs() / rd.energy_after.abs();
            let dx = sd
                .positions()
                .iter()
                .zip(sm.positions())
                .map(|(a, b)| {
                    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
                })
                .fold(0.0, f64::max);
            eprintln!(
                "step {k}: energy dense {:.6e} mg {:.6e} (rel diff {de:.1e}), max |x_mg - x_d| = {:.1e} edge lengths, mg {} cycles residual {:.1e}",
                rd.energy_after,
                rm.energy_after,
                dx / l0,
                rm.linear.iterations,
                rm.linear.relative_residual
            );
            assert!(de <= 1e-3, "step {k}: {de:e}");
        }
    }
}
