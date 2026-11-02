//! WP6 measurements (`SPEC.md` §14.3, §14.5; `NOTES-WP6.md`): hierarchical product error,
//! GMRES iterations and solve times with the sparse preconditioner, barycenter modes, the step
//! breakdown, dense vs hierarchical cost. One CSV line per run on stdout.
//!
//! ```text
//! wp6_precond solve <mesh> [--threads N] [--chi X] [--pairs V|F] [--side R|L] [--form A|C]
//!                          [--mode percomp|global] [--eps-multi X] [--schur]
//! wp6_precond error <mesh> [--threads N]            (χ = 1, 0.5, 0.25 against χ = 0)
//! wp6_precond breakdown <mesh> [--threads N]        (one full direction: BVH … GMRES, Schur rows)
//! wp6_precond dense <mesh> [--threads N]            (dense vs hierarchical direction cost)
//! wp6_precond modes <mesh> [--threads N]            (Global vs PerComponent (+ shift) directions)
//! ```
//!
//! `<mesh>`: `icosphere:L`, `torus:N:M`, `stretched:N:M` (a torus stretched ×2.5 in x with bumps),
//! or an OBJ path.

use std::sync::Arc;
use std::time::Instant;

use ddg_core::generators::{icosphere, perturb, torus};
use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::BhOptions;
use ddg_repulsive_surfaces::energy::barnes_hut::BhTree;
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::geometry::{FaceGeometry, mean_edge_length, vertex_areas};
use ddg_repulsive_surfaces::metric::hierarchical::{
    ClusterTreeOptions, FaceFractionalTree, FaceMetricTree, VertexFractionalTree,
};
use ddg_repulsive_surfaces::metric::{FractionalPairs, SobolevWeights};
use ddg_repulsive_surfaces::solve::iterative::{
    HierarchicalOptions, MetricExtras, dense_metric_with_extras, hierarchical_saddle,
};
use ddg_repulsive_surfaces::solve::precond::PreconditionerForm;
use ddg_repulsive_surfaces::solve::{
    DenseSaddle, GmresSide, SaddleSolver, SchurComplement,
};
use ddg_repulsive_surfaces::{Exponents, Surface, read_obj};

fn load(spec: &str) -> (String, Surface, Vec<[f64; 3]>) {
    let parts: Vec<&str> = spec.split(':').collect();
    let (name, s, mut p) = match parts[0] {
        "icosphere" => {
            let m = icosphere(parts[1].parse().unwrap()).unwrap();
            (spec.to_string(), Surface::new(m.positions.len(), m.faces).unwrap(), m.positions)
        }
        "torus" | "stretched" => {
            let (n, k) = (parts[1].parse().unwrap(), parts[2].parse().unwrap());
            let m = torus(n, k, 1.5, 0.6).unwrap();
            let mut p = m.positions;
            if parts[0] == "stretched" {
                for q in &mut p {
                    q[2] += 0.12 * (4.0 * q[0]).sin() * (3.0 * q[1]).cos();
                    q[0] *= 2.5;
                }
            }
            (spec.to_string(), Surface::new(p.len(), m.faces).unwrap(), p)
        }
        _ => {
            let (s, p) = read_obj(spec).unwrap();
            let name = std::path::Path::new(spec).file_stem().unwrap().to_string_lossy().into_owned();
            (name, s, p)
        }
    };
    if parts[0] != "icosphere" && parts[0] != "torus" && parts[0] != "stretched" {
        return (name, s, p);
    }
    let l = mean_edge_length(&s, &p);
    perturb(&mut p, [0.05 * l; 3], 7);
    (name, s, p)
}

struct Args {
    cmd: String,
    mesh: String,
    threads: usize,
    chi: f64,
    pairs: FractionalPairs,
    side: GmresSide,
    form: PreconditionerForm,
    mode: BarycenterMode,
    eps_multi: f64,
    schur: bool,
    willmore: f64,
    wprec: bool,
}

fn args() -> Args {
    let a: Vec<String> = std::env::args().collect();
    let mut out = Args {
        cmd: a[1].clone(),
        mesh: a[2].clone(),
        threads: 4,
        chi: 0.5,
        pairs: FractionalPairs::Vertices,
        side: GmresSide::Right,
        form: PreconditionerForm::Augmented,
        mode: BarycenterMode::PerComponent,
        eps_multi: 1e-2,
        schur: false,
        willmore: 0.0,
        wprec: true,
    };
    let mut i = 3;
    while i < a.len() {
        let v = a.get(i + 1).cloned().unwrap_or_default();
        match a[i].as_str() {
            "--threads" => out.threads = v.parse().unwrap(),
            "--chi" => out.chi = v.parse().unwrap(),
            "--pairs" => out.pairs = if v == "F" { FractionalPairs::Faces } else { FractionalPairs::Vertices },
            "--side" => out.side = if v == "L" { GmresSide::Left } else { GmresSide::Right },
            "--form" => out.form = if v == "C" { PreconditionerForm::Congruence } else { PreconditionerForm::Augmented },
            "--mode" => out.mode = if v == "global" { BarycenterMode::Global } else { BarycenterMode::PerComponent },
            "--eps-multi" => out.eps_multi = v.parse().unwrap(),
            "--willmore" => out.willmore = v.parse().unwrap(),
            "--wprec" => out.wprec = v == "1",
            "--schur" => {
                out.schur = true;
                i += 1;
                continue;
            }
            x => panic!("unknown option {x}"),
        }
        i += 2;
    }
    out
}

fn options(a: &Args) -> HierarchicalOptions {
    let mut o = HierarchicalOptions::default();
    o.cluster.chi = a.chi;
    o.preconditioner.fractional = a.pairs;
    o.preconditioner.form = a.form;
    o.preconditioner.eps_multi_rel = a.eps_multi;
    o.gmres.side = a.side;
    o.willmore_in_preconditioner = a.wprec;
    o
}

fn rel(a: &[f64], b: &[f64]) -> f64 {
    let d: f64 = a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum::<f64>().sqrt();
    d / b.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn flat(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

fn gradient(pool: &rayon::ThreadPool, s: &Surface, p: &[[f64; 3]]) -> (Vec<f64>, f64, f64, BhTree) {
    let t0 = Instant::now();
    let bh = BhTree::build(pool, s, p, BhOptions::default()).unwrap();
    let t_build = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let ev = bh.energy_and_differential(pool, &Exponents::default()).unwrap();
    (flat(&ev.differential), t_build, t1.elapsed().as_secs_f64(), bh)
}

fn constraints(schur: bool) -> Vec<Constraint> {
    let mut c = vec![Constraint::Barycenter { component: None, target: None }];
    if schur {
        c.push(Constraint::Area { component: None, target: None });
        c.push(Constraint::Volume { component: None, target: None });
    }
    c
}

fn solve(a: &Args) {
    let (name, s, p) = load(&a.mesh);
    let pool = Arc::new(make_pool(Some(a.threads)).unwrap());
    let e = Exponents::default();
    let (g, _, _, _) = gradient(&pool, &s, &p);
    let set = ConstraintSet::new(&s, &p, &constraints(a.schur), a.mode).unwrap();
    let q = set.positional_rows();
    let cq = set.jacobian(&p, &q).unwrap();
    let o = options(a);
    let t0 = Instant::now();
    let weights = SobolevWeights { willmore_weight: a.willmore, ..Default::default() };
    let mut it = match hierarchical_saddle(pool.clone(), &s, &p, &e, &weights, &MetricExtras::default(), &cq, &o, None) {
        Ok(x) => x,
        Err(err) => {
            println!("solve,{name},{},{},error,{err}", s.n_faces(), s.n_vertices());
            return;
        }
    };
    let t_build = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let nl = set.nonlinear_rows();
    let (iters, res, conv) = if nl.is_empty() {
        match it.solve(&g, &vec![0.0; cq.nrows()]) {
            Ok(sol) => (vec![sol.record.iterations], sol.record.relative_residual, sol.record.converged),
            Err(err) => {
                println!("solve,{name},{},{},diverged,{err}", s.n_faces(), s.n_vertices());
                return;
            }
        }
    } else {
        let cn = set.jacobian(&p, &nl).unwrap();
        let schur = SchurComplement::new(&mut it, nl, cn).unwrap();
        let dir = schur.direction(&mut it, &g).unwrap();
        let mut its: Vec<usize> = schur.records().iter().map(|r| r.iterations).collect();
        its.push(dir.record.iterations);
        let res = schur.records().iter().map(|r| r.relative_residual).fold(dir.record.relative_residual, f64::max);
        (its, res, dir.record.converged)
    };
    let t_solve = t1.elapsed().as_secs_f64();
    let (nop, top) = it.operator().usage();
    let (npc, tpc) = it.preconditioner().usage();
    let (eps, multi) = it.preconditioner().epsilon();
    // mesh,faces,vertices,components,chi,pairs,side,form,mode,multi,eps,iterations,residual,converged,build_s,solve_s,products,product_s,precond,precond_s,near_entries
    println!(
        "solve,{name},{},{},{},{},{:?},{:?},{:?},{:?},{},{:.3e},{},{:.2e},{},{:.3},{:.3},{},{:.3},{},{:.3},{},{},{}",
        s.n_faces(),
        s.n_vertices(),
        s.n_components(),
        a.chi,
        a.pairs,
        a.side,
        a.form,
        a.mode,
        multi,
        eps,
        iters.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("+"),
        res,
        conv,
        t_build,
        t_solve,
        nop,
        top,
        npc,
        tpc,
        it.operator().tree().stats().near_entries,
        a.willmore,
        a.wprec,
    );
}

fn error(a: &Args) {
    let (name, s, p) = load(&a.mesh);
    let pool = make_pool(Some(a.threads)).unwrap();
    let e = Exponents::default();
    let n = 3 * s.n_vertices();
    let smooth: Vec<f64> = p.iter().flat_map(|q| [q[0].sin() + q[2], q[1] * q[0], (2.0 * q[2]).cos()]).collect();
    let mut r = ddg_core::generators::SplitMix64::new(3);
    let random: Vec<f64> = (0..n).map(|_| r.uniform(-1.0, 1.0)).collect();
    let geom = FaceGeometry::new(&s, &p).unwrap();
    let mass = vertex_areas(&s, &geom.areas);
    let tau = 2.0 - e.s();
    let apply_all = |chi: f64, x: &[f64]| -> ([Vec<f64>; 4], [usize; 3], f64) {
        let o = ClusterTreeOptions { chi, ..Default::default() };
        let t0 = Instant::now();
        let tree = FaceMetricTree::build(&pool, &s, &p, &e, &o).unwrap();
        let ft = FaceFractionalTree::build(&pool, tree.face_tree(), tree.operators(), tau, &o).unwrap();
        let vt = VertexFractionalTree::build(&pool, &p, &mass, tau, &o).unwrap();
        let tb = t0.elapsed().as_secs_f64();
        let mut out: [Vec<f64>; 4] = std::array::from_fn(|_| vec![0.0; n]);
        let [o0, o1, o2, o3] = &mut out;
        tree.apply_high(&pool, x, o0).unwrap();
        tree.apply_low(&pool, x, o1).unwrap();
        ft.apply(&pool, x, o2).unwrap();
        vt.apply(&pool, x, o3).unwrap();
        (out, [tree.stats().near_entries, ft.stats().near_entries, vt.stats().near_entries], tb)
    };
    for (kind, x) in [("random", &random), ("smooth", &smooth)] {
        let (exact, _, _) = apply_all(0.0, x);
        for chi in [1.0, 0.5, 0.25] {
            let (approx, near, tb) = apply_all(chi, x);
            println!(
                "error,{name},{},{kind},{chi},{:.2e},{:.2e},{:.2e},{:.2e},{},{},{},{:.3}",
                s.n_faces(),
                rel(&approx[0], &exact[0]),
                rel(&approx[1], &exact[1]),
                rel(&approx[2], &exact[2]),
                rel(&approx[3], &exact[3]),
                near[0],
                near[1],
                near[2],
                tb
            );
        }
    }
}

fn breakdown(a: &Args) {
    let (name, s, p) = load(&a.mesh);
    let pool = Arc::new(make_pool(Some(a.threads)).unwrap());
    let e = Exponents::default();
    let w = SobolevWeights::default();
    for rep in 0..3 {
        let (g, t_bh_build, t_bh_eval, bh) = gradient(&pool, &s, &p);
        let set = ConstraintSet::new(&s, &p, &constraints(true), BarycenterMode::PerComponent).unwrap();
        let q = set.positional_rows();
        let nl = set.nonlinear_rows();
        let cq = set.jacobian(&p, &q).unwrap();
        let cn = set.jacobian(&p, &nl).unwrap();
        let mut o = options(a);
        o.cluster.leaf_size = BhOptions::default().leaf_size;
        let t0 = Instant::now();
        let mut it = hierarchical_saddle(pool.clone(), &s, &p, &e, &w, &MetricExtras::default(), &cq, &o, Some(bh.face_tree().clone())).unwrap();
        let t_setup = t0.elapsed().as_secs_f64();
        let bct_s = it.operator().tree().build_seconds();
        let factor_s = it.preconditioner().factor_seconds();
        let frac_s = it.preconditioner().fractional_seconds();
        let t1 = Instant::now();
        let schur = SchurComplement::new(&mut it, nl.clone(), cn).unwrap();
        let dir = schur.direction(&mut it, &g).unwrap();
        let t_solve = t1.elapsed().as_secs_f64();
        let mut its: Vec<usize> = schur.records().iter().map(|r| r.iterations).collect();
        its.push(dir.record.iterations);
        let (nop, top) = it.operator().usage();
        let (npc, tpc) = it.preconditioner().usage();
        // name,faces,threads,rep,bh_build,bh_energy_grad,setup(face ops+bct+factor+frac),bct,factor,frac,solve_total,iterations,products,product_s,precond,precond_s,total
        println!(
            "breakdown,{name},{},{},{rep},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{},{},{:.3},{},{:.3},{:.3}",
            s.n_faces(),
            a.threads,
            t_bh_build,
            t_bh_eval,
            t_setup,
            bct_s,
            factor_s,
            frac_s,
            t_solve,
            its.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("+"),
            nop,
            top,
            npc,
            tpc,
            t_bh_build + t_bh_eval + t_setup + t_solve
        );
    }
}

fn dense(a: &Args) {
    let (name, s, p) = load(&a.mesh);
    let pool = Arc::new(make_pool(Some(a.threads)).unwrap());
    let e = Exponents::default();
    let w = SobolevWeights::default();
    let (g, _, _, _) = gradient(&pool, &s, &p);
    let set = ConstraintSet::new(&s, &p, &constraints(false), BarycenterMode::PerComponent).unwrap();
    let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
    let b = vec![0.0; cq.nrows()];
    let t0 = Instant::now();
    let am = dense_metric_with_extras(&pool, &s, &p, &e, &w, &MetricExtras::default()).unwrap();
    let t_assemble = t0.elapsed().as_secs_f64();
    let mut d = DenseSaddle::new(&pool, &am, &cq).unwrap();
    let xd = d.solve(&g, &b).unwrap().x;
    let t_dense = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let mut h = hierarchical_saddle(pool.clone(), &s, &p, &e, &w, &MetricExtras::default(), &cq, &options(a), None).unwrap();
    let sol = h.solve(&g, &b).unwrap();
    let t_h = t1.elapsed().as_secs_f64();
    println!(
        "dense,{name},{},{},{:.3},{:.3},{:.3},{},{:.2e}",
        s.n_faces(),
        a.threads,
        t_assemble,
        t_dense,
        t_h,
        sol.record.iterations,
        rel(&sol.x, &xd)
    );
}

fn modes(a: &Args) {
    let (name, s, p) = load(&a.mesh);
    let pool = Arc::new(make_pool(Some(a.threads)).unwrap());
    let e = Exponents::default();
    let w = SobolevWeights::default();
    let (g, _, _, _) = gradient(&pool, &s, &p);
    let geom = FaceGeometry::new(&s, &p).unwrap();
    let mass = vertex_areas(&s, &geom.areas);
    let mut dirs = Vec::new();
    for mode in [BarycenterMode::Global, BarycenterMode::PerComponent] {
        let set = ConstraintSet::new(&s, &p, &constraints(false), mode).unwrap();
        let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
        let mut o = options(a);
        o.gmres.tolerance = 1e-6;
        o.gmres.max_iterations = 400;
        let mut h = hierarchical_saddle(pool.clone(), &s, &p, &e, &w, &MetricExtras::default(), &cq, &o, None).unwrap();
        let sol = h.solve(&g, &vec![0.0; cq.nrows()]).unwrap();
        let mut o4 = options(a);
        o4.gmres.tolerance = 1e-4;
        let mut h4 = hierarchical_saddle(pool.clone(), &s, &p, &e, &w, &MetricExtras::default(), &cq, &o4, None).unwrap();
        let it4 = h4.solve(&g, &vec![0.0; cq.nrows()]).unwrap().record.iterations;
        let (eps, multi) = h.preconditioner().epsilon();
        eprintln!("{name} {mode:?}: eps {eps:.2e} multi {multi}, iterations to 1e-4 {it4}, to 1e-6 {}", sol.record.iterations);
        dirs.push((sol.x, it4, sol.record.iterations));
    }
    // The per-component L² mean-motion shift (Sec. 7.1.1): add Σ_{v∈c} g_v / A_c on component c.
    let vc = s.vertex_components();
    let nc = s.n_components();
    let mut sum = vec![[0.0; 3]; nc];
    let mut area = vec![0.0; nc];
    for v in 0..s.n_vertices() {
        let c = vc[v] as usize;
        area[c] += mass[v];
        for k in 0..3 {
            sum[c][k] += g[3 * v + k];
        }
    }
    let mut shifted = dirs[1].0.clone();
    for v in 0..s.n_vertices() {
        let c = vc[v] as usize;
        for k in 0..3 {
            shifted[3 * v + k] += sum[c][k] / area[c];
        }
    }
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(a, b)| a * b).sum::<f64>();
    let nrm = |x: &[f64]| dot(x, x).sqrt();
    let xg = &dirs[0].0;
    let net: f64 = sum.iter().map(|s| (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt()).fold(0.0, f64::max);
    // name,faces,components,it_global(1e-4),it_percomp(1e-4),|xg - xpc|/|xg|,|xg - xpc_shift|/|xg|,cos(xg, xpc_shift),max component net force / |g|
    println!(
        "modes,{name},{},{},{},{},{:.3e},{:.3e},{:.6},{:.3e}",
        s.n_faces(),
        nc,
        dirs[0].1,
        dirs[1].1,
        rel(&dirs[1].0, xg),
        rel(&shifted, xg),
        dot(xg, &shifted) / (nrm(xg) * nrm(&shifted)),
        net / nrm(&g)
    );
}

fn main() {
    let a = args();
    match a.cmd.as_str() {
        "solve" => solve(&a),
        "error" => error(&a),
        "breakdown" => breakdown(&a),
        "dense" => dense(&a),
        "modes" => modes(&a),
        c => panic!("unknown command {c}"),
    }
}
