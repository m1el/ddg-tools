//! Diagnostics of the hierarchical saddle solve (`NOTES-SOLVERFIX.md`): the spectrum of the
//! preconditioned saddle operator, GMRES residual histories and the hierarchical-vs-dense product.
//!
//! ```text
//! solverfix spectrum <mesh> [--form A|C] [--eps X] [--chi X]
//! solverfix gmres    <mesh> [--form A|C] [--eps X] [--chi X] [--restart N] [--max N] [--schur]
//! solverfix product  <mesh>
//! ```
//! `<mesh>`: `h:L:pins|bary` (the H of Fig. 22 at subdivision level L with its 24 corner pins or a
//! barycenter), `torus:N:M:pins|bary` (3 pins), `spindle:NM:NN:R` (the subcritical spindle torus,
//! barycenter, exponents (3.5, 7), `low_order: false`), `cobble:N:L` (N spheres, global barycenter,
//! per-component areas), `obj:PATH` (barycenter).
#![allow(clippy::too_many_lines)]

#[path = "apps/apps.rs"]
#[allow(unused_imports)]
mod apps;

use std::sync::Arc;

use ddg_core::generators::torus;
use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::BhOptions;
use ddg_repulsive_surfaces::energy::barnes_hut::BhTree;
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::generate::primitives;
use ddg_repulsive_surfaces::metric::SobolevWeights;
use ddg_repulsive_surfaces::solve::gmres::{LinearOperator, Preconditioner, gmres};
use ddg_repulsive_surfaces::solve::iterative::{
    HierarchicalOptions, MetricExtras, dense_metric_with_extras, hierarchical_saddle,
};
use ddg_repulsive_surfaces::solve::precond::{FractionalScale, PreconditionerForm, SparsePreconditioner};
use ddg_repulsive_surfaces::solve::{DenseSaddleOperator, SchurComplement};
use ddg_repulsive_surfaces::{Exponents, Surface, read_obj};

type P = [f64; 3];

struct Case {
    surface: Surface,
    positions: Vec<P>,
    constraints: Vec<Constraint>,
    mode: BarycenterMode,
    exponents: Exponents,
    weights: SobolevWeights,
}

fn case(spec: &str) -> Case {
    let p: Vec<&str> = spec.split(':').collect();
    let bary = || Constraint::Barycenter { component: None, target: None };
    let mut c = Case {
        surface: Surface::new(3, vec![[0, 1, 2]]).unwrap(),
        positions: vec![],
        constraints: vec![],
        mode: BarycenterMode::Global,
        exponents: Exponents::default(),
        weights: SobolevWeights::default(),
    };
    match p[0] {
        "h" => {
            let (m, _) = apps::h_mesh(p[1].parse().unwrap());
            let (s, x) = apps::mesh_pair(&m).unwrap();
            c.surface = s;
            c.positions = x;
            c.constraints.push(if p[2] == "pins" {
                Constraint::Pins { vertices: (0..24).collect(), targets: None }
            } else {
                bary()
            });
            c.constraints.push(Constraint::Volume { component: None, target: None });
        }
        "torus" => {
            let m = torus(p[1].parse().unwrap(), p[2].parse().unwrap(), 1.5, 0.6).unwrap();
            c.surface = Surface::new(m.positions.len(), m.faces).unwrap();
            c.positions = m.positions;
            c.constraints.push(if p[3] == "pins" {
                let n = c.positions.len() as u32;
                Constraint::Pins { vertices: vec![0, n / 3, 2 * n / 3], targets: None }
            } else {
                bary()
            });
            c.constraints.push(Constraint::Volume { component: None, target: None });
        }
        "spindle" => {
            let (nm, nn): (usize, usize) = (p[1].parse().unwrap(), p[2].parse().unwrap());
            let r: f64 = p[3].parse().unwrap();
            let mut x = Vec::new();
            for i in 0..nm {
                for j in 0..nn {
                    let (u, v) = (std::f64::consts::TAU * i as f64 / nm as f64, std::f64::consts::TAU * j as f64 / nn as f64);
                    x.push([(1.0 + r * v.cos()) * u.cos(), (1.0 + r * v.cos()) * u.sin(), r * v.sin()]);
                }
            }
            let id = |i: usize, j: usize| ((i % nm) * nn + (j % nn)) as u32;
            let mut f = Vec::new();
            for i in 0..nm {
                for j in 0..nn {
                    f.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
                    f.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
                }
            }
            c.surface = Surface::new(x.len(), f).unwrap();
            c.positions = x;
            c.constraints = vec![bary(), Constraint::Area { component: None, target: None }];
            c.exponents = Exponents::new(3.5, 7.0).unwrap();
            c.weights.low_order = false;
            c.mode = BarycenterMode::PerComponent;
        }
        "cobble" => {
            let n: usize = p[1].parse().unwrap();
            let level: u32 = p[2].parse().unwrap();
            let cb = primitives::cobblestones(n, 0.5, 0.15, [3.2, 3.2, 0.9], level, 11).unwrap();
            let (s, x) = apps::mesh_pair(&cb.spheres).unwrap();
            c.constraints.push(bary());
            for k in 0..s.n_components() {
                c.constraints.push(Constraint::Area { component: Some(k as u32), target: None });
            }
            if p.get(3) == Some(&"pc") {
                c.mode = BarycenterMode::PerComponent;
            }
            c.surface = s;
            c.positions = x;
        }
        "obj" => {
            let (s, x) = read_obj(p[1]).unwrap();
            c.surface = s;
            c.positions = x;
            c.constraints.push(bary());
        }
        _ => panic!("unknown mesh {spec}"),
    }
    c
}

struct Args {
    cmd: String,
    mesh: String,
    form: PreconditionerForm,
    eps: Option<f64>,
    chi: f64,
    restart: usize,
    max: usize,
    schur: bool,
    scale: FractionalScale,
    coarse: bool,
}

fn args() -> Args {
    let a: Vec<String> = std::env::args().collect();
    let mut o = Args {
        cmd: a[1].clone(),
        mesh: a[2].clone(),
        form: PreconditionerForm::Augmented,
        eps: None,
        chi: 0.5,
        restart: 30,
        max: 200,
        schur: false,
        scale: FractionalScale::Symbol,
        coarse: true,
    };
    let mut i = 3;
    while i < a.len() {
        let v = a.get(i + 1).cloned().unwrap_or_default();
        match a[i].as_str() {
            "--form" => o.form = if v == "C" { PreconditionerForm::Congruence } else { PreconditionerForm::Augmented },
            "--eps" => o.eps = Some(v.parse().unwrap()),
            "--chi" => o.chi = v.parse().unwrap(),
            "--restart" => o.restart = v.parse().unwrap(),
            "--max" => o.max = v.parse().unwrap(),
            "--scale" => o.scale = if v == "U" { FractionalScale::Unit } else { FractionalScale::Symbol },
            "--nocoarse" => {
                o.coarse = false;
                i += 1;
                continue;
            }
            "--schur" => {
                o.schur = true;
                i += 1;
                continue;
            }
            x => panic!("unknown option {x}"),
        }
        i += 2;
    }
    o
}

fn hopts(a: &Args) -> HierarchicalOptions {
    let mut o = HierarchicalOptions::default();
    o.cluster.chi = a.chi;
    o.preconditioner.form = a.form;
    o.preconditioner.scale = a.scale;
    o.preconditioner.coarse = a.coarse;
    if let Some(e) = a.eps {
        o.preconditioner.eps_rel = e;
        o.preconditioner.eps_multi_rel = e;
    }
    o.gmres.restart = a.restart;
    o.gmres.max_iterations = a.max;
    o
}

fn flat(v: &[P]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

/// The cobblestones flow of `apps` (14 spheres, slab, area x2.5 over 300 steps, dense, no
/// remeshing) with one line per step: the CCD cap, the projection guard and the slab clearance.
fn cobflow(steps: usize) {
    use ddg_repulsive_surfaces::step::{FlowOptions, Schedule, Session, StopRule};
    let c = primitives::cobblestones(14, 0.5, 0.15, [3.2, 3.2, 0.9], 2, 11).unwrap();
    let (s, x) = apps::mesh_pair(&c.spheres).unwrap();
    let mut o = FlowOptions { max_iterations: steps, stop: StopRule::Off, ..FlowOptions::default() };
    o.constraints = vec![Constraint::Barycenter { component: None, target: None }];
    for k in 0..s.n_components() {
        o.constraints.push(Constraint::Area { component: Some(k as u32), target: None });
        o.schedules.push(Schedule { constraint: 1 + k, factor: 2.5, offset: 0.0, steps: 300 });
    }
    let (pen, ob) = apps::mesh_obstacle(&c.slab, 1.0).unwrap();
    o.penalties.push(pen);
    o.ccd_obstacles.push(ob);
    o.backend = ddg_repulsive_surfaces::Backend::Dense;
    o.remesh = None;
    o.threads = Some(8);
    let mut sess = Session::new(s, x, o).unwrap();
    for _ in 0..steps {
        let r = match sess.step() {
            Ok(r) => r,
            Err(e) => {
                println!("error {e}");
                break;
            }
        };
        let zmax = sess.positions().iter().map(|p| p[2].abs()).fold(0.0, f64::max);
        let ccd = r.ccd.as_ref();
        println!(
            "step {} {:?} slope {:.3e} tau0 {:.3e} tau {:.3e} tau_ccd {:?} close {} touching {} contact {:?} blocked {} frac {:.3} proj {:?} viol {:.3e} zmax {:.4}",
            r.iteration,
            r.status,
            r.sobolev_norm_sq,
            r.tau0,
            r.tau,
            r.tau_ccd,
            ccd.map_or(0, |c| c.close_pairs),
            ccd.map_or(0, |c| c.touching_pairs),
            ccd.and_then(|c| c.contact.as_ref().map(|k| k.primitive)),
            r.projection_blocked,
            r.projection_fraction,
            r.projection_status,
            r.max_violation,
            zmax
        );
    }
}

fn main() {
    let a = args();
    if a.cmd == "cobflow" {
        cobflow(a.max);
        return;
    }
    let c = case(&a.mesh);
    let pool = Arc::new(make_pool(Some(8)).unwrap());
    let set = ConstraintSet::new(&c.surface, &c.positions, &c.constraints, c.mode).unwrap();
    let q = set.positional_rows();
    let cq = set.jacobian(&c.positions, &q).unwrap();
    let nl = set.nonlinear_rows();
    let cn = set.jacobian(&c.positions, &nl).unwrap();
    let bh = BhTree::build(&pool, &c.surface, &c.positions, BhOptions::default()).unwrap();
    let g = flat(&bh.energy_and_differential(&pool, &c.exponents).unwrap().differential);
    let ex = MetricExtras::default();
    println!(
        "mesh {} faces {} vertices {} components {} k_Q {} k_N {} s {:.4}",
        a.mesh,
        c.surface.n_faces(),
        c.surface.n_vertices(),
        c.surface.n_components(),
        cq.nrows(),
        nl.len(),
        c.exponents.s()
    );
    match a.cmd.as_str() {
        "spectrum" => {
            let am = dense_metric_with_extras(&pool, &c.surface, &c.positions, &c.exponents, &c.weights, &ex).unwrap();
            let op = DenseSaddleOperator::new(&am, &cq).unwrap();
            let mut o = hopts(&a);
            o.preconditioner.cluster.chi = 0.0;
            o.preconditioner.cluster = o.cluster;
            let mut prec = SparsePreconditioner::new(
                pool.clone(),
                &c.surface,
                &c.positions,
                &c.exponents,
                &cq,
                &o.preconditioner,
                None,
                None,
            )
            .unwrap();
            let attached = prec.attach_coarse(&op, &c.surface, &cq).unwrap();
            println!("eps {:?} scale {:.4e} coarse {attached}", prec.epsilon(), prec.fractional_scale());
            let n = op.dim();
            let mut m = faer::Mat::<f64>::zeros(n, n);
            let mut e = vec![0.0; n];
            let mut t = vec![0.0; n];
            let mut z = vec![0.0; n];
            for j in 0..n {
                e.iter_mut().for_each(|v| *v = 0.0);
                e[j] = 1.0;
                prec.apply(&e, &mut t).unwrap();
                op.apply(&t, &mut z).unwrap();
                for i in 0..n {
                    m[(i, j)] = z[i];
                }
            }
            let ev = m.eigenvalues().unwrap();
            let mut re: Vec<f64> = ev.iter().map(|z| z.re).collect();
            re.sort_by(f64::total_cmp);
            let neg = ev.iter().filter(|z| z.re < 0.0).count();
            let maxim = ev.iter().map(|z| z.im.abs()).fold(0.0, f64::max);
            let absmin = ev.iter().map(|z| z.norm()).fold(f64::MAX, f64::min);
            let absmax = ev.iter().map(|z| z.norm()).fold(0.0, f64::max);
            let qs: Vec<String> = [0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0]
                .iter()
                .map(|f| format!("{:.3e}", re[((n - 1) as f64 * f) as usize]))
                .collect();
            println!(
                "spectrum n {n} negative_re {neg} max|im| {maxim:.3e} |λ| [{absmin:.3e}, {absmax:.3e}] re quantiles {}",
                qs.join(" ")
            );
            let counts: Vec<String> = [3.0, 10.0, 100.0, 1e3, 1e4]
                .iter()
                .map(|t| format!(">{t:e}: {}", ev.iter().filter(|z| z.norm() > *t).count()))
                .collect();
            println!("counts {}", counts.join(", "));
            let mut lst: Vec<_> = ev.iter().filter(|z| z.re < 0.0 || z.norm() < 1e-2 * re[n / 2].abs()).collect();
            lst.sort_by(|a, b| a.re.total_cmp(&b.re));
            for z in lst.iter().take(40) {
                println!("  outlier {:.4e} {:+.4e}i", z.re, z.im);
            }
        }
        "gmres" => {
            let o = hopts(&a);
            let mut solver =
                hierarchical_saddle(pool.clone(), &c.surface, &c.positions, &c.exponents, &c.weights, &ex, &cq, &o, None)
                    .unwrap();
            println!(
                "eps {:?} coarse {:?}",
                solver.preconditioner().epsilon(),
                solver.preconditioner().coarse().map(|c| (c.dim(), c.conditioning()))
            );
            let mut rhs = g.clone();
            rhs.extend(std::iter::repeat_n(0.0, cq.nrows()));
            let mut go = o.gmres;
            go.max_iterations = a.max;
            let r = gmres(solver.operator(), solver.preconditioner(), &rhs, &go);
            match r {
                Ok((_, rep)) => {
                    let h: Vec<String> = rep.history.iter().step_by(10).map(|x| format!("{x:.2e}")).collect();
                    println!(
                        "gmres iterations {} true {:.3e} converged {} cycles {:?}\n history/10 {}",
                        rep.iterations,
                        rep.relative_residual,
                        rep.converged,
                        rep.cycle_residuals.iter().map(|x| format!("{x:.2e}")).collect::<Vec<_>>(),
                        h.join(" ")
                    );
                }
                Err(e) => println!("gmres error {e}"),
            }
            if a.schur && !nl.is_empty() {
                let schur = SchurComplement::new(&mut solver, nl.clone(), cn.clone());
                match schur {
                    Ok(s) => match s.direction(&mut solver, &g) {
                        Ok(d) => println!(
                            "schur ok: records {:?} direction its {} res {:.2e}",
                            s.records().iter().map(|r| (r.iterations, r.relative_residual)).collect::<Vec<_>>(),
                            d.record.iterations,
                            d.record.relative_residual
                        ),
                        Err(e) => println!("schur direction error {e}"),
                    },
                    Err(e) => println!("schur error {e}"),
                }
            }
        }
        "product" => {
            let am = dense_metric_with_extras(&pool, &c.surface, &c.positions, &c.exponents, &c.weights, &ex).unwrap();
            let dop = DenseSaddleOperator::new(&am, &cq).unwrap();
            let mut o = hopts(&a);
            o.cluster.chi = 0.0;
            let solver =
                hierarchical_saddle(pool.clone(), &c.surface, &c.positions, &c.exponents, &c.weights, &ex, &cq, &o, None)
                    .unwrap();
            let n = dop.dim();
            let mut seed = 12345u64;
            let x: Vec<f64> = (0..n)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
                })
                .collect();
            let (mut y1, mut y2) = (vec![0.0; n], vec![0.0; n]);
            dop.apply(&x, &mut y1).unwrap();
            solver.operator().apply(&x, &mut y2).unwrap();
            let d: f64 = y1.iter().zip(&y2).map(|(p, q)| (p - q) * (p - q)).sum::<f64>().sqrt();
            let nn: f64 = y1.iter().map(|p| p * p).sum::<f64>().sqrt();
            println!("product chi=0 hierarchical vs dense: relative {:.3e}", d / nn);
        }
        x => panic!("unknown command {x}"),
    }
}
