//! Sanity comparison of the optimizers (`SPEC.md` §14.2, WP9): one method on one OBJ, energy
//! against wall time.
//!
//! ```text
//! cargo run --release -p ddg-repulsive-surfaces --example compare_methods -- OBJ METHOD [STEPS] [THREADS]
//! METHOD = hs | l2 | h1 | h2 | aqp | h1-lbfgs | h1-lbfgs-ref | bcqn | bcqn-ref
//! ```
//!
//! `hs`, `l2`, `h1`, `h2` use hard barycenter + area constraints; the others the soft problem
//! (`soft_problem`: area penalty, weight 100). Prints `step seconds E_tp F area` per step and a
//! `summary` line.

use std::time::Instant;

use ddg_repulsive_surfaces::constraints::Constraint;
use ddg_repulsive_surfaces::optimize::soft_problem;
use ddg_repulsive_surfaces::{
    BcqnVariant, FaceGeometry, FlowOptions, L2Kind, LbfgsVariant, Method, Metric, Session, StopRule, read_obj,
};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (obj, name) = (&a[1], a[2].as_str());
    let steps: usize = a.get(3).map_or(50, |s| s.parse().expect("STEPS"));
    let threads: usize = a.get(4).map_or(4, |s| s.parse().expect("THREADS"));
    let (surface, positions) = read_obj(obj).expect("read OBJ");
    let hard = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        stop: StopRule::Off,
        threads: Some(threads),
        ..FlowOptions::default()
    };
    let mut o = hard.clone();
    match name {
        "hs" => {}
        "l2" => o.metric = Metric::L2(L2Kind::Mass),
        "h1" => o.metric = Metric::H1,
        "h2" => o.metric = Metric::H2,
        _ => {
            o = soft_problem(&surface, &positions, &hard).expect("soft problem");
            o.method = match name {
                "aqp" => Method::Aqp { kappa: 100.0 },
                "h1-lbfgs" => Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Standard },
                "h1-lbfgs-ref" => Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Reference },
                "bcqn" => Method::Bcqn { memory: 20, variant: BcqnVariant::Standard },
                "bcqn-ref" => Method::Bcqn { memory: 20, variant: BcqnVariant::Reference },
                _ => panic!("unknown method {name}"),
            };
        }
    }
    let faces0 = surface.n_faces();
    let mut ses = Session::new(surface, positions, o).expect("session");
    let t0 = Instant::now();
    let (mut e0, mut e_last, mut a0) = (f64::NAN, f64::NAN, f64::NAN);
    for _ in 0..steps {
        let r = ses.step().expect("step");
        let area: f64 = FaceGeometry::new(ses.surface(), ses.positions()).unwrap().areas.iter().sum();
        if e0.is_nan() {
            e0 = r.energy_before;
            a0 = area;
        }
        e_last = r.energy_term_after;
        println!("{} {:.3} {:.6e} {:.6e} {:.6e}", r.iteration, t0.elapsed().as_secs_f64(), e_last, r.energy_after, area);
    }
    let area: f64 = FaceGeometry::new(ses.surface(), ses.positions()).unwrap().areas.iter().sum();
    println!(
        "summary method={name} faces0={faces0} faces={} steps={steps} seconds={:.2} F0={:.5e} E_tp={:.5e} area_ratio={:.4} diag={:?}",
        ses.surface().n_faces(),
        t0.elapsed().as_secs_f64(),
        e0,
        e_last,
        area / a0,
        ses.method_diagnostics()
    );
}
