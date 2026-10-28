//! WP16b probes (`SPEC.md` §17p): the geometry of the failed runs of the Sec. 7 evaluation.
//! Env-driven and `#[ignore]`d, e.g.
//!
//! ```text
//! RC_CURVE=final.obj cargo test --release -p ddg-repulsive-curves --test wp16b_probe \
//!   fold_probe -- --ignored --nocapture
//! ```
//!
//! * `fold_probe`: at `RC_CURVE`, the closest disjoint edge pair (or `RC_I`, `RC_J`): lengths of
//!   `I`, the joining edge `K` (second neighbours) and `J` relative to the mean edge length, the
//!   closest-point parameters, the turning angles at the ends of `K`, the paper's ordered terms
//!   `e^p` (Eq. 18) and the accurate pair integrals `e^n` (§3.4) of both orders, and the same
//!   along a pure separation of the pair (`J` moved along the normal of closest approach).

use ddg_repulsive_curves::energy::near_pair::{
    InnerTables, closest_points, joined_by_edge, ordered_near_energy, paper_pair,
};
use ddg_repulsive_curves::step::line_search::closest_pair_within;
use ddg_repulsive_curves::{Exponents, read_curve_obj};

type P = [f64; 3];

fn env(name: &str) -> Option<String> {
    std::env::var(name).ok().filter(|s| !s.is_empty())
}

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn norm(a: P) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn dot(a: P, b: P) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn terms(t: &InnerTables, e: Exponents, p: [P; 2], q: [P; 2]) -> (f64, f64, f64, f64) {
    let cl = closest_points(p[0], p[1], q[0], q[1]);
    let sw = |g: &[[f64; 3]; 4]| [g[2], g[3], g[0], g[1]];
    let (gs, gd) = (sw(&cl.grad_s), sw(&cl.grad_d));
    let (n_ij, _) = ordered_near_energy(t, p, q, (cl.t, &cl.grad_t), (cl.d, &cl.grad_d), 4, false);
    let (n_ji, _) = ordered_near_energy(t, q, p, (cl.s, &gs), (cl.d, &gd), 4, false);
    (paper_pair(p, q, e).0, paper_pair(q, p, e).0, n_ij, n_ji)
}

#[test]
#[ignore]
fn fold_probe() {
    let Some(path) = env("RC_CURVE") else {
        eprintln!("RC_CURVE not set; skipped");
        return;
    };
    let (net, pos) = read_curve_obj(&path).unwrap();
    let e = Exponents::new(2.0, 4.5).unwrap();
    let tables = InnerTables::new(e);
    let edges = net.edges();
    let len = |k: usize| norm(sub(pos[edges[k][1] as usize], pos[edges[k][0] as usize]));
    let lbar = (0..edges.len()).map(len).sum::<f64>() / edges.len() as f64;
    let (i, j) = match (env("RC_I"), env("RC_J")) {
        (Some(a), Some(b)) => (a.parse::<usize>().unwrap(), b.parse::<usize>().unwrap()),
        _ => {
            let (a, b, _) = closest_pair_within(&net, &pos, lbar).unwrap();
            (a as usize, b as usize)
        }
    };
    let ends = |k: usize| [pos[edges[k][0] as usize], pos[edges[k][1] as usize]];
    let (p, q) = (ends(i), ends(j));
    let cl = closest_points(p[0], p[1], q[0], q[1]);
    println!(
        "PAIR {i} {j} second_neighbours {} d/lbar {:.3e} s {:.4} t {:.4} lI {:.4} lJ {:.4}",
        joined_by_edge(&net, i, j),
        cl.d / lbar,
        cl.s,
        cl.t,
        len(i) / lbar,
        len(j) / lbar
    );
    println!("EDGES I {:?} J {:?}", edges[i], edges[j]);
    // The joining edge.
    for (k, ek) in edges.iter().enumerate() {
        let touches = |x: usize| edges[x].iter().any(|v| ek.contains(v));
        if k != i && k != j && touches(i) && touches(j) {
            let tk = sub(pos[ek[1] as usize], pos[ek[0] as usize]);
            let ti = sub(p[1], p[0]);
            let tj = sub(q[1], q[0]);
            let ang = |a: P, b: P| (dot(a, b) / (norm(a) * norm(b))).clamp(-1.0, 1.0).acos();
            println!(
                "K {k} {:?} lK {:.4} turn(I,K) {:.1} deg turn(K,J) {:.1} deg angle(I,J) {:.1} deg",
                ek,
                len(k) / lbar,
                ang(ti, tk).to_degrees(),
                ang(tk, tj).to_degrees(),
                ang(ti, tj).to_degrees()
            );
        }
    }
    let (pij, pji, nij, nji) = terms(&tables, e, p, q);
    println!("TERMS e^p_IJ {pij:.6e} e^p_JI {pji:.6e} e^n_IJ {nij:.6e} e^n_JI {nji:.6e}");
    // Separation along the normal of closest approach.
    let pi = [
        p[0][0] + cl.s * (p[1][0] - p[0][0]),
        p[0][1] + cl.s * (p[1][1] - p[0][1]),
        p[0][2] + cl.s * (p[1][2] - p[0][2]),
    ];
    let qj = [
        q[0][0] + cl.t * (q[1][0] - q[0][0]),
        q[0][1] + cl.t * (q[1][1] - q[0][1]),
        q[0][2] + cl.t * (q[1][2] - q[0][2]),
    ];
    let n = sub(qj, pi);
    let nn = norm(n);
    let n = [n[0] / nn, n[1] / nn, n[2] / nn];
    for k in [0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3] {
        let sh = k * lbar;
        let q2 = q.map(|x| [x[0] + sh * n[0], x[1] + sh * n[1], x[2] + sh * n[2]]);
        let c2 = closest_points(p[0], p[1], q2[0], q2[1]);
        let (a, b, c, d) = terms(&tables, e, p, q2);
        println!(
            "SEP shift/lbar {k:.0e} d/lbar {:.3e} s {:.3} t {:.3} e^p {:.5e} e^n {:.5e} ratio {:.4}",
            c2.d / lbar,
            c2.s,
            c2.t,
            a + b,
            c + d,
            (c + d) / (a + b)
        );
    }
}

/// At `RC_CURVE` with the rest lengths and barycenter of `RC_REST` (the input knot): the edge
/// multipliers `λ_I` of Eq. 22 (exact energy, dense solve), the residual violations `Φ_I`, the
/// predicted change `−λ·Φ` of the objective under the projection of the curve onto its
/// constraints, the largest `|λ_I|` with the edge length, and the vertices with the largest
/// `|g̃_v|`.
#[test]
#[ignore]
fn multiplier_probe() {
    use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
    use ddg_repulsive_curves::metric::{GramQuadrature, sobolev_gram};
    use ddg_repulsive_curves::{Constraint, ConstraintSet, DenseSaddle, projected_gradient};
    use std::sync::Arc;
    let (Some(path), Some(rest)) = (env("RC_CURVE"), env("RC_REST")) else {
        eprintln!("RC_CURVE, RC_REST not set; skipped");
        return;
    };
    let (net, pos) = read_curve_obj(&path).unwrap();
    let (_, pos0) = read_curve_obj(&rest).unwrap();
    let pool = Arc::new(make_pool(Some(4)).unwrap());
    let e = Exponents::new(2.0, 4.5).unwrap();
    let (en, diff) = energy_and_differential(&pool, &net, &pos, e);
    let geom0 = net.edge_geometry(&pos0);
    let (mut s, mut l) = ([0.0f64; 3], 0.0);
    for (m, &len) in geom0.midpoint.iter().zip(&geom0.length) {
        for c in 0..3 {
            s[c] += len * m[c];
        }
        l += len;
    }
    let cons = vec![
        Constraint::Barycenter {
            target: Some(s.map(|x| x / l)),
        },
        Constraint::EdgeLengths {
            edges: None,
            targets: Some(geom0.length.clone()),
        },
    ];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    set.analyze(&net, &pos).unwrap();
    let c = set.jacobian(&net, &pos).unwrap();
    let phi = set.values(&net, &pos).unwrap();
    let g = sobolev_gram(&pool, &net, &pos, e, GramQuadrature::Midpoint);
    let mut solver = DenseSaddle::new(pool.clone(), &g.a, &c, set.layout()).unwrap();
    let pg = projected_gradient(&mut solver, set.layout(), &diff).unwrap();
    let lam = &pg.lambda;
    let pred: f64 = lam.iter().zip(&phi).map(|(a, b)| a * b).sum();
    let geom = net.edge_geometry(&pos);
    let lbar = geom.length.iter().sum::<f64>() / geom.length.len() as f64;
    let maxphi = phi.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    println!(
        "MULT E {en:.8} |g| {:.4e} dE.g {:.4e} lbar {lbar:.4e} max|phi| {maxphi:.3e} \
         lambda.phi {pred:.4e} (F(P(curve)) - F(curve) ~ -lambda.phi or +, by the sign convention)",
        pg.gradient_norm, pg.sobolev_norm_sq
    );
    let mut idx: Vec<usize> = (3..lam.len()).collect();
    idx.sort_by(|&a, &b| lam[b].abs().total_cmp(&lam[a].abs()));
    for &r in idx.iter().take(12) {
        let edge = r - 3;
        println!(
            "LAMBDA edge {edge} {:?} lambda {:.4e} l/lbar {:.4} phi {:.3e}",
            net.edges()[edge],
            lam[r],
            geom.length[edge] / lbar,
            phi[r]
        );
    }
    let mut lams: Vec<f64> = lam[3..].iter().map(|x| x.abs()).collect();
    lams.sort_by(f64::total_cmp);
    println!(
        "LAMBDA |.| median {:.3e} p90 {:.3e} p99 {:.3e} max {:.3e}",
        lams[lams.len() / 2],
        lams[lams.len() * 9 / 10],
        lams[lams.len() * 99 / 100],
        lams[lams.len() - 1]
    );
    // Correlation with edge length: mean |λ| of the 20 shortest edges vs all.
    let mut by_len: Vec<usize> = (0..geom.length.len()).collect();
    by_len.sort_by(|&a, &b| geom.length[a].total_cmp(&geom.length[b]));
    let short: f64 = by_len[..20].iter().map(|&k| lam[k + 3].abs()).sum::<f64>() / 20.0;
    println!(
        "LAMBDA mean |.| of the 20 shortest edges {short:.3e} (shortest {:.4} lbar), of all {:.3e}",
        geom.length[by_len[0]] / lbar,
        lams.iter().sum::<f64>() / lams.len() as f64
    );
    let mut vs: Vec<(usize, f64)> = pg
        .direction
        .iter()
        .enumerate()
        .map(|(v, x)| (v, norm(*x)))
        .collect();
    vs.sort_by(|a, b| b.1.total_cmp(&a.1));
    for &(v, x) in vs.iter().take(6) {
        let inc: Vec<String> = net
            .incident_edges(v as u32)
            .iter()
            .map(|&k| format!("{k}:{:.3}", geom.length[k as usize] / lbar))
            .collect();
        println!(
            "GTILDE vertex {v} |g| {x:.4e} incident edges (l/lbar) {}",
            inc.join(" ")
        );
    }
}

/// At `RC_CURVE` (targets of `RC_REST`), the trials `γ − τ g̃` (exact energy, dense solve, the
/// direction of Eq. 22) projected with and without the Newton polish (`Absolute(1e-4)`, 3 chord
/// iterations, early abort on): status, iterations, the relative-violation history, `F − F(γ)`,
/// the largest displacement of the projection and where it is, for `τ` from `RC_TAU0` (default
/// 0.4) halved `RC_N` (default 16) times, then `τ = 0`.
#[test]
#[ignore]
fn projection_probe() {
    use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
    use ddg_repulsive_curves::metric::{GramQuadrature, sobolev_gram};
    use ddg_repulsive_curves::solve::{ProjectionOptions, Tolerance};
    use ddg_repulsive_curves::step::line_search::{ProjectionRules, project_with_path_polished};
    use ddg_repulsive_curves::{
        Constraint, ConstraintSet, DenseSaddle, EnergyEval, energy, projected_gradient,
    };
    use std::sync::Arc;
    let (Some(path), Some(rest)) = (env("RC_CURVE"), env("RC_REST")) else {
        eprintln!("RC_CURVE, RC_REST not set; skipped");
        return;
    };
    let (net, pos) = read_curve_obj(&path).unwrap();
    let (_, pos0) = read_curve_obj(&rest).unwrap();
    let pool = Arc::new(make_pool(Some(4)).unwrap());
    let e = Exponents::new(2.0, 4.5).unwrap();
    let (f0, diff) = energy_and_differential(&pool, &net, &pos, e);
    let geom0 = net.edge_geometry(&pos0);
    let lbar0 = geom0.length.iter().sum::<f64>() / geom0.length.len() as f64;
    let (mut s, mut l) = ([0.0f64; 3], 0.0);
    for (m, &len) in geom0.midpoint.iter().zip(&geom0.length) {
        for c in 0..3 {
            s[c] += len * m[c];
        }
        l += len;
    }
    let cons = vec![
        Constraint::Barycenter {
            target: Some(s.map(|x| x / l)),
        },
        Constraint::EdgeLengths {
            edges: None,
            targets: Some(geom0.length.clone()),
        },
    ];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    set.analyze(&net, &pos).unwrap();
    let c = set.jacobian(&net, &pos).unwrap();
    let g = sobolev_gram(&pool, &net, &pos, e, GramQuadrature::Midpoint);
    let mut solver = DenseSaddle::new(pool.clone(), &g.a, &c, set.layout()).unwrap();
    let pg = projected_gradient(&mut solver, set.layout(), &diff).unwrap();
    let opts = ProjectionOptions {
        tolerance: Tolerance::Absolute(1e-4),
        max_iterations: 3,
    };
    let tau0: f64 = env("RC_TAU0").and_then(|s| s.parse().ok()).unwrap_or(0.4);
    let n: usize = env("RC_N").and_then(|s| s.parse().ok()).unwrap_or(16);
    let lens = net.edge_geometry(&pos).length;
    println!("BASE F {f0:.9} slope {:.4e}", pg.sobolev_norm_sq);
    let mut taus: Vec<f64> = (0..=n).map(|k| tau0 * 0.5f64.powi(k as i32)).collect();
    taus.push(0.0);
    for tau in taus {
        let trial: Vec<P> = pos
            .iter()
            .zip(&pg.direction)
            .map(|(x, d)| [x[0] - tau * d[0], x[1] - tau * d[1], x[2] - tau * d[2]])
            .collect();
        for (name, rules) in [
            (
                "polish",
                ProjectionRules {
                    newton_polish: Some(100.0),
                    early_abort: true,
                },
            ),
            (
                "chord ",
                ProjectionRules {
                    newton_polish: None,
                    early_abort: true,
                },
            ),
        ] {
            let (proj, _) = project_with_path_polished(
                &mut solver,
                &set,
                &net,
                trial.clone(),
                &opts,
                lbar0,
                rules,
            )
            .unwrap();
            let f = energy(&net, &proj.positions, e, &EnergyEval::Exact, Some(4))
                .map(|r| r.energy)
                .unwrap_or(f64::NAN);
            let (mut vmax, mut dmax) = (0, 0.0f64);
            for (v, (a, b)) in proj.positions.iter().zip(&trial).enumerate() {
                let d = norm(sub(*a, *b));
                if d > dmax {
                    (vmax, dmax) = (v, d);
                }
            }
            let inc: Vec<String> = net
                .incident_edges(vmax as u32)
                .iter()
                .map(|&k| format!("{k}:{:.3}", lens[k as usize] / lbar0))
                .collect();
            let hist: Vec<String> = proj.history.iter().map(|h| format!("{h:.2e}")).collect();
            println!(
                "TRIAL tau {tau:.4e} {name} {:?} it {} hist [{}] F-F0 {:+.4e} armijo_rhs {:.3e} \
                 max_disp/lbar {:.3e} at v{vmax} ({})",
                proj.status,
                proj.iterations,
                hist.join(" "),
                f - f0,
                -1e-4 * tau * pg.sobolev_norm_sq,
                dmax / lbar0,
                inc.join(" ")
            );
        }
    }
}
