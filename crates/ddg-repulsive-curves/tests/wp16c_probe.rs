//! WP16c probes (`SPEC.md` §17r): why the one-sided second-neighbour rule (§18 Q10) slows the
//! hard-constraint flows on some knots. Env-driven and `#[ignore]`d, e.g.
//!
//! ```text
//! RC_CURVE=step5.obj RC_REST=knot.obj cargo test --release -p ddg-repulsive-curves \
//!   --test wp16c_probe sn_step_probe -- --ignored --nocapture
//! ```
//!
//! * `sn_step_probe`: at `RC_CURVE` (constraint targets of `RC_REST`: barycenter and edge
//!   lengths, as `eval_knots`), every second-neighbour pair `I, K, J` of the near-pair band with
//!   its geometry (turning angles at the ends of `K`, `r = d/ℓ_IJ`, `s*`, `t*`, lengths), the
//!   summed terms `E^n`, `E^p` and its one-sided correction; the projected gradient `g̃` of the
//!   corrected energy (Q10 default) next to the one without second neighbours (WP16), the share of
//!   `g̃` and of the rate `dF·g̃` that the second-neighbour term contributes; and the objective
//!   along `γ − τ g̃` (unprojected) for both energies, as the ratio
//!   `(F(γ − τ g̃) − F(γ)) / (−τ dF·g̃)` (1 for a linear model, Armijo needs ≥ 1e-4 … here
//!   the guarded search's `c`), split into the paper's part and the corrections.
//!   `RC_METRIC=hs|h1` (default `hs`), `RC_DIR=new|old` (the direction of the line scan:
//!   default `new`), `RC_TAUS` (comma separated; default `1e-4 · 2^k`, `k = 0..12`).

use std::sync::Arc;

use ddg_repulsive_curves::ccd::{CcdOptions, Contact, ContactPrimitive, first_contact};
use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::energy::near_pair::{
    InnerTables, NearPairOptions, closest_points, correction, joining_ends, near_pairs,
    ordered_near_energy, pair_correction_of_kind, paper_pair,
};
use ddg_repulsive_curves::metric::sparse_classic::classic_metric;
use ddg_repulsive_curves::metric::{GramQuadrature, Metric, sobolev_dense::sobolev_gram};
use ddg_repulsive_curves::{
    Constraint, ConstraintSet, DenseSaddle, EnergyEval, Exponents, SaddleSolver, SparseSaddle,
    energy, projected_gradient, read_curve_obj,
};

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
fn angle(a: P, b: P) -> f64 {
    (dot(a, b) / (norm(a) * norm(b)))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}
fn vnorm(v: &[[f64; 3]]) -> f64 {
    v.iter().map(|x| dot(*x, *x)).sum::<f64>().sqrt()
}
fn vdot(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter().zip(b).map(|(x, y)| dot(*x, *y)).sum()
}
fn vadd(a: &[[f64; 3]], b: Option<&Vec<[f64; 3]>>) -> Vec<[f64; 3]> {
    match b {
        None => a.to_vec(),
        Some(b) => a
            .iter()
            .zip(b)
            .map(|(x, y)| [x[0] + y[0], x[1] + y[1], x[2] + y[2]])
            .collect(),
    }
}

#[test]
#[ignore]
fn sn_step_probe() {
    let Some(path) = env("RC_CURVE") else {
        eprintln!("RC_CURVE not set; skipped");
        return;
    };
    let (net, pos) = read_curve_obj(&path).unwrap();
    let threads: usize = env("RC_THREADS").and_then(|s| s.parse().ok()).unwrap_or(4);
    let pool = Arc::new(make_pool(Some(threads)).unwrap());
    let e = Exponents::new(2.0, 4.5).unwrap();
    let new = NearPairOptions::default();
    let old = NearPairOptions {
        second_neighbours: false,
        ..new
    };
    let edges = net.edges();
    let geom = net.edge_geometry(&pos);
    let lbar = geom.length.iter().sum::<f64>() / geom.length.len() as f64;
    let all = |_: usize, _: usize| true;

    // ---- The second-neighbour pairs of the band.
    let tables = InnerTables::new(e);
    let ends = |k: usize| [pos[edges[k][0] as usize], pos[edges[k][1] as usize]];
    let (mut in_band, mut active) = (0usize, 0usize);
    let mut active_vertices: Vec<u32> = Vec::new();
    for (i, j) in near_pairs(&net, &pos, &new) {
        let (i, j) = (i as usize, j as usize);
        let Some((a, b)) = joining_ends(&net, i, j) else {
            continue;
        };
        in_band += 1;
        let (p, q) = (ends(i), ends(j));
        let pc = pair_correction_of_kind(&tables, e, &new, p, q, Some((a, b)), true, true, false);
        if pc.value == 0.0 {
            continue;
        }
        active += 1;
        let cl = closest_points(p[0], p[1], q[0], q[1]);
        let sw = |g: &[[f64; 3]; 4]| [g[2], g[3], g[0], g[1]];
        let (gs, gd) = (sw(&cl.grad_s), sw(&cl.grad_d));
        let (n_ij, _) = ordered_near_energy(
            &tables,
            p,
            q,
            (cl.t, &cl.grad_t),
            (cl.d, &cl.grad_d),
            4,
            false,
        );
        let (n_ji, _) = ordered_near_energy(&tables, q, p, (cl.s, &gs), (cl.d, &gd), 4, false);
        let ep = paper_pair(p, q, e).0 + paper_pair(q, p, e).0;
        let en = n_ij + n_ji;
        // K runs from p_a (end of I) to q_b (end of J).
        let (pa, pfar, qb, qfar) = (p[a], p[1 - a], q[b], q[1 - b]);
        let tk = sub(qb, pa);
        let turn_ik = angle(sub(pa, pfar), tk);
        let turn_kj = angle(tk, sub(qfar, qb));
        let ang_ij = angle(sub(pa, pfar), sub(qfar, qb));
        // Parameters measured from the ends at K (0 = at K).
        let s_k = if a == 0 { cl.s } else { 1.0 - cl.s };
        let t_k = if b == 0 { cl.t } else { 1.0 - cl.t };
        println!(
            "SN I {i} J {j} verts I {:?} J {:?} r {:.3} w {:.3} d/lbar {:.3e} s_fromK {:.3} t_fromK {:.3} \
             lI {:.3} lK {:.3} lJ {:.3} turnIK {:.1} turnKJ {:.1} angIJ {:.1} En {:.4e} Ep {:.4e} \
             En/Ep {:.3e} corr {:.4e}",
            edges[i],
            edges[j],
            pc.ratio,
            pc.weight,
            cl.d / lbar,
            s_k,
            t_k,
            geom.length[i] / lbar,
            norm(tk) / lbar,
            geom.length[j] / lbar,
            turn_ik,
            turn_kj,
            ang_ij,
            en,
            ep,
            en / ep,
            pc.value
        );
        active_vertices.extend(edges[i]);
        active_vertices.extend(edges[j]);
    }

    // ---- Energies and differentials.
    let (e_paper, d_paper) = energy_and_differential(&pool, &net, &pos, e);
    let c_new = correction(&pool, &net, &pos, e, &new, &all, true);
    let c_old = correction(&pool, &net, &pos, e, &old, &all, true);
    let (en_new, dn_new) = c_new
        .as_ref()
        .map_or((0.0, None), |c| (c.energy, c.differential.as_ref()));
    let (en_old, dn_old) = c_old
        .as_ref()
        .map_or((0.0, None), |c| (c.energy, c.differential.as_ref()));
    let d_new = vadd(&d_paper, dn_new);
    let d_old = vadd(&d_paper, dn_old);
    let d_sn: Vec<[f64; 3]> = d_new.iter().zip(&d_old).map(|(x, y)| sub(*x, *y)).collect();
    println!(
        "ENERGY paper {e_paper:.6} corr_first {en_old:.6} corr_second {:.6} F_new {:.6} F_old {:.6} \
         sn_in_band {in_band} sn_active {active} |dE_paper| {:.4e} |dCorr_first| {:.4e} |dCorr_second| {:.4e}",
        en_new - en_old,
        e_paper + en_new,
        e_paper + en_old,
        vnorm(&d_paper),
        vnorm(
            &vadd(&d_old, None)
                .iter()
                .zip(&d_paper)
                .map(|(x, y)| sub(*x, *y))
                .collect::<Vec<_>>()
        ),
        vnorm(&d_sn)
    );

    // ---- Finite differences of the corrected energy (Q10) and of its second-neighbour part
    // along the second-neighbour differential and a pseudo-random field (RC_FD=1).
    if env("RC_FD").as_deref() == Some("1") {
        let f_new = |p: &[[f64; 3]]| -> (f64, f64) {
            let ep = energy(&net, p, e, &EnergyEval::Exact, Some(threads))
                .unwrap()
                .energy;
            let cn = correction(&pool, &net, p, e, &new, &all, false).map_or(0.0, |c| c.energy);
            let co = correction(&pool, &net, p, e, &old, &all, false).map_or(0.0, |c| c.energy);
            (ep + cn, cn - co)
        };
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let random: Vec<[f64; 3]> = (0..pos.len()).map(|_| [rnd(), rnd(), rnd()]).collect();
        let sn_dir: Vec<[f64; 3]> = if vnorm(&d_sn) > 0.0 {
            d_sn.clone()
        } else {
            random.clone()
        };
        for (name, v) in [("sn-differential", &sn_dir), ("random", &random)] {
            let nv = vnorm(v);
            let v: Vec<[f64; 3]> = v
                .iter()
                .map(|x| [x[0] / nv, x[1] / nv, x[2] / nv])
                .collect();
            let pred = vdot(&d_new, &v);
            let pred_sn = vdot(&d_sn, &v);
            for h in [1e-4, 1e-5, 1e-6] {
                let hh = h * lbar;
                let shift = |s: f64| -> Vec<[f64; 3]> {
                    pos.iter()
                        .zip(&v)
                        .map(|(x, d)| [x[0] + s * d[0], x[1] + s * d[1], x[2] + s * d[2]])
                        .collect()
                };
                let (fp, sp) = f_new(&shift(hh));
                let (fm, sm) = f_new(&shift(-hh));
                let fd = (fp - fm) / (2.0 * hh);
                let fd_sn = (sp - sm) / (2.0 * hh);
                println!(
                    "FD {name} h/lbar {h:.0e} dF.v {pred:.8e} fd {fd:.8e} rel {:.2e} | dSN.v {pred_sn:.8e} fd {fd_sn:.8e} rel {:.2e}",
                    (fd - pred).abs() / pred.abs().max(1e-300),
                    (fd_sn - pred_sn).abs() / pred_sn.abs().max(1e-300)
                );
            }
        }
    }

    let Some(rest) = env("RC_REST") else {
        return;
    };
    // ---- The projected gradients (constraints of eval_knots).
    let (_, pos0) = read_curve_obj(&rest).unwrap();
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
    let metric = env("RC_METRIC").unwrap_or_else(|| "hs".into());
    let mut solver: Box<dyn SaddleSolver> = if metric == "h1" {
        let m = classic_metric(&net, &pos, &Metric::H1).unwrap();
        Box::new(SparseSaddle::new(&m, &c, set.layout()).unwrap())
    } else {
        let g = sobolev_gram(&pool, &net, &pos, e, GramQuadrature::Midpoint);
        Box::new(DenseSaddle::new(pool.clone(), &g.a, &c, set.layout()).unwrap())
    };
    let g_new = projected_gradient(solver.as_mut(), set.layout(), &d_new).unwrap();
    let g_old = projected_gradient(solver.as_mut(), set.layout(), &d_old).unwrap();
    let g_sn = projected_gradient(solver.as_mut(), set.layout(), &d_sn).unwrap();
    let (gn, go) = (&g_new.direction, &g_old.direction);
    let diff: Vec<[f64; 3]> = gn.iter().zip(go).map(|(x, y)| sub(*x, *y)).collect();
    println!(
        "GRADIENT metric {metric} |g_new| {:.4e} |g_old| {:.4e} |g_new - g_old|/|g_new| {:.3} rate_new {:.4e} \
         rate_old {:.4e} dSN.g_new {:.4e} (share of rate {:.3}) dSN.g_old {:.4e} |g_sn| {:.4e} \
         g_sn.A.g_sn {:.4e} cos(g_new,g_old) {:.4}",
        vnorm(gn),
        vnorm(go),
        vnorm(&diff) / vnorm(gn),
        g_new.sobolev_norm_sq,
        g_old.sobolev_norm_sq,
        vdot(&d_sn, gn),
        vdot(&d_sn, gn) / g_new.sobolev_norm_sq,
        vdot(&d_sn, go),
        vnorm(&g_sn.direction),
        g_sn.sobolev_norm_sq,
        vdot(gn, go) / (vnorm(gn) * vnorm(go))
    );
    // Where the direction is large.
    let mut idx: Vec<usize> = (0..pos.len()).collect();
    idx.sort_by(|&a, &b| norm(gn[b]).total_cmp(&norm(gn[a])));
    let near_sn = |v: usize| {
        active_vertices
            .iter()
            .map(|&w| (w as i64 - v as i64).unsigned_abs())
            .min()
            .unwrap_or(u64::MAX)
    };
    for &v in idx.iter().take(8) {
        println!(
            "TOPV v {v} |g_new|/lbar {:.4e} |g_old|/lbar {:.4e} |g_sn| {:.4e} index distance to an active SN pair {}",
            norm(gn[v]) / lbar,
            norm(go[v]) / lbar,
            norm(g_sn.direction[v]),
            near_sn(v)
        );
    }

    // ---- Line scan along −g̃ (no projection).
    let dir = if env("RC_DIR").as_deref() == Some("old") {
        go
    } else {
        gn
    };
    // The CCD of the guarded search along −g̃ (gap 0.1: a pair's distance down to a tenth).
    let vel: Vec<[f64; 3]> = dir.iter().map(|g| [-g[0], -g[1], -g[2]]).collect();
    let rep = first_contact(&pool, &net, &pos, &vel, 10.0, &[], &CcdOptions::default()).unwrap();
    match rep.contact {
        Some(Contact {
            tau,
            primitive: ContactPrimitive::EdgeEdge { a, b },
        }) => {
            let (a, b) = (a as usize, b as usize);
            let (p, q) = (ends(a), ends(b));
            let cl = closest_points(p[0], p[1], q[0], q[1]);
            let joined = joining_ends(&net, a, b);
            let mut extra = String::new();
            if let Some((ja, jb)) = joined {
                let pc = pair_correction_of_kind(
                    &tables,
                    e,
                    &new,
                    p,
                    q,
                    Some((ja, jb)),
                    true,
                    true,
                    false,
                );
                let (pa, pfar, qb, qfar) = (p[ja], p[1 - ja], q[jb], q[1 - jb]);
                let tk = sub(qb, pa);
                extra = format!(
                    " turnIK {:.1} turnKJ {:.1} corr {:.4e} r {:.3}",
                    angle(sub(pa, pfar), tk),
                    angle(tk, sub(qfar, qb)),
                    pc.value,
                    pc.ratio
                );
            }
            println!(
                "CCD tau_max {tau:.4e} edges {a} {b} verts {:?} {:?} second_neighbours {} d/lbar {:.3e} s {:.3} t {:.3} |g| at the ends/lbar {:.2} {:.2} {:.2} {:.2}{extra}",
                edges[a],
                edges[b],
                joined.is_some(),
                cl.d / lbar,
                cl.s,
                cl.t,
                norm(dir[edges[a][0] as usize]) / lbar,
                norm(dir[edges[a][1] as usize]) / lbar,
                norm(dir[edges[b][0] as usize]) / lbar,
                norm(dir[edges[b][1] as usize]) / lbar,
            );
        }
        Some(c) => println!("CCD tau_max {:.4e} {:?}", c.tau, c.primitive),
        None => println!("CCD none up to 10"),
    }
    let rate_new = vdot(&d_new, dir);
    let rate_old = vdot(&d_old, dir);
    let taus: Vec<f64> = match env("RC_TAUS") {
        Some(t) => t.split(',').map(|x| x.parse().unwrap()).collect(),
        None => (0..13).map(|k| 1e-4 * 2f64.powi(k)).collect(),
    };
    let f0_new = e_paper + en_new;
    let f0_old = e_paper + en_old;
    for tau in taus {
        let p: Vec<[f64; 3]> = pos
            .iter()
            .zip(dir)
            .map(|(x, g)| [x[0] - tau * g[0], x[1] - tau * g[1], x[2] - tau * g[2]])
            .collect();
        if net.validate_positions(&p).is_err() {
            println!("LINE tau {tau:.3e} invalid");
            continue;
        }
        let ep = energy(&net, &p, e, &EnergyEval::Exact, Some(threads))
            .unwrap()
            .energy;
        let cn = correction(&pool, &net, &p, e, &new, &all, false).map_or(0.0, |c| c.energy);
        let co = correction(&pool, &net, &p, e, &old, &all, false).map_or(0.0, |c| c.energy);
        let n_act = near_pairs(&net, &p, &new)
            .into_iter()
            .filter(|&(i, j)| {
                let (i, j) = (i as usize, j as usize);
                joining_ends(&net, i, j).is_some_and(|ab| {
                    pair_correction_of_kind(
                        &tables,
                        e,
                        &new,
                        ends_at(&p, edges, i),
                        ends_at(&p, edges, j),
                        Some(ab),
                        true,
                        true,
                        false,
                    )
                    .value
                        != 0.0
                })
            })
            .count();
        let (fn_, fo) = (ep + cn, ep + co);
        println!(
            "LINE tau {tau:.3e} dF_new {:.4e} ratio_new {:.4} dF_old {:.4e} ratio_old {:.4} d_paper {:.4e} d_corr_first {:.4e} d_corr_second {:.4e} sn_active {n_act}",
            fn_ - f0_new,
            (fn_ - f0_new) / (-tau * rate_new),
            fo - f0_old,
            (fo - f0_old) / (-tau * rate_old),
            ep - e_paper,
            co - en_old,
            (cn - co) - (en_new - en_old),
        );
    }
}

fn ends_at(p: &[[f64; 3]], edges: &[[u32; 2]], k: usize) -> [P; 2] {
    [p[edges[k][0] as usize], p[edges[k][1] as usize]]
}

/// For the checkpoints `RC_CKPT/step<k>.obj`, `k = RC_FROM ..= RC_TO`: every second-neighbour
/// pair of the near-pair band with `E^n/E^p ≥ RC_MINRATIO` (default 0.5; its geometry and
/// correction), the closest disjoint pair, and the CCD contact (gap 0.1) of the straight motion
/// from step `k` to step `k + 1` (`τ ∈ [0, 1]`, 1 = the whole step), with the pair.
#[test]
#[ignore]
fn track_probe() {
    use ddg_repulsive_curves::step::line_search::closest_pair_within;
    let Some(dir) = env("RC_CKPT") else {
        eprintln!("RC_CKPT not set; skipped");
        return;
    };
    let from: usize = env("RC_FROM").and_then(|s| s.parse().ok()).unwrap_or(0);
    let to: usize = env("RC_TO").and_then(|s| s.parse().ok()).unwrap_or(30);
    let min_ratio: f64 = env("RC_MINRATIO")
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.5);
    let pool = make_pool(Some(4)).unwrap();
    let e = Exponents::new(2.0, 4.5).unwrap();
    let tables = InnerTables::new(e);
    let new = NearPairOptions::default();
    for k in from..=to {
        let Ok((net, pos)) = read_curve_obj(format!("{dir}/step{k}.obj")) else {
            break;
        };
        let edges = net.edges();
        let geom = net.edge_geometry(&pos);
        let lbar = geom.length.iter().sum::<f64>() / geom.length.len() as f64;
        let ends = |k: usize| [pos[edges[k][0] as usize], pos[edges[k][1] as usize]];
        let mut line = format!("STEP {k}");
        for (i, j) in near_pairs(&net, &pos, &new) {
            let (i, j) = (i as usize, j as usize);
            let Some((a, b)) = joining_ends(&net, i, j) else {
                continue;
            };
            let (p, q) = (ends(i), ends(j));
            let cl = closest_points(p[0], p[1], q[0], q[1]);
            let sw = |g: &[[f64; 3]; 4]| [g[2], g[3], g[0], g[1]];
            let (gs, gd) = (sw(&cl.grad_s), sw(&cl.grad_d));
            let (n_ij, _) = ordered_near_energy(
                &tables,
                p,
                q,
                (cl.t, &cl.grad_t),
                (cl.d, &cl.grad_d),
                4,
                false,
            );
            let (n_ji, _) = ordered_near_energy(&tables, q, p, (cl.s, &gs), (cl.d, &gd), 4, false);
            let ep = paper_pair(p, q, e).0 + paper_pair(q, p, e).0;
            let rho = (n_ij + n_ji) / ep;
            if rho < min_ratio {
                continue;
            }
            let pc =
                pair_correction_of_kind(&tables, e, &new, p, q, Some((a, b)), true, true, false);
            let (pa, pfar, qb, qfar) = (p[a], p[1 - a], q[b], q[1 - b]);
            let tk = sub(qb, pa);
            line.push_str(&format!(
                " | SN {i}-{j} r {:.3} d {:.3} turn {:.0}/{:.0} En/Ep {:.2} corr {:.3e}",
                pc.ratio,
                cl.d / lbar,
                angle(sub(pa, pfar), tk),
                angle(tk, sub(qfar, qb)),
                rho,
                pc.value
            ));
        }
        if let Some((a, b, d)) = closest_pair_within(&net, &pos, lbar) {
            line.push_str(&format!(
                " | closest {a}-{b} d {:.3e} sn {}",
                d / lbar,
                joining_ends(&net, a as usize, b as usize).is_some()
            ));
        }
        if let Ok((_, next)) = read_curve_obj(format!("{dir}/step{}.obj", k + 1)) {
            let vel: Vec<[f64; 3]> = pos.iter().zip(&next).map(|(x, y)| sub(*y, *x)).collect();
            let maxv = vel.iter().map(|v| norm(*v)).fold(0.0f64, f64::max);
            let rep =
                first_contact(&pool, &net, &pos, &vel, 1.0, &[], &CcdOptions::default()).unwrap();
            line.push_str(&format!(" | max move {:.3}", maxv / lbar));
            if let Some(Contact {
                tau,
                primitive: ContactPrimitive::EdgeEdge { a, b },
            }) = rep.contact
            {
                let (p, q) = (ends(a as usize), ends(b as usize));
                let cl = closest_points(p[0], p[1], q[0], q[1]);
                line.push_str(&format!(
                    " | step contact at {tau:.3} {a}-{b} d0 {:.3e} sn {}",
                    cl.d / lbar,
                    joining_ends(&net, a as usize, b as usize).is_some()
                ));
            }
        }
        println!("{line}");
    }
}

/// For every OBJ listed in `RC_LIST`: the second-neighbour pairs of the near-pair band, how many
/// have a nonzero one-sided correction (Q10 default), their summed correction, and the pair with
/// the largest `E^n/E^p` (geometry: turning angles at the ends of `K`, the angle between `I`
/// oriented toward `K` and `J` oriented away from it — 0° for a Z-fold, 180° for a hairpin —
/// lengths relative to the mean edge length, `r`). One `SNC` line per curve.
#[test]
#[ignore]
fn sn_fold_census() {
    let Some(list) = env("RC_LIST") else {
        eprintln!("RC_LIST not set; skipped");
        return;
    };
    let e = Exponents::new(2.0, 4.5).unwrap();
    let tables = InnerTables::new(e);
    let new = NearPairOptions::default();
    for path in std::fs::read_to_string(&list).unwrap().lines() {
        let path = path.trim();
        if path.is_empty() {
            continue;
        }
        let Ok((net, pos)) = read_curve_obj(path) else {
            println!("SNC {path} unreadable");
            continue;
        };
        let edges = net.edges();
        let geom = net.edge_geometry(&pos);
        let lbar = geom.length.iter().sum::<f64>() / geom.length.len() as f64;
        let ends = |k: usize| [pos[edges[k][0] as usize], pos[edges[k][1] as usize]];
        let (mut band, mut active, mut total) = (0usize, 0usize, 0.0f64);
        let mut best: Option<(f64, String)> = None;
        for (i, j) in near_pairs(&net, &pos, &new) {
            let (i, j) = (i as usize, j as usize);
            let Some((a, b)) = joining_ends(&net, i, j) else {
                continue;
            };
            band += 1;
            let (p, q) = (ends(i), ends(j));
            let pc =
                pair_correction_of_kind(&tables, e, &new, p, q, Some((a, b)), true, true, false);
            if pc.value == 0.0 {
                continue;
            }
            active += 1;
            total += pc.value;
            let cl = closest_points(p[0], p[1], q[0], q[1]);
            let sw = |g: &[[f64; 3]; 4]| [g[2], g[3], g[0], g[1]];
            let (gs, gd) = (sw(&cl.grad_s), sw(&cl.grad_d));
            let (n_ij, _) = ordered_near_energy(
                &tables,
                p,
                q,
                (cl.t, &cl.grad_t),
                (cl.d, &cl.grad_d),
                4,
                false,
            );
            let (n_ji, _) = ordered_near_energy(&tables, q, p, (cl.s, &gs), (cl.d, &gd), 4, false);
            let ep = paper_pair(p, q, e).0 + paper_pair(q, p, e).0;
            let rho = (n_ij + n_ji) / ep;
            let (pa, pfar, qb, qfar) = (p[a], p[1 - a], q[b], q[1 - b]);
            let tk = sub(qb, pa);
            let s = format!(
                "pair {i}-{j} En/Ep {rho:.3e} corr {:.3e} r {:.3} d {:.3e} turn {:.1}/{:.1} angIJ {:.1} lI {:.3} lK {:.3} lJ {:.3}",
                pc.value,
                pc.ratio,
                cl.d / lbar,
                angle(sub(pa, pfar), tk),
                angle(tk, sub(qfar, qb)),
                angle(sub(pa, pfar), sub(qfar, qb)),
                geom.length[i] / lbar,
                norm(tk) / lbar,
                geom.length[j] / lbar
            );
            if best.as_ref().is_none_or(|(r, _)| rho > *r) {
                best = Some((rho, s));
            }
        }
        println!(
            "SNC {path} band {band} active {active} corr {total:.4e} max {}",
            best.map_or("none".into(), |b| b.1)
        );
    }
}

/// At `RC_CURVE`: every pair of the near-pair band (Q10 default options) with its kind
/// (second neighbours or not), `r`, `d/ℓ̄`, the edge lengths, and for second neighbours `E^n/E^p`
/// and the turning angles at the ends of `K`. One `BAND` line per pair.
#[test]
#[ignore]
fn band_probe() {
    let Some(path) = env("RC_CURVE") else {
        eprintln!("RC_CURVE not set; skipped");
        return;
    };
    let (net, pos) = read_curve_obj(&path).unwrap();
    let e = Exponents::new(2.0, 4.5).unwrap();
    let tables = InnerTables::new(e);
    let new = NearPairOptions::default();
    let edges = net.edges();
    let geom = net.edge_geometry(&pos);
    let lbar = geom.length.iter().sum::<f64>() / geom.length.len() as f64;
    let ends = |k: usize| [pos[edges[k][0] as usize], pos[edges[k][1] as usize]];
    for (i, j) in near_pairs(&net, &pos, &new) {
        let (i, j) = (i as usize, j as usize);
        let (p, q) = (ends(i), ends(j));
        let cl = closest_points(p[0], p[1], q[0], q[1]);
        let lm = (0.5 * (geom.length[i].powi(2) + geom.length[j].powi(2))).sqrt();
        let mut extra = String::new();
        if let Some((a, b)) = joining_ends(&net, i, j) {
            let sw = |g: &[[f64; 3]; 4]| [g[2], g[3], g[0], g[1]];
            let (gs, gd) = (sw(&cl.grad_s), sw(&cl.grad_d));
            let (n_ij, _) = ordered_near_energy(
                &tables,
                p,
                q,
                (cl.t, &cl.grad_t),
                (cl.d, &cl.grad_d),
                4,
                false,
            );
            let (n_ji, _) = ordered_near_energy(&tables, q, p, (cl.s, &gs), (cl.d, &gd), 4, false);
            let ep = paper_pair(p, q, e).0 + paper_pair(q, p, e).0;
            let (pa, pfar, qb, qfar) = (p[a], p[1 - a], q[b], q[1 - b]);
            let tk = sub(qb, pa);
            extra = format!(
                " sn En/Ep {:.3} turn {:.0}/{:.0} lK {:.3}",
                (n_ij + n_ji) / ep,
                angle(sub(pa, pfar), tk),
                angle(tk, sub(qfar, qb)),
                norm(tk) / lbar
            );
        }
        println!(
            "BAND {i}-{j} r {:.3} d {:.3e} lI {:.3} lJ {:.3}{extra}",
            cl.d / lm,
            cl.d / lbar,
            geom.length[i] / lbar,
            geom.length[j] / lbar
        );
    }
}
