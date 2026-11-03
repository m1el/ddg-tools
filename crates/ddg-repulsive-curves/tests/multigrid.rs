//! Tests of work package WP6 (`SPEC.md` §8.4, §13.4, D17, D25): coarsening and prolongation,
//! the levels with their constraints, the V-cycle with the full-multigrid start and the FCG
//! fallback, the multigrid saddle solver against the dense one (projected gradient and chord
//! projection), iterations under refinement, determinism, and the solver in a session.
//!
//! The release-only studies are `#[ignore]`d; run them alone:
//! `cargo test -p ddg-repulsive-curves --release --test multigrid -- --ignored --nocapture <name>`.

#![allow(clippy::needless_range_loop)]

use std::sync::Arc;
use std::time::Instant;

use ddg_core::sparse::SparseMatrix;
use ddg_repulsive_curves::constraints::{
    Constraint, ConstraintSet, TangentSite, default_constraints,
};
use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::metric::{DofLayout, GramQuadrature, Metric};
use ddg_repulsive_curves::multigrid::vcycle::{self, CycleOptions, Hierarchy};
use ddg_repulsive_curves::multigrid::{
    CurveHierarchy, HierarchyInput, MgOptions, MultigridSaddle, MultigridSolvers, ProlongationRow,
    coarsen, forced_black, mark_black,
};
use ddg_repulsive_curves::optimize::mean_free_rows;
use ddg_repulsive_curves::solve::{
    ProjectionOptions, ProjectionStatus, SaddleSolver, Tolerance, dense_sobolev_solver,
    project_constraints, projected_gradient,
};
use ddg_repulsive_curves::spatial::bct::BctOptions;
use ddg_repulsive_curves::subdivide::subdivide;
use ddg_repulsive_curves::{
    Axis, Backend, BhOptions, CurveNetwork, EnergyEval, Exponents, FlowOptions, LineSearch,
    RcError, Session, Solver, StartRule, read_curve_obj,
};

type P = [f64; 3];

// ---------------------------------------------------------------- helpers

fn data(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/repulsive-curves")
        .join(rel)
}

fn have_data() -> bool {
    let ok = data("FreedmanHeWang/freedman-he-wang-600.obj").exists();
    if !ok {
        eprintln!("data/repulsive-curves not present; skipping");
    }
    ok
}

fn freedman(size: usize) -> (CurveNetwork, Vec<P>) {
    read_curve_obj(data(&format!("FreedmanHeWang/freedman-he-wang-{size}.obj"))).unwrap()
}

/// Freedman 2400 refined `times` times by midpoint subdivision (`times` = 0, 1, 2).
fn freedman_refined(base: usize, times: usize) -> (CurveNetwork, Vec<P>) {
    let (mut net, mut pos) = freedman(base);
    for _ in 0..times {
        let (n, p, _) = subdivide(&net, &pos).unwrap();
        net = n;
        pos = p;
    }
    (net, pos)
}

fn closed_edges(n: usize) -> Vec<[u32; 2]> {
    (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect()
}

fn trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [
                t.sin() + 2.0 * (2.0 * t).sin(),
                t.cos() - 2.0 * (2.0 * t).cos(),
                -(3.0 * t).sin(),
            ]
        })
        .collect();
    (CurveNetwork::new(n, closed_edges(n)).unwrap(), pos)
}

/// A network with junctions, open ends and several components: a tripod (three arms from a
/// degree-3 centre, open ends), a loop of `m` vertices, a small loop of 5 vertices (frozen) and
/// an open polyline.
fn mixed_network(m: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos: Vec<P> = vec![[0.0, 0.0, 0.0]];
    let mut edges = Vec::new();
    for arm in 0..3 {
        let phi = std::f64::consts::TAU * arm as f64 / 3.0;
        let mut prev = 0u32;
        for k in 1..=m {
            let r = k as f64 / m as f64;
            pos.push([
                r * phi.cos(),
                r * phi.sin(),
                0.2 * r * r * (arm as f64 - 1.0),
            ]);
            let v = pos.len() as u32 - 1;
            edges.push([prev, v]);
            prev = v;
        }
    }
    let mut ring = |n: usize, c: P, r: f64, z: f64| {
        let base = pos.len() as u32;
        for k in 0..n {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            pos.push([
                c[0] + r * t.cos(),
                c[1] + r * t.sin(),
                z + 0.1 * (2.0 * t).sin(),
            ]);
        }
        for k in 0..n as u32 {
            edges.push([base + k, base + (k + 1) % n as u32]);
        }
    };
    ring(m, [0.0, 0.0, 0.6], 0.7, 0.6);
    ring(5, [2.5, 0.0, 0.0], 0.3, 0.0);
    let base = pos.len() as u32;
    for k in 0..m {
        let t = k as f64 / (m - 1) as f64;
        pos.push([-2.0 + 0.3 * (5.0 * t).sin(), -1.5 + 3.0 * t, -0.5]);
    }
    for k in 0..m as u32 - 1 {
        edges.push([base + k, base + k + 1]);
    }
    (CurveNetwork::new(pos.len(), edges).unwrap(), pos)
}

/// Two vertices joined by three paths of `m` interior vertices each (a "theta" graph; with
/// `m = 1` two of its white vertices would give two coarse edges between the same pair).
fn theta_graph(m: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos: Vec<P> = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]];
    let mut edges = Vec::new();
    for path in 0..3 {
        let phi = std::f64::consts::TAU * path as f64 / 3.0;
        let mut prev = 0u32;
        for k in 1..=m {
            let t = k as f64 / (m + 1) as f64;
            let r = (std::f64::consts::PI * t).sin();
            pos.push([r * phi.cos(), r * phi.sin(), 2.0 * t]);
            let v = pos.len() as u32 - 1;
            edges.push([prev, v]);
            prev = v;
        }
        edges.push([prev, 1]);
    }
    (CurveNetwork::new(pos.len(), edges).unwrap(), pos)
}

/// Deterministic pseudo-random numbers in `[-0.5, 0.5)`.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    }
}

fn rand_vec(n: usize, seed: u64) -> Vec<f64> {
    let mut r = Lcg(seed);
    (0..n).map(|_| r.next()).collect()
}

/// `[1.2e-3, ...]`.
fn es(x: &[f64]) -> String {
    let v: Vec<String> = x.iter().map(|r| format!("{r:.1e}")).collect();
    format!("[{}]", v.join(", "))
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn rel_err(a: &[f64], b: &[f64]) -> f64 {
    let d: Vec<f64> = a.iter().zip(b).map(|(x, y)| x - y).collect();
    norm(&d) / norm(b)
}

fn pool(n: usize) -> Arc<rayon::ThreadPool> {
    Arc::new(make_pool(Some(n)).unwrap())
}

fn mul(c: &SparseMatrix, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; c.nrows()];
    c.mul_vec(x, &mut y).unwrap();
    y
}

/// The pieces of one step: constraint set, Jacobian, differential.
struct Problem {
    net: CurveNetwork,
    pos: Vec<P>,
    ex: Exponents,
    set: ConstraintSet,
    c: SparseMatrix,
    diff: Vec<P>,
}

fn problem(net: CurveNetwork, pos: Vec<P>, ex: Exponents, cons: &[Constraint]) -> Problem {
    problem_planar(net, pos, ex, cons, None)
}

fn problem_planar(
    net: CurveNetwork,
    pos: Vec<P>,
    ex: Exponents,
    cons: &[Constraint],
    planar: Option<Axis>,
) -> Problem {
    let set = ConstraintSet::new(&net, &pos, cons, planar).unwrap();
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

fn mg_options(coarsest: usize, theta: f64) -> MgOptions {
    MgOptions {
        coarsest_vertices: coarsest,
        bct: BctOptions {
            theta_x: theta,
            theta_t: theta,
            ..BctOptions::default()
        },
        ..MgOptions::default()
    }
}

fn mg_solver(pr: &Problem, threads: usize, o: MgOptions) -> MultigridSaddle {
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
    MultigridSaddle::new(&p, &p, &input, o).unwrap()
}

fn level_sizes(s: &MultigridSaddle) -> Vec<usize> {
    s.hierarchy()
        .levels()
        .iter()
        .map(|l| l.network().n_vertices())
        .collect()
}

// ---------------------------------------------------------------- coarsening and prolongation

/// `J 𝟙 = 𝟙`, `Π` of a constant is the constant, the restriction is the transpose, and the
/// matrix form agrees with the matrix-free products (`rc-prolongation`).
#[test]
fn prolongation_preserves_constants_and_restriction_is_its_transpose() {
    for (net, _) in [trefoil(64), mixed_network(9), theta_graph(4)] {
        let n = net.n_vertices();
        let cg = coarsen(&net, &vec![false; n]).unwrap().unwrap();
        let nc = cg.n_coarse();
        assert!(nc < n);
        for d in [1, 2, 3] {
            let ones = vec![1.0; nc * d];
            let p = cg.prolong(&ones, d).unwrap();
            assert_eq!(p, vec![1.0; n * d], "J 1 = 1");
            let x = rand_vec(nc * d, 3);
            let y = rand_vec(n * d, 4);
            let jx = cg.prolong(&x, d).unwrap();
            let jty = cg.restrict(&y, d).unwrap();
            let a: f64 = jx.iter().zip(&y).map(|(p, q)| p * q).sum();
            let b: f64 = x.iter().zip(&jty).map(|(p, q)| p * q).sum();
            assert!((a - b).abs() <= 1e-13 * a.abs().max(1.0), "{a} vs {b}");
        }
        let j = cg.prolongation_matrix().unwrap();
        assert_eq!((j.nrows(), j.ncols()), (n, nc));
        let x = rand_vec(nc, 5);
        let mut jx = vec![0.0; n];
        j.mul_vec(&x, &mut jx).unwrap();
        assert_eq!(jx, cg.prolong(&x, 1).unwrap());
        let rows = mul(&j, &vec![1.0; nc]);
        assert!(rows.iter().all(|&r| r == 1.0));
        // Wrong lengths are errors.
        assert!(matches!(
            cg.prolong(&[1.0; 2], 3),
            Err(RcError::InvalidOption { .. })
        ));
        assert!(matches!(
            cg.restrict(&ones_len(n * 3 + 1), 3),
            Err(RcError::InvalidOption { .. })
        ));
        assert!(matches!(
            cg.prolong(&[], 0),
            Err(RcError::InvalidOption { .. })
        ));
    }
}

fn ones_len(n: usize) -> Vec<f64> {
    vec![1.0; n]
}

/// The invariants of `rc-coarsen` on one coarsening.
fn check_coarsening(net: &CurveNetwork, pos: &[P], forced: &[bool]) -> usize {
    let n = net.n_vertices();
    let Some(cg) = coarsen(net, forced).unwrap() else {
        return 0;
    };
    let f2c = cg.fine_to_coarse();
    let black = |v: u32| f2c[v as usize].is_some();
    let edges = net.edges();
    let mut whites = 0;
    for v in 0..n as u32 {
        if forced[v as usize] || net.degree(v) != 2 {
            assert!(
                black(v),
                "vertex {v} (deg {}) must stay black",
                net.degree(v)
            );
        }
        match cg.rows()[v as usize] {
            ProlongationRow::Black(c) => {
                assert!(black(v));
                assert_eq!(cg.coarse_to_fine()[c as usize], v);
            }
            ProlongationRow::White(a, b) => {
                whites += 1;
                assert!(!black(v));
                assert_eq!(net.degree(v), 2);
                assert_ne!(a, b);
                // Both neighbours are black and they are the coarse a, b.
                let nb: Vec<u32> = net
                    .incident_edges(v)
                    .iter()
                    .map(|&e| {
                        let [p, q] = edges[e as usize];
                        if p == v { q } else { p }
                    })
                    .collect();
                for &w in &nb {
                    assert!(black(w), "white {v} has a white neighbour {w}");
                }
                let mut got = [
                    cg.coarse_to_fine()[a as usize],
                    cg.coarse_to_fine()[b as usize],
                ];
                got.sort();
                let mut want = [nb[0], nb[1]];
                want.sort();
                assert_eq!(got, want);
            }
        }
    }
    // No two adjacent whites (already implied), the coarse network is valid (CurveNetwork::new
    // rejects duplicate edges and self-loops), fine positions are kept.
    let cnet = cg.coarse();
    assert_eq!(cnet.n_vertices(), n - whites);
    assert_eq!(cnet.n_edges(), net.n_edges() - whites);
    assert_eq!(cnet.components().len(), net.components().len());
    let cpos = cg.coarse_positions(pos).unwrap();
    for (c, &v) in cg.coarse_to_fine().iter().enumerate() {
        assert_eq!(cpos[c], pos[v as usize]);
        assert_eq!(cnet.degree(c as u32), net.degree(v));
    }
    // Every fine edge maps to a coarse edge that contains it.
    for (e, &ce) in cg.edge_to_coarse().iter().enumerate() {
        assert!(cg.coarse_edge_parts()[ce as usize].contains(&(e as u32)));
    }
    // Small loops stay all black.
    for comp in net.components() {
        if comp.closed && comp.vertices.len() < 6 {
            assert!(comp.vertices.iter().all(|&v| black(v)));
        }
    }
    whites
}

#[test]
fn coarsening_invariants_on_junctions_open_ends_and_components() {
    let (net, pos) = mixed_network(11);
    let n = net.n_vertices();
    let mut forced = vec![false; n];
    forced[5] = true;
    forced[40] = true;
    // Repeated coarsening down to nothing removable.
    let (mut net, mut pos, mut forced) = (net, pos, forced);
    let mut rounds = 0;
    loop {
        let w = check_coarsening(&net, &pos, &forced);
        if w == 0 {
            break;
        }
        let cg = coarsen(&net, &forced).unwrap().unwrap();
        let f: Vec<bool> = cg
            .coarse_to_fine()
            .iter()
            .map(|&v| forced[v as usize])
            .collect();
        pos = cg.coarse_positions(&pos).unwrap();
        net = cg.coarse().clone();
        forced = f;
        rounds += 1;
        assert!(rounds < 20);
    }
    assert!(rounds >= 2, "{rounds}");
    // Forced vertices survive every round.
    assert_eq!(forced.iter().filter(|&&f| f).count(), 2);

    // Loops: the first vertex is black, whites alternate, loops of < 6 vertices are frozen.
    for m in [4, 5, 6, 7, 8, 13, 64] {
        let (net, pos) = trefoil(m);
        if m == 6 {
            // A single loop of 6 coarsens to a triangle: no interacting pair, an error (the
            // hierarchy stops coarsening there).
            assert!(matches!(
                coarsen(&net, &[false; 6]),
                Err(RcError::NoInteractingPairs { .. })
            ));
            continue;
        }
        let w = check_coarsening(&net, &pos, &vec![false; m]);
        if m < 6 {
            assert_eq!(w, 0, "loop of {m}");
        } else {
            assert!(w >= m / 3, "loop of {m}: {w} whites");
        }
    }
    let b = mark_black(&trefoil(64).0, &[false; 64]).unwrap();
    assert!(b[0]);
    assert_eq!(b.iter().filter(|&&x| !x).count(), 32);

    // Theta graphs: with one interior vertex per path, the duplicate-edge check keeps all but
    // one path's vertex black; longer paths coarsen normally.
    for m in [1, 2, 3, 6] {
        let (net, pos) = theta_graph(m);
        let w = check_coarsening(&net, &pos, &vec![false; net.n_vertices()]);
        eprintln!("theta graph m = {m}: {w} whites");
        assert!(w >= 1);
    }

    // Bad input is an error.
    assert!(matches!(
        mark_black(&trefoil(10).0, &[false; 3]),
        Err(RcError::InvalidOption { .. })
    ));
}

/// Coarse constraints (§8.4): points and vertex tangents survive on their (black) vertex, edge
/// tangents move to the containing coarse edge, edge-length subsets keep the coarse edges whose
/// fine edges were all constrained, surface rows stay on surviving vertices.
#[test]
fn coarse_constraints_follow_the_vertices_and_edges() {
    let (net, pos) = trefoil(40);
    let cons = vec![
        Constraint::Barycenter { target: None },
        Constraint::EdgeLengths {
            edges: Some((0..20).collect()),
            targets: None,
        },
        Constraint::Point {
            vertex: 7,
            target: None,
        },
        Constraint::Tangent {
            site: TangentSite::Vertex(13),
            direction: None,
        },
        Constraint::Tangent {
            site: TangentSite::Edge(30),
            direction: None,
        },
    ];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    let forced = forced_black(net.n_vertices(), Some(&set));
    assert!(forced[7] && forced[13]);
    assert_eq!(forced.iter().filter(|&&f| f).count(), 2);
    let cg = coarsen(&net, &forced).unwrap().unwrap();
    let f2c = cg.fine_to_coarse();
    let coarse = cg.coarse_constraints(&set);
    assert_eq!(coarse.len(), 5);
    match &coarse[1] {
        Constraint::EdgeLengths {
            edges: Some(list), ..
        } => {
            for &ce in list {
                for &f in &cg.coarse_edge_parts()[ce as usize] {
                    assert!(f < 20);
                }
            }
            // Every coarse edge made only of listed fine edges is listed.
            let all: Vec<u32> = (0..cg.coarse().n_edges() as u32)
                .filter(|&ce| cg.coarse_edge_parts()[ce as usize].iter().all(|&f| f < 20))
                .collect();
            assert_eq!(list, &all);
        }
        c => panic!("{c:?}"),
    }
    match &coarse[2] {
        Constraint::Point { vertex, .. } => assert_eq!(Some(*vertex), f2c[7]),
        c => panic!("{c:?}"),
    }
    match &coarse[3] {
        Constraint::Tangent {
            site: TangentSite::Vertex(v),
            ..
        } => assert_eq!(Some(*v), f2c[13]),
        c => panic!("{c:?}"),
    }
    match &coarse[4] {
        Constraint::Tangent {
            site: TangentSite::Edge(e),
            ..
        } => assert_eq!(*e, cg.edge_to_coarse()[30]),
        c => panic!("{c:?}"),
    }
    let cpos = cg.coarse_positions(&pos).unwrap();
    let cset = ConstraintSet::new(cg.coarse(), &cpos, &coarse, None).unwrap();
    assert!(cset.n_active_rows() > 0);
}

// ---------------------------------------------------------------- solves against dense

/// The projected gradient (Eq. 22 through Eq. 30) by multigrid against the dense LU, with the
/// exact product (θ = 0, so only the multigrid tolerance is left) and with the default θ.
fn compare_gradient(pr: &Problem, o: MgOptions, bound: f64) -> (Vec<usize>, usize, f64) {
    let p = pool(16);
    let mut dense = dense_sobolev_solver(
        p.clone(),
        &pr.net,
        &pr.pos,
        pr.ex,
        GramQuadrature::Midpoint,
        &pr.set,
    )
    .unwrap();
    let layout = pr.set.layout();
    let gd = projected_gradient(&mut dense, layout, &pr.diff).unwrap();
    let xd = layout.gather(&gd.direction);
    let mut mg = mg_solver(pr, 16, o);
    let gm = projected_gradient(&mut mg, layout, &pr.diff).unwrap();
    let xm = layout.gather(&gm.direction);
    let e = rel_err(&xm, &xd);
    let cx = mul(&pr.c, &xm);
    let rep = &mg.last_reports()[0];
    eprintln!(
        "{} vertices, levels {:?}, theta {}: {} cycles (fallback {}), residuals {:?}, saddle residual {:.2e}, |g_mg - g_d|/|g_d| = {e:.2e}, |C g|/|g| = {:.2e}",
        pr.net.n_vertices(),
        level_sizes(&mg),
        o.bct.theta_x,
        gm.linear.iterations,
        gm.linear.fallback_used,
        rep.residuals,
        gm.linear.relative_residual,
        norm(&cx) / norm(&xm)
    );
    assert!(rep.converged, "{rep:?}");
    assert!(rep.final_residual() <= o.tolerance);
    assert!(norm(&cx) <= 1e-10 * norm(&xm), "C g = {:e}", norm(&cx));
    assert!(e <= bound, "{e:e} > {bound:e}");
    assert_eq!(gm.linear.levels, mg.hierarchy().levels().len());
    (level_sizes(&mg), gm.linear.iterations, e)
}

#[test]
fn two_level_and_multilevel_solves_reproduce_the_dense_projected_gradient() {
    let (net, pos, ex) = if have_data() {
        let (n, x) = freedman(600);
        (n, x, Exponents::new(2.0, 4.5).unwrap())
    } else {
        let (n, x) = trefoil(600);
        (n, x, Exponents::default())
    };
    let pr = problem(net, pos, ex, &default_constraints());
    // Two levels (592 → 296), exact product.
    let (sizes, _, _) = compare_gradient(&pr, mg_options(512, 0.0), 1e-2);
    assert_eq!(sizes.len(), 2);
    // Many levels (down to ≤ 40 vertices), exact product.
    let (sizes, _, _) = compare_gradient(&pr, mg_options(40, 0.0), 1e-2);
    assert!(sizes.len() >= 4, "{sizes:?}");
    // Tighter tolerance: the error follows it.
    let o = MgOptions {
        tolerance: 1e-8,
        ..mg_options(40, 0.0)
    };
    compare_gradient(&pr, o, 1e-6);
    // The default block cluster tree (θ = 0.25): its own error (§17e: 0.15 % here) on top.
    compare_gradient(&pr, mg_options(40, 0.25), 2e-2);
}

/// Junctions, open ends, several components, pinned vertices, tangents, planar mode and the
/// mean-free rows of a soft problem: the multigrid solve agrees with the dense one.
#[test]
fn multigrid_handles_junctions_components_pins_tangents_planar_and_soft_rows() {
    let ex = Exponents::default();
    let tight = MgOptions {
        tolerance: 1e-9,
        max_cycles: 60,
        ..mg_options(16, 0.0)
    };
    // Mixed network with pins and a vertex tangent.
    let (net, pos) = mixed_network(17);
    let cons = vec![
        Constraint::Barycenter { target: None },
        Constraint::EdgeLengths {
            edges: None,
            targets: None,
        },
        Constraint::Point {
            vertex: 9,
            target: None,
        },
        Constraint::Tangent {
            site: TangentSite::Vertex(30),
            direction: None,
        },
        Constraint::Tangent {
            site: TangentSite::Edge(60),
            direction: None,
        },
    ];
    let pr = problem(net, pos, ex, &cons);
    let (sizes, _, _) = compare_gradient(&pr, tight, 1e-6);
    assert!(sizes.len() >= 3, "{sizes:?}");

    // Planar mode (two coordinates per vertex).
    let n = 200;
    let pos: Vec<P> = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            let r = 1.0 + 0.3 * (3.0 * t).cos();
            [r * t.cos(), r * t.sin(), 0.5]
        })
        .collect();
    let net = CurveNetwork::new(n, closed_edges(n)).unwrap();
    let pr = problem_planar(net, pos, ex, &default_constraints(), Some(Axis::Z));
    assert_eq!(pr.set.layout().dim(), 2);
    compare_gradient(&pr, tight, 1e-6);

    // Soft problem: the mean-free rows, coarse levels with their own mean-free rows.
    let (net, pos) = trefoil(300);
    let layout = DofLayout::new(None);
    let c = mean_free_rows(&net, &pos, layout).unwrap();
    let p = pool(8);
    let (_, diff) = energy_and_differential(&p, &net, &pos, ex);
    let input = HierarchyInput {
        network: &net,
        positions: &pos,
        exponents: ex,
        quadrature: GramQuadrature::Midpoint,
        layout,
        jacobian: &c,
        constraints: None,
    };
    let mut mg = MultigridSaddle::new(&p, &p, &input, tight).unwrap();
    assert!(mg.hierarchy().levels().len() >= 3);
    let gm = projected_gradient(&mut mg, layout, &diff).unwrap();
    let gram = ddg_repulsive_curves::gram(
        &net,
        &pos,
        ex,
        &Metric::Sobolev {
            quadrature: GramQuadrature::Midpoint,
        },
        Some(8),
    )
    .unwrap();
    let mut dense = ddg_repulsive_curves::DenseSaddle::new(p.clone(), &gram.a, &c, layout).unwrap();
    let gd = projected_gradient(&mut dense, layout, &diff).unwrap();
    let e = rel_err(&layout.gather(&gm.direction), &layout.gather(&gd.direction));
    eprintln!("soft rows: {} cycles, error {e:.2e}", gm.linear.iterations);
    assert!(e <= 1e-6, "{e:e}");
}

/// A general right side `[a; b]` (Eq. 33) and the multipliers, against the dense LU.
#[test]
fn general_right_side_and_multipliers_match_the_dense_solve() {
    let (net, pos) = trefoil(240);
    let pr = problem(net, pos, Exponents::default(), &default_constraints());
    let p = pool(8);
    let mut dense = dense_sobolev_solver(
        p.clone(),
        &pr.net,
        &pr.pos,
        pr.ex,
        GramQuadrature::Midpoint,
        &pr.set,
    )
    .unwrap();
    let o = MgOptions {
        tolerance: 1e-10,
        max_cycles: 60,
        ..mg_options(30, 0.0)
    };
    let mut mg = mg_solver(&pr, 8, o);
    let a = rand_vec(mg.n_dofs(), 11);
    let b: Vec<f64> = rand_vec(mg.n_constraints(), 12)
        .iter()
        .map(|v| 1e-3 * v)
        .collect();
    let sd = dense.solve(&a, &b).unwrap();
    let sm = mg.solve(&a, &b).unwrap();
    eprintln!(
        "Eq. 33: {} cycles, saddle residual {:.2e}, x error {:.2e}, lambda error {:.2e}",
        sm.report.iterations,
        sm.report.relative_residual,
        rel_err(&sm.x, &sd.x),
        rel_err(&sm.lambda, &sd.lambda)
    );
    assert!(rel_err(&sm.x, &sd.x) <= 1e-7);
    assert!(rel_err(&sm.lambda, &sd.lambda) <= 1e-6);
    // C x = b to rounding.
    let cx = mul(&pr.c, &sm.x);
    let d: Vec<f64> = cx.iter().zip(&b).map(|(p, q)| p - q).collect();
    assert!(
        norm(&d) <= 1e-12 * norm(&b).max(norm(&sm.x)),
        "{:e}",
        norm(&d)
    );
    assert!(sm.report.relative_residual <= 1e-8);
    assert_eq!(sm.report.solver, "multigrid");
    // Wrong lengths are errors.
    assert!(matches!(
        mg.solve(&a[1..], &b),
        Err(RcError::InvalidOption { .. })
    ));
}

/// The chord projection of Eq. 23 through multigrid (one solve per iteration) converges like the
/// dense chord method on a perturbed curve.
#[test]
fn constraint_projection_through_multigrid_converges_like_the_dense_chord_method() {
    let (net, pos, ex) = if have_data() {
        let (n, x) = freedman(600);
        (n, x, Exponents::new(2.0, 4.5).unwrap())
    } else {
        let (n, x) = trefoil(600);
        (n, x, Exponents::default())
    };
    let pr = problem(net, pos, ex, &default_constraints());
    let layout = pr.set.layout();
    let p = pool(16);
    let mut dense = dense_sobolev_solver(
        p.clone(),
        &pr.net,
        &pr.pos,
        pr.ex,
        GramQuadrature::Midpoint,
        &pr.set,
    )
    .unwrap();
    let g = projected_gradient(&mut dense, layout, &pr.diff).unwrap();
    let l0 = pr.net.mean_edge_length(&pr.pos);
    // A trial step that moves the fastest vertex by half an edge length, plus a random
    // perturbation of 1 % of an edge length.
    let gmax = g
        .direction
        .iter()
        .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
        .fold(0.0, f64::max);
    let tau = 0.5 * l0 / gmax;
    let noise = rand_vec(3 * pr.pos.len(), 77);
    let trial: Vec<P> = pr
        .pos
        .iter()
        .zip(&g.direction)
        .enumerate()
        .map(|(v, (x, d))| {
            let mut y = [0.0; 3];
            for k in 0..3 {
                y[k] = x[k] - tau * d[k] + 0.01 * l0 * noise[3 * v + k];
            }
            y
        })
        .collect();
    let opts = ProjectionOptions {
        tolerance: Tolerance::MeanEdgeLength(1e-4),
        max_iterations: 8,
    };
    let pd = project_constraints(&mut dense, &pr.set, &pr.net, &trial, &opts, l0).unwrap();
    for (name, o) in [
        ("multigrid exact product", mg_options(40, 0.0)),
        ("multigrid default", mg_options(40, 0.25)),
    ] {
        let mut mg = mg_solver(&pr, 16, o);
        let pm = project_constraints(&mut mg, &pr.set, &pr.net, &trial, &opts, l0).unwrap();
        eprintln!(
            "{name}: dense {:?} {} its history {}; multigrid {:?} {} its history {} (cycles {:?})",
            pd.status,
            pd.iterations,
            es(&pd.history),
            pm.status,
            pm.iterations,
            es(&pm.history),
            pm.linear.iter().map(|l| l.iterations).collect::<Vec<_>>()
        );
        assert_eq!(pd.status, ProjectionStatus::Converged);
        assert_eq!(pm.status, ProjectionStatus::Converged);
        assert!(pm.iterations <= pd.iterations + 1);
        // Linear convergence of the chord method: every iterate at least 10x better.
        for w in pm.history.windows(2) {
            assert!(w[1] <= 0.1 * w[0], "{:?}", pm.history);
        }
        // The projected curves agree to the tolerance.
        let a: Vec<f64> = pm.positions.iter().flatten().copied().collect();
        let b: Vec<f64> = pd.positions.iter().flatten().copied().collect();
        let d = a
            .iter()
            .zip(&b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        assert!(d <= 1e-2 * l0, "{d:e}");
    }
}

// ---------------------------------------------------------------- refinement, fallback

/// Iterations under refinement on a small family (debug-friendly): trefoils of 256–2048
/// vertices, coarsest 32. The Freedman table is the ignored test below.
#[test]
fn iterations_stay_flat_under_refinement_small() {
    let mut counts = Vec::new();
    for n in [256, 512, 1024] {
        let (net, pos) = trefoil(n);
        let pr = problem(net, pos, Exponents::default(), &default_constraints());
        let mut mg = mg_solver(&pr, 16, mg_options(32, 0.25));
        let g = projected_gradient(&mut mg, pr.set.layout(), &pr.diff).unwrap();
        let rep = &mg.last_reports()[0];
        eprintln!(
            "trefoil {n}: levels {:?} cycles {} residuals {}",
            level_sizes(&mg),
            g.linear.iterations,
            es(&rep.residuals)
        );
        assert!(rep.converged);
        counts.push(g.linear.iterations);
    }
    let (lo, hi) = (*counts.iter().min().unwrap(), *counts.iter().max().unwrap());
    assert!(hi <= lo + 3, "{counts:?}");
}

/// A hierarchy whose coarse correction is useless (the coarsest solve returns zero) and whose
/// smoother is a single CG step: the V-cycles stagnate and flexible CG takes over.
struct Crippled<'a>(&'a CurveHierarchy);

impl Hierarchy for Crippled<'_> {
    fn n_levels(&self) -> usize {
        self.0.n_levels()
    }
    fn len(&self, level: usize) -> usize {
        self.0.len(level)
    }
    fn apply(&self, level: usize, x: &[f64]) -> Result<Vec<f64>, RcError> {
        self.0.apply(level, x)
    }
    fn prolong(&self, level: usize, coarse: &[f64]) -> Result<Vec<f64>, RcError> {
        self.0.prolong(level, coarse)
    }
    fn restrict(&self, level: usize, fine: &[f64]) -> Result<Vec<f64>, RcError> {
        self.0.restrict(level, fine)
    }
    fn coarse_solve(&self, b: &[f64]) -> Result<Vec<f64>, RcError> {
        Ok(vec![0.0; b.len()])
    }
}

#[test]
fn stagnating_cycles_fall_back_to_flexible_cg() {
    let (net, pos) = trefoil(400);
    let pr = problem(net, pos, Exponents::default(), &default_constraints());
    let mg = mg_solver(&pr, 8, mg_options(50, 0.0));
    let h = mg.hierarchy();
    let lv = &h.levels()[0];
    let b = lv
        .projector()
        .project(&pr.set.layout().gather(&pr.diff))
        .unwrap();
    let o = CycleOptions {
        smoothing_steps: 1,
        initial_smoothing_steps: 1,
        tolerance: 1e-3,
        max_cycles: 200,
        flexible_cg: false,
    };
    let (x, rep) = vcycle::solve(&Crippled(h), &b, &o).unwrap();
    eprintln!(
        "crippled: cycles {} fcg {} fallback {} converged {} final {:.2e}",
        rep.cycles,
        rep.fcg_iterations,
        rep.fallback_used,
        rep.converged,
        rep.final_residual()
    );
    assert!(rep.fallback_used);
    assert!(rep.cycles >= 3);
    assert!(rep.fcg_iterations > 0);
    // FCG on the crippled preconditioner still converges (it is CG with a smoother as
    // preconditioner), within the budget.
    assert!(rep.converged, "{rep:?}");
    let hx = h.apply(0, &x).unwrap();
    let r: Vec<f64> = b.iter().zip(&hx).map(|(p, q)| p - q).collect();
    assert!(norm(&r) <= 1e-3 * norm(&b) * 1.0001);
    // The same problem with the real hierarchy needs no fallback.
    let (_, rep) = vcycle::solve(h, &b, &o).unwrap();
    assert!(!rep.fallback_used && rep.converged, "{rep:?}");

    // A natural case: the tangled Freedman 1200 with weak stationary cycles (2 + 2 smoothing,
    // seven levels) stagnates; FCG finishes the solve.
    if have_data() {
        let (net, pos) = freedman(1200);
        let pr = problem(
            net,
            pos,
            Exponents::new(2.0, 4.5).unwrap(),
            &default_constraints(),
        );
        let o = MgOptions {
            smoothing_steps: 2,
            initial_smoothing_steps: 4,
            flexible_cg: false,
            ..mg_options(40, 0.25)
        };
        let mut mg = mg_solver(&pr, 16, o);
        let g = projected_gradient(&mut mg, pr.set.layout(), &pr.diff).unwrap();
        let rep = &mg.last_reports()[0];
        eprintln!(
            "freedman 1200, 2 + 2 stationary: cycles {} fcg {} residuals {}",
            rep.cycles,
            rep.fcg_iterations,
            es(&rep.residuals)
        );
        assert!(g.linear.fallback_used && rep.fallback_used);
        assert!(rep.converged, "{rep:?}");
    }
}

// ---------------------------------------------------------------- determinism, errors

#[test]
fn bitwise_identical_at_1_7_and_64_threads() {
    let (net, pos) = mixed_network(60);
    let cons = vec![
        Constraint::Barycenter { target: None },
        Constraint::EdgeLengths {
            edges: None,
            targets: None,
        },
        Constraint::Point {
            vertex: 9,
            target: None,
        },
    ];
    let pr = problem(net, pos, Exponents::default(), &cons);
    let a = rand_vec(pr.c.ncols(), 1);
    let b: Vec<f64> = rand_vec(pr.c.nrows(), 2).iter().map(|v| 1e-3 * v).collect();
    let mut results = Vec::new();
    for t in [1, 7, 64] {
        let mut mg = mg_solver(&pr, t, mg_options(40, 0.25));
        let s = mg.solve(&a, &b).unwrap();
        let g = projected_gradient(&mut mg, pr.set.layout(), &pr.diff).unwrap();
        results.push((s.x, s.lambda, g.direction, s.report.iterations));
    }
    assert!(results[0].3 > 0);
    for r in &results[1..] {
        assert!(r.0 == results[0].0, "x differs");
        assert!(r.1 == results[0].1, "lambda differs");
        assert!(r.2 == results[0].2, "gradient differs");
        assert_eq!(r.3, results[0].3);
    }
}

/// The warm starts of the projection solves (`SPEC.md` §17t, `MgOptions::recycle_projections`)
/// meet the accuracy of a cold solve: the same linearized constraints to rounding, a solution
/// within the solver tolerance of the cold one, fewer cycles, and bitwise the same at 1 and 7
/// threads.
#[test]
fn recycled_projection_solves_meet_the_cold_accuracy_with_fewer_cycles() {
    let (net, pos) = trefoil(600);
    let pr = problem(net, pos, Exponents::default(), &default_constraints());
    let cold_opts = MgOptions {
        recycle_projections: false,
        ..mg_options(60, 0.25)
    };
    let warm_opts = mg_options(60, 0.25);
    assert!(warm_opts.recycle_projections);
    // Projection right sides of the kind a line search produces: b₂ ≈ b₁/4 plus a little.
    let b1: Vec<f64> = rand_vec(pr.c.nrows(), 21)
        .iter()
        .map(|v| 1e-3 * v)
        .collect();
    let noise = rand_vec(pr.c.nrows(), 22);
    let b2: Vec<f64> = b1
        .iter()
        .zip(&noise)
        .map(|(p, q)| 0.25 * p + 1e-5 * q)
        .collect();
    let zero = vec![0.0; pr.c.ncols()];
    let mut results = Vec::new();
    for threads in [1, 7] {
        let mut cold = mg_solver(&pr, threads, cold_opts);
        let mut warm = mg_solver(&pr, threads, warm_opts);
        let c1 = cold.solve(&zero, &b1).unwrap();
        let c2 = cold.solve(&zero, &b2).unwrap();
        let w1 = warm.solve(&zero, &b1).unwrap();
        let w2 = warm.solve(&zero, &b2).unwrap();
        // The first solve has nothing to start from: identical to the cold one.
        assert!(w1.x == c1.x, "first projection solve differs");
        let err = rel_err(&w2.x, &c2.x);
        eprintln!(
            "recycled projection: cycles {} (cold {}), x vs cold {:.2e}",
            w2.report.iterations, c2.report.iterations, err
        );
        assert!(err > 0.0, "the kept solve was not used");
        assert!(err <= 1e-2, "warm x differs from cold by {err:e}");
        assert!(w2.report.iterations < c2.report.iterations);
        let cx = mul(&pr.c, &w2.x);
        let d: Vec<f64> = cx.iter().zip(&b2).map(|(p, q)| p - q).collect();
        assert!(
            norm(&d) <= 1e-12 * norm(&b2).max(norm(&w2.x)),
            "C x - b {:e}",
            norm(&d)
        );
        // A direction solve (a ≠ 0) is not recycled: identical to the cold one.
        let gw = projected_gradient(&mut warm, pr.set.layout(), &pr.diff).unwrap();
        let gc = projected_gradient(&mut cold, pr.set.layout(), &pr.diff).unwrap();
        assert!(gw.direction == gc.direction);
        results.push(w2.x);
    }
    assert!(
        results[0] == results[1],
        "warm starts depend on the thread count"
    );
}

#[test]
fn bad_options_and_inputs_are_errors_not_panics() {
    let (net, pos) = trefoil(80);
    let pr = problem(net, pos, Exponents::default(), &default_constraints());
    let p = pool(2);
    let input = HierarchyInput {
        network: &pr.net,
        positions: &pr.pos,
        exponents: pr.ex,
        quadrature: GramQuadrature::Midpoint,
        layout: pr.set.layout(),
        jacobian: &pr.c,
        constraints: Some(&pr.set),
    };
    for o in [
        MgOptions {
            tolerance: 0.0,
            ..MgOptions::default()
        },
        MgOptions {
            tolerance: f64::NAN,
            ..MgOptions::default()
        },
        MgOptions {
            max_levels: 0,
            ..MgOptions::default()
        },
        MgOptions {
            coarsest_vertices: 0,
            ..MgOptions::default()
        },
        MgOptions {
            bct: BctOptions {
                theta_x: -1.0,
                ..BctOptions::default()
            },
            ..MgOptions::default()
        },
    ] {
        assert!(matches!(
            MultigridSaddle::new(&p, &p, &input, o),
            Err(RcError::InvalidOption { .. })
        ));
        assert!(MultigridSolvers::new(o).is_err());
    }
    // Wrong position count, wrong Jacobian width.
    let short = HierarchyInput {
        positions: &pr.pos[1..],
        ..input
    };
    assert!(matches!(
        MultigridSaddle::new(&p, &p, &short, MgOptions::default()),
        Err(RcError::PositionCountMismatch { .. })
    ));
    let c2 = mean_free_rows(&trefoil(20).0, &trefoil(20).1, DofLayout::new(None)).unwrap();
    let bad = HierarchyInput {
        jacobian: &c2,
        ..input
    };
    assert!(matches!(
        MultigridSaddle::new(&p, &p, &bad, MgOptions::default()),
        Err(RcError::InvalidOption { .. })
    ));
    // Redundant rows on the finest level are the caller's error.
    let mut dup = pr.c.triplets().collect::<Vec<_>>();
    let k = pr.c.nrows() as u32;
    dup.extend(
        pr.c.triplets()
            .filter(|t| t.0 == 0)
            .map(|(_, c, v)| (k, c, v))
            .collect::<Vec<_>>(),
    );
    let cdup = SparseMatrix::from_triplets(pr.c.nrows() + 1, pr.c.ncols(), &dup).unwrap();
    let red = HierarchyInput {
        jacobian: &cdup,
        constraints: None,
        ..input
    };
    assert!(matches!(
        MultigridSaddle::new(&p, &p, &red, MgOptions::default()),
        Err(RcError::RedundantConstraints { .. })
    ));
    // A wrong right side of the generic solver.
    let mg = mg_solver(&pr, 2, mg_options(20, 0.0));
    let o = CycleOptions {
        smoothing_steps: 3,
        initial_smoothing_steps: 3,
        tolerance: 1e-3,
        max_cycles: 5,
        flexible_cg: false,
    };
    assert!(matches!(
        vcycle::solve(mg.hierarchy(), &[1.0; 3], &o),
        Err(RcError::InvalidOption { .. })
    ));
    let bad_tol = CycleOptions {
        tolerance: -1.0,
        ..o
    };
    assert!(matches!(
        vcycle::solve(mg.hierarchy(), &vec![0.0; mg.n_dofs()], &bad_tol),
        Err(RcError::InvalidOption { .. })
    ));
    // Zero right side: zero solution, no cycle.
    let (x, rep) = vcycle::solve(mg.hierarchy(), &vec![0.0; mg.n_dofs()], &o).unwrap();
    assert!(x.iter().all(|&v| v == 0.0) && rep.converged && rep.total_cycles() == 0);
    // Out-of-range levels.
    assert!(mg.hierarchy().apply_metric(99, &[]).is_err());
    assert!(mg.hierarchy().apply(0, &[1.0]).is_err());
}

// ---------------------------------------------------------------- studies (release, ignored)

/// Iterations under refinement and timings on the Freedman–He–Wang unknot, (2, 4.5),
/// barycenter + edge lengths: 600, 1200, 2400 and 2400 refined twice and four times. Per curve,
/// thread count and solver (dense LU, multigrid with the defaults, multigrid with the 3 + 3
/// stationary cycles of `SPEC.md` §10): the build (Gram + LU, or the hierarchy), the projected
/// gradient, and one projection of a trial step of half an edge length at the fastest vertex
/// plus a random perturbation of 1 % of an edge length.
/// Env: `RC_THREADS` (list, default `1`), `RC_SIZES` (list of `base` or `base:times`),
/// `RC_DENSE_MAX` (largest |V| for the dense LU, default 4736), `RC_CONFIGS` (`dense,mg,spec`).
/// `--ignored --nocapture multigrid_study`.
#[test]
#[ignore]
fn multigrid_study() {
    if !have_data() {
        return;
    }
    let var = |k: &str, d: &str| std::env::var(k).unwrap_or_else(|_| d.to_string());
    let threads: Vec<usize> = var("RC_THREADS", "1")
        .split(',')
        .map(|t| t.parse().unwrap())
        .collect();
    let sizes: Vec<(usize, usize)> = var("RC_SIZES", "600,1200,2400,2400:1,2400:2")
        .split(',')
        .map(|t| {
            let mut it = t.split(':');
            let b = it.next().unwrap().parse().unwrap();
            let k = it.next().map_or(0, |k| k.parse().unwrap());
            (b, k)
        })
        .collect();
    let dense_max: usize = var("RC_DENSE_MAX", "4736").parse().unwrap();
    let configs = var("RC_CONFIGS", "dense,mg,spec");
    let ex = Exponents::new(2.0, 4.5).unwrap();
    let spec = MgOptions {
        smoothing_steps: 3,
        initial_smoothing_steps: 6,
        flexible_cg: false,
        ..MgOptions::default()
    };
    eprintln!(
        "| curve | |V| | threads | solver | levels | build s | gradient s | gradient cycles | projection its | projection s | projection cycles | build+gradient+projection s | g error vs dense |"
    );
    eprintln!("|---|---|---|---|---|---|---|---|---|---|---|---|---|");
    for &(base, times) in &sizes {
        let (net, pos) = freedman_refined(base, times);
        let pr = problem(net, pos, ex, &default_constraints());
        let name = if times > 0 {
            format!("freedman {base} x{}", 1 << times)
        } else {
            format!("freedman {base}")
        };
        let layout = pr.set.layout();
        let l0 = pr.net.mean_edge_length(&pr.pos);
        let mut reference: Option<Vec<f64>> = None;
        let mut trial: Option<Vec<P>> = None;
        for &t in &threads {
            for config in configs.split(',') {
                let n = pr.net.n_vertices();
                if config == "dense" && n > dense_max {
                    continue;
                }
                let p = pool(t);
                let t0 = Instant::now();
                let mut solver: Box<dyn SaddleSolver> = match config {
                    "dense" => Box::new(
                        dense_sobolev_solver(
                            p.clone(),
                            &pr.net,
                            &pr.pos,
                            pr.ex,
                            GramQuadrature::Midpoint,
                            &pr.set,
                        )
                        .unwrap(),
                    ),
                    "spec" => Box::new(mg_solver(&pr, t, spec)),
                    _ => Box::new(mg_solver(&pr, t, MgOptions::default())),
                };
                let build = t0.elapsed().as_secs_f64();
                let t1 = Instant::now();
                let g = projected_gradient(solver.as_mut(), layout, &pr.diff).unwrap();
                let grad = t1.elapsed().as_secs_f64();
                let x = layout.gather(&g.direction);
                let err = match &reference {
                    Some(r) => rel_err(&x, r),
                    None => {
                        reference = Some(x.clone());
                        0.0
                    }
                };
                let tr = trial.get_or_insert_with(|| {
                    let gmax = g
                        .direction
                        .iter()
                        .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
                        .fold(0.0, f64::max);
                    let tau = 0.5 * l0 / gmax;
                    let noise = rand_vec(3 * pr.pos.len(), 77);
                    pr.pos
                        .iter()
                        .zip(&g.direction)
                        .enumerate()
                        .map(|(v, (x, d))| {
                            let mut y = [0.0; 3];
                            for k in 0..3 {
                                y[k] = x[k] - tau * d[k] + 0.01 * l0 * noise[3 * v + k];
                            }
                            y
                        })
                        .collect()
                });
                let t2 = Instant::now();
                let pj = project_constraints(
                    solver.as_mut(),
                    &pr.set,
                    &pr.net,
                    tr,
                    &ProjectionOptions::default(),
                    l0,
                )
                .unwrap();
                let proj = t2.elapsed().as_secs_f64();
                eprintln!(
                    "| {name} | {n} | {t} | {config} | {} | {build:.3} | {grad:.3} | {} | {} {:?} | {proj:.3} | {:?} | {:.3} | {err:.1e} |",
                    g.linear.levels,
                    g.linear.iterations,
                    pj.iterations,
                    pj.status,
                    pj.linear.iter().map(|l| l.iterations).collect::<Vec<_>>(),
                    build + grad + proj,
                );
            }
        }
    }
}

/// Probe (ignored): Galerkin consistency `xᵀ Jᵀ Ā_f J x / xᵀ Ā_c x` for coarse vectors of
/// increasing frequency, and cycles against ν and the coarsest size.
#[test]
#[ignore]
fn probe_galerkin_and_smoothing() {
    let ex = Exponents::new(2.0, 4.5).unwrap();
    let (net, pos) = if have_data() {
        freedman(1200)
    } else {
        trefoil(1200)
    };
    let pr = problem(net, pos, ex, &default_constraints());
    let mg = mg_solver(&pr, 32, mg_options(40, 0.0));
    let h = mg.hierarchy();
    for l in 0..h.levels().len() - 1 {
        let cg = h.levels()[l].coarsening().unwrap();
        let nc = cg.n_coarse();
        let mut line = format!("level {l} -> {}: ", l + 1);
        for f in [1.0, 4.0, 16.0, nc as f64 / 8.0, nc as f64 / 3.0] {
            // x = (sin(2π f t), cos(2π f t), 0) by coarse vertex index (loop order).
            let x: Vec<f64> = (0..nc)
                .flat_map(|i| {
                    let t = std::f64::consts::TAU * f * i as f64 / nc as f64;
                    [t.sin(), t.cos(), 0.0]
                })
                .collect();
            let jx = cg.prolong(&x, 3).unwrap();
            let af = h.apply_metric(l, &jx).unwrap();
            let ac = h.apply_metric(l + 1, &x).unwrap();
            let q1: f64 = jx.iter().zip(&af).map(|(a, b)| a * b).sum();
            let q2: f64 = x.iter().zip(&ac).map(|(a, b)| a * b).sum();
            line += &format!("f {f:.0}: {:.3}  ", q1 / q2);
        }
        eprintln!("{line}");
    }
    for (base, times) in [(600, 0), (1200, 0), (2400, 0), (2400, 1), (2400, 2)] {
        if !have_data() {
            break;
        }
        let (net, pos) = freedman_refined(base, times);
        let pr = problem(net, pos, ex, &default_constraints());
        for (nu, nu0, fcg) in [
            (3, 6, false),
            (6, 12, true),
            (8, 16, true),
            (12, 12, true),
            (12, 24, true),
            (12, 24, false),
            (16, 24, true),
        ] {
            let o = MgOptions {
                smoothing_steps: nu,
                initial_smoothing_steps: nu0,
                flexible_cg: fcg,
                ..MgOptions::default()
            };
            let threads: usize = std::env::var("RC_THREADS")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(1);
            let mut mg = mg_solver(&pr, threads, o);
            let t = Instant::now();
            let mut cyc = Vec::new();
            for _ in 0..3 {
                let g = projected_gradient(&mut mg, pr.set.layout(), &pr.diff).unwrap();
                cyc.push(g.linear.iterations);
            }
            eprintln!(
                "|V| {} nu {nu} nu0 {nu0} fcg {fcg}: levels {} cycles {:?} solve {:.3} s residuals {}",
                pr.net.n_vertices(),
                level_sizes(&mg).len(),
                cyc,
                t.elapsed().as_secs_f64() / 3.0,
                es(&mg.last_reports()[0].residuals)
            );
        }
    }
}

/// Probe (ignored): cost of the pieces on Freedman 2400 ×4 at 1 and 16 threads.
#[test]
#[ignore]
fn probe_costs() {
    if !have_data() {
        return;
    }
    let ex = Exponents::new(2.0, 4.5).unwrap();
    let (net, pos) = freedman_refined(2400, 2);
    let pr = problem(net, pos, ex, &default_constraints());
    for t in [1, 16] {
        let t0 = Instant::now();
        let mg = mg_solver(&pr, t, mg_options(512, 0.25));
        let build = t0.elapsed().as_secs_f64();
        let h = mg.hierarchy();
        let x = rand_vec(h.len(0), 3);
        let time = |f: &dyn Fn()| {
            let t = Instant::now();
            for _ in 0..20 {
                f();
            }
            t.elapsed().as_secs_f64() / 20.0 * 1e3
        };
        let ta = time(&|| {
            h.apply_metric(0, &x).unwrap();
        });
        let tp = time(&|| {
            h.levels()[0].projector().project(&x).unwrap();
        });
        let th = time(&|| {
            h.apply(0, &x).unwrap();
        });
        let xc = rand_vec(h.len(h.n_levels() - 1), 4);
        let tc = time(&|| {
            h.coarse_solve(&xc).unwrap();
        });
        eprintln!(
            "threads {t}: build {build:.3} s, levels {:?}, A {ta:.2} ms, P {tp:.2} ms, H {th:.2} ms, coarse solve {tc:.2} ms, near entries {:?}",
            level_sizes(&mg),
            h.levels()[0].hierarchical_metric().map(|m| m.stats())
        );
    }
}

// ---------------------------------------------------------------- the solver in a session

fn check_constraints(s: &Session) -> f64 {
    let phi = s.constraints().values(s.network(), s.positions()).unwrap();
    let tol = s
        .constraints()
        .row_tolerances(
            s.network(),
            s.positions(),
            s.options().projection.tolerance,
            s.mean_edge_length0(),
        )
        .unwrap();
    let mut worst: f64 = 0.0;
    for (k, (p, t)) in phi.iter().zip(&tol).enumerate() {
        assert!(p.abs() <= *t, "row {k}: violation {p} > {t}");
        worst = worst.max(p.abs() / t);
    }
    worst
}

fn session_options(solver: Solver, energy: EnergyEval, threads: usize) -> FlowOptions {
    FlowOptions {
        backend: Backend::Explicit { energy, solver },
        threads: Some(threads),
        ..FlowOptions::default()
    }
}

/// The multigrid factory in a session (`Solver::Custom`): ten steps on a trefoil of 300
/// vertices (four levels) follow the dense flow closely, every step reports the multigrid solve,
/// the constraints hold.
#[test]
fn session_with_the_multigrid_factory_follows_the_dense_flow() {
    let (net, pos) = trefoil(300);
    let o = MgOptions {
        coarsest_vertices: 40,
        ..MgOptions::default()
    };
    let factory = Arc::new(MultigridSolvers::new(o).unwrap());
    assert_eq!(
        ddg_repulsive_curves::SolverFactory::name(&*factory),
        "multigrid"
    );
    let mut mg = Session::new(
        net.clone(),
        pos.clone(),
        session_options(Solver::Custom(factory), EnergyEval::Exact, 8),
    )
    .unwrap();
    let mut dense = Session::new(
        net,
        pos,
        session_options(Solver::Dense, EnergyEval::Exact, 8),
    )
    .unwrap();
    let e0 = mg.energy();
    for _ in 0..10 {
        let rm = mg.step().unwrap();
        let rd = dense.step().unwrap();
        assert_eq!(rm.linear.solver, "multigrid");
        assert!(rm.linear.levels >= 3, "{:?}", rm.linear);
        assert!(rm.linear.relative_residual <= 1e-2, "{:?}", rm.linear);
        eprintln!(
            "step {}: E mg {:.6} dense {:.6}, |g| mg {:.4e} dense {:.4e}, cycles {}, levels {}",
            rm.iteration,
            rm.energy_after,
            rd.energy_after,
            rm.gradient_norm,
            rd.gradient_norm,
            rm.linear.iterations,
            rm.linear.levels
        );
        // Same base curve at the first step: same direction up to the solver tolerances.
        if rm.iteration == 1 {
            assert!(
                (rm.gradient_norm - rd.gradient_norm).abs() <= 1e-2 * rd.gradient_norm,
                "{} vs {}",
                rm.gradient_norm,
                rd.gradient_norm
            );
        }
        assert!(rm.energy_after <= rm.energy_before);
    }
    assert!(mg.energy() < e0);
    assert!((mg.energy() - dense.energy()).abs() <= 1e-3 * dense.energy());
    check_constraints(&mg);
}

/// Freedman–He–Wang 600, (2, 4.5), Barnes–Hut energy + multigrid (or `RC_SOLVER=dense`,
/// `RC_ENERGY=exact`, `RC_START=paper|grow`, `RC_THREADS`), to the near-circle state and the
/// stop rule. `--ignored --nocapture freedman_600_multigrid_flow`.
#[test]
#[ignore]
fn freedman_600_multigrid_flow() {
    if !have_data() {
        return;
    }
    let var = |k: &str, d: &str| std::env::var(k).unwrap_or_else(|_| d.to_string());
    let threads: usize = var("RC_THREADS", "16").parse().unwrap();
    let steps: usize = var("RC_STEPS", "1500").parse().unwrap();
    let size: usize = var("RC_SIZE", "600").parse().unwrap();
    let (net, pos) = freedman(size);
    let energy = match var("RC_ENERGY", "bh").as_str() {
        "exact" => EnergyEval::Exact,
        _ => EnergyEval::BarnesHut(BhOptions::default()),
    };
    let solver = match var("RC_SOLVER", "mg").as_str() {
        "dense" => Solver::Dense,
        _ => Solver::Custom(Arc::new(
            MultigridSolvers::new(MgOptions::default()).unwrap(),
        )),
    };
    let start = match var("RC_START", "grow").as_str() {
        "paper" => StartRule::TwoThirdsTauMax,
        _ => StartRule::grow_previous(),
    };
    let opts = FlowOptions {
        exponents: Exponents::new(2.0, 4.5).unwrap(),
        line_search: LineSearch::guarded_with_start(start),
        max_iterations: steps,
        ..session_options(solver, energy, threads)
    };
    eprintln!("options: {:?}", opts.backend);
    let t0 = Instant::now();
    let mut s = Session::new(net, pos, opts).unwrap();
    eprintln!(
        "start: |V| {} E {:.4}",
        s.network().n_vertices(),
        s.energy()
    );
    let mut round_at = None;
    let mut split = [0.0f64; 7];
    let mut cycles = 0usize;
    let mut solves = 0usize;
    loop {
        let r = s.step().unwrap();
        split[0] += r.seconds.energy;
        split[1] += r.seconds.metric;
        split[2] += r.seconds.solve;
        split[3] += r.seconds.line_search_energy;
        split[4] += r.seconds.projection;
        split[5] += r.seconds.ccd;
        split[6] += r.seconds.total;
        cycles += r.linear.iterations;
        solves += 1;
        let (sp, pl, tu) = roundness(s.positions());
        if r.iteration.is_multiple_of(25) || s.status().is_terminal() {
            eprintln!(
                "step {:4} {:?} E {:.5} |g| {:.2e} tau {:.2e} halv {} proj {} cycles {} spread {sp:.4} planarity {pl:.4} turning/2pi {tu:.4} t {:.1}s",
                r.iteration,
                r.status,
                r.energy_after,
                r.gradient_norm,
                r.tau,
                r.halvings,
                r.projection_iterations,
                r.linear.iterations,
                t0.elapsed().as_secs_f64()
            );
        }
        if round_at.is_none() && sp < 0.02 && pl < 0.02 && tu < 1.01 {
            round_at = Some((r.iteration, t0.elapsed().as_secs_f64(), s.energy()));
            eprintln!(
                "near-circle at step {} ({:.1} s)",
                r.iteration,
                t0.elapsed().as_secs_f64()
            );
        }
        if s.status().is_terminal() {
            break;
        }
    }
    let n = s.iteration() as f64;
    let exact_final = ddg_repulsive_curves::energy(
        s.network(),
        s.positions(),
        Exponents::new(2.0, 4.5).unwrap(),
        &EnergyEval::Exact,
        Some(threads),
    )
    .unwrap()
    .energy;
    eprintln!(
        "end: status {:?} steps {} energy {:.5} (exact {exact_final:.5}) time {:.1}s, {:.3} s/step; per step: energy {:.3} metric {:.3} solve {:.3} ls-energy {:.3} projection {:.3} ccd {:.3}; gradient cycles/step {:.2}",
        s.status(),
        s.iteration(),
        s.energy(),
        t0.elapsed().as_secs_f64(),
        split[6] / n,
        split[0] / n,
        split[1] / n,
        split[2] / n,
        split[3] / n,
        split[4] / n,
        split[5] / n,
        cycles as f64 / solves as f64,
    );
    let worst = check_constraints(&s);
    eprintln!("largest violation / tolerance {worst:.3}");
    let (it, secs, e) = round_at.expect("near-circle within the budget");
    eprintln!("RESULT near-circle iterations {it} seconds {secs:.1} energy {e:.5}");
}

/// Roundness of a closed loop: `(max r − min r)/mean r` about the length-weighted barycenter,
/// planarity, total turning / 2π (as in `tests/flow.rs`).
fn roundness(pos: &[P]) -> (f64, f64, f64) {
    let n = pos.len() as f64;
    let n3 = |a: P| (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    let mut c = [0.0; 3];
    let mut total = 0.0;
    for i in 0..pos.len() {
        let (a, b) = (pos[i], pos[(i + 1) % pos.len()]);
        let l = n3([b[0] - a[0], b[1] - a[1], b[2] - a[2]]);
        total += l;
        for k in 0..3 {
            c[k] += l * 0.5 * (a[k] + b[k]);
        }
    }
    for k in 0..3 {
        c[k] /= total;
    }
    let r: Vec<f64> = pos
        .iter()
        .map(|p| n3([p[0] - c[0], p[1] - c[1], p[2] - c[2]]))
        .collect();
    let mean = r.iter().sum::<f64>() / n;
    let spread =
        (r.iter().cloned().fold(0.0, f64::max) - r.iter().cloned().fold(f64::MAX, f64::min)) / mean;
    let mut cov = [[0.0; 3]; 3];
    for p in pos {
        let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += d[i] * d[j] / n;
            }
        }
    }
    let tr = cov[0][0] + cov[1][1] + cov[2][2];
    let mut v = [0.3, 0.5, 0.8];
    let mut lam = 0.0;
    for _ in 0..500 {
        let w: Vec<f64> = (0..3)
            .map(|i| tr * v[i] - (0..3).map(|j| cov[i][j] * v[j]).sum::<f64>())
            .collect();
        let nw = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        lam = nw;
        v = [w[0] / nw, w[1] / nw, w[2] / nw];
    }
    let smallest = (tr - lam).max(0.0).sqrt();
    let m = pos.len();
    let mut turning = 0.0;
    for i in 0..m {
        let a = pos[(i + m - 1) % m];
        let b = pos[i];
        let cc = pos[(i + 1) % m];
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [cc[0] - b[0], cc[1] - b[1], cc[2] - b[2]];
        let cosang = (u[0] * w[0] + u[1] * w[1] + u[2] * w[2]) / (n3(u) * n3(w));
        turning += cosang.clamp(-1.0, 1.0).acos();
    }
    (spread, smallest / mean, turning / std::f64::consts::TAU)
}

/// Three session steps with the multigrid factory give bitwise identical curves at 1, 7 and 64
/// threads (D25), with the Barnes–Hut energy.
#[test]
fn session_with_multigrid_is_bitwise_identical_across_threads() {
    let (net, pos) = trefoil(300);
    let o = MgOptions {
        coarsest_vertices: 40,
        ..MgOptions::default()
    };
    let mut curves = Vec::new();
    for t in [1, 7, 64] {
        let factory = Arc::new(MultigridSolvers::new(o).unwrap());
        let mut s = Session::new(
            net.clone(),
            pos.clone(),
            session_options(
                Solver::Custom(factory),
                EnergyEval::BarnesHut(BhOptions::default()),
                t,
            ),
        )
        .unwrap();
        for _ in 0..3 {
            s.step().unwrap();
        }
        curves.push(s.positions().to_vec());
    }
    assert!(curves[1] == curves[0] && curves[2] == curves[0]);
}
