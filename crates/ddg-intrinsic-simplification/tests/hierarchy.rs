//! Tests of the intrinsic hierarchy, the multigrid on it and the conformalized mean curvature
//! flow (SPEC §9).

mod common;

use common::{Shape, load_data};
use ddg_core::generators::{self, SplitMix64};
use ddg_core::{SparseMatrix, TriMesh};
use ddg_intrinsic_simplification::cmcf::{CmcfError, CmcfFlow, CmcfMultigridOptions};
use ddg_intrinsic_simplification::hierarchy::{
    Hierarchy, HierarchyError, HierarchyOptions, build_hierarchy,
};
use ddg_intrinsic_simplification::multigrid::MultigridOptions;
use ddg_intrinsic_simplification::{SimplifyOptions, Target, simplify};

fn shapes() -> Vec<Shape> {
    let g = |name: &str, m: generators::GeneratedMesh| Shape::new(name, m.positions, m.faces);
    vec![
        g("icosphere 3", generators::icosphere(3).unwrap()),
        g("torus 24x12", generators::torus(24, 12, 1.0, 0.4).unwrap()),
        g(
            "bumpy grid 16x16",
            generators::bumpy_grid_disk(16, 16, 0.25).unwrap(),
        ),
        g(
            "annulus 24x6",
            generators::annulus(24, 6, 0.5, 1.0).unwrap(),
        ),
    ]
}

fn small_options(fresh: bool) -> HierarchyOptions {
    HierarchyOptions {
        ratio: 0.25,
        n_min: 30,
        fresh_levels: fresh,
        ..HierarchyOptions::default()
    }
}

fn row_sums(p: &SparseMatrix) -> Vec<f64> {
    let mut s = vec![0.0; p.nrows()];
    for (r, _, v) in p.triplets() {
        s[r as usize] += v;
    }
    s
}

fn random_vec(n: usize, seed: u64) -> Vec<f64> {
    let mut r = SplitMix64::new(seed);
    (0..n).map(|_| r.uniform(-1.0, 1.0)).collect()
}

fn quad(a: &SparseMatrix, x: &[f64]) -> f64 {
    let mut y = vec![0.0; x.len()];
    a.mul_vec(x, &mut y).unwrap();
    x.iter().zip(&y).map(|(a, b)| a * b).sum()
}

fn check_prolongations(name: &str, h: &Hierarchy) {
    assert!(h.n_levels() >= 2, "{name}: only {} level", h.n_levels());
    for (k, p) in h.prolongations().iter().enumerate() {
        let (nk, nk1) = (h.levels()[k].n_vertices(), h.levels()[k + 1].n_vertices());
        assert_eq!((p.nrows(), p.ncols()), (nk, nk1), "{name}: P_{k} shape");
        assert!(nk1 < nk, "{name}: level {} is not coarser", k + 1);
        let sums = row_sums(p);
        for (i, s) in sums.iter().enumerate() {
            assert!((s - 1.0).abs() < 1e-12, "{name}: P_{k} row {i} sums to {s}");
        }
        // Constants are reproduced.
        let mut y = vec![0.0; nk];
        p.mul_vec(&vec![1.0; nk1], &mut y).unwrap();
        assert!(
            y.iter().all(|v| (v - 1.0).abs() < 1e-12),
            "{name}: P_{k} 1 != 1"
        );
        assert!(
            p.values().iter().all(|&v| v >= -1e-15),
            "{name}: P_{k} has a negative weight"
        );
        // A vertex that survives has one entry, 1, in its own column.
        let coarse: std::collections::HashMap<u32, usize> = h.levels()[k + 1]
            .vertices
            .iter()
            .enumerate()
            .map(|(c, &v)| (v, c))
            .collect();
        let mut count = vec![0usize; nk];
        for (r, _, _) in p.triplets() {
            count[r as usize] += 1;
        }
        let mut kept = 0;
        for (i, &v) in h.levels()[k].vertices.iter().enumerate() {
            if let Some(&c) = coarse.get(&v) {
                kept += 1;
                assert_eq!(
                    count[i], 1,
                    "{name}: kept vertex {v} has {} entries",
                    count[i]
                );
                assert!((p.get(i, c) - 1.0).abs() < 1e-15, "{name}: kept vertex {v}");
            }
        }
        assert_eq!(kept, nk1, "{name}: P_{k} kept vertices");
    }
}

#[test]
fn prolongations_have_unit_rows_and_reproduce_constants() {
    for s in shapes() {
        for fresh in [false, true] {
            let mesh = s.coarse(&[]);
            let h = build_hierarchy(&mesh, &small_options(fresh))
                .unwrap_or_else(|e| panic!("{} fresh={fresh}: {e}", s.name));
            let name = format!("{} fresh={fresh}", s.name);
            check_prolongations(&name, &h);
            // The counts follow the reference's rule.
            let counts = &h.report().counts;
            assert_eq!(counts[0], s.positions.len());
            for w in counts.windows(2) {
                assert_eq!(w[1], ((w[0] as f64) * 0.25).round() as usize, "{name}");
            }
            let last = *counts.last().unwrap();
            assert!(last <= 30 || last < counts[0], "{name}: {counts:?}");
        }
    }
}

#[test]
fn level_operators_match_the_coarse_mesh_and_are_psd() {
    for s in shapes() {
        let mesh = s.coarse(&[]);
        let h = build_hierarchy(&mesh, &small_options(false)).unwrap();
        for (k, level) in h.levels().iter().enumerate().skip(1) {
            // The same run stopped at this count has the same coarse mesh.
            let mut m = mesh.clone();
            let n = level.n_vertices();
            let opts = SimplifyOptions {
                target: Target::Vertices(n),
                ..SimplifyOptions::default()
            };
            simplify(&mut m, &opts).unwrap();
            assert_eq!(m.vertices(), &level.vertices[..], "{} level {k}", s.name);
            let lm = m.cotan_laplacian();
            let diff = lm
                .triplets()
                .map(|(i, j, v)| (v - level.laplacian.get(i as usize, j as usize)).abs())
                .fold(0.0f64, f64::max);
            let scale = lm.values().iter().fold(0.0f64, |a, v| a.max(v.abs()));
            assert!(
                diff <= 1e-9 * scale,
                "{} level {k}: Laplacian differs by {diff}",
                s.name
            );
            let mass = m.lumped_mass();
            for (a, b) in mass.iter().zip(&level.mass) {
                assert!(
                    (a - b).abs() <= 1e-12 * a.abs().max(1e-300),
                    "{} mass",
                    s.name
                );
            }
        }
        for (k, level) in h.levels().iter().enumerate() {
            let l = &level.laplacian;
            assert!(l.symmetry_defect().unwrap() < 1e-9, "{} level {k}", s.name);
            let ones = vec![1.0; level.n_vertices()];
            let mut y = vec![0.0; ones.len()];
            l.mul_vec(&ones, &mut y).unwrap();
            assert!(
                y.iter().all(|v| v.abs() < 1e-9),
                "{} level {k}: L 1 != 0",
                s.name
            );
            for seed in 0..5 {
                let x = random_vec(level.n_vertices(), seed);
                assert!(quad(l, &x) >= -1e-9, "{} level {k}: L is not PSD", s.name);
            }
            assert!(level.mass.iter().all(|&m| m > 0.0));
        }
    }
}

#[test]
fn galerkin_levels_are_spd() {
    for s in shapes() {
        let mesh = s.coarse(&[]);
        let h = build_hierarchy(&mesh, &small_options(false)).unwrap();
        let l0 = &h.levels()[0];
        let mut a = l0.laplacian.clone();
        a.add_to_diagonal(&l0.mass).unwrap();
        let mg = h.multigrid(&a, &[], MultigridOptions::default()).unwrap();
        assert_eq!(mg.level_matrices().len(), h.n_levels());
        for (k, m) in mg.level_matrices().iter().enumerate() {
            assert!(
                m.symmetry_defect().unwrap()
                    < 1e-9 * m.values().iter().fold(1e-300f64, |a, v| a.max(v.abs()))
            );
            assert!(m.diagonal().iter().all(|&d| d > 0.0));
            for seed in 0..5 {
                let x = random_vec(m.nrows(), 100 + seed);
                assert!(
                    quad(m, &x) > 0.0,
                    "{} Galerkin level {k} is not SPD",
                    s.name
                );
            }
        }
    }
}

fn direct_solve(a: &SparseMatrix, b: &[f64], known: &[u32]) -> Vec<f64> {
    use ddg_core::sparse::{Method, Parallelism, Symbolic};
    let n = a.nrows();
    let mut is_k = vec![false; n];
    for &k in known {
        is_k[k as usize] = true;
    }
    let unk: Vec<u32> = (0..n as u32).filter(|&i| !is_k[i as usize]).collect();
    let auu = a.submatrix(&unk, &unk).unwrap();
    let mut rhs: Vec<f64> = unk.iter().map(|&i| b[i as usize]).collect();
    let f = Symbolic::analyze(&auu, Method::Llt)
        .unwrap()
        .factor(&auu, Parallelism::Sequential)
        .unwrap();
    f.solve_in_place(&mut rhs, Parallelism::Sequential).unwrap();
    let mut x = vec![0.0; n];
    for (&i, &v) in unk.iter().zip(&rhs) {
        x[i as usize] = v;
    }
    x
}

fn rel_diff(x: &[f64], y: &[f64]) -> f64 {
    let d: f64 = x
        .iter()
        .zip(y)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    d / y.iter().map(|v| v * v).sum::<f64>().sqrt()
}

#[test]
fn multigrid_on_the_intrinsic_hierarchy_equals_the_direct_solve_on_spot() {
    let Some(spot) = load_data("intrinsic-simplification/spot.obj") else {
        return;
    };
    let mesh = spot.coarse(&[]);
    for fresh in [false, true] {
        let h = build_hierarchy(
            &mesh,
            &HierarchyOptions {
                fresh_levels: fresh,
                ..HierarchyOptions::default()
            },
        )
        .unwrap();
        assert_eq!(h.report().counts, vec![2930, 733, 183], "fresh={fresh}");
        let l0 = &h.levels()[0];
        let n = l0.n_vertices();
        let total: f64 = l0.mass.iter().sum();
        // L + eps M.
        let mut a = l0.laplacian.clone();
        a.add_to_diagonal(&l0.mass.iter().map(|m| m / total).collect::<Vec<_>>())
            .unwrap();
        let b = random_vec(n, 1);
        let xd = direct_solve(&a, &b, &[]);
        let opts = MultigridOptions {
            tolerance: 1e-11,
            max_cycles: 200,
            ..MultigridOptions::default()
        };
        let mg = h.multigrid(&a, &[], opts.clone()).unwrap();
        let (x, rep) = mg.solve(&b, &[], None).unwrap();
        assert!(rep.converged, "fresh={fresh}: {rep:?}");
        assert!(rel_diff(&x, &xd) < 1e-9, "V-cycles: {}", rel_diff(&x, &xd));
        assert!(
            rep.convergence_factor < 0.5,
            "factor {}",
            rep.convergence_factor
        );
        let (x, rep) = mg.pcg(&b, &[], None, 100).unwrap();
        assert!(rep.converged);
        assert!(rel_diff(&x, &xd) < 1e-9, "PCG: {}", rel_diff(&x, &xd));
        // The reference's Dirichlet example: L z = M 1, z = 0 at 0, 1, 2.
        let known = [0u32, 1, 2];
        let mut bd = l0.mass.clone();
        for &k in &known {
            bd[k as usize] = 0.0;
        }
        let xd = direct_solve(&l0.laplacian, &bd, &known);
        let mg = h.multigrid(&l0.laplacian, &known, opts).unwrap();
        let z0 = random_vec(n, 5);
        let (x, rep) = mg.solve(&bd, &[0.0; 3], Some(&z0)).unwrap();
        assert!(rep.converged, "Dirichlet fresh={fresh}: {rep:?}");
        assert!(rel_diff(&x, &xd) < 1e-9, "Dirichlet: {}", rel_diff(&x, &xd));
        assert!(rep.convergence_factor < 0.5);
    }
}

#[test]
fn bad_hierarchy_options_are_errors() {
    let mesh = shapes()[0].coarse(&[]);
    for ratio in [0.0, 1.0, -0.5, f64::NAN, 1.5] {
        let e = build_hierarchy(
            &mesh,
            &HierarchyOptions {
                ratio,
                ..HierarchyOptions::default()
            },
        )
        .unwrap_err();
        assert!(matches!(e, HierarchyError::BadOptions { .. }), "{e}");
        assert!(e.to_string().contains("ratio"));
    }
    // A mesh already below n_min: one level, no prolongation.
    let h = build_hierarchy(
        &mesh,
        &HierarchyOptions {
            n_min: 100_000,
            ..HierarchyOptions::default()
        },
    )
    .unwrap();
    assert_eq!(h.n_levels(), 1);
    assert!(h.prolongations().is_empty());
}

fn flow_pair(s: &Shape, h: f64) -> (TriMesh, Hierarchy, CmcfFlow) {
    let mesh = TriMesh::new(s.positions.len(), s.faces.clone()).unwrap();
    let hier = build_hierarchy(
        &s.coarse(&[]),
        &HierarchyOptions {
            n_min: 100,
            ..HierarchyOptions::default()
        },
    )
    .unwrap();
    let flow = CmcfFlow::new(&mesh, &s.positions, h).unwrap();
    (mesh, hier, flow)
}

fn pos_diff(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let (mut d, mut s) = (0.0, 0.0);
    for (p, q) in a.iter().zip(b) {
        for c in 0..3 {
            d += (p[c] - q[c]).powi(2);
            s += q[c].powi(2);
        }
    }
    (d / s).sqrt()
}

#[test]
fn cmcf_multigrid_steps_agree_with_direct_steps() {
    let mut shapes = vec![Shape::new(
        "icosphere 4",
        generators::icosphere(4).unwrap().positions,
        generators::icosphere(4).unwrap().faces,
    )];
    if let Some(spot) = load_data("intrinsic-simplification/spot.obj") {
        shapes.push(spot);
    }
    for s in shapes {
        // Unit-area scale so that h is an area fraction.
        let g0 = {
            let mesh = TriMesh::new(s.positions.len(), s.faces.clone()).unwrap();
            ddg_core::Geometry::from_positions(&mesh, &s.positions)
                .unwrap()
                .total_area()
        };
        let scale = 1.0 / g0.sqrt();
        let positions: Vec<[f64; 3]> = s.positions.iter().map(|p| p.map(|c| c * scale)).collect();
        let s = Shape::new(&s.name, positions, s.faces.clone());
        let (_mesh, hier, mut direct) = flow_pair(&s, 1e-3);
        let mut mg = direct.clone();
        let opts = CmcfMultigridOptions {
            multigrid: MultigridOptions {
                tolerance: 1e-12,
                ..MultigridOptions::default()
            },
            ..CmcfMultigridOptions::default()
        };
        for step in 1..=5 {
            // One step from the same state.
            let mut probe = direct.clone();
            let rd = direct.step_direct().unwrap();
            let rm = probe.step_multigrid(&hier, &opts).unwrap();
            assert!(rm.converged && rm.iterations.len() == 3, "{rm:?}");
            let d = pos_diff(probe.positions(), direct.positions());
            assert!(d < 1e-8, "{} step {step}: one-step difference {d}", s.name);
            assert!(rd.relative_residual < 1e-12 && rm.relative_residual < 1e-10);
            // The multigrid trajectory.
            mg.step_multigrid(&hier, &opts).unwrap();
            let d = pos_diff(mg.positions(), direct.positions());
            assert!(
                d < 1e-8,
                "{} step {step}: trajectory difference {d}",
                s.name
            );
        }
        assert_eq!(mg.steps_done(), 5);
    }
}

#[test]
fn cmcf_shrinks_a_sphere_like_the_continuum_flow() {
    // Continuum: r(t)^2 = r0^2 - 4t for the unit sphere; backward Euler with lumped masses is
    // first order in h.
    let g = generators::icosphere(4).unwrap();
    let mesh = TriMesh::new(g.positions.len(), g.faces.clone()).unwrap();
    let mut flow = CmcfFlow::new(&mesh, &g.positions, 0.002).unwrap();
    for _ in 0..10 {
        flow.step_direct().unwrap();
    }
    let r2: f64 = flow
        .positions()
        .iter()
        .map(|p| p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
        .sum::<f64>()
        / g.positions.len() as f64;
    let expected = 1.0 - 4.0 * 0.02;
    assert!(
        (r2 - expected).abs() < 0.03 * expected,
        "r^2 = {r2}, continuum {expected}"
    );
    // Rescaling restores the input's area and centre.
    flow.normalize().unwrap();
    let g1 = ddg_core::Geometry::from_positions(flow.mesh(), flow.positions()).unwrap();
    let g0 = ddg_core::Geometry::from_positions(&mesh, &g.positions).unwrap();
    assert!((g1.total_area() - g0.total_area()).abs() < 1e-12 * g0.total_area());
}

#[test]
fn cmcf_bad_input_is_an_error() {
    let g = generators::icosphere(1).unwrap();
    let mesh = TriMesh::new(g.positions.len(), g.faces.clone()).unwrap();
    for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let e = CmcfFlow::new(&mesh, &g.positions, h).err().unwrap();
        assert!(matches!(e, CmcfError::BadInput { .. }), "{e}");
    }
    let mut bad = g.positions.clone();
    bad[0] = [f64::NAN; 3];
    assert!(CmcfFlow::new(&mesh, &bad, 0.01).is_err());
    assert!(CmcfFlow::new(&mesh, &g.positions[..5], 0.01).is_err());
    // A hierarchy of another mesh.
    let other = Shape::new(
        "icosphere 3",
        generators::icosphere(3).unwrap().positions,
        generators::icosphere(3).unwrap().faces,
    );
    let hier = build_hierarchy(&other.coarse(&[]), &small_options(false)).unwrap();
    let mut flow = CmcfFlow::new(&mesh, &g.positions, 0.01).unwrap();
    let e = flow
        .step_multigrid(&hier, &CmcfMultigridOptions::default())
        .unwrap_err();
    assert!(matches!(e, CmcfError::HierarchyMismatch { .. }), "{e}");
}
