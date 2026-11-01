//! Tests of the Galerkin multigrid solver (SPEC §9.2) on hierarchies of regular subdivision,
//! where the prolongation is exact linear interpolation.

use std::path::PathBuf;

use ddg_core::generators::{self, SplitMix64};
use ddg_core::io::read_obj;
use ddg_core::operators::{cotan_laplacian, lumped_mass_normalized};
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::{Geometry, SparseMatrix, TriMesh};
use ddg_intrinsic_simplification::multigrid::{
    Multigrid, MultigridError, MultigridOptions, midpoint_subdivision_prolongation,
};

/// One 4-to-1 split of a mesh with midpoint positions; numbering of
/// `ddg_intrinsic_apps::fem::regular_subdivision` (new vertex `n + e` on edge `e`).
fn subdivide(mesh: &TriMesh, pos: &[[f64; 3]]) -> (TriMesh, Vec<[f64; 3]>) {
    let n = mesh.n_vertices() as u32;
    let mut p = pos.to_vec();
    for &[a, b] in mesh.edges() {
        let (pa, pb) = (pos[a as usize], pos[b as usize]);
        p.push([0, 1, 2].map(|k| 0.5 * (pa[k] + pb[k])));
    }
    let mut faces = Vec::new();
    for (f, &[a, b, c]) in mesh.faces().iter().enumerate() {
        let [mab, mbc, mca] = [0, 1, 2].map(|k| n + mesh.edge(3 * f + k) as u32);
        faces.push([a, mab, mca]);
        faces.push([mab, b, mbc]);
        faces.push([mca, mbc, c]);
        faces.push([mab, mbc, mca]);
    }
    (TriMesh::new(p.len(), faces).unwrap(), p)
}

struct Hierarchy {
    /// Finest mesh and positions.
    mesh: TriMesh,
    pos: Vec<[f64; 3]>,
    /// `P_0` (fine <- next) ... finest first.
    p: Vec<SparseMatrix>,
}

fn hierarchy(mesh: TriMesh, pos: Vec<[f64; 3]>, splits: usize) -> Hierarchy {
    let mut levels = vec![(mesh, pos)];
    for _ in 0..splits {
        let (m, p) = {
            let (m, p) = levels.last().unwrap();
            subdivide(m, p)
        };
        levels.push((m, p));
    }
    // P between consecutive levels, finest first.
    let mut p = Vec::new();
    for k in (0..splits).rev() {
        p.push(midpoint_subdivision_prolongation(&levels[k].0).unwrap());
    }
    let (mesh, pos) = levels.pop().unwrap();
    Hierarchy { mesh, pos, p }
}

fn system(h: &Hierarchy, mass_weight: f64) -> (SparseMatrix, Vec<f64>) {
    let g = Geometry::from_positions(&h.mesh, &h.pos).unwrap();
    let mut a = cotan_laplacian(&h.mesh, &g).unwrap();
    let m = lumped_mass_normalized(&h.mesh, &g).unwrap();
    a.add_to_diagonal(&m.iter().map(|x| x * mass_weight).collect::<Vec<_>>())
        .unwrap();
    (a, m)
}

fn direct(a: &SparseMatrix, b: &[f64], known: &[u32], vals: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let mut is_k = vec![false; n];
    for &k in known {
        is_k[k as usize] = true;
    }
    let unk: Vec<u32> = (0..n as u32).filter(|&i| !is_k[i as usize]).collect();
    let auu = a.submatrix(&unk, &unk).unwrap();
    let auk = a.submatrix(&unk, known).unwrap();
    let mut rhs: Vec<f64> = unk.iter().map(|&i| b[i as usize]).collect();
    let mut t = vec![0.0; unk.len()];
    auk.mul_vec(vals, &mut t).unwrap();
    for (r, ti) in rhs.iter_mut().zip(&t) {
        *r -= ti;
    }
    let f = Symbolic::analyze(&auu, Method::Llt)
        .unwrap()
        .factor(&auu, Parallelism::Sequential)
        .unwrap();
    f.solve_in_place(&mut rhs, Parallelism::Sequential).unwrap();
    let mut x = vec![0.0; n];
    for (&i, &v) in unk.iter().zip(&rhs) {
        x[i as usize] = v;
    }
    for (&i, &v) in known.iter().zip(vals) {
        x[i as usize] = v;
    }
    x
}

fn rel_diff(x: &[f64], y: &[f64]) -> f64 {
    let d: f64 = x.iter().zip(y).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();
    d / y.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn random_rhs(n: usize, seed: u64) -> Vec<f64> {
    let mut r = SplitMix64::new(seed);
    (0..n).map(|_| r.uniform(-1.0, 1.0)).collect()
}

fn flat_hierarchy() -> Hierarchy {
    let g = generators::grid_disk(4, 4).unwrap();
    hierarchy(g.topology().unwrap(), g.positions, 4)
}

#[test]
fn flat_grid_matches_direct_solve_and_converges_fast() {
    let h = flat_hierarchy();
    assert_eq!(h.mesh.n_vertices(), 65 * 65);
    let (a, _) = system(&h, 1.0);
    let b = random_rhs(a.nrows(), 1);
    let mg = Multigrid::new(&a, &h.p, &[], MultigridOptions::default()).unwrap();
    let (x, rep) = mg.solve(&b, &[], None).unwrap();
    let xd = direct(&a, &b, &[], &[]);
    assert!(rep.converged);
    assert!(rel_diff(&x, &xd) < 1e-9, "rel diff {}", rel_diff(&x, &xd));
    assert!(rep.convergence_factor < 0.2, "factor {}", rep.convergence_factor);
    println!(
        "flat grid: levels {:?}, {} cycles, factor {:.3}, setup {:.3}s solve {:.3}s",
        rep.level_sizes, rep.cycles, rep.convergence_factor, mg.info().setup_seconds, rep.solve_seconds
    );
    // Residual decreases at every cycle.
    assert!(rep.residuals.windows(2).all(|w| w[1] < w[0]));
}

fn spot_hierarchy(splits: usize) -> Option<Hierarchy> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/models/spot/spot/spot_triangulated.obj");
    if !path.is_file() {
        eprintln!("SKIPPED: {} not found; nothing was checked for spot", path.display());
        return None;
    }
    let obj = read_obj(&path).unwrap();
    let faces = obj.triangles().expect("spot_triangulated.obj is triangles");
    let mesh = TriMesh::new(obj.positions.len(), faces).unwrap();
    Some(hierarchy(mesh, obj.positions, splits))
}

#[test]
fn spot_matches_direct_solve_and_converges_fast() {
    let Some(h) = spot_hierarchy(2) else { return };
    let (a, _) = system(&h, 1.0);
    let b = random_rhs(a.nrows(), 2);
    let mg = Multigrid::new(&a, &h.p, &[], MultigridOptions::default()).unwrap();
    let (x, rep) = mg.solve(&b, &[], None).unwrap();
    let xd = direct(&a, &b, &[], &[]);
    assert!(rep.converged);
    assert!(rel_diff(&x, &xd) < 1e-9, "rel diff {}", rel_diff(&x, &xd));
    // Spot's mesh is irregular (edge lengths vary 27-fold, anisotropic cells), and subdividing
    // it keeps that: with 2 + 2 point Gauss-Seidel sweeps the factor is about 0.6 (measured),
    // with 8 + 8 about 0.1; 2 + 2 sweeps give < 0.2 on regular meshes (the grid and sphere
    // tests). Conjugate gradients with the V-cycle needs about 16 iterations.
    assert!(rep.convergence_factor < 0.7, "factor {}", rep.convergence_factor);
    let strong = MultigridOptions {
        pre_sweeps: 8,
        post_sweeps: 8,
        ..Default::default()
    };
    let mg8 = Multigrid::new(&a, &h.p, &[], strong).unwrap();
    let (_, rep8) = mg8.solve(&b, &[], None).unwrap();
    assert!(rep8.convergence_factor < 0.2, "8/8 factor {}", rep8.convergence_factor);
    let (xp, repp) = mg.pcg(&b, &[], None, 100).unwrap();
    assert!(repp.converged && repp.cycles <= 25, "pcg {} its", repp.cycles);
    assert!(rel_diff(&xp, &xd) < 1e-8);
    println!(
        "spot: levels {:?}, {} cycles, factor {:.3}, setup {:.3}s solve {:.3}s",
        rep.level_sizes, rep.cycles, rep.convergence_factor, mg.info().setup_seconds, rep.solve_seconds
    );
}

#[test]
fn icosphere_matches_direct_solve_and_converges_fast() {
    let g = generators::icosphere(2).unwrap();
    let h = hierarchy(g.topology().unwrap(), g.positions.clone(), 3);
    let (a, _) = system(&h, 1.0);
    let b = random_rhs(a.nrows(), 7);
    let mg = Multigrid::new(&a, &h.p, &[], MultigridOptions::default()).unwrap();
    let (x, rep) = mg.solve(&b, &[], None).unwrap();
    let xd = direct(&a, &b, &[], &[]);
    assert!(rep.converged);
    assert!(rel_diff(&x, &xd) < 1e-9);
    assert!(rep.convergence_factor < 0.2, "factor {}", rep.convergence_factor);
    // With fixed vertices (Laplacian alone, no mass term).
    let (l, m) = system(&h, 0.0);
    let known = [0u32, 17, 5000];
    let vals = [1.0, 0.0, -1.0];
    let mg = Multigrid::new(&l, &h.p, &known, MultigridOptions::default()).unwrap();
    let (x, rep) = mg.solve(&m, &vals, None).unwrap();
    assert!(rep.converged);
    assert!(rel_diff(&x, &direct(&l, &m, &known, &vals)) < 1e-9);
    assert!(rep.convergence_factor < 0.2, "dirichlet factor {}", rep.convergence_factor);
    println!("icosphere: {:?} factor {:.3}", rep.level_sizes, rep.convergence_factor);
}

#[test]
fn dirichlet_three_fixed_vertices_like_reference_example_04() {
    // L z = M 1 with three vertices fixed (the reference's 04 example uses a Laplace-type system
    // with three fixed vertices).
    let Some(h) = spot_hierarchy(2) else { return };
    let g = Geometry::from_positions(&h.mesh, &h.pos).unwrap();
    let l = cotan_laplacian(&h.mesh, &g).unwrap();
    let m = lumped_mass_normalized(&h.mesh, &g).unwrap();
    let n = l.nrows() as u32;
    let known = [3u32, n / 2, n - 5];
    let vals = [0.0, 1.0, -0.5];
    let mg = Multigrid::new(&l, &h.p, &known, MultigridOptions::default()).unwrap();
    let (x, rep) = mg.solve(&m, &vals, None).unwrap();
    let xd = direct(&l, &m, &known, &vals);
    assert!(rep.converged, "{rep:?}");
    assert!(rel_diff(&x, &xd) < 1e-9, "rel diff {}", rel_diff(&x, &xd));
    for (&k, &v) in known.iter().zip(&vals) {
        assert_eq!(x[k as usize], v);
    }
    assert!(rep.convergence_factor < 0.7, "factor {}", rep.convergence_factor);
    assert_eq!(mg.info().n_known, 3);
    assert_eq!(mg.info().level_sizes[0], l.nrows() - 3);
}

#[test]
fn flat_grid_dirichlet_boundary_and_dropped_columns_cascade() {
    let h = flat_hierarchy();
    let (a, _) = system(&h, 0.0);
    let n = a.nrows();
    // Fix every vertex in the support of coarse column 7 of P_0, so that column is dropped.
    let p0 = &h.p[0];
    let support: Vec<u32> = (0..n as u32)
        .filter(|&i| p0.get(i as usize, 7).abs() > 1e-15)
        .collect();
    assert!(support.len() >= 2);
    let vals: Vec<f64> = (0..support.len()).map(|i| i as f64 * 0.1).collect();
    let b = random_rhs(n, 3);
    let mg = Multigrid::new(&a, &h.p, &support, MultigridOptions::default()).unwrap();
    assert!(mg.info().dropped_columns[0] >= 1, "{:?}", mg.info());
    let (x, rep) = mg.solve(&b, &vals, None).unwrap();
    let xd = direct(&a, &b, &support, &vals);
    assert!(rep.converged);
    assert!(rel_diff(&x, &xd) < 1e-8, "rel diff {}", rel_diff(&x, &xd));
    // Dropping shrinks the coarse levels consistently.
    let sizes = &mg.info().level_sizes;
    for (k, p) in mg.prolongations().iter().enumerate() {
        assert_eq!(p.nrows(), sizes[k]);
        assert_eq!(p.ncols(), sizes[k + 1]);
    }
}

#[test]
fn every_level_is_symmetric_positive_definite() {
    let h = flat_hierarchy();
    let (a, _) = system(&h, 0.5);
    let known = [0u32, 100, 2000];
    let mg = Multigrid::new(&a, &h.p, &known, MultigridOptions::default()).unwrap();
    let mut rng = SplitMix64::new(9);
    for (lv, m) in mg.level_matrices().iter().enumerate() {
        let scale = m.values().iter().fold(0.0f64, |s, v| s.max(v.abs()));
        assert!(m.symmetry_defect().unwrap() < 1e-12 * scale, "level {lv} symmetry");
        Symbolic::analyze(m, Method::Llt)
            .unwrap()
            .factor(m, Parallelism::Sequential)
            .unwrap_or_else(|e| panic!("level {lv} is not positive definite: {e}"));
        for _ in 0..3 {
            let x: Vec<f64> = (0..m.nrows()).map(|_| rng.uniform(-1.0, 1.0)).collect();
            let mut y = vec![0.0; x.len()];
            m.mul_vec(&x, &mut y).unwrap();
            let q: f64 = x.iter().zip(&y).map(|(a, b)| a * b).sum();
            assert!(q > 0.0, "level {lv}: x'Ax = {q}");
        }
    }
    assert_eq!(mg.info().coarsest_shift, 0.0);
}

#[test]
fn vcycle_is_a_symmetric_operator_and_pcg_converges() {
    let h = flat_hierarchy();
    let (a, _) = system(&h, 1.0);
    let b = random_rhs(a.nrows(), 4);
    let mg = Multigrid::new(&a, &h.p, &[], MultigridOptions::default()).unwrap();
    let n = a.nrows();
    let (u, v) = (random_rhs(n, 5), random_rhs(n, 6));
    let (mut mu, mut mv) = (vec![0.0; n], vec![0.0; n]);
    mg.apply_preconditioner(&u, &mut mu).unwrap();
    mg.apply_preconditioner(&v, &mut mv).unwrap();
    let uv: f64 = u.iter().zip(&mv).map(|(a, b)| a * b).sum();
    let vu: f64 = v.iter().zip(&mu).map(|(a, b)| a * b).sum();
    assert!((uv - vu).abs() < 1e-10 * uv.abs().max(1.0), "{uv} vs {vu}");
    let uu: f64 = u.iter().zip(&mu).map(|(a, b)| a * b).sum();
    assert!(uu > 0.0);
    let (x, rep) = mg.pcg(&b, &[], None, 100).unwrap();
    let xd = direct(&a, &b, &[], &[]);
    assert!(rep.converged, "{rep:?}");
    assert!(rel_diff(&x, &xd) < 1e-8);
    let (_, plain) = mg.solve(&b, &[], None).unwrap();
    println!("pcg {} iterations, V-cycles {}", rep.cycles, plain.cycles);
    assert!(rep.cycles <= plain.cycles);
}

#[test]
fn singular_coarsest_gets_a_reported_shift() {
    // Pure Neumann Laplacian (singular) with one fixed vertex is SPD on the fine level, but the
    // hierarchy with the fixed vertex's column dropped keeps SPD; to make the coarsest singular
    // use a single level without fixed vertices.
    let g = generators::grid_disk(6, 6).unwrap();
    let m = g.topology().unwrap();
    let geo = Geometry::from_positions(&m, &g.positions).unwrap();
    let l = cotan_laplacian(&m, &geo).unwrap();
    let mg = Multigrid::new(&l, &[], &[], MultigridOptions::default()).unwrap();
    assert!(mg.info().coarsest_shift > 0.0, "{:?}", mg.info());
    let strict = MultigridOptions {
        coarsest_shift: 0.0,
        ..Default::default()
    };
    assert!(matches!(
        Multigrid::new(&l, &[], &[], strict),
        Err(MultigridError::Coarsest { .. })
    ));
}

#[test]
fn bad_input_is_an_error_not_a_panic() {
    let h = flat_hierarchy();
    let (a, _) = system(&h, 1.0);
    let opts = MultigridOptions::default;
    // Not square.
    let rect = SparseMatrix::from_triplets(3, 4, &[(0, 0, 1.0)]).unwrap();
    assert!(matches!(
        Multigrid::new(&rect, &[], &[], opts()),
        Err(MultigridError::NotSquare { rows: 3, cols: 4 })
    ));
    // P with the wrong number of rows.
    let bad_p = SparseMatrix::from_triplets(a.nrows() - 1, 3, &[(0, 0, 1.0)]).unwrap();
    assert!(matches!(
        Multigrid::new(&a, &[bad_p], &[], opts()),
        Err(MultigridError::SizeMismatch { .. })
    ));
    // Second P not matching the first's columns.
    let mut ps = h.p.clone();
    ps.swap(0, 1);
    assert!(matches!(
        Multigrid::new(&a, &ps, &[], opts()),
        Err(MultigridError::SizeMismatch { .. })
    ));
    // A level with zero columns.
    let zero = SparseMatrix::from_triplets(a.nrows(), 0, &[]).unwrap();
    assert!(matches!(
        Multigrid::new(&a, &[zero], &[], opts()),
        Err(MultigridError::EmptyLevel { level: 0, .. })
    ));
    // Fixed vertices out of range, repeated, all.
    assert!(matches!(
        Multigrid::new(&a, &h.p, &[a.nrows() as u32], opts()),
        Err(MultigridError::BadKnown { .. })
    ));
    assert!(matches!(
        Multigrid::new(&a, &h.p, &[5, 5], opts()),
        Err(MultigridError::BadKnown { .. })
    ));
    let all: Vec<u32> = (0..a.nrows() as u32).collect();
    assert!(matches!(
        Multigrid::new(&a, &h.p, &all, opts()),
        Err(MultigridError::BadKnown { .. })
    ));
    // Every free vertex of the support removed: all columns dropped.
    let tiny = SparseMatrix::from_triplets(2, 2, &[(0, 0, 2.0), (1, 1, 2.0)]).unwrap();
    let p = SparseMatrix::from_triplets(2, 1, &[(0, 0, 1.0)]).unwrap();
    assert!(matches!(
        Multigrid::new(&tiny, &[p], &[0], opts()),
        Err(MultigridError::EmptyLevel { .. })
    ));
    // Not symmetric; non-positive diagonal; bad options.
    let ns = SparseMatrix::from_triplets(2, 2, &[(0, 0, 2.0), (1, 1, 2.0), (0, 1, 1.0)]).unwrap();
    assert!(matches!(
        Multigrid::new(&ns, &[], &[], opts()),
        Err(MultigridError::NotSymmetric { .. })
    ));
    let neg = SparseMatrix::from_triplets(2, 2, &[(0, 0, -2.0), (1, 1, 2.0)]).unwrap();
    assert!(matches!(
        Multigrid::new(&neg, &[], &[], opts()),
        Err(MultigridError::BadDiagonal { level: 0, row: 0, .. })
    ));
    let bad = MultigridOptions {
        tolerance: f64::NAN,
        ..opts()
    };
    assert!(matches!(
        Multigrid::new(&a, &h.p, &[], bad),
        Err(MultigridError::BadOption { .. })
    ));
    // Solves.
    let mg = Multigrid::new(&a, &h.p, &[1, 2], opts()).unwrap();
    let b = vec![0.0; a.nrows()];
    assert!(matches!(
        mg.solve(&b[1..], &[0.0, 0.0], None),
        Err(MultigridError::SizeMismatch { .. })
    ));
    assert!(matches!(
        mg.solve(&b, &[0.0], None),
        Err(MultigridError::SizeMismatch { .. })
    ));
    assert!(matches!(
        mg.solve(&b, &[0.0, 0.0], Some(&b[1..])),
        Err(MultigridError::SizeMismatch { .. })
    ));
    let mut nan = b.clone();
    nan[4] = f64::NAN;
    assert!(matches!(
        mg.solve(&nan, &[0.0, 0.0], None),
        Err(MultigridError::NonFinite { .. })
    ));
    let mut x = vec![0.0; 3];
    assert!(mg.vcycle(&b, &mut x).is_err());
    // Messages are readable.
    let msg = Multigrid::new(&rect, &[], &[], opts()).unwrap_err().to_string();
    assert!(msg.contains("3 x 4"));
}

#[test]
fn zero_rhs_and_hot_start() {
    let h = flat_hierarchy();
    let (a, _) = system(&h, 1.0);
    let mg = Multigrid::new(&a, &h.p, &[], MultigridOptions::default()).unwrap();
    let b = vec![0.0; a.nrows()];
    let (x, rep) = mg.solve(&b, &[], None).unwrap();
    assert!(rep.converged && rep.cycles == 0);
    assert!(x.iter().all(|&v| v == 0.0));
    // Starting from the exact solution takes no cycle.
    let b = random_rhs(a.nrows(), 8);
    let xd = direct(&a, &b, &[], &[]);
    let (_, rep) = mg.solve(&b, &[], Some(&xd)).unwrap();
    assert!(rep.cycles <= 1, "{}", rep.cycles);
}
