//! Meshes with several connected components: the solution of the canonical problem on the whole
//! mesh (masses divided by the total area, one additive constant of `u` per closed component)
//! against separate solves of the components with the equivalent `λ`, against the generic
//! convex solver, under relabeling, with steering per component, and for every entry point.

mod common;

use common::{
    Features, Reference, TestMesh, max_abs_diff, random_features, random_mesh, reference_solve,
};
use ddg_core::generators::SplitMix64;
use ddg_core::sparse::ThreadPool;
use ddg_optimal_cones::{
    BoundaryTreatment, Cone, ConeError, MethodUsed, Parallelism, Solution, SolveMethod,
    SolveOptions, Status, certify, certify_with_boundary, evaluate_cones, factor_once_support,
    solve, solve_with_boundary,
};

/// The disjoint union of meshes, in this order (so component `c` is `parts[c]`), with the vertex
/// offset of every part.
fn union(parts: &[&TestMesh]) -> (TestMesh, Vec<usize>) {
    let mut positions = Vec::new();
    let mut faces = Vec::new();
    let mut offsets = Vec::new();
    for part in parts {
        let offset = positions.len();
        offsets.push(offset);
        positions.extend_from_slice(&part.positions);
        faces.extend(
            part.mesh
                .faces()
                .iter()
                .map(|f| f.map(|v| v + offset as u32)),
        );
    }
    let name = parts
        .iter()
        .map(|p| p.name.as_str())
        .collect::<Vec<_>>()
        .join("+");
    (TestMesh::from_positions(&name, positions, faces), offsets)
}

/// The features of the union, restricted to one part and renumbered, with `λ` replaced.
fn restrict(f: &Features, offset: usize, n: usize, lambda: f64) -> Features {
    let slice = |a: &Option<Vec<f64>>| a.as_ref().map(|a| a[offset..offset + n].to_vec());
    Features {
        lambda,
        distortion_weights: slice(&f.distortion_weights),
        cone_costs: slice(&f.cone_costs),
        lower_bounds: slice(&f.lower_bounds),
        upper_bounds: slice(&f.upper_bounds),
        prescribed_cones: f
            .prescribed_cones
            .iter()
            .filter(|c| (offset..offset + n).contains(&(c.vertex as usize)))
            .map(|c| Cone {
                vertex: c.vertex - offset as u32,
                angle: c.angle,
            })
            .collect(),
    }
}

/// The features of the union from features of the parts (each for its own part).
fn concat(parts: &[&TestMesh], features: &[Features], lambda: f64) -> Features {
    let mut out = Features::plain(lambda);
    let mut offset = 0u32;
    let cat = |get: &dyn Fn(&Features) -> &Option<Vec<f64>>, default: f64| {
        if features.iter().all(|f| get(f).is_none()) {
            return None;
        }
        let mut v = Vec::new();
        for (p, f) in parts.iter().zip(features) {
            match get(f) {
                Some(a) => v.extend_from_slice(a),
                None => v.extend(std::iter::repeat_n(default, p.mesh.n_vertices())),
            }
        }
        Some(v)
    };
    out.distortion_weights = cat(&|f| &f.distortion_weights, 1.0);
    out.cone_costs = cat(&|f| &f.cone_costs, 1.0);
    out.lower_bounds = cat(&|f| &f.lower_bounds, f64::NEG_INFINITY);
    out.upper_bounds = cat(&|f| &f.upper_bounds, f64::INFINITY);
    for (p, f) in parts.iter().zip(features) {
        out.prescribed_cones
            .extend(f.prescribed_cones.iter().map(|c| Cone {
                vertex: c.vertex + offset,
                angle: c.angle,
            }));
        offset += p.mesh.n_vertices() as u32;
    }
    out
}

fn area(mesh: &TestMesh) -> f64 {
    mesh.geometry.total_area()
}

/// Solves the union and every part alone with `λ · A / A_c`, and checks that they agree: `µ`
/// and `u` per vertex, the adjoint up to the factor `A_c / A`, the objective as the area
/// weighted sum, and the report of every component.
fn check_union_against_parts(
    parts: &[&TestMesh],
    features: &Features,
    boundary: BoundaryTreatment,
    options: &SolveOptions<'_>,
) -> Option<Solution> {
    let (whole, offsets) = union(parts);
    let total_area = area(&whole);
    let solved = match solve_with_boundary(&features.problem(&whole), boundary, options) {
        Ok(s) => s,
        Err(ConeError::Infeasible { .. }) => return None,
        Err(e) => panic!("{}: {e}", whole.name),
    };
    assert_eq!(solved.report.status, Status::Converged, "{}", whole.name);
    assert!(solved.report.certificate.holds(1e-9));
    assert_eq!(solved.report.components.len(), parts.len());
    // The certificate of the whole problem is reproduced from the triple.
    let again = certify_with_boundary(
        &features.problem(&whole),
        boundary,
        &solved.log_scale_optimal,
        &solved.adjoint,
        &solved.measure,
    )
    .unwrap();
    assert_eq!(again, solved.report.certificate);

    let mut objective = 0.0;
    for (c, part) in parts.iter().enumerate() {
        let n = part.mesh.n_vertices();
        let fraction = area(part) / total_area;
        let report = &solved.report.components[c];
        assert_eq!(report.component, c as u32);
        assert_eq!(report.n_vertices, n);
        assert!((report.area_fraction - fraction).abs() < 1e-13);
        let lambda_c = features.lambda * total_area / area(part);
        assert!((report.equivalent_lambda - lambda_c).abs() <= 1e-12 * lambda_c.max(1.0));
        assert_eq!(report.report.status, Status::Converged);
        assert!(report.report.certificate.holds(1e-9));
        let alone = solve_with_boundary(
            &restrict(features, offsets[c], n, lambda_c).problem(part),
            boundary,
            options,
        )
        .unwrap();
        assert_eq!(alone.report.status, Status::Converged);
        let range = offsets[c]..offsets[c] + n;
        let tag = format!("{} component {c}", whole.name);
        let dmu = max_abs_diff(&solved.measure[range.clone()], &alone.measure);
        let du = max_abs_diff(
            &solved.log_scale_optimal[range.clone()],
            &alone.log_scale_optimal,
        );
        assert!(dmu < 1e-7, "{tag}: |Δµ| = {dmu:e}");
        assert!(du < 1e-7, "{tag}: |Δu| = {du:e}");
        // The adjoint of the part alone is A / A_c times the adjoint of the whole problem
        // (unless its constant is free: then compare up to that constant).
        if report.adjoint_shift == (0.0, 0.0) && alone.adjoint_shift == (0.0, 0.0) {
            let scaled: Vec<f64> = alone.adjoint.iter().map(|p| p * fraction).collect();
            let dphi = max_abs_diff(&solved.adjoint[range.clone()], &scaled);
            let size = scaled.iter().fold(1e-3f64, |m, p| m.max(p.abs()));
            assert!(dphi < 1e-6 * size, "{tag}: |Δφ| = {dphi:e}");
        }
        let share = fraction * alone.objective;
        assert!(
            (report.objective - share).abs() <= 1e-9 * share.abs().max(1e-6),
            "{tag}: objective {} against {share}",
            report.objective
        );
        objective += report.objective;
    }
    assert!((solved.objective - objective).abs() <= 1e-12 * objective.abs().max(1.0));
    Some(solved)
}

fn sphere(seed: u64) -> TestMesh {
    random_mesh(3, seed)
}

fn disk(seed: u64) -> TestMesh {
    random_mesh(0, seed)
}

#[test]
fn union_of_two_closed_meshes_equals_the_separate_solves() {
    let (a, b) = (sphere(1), random_mesh(4, 2)); // a sphere and a torus
    for lambda in [0.02, 0.1, 0.5] {
        check_union_against_parts(
            &[&a, &b],
            &Features::plain(lambda),
            BoundaryTreatment::Fixed,
            &SolveOptions::default(),
        )
        .unwrap();
    }
}

#[test]
fn union_of_a_closed_mesh_and_a_mesh_with_boundary() {
    let (a, b) = (disk(3), sphere(4));
    for lambda in [0.005, 0.05, 0.5] {
        let s = check_union_against_parts(
            &[&a, &b],
            &Features::plain(lambda),
            BoundaryTreatment::Fixed,
            &SolveOptions::default(),
        )
        .unwrap();
        assert!(!s.report.components[0].closed && s.report.components[1].closed);
        // Gauss–Bonnet on the closed component alone: its angles sum to 4π.
        let n0 = a.mesh.n_vertices();
        let sum: f64 = s.measure[n0..].iter().sum();
        assert!((sum - 4.0 * std::f64::consts::PI).abs() < 1e-9, "{sum}");
    }
}

#[test]
fn three_components_and_every_boundary_treatment() {
    let (a, b, c) = (disk(5), random_mesh(1, 6), sphere(7)); // disk, annulus, sphere
    for boundary in [
        BoundaryTreatment::Fixed,
        BoundaryTreatment::Free,
        BoundaryTreatment::Convex,
        BoundaryTreatment::Penalized,
    ] {
        check_union_against_parts(
            &[&a, &b, &c],
            &Features::plain(0.05),
            boundary,
            &SolveOptions::default(),
        );
    }
}

#[test]
fn the_normalization_is_that_of_the_whole_mesh() {
    // Solving a component alone with the same λ is a different problem: on a small component
    // it puts cones where the whole problem does not (or the other way round). The whole problem
    // is the one with λ · A / A_c per component.
    let big = random_mesh(0, 11);
    let small_positions: Vec<[f64; 3]> = big
        .positions
        .iter()
        .map(|p| [p[0] * 0.2 + 9.0, p[1] * 0.2, p[2] * 0.2])
        .collect();
    let small = TestMesh::from_positions("small", small_positions, big.mesh.faces().to_vec());
    let (whole, _) = union(&[&big, &small]);
    let lambda = 0.004;
    let s = solve(
        &Features::plain(lambda).problem(&whole),
        &SolveOptions::default(),
    )
    .unwrap();
    let alone_same = solve(
        &Features::plain(lambda).problem(&small),
        &SolveOptions::default(),
    )
    .unwrap();
    let n0 = big.mesh.n_vertices();
    let in_whole = s.measure[n0..].iter().filter(|m| **m != 0.0).count();
    let alone = alone_same.measure.iter().filter(|m| **m != 0.0).count();
    // The small copy is scaled: alone it is the big mesh's problem (scale invariance), so it
    // has the big mesh's cones; inside the union its distortion counts 25 times less.
    let in_big = s.measure[..n0].iter().filter(|m| **m != 0.0).count();
    eprintln!("cones: big {in_big}, small alone {alone}, small in the union {in_whole}");
    assert_eq!(alone, in_big, "scale invariance");
    assert!(
        in_whole < alone,
        "{in_whole} cones in the union, {alone} alone"
    );
}

#[test]
fn agrees_with_the_generic_convex_solver() {
    let mut checked = 0;
    for seed in 0..40u64 {
        let (a, b) = (
            random_mesh(seed as usize, 100 + seed),
            random_mesh(seed as usize + 3, 200 + seed),
        );
        let (whole, _) = union(&[&a, &b]);
        // Positive weights: a unique minimizer.
        let f = random_features(&whole, 300 + seed, common::FeatureMask::UNIQUE);
        let ours = solve(&f.problem(&whole), &SolveOptions::default());
        match (ours, reference_solve(&whole, &f)) {
            (
                Ok(s),
                Reference::Solved {
                    objective,
                    log_scale,
                    measure,
                },
            ) => {
                assert_eq!(s.report.status, Status::Converged, "{}", whole.name);
                let rel = (s.objective - objective).abs() / objective.abs().max(1e-3);
                assert!(
                    rel < 1e-7,
                    "{}: objective {} against {objective}",
                    whole.name,
                    s.objective
                );
                let dmu = max_abs_diff(&s.measure, &measure);
                let du = max_abs_diff(&s.log_scale_optimal, &log_scale);
                assert!(
                    dmu < 1e-4 && du < 1e-4,
                    "{}: |Δµ| {dmu:e}, |Δu| {du:e}",
                    whole.name
                );
                checked += 1;
            }
            (Err(ConeError::Infeasible { component, .. }), Reference::Infeasible) => {
                assert!(component.is_some());
            }
            (_, Reference::Unreliable(why)) => eprintln!("{}: Clarabel {why}", whole.name),
            (ours, theirs) => panic!("{}: {ours:?} against {theirs:?}", whole.name),
        }
    }
    eprintln!("{checked} of 40 disconnected problems compared with Clarabel");
    assert!(checked >= 20, "only {checked} cases compared");
}

#[test]
fn steering_works_per_component() {
    let (a, b) = (disk(21), sphere(22));
    let na = a.mesh.n_vertices();
    let nb = b.mesh.n_vertices();
    // Disk: a free site, a forbidden region, positive angles only. Sphere: a prescribed cone,
    // bounds, and weights.
    let interior_a = a.interior();
    let mut costs_a = vec![1.0; na];
    costs_a[interior_a[0]] = 0.0;
    for &v in interior_a.iter().skip(1).take(3) {
        costs_a[v] = f64::INFINITY;
    }
    let fa = Features {
        cone_costs: Some(costs_a),
        lower_bounds: Some(vec![0.0; na]),
        ..Features::plain(1.0)
    };
    let fb = Features {
        upper_bounds: Some(vec![2.0; nb]),
        distortion_weights: Some((0..nb).map(|v| 1.0 + (v % 3) as f64).collect()),
        prescribed_cones: vec![Cone {
            vertex: 2,
            angle: 1.5,
        }],
        ..Features::plain(1.0)
    };
    let f = concat(&[&a, &b], &[fa, fb], 0.02);
    let s = check_union_against_parts(
        &[&a, &b],
        &f,
        BoundaryTreatment::Fixed,
        &SolveOptions::default(),
    )
    .unwrap();
    assert_eq!(s.measure[na + 2], 1.5);
    for &v in interior_a.iter().skip(1).take(3) {
        assert_eq!(s.measure[v], 0.0);
    }
    assert!(s.measure[..na].iter().all(|&m| m >= 0.0));
    assert!(s.measure[na..].iter().all(|&m| m <= 2.0));
    // Infeasible steering on one closed component names it.
    let g = Features {
        cone_costs: Some(
            (0..na + nb)
                .map(|v| if v < na { 1.0 } else { f64::INFINITY })
                .collect(),
        ),
        ..Features::plain(0.1)
    };
    let (whole, _) = union(&[&a, &b]);
    match solve(&g.problem(&whole), &SolveOptions::default()) {
        Err(
            e @ ConeError::Infeasible {
                component: Some(1),
                euler_characteristic: 2,
                ..
            },
        ) => {
            assert!(e.to_string().contains("closed component 1"), "{e}");
        }
        other => panic!("{other:?}"),
    }
    let unweighted = Features {
        distortion_weights: Some(
            (0..na + nb)
                .map(|v| if v < na { 1.0 } else { 0.0 })
                .collect(),
        ),
        ..Features::plain(0.1)
    };
    assert_eq!(
        solve(&unweighted.problem(&whole), &SolveOptions::default()).unwrap_err(),
        ConeError::NoDistortionWeight { component: Some(1) }
    );
}

#[test]
fn relabeling_vertices_and_faces_gives_the_same_result() {
    let (a, b, c) = (disk(31), sphere(32), random_mesh(2, 33));
    let (whole, _) = union(&[&a, &b, &c]);
    let nv = whole.mesh.n_vertices();
    let mut rng = SplitMix64::new(77);
    // A random permutation of the vertices and of the faces, and a rotation of the corners.
    let mut perm: Vec<usize> = (0..nv).collect();
    for i in (1..nv).rev() {
        perm.swap(i, (rng.next_u64() as usize) % (i + 1));
    }
    let mut faces: Vec<[u32; 3]> = whole
        .mesh
        .faces()
        .iter()
        .map(|f| f.map(|v| perm[v as usize] as u32))
        .collect();
    for i in (1..faces.len()).rev() {
        faces.swap(i, (rng.next_u64() as usize) % (i + 1));
    }
    for (i, f) in faces.iter_mut().enumerate() {
        f.rotate_left(i % 3);
    }
    let mut positions = vec![[0.0; 3]; nv];
    for v in 0..nv {
        positions[perm[v]] = whole.positions[v];
    }
    let relabeled = TestMesh::from_positions("relabeled", positions, faces);
    for lambda in [0.01, 0.2] {
        let s = solve(
            &Features::plain(lambda).problem(&whole),
            &SolveOptions::default(),
        )
        .unwrap();
        let t = solve(
            &Features::plain(lambda).problem(&relabeled),
            &SolveOptions::default(),
        )
        .unwrap();
        assert_eq!(t.report.status, Status::Converged);
        let back = |x: &[f64]| -> Vec<f64> { (0..nv).map(|v| x[perm[v]]).collect() };
        assert!(max_abs_diff(&s.measure, &back(&t.measure)) < 1e-8);
        assert!(max_abs_diff(&s.log_scale_optimal, &back(&t.log_scale_optimal)) < 1e-8);
        assert!((s.objective - t.objective).abs() < 1e-11 * s.objective.max(1.0));
        // The components are numbered by their lowest face, so their order may differ.
        let mut x: Vec<usize> = s.report.components.iter().map(|c| c.n_vertices).collect();
        let mut y: Vec<usize> = t.report.components.iter().map(|c| c.n_vertices).collect();
        x.sort_unstable();
        y.sort_unstable();
        assert_eq!(x, y);
    }
}

#[test]
fn methods_threads_and_entry_points_agree() {
    let (a, b, c) = (disk(41), sphere(42), random_mesh(4, 43));
    let (whole, offsets) = union(&[&a, &b, &c]);
    let f = Features::plain(0.03);
    let problem = f.problem(&whole);
    assert_eq!(
        factor_once_support(&problem, BoundaryTreatment::Fixed).unwrap(),
        None
    );
    let base = solve(&problem, &SolveOptions::default()).unwrap();
    assert_eq!(base.report.method.used, MethodUsed::FactorOnce);
    let reference = solve(
        &problem,
        &SolveOptions {
            method: SolveMethod::Reference,
            ..SolveOptions::default()
        },
    )
    .unwrap();
    assert_eq!(reference.report.method.used, MethodUsed::Reference);
    assert!(max_abs_diff(&base.measure, &reference.measure) < 1e-8);
    assert!(
        reference
            .report
            .components
            .iter()
            .all(|c| !c.report.stages.is_empty())
    );
    assert_eq!(
        reference.report.stages.len(),
        reference
            .report
            .components
            .iter()
            .map(|c| c.report.stages.len())
            .sum::<usize>()
    );
    // On a pool every component runs on one thread: the same bits as sequentially.
    let pool = ThreadPool::new(4).unwrap();
    let pooled = solve(
        &problem,
        &SolveOptions {
            parallelism: Parallelism::Pool(&pool),
            ..SolveOptions::default()
        },
    )
    .unwrap();
    assert_eq!(pooled.measure, base.measure);
    assert_eq!(pooled.log_scale_optimal, base.log_scale_optimal);
    assert_eq!(pooled.adjoint, base.adjoint);
    // The certificate of the returned triple.
    let cert = certify(
        &problem,
        &base.log_scale_optimal,
        &base.adjoint,
        &base.measure,
    )
    .unwrap();
    assert_eq!(cert, base.report.certificate);
    // A zero distortion weight in one component makes the factor-once method decline, naming a
    // vertex of the whole mesh.
    let mut w = vec![1.0; whole.mesh.n_vertices()];
    let v = offsets[1] + 3;
    w[v] = 0.0;
    let weighted = Features {
        distortion_weights: Some(w),
        ..Features::plain(0.03)
    };
    assert_eq!(
        factor_once_support(&weighted.problem(&whole), BoundaryTreatment::Fixed).unwrap(),
        Some(ddg_optimal_cones::Unsupported::ZeroDistortionWeight { vertex: v as u32 })
    );
    // evaluate_cones on the union: u per component, the adjoint of the whole problem.
    let e = evaluate_cones(&problem, &base.cones, Parallelism::Sequential).unwrap();
    assert!(max_abs_diff(&e.log_scale, &base.log_scale_optimal) < 1e-8);
    assert!((e.distortion - base.distortion_optimal).abs() < 1e-10);
    assert_eq!(e.component_defects.len(), 3);
    assert!(e.gauss_bonnet_defect.abs() < 1e-9);
    let gap = base.measure.iter().zip(&e.adjoint).count();
    assert_eq!(gap, whole.mesh.n_vertices());
    // Cones missing on the sphere: its defect is reported, the others are zero.
    let without: Vec<Cone> = base
        .cones
        .iter()
        .filter(|c| (c.vertex as usize) < offsets[1] || (c.vertex as usize) >= offsets[2])
        .copied()
        .collect();
    let e = evaluate_cones(&problem, &without, Parallelism::Sequential).unwrap();
    assert!((e.component_defects[1] + 4.0 * std::f64::consts::PI).abs() < 1e-9);
    assert_eq!(e.component_defects[0], 0.0);
    assert_eq!(e.gauss_bonnet_defect, e.component_defects[1]);
    // Merging works per component; the extracted cones' scale factor per component.
    let merged = solve(
        &problem,
        &SolveOptions {
            merge_cone_clusters: true,
            cone_threshold: 1e-3,
            ..SolveOptions::default()
        },
    )
    .unwrap();
    let check = evaluate_cones(&problem, &merged.cones, Parallelism::Sequential).unwrap();
    assert!(max_abs_diff(&check.log_scale, &merged.log_scale) < 1e-8);
    assert!((check.distortion - merged.distortion).abs() < 1e-10);
}
