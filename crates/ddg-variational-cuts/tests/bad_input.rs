//! Bad input is reported through `CutError` with a message that names the vertex, edge or face
//! concerned; nothing panics.

mod common;

use common::{Shape, bumpy, sphere};
use ddg_core::generators::grid_disk;
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CurvatureSource, Cut, CutError, CutProblem, Energy, EvalOptions, Evaluator, FluxRecovery,
    PointwiseEnergy, RegionWeights, ScaleWeight, UnitAreaScale,
};
use proptest::prelude::*;

fn level_set(p: [f64; 3]) -> f64 {
    0.3 * p[0] + 0.5 * p[1] + 0.4 * p[2] - 0.35
}

/// The error of an evaluation that must fail.
fn error(shape: &Shape, phi: &[f64], problem: &CutProblem, options: &EvalOptions<'_>) -> CutError {
    shape
        .evaluator()
        .evaluate(phi, problem, options)
        .expect_err("the evaluation must fail")
}

fn dirichlet() -> CutProblem {
    CutProblem::new(Energy::Dirichlet, 0.1)
}

#[test]
fn geometry_of_another_mesh() {
    let (a, b) = (sphere(1), bumpy(3, 0.1));
    let e = Evaluator::new(&a.mesh, &b.geometry).unwrap_err();
    assert!(matches!(e, CutError::Geometry(_)));
    assert!(
        e.to_string()
            .starts_with("geometry does not belong to this mesh")
    );
}

#[test]
fn level_set_function() {
    let shape = sphere(1);
    let n = shape.positions.len();
    let phi = shape.sample(level_set);
    let options = EvalOptions::default();

    let e = error(&shape, &phi[..n - 1], &dirichlet(), &options);
    assert_eq!(
        e.to_string(),
        format!(
            "level-set function has {} entries, but there are {n} vertices",
            n - 1
        )
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut broken = phi.clone();
        broken[17] = bad;
        let e = error(&shape, &broken, &dirichlet(), &options);
        assert!(matches!(e, CutError::NonFiniteLevelSet { vertex: 17, .. }));
        assert_eq!(
            e.to_string(),
            format!("level-set value of vertex 17 is {bad}; it must be finite")
        );
        assert!(Cut::from_level_set(&shape.mesh, &broken, 1e-3).is_err());
    }
    // The helper for positions checks its input too.
    let cut = Cut::from_level_set(&shape.mesh, &phi, 1e-3).unwrap();
    let e = cut
        .point_positions(&shape.mesh, &shape.positions[1..])
        .unwrap_err();
    assert!(matches!(
        e,
        CutError::WrongLength {
            what: "positions",
            ..
        }
    ));
    let other = bumpy(3, 0.1);
    assert!(cut.point_positions(&other.mesh, &other.positions).is_err());
}

#[test]
fn problem_description() {
    let shape = bumpy(4, 0.2);
    let n = shape.positions.len();
    let nf = shape.mesh.n_faces();
    let phi = shape.sample(level_set);
    let options = EvalOptions::default();
    let check = |problem: CutProblem, message: &str| {
        assert_eq!(error(&shape, &phi, &problem, &options).to_string(), message);
    };

    check(
        CutProblem {
            regions: vec![RegionWeights::default(); 3],
            ..dirichlet()
        },
        "the problem describes 3 regions; this version supports exactly 2",
    );
    for (bad, shown) in [(-1.0, "-1"), (f64::NAN, "NaN"), (f64::INFINITY, "inf")] {
        let mut problem = dirichlet();
        problem.regions[1].distortion = bad;
        check(
            problem,
            &format!(
                "distortion weight of region 1 is {shown}; it must be finite and non-negative"
            ),
        );
        check(
            CutProblem::new(Energy::Dirichlet, bad),
            &format!("length weight is {shown}; it must be finite and non-negative"),
        );
        let mut eta = vec![0.5; n];
        eta[9] = bad;
        check(
            CutProblem {
                penalty: Some(eta),
                ..dirichlet()
            },
            &format!("penalty of vertex 9 is {shown}; it must be finite and non-negative"),
        );
    }
    let mut problem = dirichlet();
    problem.regions[0].area = f64::NAN;
    check(problem, "area weight of region 0 is NaN; it must be finite");
    // A negative area weight is legal.
    let mut problem = dirichlet();
    problem.regions[0].area = -3.0;
    assert!(shape.evaluator().evaluate(&phi, &problem, &options).is_ok());

    check(
        CutProblem {
            penalty: Some(vec![0.0; n + 2]),
            ..dirichlet()
        },
        &format!("penalty has {} entries, but there are {n} vertices", n + 2),
    );
    for bad in [0.0, -2.0, f64::NAN, f64::INFINITY] {
        let mut a = vec![1.0; n];
        a[11] = bad;
        check(
            CutProblem::new(Energy::ScaledDirichlet(ScaleWeight::PerVertex(a)), 0.0),
            &format!(
                "scaled Dirichlet weight of vertex 11 is {bad}; it must be finite and positive"
            ),
        );
    }
    check(
        CutProblem::new(
            Energy::ScaledDirichlet(ScaleWeight::PerVertex(vec![1.0; 3])),
            0.0,
        ),
        &format!("scaled Dirichlet weight has 3 entries, but there are {n} vertices"),
    );
    let feature = |r_eps: f64, smoothing: f64| {
        CutProblem::new(
            Energy::ScaledDirichlet(ScaleWeight::FeatureSize { r_eps, smoothing }),
            0.0,
        )
    };
    check(
        feature(0.0, 0.02),
        "feature size: r_eps must be positive and finite, got 0",
    );
    check(
        feature(f64::INFINITY, 0.02),
        "feature size: r_eps must be positive and finite, got inf",
    );
    check(
        feature(0.1, -0.02),
        "feature size: the smoothing constant must be finite and non-negative, got -0.02",
    );
    assert!(
        shape
            .evaluator()
            .feature_size_weight(-1.0, 0.02, ddg_variational_cuts::Parallelism::Sequential)
            .is_err()
    );

    check(
        CutProblem {
            curvature: CurvatureSource::Uniform(f64::NAN),
            ..dirichlet()
        },
        "uniform curvature override is NaN; it must be finite",
    );
    let mut k = vec![1.0; nf];
    k[5] = f64::NEG_INFINITY;
    check(
        CutProblem {
            curvature: CurvatureSource::PerFace(k),
            ..dirichlet()
        },
        "curvature override of face 5 is -inf; it must be finite",
    );
    check(
        CutProblem {
            curvature: CurvatureSource::PerFace(vec![1.0; nf - 1]),
            ..dirichlet()
        },
        &format!(
            "curvature override has {} entries, but there are {nf} faces",
            nf - 1
        ),
    );
}

#[test]
fn options() {
    let shape = bumpy(4, 0.2);
    let phi = shape.sample(level_set);
    let check = |options: EvalOptions<'_>, fragment: &str| {
        let e = error(&shape, &phi, &dirichlet(), &options);
        assert!(matches!(e, CutError::InvalidOptions { .. }));
        assert!(e.to_string().contains(fragment), "{e}");
    };
    for bad in [0.5, -1e-3, 0.7, f64::NAN] {
        check(
            EvalOptions {
                min_crossing_offset: bad,
                ..EvalOptions::default()
            },
            "min_crossing_offset must lie in [0, 0.5)",
        );
        assert!(Cut::from_level_set(&shape.mesh, &phi, bad).is_err());
    }
    for bad in [0.0, -5.0, f64::NAN] {
        check(
            EvalOptions {
                sliver_aspect_ratio: bad,
                ..EvalOptions::default()
            },
            "sliver_aspect_ratio must be positive",
        );
    }
    for bad in [-0.5, f64::NAN, f64::INFINITY] {
        check(
            EvalOptions {
                flux_recovery: FluxRecovery::Smoothed { length: bad },
                ..EvalOptions::default()
            },
            "smoothing length of the flux recovery must be finite and non-negative",
        );
    }
    // Legal extremes: no sliver at all, every vertex a sliver, many passes.
    for (ratio, passes) in [(f64::INFINITY, 3), (0.5, 4), (5.0, 0)] {
        let options = EvalOptions {
            sliver_aspect_ratio: ratio,
            smoothing_passes: passes,
            ..EvalOptions::paper()
        };
        let e = shape
            .evaluator()
            .evaluate(&phi, &dirichlet(), &options)
            .unwrap();
        let expect = if ratio.is_infinite() {
            0
        } else {
            e.cut.points.len()
        };
        if ratio != 5.0 {
            assert_eq!(e.report.n_sliver_points, expect);
        }
        assert!(e.speed.iter().all(|s| s.is_finite()));
    }
}

#[test]
fn pointwise_energy_that_is_not_finite() {
    let shape = bumpy(4, 0.2);
    let phi = shape.sample(level_set);
    let energy = PointwiseEnergy::new(|u| if u < 0.0 { f64::NAN } else { u }, |_| 1.0);
    let e = error(
        &shape,
        &phi,
        &CutProblem::new(Energy::Pointwise(energy), 0.0),
        &EvalOptions::default(),
    );
    match &e {
        CutError::NonFiniteEnergyDensity {
            what,
            value,
            log_scale,
            ..
        } => {
            assert_eq!(*what, "b");
            assert!(value.is_nan() && *log_scale < 0.0);
        }
        other => panic!("unexpected {other:?}"),
    }
    assert!(e.to_string().starts_with("pointwise energy: b("));
    let energy = PointwiseEnergy::new(|u| u, |_| f64::INFINITY);
    let e = error(
        &shape,
        &phi,
        &CutProblem::new(Energy::Pointwise(energy), 0.0),
        &EvalOptions::default(),
    );
    assert!(matches!(
        e,
        CutError::NonFiniteEnergyDensity { what: "b'", .. }
    ));
}

#[test]
fn crossing_on_a_vertex_without_a_minimal_offset() {
    // With min_crossing_offset = 0 nothing is moved, so a zero of φ at a vertex gives a
    // triangle without area; the error names the face and the two edges.
    let shape = Shape::from_generated("grid", grid_disk(2, 2).unwrap());
    let phi = shape.sample(|p| p[0] - 0.5);
    assert!(phi.contains(&0.0));
    let zero = EvalOptions {
        min_crossing_offset: 0.0,
        ..EvalOptions::default()
    };
    let e = error(&shape, &phi, &dirichlet(), &zero);
    match &e {
        CutError::DegenerateCutTriangle { face, vertices, .. } => {
            assert_eq!(shape.mesh.faces()[*face as usize], *vertices);
        }
        other => panic!("unexpected {other:?}"),
    }
    let message = e.to_string();
    assert!(
        message.contains("is cut into a triangle without positive area"),
        "{message}"
    );
    assert!(message.contains("at the fractions"), "{message}");
    // With the default offset the same input is fine, and the moved crossings are counted.
    let e = shape
        .evaluator()
        .evaluate(&phi, &dirichlet(), &EvalOptions::default())
        .unwrap();
    assert!(e.report.n_clamped_crossings > 0);
    assert_eq!(e.report.smallest_crossing_offset, 0.0);
    assert!(e.speed.iter().all(|s| s.is_finite()));
    // An extremely small offset still evaluates, but the mesh with geometry may be refused:
    // its edge lengths no longer satisfy the triangle inequality in floating point.
    let tiny = EvalOptions {
        min_crossing_offset: 1e-17,
        ..EvalOptions::default()
    };
    match shape.evaluator().evaluate(&phi, &dirichlet(), &tiny) {
        Ok(e) => {
            if let Err(error) = e.cut_mesh() {
                assert!(
                    matches!(error, CutError::DegenerateCutTriangle { .. }),
                    "{error}"
                );
            }
        }
        Err(error) => assert!(
            matches!(
                error,
                CutError::DegenerateCutTriangle { .. } | CutError::Sparse(_)
            ),
            "{error}"
        ),
    }
}

#[test]
fn helpers_reject_bad_arguments() {
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let e = UnitAreaScale::from_total_area(bad).unwrap_err();
        assert!(matches!(e, CutError::InvalidTotalArea { .. }), "{bad}");
    }
    let scale = UnitAreaScale::from_total_area(4.0).unwrap();
    assert_eq!(
        (scale.length_scale, scale.length_to_unit_area(3.0)),
        (0.5, 1.5)
    );
    assert_eq!(scale.length_from_unit_area(1.5), 3.0);
    assert_eq!(
        (scale.area_to_unit_area(2.0), scale.area_from_unit_area(0.5)),
        (0.5, 2.0)
    );

    let shape = sphere(1);
    let phi = shape.sample(level_set);
    let e = shape
        .evaluator()
        .evaluate(&phi, &dirichlet(), &EvalOptions::default())
        .unwrap();
    let n = e.triangles.n_vertices();
    let mut y = vec![0.0; n];
    let error = e
        .triangles
        .apply_laplacian(None, &vec![0.0; n - 1], &mut y)
        .unwrap_err();
    assert_eq!(
        error.to_string(),
        format!(
            "x has {} entries, but there are {n} cut-mesh vertices",
            n - 1
        )
    );
    assert!(
        e.triangles
            .apply_laplacian(Some(&[1.0]), &vec![0.0; n], &mut y)
            .is_err()
    );
    let cut_mesh = e.cut_mesh().unwrap();
    assert!(
        cut_mesh
            .positions(&e.cut, &shape.mesh, &shape.positions[2..])
            .is_err()
    );
}

#[test]
fn a_mesh_without_interior_vertices() {
    // A single triangle and a strip: there is nothing to solve, and nothing fails.
    let mesh = TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap();
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.3],
        [0.0, 1.0, 0.0],
    ];
    let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
    let mut evaluator = Evaluator::new(&mesh, &geometry).unwrap();
    assert!(evaluator.symbolic().is_none() && evaluator.vertex_of_unknown().is_empty());
    for energy in [Energy::Dirichlet, Energy::Hencky] {
        let e = evaluator
            .evaluate(
                &[-1.0, 2.0, -0.5, -1.5],
                &CutProblem::new(energy, 1.0),
                &EvalOptions::default(),
            )
            .unwrap();
        assert_eq!(e.cut.points.len(), 2);
        assert!(e.log_scale.iter().all(|&u| u == 0.0));
        assert!(e.speed.iter().all(|&s| s == 0.0));
        assert!(e.energy.total > 0.0 && e.energy.total == e.energy.length);
    }
    assert!(evaluator.factorization().is_none());
}

/// Level-set values that stress the crossing computation.
fn level_value() -> impl Strategy<Value = f64> {
    prop_oneof![
        4 => -1.0..1.0f64,
        1 => Just(0.0),
        1 => Just(-0.0),
        1 => prop_oneof![Just(1e-300), Just(-1e-300), Just(1e300), Just(-1e300), Just(f64::MAX), Just(f64::MIN_POSITIVE)],
        1 => (-1.0..1.0f64).prop_map(|x| x * 1e-12),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 48, failure_persistence: None, ..ProptestConfig::default() })]

    /// Any finite level-set function evaluates without panic and without a non-finite number,
    /// for every energy, on a closed mesh and on one with boundary.
    #[test]
    fn any_finite_level_set_function_evaluates(
        values in proptest::collection::vec(level_value(), 42),
        which in 0..2usize,
        energy in 0..3usize,
        paper in proptest::bool::ANY,
    ) {
        let shape = if which == 0 { sphere(1) } else { bumpy(5, 0.3) };
        let n = shape.positions.len();
        let phi: Vec<f64> = (0..n).map(|v| values[v % values.len()]).collect();
        let energy = match energy {
            0 => Energy::Dirichlet,
            1 => Energy::Hencky,
            _ => Energy::ScaledDirichlet(ScaleWeight::PerVertex(
                shape.sample(|p| 1.0 + 0.5 * p[0] * p[0]),
            )),
        };
        let options = if paper { EvalOptions::paper() } else { EvalOptions::default() };
        let e = shape
            .evaluator()
            .evaluate(&phi, &CutProblem::new(energy, 0.2), &options)
            .unwrap();
        prop_assert!(e.energy.total.is_finite());
        prop_assert!(e.log_scale.iter().all(|x| x.is_finite()));
        prop_assert!(e.normal_derivative.iter().all(|x| x.is_finite()));
        prop_assert!(e.gradient.iter().all(|x| x.is_finite()));
        prop_assert!(e.speed.iter().chain(&e.speed_raw).all(|x| x.is_finite()));
        prop_assert!(e.cut.points.iter().all(|p| p.t >= 1e-3 && p.t <= 1.0 - 1e-3));
        let area: f64 = e.energy.area.iter().chain(&e.report.excluded_area).sum();
        prop_assert!((area - 1.0).abs() < 1e-12);
        // The mesh with geometry can be built whenever the minimal offset is the default.
        let cut_mesh = e.cut_mesh().unwrap();
        prop_assert_eq!(cut_mesh.mesh.n_faces(), e.triangles.n_faces());
    }
}
