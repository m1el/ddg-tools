//! SPEC §6 item 5, the static part: the evaluation does not depend on the numbering of the
//! vertices and faces, on the length unit, or on the unit of the level-set function; and
//! equivalent descriptions of one energy give one result.

mod common;

use common::{Shape, bumpy, sphere};
use ddg_core::generators::SplitMix64;
use ddg_variational_cuts::{
    CutProblem, Energy, EvalOptions, Evaluation, PointwiseEnergy, RegionWeights, ScaleWeight,
};

fn level_set(p: [f64; 3]) -> f64 {
    0.42 * p[0] - 0.33 * p[1] + 0.61 * p[2] - 0.21 + 0.17 * p[0] * p[1]
}

fn weight(p: [f64; 3]) -> f64 {
    0.7 + 0.3 * p[0] * p[0] + 0.2 * (3.0 * p[1]).sin().powi(2)
}

fn penalty(p: [f64; 3]) -> f64 {
    0.2 + 0.1 * (2.0 * p[0] + p[2]).cos().powi(2)
}

/// The problems of the test, as functions of the shape (per-vertex data follows the vertices).
fn problems(shape: &Shape) -> Vec<CutProblem> {
    let regions = vec![
        RegionWeights {
            distortion: 0.8,
            area: 0.4,
        },
        RegionWeights {
            distortion: 1.3,
            area: -0.1,
        },
    ];
    [
        Energy::Dirichlet,
        Energy::Hencky,
        Energy::ScaledDirichlet(ScaleWeight::PerVertex(shape.sample(weight))),
        Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
            r_eps: 0.25,
            smoothing: 0.02,
        }),
    ]
    .into_iter()
    .map(|energy| CutProblem {
        regions: regions.clone(),
        penalty: Some(shape.sample(penalty)),
        ..CutProblem::new(energy, 0.3)
    })
    .collect()
}

fn evaluate(shape: &Shape, phi: &[f64], problem: &CutProblem) -> Evaluation {
    shape
        .evaluator()
        .evaluate(phi, problem, &EvalOptions::default())
        .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0)
}

/// Compares everything that has a meaning independent of the numbering. `map[v]` is the vertex
/// of `b`'s mesh that vertex `v` of `a`'s mesh became.
fn assert_same(a: (&Shape, &Evaluation), b: (&Shape, &Evaluation), map: &[usize], tolerance: f64) {
    let ((sa, ea), (sb, eb)) = (a, b);
    assert!(close(ea.energy.total, eb.energy.total, tolerance));
    assert!(close(ea.energy.length, eb.energy.length, tolerance));
    assert!(close(ea.energy.penalty, eb.energy.penalty, tolerance));
    for r in 0..2 {
        assert!(close(
            ea.energy.distortion[r],
            eb.energy.distortion[r],
            tolerance
        ));
        assert!(close(ea.energy.area[r], eb.energy.area[r], tolerance));
    }
    assert_eq!(ea.report.n_cut_points, eb.report.n_cut_points);
    assert_eq!(ea.report.n_sliver_points, eb.report.n_sliver_points);
    let scale = ea.log_scale.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    for v in 0..sa.mesh.n_vertices() {
        assert_eq!(ea.cut.vertex_region[v], eb.cut.vertex_region[map[v]]);
        assert!((ea.log_scale[v] - eb.log_scale[map[v]]).abs() <= tolerance * scale);
        if let (Some(pa), Some(pb)) = (&ea.adjoint, &eb.adjoint) {
            let scale = pa.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            assert!((pa[v] - pb[map[v]]).abs() <= tolerance * scale);
        }
    }
    // Cut points are identified by their edge.
    let largest = ea.speed.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    for (k, point) in ea.cut.points.iter().enumerate() {
        let [x, y] = sa.mesh.edge_vertices(point.edge as usize);
        let edge = sb.mesh.find_edge(map[x], map[y]).unwrap();
        let other = eb.cut.edge_point[edge] as usize;
        // The parameter is measured from the lower-numbered end.
        let t = if map[x] < map[y] {
            point.t
        } else {
            1.0 - point.t
        };
        assert!((eb.cut.points[other].t - t).abs() < 1e-12);
        assert!((ea.speed[k] - eb.speed[other]).abs() <= tolerance * largest);
        assert!((ea.speed_raw[k] - eb.speed_raw[other]).abs() <= tolerance * largest);
        // Per side: the copy in the same region.
        let (va, vb) = (
            ea.triangles.point_vertices_by_region(k),
            eb.triangles.point_vertices_by_region(other),
        );
        for side in 0..2 {
            let (da, db) = (
                ea.normal_derivative[va[side]],
                eb.normal_derivative[vb[side]],
            );
            assert!(close(da, db, tolerance), "{da} vs {db}");
            assert!(close(
                ea.dual_length[va[side]],
                eb.dual_length[vb[side]],
                tolerance
            ));
            assert!(close(
                ea.gradient[va[side]],
                eb.gradient[vb[side]],
                tolerance * largest.max(1.0)
            ));
        }
    }
}

#[test]
fn invariant_under_relabeling() {
    for shape in [sphere(2), bumpy(9, 0.2)] {
        let n = shape.positions.len();
        let mut rng = SplitMix64::new(7);
        // A random permutation of the vertices, a random order and rotation of the faces.
        let mut map: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            map.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        let mut positions = vec![[0.0; 3]; n];
        for v in 0..n {
            positions[map[v]] = shape.positions[v];
        }
        let mut faces: Vec<[u32; 3]> = shape
            .mesh
            .faces()
            .iter()
            .map(|f| {
                let g = f.map(|v| map[v as usize] as u32);
                let k = (rng.next_u64() % 3) as usize;
                [g[k], g[(k + 1) % 3], g[(k + 2) % 3]]
            })
            .collect();
        for i in (1..faces.len()).rev() {
            faces.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        let relabeled = Shape::new("relabeled", positions, faces);
        let (phi_a, phi_b) = (shape.sample(level_set), relabeled.sample(level_set));
        for (pa, pb) in problems(&shape).iter().zip(&problems(&relabeled)) {
            let (ea, eb) = (
                evaluate(&shape, &phi_a, pa),
                evaluate(&relabeled, &phi_b, pb),
            );
            assert_same((&shape, &ea), (&relabeled, &eb), &map, 1e-9);
        }
    }
}

#[test]
fn invariant_under_uniform_scaling() {
    for shape in [sphere(2), bumpy(9, 0.2)] {
        let identity: Vec<usize> = (0..shape.positions.len()).collect();
        let phi = shape.sample(level_set);
        for factor in [8.0, 0.37, 1e4, 1e-3] {
            let scaled = Shape::new(
                "scaled",
                shape
                    .positions
                    .iter()
                    .map(|p| p.map(|x| factor * x))
                    .collect(),
                shape.mesh.faces().to_vec(),
            );
            // The level-set function in the new unit; any positive multiple gives the same cut.
            let phi_scaled: Vec<f64> = phi.iter().map(|x| factor * x).collect();
            for problem in problems(&shape) {
                let a = evaluate(&shape, &phi, &problem);
                let b = evaluate(&scaled, &phi_scaled, &problem);
                // Everything in an evaluation is on the unit-area surface.
                assert_same((&shape, &a), (&scaled, &b), &identity, 1e-9);
                assert!(close(
                    b.scale.total_area,
                    factor * factor * a.scale.total_area,
                    1e-12
                ));
                assert!(close(
                    b.scale.length_scale * factor,
                    a.scale.length_scale,
                    1e-12
                ));
                // The cut mesh that is handed out is in the caller's unit.
                let (ma, mb) = (a.cut_mesh().unwrap(), b.cut_mesh().unwrap());
                assert_eq!(ma.mesh, mb.mesh);
                for (la, lb) in ma
                    .geometry
                    .edge_lengths()
                    .iter()
                    .zip(mb.geometry.edge_lengths())
                {
                    assert!(close(factor * la, *lb, 1e-12));
                }
                // The same function in another unit of its own is the same cut.
                let c = evaluate(&shape, &phi_scaled, &problem);
                assert_eq!(c.cut.vertex_region, a.cut.vertex_region);
                assert!(close(c.energy.total, a.energy.total, 1e-11));
            }
        }
    }
}

#[test]
fn equivalent_descriptions_of_an_energy_agree() {
    let shape = bumpy(10, 0.25);
    let phi = shape.sample(level_set);
    let n = shape.positions.len();
    let run = |energy: Energy| evaluate(&shape, &phi, &CutProblem::new(energy, 0.0));
    // b(u) = u² given as a pointwise energy is the Hencky energy.
    let hencky = run(Energy::Hencky);
    let pointwise = run(Energy::Pointwise(PointwiseEnergy::new(
        |u| u * u,
        |u| 2.0 * u,
    )));
    assert_eq!(hencky.energy, pointwise.energy);
    assert_eq!(hencky.speed, pointwise.speed);
    assert_eq!(hencky.adjoint, pointwise.adjoint);
    // Scaled Dirichlet with a constant weight is that constant times the Dirichlet energy.
    let dirichlet = run(Energy::Dirichlet);
    let scaled = run(Energy::ScaledDirichlet(ScaleWeight::PerVertex(vec![
        2.5;
        n
    ])));
    assert!(close(
        scaled.energy.total,
        2.5 * dirichlet.energy.total,
        1e-11
    ));
    let largest = dirichlet.speed.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    for (s, d) in scaled.speed.iter().zip(&dirichlet.speed) {
        assert!((s - 2.5 * d).abs() < 1e-9 * largest);
    }
    // A pointwise energy with a constant density is the area: b = 1 gives J_r = Area_r, an
    // adjoint that vanishes, and the gradient density b(0) = 1 on both sides.
    let constant = run(Energy::Pointwise(PointwiseEnergy::new(|_| 1.0, |_| 0.0)));
    for r in 0..2 {
        assert!(close(
            constant.energy.distortion[r],
            constant.energy.area[r],
            1e-12
        ));
    }
    assert!(constant.adjoint.as_ref().unwrap().iter().all(|&p| p == 0.0));
    assert!(constant.speed.iter().all(|&s| s == 0.0));
    // The same through the area weights: β = (1, 0) moves the cut like b = 1 in region 0 only.
    let area = evaluate(
        &shape,
        &phi,
        &CutProblem {
            regions: vec![
                RegionWeights {
                    distortion: 0.0,
                    area: 1.0,
                },
                RegionWeights {
                    distortion: 0.0,
                    area: 0.0,
                },
            ],
            ..CutProblem::new(Energy::Dirichlet, 0.0)
        },
    );
    assert!(area.speed.iter().all(|&s| s == 1.0));
    assert!(close(area.energy.total, area.energy.area[0], 1e-14));
}

#[test]
fn evaluations_are_deterministic_and_reusable() {
    // One evaluator serves many cuts and problems; the order of the calls does not matter.
    let shape = sphere(2);
    let mut evaluator = shape.evaluator();
    let phi = shape.sample(level_set);
    let other: Vec<f64> = phi.iter().map(|x| x + 0.4).collect();
    let all = problems(&shape);
    let first: Vec<Evaluation> = all
        .iter()
        .map(|p| {
            evaluator
                .evaluate(&phi, p, &EvalOptions::default())
                .unwrap()
        })
        .collect();
    for p in all.iter().rev() {
        evaluator
            .evaluate(&other, p, &EvalOptions::default())
            .unwrap();
    }
    for (p, a) in all.iter().zip(&first) {
        let b = evaluator
            .evaluate(&phi, p, &EvalOptions::default())
            .unwrap();
        assert_eq!(a.log_scale, b.log_scale);
        assert_eq!(a.speed, b.speed);
        assert_eq!(a.energy, b.energy);
    }
}
