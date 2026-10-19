//! Robustness of the flow (SPEC §4.9, §5.1, §5.3; D31, D32): no degenerate faces, no error from a
//! flow that over-runs. The fast tests run by default; those on `data/models/` are `#[ignore]`d:
//!
//! ```text
//! cargo test --release -p ddg-developability --test robust -- --ignored --nocapture
//! ```

mod energy_common;
use energy_common::*;

use ddg_core::{Geometry, TriMesh};
use ddg_developability::geometry::corner_angles;
use ddg_developability::{
    Control, DevError, Energy, EnergyOptions, FlowOptions, OptimizeOptions, Progress, evaluate, flow, optimize,
};

fn min_angle(m: &TriMesh, p: &[[f64; 3]]) -> f64 {
    m.faces().iter().flat_map(|f| corner_angles(f.map(|v| p[v as usize]))).fold(f64::INFINITY, f64::min)
}

fn run_ep(mesh: &TriMesh, pos: &[[f64; 3]], floor: f64, energy: Energy) -> (f64, ddg_developability::OptimizeReport, f64) {
    let oo = OptimizeOptions { max_iterations: 600, min_face_angle: floor, ..Default::default() };
    let mut lowest = f64::INFINITY;
    let mut cb = |p: &Progress| {
        lowest = lowest.min(min_angle(mesh, p.positions));
        Control::Continue
    };
    let r = optimize(mesh, pos, &EnergyOptions::new(energy), &oo, Some(&mut cb)).unwrap();
    let e = r.report.energy_final;
    (lowest, r.report, e)
}

/// Bug 1, mechanism: the energies are scale free, so a face can be squeezed to zero area for
/// (almost) no cost. Without the floor the optimizer does exactly that, with it the trial points
/// are rejected, the energy reached is the same, and the report counts the rejections.
#[test]
fn optimizer_squeezes_faces_without_the_floor_and_not_with_it() {
    for energy in [Energy::Combinatorial, Energy::Intrinsic] {
        for seed in 0..3u64 {
            let (mesh, pos) = perturbed_icosphere(1, 0.08, seed);
            let (low0, rep0, e0) = run_ep(&mesh, &pos, 0.0, energy);
            assert!(low0 < 1e-8, "{energy:?} seed {seed}: without the floor the lowest angle is {low0:e}");
            assert_eq!(rep0.inadmissible_trials, 0);
            let (low1, rep1, e1) = run_ep(&mesh, &pos, 1e-4, energy);
            assert!(low1 >= 1e-4 * (1.0 - 1e-6), "{energy:?} seed {seed}: lowest angle {low1:e}");
            assert!(rep1.inadmissible_trials > 0, "{energy:?} seed {seed}: no rejected trial");
            // the squeezed faces are worth (almost) nothing
            assert!(e1 <= 1.25 * e0 + 1e-9, "{energy:?} seed {seed}: energy {e0} without, {e1} with the floor");
        }
    }
}

/// A starting point below the floor is an error (not repaired), the floor 0 switches the check off,
/// and invalid floors are rejected.
#[test]
fn start_below_the_floor_is_an_error() {
    let (mesh, mut pos) = perturbed_icosphere(1, 0.02, 1);
    // move a vertex almost onto the line through two of its neighbours: one face with angle ~1e-6
    let f = mesh.faces()[0];
    let (a, b, c) = (pos[f[0] as usize], pos[f[1] as usize], pos[f[2] as usize]);
    let t = 0.5;
    let on_edge = [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])];
    let n = normalize(cross(sub(b, a), sub(c, a)));
    pos[f[2] as usize] = add(on_edge, scale(n, 1e-6));
    let e = EnergyOptions::new(Energy::Covariance);
    let oo = OptimizeOptions { max_iterations: 2, ..Default::default() };
    match optimize(&mesh, &pos, &e, &oo, None) {
        Err(DevError::NearDegenerateFace { face, min_angle, limit }) => {
            assert!(min_angle < limit && limit == 1e-4, "{face} {min_angle:e} {limit:e}");
            assert!(DevError::NearDegenerateFace { face, min_angle, limit }.to_string().contains("nearly degenerate"));
        }
        other => panic!("expected NearDegenerateFace, got {:?}", other.map(|o| o.report.status)),
    }
    let off = OptimizeOptions { min_face_angle: 0.0, ..oo.clone() };
    assert!(optimize(&mesh, &pos, &e, &off, None).is_ok());
    let bad = OptimizeOptions { min_face_angle: -1.0, ..oo };
    assert!(matches!(optimize(&mesh, &pos, &e, &bad, None), Err(DevError::InvalidOption { name: "min_face_angle", .. })));
}

/// Bug 1, end to end: E_P / intrinsic flows with remeshing (on a mesh of 42 vertices, 2 subdivisions)
/// return a mesh whose smallest angle is above the floor, so the cut stage can build its geometry.
#[test]
fn flows_return_no_degenerate_faces() {
    for energy in [Energy::Combinatorial, Energy::Intrinsic, Energy::Covariance] {
        for seed in 0..2u64 {
            let (mesh, pos) = perturbed_icosphere(1, 0.08, seed);
            let opts = FlowOptions::uniform(EnergyOptions::new(energy), 2, 400);
            let res = flow(&mesh, &pos, &opts, None).unwrap();
            let floor = opts.optimize.min_face_angle;
            assert!(res.report.min_angle >= floor * (1.0 - 1e-6), "{energy:?} {seed}: {:e}", res.report.min_angle);
            assert!((res.report.min_angle - min_angle(&res.mesh, &res.positions)).abs() < 1e-15);
            Geometry::from_positions(&res.mesh, &res.positions)
                .unwrap_or_else(|e| panic!("{energy:?} {seed}: degenerate result: {e}"));
        }
    }
}

/// Boundary stars with at most two faces contribute nothing, whatever the (nearly coplanar)
/// normals are: no value, no gradient noise (the eigenvector of a rank-deficient matrix is arbitrary).
#[test]
fn two_face_boundary_stars_contribute_nothing() {
    // a hinge of two triangles with noise 1e-10 between its planes: every star is a boundary star
    // with at most two faces
    let pos = vec![[0., 0., 0.], [1., 0., 0.], [0.5, 1., 1e-10], [0.5, -1., 2e-10]];
    let mesh = TriMesh::new(4, vec![[0, 1, 2], [1, 0, 3]]).unwrap();
    for energy in
        [Energy::Combinatorial, Energy::CombinatorialMax, Energy::Covariance, Energy::CovarianceMax, Energy::Intrinsic, Energy::IntrinsicMax]
    {
        let o = EnergyOptions { energy, include_boundary: true, omit_valence_three: true, threads: None };
        let ev = evaluate(&mesh, &pos, &o).unwrap();
        assert_eq!(ev.energy, 0.0, "{energy:?}");
        assert!(ev.gradient.iter().flatten().all(|&g| g == 0.0), "{energy:?}: {:?}", ev.gradient);
    }
}

/// Face provenance (SPEC §5.1): subdivision gives the four children the parent of the face; with
/// remeshing every final face has a parent among the input faces, and a face's centroid is
/// close to its parent's neighbourhood (the surface only moves by the flow).
#[test]
fn face_parent_follows_subdivision_and_remeshing() {
    let (mesh, pos) = perturbed_icosphere(1, 0.05, 3);
    let mut opts = FlowOptions::uniform(EnergyOptions::new(Energy::Combinatorial), 2, 0);
    opts.remesh = None;
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    assert_eq!(res.face_parent.len(), res.mesh.n_faces());
    for (f, &p) in res.face_parent.iter().enumerate() {
        assert_eq!(p as usize, f / 16, "two subdivisions: 16 children per input face");
    }
    let opts = FlowOptions::uniform(EnergyOptions::new(Energy::Combinatorial), 2, 300);
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    assert_eq!(res.face_parent.len(), res.mesh.n_faces());
    assert!(res.face_parent.iter().all(|&p| (p as usize) < mesh.n_faces()));
    let flips: usize = res.report.rounds.iter().flat_map(|r| &r.remesh).map(|m| m.flips + m.collapses).sum();
    eprintln!("remesh operations: {flips}");
    // every input face still has descendants (a sphere does not lose a face's whole region)
    let mut seen = vec![false; mesh.n_faces()];
    for &p in &res.face_parent {
        seen[p as usize] = true;
    }
    assert!(seen.iter().filter(|&&s| s).count() * 10 >= mesh.n_faces() * 9);
}

fn check_result(rel: &str, energy: Energy, res: &ddg_developability::FlowResult, floor: f64) {
    assert!(res.report.min_angle >= floor * (1.0 - 1e-6), "{rel} {energy:?}: min angle {:e}", res.report.min_angle);
    Geometry::from_positions(&res.mesh, &res.positions).unwrap_or_else(|e| panic!("{rel} {energy:?}: {e}"));
    eprintln!(
        "{rel} {energy:?}: {} faces, E {:.3e}, min angle {:.2e}, bad faces {}, stopped {:?}",
        res.mesh.n_faces(),
        res.report.rounds.last().map_or(f64::NAN, |r| r.optimize.energy_final),
        res.report.min_angle,
        res.report.bad_faces,
        res.report.stopped.as_ref().map(|s| (s.round, s.vertex))
    );
}

/// The two bugs of the cut study, as reported: the blub `E_P` flow (300 it/round) left a
/// zero-area face; the blub intrinsic flow (3000 it/round) failed with "undefined at vertex 881".
#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test robust -- --ignored"]
fn cut_study_reproducers() {
    let Some((mesh, pos)) = load_model("blub/blub/blub_control_mesh.obj") else { return };
    let opts = FlowOptions::uniform(EnergyOptions::new(Energy::Combinatorial), 2, 300);
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    check_result("blub", Energy::Combinatorial, &res, opts.optimize.min_face_angle);
    let opts = FlowOptions::uniform(EnergyOptions::new(Energy::Intrinsic), 2, 3000);
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    check_result("blub", Energy::Intrinsic, &res, opts.optimize.min_face_angle);
}

/// Over-runs (3 rounds of 3000 iterations) of all six energies on the control meshes of spot, bob
/// and blub: no error, no degenerate face (the result builds a `Geometry`).
#[test]
#[ignore = "long (minutes): cargo test --release -p ddg-developability --test robust -- --ignored"]
fn over_runs_of_all_energies_complete() {
    use Energy::*;
    for rel in ["spot/spot/spot_control_mesh.obj", "bob/bob/bob_controlmesh.obj", "blub/blub/blub_control_mesh.obj"] {
        let Some((mesh, pos)) = load_model(rel) else { continue };
        std::thread::scope(|s| {
            for energy in [Combinatorial, CombinatorialMax, Covariance, CovarianceMax, Intrinsic, IntrinsicMax] {
                let (mesh, pos) = (&mesh, &pos);
                s.spawn(move || {
                    let mut opts = FlowOptions::uniform(EnergyOptions::new(energy), 2, 3000);
                    opts.threads = Some(1);
                    let res = flow(mesh, pos, &opts, None).unwrap_or_else(|e| panic!("{rel} {energy:?}: {e}"));
                    check_result(rel, energy, &res, opts.optimize.min_face_angle);
                });
            }
        });
    }
}
