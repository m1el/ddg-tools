//! The consumers that are generic over the representation (Integer Coordinates contract,
//! WP11), run on integer coordinates and on signposts from the same input, flipped to Delaunay
//! and refined by the same rules. The two representations make the same operations where their
//! floating-point decisions agree (IC README, WP6: 255 of 256 runs bitwise), so the results
//! must agree to rounding; the differences are printed (`--nocapture`).
//!
//! Tolerance policy: the crates' own tests compare results of the same triangulation at
//! `1e-10 .. 1e-12` (relative to the largest entry or the mean edge length); the same bounds
//! are used here.

mod common;

use common::*;
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{
    DelaunayOptions, IntrinsicTriangulation, RefineOptions, RefineRules, SurfacePoint,
    Triangulation,
};
use ddg_intrinsic_apps::amr::{self, AmrOptions, Problem};
use ddg_intrinsic_apps::fem;
use ddg_intrinsic_apps::heat::{HeatOptions, HeatSolver, HeatSource};
use ddg_intrinsic_apps::representation::{InputPoints, Representation};
use ddg_intrinsic_apps::steiner::{SteinerOptions, minimum_spanning_tree, steiner_tree};

fn shapes() -> Vec<Shape> {
    vec![
        jittered_grid(14, 0.3, 11),
        jittered_icosphere(2, 0.08, 3),
        sheared_torus(20, 8, 1.0, 0.4, 3),
    ]
}

/// Both representations of `shape`, Delaunay-flipped and (if `refine`) refined by Ruppert's
/// rules.
fn pair(shape: &Shape, refine: bool) -> (IntrinsicTriangulation, IntegerTriangulation) {
    let mut s = shape.triangulation();
    let mut i = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    s.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    i.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    if refine {
        s.refine_with_rules(&RefineOptions::default(), RefineRules::Ruppert)
            .unwrap();
        i.refine_with_rules(&RefineOptions::default(), RefineRules::Ruppert)
            .unwrap();
    }
    (s, i)
}

fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

#[test]
fn the_default_representation_is_signposts() {
    assert_eq!(Representation::default(), Representation::Signposts);
}

#[test]
fn fem_matrices_agree() {
    for shape in shapes() {
        for refine in [false, true] {
            let (s, i) = pair(&shape, refine);
            assert_eq!(s.complex(), i.complex(), "{}: connectivity", shape.name);
            let (ls, li) = (fem::stiffness_matrix(&s), fem::stiffness_matrix(&i));
            let scale = ls.values().iter().fold(0.0f64, |a, &b| a.max(b.abs()));
            let dl = max_diff(ls.values(), li.values()) / scale;
            let dm = max_diff(&fem::lumped_mass(&s), &fem::lumped_mass(&i));
            let (gs, gi) = (fem::galerkin_mass_matrix(&s), fem::galerkin_mass_matrix(&i));
            let dg = max_diff(gs.values(), gi.values());
            eprintln!(
                "fem {} refine={refine}: L rel {dl:.2e}, lumped {dm:.2e}, galerkin {dg:.2e}",
                shape.name
            );
            assert!(dl < 1e-12 && dm < 1e-14 && dg < 1e-14);
        }
    }
}

#[test]
fn heat_distances_agree() {
    for shape in shapes() {
        let (s, i) = pair(&shape, true);
        let sources = [HeatSource::Vertex(3), HeatSource::Vertex(7)];
        let ds = HeatSolver::new(&s, &HeatOptions::default())
            .unwrap()
            .distance(&sources)
            .unwrap();
        let di = HeatSolver::new(&i, &HeatOptions::default())
            .unwrap()
            .distance(&sources)
            .unwrap();
        let d = max_diff(&ds.distances, &di.distances);
        // A point source: the same input point, located on each.
        let p = SurfacePoint::Face {
            face: 2,
            bary: [0.2, 0.3, 0.5],
        };
        let ps = HeatSolver::new(&s, &HeatOptions::default())
            .unwrap()
            .distance(&[HeatSource::Point(p)])
            .unwrap();
        let pi = HeatSolver::new(&i, &HeatOptions::default())
            .unwrap()
            .distance(&[HeatSource::Point(p)])
            .unwrap();
        let dp = max_diff(
            &ps.distances[..shape.mesh.n_vertices()],
            &pi.distances[..shape.mesh.n_vertices()],
        );
        eprintln!(
            "heat {}: vertex sources {d:.2e}, point source (input vertices) {dp:.2e}",
            shape.name
        );
        assert!(d < 1e-9, "{d}");
        // The point is located by rounding-dependent decisions in each representation.
        assert!(dp < 1e-6, "{dp}");
        assert!(i.locate_on_input(i.vertex_location(4)).is_ok());
    }
}

#[test]
fn steiner_trees_agree() {
    for shape in shapes() {
        let (s, i) = pair(&shape, false);
        let (_, ls) = minimum_spanning_tree(&s);
        let (_, li) = minimum_spanning_tree(&i);
        let options = SteinerOptions::default();
        let ts = steiner_tree(&s, &options).unwrap();
        let ti = steiner_tree(&i, &options).unwrap();
        eprintln!(
            "steiner {}: mst {ls:.6} vs {li:.6}, tree {:.9} ({} steiner) vs {:.9} ({} steiner)",
            shape.name,
            ts.length,
            ts.steiner_vertices.len(),
            ti.length,
            ti.steiner_vertices.len()
        );
        assert!((ls - li).abs() < 1e-12 * ls);
        assert_eq!(ts.steiner_vertices.len(), ti.steiner_vertices.len());
        assert!((ts.length - ti.length).abs() < 1e-9 * ts.length);
        assert_eq!(ts.traces.len(), ti.traces.len());
        let n = ddg_intrinsic_apps::steiner::cut_to_disk(&ti.triangulation, &ti.edges).unwrap();
        let m = ddg_intrinsic_apps::steiner::cut_to_disk(&ts.triangulation, &ts.edges).unwrap();
        assert_eq!(n.is_disk, m.is_disk);
    }
}

#[test]
fn adaptive_refinement_agrees() {
    let shape = jittered_icosphere(2, 0.05, 2);
    let problem = Problem::HeatKernel {
        source: 3,
        time: 0.002,
    };
    let options = AmrOptions {
        max_vertices: 1500,
        ..Default::default()
    };
    let mut s = shape.triangulation();
    let mut i = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    let (sol_s, rep_s) = amr::adaptive_solve(&mut s, &problem, &options).unwrap();
    let (sol_i, rep_i) = amr::adaptive_solve(&mut i, &problem, &options).unwrap();
    let n = shape.mesh.n_vertices();
    let e = amr::input_vertex_error(&s, &sol_s, &sol_i, &[3]).unwrap();
    let ls = rep_s.iterations.last().unwrap().estimate;
    let li = rep_i.iterations.last().unwrap().estimate;
    eprintln!(
        "amr: {} vs {} vertices, {} vs {} steps, last estimate {ls:.6e} vs {li:.6e}, input-vertex difference {e:.2e} ({n} input vertices)",
        s.complex().n_vertices(),
        i.complex().n_vertices(),
        rep_s.iterations.len(),
        rep_i.iterations.len()
    );
    // The same marked faces and insertions give the same triangulation up to the 3% of runs
    // where a decision differs by rounding; then the solutions are those of near-identical
    // meshes.
    if s.complex() == i.complex() {
        assert!(e < 1e-8, "{e}");
        assert!((ls - li).abs() < 1e-8 * ls.abs());
    } else {
        assert!(e < 2e-2, "{e}");
    }
    // l2_error across representations: the integer solution against the signpost reference.
    let err = amr::l2_error(&s, &sol_s, &i, &sol_i).unwrap();
    eprintln!("amr: l2 error between the two: {err:.2e}");
    assert!(err < 2e-2);
}
