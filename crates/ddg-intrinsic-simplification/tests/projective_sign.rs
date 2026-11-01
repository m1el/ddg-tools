//! The sign of the exponent in Eq. 10 (SPEC S25): both are valid tracking maps; the coarse mesh
//! (lengths, connectivity) does not depend on it, only where the tracked points land.

use ddg_core::generators::icosphere;
use ddg_intrinsic_simplification::formulas::{ProjectiveSign, projective_interpolation_signed};
use ddg_intrinsic_simplification::{
    CoarseMesh, CoarseMeshOptions, CoarsePoint, SimplifyOptions, Target, simplify,
};

fn run(sign: ProjectiveSign) -> CoarseMesh {
    let g = icosphere(3).unwrap();
    let mut mesh = CoarseMesh::from_positions(
        &g.positions,
        &g.faces,
        &[],
        CoarseMeshOptions {
            projective_sign: sign,
            ..CoarseMeshOptions::default()
        },
    )
    .unwrap();
    let opts = SimplifyOptions {
        target: Target::Vertices(60),
        ..SimplifyOptions::default()
    };
    simplify(&mut mesh, &opts).unwrap();
    mesh
}

#[test]
fn the_two_signs_are_inverse_maps_of_a_face() {
    let b = [0.2, 0.5, 0.3];
    let at = [true, false, false];
    for u in [-1.5, 0.3, 2.0] {
        let p = projective_interpolation_signed(b, at, u, ProjectiveSign::Paper);
        let q = projective_interpolation_signed(p, at, u, ProjectiveSign::Cetm);
        for k in 0..3 {
            assert!((q[k] - b[k]).abs() < 1e-14, "u = {u}: {q:?} vs {b:?}");
        }
        assert!((p.iter().sum::<f64>() - 1.0).abs() < 1e-14);
    }
}

#[test]
fn both_signs_give_the_same_coarse_mesh_and_valid_different_maps() {
    let a = run(ProjectiveSign::Paper);
    let b = run(ProjectiveSign::Cetm);
    a.certify(false).unwrap();
    b.certify(false).unwrap();
    assert_eq!(a.edge_lengths(), b.edge_lengths());
    assert_eq!(a.faces(), b.faces());
    let (oa, ob) = (a.compact().unwrap(), b.compact().unwrap());
    let mut differ = 0;
    for (pa, pb) in oa.points.iter().zip(&ob.points) {
        for p in [pa, pb] {
            if let CoarsePoint::Face { bary, .. } = p {
                assert!(bary.iter().all(|&x| x >= 0.0 && x.is_finite()));
                assert!((bary.iter().sum::<f64>() - 1.0).abs() < 1e-9);
            }
        }
        if pa != pb {
            differ += 1;
        }
    }
    assert!(differ > 0, "the two signs gave identical maps");
}
