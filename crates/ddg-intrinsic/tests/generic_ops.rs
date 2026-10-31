//! The representation-generic operators equal the `IntrinsicTriangulation` methods bit for bit,
//! on fresh, flipped and refined triangulations.

mod common;

use common::standard_shapes;
use ddg_core::generators::SplitMix64;
use ddg_intrinsic::{DelaunayOptions, IntrinsicTriangulation, RefineOptions, generic_ops};

fn same(tri: &IntrinsicTriangulation, what: &str) {
    let a = tri.cotan_laplacian();
    let b = generic_ops::cotan_laplacian(tri);
    assert_eq!(a.col_ptr(), b.col_ptr(), "{what}");
    assert_eq!(a.row_idx(), b.row_idx(), "{what}");
    let bits = |x: &[f64]| x.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(a.values()), bits(b.values()), "{what}");
    assert_eq!(
        bits(&tri.lumped_mass()),
        bits(&generic_ops::lumped_mass(tri)),
        "{what}"
    );
    assert_eq!(
        bits(&tri.lumped_mass_normalized()),
        bits(&generic_ops::lumped_mass_normalized(tri)),
        "{what}"
    );
}

#[test]
fn generic_operators_are_bitwise_the_methods() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        same(&tri, &format!("{} fresh", shape.name));
        let mut rng = SplitMix64::new(5 + s as u64);
        for _ in 0..tri.n_edges() {
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            let _ = tri.flip_edge(e);
        }
        same(&tri, &format!("{} flipped", shape.name));
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        tri.refine(&RefineOptions::default()).unwrap();
        same(&tri, &format!("{} refined", shape.name));
    }
}
