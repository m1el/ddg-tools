//! The index sum of `singularities` when a transported field difference is exactly `π`
//! (id `app-vector-singularity-index`; lean/proofs `TriMesh.sum_faceIndex`).
//!
//! With `ω ∈ (−π, π]` both halfedges of an edge get `ω = +π` when `u_j = −r_ij u_i` exactly
//! (in exact arithmetic), so the twin terms do not cancel and the sum is `nχ + #(such edges)`.
//! In floating point the twin's `ω` is `±π` depending on the sign of a rounding error.

mod common;

use common::*;
use ddg_intrinsic_apps::vector::*;

/// For every halfedge `h = (i → j)` of an icosphere, a generic unit field except
/// `u_j = −r_h u_i`, which puts `ω_h` at exactly `π` (the quotient is `(−|r|², +0)`).
/// Counts the halfedges for which the reported index sum differs from `nχ = 2`.
#[test]
#[ignore = "documents the ω = π tie of `singularities`; lean/proofs TriMesh.sum_faceIndex: index_sum = nχ + #ties"]
fn an_exact_pi_tie_changes_the_index_sum() {
    let shape = icosphere(1);
    let tri = shape.triangulation();
    let n_sym = 1u32;
    let chi = tri.euler_characteristic();
    let mut wrong = Vec::new();
    {
        let values: Vec<C64> = (0..tri.n_vertices())
            .map(|v| C64::from_angle(0.7 + 1.3 * v as f64 + 0.11 * (v * v) as f64))
            .collect();
        let s = singularities(&tri, &TangentField { n_sym, values }).unwrap();
        assert_eq!(s.index_sum, n_sym as i64 * chi, "the generic field itself");
    }
    for h in 0..tri.n_halfedges() {
        let (i, j) = (tri.complex().from_vertex(h), tri.complex().to_vertex(h));
        // A generic field (no near-ties of its own), then the tie on `h`.
        let mut values: Vec<C64> = (0..tri.n_vertices())
            .map(|v| C64::from_angle(0.7 + 1.3 * v as f64 + 0.11 * (v * v) as f64))
            .collect();
        let r = C64::from_angle(n_sym as f64 * transport_angle(&tri, h));
        values[j] = -(r * values[i]);
        let field = TangentField { n_sym, values };
        let s = singularities(&tri, &field).unwrap();
        assert!(s.max_rounding < 1e-9, "h {h}: {}", s.max_rounding);
        if s.index_sum != n_sym as i64 * chi {
            wrong.push((h, s.index_sum));
        }
    }
    eprintln!("{} of {} halfedges: index sum ≠ nχ = {}: {:?}", wrong.len(), tri.n_halfedges(), n_sym as i64 * chi, &wrong[..wrong.len().min(10)]);
    assert!(wrong.is_empty(), "{} halfedges with a wrong index sum", wrong.len());
}
