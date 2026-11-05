//! The index sum of `singularities` when the transported values of an edge are exactly opposite
//! (id `app-vector-singularity-index`; lean/proofs `TriMesh.sum_round_faceIndexFixed`).
//!
//! With `ω ∈ (−π, π]` evaluated on both halfedges, a tie `u_j = −r_ij u_i` gives `+π` on both
//! (in exact arithmetic; in floating point `±π` by the sign of a rounding error), so the twin
//! terms do not cancel and the sum is `nχ + 1`. `singularities` gives the halfedge with the
//! higher index `−π` at such a tie, so the sum is `nχ` always.

mod common;

use common::*;
use ddg_intrinsic_apps::vector::*;

fn generic_field(n_vertices: usize) -> Vec<C64> {
    (0..n_vertices)
        .map(|v| C64::from_angle(0.7 + 1.3 * v as f64 + 0.11 * (v * v) as f64))
        .collect()
}

/// For every halfedge `h = (i → j)` of an icosphere and `n = 1, 2`: a generic unit field except
/// `u_j = −r_h u_i`, which puts `ω_h` at exactly `π` (the quotient is `(−|r|², +0)`). Before the
/// fix 25 of the 240 cases with `n = 1` reported `index_sum = 3`.
#[test]
fn an_exact_pi_tie_keeps_the_index_sum() {
    let shape = icosphere(1);
    let tri = shape.triangulation();
    let chi = tri.euler_characteristic();
    for n_sym in [1u32, 2] {
        let base = singularities(&tri, &TangentField { n_sym, values: generic_field(tri.n_vertices()) }).unwrap();
        assert_eq!(base.index_sum, n_sym as i64 * chi, "the generic field itself, n {n_sym}");
        let mut wrong = Vec::new();
        for h in 0..tri.n_halfedges() {
            let (i, j) = (tri.complex().from_vertex(h), tri.complex().to_vertex(h));
            let mut values = generic_field(tri.n_vertices());
            let r = C64::from_angle(n_sym as f64 * transport_angle(&tri, h));
            values[j] = -(r * values[i]);
            let s = singularities(&tri, &TangentField { n_sym, values }).unwrap();
            assert!(s.max_rounding < 1e-9, "n {n_sym} h {h}: {}", s.max_rounding);
            assert_eq!(s.n_skipped, 0);
            if s.index_sum != n_sym as i64 * chi {
                wrong.push((h, s.index_sum));
            }
        }
        assert!(
            wrong.is_empty(),
            "n {n_sym}: {} of {} halfedges give index_sum ≠ nχ = {}: {:?}",
            wrong.len(),
            tri.n_halfedges(),
            n_sym as i64 * chi,
            wrong
        );
    }
}
