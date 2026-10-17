//! Tangent frames, transport and the connection matrices of an intrinsic triangulation.

use std::f64::consts::{PI, TAU};

use ddg_core::SparseError;
use ddg_intrinsic::IntrinsicTriangulation;
use ddg_intrinsic::formulas::triangle_layout;

use super::complex::{C64, ComplexTriplets};
use super::knoeppel::local_matrices;

/// The factor that turns a stored signpost of a halfedge leaving `v` into its angle in the
/// tangent frame of `v`: `1` at an interior vertex (signposts are already normalized to
/// `2π`), `Θ/π` at a boundary vertex, whose tangent angles are not normalized (SPEC N41).
#[inline]
pub(crate) fn frame_factor(tri: &IntrinsicTriangulation, v: usize) -> f64 {
    if tri.is_boundary_vertex(v) {
        tri.vertex_angle_sum(v) / PI
    } else {
        1.0
    }
}

/// The scale `s_v` from corner angles to tangent angles at `v`: `2π/Θ` inside, `1` on the
/// boundary (Knöppel et al. 2013, Eq. 11).
#[inline]
pub(crate) fn tangent_scale(tri: &IntrinsicTriangulation, v: usize) -> f64 {
    if tri.is_boundary_vertex(v) {
        1.0
    } else {
        TAU / tri.vertex_angle_sum(v)
    }
}

/// The direction of halfedge `h` at its tail, in the tangent frame of the tail.
#[inline]
pub(crate) fn tail_angle(tri: &IntrinsicTriangulation, h: usize) -> f64 {
    tri.signpost(h) * frame_factor(tri, tri.from_vertex(h))
}

/// The direction, at the head of halfedge `h`, that points back along `h`, in the tangent
/// frame of the head: the twin's direction, or for a boundary halfedge the direction of the
/// boundary edge that arrives at the head (`Θ` in the unnormalized boundary frame).
#[inline]
pub(crate) fn head_back_angle(tri: &IntrinsicTriangulation, h: usize) -> f64 {
    match tri.twin(h) {
        Some(t) => tail_angle(tri, t),
        None => tri.vertex_angle_sum(tri.to_vertex(h)),
    }
}

/// The angle `ρ = (φ_ji + π) − φ_ij` by which a tangent vector turns, in the coordinates of the
/// frames, when it is transported from the tail `i` of halfedge `h` to its head `j` along the
/// edge (Sharp, Soliman, Crane 2019, "The Vector Heat Method", Eq. 10). A vector with angle
/// `α` at `i` has the angle `α + ρ` at `j`.
///
/// [formula: app-vector-transport]
pub fn transport_angle(tri: &IntrinsicTriangulation, h: usize) -> f64 {
    head_back_angle(tri, h) + PI - tail_angle(tri, h)
}

/// The transport coefficient `r_ij = e^{i n ρ_ij}` of halfedge `h` for `n`-vectors.
#[inline]
pub(crate) fn transport(tri: &IntrinsicTriangulation, h: usize, n_sym: u32) -> C64 {
    C64::from_angle(n_sym as f64 * transport_angle(tri, h))
}

/// The curvature that the tangent frames put into face `f`: `Σ_c (s_c − 1) θ_c` over its
/// corners. Transport around the face multiplies an `n`-vector by `e^{i n Ω_f}`. The sum over
/// all faces is `Σ_v s_v Θ_v − Σ_v Θ_v`, which on a closed surface is `2π χ`.
///
/// [formula: app-vector-face-holonomy]
pub fn face_holonomy(tri: &IntrinsicTriangulation, f: usize) -> f64 {
    (0..3)
        .map(|k| {
            let h = 3 * f + k;
            (tangent_scale(tri, tri.from_vertex(h)) - 1.0) * tri.corner_angle(h)
        })
        .sum()
}

/// The finite-difference connection Laplacian (VHM Sec. 5.3): `L_ii = Σ_j w_ij`,
/// `L_ij = −w_ij r_ji`, so that `u* L u = Σ_ij w_ij |u_j − r_ij u_i|²` over the edges and
/// `(L u)_i = Σ_j w_ij (u_i − r_ji u_j)`. Edges between the same vertices add up; an edge from
/// a vertex to itself adds `2 w (1 − cos nρ)` to its diagonal.
///
/// [formula: app-vector-fd-connection-laplacian]
pub(crate) fn fd_laplacian(tri: &IntrinsicTriangulation, n_sym: u32) -> ComplexTriplets {
    let mut out = ComplexTriplets::new(tri.n_vertices());
    for v in 0..tri.n_vertices() {
        out.add(v, v, C64::default());
    }
    for h in 0..tri.n_halfedges() {
        let (i, j) = (tri.from_vertex(h), tri.to_vertex(h));
        let f = h / 3;
        let k = h % 3;
        let l = tri.face_lengths(f);
        let cot = ddg_core::geometry::corner_cotan(l[k], l[(k + 1) % 3], l[(k + 2) % 3])
            .expect("every face keeps the triangle inequalities");
        let w = 0.5 * cot;
        let r = transport(tri, h, n_sym);
        out.add(i, i, C64::new(w, 0.0));
        out.add(j, j, C64::new(w, 0.0));
        out.add(i, j, r.conj().scale(-w));
        out.add(j, i, r.scale(-w));
    }
    out
}

/// The lumped mass of the finite-difference discretization, on the diagonal.
pub(crate) fn lumped_mass(tri: &IntrinsicTriangulation) -> ComplexTriplets {
    let mut out = ComplexTriplets::new(tri.n_vertices());
    for (v, m) in tri.lumped_mass().into_iter().enumerate() {
        out.add(v, v, C64::new(m, 0.0));
    }
    out
}

/// The finite element mass and stiffness matrices of Knöppel et al. (2013, Sec. 6.1.1).
pub(crate) fn fe_matrices(
    tri: &IntrinsicTriangulation,
    n_sym: u32,
) -> (ComplexTriplets, ComplexTriplets) {
    let n = tri.n_vertices();
    let mut mass = ComplexTriplets::new(n);
    let mut stiff = ComplexTriplets::new(n);
    for f in 0..tri.n_faces() {
        let p = triangle_layout(tri.face_lengths(f)).expect("faces keep the triangle inequalities");
        let area = tri.face_area(f);
        let omega = n_sym as f64 * face_holonomy(tri, f);
        let r = [0, 1, 2].map(|k| transport(tri, 3 * f + k, n_sym));
        let (m, a) = local_matrices(p, area, omega, r);
        let verts = tri.faces()[f];
        for x in 0..3 {
            for y in 0..3 {
                mass.add(verts[x] as usize, verts[y] as usize, m[x][y]);
                stiff.add(verts[x] as usize, verts[y] as usize, a[x][y]);
            }
        }
    }
    (mass, stiff)
}

/// Builds the real forms of a pair of complex triplet lists.
pub(crate) fn build_pair(
    a: &ComplexTriplets,
    m: &ComplexTriplets,
) -> Result<(ddg_core::SparseMatrix, ddg_core::SparseMatrix), SparseError> {
    Ok((a.build()?, m.build()?))
}
