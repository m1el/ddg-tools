//! Tangent frames, transport and the connection matrices of an intrinsic triangulation, on any
//! representation that has input-aligned frames ([`TangentFrames`]).

use std::borrow::Cow;
use std::f64::consts::{PI, TAU};

use ddg_core::SparseError;
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::formulas::triangle_layout;
use ddg_intrinsic::{IntrinsicTriangulation, Triangulation};

use super::VectorError;
use super::complex::{C64, ComplexTriplets};
use super::knoeppel::local_matrices;

/// A representation of an intrinsic triangulation whose vertices carry **input-aligned tangent
/// frames**: the direction of every halfedge at its tail, normalized as the signposts store it
/// (`[0, 2π)` at an interior vertex, `[0, π]` at a boundary vertex, SPEC N41), in the frame of
/// the input mesh at an input vertex and of the input face or edge at an inserted vertex.
///
/// * Signposts ([`IntrinsicTriangulation`]): the stored signposts, borrowed.
/// * Integer coordinates ([`IntegerTriangulation`]): computed from the roundabouts and the
///   straight segments of the edges over the input mesh
///   ([`IntegerTriangulation::frame_directions`], Integer Coordinates SPEC §4.12, I26); on the
///   same triangulation they equal the signposts to rounding.
pub trait TangentFrames: Triangulation {
    /// The direction of every halfedge at its tail, normalized (see the trait).
    ///
    /// # Errors
    /// [`VectorError::Frames`] when the representation cannot compute them (integers that are
    /// not those of a realization).
    fn tangent_frames(&self) -> Result<Cow<'_, [f64]>, VectorError>;
}

impl TangentFrames for IntrinsicTriangulation {
    fn tangent_frames(&self) -> Result<Cow<'_, [f64]>, VectorError> {
        Ok(Cow::Borrowed(self.signposts()))
    }
}

impl TangentFrames for IntegerTriangulation {
    fn tangent_frames(&self) -> Result<Cow<'_, [f64]>, VectorError> {
        self.frame_directions()
            .map(|(directions, _)| Cow::Owned(directions))
            .map_err(|e| VectorError::Frames {
                reason: e.to_string(),
            })
    }
}

/// A triangulation with its frames, read once.
pub(crate) struct Framed<'a, T: ?Sized> {
    pub(crate) tri: &'a T,
    directions: Cow<'a, [f64]>,
}

impl<'a, T: TangentFrames + ?Sized> Framed<'a, T> {
    pub(crate) fn new(tri: &'a T) -> Result<Framed<'a, T>, VectorError> {
        let directions = tri.tangent_frames()?;
        if directions.len() != tri.complex().n_halfedges() {
            return Err(VectorError::Frames {
                reason: format!(
                    "{} directions for {} halfedges",
                    directions.len(),
                    tri.complex().n_halfedges()
                ),
            });
        }
        Ok(Framed { tri, directions })
    }

    /// The normalized direction of halfedge `h` (the signpost).
    #[inline]
    pub(crate) fn signpost(&self, h: usize) -> f64 {
        self.directions[h]
    }

    /// The direction of halfedge `h` at its tail, in the tangent frame of the tail.
    #[inline]
    pub(crate) fn tail_angle(&self, h: usize) -> f64 {
        self.signpost(h) * frame_factor(self.tri, self.tri.complex().from_vertex(h))
    }

    /// The direction, at the head of halfedge `h`, that points back along `h`, in the tangent
    /// frame of the head: the twin's direction, or for a boundary halfedge the direction of the
    /// boundary edge that arrives at the head (`Θ` in the unnormalized boundary frame).
    #[inline]
    pub(crate) fn head_back_angle(&self, h: usize) -> f64 {
        let c = self.tri.complex();
        match c.twin(h) {
            Some(t) => self.tail_angle(t),
            None => self.tri.vertex_angle_sum(c.to_vertex(h)),
        }
    }

    /// [`transport_angle`] of halfedge `h`.
    ///
    /// [formula: app-vector-transport]
    #[inline]
    pub(crate) fn transport_angle(&self, h: usize) -> f64 {
        self.head_back_angle(h) + PI - self.tail_angle(h)
    }

    /// The transport coefficient `r_ij = e^{i n ρ_ij}` of halfedge `h` for `n`-vectors.
    #[inline]
    pub(crate) fn transport(&self, h: usize, n_sym: u32) -> C64 {
        C64::from_angle(n_sym as f64 * self.transport_angle(h))
    }
}

/// The factor that turns a normalized direction of a halfedge leaving `v` into its angle in the
/// tangent frame of `v`: `1` at an interior vertex (directions are already normalized to `2π`),
/// `Θ/π` at a boundary vertex, whose tangent angles are not normalized (SPEC N41).
#[inline]
pub(crate) fn frame_factor<T: Triangulation + ?Sized>(tri: &T, v: usize) -> f64 {
    if tri.complex().is_boundary_vertex(v) {
        tri.vertex_angle_sum(v) / PI
    } else {
        1.0
    }
}

/// The scale `s_v` from corner angles to tangent angles at `v`: `2π/Θ` inside, `1` on the
/// boundary (Knöppel et al. 2013, Eq. 11).
#[inline]
pub(crate) fn tangent_scale<T: Triangulation + ?Sized>(tri: &T, v: usize) -> f64 {
    if tri.complex().is_boundary_vertex(v) {
        1.0
    } else {
        TAU / tri.vertex_angle_sum(v)
    }
}

/// The angle `ρ = (φ_ji + π) − φ_ij` by which a tangent vector turns, in the coordinates of the
/// frames, when it is transported from the tail `i` of halfedge `h` to its head `j` along the
/// edge (Sharp, Soliman, Crane 2019, "The Vector Heat Method", Eq. 10). A vector with angle
/// `α` at `i` has the angle `α + ρ` at `j`. On signposts; for any representation,
/// [`transport_angles`].
///
/// [formula: app-vector-transport]
pub fn transport_angle(tri: &IntrinsicTriangulation, h: usize) -> f64 {
    Framed {
        tri,
        directions: Cow::Borrowed(tri.signposts()),
    }
    .transport_angle(h)
}

/// [`transport_angle`] of every halfedge, on any representation with input-aligned frames.
///
/// # Errors
/// [`VectorError::Frames`] (see [`TangentFrames`]).
pub fn transport_angles<T: TangentFrames + ?Sized>(tri: &T) -> Result<Vec<f64>, VectorError> {
    let framed = Framed::new(tri)?;
    Ok((0..tri.complex().n_halfedges())
        .map(|h| framed.transport_angle(h))
        .collect())
}

/// The curvature that the tangent frames put into face `f`: `Σ_c (s_c − 1) θ_c` over its
/// corners. Transport around the face multiplies an `n`-vector by `e^{i n Ω_f}`. The sum over
/// all faces is `Σ_v s_v Θ_v − Σ_v Θ_v`, which on a closed surface is `2π χ`. It depends on the
/// lengths alone, so it is the same in every representation.
///
/// [formula: app-vector-face-holonomy]
pub fn face_holonomy<T: Triangulation + ?Sized>(tri: &T, f: usize) -> f64 {
    let c = tri.complex();
    (0..3)
        .map(|k| {
            let h = 3 * f + k;
            (tangent_scale(tri, c.from_vertex(h)) - 1.0) * tri.corner_angle(h)
        })
        .sum()
}

/// The finite-difference connection Laplacian (VHM Sec. 5.3): `L_ii = Σ_j w_ij`,
/// `L_ij = −w_ij r_ji`, so that `u* L u = Σ_ij w_ij |u_j − r_ij u_i|²` over the edges and
/// `(L u)_i = Σ_j w_ij (u_i − r_ji u_j)`. Edges between the same vertices add up; an edge from
/// a vertex to itself adds `2 w (1 − cos nρ)` to its diagonal.
///
/// [formula: app-vector-fd-connection-laplacian]
pub(crate) fn fd_laplacian<T: TangentFrames + ?Sized>(
    framed: &Framed<'_, T>,
    n_sym: u32,
) -> ComplexTriplets {
    let tri = framed.tri;
    let c = tri.complex();
    let mut out = ComplexTriplets::new(c.n_vertices());
    for v in 0..c.n_vertices() {
        out.add(v, v, C64::default());
    }
    for h in 0..c.n_halfedges() {
        let (i, j) = (c.from_vertex(h), c.to_vertex(h));
        let f = h / 3;
        let k = h % 3;
        let l = tri.face_lengths(f);
        let cot = ddg_core::geometry::corner_cotan(l[k], l[(k + 1) % 3], l[(k + 2) % 3])
            .expect("every face keeps the triangle inequalities");
        let w = 0.5 * cot;
        let r = framed.transport(h, n_sym);
        out.add(i, i, C64::new(w, 0.0));
        out.add(j, j, C64::new(w, 0.0));
        out.add(i, j, r.conj().scale(-w));
        out.add(j, i, r.scale(-w));
    }
    out
}

/// The lumped mass of the finite-difference discretization, on the diagonal.
pub(crate) fn lumped_mass<T: Triangulation + ?Sized>(tri: &T) -> ComplexTriplets {
    let mut out = ComplexTriplets::new(tri.complex().n_vertices());
    for (v, m) in ddg_intrinsic::generic_ops::lumped_mass(tri)
        .into_iter()
        .enumerate()
    {
        out.add(v, v, C64::new(m, 0.0));
    }
    out
}

/// The finite element mass and stiffness matrices of Knöppel et al. (2013, Sec. 6.1.1).
pub(crate) fn fe_matrices<T: TangentFrames + ?Sized>(
    framed: &Framed<'_, T>,
    n_sym: u32,
) -> (ComplexTriplets, ComplexTriplets) {
    let tri = framed.tri;
    let c = tri.complex();
    let n = c.n_vertices();
    let mut mass = ComplexTriplets::new(n);
    let mut stiff = ComplexTriplets::new(n);
    for f in 0..c.n_faces() {
        let p = triangle_layout(tri.face_lengths(f)).expect("faces keep the triangle inequalities");
        let area = tri.face_area(f);
        let omega = n_sym as f64 * face_holonomy(tri, f);
        let r = [0, 1, 2].map(|k| framed.transport(3 * f + k, n_sym));
        let (m, a) = local_matrices(p, area, omega, r);
        let verts = c.faces()[f];
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
