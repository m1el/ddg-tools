//! Attribute transfer between the input mesh `T⁰` and the intrinsic triangulation `T¹`
//! (paper Sec. 4.3, Eq. 11–14; SPEC §6; work package WP7).
//!
//! Both triangulations are refined by the common subdivision `S`. A function `f` on the
//! vertices of `T¹` is the piecewise linear function `P₁ f` on `S`, and a function `g` on the
//! vertices of `T⁰` is `P₀ g` ([`IntegerTriangulation::interpolation_matrices`]). The `L²`
//! transfer of `f` to the input mesh is the `f̂` that minimizes (Eq. 13)
//!
//! ```text
//!     ‖f − f̂‖²  =  (P₁ f − P₀ f̂)ᵀ M_S (P₁ f − P₀ f̂),
//! ```
//!
//! that is, the solution of the normal equations (Eq. 14, `ic-transfer-normal-equations`)
//!
//! ```text
//!     P₀ᵀ M_S P₀ f̂ = P₀ᵀ M_S P₁ f.
//! ```
//!
//! The matrix on the left is the Galerkin mass matrix `M₀` of the input mesh, exactly: on every
//! input triangle the pieces of `S` are linear (`ic-transfer-mass-identity`; as integrals,
//! `ic-transfer-mass-integral`: a piece contributes `∫_piece λ_a λ_b`, and pieces that tile the
//! input triangle sum to its P1 mass integral; the tiling itself is a hypothesis of the proof). [`Transfer`] still forms `P₀ᵀ M_S P₀` from the pieces instead of using the
//! mesh's own `M₀`, which is the check of the identity (`tests/transfer.rs` compares the two);
//! the factorization of it is kept for repeated transfers. The reverse direction `T⁰ → T¹`
//! solves `M₁ ĝ = P₁ᵀ M_S P₀ g` with `M₁ = P₁ᵀ M_S P₁`, the Galerkin mass matrix of `T¹`.
//!
//! # `M_S`: the consistent mass matrix, unnormalized
//!
//! `M_S` is the Galerkin (consistent) P1 mass matrix of the pieces: a triangle of area `A` with
//! corners `r, s` contributes `A/12 · (1 + δ_rs)` to entry `(r, s)`. Each piece is a convex
//! polygon of one input triangle; it is cut into a fan ([`Overlay::triangles`]) and the area of a
//! fan triangle is `|det β| · Area(T)`, where `β` holds the barycentric coordinates of its
//! corners in the input triangle `T` (from the place of the overlay vertices on the input mesh)
//! and `Area(T)` comes from the lengths of `T`. This is exact planar geometry and needs no
//! positions in space.
//!
//! **Normalization.** `M_S` is *not* divided by the total area: `Σ M_S = Area(T⁰)`. The area
//! normalization of the workspace's `CLAUDE.md` (`Σ m_i = 1`) belongs to the cones problem, whose
//! `λ` should be dimensionless. A transfer solves `M₀ f̂ = …` with the same `M_S` on both sides,
//! so the scale cancels, and the `L²` norms reported here are in the units of the mesh
//! (`‖1‖² = Area`). Divide by [`Transfer::total_area`] for a relative error.
//!
//! # Copying, for the comparison of Fig. 14
//!
//! [`Transfer::copy_to_input`] keeps the values at the input vertices (they are the first
//! `|V⁰|` vertices of `T¹`; the others are dropped). [`Transfer::copy_to_intrinsic`] evaluates
//! the function of `T⁰` at the places of the vertices of `T¹` (`P₀ g` at the first `|V¹|`
//! vertices of `S`, which is `g` itself at the input vertices).
//!
//! # Formula ids
//!
//! `ic-transfer-normal-equations` (Lean `transfer_minimizer_iff`, `transfer_normal_unique`) is
//! tagged on [`Transfer::to_input`]; `ic-transfer-mass-identity` and `ic-transfer-mass-integral`
//! (Lean `integral_piece_mass`, `transfer_mass_tiling`) on [`Transfer::new`]; the area of a fan
//! triangle, `|det β| · Area(T)` (Lean `area_piece`), on `triangle_areas`.

use std::time::{Duration, Instant};

use ddg_core::SparseMatrix;
use ddg_core::sparse::{
    Factorization, Method, Parallelism, RefineOptions, Symbolic, relative_residual,
};
use ddg_intrinsic::{Overlay, SurfacePoint};

use crate::coords::IntegerTriangulation;
use crate::error::IntegerError;
use crate::subdivision::CommonSubdivision;

fn bad(reason: impl Into<String>) -> IntegerError {
    IntegerError::Transfer {
        reason: reason.into(),
    }
}

fn sparse(what: &str) -> impl Fn(ddg_core::sparse::SparseError) -> IntegerError + '_ {
    move |e| bad(format!("{what}: {e:?}"))
}

/// Which way a transfer goes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// `T¹ → T⁰` (Eq. 13–14): `M₀ f̂ = P₀ᵀ M_S P₁ f`.
    ToInput,
    /// `T⁰ → T¹`: `M₁ ĝ = P₁ᵀ M_S P₀ g`.
    ToIntrinsic,
}

/// What building a [`Transfer`] did, with the wall-clock times of its steps (single thread).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransferReport {
    /// Vertices of `S`.
    pub n_s_vertices: usize,
    /// Pieces of `S` (convex polygons).
    pub n_pieces: usize,
    /// Triangles of the fan triangulation of the pieces.
    pub n_triangles: usize,
    /// `Σ M_S`, the total area of the surface.
    pub total_area: f64,
    /// Time to assemble `M_S` from the pieces.
    pub t_mass_s: Duration,
    /// Time to build `P₀`, `P₁` and their transposes.
    pub t_interpolation: Duration,
    /// Time to form `M₀ = P₀ᵀ M_S P₀` and `M₁ = P₁ᵀ M_S P₁`.
    pub t_congruence: Duration,
    /// Time of the two analyses and Cholesky factorizations.
    pub t_factor: Duration,
}

/// What one solve of a transfer did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolveReport {
    /// The relative residual `‖A x − b‖ / ‖b‖` of the normal equations (`‖A x‖` if `b = 0`).
    pub residual: f64,
    /// Steps of iterative refinement that were applied (0 when the plain solve was accurate).
    pub refinement_steps: usize,
    /// The `L²` distance `‖P_from f − P_to x‖_{M_S}` of the result, the square root of the energy
    /// of Eq. 13, in the units of the mesh.
    pub l2_distance: f64,
}

/// Target accuracy of the normal equations; one or two refinement steps reach it.
const SOLVE_TOLERANCE: f64 = 1e-14;

/// The transfer operators of one state: `M_S`, `P₀`, `P₁`, `M₀`, `M₁` and the factorizations of
/// the last two.
#[derive(Debug)]
pub struct Transfer {
    fingerprint: u64,
    n0: usize,
    n1: usize,
    overlay: Overlay,
    ms: SparseMatrix,
    p0: SparseMatrix,
    p1: SparseMatrix,
    p0t: SparseMatrix,
    p1t: SparseMatrix,
    m0: SparseMatrix,
    m1: SparseMatrix,
    fact0: Factorization,
    fact1: Factorization,
    report: TransferReport,
}

/// The barycentric coordinates in input face `g` of the overlay vertex whose place on the input
/// mesh is `p`.
fn bary_in_face(
    mesh: &ddg_core::TriMesh,
    g: usize,
    p: &SurfacePoint,
) -> Result<[f64; 3], IntegerError> {
    let face = mesh.faces()[g];
    let slot = |v: usize| face.iter().position(|&w| w as usize == v);
    let mut b = [0.0; 3];
    match *p {
        SurfacePoint::Vertex(v) => {
            let k = slot(v as usize)
                .ok_or_else(|| bad(format!("input vertex {v} is not a corner of face {g}")))?;
            b[k] = 1.0;
        }
        SurfacePoint::Edge { edge, t } => {
            if edge as usize >= mesh.n_edges() {
                return Err(IntegerError::IndexOutOfRange {
                    what: "input edge",
                    index: edge as usize,
                    bound: mesh.n_edges(),
                });
            }
            let [a, c] = mesh.edge_vertices(edge as usize);
            let (ka, kc) = (slot(a), slot(c));
            match (ka, kc) {
                (Some(ka), Some(kc)) => {
                    b[ka] = 1.0 - t;
                    b[kc] = t;
                }
                _ => return Err(bad(format!("input edge {edge} is not an edge of face {g}"))),
            }
        }
        SurfacePoint::Face { face: f, bary } => {
            if f as usize != g {
                return Err(bad(format!("a vertex lies in input face {f}, not in {g}")));
            }
            b = bary;
        }
    }
    Ok(b)
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// The areas of the fan triangles of the pieces, from the barycentric coordinates of their
/// corners in the input face and the area of that face.
///
/// [formula: ic-transfer-mass-integral]
fn triangle_areas(
    tri: &IntegerTriangulation,
    overlay: &Overlay,
    triangles: &[[u32; 3]],
    parents: &[u32],
) -> Result<Vec<f64>, IntegerError> {
    let mesh = tri.input().mesh();
    let face_areas = tri.input().geometry().face_areas();
    let mut areas = Vec::with_capacity(triangles.len());
    for (t, parent) in triangles.iter().zip(parents) {
        let g = overlay.face_input[*parent as usize] as usize;
        if g >= mesh.n_faces() {
            return Err(IntegerError::IndexOutOfRange {
                what: "input face",
                index: g,
                bound: mesh.n_faces(),
            });
        }
        let mut beta = [[0.0; 3]; 3];
        for (row, &s) in beta.iter_mut().zip(t) {
            let vertex = overlay
                .vertices
                .get(s as usize)
                .ok_or(IntegerError::IndexOutOfRange {
                    what: "overlay vertex",
                    index: s as usize,
                    bound: overlay.vertices.len(),
                })?;
            *row = bary_in_face(mesh, g, &vertex.on_input)?;
        }
        areas.push(det3(&beta).abs() * face_areas[g]);
    }
    Ok(areas)
}

/// `Aᵀ M B` for sparse `A` (`n × a`), `M` (`n × n`), `B` (`n × b`); entries are summed per
/// position. Rows of `A` and `B` have at most three entries, so the cost is `9 nnz(M)`.
fn congruence(
    a: &SparseMatrix,
    m: &SparseMatrix,
    b: &SparseMatrix,
) -> Result<SparseMatrix, IntegerError> {
    let rows = |p: &SparseMatrix| {
        let mut r: Vec<Vec<(u32, f64)>> = vec![Vec::new(); p.nrows()];
        for (i, j, v) in p.triplets() {
            r[i as usize].push((j, v));
        }
        r
    };
    let (ra, rb) = (rows(a), rows(b));
    let mut triplets = Vec::with_capacity(9 * m.nnz());
    for (i, j, v) in m.triplets() {
        for &(p, x) in &ra[i as usize] {
            for &(q, y) in &rb[j as usize] {
                triplets.push((p, q, v * x * y));
            }
        }
    }
    SparseMatrix::from_triplets(a.ncols(), b.ncols(), &triplets).map_err(sparse("congruence"))
}

fn check_values(what: &str, x: &[f64], expected: usize) -> Result<(), IntegerError> {
    if x.len() != expected {
        return Err(bad(format!(
            "{what} has {} values, expected {expected}",
            x.len()
        )));
    }
    if let Some(i) = x.iter().position(|v| !v.is_finite()) {
        return Err(bad(format!("{what}[{i}] = {} is not finite", x[i])));
    }
    Ok(())
}

impl Transfer {
    /// Builds the operators of `tri` from its common subdivision: `M_S`, `P₀`, `P₁`,
    /// `M₀ = P₀ᵀ M_S P₀`, `M₁ = P₁ᵀ M_S P₁`, and their Cholesky factorizations.
    ///
    /// `M₀` is formed from the pieces, not from the input mesh; that it equals the input mesh's
    /// Galerkin mass matrix is `ic-transfer-mass-identity`, checked by the tests.
    /// [formula: ic-transfer-mass-identity]
    /// [formula: ic-transfer-mass-integral]
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] (`InvalidOptions`) if `subdivision` is not the subdivision
    /// of `tri`'s current state (its fingerprint differs, SPEC I19);
    /// [`IntegerError::Transfer`] if the pieces are inconsistent with the input mesh or a mass
    /// matrix is not positive definite (an input vertex in no face, a zero-area face).
    pub fn new(
        tri: &IntegerTriangulation,
        subdivision: &CommonSubdivision,
    ) -> Result<Transfer, IntegerError> {
        let t = Instant::now();
        let (p0, p1) = tri.interpolation_matrices(subdivision)?; // checks the fingerprint
        let t_interpolation_a = t.elapsed();

        let t = Instant::now();
        let overlay = &subdivision.overlay;
        let (triangles, parents) = overlay.triangles();
        let areas = triangle_areas(tri, overlay, &triangles, &parents)?;
        let mut trip = Vec::with_capacity(9 * triangles.len());
        let mut total_area = 0.0;
        for (tr, &a) in triangles.iter().zip(&areas) {
            total_area += a;
            for (r, &i) in tr.iter().enumerate() {
                for (s, &j) in tr.iter().enumerate() {
                    trip.push((i, j, a / 12.0 * if r == s { 2.0 } else { 1.0 }));
                }
            }
        }
        let ns = overlay.vertices.len();
        let ms = SparseMatrix::from_triplets(ns, ns, &trip).map_err(sparse("M_S"))?;
        let t_mass_s = t.elapsed();

        let t = Instant::now();
        let (p0t, p1t) = (p0.transpose(), p1.transpose());
        let t_interpolation = t_interpolation_a + t.elapsed();

        let t = Instant::now();
        let m0 = congruence(&p0, &ms, &p0)?;
        let m1 = congruence(&p1, &ms, &p1)?;
        let t_congruence = t.elapsed();

        let t = Instant::now();
        let factor = |m: &SparseMatrix, name: &str| -> Result<Factorization, IntegerError> {
            let symbolic = Symbolic::analyze(m, Method::Llt).map_err(sparse(name))?;
            symbolic
                .factor(m, Parallelism::Sequential)
                .map_err(sparse(name))
        };
        let fact0 = factor(&m0, "factorization of M0 = P0' M_S P0")?;
        let fact1 = factor(&m1, "factorization of M1 = P1' M_S P1")?;
        let t_factor = t.elapsed();

        Ok(Transfer {
            fingerprint: subdivision.fingerprint,
            n0: p0.ncols(),
            n1: p1.ncols(),
            overlay: overlay.clone(),
            ms,
            p0,
            p1,
            p0t,
            p1t,
            m0,
            m1,
            fact0,
            fact1,
            report: TransferReport {
                n_s_vertices: ns,
                n_pieces: overlay.n_faces(),
                n_triangles: triangles.len(),
                total_area,
                t_mass_s,
                t_interpolation,
                t_congruence,
                t_factor,
            },
        })
    }

    /// Builds the common subdivision of `tri` and the operators of it.
    ///
    /// # Errors
    /// Those of [`IntegerTriangulation::common_subdivision`] and [`Transfer::new`].
    pub fn from_state(tri: &IntegerTriangulation) -> Result<Transfer, IntegerError> {
        Transfer::new(tri, &tri.common_subdivision()?)
    }

    /// The fingerprint of the state the operators were made from.
    pub fn fingerprint(&self) -> u64 {
        self.fingerprint
    }

    /// Whether `tri` is still the state these operators were made from (SPEC I19). After a flip,
    /// split, insertion or removal they are stale and [`Transfer::new`] must be called again.
    pub fn is_current_for(&self, tri: &IntegerTriangulation) -> bool {
        tri.fingerprint() == self.fingerprint
    }

    /// What building did (sizes, step times).
    pub fn report(&self) -> &TransferReport {
        &self.report
    }

    /// `Σ M_S`: the total area of the surface.
    pub fn total_area(&self) -> f64 {
        self.report.total_area
    }

    /// Number of vertices of `T⁰`.
    pub fn n_input_vertices(&self) -> usize {
        self.n0
    }

    /// Number of vertices of `T¹`.
    pub fn n_intrinsic_vertices(&self) -> usize {
        self.n1
    }

    /// The common subdivision the operators were made from.
    pub fn overlay(&self) -> &Overlay {
        &self.overlay
    }

    /// `M_S`, the unnormalized Galerkin mass matrix of the common subdivision.
    pub fn mass_subdivision(&self) -> &SparseMatrix {
        &self.ms
    }

    /// `P₀`, `|V^S| × |V⁰|`.
    pub fn p0(&self) -> &SparseMatrix {
        &self.p0
    }

    /// `P₁`, `|V^S| × |V¹|`.
    pub fn p1(&self) -> &SparseMatrix {
        &self.p1
    }

    /// `M₀ = P₀ᵀ M_S P₀`: equal to the Galerkin mass matrix of the input mesh.
    pub fn mass_input(&self) -> &SparseMatrix {
        &self.m0
    }

    /// `M₁ = P₁ᵀ M_S P₁`: the Galerkin mass matrix of `T¹`.
    pub fn mass_intrinsic(&self) -> &SparseMatrix {
        &self.m1
    }

    /// Transfers `f` on `V¹` to `V⁰` by `L²` projection (Eq. 13–14).
    /// [formula: ic-transfer-normal-equations]
    ///
    /// # Errors
    /// [`IntegerError::Transfer`] if `f` does not have `|V¹|` finite values.
    pub fn to_input(&self, f: &[f64]) -> Result<Vec<f64>, IntegerError> {
        Ok(self.to_input_report(f)?.0)
    }

    /// [`Transfer::to_input`] with the residual of the normal equations and the `L²` distance.
    ///
    /// # Errors
    /// As [`Transfer::to_input`].
    pub fn to_input_report(&self, f: &[f64]) -> Result<(Vec<f64>, SolveReport), IntegerError> {
        self.solve(Direction::ToInput, f)
    }

    /// Transfers `g` on `V⁰` to `V¹` by `L²` projection: `M₁ ĝ = P₁ᵀ M_S P₀ g`.
    ///
    /// # Errors
    /// [`IntegerError::Transfer`] if `g` does not have `|V⁰|` finite values.
    pub fn to_intrinsic(&self, g: &[f64]) -> Result<Vec<f64>, IntegerError> {
        Ok(self.to_intrinsic_report(g)?.0)
    }

    /// [`Transfer::to_intrinsic`] with the residual and the `L²` distance.
    ///
    /// # Errors
    /// As [`Transfer::to_intrinsic`].
    pub fn to_intrinsic_report(&self, g: &[f64]) -> Result<(Vec<f64>, SolveReport), IntegerError> {
        self.solve(Direction::ToIntrinsic, g)
    }

    fn solve(&self, dir: Direction, x: &[f64]) -> Result<(Vec<f64>, SolveReport), IntegerError> {
        let (from, to, _from_t, to_t, matrix, fact, n_from, n_to, name) = match dir {
            Direction::ToInput => (
                &self.p1,
                &self.p0,
                &self.p1t,
                &self.p0t,
                &self.m0,
                &self.fact0,
                self.n1,
                self.n0,
                "f on V1",
            ),
            Direction::ToIntrinsic => (
                &self.p0,
                &self.p1,
                &self.p0t,
                &self.p1t,
                &self.m1,
                &self.fact1,
                self.n0,
                self.n1,
                "g on V0",
            ),
        };
        check_values(name, x, n_from)?;
        let ns = self.ms.nrows();
        let mut on_s = vec![0.0; ns];
        let mut mass_on_s = vec![0.0; ns];
        from.mul_vec(x, &mut on_s).map_err(sparse("P x"))?;
        self.ms
            .mul_vec(&on_s, &mut mass_on_s)
            .map_err(sparse("M_S P x"))?;
        let mut rhs = vec![0.0; n_to];
        to_t.mul_vec(&mass_on_s, &mut rhs)
            .map_err(sparse("P' M_S P x"))?;
        let mut out = vec![0.0; n_to];
        let refine = fact
            .solve_refined(
                matrix,
                &rhs,
                &mut out,
                RefineOptions {
                    max_steps: 3,
                    tolerance: SOLVE_TOLERANCE,
                },
                Parallelism::Sequential,
            )
            .map_err(sparse("solve"))?;
        // Eq. 13 at the minimizer: sqrt of (P_from x - P_to out)' M_S (...).
        let mut diff = vec![0.0; ns];
        to.mul_vec(&out, &mut diff).map_err(sparse("P out"))?;
        for (d, a) in diff.iter_mut().zip(&on_s) {
            *d = a - *d;
        }
        let l2_distance = self.mass_norm(&diff)?;
        let residual = relative_residual(matrix, &out, &rhs).map_err(sparse("residual"))?;
        Ok((
            out,
            SolveReport {
                residual,
                refinement_steps: refine.steps,
                l2_distance,
            },
        ))
    }

    /// Copies `f` on `V¹` to `V⁰`: the values at the input vertices, which are the first `|V⁰|`
    /// vertices of `T¹`.
    ///
    /// # Errors
    /// [`IntegerError::Transfer`] if `f` does not have `|V¹|` finite values.
    pub fn copy_to_input(&self, f: &[f64]) -> Result<Vec<f64>, IntegerError> {
        check_values("f on V1", f, self.n1)?;
        if self.n0 > self.n1 {
            return Err(bad("T1 has fewer vertices than T0"));
        }
        Ok(f[..self.n0].to_vec())
    }

    /// Copies `g` on `V⁰` to `V¹`: the function of `T⁰` evaluated at the place of every vertex
    /// of `T¹` (its own value at an input vertex, the linear interpolation along an input edge or
    /// in an input face at an inserted vertex).
    ///
    /// # Errors
    /// [`IntegerError::Transfer`] if `g` does not have `|V⁰|` finite values.
    pub fn copy_to_intrinsic(&self, g: &[f64]) -> Result<Vec<f64>, IntegerError> {
        check_values("g on V0", g, self.n0)?;
        let mut on_s = vec![0.0; self.ms.nrows()];
        self.p0.mul_vec(g, &mut on_s).map_err(sparse("P0 g"))?;
        on_s.truncate(self.n1);
        Ok(on_s)
    }

    /// `sqrt(uᵀ M_S u)` for a function `u` given at the vertices of `S`.
    ///
    /// # Errors
    /// [`IntegerError::Transfer`] for a wrong length or a non-finite value.
    pub fn mass_norm(&self, u: &[f64]) -> Result<f64, IntegerError> {
        check_values("u on S", u, self.ms.nrows())?;
        let mut mu = vec![0.0; u.len()];
        self.ms.mul_vec(u, &mut mu).map_err(sparse("M_S u"))?;
        let q: f64 = u.iter().zip(&mu).map(|(a, b)| a * b).sum();
        Ok(q.max(0.0).sqrt())
    }

    /// The positions in space of the vertices of `S`, `P₀ x` for the positions `x` of the input
    /// vertices. Exact, since the pieces are flat.
    ///
    /// # Errors
    /// [`IntegerError::Transfer`] for a wrong length or a non-finite coordinate.
    pub fn subdivision_positions(
        &self,
        input_positions: &[[f64; 3]],
    ) -> Result<Vec<[f64; 3]>, IntegerError> {
        if input_positions.len() != self.n0 {
            return Err(bad(format!(
                "{} positions for {} input vertices",
                input_positions.len(),
                self.n0
            )));
        }
        let mut out = vec![[0.0; 3]; self.ms.nrows()];
        for c in 0..3 {
            let x: Vec<f64> = input_positions.iter().map(|p| p[c]).collect();
            check_values("positions", &x, self.n0)?;
            let mut y = vec![0.0; out.len()];
            self.p0.mul_vec(&x, &mut y).map_err(sparse("P0 x"))?;
            for (o, v) in out.iter_mut().zip(&y) {
                o[c] = *v;
            }
        }
        Ok(out)
    }

    /// The `L²` error `‖u − exact‖` over the surface of the function `u` given at the vertices of
    /// `S` (piecewise linear on the fan triangles of the pieces), against `exact(point)`, by a
    /// degree 4 quadrature (6 points) on every fan triangle of `S`; `positions` are those of
    /// [`Transfer::subdivision_positions`]. In the units of the mesh.
    ///
    /// The quadrature is exact for the polynomial part of the integrand only; `exact` should be
    /// smooth on the pieces.
    ///
    /// # Errors
    /// [`IntegerError::Transfer`] for wrong lengths or non-finite values.
    pub fn l2_error_to_function(
        &self,
        positions: &[[f64; 3]],
        u: &[f64],
        exact: impl Fn([f64; 3]) -> f64,
    ) -> Result<f64, IntegerError> {
        let ns = self.ms.nrows();
        check_values("u on S", u, ns)?;
        if positions.len() != ns {
            return Err(bad(format!(
                "{} positions for {ns} vertices of S",
                positions.len()
            )));
        }
        // Dunavant degree 4: (barycentric a, a, 1-2a) with weights (fractions of the area).
        const A1: f64 = 0.445_948_490_915_965;
        const A2: f64 = 0.091_576_213_509_771;
        const W1: f64 = 0.223_381_589_678_011;
        const W2: f64 = 0.109_951_743_655_322;
        let mut rule = Vec::with_capacity(6);
        for (a, w) in [(A1, W1), (A2, W2)] {
            let c = 1.0 - 2.0 * a;
            rule.push(([a, a, c], w));
            rule.push(([a, c, a], w));
            rule.push(([c, a, a], w));
        }
        let (triangles, _) = self.overlay.triangles();
        let mut sum = 0.0;
        for t in &triangles {
            let [p, q, r] = t.map(|i| positions[i as usize]);
            let (e1, e2) = (sub(q, p), sub(r, p));
            let area = 0.5 * norm(cross(e1, e2));
            let mut local = 0.0;
            for (b, w) in &rule {
                let x = [
                    b[0] * p[0] + b[1] * q[0] + b[2] * r[0],
                    b[0] * p[1] + b[1] * q[1] + b[2] * r[1],
                    b[0] * p[2] + b[1] * q[2] + b[2] * r[2],
                ];
                let value =
                    b[0] * u[t[0] as usize] + b[1] * u[t[1] as usize] + b[2] * u[t[2] as usize];
                let d = value - exact(x);
                local += w * d * d;
            }
            sum += area * local;
        }
        Ok(sum.sqrt())
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
