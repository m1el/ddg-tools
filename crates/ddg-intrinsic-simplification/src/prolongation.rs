//! Prolongation from the coarse to the fine mesh (paper Sec. 7.2 and 7.3; SPEC §7.2, §7.3, S18).
//!
//! * [`CoarseOutput::scalar_prolongation`] (`P`, Sec. 7.2): row `i` is the barycentric
//!   interpolation at the image of fine vertex `i`; rows sum to one.
//! * [`CoarseOutput::vector_prolongation`] (`P^vec`, Sec. 7.3): complex entries `b_j ψ_i z_j`,
//!   acting on coarse vectors in the **unnormalized** frames of the coarse vertices
//!   ([`unnormalize`], [`normalize`]) and producing fine vectors in the **normalized** frames of
//!   the fine mesh ([`FineFrames`]: those of [`IntrinsicTriangulation::new`] on the input). An
//!   `n`-direction field (Knöppel et al. 2013 encoding, `u = |u| e^{inα}`) takes the `n`-th
//!   powers of `ψ` and `z`.
//!
//! # Vector prolongation, as built
//!
//! **Coarse interpolation (Sec. 7.3.1).** In the frame of coarse face `ijk` the oriented side
//! `ij` points along `φ = 0`; the sides have the angles `β_ij = 0`, `β_jk = π − θ_j`,
//! `β_ki = (π − θ_j) + (π − θ_k)` (the reference's `fsv`). The unit vector of a halfedge in
//! the unnormalized frame at its tail is `w = e^{i φ/s}` (`φ` the angular coordinate, `s` the
//! scale `2π/Θ` inside, which is the paper's `(e/ℓ)^{Θ/2π}`, and `π/Θ` on the boundary, the
//! generalization to the half-plane frames of S4); `z = e^{iβ}/w` rotates a vector from the
//! frame of the corner to the frame of the face. [formula: iem-vector-interpolation]
//!
//! **Fine sampling (Sec. 7.3.2).** For a fine vertex `i` whose image `x_i` lies in coarse face
//! `fc`, and a fine neighbour `j` whose image lies in the same face, `ψ_ij = e^{iφ_ij} /
//! w̄_ij`: the fine edge's direction in the normalized frame of `i` over the unit vector from
//! `x_i` to `x_j` in the face (laid out in the plane from the coarse lengths). `ψ_i` is the
//! average of the `ψ_ij`, normalized to unit modulus. The fallbacks of S18, in order:
//!
//! 1. the neighbours in the same face (a neighbour that is a coarse vertex counts as lying in
//!    every face around it, at its corner; the reference ignores such neighbours, R29);
//! 2. neighbours in a face that is edge-adjacent to `fc`, laid out across the shared side by
//!    the lengths of the two faces (the reference's addition R29, without its sign bug R31);
//! 3. the average of the rows of the neighbours that have one (the paper's "average of known
//!    interpolated values"), each row carried to the frame of `i` by the fine transport
//!    `e^{i n (φ_ij − φ_ji − π)}`, repeated breadth-first until every row is filled. A vertex
//!    in a component without any known row, which cannot happen on a connected mesh with at
//!    least one usable vertex, gets the face's own corners with `ψ = 1` and is counted.
//!
//! A fine vertex that is a coarse vertex `v̄` has one entry, at its own column: the average over
//! the faces `f` around `v̄` (corner slot `s`) and the neighbours `j` in `f` of `ψ_ij z_s`,
//! normalized (the reference uses `1`, R29), then the same fallback 3. No row is ever empty
//! ([`VectorProlongationReport`] counts the rows of each kind).
//!
//! The input to a prolongation is *unnormalized*: `u_i = |u⁰_i| e^{i arg(u⁰_i)/s_i}` for the
//! normalized `u⁰_i` stored at the vertex (the paper writes `u⁰ e^{iΘφ/2π}`; the reference
//! skips this step, R28). The map is nonlinear, so it is applied to the field by
//! [`unnormalize`] and not folded into the matrix; the matrix is complex-linear in the
//! unnormalized vectors.

#![allow(clippy::needless_range_loop)]

use std::f64::consts::{PI, TAU};
use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, Sub};

use ddg_core::{Geometry, SparseMatrix, TriMesh};
use ddg_intrinsic::IntrinsicTriangulation;
use ddg_intrinsic::formulas::triangle_layout;

use crate::output::{CoarseOutput, CoarsePoint};

/// A complex number (a vector of a tangent plane).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct C64 {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl C64 {
    /// `re + i im`.
    #[inline]
    pub const fn new(re: f64, im: f64) -> C64 {
        C64 { re, im }
    }
    /// `e^{iθ}`.
    #[inline]
    pub fn from_angle(theta: f64) -> C64 {
        let (s, c) = theta.sin_cos();
        C64 { re: c, im: s }
    }
    /// The modulus.
    #[inline]
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    /// The argument in `(-π, π]`.
    #[inline]
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
    /// The conjugate.
    #[inline]
    pub fn conj(self) -> C64 {
        C64 {
            re: self.re,
            im: -self.im,
        }
    }
    /// Multiplication by a real.
    #[inline]
    pub fn scale(self, s: f64) -> C64 {
        C64 {
            re: self.re * s,
            im: self.im * s,
        }
    }
}

impl Add for C64 {
    type Output = C64;
    #[inline]
    fn add(self, o: C64) -> C64 {
        C64::new(self.re + o.re, self.im + o.im)
    }
}
impl Sub for C64 {
    type Output = C64;
    #[inline]
    fn sub(self, o: C64) -> C64 {
        C64::new(self.re - o.re, self.im - o.im)
    }
}
impl Mul for C64 {
    type Output = C64;
    #[inline]
    fn mul(self, o: C64) -> C64 {
        C64::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}
impl Div for C64 {
    type Output = C64;
    #[inline]
    fn div(self, o: C64) -> C64 {
        let n = o.re * o.re + o.im * o.im;
        let p = self * o.conj();
        C64::new(p.re / n, p.im / n)
    }
}
impl AddAssign for C64 {
    #[inline]
    fn add_assign(&mut self, o: C64) {
        self.re += o.re;
        self.im += o.im;
    }
}

/// The scale `s` from corner angles to normalized angles at a vertex of angle sum `theta`:
/// `2π/Θ` inside, `π/Θ` on the boundary (SPEC S4).
#[inline]
fn scale_of(theta: f64, boundary: bool) -> f64 {
    (if boundary { PI } else { TAU }) / theta
}

/// The unnormalized vector of the normalized `u0` stored at a vertex of angle sum `theta`:
/// `|u0| e^{i arg(u0)/s}` (paper Sec. 7.3: `u_i = u⁰_i e^{iΘ_i φ_i/2π}` with `φ_i = arg u⁰_i`
/// in the paper's loose notation), `arg` taken in `[0, 2π)` inside and in `(-π, π]` on the
/// boundary, where a vector in the half-plane has `arg ∈ [0, π]`. For `n`-direction fields
/// (angle `nα`) the same formula applies to the encoded value: `nα/s`.
pub fn unnormalize(u0: C64, theta: f64, boundary: bool) -> C64 {
    let s = scale_of(theta, boundary);
    let a = if boundary {
        u0.arg()
    } else {
        u0.arg().rem_euclid(TAU)
    };
    C64::from_angle(a / s).scale(u0.abs())
}

/// The inverse of [`unnormalize`]: `|u| e^{i s arg(u)}`, `arg` in `[0, 2π)` (the unnormalized
/// angles of a boundary vertex run over `[0, Θ]`, so a cone angle above `2π` is ambiguous there).
pub fn normalize(u: C64, theta: f64, boundary: bool) -> C64 {
    let s = scale_of(theta, boundary);
    // An unnormalized boundary angle runs over [0, Θ], which can pass π: no principal value.
    let a = u.arg().rem_euclid(TAU);
    C64::from_angle(a * s).scale(u.abs())
}

/// Errors of the prolongation.
#[derive(Debug, Clone, PartialEq)]
pub enum ProlongationError {
    /// The fine frames could not be built.
    FineFrames(String),
    /// The symmetry `n` of an `n`-direction field must be at least 1.
    InvalidSymmetry(u32),
    /// The fine frames and the coarse output are not of the same input mesh.
    Mismatch {
        /// What differs.
        what: &'static str,
        /// The value on the coarse output.
        coarse: usize,
        /// The value on the other argument.
        other: usize,
    },
    /// A face of the coarse mesh is not a valid triangle by its lengths.
    InvalidFace {
        /// The coarse face.
        face: usize,
        /// Its side lengths.
        lengths: [f64; 3],
    },
    /// The sparse matrix could not be built.
    Sparse(String),
}

impl fmt::Display for ProlongationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ProlongationError::FineFrames(s) => write!(f, "fine frames: {s}"),
            ProlongationError::InvalidSymmetry(n) => {
                write!(
                    f,
                    "the symmetry n of a direction field must be at least 1, got {n}"
                )
            }
            ProlongationError::Mismatch {
                what,
                coarse,
                other,
            } => write!(
                f,
                "the coarse output has {coarse} {what} but the other argument has {other}"
            ),
            ProlongationError::InvalidFace { face, lengths } => write!(
                f,
                "coarse face {face} with side lengths {lengths:?} is not a triangle"
            ),
            ProlongationError::Sparse(s) => write!(f, "sparse matrix: {s}"),
        }
    }
}

impl std::error::Error for ProlongationError {}

/// One fine neighbour `j` of fine vertex `i`, with both directions of the edge in the
/// normalized frames.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FineNeighbour {
    /// The neighbour `j`.
    pub vertex: u32,
    /// The direction of the edge in the normalized frame of `i`.
    pub phi_ij: f64,
    /// The direction of the edge in the normalized frame of `j` (`π` for the virtual
    /// direction of the incoming boundary edge).
    pub phi_ji: f64,
}

/// What the vector prolongation needs from the fine mesh: the neighbours of every vertex with
/// the directions of the edges in the normalized frames of the input mesh (the signposts of
/// [`IntrinsicTriangulation::new`], SPEC §7.3).
#[derive(Debug, Clone, PartialEq)]
pub struct FineFrames {
    ptr: Vec<usize>,
    neighbours: Vec<FineNeighbour>,
    /// Per fine vertex: the angle sum `Θ`.
    pub angle_sum: Vec<f64>,
    /// Per fine vertex: whether it is on the boundary.
    pub boundary: Vec<bool>,
}

impl FineFrames {
    /// The frames of the input mesh (loops are skipped; a boundary vertex also has the
    /// neighbour along its incoming boundary edge, at the direction `π`).
    ///
    /// # Errors
    /// [`ProlongationError::FineFrames`] if the geometry does not belong to the mesh.
    pub fn from_mesh(mesh: &TriMesh, geometry: &Geometry) -> Result<FineFrames, ProlongationError> {
        let tri = IntrinsicTriangulation::new(mesh, geometry)
            .map_err(|e| ProlongationError::FineFrames(e.to_string()))?;
        Ok(FineFrames::from_triangulation(&tri))
    }

    /// The frames of an intrinsic triangulation (normally the input mesh itself, but any
    /// triangulation of the input works when the coarse mesh is tracked from it).
    pub fn from_triangulation(tri: &IntrinsicTriangulation) -> FineFrames {
        let c = tri.complex();
        let n = c.n_vertices();
        let mut lists: Vec<Vec<FineNeighbour>> = vec![Vec::new(); n];
        for h in 0..c.n_halfedges() {
            let (i, j) = (c.from_vertex(h), c.to_vertex(h));
            if i == j {
                continue;
            }
            let here = tri.signpost(h);
            match c.twin(h) {
                Some(t) => lists[i].push(FineNeighbour {
                    vertex: j as u32,
                    phi_ij: here,
                    phi_ji: tri.signpost(t),
                }),
                None => {
                    // Boundary halfedge i -> j: j sees i along its incoming boundary edge,
                    // the virtual direction pi.
                    lists[i].push(FineNeighbour {
                        vertex: j as u32,
                        phi_ij: here,
                        phi_ji: PI,
                    });
                    lists[j].push(FineNeighbour {
                        vertex: i as u32,
                        phi_ij: PI,
                        phi_ji: here,
                    });
                }
            }
        }
        let mut ptr = vec![0usize];
        let mut neighbours = Vec::new();
        for l in lists {
            neighbours.extend(l);
            ptr.push(neighbours.len());
        }
        FineFrames {
            ptr,
            neighbours,
            angle_sum: (0..n).map(|v| tri.vertex_angle_sum(v)).collect(),
            boundary: (0..n).map(|v| c.is_boundary_vertex(v)).collect(),
        }
    }

    /// The number of fine vertices.
    pub fn n_vertices(&self) -> usize {
        self.angle_sum.len()
    }

    /// The neighbours of fine vertex `i`.
    pub fn neighbours(&self, i: usize) -> &[FineNeighbour] {
        &self.neighbours[self.ptr[i]..self.ptr[i + 1]]
    }
}

/// A sparse complex matrix, rows compressed (`P^vec`).
#[derive(Debug, Clone, PartialEq)]
pub struct ComplexSparse {
    /// The number of rows.
    pub nrows: usize,
    /// The number of columns.
    pub ncols: usize,
    /// Row `i` is `col_idx[row_ptr[i]..row_ptr[i + 1]]` with the values at the same places.
    pub row_ptr: Vec<usize>,
    /// Column indices, strictly increasing within a row.
    pub col_idx: Vec<u32>,
    /// Values.
    pub values: Vec<C64>,
}

impl ComplexSparse {
    /// The columns and values of row `i`.
    pub fn row(&self, i: usize) -> (&[u32], &[C64]) {
        let r = self.row_ptr[i]..self.row_ptr[i + 1];
        (&self.col_idx[r.clone()], &self.values[r])
    }

    /// The number of stored entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// The number of rows without an entry.
    pub fn n_empty_rows(&self) -> usize {
        (0..self.nrows)
            .filter(|&i| self.row_ptr[i] == self.row_ptr[i + 1])
            .count()
    }

    /// `y = A x`.
    ///
    /// # Errors
    /// [`ProlongationError::Mismatch`] if `x` does not have `ncols` entries.
    pub fn mul_vec(&self, x: &[C64]) -> Result<Vec<C64>, ProlongationError> {
        if x.len() != self.ncols {
            return Err(ProlongationError::Mismatch {
                what: "columns",
                coarse: self.ncols,
                other: x.len(),
            });
        }
        Ok((0..self.nrows)
            .map(|i| {
                let (c, v) = self.row(i);
                c.iter()
                    .zip(v)
                    .fold(C64::default(), |s, (&j, &a)| s + a * x[j as usize])
            })
            .collect())
    }
}

/// How the rows of the vector prolongation were obtained (SPEC S18). The counts sum to the
/// number of fine vertices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct VectorProlongationReport {
    /// Removed fine vertices sampled from neighbours in the same coarse face (paper).
    pub n_same_face: usize,
    /// Removed fine vertices sampled from neighbours in an edge-adjacent face (R29).
    pub n_adjacent_face: usize,
    /// Fine vertices that are coarse vertices, sampled from their neighbours' images.
    pub n_kept_from_images: usize,
    /// Removed fine vertices filled by the breadth-first average of neighbouring rows.
    pub n_averaged: usize,
    /// Kept fine vertices filled by the breadth-first average of neighbouring rows.
    pub n_kept_averaged: usize,
    /// Vertices with no known neighbour anywhere (own corners, `ψ = 1`).
    pub n_default: usize,
    /// Vertices whose neighbours' `ψ_ij` cancelled to zero in a tier (passed to the next).
    pub n_cancelled: usize,
    /// Neighbour images that coincide with the vertex's image (skipped).
    pub n_coincident_neighbours: usize,
    /// Rounds of the breadth-first average.
    pub averaging_rounds: usize,
}

type Row = Vec<(u32, C64)>;

/// Below this modulus the average of unit `ψ_ij` counts as cancelled.
const MIN_AVERAGE: f64 = 1e-8;

impl CoarseOutput {
    /// The scalar prolongation `P` of the fine vertices, `|V| × |V̄|` (SPEC §7.2): row `i` has
    /// a one at the column of `i` if it is a coarse vertex and the barycentric coordinates at
    /// the corners of its coarse face otherwise (summed if a corner repeats).
    ///
    /// # Errors
    /// [`ProlongationError::Sparse`] if the matrix is too large.
    pub fn scalar_prolongation(&self) -> Result<SparseMatrix, ProlongationError> {
        self.point_prolongation_rows(self.n_fine_vertices)
    }

    /// The same for every tracked point (the fine vertices, then the extra points): the
    /// interpolation of a coarse function at each point.
    ///
    /// # Errors
    /// As [`scalar_prolongation`](Self::scalar_prolongation).
    pub fn point_prolongation(&self) -> Result<SparseMatrix, ProlongationError> {
        self.point_prolongation_rows(self.points.len())
    }

    fn point_prolongation_rows(&self, n: usize) -> Result<SparseMatrix, ProlongationError> {
        let mut t: Vec<(u32, u32, f64)> = Vec::with_capacity(3 * n);
        for (i, p) in self.points.iter().take(n).enumerate() {
            match *p {
                CoarsePoint::Vertex(v) => t.push((i as u32, v, 1.0)),
                CoarsePoint::Face { face, bary } => {
                    for k in 0..3 {
                        t.push((i as u32, self.faces[face as usize][k], bary[k]));
                    }
                }
            }
        }
        SparseMatrix::from_triplets(n, self.vertices.len(), &t)
            .map_err(|e| ProlongationError::Sparse(e.to_string()))
    }

    /// The unnormalized vectors of a field given in the normalized frames of the coarse
    /// vertices (indexed by coarse id); see [`unnormalize`].
    ///
    /// # Errors
    /// [`ProlongationError::Mismatch`] if `u0` does not have one entry per coarse vertex.
    pub fn unnormalize_field(&self, u0: &[C64]) -> Result<Vec<C64>, ProlongationError> {
        self.map_field(u0, unnormalize)
    }

    /// The inverse of [`unnormalize_field`](Self::unnormalize_field).
    ///
    /// # Errors
    /// As [`unnormalize_field`](Self::unnormalize_field).
    pub fn normalize_field(&self, u: &[C64]) -> Result<Vec<C64>, ProlongationError> {
        self.map_field(u, normalize)
    }

    fn map_field(
        &self,
        u: &[C64],
        f: fn(C64, f64, bool) -> C64,
    ) -> Result<Vec<C64>, ProlongationError> {
        if u.len() != self.vertices.len() {
            return Err(ProlongationError::Mismatch {
                what: "vertices",
                coarse: self.vertices.len(),
                other: u.len(),
            });
        }
        Ok(u.iter()
            .enumerate()
            .map(|(v, &x)| f(x, self.angle_sum[v], self.boundary[v]))
            .collect())
    }

    /// The vector prolongation `P^vec` (`|V| × |V̄|`, complex) and the report of how its rows
    /// were found (module documentation, SPEC §7.3, S18). `n_sym = 1` for vectors, `n` for
    /// `n`-direction fields. Applied to the **unnormalized** coarse vectors
    /// ([`unnormalize_field`](Self::unnormalize_field)); the result is in the normalized
    /// frames of the fine mesh.
    ///
    /// [formula: iem-vector-interpolation]
    ///
    /// # Errors
    /// [`ProlongationError::InvalidSymmetry`], [`ProlongationError::Mismatch`] (the fine
    /// frames are of another mesh), [`ProlongationError::InvalidFace`].
    pub fn vector_prolongation(
        &self,
        fine: &FineFrames,
        n_sym: u32,
    ) -> Result<(ComplexSparse, VectorProlongationReport), ProlongationError> {
        if n_sym == 0 {
            return Err(ProlongationError::InvalidSymmetry(n_sym));
        }
        let nfine = self.n_fine_vertices;
        if fine.n_vertices() != nfine {
            return Err(ProlongationError::Mismatch {
                what: "fine vertices",
                coarse: nfine,
                other: fine.n_vertices(),
            });
        }
        let n = n_sym as f64;
        let nf = self.faces.len();
        let mut layouts: Vec<[[f64; 2]; 3]> = Vec::with_capacity(nf);
        for f in 0..nf {
            layouts.push(triangle_layout(self.face_lengths[f]).ok_or(
                ProlongationError::InvalidFace {
                    face: f,
                    lengths: self.face_lengths[f],
                },
            )?);
        }
        // z as angles: beta - phi / s.
        let z_angle: Vec<[f64; 3]> = (0..nf)
            .map(|f| {
                let th = self.corner_angles[f];
                let beta = [0.0, PI - th[1], PI + th[0]];
                [0, 1, 2].map(|k| {
                    let v = self.faces[f][k] as usize;
                    let s = scale_of(self.angle_sum[v], self.boundary[v]);
                    beta[k] - self.phi[3 * f + k] / s
                })
            })
            .collect();
        let (slot_ptr, slots) = self.vertex_slots();
        let mut report = VectorProlongationReport::default();
        let mut rows: Vec<Option<Row>> = vec![None; nfine];
        let mut cand: Vec<[f64; 3]> = Vec::new();

        for i in 0..nfine {
            match self.points[i] {
                CoarsePoint::Face { face, bary } => {
                    let fc = face as usize;
                    let pi = place(&layouts[fc], bary);
                    let eps = 1e-12 * self.face_lengths[fc].iter().cloned().fold(0.0, f64::max);
                    // Tier 1: neighbours in the same face.
                    let mut sum = C64::default();
                    let mut cnt = 0usize;
                    for nb in fine.neighbours(i) {
                        self.candidates(nb.vertex as usize, fc, &mut cand);
                        for b in &cand {
                            let d = place(&layouts[fc], *b);
                            accumulate(&mut sum, &mut cnt, pi, d, nb.phi_ij, eps, &mut report);
                        }
                    }
                    let mut tier = 0;
                    if cnt > 0 && sum.abs() >= MIN_AVERAGE * cnt as f64 {
                        tier = 1;
                    } else {
                        if cnt > 0 {
                            report.n_cancelled += 1;
                        }
                        // Tier 2: neighbours in a face across a side of fc.
                        sum = C64::default();
                        cnt = 0;
                        for k in 0..3 {
                            let Some((g, gpos)) = self.unfold(&layouts[fc], fc, k) else {
                                continue;
                            };
                            for nb in fine.neighbours(i) {
                                self.candidates(nb.vertex as usize, g, &mut cand);
                                for b in &cand {
                                    let d = place(&gpos, *b);
                                    accumulate(
                                        &mut sum,
                                        &mut cnt,
                                        pi,
                                        d,
                                        nb.phi_ij,
                                        eps,
                                        &mut report,
                                    );
                                }
                            }
                        }
                        if cnt > 0 && sum.abs() >= MIN_AVERAGE * cnt as f64 {
                            tier = 2;
                        } else if cnt > 0 {
                            report.n_cancelled += 1;
                        }
                    }
                    if tier > 0 {
                        if tier == 1 {
                            report.n_same_face += 1;
                        } else {
                            report.n_adjacent_face += 1;
                        }
                        let psi = sum.arg();
                        rows[i] = Some(
                            (0..3)
                                .map(|k| {
                                    (
                                        self.faces[fc][k],
                                        C64::from_angle(n * (psi + z_angle[fc][k])).scale(bary[k]),
                                    )
                                })
                                .collect(),
                        );
                    }
                }
                CoarsePoint::Vertex(v) => {
                    let mut sum = C64::default();
                    let mut cnt = 0usize;
                    for &(f, s) in &slots[slot_ptr[v as usize]..slot_ptr[v as usize + 1]] {
                        let (f, s) = (f as usize, s as usize);
                        let pi = layouts[f][s];
                        let eps = 1e-12 * self.face_lengths[f].iter().cloned().fold(0.0, f64::max);
                        for nb in fine.neighbours(i) {
                            self.candidates(nb.vertex as usize, f, &mut cand);
                            for b in &cand {
                                let pj = place(&layouts[f], *b);
                                let d = [pj[0] - pi[0], pj[1] - pi[1]];
                                if d[0].hypot(d[1]) <= eps {
                                    report.n_coincident_neighbours += 1;
                                    continue;
                                }
                                // psi z for this face and neighbour, a unit complex number.
                                let a = nb.phi_ij - d[1].atan2(d[0]) + z_angle[f][s];
                                sum += C64::from_angle(a);
                                cnt += 1;
                            }
                        }
                    }
                    if cnt > 0 && sum.abs() >= MIN_AVERAGE * cnt as f64 {
                        report.n_kept_from_images += 1;
                        rows[i] = Some(vec![(v, C64::from_angle(n * sum.arg()))]);
                    } else if cnt > 0 {
                        report.n_cancelled += 1;
                    }
                }
            }
        }

        // Tier 3: breadth-first average of the neighbours' rows.
        loop {
            let mut filled: Vec<(usize, Row)> = Vec::new();
            let mut unknown = 0usize;
            for i in 0..nfine {
                if rows[i].is_some() {
                    continue;
                }
                unknown += 1;
                let mut acc: Vec<(u32, C64)> = Vec::new();
                let mut m = 0usize;
                for nb in fine.neighbours(i) {
                    if let Some(r) = &rows[nb.vertex as usize] {
                        let rot = C64::from_angle(n * (nb.phi_ij - nb.phi_ji - PI));
                        acc.extend(r.iter().map(|&(c, a)| (c, a * rot)));
                        m += 1;
                    }
                }
                if m > 0 {
                    let inv = 1.0 / m as f64;
                    for e in &mut acc {
                        e.1 = e.1.scale(inv);
                    }
                    filled.push((i, acc));
                }
            }
            if unknown == 0 {
                break;
            }
            if filled.is_empty() {
                // Nothing known in these components: the own corners with psi = 1.
                for i in 0..nfine {
                    if rows[i].is_none() {
                        report.n_default += 1;
                        rows[i] = Some(match self.points[i] {
                            CoarsePoint::Vertex(v) => vec![(v, C64::new(1.0, 0.0))],
                            CoarsePoint::Face { face, bary } => (0..3)
                                .map(|k| {
                                    (
                                        self.faces[face as usize][k],
                                        C64::from_angle(n * z_angle[face as usize][k])
                                            .scale(bary[k]),
                                    )
                                })
                                .collect(),
                        });
                    }
                }
                break;
            }
            report.averaging_rounds += 1;
            for (i, r) in filled {
                match self.points[i] {
                    CoarsePoint::Vertex(_) => report.n_kept_averaged += 1,
                    CoarsePoint::Face { .. } => report.n_averaged += 1,
                }
                rows[i] = Some(r);
            }
        }

        // Assemble: sort each row by column, sum repeated corners.
        let mut row_ptr = vec![0usize];
        let mut col_idx = Vec::new();
        let mut values = Vec::new();
        for r in rows {
            let mut r = r.expect("every row is filled above");
            r.sort_by_key(|e| e.0);
            let start = col_idx.len();
            for (c, a) in r {
                if col_idx.len() > start && *col_idx.last().unwrap() == c {
                    *values.last_mut().unwrap() += a;
                } else {
                    col_idx.push(c);
                    values.push(a);
                }
            }
            row_ptr.push(col_idx.len());
        }
        Ok((
            ComplexSparse {
                nrows: nfine,
                ncols: self.vertices.len(),
                row_ptr,
                col_idx,
                values,
            },
            report,
        ))
    }

    /// The barycentric coordinates of the image of fine point `j` in face `g`: the point's own
    /// coordinates if it is in `g`, a unit vector for each corner of `g` that is the vertex
    /// the point is. Empty otherwise.
    fn candidates(&self, j: usize, g: usize, out: &mut Vec<[f64; 3]>) {
        out.clear();
        match self.points[j] {
            CoarsePoint::Face { face, bary } => {
                if face as usize == g {
                    out.push(bary);
                }
            }
            CoarsePoint::Vertex(v) => {
                for k in 0..3 {
                    if self.faces[g][k] == v {
                        let mut b = [0.0; 3];
                        b[k] = 1.0;
                        out.push(b);
                    }
                }
            }
        }
    }

    /// The face across side `k` of face `f` (laid out as `lay`), if any, with the positions of
    /// its slots in the frame of `f`: the shared side is traversed the other way, and the
    /// apex is on the far side of it at the distances of the face's lengths.
    fn unfold(&self, lay: &[[f64; 2]; 3], f: usize, k: usize) -> Option<(usize, [[f64; 2]; 3])> {
        let t = self.twin[3 * f + k]? as usize;
        let (g, s) = (t / 3, t % 3);
        let a = lay[(k + 1) % 3];
        let b = lay[k];
        let lg = self.face_lengths[g];
        let d = lg[s];
        let (ra, rb) = (lg[(s + 2) % 3], lg[(s + 1) % 3]);
        let ab = [b[0] - a[0], b[1] - a[1]];
        let len = ab[0].hypot(ab[1]);
        if !(len > 0.0 && d > 0.0) {
            return None;
        }
        let u = [ab[0] / len, ab[1] / len];
        let x = (d * d + ra * ra - rb * rb) / (2.0 * d);
        let h = (ra * ra - x * x).max(0.0).sqrt();
        let q = [a[0] + x * u[0] - h * u[1], a[1] + x * u[1] + h * u[0]];
        let mut pos = [[0.0; 2]; 3];
        pos[s] = a;
        pos[(s + 1) % 3] = b;
        pos[(s + 2) % 3] = q;
        Some((g, pos))
    }
}

/// The point with barycentric coordinates `b` in the layout `lay`.
#[inline]
fn place(lay: &[[f64; 2]; 3], b: [f64; 3]) -> [f64; 2] {
    [
        b[0] * lay[0][0] + b[1] * lay[1][0] + b[2] * lay[2][0],
        b[0] * lay[0][1] + b[1] * lay[1][1] + b[2] * lay[2][1],
    ]
}

/// Adds `ψ_ij = e^{iφ_ij} / w̄_ij` for the neighbour at `pj` to `sum` (skipped, and counted,
/// if its image coincides with `pi`).
#[inline]
fn accumulate(
    sum: &mut C64,
    cnt: &mut usize,
    pi: [f64; 2],
    pj: [f64; 2],
    phi_ij: f64,
    eps: f64,
    report: &mut VectorProlongationReport,
) {
    let d = [pj[0] - pi[0], pj[1] - pi[1]];
    if d[0].hypot(d[1]) <= eps {
        report.n_coincident_neighbours += 1;
        return;
    }
    *sum += C64::from_angle(phi_ij - d[1].atan2(d[0]));
    *cnt += 1;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CoarseMesh, CoarseMeshOptions};

    /// Unfolding the face across a side of a flat mesh gives the planar positions up to the
    /// rigid motion of the first face.
    #[test]
    fn unfold_is_the_planar_layout_of_the_neighbour() {
        let n = 5usize;
        let id = |i: usize, j: usize| (j * (n + 1) + i) as u32;
        let mut positions = Vec::new();
        for j in 0..=n {
            for i in 0..=n {
                let w = ((i * 7 + j * 13) % 5) as f64 * 0.05;
                positions.push([i as f64 + w, j as f64 - w * 0.5, 0.0]);
            }
        }
        let mut faces = Vec::new();
        for j in 0..n {
            for i in 0..n {
                let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i, j + 1), id(i + 1, j + 1));
                faces.push([a, b, d]);
                faces.push([a, d, c]);
            }
        }
        let mesh =
            CoarseMesh::from_positions(&positions, &faces, &[], CoarseMeshOptions::default())
                .unwrap();
        let out = mesh.compact().unwrap();
        let mut checked = 0;
        for f in 0..out.n_faces() {
            let lay = triangle_layout(out.face_lengths[f]).unwrap();
            // The rigid motion from the layout of f to the plane: corner 0 and the direction of
            // side 0.
            let p = |v: u32| positions[out.vertices[v as usize] as usize];
            let (o, e) = (p(out.faces[f][0]), p(out.faces[f][1]));
            let ang = (e[1] - o[1]).atan2(e[0] - o[0]);
            let (s, c) = ang.sin_cos();
            let to_plane = |q: [f64; 2]| [o[0] + c * q[0] - s * q[1], o[1] + s * q[0] + c * q[1]];
            for k in 0..3 {
                let Some((g, pos)) = out.unfold(&lay, f, k) else {
                    continue;
                };
                for m in 0..3 {
                    let want = p(out.faces[g][m]);
                    let got = to_plane(pos[m]);
                    assert!(
                        (want[0] - got[0]).abs() < 1e-12 && (want[1] - got[1]).abs() < 1e-12,
                        "face {f} side {k}: slot {m} of face {g}: {got:?} vs {want:?}"
                    );
                }
                checked += 1;
            }
        }
        assert!(checked > 40);
    }
}
