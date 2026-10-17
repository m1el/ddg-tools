//! Complex numbers and complex Hermitian matrices stored as real symmetric ones.
//!
//! A complex vector `u ∈ Cⁿ` is stored as `2n` reals, `(Re u_i, Im u_i)` at `2i, 2i + 1`. A
//! complex matrix entry `a + ib` at `(i, j)` becomes the block `[[a, −b], [b, a]]` at rows
//! `2i, 2i + 1` and columns `2j, 2j + 1`: then the real matrix acts on the stored vector as the
//! complex matrix acts on `u`, a Hermitian matrix becomes a symmetric one, and
//! `Re(u* A v) = uᵀ A v` in stored form. Every eigenvalue of a Hermitian matrix appears twice
//! (for `u` and `i u`), and the real matrix is positive (semi)definite exactly when the complex
//! one is.

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

use ddg_core::{SparseError, SparseMatrix};

/// A complex number.
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

    /// The complex conjugate.
    #[inline]
    pub fn conj(self) -> C64 {
        C64 {
            re: self.re,
            im: -self.im,
        }
    }

    /// `|z|²`.
    #[inline]
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// `|z|`.
    #[inline]
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    /// The argument in `(−π, π]`.
    #[inline]
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }

    /// `s z` for a real `s`.
    #[inline]
    pub fn scale(self, s: f64) -> C64 {
        C64 {
            re: s * self.re,
            im: s * self.im,
        }
    }

    /// `zⁿ` for `n ≥ 0`.
    pub fn powi(self, n: u32) -> C64 {
        let mut out = C64::new(1.0, 0.0);
        for _ in 0..n {
            out = out * self;
        }
        out
    }
}

impl Add for C64 {
    type Output = C64;
    #[inline]
    fn add(self, o: C64) -> C64 {
        C64::new(self.re + o.re, self.im + o.im)
    }
}

impl AddAssign for C64 {
    #[inline]
    fn add_assign(&mut self, o: C64) {
        self.re += o.re;
        self.im += o.im;
    }
}

impl Sub for C64 {
    type Output = C64;
    #[inline]
    fn sub(self, o: C64) -> C64 {
        C64::new(self.re - o.re, self.im - o.im)
    }
}

impl Neg for C64 {
    type Output = C64;
    #[inline]
    fn neg(self) -> C64 {
        C64::new(-self.re, -self.im)
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
        let d = o.norm_sqr();
        let n = self * o.conj();
        C64::new(n.re / d, n.im / d)
    }
}

/// Collects the entries of a complex `n × n` matrix and builds its real `2n × 2n` form.
#[derive(Debug, Clone, Default)]
pub(crate) struct ComplexTriplets {
    n: usize,
    triplets: Vec<(u32, u32, f64)>,
}

impl ComplexTriplets {
    pub(crate) fn new(n: usize) -> Self {
        ComplexTriplets {
            n,
            triplets: Vec::new(),
        }
    }

    /// Adds `z` to entry `(i, j)`.
    pub(crate) fn add(&mut self, i: usize, j: usize, z: C64) {
        let (r, c) = (2 * i as u32, 2 * j as u32);
        self.triplets.push((r, c, z.re));
        self.triplets.push((r, c + 1, -z.im));
        self.triplets.push((r + 1, c, z.im));
        self.triplets.push((r + 1, c + 1, z.re));
    }

    pub(crate) fn build(&self) -> Result<SparseMatrix, SparseError> {
        SparseMatrix::from_triplets(2 * self.n, 2 * self.n, &self.triplets)
    }
}

/// The stored form of complex values.
pub(crate) fn to_real(values: &[C64]) -> Vec<f64> {
    values.iter().flat_map(|z| [z.re, z.im]).collect()
}

/// Complex values from their stored form.
pub(crate) fn from_real(values: &[f64]) -> Vec<C64> {
    values
        .chunks_exact(2)
        .map(|p| C64::new(p[0], p[1]))
        .collect()
}
