//! The smooth perturbed torus of Sec. 8.1 (Fig. 10): `f = f₀ + Φ∘f₀` with `f₀` the standard
//! torus and `Φ: R³ → R³` a seeded trigonometric polynomial. The paper does not give `Φ`; ours is
//! `Φ(x) = amplitude Σ_k (A_k cos(k·x) + B_k sin(k·x)) / (1 + |k|²)` over the integer vectors
//! `k ≠ 0` with `|k|∞ ≤ degree` (one of `±k`), with `A_k, B_k ∈ [−1, 1]³` from a SplitMix64 stream
//! (SPEC D31-style analogue; WP16 records it).

use std::f64::consts::TAU;

use super::{Mesh, Rng, bad};
use crate::error::RsError;

type Term = ([f64; 3], [f64; 3], [f64; 3], f64);

/// A smooth embedded torus given parametrically.
#[derive(Debug, Clone)]
pub struct PerturbedTorus {
    major: f64,
    minor: f64,
    terms: Vec<Term>,
}

impl PerturbedTorus {
    /// The torus with radii `major > minor > 0` perturbed by a random polynomial of the given
    /// `degree ≥ 1`, `amplitude ≥ 0`.
    ///
    /// # Errors
    /// Invalid radii, degree or amplitude.
    pub fn new(major: f64, minor: f64, degree: u32, amplitude: f64, seed: u64) -> Result<Self, RsError> {
        if !(minor > 0.0 && major > minor) {
            return Err(bad("radii", format!("{major} {minor}"), "0 < minor < major"));
        }
        if degree == 0 || !(amplitude.is_finite() && amplitude >= 0.0) {
            return Err(bad("degree/amplitude", format!("{degree} {amplitude}"), "degree >= 1, amplitude >= 0"));
        }
        let mut rng = Rng::new(seed);
        let d = degree as i32;
        let mut terms = Vec::new();
        for a in -d..=d {
            for b in -d..=d {
                for c in -d..=d {
                    let first = [a, b, c].into_iter().find(|&x| x != 0);
                    if first.is_none_or(|x| x < 0) {
                        continue; // k = 0 or the negative of a kept vector
                    }
                    let k = [a as f64, b as f64, c as f64];
                    let w = amplitude / (1.0 + k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
                    let av = [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)];
                    let bv = [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)];
                    terms.push((k, av, bv, w));
                }
            }
        }
        Ok(PerturbedTorus { major, minor, terms })
    }

    /// The unperturbed torus `f₀(u, v)`.
    pub fn base(&self, u: f64, v: f64) -> [f64; 3] {
        let ring = self.major + self.minor * v.cos();
        [ring * u.cos(), ring * u.sin(), self.minor * v.sin()]
    }

    /// `Φ(x)`.
    pub fn phi(&self, x: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0; 3];
        for (k, a, b, w) in &self.terms {
            let t = k[0] * x[0] + k[1] * x[1] + k[2] * x[2];
            let (s, c) = t.sin_cos();
            for i in 0..3 {
                out[i] += w * (a[i] * c + b[i] * s);
            }
        }
        out
    }

    /// `f(u, v) = f₀(u, v) + Φ(f₀(u, v))`, `2π`-periodic in both parameters.
    pub fn point(&self, u: f64, v: f64) -> [f64; 3] {
        let p = self.base(u, v);
        let d = self.phi(p);
        [p[0] + d[0], p[1] + d[1], p[2] + d[2]]
    }

    /// Samples the surface on an `nu × nv` grid of the parameter domain `[0, 2π)²`. The
    /// parameters are warped by `x ↦ x − s sin(2π x)/(2π)` on `[0, 1]` (`s = squeeze_u`,
    /// `squeeze_v` in `[0, 1)`; 0: uniform), which squeezes the grid towards `x = 0` and `1`
    /// and stretches it at `½`; `nu ≠ nv` gives anisotropic cells. `2 nu nv` faces, outward.
    ///
    /// # Errors
    /// `nu, nv < 3` or a squeeze outside `[0, 1)`.
    pub fn sample(&self, nu: usize, nv: usize, squeeze_u: f64, squeeze_v: f64) -> Result<Mesh, RsError> {
        if nu < 3 || nv < 3 {
            return Err(bad("grid", format!("{nu}x{nv}"), "at least 3 per direction"));
        }
        if !((0.0..1.0).contains(&squeeze_u) && (0.0..1.0).contains(&squeeze_v)) {
            return Err(bad("squeeze", format!("{squeeze_u} {squeeze_v}"), "in [0, 1)"));
        }
        let warp = |x: f64, s: f64| TAU * (x - s * (TAU * x).sin() / TAU);
        let mut positions = Vec::with_capacity(nu * nv);
        for i in 0..nu {
            for j in 0..nv {
                positions.push(self.point(
                    warp(i as f64 / nu as f64, squeeze_u),
                    warp(j as f64 / nv as f64, squeeze_v),
                ));
            }
        }
        let id = |i: usize, j: usize| ((i % nu) * nv + j % nv) as u32;
        let mut faces = Vec::with_capacity(2 * nu * nv);
        for i in 0..nu {
            for j in 0..nv {
                faces.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
                faces.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
            }
        }
        Mesh { positions, faces }.validated()
    }
}
