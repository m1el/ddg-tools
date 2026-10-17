//! The local matrices of the finite element connection Laplacian of Knöppel, Crane, Pinkall,
//! Schröder (2013, Sec. 6.1.1 and App. D): piecewise linear basis sections on a triangle that
//! carries a connection of constant curvature.
//!
//! # The model (derived here; SPEC N43)
//!
//! Triangle `T` with corners `p₀, p₁, p₂` (counter-clockwise, in the plane of its layout),
//! area `A`, holonomy `Ω` (transport counter-clockwise around `T` multiplies by `e^{iΩ}`),
//! curvature density `κ = Ω / A`. In the gauge centred at a corner `p` the connection is
//! `∇ψ = dψ − i A_p ψ` with `A_p(w) = (κ/2) (x − p) × w`, which vanishes along rays from `p`;
//! the section `ξ_m` that is parallel along rays from `p_m` with `ξ_m(p_m) = X_m` (the frame
//! vector of the vertex) is constant in the gauge of `p_m`, and the basis section is
//! `ψ_m = b_m ξ_m` (`b_m` the barycentric coordinate). Transport around the triangle
//! `p_j → x → p_k → p_j` (clockwise, area fraction `b_i`) gives, for `(i, j, k)` counter-clockwise,
//!
//! ```text
//! ⟨ψ_j, ψ_k⟩(x)   = r̄_jk e^{iΩ b_i} b_j b_k
//! ⟨∇ψ_j, ∇ψ_k⟩(x) = r̄_jk e^{iΩ b_i} [∇b_j·∇b_k + i (b_j a_j·∇b_k − b_k a_k·∇b_j) + b_j b_k a_j·a_k]
//! ```
//!
//! with `a_m(x) = (κ/2) J (x − p_m)`, `J` the rotation by `+90°`, and `r_jk` the transport
//! from `j` to `k` along the side. Both are polynomials in the barycentric coordinates times
//! `e^{iΩ b_i}`; they are integrated exactly, monomial by monomial
//! ([`monomial_integral`]). For `Ω → 0` the mass entries tend to `r̄_jk A/12` and the
//! stiffness entries to `−½ r̄_jk cot θ_i`: the cotan Laplacian. Knöppel et al.'s closed forms
//! (their Eq. 17 and App. D.3) are the same integrals; `tests` compares the mass entry and the
//! diagonal with them and the whole local matrices with a numerical quadrature of the model.

use super::complex::C64;

/// Polynomials of total degree at most 4 in the barycentric coordinates `(b₀, b₁, b₂)`.
#[derive(Clone, Copy)]
struct Poly([f64; 125]);

const fn index(e: [usize; 3]) -> usize {
    e[0] * 25 + e[1] * 5 + e[2]
}

impl Poly {
    fn zero() -> Poly {
        Poly([0.0; 125])
    }

    fn constant(c: f64) -> Poly {
        let mut p = Poly::zero();
        p.0[0] = c;
        p
    }

    fn linear(c: [f64; 3]) -> Poly {
        let mut p = Poly::zero();
        p.0[index([1, 0, 0])] = c[0];
        p.0[index([0, 1, 0])] = c[1];
        p.0[index([0, 0, 1])] = c[2];
        p
    }

    fn monomials(&self) -> impl Iterator<Item = ([usize; 3], f64)> + '_ {
        (0..125).filter_map(move |i| {
            let c = self.0[i];
            (c != 0.0).then_some(([i / 25, (i / 5) % 5, i % 5], c))
        })
    }

    fn mul(&self, other: &Poly) -> Poly {
        let mut out = Poly::zero();
        for (e, a) in self.monomials() {
            for (f, b) in other.monomials() {
                let g = [e[0] + f[0], e[1] + f[1], e[2] + f[2]];
                assert!(g.iter().sum::<usize>() <= 4, "degree above 4");
                out.0[index(g)] += a * b;
            }
        }
        out
    }

    fn add(&self, other: &Poly) -> Poly {
        let mut out = *self;
        for i in 0..125 {
            out.0[i] += other.0[i];
        }
        out
    }

    fn scale(&self, s: f64) -> Poly {
        let mut out = *self;
        for v in &mut out.0 {
            *v *= s;
        }
        out
    }

    /// `∫_T self · e^{iΩ b_c} dA`.
    fn integrate(&self, c: usize, omega: f64, area: f64) -> C64 {
        let mut sum = C64::default();
        for (e, a) in self.monomials() {
            sum += monomial_integral(e, c, omega, area).scale(a);
        }
        sum
    }
}

fn factorial(n: usize) -> f64 {
    (1..=n).map(|k| k as f64).product()
}

/// `∫₀¹ s^p (1 − s)^m e^{iΩ s} ds`.
///
/// For `|Ω| < 3` by the power series of the exponential (`Σ_k (iΩ)^k/k! · B(p + k + 1, m + 1)`,
/// whose terms grow to at most `e^{|Ω|}` times the sum); otherwise by expanding `(1 − s)^m` and
/// the recurrence `K_k = (e^{iΩ} − k K_{k−1}) / (iΩ)` for `K_k = ∫₀¹ s^k e^{iΩ s} ds`, which
/// is stable for `k < |Ω|` and amplifies rounding by at most `k!/|Ω|^k` beyond.
fn beta_exp(p: usize, m: usize, omega: f64) -> C64 {
    if omega.abs() < 3.0 {
        // Term k: (iΩ)^k p_k! m! / ((p + k + m + 1)! ) · (p + k)!/(k! p_k!) ...: written as a
        // running product.
        // t_k = (iΩ)^k / k! · (p + k)! m! / (p + k + m + 1)!
        let mut t = C64::new(factorial(p) * factorial(m) / factorial(p + m + 1), 0.0);
        let mut sum = t;
        for k in 0..200 {
            // t_{k+1} / t_k = iΩ (p + k + 1) / ((k + 1)(p + k + m + 2)).
            let ratio = omega * (p + k + 1) as f64 / ((k + 1) as f64 * (p + k + m + 2) as f64);
            t = C64::new(0.0, ratio) * t;
            sum += t;
            if t.abs() < 1e-18 * sum.abs().max(1e-300) {
                break;
            }
        }
        sum
    } else {
        let e = C64::from_angle(omega);
        let i_omega = C64::new(0.0, omega);
        let mut k_vals = Vec::with_capacity(p + m + 1);
        k_vals.push((e - C64::new(1.0, 0.0)) / i_omega);
        for k in 1..=p + m {
            let prev = k_vals[k - 1];
            k_vals.push((e - prev.scale(k as f64)) / i_omega);
        }
        let mut sum = C64::default();
        let mut binom = 1.0;
        for l in 0..=m {
            let sign = if l % 2 == 0 { 1.0 } else { -1.0 };
            sum += k_vals[p + l].scale(sign * binom);
            binom = binom * (m - l) as f64 / (l + 1) as f64;
        }
        sum
    }
}

/// `∫_T b₀^{e₀} b₁^{e₁} b₂^{e₂} e^{iΩ b_c} dA` over a triangle of area `area`.
///
/// With `s = b_c` and `q, r` the other two exponents the integral over the slices `b_c = s`
/// gives `2A · q! r! / (q + r + 1)! · ∫₀¹ s^{e_c} (1 − s)^{q + r + 1} e^{iΩ s} ds`.
pub(crate) fn monomial_integral(e: [usize; 3], c: usize, omega: f64, area: f64) -> C64 {
    let (q, r) = (e[(c + 1) % 3], e[(c + 2) % 3]);
    let weight = 2.0 * area * factorial(q) * factorial(r) / factorial(q + r + 1);
    beta_exp(e[c], q + r + 1, omega).scale(weight)
}

/// The local mass and stiffness matrices of one triangle.
///
/// `p` is the layout of the triangle (counter-clockwise), `area` its area, `omega` the
/// holonomy of the bundle over it, `r[a]` the transport coefficient from corner `a` to corner
/// `a + 1` (for `n`-vectors: the `n`-th power). Returns `(M, A)` with
/// `M[a][b] = ∫ ⟨ψ_a, ψ_b⟩` and `A[a][b] = ∫ ⟨∇ψ_a, ∇ψ_b⟩` (antilinear in the first slot), so
/// that `‖Σ u_a ψ_a‖² = u* M u` and the Dirichlet energy is `u* A u`.
///
/// [formula: app-vector-fem-connection-laplacian]
pub(crate) fn local_matrices(
    p: [[f64; 2]; 3],
    area: f64,
    omega: f64,
    r: [C64; 3],
) -> ([[C64; 3]; 3], [[C64; 3]; 3]) {
    let kappa = omega / area;
    let dot = |u: [f64; 2], v: [f64; 2]| u[0] * v[0] + u[1] * v[1];
    let sub = |u: [f64; 2], v: [f64; 2]| [u[0] - v[0], u[1] - v[1]];
    // d[m] = p_{m+2} − p_{m+1}; ∇b_m = J d[m] / (2A).
    let d: [[f64; 2]; 3] = [0, 1, 2].map(|m| sub(p[(m + 2) % 3], p[(m + 1) % 3]));
    let grad_dot = |a: usize, b: usize| dot(d[a], d[b]) / (4.0 * area * area);
    // (x − p_j) · v as a linear form in b.
    let rel_dot = |j: usize, v: [f64; 2]| [0, 1, 2].map(|m| dot(sub(p[m], p[j]), v));
    // a_j · ∇b_k = (κ/2) (x − p_j) · d[k] / (2A).
    let a_dot_grad = |j: usize, k: usize| {
        Poly::linear(rel_dot(j, d[k])).scale(kappa / (4.0 * area))
    };
    // a_j · a_k = (κ²/4) (x − p_j)·(x − p_k).
    let a_dot_a = |j: usize, k: usize| {
        let mut out = Poly::zero();
        for m in 0..3 {
            for l in 0..3 {
                let c = dot(sub(p[m], p[j]), sub(p[l], p[k]));
                let mut e = [0usize; 3];
                e[m] += 1;
                e[l] += 1;
                out.0[index(e)] += c;
            }
        }
        out.scale(0.25 * kappa * kappa)
    };
    let b = |m: usize| {
        let mut c = [0.0; 3];
        c[m] = 1.0;
        Poly::linear(c)
    };

    let mut mass = [[C64::default(); 3]; 3];
    let mut stiff = [[C64::default(); 3]; 3];
    for a in 0..3 {
        mass[a][a] = C64::new(area / 6.0, 0.0);
        let extra = b(a).mul(&b(a)).mul(&a_dot_a(a, a)).integrate(0, 0.0, area);
        stiff[a][a] = C64::new(area * grad_dot(a, a) + extra.re, 0.0);
    }
    for c in 0..3 {
        // (c, j, k) counter-clockwise.
        let (j, k) = ((c + 1) % 3, (c + 2) % 3);
        let rbar = r[j].conj();
        let m = b(j).mul(&b(k)).integrate(c, omega, area);
        mass[j][k] = rbar * m;
        mass[k][j] = mass[j][k].conj();
        let real = Poly::constant(grad_dot(j, k)).add(&b(j).mul(&b(k)).mul(&a_dot_a(j, k)));
        let imag = b(j)
            .mul(&a_dot_grad(j, k))
            .add(&b(k).mul(&a_dot_grad(k, j)).scale(-1.0));
        let integral = real.integrate(c, omega, area) + C64::new(0.0, 1.0) * imag.integrate(c, omega, area);
        stiff[j][k] = rbar * integral;
        stiff[k][j] = stiff[j][k].conj();
    }
    (mass, stiff)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn layout(l: [f64; 3]) -> ([[f64; 2]; 3], f64) {
        let p = ddg_intrinsic::formulas::triangle_layout(l).unwrap();
        let area = 0.5 * (p[1][0] * p[2][1] - p[1][1] * p[2][0]);
        (p, area)
    }

    #[test]
    fn both_evaluations_of_the_one_dimensional_integral_agree() {
        for p in 0..5 {
            for m in 0..6 {
                for omega in [2.9999999, 3.0, -3.0, 0.7, 5.5, 12.0] {
                    // Quadrature: composite Simpson with many intervals.
                    let n = 20000;
                    let h = 1.0 / n as f64;
                    let mut q = C64::default();
                    for i in 0..=n {
                        let s = i as f64 * h;
                        let w = if i == 0 || i == n {
                            1.0
                        } else if i % 2 == 1 {
                            4.0
                        } else {
                            2.0
                        };
                        let f = C64::from_angle(omega * s)
                            .scale(s.powi(p as i32) * (1.0 - s).powi(m as i32));
                        q += f.scale(w * h / 3.0);
                    }
                    let v = beta_exp(p, m, omega);
                    assert!((v - q).abs() < 1e-12, "p {p} m {m} Ω {omega}: {v:?} vs {q:?}");
                }
            }
        }
    }

    /// Knöppel et al. Eq. 17: the off-diagonal mass entry, and its flat limit.
    #[test]
    fn the_mass_matrix_is_knoeppels() {
        let (p, area) = layout([1.0, 0.8, 0.7]);
        for omega in [0.05, 0.3, -1.1, 2.5, std::f64::consts::PI] {
            let r = [0.3f64, -0.2, 1.4].map(C64::from_angle);
            let (m, _) = local_matrices(p, area, omega, r);
            let e = C64::from_angle(omega);
            let o = C64::new(0.0, omega);
            let num = e.scale(6.0) - C64::new(6.0, 0.0) - o.scale(6.0)
                + C64::new(3.0 * omega * omega, 0.0)
                + C64::new(0.0, omega.powi(3));
            let closed = num.scale(area / (3.0 * omega.powi(4)));
            for c in 0..3 {
                let (j, k) = ((c + 1) % 3, (c + 2) % 3);
                let expected = r[j].conj() * closed;
                // The closed form cancels for small Ω; compare with a matching tolerance.
                let tol = if omega.abs() < 0.1 { 1e-8 } else { 1e-13 };
                assert!((m[j][k] - expected).abs() < tol, "{:?} {:?}", m[j][k], expected);
            }
        }
        let (m, _) = local_matrices(p, area, 0.0, [C64::new(1.0, 0.0); 3]);
        assert!((m[0][1].re - area / 12.0).abs() < 1e-15);
    }

    /// For Ω = 0 the stiffness matrix is the cotan matrix times the transport.
    #[test]
    fn the_flat_limit_is_the_cotan_matrix() {
        let l = [1.0, 0.8, 0.7];
        let (p, area) = layout(l);
        let r = [0.3f64, -0.2, 1.4].map(C64::from_angle);
        let (_, a) = local_matrices(p, area, 0.0, r);
        for c in 0..3 {
            let (j, k) = ((c + 1) % 3, (c + 2) % 3);
            let theta = ddg_core::geometry::corner_cotan(l[j], l[c], l[k]).unwrap();
            let expected = r[j].conj().scale(-0.5 * theta);
            assert!((a[j][k] - expected).abs() < 1e-14);
        }
    }

    /// The diagonal of Knöppel et al. App. D.3, `(1/4A) [|p_jk|² + Ω² (3g₁₁ − 3g₁₂ + g₂₂)/90]`
    /// written with the edge vectors at the corner, compared numerically: it must equal
    /// `A |∇b|² + (κ²/4) ∫ b² |x − p|²`.
    #[test]
    fn the_diagonal_has_the_curvature_term() {
        let (p, area) = layout([1.0, 0.8, 0.7]);
        let omega = 1.3;
        let (_, a) = local_matrices(p, area, omega, [C64::new(1.0, 0.0); 3]);
        for i in 0..3 {
            // Quadrature of b_i² |x − p_i|² over the triangle (Duffy, Gauss–Legendre is
            // overkill here: the integrand is a polynomial of degree 4, use the monomial rule
            // with an independent expansion).
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            let u = [p[j][0] - p[i][0], p[j][1] - p[i][1]];
            let v = [p[k][0] - p[i][0], p[k][1] - p[i][1]];
            let (g11, g12, g22) = (
                u[0] * u[0] + u[1] * u[1],
                u[0] * v[0] + u[1] * v[1],
                v[0] * v[0] + v[1] * v[1],
            );
            // x − p_i = b_j u + b_k v; ∫ b_i² b_j² = ∫ b_i² b_k² = 2A·2!2!/6! = A/90,
            // ∫ b_i² b_j b_k = 2A·2!/6! = A/180.
            let int = area / 90.0 * g11 + area / 90.0 * g22 + 2.0 * g12 * area / 180.0;
            let kappa = omega / area;
            let d = [p[k][0] - p[j][0], p[k][1] - p[j][1]];
            let expected = (d[0] * d[0] + d[1] * d[1]) / (4.0 * area) + 0.25 * kappa * kappa * int;
            assert!((a[i][i].re - expected).abs() < 1e-13);
            assert_eq!(a[i][i].im, 0.0);
        }
    }

    /// The whole local matrices against a direct quadrature of the model: the basis sections
    /// in one gauge (centred at the origin), their covariant derivatives by central
    /// differences.
    #[test]
    fn the_local_matrices_match_a_quadrature_of_the_model() {
        let (p, area) = layout([1.0, 0.9, 0.6]);
        let omega = 2.2;
        let kappa = omega / area;
        // Frame vectors X_m of the three corners in the gauge at the origin: arbitrary.
        let x_frame = [0.4f64, 2.0, -1.0].map(C64::from_angle);
        // In the gauge A_o(w) = (κ/2) x × w, ξ_m(x) = X_m e^{-iχ_m(x)} with
        // χ_m(x) = (κ/2) ((0 − p_m.x) y − (0 − p_m.y) x) ... evaluated against p_m.
        let chi = |m: usize, x: [f64; 2]| {
            0.5 * kappa * ((-p[m][0]) * x[1] - (-p[m][1]) * x[0])
                - 0.5 * kappa * ((-p[m][0]) * p[m][1] - (-p[m][1]) * p[m][0])
        };
        let bary = |x: [f64; 2]| ddg_intrinsic::formulas::barycentric(&p, x);
        let psi = |m: usize, x: [f64; 2]| x_frame[m] * C64::from_angle(-chi(m, x)).scale(bary(x)[m]);
        // Transport from corner a to a + 1 along the side, in the gauge at the origin:
        // ∫ A_o along the straight side from p_a to p_b is (κ/2) p_a × p_b.
        let r: [C64; 3] = [0, 1, 2].map(|a| {
            let bb = (a + 1) % 3;
            let phase = 0.5 * kappa * (p[a][0] * p[bb][1] - p[a][1] * p[bb][0]);
            // Transported X_a in gauge-o coordinates: X_a e^{i phase}; its coordinate in the
            // frame X_b is that over X_b.
            x_frame[a] * C64::from_angle(phase) / x_frame[bb]
        });
        let (m_local, a_local) = local_matrices(p, area, omega, r);
        // Quadrature on a fine barycentric grid (midpoints of sub-triangles).
        let n = 200;
        let mut mass = [[C64::default(); 3]; 3];
        let mut stiff = [[C64::default(); 3]; 3];
        let h = 1e-6;
        let cell = area / (n * n) as f64;
        for i in 0..n {
            for j in 0..n - i {
                for up in [true, false] {
                    if !up && j + i + 1 >= n {
                        continue;
                    }
                    let (s, t) = if up {
                        ((i as f64 + 1.0 / 3.0) / n as f64, (j as f64 + 1.0 / 3.0) / n as f64)
                    } else {
                        ((i as f64 + 2.0 / 3.0) / n as f64, (j as f64 + 2.0 / 3.0) / n as f64)
                    };
                    let x = [
                        p[0][0] + s * (p[1][0] - p[0][0]) + t * (p[2][0] - p[0][0]),
                        p[0][1] + s * (p[1][1] - p[0][1]) + t * (p[2][1] - p[0][1]),
                    ];
                    let cov = |m: usize, dir: usize| {
                        let mut xp = x;
                        let mut xm = x;
                        xp[dir] += h;
                        xm[dir] -= h;
                        let deriv = (psi(m, xp) - psi(m, xm)).scale(0.5 / h);
                        // ∇ψ = dψ − i A_o ψ, A_o = (κ/2)(−y, x).
                        let a_o = if dir == 0 { -0.5 * kappa * x[1] } else { 0.5 * kappa * x[0] };
                        deriv - C64::new(0.0, a_o) * psi(m, x)
                    };
                    for a in 0..3 {
                        for b in 0..3 {
                            mass[a][b] += (psi(a, x).conj() * psi(b, x)).scale(cell);
                            for dir in 0..2 {
                                stiff[a][b] += (cov(a, dir).conj() * cov(b, dir)).scale(cell);
                            }
                        }
                    }
                }
            }
        }
        for a in 0..3 {
            for b in 0..3 {
                assert!((mass[a][b] - m_local[a][b]).abs() < 2e-4 * area, "M {a}{b}");
                assert!((stiff[a][b] - a_local[a][b]).abs() < 2e-3, "A {a}{b}: {:?} {:?}", stiff[a][b], a_local[a][b]);
            }
        }
    }
}
