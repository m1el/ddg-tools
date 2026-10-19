//! Tests of the symmetric eigensolvers against a double-double reference and known spectra.
#![allow(clippy::needless_range_loop)]
use ddg_developability::eig::{EigError, sym_eig2, sym_eig3};

const EPS: f64 = f64::EPSILON;

// ---------- small seeded generator (splitmix64) ----------
struct Rng(u64);
impl Rng {
    fn u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    /// uniform in [-1, 1)
    fn unif(&mut self) -> f64 {
        ((self.u64() >> 11) as f64) / ((1u64 << 52) as f64) - 1.0
    }
    fn unit3(&mut self) -> [f64; 3] {
        loop {
            let v = [self.unif(), self.unif(), self.unif()];
            let n = norm3(v);
            if n > 0.1 && n <= 1.0 {
                return [v[0] / n, v[1] / n, v[2] / n];
            }
        }
    }
}

fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
fn frob(a: &[[f64; 3]; 3]) -> f64 {
    a.iter().flatten().map(|x| x * x).sum::<f64>().sqrt()
}
fn matvec(a: &[[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    let mut r = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i] += a[i][j] * x[j];
        }
    }
    r
}

/// Check residuals, orthonormality, handedness, ordering. Returns (max residual / (eps*norm), max orth / eps).
fn check3(a: [[f64; 3]; 3]) -> (f64, f64) {
    let e = sym_eig3(a).unwrap();
    let nrm = frob(&a);
    let mut res: f64 = 0.0;
    let mut orth: f64 = 0.0;
    for k in 0..3 {
        let x = e.vectors[k];
        let ax = matvec(&a, x);
        let r = [ax[0] - e.values[k] * x[0], ax[1] - e.values[k] * x[1], ax[2] - e.values[k] * x[2]];
        res = res.max(norm3(r));
        for l in 0..3 {
            let y = e.vectors[l];
            let d = x[0] * y[0] + x[1] * y[1] + x[2] * y[2] - if k == l { 1.0 } else { 0.0 };
            orth = orth.max(d.abs());
        }
    }
    assert!(e.values[0] <= e.values[1] && e.values[1] <= e.values[2], "not ascending {:?}", e.values);
    let (u, v, w) = (e.vectors[0], e.vectors[1], e.vectors[2]);
    let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    let hd = norm3([cr[0] - w[0], cr[1] - w[1], cr[2] - w[2]]);
    assert!(hd < 1e-14, "not right-handed: {hd}");
    let scale = if nrm > 0.0 { EPS * nrm } else { 1.0 };
    (res / scale, orth / EPS)
}

// ---------- double-double reference ----------
#[derive(Clone, Copy)]
struct Dd(f64, f64);
fn two_sum(a: f64, b: f64) -> Dd {
    let s = a + b;
    let bb = s - a;
    Dd(s, (a - (s - bb)) + (b - bb))
}
fn quick(a: f64, b: f64) -> Dd {
    let s = a + b;
    Dd(s, b - (s - a))
}
impl Dd {
    fn add(self, o: Dd) -> Dd {
        let s = two_sum(self.0, o.0);
        let t = two_sum(self.1, o.1);
        let u = quick(s.0, s.1 + t.0);
        quick(u.0, u.1 + t.1)
    }
    fn neg(self) -> Dd {
        Dd(-self.0, -self.1)
    }
    fn sub(self, o: Dd) -> Dd {
        self.add(o.neg())
    }
    fn mul(self, o: Dd) -> Dd {
        let p = self.0 * o.0;
        let e = self.0.mul_add(o.0, -p);
        quick(p, e + (self.0 * o.1 + self.1 * o.0))
    }
    fn div(self, o: Dd) -> Dd {
        let q1 = self.0 / o.0;
        let r = self.sub(o.mul(Dd(q1, 0.0)));
        let q2 = r.0 / o.0;
        let r = r.sub(o.mul(Dd(q2, 0.0)));
        let q3 = r.0 / o.0;
        quick(q1, q2).add(Dd(q3, 0.0))
    }
    fn sqrt(self) -> Dd {
        let s = self.0.sqrt();
        let s0 = Dd(s, 0.0);
        let d = self.sub(s0.mul(s0));
        s0.add(Dd(d.0 / (2.0 * s), 0.0))
    }
    fn abs(self) -> Dd {
        if self.0 < 0.0 { self.neg() } else { self }
    }
    fn f(self) -> f64 {
        self.0 + self.1
    }
}

/// Jacobi in double-double: (values, vectors as columns v[i][j]) unsorted.
fn dd_eig3(a: [[f64; 3]; 3]) -> ([Dd; 3], [[Dd; 3]; 3]) {
    let z = Dd(0.0, 0.0);
    let o = Dd(1.0, 0.0);
    let mut m = [[z; 3]; 3];
    let mut v = [[z; 3]; 3];
    for i in 0..3 {
        v[i][i] = o;
        for j in 0..3 {
            m[i][j] = Dd(a[i.min(j)][i.max(j)], 0.0);
        }
    }
    for _ in 0..30 {
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if m[p][q].0.abs() <= 1e-34 * m[p][p].0.abs().max(m[q][q].0.abs()) {
                m[p][q] = z;
                m[q][p] = z;
                continue;
            }
            let theta = m[q][q].sub(m[p][p]).div(m[p][q].add(m[p][q]));
            let sg = if theta.0 >= 0.0 { o } else { o.neg() };
            let t = sg.div(theta.abs().add(theta.mul(theta).add(o).sqrt()));
            let c = o.div(t.mul(t).add(o).sqrt());
            let s = t.mul(c);
            // J^T M J with J = [[c, s], [-s, c]] in the (p,q) plane
            let mut j = [[z; 3]; 3];
            for i in 0..3 {
                j[i][i] = o;
            }
            j[p][p] = c;
            j[q][q] = c;
            j[p][q] = s;
            j[q][p] = s.neg();
            let mul = |x: &[[Dd; 3]; 3], y: &[[Dd; 3]; 3]| {
                let mut r = [[z; 3]; 3];
                for i in 0..3 {
                    for k in 0..3 {
                        for l in 0..3 {
                            r[i][k] = r[i][k].add(x[i][l].mul(y[l][k]));
                        }
                    }
                }
                r
            };
            let mut jt = [[z; 3]; 3];
            for i in 0..3 {
                for k in 0..3 {
                    jt[i][k] = j[k][i];
                }
            }
            m = mul(&jt, &mul(&m, &j));
            m[p][q] = z;
            m[q][p] = z;
            v = mul(&v, &j);
        }
    }
    ([m[0][0], m[1][1], m[2][2]], v)
}

/// Sorted reference eigenvalues and eigenvectors (as rows).
fn reference(a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let (d, v) = dd_eig3(a);
    let mut idx = [0, 1, 2];
    idx.sort_by(|&i, &j| d[i].f().partial_cmp(&d[j].f()).unwrap());
    let mut vals = [0.0; 3];
    let mut vecs = [[0.0; 3]; 3];
    for (k, &i) in idx.iter().enumerate() {
        vals[k] = d[i].f();
        for r in 0..3 {
            vecs[k][r] = v[r][i].f();
        }
    }
    (vals, vecs)
}

fn angle_err(x: [f64; 3], y: [f64; 3]) -> f64 {
    let d = x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let c = [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]];
    norm3(c).atan2(d.abs())
}

// ---------- matrix builders ----------
fn from_spectrum(q: [[f64; 3]; 3], d: [f64; 3]) -> [[f64; 3]; 3] {
    // A = sum d_k q_k q_k^T, q_k = q[k]
    let mut a = [[0.0; 3]; 3];
    for k in 0..3 {
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += d[k] * q[k][i] * q[k][j];
            }
        }
    }
    a
}
fn random_basis(r: &mut Rng) -> [[f64; 3]; 3] {
    let u = r.unit3();
    let mut w = r.unit3();
    let dot = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
    for i in 0..3 {
        w[i] -= dot * u[i];
    }
    let n = norm3(w);
    let v = [w[0] / n, w[1] / n, w[2] / n];
    let c = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    [u, v, c]
}
/// Sum theta_k N_k N_k^T, normals nearly in a plane (out-of-plane perturbation `pert`).
fn planar_normals_matrix(r: &mut Rng, n: usize, pert: f64) -> [[f64; 3]; 3] {
    let basis = random_basis(r);
    let mut a = [[0.0; 3]; 3];
    let mut ws: Vec<f64> = (0..n).map(|_| 0.05 + (r.unif() + 1.0)).collect();
    let tot: f64 = ws.iter().sum();
    for w in &mut ws {
        *w *= 2.0 * std::f64::consts::PI / tot;
    }
    for w in ws {
        let phi = std::f64::consts::PI * r.unif();
        let z = pert * r.unif();
        let loc = [phi.cos(), phi.sin(), z];
        let l = norm3(loc);
        let mut nn = [0.0; 3];
        for i in 0..3 {
            nn[i] = (loc[0] * basis[0][i] + loc[1] * basis[1][i] + loc[2] * basis[2][i]) / l;
        }
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += w * nn[i] * nn[j];
            }
        }
    }
    a
}

#[test]
fn random_symmetric() {
    let mut r = Rng(1);
    let (mut wr, mut wo) = (0.0f64, 0.0f64);
    for i in 0..20000 {
        let s = [1.0, 1e-3, 1e5][i % 3];
        let mut a = [[0.0; 3]; 3];
        for p in 0..3 {
            for q in p..3 {
                a[p][q] = s * r.unif();
                a[q][p] = a[p][q];
            }
        }
        let (res, orth) = check3(a);
        wr = wr.max(res);
        wo = wo.max(orth);
    }
    println!("random: worst residual {wr:.2} eps*||A||, orthonormality {wo:.2} eps");
    assert!(wr < 10.0 && wo < 10.0);
}

#[test]
fn known_spectrum_tiny_lambda_min() {
    // Q diag Qt is rounded to f64, so compare against the dd eigenvalues of the rounded matrix.
    let mut r = Rng(2);
    let mut worst: f64 = 0.0;
    for e in 0..=18 {
        for _ in 0..200 {
            let lam_min = if e == 0 { 0.0 } else { 10f64.powi(-e) };
            let q = random_basis(&mut r);
            let d = [lam_min, 1.0 + r.unif().abs(), 5.0 + r.unif().abs()];
            let a = from_spectrum(q, d);
            let (rv, _) = reference(a);
            let got = sym_eig3(a).unwrap();
            let err = (got.values[0] - rv[0]).abs() / (EPS * frob(&a));
            worst = worst.max(err);
            for k in 1..3 {
                assert!((got.values[k] - rv[k]).abs() <= 8.0 * EPS * frob(&a));
            }
            check3(a);
        }
    }
    println!("tiny lambda_min: worst abs error {worst:.3} eps*||A||");
    assert!(worst < 4.0);
}

#[test]
fn angle_weighted_normals() {
    let mut r = Rng(3);
    let mut worst_val: f64 = 0.0;
    let mut worst_vec: f64 = 0.0; // angle * gap / (eps ||A||)
    let mut worst_res: f64 = 0.0;
    for e in 2..=9 {
        let pert = 10f64.powi(-e);
        for n in 3..=12 {
            for _ in 0..100 {
                let a = planar_normals_matrix(&mut r, n, pert);
                let (rv, rvec) = reference(a);
                let got = sym_eig3(a).unwrap();
                let nrm = frob(&a);
                worst_val = worst_val.max((got.values[0] - rv[0]).abs() / (EPS * nrm));
                let (res, _) = check3(a);
                worst_res = worst_res.max(res);
                let gap = rv[1] - rv[0];
                let ang = angle_err(got.vectors[0], rvec[0]);
                // Davis-Kahan: angle <= c eps ||A|| / gap
                worst_vec = worst_vec.max(ang * gap / (EPS * nrm));
            }
        }
    }
    println!(
        "normals: lambda_min abs err {worst_val:.3} eps*||A||, residual {worst_res:.2}, eigvec angle*gap {worst_vec:.3} eps*||A||"
    );
    assert!(worst_val < 4.0 && worst_res < 10.0 && worst_vec < 20.0);
}

#[test]
fn near_degenerate_pair_eigenvector() {
    let mut r = Rng(4);
    let mut worst: f64 = 0.0;
    for e in 2..=14 {
        for _ in 0..200 {
            let gap = 10f64.powi(-e);
            let q = random_basis(&mut r);
            let a = from_spectrum(q, [1e-3, 1e-3 + gap, 3.0]);
            let (rv, rvec) = reference(a);
            let got = sym_eig3(a).unwrap();
            let g = rv[1] - rv[0];
            let ang = angle_err(got.vectors[0], rvec[0]);
            worst = worst.max(ang * g / (EPS * frob(&a)));
            assert!(ang <= std::f64::consts::FRAC_PI_2);
        }
    }
    println!("near-degenerate: angle*gap {worst:.3} eps*||A||");
    assert!(worst < 20.0);
}

#[test]
fn exact_cases() {
    // zero
    let e = sym_eig3([[0.0; 3]; 3]).unwrap();
    assert_eq!(e.values, [0.0; 3]);
    check3([[0.0; 3]; 3]);
    // identity (and scaled)
    let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(sym_eig3(id).unwrap().values, [1.0; 3]);
    check3(id);
    // diagonal, unsorted
    let e = sym_eig3([[3.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 2.0]]).unwrap();
    assert_eq!(e.values, [-1.0, 2.0, 3.0]);
    check3([[3.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 2.0]]);
    // rank 1 (integer entries, exact): u u^T with u = (1,2,3)
    let u = [1.0, 2.0, 3.0];
    let mut a = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = u[i] * u[j];
        }
    }
    let e = sym_eig3(a).unwrap();
    assert!(e.values[0].abs() <= 4.0 * EPS * 14.0 && e.values[1].abs() <= 4.0 * EPS * 14.0);
    assert!((e.values[2] - 14.0).abs() <= 4.0 * EPS * 14.0);
    assert!(angle_err(e.vectors[2], [1.0 / 14f64.sqrt(), 2.0 / 14f64.sqrt(), 3.0 / 14f64.sqrt()]) < 1e-15);
    check3(a);
    // rank 2 with exact zero: M^T M, M integer 2x3
    let m = [[1.0, 2.0, 0.0], [3.0, -1.0, 4.0]];
    let mut a = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = m[0][i] * m[0][j] + m[1][i] * m[1][j];
        }
    }
    let e = sym_eig3(a).unwrap();
    assert!(e.values[0].abs() <= 4.0 * EPS * frob(&a), "{:e}", e.values[0]);
    // null vector is m0 x m1
    let nv = [2.0 * 4.0, -(1.0 * 4.0), -1.0 - 2.0 * 3.0];
    let n = norm3(nv);
    assert!(angle_err(e.vectors[0], [nv[0] / n, nv[1] / n, nv[2] / n]) < 1e-14);
    check3(a);
    // two equal eigenvalues and three equal (rotated)
    let mut r = Rng(5);
    for _ in 0..200 {
        let q = random_basis(&mut r);
        for d in [[1.0, 1.0, 4.0], [0.5, 4.0, 4.0], [2.0, 2.0, 2.0], [0.0, 0.0, 3.0], [0.0, 3.0, 3.0]] {
            let a = from_spectrum(q, d);
            let e = sym_eig3(a).unwrap();
            for k in 0..3 {
                assert!((e.values[k] - d[k]).abs() <= 8.0 * EPS * 4.0);
            }
            let (res, orth) = check3(a);
            assert!(res < 10.0 && orth < 10.0);
        }
    }
}

#[test]
fn extreme_scales_and_bad_input() {
    let mut r = Rng(6);
    for &s in &[1e-300, 1e-310, 1e300, 1e307, 1e-160, 1e160] {
        for _ in 0..200 {
            let q = random_basis(&mut r);
            let a0 = from_spectrum(q, [1e-6, 1.0, 2.0]);
            let mut a = a0;
            for row in &mut a {
                for x in row {
                    *x *= s;
                }
            }
            let e = sym_eig3(a).unwrap();
            let e0 = sym_eig3(a0).unwrap();
            for k in 0..3 {
                assert!(e.values[k].is_finite());
                // subnormal inputs lose bits by themselves: allow absolute slack of that size
                let tol = 8.0 * EPS * 2.0 * s + if s < 1e-300 { 1e-322 } else { 0.0 };
                assert!((e.values[k] - e0.values[k] * s).abs() <= tol, "s={s:e} k={k}");
            }
            for k in 0..3 {
                let n = norm3(e.vectors[k]);
                assert!((n - 1.0).abs() < 1e-14);
            }
        }
    }
    let nan = f64::NAN;
    assert_eq!(sym_eig3([[nan, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), Err(EigError::NonFinite));
    assert_eq!(sym_eig3([[1.0, f64::INFINITY, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), Err(EigError::NonFinite));
    // lower triangle is ignored
    let a = [[2.0, 1.0, 0.0], [nan, 2.0, 0.0], [nan, nan, 5.0]];
    let e = sym_eig3(a).unwrap();
    assert!((e.values[0] - 1.0).abs() < 1e-15 && (e.values[1] - 3.0).abs() < 1e-15);
    assert_eq!(sym_eig2([[1.0, nan], [nan, 1.0]]), Err(EigError::NonFinite));
}

#[test]
fn eig2() {
    let mut r = Rng(7);
    for i in 0..20000 {
        let s = [1.0, 1e-200, 1e200][i % 3];
        let (a, b, c) = (r.unif() * s, r.unif() * s, r.unif() * s);
        let m = [[a, b], [b, c]];
        let e = sym_eig2(m).unwrap();
        let nrm = (a * a + 2.0 * b * b + c * c).sqrt();
        assert!(e.values[0] <= e.values[1]);
        for k in 0..2 {
            let x = e.vectors[k];
            let ax = [a * x[0] + b * x[1], b * x[0] + c * x[1]];
            let res = ((ax[0] - e.values[k] * x[0]).powi(2) + (ax[1] - e.values[k] * x[1]).powi(2)).sqrt();
            // squares of 1e-200 stay finite; 1e200 squares are only used via sqrt of scaled sums below
            assert!(res <= 4.0 * EPS * nrm, "res {res:e} nrm {nrm:e}");
            assert!((x[0] * x[0] + x[1] * x[1] - 1.0).abs() < 4.0 * EPS);
        }
        assert!((e.vectors[0][0] * e.vectors[1][1] - e.vectors[0][1] * e.vectors[1][0] - 1.0).abs() < 4.0 * EPS);
    }
    // exact cases
    assert_eq!(sym_eig2([[0.0; 2]; 2]).unwrap().values, [0.0, 0.0]);
    assert_eq!(sym_eig2([[2.0, 0.0], [0.0, 2.0]]).unwrap().values, [2.0, 2.0]);
    assert_eq!(sym_eig2([[3.0, 0.0], [0.0, 1.0]]).unwrap().values, [1.0, 3.0]);
    let e = sym_eig2([[1.0, 2.0], [2.0, 4.0]]).unwrap(); // rank 1
    assert!(e.values[0].abs() < 4.0 * EPS && (e.values[1] - 5.0).abs() < 8.0 * EPS);
    // nearly equal eigenvalues: a = 1, c = 1 + 2^-40, b = 1e-20: eigenvalues are accurate
    let (a, c, b) = (1.0, 1.0 + 2f64.powi(-40), 1e-20);
    let e = sym_eig2([[a, b], [b, c]]).unwrap();
    assert_eq!(e.values[0], 1.0);
    assert!((e.values[1] - c).abs() <= EPS);
    // tiny eigenvalue next to a big one: [[1, 1e-9],[1e-9, 1e-18]] -> det = 0 (exactly 1e-18 - 1e-18)
    let e = sym_eig2([[1.0, 1e-9], [1e-9, 1e-18]]).unwrap();
    assert!(e.values[0].abs() < 4.0 * EPS * 1e-18 + 1e-33, "{:e}", e.values[0]);
}
