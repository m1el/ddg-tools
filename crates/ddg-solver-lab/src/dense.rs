//! Dense active-set method for
//!
//! ```text
//! minimize ½ µᵀ H µ + qᵀ µ + λ ‖µ‖₁     (optionally subject to Σ µ = c)
//! ```
//!
//! over a growing set of variables, warm started. `H` is positive semidefinite. With the
//! constraint, variable 0 is the *pin*: its row and column of `H` and its `q` are zero (it does
//! not influence `u`), and it absorbs the sum.
//!
//! The method is the classical one for ℓ1 problems (Osborne et al. 2000): on the current
//! support with fixed signs solve the stationarity system with a Cholesky factor of `H_BB`
//! (maintained by appending and deleting rows in O(a²)), walk from the current point towards
//! that solution and stop where a component reaches zero (it leaves), and when the solution is
//! reached add the most violated variable. The objective decreases at every step.

#[derive(Debug, Default, Clone)]
pub struct DenseStats {
    pub iterations: usize,
    pub adds: usize,
    pub drops: usize,
    pub failed_adds: usize,
    pub peak_active: usize,
}

pub struct Dense {
    pub lambda: f64,
    pub sum: Option<f64>,
    dim: usize,
    cap: usize,
    h: Vec<f64>, // cap × cap column major
    q: Vec<f64>,
    pub mu: Vec<f64>,
    sign: Vec<f64>,
    blocked: Vec<bool>,
    act: Vec<usize>, // B, in Cholesky order (pin excluded)
    chol: Vec<f64>,  // row-major lower triangle, leading dimension chol_ld
    chol_ld: usize,
    pin_active: bool,
    pub nu: f64,
    pub stats: DenseStats,
}

impl Dense {
    pub fn new(lambda: f64, sum: Option<f64>) -> Dense {
        let mut d = Dense {
            lambda,
            sum,
            dim: 0,
            cap: 0,
            h: vec![],
            q: vec![],
            mu: vec![],
            sign: vec![],
            blocked: vec![],
            act: vec![],
            chol: vec![],
            chol_ld: 0,
            pin_active: false,
            nu: 0.0,
            stats: DenseStats::default(),
        };
        if let Some(c) = sum {
            d.extend(&[0.0; 0], 1, &[0.0]);
            d.mu[0] = c;
            d.sign[0] = if c >= 0.0 { 1.0 } else { -1.0 };
            d.pin_active = true;
        }
        d
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Number of nonzero variables (pin included).
    pub fn n_active(&self) -> usize {
        self.act.len() + self.pin_active as usize
    }

    /// Indices (in W) of the active variables.
    pub fn active(&self) -> Vec<usize> {
        let mut v: Vec<usize> = self.act.clone();
        if self.pin_active {
            v.push(0);
        }
        v
    }

    pub fn is_active(&self, j: usize) -> bool {
        self.sign[j] != 0.0
    }

    /// Appends `k` variables. `cross` is `(dim + k) × k` column major: column `t` holds
    /// `H[·, dim + t]` for all rows of the enlarged matrix (rows `dim..` are the new block).
    pub fn extend(&mut self, cross: &[f64], k: usize, q_new: &[f64]) {
        let old = self.dim;
        let new = old + k;
        if new > self.cap {
            let cap = (new * 3 / 2).max(64);
            let mut h = vec![0.0; cap * cap];
            for j in 0..old {
                h[j * cap..j * cap + old].copy_from_slice(&self.h[j * self.cap..j * self.cap + old]);
            }
            self.h = h;
            self.cap = cap;
        }
        for t in 0..k {
            for i in 0..new {
                let v = if cross.is_empty() { 0.0 } else { cross[t * new + i] };
                self.h[(old + t) * self.cap + i] = v;
                self.h[i * self.cap + old + t] = v;
            }
        }
        self.q.extend_from_slice(q_new);
        self.mu.extend(std::iter::repeat_n(0.0, k));
        self.sign.extend(std::iter::repeat_n(0.0, k));
        self.blocked.extend(std::iter::repeat_n(false, k));
        self.dim = new;
    }

    #[inline]
    fn hh(&self, i: usize, j: usize) -> f64 {
        self.h[j * self.cap + i]
    }

    fn chol_append(&mut self, j: usize) -> bool {
        let a = self.act.len();
        if a + 1 > self.chol_ld {
            let ld = ((a + 1) * 3 / 2).max(64);
            let mut c = vec![0.0; ld * ld];
            for i in 0..a {
                c[i * ld..i * ld + i + 1].copy_from_slice(&self.chol[i * self.chol_ld..i * self.chol_ld + i + 1]);
            }
            self.chol = c;
            self.chol_ld = ld;
        }
        let ld = self.chol_ld;
        let mut y = vec![0.0; a];
        for i in 0..a {
            let mut s = self.hh(self.act[i], j);
            for t in 0..i {
                s -= self.chol[i * ld + t] * y[t];
            }
            y[i] = s / self.chol[i * ld + i];
        }
        let d2 = self.hh(j, j) - y.iter().map(|v| v * v).sum::<f64>();
        if !(d2 > 1e-13 * self.hh(j, j).abs()) || !d2.is_finite() {
            return false;
        }
        self.chol[a * ld..a * ld + a].copy_from_slice(&y);
        self.chol[a * ld + a] = d2.sqrt();
        self.act.push(j);
        true
    }

    fn chol_delete(&mut self, k: usize) {
        let a = self.act.len();
        let ld = self.chol_ld;
        // remove row k
        for i in k..a - 1 {
            for t in 0..=i + 1 {
                self.chol[i * ld + t] = self.chol[(i + 1) * ld + t];
            }
        }
        // rows are now 0..a-1; row i >= k has entries in columns 0..=i+1
        for i in k..a - 1 {
            let x = self.chol[i * ld + i];
            let y = self.chol[i * ld + i + 1];
            let r = x.hypot(y);
            let (c, s) = (x / r, y / r);
            for t in i..a - 1 {
                let xt = self.chol[t * ld + i];
                let yt = self.chol[t * ld + i + 1];
                self.chol[t * ld + i] = c * xt + s * yt;
                self.chol[t * ld + i + 1] = -s * xt + c * yt;
            }
            self.chol[i * ld + i + 1] = 0.0;
        }
        // clear the stale last row/col
        for t in 0..a {
            self.chol[(a - 1) * ld + t] = 0.0;
        }
        for t in 0..a - 1 {
            self.chol[t * ld + a - 1] = 0.0;
        }
        self.act.remove(k);
    }

    fn chol_solve(&self, b: &[f64]) -> Vec<f64> {
        let a = self.act.len();
        let ld = self.chol_ld;
        let mut x = b.to_vec();
        for i in 0..a {
            let mut s = x[i];
            for t in 0..i {
                s -= self.chol[i * ld + t] * x[t];
            }
            x[i] = s / self.chol[i * ld + i];
        }
        for i in (0..a).rev() {
            let mut s = x[i];
            for t in i + 1..a {
                s -= self.chol[t * ld + i] * x[t];
            }
            x[i] = s / self.chol[i * ld + i];
        }
        x
    }

    fn deactivate(&mut self, j: usize) {
        if j == 0 && self.sum.is_some() {
            self.pin_active = false;
        } else {
            let pos = self.act.iter().position(|&x| x == j).unwrap();
            self.chol_delete(pos);
        }
        self.sign[j] = 0.0;
        self.mu[j] = 0.0;
    }

    /// Runs the active-set iteration at the current `lambda` until no variable of the working
    /// set violates its optimality condition by more than `tol · λ`.
    pub fn solve(&mut self, tol: f64, max_iter: usize) -> bool {
        let lam = self.lambda;
        for _ in 0..max_iter {
            self.stats.iterations += 1;
            let b: Vec<usize> = self.act.clone();
            let rhs: Vec<f64> = b.iter().map(|&j| -self.q[j] - lam * self.sign[j]).collect();
            let mut star = vec![0.0; self.dim];
            let nu;
            match self.sum {
                None => {
                    let x = self.chol_solve(&rhs);
                    for (t, &j) in b.iter().enumerate() {
                        star[j] = x[t];
                    }
                    nu = 0.0;
                }
                Some(c) => {
                    if self.pin_active {
                        nu = -lam * self.sign[0];
                        let r2: Vec<f64> = rhs.iter().map(|v| v - nu).collect();
                        let x = self.chol_solve(&r2);
                        let mut sum = 0.0;
                        for (t, &j) in b.iter().enumerate() {
                            star[j] = x[t];
                            sum += x[t];
                        }
                        star[0] = c - sum;
                    } else {
                        let x = self.chol_solve(&rhs);
                        let y = self.chol_solve(&vec![1.0; b.len()]);
                        let sx: f64 = x.iter().sum();
                        let sy: f64 = y.iter().sum();
                        nu = (sx - c) / sy;
                        for (t, &j) in b.iter().enumerate() {
                            star[j] = x[t] - nu * y[t];
                        }
                    }
                }
            }
            // line search
            let mut t = 1.0f64;
            let mut drop = None;
            let mut actl = b.clone();
            if self.pin_active {
                actl.push(0);
            }
            for &i in &actl {
                if self.sign[i] * star[i] < 0.0 {
                    let ti = self.mu[i] / (self.mu[i] - star[i]);
                    if ti < t {
                        t = ti;
                        drop = Some(i);
                    }
                }
            }
            if let Some(d) = drop {
                for &i in &actl {
                    self.mu[i] += t * (star[i] - self.mu[i]);
                }
                self.deactivate(d);
                self.stats.drops += 1;
                if t < 1e-12 {
                    self.blocked[d] = true;
                } else {
                    self.blocked.iter_mut().for_each(|x| *x = false);
                }
                continue;
            }
            for &i in &actl {
                self.mu[i] = star[i];
            }
            self.nu = nu;
            // violations
            let mut best: Option<(usize, f64)> = None;
            let mut g = self.q.clone();
            for &i in &b {
                let m = self.mu[i];
                if m != 0.0 {
                    for j in 0..self.dim {
                        g[j] += self.hh(j, i) * m;
                    }
                }
            }
            for j in 0..self.dim {
                if self.sign[j] != 0.0 || self.blocked[j] {
                    continue;
                }
                let ex = (g[j] + nu).abs() - lam;
                if ex > tol * lam && best.is_none_or(|(_, e)| ex > e) {
                    best = Some((j, ex));
                }
            }
            match best {
                None => return true,
                Some((j, _)) => {
                    let s = -(g[j] + nu).signum();
                    if self.sum.is_some() && j == 0 {
                        self.pin_active = true;
                        self.sign[0] = s;
                        self.stats.adds += 1;
                    } else if self.chol_append(j) {
                        self.sign[j] = s;
                        self.stats.adds += 1;
                    } else {
                        self.blocked[j] = true;
                        self.stats.failed_adds += 1;
                    }
                    self.stats.peak_active = self.stats.peak_active.max(self.n_active());
                }
            }
        }
        false
    }
}
