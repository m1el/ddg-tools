//! Oriented point clouds and their validation (SPEC.md §2.1).

use crate::PatError;

/// How [`PointCloud::new`] treats normals that are not of unit length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NormalPolicy {
    /// Reject a normal whose length differs from 1 by more than `tol` (default `1e-6`).
    RequireUnit {
        /// Allowed deviation of the length from 1.
        tol: f64,
    },
    /// Divide every normal by its length (explicit opt-in; counted in the report).
    Normalize,
}

impl Default for NormalPolicy {
    fn default() -> Self {
        NormalPolicy::RequireUnit { tol: 1e-6 }
    }
}

/// What validation found.
#[derive(Debug, Clone, PartialEq)]
pub struct CloudReport {
    /// Number of points.
    pub n: usize,
    /// Normals changed by [`NormalPolicy::Normalize`] (bitwise different after division).
    pub normalized: usize,
    /// Points whose position equals that of a point with a smaller index.
    pub exact_duplicates: usize,
    /// Axis-aligned bounding box `[min, max]`.
    pub bbox: [[f64; 3]; 2],
}

/// A validated oriented point cloud: finite positions, finite non-zero (unit) normals, `n ≥ 2`.
#[derive(Debug, Clone, PartialEq)]
pub struct PointCloud {
    positions: Vec<[f64; 3]>,
    normals: Vec<[f64; 3]>,
}

fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

impl PointCloud {
    /// Validate and build. Errors name the first offending point.
    pub fn new(
        positions: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        policy: NormalPolicy,
    ) -> Result<(Self, CloudReport), PatError> {
        if positions.len() != normals.len() {
            return Err(PatError::LengthMismatch { positions: positions.len(), normals: normals.len() });
        }
        let n = positions.len();
        if n < 2 {
            return Err(PatError::TooFewPoints { n, need: 2 });
        }
        if let NormalPolicy::RequireUnit { tol } = policy
            && !(tol.is_finite() && tol >= 0.0)
        {
            return Err(PatError::InvalidOption { what: format!("NormalPolicy::RequireUnit tol = {tol}") });
        }
        let mut normals = normals;
        let mut normalized = 0;
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for (i, (p, nr)) in positions.iter().zip(normals.iter_mut()).enumerate() {
            if !p.iter().all(|c| c.is_finite()) {
                return Err(PatError::NonFinite { index: i, what: "position" });
            }
            if !nr.iter().all(|c| c.is_finite()) {
                return Err(PatError::NonFinite { index: i, what: "normal" });
            }
            let len = norm3(*nr);
            if len == 0.0 {
                return Err(PatError::ZeroNormal { index: i });
            }
            match policy {
                NormalPolicy::RequireUnit { tol } => {
                    if (len - 1.0).abs() > tol {
                        return Err(PatError::NonUnitNormal { index: i, norm: len });
                    }
                }
                NormalPolicy::Normalize => {
                    let m = [nr[0] / len, nr[1] / len, nr[2] / len];
                    if m.iter().zip(nr.iter()).any(|(a, b)| a.to_bits() != b.to_bits()) {
                        normalized += 1;
                    }
                    *nr = m;
                }
            }
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        let exact_duplicates = count_duplicates(&positions);
        let report = CloudReport { n, normalized, exact_duplicates, bbox: [lo, hi] };
        Ok((PointCloud { positions, normals }, report))
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Always false (a cloud has at least two points).
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Positions.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// Normals (unit length within the policy's tolerance).
    pub fn normals(&self) -> &[[f64; 3]] {
        &self.normals
    }

    /// 64-bit FNV-1a hash of the positions and normals (bit patterns), used to tie saved tori to
    /// their cloud.
    pub fn hash(&self) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut eat = |x: f64| {
            for b in x.to_bits().to_le_bytes() {
                h ^= b as u64;
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        for p in &self.positions {
            p.iter().for_each(|&c| eat(c));
        }
        for p in &self.normals {
            p.iter().for_each(|&c| eat(c));
        }
        h
    }
}

fn key(p: &[f64; 3]) -> [u64; 3] {
    // +0.0 so that -0.0 and 0.0 compare equal
    [(p[0] + 0.0).to_bits(), (p[1] + 0.0).to_bits(), (p[2] + 0.0).to_bits()]
}

fn count_duplicates(p: &[[f64; 3]]) -> usize {
    let mut keys: Vec<[u64; 3]> = p.iter().map(key).collect();
    keys.sort_unstable();
    keys.windows(2).filter(|w| w[0] == w[1]).count()
}
