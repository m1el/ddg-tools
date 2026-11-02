//! The vertex map of a remeshing and the transport of per-vertex data through it
//! (`SPEC.md` §10 step 3, D17).

use crate::error::RsError;

/// New vertex → old vertices with weights. A vertex that was not touched maps to itself with
/// weight 1; a split midpoint to `½ a + ½ b` (of the endpoints' own entries); a collapse
/// survivor to the weights of its position rule (`½, ½` for the midpoint; `1, 0` when it sits
/// at an endpoint, the zero entry kept so that every merged old vertex is listed). The
/// weights of each entry sum to 1 (up to round-off).
///
/// Smoothing moves vertices tangentially without changing the map: per-vertex data follow the
/// vertex, not the point in space.
#[derive(Debug, Clone, PartialEq)]
pub struct VertexMap {
    offsets: Vec<u32>,
    sources: Vec<u32>,
    weights: Vec<f64>,
    n_old: usize,
}

impl VertexMap {
    pub(crate) fn from_parts(
        offsets: Vec<u32>,
        sources: Vec<u32>,
        weights: Vec<f64>,
        n_old: usize,
    ) -> Self {
        VertexMap {
            offsets,
            sources,
            weights,
            n_old,
        }
    }

    /// The identity map on `n` vertices.
    pub fn identity(n: usize) -> Self {
        VertexMap {
            offsets: (0..=n as u32).collect(),
            sources: (0..n as u32).collect(),
            weights: vec![1.0; n],
            n_old: n,
        }
    }

    /// Number of new vertices.
    pub fn n_new(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Number of old vertices (the length per-vertex input data must have).
    pub fn n_old(&self) -> usize {
        self.n_old
    }

    /// The old vertices of new vertex `v` and their weights.
    pub fn entry(&self, v: u32) -> (&[u32], &[f64]) {
        let (lo, hi) = (
            self.offsets[v as usize] as usize,
            self.offsets[v as usize + 1] as usize,
        );
        (&self.sources[lo..hi], &self.weights[lo..hi])
    }

    fn check_len(&self, len: usize) -> Result<(), RsError> {
        if len != self.n_old {
            return Err(RsError::PositionCountMismatch {
                expected: self.n_old,
                found: len,
            });
        }
        Ok(())
    }

    /// Transports a scalar per-vertex field: `new[v] = Σ w_i old[i]` (split and collapse
    /// midpoints interpolate linearly).
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`] if `old` does not have one value per old vertex.
    pub fn transport(&self, old: &[f64]) -> Result<Vec<f64>, RsError> {
        self.check_len(old.len())?;
        Ok((0..self.n_new() as u32)
            .map(|v| {
                let (s, w) = self.entry(v);
                s.iter().zip(w).map(|(&i, &w)| w * old[i as usize]).sum()
            })
            .collect())
    }

    /// Transports a 3-vector per-vertex field (e.g. squared-error targets, `SPEC.md` §5.2,
    /// D17), componentwise as [`transport`](Self::transport).
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`] if `old` does not have one value per old vertex.
    pub fn transport_vec3(&self, old: &[[f64; 3]]) -> Result<Vec<[f64; 3]>, RsError> {
        self.check_len(old.len())?;
        Ok((0..self.n_new() as u32)
            .map(|v| {
                let (s, w) = self.entry(v);
                let mut r = [0.0; 3];
                for (&i, &w) in s.iter().zip(w) {
                    for k in 0..3 {
                        r[k] += w * old[i as usize][k];
                    }
                }
                r
            })
            .collect())
    }

    /// Transports a vertex set (e.g. attractor membership, `SPEC.md` §10 step 3): new vertex
    /// `v` is in the set iff the weights of its old vertices that are in the set sum to at
    /// least `threshold` (`½`: a split midpoint joins if both ends are members, a collapse
    /// survivor if the point it sits at is). Returns the new set, ascending.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a member that is not an old vertex.
    pub fn transport_set(&self, old: &[u32], threshold: f64) -> Result<Vec<u32>, RsError> {
        let mut member = vec![0.0; self.n_old];
        for &i in old {
            if i as usize >= self.n_old {
                return Err(RsError::InvalidOption {
                    name: "vertex set",
                    value: i.to_string(),
                    requirement: "every member must be an old vertex",
                });
            }
            member[i as usize] = 1.0;
        }
        let m = self.transport(&member)?;
        Ok((0..m.len() as u32)
            .filter(|&v| m[v as usize] >= threshold)
            .collect())
    }

    /// The composition `self ∘ earlier`: the map from the vertices before `earlier` to the
    /// vertices after `self` (for several passes in a row).
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`] if `earlier` does not end where `self` starts.
    pub fn compose(&self, earlier: &VertexMap) -> Result<VertexMap, RsError> {
        self.check_len(earlier.n_new())?;
        let mut offsets = vec![0u32];
        let mut sources: Vec<u32> = Vec::new();
        let mut weights: Vec<f64> = Vec::new();
        for v in 0..self.n_new() as u32 {
            let start = sources.len();
            let (s, w) = self.entry(v);
            for (&mid, &w1) in s.iter().zip(w) {
                let (s2, w2) = earlier.entry(mid);
                for (&o, &w2) in s2.iter().zip(w2) {
                    if let Some(k) = sources[start..].iter().position(|&x| x == o) {
                        weights[start + k] += w1 * w2;
                    } else {
                        sources.push(o);
                        weights.push(w1 * w2);
                    }
                }
            }
            offsets.push(sources.len() as u32);
        }
        Ok(VertexMap {
            offsets,
            sources,
            weights,
            n_old: earlier.n_old,
        })
    }
}
