//! The ICE error metric (paper Sec. 5; SPEC §4): channels of mass with error vectors, the
//! redistribution weights of Eq. 9, the mass and centre-of-mass update of Eq. 5 and 8 (or the
//! memoryless metric of Sec. 5.1), and the cost of Eq. 6 summed over all neighbours and
//! channels.
//!
//! The closed-form pieces are plain functions ([`redistribution_weights`], [`transport_update`],
//! [`memoryless_cost`]) so that the Lean oracle can check them; [`ChannelState`] holds the
//! per-vertex masses and vectors and evaluates or commits the cost of a removal from the
//! [`PreparedRemoval`] of the core (SPEC §3.6: the tentative and the committed removal give the
//! same `PreparedRemoval`, so the key is the commit cost for an unchanged neighbourhood, S5).

use std::f64::consts::PI;

use crate::error::IemError;
use crate::mesh::CoarseMesh;
use crate::removal::{PreparedRemoval, RemovalKind};

/// The tolerance below which a vertex counts as already flat (SPEC §4.2, S11): `|K_i| ≤ 1e-12`
/// (or `|κ_i|`), and its mass goes to its neighbours in equal shares.
pub const FLAT_TOLERANCE: f64 = 1e-12;

/// Where the mass of a channel comes from (paper Sec. 5.4, 5.4.1, 8.6.1; SPEC §4.1).
#[derive(Debug, Clone, PartialEq)]
pub enum ChannelSource {
    /// `K⁺ = max(K, 0)` at interior vertices, `κ⁺` at boundary vertices.
    CurvaturePositive,
    /// `K⁻ = −min(K, 0)` at interior vertices, `κ⁻` at boundary vertices.
    CurvatureNegative,
    /// The lumped vertex area `⅓ Σ A_f` divided by the total area (so the channel totals 1;
    /// as the reference, R2).
    Area,
    /// A value per fine vertex (`n_fine_vertices` entries, finite). A signed array is split into
    /// two channels, `max(x, 0)` and `−min(x, 0)` (Sec. 5.4.1), both with this channel's weight;
    /// an array without negative entries gives one channel.
    PerVertex(Vec<f64>),
}

/// A channel and its weight `w^p ≥ 0`.
#[derive(Debug, Clone, PartialEq)]
pub struct Channel {
    /// Where its mass comes from.
    pub source: ChannelSource,
    /// Its weight in the cost (Sec. 5.4.1), finite and `≥ 0`.
    pub weight: f64,
}

/// The channels of the metric and their weights (SPEC §4.1).
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelWeights {
    /// The channels. Channels of weight 0 are not tracked (they cannot change a cost).
    pub channels: Vec<Channel>,
}

impl Default for ChannelWeights {
    /// The default of S10: the two curvature channels with weight 1.
    fn default() -> Self {
        ChannelWeights::curvature()
    }
}

impl ChannelWeights {
    /// `K⁺` and `K⁻` with weight 1 (S10, the paper's ICE metric; the reference's `w = 0`).
    pub fn curvature() -> Self {
        ChannelWeights {
            channels: vec![
                Channel {
                    source: ChannelSource::CurvaturePositive,
                    weight: 1.0,
                },
                Channel {
                    source: ChannelSource::CurvatureNegative,
                    weight: 1.0,
                },
            ],
        }
    }

    /// The reference's single knob `w ∈ [0, 1]` (Fig. 11; R2): the curvature channels with
    /// weight `1 − w`, the area channel with weight `w`.
    pub fn mixture(w: f64) -> Self {
        ChannelWeights {
            channels: vec![
                Channel {
                    source: ChannelSource::CurvaturePositive,
                    weight: 1.0 - w,
                },
                Channel {
                    source: ChannelSource::CurvatureNegative,
                    weight: 1.0 - w,
                },
                Channel {
                    source: ChannelSource::Area,
                    weight: w,
                },
            ],
        }
    }

    /// These channels plus one more.
    pub fn with(mut self, source: ChannelSource, weight: f64) -> Self {
        self.channels.push(Channel { source, weight });
        self
    }
}

/// Which share each neighbour gets when the removed vertex is already flat (SPEC §4.2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FlatRule {
    /// `|K_i| ≤ tolerance` gives equal shares (S11; default [`FLAT_TOLERANCE`]).
    Threshold(f64),
    /// The reference (R12): equal shares only when `Σ_l |K̃_l − K_l|` is exactly zero, by adding
    /// `1e-6` to every change.
    ReferenceEpsilon,
}

/// How a boundary vertex's angular coordinates are normalized in the transport (SPEC S4, R3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryFrames {
    /// The crate's `φ = (π/Θ) θ` (IT-SPEC N41): the incoming boundary edge at `π`.
    Paper,
    /// The reference's `A = (2π/(Θ + π)) θ` (R3), computed from the crate's coordinates as
    /// `A = φ · 2Θ/(Θ + π)`.
    Reference,
}

/// Which angular coordinates the transport of a removal uses (SPEC §4.3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Frames {
    /// The frames of the flattened star at both ends (SPEC §4.3; the default).
    AfterFlattening,
    /// The frames before flattening at both ends (the reference's tentative cost, R16).
    BeforeFlattening,
    /// Before flattening at the removed vertex, after at the neighbour (the reference's commit,
    /// R16 and R19).
    ReferenceCommit,
}

/// How the cost of a removal is evaluated (SPEC §4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MetricOptions {
    /// The memory-based metric (Sec. 5.2, 5.3; default) or the memoryless one (Sec. 5.1:
    /// `t = 0` always, the cost is `Σ α_ij m_i ℓ_ij`).
    pub memory: bool,
    /// The rule for an already flat vertex.
    pub flat: FlatRule,
    /// Sum the cost only over the (neighbour, channel) pairs that receive mass (the reference,
    /// R11) instead of over all of `N_i` (paper Eq. 6, S11; default `false`).
    pub receiving_only: bool,
    /// The normalization of boundary frames.
    pub boundary_frames: BoundaryFrames,
}

impl Default for MetricOptions {
    fn default() -> Self {
        MetricOptions {
            memory: true,
            flat: FlatRule::Threshold(FLAT_TOLERANCE),
            receiving_only: false,
            boundary_frames: BoundaryFrames::Paper,
        }
    }
}

/// The redistribution weights of Eq. 9 (SPEC §4.2): `α_ij = |K̃_j − K_j| / Σ_l |K̃_l − K_l|`,
/// from the curvature `curvature` (`K_i` or `κ_i`) of the removed vertex and each neighbour's
/// curvature before and after flattening; equal shares `1/|N_i|` when `|K_i| ≤ 1e-12` (S11).
/// When `|K_i|` is larger the denominator is at least `|K_i|` (`iem-flatten-conserves-curvature`);
/// should rounding make it zero or not finite, the shares are equal too. Empty for no
/// neighbours.
///
/// [formula: iem-redistribution-weights]
pub fn redistribution_weights(curvature: f64, before: &[f64], after: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(before.len());
    let deltas = before.iter().zip(after).map(|(b, a)| (a - b).abs());
    weights_into(
        curvature,
        deltas,
        FlatRule::Threshold(FLAT_TOLERANCE),
        &mut out,
    );
    out
}

/// [`redistribution_weights`] into a buffer, from the changes `|K̃_j − K_j|`, with a flat rule.
///
/// [formula: iem-redistribution-weights]
fn weights_into(
    curvature: f64,
    deltas: impl Iterator<Item = f64>,
    rule: FlatRule,
    out: &mut Vec<f64>,
) {
    out.clear();
    out.extend(deltas);
    let n = out.len();
    if n == 0 {
        return;
    }
    let equal = |out: &mut Vec<f64>| out.iter_mut().for_each(|a| *a = 1.0 / n as f64);
    match rule {
        FlatRule::Threshold(tolerance) => {
            let sum: f64 = out.iter().sum();
            if curvature.abs() <= tolerance || !(sum > 0.0 && sum.is_finite()) {
                equal(out);
            } else {
                out.iter_mut().for_each(|a| *a /= sum);
            }
        }
        FlatRule::ReferenceEpsilon => {
            let mut sum: f64 = out.iter().sum();
            if sum == 0.0 {
                out.iter_mut().for_each(|a| *a += 1e-6);
                sum = out.iter().sum();
            }
            if !(sum > 0.0 && sum.is_finite()) {
                equal(out);
            } else {
                out.iter_mut().for_each(|a| *a /= sum);
            }
        }
    }
}

/// The mass and error vector at a neighbour `j` after it receives the share `alpha` of the
/// mass `m_i` of the removed vertex (paper Eq. 5 and Eq. 8):
/// `m̃_j = m_j + α m_i`, `t̃_j = (α m_i (R_ij t_i + e_ji) + m_j t_j) / m̃_j`, with `R_ij` the
/// rotation `rotation = (cos, sin)` from `i`'s frame to `j`'s and `e_ji` the edge vector from
/// `j` to `i` in `j`'s frame. With `rotation = (1, 0)` this is Eq. 7 in the plane: if `t_i`,
/// `t_j` point from `x_i`, `x_j` to the centres of mass of `m_i`, `m_j`, then `t̃_j` points from
/// `x_j` to the centre of mass of both. If `m̃_j = 0` the vector stays `t_j`.
///
/// [formula: iem-center-of-mass-update]
#[inline]
pub fn transport_update(
    alpha: f64,
    m_i: f64,
    t_i: [f64; 2],
    rotation: [f64; 2],
    edge: [f64; 2],
    m_j: f64,
    t_j: [f64; 2],
) -> (f64, [f64; 2]) {
    let a = alpha * m_i;
    let m = a + m_j;
    if m == 0.0 {
        return (m, t_j);
    }
    let [c, s] = rotation;
    let x = c * t_i[0] - s * t_i[1] + edge[0];
    let y = s * t_i[0] + c * t_i[1] + edge[1];
    (m, [(a * x + m_j * t_j[0]) / m, (a * y + m_j * t_j[1]) / m])
}

/// The memoryless cost of one channel (paper Eq. 6, second form; Sec. 5.1):
/// `Σ_j α_ij m_i ℓ_ij`. It equals `Σ_j m̃_j ‖t̃_j‖` when every `t` is zero, and it is the
/// 1-Wasserstein distance of moving the mass `m_i` at `i` to its neighbours over the flat star
/// (a single source: the coupling is forced; the distances are the `ℓ_ij`).
///
/// [formula: iem-memoryless-cost]
pub fn memoryless_cost(m_i: f64, alphas: &[f64], lengths: &[f64]) -> f64 {
    alphas.iter().zip(lengths).map(|(a, l)| a * m_i * l).sum()
}

/// The masses and error vectors of every channel at every fine vertex (SPEC §2.2, §4).
/// Masses of removed vertices are zero (they live at the neighbours). Vectors are in the
/// vertex's normalized frame (S4; with [`BoundaryFrames::Reference`], the reference's frame at
/// boundary vertices).
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelState {
    k: usize,
    weights: Vec<f64>,
    labels: Vec<String>,
    mass: Vec<f64>,
    vector: Vec<[f64; 2]>,
}

/// One neighbour's new mass and vector in one channel (a commit's write).
pub(crate) type Update = (u32, usize, f64, [f64; 2]);

impl ChannelState {
    /// The initial state on `mesh` (paper Alg. 1 line 5): the masses of the channels at its
    /// vertices from its current curvature and areas, every vector zero. Channels of weight 0
    /// are left out. With `fold_weights` (the reference, R2) each weight is multiplied into the
    /// masses and the weights become 1 (the cost is the same up to rounding: `t̃` does not
    /// depend on a common factor of the masses).
    ///
    /// # Errors
    /// [`IemError::InvalidOptions`]: a weight negative or not finite, a `PerVertex` array of
    /// the wrong length or with a value that is not finite.
    pub fn new(
        mesh: &CoarseMesh,
        weights: &ChannelWeights,
        fold_weights: bool,
    ) -> Result<ChannelState, IemError> {
        let n_input = mesh.n_fine_vertices();
        // Masses are indexed by fine id, the vertices inserted by pre-refinement included
        // (user-supplied per-vertex values are 0 there).
        let n_fine = mesh.fine_id_bound();
        let bad = |what: String| Err(IemError::InvalidOptions { what });
        // Per channel: label, weight, values over the fine ids.
        let mut columns: Vec<(String, f64, Vec<f64>)> = Vec::new();
        let curvature = |sign: f64| {
            let mut x = vec![0.0; n_fine];
            for v in 0..mesh.n_vertices() {
                x[mesh.fine_vertex(v) as usize] = (sign * mesh.curvature(v)).max(0.0);
            }
            x
        };
        for (index, ch) in weights.channels.iter().enumerate() {
            if !(ch.weight.is_finite() && ch.weight >= 0.0) {
                return bad(format!(
                    "channel {index} ({:?}) has weight {}, which is not finite and non-negative",
                    label_of(&ch.source),
                    ch.weight
                ));
            }
            match &ch.source {
                ChannelSource::CurvaturePositive => {
                    columns.push(("K+".into(), ch.weight, curvature(1.0)));
                }
                ChannelSource::CurvatureNegative => {
                    columns.push(("K-".into(), ch.weight, curvature(-1.0)));
                }
                ChannelSource::Area => {
                    let lumped = mesh.lumped_mass();
                    let total: f64 = lumped.iter().sum();
                    let mut x = vec![0.0; n_fine];
                    for (v, m) in lumped.into_iter().enumerate() {
                        x[mesh.fine_vertex(v) as usize] = m / total;
                    }
                    columns.push(("area".into(), ch.weight, x));
                }
                ChannelSource::PerVertex(values) => {
                    if values.len() != n_input {
                        return bad(format!(
                            "channel {index} has {} values, but the mesh has {n_input} fine \
                             vertices",
                            values.len()
                        ));
                    }
                    if let Some(v) = values.iter().position(|x| !x.is_finite()) {
                        return bad(format!(
                            "channel {index}: the value at fine vertex {v} is {}",
                            values[v]
                        ));
                    }
                    let mut pos = vec![0.0; n_fine];
                    let mut neg = vec![0.0; n_fine];
                    for v in 0..mesh.n_vertices() {
                        let f = mesh.fine_vertex(v) as usize;
                        let x = values.get(f).copied().unwrap_or(0.0);
                        pos[f] = x.max(0.0);
                        neg[f] = (-x).max(0.0);
                    }
                    let signed = neg.iter().any(|&x| x > 0.0);
                    columns.push((format!("user{index}+"), ch.weight, pos));
                    if signed {
                        columns.push((format!("user{index}-"), ch.weight, neg));
                    }
                }
            }
        }
        columns.retain(|c| c.1 > 0.0);
        let k = columns.len();
        let mut mass = vec![0.0; n_fine * k];
        for (p, (_, w, x)) in columns.iter().enumerate() {
            for (v, &m) in x.iter().enumerate() {
                mass[v * k + p] = if fold_weights { w * m } else { m };
            }
        }
        Ok(ChannelState {
            k,
            weights: columns
                .iter()
                .map(|c| if fold_weights { 1.0 } else { c.1 })
                .collect(),
            labels: columns.into_iter().map(|c| c.0).collect(),
            mass,
            vector: vec![[0.0; 2]; n_fine * k],
        })
    }

    /// The number of tracked channels (after splitting signed arrays and leaving out channels
    /// of weight 0).
    pub fn n_channels(&self) -> usize {
        self.k
    }

    /// The weight of channel `p`.
    pub fn weight(&self, p: usize) -> f64 {
        self.weights[p]
    }

    /// The label of channel `p` (`K+`, `K-`, `area`, `user<i>+`, `user<i>-`).
    pub fn label(&self, p: usize) -> &str {
        &self.labels[p]
    }

    /// The mass of channel `p` at fine vertex `v`.
    pub fn mass(&self, v: u32, p: usize) -> f64 {
        self.mass[v as usize * self.k + p]
    }

    /// The error vector of channel `p` at fine vertex `v`, in its frame.
    pub fn vector(&self, v: u32, p: usize) -> [f64; 2] {
        self.vector[v as usize * self.k + p]
    }

    /// The total mass of channel `p` (conserved by removals).
    pub fn total_mass(&self, p: usize) -> f64 {
        self.mass.iter().skip(p).step_by(self.k.max(1)).sum()
    }

    /// The cost `C_i` of the removal `prepared` (paper Eq. 6, Sec. 5.4.1; SPEC §4.4), in the
    /// frames `frames`. Does not change the state.
    pub fn removal_cost(
        &self,
        prepared: &PreparedRemoval,
        options: &MetricOptions,
        frames: Frames,
    ) -> f64 {
        let mut alphas = Vec::new();
        self.evaluate(prepared, options, frames, &mut alphas, None)
    }

    /// Commits the transport of the removal `prepared` (Eq. 5 and 8, or Eq. 5 alone for the
    /// memoryless metric): every neighbour's masses and vectors are updated, the removed
    /// vertex's are dropped. Returns the cost, as [`removal_cost`](Self::removal_cost).
    pub fn apply_removal(
        &mut self,
        prepared: &PreparedRemoval,
        options: &MetricOptions,
        frames: Frames,
    ) -> f64 {
        let mut alphas = Vec::new();
        let mut updates = Vec::new();
        self.apply_with(prepared, options, frames, &mut alphas, &mut updates)
    }

    /// [`apply_removal`](Self::apply_removal) with caller-owned scratch buffers.
    pub(crate) fn apply_with(
        &mut self,
        prepared: &PreparedRemoval,
        options: &MetricOptions,
        frames: Frames,
        alphas: &mut Vec<f64>,
        updates: &mut Vec<Update>,
    ) -> f64 {
        updates.clear();
        let cost = self.evaluate(prepared, options, frames, alphas, Some(updates));
        self.commit(updates, prepared.vertex);
        cost
    }

    /// Writes the updates of [`evaluate`](Self::evaluate) and drops the removed vertex's
    /// masses and vectors.
    pub(crate) fn commit(&mut self, updates: &[Update], removed: u32) {
        for &(j, p, m, t) in updates {
            let at = j as usize * self.k + p;
            self.mass[at] = m;
            self.vector[at] = t;
        }
        let i = removed as usize * self.k;
        for p in 0..self.k {
            self.mass[i + p] = 0.0;
            self.vector[i + p] = [0.0; 2];
        }
    }

    /// A state without channels (a placeholder).
    pub(crate) fn empty() -> ChannelState {
        ChannelState {
            k: 0,
            weights: Vec::new(),
            labels: Vec::new(),
            mass: Vec::new(),
            vector: Vec::new(),
        }
    }

    /// The cost of `prepared` with a caller-owned buffer for the weights; with `updates`, the
    /// new masses and vectors of the neighbours are pushed there.
    pub(crate) fn evaluate(
        &self,
        prepared: &PreparedRemoval,
        options: &MetricOptions,
        frames: Frames,
        alphas: &mut Vec<f64>,
        mut updates: Option<&mut Vec<Update>>,
    ) -> f64 {
        let nb = &prepared.neighbours;
        let deltas = nb
            .iter()
            .map(|n| (n.curvature_after - n.curvature_before).abs());
        weights_into(prepared.curvature, deltas, options.flat, alphas);
        let k = self.k;
        let i = prepared.vertex as usize * k;
        let i_boundary = prepared.kind != RemovalKind::Interior;
        // The scale of the reference's boundary frames, A = φ · 2Θ/(Θ + π) (R3).
        let reference_frames = options.boundary_frames == BoundaryFrames::Reference;
        let scale = |boundary: bool, curvature: f64| {
            if reference_frames && boundary {
                let theta = PI - curvature;
                2.0 * theta / (theta + PI)
            } else {
                1.0
            }
        };
        let mut cost = 0.0;
        if !options.memory {
            // Eq. 6, second form: Σ_p w^p Σ_j α_ij m_i^p ℓ_ij.
            let mut moved = 0.0;
            for (a, n) in alphas.iter().zip(nb) {
                moved += a * n.length;
            }
            for p in 0..k {
                let m_i = self.mass[i + p];
                cost += self.weights[p] * m_i * moved;
            }
            if let Some(out) = updates.as_deref_mut() {
                for (a, n) in alphas.iter().zip(nb) {
                    let j = n.vertex as usize * k;
                    for p in 0..k {
                        let am = a * self.mass[i + p];
                        if am != 0.0 {
                            out.push((n.vertex, p, self.mass[j + p] + am, [0.0; 2]));
                        }
                    }
                }
            }
            return cost;
        }
        for (a, n) in alphas.iter().zip(nb) {
            let (phi_i, phi_j) = match frames {
                Frames::AfterFlattening => (
                    n.phi_ij * scale(i_boundary, 0.0),
                    n.phi_ji * scale(n.boundary, n.curvature_after),
                ),
                Frames::BeforeFlattening => (
                    n.phi_ij_before * scale(i_boundary, prepared.curvature),
                    n.phi_ji_before * scale(n.boundary, n.curvature_before),
                ),
                Frames::ReferenceCommit => (
                    n.phi_ij_before * scale(i_boundary, prepared.curvature),
                    n.phi_ji * scale(n.boundary, n.curvature_after),
                ),
            };
            let angle = phi_j + PI - phi_i;
            let rotation = [angle.cos(), angle.sin()];
            let edge = [n.length * phi_j.cos(), n.length * phi_j.sin()];
            let j = n.vertex as usize * k;
            for p in 0..k {
                let m_i = self.mass[i + p];
                let m_j = self.mass[j + p];
                let t_j = self.vector[j + p];
                if a * m_i == 0.0 {
                    if !options.receiving_only {
                        cost += self.weights[p] * m_j * t_j[0].hypot(t_j[1]);
                    }
                    continue;
                }
                let (m, t) =
                    transport_update(*a, m_i, self.vector[i + p], rotation, edge, m_j, t_j);
                cost += self.weights[p] * m * t[0].hypot(t[1]);
                if let Some(out) = updates.as_deref_mut() {
                    out.push((n.vertex, p, m, t));
                }
            }
        }
        cost
    }
}

fn label_of(source: &ChannelSource) -> &'static str {
    match source {
        ChannelSource::CurvaturePositive => "CurvaturePositive",
        ChannelSource::CurvatureNegative => "CurvatureNegative",
        ChannelSource::Area => "Area",
        ChannelSource::PerVertex(_) => "PerVertex",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weights_sum_to_one_and_flat_rule() {
        let a = redistribution_weights(0.3, &[0.1, -0.2, 0.0], &[0.2, -0.1, 0.1]);
        assert!((a.iter().sum::<f64>() - 1.0).abs() < 1e-15);
        assert!(a.iter().all(|&x| x >= 0.0));
        let flat = redistribution_weights(1e-13, &[0.1, 0.2], &[0.3, 0.2]);
        assert_eq!(flat, vec![0.5, 0.5]);
        assert!(redistribution_weights(1.0, &[], &[]).is_empty());
    }

    #[test]
    fn plane_centre_of_mass() {
        // Masses 2 at x_i = (1, 0) (centre of mass at (1, 1): t_i = (0, 1)) and 1 at x_j = 0
        // (t_j = (0, -1)), all of m_i to j: the combined centre of mass is (2/3, 1/3).
        let (m, t) = transport_update(
            1.0,
            2.0,
            [0.0, 1.0],
            [1.0, 0.0],
            [1.0, 0.0],
            1.0,
            [0.0, -1.0],
        );
        assert_eq!(m, 3.0);
        assert!((t[0] - 2.0 / 3.0).abs() < 1e-15 && (t[1] - 1.0 / 3.0).abs() < 1e-15);
        // Memoryless equals Σ m̃ ‖t̃‖ with zero vectors.
        let alphas = [0.25, 0.75];
        let lengths = [2.0, 0.5];
        let direct: f64 = alphas
            .iter()
            .zip(&lengths)
            .map(|(&a, &l)| {
                let (m, t) =
                    transport_update(a, 4.0, [0.0; 2], [0.0, 1.0], [0.0, l], 0.0, [0.0; 2]);
                m * t[0].hypot(t[1])
            })
            .sum();
        assert!((memoryless_cost(4.0, &alphas, &lengths) - direct).abs() < 1e-15);
    }
}
