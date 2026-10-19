//! Combinatorial energies `E_P` (Eq. 1–2) and `π_max` (Eq. 5), SPEC §3.2–§3.3.

use super::{Scratch, Star, StarOut};
use crate::geometry::v3::*;

/// The edge-connected bipartitions of a star of `k` faces (SPEC §3.2, D9), in enumeration order,
/// as pairs of arcs `(start, len)` of the star's face sequence (indices mod `k`):
///
/// * interior (`closed`): `(a, b)`, `0 ≤ a < b ≤ k − 1` lexicographic, `F_1 = arc(a+1, b−a)`,
///   `F_2 = arc(b+1, k−(b−a))`; `k(k−1)/2` of them;
/// * boundary: `a = 0 … k−2`, `F_1 = arc(0, a+1)`, `F_2 = arc(a+1, k−a−1)`; `k − 1` of them.
///
/// Both parts are non-empty and edge-connected.
///
/// [formula: dev-star-bipartitions]
pub fn bipartitions(k: usize, closed: bool) -> impl Iterator<Item = [(usize, usize); 2]> {
    let interior = closed.then(|| {
        (0..k).flat_map(move |a| {
            (a + 1..k).map(move |b| [((a + 1) % k, b - a), ((b + 1) % k, k - (b - a))])
        })
    });
    let boundary =
        (!closed).then(|| (0..k.saturating_sub(1)).map(move |a| [(0, a + 1), (a + 1, k - a - 1)]));
    interior.into_iter().flatten().chain(boundary.into_iter().flatten())
}

/// Mean `N̄` and width `(1/n) Σ |N_σ − N̄|²` of the arc `(start, len)` of `normals` (mod `k`).
#[inline]
fn arc_width(normals: &[[f64; 3]], (start, len): (usize, usize)) -> ([f64; 3], f64) {
    let k = normals.len();
    let mut m = [0.0; 3];
    for i in 0..len {
        m = add(m, normals[(start + i) % k]);
    }
    let m = scale(m, 1.0 / len as f64);
    let mut w = 0.0;
    for i in 0..len {
        w += norm2(sub(normals[(start + i) % k], m));
    }
    (m, w / len as f64)
}

/// Combinatorial energy `E_P,v = min_P π(P)` with
/// `π(P) = Σ_p (1/n_p) Σ_{σ∈F_p} |N_σ − N̄_p|²` (the first form of Eq. 1; equals
/// `Σ_p (1/n_p²) Σ_{unordered pairs} |N_σ1 − N_σ2|²` and `Σ_p (1 − |N̄_p|²)`, D1), and its
/// gradient for the minimizing partition (lowest enumeration index on ties):
/// `y_c = (2/n_p)(N_σ − N̄_p)` for `σ ∈ F_p`, `s_c = 0`.
///
/// Margin: second-best minus best `π` (`+∞` with fewer than two partitions). A boundary star with
/// one face has no bipartition: `e_v = 0`.
///
/// [formula: dev-combinatorial-width]
/// [formula: dev-combinatorial-gradient]
pub fn combinatorial(star: &Star, y: &mut [[f64; 3]], s: &mut [f64]) -> StarOut {
    let nrm = star.normals;
    let k = nrm.len();
    s.fill(0.0);
    y.fill([0.0; 3]);
    let mut best = f64::INFINITY;
    let mut second = f64::INFINITY;
    let mut best_p = None;
    for parts in bipartitions(k, star.closed) {
        let w = arc_width(nrm, parts[0]).1 + arc_width(nrm, parts[1]).1;
        if w < best {
            second = best;
            best = w;
            best_p = Some(parts);
        } else if w < second {
            second = w;
        }
    }
    let Some(parts) = best_p else {
        return StarOut { value: 0.0, margin: f64::INFINITY, z: [0.0; 3] };
    };
    for part in parts {
        let (m, _) = arc_width(nrm, part);
        let c = 2.0 / part.1 as f64;
        for i in 0..part.1 {
            let f = (part.0 + i) % k;
            y[f] = scale(sub(nrm[f], m), c);
        }
    }
    StarOut { value: best, margin: second - best, z: [0.0; 3] }
}

/// Active pair of an arc: the largest `|N_σ1 − N_σ2|²` within it (first found on ties).
#[derive(Clone, Copy)]
struct ArcMax {
    value: f64,
    pair: Option<(u32, u32)>,
}

/// Combinatorial max energy `e_v = min_P max_p max_{σ1,σ2∈F_p} |N_σ1 − N_σ2|²` (Eq. 5) and the
/// subgradient of the active pair `(σ1, σ2)` of the minimizing partition:
/// `y_σ1 = 2(N_σ1 − N_σ2)`, `y_σ2 = −2(N_σ1 − N_σ2)`, other `y = 0`, `s = 0`.
///
/// Arc maxima are tabulated by extending arcs, O(k³) per star. Ties: lowest partition index; within
/// an arc the first pair found; between the parts `F_1`. Margin: the smaller of (second-best minus
/// best partition value) and (active pair minus the next largest pair of the chosen partition);
/// `+∞` where there is no alternative.
///
/// [formula: dev-combinatorial-max]
pub fn combinatorial_max(
    star: &Star,
    y: &mut [[f64; 3]],
    s: &mut [f64],
    scratch: &mut Scratch,
) -> StarOut {
    let nrm = star.normals;
    let k = nrm.len();
    s.fill(0.0);
    y.fill([0.0; 3]);
    if k < 2 {
        return StarOut { value: 0.0, margin: f64::INFINITY, z: [0.0; 3] };
    }
    // Pairwise squared distances.
    let d = &mut scratch.f1;
    d.clear();
    d.resize(k * k, 0.0);
    for i in 0..k {
        for j in i + 1..k {
            let v = norm2(sub(nrm[i], nrm[j]));
            d[i * k + j] = v;
            d[j * k + i] = v;
        }
    }
    // mx[s * k + (len − 1)] = active pair of arc(s, len), len = 1 .. k−1.
    let mx = &mut scratch.arcs;
    mx.clear();
    mx.resize(k * k, ArcMaxSlot::default());
    for st in 0..k {
        let mut cur = ArcMax { value: 0.0, pair: None };
        mx[st * k] = cur.into();
        for len in 2..k {
            let new = (st + len - 1) % k;
            for m in 0..len - 1 {
                let old = (st + m) % k;
                let v = d[new * k + old];
                if v > cur.value {
                    cur = ArcMax { value: v, pair: Some((old as u32, new as u32)) };
                }
            }
            mx[st * k + len - 1] = cur.into();
        }
    }
    let mut best = f64::INFINITY;
    let mut second = f64::INFINITY;
    let mut chosen: Option<([(usize, usize); 2], ArcMax)> = None;
    for parts in bipartitions(k, star.closed) {
        let a: ArcMax = mx[parts[0].0 * k + parts[0].1 - 1].into();
        let b: ArcMax = mx[parts[1].0 * k + parts[1].1 - 1].into();
        let act = if b.value > a.value { b } else { a };
        if act.value < best {
            second = best;
            best = act.value;
            chosen = Some((parts, act));
        } else if act.value < second {
            second = act.value;
        }
    }
    let Some((parts, act)) = chosen else {
        return StarOut { value: 0.0, margin: f64::INFINITY, z: [0.0; 3] };
    };
    let mut margin = second - best;
    if let Some((p1, p2)) = act.pair {
        let (p1, p2) = (p1 as usize, p2 as usize);
        let g = scale(sub(nrm[p1], nrm[p2]), 2.0);
        y[p1] = g;
        y[p2] = scale(g, -1.0);
        // Next largest pair within the chosen partition.
        let mut next = f64::NEG_INFINITY;
        for (st, len) in parts {
            for i in 0..len {
                for j in i + 1..len {
                    let (a, b) = ((st + i) % k, (st + j) % k);
                    if (a.min(b), a.max(b)) == (p1.min(p2), p1.max(p2)) {
                        continue;
                    }
                    next = next.max(d[a * k + b]);
                }
            }
        }
        margin = margin.min(act.value - next);
    }
    StarOut { value: best, margin, z: [0.0; 3] }
}

/// Storage form of [`ArcMax`] in the scratch table.
#[derive(Clone, Copy, Default)]
pub(crate) struct ArcMaxSlot {
    value: f64,
    a: u32,
    b: u32,
}

const NO_PAIR: u32 = u32::MAX;

impl From<ArcMax> for ArcMaxSlot {
    fn from(m: ArcMax) -> Self {
        let (a, b) = m.pair.unwrap_or((NO_PAIR, NO_PAIR));
        ArcMaxSlot { value: m.value, a, b }
    }
}

impl From<ArcMaxSlot> for ArcMax {
    fn from(m: ArcMaxSlot) -> Self {
        ArcMax { value: m.value, pair: (m.a != NO_PAIR).then_some((m.a, m.b)) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bipartition_counts() {
        for k in 1..10 {
            assert_eq!(bipartitions(k, true).count(), k * (k.max(1) - 1) / 2, "k = {k}");
            assert_eq!(bipartitions(k, false).count(), k - 1, "k = {k}");
        }
        // every part non-empty, parts complementary
        for k in 2..9 {
            for closed in [true, false] {
                for [p, q] in bipartitions(k, closed) {
                    assert!(p.1 >= 1 && q.1 >= 1 && p.1 + q.1 == k);
                    assert_eq!((p.0 + p.1) % k, q.0);
                    assert_eq!((q.0 + q.1) % k, p.0);
                }
            }
        }
    }
}
