//! The integer formulas of the paper on plain integers: Eq. 1, 2 (corner coordinates and
//! emanating counts), the validity of a face, the reversal of a crossing, one step of Supp.
//! Alg. 1, the index of Supp. Alg. 3, Eq. 3 and 4 (the flip) and Eq. 6.
//!
//! Every function here is the Rust counterpart of a definition of the Lean oracle
//! (`lean/oracle/DdgOracle/IntegerCoords/Formulas.lean`), which `lean/proofs` proves theorems
//! about, and is compared with it by exact integer equality (`tests/oracle.rs`).
//!
//! # Face convention (SPEC §2.1, the oracle's)
//!
//! A face has corners `0, 1, 2` and sides `0, 1, 2`; side `k` is the halfedge `3f + k`, from
//! corner `k` to corner `k + 1`, and its crossings are numbered `0, 1, …` from corner `k`. Corner
//! `k` lies between side `k + 2` (which ends at it) and side `k` (which starts at it); the side
//! opposite corner `k` is side `k + 1`.
//!
//! # Arithmetic
//!
//! Normal coordinates are `i64` (SPEC I2). Eq. 1 and 2 are evaluated in `i128`, where no sum of
//! three `i64` can overflow, and their results always fit an `i64` again (a corner coordinate is at
//! most the smaller of the two sides at its corner, an emanating count at most the opposite side).
//! Eq. 1 is **never** a floor division (SPEC I14): its numerator is never negative, and an odd
//! numerator is reported (`None`), not rounded.

/// `n⁺ = max(n, 0)`: the number of transversal crossings of an edge with `n` as its normal
/// coordinate.
///
/// [formula: ic-normal-sign]
#[inline]
pub fn n_plus(n: i64) -> i64 {
    n.max(0)
}

/// `n⁻ = −min(n, 0)`: `1` for an edge that runs along an input edge (`n = −1`), `0` otherwise.
/// Every valid state has `n ≥ −1`; for `n = i64::MIN`, which no state has, the result
/// saturates at `i64::MAX`.
///
/// [formula: ic-normal-sign]
#[inline]
pub fn n_minus(n: i64) -> i64 {
    n.min(0).checked_neg().unwrap_or(i64::MAX)
}

/// Eq. 2, `e^{ij}_k = max(0, n⁺_ij − n⁺_jk − n⁺_ki)`: the number of arcs from corner `k` to the
/// opposite side `ij`. Arguments: the opposite side, then the two sides at the corner (in either
/// order).
///
/// [formula: ic-emanating]
#[inline]
pub fn emanating(opposite: i64, a: i64, b: i64) -> i64 {
    let e = n_plus(opposite) as i128 - n_plus(a) as i128 - n_plus(b) as i128;
    // 0 ≤ e ≤ n⁺(opposite), so it fits.
    e.max(0) as i64
}

/// Twice Eq. 1 at corner `k` of the face `ijk`:
/// `max(0, n⁺_jk + n⁺_ki − n⁺_ij) − e^{jk}_i − e^{ki}_j`. Arguments: the opposite side `ij`, then
/// the two sides `jk`, `ki` at the corner. Never negative (Lean `cornerNumerator_nonneg`); it can
/// exceed `i64::MAX` when both sides at the corner are close to it, hence `i128`.
///
/// [formula: ic-corner-coordinates]
#[inline]
pub fn corner_numerator(opposite: i64, a: i64, b: i64) -> i128 {
    let (o, x, y) = (
        n_plus(opposite) as i128,
        n_plus(a) as i128,
        n_plus(b) as i128,
    );
    (x + y - o).max(0) - emanating(a, opposite, b) as i128 - emanating(b, opposite, a) as i128
}

/// Eq. 1, `c^{ij}_k`: the number of arcs that cut corner `k`, or `None` when the numerator is odd
/// (the face is not valid, SPEC §2.4 item 1) or negative (impossible).
///
/// [formula: ic-corner-coordinates]
#[inline]
pub fn corner_coordinate(opposite: i64, a: i64, b: i64) -> Option<i64> {
    let m = corner_numerator(opposite, a, b);
    if m < 0 || m % 2 != 0 {
        None
    } else {
        // m / 2 ≤ min(n⁺a, n⁺b): fits.
        Some((m / 2) as i64)
    }
}

/// Whether the face with side counts `n` (side `k` from corner `k` to `k + 1`) is valid: the
/// three numerators of Eq. 1 are even (SPEC §2.4 item 1).
///
/// [formula: ic-face-validity]
pub fn face_valid(n: [i64; 3]) -> bool {
    (0..3).all(|k| corner_numerator(n[(k + 1) % 3], n[k], n[(k + 2) % 3]) % 2 == 0)
}

/// The numerator of Eq. 1 at corner `k` of the face with side counts `n` (opposite side `k + 1`,
/// sides `k` and `k + 2` at the corner).
///
/// [formula: ic-corner-coordinates]
#[inline]
pub fn face_corner_numerator(n: [i64; 3], k: usize) -> i128 {
    corner_numerator(n[(k + 1) % 3], n[k % 3], n[(k + 2) % 3])
}

/// Eq. 2 at corner `k` of the face with side counts `n`: the arcs from corner `k` to side `k + 1`.
///
/// [formula: ic-emanating]
#[inline]
pub fn face_emanating(n: [i64; 3], k: usize) -> i64 {
    emanating(n[(k + 1) % 3], n[k % 3], n[(k + 2) % 3])
}

/// The decoded corner data of a valid face: what Supp. Alg. 1, 3, 4 and 9 read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceCoordinates {
    /// `n⁺` of side `k`.
    pub n: [i64; 3],
    /// `c` of corner `k` (Eq. 1): arcs that cut it, from side `k + 2` to side `k`.
    pub c: [i64; 3],
    /// `e` of corner `k` (Eq. 2): arcs from corner `k` to side `k + 1`. At most one is positive
    /// (SPEC §2.4 item 2, Lean `emanating_exclusive`).
    pub e: [i64; 3],
}

impl FaceCoordinates {
    /// Decodes a face from the normal coordinates of its sides; `Err(k)` names the first corner
    /// whose numerator of Eq. 1 is odd.
    ///
    /// [formula: ic-corner-coordinates]
    pub fn of(n: [i64; 3]) -> Result<FaceCoordinates, usize> {
        let mut c = [0i64; 3];
        for (k, ck) in c.iter_mut().enumerate() {
            *ck = corner_coordinate(n[(k + 1) % 3], n[k], n[(k + 2) % 3]).ok_or(k)?;
        }
        Ok(FaceCoordinates {
            n: n.map(n_plus),
            c,
            e: [0, 1, 2].map(|k| face_emanating(n, k)),
        })
    }

    /// The first crossing (from the side's tail) of side `k + 1` that an arc emanating from
    /// corner `k` reaches: the arcs cutting corner `k + 1` come first. The emanating arcs are the
    /// crossings `first … first + e_k − 1`.
    #[inline]
    pub fn first_emanating(&self, k: usize) -> i64 {
        self.c[(k + 1) % 3]
    }
}

/// The reversal of a combinatorial crossing `(ij→, p) ↦ (ji→, n − p − 1)`: the position part.
///
/// [formula: ic-reverse-crossing]
#[inline]
pub fn reverse_crossing(n: i64, p: i64) -> i64 {
    n - p - 1
}

/// Where a curve that crosses an edge goes next in the face it enters (Supp. Alg. 1, Fig. 7).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Step {
    /// It leaves through `ik→` at this crossing index (Fig. 7 case 3, "goes right").
    ToIk(i64),
    /// It leaves through `kj→` at this crossing index (Fig. 7 case 1, "goes left").
    ToKj(i64),
    /// It ends at the opposite vertex `k` (Fig. 7 case 2).
    Vertex,
}

/// Why [`trace_step`] has no answer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepError {
    /// The face `jik` is not valid (an odd numerator of Eq. 1).
    InvalidFace,
    /// `p` is not in `[0, n⁺_ij)`.
    OutOfRange,
}

/// One step of Supp. Alg. 1: the curve crosses `ij→` at index `p` (from `i`) and enters the face
/// `jik` of the twin, whose opposite vertex is `k`. With `c^{kj}_i` the corner `i` and `c^{ik}_j`
/// the corner `j` of `jik`: `p < c^{kj}_i` goes to `ik→` at `p`; `p ≥ n_ij − c^{ik}_j` goes to
/// `kj→` at `n_kj + p − n_ij`; otherwise the curve ends at `k`.
///
/// [formula: ic-trace-step]
pub fn trace_step(n_ij: i64, n_ik: i64, n_kj: i64, p: i64) -> Result<Step, StepError> {
    if !face_valid([n_ij, n_ik, n_kj]) {
        return Err(StepError::InvalidFace);
    }
    let (nij, nkj) = (n_plus(n_ij), n_plus(n_kj));
    if p < 0 || p >= nij {
        return Err(StepError::OutOfRange);
    }
    let ci = corner_coordinate(n_kj, n_ij, n_ik).ok_or(StepError::InvalidFace)?;
    let cj = corner_coordinate(n_ik, n_ij, n_kj).ok_or(StepError::InvalidFace)?;
    Ok(if p < ci {
        Step::ToIk(p)
    } else if nij - cj <= p {
        // p − (n_ij − c_j) < c_j ≤ n_kj, so this is in [n_kj − c_j, n_kj).
        Step::ToKj(nkj - (nij - p))
    } else {
        Step::Vertex
    })
}

/// Supp. Alg. 3 line 9: the curve leaves `a ∈ V⁰` and first crosses `ij→` at `p`, where `ija` is
/// the face of `ij→` and `ai→ = Next(Next(ij→))` (the paper calls it `aj`, SPEC §3.1). Its local
/// index around `a` is `r_ai + n⁻_ai + (p − c^{ja}_i)`, reduced modulo `deg₀(a)` at an interior
/// vertex (`deg = 0`: a boundary vertex, whose numbering does not wrap). `None` when the face
/// `ija` is not valid or a value leaves the `i64` range.
///
/// Alg. 4 inverts it: the crossing of the `q`-th input halfedge after `r_ai` is
/// `p = c^{ja}_i + (k − r_ai − n⁻_ai)`, which the paper writes `k + n_ai − r_ai` (both agree,
/// since `c^{ja}_i = n⁺_ai` whenever `a` emanates, SPEC §3.1).
///
/// [formula: ic-identify-edge]
pub fn identify_index(r_ai: i64, n_ai: i64, n_ij: i64, n_ja: i64, p: i64, deg: u32) -> Option<i64> {
    let c = corner_coordinate(n_ja, n_ij, n_ai)?;
    let k = r_ai
        .checked_add(n_minus(n_ai))?
        .checked_add(p.checked_sub(c)?)?;
    Some(if deg == 0 {
        k
    } else {
        k.rem_euclid(deg as i64)
    })
}

/// Eq. 6: `n_mk = max(n_ki, n_jk, 0)` for the new edge from the split point `m` of an edge `ij`
/// with `n_ij = −1` to the opposite corner `k` (the halves keep `n = −1`); the `0` matters when
/// both other sides run along input edges.
///
/// [formula: ic-split-edge-normal]
#[inline]
pub fn split_edge_normal(n_ki: i64, n_jk: i64) -> i64 {
    n_ki.max(n_jk).max(0)
}

/// Twice Eq. 3 without `n⁻_ij` (the "doubled form", so that nothing is halved before it is
/// known to be even, SPEC I14): edge `ij` with the faces `ijk` and `jil`, arguments
/// `n_ij, n_jk, n_ki, n_il, n_lj`,
///
/// ```text
/// 2c^{ji}_l + 2c^{ij}_k + |c^{il}_j − c^{ki}_j| + |c^{lj}_i − c^{jk}_i| − e^{ji}_l − e^{ij}_k
///   + 2(e^{lj}_i + e^{jk}_i + e^{il}_j + e^{ki}_j)
/// ```
///
/// with every corner quantity that of one of the two old faces: `c^{ji}_l` is corner `l` of
/// `jil` and `c^{ij}_k` corner `k` of `ijk` (the paper prints `c^{jk}_l + c^{il}_k`, which name
/// no corner; SPEC §4.2), `c^{il}_j`, `e^{il}_j` are corner `j` of `jil`, `c^{ki}_j`, `e^{ki}_j`
/// corner `j` of `ijk`, and so on. `None` when one of the two faces is not valid (an odd
/// numerator of Eq. 1). Evaluated in `i128`, where it cannot overflow for `i64` arguments.
///
/// [formula: ic-flip-normal]
pub fn flip_doubled(n_ij: i64, n_jk: i64, n_ki: i64, n_il: i64, n_lj: i64) -> Option<i128> {
    let c = |opp, a, b| corner_coordinate(opp, a, b).map(i128::from);
    let e = |opp, a, b| i128::from(emanating(opp, a, b));
    let cl = c(n_ij, n_il, n_lj)?; // c^{ji}_l, corner l of jil
    let ck = c(n_ij, n_jk, n_ki)?; // c^{ij}_k, corner k of ijk
    let cj_l = c(n_il, n_ij, n_lj)?; // c^{il}_j, corner j of jil
    let cj_r = c(n_ki, n_ij, n_jk)?; // c^{ki}_j, corner j of ijk
    let ci_l = c(n_lj, n_ij, n_il)?; // c^{lj}_i, corner i of jil
    let ci_r = c(n_jk, n_ij, n_ki)?; // c^{jk}_i, corner i of ijk
    let el = e(n_ij, n_il, n_lj); // e^{ji}_l
    let ek = e(n_ij, n_jk, n_ki); // e^{ij}_k
    let ei_l = e(n_lj, n_ij, n_il); // e^{lj}_i
    let ei_r = e(n_jk, n_ij, n_ki); // e^{jk}_i
    let ej_l = e(n_il, n_ij, n_lj); // e^{il}_j
    let ej_r = e(n_ki, n_ij, n_jk); // e^{ki}_j
    Some(
        2 * cl + 2 * ck + (cj_l - cj_r).abs() + (ci_l - ci_r).abs() - el - ek
            + 2 * (ei_l + ei_r + ej_l + ej_r),
    )
}

/// Why [`flip_normal`] has no answer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlipNormalError {
    /// One of the two faces `ijk`, `jil` is not valid (an odd numerator of Eq. 1).
    InvalidFace,
    /// The doubled form of Eq. 3 is odd (impossible for two valid faces, Lean `flip_count`;
    /// reported instead of rounded, SPEC I14).
    OddDoubled(i128),
    /// The result does not fit an `i64` (SPEC I2).
    Overflow,
}

/// Eq. 3: the normal coordinate of the new edge `kl` after flipping `ij` (faces `ijk`, `jil`;
/// arguments `n_ij, n_jk, n_ki, n_il, n_lj`): [`flip_doubled`] halved, plus `n⁻_ij`. The
/// result can be below `−1` for face pairs that are combinatorially valid but have two or more
/// arcs from `k` to `l` (SPEC §4.2); the caller decides what that means. Applying it to the new
/// edge, `flip_normal(n_kl, n_ki, n_il, n_lj, n_jk)`, returns `n_ij` (`ic-flip-involution`).
///
/// [formula: ic-flip-normal]
pub fn flip_normal(
    n_ij: i64,
    n_jk: i64,
    n_ki: i64,
    n_il: i64,
    n_lj: i64,
) -> Result<i64, FlipNormalError> {
    let d = flip_doubled(n_ij, n_jk, n_ki, n_il, n_lj).ok_or(FlipNormalError::InvalidFace)?;
    if d % 2 != 0 {
        return Err(FlipNormalError::OddDoubled(d));
    }
    i64::try_from(d / 2 + i128::from(n_minus(n_ij))).map_err(|_| FlipNormalError::Overflow)
}

/// Eq. 4, `r_kl = (r_ki + e^{il}_k + n⁻_ki) mod deg₀(k)`: the roundabout of the new halfedge
/// `k → l` after a flip, from the roundabout of `k → i` (the halfedge just clockwise of it),
/// the emanating count `e^{il}_k` of corner `k` in the **new** face `kil`, and `n⁻_ki`.
/// `modulus` is [`LocalIndices::modulus`](crate::LocalIndices::modulus): `deg₀(k)` at an
/// interior vertex, `0` (no reduction, SPEC I3) at a boundary vertex. `None` when the result
/// leaves the `i64` range (only possible without reduction). The same with `(l, j, k)` for
/// `r_lk`.
///
/// [formula: ic-roundabout-update]
pub fn roundabout_update(r_ki: i64, e_new: i64, n_minus_ki: i64, modulus: u32) -> Option<i64> {
    let s = i128::from(r_ki) + i128::from(e_new) + i128::from(n_minus_ki);
    let s = if modulus == 0 {
        s
    } else {
        s.rem_euclid(i128::from(modulus))
    };
    i64::try_from(s).ok()
}

/// What splitting one side of a face at a point between its crossings does to the face (SPEC
/// I20): see [`split_side`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SideSplit {
    /// The number of arcs of the face that cross the new edge from the split point `m` to the
    /// opposite corner: its normal coordinate (never `−1`: no input edge ends at `m`).
    pub to_opposite: i64,
    /// The number of arcs from the opposite corner that end on the split side before `m`
    /// (counted from the side's tail): the input halfedges at the opposite corner between its
    /// halfedge to the side's tail and the new edge, for Eq. 4's analogue
    /// `r_{om} = r_{o,tail} + n⁻ + emanating_before`.
    pub emanating_before: i64,
}

/// The integers of a split of side `s` (from corner `s` to corner `s + 1`) of a face with side
/// counts `n`, at a point with `before` crossings of side `s` between corner `s` and it
/// (`0 ≤ before ≤ n⁺_s`, the point not on a crossing). With `a, b, c` the corners `s, s + 1,
/// s + 2` and `m` the point, the straight new edge `mc` crosses exactly the arcs whose ends lie
/// on different sides of the chord `mc` on the face's boundary (`ic-chord-interleave`):
///
/// ```text
/// n_mc = c_c + e_a + e_b + max(0, c_a − before) + max(0, before − (n⁺_ab − c_b))
/// ```
///
/// (the arcs around `c`, from `a` to `bc`, from `b` to `ca`, the arcs around `a` that reach `ab`
/// beyond `m`, and those around `b` that reach it before `m`; the arcs from `c` end at an end
/// of the chord). For a side along an input edge (`n_ab = −1`, `before = 0`) this is Eq. 6,
/// `max(n_ca, n_bc, 0)`. `None` for a face that is not valid or `before` out of range.
///
/// [formula: ic-split-crossed-edge]
pub fn split_side(n: [i64; 3], s: usize, before: i64) -> Option<SideSplit> {
    let f = FaceCoordinates::of(n).ok()?;
    let (a, b, c) = (s % 3, (s + 1) % 3, (s + 2) % 3);
    let n_ab = f.n[a];
    if !(0..=n_ab).contains(&before) {
        return None;
    }
    let i = |x: i64| i128::from(x);
    let to_opposite = i(f.c[c])
        + i(f.e[a])
        + i(f.e[b])
        + (i(f.c[a]) - i(before)).max(0)
        + (i(before) - (i(n_ab) - i(f.c[b]))).max(0);
    let emanating_before = (before - f.c[a]).clamp(0, f.e[c]);
    Some(SideSplit {
        to_opposite: i64::try_from(to_opposite).ok()?,
        emanating_before,
    })
}

/// What splitting one side of a face **at** its crossing number `j` does to the face (SPEC
/// I28): see [`split_side_at_crossing`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SideSplitAtCrossing {
    /// The normal coordinate of the new edge `mc` from the split point `m` to the opposite
    /// corner `c`: `−1` when the arc through `m` is one of the arcs from `c` (the edge then runs
    /// along that input edge), else the number of other arcs it crosses.
    pub to_opposite: i64,
    /// The number of arcs from `c` that end on the split side before `m`, as
    /// [`SideSplit::emanating_before`] (the arc through `m` itself not counted).
    pub emanating_before: i64,
}

/// The integers of a split of side `s` (from corner `s`) of a face with side counts `n` **at**
/// its crossing number `j` (`0 ≤ j < n⁺_s`, counted from corner `s`): the new vertex `m` lies on
/// the arc `A` of the input curve that crosses the side there, so `A` now starts at `m`. With
/// `a, b, c` the corners `s, s + 1, s + 2` and `c_*`, `e_*` the corner and emanating counts:
///
/// ```text
/// j <  c_a             (A cuts a):   n_mc = split_side(n, s, j + 1).to_opposite
/// c_a ≤ j < c_a + e_c  (A is from c): n_mc = −1,  e_before = j − c_a
/// j ≥  c_a + e_c       (A cuts b):   n_mc = split_side(n, s, j).to_opposite
/// ```
///
/// The curve `A` shares its end `m` with the new edge, so it does not cross it: of the two sides
/// of `A` on which `m` may be put (`before = j` or `j + 1`, [`split_side`]) it takes the one
/// where `A` separates `m` from the rest of the face less than from `c`, and for an arc ending
/// at `c` the new edge is `A` itself. `None` for a face that is not valid or `j` out of range.
///
/// Formula id `ic-split-at-crossing` (SPEC I28; proved in Lean,
/// `lean/proofs/DdgProofs/IntegerCoords/SplitAtCrossing.lean`, `Complex/SplitAtCrossing.lean`).
///
/// [formula: ic-split-at-crossing]
pub fn split_side_at_crossing(n: [i64; 3], s: usize, j: i64) -> Option<SideSplitAtCrossing> {
    let f = FaceCoordinates::of(n).ok()?;
    let (a, c) = (s % 3, (s + 2) % 3);
    if !(0..f.n[a]).contains(&j) {
        return None;
    }
    let fan = j >= f.c[a] && j < f.c[a] + f.e[c];
    let before = if j < f.c[a] { j + 1 } else { j };
    let x = split_side(n, s, before)?;
    Some(SideSplitAtCrossing {
        to_opposite: if fan { -1 } else { x.to_opposite },
        emanating_before: x.emanating_before,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eq1_never_rounds_and_handles_extremes() {
        // One arc from corner 2 to side 0 (opposite it), nothing else: valid.
        assert_eq!(FaceCoordinates::of([1, 0, 0]).unwrap().e, [0, 0, 1]);
        // Odd sum without an emanating corner: invalid.
        assert_eq!(FaceCoordinates::of([1, 1, 1]), Err(0));
        let m = i64::MAX;
        assert_eq!(corner_numerator(0, m, m), 2 * m as i128);
        assert_eq!(corner_coordinate(0, m, m), Some(m));
        assert_eq!(emanating(m, -1, -1), m);
        assert_eq!(n_minus(i64::MIN), i64::MAX);
    }

    #[test]
    fn identify_inverts_the_extraction_offset() {
        // Face ija with a emanating 3 arcs to ij: n_ij = n_ai + n_ja + 3.
        let (n_ai, n_ja) = (2, 1);
        let n_ij = n_ai + n_ja + 3;
        let c_i = corner_coordinate(n_ja, n_ij, n_ai).unwrap();
        assert_eq!(c_i, n_ai);
        for q in 0..3 {
            assert_eq!(
                identify_index(5, n_ai, n_ij, n_ja, c_i + q, 7),
                Some((5 + q) % 7)
            );
        }
    }
}
