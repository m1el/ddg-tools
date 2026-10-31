import DdgOracle.IntegerCoords.Formulas
import Mathlib.Tactic.Linarith
import Mathlib.Algebra.Order.Group.Int

/-!
# Integer coordinates: arithmetic of the formulas (tier A)

`crates/ddg-integer-coords/SPEC.md` §2.2–2.4, §4.2, §4.5, §9.2; the paper Sec. 2.3, 2.5, 3.2, 3.4.
All statements are about the oracle's integer definitions (`lean/oracle/DdgOracle/IntegerCoords/
Formulas.lean`, namespace `Ddg.IntegerCoords`), i.e. about the code the oracle runs.

Face convention: side `k` runs from corner `k` to corner `k + 1`; corner `k` is opposite side
`k + 1`. A face `(n₀, n₁, n₂)` has the corner coordinates
`c_k = cornerCoord n_{k+1} n_k n_{k+2}` and the emanating counts `e_k = emanating n_{k+1} n_k n_{k+2}`.

* **`ic-normal-sign`** — `n⁺ − n⁻ = n`, `n⁺ n⁻ = 0`, both nonnegative, `n⁻ = 1` exactly on
  `n = −1` (`nPlus_sub_nMinus`, `nPlus_mul_nMinus`, `nMinus_eq_one_iff`).
* **`ic-emanating`** (arithmetic half) — `e ≥ 0`, and at most one corner of a face has `e > 0`
  (`emanating_exclusive`): item 2 of SPEC §2.4 is implied by Eq. 2 itself.
* **`ic-corner-coordinates`** (arithmetic half) — the numerator of Eq. 1 is never negative
  (`cornerNumerator_nonneg`), so the `/ 2` of the oracle never sees a negative integer (SPEC I14);
  on counts of the form `n_k = c_k + c_{k+1} + e_{k+2}` with at most one `e > 0` and `c = 0` at an
  emanating corner, Eq. 1 and 2 return `c` and `e` (`corner_eval`).
* **`ic-face-validity`** — the face is valid (the three numerators even) if and only if
  `n⁺₀ + n⁺₁ + n⁺₂` is even or some corner emanates (`faceValid_iff`), and then it has the form above
  (`exists_params_of_faceValid`).
* **`ic-reverse-crossing`** — reversal is an involution on `0 ≤ p < n` (`reverseCrossing_reverseCrossing`,
  `reverseCrossing_mem`).
* **`ic-flip-involution`** — Eq. 3 applied to the new edge of a flip returns the old normal
  coordinate, for any two valid faces on the edge (`flipNormal_flipNormal`); the doubled value is even
  and the result is at least `−1` (`flipDoubled_even`, `flipNormal_ge`).
* **`ic-split-edge-normal`** — Eq. 6 counts every arc of a face with a side carrying no
  crossing (`splitEdgeNormal_eq`; the face-split form is `DdgProofs/IntegerCoords/FaceSplit.lean`).
-/

namespace DdgProofs

namespace IntegerCoords

open Ddg.IntegerCoords

/-! ## Signs -/

/-- `n⁺ − n⁻ = n`. [formula: ic-normal-sign] -/
theorem nPlus_sub_nMinus (n : ℤ) : nPlus n - nMinus n = n := by
  simp only [nPlus, nMinus]; omega

/-- `n⁺ n⁻ = 0`. [formula: ic-normal-sign] -/
theorem nPlus_mul_nMinus (n : ℤ) : nPlus n * nMinus n = 0 := by
  rcases le_or_gt 0 n with h | h
  · have : nMinus n = 0 := by simp only [nMinus]; omega
    rw [this, mul_zero]
  · have : nPlus n = 0 := by simp only [nPlus]; omega
    rw [this, zero_mul]

/-- `n⁺ ≥ 0` and `n⁻ ≥ 0`. [formula: ic-normal-sign] -/
theorem nPlus_nonneg (n : ℤ) : 0 ≤ nPlus n ∧ 0 ≤ nMinus n := by
  simp only [nPlus, nMinus]; omega

/-- `n⁻ = 1` exactly on a shared edge `n = −1`. [formula: ic-normal-sign] -/
theorem nMinus_eq_one_iff (n : ℤ) : nMinus n = 1 ↔ n = -1 := by
  simp only [nMinus]; omega

/-- [formula: ic-normal-sign] -/
theorem nMinus_eq_zero_iff (n : ℤ) : nMinus n = 0 ↔ 0 ≤ n := by
  simp only [nMinus]; omega

/-! ## Emanating counts and corner coordinates -/

/-- `e ≥ 0`. [formula: ic-emanating] -/
theorem emanating_nonneg (a b c : ℤ) : 0 ≤ emanating a b c := by
  simp only [emanating]; omega

/-- At most one corner of a face emanates: if `e^{ij}_k > 0` then `e^{jk}_i = e^{ki}_j = 0`.
[formula: ic-emanating] -/
theorem emanating_exclusive {a b c : ℤ} (h : 0 < emanating a b c) :
    emanating b c a = 0 ∧ emanating c a b = 0 := by
  simp only [emanating, nPlus] at *; omega

/-- `emanating` is symmetric in the two sides at the corner. [formula: ic-emanating] -/
theorem emanating_comm (a b c : ℤ) : emanating a b c = emanating a c b := by
  simp only [emanating]; omega

/-- The numerator of Eq. 1 is never negative. [formula: ic-corner-coordinates] -/
theorem cornerNumerator_nonneg (a b c : ℤ) : 0 ≤ cornerNumerator a b c := by
  simp only [cornerNumerator, emanating, nPlus]; omega

/-- `cornerNumerator` is symmetric in the two sides at the corner. [formula: ic-corner-coordinates] -/
theorem cornerNumerator_comm (a b c : ℤ) : cornerNumerator a b c = cornerNumerator a c b := by
  simp only [cornerNumerator, emanating, nPlus]; omega

/-- An emanating corner is cut by no arc: `e^{ij}_k > 0 → c^{ij}_k = 0`. [formula: ic-corner-coordinates] -/
theorem cornerNumerator_eq_zero_of_emanating {a b c : ℤ} (h : 0 < emanating a b c) :
    cornerNumerator a b c = 0 := by
  simp only [cornerNumerator, emanating, nPlus] at *; omega

/-- The checked form of Eq. 1 agrees with the plain one exactly when the numerator is even
(it is never negative). [formula: ic-corner-coordinates] -/
theorem cornerCoord?_eq (a b c : ℤ) :
    cornerCoord? a b c = if cornerNumerator a b c % 2 = 0 then some (cornerCoord a b c) else none := by
  have h := cornerNumerator_nonneg a b c
  simp only [cornerCoord?, cornerCoord]
  by_cases h2 : cornerNumerator a b c % 2 = 0
  · simp [h2, not_lt.mpr h]
  · simp [h2, not_lt.mpr h]

/-- **Eq. 1 and 2 on a face given by its arcs.** If the side counts (as `n⁺`) are
`n₀ = c₀ + c₁ + e₂`, `n₁ = c₁ + c₂ + e₀`, `n₂ = c₂ + c₀ + e₁` (side `k` carries the arcs that cut
its two corners and those from the opposite corner) with at most one `e > 0` and `c = 0` at an
emanating corner, then Eq. 1 at corner `0` is `2c₀` before halving and Eq. 2 at corner `0` is `e₀`;
by symmetry the same at every corner. [formula: ic-corner-coordinates] [formula: ic-emanating] -/
theorem corner_eval {n0 n1 n2 c0 c1 c2 e0 e1 e2 : ℤ} (hc0 : 0 ≤ c0) (hc1 : 0 ≤ c1) (hc2 : 0 ≤ c2)
    (he0 : 0 ≤ e0) (he1 : 0 ≤ e1) (he2 : 0 ≤ e2)
    (hx : (e0 = 0 ∧ e1 = 0) ∨ (e1 = 0 ∧ e2 = 0) ∨ (e2 = 0 ∧ e0 = 0))
    (hz0 : e0 = 0 ∨ c0 = 0) (hz1 : e1 = 0 ∨ c1 = 0) (hz2 : e2 = 0 ∨ c2 = 0)
    (h0 : nPlus n0 = c0 + c1 + e2) (h1 : nPlus n1 = c1 + c2 + e0) (h2 : nPlus n2 = c2 + c0 + e1) :
    cornerNumerator n1 n0 n2 = 2 * c0 ∧ emanating n1 n0 n2 = e0 ∧ cornerCoord n1 n0 n2 = c0 := by
  have hn : cornerNumerator n1 n0 n2 = 2 * c0 := by
    simp only [cornerNumerator, emanating]; rw [h0, h1, h2]; omega
  refine ⟨hn, ?_, ?_⟩
  · simp only [emanating]; rw [h0, h1, h2]; omega
  · simp only [cornerCoord]; rw [hn]; omega

/-- **Validity of a face** (SPEC §2.4 item 1): the three numerators of Eq. 1 are even if and only
if `n⁺₀ + n⁺₁ + n⁺₂` is even or one corner emanates. [formula: ic-face-validity] -/
theorem faceValid_iff (n0 n1 n2 : ℤ) :
    faceValid n0 n1 n2 = true ↔
      (nPlus n0 + nPlus n1 + nPlus n2) % 2 = 0 ∨ 0 < emanating n1 n0 n2 ∨
        0 < emanating n2 n1 n0 ∨ 0 < emanating n0 n2 n1 := by
  simp only [faceValid, cornerNumerator, emanating, nPlus, Bool.and_eq_true, beq_iff_eq]
  omega

/-- **A valid face has the form of `corner_eval`**, with `c` and `e` given by Eq. 1 and 2.
[formula: ic-face-validity] [formula: ic-corner-coordinates] -/
theorem exists_params_of_faceValid {n0 n1 n2 : ℤ} (h : faceValid n0 n1 n2 = true) :
    nPlus n0 = cornerCoord n1 n0 n2 + cornerCoord n2 n1 n0 + emanating n0 n2 n1 ∧
      nPlus n1 = cornerCoord n2 n1 n0 + cornerCoord n0 n2 n1 + emanating n1 n0 n2 ∧
      nPlus n2 = cornerCoord n0 n2 n1 + cornerCoord n1 n0 n2 + emanating n2 n1 n0 := by
  simp only [faceValid, cornerCoord, cornerNumerator, emanating, nPlus, Bool.and_eq_true,
    beq_iff_eq] at *
  omega

/-! ## Crossings -/

/-- Reversal of a combinatorial crossing is an involution. [formula: ic-reverse-crossing] -/
theorem reverseCrossing_reverseCrossing (n p : ℤ) : reverseCrossing n (reverseCrossing n p) = p := by
  simp only [reverseCrossing]; omega

/-- Reversal maps the crossings `0 ≤ p < n` of a halfedge onto those of its twin, reversing the
order. [formula: ic-reverse-crossing] -/
theorem reverseCrossing_mem {n p : ℤ} : (0 ≤ reverseCrossing n p ∧ reverseCrossing n p < n) ↔
    (0 ≤ p ∧ p < n) := by
  simp only [reverseCrossing]; omega

/-- [formula: ic-reverse-crossing] -/
theorem reverseCrossing_lt_iff (n p q : ℤ) : reverseCrossing n p < reverseCrossing n q ↔ q < p := by
  simp only [reverseCrossing]; omega

/-! ## Edge split along a shared edge -/

/-- **Eq. 6** counts every arc of the face `ijk` when the side `ij` carries no crossing: with
`n⁺_ij = 0` the face has `c_i = c_j = e_k = 0`, its arcs are the `c_k` arcs around `k` and the
`e_i + e_j` arcs from `i` or `j` (at most one of them nonzero), and `max(n_ki, n_jk, 0)` is their
number. Arguments: `n_ij`, `n_jk`, `n_ki` as `n₀, n₁, n₂`. [formula: ic-split-edge-normal] -/
theorem splitEdgeNormal_eq {n0 n1 n2 : ℤ} (hv : faceValid n0 n1 n2 = true) (h0 : nPlus n0 = 0) :
    splitEdgeNormal n2 n1 = cornerCoord n0 n2 n1 + emanating n1 n0 n2 + emanating n2 n1 n0 ∧
      cornerCoord n1 n0 n2 = 0 ∧ cornerCoord n2 n1 n0 = 0 ∧ emanating n0 n2 n1 = 0 := by
  simp only [faceValid, splitEdgeNormal, cornerCoord, cornerNumerator, emanating, nPlus,
    Bool.and_eq_true, beq_iff_eq] at *
  omega

end IntegerCoords

end DdgProofs
