import DdgProofs.PointsAsTori.Basic

/-!
# Points as Tori: neighborhood scale, local frame and features (SPEC §4.1)

For a neighborhood `P 0 = p_i, P 1, …, P k` with normals `N`:

```
σ_i = median_m |P m − P 0|                         (all k+1 values, m = 0 included)
n̂   = n / max(|n|, 1e-10)
u   = e_y if |c_x| > 0.9 else e_x,   ŝ = (c × u)/|c × u|,   t̂ = c × ŝ      (c = n̂_i)
f_m = (d·ŝ, d·t̂, d·c, n̂_m·ŝ, n̂_m·t̂, n̂_m·c),   d = (P m − P 0)/σ_i
```

* `median` is the middle element of the sorted list, or the mean of the two middle elements for
  an even length; it does not depend on the sorting algorithm (`median_eq_of_sorted`) and is
  positively homogeneous (`median_map_mul`). So `σ` is translation invariant and scales with the
  cloud (`sigma_affine`).
* The frame rule gives an orthonormal right-handed frame for every unit `c`, on both branches of
  the switch (`refFrame`).
* The features are invariant under `p ↦ a p + τ` (`a > 0`) of the positions, for **any** frame
  rule depending only on `c` (`features_affine`). (They are not rotation invariant: the switch
  depends on the coordinates of `c`.)
-/

noncomputable section

namespace DDG

namespace PAT

open Matrix Real
open RC (V3 vnorm vnorm_nonneg vnorm_sq vnorm_mul_self vnorm_eq_zero_iff vnorm_pos_iff vnorm_smul
  dot_self_nonneg)

/-! ## The median -/

/-- Sorting (any sorting algorithm gives the same list, `sortR_eq`). -/
def sortR (l : List ℝ) : List ℝ := l.insertionSort (· ≤ ·)

/-- The median: the middle element of the sorted list for an odd length, the mean of the two
middle elements for an even length. -/
def median (l : List ℝ) : ℝ :=
  if l.length % 2 = 0 then
    ((sortR l).getD (l.length / 2 - 1) 0 + (sortR l).getD (l.length / 2) 0) / 2
  else (sortR l).getD (l.length / 2) 0

/-- The sorted list is the unique sorted permutation. -/
theorem sortR_eq {l s : List ℝ} (hs : s.Pairwise (· ≤ ·)) (hp : s.Perm l) : sortR l = s :=
  List.Perm.eq_of_pairwise' (List.pairwise_insertionSort _ _) hs
    ((List.perm_insertionSort _ l).trans hp.symm)

/-- **The median of an even list is the mean of the middle two** of *any* sorted permutation
`s` (for `k = 64` neighbors, the 32nd and 33rd smallest distances).

[formula: pat-median-scale] -/
theorem median_eq_of_sorted {l s : List ℝ} (hs : s.Pairwise (· ≤ ·)) (hp : s.Perm l)
    {m : ℕ} (hm : l.length = 2 * m) (hm0 : 0 < m) :
    median l = (s[m - 1]'(by rw [hp.length_eq]; omega) + s[m]'(by rw [hp.length_eq]; omega)) / 2 := by
  have hlen : s.length = 2 * m := by rw [hp.length_eq, hm]
  simp only [median, sortR_eq hs hp, hm, Nat.mul_mod_right, ite_true,
    Nat.mul_div_cancel_left m (by norm_num : 0 < 2)]
  rw [List.getD_eq_getElem _ _ (by omega), List.getD_eq_getElem _ _ (by omega)]

theorem getD_map_mul (c : ℝ) (l : List ℝ) (i : ℕ) :
    (l.map (c * ·)).getD i 0 = c * l.getD i 0 := by
  simp only [List.getD_eq_getElem?_getD, List.getElem?_map]
  cases l[i]? <;> simp

/-- **The median is positively homogeneous**: `median(c·l) = c·median(l)` for `c > 0`.

[formula: pat-median-scale] -/
theorem median_map_mul {c : ℝ} (hc : 0 < c) (l : List ℝ) :
    median (l.map (c * ·)) = c * median l := by
  have hsort : sortR (l.map (c * ·)) = (sortR l).map (c * ·) := by
    unfold sortR
    refine (List.map_insertionSort (r := (· ≤ ·)) (s := (· ≤ ·)) (c * ·) l
      (fun a _ b _ => ?_)).symm
    exact ⟨fun h => mul_le_mul_of_nonneg_left h hc.le, fun h => le_of_mul_le_mul_left h hc⟩
  simp only [median, List.length_map, hsort, getD_map_mul]
  split_ifs <;> ring

/-- The scale `σ = median_m |P m − P 0|` of a neighborhood (center first). -/
def sigma {k : ℕ} (P : Fin (k + 1) → V3) : ℝ := median (List.ofFn fun m => vnorm (P m - P 0))

/-- **`σ` is translation invariant and scales with the cloud**: `σ(aP + τ) = a σ(P)`, `a > 0`.

[formula: pat-median-scale] -/
theorem sigma_affine {k : ℕ} (P : Fin (k + 1) → V3) {a : ℝ} (ha : 0 < a) (τ : V3) :
    sigma (fun m => a • P m + τ) = a * sigma P := by
  unfold sigma
  have : (List.ofFn fun m => vnorm ((a • P m + τ) - (a • P 0 + τ))) =
      (List.ofFn fun m => vnorm (P m - P 0)).map (a * ·) := by
    rw [List.map_ofFn]
    congr 1
    funext m
    simp only [Function.comp, add_sub_add_right_eq_sub, ← smul_sub, vnorm_smul, abs_of_pos ha]
  rw [this, median_map_mul ha]

/-! ## The local frame -/

/-- The helper axis: `e_y` if `|c_x| > 0.9`, else `e_x`. -/
def pickU (c : V3) : V3 := if 0.9 < |c 0| then ![0, 1, 0] else ![1, 0, 0]

/-- `ŝ = (c × u)/|c × u|`. -/
def frameS (c : V3) : V3 := (vnorm (c ⨯₃ pickU c))⁻¹ • (c ⨯₃ pickU c)

/-- `t̂ = c × ŝ`. -/
def frameT (c : V3) : V3 := c ⨯₃ frameS c

theorem cross_pickU_dot (c : V3) (hc : c ⬝ᵥ c = 1) :
    0 < (c ⨯₃ pickU c) ⬝ᵥ (c ⨯₃ pickU c) := by
  have h1 : c 0 ^ 2 + c 1 ^ 2 + c 2 ^ 2 = 1 := by
    simpa [dotProduct, Fin.sum_univ_three, sq] using hc
  unfold pickU
  split_ifs with h
  · have h0 : 0.81 < c 0 ^ 2 := by
      have := sq_lt_sq' (by linarith [abs_nonneg (c 0)]) h
      rw [sq_abs] at this; norm_num at this ⊢; linarith
    simp only [cross_apply, dotProduct, Fin.sum_univ_three]
    simp
    nlinarith [sq_nonneg (c 2)]
  · push Not at h
    have h0 : c 0 ^ 2 ≤ 0.81 := by
      have := sq_le_sq' (by linarith [abs_nonneg (c 0), neg_abs_le (c 0)]) h
      rw [sq_abs] at this; norm_num at this ⊢; linarith
    simp only [cross_apply, dotProduct, Fin.sum_univ_three]
    simp
    nlinarith [sq_nonneg (c 2), sq_nonneg (c 1)]

/-- Any unit `c` and unit `s ⊥ c` complete to a right-handed frame `(s, c × s, c)`. -/
def frameOf (c s : V3) (hc : c ⬝ᵥ c = 1) (hs : s ⬝ᵥ s = 1) (hsc : s ⬝ᵥ c = 0) : Frame where
  s := s
  t := c ⨯₃ s
  n := c
  ss := hs
  tt := by rw [cross_dot_cross, hc, hs, dotProduct_comm c s, hsc]; ring
  nn := hc
  st := by rw [dot_cross_self]
  sn := hsc
  tn := by rw [dotProduct_comm, dot_self_cross]
  cross := by rw [cross_cross_eq_smul_sub_smul', hs, dotProduct_comm c s, hsc]; simp

/-- **The reference frame rule gives an orthonormal right-handed frame `(ŝ, t̂, c)`** for every
unit `c`, on both branches of the switch `|c_x| > 0.9`: `|c × u| > 0`, `ŝ·ŝ = t̂·t̂ = 1`,
`ŝ·t̂ = ŝ·c = t̂·c = 0` and `ŝ × t̂ = c` (the fields of `Frame`).

[formula: pat-local-frame] -/
theorem refFrame_spec (c : V3) (hc : c ⬝ᵥ c = 1) :
    0 < vnorm (c ⨯₃ pickU c) ∧ frameS c ⬝ᵥ frameS c = 1 ∧ frameS c ⬝ᵥ c = 0 ∧
      frameT c ⬝ᵥ frameT c = 1 ∧ frameS c ⬝ᵥ frameT c = 0 ∧ frameT c ⬝ᵥ c = 0 ∧
      frameS c ⨯₃ frameT c = c := by
  have hpos := cross_pickU_dot c hc
  have hv : 0 < vnorm (c ⨯₃ pickU c) := Real.sqrt_pos.mpr hpos
  have hss : frameS c ⬝ᵥ frameS c = 1 := by
    have := RC.unit_dot_self (v := c ⨯₃ pickU c) (fun h => by rw [h] at hpos; simp at hpos)
    exact this
  have hsc : frameS c ⬝ᵥ c = 0 := by
    rw [dotProduct_comm]
    simp only [frameS, dotProduct_smul, dot_self_cross, smul_zero]
  let F := frameOf c (frameS c) hc hss hsc
  exact ⟨hv, hss, hsc, F.tt, F.st, F.tn, F.cross⟩

/-- The reference frame as a `Frame`. -/
def refFrame (c : V3) (hc : c ⬝ᵥ c = 1) : Frame :=
  frameOf c (frameS c) hc (refFrame_spec c hc).2.1 (refFrame_spec c hc).2.2.1

/-! ## Features -/

/-- The guarded unit normal `n / max(|n|, 1e-10)`. -/
def unitGuard (v : V3) : V3 := (max (vnorm v) 1e-10)⁻¹ • v

/-- The feature vector of token `m` for a frame rule `frame : c ↦ (ŝ, t̂)`. -/
def features {k : ℕ} (frame : V3 → V3 × V3) (P N : Fin (k + 1) → V3) (m : Fin (k + 1)) :
    Fin 6 → ℝ :=
  let c := unitGuard (N 0)
  let d := (sigma P)⁻¹ • (P m - P 0)
  let nm := unitGuard (N m)
  ![d ⬝ᵥ (frame c).1, d ⬝ᵥ (frame c).2, d ⬝ᵥ c, nm ⬝ᵥ (frame c).1, nm ⬝ᵥ (frame c).2, nm ⬝ᵥ c]

/-- The reference rule `c ↦ (ŝ, t̂)`. -/
def refRule (c : V3) : V3 × V3 := (frameS c, frameT c)

/-- **The features are invariant under translation and uniform scaling of the cloud**
(`p ↦ a p + τ`, `a > 0`, normals unchanged), for any frame rule that depends only on `c`, in
particular the reference rule `refRule`. With `σ = 0` (all points equal) the features are `0`
on both sides (Lean's `0⁻¹ = 0`), so no hypothesis is needed.

[formula: pat-features] [formula: pat-feature-invariance] -/
theorem features_affine {k : ℕ} (frame : V3 → V3 × V3) (P N : Fin (k + 1) → V3) {a : ℝ}
    (ha : 0 < a) (τ : V3) (m : Fin (k + 1)) :
    features frame (fun j => a • P j + τ) N m = features frame P N m := by
  have hd : (sigma fun j => a • P j + τ)⁻¹ • ((a • P m + τ) - (a • P 0 + τ)) =
      (sigma P)⁻¹ • (P m - P 0) := by
    rw [sigma_affine P ha τ, add_sub_add_right_eq_sub, ← smul_sub, smul_smul]
    congr 1
    rw [mul_inv]
    field_simp
  simp only [features, hd]

end PAT

end DDG
