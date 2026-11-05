import DdgProofs.Cuts.Quantile

/-!
# The feedback loop of the cone solver: the scale mismatch

`crates/ddg-optimal-cones/src/feedback.rs`, `scale_mismatch`. At a vertex `i` the solver's log
scale factor `u_i` is compared with the realized one: `R_i` is the area-weighted mean of the
realized log area factors `r_t` of the faces at `i` (`faceMean` of the pairs `(r_t, A_t)`), and

```
raw_i = u_i − R_i,      c = weighted median of raw over the measured vertices (weights m_j),
δ_i   = raw_i − c.
```

The weighted median is the weighted quantile at `q = ½` of the cut flow (`cut-weighted-quantile`,
`DdgProofs/Cuts/Quantile.lean`): the Rust function is the same scan.

* **`cut-weighted-quantile`** (used again): `isWeightedQuantile_shift`, `weightedQuantile_shift`
  — shifting every value by `a` shifts the quantile by `a` (nonnegative weights with a positive
  sum, `0 < q ≤ 1`).
* **`cone-feedback-mismatch`**:
  - `faceMean_shift` — adding `b` to every face value adds `b` to the mean (positive area sum);
  - `mismatch_shift` — `δ` does not change when a constant `a` is added to `u` (the free
    constant on a closed mesh) and a constant `b` to every realized face value (a uniform
    scale of the plane by `e^b`);
  - `mismatch_eq_zero_of_const` — if `u_j − R_j` is one constant at every measured vertex (the
    flattening realizes the solver's scale factor up to such constants), `δ = 0` there: such a
    flattening never reweights a vertex.

What is not proved: that `δ` measures how good the linearization is (that is the measurement in
the crate README), and anything about the threshold rule, which is a design choice
(`lean/proofs/CITED.md`, `cone-feedback-rule`).
-/

noncomputable section

namespace DDG

namespace ConeFeedback

open Cuts

/-! ## Shifted values -/

/-- Every value of a list of (value, weight) pairs shifted by `a`, the weights kept. -/
def shiftValues (a : ℝ) (xs : List (ℝ × ℝ)) : List (ℝ × ℝ) := xs.map fun p => (p.1 + a, p.2)

theorem shiftValues_cons (a : ℝ) (p : ℝ × ℝ) (l : List (ℝ × ℝ)) :
    shiftValues a (p :: l) = (p.1 + a, p.2) :: shiftValues a l := rfl

theorem totalWeight_shiftValues (a : ℝ) (xs : List (ℝ × ℝ)) :
    totalWeight (shiftValues a xs) = totalWeight xs := by
  induction xs with
  | nil => simp [shiftValues, totalWeight]
  | cons p l ih => simp only [shiftValues_cons, totalWeight_cons, ih]

theorem weightLE_shiftValues (a : ℝ) (xs : List (ℝ × ℝ)) (v : ℝ) :
    weightLE (shiftValues a xs) (v + a) = weightLE xs v := by
  induction xs with
  | nil => simp [shiftValues, weightLE]
  | cons p l ih => simp only [shiftValues_cons, weightLE_cons, ih, add_le_add_iff_right]

theorem mem_shiftValues {a : ℝ} {xs : List (ℝ × ℝ)} {p : ℝ × ℝ} (hp : p ∈ shiftValues a xs) :
    ∃ p₀ ∈ xs, p = (p₀.1 + a, p₀.2) := by
  obtain ⟨p₀, hp₀, rfl⟩ := List.mem_map.mp hp
  exact ⟨p₀, hp₀, rfl⟩

/-- **A shifted quantile.** If `v` is the `q`-quantile of the pairs, `v + a` is the `q`-quantile
of the pairs with every value shifted by `a`.

[formula: cut-weighted-quantile] -/
theorem isWeightedQuantile_shift {xs : List (ℝ × ℝ)} {q v : ℝ} (h : IsWeightedQuantile xs q v)
    (a : ℝ) : IsWeightedQuantile (shiftValues a xs) q (v + a) := by
  obtain ⟨⟨p, hp, hv⟩, hreach, hleast⟩ := h
  refine ⟨⟨(p.1 + a, p.2), List.mem_map.mpr ⟨p, hp, rfl⟩, by simp [hv]⟩, ?_, ?_⟩
  · rw [totalWeight_shiftValues, weightLE_shiftValues]
    exact hreach
  · intro p' hp' hlt
    obtain ⟨p₀, hp₀, rfl⟩ := mem_shiftValues hp'
    have hlt₀ : p₀.1 < v := by simpa using hlt
    show weightLE (shiftValues a xs) (p₀.1 + a) < q * totalWeight (shiftValues a xs)
    rw [totalWeight_shiftValues, weightLE_shiftValues]
    exact hleast p₀ hp₀ hlt₀

/-- **The weighted quantile is shift-equivariant**: for nonnegative weights with a positive sum
and `0 < q ≤ 1`, shifting every value by `a` shifts the quantile by `a`.

[formula: cut-weighted-quantile] -/
theorem weightedQuantile_shift {xs : List (ℝ × ℝ)} (hw : ∀ p ∈ xs, 0 ≤ p.2)
    (hW : 0 < totalWeight xs) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) (a : ℝ) :
    weightedQuantile (shiftValues a xs) q = weightedQuantile xs q + a := by
  have hw' : ∀ p ∈ shiftValues a xs, 0 ≤ p.2 := by
    intro p hp
    obtain ⟨p₀, hp₀, rfl⟩ := mem_shiftValues hp
    exact hw p₀ hp₀
  have hW' : 0 < totalWeight (shiftValues a xs) := by
    rw [totalWeight_shiftValues]
    exact hW
  exact (weightedQuantile_spec hw' hW' hq0 hq1).unique
    (isWeightedQuantile_shift (weightedQuantile_spec hw hW hq0 hq1) a)

/-! ## The realized value at a vertex -/

/-- `R_i`: the mean of the realized face values `r_t` weighted by the face areas `A_t`, for the
list of the pairs `(r_t, A_t)` of the faces at the vertex. -/
def faceMean (fs : List (ℝ × ℝ)) : ℝ := (fs.map fun p => p.2 * p.1).sum / totalWeight fs

theorem sum_mul_shiftValues (b : ℝ) (fs : List (ℝ × ℝ)) :
    ((shiftValues b fs).map fun p => p.2 * p.1).sum =
      (fs.map fun p => p.2 * p.1).sum + b * totalWeight fs := by
  induction fs with
  | nil => simp [shiftValues, totalWeight]
  | cons p l ih =>
    simp only [shiftValues_cons, List.map_cons, List.sum_cons, ih, totalWeight_cons]
    ring

/-- Adding `b` to every face value adds `b` to the area-weighted mean.

[formula: cone-feedback-mismatch] -/
theorem faceMean_shift (b : ℝ) {fs : List (ℝ × ℝ)} (h : totalWeight fs ≠ 0) :
    faceMean (shiftValues b fs) = faceMean fs + b := by
  unfold faceMean
  rw [sum_mul_shiftValues, totalWeight_shiftValues, add_div, mul_div_assoc, div_self h,
    mul_one]

/-! ## The mismatch -/

variable {ι : Type*}

/-- `u_i − R_i`: `nb i` lists the pairs `(r_t, A_t)` of the measured faces at `i`. -/
def raw (u : ι → ℝ) (nb : ι → List (ℝ × ℝ)) (i : ι) : ℝ := u i - faceMean (nb i)

/-- The pairs `(u_j − R_j, m_j)` of the measured vertices `vs` of a component. -/
def pairs (u : ι → ℝ) (nb : ι → List (ℝ × ℝ)) (m : ι → ℝ) (vs : List ι) : List (ℝ × ℝ) :=
  vs.map fun j => (raw u nb j, m j)

/-- The offset `c`: the weighted median of `u − R` over the measured vertices. -/
def offset (u : ι → ℝ) (nb : ι → List (ℝ × ℝ)) (m : ι → ℝ) (vs : List ι) : ℝ :=
  weightedQuantile (pairs u nb m vs) (1 / 2)

/-- `δ_i = (u_i − R_i) − c`. -/
def mismatch (u : ι → ℝ) (nb : ι → List (ℝ × ℝ)) (m : ι → ℝ) (vs : List ι) (i : ι) : ℝ :=
  raw u nb i - offset u nb m vs

section Shift

variable {u u' : ι → ℝ} {nb nb' : ι → List (ℝ × ℝ)} {a b : ℝ}

theorem raw_shift (hu : ∀ j, u' j = u j + a) (hnb : ∀ j, nb' j = shiftValues b (nb j)) {i : ι}
    (h : totalWeight (nb i) ≠ 0) : raw u' nb' i = raw u nb i + (a - b) := by
  simp only [raw]
  rw [hu i, hnb i, faceMean_shift b h]
  ring

theorem pairs_shift (m : ι → ℝ) (vs : List ι) (hu : ∀ j, u' j = u j + a)
    (hnb : ∀ j, nb' j = shiftValues b (nb j)) (h : ∀ j ∈ vs, totalWeight (nb j) ≠ 0) :
    pairs u' nb' m vs = shiftValues (a - b) (pairs u nb m vs) := by
  simp only [pairs, shiftValues, List.map_map]
  apply List.map_congr_left
  intro j hj
  simp only [Function.comp_apply, raw_shift hu hnb (h j hj)]

theorem offset_shift {m : ι → ℝ} {vs : List ι} (hu : ∀ j, u' j = u j + a)
    (hnb : ∀ j, nb' j = shiftValues b (nb j)) (h : ∀ j ∈ vs, totalWeight (nb j) ≠ 0)
    (hm : ∀ j ∈ vs, 0 ≤ m j) (hM : 0 < totalWeight (pairs u nb m vs)) :
    offset u' nb' m vs = offset u nb m vs + (a - b) := by
  have hw : ∀ p ∈ pairs u nb m vs, 0 ≤ p.2 := by
    intro p hp
    obtain ⟨j, hj, rfl⟩ := List.mem_map.mp hp
    exact hm j hj
  unfold offset
  rw [pairs_shift m vs hu hnb h]
  exact weightedQuantile_shift hw hM (by norm_num) (by norm_num) (a - b)

/-- **The mismatch does not depend on the constants.** Adding `a` to the solver's `u` (its free
constant on a closed mesh) and `b` to every realized face value (scaling the plane by `e^b`)
leaves `δ_i` unchanged, at every vertex whose measured faces have a positive area, provided
the masses are nonnegative with a positive sum over the measured vertices `vs` and every
measured vertex has measured faces of positive area.

[formula: cone-feedback-mismatch] -/
theorem mismatch_shift {m : ι → ℝ} {vs : List ι} (hu : ∀ j, u' j = u j + a)
    (hnb : ∀ j, nb' j = shiftValues b (nb j)) (h : ∀ j ∈ vs, totalWeight (nb j) ≠ 0)
    (hm : ∀ j ∈ vs, 0 ≤ m j) (hM : 0 < totalWeight (pairs u nb m vs)) {i : ι}
    (hi : totalWeight (nb i) ≠ 0) : mismatch u' nb' m vs i = mismatch u nb m vs i := by
  unfold mismatch
  rw [raw_shift hu hnb hi, offset_shift hu hnb h hm hM]
  ring

end Shift

/-- **A flattening that realizes the solver's scale factor up to constants has no mismatch.**
If `u_j − R_j = κ` at every measured vertex `j ∈ vs` (nonnegative masses with a positive sum),
then `δ_i = 0` at every measured vertex: the rule never reweights one.

[formula: cone-feedback-mismatch] -/
theorem mismatch_eq_zero_of_const {u : ι → ℝ} {nb : ι → List (ℝ × ℝ)} {m : ι → ℝ}
    {vs : List ι} {κ : ℝ} (hκ : ∀ j ∈ vs, raw u nb j = κ) (hm : ∀ j ∈ vs, 0 ≤ m j)
    (hM : 0 < totalWeight (pairs u nb m vs)) {i : ι} (hi : i ∈ vs) :
    mismatch u nb m vs i = 0 := by
  have hw : ∀ p ∈ pairs u nb m vs, 0 ≤ p.2 := by
    intro p hp
    obtain ⟨j, hj, rfl⟩ := List.mem_map.mp hp
    exact hm j hj
  have hq := weightedQuantile_spec hw hM (q := 1 / 2) (by norm_num) (by norm_num)
  obtain ⟨p, hp, hv⟩ := hq.mem
  obtain ⟨j, hj, rfl⟩ := List.mem_map.mp hp
  unfold mismatch offset
  rw [hv]
  show raw u nb i - raw u nb j = 0
  rw [hκ i hi, hκ j hj, sub_self]

end ConeFeedback

end DDG
