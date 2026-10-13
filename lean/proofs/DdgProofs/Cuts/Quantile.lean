import DdgProofs.Cuts.Ripple

/-!
# Cut flow: the weighted quantile, the robust step rule and the local limit

`crates/ddg-variational-cuts/SPEC.md` §1 ("the roughening rate of the state"), §4 ("Step size"),
D22, D23.

* **`cut-weighted-quantile`** — `quantile_q(x; ℓ)`: with the values in ascending order, the
  first value at which the weights summed so far reach `q · Σℓ`. `weightedQuantile xs q` is this
  for the list `xs` of the pairs `(x_c, ℓ_c)`: an insertion sort by value (stable; the order of
  equal values does not affect the result) followed by the scan `quantileScan`. For nonnegative
  weights with a positive sum and `0 < q ≤ 1` (`oracle_weightedQuantile_spec`,
  `IsWeightedQuantile`): the result is one of the values; the weight of the values not exceeding
  it is at least `q · Σℓ`; every smaller value has less than that below or at it — so it is the
  least value with this property, and less than `q · Σℓ` lies strictly below it
  (`oracle_weightLT_weightedQuantile_lt`). Hence it lies between the smallest and the largest
  value, does not decrease in `q`, and (positive weights) is the largest value for `q = 1`
  (`oracle_weightedQuantile_bounds`, `oracle_weightedQuantile_mono`,
  `oracle_weightedQuantile_one`); an empty list gives `0` (`oracle_weightedQuantile_nil`).
* **`cut-roughening-rate`** — the rate of the state is `ω = weightedQuantile 0.9 (ω_c) (ℓ_c)`: it
  is one of the local rates, at least nine tenths of the length of the cut has `ω_c ≤ ω`, and
  at such a cut point every ripple shorter than `2πw/ω` decays, that is shorter than `ℓ` when
  `α_L = ωℓ/(2π)`, `ω > 0` (`roughening_quantile_spec`, `rippleRate_neg_of_le_of_lt`,
  `rippleRate_neg_of_le`). Ripples longer than `ℓ` grow where `ω_c ≥ ω`, which is more than a
  tenth of the cut; where `ω_c < ω` they decay down to the wavelength `ℓ ω/ω_c`
  (`rippleRate_pos_of_ge`, `roughening_quantile_above`, `rippleRate_neg_of_lt_longer`).
  `ω ≥ 0` if all `ω_c ≥ 0` (the Dirichlet energy, `roughening_quantile_nonneg`); for `ω ≤ 0` the
  weight rule gives no positive weight, and where `ω_c < 0` every ripple decays
  (`lengthWeightForRipples_nonpos`, `rippleRate_neg_of_neg`).
* **`cut-step-cap`** — the oracle's `robustStep` (`= adaptiveStep` with `h = 1` and the speed
  `m̂`) and `adaptiveStep` skip a term whose denominator is not positive, satisfy the limits that
  apply and are the largest such step (`robustStep_spec`, `adaptiveStep_skip_*`,
  `le_adaptiveStep`); with `ω ≤ 0` they are the capped rule (`robustStep_capped`,
  `adaptiveStep_skip_omega`). With the robust statistics every cut point with
  `speedInEdges ρ_c h_c ≤ m̂` moves by at most `max_motion · h_c` in the step
  (`robustStep_motion_le'`); with the extremes, `m̂ = max|ρ_c|/min h`, every cut point does
  (`robustStep_extremes`, `extremesStep_motion_le`). The larger of two `m̂` gives the smaller of
  the two steps, and the step respects every speed up to `m̂` (`robustStep_max`,
  `robustStep_mul_le_of_le`).
* **`cut-step-local`** — the oracle's `localLimitFactor max_motion h_c τ v_c`
  (`min(1, max_motion · h_c/(τ v_c))` where `τ v_c > 0`, else `1`) lies in `(0, 1]`, is `1` where
  `τ v_c ≤ max_motion · h_c`, and with `v_c ≥ |r_c|` the limited rate moves the cut point by at
  most `max_motion · h_c`; a limited point takes the step `max_motion · h_c/v_c`; the limit keeps
  the sign of the rate, never increases its size, and keeps a vanishing rate
  (`localLimitFactor_mem_Ioc`, `localLimitFactor_eq_one`, `localLimit_motion_le`,
  `localLimit_step_eq`, `localLimit_sign`); a rate limited point by point is still a descent
  direction (`firstVariation_localLimit`); the oracle's
  `localSpeed r_c ω_c ω = max(|r_c|, ½(ω_c − ω)₊)` is `≥ |r_c|` (`abs_le_localSpeed`,
  `localLimit_motion_le_adaptive`), and `|s_c| ≤ ω_c/2` when `Q_0, Q_1 ≥ 0`
  (`abs_cutSpeed_le_half_roughRate`).

The tagged theorems are about the oracle's definitions `Ddg.weightedQuantile` (with its
`sortAscending`, `insertAscending`, `sumLeft`, `firstReaching`), `robustStep`, `adaptiveStep`,
`stepCap`, `speedInEdges`, `localSpeed`, `localLimitFactor`, `roughRate`, `cutSpeed`,
`lengthWeightForRipples`. This package's own `quantileScan`, `Cuts.weightedQuantile` (Mathlib's
insertion sort and the same scan), `totalWeight`, `weightLE`, `weightLT`, `IsWeightedQuantile` and
`localFactor` are the vehicle of the proofs: `IsWeightedQuantile` characterizes the quantile
without an order, both quantiles satisfy it and are therefore equal
(`oracle_weightedQuantile_eq`). `sortAscending_reverses_ties` records that the oracle's sort
puts equal values in the reverse of their original order.
-/

noncomputable section

namespace DdgProofs

namespace Cuts

open scoped DdgProofs.Oracle
open Real

/-! ## Weights of a list of (value, weight) pairs -/

/-- The sum of the weights. -/
def totalWeight (l : List (ℝ × ℝ)) : ℝ := (l.map Prod.snd).sum

/-- The sum of the weights of the entries whose value does not exceed `v`. -/
def weightLE (l : List (ℝ × ℝ)) (v : ℝ) : ℝ := (l.map fun p => if p.1 ≤ v then p.2 else 0).sum

theorem totalWeight_nil : totalWeight [] = 0 := by simp [totalWeight]

theorem totalWeight_cons (p : ℝ × ℝ) (l : List (ℝ × ℝ)) :
    totalWeight (p :: l) = p.2 + totalWeight l := by simp [totalWeight]

theorem totalWeight_append (l₁ l₂ : List (ℝ × ℝ)) :
    totalWeight (l₁ ++ l₂) = totalWeight l₁ + totalWeight l₂ := by simp [totalWeight]

theorem weightLE_nil (v : ℝ) : weightLE [] v = 0 := by simp [weightLE]

theorem weightLE_cons (p : ℝ × ℝ) (l : List (ℝ × ℝ)) (v : ℝ) :
    weightLE (p :: l) v = (if p.1 ≤ v then p.2 else 0) + weightLE l v := by simp [weightLE]

theorem weightLE_append (l₁ l₂ : List (ℝ × ℝ)) (v : ℝ) :
    weightLE (l₁ ++ l₂) v = weightLE l₁ v + weightLE l₂ v := by simp [weightLE]

theorem totalWeight_perm {l l' : List (ℝ × ℝ)} (h : l.Perm l') : totalWeight l = totalWeight l' :=
  (h.map _).sum_eq

theorem weightLE_perm {l l' : List (ℝ × ℝ)} (h : l.Perm l') (v : ℝ) :
    weightLE l v = weightLE l' v :=
  (h.map _).sum_eq

theorem totalWeight_nonneg {l : List (ℝ × ℝ)} (hw : ∀ p ∈ l, 0 ≤ p.2) : 0 ≤ totalWeight l := by
  induction l with
  | nil => simp [totalWeight_nil]
  | cons p l ih =>
    rw [totalWeight_cons]
    have h1 := hw p (List.mem_cons.mpr (Or.inl rfl))
    have h2 := ih fun q hq => hw q (List.mem_cons.mpr (Or.inr hq))
    linarith

theorem totalWeight_pos {l : List (ℝ × ℝ)} (hne : l ≠ []) (hw : ∀ p ∈ l, 0 < p.2) :
    0 < totalWeight l := by
  cases l with
  | nil => exact absurd rfl hne
  | cons p l =>
    rw [totalWeight_cons]
    have h1 := hw p (List.mem_cons.mpr (Or.inl rfl))
    have h2 := totalWeight_nonneg fun q hq => (hw q (List.mem_cons.mpr (Or.inr hq))).le
    linarith

theorem weightLE_nonneg {l : List (ℝ × ℝ)} (hw : ∀ p ∈ l, 0 ≤ p.2) (v : ℝ) :
    0 ≤ weightLE l v := by
  induction l with
  | nil => simp [weightLE_nil]
  | cons p l ih =>
    rw [weightLE_cons]
    have h1 := hw p (List.mem_cons.mpr (Or.inl rfl))
    have h2 := ih fun q hq => hw q (List.mem_cons.mpr (Or.inr hq))
    split_ifs <;> linarith

theorem weightLE_le_totalWeight {l : List (ℝ × ℝ)} (hw : ∀ p ∈ l, 0 ≤ p.2) (v : ℝ) :
    weightLE l v ≤ totalWeight l := by
  induction l with
  | nil => simp [weightLE_nil, totalWeight_nil]
  | cons p l ih =>
    rw [weightLE_cons, totalWeight_cons]
    have h1 := hw p (List.mem_cons.mpr (Or.inl rfl))
    have h2 := ih fun q hq => hw q (List.mem_cons.mpr (Or.inr hq))
    split_ifs <;> linarith

theorem weightLE_eq_totalWeight {l : List (ℝ × ℝ)} {v : ℝ} (h : ∀ p ∈ l, p.1 ≤ v) :
    weightLE l v = totalWeight l := by
  induction l with
  | nil => simp [weightLE_nil, totalWeight_nil]
  | cons p l ih =>
    rw [weightLE_cons, totalWeight_cons, ite_eq_left (h p (List.mem_cons.mpr (Or.inl rfl))),
      ih fun q hq => h q (List.mem_cons.mpr (Or.inr hq))]

theorem weightLE_eq_zero {l : List (ℝ × ℝ)} {v : ℝ} (h : ∀ p ∈ l, v < p.1) : weightLE l v = 0 := by
  induction l with
  | nil => simp [weightLE_nil]
  | cons p l ih =>
    rw [weightLE_cons, ite_eq_right (not_le.mpr (h p (List.mem_cons.mpr (Or.inl rfl)))),
      ih fun q hq => h q (List.mem_cons.mpr (Or.inr hq)), add_zero]

/-- With positive weights, if the entries up to `v` carry all the weight, every value is `≤ v`. -/
theorem forall_le_of_weightLE_eq {l : List (ℝ × ℝ)} (hw : ∀ p ∈ l, 0 < p.2) {v : ℝ}
    (h : weightLE l v = totalWeight l) : ∀ p ∈ l, p.1 ≤ v := by
  induction l with
  | nil => intro p hp; simp at hp
  | cons a l ih =>
    rw [weightLE_cons, totalWeight_cons] at h
    have ha := hw a (List.mem_cons.mpr (Or.inl rfl))
    have hwl : ∀ q ∈ l, 0 < q.2 := fun q hq => hw q (List.mem_cons.mpr (Or.inr hq))
    have hle := weightLE_le_totalWeight (fun q hq => (hwl q hq).le) v
    by_cases hav : a.1 ≤ v
    · rw [ite_eq_left hav] at h
      intro p hp
      rcases List.mem_cons.mp hp with rfl | hp
      · exact hav
      · exact ih hwl (by linarith) p hp
    · rw [ite_eq_right hav] at h
      linarith

/-! ## The scan and the quantile -/

/-- The scan of SPEC §1 over a list of (value, weight) pairs in ascending order of the values:
the first value at which the weights summed so far (`seen` before the list) reach `target`; if
none does, the last value (`last` for an empty list). -/
def quantileScan : List (ℝ × ℝ) → ℝ → ℝ → ℝ → ℝ
  | [], _, _, last => last
  | p :: rest, seen, target, _ =>
    if target ≤ seen + p.2 then p.1 else quantileScan rest (seen + p.2) target p.1

theorem quantileScan_nil (seen target last : ℝ) : quantileScan [] seen target last = last := rfl

theorem quantileScan_cons (p : ℝ × ℝ) (rest : List (ℝ × ℝ)) (seen target last : ℝ) :
    quantileScan (p :: rest) seen target last =
      if target ≤ seen + p.2 then p.1 else quantileScan rest (seen + p.2) target p.1 := rfl

/-- The scan returns one of the values of the list, or the fallback. -/
theorem quantileScan_mem_or : ∀ (l : List (ℝ × ℝ)) (seen target last : ℝ),
    quantileScan l seen target last = last ∨ ∃ p ∈ l, quantileScan l seen target last = p.1
  | [], _, _, _ => Or.inl rfl
  | p :: rest, seen, target, last => by
    right
    rw [quantileScan_cons]
    by_cases h : target ≤ seen + p.2
    · exact ⟨p, List.mem_cons.mpr (Or.inl rfl), by rw [ite_eq_left h]⟩
    · rw [ite_eq_right h]
      rcases quantileScan_mem_or rest (seen + p.2) target p.1 with h1 | ⟨p', hp', h1⟩
      · exact ⟨p, List.mem_cons.mpr (Or.inl rfl), h1⟩
      · exact ⟨p', List.mem_cons.mpr (Or.inr hp'), h1⟩

/-- Where the scan stops: the list splits as `pre ++ p :: post` with the result `p.1`, the
weights of `pre` not reaching the target and those of `pre ++ [p]` reaching it. -/
theorem quantileScan_decomp : ∀ (l : List (ℝ × ℝ)) (seen target last : ℝ), seen < target →
    target ≤ seen + totalWeight l →
    ∃ (pre : List (ℝ × ℝ)) (p : ℝ × ℝ) (post : List (ℝ × ℝ)), l = pre ++ p :: post ∧
      quantileScan l seen target last = p.1 ∧ seen + totalWeight pre < target ∧
      target ≤ seen + totalWeight pre + p.2
  | [], seen, target, _, h1, h2 => by
    rw [totalWeight_nil, add_zero] at h2
    exact absurd h1 (not_lt.mpr h2)
  | p :: rest, seen, target, last, h1, h2 => by
    by_cases h : target ≤ seen + p.2
    · refine ⟨[], p, rest, rfl, by rw [quantileScan_cons, ite_eq_left h], ?_, ?_⟩
      · rw [totalWeight_nil, add_zero]
        exact h1
      · rw [totalWeight_nil, add_zero]
        exact h
    · rw [totalWeight_cons] at h2
      obtain ⟨pre, p', post, hl, hs, hlt, hle⟩ :=
        quantileScan_decomp rest (seen + p.2) target p.1 (not_le.mp h) (by linarith)
      refine ⟨p :: pre, p', post, by rw [hl]; rfl, ?_, ?_, ?_⟩
      · rw [quantileScan_cons, ite_eq_right h]
        exact hs
      · rw [totalWeight_cons]
        linarith
      · rw [totalWeight_cons]
        linarith

/-- **The scan on an ascending list.** For nonnegative weights and `0 < target ≤ Σ weights`: the
result is a value of the list, the weight of the values not exceeding it reaches the target, and
no smaller value does. -/
theorem quantileScan_sorted {l : List (ℝ × ℝ)} (hs : l.Pairwise fun a b => a.1 ≤ b.1)
    (hw : ∀ p ∈ l, 0 ≤ p.2) {target : ℝ} (h0 : 0 < target) (h1 : target ≤ totalWeight l)
    (last : ℝ) :
    (∃ p ∈ l, quantileScan l 0 target last = p.1) ∧
      target ≤ weightLE l (quantileScan l 0 target last) ∧
      ∀ p ∈ l, p.1 < quantileScan l 0 target last → weightLE l p.1 < target := by
  obtain ⟨pre, p, post, hl, hscan, hlt, hle⟩ :=
    quantileScan_decomp l 0 target last h0 (by rw [zero_add]; exact h1)
  rw [zero_add] at hlt hle
  subst hl
  rw [hscan]
  rw [List.pairwise_append] at hs
  obtain ⟨-, hpost, hcross⟩ := hs
  rw [List.pairwise_cons] at hpost
  have hpre_le : ∀ a ∈ pre, a.1 ≤ p.1 := fun a ha =>
    hcross a ha p (List.mem_cons.mpr (Or.inl rfl))
  have hpost_ge : ∀ b ∈ post, p.1 ≤ b.1 := hpost.1
  have hwpre : ∀ a ∈ pre, 0 ≤ a.2 := fun a ha => hw a (List.mem_append.mpr (Or.inl ha))
  have hwpost : ∀ b ∈ post, 0 ≤ b.2 := fun b hb =>
    hw b (List.mem_append.mpr (Or.inr (List.mem_cons.mpr (Or.inr hb))))
  refine ⟨⟨p, List.mem_append.mpr (Or.inr (List.mem_cons.mpr (Or.inl rfl))), rfl⟩, ?_, ?_⟩
  · rw [weightLE_append, weightLE_cons, ite_eq_left (le_refl p.1), weightLE_eq_totalWeight hpre_le]
    have := weightLE_nonneg hwpost p.1
    linarith
  · intro a _ hlt'
    rw [weightLE_append, weightLE_cons, ite_eq_right (not_le.mpr hlt'),
      weightLE_eq_zero fun b hb => lt_of_lt_of_le hlt' (hpost_ge b hb)]
    have := weightLE_le_totalWeight hwpre a.1
    linarith

/-- `v` is the `q`-quantile of the (value, weight) pairs `xs` in the sense of SPEC §1, stated
without reference to an order: `v` is one of the values, the weight of the values not exceeding
it is at least `q · Σℓ`, and every value below it has less than that at or below itself. -/
structure IsWeightedQuantile (xs : List (ℝ × ℝ)) (q v : ℝ) : Prop where
  mem : ∃ p ∈ xs, v = p.1
  reach : q * totalWeight xs ≤ weightLE xs v
  least : ∀ p ∈ xs, p.1 < v → weightLE xs p.1 < q * totalWeight xs

/-- There is at most one such value. -/
theorem IsWeightedQuantile.unique {xs : List (ℝ × ℝ)} {q v v' : ℝ}
    (h : IsWeightedQuantile xs q v) (h' : IsWeightedQuantile xs q v') : v = v' := by
  obtain ⟨p, hp, rfl⟩ := h.mem
  obtain ⟨p', hp', rfl⟩ := h'.mem
  apply le_antisymm
  · by_contra hcon
    have := h.least p' hp' (not_le.mp hcon)
    linarith [h'.reach]
  · by_contra hcon
    have := h'.least p hp (not_le.mp hcon)
    linarith [h.reach]

/-- The total weight of the values strictly below `v`. -/
def weightLT (l : List (ℝ × ℝ)) (v : ℝ) : ℝ := (l.map fun p => if p.1 < v then p.2 else 0).sum

/-- **Less than the fraction `q` of the weight lies strictly below the quantile**, so more than
the fraction `1 − q` lies at or above it. -/
theorem IsWeightedQuantile.weightLT_lt {xs : List (ℝ × ℝ)} {q v : ℝ}
    (h : IsWeightedQuantile xs q v) (hq : 0 < q) (hW : 0 < totalWeight xs) :
    weightLT xs v < q * totalWeight xs := by
  classical
  by_cases hne : ∃ p ∈ xs, p.1 < v
  · have hs : (xs.toFinset.filter fun p => p.1 < v).Nonempty := by
      obtain ⟨p, hp, hlt⟩ := hne
      exact ⟨p, Finset.mem_filter.mpr ⟨List.mem_toFinset.mpr hp, hlt⟩⟩
    obtain ⟨m, hm, hmax⟩ := Finset.exists_max_image _ (fun p : ℝ × ℝ => p.1) hs
    rw [Finset.mem_filter, List.mem_toFinset] at hm
    have e : weightLT xs v = weightLE xs m.1 := by
      unfold weightLT weightLE
      congr 1
      apply List.map_congr_left
      intro p hp
      by_cases hlt : p.1 < v
      · have hle : p.1 ≤ m.1 := hmax p (Finset.mem_filter.mpr ⟨List.mem_toFinset.mpr hp, hlt⟩)
        rw [ite_eq_left hlt, ite_eq_left hle]
      · have hnle : ¬ p.1 ≤ m.1 := fun hle => hlt (lt_of_le_of_lt hle hm.2)
        rw [ite_eq_right hlt, ite_eq_right hnle]
    rw [e]
    exact h.least m hm.1 hm.2
  · have e : weightLT xs v = 0 := by
      unfold weightLT
      apply List.sum_eq_zero
      intro x hx
      obtain ⟨p, hp, rfl⟩ := List.mem_map.mp hx
      exact ite_eq_right fun hlt => hne ⟨p, hp, hlt⟩
    rw [e]
    exact mul_pos hq hW

/-- **The scan of any ascending arrangement of the pairs returns the quantile**: for
nonnegative weights with a positive sum and `0 < q ≤ 1`. (So the result does not depend on how
equal values are ordered.) -/
theorem isWeightedQuantile_scan {xs ys : List (ℝ × ℝ)} (hperm : ys.Perm xs)
    (hsorted : ys.Pairwise fun a b => a.1 ≤ b.1) (hw : ∀ p ∈ xs, 0 ≤ p.2)
    (hW : 0 < totalWeight xs) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) :
    IsWeightedQuantile xs q (quantileScan ys 0 (q * totalWeight ys) 0) := by
  have hWy := totalWeight_perm hperm
  have hwy : ∀ p ∈ ys, 0 ≤ p.2 := fun p hp => hw p (hperm.mem_iff.mp hp)
  obtain ⟨⟨p, hp, hpe⟩, h2, h3⟩ := quantileScan_sorted hsorted hwy
    (target := q * totalWeight ys) (by rw [hWy]; exact mul_pos hq0 hW)
    (by rw [hWy]; nlinarith) 0
  obtain ⟨v, hv⟩ : ∃ v, v = quantileScan ys 0 (q * totalWeight ys) 0 := ⟨_, rfl⟩
  rw [← hv] at hpe h2 h3 ⊢
  refine ⟨⟨p, hperm.mem_iff.mp hp, hpe⟩, ?_, ?_⟩
  · rw [weightLE_perm hperm, hWy] at h2
    exact h2
  · intro a ha hlt
    have := h3 a (hperm.mem_iff.mpr ha) hlt
    rw [weightLE_perm hperm, hWy] at this
    exact this

/-- The order of the scan: ascending values. -/
def valueLE (a b : ℝ × ℝ) : Prop := a.1 ≤ b.1

instance : DecidableRel valueLE := fun a b => inferInstanceAs (Decidable (a.1 ≤ b.1))

instance : Std.Total valueLE := ⟨fun a b => le_total a.1 b.1⟩

instance : IsTrans (ℝ × ℝ) valueLE := ⟨fun _ _ _ h1 h2 => le_trans h1 h2⟩

/-- **The weighted quantile of SPEC §1**, `quantile_q(x; ℓ)`, for the list `xs` of the pairs
`(x_c, ℓ_c)` in the order of the cut points: the values are put in ascending order (insertion
sort, which keeps the order of equal values), `Σℓ` is formed, and the scan returns the first
value at which the weights summed so far reach `q · Σℓ`. -/
def weightedQuantile (xs : List (ℝ × ℝ)) (q : ℝ) : ℝ :=
  quantileScan (xs.insertionSort valueLE) 0 (q * totalWeight (xs.insertionSort valueLE)) 0

/-- **A cut without points gives `0`.** -/
theorem weightedQuantile_nil (q : ℝ) : weightedQuantile [] q = 0 := by
  simp [weightedQuantile, quantileScan_nil]

/-- The quantile is `0` (the fallback) or one of the values, whatever the weights and `q`. -/
theorem weightedQuantile_mem_or (xs : List (ℝ × ℝ)) (q : ℝ) :
    weightedQuantile xs q = 0 ∨ ∃ p ∈ xs, weightedQuantile xs q = p.1 := by
  rcases quantileScan_mem_or (xs.insertionSort valueLE) 0
    (q * totalWeight (xs.insertionSort valueLE)) 0 with h | ⟨p, hp, h⟩
  · exact Or.inl h
  · exact Or.inr ⟨p, (List.perm_insertionSort valueLE xs).mem_iff.mp hp, h⟩

/-- **The weighted quantile, characterized.** For nonnegative weights with a positive sum and
`0 < q ≤ 1`: the quantile is one of the given values; the total weight of the values not
exceeding it is at least `q` times the total weight; and every given value below it has less
than that at or below itself. -/
theorem weightedQuantile_spec {xs : List (ℝ × ℝ)} (hw : ∀ p ∈ xs, 0 ≤ p.2)
    (hW : 0 < totalWeight xs) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) :
    IsWeightedQuantile xs q (weightedQuantile xs q) :=
  isWeightedQuantile_scan (List.perm_insertionSort valueLE xs)
    (List.pairwise_insertionSort valueLE xs) hw hW hq0 hq1

/-- It lies between the smallest and the largest value. -/
theorem weightedQuantile_bounds {xs : List (ℝ × ℝ)} (hw : ∀ p ∈ xs, 0 ≤ p.2)
    (hW : 0 < totalWeight xs) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) {lo hi : ℝ}
    (hlo : ∀ p ∈ xs, lo ≤ p.1) (hhi : ∀ p ∈ xs, p.1 ≤ hi) :
    lo ≤ weightedQuantile xs q ∧ weightedQuantile xs q ≤ hi := by
  obtain ⟨p, hp, he⟩ := (weightedQuantile_spec hw hW hq0 hq1).mem
  rw [he]
  exact ⟨hlo p hp, hhi p hp⟩

/-- The quantile is the least value at or below which at least the fraction `q` of the weight
lies. -/
theorem weightedQuantile_le {xs : List (ℝ × ℝ)} (hw : ∀ p ∈ xs, 0 ≤ p.2)
    (hW : 0 < totalWeight xs) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) {p : ℝ × ℝ} (hp : p ∈ xs)
    (h : q * totalWeight xs ≤ weightLE xs p.1) : weightedQuantile xs q ≤ p.1 := by
  by_contra hcon
  have := (weightedQuantile_spec hw hW hq0 hq1).least p hp (not_le.mp hcon)
  linarith

/-- The quantile does not decrease in `q`. -/
theorem weightedQuantile_mono {xs : List (ℝ × ℝ)} (hw : ∀ p ∈ xs, 0 ≤ p.2)
    (hW : 0 < totalWeight xs) {q q' : ℝ} (hq0 : 0 < q) (hqq : q ≤ q') (hq1 : q' ≤ 1) :
    weightedQuantile xs q ≤ weightedQuantile xs q' := by
  have hq0' : 0 < q' := lt_of_lt_of_le hq0 hqq
  have hs := weightedQuantile_spec hw hW hq0' hq1
  obtain ⟨p, hp, he⟩ := hs.mem
  have h2 := hs.reach
  rw [he] at h2 ⊢
  exact weightedQuantile_le hw hW hq0 (hqq.trans hq1) hp
    ((mul_le_mul_of_nonneg_right hqq hW.le).trans h2)

/-- With `q = 1` and positive weights the quantile is the largest value. -/
theorem weightedQuantile_one {xs : List (ℝ × ℝ)} (hne : xs ≠ []) (hw : ∀ p ∈ xs, 0 < p.2) :
    (∃ p ∈ xs, weightedQuantile xs 1 = p.1) ∧ ∀ p ∈ xs, p.1 ≤ weightedQuantile xs 1 := by
  have hs := weightedQuantile_spec (fun p hp => (hw p hp).le) (totalWeight_pos hne hw) one_pos
    le_rfl
  refine ⟨hs.mem, ?_⟩
  have h2 := hs.reach
  rw [one_mul] at h2
  have hle := weightLE_le_totalWeight (fun p hp => (hw p hp).le) (weightedQuantile xs 1)
  exact forall_le_of_weightLE_eq hw (le_antisymm hle h2)

/-- A quantile of nonnegative values is nonnegative (any weights, any `q`, also the empty
list). -/
theorem weightedQuantile_nonneg {xs : List (ℝ × ℝ)} (h : ∀ p ∈ xs, 0 ≤ p.1) (q : ℝ) :
    0 ≤ weightedQuantile xs q := by
  rcases weightedQuantile_mem_or xs q with h0 | ⟨p, hp, he⟩
  · rw [h0]
  · rw [he]
    exact h p hp

/-! ## The oracle's quantile

`Ddg.weightedQuantile q xs ls` sorts the pairs `xs.zip ls` with its own insertion sort
(`sortAscending`), sums the weights from the left (`sumLeft`) and scans (`firstReaching`). -/

/-- The oracle's `sumLeft` at `ℝ` is the sum of the list. -/
theorem sumLeft_eq_sum (l : List ℝ) : Ddg.sumLeft l = l.sum := by
  have h : ∀ (l : List ℝ) (a : ℝ), l.foldl (fun a x => a + x) a = a + l.sum := by
    intro l
    induction l with
    | nil => intro a; simp
    | cons x t ih =>
      intro a
      rw [List.foldl_cons, ih, List.sum_cons]
      ring
  unfold Ddg.sumLeft
  simp only [Oracle.zero_eq]
  exact (h l 0).trans (zero_add _)

theorem insertAscending_perm (e : ℝ × ℝ) :
    ∀ l : List (ℝ × ℝ), (Ddg.insertAscending e l).Perm (e :: l)
  | [] => by simp [Ddg.insertAscending]
  | f :: fs => by
    by_cases h : e.1 < f.1
    · simp [Ddg.insertAscending, h]
    · simp only [Ddg.insertAscending, h, ite_false]
      exact ((insertAscending_perm e fs).cons f).trans (List.Perm.swap e f fs)

theorem insertAscending_pairwise (e : ℝ × ℝ) :
    ∀ l : List (ℝ × ℝ), (l.Pairwise fun a b => a.1 ≤ b.1) →
      (Ddg.insertAscending e l).Pairwise fun a b => a.1 ≤ b.1
  | [], _ => by simp [Ddg.insertAscending]
  | f :: fs, hl => by
    rw [List.pairwise_cons] at hl
    by_cases h : e.1 < f.1
    · simp only [Ddg.insertAscending, h, ite_true]
      rw [List.pairwise_cons, List.pairwise_cons]
      refine ⟨fun b hb => ?_, hl⟩
      rcases List.mem_cons.mp hb with rfl | hb
      · exact h.le
      · exact h.le.trans (hl.1 b hb)
    · simp only [Ddg.insertAscending, h, ite_false]
      rw [List.pairwise_cons]
      refine ⟨fun b hb => ?_, insertAscending_pairwise e fs hl.2⟩
      rcases List.mem_cons.mp ((insertAscending_perm e fs).mem_iff.mp hb) with rfl | hb
      · exact not_lt.mp h
      · exact hl.1 b hb

/-- The oracle's `sortAscending` returns a permutation of its input … -/
theorem sortAscending_perm : ∀ l : List (ℝ × ℝ), (Ddg.sortAscending l).Perm l
  | [] => by simp [Ddg.sortAscending]
  | e :: es => by
    simp only [Ddg.sortAscending]
    exact (insertAscending_perm e _).trans ((sortAscending_perm es).cons e)

/-- … in ascending order of the values. -/
theorem sortAscending_pairwise :
    ∀ l : List (ℝ × ℝ), (Ddg.sortAscending l).Pairwise fun a b => a.1 ≤ b.1
  | [] => by simp [Ddg.sortAscending]
  | e :: es => by
    simp only [Ddg.sortAscending]
    exact insertAscending_pairwise e _ (sortAscending_pairwise es)

/-- **The oracle's sort puts equal values in the reverse of their original order**: a new entry
is inserted after the entries with an equal key, and the entries are inserted from the last to
the first. (Its docstring and SPEC §1 say "equal values in the order of the cut points". The
value of the quantile is not affected, `isWeightedQuantile_scan`; the order in which equal
values enter the floating-point sums is.) -/
theorem sortAscending_reverses_ties :
    Ddg.sortAscending [((1 : ℝ), (10 : ℝ)), (1, 20)] = [(1, 20), (1, 10)] := by
  simp [Ddg.sortAscending, Ddg.insertAscending]

/-- The oracle's `firstReaching` on a nonempty list is the scan `quantileScan`. -/
theorem firstReaching_eq_quantileScan (target : ℝ) :
    ∀ (l : List (ℝ × ℝ)) (acc last : ℝ), l ≠ [] →
      Ddg.firstReaching target acc l = quantileScan l acc target last
  | [], _, _, h => absurd rfl h
  | [e], acc, last, _ => by
    simp [Ddg.firstReaching, quantileScan_cons, quantileScan_nil]
  | e :: e' :: es, acc, last, _ => by
    rw [quantileScan_cons]
    by_cases h : target ≤ acc + e.2
    · simp [Ddg.firstReaching, h]
    · simp only [Ddg.firstReaching, h, ite_false]
      exact firstReaching_eq_quantileScan target (e' :: es) (acc + e.2) e.1 (by simp)

/-- The oracle's `weightedQuantile` is the scan of its sorted list. -/
theorem oracle_weightedQuantile_eq_scan (q : ℝ) (xs ls : List ℝ) :
    Ddg.weightedQuantile q xs ls =
      quantileScan (Ddg.sortAscending (xs.zip ls)) 0
        (q * totalWeight (Ddg.sortAscending (xs.zip ls))) 0 := by
  have e : Ddg.weightedQuantile q xs ls =
      Ddg.firstReaching (q * Ddg.sumLeft ((Ddg.sortAscending (xs.zip ls)).map fun e => e.2)) 0
        (Ddg.sortAscending (xs.zip ls)) := by
    simp only [Ddg.weightedQuantile, Oracle.zero_eq]
  rw [e, sumLeft_eq_sum]
  by_cases hnil : Ddg.sortAscending (xs.zip ls) = []
  · rw [hnil]
    simp [Ddg.firstReaching, quantileScan_nil]
  · rw [firstReaching_eq_quantileScan _ _ _ 0 hnil]
    rfl

/-- **The oracle's `weightedQuantile` is the quantile of SPEC §1.** For the pairs `xs.zip ls`
with nonnegative weights of positive sum and `0 < q ≤ 1`: the result is one of the given values;
the total weight of the values not exceeding it is at least `q` times the total weight; every
given value below it has less than that at or below itself.

[formula: cut-weighted-quantile] -/
theorem oracle_weightedQuantile_spec {xs ls : List ℝ} (hw : ∀ p ∈ xs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (xs.zip ls)) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) :
    IsWeightedQuantile (xs.zip ls) q (Ddg.weightedQuantile q xs ls) := by
  rw [oracle_weightedQuantile_eq_scan]
  exact isWeightedQuantile_scan (sortAscending_perm _) (sortAscending_pairwise _) hw hW hq0 hq1

/-- The oracle's quantile equals this package's `weightedQuantile` of the zipped pairs (both are
the unique value of `IsWeightedQuantile`). -/
theorem oracle_weightedQuantile_eq {xs ls : List ℝ} (hw : ∀ p ∈ xs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (xs.zip ls)) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) :
    Ddg.weightedQuantile q xs ls = weightedQuantile (xs.zip ls) q :=
  (oracle_weightedQuantile_spec hw hW hq0 hq1).unique (weightedQuantile_spec hw hW hq0 hq1)

/-- **The quantile is one of the given values and lies between the smallest and the largest.**

[formula: cut-weighted-quantile] -/
theorem oracle_weightedQuantile_bounds {xs ls : List ℝ} (hw : ∀ p ∈ xs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (xs.zip ls)) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) {lo hi : ℝ}
    (hlo : ∀ p ∈ xs.zip ls, lo ≤ p.1) (hhi : ∀ p ∈ xs.zip ls, p.1 ≤ hi) :
    (∃ p ∈ xs.zip ls, Ddg.weightedQuantile q xs ls = p.1) ∧
      lo ≤ Ddg.weightedQuantile q xs ls ∧ Ddg.weightedQuantile q xs ls ≤ hi := by
  obtain ⟨p, hp, he⟩ := (oracle_weightedQuantile_spec hw hW hq0 hq1).mem
  exact ⟨⟨p, hp, he⟩, he ▸ hlo p hp, he ▸ hhi p hp⟩

/-- **At least the fraction `q` of the weight lies at or below the quantile.**

[formula: cut-weighted-quantile] -/
theorem oracle_le_weightLE_weightedQuantile {xs ls : List ℝ} (hw : ∀ p ∈ xs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (xs.zip ls)) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) :
    q * totalWeight (xs.zip ls) ≤ weightLE (xs.zip ls) (Ddg.weightedQuantile q xs ls) :=
  (oracle_weightedQuantile_spec hw hW hq0 hq1).reach

/-- **Less than the fraction `q` of the weight lies strictly below the quantile** (more than
`1 − q` of it lies at or above).

[formula: cut-weighted-quantile] -/
theorem oracle_weightLT_weightedQuantile_lt {xs ls : List ℝ} (hw : ∀ p ∈ xs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (xs.zip ls)) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) :
    weightLT (xs.zip ls) (Ddg.weightedQuantile q xs ls) < q * totalWeight (xs.zip ls) :=
  (oracle_weightedQuantile_spec hw hW hq0 hq1).weightLT_lt hq0 hW

/-- **The quantile does not decrease in `q`.**

[formula: cut-weighted-quantile] -/
theorem oracle_weightedQuantile_mono {xs ls : List ℝ} (hw : ∀ p ∈ xs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (xs.zip ls)) {q q' : ℝ} (hq0 : 0 < q) (hqq : q ≤ q') (hq1 : q' ≤ 1) :
    Ddg.weightedQuantile q xs ls ≤ Ddg.weightedQuantile q' xs ls := by
  rw [oracle_weightedQuantile_eq hw hW hq0 (hqq.trans hq1),
    oracle_weightedQuantile_eq hw hW (lt_of_lt_of_le hq0 hqq) hq1]
  exact weightedQuantile_mono hw hW hq0 hqq hq1

/-- **With `q = 1` and positive weights the quantile is the largest value.**

[formula: cut-weighted-quantile] -/
theorem oracle_weightedQuantile_one {xs ls : List ℝ} (hne : xs.zip ls ≠ [])
    (hw : ∀ p ∈ xs.zip ls, 0 < p.2) :
    (∃ p ∈ xs.zip ls, Ddg.weightedQuantile 1 xs ls = p.1) ∧
      ∀ p ∈ xs.zip ls, p.1 ≤ Ddg.weightedQuantile 1 xs ls := by
  rw [oracle_weightedQuantile_eq (fun p hp => (hw p hp).le) (totalWeight_pos hne hw) one_pos
    le_rfl]
  exact weightedQuantile_one hne hw

/-- **An empty list gives `0`.**

[formula: cut-weighted-quantile] -/
theorem oracle_weightedQuantile_nil (q : ℝ) (ls : List ℝ) : Ddg.weightedQuantile q [] ls = 0 := by
  rw [oracle_weightedQuantile_eq_scan]
  simp [Ddg.sortAscending, quantileScan_nil]

/-- If all values are equal the quantile is that value. -/
theorem oracle_weightedQuantile_const {xs ls : List ℝ} (hw : ∀ p ∈ xs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (xs.zip ls)) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) {c : ℝ}
    (hc : ∀ p ∈ xs.zip ls, p.1 = c) : Ddg.weightedQuantile q xs ls = c := by
  obtain ⟨p, hp, he⟩ := (oracle_weightedQuantile_spec hw hW hq0 hq1).mem
  rw [he, hc p hp]

/-- A quantile of nonnegative values is nonnegative (nonnegative weights of positive sum, or
the empty list). -/
theorem oracle_weightedQuantile_nonneg {xs ls : List ℝ} (hw : ∀ p ∈ xs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (xs.zip ls)) {q : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1)
    (h : ∀ p ∈ xs.zip ls, 0 ≤ p.1) : 0 ≤ Ddg.weightedQuantile q xs ls := by
  rw [oracle_weightedQuantile_eq hw hW hq0 hq1]
  exact weightedQuantile_nonneg h q

/-! ## The roughening rate of a state (SPEC §1, D22) -/

/-- **The roughening rate of the state**, `ω = quantile_0.9(ω_c; ℓ_c)` (the oracle's
`weightedQuantile` of the local rates and the lengths): it is one of the local rates, and at
least nine tenths of the length of the cut have a local rate `ω_c ≤ ω`.

[formula: cut-roughening-rate] -/
theorem roughening_quantile_spec {ωs ls : List ℝ} (hw : ∀ p ∈ ωs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (ωs.zip ls)) :
    (∃ p ∈ ωs.zip ls, Ddg.weightedQuantile 0.9 ωs ls = p.1) ∧
      0.9 * totalWeight (ωs.zip ls) ≤ weightLE (ωs.zip ls) (Ddg.weightedQuantile 0.9 ωs ls) :=
  ⟨(oracle_weightedQuantile_spec hw hW (by norm_num) (by norm_num)).mem,
    (oracle_weightedQuantile_spec hw hW (by norm_num) (by norm_num)).reach⟩

/-- For the Dirichlet energy all local rates are nonnegative, and so is the rate of the state. -/
theorem roughening_quantile_nonneg {ωs ls : List ℝ} (hw : ∀ p ∈ ωs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (ωs.zip ls)) (h : ∀ p ∈ ωs.zip ls, 0 ≤ p.1) :
    0 ≤ Ddg.weightedQuantile 0.9 ωs ls :=
  oracle_weightedQuantile_nonneg hw hW (by norm_num) (by norm_num) h

/-- The rate of a ripple increases with the roughening rate. -/
theorem rippleRate_mono_omega {ω ω' : ℝ} (h : ω ≤ ω') (w : ℝ) {k : ℝ} (hk : 0 ≤ k) :
    rippleRate ω w k ≤ rippleRate ω' w k := by
  unfold rippleRate
  nlinarith [mul_le_mul_of_nonneg_right h hk]

/-- **Where the local rate is at most that of the state, ripples shorter than `ℓ` decay**: with
`α_L = lengthWeightForRipples ω ℓ`, `ω > 0`, a cut point with `ω_c ≤ ω` has a negative rate for
every `k > 2π/ℓ`. By `roughening_quantile_spec` this is nine tenths of the cut.

[formula: cut-roughening-rate] -/
theorem rippleRate_neg_of_le {ωc ω ℓ k : ℝ} (hc : ωc ≤ ω) (hω : 0 < ω) (hℓ : 0 < ℓ)
    (hk : 2 * π / ℓ < k) : rippleRate ωc (Ddg.lengthWeightForRipples ω ℓ) k < 0 := by
  have hk0 : 0 < k := lt_trans (by positivity) hk
  exact lt_of_le_of_lt (rippleRate_mono_omega hc _ hk0.le) (rippleRate_neg_of_shorter hω hℓ hk)

/-- **Where the local rate is at most that of the state, ripples shorter than `2πw/ω` decay**
(SPEC §1): for `ω_c ≤ ω`, `w > 0`, `k > 0` and `k > ω/w` the local rate `ω_c k − w k²` is
negative.

[formula: cut-roughening-rate] -/
theorem rippleRate_neg_of_le_of_lt {ωc ω w k : ℝ} (hc : ωc ≤ ω) (hw : 0 < w) (hk0 : 0 < k)
    (hk : ω / w < k) : rippleRate ωc w k < 0 := by
  have h1 : ω < w * k := by
    rw [div_lt_iff₀ hw] at hk
    linarith
  have h2 : rippleRate ω w k < 0 := by
    unfold rippleRate
    nlinarith
  exact lt_of_le_of_lt (rippleRate_mono_omega hc w hk0.le) h2

/-- **Where the local rate is at least that of the state, ripples longer than `ℓ` grow**: with
`α_L = lengthWeightForRipples ω ℓ`, `ω > 0`, a cut point with `ω_c ≥ ω` has a positive rate for
every `0 < k < 2π/ℓ`. By `roughening_quantile_above` this is more than a tenth of the cut (not
nine tenths, unless the rates are equal).

[formula: cut-roughening-rate] -/
theorem rippleRate_pos_of_ge {ωc ω ℓ k : ℝ} (hc : ω ≤ ωc) (hω : 0 < ω) (hℓ : 0 < ℓ)
    (hk0 : 0 < k) (hk : k < 2 * π / ℓ) : 0 < rippleRate ωc (Ddg.lengthWeightForRipples ω ℓ) k :=
  lt_of_lt_of_le (rippleRate_pos_of_longer hω hℓ hk0 hk) (rippleRate_mono_omega hc _ hk0.le)

/-- Where the local rate is below that of the state, ripples longer than `ℓ` decay too, down to
the wavelength `ℓ ω/ω_c`: the rate is negative for every `k > (ω_c/ω) · 2π/ℓ`, `k > 0`. -/
theorem rippleRate_neg_of_lt_longer {ωc ω ℓ k : ℝ} (hω : 0 < ω) (hℓ : 0 < ℓ) (hk0 : 0 < k)
    (hk : ωc / ω * (2 * π / ℓ) < k) : rippleRate ωc (Ddg.lengthWeightForRipples ω ℓ) k < 0 := by
  have hw : 0 < Ddg.lengthWeightForRipples ω ℓ := by
    rw [lengthWeightForRipples_eq]
    positivity
  refine rippleRate_neg_of_le_of_lt le_rfl hw hk0 (lt_of_eq_of_lt ?_ hk)
  rw [lengthWeightForRipples_eq]
  have := Real.pi_ne_zero
  field_simp

/-- **Less than nine tenths of the length of the cut have a local rate below that of the
state**: more than a tenth has `ω_c ≥ ω`. -/
theorem roughening_quantile_above {ωs ls : List ℝ} (hw : ∀ p ∈ ωs.zip ls, 0 ≤ p.2)
    (hW : 0 < totalWeight (ωs.zip ls)) :
    weightLT (ωs.zip ls) (Ddg.weightedQuantile 0.9 ωs ls) < 0.9 * totalWeight (ωs.zip ls) :=
  oracle_weightLT_weightedQuantile_lt hw hW (by norm_num) (by norm_num)

/-- **A state with `ω ≤ 0` gets no weight from the rule** `α_L = ωℓ/(2π)`: the result is not
positive. -/
theorem lengthWeightForRipples_nonpos {ω ℓ : ℝ} (hω : ω ≤ 0) (hℓ : 0 ≤ ℓ) :
    Ddg.lengthWeightForRipples ω ℓ ≤ 0 := by
  rw [lengthWeightForRipples_eq]
  exact div_nonpos_of_nonpos_of_nonneg (mul_nonpos_of_nonpos_of_nonneg hω hℓ) (by positivity)

/-- **Where the local rate is negative the distortion term damps every ripple** (`w ≥ 0`). -/
theorem rippleRate_neg_of_neg {ω w k : ℝ} (hω : ω < 0) (hw : 0 ≤ w) (hk : 0 < k) :
    rippleRate ω w k < 0 := by
  unfold rippleRate
  have h1 : ω * k < 0 := mul_neg_of_neg_of_pos hω hk
  have h2 : 0 ≤ w * k ^ 2 := by positivity
  linarith

/-- **The speed of a cut point is at most half its local roughening rate** when the quadratic
parts are nonnegative (the Dirichlet energy with default weights): `|Q_0 − Q_1| ≤ Q_0 + Q_1`. -/
theorem abs_cutSpeed_le_half_roughRate {Q0 Q1 : ℝ} (h0 : 0 ≤ Q0) (h1 : 0 ≤ Q1) :
    |Ddg.cutSpeed Q0 Q1| ≤ Ddg.roughRate Q0 Q1 / 2 := by
  rw [cutSpeed_eq, roughRate_eq, abs_le]
  constructor <;> linarith

/-! ## The step rule (SPEC §4) -/

/-- **A term whose denominator is not positive is skipped**: with `ω ≤ 0` and `m̂ ≤ 0` the step
is `τ_max` … -/
theorem adaptiveStep_skip_both (τmax ts mm h : ℝ) {ω ρ : ℝ} (hω : ω ≤ 0) (hρ : ρ ≤ 0) :
    Ddg.adaptiveStep τmax ts ω mm h ρ = τmax := by
  rw [adaptiveStep_eq, ite_eq_right (not_lt.mpr hρ), ite_eq_right (not_lt.mpr hω)]

/-- … with `ω ≤ 0` and a positive speed it is the capped rule
`min(τ_max, max_motion · h / m̂)` (the oracle's `stepCap`) …

[formula: cut-step-cap] -/
theorem adaptiveStep_skip_omega (τmax ts mm h : ℝ) {ω ρ : ℝ} (hω : ω ≤ 0) (hρ : 0 < ρ) :
    Ddg.adaptiveStep τmax ts ω mm h ρ = Ddg.stepCap τmax mm h ρ := by
  rw [adaptiveStep_eq, ite_eq_left hρ, ite_eq_right (not_lt.mpr hω), stepCap_eq]

/-- … with a positive `ω` and no speed it is `min(τ_max, time_scale/ω)` … -/
theorem adaptiveStep_skip_rho (τmax ts mm h : ℝ) {ω ρ : ℝ} (hω : 0 < ω) (hρ : ρ ≤ 0) :
    Ddg.adaptiveStep τmax ts ω mm h ρ = min τmax (ts / ω) := by
  rw [adaptiveStep_eq, ite_eq_right (not_lt.mpr hρ), ite_eq_left hω]

/-- … and with both positive `min(τ_max, time_scale/ω, max_motion · h / m̂)`. -/
theorem adaptiveStep_full (τmax ts mm h : ℝ) {ω ρ : ℝ} (hω : 0 < ω) (hρ : 0 < ρ) :
    Ddg.adaptiveStep τmax ts ω mm h ρ = min (min τmax (ts / ω)) (mm * h / ρ) := by
  rw [adaptiveStep_eq, ite_eq_left hρ, ite_eq_left hω]

/-- **The adaptive step is the largest step that satisfies the limits that apply**: every `τ`
with `τ ≤ τ_max`, `τ ω ≤ time_scale` if `ω > 0`, and `τ m̂ ≤ max_motion · h` if `m̂ > 0`, is at
most `adaptiveStep`.

[formula: cut-step-cap] -/
theorem le_adaptiveStep {τmax ts ω mm h ρ τ : ℝ} (h1 : τ ≤ τmax) (h2 : 0 < ω → τ * ω ≤ ts)
    (h3 : 0 < ρ → τ * ρ ≤ mm * h) : τ ≤ Ddg.adaptiveStep τmax ts ω mm h ρ := by
  rw [adaptiveStep_eq]
  split_ifs with hρ hω hω
  · exact le_min (le_min h1 ((le_div_iff₀ hω).mpr (h2 hω))) ((le_div_iff₀ hρ).mpr (h3 hρ))
  · exact le_min h1 ((le_div_iff₀ hρ).mpr (h3 hρ))
  · exact le_min h1 ((le_div_iff₀ hω).mpr (h2 hω))
  · exact h1

/-- The adaptive step is nonnegative for nonnegative parameters. -/
theorem adaptiveStep_nonneg {τmax ts mm h : ℝ} (hτ : 0 ≤ τmax) (hts : 0 ≤ ts) (hmh : 0 ≤ mm * h)
    (ω ρ : ℝ) : 0 ≤ Ddg.adaptiveStep τmax ts ω mm h ρ :=
  le_adaptiveStep hτ (fun _ => by rw [zero_mul]; exact hts) (fun _ => by rw [zero_mul]; exact hmh)

/-- **The robust step** (`h = 1`, speed `m̂` in crossed edges per unit time): a cut point whose
own speed `m_c = |ρ_c|/h_c` is at most `m̂` moves by at most `max_motion · h_c` in the step. By
`le_weightLE_weightedQuantile` that is at least nine tenths of the length of the cut when
`m̂ = quantile_0.9(m_c; ℓ_c)`; the others are limited one by one (`limited_motion_le`).

[formula: cut-step-cap] -/
theorem robustStep_motion_le {τmax ts mm ω mhat ρc hc : ℝ} (hτ : 0 ≤ τmax) (hts : 0 ≤ ts)
    (hmm : 0 ≤ mm) (hhc : 0 < hc) (hm : |ρc| / hc ≤ mhat) :
    Ddg.adaptiveStep τmax ts ω mm 1 mhat * |ρc| ≤ mm * hc := by
  have hmh : 0 ≤ mm * 1 := by rw [mul_one]; exact hmm
  have hτ0 := adaptiveStep_nonneg hτ hts hmh ω mhat
  have hρ : |ρc| ≤ mhat * hc := (div_le_iff₀ hhc).mp hm
  by_cases hpos : 0 < mhat
  · have h1 := adaptiveStep_mul_rho_le τmax ts ω mm 1 hpos
    calc Ddg.adaptiveStep τmax ts ω mm 1 mhat * |ρc|
        ≤ Ddg.adaptiveStep τmax ts ω mm 1 mhat * (mhat * hc) :=
          mul_le_mul_of_nonneg_left hρ hτ0
      _ = Ddg.adaptiveStep τmax ts ω mm 1 mhat * mhat * hc := by ring
      _ ≤ mm * 1 * hc := mul_le_mul_of_nonneg_right h1 hhc.le
      _ = mm * hc := by ring
  · have h0 : |ρc| ≤ 0 := hρ.trans (mul_nonpos_of_nonpos_of_nonneg (not_lt.mp hpos) hhc.le)
    have h0' : |ρc| = 0 := le_antisymm h0 (abs_nonneg _)
    rw [h0', mul_zero]
    exact mul_nonneg hmm hhc.le

/-- With the length step on the mesh `m̂` is the larger of two quantiles; the step then respects
both: `τ m ≤ max_motion · h` for every speed `0 ≤ m ≤ m̂`. -/
theorem adaptiveStep_mul_le_of_le {τmax ts mm h ω ρ m : ℝ} (hτ : 0 ≤ τmax) (hts : 0 ≤ ts)
    (hmh : 0 ≤ mm * h) (hρ : 0 < ρ) (hm : m ≤ ρ) :
    Ddg.adaptiveStep τmax ts ω mm h ρ * m ≤ mm * h :=
  (mul_le_mul_of_nonneg_left hm (adaptiveStep_nonneg hτ hts hmh ω ρ)).trans
    (adaptiveStep_mul_rho_le τmax ts ω mm h hρ)

/-- **The extremes** (largest speed, shortest crossed edge): every cut point moves by at most
`max_motion · h_c`.

[formula: cut-step-cap] -/
theorem extremesStep_motion_le {τmax ts mm ω hmin ρmax ρc hc : ℝ} (hτ : 0 ≤ τmax)
    (hts : 0 ≤ ts) (hmm : 0 ≤ mm) (hh : 0 ≤ hmin) (hhc : hmin ≤ hc) (hρ : |ρc| ≤ ρmax) :
    Ddg.adaptiveStep τmax ts ω mm hmin ρmax * |ρc| ≤ mm * hc := by
  have hmh : 0 ≤ mm * hmin := mul_nonneg hmm hh
  have h2 : mm * hmin ≤ mm * hc := mul_le_mul_of_nonneg_left hhc hmm
  by_cases hpos : 0 < ρmax
  · exact (adaptiveStep_mul_le_of_le hτ hts hmh hpos hρ).trans h2
  · have h0' : |ρc| = 0 := le_antisymm (hρ.trans (not_lt.mp hpos)) (abs_nonneg _)
    rw [h0', mul_zero]
    exact hmh.trans h2

/-! ## The local limit (SPEC §4, D22) -/

/-- The factor of the local limit for the allowed motion `max_motion · h_c` and the motion
`τ v_c` that the cut point would make: `allowed/moved` if that is less than one, else `1`. For
`moved > 0` it is `min(1, allowed/moved)` (`localFactor_eq_min`), the expression of SPEC §4;
for `moved = 0` it is `1`. -/
def localFactor (allowed moved : ℝ) : ℝ := if allowed < moved then allowed / moved else 1

theorem localFactor_eq_min {allowed moved : ℝ} (hm : 0 < moved) :
    localFactor allowed moved = min 1 (allowed / moved) := by
  unfold localFactor
  by_cases h : allowed < moved
  · rw [ite_eq_left h, min_eq_right ((div_le_one hm).mpr h.le)]
  · rw [ite_eq_right h, min_eq_left ((one_le_div hm).mpr (not_lt.mp h))]

/-- **A cut point within its limit is unchanged**: if `τ v_c ≤ max_motion · h_c` the factor is
`1`.

[formula: cut-step-local] -/
theorem localFactor_eq_one {allowed moved : ℝ} (h : moved ≤ allowed) :
    localFactor allowed moved = 1 := by
  unfold localFactor
  rw [ite_eq_right (not_lt.mpr h)]

/-- **The factor lies in `(0, 1]`** (for a positive allowed motion).

[formula: cut-step-local] -/
theorem localFactor_mem_Ioc {allowed : ℝ} (ha : 0 < allowed) (moved : ℝ) :
    0 < localFactor allowed moved ∧ localFactor allowed moved ≤ 1 := by
  unfold localFactor
  by_cases h : allowed < moved
  · rw [ite_eq_left h]
    have hm : 0 < moved := ha.trans h
    exact ⟨div_pos ha hm, (div_le_one hm).mpr h.le⟩
  · rw [ite_eq_right h]
    exact ⟨one_pos, le_rfl⟩

theorem localFactor_nonneg {allowed : ℝ} (ha : 0 ≤ allowed) (moved : ℝ) :
    0 ≤ localFactor allowed moved := by
  unfold localFactor
  by_cases h : allowed < moved
  · rw [ite_eq_left h]
    exact div_nonneg ha (ha.trans h.le)
  · rw [ite_eq_right h]
    exact zero_le_one

/-- The factor times the motion is at most the allowed motion. -/
theorem localFactor_mul_le {allowed : ℝ} (ha : 0 ≤ allowed) (moved : ℝ) :
    localFactor allowed moved * moved ≤ allowed := by
  unfold localFactor
  by_cases h : allowed < moved
  · rw [ite_eq_left h, div_mul_cancel₀ _ (ha.trans_lt h).ne']
  · rw [ite_eq_right h, one_mul]
    exact not_lt.mp h

/-- **No cut point moves by more than `max_motion · h_c`.** With the limited rate
`r_c · localFactor(max_motion · h_c, τ v_c)` and `v_c ≥ |r_c|`, `τ ≥ 0`:
`τ |r_c| ≤ max_motion · h_c`.

[formula: cut-step-local] -/
theorem limited_motion_le {mm h τ r v : ℝ} (hmh : 0 ≤ mm * h) (hτ : 0 ≤ τ) (hv : |r| ≤ v) :
    τ * |r * localFactor (mm * h) (τ * v)| ≤ mm * h := by
  have hf := localFactor_nonneg hmh (τ * v)
  rw [abs_mul, abs_of_nonneg hf]
  calc τ * (|r| * localFactor (mm * h) (τ * v))
      = localFactor (mm * h) (τ * v) * (τ * |r|) := by ring
    _ ≤ localFactor (mm * h) (τ * v) * (τ * v) :=
        mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hv hτ) hf
    _ ≤ mm * h := localFactor_mul_le hmh _

/-- **The limit keeps the sign of the rate, does not increase its size, and keeps a vanishing
rate** (so rest states are unchanged).

[formula: cut-step-local] -/
theorem limited_sign {allowed : ℝ} (ha : 0 < allowed) (moved r : ℝ) :
    (0 < r * localFactor allowed moved ↔ 0 < r) ∧ (r * localFactor allowed moved < 0 ↔ r < 0) ∧
      (r * localFactor allowed moved = 0 ↔ r = 0) ∧ |r * localFactor allowed moved| ≤ |r| := by
  obtain ⟨hf0, hf1⟩ := localFactor_mem_Ioc ha moved
  refine ⟨mul_pos_iff_of_pos_right hf0, ?_, ?_, ?_⟩
  · constructor
    · intro h
      by_contra hr
      have := mul_nonneg (not_lt.mp hr) hf0.le
      linarith
    · intro h
      exact mul_neg_of_neg_of_pos h hf0
  · rw [mul_eq_zero]
    exact ⟨fun h => h.resolve_right hf0.ne', Or.inl⟩
  · rw [abs_mul, abs_of_pos hf0]
    exact mul_le_of_le_one_right (abs_nonneg r) hf1

/-- The speed that the adaptive rule takes for a cut point,
`v_c = max(|r_c|, ½ (ω_c − ω)₊)`, is at least `|r_c|` (so `limited_motion_le` applies) and at
least half the excess of the local roughening rate. -/
theorem localSpeed_ge (r ωc ω : ℝ) :
    |r| ≤ max |r| (1 / 2 * max (ωc - ω) 0) ∧ 1 / 2 * (ωc - ω) ≤ max |r| (1 / 2 * max (ωc - ω) 0) := by
  refine ⟨le_max_left _ _, le_trans ?_ (le_max_right _ _)⟩
  have := le_max_left (ωc - ω) 0
  linarith

/-- **The limited rate is still a descent direction.** If the first variation of the energy is
`dE(σ) = Σ_c g_c σ_c ℓ_c` and every cut point moves with `σ_c = −f_c g_c`, `f_c ≥ 0` (the flow
`φ_t = g` limited point by point), then `dE = −Σ_c f_c g_c² ℓ_c ≤ 0`.

[formula: cut-step-local] -/
theorem firstVariation_limited {C : Type*} [Fintype C] (g f : C → ℝ) {ℓ : C → ℝ}
    (hℓ : ∀ c, 0 ≤ ℓ c) (hf : ∀ c, 0 ≤ f c) :
    firstVariation g ℓ (fun c => -(g c * f c)) = -∑ c, f c * g c ^ 2 * ℓ c ∧
      firstVariation g ℓ (fun c => -(g c * f c)) ≤ 0 := by
  have h1 : firstVariation g ℓ (fun c => -(g c * f c)) = -∑ c, f c * g c ^ 2 * ℓ c := by
    unfold firstVariation
    rw [← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl fun c _ => by ring
  refine ⟨h1, ?_⟩
  rw [h1, neg_nonpos]
  exact Finset.sum_nonneg fun c _ => mul_nonneg (mul_nonneg (hf c) (sq_nonneg _)) (hℓ c)

/-! ## The oracle's step rule and local limit -/

/-- The oracle's `speedInEdges` at `ℝ`: `m_c = |ρ_c|/h_c`. -/
theorem speedInEdges_eq (ρ hc : ℝ) : Ddg.speedInEdges ρ hc = |ρ| / hc := by
  simp only [Ddg.speedInEdges, Oracle.abs_eq]

/-- The oracle's `robustStep` is its `adaptiveStep` with the edge length `1` and the speed
`m̂` in crossed edges per unit time. -/
theorem robustStep_eq_adaptiveStep (τmax ts ω mm mhat : ℝ) :
    Ddg.robustStep τmax ts ω mm mhat = Ddg.adaptiveStep τmax ts ω mm 1 mhat := by
  rw [adaptiveStep_eq]
  simp only [Ddg.robustStep, Oracle.zero_eq, mul_one]

/-- **The robust step**: `τ ≤ τ_max`; `τ ω ≤ time_scale` if `ω > 0`; `τ m̂ ≤ max_motion` if
`m̂ > 0`; and it is the largest `τ` with these properties (a term whose denominator is not
positive imposes nothing: it is skipped).

[formula: cut-step-cap] -/
theorem robustStep_spec (τmax ts ω mm mhat : ℝ) :
    Ddg.robustStep τmax ts ω mm mhat ≤ τmax ∧
      (0 < ω → Ddg.robustStep τmax ts ω mm mhat * ω ≤ ts) ∧
      (0 < mhat → Ddg.robustStep τmax ts ω mm mhat * mhat ≤ mm) ∧
      ∀ τ : ℝ, τ ≤ τmax → (0 < ω → τ * ω ≤ ts) → (0 < mhat → τ * mhat ≤ mm) →
        τ ≤ Ddg.robustStep τmax ts ω mm mhat := by
  rw [robustStep_eq_adaptiveStep]
  refine ⟨adaptiveStep_le_max _ _ _ _ _ _, fun hω => adaptiveStep_mul_omega_le _ _ _ _ _ hω,
    fun hm => ?_, fun τ h1 h2 h3 => le_adaptiveStep h1 h2 fun hm => ?_⟩
  · have := adaptiveStep_mul_rho_le τmax ts ω mm 1 hm
    rwa [mul_one] at this
  · rw [mul_one]
    exact h3 hm

/-- **The extremes** are the same rule with `m̂ = max|ρ_c| / min h`: for a positive shortest
crossed edge, `robustStep … (ρ_max/h_min) = adaptiveStep … h_min ρ_max`
(`min(τ_max, time_scale/ω, max_motion · h_min/ρ_max)`). -/
theorem robustStep_extremes (τmax ts ω mm ρmax : ℝ) {hmin : ℝ} (hh : 0 < hmin) :
    Ddg.robustStep τmax ts ω mm (ρmax / hmin) = Ddg.adaptiveStep τmax ts ω mm hmin ρmax := by
  rw [robustStep_eq_adaptiveStep, adaptiveStep_eq, adaptiveStep_eq]
  have e1 : (0 < ρmax / hmin) ↔ (0 < ρmax) := div_pos_iff_of_pos_right hh
  have e2 : mm * 1 / (ρmax / hmin) = mm * hmin / ρmax := by
    rw [mul_one, div_div_eq_mul_div]
  by_cases hρ : 0 < ρmax
  · rw [ite_eq_left (e1.mpr hρ), ite_eq_left hρ, e2]
  · rw [ite_eq_right (fun h => hρ (e1.mp h)), ite_eq_right hρ]

/-- With the length step on the mesh `m̂` is the larger of two quantiles; the robust step then
respects both: `τ m ≤ max_motion` for every speed `m ≤ m̂` (`m̂ > 0`). -/
theorem robustStep_mul_le_of_le {τmax ts mm ω mhat m : ℝ} (hτ : 0 ≤ τmax) (hts : 0 ≤ ts)
    (hmm : 0 ≤ mm) (hpos : 0 < mhat) (hm : m ≤ mhat) :
    Ddg.robustStep τmax ts ω mm mhat * m ≤ mm := by
  rw [robustStep_eq_adaptiveStep]
  have hmh : 0 ≤ mm * 1 := by rw [mul_one]; exact hmm
  have := adaptiveStep_mul_le_of_le (ω := ω) hτ hts hmh hpos hm
  rwa [mul_one] at this

/-- **The larger of two `m̂` gives the smaller of the two steps** (`max_motion ≥ 0`): with the
length step on the mesh, the speed is taken with and without the length force. -/
theorem robustStep_max (τmax ts ω : ℝ) {mm : ℝ} (hmm : 0 ≤ mm) (m₁ m₂ : ℝ) :
    Ddg.robustStep τmax ts ω mm (max m₁ m₂) =
      min (Ddg.robustStep τmax ts ω mm m₁) (Ddg.robustStep τmax ts ω mm m₂) := by
  obtain ⟨a1, a2, a3, a4⟩ := robustStep_spec τmax ts ω mm (max m₁ m₂)
  obtain ⟨b1, b2, b3, b4⟩ := robustStep_spec τmax ts ω mm m₁
  obtain ⟨c1, c2, c3, c4⟩ := robustStep_spec τmax ts ω mm m₂
  have key : ∀ m : ℝ, m ≤ max m₁ m₂ → 0 < m →
      Ddg.robustStep τmax ts ω mm (max m₁ m₂) * m ≤ mm := by
    intro m hm hpos
    have hmax : 0 < max m₁ m₂ := lt_of_lt_of_le hpos hm
    by_cases hneg : Ddg.robustStep τmax ts ω mm (max m₁ m₂) < 0
    · exact (mul_neg_of_neg_of_pos hneg hpos).le.trans hmm
    · exact (mul_le_mul_of_nonneg_left hm (not_lt.mp hneg)).trans (a3 hmax)
  apply le_antisymm
  · exact le_min (b4 _ a1 a2 (key m₁ (le_max_left _ _))) (c4 _ a1 a2 (key m₂ (le_max_right _ _)))
  · refine a4 _ ((min_le_left _ _).trans b1)
      (fun hω => (mul_le_mul_of_nonneg_right (min_le_left _ _) hω.le).trans (b2 hω)) fun hmax => ?_
    rcases le_total m₁ m₂ with h | h
    · rw [max_eq_right h] at hmax ⊢
      exact (mul_le_mul_of_nonneg_right (min_le_right _ _) hmax.le).trans (c3 hmax)
    · rw [max_eq_left h] at hmax ⊢
      exact (mul_le_mul_of_nonneg_right (min_le_left _ _) hmax.le).trans (b3 hmax)

/-- The capped rule is the robust step without the term `time_scale/ω` (`ω ≤ 0`):
`min(τ_max, max_motion/m̂)` for `m̂ > 0`, `τ_max` otherwise. -/
theorem robustStep_capped (τmax ts mm : ℝ) {ω : ℝ} (hω : ω ≤ 0) (mhat : ℝ) :
    Ddg.robustStep τmax ts ω mm mhat = if 0 < mhat then min τmax (mm / mhat) else τmax := by
  rw [robustStep_eq_adaptiveStep, adaptiveStep_eq, ite_eq_right (not_lt.mpr hω), mul_one]

/-- **With the robust step a cut point whose own speed is at most `m̂` needs no limit**: if
`speedInEdges ρ_c h_c ≤ m̂` then `τ |ρ_c| ≤ max_motion · h_c`. With
`m̂ = weightedQuantile 0.9 (m_c) (ℓ_c)` that is at least nine tenths of the length of the cut
(`oracle_le_weightLE_weightedQuantile`).

[formula: cut-step-cap] -/
theorem robustStep_motion_le' {τmax ts mm ω mhat ρc hc : ℝ} (hτ : 0 ≤ τmax) (hts : 0 ≤ ts)
    (hmm : 0 ≤ mm) (hhc : 0 < hc) (hm : Ddg.speedInEdges ρc hc ≤ mhat) :
    Ddg.robustStep τmax ts ω mm mhat * |ρc| ≤ mm * hc := by
  rw [robustStep_eq_adaptiveStep]
  rw [speedInEdges_eq] at hm
  exact robustStep_motion_le hτ hts hmm hhc hm

/-- The oracle's `localSpeed` at `ℝ`: `max(|r|, ½ (ω_c − ω)₊)`. -/
theorem localSpeed_eq (r ωc ω : ℝ) :
    Ddg.localSpeed r ωc ω = max |r| (max 0 (ωc - ω) / 2) := by
  simp only [Ddg.localSpeed, Oracle.abs_eq, Oracle.zero_eq, Oracle.two_eq]

/-- **The speed of the local limit is at least `|r_c|`** and at least half the excess of the
local roughening rate over that of the state.

[formula: cut-step-local] -/
theorem abs_le_localSpeed (r ωc ω : ℝ) :
    |r| ≤ Ddg.localSpeed r ωc ω ∧ (ωc - ω) / 2 ≤ Ddg.localSpeed r ωc ω := by
  rw [localSpeed_eq]
  refine ⟨le_max_left _ _, le_trans ?_ (le_max_right _ _)⟩
  have := le_max_right 0 (ωc - ω)
  linarith

/-- With no excess (`ω_c ≤ ω`) the speed of the local limit is `|r_c|`: the capped rule. -/
theorem localSpeed_of_le (r : ℝ) {ωc ω : ℝ} (h : ωc ≤ ω) : Ddg.localSpeed r ωc ω = |r| := by
  rw [localSpeed_eq, max_eq_left (by linarith : ωc - ω ≤ 0), zero_div,
    max_eq_left (abs_nonneg r)]

/-- The oracle's `localLimitFactor` at `ℝ`. -/
theorem localLimitFactor_def (mm hc τ v : ℝ) :
    Ddg.localLimitFactor mm hc τ v =
      if 0 < τ * v then min 1 (mm * hc / (τ * v)) else 1 := by
  simp only [Ddg.localLimitFactor, Oracle.zero_eq, Oracle.one_eq]

/-- For a nonnegative allowed motion the oracle's factor is `localFactor`. -/
theorem localLimitFactor_eq {mm hc : ℝ} (ha : 0 ≤ mm * hc) (τ v : ℝ) :
    Ddg.localLimitFactor mm hc τ v = localFactor (mm * hc) (τ * v) := by
  rw [localLimitFactor_def]
  by_cases h : 0 < τ * v
  · rw [ite_eq_left h, localFactor_eq_min h]
  · rw [ite_eq_right h, localFactor_eq_one ((not_lt.mp h).trans ha)]

/-- **The factor of the local limit lies in `(0, 1]`** (`max_motion · h_c > 0`).

[formula: cut-step-local] -/
theorem localLimitFactor_mem_Ioc {mm hc : ℝ} (ha : 0 < mm * hc) (τ v : ℝ) :
    0 < Ddg.localLimitFactor mm hc τ v ∧ Ddg.localLimitFactor mm hc τ v ≤ 1 := by
  rw [localLimitFactor_eq ha.le]
  exact localFactor_mem_Ioc ha _

/-- **A cut point with `τ v_c ≤ max_motion · h_c` is unchanged**: the factor is `1`.

[formula: cut-step-local] -/
theorem localLimitFactor_eq_one {mm hc τ v : ℝ} (h : τ * v ≤ mm * hc) :
    Ddg.localLimitFactor mm hc τ v = 1 := by
  rw [localLimitFactor_def]
  by_cases hpos : 0 < τ * v
  · rw [ite_eq_left hpos, min_eq_left ((one_le_div hpos).mpr h)]
  · rw [ite_eq_right hpos]

/-- **A limited cut point takes the step `max_motion · h_c / v_c`** (D22): if
`τ v_c > max_motion · h_c ≥ 0` and `τ > 0`, then `τ · factor = max_motion · h_c / v_c`.

[formula: cut-step-local] -/
theorem localLimit_step_eq {mm hc τ v : ℝ} (ha : 0 ≤ mm * hc) (hτ : 0 < τ)
    (h : mm * hc < τ * v) : τ * Ddg.localLimitFactor mm hc τ v = mm * hc / v := by
  have hpos : 0 < τ * v := lt_of_le_of_lt ha h
  rw [localLimitFactor_def, ite_eq_left hpos, min_eq_right ((div_le_one hpos).mpr h.le),
    ← mul_div_assoc, mul_div_mul_left _ _ hτ.ne']

/-- **No cut point moves by more than `max_motion · h_c`.** With the limited rate
`r_c · localLimitFactor(max_motion, h_c, τ, v_c)`, `v_c ≥ |r_c|`, `τ ≥ 0`,
`max_motion · h_c ≥ 0`: `τ |r_c| ≤ max_motion · h_c`.

[formula: cut-step-local] -/
theorem localLimit_motion_le {mm hc τ r v : ℝ} (ha : 0 ≤ mm * hc) (hτ : 0 ≤ τ) (hv : |r| ≤ v) :
    τ * |r * Ddg.localLimitFactor mm hc τ v| ≤ mm * hc := by
  rw [localLimitFactor_eq ha]
  exact limited_motion_le ha hτ hv

/-- The same with the oracle's speed of the adaptive rule, `v_c = localSpeed r_c ω_c ω`.

[formula: cut-step-local] -/
theorem localLimit_motion_le_adaptive {mm hc τ r : ℝ} (ha : 0 ≤ mm * hc) (hτ : 0 ≤ τ)
    (ωc ω : ℝ) :
    τ * |r * Ddg.localLimitFactor mm hc τ (Ddg.localSpeed r ωc ω)| ≤ mm * hc :=
  localLimit_motion_le ha hτ (abs_le_localSpeed r ωc ω).1

/-- **The local limit keeps the sign of the rate, does not increase its size, and keeps a
vanishing rate**: rest states are unchanged.

[formula: cut-step-local] -/
theorem localLimit_sign {mm hc : ℝ} (ha : 0 < mm * hc) (τ v r : ℝ) :
    (0 < r * Ddg.localLimitFactor mm hc τ v ↔ 0 < r) ∧
      (r * Ddg.localLimitFactor mm hc τ v < 0 ↔ r < 0) ∧
      (r * Ddg.localLimitFactor mm hc τ v = 0 ↔ r = 0) ∧
      |r * Ddg.localLimitFactor mm hc τ v| ≤ |r| := by
  rw [localLimitFactor_eq ha.le]
  exact limited_sign ha _ _

/-- **The limited rate is still a descent direction**: with `dE(σ) = Σ_c g_c σ_c ℓ_c` and
`σ_c = −g_c · localLimitFactor(…)_c`, `dE ≤ 0`.

[formula: cut-step-local] -/
theorem firstVariation_localLimit {C : Type*} [Fintype C] (g : C → ℝ) {ℓ : C → ℝ}
    (hℓ : ∀ c, 0 ≤ ℓ c) (mm τ : ℝ) (hc v : C → ℝ) (ha : ∀ c, 0 ≤ mm * hc c) :
    firstVariation g ℓ (fun c => -(g c * Ddg.localLimitFactor mm (hc c) τ (v c))) ≤ 0 :=
  (firstVariation_limited g (fun c => Ddg.localLimitFactor mm (hc c) τ (v c)) hℓ fun c => by
    rw [localLimitFactor_eq (ha c)]
    exact localFactor_nonneg (ha c) _).2

end Cuts

end DdgProofs
