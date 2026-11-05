import DdgProofs.Developability.Vec

/-!
# Developability: the combinatorial energies `E_P` and `π_max` (paper Eq. 1, 2, 5)

`crates/ddg-developability/SPEC.md` §3.2, §3.3; App. B.2 of the paper. The oracle's
`interiorBipartitions`, `arcPart`, `arcComplement`, `boundarySplits`, `partMean`, `partWidth`,
`partWidthPairs`, `partWidthMean`, `partWidthSensitivity`, `partMaxWidth`, `pairSensitivity`
at `ℝ`.

**Bipartitions** (`dev-star-bipartitions`). The faces of an interior star of valence `k` form a
cycle `σ_0 … σ_{k−1}`. A set of faces is *edge-connected* in the cycle iff it is a cyclic interval
`{s, s+1, …, s+l−1} mod k` (or all faces); this characterization is taken as the definition here
(it is a fact about the cycle graph, not proved). Proved:

* `mem_interiorBipartitions`, `length_interiorBipartitions` — the enumeration is the pairs
  `0 ≤ a < b ≤ k − 1`, `k(k−1)/2` of them; `length_boundarySplits` — `k − 1` splits of a path;
* `arcPart_arcComplement` — for each pair the two parts `F_1 = {a+1, …, b}`,
  `F_2 = {b+1, …, k−1, 0, …, a}` are non-empty, disjoint, and cover the `k` faces;
* `arc_bipartition_injective` — different pairs give different unordered partitions `{F_1, F_2}`;
* `cycInterval_mem_bipartitions` — every proper non-empty cyclic interval is `F_1` or `F_2` of
  some pair (and then its complement is the other part).

Together: the enumeration lists every partition of the cycle into two non-empty cyclic intervals
exactly once (D9).

**Widths.**

* `partWidth_eq_pairs_eq_mean` — Eq. 1, three forms equal for unit normals (`n ≥ 1`), the second
  over **unordered** pairs (D1); `sum_ordered_pairs` — over ordered pairs the sum is twice that;
* `hasDerivAt_partWidth` — along differentiable normals `N_σ(t)`, `dπ = Σ_σ (2/n)(N_σ − N̄)·N_σ'`
  (any vectors, not only unit ones);
* `partMaxWidth_spec` — `π_max` of a part is the maximum of `|N_σ1 − N_σ2|²` over its pairs
  (Eq. 5: an upper bound of every pair term, attained by some pair);
  `hasDerivAt_pairTerm` — the derivative of the active pair term is
  `2(N_σ1 − N_σ2)·(N_σ1' − N_σ2')`, the sensitivities `pairSensitivity`.
-/

noncomputable section

namespace DDG

namespace Dev

open scoped DDG.Oracle
open Finset

/-! ## Bipartitions of a cycle and of a path -/

section bipartitions

open Ddg.Dev

theorem mem_interiorBipartitions (k a b : ℕ) :
    (a, b) ∈ interiorBipartitions k ↔ a < b ∧ b < k := by
  simp only [interiorBipartitions, List.mem_flatMap, List.mem_range, List.mem_map,
    List.mem_filter, decide_eq_true_eq, Prod.mk.injEq]
  constructor
  · rintro ⟨a', ha', b', ⟨hb', hab'⟩, rfl, rfl⟩; exact ⟨hab', hb'⟩
  · rintro ⟨hab, hb⟩; exact ⟨a, by omega, b, ⟨hb, hab⟩, rfl, rfl⟩

private theorem length_filter_range_lt (a k : ℕ) :
    ((List.range k).filter (fun b => decide (a < b))).length = k - (a + 1) := by
  induction k with
  | zero => simp
  | succ k ih =>
    rw [List.range_succ, List.filter_append, List.length_append, ih]
    by_cases h : a < k <;> simp [h] <;> omega

private theorem list_sum_map_range (k : ℕ) (f : ℕ → ℕ) :
    ((List.range k).map f).sum = ∑ i ∈ range k, f i := by
  induction k with
  | zero => simp
  | succ k ih => rw [List.range_succ, List.map_append, List.sum_append, ih, sum_range_succ]; simp

/-- The number of bipartitions of an interior star of valence `k`: `k(k−1)/2`.

[formula: dev-star-bipartitions] -/
theorem length_interiorBipartitions (k : ℕ) :
    (interiorBipartitions k).length = k * (k - 1) / 2 := by
  have hL : (interiorBipartitions k).length = ∑ a ∈ range k, (k - (a + 1)) := by
    rw [interiorBipartitions, List.length_flatMap, ← list_sum_map_range]
    congr 1
    apply List.map_congr_left
    intro a _
    simp [length_filter_range_lt]
  have hr : ∑ a ∈ range k, (k - (a + 1)) = ∑ a ∈ range k, a := by
    rw [← sum_range_reflect (fun a => a) k]
    apply sum_congr rfl
    intro a ha
    rw [mem_range] at ha
    omega
  rw [hL, hr]
  exact (Nat.div_eq_of_eq_mul_left (by norm_num) (sum_range_id_mul_two k).symm).symm

/-- A boundary star with `k` faces has the `k − 1` splits `a = 0 … k − 2`. -/
theorem length_boundarySplits (k : ℕ) : (boundarySplits k).length = k - 1 := by
  simp [boundarySplits]

theorem mem_boundarySplits (k a : ℕ) : a ∈ boundarySplits k ↔ a + 2 ≤ k := by
  simp [boundarySplits]; omega

theorem mem_arcPart (a b x : ℕ) : x ∈ arcPart a b ↔ a < x ∧ x ≤ b := by
  simp only [arcPart, List.mem_map, List.mem_range]
  constructor
  · rintro ⟨i, hi, rfl⟩; omega
  · rintro ⟨h1, h2⟩; exact ⟨x - (a + 1), by omega, by omega⟩

private theorem mod_of_lt_two {n k : ℕ} (h1 : k ≤ n) (h2 : n < 2 * k) : n % k = n - k := by
  rw [Nat.mod_eq_sub_mod h1, Nat.mod_eq_of_lt (by omega)]

theorem mem_arcComplement {k a b : ℕ} (hab : a < b) (hb : b < k) (x : ℕ) :
    x ∈ arcComplement k a b ↔ x < k ∧ ¬(a < x ∧ x ≤ b) := by
  simp only [arcComplement, List.mem_map, List.mem_range]
  constructor
  · rintro ⟨i, hi, rfl⟩
    by_cases h : b + 1 + i < k
    · rw [Nat.mod_eq_of_lt h]; omega
    · rw [mod_of_lt_two (by omega) (by omega)]; omega
  · rintro ⟨hx, hn⟩
    by_cases h : b < x
    · exact ⟨x - b - 1, by omega, by rw [Nat.mod_eq_of_lt (by omega)]; omega⟩
    · refine ⟨x + k - b - 1, by omega, ?_⟩
      rw [mod_of_lt_two (by omega) (by omega)]; omega

/-- The part `F_1` as a set of face indices. -/
def part1 (a b : ℕ) : Finset ℕ := (arcPart a b).toFinset

/-- The part `F_2` as a set of face indices. -/
def part2 (k a b : ℕ) : Finset ℕ := (arcComplement k a b).toFinset

theorem mem_part1 (a b x : ℕ) : x ∈ part1 a b ↔ a < x ∧ x ≤ b := by
  rw [part1, List.mem_toFinset, mem_arcPart]

theorem mem_part2 {k a b : ℕ} (hab : a < b) (hb : b < k) (x : ℕ) :
    x ∈ part2 k a b ↔ x < k ∧ ¬(a < x ∧ x ≤ b) := by
  rw [part2, List.mem_toFinset, mem_arcComplement hab hb]

/-- **The two parts of a bipartition** `(a, b)`, `a < b < k`: `F_1 = {a+1, …, b}` and
`F_2 = {b+1, …, k−1, 0, …, a}` are non-empty, disjoint, and together the `k` faces; `F_1` never
contains face `0`, `F_2` always does. -/
theorem arcPart_arcComplement {k a b : ℕ} (hab : a < b) (hb : b < k) :
    (part1 a b).Nonempty ∧ (part2 k a b).Nonempty ∧ Disjoint (part1 a b) (part2 k a b) ∧
      part1 a b ∪ part2 k a b = range k ∧ 0 ∉ part1 a b ∧ 0 ∈ part2 k a b := by
  refine ⟨⟨b, (mem_part1 _ _ _).2 ⟨hab, le_rfl⟩⟩, ⟨0, (mem_part2 hab hb _).2 (by omega)⟩, ?_, ?_,
    by rw [mem_part1]; omega, by rw [mem_part2 hab hb]; omega⟩
  · rw [Finset.disjoint_left]
    intro x h1 h2
    rw [mem_part1] at h1; rw [mem_part2 hab hb] at h2
    exact h2.2 h1
  · ext x
    rw [mem_union, mem_part1, mem_part2 hab hb, mem_range]
    omega

/-- **Different pairs give different unordered partitions**: if `{F_1, F_2}` of `(a, b)` equals
`{F_1, F_2}` of `(a', b')` (in either order), then `(a, b) = (a', b')`. -/
theorem arc_bipartition_injective {k a b a' b' : ℕ} (hab : a < b) (_hb : b < k)
    (hab' : a' < b') (hb' : b' < k)
    (h : (part1 a b = part1 a' b' ∧ part2 k a b = part2 k a' b') ∨
      (part1 a b = part2 k a' b' ∧ part2 k a b = part1 a' b')) :
    a = a' ∧ b = b' := by
  rcases h with ⟨h1, -⟩ | ⟨h1, -⟩
  · have e : ∀ x, (a < x ∧ x ≤ b) ↔ (a' < x ∧ x ≤ b') := by
      intro x; rw [← mem_part1, ← mem_part1, h1]
    have e1 := e (a + 1); have e2 := e (a' + 1); have e3 := e b; have e4 := e b'
    omega
  · have : (0 : ℕ) ∈ part1 a b := by rw [h1, mem_part2 hab' hb']; omega
    rw [mem_part1] at this; omega

/-- The cyclic interval `{s, s+1, …, s+l−1} mod k`. -/
def cycInterval (k s l : ℕ) : Finset ℕ := ((List.range l).map fun i => (s + i) % k).toFinset

theorem mem_cycInterval {k s l : ℕ} (hs : s < k) (hl : l ≤ k) (x : ℕ) :
    x ∈ cycInterval k s l ↔ x < k ∧ (s ≤ x ∧ x < s + l ∨ x + k < s + l) := by
  simp only [cycInterval, List.mem_toFinset, List.mem_map, List.mem_range]
  constructor
  · rintro ⟨i, hi, rfl⟩
    by_cases h : s + i < k
    · rw [Nat.mod_eq_of_lt h]; omega
    · rw [mod_of_lt_two (by omega) (by omega)]; omega
  · rintro ⟨hx, h | h⟩
    · exact ⟨x - s, by omega, by rw [Nat.mod_eq_of_lt (by omega)]; omega⟩
    · refine ⟨x + k - s, by omega, ?_⟩
      rw [mod_of_lt_two (by omega) (by omega)]; omega

/-- **Every proper non-empty cyclic interval is a part of an enumerated bipartition**: for
`s < k` and `1 ≤ l ≤ k − 1` there is a pair `a < b < k` whose `F_1` or `F_2` is the interval
`{s, …, s+l−1} mod k` (the other part is then its complement, `arcPart_arcComplement`). -/
theorem cycInterval_mem_bipartitions {k s l : ℕ} (hs : s < k) (hl1 : 1 ≤ l) (hl : l < k) :
    ∃ a b, a < b ∧ b < k ∧ (cycInterval k s l = part1 a b ∨ cycInterval k s l = part2 k a b) := by
  rcases Nat.eq_zero_or_pos s with rfl | hs0
  · refine ⟨l - 1, k - 1, by omega, by omega, Or.inr ?_⟩
    ext x
    rw [mem_cycInterval hs hl.le, mem_part2 (by omega) (by omega)]
    omega
  · by_cases hw : s + l ≤ k
    · refine ⟨s - 1, s + l - 1, by omega, by omega, Or.inl ?_⟩
      ext x
      rw [mem_cycInterval hs hl.le, mem_part1]
      omega
    · refine ⟨s + l - k - 1, s - 1, by omega, by omega, Or.inr ?_⟩
      ext x
      rw [mem_cycInterval hs hl.le, mem_part2 (by omega) (by omega)]
      omega

end bipartitions

/-! ## Sums -/

theorem sumTo_eq_range (n : ℕ) (f : ℕ → ℝ) : Ddg.sumTo n f = ∑ i ∈ range n, f i := by
  induction n with
  | zero => simp [Ddg.sumTo]
  | succ n ih => rw [sum_range_succ, ← ih]; rfl

/-- `Σ (a_i − m)² = Σ a_i² − 2 m Σ a_i + n m²`. -/
theorem sum_sub_sq (n : ℕ) (a : ℕ → ℝ) (m : ℝ) :
    ∑ i ∈ range n, (a i - m) ^ 2 =
      ∑ i ∈ range n, a i ^ 2 - 2 * m * ∑ i ∈ range n, a i + n * m ^ 2 := by
  induction n with
  | zero => simp
  | succ n ih => rw [sum_range_succ, ih, sum_range_succ, sum_range_succ]; push_cast; ring

/-- Unordered pairs: `Σ_{j<i<n} (a_i − a_j)² = n Σ a_i² − (Σ a_i)²`. -/
theorem sum_pairs_sq (n : ℕ) (a : ℕ → ℝ) :
    ∑ i ∈ range n, ∑ j ∈ range i, (a i - a j) ^ 2 =
      n * ∑ i ∈ range n, a i ^ 2 - (∑ i ∈ range n, a i) ^ 2 := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [sum_range_succ, ih]
    have h : ∑ j ∈ range n, (a n - a j) ^ 2 =
        n * a n ^ 2 - 2 * a n * ∑ j ∈ range n, a j + ∑ j ∈ range n, a j ^ 2 := by
      have := sum_sub_sq n a (a n)
      rw [show ∑ j ∈ range n, (a n - a j) ^ 2 = ∑ j ∈ range n, (a j - a n) ^ 2 from
        sum_congr rfl (fun j _ => by ring), this]
      ring
    rw [h, sum_range_succ, sum_range_succ]
    push_cast; ring

/-- Ordered pairs: `Σ_{i,j<n} (a_i − a_j)² = 2 (n Σ a_i² − (Σ a_i)²)`. -/
theorem sum_ordered_pairs_sq (n : ℕ) (a : ℕ → ℝ) :
    ∑ i ∈ range n, ∑ j ∈ range n, (a i - a j) ^ 2 =
      2 * (n * ∑ i ∈ range n, a i ^ 2 - (∑ i ∈ range n, a i) ^ 2) := by
  have : ∀ i, ∑ j ∈ range n, (a i - a j) ^ 2 =
      n * a i ^ 2 - 2 * a i * ∑ j ∈ range n, a j + ∑ j ∈ range n, a j ^ 2 := by
    intro i
    rw [show ∑ j ∈ range n, (a i - a j) ^ 2 = ∑ j ∈ range n, (a j - a i) ^ 2 from
      sum_congr rfl (fun j _ => by ring), sum_sub_sq]
    ring
  simp only [this, sum_add_distrib, sum_sub_distrib, ← mul_sum, ← sum_mul, sum_const, card_range,
    nsmul_eq_mul]
  ring

/-! ## The width of a part, Eq. 1 -/

section width

open Ddg.Dev

theorem norm2_sub_eq (p q : V3) :
    norm2 (sub p q) = (p.x - q.x) ^ 2 + (p.y - q.y) ^ 2 + (p.z - q.z) ^ 2 := by
  simp only [norm2_def, sub_x, sub_y, sub_z]; ring

theorem norm2_eq (p : V3) : norm2 p = p.x ^ 2 + p.y ^ 2 + p.z ^ 2 := by
  simp only [norm2_def]; ring

theorem partMean_x (n : ℕ) (N : ℕ → V3) :
    (partMean n N).x = (∑ i ∈ range n, (N i).x) / n := by
  simp [partMean, vsumTo, sumTo_eq_range, div_eq_inv_mul]

theorem partMean_y (n : ℕ) (N : ℕ → V3) :
    (partMean n N).y = (∑ i ∈ range n, (N i).y) / n := by
  simp [partMean, vsumTo, sumTo_eq_range, div_eq_inv_mul]

theorem partMean_z (n : ℕ) (N : ℕ → V3) :
    (partMean n N).z = (∑ i ∈ range n, (N i).z) / n := by
  simp [partMean, vsumTo, sumTo_eq_range, div_eq_inv_mul]

/-- First form of Eq. 1 for any vectors: `(1/n) Σ |N_σ − N̄|² = (1/n) Σ |N_σ|² − |N̄|²`. -/
theorem partWidth_eq (n : ℕ) (hn : 0 < n) (N : ℕ → V3) :
    partWidth n N = (∑ i ∈ range n, norm2 (N i)) / n - norm2 (partMean n N) := by
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  simp only [partWidth, sumTo_eq_range, scalar_ofNat, norm2_eq, sub_x, sub_y, sub_z,
    sum_add_distrib, partMean_x, partMean_y, partMean_z]
  rw [sum_sub_sq, sum_sub_sq, sum_sub_sq]
  field_simp
  ring

/-- Second form of Eq. 1 for any vectors: `(1/n²) Σ_{σ1<σ2} |N_σ1 − N_σ2|² = (1/n) Σ |N_σ|² − |N̄|²`. -/
theorem partWidthPairs_eq (n : ℕ) (hn : 0 < n) (N : ℕ → V3) :
    partWidthPairs n N = (∑ i ∈ range n, norm2 (N i)) / n - norm2 (partMean n N) := by
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  simp only [partWidthPairs, sumTo_eq_range, scalar_ofNat, norm2_eq, sub_x, sub_y, sub_z,
    sum_add_distrib, partMean_x, partMean_y, partMean_z]
  rw [sum_pairs_sq, sum_pairs_sq, sum_pairs_sq]
  field_simp
  ring

/-- **Eq. 1, the three forms** (SPEC §3.2, D1): for a part of `n ≥ 1` unit normals,
`(1/n) Σ_σ |N_σ − N̄|² = (1/n²) Σ_{σ1<σ2} |N_σ1 − N_σ2|² = 1 − |N̄|²`, the middle sum over
**unordered** pairs (over ordered pairs it is twice as large, `sum_ordered_pairs`).

[formula: dev-combinatorial-width] -/
theorem partWidth_eq_pairs_eq_mean (n : ℕ) (hn : 0 < n) (N : ℕ → V3)
    (hunit : ∀ i < n, norm2 (N i) = 1) :
    partWidth n N = partWidthPairs n N ∧ partWidthPairs n N = partWidthMean n N := by
  have hs : ∑ i ∈ range n, norm2 (N i) = n := by
    rw [sum_congr rfl (fun i hi => hunit i (mem_range.mp hi))]; simp
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  refine ⟨by rw [partWidth_eq n hn, partWidthPairs_eq n hn], ?_⟩
  rw [partWidthPairs_eq n hn, hs, div_self hn', partWidthMean]
  simp

/-- D1: the second line of Eq. 1 summed over **ordered** pairs is twice the unordered sum. -/
theorem sum_ordered_pairs (n : ℕ) (N : ℕ → V3) :
    ∑ i ∈ range n, ∑ j ∈ range n, norm2 (sub (N i) (N j)) =
      2 * ∑ i ∈ range n, ∑ j ∈ range i, norm2 (sub (N i) (N j)) := by
  simp only [norm2_sub_eq, sum_add_distrib, sum_ordered_pairs_sq, sum_pairs_sq]
  ring

/-- **Gradient of the width of a part** (SPEC §3.2): along differentiable vectors `N_σ(t)`
(`n ≥ 1`), `d/dt π = Σ_σ y_σ · N_σ'` with `y_σ = (2/n)(N_σ − N̄)` (`partWidthSensitivity`).
Composed with Eq. 8 this is the gradient of App. B.2.

[formula: dev-combinatorial-gradient] -/
theorem hasDerivAt_partWidth (n : ℕ) (hn : 0 < n) (N : ℕ → ℝ → V3) (N' : ℕ → V3) (t : ℝ)
    (hN : ∀ i < n, HasDerivAt3 (N i) (N' i) t) :
    HasDerivAt (fun s => partWidth n (fun i => N i s))
      (Ddg.sumTo n (fun i => dot (partWidthSensitivity n (fun i => N i t) i) (N' i))) t := by
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  simp only [partWidth_eq n hn]
  -- the sum of squared norms
  have h1 : HasDerivAt (fun s => ∑ i ∈ range n, norm2 (N i s))
      (∑ i ∈ range n, 2 * dot (N i t) (N' i)) t :=
    HasDerivAt.fun_sum (fun i hi => (hN i (mem_range.mp hi)).norm2)
  -- the mean
  have hM : HasDerivAt3 (fun s => partMean n (fun i => N i s))
      ⟨(∑ i ∈ range n, (N' i).x) / n, (∑ i ∈ range n, (N' i).y) / n,
        (∑ i ∈ range n, (N' i).z) / n⟩ t := by
    refine ⟨?_, ?_, ?_⟩ <;> simp only [partMean_x, partMean_y, partMean_z]
    · exact (HasDerivAt.fun_sum (fun i hi => (hN i (mem_range.mp hi)).x)).div_const _
    · exact (HasDerivAt.fun_sum (fun i hi => (hN i (mem_range.mp hi)).y)).div_const _
    · exact (HasDerivAt.fun_sum (fun i hi => (hN i (mem_range.mp hi)).z)).div_const _
  convert (h1.div_const (n : ℝ)).sub hM.norm2 using 1
  rw [sumTo_eq_range]
  set M := partMean n (fun i => N i t)
  have key : ∀ i, dot (partWidthSensitivity n (fun i => N i t) i) (N' i) =
      (2 / n) * dot (N i t) (N' i) - (2 / n) * dot M (N' i) := by
    intro i
    simp only [partWidthSensitivity, dot_def, smul_x, smul_y, smul_z, sub_x, sub_y, sub_z,
      scalar_ofNat, DDG.Oracle.two_eq, M]
    ring
  have key2 : ∑ i ∈ range n, dot M (N' i) =
      M.x * ∑ i ∈ range n, (N' i).x + M.y * ∑ i ∈ range n, (N' i).y +
        M.z * ∑ i ∈ range n, (N' i).z := by
    simp only [dot_def, sum_add_distrib, ← mul_sum]
  rw [sum_congr rfl (fun i _ => key i), sum_sub_distrib, ← mul_sum, ← mul_sum, key2, ← mul_sum]
  simp only [dot_def]
  field_simp

end width

/-! ## The max width, Eq. 5 -/

section maxwidth

open Ddg.Dev

theorem maxTo1_succ (n : ℕ) (f : ℕ → ℝ) : maxTo1 (n + 1) f = max (maxTo1 n f) (f (n + 1)) := rfl

theorem le_maxTo1 (f : ℕ → ℝ) : ∀ {n m : ℕ}, m ≤ n → f m ≤ maxTo1 n f
  | 0, m, h => by rw [Nat.le_zero.mp h]; exact le_rfl
  | n + 1, m, h => by
    rw [maxTo1_succ]
    rcases Nat.lt_or_eq_of_le h with h | rfl
    · exact le_max_of_le_left (le_maxTo1 f (by omega))
    · exact le_max_right _ _

theorem exists_maxTo1 (f : ℕ → ℝ) : ∀ n : ℕ, ∃ m ≤ n, maxTo1 n f = f m
  | 0 => ⟨0, le_rfl, rfl⟩
  | n + 1 => by
    rw [maxTo1_succ]
    obtain ⟨m, hm, he⟩ := exists_maxTo1 f n
    rcases le_total (maxTo1 n f) (f (n + 1)) with h | h
    · exact ⟨n + 1, le_rfl, max_eq_right h⟩
    · exact ⟨m, by omega, by rw [max_eq_left h, he]⟩

theorem norm2_sub_comm (p q : V3) : norm2 (sub p q) = norm2 (sub q p) := by
  rw [norm2_sub_eq, norm2_sub_eq]; ring

/-- **Eq. 5, the width of a part for `π_max`** (SPEC §3.3): for the normals `N_0 … N_n` of a part,
`partMaxWidth` is an upper bound of `|N_σ1 − N_σ2|²` for every pair and equals it for some pair:
it is `max_{σ1,σ2 ∈ F_p} |N_σ1 − N_σ2|²` (a one-face part has width `0`).

[formula: dev-combinatorial-max] -/
theorem partMaxWidth_spec (n : ℕ) (N : ℕ → V3) :
    (∀ i ≤ n, ∀ j ≤ n, norm2 (sub (N i) (N j)) ≤ partMaxWidth n N) ∧
      ∃ i ≤ n, ∃ j ≤ i, partMaxWidth n N = norm2 (sub (N i) (N j)) := by
  constructor
  · intro i hi j hj
    rcases le_total j i with h | h
    · exact (le_maxTo1 (fun j => norm2 (sub (N i) (N j))) h).trans
        (le_maxTo1 (fun i => maxTo1 i (fun j => norm2 (sub (N i) (N j)))) hi)
    · rw [norm2_sub_comm]
      exact (le_maxTo1 (fun i' => norm2 (sub (N j) (N i'))) h).trans
        (le_maxTo1 (fun i => maxTo1 i (fun j => norm2 (sub (N i) (N j)))) hj)
  · obtain ⟨i, hi, he⟩ := exists_maxTo1 (fun i => maxTo1 i (fun j => norm2 (sub (N i) (N j)))) n
    obtain ⟨j, hj, he'⟩ := exists_maxTo1 (fun j => norm2 (sub (N i) (N j))) i
    exact ⟨i, hi, j, hj, by rw [partMaxWidth, he, he']⟩

/-- The derivative of the active pair term of `π_max`: along differentiable normals,
`d/dt |N_1 − N_2|² = y·N_1' − y·N_2'` with `y = 2(N_1 − N_2)` (`pairSensitivity`). -/
theorem hasDerivAt_pairTerm {N1 N2 : ℝ → V3} {N1' N2' : V3} {t : ℝ}
    (h1 : HasDerivAt3 N1 N1' t) (h2 : HasDerivAt3 N2 N2' t) :
    HasDerivAt (fun s => norm2 (sub (N1 s) (N2 s)))
      (dot (pairSensitivity (N1 t) (N2 t)) N1' - dot (pairSensitivity (N1 t) (N2 t)) N2') t := by
  convert (h1.sub h2).norm2 using 1
  simp only [pairSensitivity, dot_def, smul_x, smul_y, smul_z, sub_x, sub_y, sub_z,
    DDG.Oracle.two_eq]
  ring

end maxwidth

end Dev

end DDG
