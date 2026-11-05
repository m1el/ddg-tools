import DdgProofs.RepulsiveCurves.Constraints

/-!
# Repulsive Curves: length schedules and subdivision (SPEC §5.4, D13)

* `rc-length-schedule`: the total target `L*_k = L⁰ + sign(f − 1) · min(k · rate · L⁰, |f − 1| L⁰)`
  starts at `L⁰`, stays between `L⁰` and `f L⁰`, moves by at most `rate · L⁰` per step and reaches
  `f L⁰` once `k · rate ≥ |f − 1|`; the edge targets `ℓ⁰_I L*_k / L⁰` keep every ratio and sum to
  `L*_k` (the reference's equal absolute increments do not keep the ratios, D13).
* `rc-subdivision-constraints`: splitting every edge at its midpoint (`V ⊕ E` vertices, the midpoint
  of `I` is `inr I`; children `inl I = (i₁, m)` and `inr I = (m, i₂)`, the SPEC's numbering
  `|V| + I`, `I`, `|E| + I`) keeps the old positions (`Point`), halves the edge lengths
  (`EdgeLengths` with halved targets), keeps the tangents of the children (`Tangent(Edge)`), the
  total length (`TotalLength`), the length-weighted barycenter (`Barycenter`), and the unit
  directions from a vertex to its new neighbours (`Tangent(Vertex)`). New `Surface` rows on
  midpoints are not satisfied in general (the next re-target projection fixes them).
-/

noncomputable section

namespace DDG

namespace RC

open Matrix Real Finset

/-! ## Length schedules -/

/-- The total target length of step `k`: `L*_k = L⁰ + sign(f − 1) · min(k · rate · L⁰, |f − 1| L⁰)`.

[formula: rc-length-schedule] -/
def scheduleTarget (L0 f rate : ℝ) (k : ℕ) : ℝ :=
  L0 + Real.sign (f - 1) * min (k * rate * L0) (|f - 1| * L0)

/-- **The schedule of the total target** (`L⁰ > 0`, `rate > 0`): `L*_0 = L⁰`; `L*_k = f L⁰` once
`k · rate ≥ |f − 1|`; `L*_k` lies between `L⁰` and `f L⁰`; `|L*_{k+1} − L*_k| ≤ rate · L⁰`.

[formula: rc-length-schedule] -/
theorem schedule_props {L0 f rate : ℝ} (hL : 0 < L0) (hr : 0 < rate) :
    scheduleTarget L0 f rate 0 = L0 ∧
      (∀ k : ℕ, |f - 1| ≤ k * rate → scheduleTarget L0 f rate k = f * L0) ∧
      (∀ k, min L0 (f * L0) ≤ scheduleTarget L0 f rate k ∧
        scheduleTarget L0 f rate k ≤ max L0 (f * L0)) ∧
      (∀ k, |scheduleTarget L0 f rate (k + 1) - scheduleTarget L0 f rate k| ≤ rate * L0) := by
  have hmin0 : ∀ k : ℕ, 0 ≤ min ((k : ℝ) * rate * L0) (|f - 1| * L0) := fun k =>
    le_min (by positivity) (by positivity)
  have hminle : ∀ k : ℕ, min ((k : ℝ) * rate * L0) (|f - 1| * L0) ≤ |f - 1| * L0 := fun k =>
    min_le_right _ _
  refine ⟨?_, fun k hk => ?_, fun k => ?_, fun k => ?_⟩
  · simp [scheduleTarget, abs_nonneg, hL.le, mul_nonneg]
  · have : min ((k : ℝ) * rate * L0) (|f - 1| * L0) = |f - 1| * L0 :=
      min_eq_right (by nlinarith)
    rw [scheduleTarget, this]
    rcases lt_trichotomy (f - 1) 0 with h | h | h
    · rw [Real.sign_of_neg h, abs_of_neg h]; ring
    · rw [h, Real.sign_zero]; have : f = 1 := by linarith
      rw [this]; ring
    · rw [Real.sign_of_pos h, abs_of_pos h]; ring
  · have h0 := hmin0 k
    have h1 := hminle k
    unfold scheduleTarget
    rcases lt_trichotomy (f - 1) 0 with h | h | h
    · rw [Real.sign_of_neg h]
      rw [abs_of_neg h] at h1 h0 ⊢
      constructor
      · exact (min_le_right _ _).trans (by linarith)
      · exact (by linarith : _ ≤ L0).trans (le_max_left _ _)
    · rw [h, Real.sign_zero, zero_mul, add_zero]
      exact ⟨min_le_left _ _, le_max_left _ _⟩
    · rw [Real.sign_of_pos h]
      rw [abs_of_pos h] at h1 h0 ⊢
      constructor
      · exact (min_le_left _ _).trans (by linarith)
      · exact (by linarith : _ ≤ f * L0).trans (le_max_right _ _)
  · unfold scheduleTarget
    have hstep : |min (((k + 1 : ℕ) : ℝ) * rate * L0) (|f - 1| * L0)
        - min ((k : ℝ) * rate * L0) (|f - 1| * L0)| ≤ rate * L0 := by
      refine (abs_min_sub_min_le_max _ _ _ _).trans ?_
      rw [sub_self, abs_zero]
      refine max_le ?_ (by positivity)
      push_cast
      rw [show ((k : ℝ) + 1) * rate * L0 - k * rate * L0 = rate * L0 by ring,
        abs_of_pos (by positivity)]
    have hs : |Real.sign (f - 1)| ≤ 1 := by
      rcases lt_trichotomy (f - 1) 0 with h | h | h
      · rw [Real.sign_of_neg h]; norm_num
      · rw [h, Real.sign_zero]; norm_num
      · rw [Real.sign_of_pos h]; norm_num
    rw [add_sub_add_left_eq_sub, ← mul_sub, abs_mul]
    calc |Real.sign (f - 1)| * |min (((k + 1 : ℕ) : ℝ) * rate * L0) (|f - 1| * L0)
          - min ((k : ℝ) * rate * L0) (|f - 1| * L0)| ≤ 1 * (rate * L0) :=
          mul_le_mul hs hstep (abs_nonneg _) zero_le_one
      _ = rate * L0 := one_mul _

/-- **Edge targets with a common factor** (D13): `ℓ*_I = ℓ⁰_I L*/L⁰` keeps every ratio
`ℓ*_I/ℓ*_J = ℓ⁰_I/ℓ⁰_J` and, when `Σ ℓ⁰ = L⁰`, sums to `L*`.

[formula: rc-length-schedule] -/
theorem edge_targets {E : Type*} [Fintype E] (ℓ0 : E → ℝ) {L0 Lstar : ℝ} (hL : L0 ≠ 0) :
    (∀ I J, ℓ0 I * Lstar / L0 * ℓ0 J = ℓ0 J * Lstar / L0 * ℓ0 I) ∧
      (∑ I, ℓ0 I = L0 → ∑ I, ℓ0 I * Lstar / L0 = Lstar) := by
  refine ⟨fun I J => by ring, fun h => ?_⟩
  rw [← Finset.sum_div, ← Finset.sum_mul, h]
  field_simp

/-! ## Subdivision -/

section subdivision

variable {V E : Type*} (γ : V → V3) (ep : E → Fin 2 → V)

/-- Positions after splitting every edge at its midpoint: old vertices keep theirs, the new vertex
`inr I` sits at `x_I`. -/
def subPos : V ⊕ E → V3 := Sum.elim γ (mid γ ep)

/-- Edges after the split: `inl I = (i₁, m_I)` (keeps index `I`), `inr I = (m_I, i₂)`. -/
def subEp : E ⊕ E → Fin 2 → V ⊕ E
  | Sum.inl I => ![Sum.inl (ep I 0), Sum.inr I]
  | Sum.inr I => ![Sum.inr I, Sum.inl (ep I 1)]

theorem subEdgeVec (I : E) :
    edgeVec (subPos γ ep) (subEp ep) (Sum.inl I) = (1 / 2 : ℝ) • edgeVec γ ep I ∧
      edgeVec (subPos γ ep) (subEp ep) (Sum.inr I) = (1 / 2 : ℝ) • edgeVec γ ep I := by
  constructor <;> (ext i; simp [edgeVec, subPos, subEp, mid]; ring)

/-- **Children of a split edge** (SPEC §5.4): each child has half the length (so `EdgeLengths`
with halved targets still hold) and the same unit tangent (`Tangent(Edge(I))` holds on both).

[formula: rc-subdivision-constraints] -/
theorem sub_children (I : E) :
    len (subPos γ ep) (subEp ep) (Sum.inl I) = len γ ep I / 2 ∧
      len (subPos γ ep) (subEp ep) (Sum.inr I) = len γ ep I / 2 ∧
      tangent (subPos γ ep) (subEp ep) (Sum.inl I) = tangent γ ep I ∧
      tangent (subPos γ ep) (subEp ep) (Sum.inr I) = tangent γ ep I := by
  obtain ⟨h1, h2⟩ := subEdgeVec γ ep I
  have hl : ∀ K, edgeVec (subPos γ ep) (subEp ep) K = (1 / 2 : ℝ) • edgeVec γ ep I →
      len (subPos γ ep) (subEp ep) K = len γ ep I / 2 := fun K hK => by
    rw [len, hK, vnorm_smul, abs_of_pos (by norm_num), len]; ring
  have ht : ∀ K, edgeVec (subPos γ ep) (subEp ep) K = (1 / 2 : ℝ) • edgeVec γ ep I →
      tangent (subPos γ ep) (subEp ep) K = tangent γ ep I := fun K hK => by
    rw [tangent, hl K hK, hK, smul_smul, tangent]
    congr 1
    by_cases h : len γ ep I = 0
    · simp [h]
    · field_simp
  exact ⟨hl _ h1, hl _ h2, ht _ h1, ht _ h2⟩

/-- **Old vertices keep their positions** (`Point` constraints still hold).

[formula: rc-subdivision-constraints] -/
theorem sub_old_positions (v : V) : subPos γ ep (Sum.inl v) = γ v := rfl

variable [Fintype E]

/-- **The total length and the length-weighted barycenter are unchanged by the split**
(`TotalLength`, `Barycenter` still hold).

[formula: rc-subdivision-constraints] -/
theorem sub_total_bary :
    totalLen (subPos γ ep) (subEp ep) = totalLen γ ep ∧
      bary (subPos γ ep) (subEp ep) = bary γ ep := by
  have hT : totalLen (subPos γ ep) (subEp ep) = totalLen γ ep := by
    simp only [totalLen, Fintype.sum_sum_type, (sub_children γ ep _).1, (sub_children γ ep _).2.1,
      ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun I _ => by ring
  refine ⟨hT, ?_⟩
  rw [bary, bary, hT]
  congr 1
  simp only [Fintype.sum_sum_type, (sub_children γ ep _).1, (sub_children γ ep _).2.1,
    ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun I _ => ?_
  ext i
  simp [mid, subPos, subEp]
  ring

omit [Fintype E] in
/-- **The direction from a vertex to its new neighbour is unchanged** (`Tangent(Vertex(v))` still
holds): for an edge `I` with `i₂ = v`, the unit vector from the midpoint to `v` equals the unit
vector from `i₁` to `v`; likewise for `i₁ = v`.

[formula: rc-subdivision-constraints] -/
theorem sub_vertex_direction (I : E) :
    (vnorm (subPos γ ep (Sum.inl (ep I 1)) - subPos γ ep (Sum.inr I)))⁻¹ •
        (subPos γ ep (Sum.inl (ep I 1)) - subPos γ ep (Sum.inr I)) =
      (vnorm (γ (ep I 1) - γ (ep I 0)))⁻¹ • (γ (ep I 1) - γ (ep I 0)) ∧
    (vnorm (subPos γ ep (Sum.inr I) - subPos γ ep (Sum.inl (ep I 0))))⁻¹ •
        (subPos γ ep (Sum.inr I) - subPos γ ep (Sum.inl (ep I 0))) =
      (vnorm (γ (ep I 1) - γ (ep I 0)))⁻¹ • (γ (ep I 1) - γ (ep I 0)) := by
  have e1 : subPos γ ep (Sum.inl (ep I 1)) - subPos γ ep (Sum.inr I) =
      (1 / 2 : ℝ) • (γ (ep I 1) - γ (ep I 0)) := by
    ext i; simp [subPos, mid]; ring
  have e2 : subPos γ ep (Sum.inr I) - subPos γ ep (Sum.inl (ep I 0)) =
      (1 / 2 : ℝ) • (γ (ep I 1) - γ (ep I 0)) := by
    ext i; simp [subPos, mid]; ring
  have key : ∀ d : V3, (vnorm ((1 / 2 : ℝ) • d))⁻¹ • ((1 / 2 : ℝ) • d) = (vnorm d)⁻¹ • d := by
    intro d
    rw [vnorm_smul, abs_of_pos (by norm_num), smul_smul]
    by_cases h : vnorm d = 0
    · simp [h]
    · congr 1; field_simp
  rw [e1, e2, key]
  exact ⟨rfl, rfl⟩

end subdivision

end RC

end DDG
