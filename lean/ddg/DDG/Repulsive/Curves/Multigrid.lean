import DDG.Repulsive.Curves.Saddle

/-!
# Repulsive Curves: coarsening and prolongation (paper B.3.2, B.3.3; SPEC §8.4, D17)

* `rc-coarsen`: the black/white rule of SPEC §8.4 in the form its correctness uses. Vertices are
  explored in an order `ord` (the breadth-first order of the code; any injective order); a vertex
  may be white only if it has degree 2, carries no `Point`/`Tangent(Vertex)` constraint and no
  already explored neighbour is white (`CoarsenRule`). Then no two adjacent vertices are white
  (`no_adjacent_whites`), so every neighbour of a white vertex is black and the coarse edges
  `(a, b)` that replace a white vertex join black vertices. (That loops with fewer than 6 vertices
  are frozen, so that no level gets two edges between the same vertices, is a check of the code.)
* `rc-prolongation`: `J` keeps the value at a black vertex and averages the two black neighbours at
  a white one; `J 𝟙 = 𝟙`; and `Π = P_ℓ M P_{ℓ+1}` (`M = J ⊗ I₃` or any matrix) has its range in
  `range P_ℓ = ker C_ℓ` (`projector_props`).
-/

noncomputable section

namespace DDG

namespace RC

open Matrix Finset

/-! ## Coarsening -/

section coarsen

variable {V E : Type*} (ep : E → Fin 2 → V)

/-- `v` and `w` are joined by an edge. -/
def Adj (v w : V) : Prop := ∃ I, (ep I 0 = v ∧ ep I 1 = w) ∨ (ep I 0 = w ∧ ep I 1 = v)


theorem adj_symm {v w : V} (h : Adj ep v w) : Adj ep w v := by
  obtain ⟨I, h | h⟩ := h
  · exact ⟨I, Or.inr h⟩
  · exact ⟨I, Or.inl h⟩

/-- The coarsening rule of SPEC §8.4 (D17): the exploration order is injective; a white vertex has
degree 2, carries no `Point`/`Tangent(Vertex)` constraint, and none of its already explored
neighbours is white.

[formula: rc-coarsen] -/
structure CoarsenRule (ord : V → ℕ) (white : V → Prop) (deg : V → ℕ) (pinned : V → Prop) :
    Prop where
  ord_inj : Function.Injective ord
  white_deg : ∀ v, white v → deg v = 2
  white_free : ∀ v, white v → ¬ pinned v
  white_nbr : ∀ v w, white v → Adj ep v w → ord w < ord v → ¬ white w

/-- **No two adjacent vertices are white** (for a network without self-loops), so every neighbour of
a white vertex is black, and the coarse edge `(a, b)` replacing a white vertex `v` with neighbours
`a, b` joins two black vertices.

[formula: rc-coarsen] -/
theorem no_adjacent_whites {ord : V → ℕ} {white : V → Prop} {deg : V → ℕ} {pinned : V → Prop}
    (h : CoarsenRule ep ord white deg pinned) (hloop : ∀ I, ep I 0 ≠ ep I 1) {v w : V}
    (hvw : Adj ep v w) : ¬ (white v ∧ white w) := by
  rintro ⟨hv, hw⟩
  have hne : v ≠ w := by
    rintro rfl
    obtain ⟨I, h1 | h1⟩ := hvw
    · exact hloop I (h1.1.trans h1.2.symm)
    · exact hloop I (h1.1.trans h1.2.symm)
  have hord : ord v ≠ ord w := fun h' => hne (h.ord_inj h')
  rcases lt_or_gt_of_ne hord with hlt | hlt
  · exact h.white_nbr w v hw (adj_symm ep hvw) hlt hv
  · exact h.white_nbr v w hv hvw hlt hw

/-- **The coarse edges join black vertices**: a fine edge with two black ends, or the edge `(a, b)`
replacing a white vertex with neighbours `a` and `b`; and a white vertex has degree 2 and no
constraint.

[formula: rc-coarsen] -/
theorem coarse_edge_black {ord : V → ℕ} {white : V → Prop} {deg : V → ℕ} {pinned : V → Prop}
    (h : CoarsenRule ep ord white deg pinned) (hloop : ∀ I, ep I 0 ≠ ep I 1) {v a b : V}
    (hv : white v) (ha : Adj ep v a) (hb : Adj ep v b) :
    ¬ white a ∧ ¬ white b ∧ deg v = 2 ∧ ¬ pinned v :=
  ⟨fun hw => no_adjacent_whites ep h hloop ha ⟨hv, hw⟩,
    fun hw => no_adjacent_whites ep h hloop hb ⟨hv, hw⟩, h.white_deg v hv, h.white_free v hv⟩

end coarsen

/-! ## Prolongation -/

section prolongation

variable {V C : Type*} [Fintype C] [DecidableEq C] (black : V → Prop) [DecidablePred black]
  (idx : V → C) (na nb : V → V)

/-- The prolongation `J ∈ R^{|V_f|×|V_c|}`: `1` at a black vertex (coarse index `idx v`), `½` and `½`
from the two black neighbours `na v`, `nb v` of a white one.

[formula: rc-prolongation] -/
def prolong : Matrix V C ℝ := fun v c =>
  if black v then (if idx v = c then 1 else 0)
  else (1 / 2) * (if idx (na v) = c then 1 else 0) + (1 / 2) * (if idx (nb v) = c then 1 else 0)

/-- **The prolongation keeps values at black vertices, averages at white ones, and `J 𝟙 = 𝟙`.**

[formula: rc-prolongation] -/
theorem prolong_props (u : C → ℝ) (v : V) :
    (prolong black idx na nb *ᵥ u) v =
        (if black v then u (idx v) else (u (idx (na v)) + u (idx (nb v))) / 2) ∧
      prolong black idx na nb *ᵥ (fun _ => (1 : ℝ)) = fun _ => 1 := by
  have key : ∀ (u : C → ℝ) (v : V), (prolong black idx na nb *ᵥ u) v =
      (if black v then u (idx v) else (u (idx (na v)) + u (idx (nb v))) / 2) := by
    intro u v
    simp only [prolong, mulVec, dotProduct]
    split_ifs with hb
    · simp [ite_mul, Finset.sum_ite_eq]
    · simp only [add_mul, Finset.sum_add_distrib, mul_assoc, ← Finset.mul_sum, ite_mul, one_mul,
        zero_mul, Finset.sum_ite_eq, Finset.mem_univ, ite_true]
      ring
  refine ⟨key u v, funext fun v => ?_⟩
  rw [key]
  split_ifs <;> norm_num

/-- **The projected prolongation has its range in the constraint space** (B.3.3): for
`P_ℓ = I − C_ℓ†C_ℓ` (`C_ℓ` of full row rank) and any `M` (the reference's `P_ℓ (J ⊗ I₃) P_{ℓ+1}`),
`C_ℓ P_ℓ M P' = 0`, i.e. `C_ℓ (Π y) = 0` for every coarse `y`.

[formula: rc-prolongation] -/
theorem projected_prolongation_range {nf nc m : Type*} [Fintype nf] [DecidableEq nf] [Fintype m]
    [DecidableEq m] [Fintype nc] {Cf : Matrix m nf ℝ} (hC : IsUnit (Cf * Cfᵀ).det) {k : Type*}
    [Fintype k] (M : Matrix nf k ℝ)
    (P' : Matrix k nc ℝ) :
    Cf * (projP Cf * M * P') = 0 ∧ ∀ y, Cf *ᵥ ((projP Cf * M * P') *ᵥ y) = 0 := by
  have h0 : Cf * (projP Cf * M * P') = 0 := by
    rw [← Matrix.mul_assoc, ← Matrix.mul_assoc, (projector_props hC).2.1, Matrix.zero_mul,
      Matrix.zero_mul]
  exact ⟨h0, fun y => by rw [mulVec_mulVec, h0, zero_mulVec]⟩

end prolongation

end RC

end DDG
