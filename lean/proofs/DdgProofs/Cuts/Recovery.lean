import DdgProofs.Cuts.Step

/-!
# Cut flow: recovery of the normal derivative from the integrated fluxes

`crates/ddg-variational-cuts/SPEC.md` §3.5, D12. The integrated fluxes `v_p` at the cut vertices
`C` of a region are turned into normal derivatives `d_p ≈ (∂u/∂n)_p` by one linear solve,

`ℓ_p d_p + Σ_{σ = pq at p} κ_σ (d_p − d_q) = v_p`   for every `p ∈ C`,

with the dual lengths `ℓ_p > 0` and the couplings `κ_σ = c² h_f² / |pq| ≥ 0` of the cut segments.
Here the couplings are a symmetric function `k p q ≥ 0` on pairs of cut vertices (`κ_σ` if the
segment `σ` joins `p` and `q`, `0` otherwise); a *piece* of the cut is a set `S` of cut vertices
with `k p q = 0` whenever `p ∈ S`, `q ∉ S`.

**`cut-flux-recovery`** (`Ddg.fluxDiagEntry`, `Ddg.fluxOffEntry`, `Ddg.fluxRhs`,
`Ddg.fluxCoupling`, `Ddg.faceMeanSquareEdge`):

* `fluxMatrix_isSymm`, `fluxMatrix_posDef`: the matrix `diag(ℓ) + Σ_σ κ_σ [[1,−1],[−1,1]]`
  assembled from the oracle's entries is symmetric positive definite;
  `fluxRecovery_existsUnique`: the recovered `d` exists and is unique;
* `fluxRecovery_sum_piece`, `fluxRecovery_sum`: `Σ_{p ∈ S} ℓ_p d_p = Σ_{p ∈ S} v_p` on every
  piece of the cut (the total flux is kept);
* `fluxRecovery_ge_piece`, `fluxRecovery_le_piece`: on every piece, `d` stays between the
  smallest and the largest `v_p/ℓ_p` of that piece (discrete maximum principle; the off-diagonal
  entries are `≤ 0`); `fluxRecovery_const_piece`: a constant quotient `v_p/ℓ_p` is reproduced;
* `fluxRecovery_nodal`: with all couplings `0` (`c = 0`, `fluxCoupling_zero`) it is the nodal
  quotient `d_p = v_p/ℓ_p` of the paper;
* `fluxRecovery_smul`, `fluxRecovery_add`: the map `v ↦ d` is linear (so `q = −2c·v` gives
  `e = −2c·d`, used by `gradScaledD_const`).
-/

noncomputable section

namespace DdgProofs

namespace Cuts

open Matrix
open scoped DdgProofs.Oracle

section Recovery

variable {C : Type*} [Fintype C] [DecidableEq C]

/-- A matrix with symmetric nonpositive off-diagonal entries and zero row sums (a graph
Laplacian with nonnegative edge weights) is positive semidefinite. -/
theorem lapOfOffDiag_posSemidef {o : C → C → ℝ} (ho : ∀ i j, o i j = o j i)
    (hneg : ∀ i j, o i j ≤ 0) : (lapOfOffDiag o).PosSemidef := by
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg ?_ fun u => ?_
  · exact Matrix.isHermitian_iff_isSymm.mpr (lapOfOffDiag_isSymm ho)
  · rw [star_trivial, dotProduct_lapOfOffDiag_mulVec ho]
    refine mul_nonneg (by norm_num) (Finset.sum_nonneg fun i _ => Finset.sum_nonneg fun j _ => ?_)
    exact mul_nonneg (neg_nonneg.mpr (hneg i j)) (sq_nonneg _)

/-- The oracle's entries at `ℝ`. -/
theorem fluxDiagEntry_eq (ℓ κ : ℝ) : Ddg.fluxDiagEntry ℓ κ = ℓ + κ := rfl

theorem fluxOffEntry_eq (κ : ℝ) : Ddg.fluxOffEntry κ = -κ := rfl

theorem fluxRhs_eq (v : ℝ) : Ddg.fluxRhs v = v := rfl

/-- The oracle's coupling `κ_σ = c² h_f² / |pq|` is nonnegative. -/
theorem fluxCoupling_nonneg (c : ℝ) {h2 s : ℝ} (hh : 0 ≤ h2) (hs : 0 ≤ s) :
    0 ≤ Ddg.fluxCoupling c h2 s := by
  unfold Ddg.fluxCoupling
  exact div_nonneg (mul_nonneg (mul_self_nonneg c) hh) hs

/-- With `c = 0` there is no coupling (`FluxRecovery::Nodal`). -/
theorem fluxCoupling_zero (h2 s : ℝ) : Ddg.fluxCoupling 0 h2 s = 0 := by
  unfold Ddg.fluxCoupling
  simp

/-- The oracle's `h_f² = (ℓ_ab² + ℓ_bc² + ℓ_ca²)/3` is nonnegative. -/
theorem faceMeanSquareEdge_nonneg (lab lbc lca : ℝ) : 0 ≤ Ddg.faceMeanSquareEdge lab lbc lca := by
  simp only [Ddg.faceMeanSquareEdge, Oracle.three_eq]
  have := mul_self_nonneg lab
  have := mul_self_nonneg lbc
  have := mul_self_nonneg lca
  positivity

/-- **The matrix of the flux recovery**, assembled from the oracle's entries: diagonal
`fluxDiagEntry ℓ_p (Σ_{q ≠ p} k p q) = ℓ_p + Σ_{σ ∋ p} κ_σ`, off-diagonal
`fluxOffEntry (k p q) = −κ_σ` (SPEC §3.5). -/
def fluxMatrix (ℓ : C → ℝ) (k : C → C → ℝ) : Matrix C C ℝ :=
  Matrix.of fun p q =>
    if p = q then Ddg.fluxDiagEntry (ℓ p) (∑ r, if r = p then 0 else k p r)
    else Ddg.fluxOffEntry (k p q)

/-- The flux matrix is `diag(ℓ)` plus the graph Laplacian of the couplings. -/
theorem fluxMatrix_eq (ℓ : C → ℝ) (k : C → C → ℝ) :
    fluxMatrix ℓ k = Matrix.diagonal ℓ + lapOfOffDiag (fun p q => -k p q) := by
  ext p q
  simp only [fluxMatrix, lapOfOffDiag, Matrix.of_apply, Matrix.add_apply, fluxDiagEntry_eq,
    fluxOffEntry_eq]
  by_cases h : p = q
  · subst h
    have e : ∑ r, (if r = p then 0 else -k p r) = -∑ r, (if r = p then 0 else k p r) := by
      rw [← Finset.sum_neg_distrib]
      exact Finset.sum_congr rfl fun r _ => by split_ifs <;> simp
    simp [e]
  · simp [h]

/-- **Row `p` of the system**: `ℓ_p d_p + Σ_q k p q (d_p − d_q)`, the left-hand side of SPEC
§3.5. -/
theorem fluxMatrix_mulVec_apply (ℓ : C → ℝ) (k : C → C → ℝ) (d : C → ℝ) (p : C) :
    (fluxMatrix ℓ k *ᵥ d) p = ℓ p * d p + ∑ q, k p q * (d p - d q) := by
  rw [fluxMatrix_eq, Matrix.add_mulVec, Pi.add_apply, Matrix.mulVec_diagonal,
    lapOfOffDiag_mulVec_apply]
  simp only [neg_neg]

/-- The system in the form of SPEC §3.5. -/
theorem fluxRecovery_iff_rows (ℓ : C → ℝ) (k : C → C → ℝ) (v d : C → ℝ) :
    fluxMatrix ℓ k *ᵥ d = (fun p => Ddg.fluxRhs (v p)) ↔
      ∀ p, ℓ p * d p + ∑ q, k p q * (d p - d q) = v p := by
  simp only [funext_iff, fluxMatrix_mulVec_apply, fluxRhs_eq]

/-- **The flux matrix is symmetric.**

[formula: cut-flux-recovery] -/
theorem fluxMatrix_isSymm (ℓ : C → ℝ) {k : C → C → ℝ} (hk : ∀ p q, k p q = k q p) :
    (fluxMatrix ℓ k).IsSymm := by
  rw [fluxMatrix_eq]
  exact (Matrix.isSymm_diagonal ℓ).add (lapOfOffDiag_isSymm fun p q => by rw [hk p q])

/-- **The flux matrix is positive definite** for positive dual lengths and symmetric nonnegative
couplings.

[formula: cut-flux-recovery] -/
theorem fluxMatrix_posDef {ℓ : C → ℝ} (hℓ : ∀ p, 0 < ℓ p) {k : C → C → ℝ}
    (hk : ∀ p q, k p q = k q p) (hk0 : ∀ p q, 0 ≤ k p q) : (fluxMatrix ℓ k).PosDef := by
  rw [fluxMatrix_eq]
  have hpsd : (lapOfOffDiag (fun p q => -k p q)).PosSemidef :=
    lapOfOffDiag_posSemidef (fun p q => by rw [hk p q]) (fun p q => neg_nonpos.mpr (hk0 p q))
  have := lengthStep_posDef hℓ hpsd zero_le_one
  rwa [one_smul] at this

/-- **The recovered normal derivative exists and is unique.**

[formula: cut-flux-recovery] -/
theorem fluxRecovery_existsUnique {ℓ : C → ℝ} (hℓ : ∀ p, 0 < ℓ p) {k : C → C → ℝ}
    (hk : ∀ p q, k p q = k q p) (hk0 : ∀ p q, 0 ≤ k p q) (v : C → ℝ) :
    ∃! d : C → ℝ, fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p) :=
  existsUnique_mulVec_eq ((Matrix.isUnit_iff_isUnit_det _).mp (fluxMatrix_posDef hℓ hk hk0).isUnit) _

/-- **The total flux of every piece of the cut is kept**: if no coupling leaves the set `S` of
cut vertices, `Σ_{p ∈ S} ℓ_p d_p = Σ_{p ∈ S} v_p`.

[formula: cut-flux-recovery] -/
theorem fluxRecovery_sum_piece {ℓ v d : C → ℝ} {k : C → C → ℝ} (hk : ∀ p q, k p q = k q p)
    (h : fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p)) (S : Finset C)
    (hS : ∀ p ∈ S, ∀ q, q ∉ S → k p q = 0) : ∑ p ∈ S, ℓ p * d p = ∑ p ∈ S, v p := by
  have hrow := (fluxRecovery_iff_rows ℓ k v d).mp h
  have h1 : ∀ p ∈ S, ∑ q, k p q * (d p - d q) = ∑ q ∈ S, k p q * (d p - d q) := by
    intro p hp
    symm
    refine Finset.sum_subset (Finset.subset_univ S) fun q _ hq => ?_
    rw [hS p hp q hq, zero_mul]
  have h2 : ∑ p ∈ S, ∑ q ∈ S, k p q * (d p - d q) = 0 := by
    have e : ∑ p ∈ S, ∑ q ∈ S, k p q * (d p - d q) =
        ∑ p ∈ S, ∑ q ∈ S, -(k p q * (d p - d q)) := by
      rw [Finset.sum_comm]
      exact Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun q _ => by
        rw [hk q p]
        ring
    simp only [Finset.sum_neg_distrib] at e
    linarith
  calc ∑ p ∈ S, ℓ p * d p = ∑ p ∈ S, (ℓ p * d p + ∑ q ∈ S, k p q * (d p - d q)) := by
        rw [Finset.sum_add_distrib, h2, add_zero]
    _ = ∑ p ∈ S, v p := Finset.sum_congr rfl fun p hp => by
        rw [← h1 p hp]
        exact hrow p

/-- **The total flux is kept**: `Σ_p ℓ_p d_p = Σ_p v_p`.

[formula: cut-flux-recovery] -/
theorem fluxRecovery_sum {ℓ v d : C → ℝ} {k : C → C → ℝ} (hk : ∀ p q, k p q = k q p)
    (h : fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p)) : ∑ p, ℓ p * d p = ∑ p, v p :=
  fluxRecovery_sum_piece hk h Finset.univ fun _ _ q hq => absurd (Finset.mem_univ q) hq

/-- **Discrete maximum principle, lower bound.** On a piece `S` of the cut, a lower bound of the
nodal quotients `v_p/ℓ_p` is a lower bound of the recovered `d`.

[formula: cut-flux-recovery] -/
theorem fluxRecovery_ge_piece {ℓ v d : C → ℝ} {k : C → C → ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (hk0 : ∀ p q, 0 ≤ k p q) (h : fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p))
    (S : Finset C) (hS : ∀ p ∈ S, ∀ q, q ∉ S → k p q = 0) {m : ℝ}
    (hm : ∀ p ∈ S, m ≤ v p / ℓ p) : ∀ p ∈ S, m ≤ d p := by
  intro p hp
  have hrow := (fluxRecovery_iff_rows ℓ k v d).mp h
  obtain ⟨i₀, hi₀, hmin⟩ := Finset.exists_min_image S d ⟨p, hp⟩
  have hK : ∑ q, k i₀ q * (d i₀ - d q) ≤ 0 := by
    refine Finset.sum_nonpos fun q _ => ?_
    by_cases hq : q ∈ S
    · exact mul_nonpos_of_nonneg_of_nonpos (hk0 i₀ q) (sub_nonpos.mpr (hmin q hq))
    · rw [hS i₀ hi₀ q hq, zero_mul]
  have h1 : v i₀ ≤ d i₀ * ℓ i₀ := by
    have := hrow i₀
    linarith
  have h2 : v i₀ / ℓ i₀ ≤ d i₀ := (div_le_iff₀ (hℓ i₀)).mpr h1
  exact (hm i₀ hi₀).trans (h2.trans (hmin p hp))

/-- **The recovery is linear**: scaling the fluxes scales the result. -/
theorem fluxRecovery_smul {ℓ v d : C → ℝ} {k : C → C → ℝ}
    (h : fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p)) (a : ℝ) :
    fluxMatrix ℓ k *ᵥ (a • d) = fun p => Ddg.fluxRhs (a * v p) := by
  rw [Matrix.mulVec_smul, h]
  funext p
  simp [fluxRhs_eq]

/-- **The recovery is linear**: the fluxes add. -/
theorem fluxRecovery_add {ℓ v v' d d' : C → ℝ} {k : C → C → ℝ}
    (h : fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p))
    (h' : fluxMatrix ℓ k *ᵥ d' = fun p => Ddg.fluxRhs (v' p)) :
    fluxMatrix ℓ k *ᵥ (d + d') = fun p => Ddg.fluxRhs (v p + v' p) := by
  rw [Matrix.mulVec_add, h, h']
  funext p
  simp [fluxRhs_eq]

/-- **Discrete maximum principle, upper bound.**

[formula: cut-flux-recovery] -/
theorem fluxRecovery_le_piece {ℓ v d : C → ℝ} {k : C → C → ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (hk0 : ∀ p q, 0 ≤ k p q) (h : fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p))
    (S : Finset C) (hS : ∀ p ∈ S, ∀ q, q ∉ S → k p q = 0) {m : ℝ}
    (hm : ∀ p ∈ S, v p / ℓ p ≤ m) : ∀ p ∈ S, d p ≤ m := by
  intro p hp
  have h' := fluxRecovery_smul h (-1)
  have := fluxRecovery_ge_piece hℓ hk0 h' S hS (m := -m)
    (fun q hq => by
      have := hm q hq
      rw [neg_one_mul, neg_div]
      linarith) p hp
  simp only [Pi.smul_apply, smul_eq_mul, neg_one_mul] at this
  linarith

/-- **The recovered derivative lies between the smallest and the largest nodal quotient**
(whole system).

[formula: cut-flux-recovery] -/
theorem fluxRecovery_bounds {ℓ v d : C → ℝ} {k : C → C → ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (hk0 : ∀ p q, 0 ≤ k p q) (h : fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p))
    {lo hi : ℝ} (hlo : ∀ p, lo ≤ v p / ℓ p) (hhi : ∀ p, v p / ℓ p ≤ hi) (p : C) :
    lo ≤ d p ∧ d p ≤ hi :=
  ⟨fluxRecovery_ge_piece hℓ hk0 h Finset.univ (fun _ _ q hq => absurd (Finset.mem_univ q) hq)
      (fun q _ => hlo q) p (Finset.mem_univ p),
    fluxRecovery_le_piece hℓ hk0 h Finset.univ (fun _ _ q hq => absurd (Finset.mem_univ q) hq)
      (fun q _ => hhi q) p (Finset.mem_univ p)⟩

/-- **A constant normal derivative is reproduced**: if `v_p = a ℓ_p` on a piece of the cut then
`d_p = a` there.

[formula: cut-flux-recovery] -/
theorem fluxRecovery_const_piece {ℓ v d : C → ℝ} {k : C → C → ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (hk0 : ∀ p q, 0 ≤ k p q) (h : fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p))
    (S : Finset C) (hS : ∀ p ∈ S, ∀ q, q ∉ S → k p q = 0) {a : ℝ}
    (hv : ∀ p ∈ S, v p = a * ℓ p) : ∀ p ∈ S, d p = a := by
  intro p hp
  have hquot : ∀ q ∈ S, v q / ℓ q = a := fun q hq => by
    rw [hv q hq]
    have := (hℓ q).ne'
    field_simp
  exact le_antisymm
    (fluxRecovery_le_piece hℓ hk0 h S hS (fun q hq => (hquot q hq).le) p hp)
    (fluxRecovery_ge_piece hℓ hk0 h S hS (fun q hq => (hquot q hq).ge) p hp)

/-- **Without coupling the recovery is the nodal quotient of the paper**: `d_p = v_p/ℓ_p`, the
oracle's `normalDerivative`.

[formula: cut-flux-recovery] -/
theorem fluxRecovery_nodal {ℓ v d : C → ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (h : fluxMatrix ℓ (fun _ _ => 0) *ᵥ d = fun p => Ddg.fluxRhs (v p)) (p : C) :
    d p = Ddg.normalDerivative (v p) (ℓ p) := by
  have hrow := (fluxRecovery_iff_rows ℓ (fun _ _ => 0) v d).mp h p
  simp only [zero_mul, Finset.sum_const_zero, add_zero] at hrow
  unfold Ddg.normalDerivative
  rw [eq_div_iff (hℓ p).ne', ← hrow]
  ring

end Recovery

end Cuts

end DdgProofs
