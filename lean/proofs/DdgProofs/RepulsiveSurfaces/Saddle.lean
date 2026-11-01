import DdgProofs.RepulsiveCurves.Saddle

/-!
# Repulsive Surfaces: the saddle point systems and the Schur complement (paper Eq. 14, 15, 18–21,
Sec. 5.3, 7.2; SPEC §6.1, §6.2)

Generic in the index types: `A : Matrix n n ℝ` is the (constrained) metric — `A₃`, or the
operator `A′` with the positional rows inside (`rs-nested-saddle`) — and `C : Matrix k n ℝ` the
Jacobian of the nonlinear constraints. The Schur complement is `S = −C A⁻¹ Cᵀ` (Eq. 19).

* `rs-schur-direction` (Eq. 18–20): the block inverse of `M = [A Cᵀ; C 0]` printed in Sec. 5.3
  (`M · M⁻¹ = 1`); the direction `x = A⁻¹g + A⁻¹CᵀS⁻¹CA⁻¹g` with `λ = −S⁻¹CA⁻¹g` solves
  `Ax + Cᵀλ = g`, `Cx = 0`. Scaling the rows of `C` by an invertible `D` (the reference's
  volume row `Σ a_S N_S = 3∇V`, D9) leaves `x` unchanged.
* `rs-schur-correction` (Eq. 15, 21): `h = A⁻¹CᵀS⁻¹Φ` with `μ = −S⁻¹Φ` solves `Ah + Cᵀμ = 0`,
  `Ch = −Φ`: it is the paper's Newton step, to be **added** (`f ← f + h`), and the `A`-smallest
  solution of `Ch = −Φ`. With rows scaled by `D`, the step becomes `A⁻¹CᵀS⁻¹D⁻¹Φ` (a third of
  the correction for the reference's 3× volume row, D9).

  **SPEC §6.2 has the sign wrong**: it says "`h = A′⁻¹ C_Nᵀ S⁻¹ Φ(f)` subtracted
  (`rs-schur-correction`: `C_N h = Φ`)". With `S = −C A⁻¹ Cᵀ` that `h` has `C h = −Φ`, so
  subtracting it doubles the linearized violation (`Φ − C h = 2Φ`, `schur_correction_sign`).
  Either add it (paper Eq. 21: `h = −A⁻¹CᵀS⁻¹(−Φ)`, "which we add to the updated embedding") or
  subtract `A⁻¹Cᵀ(CA⁻¹Cᵀ)⁻¹Φ = −h`.
* `rs-nested-saddle` (Sec. 7.2): with the positional rows `C_Q` inside
  `A′ = [A₃ C_Qᵀ; C_Q 0]` and `C̃_N = [C_N 0]`, the saddle `[A′ C̃_Nᵀ; C̃_N 0]` is the full system
  `[A₃ C_Qᵀ C_Nᵀ; C_Q 0 0; C_N 0 0]` (`nested_saddle_iff`); so the Schur direction with `A′`
  satisfies every constraint row (`nested_direction`).
-/

noncomputable section

namespace DdgProofs

namespace RS

open Matrix Finset

variable {n k : Type*} [Fintype n] [Fintype k] [DecidableEq n] [DecidableEq k]

/-- The Schur complement `M/A = −C A⁻¹ Cᵀ` (Eq. 19). -/
def schur (A : Matrix n n ℝ) (C : Matrix k n ℝ) : Matrix k k ℝ := -(C * A⁻¹ * Cᵀ)

/-- The block inverse of `M = [A Cᵀ; C 0]` printed in Sec. 5.3. -/
def schurInv (A : Matrix n n ℝ) (C : Matrix k n ℝ) : Matrix (n ⊕ k) (n ⊕ k) ℝ :=
  fromBlocks (A⁻¹ + A⁻¹ * Cᵀ * (schur A C)⁻¹ * C * A⁻¹) (-(A⁻¹ * Cᵀ * (schur A C)⁻¹))
    (-((schur A C)⁻¹ * C * A⁻¹)) (schur A C)⁻¹

/-- **The block inverse of Sec. 5.3**: for invertible `A` and `S = −C A⁻¹ Cᵀ`,
`[A Cᵀ; C 0] · [A⁻¹ + A⁻¹CᵀS⁻¹CA⁻¹, −A⁻¹CᵀS⁻¹; −S⁻¹CA⁻¹, S⁻¹] = 1`.

[formula: rs-schur-direction] -/
theorem schur_block_inverse {A : Matrix n n ℝ} {C : Matrix k n ℝ} (hA : IsUnit A.det)
    (hS : IsUnit (schur A C).det) :
    fromBlocks A Cᵀ C 0 * schurInv A C = 1 := by
  have hAA : A * A⁻¹ = 1 := mul_nonsing_inv A hA
  have hSS : schur A C * (schur A C)⁻¹ = 1 := mul_nonsing_inv _ hS
  have hCAC : C * A⁻¹ * Cᵀ = -schur A C := by rw [schur, neg_neg]
  rw [schurInv, fromBlocks_multiply, ← fromBlocks_one]
  congr 1
  · -- A(A⁻¹ + A⁻¹CᵀS⁻¹CA⁻¹) − CᵀS⁻¹CA⁻¹ = 1
    rw [Matrix.mul_add, hAA]
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc A A⁻¹, hAA, Matrix.one_mul, Matrix.mul_neg]
    abel
  · simp only [Matrix.mul_neg, ← Matrix.mul_assoc, hAA, Matrix.one_mul, neg_add_cancel]
  · -- C(A⁻¹ + A⁻¹CᵀS⁻¹CA⁻¹) − 0 = CA⁻¹ + (CA⁻¹Cᵀ)S⁻¹CA⁻¹ = 0
    rw [Matrix.zero_mul, add_zero, Matrix.mul_add]
    have : C * (A⁻¹ * Cᵀ * (schur A C)⁻¹ * C * A⁻¹) = -(C * A⁻¹) := by
      calc C * (A⁻¹ * Cᵀ * (schur A C)⁻¹ * C * A⁻¹)
          = (C * A⁻¹ * Cᵀ) * (schur A C)⁻¹ * C * A⁻¹ := by simp only [Matrix.mul_assoc]
        _ = -(C * A⁻¹) := by
            rw [hCAC, Matrix.neg_mul, Matrix.neg_mul, Matrix.neg_mul, hSS, Matrix.one_mul]
    rw [this, add_neg_cancel]
  · rw [Matrix.zero_mul, add_zero, Matrix.mul_neg]
    have : C * (A⁻¹ * Cᵀ * (schur A C)⁻¹) = (C * A⁻¹ * Cᵀ) * (schur A C)⁻¹ := by
      simp only [Matrix.mul_assoc]
    rw [this, hCAC, Matrix.neg_mul, neg_neg, hSS]

omit [DecidableEq n] [DecidableEq k] in
/-- Evaluating `M (Sum.elim x λ)` blockwise. -/
theorem saddle_mulVec (A : Matrix n n ℝ) (C : Matrix k n ℝ) (x : n → ℝ) (μ : k → ℝ) :
    fromBlocks A Cᵀ C 0 *ᵥ Sum.elim x μ = Sum.elim (A *ᵥ x + Cᵀ *ᵥ μ) (C *ᵥ x) := by
  rw [fromBlocks_mulVec]
  simp only [Sum.elim_comp_inl, Sum.elim_comp_inr, zero_mulVec, add_zero]

/-- **The constrained direction** (Eq. 20, SPEC §6.1): `x = A⁻¹g + A⁻¹CᵀS⁻¹CA⁻¹g` and
`λ = −S⁻¹CA⁻¹g` solve `A x + Cᵀ λ = g`, `C x = 0`.

[formula: rs-schur-direction] -/
theorem schur_direction {A : Matrix n n ℝ} {C : Matrix k n ℝ} (hA : IsUnit A.det)
    (hS : IsUnit (schur A C).det) (g : n → ℝ) :
    let x := (A⁻¹ + A⁻¹ * Cᵀ * (schur A C)⁻¹ * C * A⁻¹) *ᵥ g
    let μ := -((schur A C)⁻¹ * C * A⁻¹) *ᵥ g
    A *ᵥ x + Cᵀ *ᵥ μ = g ∧ C *ᵥ x = 0 := by
  intro x μ
  have h := congrArg (fun M => M *ᵥ Sum.elim g (0 : k → ℝ)) (schur_block_inverse hA hS)
  have hv : schurInv A C *ᵥ Sum.elim g (0 : k → ℝ) = Sum.elim x μ := by
    rw [schurInv, fromBlocks_mulVec]
    simp only [Sum.elim_comp_inl, Sum.elim_comp_inr, mulVec_zero, add_zero]
    rfl
  simp only [one_mulVec] at h
  rw [← mulVec_mulVec, hv, saddle_mulVec] at h
  have e1 := congrArg (fun v => fun i => v (Sum.inl i)) h
  have e2 := congrArg (fun v => fun i => v (Sum.inr i)) h
  simp only [Sum.elim_inl, Sum.elim_inr] at e1 e2
  exact ⟨e1, e2⟩

/-- **The Newton correction** (Eq. 15, 21, SPEC §6.2): `h = A⁻¹CᵀS⁻¹Φ` and `μ = −S⁻¹Φ` solve
`A h + Cᵀ μ = 0`, `C h = −Φ`; so `f + h` satisfies the linearized constraints
`Φ + C h = 0` (the step is **added**, paper Eq. 21).

[formula: rs-schur-correction] -/
theorem schur_correction {A : Matrix n n ℝ} {C : Matrix k n ℝ} (hA : IsUnit A.det)
    (hS : IsUnit (schur A C).det) (Φ : k → ℝ) :
    let h := (A⁻¹ * Cᵀ * (schur A C)⁻¹) *ᵥ Φ
    let μ := -(schur A C)⁻¹ *ᵥ Φ
    A *ᵥ h + Cᵀ *ᵥ μ = 0 ∧ C *ᵥ h = -Φ := by
  intro h μ
  have e := congrArg (fun M => M *ᵥ Sum.elim (0 : n → ℝ) (-Φ)) (schur_block_inverse hA hS)
  have hv : schurInv A C *ᵥ Sum.elim (0 : n → ℝ) (-Φ) = Sum.elim h μ := by
    rw [schurInv, fromBlocks_mulVec]
    simp only [Sum.elim_comp_inl, Sum.elim_comp_inr, mulVec_zero, zero_add, h, μ, neg_mulVec,
      mulVec_neg, neg_neg]
  simp only [one_mulVec] at e
  rw [← mulVec_mulVec, hv, saddle_mulVec] at e
  have e1 := congrArg (fun v => fun i => v (Sum.inl i)) e
  have e2 := congrArg (fun v => fun i => v (Sum.inr i)) e
  simp only [Sum.elim_inl, Sum.elim_inr] at e1 e2
  exact ⟨e1, e2⟩

/-- **The sign of the correction** (SPEC §6.2 error): with `h = A⁻¹CᵀS⁻¹Φ`, `Φ + C h = 0` but
`Φ − C h = 2Φ`: subtracting `h` (as SPEC §6.2 prints) doubles the linearized violation.

[formula: rs-schur-correction] -/
theorem schur_correction_sign {A : Matrix n n ℝ} {C : Matrix k n ℝ} (hA : IsUnit A.det)
    (hS : IsUnit (schur A C).det) (Φ : k → ℝ) :
    Φ + C *ᵥ ((A⁻¹ * Cᵀ * (schur A C)⁻¹) *ᵥ Φ) = 0 ∧
      Φ - C *ᵥ ((A⁻¹ * Cᵀ * (schur A C)⁻¹) *ᵥ Φ) = (2 : ℝ) • Φ := by
  have := (schur_correction hA hS Φ).2
  rw [this]
  constructor
  · abel
  · rw [sub_neg_eq_add, two_smul]

/-- **The correction is the `A`-smallest solution** of `C y = −Φ` (Sec. 3.3.2, "least-norm
solution with respect to the metric encoded by `A₃`"), for symmetric `A`.

[formula: rs-schur-correction] -/
theorem schur_correction_least_norm {A : Matrix n n ℝ} {C : Matrix k n ℝ} (hA : IsUnit A.det)
    (hAs : Aᵀ = A) (hS : IsUnit (schur A C).det) (Φ : k → ℝ) {y : n → ℝ}
    (hy : C *ᵥ y = -Φ) :
    let h := (A⁻¹ * Cᵀ * (schur A C)⁻¹) *ᵥ Φ
    y ⬝ᵥ A *ᵥ y = h ⬝ᵥ A *ᵥ h + (y - h) ⬝ᵥ A *ᵥ (y - h) := by
  intro h
  obtain ⟨h1, h2⟩ := schur_correction hA hS Φ
  exact RC.chord_least_norm hAs ⟨h1, h2⟩ hy

/-- **Row scaling** (D9): replacing `C` by `D C` for an invertible `D` leaves the direction
unchanged and turns the correction `A⁻¹CᵀS⁻¹Φ` into `A⁻¹CᵀS⁻¹D⁻¹Φ` (for the reference's
`D = 3` on the volume row, one Newton step removes a third of the volume error).

[formula: rs-schur-correction] -/
theorem schur_row_scaling {A : Matrix n n ℝ} {C : Matrix k n ℝ} {D : Matrix k k ℝ}
    (hD : IsUnit D.det) :
    (A⁻¹ * (D * C)ᵀ * (schur A (D * C))⁻¹ * (D * C) * A⁻¹ =
        A⁻¹ * Cᵀ * (schur A C)⁻¹ * C * A⁻¹) ∧
      A⁻¹ * (D * C)ᵀ * (schur A (D * C))⁻¹ = A⁻¹ * Cᵀ * (schur A C)⁻¹ * D⁻¹ := by
  have hsch : schur A (D * C) = D * schur A C * Dᵀ := by
    simp only [schur, Matrix.transpose_mul, Matrix.mul_neg, Matrix.neg_mul, Matrix.mul_assoc]
  have hDt : IsUnit Dᵀ.det := by rwa [det_transpose]
  have hinv : (schur A (D * C))⁻¹ = (Dᵀ)⁻¹ * (schur A C)⁻¹ * D⁻¹ := by
    rw [hsch, Matrix.mul_inv_rev, Matrix.mul_inv_rev, Matrix.mul_assoc]
  have hDDt : Dᵀ * (Dᵀ)⁻¹ = 1 := mul_nonsing_inv _ hDt
  have hDD : D⁻¹ * D = 1 := nonsing_inv_mul _ hD
  rw [hinv, Matrix.transpose_mul]
  constructor
  · calc A⁻¹ * (Cᵀ * Dᵀ) * ((Dᵀ)⁻¹ * (schur A C)⁻¹ * D⁻¹) * (D * C) * A⁻¹
        = A⁻¹ * Cᵀ * (Dᵀ * (Dᵀ)⁻¹) * (schur A C)⁻¹ * (D⁻¹ * D) * C * A⁻¹ := by
          simp only [Matrix.mul_assoc]
      _ = A⁻¹ * Cᵀ * (schur A C)⁻¹ * C * A⁻¹ := by rw [hDDt, hDD]; simp only [Matrix.mul_one]
  · calc A⁻¹ * (Cᵀ * Dᵀ) * ((Dᵀ)⁻¹ * (schur A C)⁻¹ * D⁻¹)
        = A⁻¹ * Cᵀ * (Dᵀ * (Dᵀ)⁻¹) * (schur A C)⁻¹ * D⁻¹ := by simp only [Matrix.mul_assoc]
      _ = A⁻¹ * Cᵀ * (schur A C)⁻¹ * D⁻¹ := by rw [hDDt, Matrix.mul_one]

/-! ## The nested saddle (Sec. 7.2) -/

section nested

variable {q : Type*} [Fintype q] [DecidableEq q]

/-- `A′ = [A₃ C_Qᵀ; C_Q 0]`: the positional rows inside the operator. -/
def innerSaddle (A : Matrix n n ℝ) (CQ : Matrix q n ℝ) : Matrix (n ⊕ q) (n ⊕ q) ℝ :=
  fromBlocks A CQᵀ CQ 0

/-- `C̃_N = [C_N 0]`: the nonlinear rows padded with zeros on the multiplier columns. -/
def padRows (CN : Matrix k n ℝ) : Matrix k (n ⊕ q) ℝ := fun r j => Sum.elim (CN r) 0 j

omit [DecidableEq k] [DecidableEq n] [DecidableEq q] in
/-- **The nested saddle is the full saddle** (Sec. 7.2, SPEC §6.1):
`[A′ C̃_Nᵀ; C̃_N 0] (x, λ_Q, λ_N) = (b, r, s)` iff `A₃x + C_Qᵀλ_Q + C_Nᵀλ_N = b`, `C_Q x = r`,
`C_N x = s`.

[formula: rs-nested-saddle] -/
theorem nested_saddle_iff (A : Matrix n n ℝ) (CQ : Matrix q n ℝ) (CN : Matrix k n ℝ)
    (x : n → ℝ) (lQ : q → ℝ) (lN : k → ℝ) (b : n → ℝ) (r : q → ℝ) (s : k → ℝ) :
    fromBlocks (innerSaddle A CQ) (padRows (q := q) CN)ᵀ (padRows CN) 0 *ᵥ
        Sum.elim (Sum.elim x lQ) lN = Sum.elim (Sum.elim b r) s ↔
      A *ᵥ x + CQᵀ *ᵥ lQ + CNᵀ *ᵥ lN = b ∧ CQ *ᵥ x = r ∧ CN *ᵥ x = s := by
  have hT : (padRows (q := q) CN)ᵀ *ᵥ lN = Sum.elim (CNᵀ *ᵥ lN) 0 := by
    ext j; cases j <;> simp [padRows, mulVec, dotProduct, transpose_apply]
  have hC : padRows (q := q) CN *ᵥ Sum.elim x lQ = CN *ᵥ x := by
    ext r; simp [padRows, mulVec, dotProduct, Fintype.sum_sum_type]
  rw [fromBlocks_mulVec]
  simp only [Sum.elim_comp_inl, Sum.elim_comp_inr, zero_mulVec, add_zero]
  rw [innerSaddle, saddle_mulVec, hT, hC]
  constructor
  · intro h
    have h1 := congrArg (fun v => fun i => v (Sum.inl (Sum.inl i))) h
    have h2 := congrArg (fun v => fun i => v (Sum.inl (Sum.inr i))) h
    have h3 := congrArg (fun v => fun i => v (Sum.inr i)) h
    simp only [Sum.elim_inl, Sum.elim_inr, Pi.add_apply] at h1 h2 h3
    refine ⟨?_, ?_, h3⟩
    · funext i; have := congrFun h1 i; simpa using this
    · funext i; have := congrFun h2 i; simpa using this
  · rintro ⟨h1, h2, h3⟩
    rw [h3]
    congr 1
    ext i
    cases i with
    | inl i => simp [← h1]
    | inr i => simp [← h2]

/-- **The Schur direction on `A′` satisfies every row** (SPEC §6.1): if `A′` and the Schur
complement of `C̃_N` are invertible, the direction of `schur_direction` (with right side
`(g, 0)`) solves the full saddle `A₃x + C_Qᵀλ_Q + C_Nᵀλ_N = g`, `C_Q x = 0`, `C_N x = 0`.

[formula: rs-nested-saddle] -/
theorem nested_direction {A : Matrix n n ℝ} {CQ : Matrix q n ℝ} {CN : Matrix k n ℝ}
    (hA : IsUnit (innerSaddle A CQ).det)
    (hS : IsUnit (schur (innerSaddle A CQ) (padRows (q := q) CN)).det) (g : n → ℝ) :
    let A' := innerSaddle A CQ
    let C := padRows (q := q) CN
    let z := (A'⁻¹ + A'⁻¹ * Cᵀ * (schur A' C)⁻¹ * C * A'⁻¹) *ᵥ Sum.elim g 0
    let μ := -((schur A' C)⁻¹ * C * A'⁻¹) *ᵥ Sum.elim g 0
    A *ᵥ (fun i => z (Sum.inl i)) + CQᵀ *ᵥ (fun j => z (Sum.inr j)) + CNᵀ *ᵥ μ = g ∧
      CQ *ᵥ (fun i => z (Sum.inl i)) = 0 ∧ CN *ᵥ (fun i => z (Sum.inl i)) = 0 := by
  intro A' C z μ
  obtain ⟨h1, h2⟩ := schur_direction hA hS (Sum.elim g 0)
  have hz : z = Sum.elim (fun i => z (Sum.inl i)) (fun j => z (Sum.inr j)) := by
    ext i; cases i <;> rfl
  have := (nested_saddle_iff A CQ CN (fun i => z (Sum.inl i)) (fun j => z (Sum.inr j)) μ g 0 0).mp
  apply this
  rw [← hz]
  change fromBlocks A' Cᵀ C 0 *ᵥ Sum.elim z μ = _
  rw [saddle_mulVec]
  rw [h1, h2]

end nested

end RS

end DdgProofs
