import DDG.Repulsive.Surfaces.Saddle

/-!
# Repulsive Surfaces: the fractional preconditioner (paper Sec. 5.2; SPEC §6.3, D14)

* `rs-preconditioner-order`: the order bookkeeping `(−Δ)⁻¹ (−Δ)^{2−s} (−Δ)⁻¹ = (−Δ)^{−s}` in
  matrix form with stiffness and mass matrices. For invertible `M`, `W` and any family of powers
  `P t` of `Δ = M⁻¹W` with `P (a + b) = P a P b` and `P 1 = Δ`,
  `W⁻¹ (M P(2 − s)) W⁻¹ = P(−s) M⁻¹` (`precond_order`): the mass factors cancel. Such a family
  exists when `Δ = Q D Q⁻¹` with a positive diagonal `D` (`P t = Q D^t Q⁻¹`,
  `powFamily_props`; for symmetric positive definite `M`, `W` this is the generalized eigenbasis).
  That `L_{2−s}` approximates `M Δ^{2−s}`, and the resulting quality of the preconditioner, is cited
  (`rs-preconditioner-approx`).
* `rs-fractional-preconditioner` (SPEC §6.3, D33, the default augmented form
  `P = L′⁻¹ [[F, C_Qᵀ],[C_Q, 0]] L′⁻¹`): equal to `L′⁻¹ + L′⁻¹ blockdiag(F − K, 0) L′⁻¹`,
  symmetric, nonsingular iff the middle factor is, and `C_Q (P [r; b])_x = b` for all `r, b`
  (`aug_precond_props`).
* The printed form `P⁻¹ = L′⁻¹ · blockdiag(L_frac, 0) · L′⁻¹` (`PreconditionerForm::Printed`,
  reference mode) is symmetric, positive semidefinite when `L_frac` is, and its output satisfies
  `C_Q z = 0` (`frac_precond_props`, untagged since D33: the id names the default form).
-/

noncomputable section

namespace DDG

namespace RS

open Matrix Finset

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **The order bookkeeping of the preconditioner** (SPEC §6.3): for invertible `M`, `W` and a
power family `P` of `Δ = M⁻¹W` (`P(a + b) = P a · P b`, `P 1 = Δ`),
`W⁻¹ (M P(2 − s)) W⁻¹ = P(−s) M⁻¹`.

[formula: rs-preconditioner-order] -/
theorem precond_order {M W : Matrix n n ℝ} (hM : IsUnit M.det) (hW : IsUnit W.det)
    (P : ℝ → Matrix n n ℝ) (hadd : ∀ a b, P (a + b) = P a * P b) (h1 : P 1 = M⁻¹ * W) (s : ℝ) :
    W⁻¹ * (M * P (2 - s)) * W⁻¹ = P (-s) * M⁻¹ := by
  have hMinv : IsUnit M⁻¹.det := by rw [det_nonsing_inv]; exact hM.ringInverse
  have hΔ : IsUnit (M⁻¹ * W).det := by rw [det_mul]; exact hMinv.mul hW
  -- `P 0 = 1`
  have h0 : P 0 = 1 := by
    have e : P 1 = P 1 * P 0 := by rw [← hadd, add_zero]
    have := congrArg (fun X => (M⁻¹ * W)⁻¹ * X) e
    simp only [h1, ← Matrix.mul_assoc, nonsing_inv_mul _ hΔ, Matrix.one_mul] at this
    exact this.symm
  -- `W⁻¹ = P(−1) M⁻¹`
  have hWinv : W⁻¹ = P (-1) * M⁻¹ := by
    refine Matrix.inv_eq_left_inv ?_
    rw [Matrix.mul_assoc, ← h1, ← hadd, neg_add_cancel, h0]
  rw [hWinv]
  calc P (-1) * M⁻¹ * (M * P (2 - s)) * (P (-1) * M⁻¹)
      = P (-1) * (M⁻¹ * M) * P (2 - s) * P (-1) * M⁻¹ := by simp only [Matrix.mul_assoc]
    _ = P (-1 + (2 - s) + -1) * M⁻¹ := by
        rw [nonsing_inv_mul _ hM, Matrix.mul_one, hadd, hadd]
    _ = P (-s) * M⁻¹ := by congr 2; ring

/-- The power family `P t = Q D^t Q⁻¹` of a diagonalizable `Δ = Q D Q⁻¹` with positive `D`. -/
def powFamily (Q : Matrix n n ℝ) (d : n → ℝ) (t : ℝ) : Matrix n n ℝ :=
  Q * diagonal (fun i => d i ^ t) * Q⁻¹

/-- **Powers of a diagonalizable matrix with positive eigenvalues** form a power family:
`P(a + b) = P a P b`, `P 1 = Q D Q⁻¹`.

[formula: rs-preconditioner-order] -/
theorem powFamily_props {Q : Matrix n n ℝ} (hQ : IsUnit Q.det) {d : n → ℝ} (hd : ∀ i, 0 < d i) :
    (∀ a b, powFamily Q d (a + b) = powFamily Q d a * powFamily Q d b) ∧
      powFamily Q d 1 = Q * diagonal d * Q⁻¹ := by
  refine ⟨fun a b => ?_, ?_⟩
  · have hQQ : Q⁻¹ * Q = 1 := nonsing_inv_mul _ hQ
    calc powFamily Q d (a + b) = Q * (diagonal (fun i => d i ^ a) * diagonal (fun i => d i ^ b)) *
          Q⁻¹ := by
          rw [powFamily, diagonal_mul_diagonal]
          congr 3; funext i; exact Real.rpow_add (hd i) a b
      _ = Q * diagonal (fun i => d i ^ a) * (Q⁻¹ * Q) * diagonal (fun i => d i ^ b) * Q⁻¹ := by
          rw [hQQ, Matrix.mul_one]; simp only [Matrix.mul_assoc]
      _ = powFamily Q d a * powFamily Q d b := by simp only [powFamily, Matrix.mul_assoc]
  · simp [powFamily, Real.rpow_one]

/-! ## The preconditioner with positional rows -/

section frac

variable {q : Type*} [Fintype q] [DecidableEq q]

/-- `P⁻¹ = L′⁻¹ · blockdiag(L_frac, 0) · L′⁻¹` (SPEC §6.3; the multiplier part is zeroed between
the two solves). The printed form (reference mode). -/
def fracPrecond (Lw : Matrix n n ℝ) (CQ : Matrix q n ℝ) (Lfrac : Matrix n n ℝ) :
    Matrix (n ⊕ q) (n ⊕ q) ℝ :=
  (fromBlocks Lw CQᵀ CQ 0)⁻¹ * fromBlocks Lfrac 0 0 0 * (fromBlocks Lw CQᵀ CQ 0)⁻¹

/-- **The preconditioner is symmetric, PSD, and lands in the positional tangent space**
(SPEC §6.3): for symmetric `W′` and `L_frac` and invertible `L′`, `P⁻¹` is symmetric; if `L_frac`
is PSD, so is `P⁻¹`; and `z = P⁻¹ r` satisfies `C_Q z_x = 0`. This is the printed form; the default form is
`aug_precond_props`. -/
theorem frac_precond_props {Lw : Matrix n n ℝ} {CQ : Matrix q n ℝ} {Lfrac : Matrix n n ℝ}
    (hLw : Lwᵀ = Lw) (hLf : Lfracᵀ = Lfrac) (hL : IsUnit (fromBlocks Lw CQᵀ CQ 0).det) :
    (fracPrecond Lw CQ Lfrac)ᵀ = fracPrecond Lw CQ Lfrac ∧
      (Lfrac.PosSemidef → ∀ r, 0 ≤ r ⬝ᵥ fracPrecond Lw CQ Lfrac *ᵥ r) ∧
      ∀ r, CQ *ᵥ (fun i => (fracPrecond Lw CQ Lfrac *ᵥ r) (Sum.inl i)) = 0 := by
  set L := fromBlocks Lw CQᵀ CQ 0 with hLdef
  have hLs : Lᵀ = L := by
    rw [hLdef, fromBlocks_transpose, hLw, transpose_transpose, transpose_zero]
  have hLinvs : L⁻¹ᵀ = L⁻¹ := by rw [transpose_nonsing_inv, hLs]
  have hD : (fromBlocks Lfrac 0 0 0 : Matrix (n ⊕ q) (n ⊕ q) ℝ)ᵀ = fromBlocks Lfrac 0 0 0 := by
    rw [fromBlocks_transpose, hLf]; simp
  refine ⟨?_, fun hpsd r => ?_, fun r => ?_⟩
  · rw [fracPrecond, ← hLdef, transpose_mul, transpose_mul, hLinvs, hD, Matrix.mul_assoc]
  · -- `rᵀ L⁻¹ D L⁻¹ r = yᵀ D y` with `y = L⁻¹ r`
    set y := L⁻¹ *ᵥ r with hy
    have e : r ⬝ᵥ fracPrecond Lw CQ Lfrac *ᵥ r =
        y ⬝ᵥ fromBlocks Lfrac 0 0 0 *ᵥ y := by
      rw [fracPrecond, ← hLdef, ← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec,
        ← hLinvs, vecMul_transpose, hLinvs]
    rw [e, fromBlocks_mulVec]
    simp only [zero_mulVec, add_zero]
    have : y ⬝ᵥ Sum.elim (Lfrac *ᵥ (y ∘ Sum.inl)) 0 =
        (y ∘ Sum.inl) ⬝ᵥ Lfrac *ᵥ (y ∘ Sum.inl) := by
      simp [dotProduct, Fintype.sum_sum_type]
    rw [this]
    exact hpsd.dotProduct_mulVec_nonneg _
  · -- `L z = D y` and the second block row of `L` is `C_Q z_x`
    have hz : L *ᵥ (fracPrecond Lw CQ Lfrac *ᵥ r) = fromBlocks Lfrac 0 0 0 *ᵥ (L⁻¹ *ᵥ r) := by
      rw [fracPrecond, ← hLdef, ← mulVec_mulVec, ← mulVec_mulVec, mulVec_mulVec,
        mul_nonsing_inv _ hL, one_mulVec]
    have := congrArg (fun v => fun j => v (Sum.inr j)) hz
    simp only [hLdef, fromBlocks_mulVec, Sum.elim_inr, zero_mulVec, add_zero] at this
    exact this

/-- **The augmented preconditioner** (D33): `P = L′⁻¹ N L′⁻¹` with `L′ = [[K, Cᵀ],[C, 0]]` and
`N = [[F, Cᵀ],[C, 0]]`. -/
def augPrecond (K : Matrix n n ℝ) (C : Matrix q n ℝ) (F : Matrix n n ℝ) :
    Matrix (n ⊕ q) (n ⊕ q) ℝ :=
  (fromBlocks K Cᵀ C 0)⁻¹ * fromBlocks F Cᵀ C 0 * (fromBlocks K Cᵀ C 0)⁻¹

/-- **The augmented fractional preconditioner** (SPEC §6.3, D33): for symmetric `K`, `F` and
invertible `L′ = [[K, Cᵀ],[C, 0]]`, with `N = [[F, Cᵀ],[C, 0]]`, `P = L′⁻¹ N L′⁻¹`
(a) equals `L′⁻¹ + L′⁻¹ blockdiag(F − K, 0) L′⁻¹`, (b) is symmetric, (c) is nonsingular iff `N`
is, and (d) reproduces the constraint rows for every right-hand side `[r; b]`:
`C (P [r; b])_x = b` (no further hypothesis).

[formula: rs-fractional-preconditioner] -/
theorem aug_precond_props {K F : Matrix n n ℝ} {C : Matrix q n ℝ}
    (hK : Kᵀ = K) (hF : Fᵀ = F) (hL : IsUnit (fromBlocks K Cᵀ C 0).det) :
    augPrecond K C F = (fromBlocks K Cᵀ C 0)⁻¹ +
        (fromBlocks K Cᵀ C 0)⁻¹ * fromBlocks (F - K) 0 0 0 * (fromBlocks K Cᵀ C 0)⁻¹ ∧
      (augPrecond K C F)ᵀ = augPrecond K C F ∧
      (IsUnit (augPrecond K C F).det ↔ IsUnit (fromBlocks F Cᵀ C 0).det) ∧
      ∀ (r : n → ℝ) (b : q → ℝ),
        C *ᵥ (fun i => (augPrecond K C F *ᵥ Sum.elim r b) (Sum.inl i)) = b := by
  set L := fromBlocks K Cᵀ C 0 with hLdef
  set N := fromBlocks F Cᵀ C 0 with hNdef
  have hLs : Lᵀ = L := by
    rw [hLdef, fromBlocks_transpose, hK, transpose_transpose, transpose_zero]
  have hNs : Nᵀ = N := by
    rw [hNdef, fromBlocks_transpose, hF, transpose_transpose, transpose_zero]
  have hLinvs : L⁻¹ᵀ = L⁻¹ := by rw [transpose_nonsing_inv, hLs]
  have hLi : IsUnit L⁻¹.det := by rw [det_nonsing_inv]; exact hL.ringInverse
  refine ⟨?_, ?_, ?_, fun r b => ?_⟩
  · have hN : N = L + fromBlocks (F - K) 0 0 0 := by
      rw [hLdef, hNdef, fromBlocks_add]; simp
    have hLLL : L⁻¹ * L * L⁻¹ = L⁻¹ := by rw [nonsing_inv_mul _ hL, Matrix.one_mul]
    rw [augPrecond, ← hLdef, ← hNdef, hN, Matrix.mul_add, Matrix.add_mul, hLLL]
  · rw [augPrecond, ← hLdef, ← hNdef, transpose_mul, transpose_mul, hLinvs, hNs, Matrix.mul_assoc]
  · rw [augPrecond, ← hLdef, ← hNdef]
    simp only [Matrix.det_mul, IsUnit.mul_iff, hLi, and_true, true_and]
  · set v : n ⊕ q → ℝ := Sum.elim r b with hv
    set y := L⁻¹ *ᵥ v with hy
    have hLy : L *ᵥ y = v := by rw [hy, mulVec_mulVec, mul_nonsing_inv _ hL, one_mulVec]
    have hz : L *ᵥ (augPrecond K C F *ᵥ v) = N *ᵥ y := by
      rw [augPrecond, ← hLdef, ← hNdef, ← mulVec_mulVec, ← mulVec_mulVec, mulVec_mulVec,
        mul_nonsing_inv _ hL, one_mulVec]
    have h1 := congrArg (fun w => fun j => w (Sum.inr j)) hz
    have h2 := congrArg (fun w => fun j => w (Sum.inr j)) hLy
    simp only [hLdef, hNdef, fromBlocks_mulVec, Sum.elim_inr, zero_mulVec, add_zero, hv] at h1 h2
    exact h1.trans h2

end frac

end RS

end DDG
