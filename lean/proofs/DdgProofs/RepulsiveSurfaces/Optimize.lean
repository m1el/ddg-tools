import Mathlib

/-!
# Repulsive Surfaces: the L-BFGS two-loop recursion (Sec. 8.2 comparison methods; SPEC D26)

`rs-lbfgs-two-loop`: with stored pairs `(s_i, y_i)` (newest first in the list), `ρ_i = 1/(y_iᵀs_i)`
and an initial operator `H₀` (here `W′⁻¹`, the solve with the constrained cotan Laplacian), the
BFGS inverse update is `H⁺ = Vᵀ H V + ρ s sᵀ`, `V = I − ρ y sᵀ`. The two-loop recursion, written
recursively (the first loop is the descent into the list, newest to oldest, the second loop the
return, oldest to newest),

```
twoLoop H₀ []            q = H₀ q
twoLoop H₀ ((s, y) :: p) q = let α = ρ sᵀq; r = twoLoop H₀ p (q − α y); r + (α − ρ yᵀr) s
```

returns `H_k q` (`twoLoop_eq`). If `H₀` is symmetric positive definite and every `y_iᵀs_i > 0`,
`H_k` is symmetric positive definite (`bfgsInv_posDef`), so `−H_k g` is a descent direction.
-/

noncomputable section

namespace DdgProofs

namespace RS

open Matrix

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- The BFGS inverse-Hessian approximation from `H₀` and the pairs (newest first):
`H⁺ = (I − ρ s yᵀ) H (I − ρ y sᵀ) + ρ s sᵀ`, `ρ = 1/(yᵀs)`. -/
def bfgsInv (H₀ : Matrix n n ℝ) : List ((n → ℝ) × (n → ℝ)) → Matrix n n ℝ
  | [] => H₀
  | (s, y) :: p =>
      (1 - (y ⬝ᵥ s)⁻¹ • vecMulVec s y) * bfgsInv H₀ p * (1 - (y ⬝ᵥ s)⁻¹ • vecMulVec y s)
        + (y ⬝ᵥ s)⁻¹ • vecMulVec s s

/-- The two-loop recursion (recursive form).

[formula: rs-lbfgs-two-loop] -/
def twoLoop (H₀ : Matrix n n ℝ) : List ((n → ℝ) × (n → ℝ)) → (n → ℝ) → (n → ℝ)
  | [], q => H₀ *ᵥ q
  | (s, y) :: p, q =>
      let α := (y ⬝ᵥ s)⁻¹ * (s ⬝ᵥ q)
      let r := twoLoop H₀ p (q - α • y)
      r + (α - (y ⬝ᵥ s)⁻¹ * (y ⬝ᵥ r)) • s

omit [DecidableEq n] in
theorem vecMulVec_mulVec' (a b q : n → ℝ) : vecMulVec a b *ᵥ q = (b ⬝ᵥ q) • a := by
  ext i
  simp [vecMulVec, mulVec, dotProduct, Finset.mul_sum, mul_comm, mul_left_comm]

/-- **The two-loop recursion applies the L-BFGS inverse** (N&W Algorithm 7.4): `twoLoop H₀ p q =
H_k q`.

[formula: rs-lbfgs-two-loop] -/
theorem twoLoop_eq (H₀ : Matrix n n ℝ) (p : List ((n → ℝ) × (n → ℝ))) (q : n → ℝ) :
    twoLoop H₀ p q = bfgsInv H₀ p *ᵥ q := by
  induction p generalizing q with
  | nil => rfl
  | cons sy p ih =>
    obtain ⟨s, y⟩ := sy
    simp only [twoLoop, bfgsInv, ih]
    rw [add_mulVec, ← mulVec_mulVec, ← mulVec_mulVec, sub_mulVec, one_mulVec, smul_mulVec,
      vecMulVec_mulVec', sub_mulVec, one_mulVec, smul_mulVec, vecMulVec_mulVec', smul_mulVec,
      vecMulVec_mulVec']
    set r := bfgsInv H₀ p *ᵥ (q - ((y ⬝ᵥ s)⁻¹ * (s ⬝ᵥ q)) • y)
    have hr : bfgsInv H₀ p *ᵥ (q - (y ⬝ᵥ s)⁻¹ • (s ⬝ᵥ q) • y) = r := by
      rw [smul_smul]
    rw [hr, smul_smul, smul_smul, sub_smul]
    abel

omit [DecidableEq n] in
theorem dot_vecMulVec' (a b x : n → ℝ) : x ⬝ᵥ vecMulVec a b *ᵥ x = (x ⬝ᵥ a) * (b ⬝ᵥ x) := by
  rw [vecMulVec_mulVec', dotProduct_smul, smul_eq_mul, mul_comm]

/-- **The update keeps symmetry and positive definiteness** (N&W Sec. 6.1): if `H₀` is symmetric
positive definite and `y_iᵀ s_i > 0` for every stored pair, `H_k` is symmetric and
`xᵀ H_k x > 0` for `x ≠ 0`; so `−H_k g` is a descent direction for `g ≠ 0`.

[formula: rs-lbfgs-two-loop] -/
theorem bfgsInv_posDef {H₀ : Matrix n n ℝ} (hsym : H₀ᵀ = H₀) (hpd : ∀ x, x ≠ 0 → 0 < x ⬝ᵥ H₀ *ᵥ x)
    (p : List ((n → ℝ) × (n → ℝ))) (hp : ∀ sy ∈ p, 0 < sy.2 ⬝ᵥ sy.1) :
    (bfgsInv H₀ p)ᵀ = bfgsInv H₀ p ∧ ∀ x, x ≠ 0 → 0 < x ⬝ᵥ bfgsInv H₀ p *ᵥ x := by
  induction p with
  | nil => exact ⟨hsym, hpd⟩
  | cons sy p ih =>
    obtain ⟨s, y⟩ := sy
    obtain ⟨ihs, ihp⟩ := ih fun sy h => hp sy (List.mem_cons_of_mem _ h)
    have hys : 0 < y ⬝ᵥ s := hp (s, y) List.mem_cons_self
    set ρ := (y ⬝ᵥ s)⁻¹ with hρ
    have hρpos : 0 < ρ := inv_pos.mpr hys
    set H := bfgsInv H₀ p
    have hT : (1 - ρ • vecMulVec s y)ᵀ = 1 - ρ • vecMulVec y s := by
      rw [transpose_sub, transpose_one, transpose_smul, transpose_vecMulVec]
    have hT' : (1 - ρ • vecMulVec y s)ᵀ = 1 - ρ • vecMulVec s y := by
      rw [transpose_sub, transpose_one, transpose_smul, transpose_vecMulVec]
    constructor
    · simp only [bfgsInv, ← hρ]
      rw [transpose_add, transpose_mul, transpose_mul, hT', ihs, hT, transpose_smul,
        transpose_vecMulVec, Matrix.mul_assoc]
    · intro x hx
      simp only [bfgsInv, ← hρ]
      -- `xᵀH⁺x = (Vx)ᵀ H (Vx) + ρ (sᵀx)²`
      set z := (1 - ρ • vecMulVec y s) *ᵥ x with hz
      have e : x ⬝ᵥ ((1 - ρ • vecMulVec s y) * H * (1 - ρ • vecMulVec y s) + ρ • vecMulVec s s) *ᵥ x
          = z ⬝ᵥ H *ᵥ z + ρ * (s ⬝ᵥ x) ^ 2 := by
        rw [add_mulVec, dotProduct_add, ← mulVec_mulVec, ← mulVec_mulVec, ← hz, ← hT',
          dotProduct_mulVec, vecMul_transpose, ← hz, smul_mulVec, dotProduct_smul, dot_vecMulVec',
          smul_eq_mul, dotProduct_comm x s]
        ring
      rw [e]
      by_cases h0 : s ⬝ᵥ x = 0
      · have hzx : z = x := by
          rw [hz, sub_mulVec, one_mulVec, smul_mulVec, vecMulVec_mulVec', h0]; simp
        rw [h0, hzx]
        simpa using ihp x hx
      · have h1 : 0 ≤ z ⬝ᵥ H *ᵥ z := by
          by_cases hz0 : z = 0
          · rw [hz0]; simp
          · exact (ihp z hz0).le
        have h2 : 0 < ρ * (s ⬝ᵥ x) ^ 2 := mul_pos hρpos (by positivity)
        linarith

end RS

end DdgProofs
