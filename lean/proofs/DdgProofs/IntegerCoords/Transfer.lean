import Mathlib.Data.Matrix.Mul
import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.Abel
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.FinCases
import Mathlib.LinearAlgebra.Matrix.Symmetric
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Algebra.Group.Fin.Basic
import Mathlib.Tactic.LinearCombination
import Mathlib.Data.Fin.VecNotation

/-!
# Integer coordinates: attribute transfer (tier A)

`crates/ddg-integer-coords/SPEC.md` §6; the paper Sec. 4.3, Eq. 11–14.

`P₀ : S × V⁰` and `P₁ : S × V¹` interpolate functions of `T⁰` and `T¹` to the vertices of the
common subdivision `S`, `M` is the Galerkin mass matrix of `S`. The transfer of `f` on `V¹` is a
minimizer `f̂` of `‖f − f̂‖² = (P₁f − P₀f̂)ᵀ M (P₁f − P₀f̂)` (Eq. 13).

* **`ic-transfer-normal-equations`** — for `M` symmetric positive semidefinite, `g` minimizes Eq. 13
  if and only if it solves the normal equations `P₀ᵀ M P₀ g = P₀ᵀ M P₁ f` (Eq. 14 with
  `P₀ᵀ M P₀` for `M₀`; `transfer_minimizer_iff`); if `P₀ᵀ M P₀` is positive definite the solution is
  unique (`transfer_normal_unique`). The energy of `g + h` is that of `g` plus
  `(P₀h)ᵀ M (P₀h)` when `g` solves them (`transfer_energy_add`).
* **`ic-transfer-mass-identity`** (`P₀ᵀ M_S P₀ = M₀`), the exact finite part. A piece of `S` inside
  the input triangle `T` is a triangle whose corners have barycentric rows `β r` with respect to
  `T`; the hat functions of `T` restricted to the piece are affine with the values `β r a` at its
  corners, so `P₀` represents them exactly, and the contribution of the piece to
  `(P₀ᵀ M_S P₀)_{ab}` is `pieceMass β a b · Area(T)` (`pieceMass`: `det β / 12 Σ_{r,s} (1 + δ_rs)
  β_ra β_sb`, the P1 mass formula of the piece, area `det β · Area(T)`, applied to the values of
  `λ_a`, `λ_b`). Proved: the piece equal to `T` gives `M₀ / Area(T)`, i.e. `(1 + δ_ab)/12`
  (`pieceMass_one`); `pieceMass` is additive under splitting a (signed) triangle at any point of its
  plane into three (`pieceMass_split`), hence under splitting an edge (`pieceMass_split_edge`) and
  over any fan of a convex polygon, independently of the apex (`pieceMass_fan_apex`). **Cited**,
  not formalized: that the pieces of the common subdivision tile `T` (a planar subdivision of a
  triangle into triangles has, as signed chains, the triangle as its sum, so the valuation
  `pieceMass` sums to `(1 + δ_ab)/12`), and the exactness of the P1 mass formula as an integral
  (`∫_T λ_a λ_b = Area/12 (1 + δ_ab)`); both are standard finite-element facts (e.g. Ciarlet,
  *The Finite Element Method for Elliptic Problems*, 1978, Sec. 4.1), the latter checked numerically
  in Rust.
-/

namespace DdgProofs

namespace IntegerCoords

open Matrix

/-! ## Normal equations (Eq. 13–14) -/

section Normal

variable {S V W : Type*} [Fintype S] [Fintype V] [Fintype W]

/-- The residual `P₁ f − P₀ g` on the common subdivision. -/
def transferResidual (P0 : Matrix S V ℝ) (P1 : Matrix S W ℝ) (f : W → ℝ) (g : V → ℝ) : S → ℝ :=
  P1 *ᵥ f - P0 *ᵥ g

/-- Eq. 13: `(P₁f − P₀g)ᵀ M (P₁f − P₀g)`. -/
def transferEnergy (M : Matrix S S ℝ) (P0 : Matrix S V ℝ) (P1 : Matrix S W ℝ) (f : W → ℝ)
    (g : V → ℝ) : ℝ :=
  transferResidual P0 P1 f g ⬝ᵥ (M *ᵥ transferResidual P0 P1 f g)

lemma quad_sub (M : Matrix S S ℝ) (hM : M.IsSymm) (r w : S → ℝ) (t : ℝ) :
    (r - t • w) ⬝ᵥ (M *ᵥ (r - t • w)) =
      r ⬝ᵥ (M *ᵥ r) - 2 * t * (w ⬝ᵥ (M *ᵥ r)) + t ^ 2 * (w ⬝ᵥ (M *ᵥ w)) := by
  have hsym : r ⬝ᵥ (M *ᵥ w) = w ⬝ᵥ (M *ᵥ r) := by
    rw [dotProduct_mulVec, ← mulVec_transpose, hM.eq, dotProduct_comm]
  simp only [mulVec_sub, mulVec_smul, sub_dotProduct, dotProduct_sub, smul_dotProduct,
    dotProduct_smul, smul_eq_mul, hsym]
  ring

omit [Fintype S] in
lemma residual_add (P0 : Matrix S V ℝ) (P1 : Matrix S W ℝ) (f : W → ℝ) (g h : V → ℝ) (t : ℝ) :
    transferResidual P0 P1 f (g + t • h) = transferResidual P0 P1 f g - t • (P0 *ᵥ h) := by
  simp only [transferResidual, mulVec_add, mulVec_smul]
  abel

/-- `(P₀ h)ᵀ M r = hᵀ (P₀ᵀ M r)`. -/
lemma dot_P0 (M : Matrix S S ℝ) (P0 : Matrix S V ℝ) (h : V → ℝ) (r : S → ℝ) :
    (P0 *ᵥ h) ⬝ᵥ (M *ᵥ r) = h ⬝ᵥ (P0ᵀ *ᵥ (M *ᵥ r)) := by
  rw [dotProduct_comm, dotProduct_mulVec, ← mulVec_transpose, dotProduct_comm]

/-- The normal equations say that `P₀ᵀ M` annihilates the residual. -/
lemma normal_iff_residual (M : Matrix S S ℝ) (P0 : Matrix S V ℝ) (P1 : Matrix S W ℝ) (f : W → ℝ)
    (g : V → ℝ) :
    (P0ᵀ * M * P0) *ᵥ g = (P0ᵀ * M * P1) *ᵥ f ↔
      P0ᵀ *ᵥ (M *ᵥ transferResidual P0 P1 f g) = 0 := by
  simp only [transferResidual, mulVec_sub, ← mulVec_mulVec, Matrix.mul_assoc]
  rw [sub_eq_zero, eq_comm]

/-- A real quadratic `a t + b t²` that is never negative has `a = 0`. -/
lemma linear_coeff_zero {a b : ℝ} (hb : 0 ≤ b) (h : ∀ t : ℝ, 0 ≤ a * t + b * t ^ 2) : a = 0 := by
  by_contra ha
  have h1 := h (-a / (b + 1))
  have hb1 : 0 < b + 1 := by linarith
  have : a * (-a / (b + 1)) + b * (-a / (b + 1)) ^ 2 = -(a ^ 2 / (b + 1) ^ 2) := by
    field_simp; ring
  rw [this] at h1
  have : 0 < a ^ 2 / (b + 1) ^ 2 := by positivity
  linarith

/-- If `g` solves the normal equations, the energy of `g + h` is that of `g` plus
`(P₀h)ᵀ M (P₀h)`. [formula: ic-transfer-normal-equations] -/
theorem transfer_energy_add (M : Matrix S S ℝ) (hM : M.IsSymm) (P0 : Matrix S V ℝ)
    (P1 : Matrix S W ℝ) (f : W → ℝ) (g : V → ℝ)
    (hn : (P0ᵀ * M * P0) *ᵥ g = (P0ᵀ * M * P1) *ᵥ f) (h : V → ℝ) :
    transferEnergy M P0 P1 f (g + h) =
      transferEnergy M P0 P1 f g + (P0 *ᵥ h) ⬝ᵥ (M *ᵥ (P0 *ᵥ h)) := by
  have e := residual_add P0 P1 f g h 1
  rw [one_smul] at e
  unfold transferEnergy
  rw [e, quad_sub M hM, dot_P0, (normal_iff_residual M P0 P1 f g).1 hn, dotProduct_zero]
  ring

/-- **Eq. 14 characterizes the minimizers of Eq. 13**: for `M` symmetric positive semidefinite,
`g` minimizes `(P₁f − P₀g)ᵀ M (P₁f − P₀g)` if and only if `P₀ᵀ M P₀ g = P₀ᵀ M P₁ f`.
[formula: ic-transfer-normal-equations] -/
theorem transfer_minimizer_iff (M : Matrix S S ℝ) (hM : M.IsSymm)
    (hpsd : ∀ x : S → ℝ, 0 ≤ x ⬝ᵥ (M *ᵥ x)) (P0 : Matrix S V ℝ) (P1 : Matrix S W ℝ)
    (f : W → ℝ) (g : V → ℝ) :
    (∀ g' : V → ℝ, transferEnergy M P0 P1 f g ≤ transferEnergy M P0 P1 f g') ↔
      (P0ᵀ * M * P0) *ᵥ g = (P0ᵀ * M * P1) *ᵥ f := by
  classical
  constructor
  · intro hmin
    rw [normal_iff_residual]
    funext v
    -- the derivative of the energy in the direction of the unit vector at `v` vanishes
    let h : V → ℝ := fun w => if w = v then 1 else 0
    have key : ∀ t : ℝ, 0 ≤ (-2 * ((P0 *ᵥ h) ⬝ᵥ (M *ᵥ transferResidual P0 P1 f g))) * t +
        ((P0 *ᵥ h) ⬝ᵥ (M *ᵥ (P0 *ᵥ h))) * t ^ 2 := by
      intro t
      have := hmin (g + t • h)
      unfold transferEnergy at this
      rw [residual_add, quad_sub M hM] at this
      linarith
    have h0 := linear_coeff_zero (hpsd _) key
    rw [dot_P0] at h0
    have : h ⬝ᵥ (P0ᵀ *ᵥ (M *ᵥ transferResidual P0 P1 f g)) =
        (P0ᵀ *ᵥ (M *ᵥ transferResidual P0 P1 f g)) v := by
      simp [h, dotProduct]
    rw [this] at h0
    simp only [Pi.zero_apply]
    linarith
  · intro hn g'
    have := transfer_energy_add M hM P0 P1 f g hn (g' - g)
    rw [add_sub_cancel] at this
    rw [this]
    linarith [hpsd (P0 *ᵥ (g' - g))]

/-- **The transfer is unique** when `P₀ᵀ M P₀` is positive definite (as the paper's `M₀` is).
[formula: ic-transfer-normal-equations] -/
theorem transfer_normal_unique (M : Matrix S S ℝ) (P0 : Matrix S V ℝ) (P1 : Matrix S W ℝ)
    (hpd : ∀ x : V → ℝ, x ≠ 0 → 0 < (P0 *ᵥ x) ⬝ᵥ (M *ᵥ (P0 *ᵥ x))) (f : W → ℝ) {g g' : V → ℝ}
    (hg : (P0ᵀ * M * P0) *ᵥ g = (P0ᵀ * M * P1) *ᵥ f)
    (hg' : (P0ᵀ * M * P0) *ᵥ g' = (P0ᵀ * M * P1) *ᵥ f) : g = g' := by
  by_contra hne
  have hx : g - g' ≠ 0 := sub_ne_zero.2 hne
  have h0 : (P0ᵀ * M * P0) *ᵥ (g - g') = 0 := by rw [mulVec_sub, hg, hg', sub_self]
  have h0' : P0ᵀ *ᵥ (M *ᵥ (P0 *ᵥ (g - g'))) = 0 := by
    rw [mulVec_mulVec, mulVec_mulVec]; exact h0
  have := hpd _ hx
  rw [dot_P0, h0', dotProduct_zero] at this
  exact lt_irrefl _ this

end Normal

/-! ## The mass identity, one input triangle at a time -/

/-- Determinant of three barycentric rows `β 0, β 1, β 2`. -/
def bdet (β : Fin 3 → Fin 3 → ℝ) : ℝ :=
  β 0 0 * (β 1 1 * β 2 2 - β 1 2 * β 2 1) - β 0 1 * (β 1 0 * β 2 2 - β 1 2 * β 2 0) +
    β 0 2 * (β 1 0 * β 2 1 - β 1 1 * β 2 0)

/-- The contribution of a piece with barycentric corner rows `β` to `(P₀ᵀ M_S P₀)_{ab} / Area(T)`:
the P1 mass formula `Area_p/12 · Σ_{r,s} (1 + δ_rs) λ_a(p_r) λ_b(p_s)` with `Area_p = det β · Area(T)`.
-/
noncomputable def pieceMass (β : Fin 3 → Fin 3 → ℝ) (a b : Fin 3) : ℝ :=
  bdet β / 12 * ∑ r : Fin 3, ∑ s : Fin 3, (if r = s then 2 else 1) * β r a * β s b

/-- **The whole triangle gives `M₀`**: with the identity rows, `pieceMass = (1 + δ_ab)/12`, i.e.
`M₀ / Area(T)`. [formula: ic-transfer-mass-identity] -/
theorem pieceMass_one (a b : Fin 3) :
    pieceMass (fun r a => if r = a then 1 else 0) a b = (if a = b then 2 else 1) / 12 := by
  fin_cases a <;> fin_cases b <;> simp [pieceMass, bdet] <;> norm_num

/-- Replace row `r` of `β` by `x`. -/
def setRow (β : Fin 3 → Fin 3 → ℝ) (r : Fin 3) (x : Fin 3 → ℝ) : Fin 3 → Fin 3 → ℝ :=
  fun r' => if r' = r then x else β r'

/-- **`pieceMass` is a valuation**: splitting a signed triangle at any point `x` of its plane (rows
summing to one) into three signed triangles leaves the sum unchanged.
[formula: ic-transfer-mass-identity] -/
theorem pieceMass_split (β : Fin 3 → Fin 3 → ℝ) (hβ : ∀ r, β r 0 + β r 1 + β r 2 = 1)
    (x : Fin 3 → ℝ) (hx : x 0 + x 1 + x 2 = 1) (a b : Fin 3) :
    pieceMass β a b = pieceMass (setRow β 0 x) a b + pieceMass (setRow β 1 x) a b +
      pieceMass (setRow β 2 x) a b := by
  have e0 : β 0 2 = 1 - β 0 0 - β 0 1 := by linarith [hβ 0]
  have e1 : β 1 2 = 1 - β 1 0 - β 1 1 := by linarith [hβ 1]
  have e2 : β 2 2 = 1 - β 2 0 - β 2 1 := by linarith [hβ 2]
  have ex : x 2 = 1 - x 0 - x 1 := by linarith
  have hc : ∀ k : Fin 3, k = 0 ∨ k = 1 ∨ k = 2 := by decide
  rcases hc a with rfl | rfl | rfl <;> rcases hc b with rfl | rfl | rfl <;>
    simp only [pieceMass, bdet, setRow, Fin.sum_univ_three, Fin.isValue, Fin.reduceEq,
      ↓reduceIte, e0, e1, e2, ex] <;> ring

/-- **Splitting an edge**: a point `x` on the edge between corners `0` and `1` (`x = (1−t) β₀ + t β₁`)
splits the piece into two with the same total. [formula: ic-transfer-mass-identity] -/
theorem pieceMass_split_edge (β : Fin 3 → Fin 3 → ℝ) (hβ : ∀ r, β r 0 + β r 1 + β r 2 = 1)
    (t : ℝ) (a b : Fin 3) :
    pieceMass β a b = pieceMass (setRow β 0 (fun c => (1 - t) * β 0 c + t * β 1 c)) a b +
      pieceMass (setRow β 1 (fun c => (1 - t) * β 0 c + t * β 1 c)) a b := by
  have h := pieceMass_split β hβ (fun c => (1 - t) * β 0 c + t * β 1 c)
    (by linear_combination (1 - t) * hβ 0 + t * hβ 1) a b
  have z : bdet (setRow β 2 (fun c => (1 - t) * β 0 c + t * β 1 c)) = 0 := by
    simp only [bdet, setRow, Fin.isValue, Fin.reduceEq, ↓reduceIte]; ring
  have z' : pieceMass (setRow β 2 (fun c => (1 - t) * β 0 c + t * β 1 c)) a b = 0 := by
    rw [pieceMass, z, zero_div, zero_mul]
  rw [h, z', add_zero]

/-- The triangle with rows `p, q, r`. -/
def tri (p q r : Fin 3 → ℝ) : Fin 3 → Fin 3 → ℝ := ![p, q, r]

/-- **A fan does not depend on its apex**: for a closed polygon `p 0, …, p (N−1)` of points of the
plane (rows summing to one), `Σ_t pieceMass (x, p t, p (t+1))` is the same for every apex `x`
(the signed decomposition of the polygon into triangles from `x`).
[formula: ic-transfer-mass-identity] -/
theorem pieceMass_fan_apex {N : ℕ} (p : Fin (N + 1) → Fin 3 → ℝ)
    (hp : ∀ t, p t 0 + p t 1 + p t 2 = 1) (x y : Fin 3 → ℝ) (hx : x 0 + x 1 + x 2 = 1)
    (hy : y 0 + y 1 + y 2 = 1) (a b : Fin 3) :
    ∑ t : Fin (N + 1), pieceMass (tri x (p t) (p (t + 1))) a b =
      ∑ t : Fin (N + 1), pieceMass (tri y (p t) (p (t + 1))) a b := by
  have hc : ∀ k : Fin 3, k = 0 ∨ k = 1 ∨ k = 2 := by decide
  -- split each triangle of the `x`-fan at `y`
  have hsplit : ∀ t, pieceMass (tri x (p t) (p (t + 1))) a b =
      pieceMass (tri y (p t) (p (t + 1))) a b + pieceMass (tri x y (p (t + 1))) a b +
        pieceMass (tri x (p t) y) a b := by
    intro t
    have hβ : ∀ r, tri x (p t) (p (t + 1)) r 0 + tri x (p t) (p (t + 1)) r 1 +
        tri x (p t) (p (t + 1)) r 2 = 1 := by
      intro r
      rcases hc r with rfl | rfl | rfl <;>
        simp only [tri, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
          Matrix.head_cons, Matrix.tail_cons, hx, hp]
    have := pieceMass_split _ hβ y hy a b
    have e0 : setRow (tri x (p t) (p (t + 1))) 0 y = tri y (p t) (p (t + 1)) := by
      funext i
      rcases hc i with rfl | rfl | rfl <;>
        simp only [setRow, tri, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
          Matrix.head_cons, Matrix.tail_cons, Fin.reduceEq, ↓reduceIte]
    have e1 : setRow (tri x (p t) (p (t + 1))) 1 y = tri x y (p (t + 1)) := by
      funext i
      rcases hc i with rfl | rfl | rfl <;>
        simp only [setRow, tri, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
          Matrix.head_cons, Matrix.tail_cons, Fin.reduceEq, ↓reduceIte]
    have e2 : setRow (tri x (p t) (p (t + 1))) 2 y = tri x (p t) y := by
      funext i
      rcases hc i with rfl | rfl | rfl <;>
        simp only [setRow, tri, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
          Matrix.head_cons, Matrix.tail_cons, Fin.reduceEq, ↓reduceIte]
    rw [this, e0, e1, e2]
  -- the two side terms telescope around the closed polygon
  have hrot : ∑ t : Fin (N + 1), pieceMass (tri x y (p (t + 1))) a b =
      ∑ t : Fin (N + 1), pieceMass (tri x y (p t)) a b :=
    Fintype.sum_equiv (Equiv.addRight 1) _ _ (fun t => by rw [Equiv.coe_addRight])
  have hopp : ∀ q, pieceMass (tri x q y) a b = -pieceMass (tri x y q) a b := by
    intro q
    simp only [pieceMass, tri, bdet, Fin.sum_univ_three, Fin.isValue, Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.reduceEq,
      ↓reduceIte]
    ring
  simp only [hsplit, Finset.sum_add_distrib, hopp, Finset.sum_neg_distrib, hrot]
  ring

end IntegerCoords

end DdgProofs
