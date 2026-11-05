import Mathlib

/-!
# Repulsive Surfaces: the line-search tests (SPEC §7.2, D15)

`rs-armijo-reference`: the reference's sufficient-decrease test
`F(f) − F(f − τx) ≥ c τ ‖x‖ cos∠(g, x)` (`ref: line_search.cpp:86,103`) is **the standard Armijo
test** `F(f − τx) ≤ F(f) − c′ τ ⟨g, x⟩` with `c′ = c/‖g‖`. Since `c′` depends on `‖g‖`, the
reference's test is not invariant under scaling `F` (it is stricter than the standard one with the
same `c` exactly when `‖g‖ < 1`, for a descent direction).
-/

noncomputable section

namespace DDG

namespace RS

variable {n : Type*} [Fintype n]

/-- The Euclidean norm `√(v·v)` on `n → ℝ`. -/
def enorm (v : n → ℝ) : ℝ := √(v ⬝ᵥ v)

/-- **The reference's Armijo test is the standard one with `c/‖g‖`** (SPEC §7.2): for `‖g‖ > 0`,
`‖x‖ > 0` and `cos∠(g, x) = ⟨g, x⟩/(‖g‖‖x‖)`,
`F₀ − F₁ ≥ c τ ‖x‖ cos∠(g, x) ⟺ F₁ ≤ F₀ − (c/‖g‖) τ ⟨g, x⟩`.

[formula: rs-armijo-reference] -/
theorem armijo_reference (F₀ F₁ c τ : ℝ) (g x : n → ℝ) (hg : 0 < enorm g) (hx : 0 < enorm x) :
    F₀ - F₁ ≥ c * τ * enorm x * ((g ⬝ᵥ x) / (enorm g * enorm x)) ↔
      F₁ ≤ F₀ - (c / enorm g) * τ * (g ⬝ᵥ x) := by
  have e : c * τ * enorm x * ((g ⬝ᵥ x) / (enorm g * enorm x)) = (c / enorm g) * τ * (g ⬝ᵥ x) := by
    field_simp
  rw [e]
  constructor <;> intro h <;> linarith

/-- **Which test is stricter** (SPEC §7.2): for a descent direction (`⟨g, x⟩ > 0`), `τ > 0` and
`c > 0`, the reference's required decrease `(c/‖g‖) τ ⟨g, x⟩` exceeds the standard one
`c τ ⟨g, x⟩` iff `‖g‖ < 1`.

[formula: rs-armijo-reference] -/
theorem armijo_reference_stricter {c τ : ℝ} (hc : 0 < c) (hτ : 0 < τ) (g x : n → ℝ)
    (hg : 0 < enorm g) (hgx : 0 < g ⬝ᵥ x) :
    c * τ * (g ⬝ᵥ x) < (c / enorm g) * τ * (g ⬝ᵥ x) ↔ enorm g < 1 := by
  have hp : 0 < c * τ * (g ⬝ᵥ x) := by positivity
  have e : (c / enorm g) * τ * (g ⬝ᵥ x) = (c * τ * (g ⬝ᵥ x)) / enorm g := by ring
  rw [e, lt_div_iff₀ hg]
  constructor
  · intro h; nlinarith
  · intro h; nlinarith

end RS

end DDG
