import Mathlib

/-!
# Piecewise-affine functions and the exactness of a Newton step that keeps its pieces

The residual maps of the cone problem are piecewise affine: the only nonlinearity is a per-vertex
function that is affine on each of finitely many pieces. A semismooth Newton step solves the
linearization on the current pieces. Hence:

* `IsPiecewiseAffine Q piece slope`: `Q` is affine on each piece, with generalized derivative
  `slope`;
* `newton_stageI_exact`, `newton_stageII_exact`: if a Newton step does not change the pieces,
  the new iterate is an **exact** zero of the residual (no invertibility is assumed).

There is no geometry in this file. (It was the abstract part of `DdgProofs.Newton`; the
declarations keep their names.)
-/

noncomputable section

namespace DDG

open Matrix

/-- `Q` is affine on each piece: whenever `x'` lies in the same piece as `x`,
`Q x' = Q x + slope x · (x' − x)`. The function `slope` is then a generalized (semismooth)
derivative of `Q`. -/
def IsPiecewiseAffine {κ : Type*} (Q : ℝ → ℝ) (piece : ℝ → κ) (slope : ℝ → ℝ) : Prop :=
  ∀ x x' : ℝ, piece x' = piece x → Q x' = Q x + slope x * (x' - x)

/-! ## Abstract exactness of a Newton step that keeps the pieces -/

section Abstract

variable {ι κ : Type*} [Fintype ι]

/-- **Stage II, abstract form.** Residual
`F(u, φ, μ) = (L u − μ + Ω,  Lᵀφ + w u,  μ − Q(φ + μ))` with `Q` piecewise affine per vertex.
If `(v, q, ν)` solves the Newton system `DF(u, φ, μ)(v, q, ν) = −F(u, φ, μ)`, where
`DF(v, q, ν) = (L v − ν,  Lᵀq + w v,  ν − slope(φ + μ)·(q + ν))`, and the pieces of `φ + μ` are the
same after the step, then `F(u + v, φ + q, μ + ν) = 0`. No invertibility is assumed. -/
theorem newton_stageII_exact (L : Matrix ι ι ℝ) (w Ω : ι → ℝ) (Q : ι → ℝ → ℝ)
    (piece : ι → ℝ → κ) (slope : ι → ℝ → ℝ)
    (hQ : ∀ i, IsPiecewiseAffine (Q i) (piece i) (slope i)) (u φ μ v q ν : ι → ℝ)
    (h1 : L *ᵥ v - ν = -(L *ᵥ u - μ + Ω))
    (h2 : Lᵀ *ᵥ q + (fun i => w i * v i) = -(Lᵀ *ᵥ φ + fun i => w i * u i))
    (h3 : ∀ i, ν i - slope i (φ i + μ i) * (q i + ν i) = -(μ i - Q i (φ i + μ i)))
    (hact : ∀ i, piece i ((φ i + q i) + (μ i + ν i)) = piece i (φ i + μ i)) :
    L *ᵥ (u + v) - (μ + ν) + Ω = 0 ∧
      Lᵀ *ᵥ (φ + q) + (fun i => w i * (u i + v i)) = 0 ∧
      ∀ i, (μ i + ν i) - Q i ((φ i + q i) + (μ i + ν i)) = 0 := by
  refine ⟨?_, ?_, fun i => ?_⟩
  · rw [Matrix.mulVec_add]
    ext i
    have := congrFun h1 i
    simp only [Pi.sub_apply, Pi.add_apply, Pi.neg_apply, Pi.zero_apply] at this ⊢
    linarith
  · rw [Matrix.mulVec_add]
    ext i
    have := congrFun h2 i
    simp only [Pi.add_apply, Pi.neg_apply, Pi.zero_apply] at this ⊢
    linarith
  · have hq := hQ i (φ i + μ i) ((φ i + q i) + (μ i + ν i)) (hact i)
    have := h3 i
    rw [hq]
    linarith

/-- **Stage I, abstract form.** Residual `F_γ(u, φ) = (L u + Ω − Q(φ),  Lᵀφ + w u)` with `Q`
piecewise affine per vertex. If `(v, q)` solves the Newton system
`DF_γ(u, φ)(v, q) = −F_γ(u, φ)`, where `DF_γ(v, q) = (L v − slope(φ)·q,  Lᵀq + w v)`, and the pieces
of `φ` are the same after the step, then `F_γ(u + v, φ + q) = 0`. -/
theorem newton_stageI_exact (L : Matrix ι ι ℝ) (w Ω : ι → ℝ) (Q : ι → ℝ → ℝ)
    (piece : ι → ℝ → κ) (slope : ι → ℝ → ℝ)
    (hQ : ∀ i, IsPiecewiseAffine (Q i) (piece i) (slope i)) (u φ v q : ι → ℝ)
    (h1 : ∀ i, (L *ᵥ v) i - slope i (φ i) * q i = -((L *ᵥ u) i + Ω i - Q i (φ i)))
    (h2 : Lᵀ *ᵥ q + (fun i => w i * v i) = -(Lᵀ *ᵥ φ + fun i => w i * u i))
    (hact : ∀ i, piece i (φ i + q i) = piece i (φ i)) :
    (∀ i, (L *ᵥ (u + v)) i + Ω i - Q i (φ i + q i) = 0) ∧
      Lᵀ *ᵥ (φ + q) + (fun i => w i * (u i + v i)) = 0 := by
  refine ⟨fun i => ?_, ?_⟩
  · have hq := hQ i (φ i) (φ i + q i) (hact i)
    have := h1 i
    rw [Matrix.mulVec_add, hq]
    simp only [Pi.add_apply]
    linarith
  · rw [Matrix.mulVec_add]
    ext i
    have := congrFun h2 i
    simp only [Pi.add_apply, Pi.neg_apply, Pi.zero_apply] at this ⊢
    linarith

end Abstract

end DDG
