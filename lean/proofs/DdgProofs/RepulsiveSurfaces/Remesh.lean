import DDG.Repulsive.Surfaces.Basic

/-!
# Repulsive Surfaces: the formulas of the remesher (paper Sec. 6; SPEC §10, D16)

* `rs-circumcenter` (`ref: remeshing.cpp:73–87`): with squared side lengths `A = |P₁ − P₂|²`,
  `B = |P₂ − P₀|²`, `C = |P₀ − P₁|²`, the barycentric weights
  `A(B + C − A) : B(C + A − B) : C(A + B − C)` (sum `W = 16·area²`) give a point equidistant from
  the three corners (`circumcenter_equidistant`), in the triangle's plane (an affine combination).
* `rs-smoothing-weights`: `Σ_{S∋v} a_S = 3 a_v`; hence the reference's displacement
  `½ Σ a_S (c_S − f_v) / (3 a_v)` is the paper's `ρ Σ a_S (c_S − f_v) / Σ a_S` with `ρ = ½`.
* `rs-delaunay-test`: for opposite angles `θ₁, θ₂ ∈ (0, π)`, `θ₁ + θ₂ > π` iff
  `cot θ₁ + cot θ₂ < 0` (the cotan weight of the edge is negative), and `= π` iff the sum is `0`.
-/

noncomputable section

namespace DDG

namespace RS

open Matrix Real Finset
open RC (V3)

/-! ## Circumcenter -/

/-- The barycentric circumcenter weights of `ref: remeshing.cpp`. -/
def ccWeights (P₀ P₁ P₂ : V3) : Fin 3 → ℝ :=
  let A := (P₁ - P₂) ⬝ᵥ (P₁ - P₂)
  let B := (P₂ - P₀) ⬝ᵥ (P₂ - P₀)
  let C := (P₀ - P₁) ⬝ᵥ (P₀ - P₁)
  ![A * (B + C - A), B * (C + A - B), C * (A + B - C)]

/-- The circumcenter `Σ w_i P_i / Σ w_i`.

[formula: rs-circumcenter] -/
def circumcenter (P₀ P₁ P₂ : V3) : V3 :=
  (∑ i, ccWeights P₀ P₁ P₂ i)⁻¹ •
    (ccWeights P₀ P₁ P₂ 0 • P₀ + ccWeights P₀ P₁ P₂ 1 • P₁ + ccWeights P₀ P₁ P₂ 2 • P₂)

/-- `|P − Q|² = |P − R|²` from `2 (P − R)·(R − Q)`-type scalar identities: with `X = c − P₀`,
equidistance from `P₀` and `P₀ + u` is `2 X·u = u·u`. -/
theorem equidist_iff (X u : V3) : X ⬝ᵥ X = (X - u) ⬝ᵥ (X - u) ↔ 2 * (X ⬝ᵥ u) = u ⬝ᵥ u := by
  rw [sub_dotProduct, dotProduct_sub, dotProduct_sub, dotProduct_comm u X]
  constructor <;> intro h <;> linarith

/-- **The barycentric weights give the circumcenter** (SPEC §10): for a nondegenerate triangle
(`W = Σ w_i ≠ 0`), `c` is equidistant from `P₀`, `P₁`, `P₂`, and the weights sum to `W`.

[formula: rs-circumcenter] -/
theorem circumcenter_equidistant (P₀ P₁ P₂ : V3) (hW : ∑ i, ccWeights P₀ P₁ P₂ i ≠ 0) :
    let c := circumcenter P₀ P₁ P₂
    (c - P₀) ⬝ᵥ (c - P₀) = (c - P₁) ⬝ᵥ (c - P₁) ∧ (c - P₀) ⬝ᵥ (c - P₀) = (c - P₂) ⬝ᵥ (c - P₂) := by
  intro c
  set u := P₁ - P₀ with hu
  set v := P₂ - P₀ with hv
  set W := ∑ i, ccWeights P₀ P₁ P₂ i with hWdef
  -- the weights in terms of `u·u`, `v·v`, `u·v`
  have hA : (P₁ - P₂) ⬝ᵥ (P₁ - P₂) = u ⬝ᵥ u + v ⬝ᵥ v - 2 * (u ⬝ᵥ v) := by
    rw [hu, hv]; simp only [dotProduct, Fin.sum_univ_three, Pi.sub_apply]; ring
  have hB : (P₂ - P₀) ⬝ᵥ (P₂ - P₀) = v ⬝ᵥ v := rfl
  have hC : (P₀ - P₁) ⬝ᵥ (P₀ - P₁) = u ⬝ᵥ u := by
    rw [hu]; simp only [dotProduct, Fin.sum_univ_three, Pi.sub_apply]; ring
  set w0 := ccWeights P₀ P₁ P₂ 0
  set w1 := ccWeights P₀ P₁ P₂ 1
  set w2 := ccWeights P₀ P₁ P₂ 2
  -- `c − P₀ = W⁻¹ (w₁ u + w₂ v)`
  have hc : c - P₀ = W⁻¹ • (w1 • u + w2 • v) := by
    have hWc : W • c = w0 • P₀ + w1 • P₁ + w2 • P₂ := by
      simp only [c, circumcenter, ← hWdef, smul_smul, mul_inv_cancel₀ hW, one_smul]
      rfl
    have hsum : W = w0 + w1 + w2 := by rw [hWdef, Fin.sum_univ_three]
    have h2 : W • (c - P₀) = w1 • u + w2 • v := by
      rw [smul_sub, hWc]
      conv_lhs => rw [hsum]
      rw [hu, hv]; simp only [add_smul, smul_sub]; abel
    rw [← h2, smul_smul, inv_mul_cancel₀ hW, one_smul]
  have hw1 : w1 = v ⬝ᵥ v * (u ⬝ᵥ u + (u ⬝ᵥ u + v ⬝ᵥ v - 2 * (u ⬝ᵥ v)) - v ⬝ᵥ v) := by
    simp only [w1, ccWeights, hA, hB, hC]; rfl
  have hw2 : w2 = u ⬝ᵥ u * ((u ⬝ᵥ u + v ⬝ᵥ v - 2 * (u ⬝ᵥ v)) + v ⬝ᵥ v - u ⬝ᵥ u) := by
    simp only [w2, ccWeights, hA, hB, hC]; rfl
  have hWv : W = (u ⬝ᵥ u + v ⬝ᵥ v - 2 * (u ⬝ᵥ v)) * (v ⬝ᵥ v + u ⬝ᵥ u -
      (u ⬝ᵥ u + v ⬝ᵥ v - 2 * (u ⬝ᵥ v))) + w1 + w2 := by
    simp only [hWdef, Fin.sum_univ_three, ccWeights, hA, hB, hC, w1, w2]; rfl
  have e1 : c - P₁ = (c - P₀) - u := by rw [hu]; abel
  have e2 : c - P₂ = (c - P₀) - v := by rw [hv]; abel
  rw [e1, e2, equidist_iff, equidist_iff, hc, smul_dotProduct, smul_dotProduct, add_dotProduct,
    add_dotProduct, smul_dotProduct, smul_dotProduct, smul_dotProduct, smul_dotProduct]
  simp only [smul_eq_mul]
  rw [dotProduct_comm v u]
  constructor
  · field_simp
    rw [hWv, hw1, hw2]; ring
  · field_simp
    rw [hWv, hw1, hw2]; ring

/-! ## Smoothing weights -/

section smoothing

variable {V F : Type*} [Fintype V] [Fintype F] [DecidableEq V]

omit [Fintype V] in
/-- **`Σ_{S∋v} a_S = 3 a_v`**, and the reference's `½ Σ_{S∋v} a_S (c_S − f_v)/(3a_v)` is the
paper's `ρ Σ_{S∋v} a_S (c_S − f_v)/Σ_{S∋v} a_S` with `ρ = ½` (SPEC §10).

[formula: rs-smoothing-weights] -/
theorem smoothing_weights (f : V → V3) (fv : F → Fin 3 → V) (v : V) (d : F → V3) :
    (∑ S, ∑ a, if fv S a = v then area f fv S else 0) = 3 * dualArea f fv v ∧
      (1 / 2 : ℝ) • ((3 * dualArea f fv v)⁻¹ • ∑ S, ∑ a, if fv S a = v then area f fv S • d S
        else 0) =
      (1 / 2 : ℝ) • ((∑ S, ∑ a, if fv S a = v then area f fv S else 0)⁻¹ •
        ∑ S, ∑ a, if fv S a = v then area f fv S • d S else 0) := by
  have h : (∑ S, ∑ a, if fv S a = v then area f fv S else 0) = 3 * dualArea f fv v := by
    rw [dualArea]; ring
  exact ⟨h, by rw [h]⟩

end smoothing

/-! ## Delaunay test -/

/-- **The Delaunay flip test** (SPEC §10): for angles `θ₁, θ₂ ∈ (0, π)` opposite an interior
edge, `θ₁ + θ₂ > π ⟺ cot θ₁ + cot θ₂ < 0` and `θ₁ + θ₂ = π ⟺ cot θ₁ + cot θ₂ = 0`; the edge's
cotan weight `½(cot θ₁ + cot θ₂)` is negative exactly when the edge is not Delaunay.

[formula: rs-delaunay-test] -/
theorem delaunay_test {θ₁ θ₂ : ℝ} (h1 : 0 < θ₁) (h1' : θ₁ < π) (h2 : 0 < θ₂) (h2' : θ₂ < π) :
    (π < θ₁ + θ₂ ↔ cos θ₁ / sin θ₁ + cos θ₂ / sin θ₂ < 0) ∧
      (θ₁ + θ₂ = π ↔ cos θ₁ / sin θ₁ + cos θ₂ / sin θ₂ = 0) := by
  have s1 : 0 < sin θ₁ := sin_pos_of_pos_of_lt_pi h1 h1'
  have s2 : 0 < sin θ₂ := sin_pos_of_pos_of_lt_pi h2 h2'
  have e : cos θ₁ / sin θ₁ + cos θ₂ / sin θ₂ = sin (θ₁ + θ₂) / (sin θ₁ * sin θ₂) := by
    rw [sin_add]; field_simp; ring
  have hp : 0 < sin θ₁ * sin θ₂ := mul_pos s1 s2
  rw [e]
  -- the sign of `sin` on `(0, 2π)`
  have hsum : 0 < θ₁ + θ₂ ∧ θ₁ + θ₂ < 2 * π := ⟨by linarith, by linarith⟩
  have neg_of : π < θ₁ + θ₂ → sin (θ₁ + θ₂) < 0 := fun h => by
    have := sin_pos_of_pos_of_lt_pi (x := θ₁ + θ₂ - π) (by linarith) (by linarith)
    rw [sin_sub_pi] at this; linarith
  have pos_of : θ₁ + θ₂ < π → 0 < sin (θ₁ + θ₂) := fun h => sin_pos_of_pos_of_lt_pi hsum.1 h
  constructor
  · rw [div_neg_iff]
    constructor
    · intro h; exact Or.inr ⟨neg_of h, hp⟩
    · rintro (⟨-, hb⟩ | ⟨h, -⟩)
      · linarith
      · by_contra hc
        rcases (not_lt.mp hc).lt_or_eq with hc | hc
        · linarith [pos_of hc]
        · rw [hc, sin_pi] at h; exact lt_irrefl _ h
  · rw [div_eq_zero_iff]
    constructor
    · intro h; left; rw [h, sin_pi]
    · rintro (h | h)
      · rcases lt_trichotomy (θ₁ + θ₂) π with hc | hc | hc
        · linarith [pos_of hc]
        · exact hc
        · linarith [neg_of hc]
      · linarith

end RS

end DDG
