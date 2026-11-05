import DDG.SDF.Blend

/-!
# Points as Tori: applications (SPEC §8) and the lifted analytic shapes (SPEC §11)

* `pat-orientation-metric`: `M(r) = √(Σ_k (κ_k ⟨r, b_k⟩)²)` is a norm when all `κ_k ≠ 0` and
  `b_1, b_2, b_3` are linearly independent (triple product `≠ 0`).
* `pat-mc-interpolation`: the vertex `t = f₀/(f₀ − f₁)` is the unique root of the linear
  interpolant on an edge with a sign change, and lies on the edge.
* `pat-boolean`: the sign rules of `min`, `max`, `max(·, −·)`, the preservation of Lipschitz bounds,
  and the gradient taken from the active argument.
* `pat-revolution`, `pat-extrusion`, `pat-extrusion-gradient`: Quilez's revolution and extrusion of
  a 2D SDF, their case forms and gradients.
-/

noncomputable section

namespace DDG

namespace PAT

open Matrix Real

/-! ## The orientation metric (SPEC §8.4) -/

/-- Triangle inequality for `vnorm`. -/
theorem vnorm_add_le (a b : V3) : vnorm (a + b) ≤ vnorm a + vnorm b := by
  have h := dot_le_vnorm_mul a b
  have hab : vnorm (a + b) ^ 2 ≤ (vnorm a + vnorm b) ^ 2 := by
    rw [vnorm_sq]
    have e : (a + b) ⬝ᵥ (a + b) = a ⬝ᵥ a + 2 * (a ⬝ᵥ b) + b ⬝ᵥ b := by
      simp only [add_dotProduct, dotProduct_add, dotProduct_comm b a]; ring
    rw [e, ← vnorm_sq a, ← vnorm_sq b]; nlinarith
  exact le_of_sq_le_sq hab (add_nonneg (vnorm_nonneg a) (vnorm_nonneg b))

/-- The coordinates `(κ₁⟨r, b₁⟩, κ₂⟨r, b₂⟩, κ₃⟨r, b₃⟩)`. -/
def orientVec (k : Fin 3 → ℝ) (b : Fin 3 → V3) (r : V3) : V3 := fun j => k j * (b j ⬝ᵥ r)

/-- The anisotropic metric of Fig. 27: `M(r) = (Σ_j κ_j² ⟨r, b_j⟩²)^{1/2}` (SPEC §8.4 with
`κ = (κmin, κmax, ε)`, `b = (v̂min, v̂max, n_i)`). -/
def orientMetric (k : Fin 3 → ℝ) (b : Fin 3 → V3) (r : V3) : ℝ := vnorm (orientVec k b r)

theorem orientMetric_eq (k : Fin 3 → ℝ) (b : Fin 3 → V3) (r : V3) :
    orientMetric k b r =
      √((k 0) ^ 2 * (r ⬝ᵥ b 0) ^ 2 + (k 1) ^ 2 * (r ⬝ᵥ b 1) ^ 2 + (k 2) ^ 2 * (r ⬝ᵥ b 2) ^ 2) := by
  unfold orientMetric vnorm orientVec
  congr 1
  simp only [dotProduct, Fin.sum_univ_three]; ring

/-- Cramer's identity `det[b₀,b₁,b₂] r = (r·b₀)(b₁ × b₂) + (r·b₁)(b₂ × b₀) + (r·b₂)(b₀ × b₁)`. -/
theorem cramer3 (b0 b1 b2 r : V3) :
    (b0 ⬝ᵥ (b1 ⨯₃ b2)) • r =
      (r ⬝ᵥ b0) • (b1 ⨯₃ b2) + (r ⬝ᵥ b1) • (b2 ⨯₃ b0) + (r ⬝ᵥ b2) • (b0 ⨯₃ b1) := by
  ext i
  fin_cases i <;> simp [cross_apply, dotProduct, Fin.sum_univ_three] <;> ring

/-- **The orientation metric is a norm** when `κ_j ≠ 0` and `b₀·(b₁ × b₂) ≠ 0`: definite,
absolutely homogeneous and subadditive. (For the fitted directions `b₀ × b₁ ∥ n*` and
`n_i·n* = 1/A > 0`, so the triple product is nonzero when `v̂min ⊥ v̂max` are unit.)

[formula: pat-orientation-metric] -/
theorem orientMetric_norm (k : Fin 3 → ℝ) (b : Fin 3 → V3) :
    (∀ r, 0 ≤ orientMetric k b r) ∧
    (∀ (t : ℝ) r, orientMetric k b (t • r) = |t| * orientMetric k b r) ∧
    (∀ r s, orientMetric k b (r + s) ≤ orientMetric k b r + orientMetric k b s) ∧
    ((∀ j, k j ≠ 0) → b 0 ⬝ᵥ (b 1 ⨯₃ b 2) ≠ 0 → ∀ r, orientMetric k b r = 0 → r = 0) := by
  refine ⟨fun r => vnorm_nonneg _, fun t r => ?_, fun r s => ?_, fun hk hb r h => ?_⟩
  · unfold orientMetric
    rw [show orientVec k b (t • r) = t • orientVec k b r by
      ext j; simp [orientVec, dotProduct_smul]; ring, vnorm_smul]
  · unfold orientMetric
    rw [show orientVec k b (r + s) = orientVec k b r + orientVec k b s by
      ext j; simp [orientVec, dotProduct_add]; ring]
    exact vnorm_add_le _ _
  · have h0 : orientVec k b r = 0 := vnorm_eq_zero_iff.mp h
    have hd : ∀ j, b j ⬝ᵥ r = 0 := by
      intro j
      have := congrFun h0 j
      simp only [orientVec, Pi.zero_apply] at this
      exact (mul_eq_zero.mp this).resolve_left (hk j)
    have hc := cramer3 (b 0) (b 1) (b 2) r
    rw [dotProduct_comm r (b 0), dotProduct_comm r (b 1), dotProduct_comm r (b 2), hd 0, hd 1,
      hd 2] at hc
    simp only [zero_smul, add_zero] at hc
    exact (smul_eq_zero.mp hc).resolve_left hb

/-! ## Marching cubes (SPEC §8.6) -/

/-- **The marching-cubes vertex** on an edge with values `f₀, f₁` of opposite signs
(`f₀ f₁ < 0`): `t = f₀/(f₀ − f₁) ∈ (0, 1)` is the unique root of the linear interpolant
`s ↦ f₀ + s(f₁ − f₀)`, and `x = x₀ + t(x₁ − x₀)` is on the edge.

[formula: pat-mc-interpolation] -/
theorem mc_interpolation {f0 f1 : ℝ} (h : f0 * f1 < 0) :
    0 < f0 / (f0 - f1) ∧ f0 / (f0 - f1) < 1 ∧ f0 + f0 / (f0 - f1) * (f1 - f0) = 0 ∧
    ∀ s, f0 + s * (f1 - f0) = 0 → s = f0 / (f0 - f1) := by
  have hne : f0 - f1 ≠ 0 := by
    intro h'
    have : f0 = f1 := by linarith
    rw [this] at h; nlinarith [sq_nonneg f1]
  rcases lt_or_gt_of_ne (show f0 ≠ 0 by rintro rfl; simp at h) with h0 | h0
  · have h1 : 0 < f1 := by nlinarith
    have hd : f0 - f1 < 0 := by linarith
    refine ⟨div_pos_of_neg_of_neg h0 hd, by rw [div_lt_one_of_neg hd]; linarith, by field_simp; ring,
      fun s hs => by field_simp; linarith⟩
  · have h1 : f1 < 0 := by nlinarith
    have hd : 0 < f0 - f1 := by linarith
    refine ⟨div_pos h0 hd, by rw [div_lt_one hd]; linarith, by field_simp; ring,
      fun s hs => by field_simp; linarith⟩

/-! ## Booleans (SPEC §8.2) -/

/-- **CSG of SDFs** (SPEC §8.2): the union `min` is negative iff one argument is, the
intersection `max` iff both are, the difference `max(φ_A, −φ_B)` iff `φ_A < 0 < φ_B`; `min` and
`max` preserve `L`-Lipschitz bounds (so exact SDFs give 1-Lipschitz results, the bound used by the
marching-cubes band and the tracer); and where one argument is strictly smaller its gradient is
the gradient of the `min`.

[formula: pat-boolean] -/
theorem boolean_rules :
    (∀ a b : ℝ, min a b < 0 ↔ a < 0 ∨ b < 0) ∧
    (∀ a b : ℝ, max a b < 0 ↔ a < 0 ∧ b < 0) ∧
    (∀ a b : ℝ, max a (-b) < 0 ↔ a < 0 ∧ 0 < b) ∧
    (∀ (fA fB : V3 → ℝ) (L : ℝ), (∀ x y, |fA x - fA y| ≤ L * vnorm (x - y)) →
      (∀ x y, |fB x - fB y| ≤ L * vnorm (x - y)) →
      (∀ x y, |min (fA x) (fB x) - min (fA y) (fB y)| ≤ L * vnorm (x - y)) ∧
      (∀ x y, |max (fA x) (fB x) - max (fA y) (fB y)| ≤ L * vnorm (x - y))) := by
  refine ⟨fun a b => min_lt_iff, fun a b => max_lt_iff, fun a b => by
    rw [max_lt_iff, neg_lt_zero], fun fA fB L hA hB => ⟨fun x y => ?_, fun x y => ?_⟩⟩
  · exact (abs_min_sub_min_le_max _ _ _ _).trans (max_le (hA x y) (hB x y))
  · exact (abs_max_sub_max_le_max _ _ _ _).trans (max_le (hA x y) (hB x y))

/-- The gradient of `min(φ_A, φ_B)` where `φ_A < φ_B` is `∇φ_A`.

[formula: pat-boolean] -/
theorem boolean_gradient {fA fB : V3 → ℝ} {DA DB : V3 →L[ℝ] ℝ} {x : V3}
    (hA : HasFDerivAt fA DA x) (hB : HasFDerivAt fB DB x) (hlt : fA x < fB x) :
    HasFDerivAt (fun y => min (fA y) (fB y)) DA x := by
  have hc : ∀ᶠ y in nhds x, fA y < fB y :=
    (hA.continuousAt.prodMk hB.continuousAt).eventually
      (isOpen_lt continuous_fst continuous_snd |>.mem_nhds hlt)
  exact hA.congr_of_eventuallyEq (hc.mono fun y hy => min_eq_left hy.le)

/-! ## Revolution and extrusion of a 2D SDF (SPEC §11) -/

/-- Revolution about `y`: `Φ(x) = φ₂(|(x₀, x₂)| − o, x₁)`. -/
def revolve (φ2 : ℝ × ℝ → ℝ) (o : ℝ) (x : V3) : ℝ := φ2 (√(x 0 ^ 2 + x 2 ^ 2) - o, x 1)

/-- **Gradient of the revolution**: with `∇φ₂ = (a, b)` at `(ρ − o, x₁)` and `ρ = |(x₀, x₂)| > 0`,
`∇Φ = (a x₀/ρ, b, a x₂/ρ)`. The revolved field is invariant under rotations about `y` (it depends
on `x` through `ρ` and `x₁` only).

[formula: pat-revolution] -/
theorem revolve_gradient {φ2 : ℝ × ℝ → ℝ} {a b o : ℝ} {x : V3}
    (hρ : x 0 ^ 2 + x 2 ^ 2 ≠ 0)
    (hφ : HasFDerivAt φ2 (a • ContinuousLinearMap.fst ℝ ℝ ℝ + b • ContinuousLinearMap.snd ℝ ℝ ℝ)
      (√(x 0 ^ 2 + x 2 ^ 2) - o, x 1)) :
    HasFDerivAt (revolve φ2 o)
      (dualOf ![a * x 0 / √(x 0 ^ 2 + x 2 ^ 2), b, a * x 2 / √(x 0 ^ 2 + x 2 ^ 2)]) x := by
  have hc : ∀ i : Fin 3, HasFDerivAt (fun y : V3 => y i)
      (ContinuousLinearMap.proj (R := ℝ) (φ := fun _ : Fin 3 => ℝ) i) x :=
    fun i => hasFDerivAt_apply (𝕜 := ℝ) (F' := fun _ : Fin 3 => ℝ) i x
  have hs := (((hc 0).pow 2).add ((hc 2).pow 2)).sqrt hρ
  have hp := (hs.sub_const o).prodMk (hc 1)
  have := hφ.comp x hp
  refine (this.congr_fderiv ?_)
  have hpos : 0 < √(x 0 ^ 2 + x 2 ^ 2) := Real.sqrt_pos.mpr (lt_of_le_of_ne (by positivity) (Ne.symm hρ))
  ext h
  simp [dotProduct, Fin.sum_univ_three]
  field_simp
  ring

/-- Extrusion along `z`: `Φ = min(max(w₁, w₂), 0) + |max(w, 0)|`, `w = (φ₂(x₀, x₁), |x₂| − h)`. -/
def extrudeVal (w1 w2 : ℝ) : ℝ := min (max w1 w2) 0 + √(max w1 0 ^ 2 + max w2 0 ^ 2)

/-- **The extrusion formula by cases**: `|w|` if `w₁, w₂ ≥ 0`; `w₁` if `w₁ ≥ 0 ≥ w₂`; `w₂` if
`w₂ ≥ 0 ≥ w₁`; `max(w₁, w₂)` if both are `≤ 0`.

[formula: pat-extrusion] -/
theorem extrude_cases (w1 w2 : ℝ) :
    (0 ≤ w1 → 0 ≤ w2 → extrudeVal w1 w2 = √(w1 ^ 2 + w2 ^ 2)) ∧
    (0 ≤ w1 → w2 ≤ 0 → extrudeVal w1 w2 = w1) ∧
    (w1 ≤ 0 → 0 ≤ w2 → extrudeVal w1 w2 = w2) ∧
    (w1 ≤ 0 → w2 ≤ 0 → extrudeVal w1 w2 = max w1 w2) := by
  unfold extrudeVal
  refine ⟨fun h1 h2 => ?_, fun h1 h2 => ?_, fun h1 h2 => ?_, fun h1 h2 => ?_⟩
  · rw [max_eq_left h1, max_eq_left h2, min_eq_right (le_max_of_le_left h1), zero_add]
  · rw [max_eq_left h1, max_eq_right h2, min_eq_right (le_max_of_le_left h1), zero_add]
    simp [Real.sqrt_sq h1]
  · rw [max_eq_right h1, max_eq_left h2, min_eq_right (le_max_of_le_right h2), zero_add]
    simp [Real.sqrt_sq h2]
  · rw [max_eq_right h1, max_eq_right h2, min_eq_left (max_le h1 h2)]
    simp

/-- **Gradients of the extrusion** in the open regions, for `w₁ = f(x)` and `w₂ = |x₂| − h` with
gradients `∇w₁ = G₁`, `∇w₂ = G₂` (for `x₂ ≠ 0`, `G₂ = sgn(x₂) e_z`): `w/|w|·(G₁, G₂)` outside
both slabs, `G₁` or `G₂` where only one is positive, and the gradient of the larger one where both
are negative.

[formula: pat-extrusion-gradient] -/
theorem extrude_gradient {f1 f2 : V3 → ℝ} {G1 G2 : V3} {x : V3}
    (h1 : HasFDerivAt f1 (dualOf G1) x) (h2 : HasFDerivAt f2 (dualOf G2) x) :
    (0 < f1 x → 0 < f2 x → HasFDerivAt (fun y => extrudeVal (f1 y) (f2 y))
      (dualOf ((√(f1 x ^ 2 + f2 x ^ 2))⁻¹ • (f1 x • G1 + f2 x • G2))) x) ∧
    (0 < f1 x → f2 x < 0 → HasFDerivAt (fun y => extrudeVal (f1 y) (f2 y)) (dualOf G1) x) ∧
    (f1 x < 0 → 0 < f2 x → HasFDerivAt (fun y => extrudeVal (f1 y) (f2 y)) (dualOf G2) x) ∧
    (f2 x < f1 x → f1 x < 0 → HasFDerivAt (fun y => extrudeVal (f1 y) (f2 y)) (dualOf G1) x) ∧
    (f1 x < f2 x → f2 x < 0 → HasFDerivAt (fun y => extrudeVal (f1 y) (f2 y)) (dualOf G2) x) := by
  have c1 := h1.continuousAt
  have c2 := h2.continuousAt
  have ev : ∀ {P : ℝ → Prop}, IsOpen {t | P t} → ∀ {F : V3 → ℝ}, ContinuousAt F x → P (F x) →
      ∀ᶠ y in nhds x, P (F y) := fun hP _ hF hx => hF.eventually (hP.mem_nhds hx)
  refine ⟨fun a b => ?_, fun a b => ?_, fun a b => ?_, fun a b => ?_, fun a b => ?_⟩
  · have hne : f1 x ^ 2 + f2 x ^ 2 ≠ 0 := by positivity
    have hs := ((h1.pow 2).add (h2.pow 2)).sqrt hne
    have e1 := ev (P := fun t => 0 < t) isOpen_Ioi c1 a
    have e2 := ev (P := fun t => 0 < t) isOpen_Ioi c2 b
    refine (hs.congr_fderiv ?_).congr_of_eventuallyEq ((e1.and e2).mono fun y hy =>
      (extrude_cases _ _).1 hy.1.le hy.2.le)
    ext v; simp [add_dotProduct, smul_dotProduct]; field_simp; try ring
  · have e1 := ev (P := fun t => 0 < t) isOpen_Ioi c1 a
    have e2 := ev (P := fun t => t < 0) isOpen_Iio c2 b
    exact h1.congr_of_eventuallyEq ((e1.and e2).mono fun y hy =>
      (extrude_cases _ _).2.1 hy.1.le hy.2.le)
  · have e1 := ev (P := fun t => t < 0) isOpen_Iio c1 a
    have e2 := ev (P := fun t => 0 < t) isOpen_Ioi c2 b
    exact h2.congr_of_eventuallyEq ((e1.and e2).mono fun y hy =>
      (extrude_cases _ _).2.2.1 hy.1.le hy.2.le)
  · have e1 := ev (P := fun t => t < 0) isOpen_Iio c1 b
    have e12 : ∀ᶠ y in nhds x, f2 y < f1 y :=
      (c2.prodMk c1).eventually (isOpen_lt continuous_fst continuous_snd |>.mem_nhds a)
    refine h1.congr_of_eventuallyEq ((e1.and e12).mono fun y hy => ?_)
    change extrudeVal (f1 y) (f2 y) = f1 y
    rw [(extrude_cases _ _).2.2.2 hy.1.le (hy.2.trans hy.1).le, max_eq_left hy.2.le]
  · have e2 := ev (P := fun t => t < 0) isOpen_Iio c2 b
    have e12 : ∀ᶠ y in nhds x, f1 y < f2 y :=
      (c1.prodMk c2).eventually (isOpen_lt continuous_fst continuous_snd |>.mem_nhds a)
    refine h2.congr_of_eventuallyEq ((e2.and e12).mono fun y hy => ?_)
    change extrudeVal (f1 y) (f2 y) = f2 y
    rw [(extrude_cases _ _).2.2.2 (hy.2.trans hy.1).le hy.1.le, max_eq_right hy.2.le]

/-! ## The sphere-tracing step (SPEC §8.5, D17) -/

/-- A step of `|φ(x) − ℓ| / L` never reaches the level set `φ = ℓ` when `φ` is `L`-Lipschitz: every
point `y` strictly closer to `x` than `|φ(x) − ℓ| / L` has `φ(y) ≠ ℓ`. (The tracer steps by
`safety · |φ − ℓ| / max(‖∇φ‖, g_min)`, `safety ≤ 1`; where `max(‖∇φ‖, g_min)` bounds the local
Lipschitz constant, this is the lemma with `L` that bound.)

[formula: pat-render-step] -/
theorem sphere_trace_step_safe {φ : V3 → ℝ} {L ℓ : ℝ} (hL : 0 < L)
    (hφ : ∀ x y, |φ x - φ y| ≤ L * vnorm (x - y)) (x y : V3)
    (hy : vnorm (y - x) < |φ x - ℓ| / L) : φ y ≠ ℓ := by
  intro h
  have h1 := hφ x y
  rw [DDG.vnorm_sub_comm x y, h] at h1
  have h2 : L * vnorm (y - x) < |φ x - ℓ| := by
    rwa [lt_div_iff₀ hL, mul_comm] at hy
  linarith

end PAT

end DDG
