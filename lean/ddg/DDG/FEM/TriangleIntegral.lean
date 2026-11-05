import DDG.Geometry.Plane

/-!
# Integrals over a triangle in the plane

The infrastructure behind the finite element formulas of the library: Lebesgue integrals over a
triangle of the plane `ℝ × ℝ` (Lebesgue measure `volume = volume.prod volume`). Points of the
plane are pairs; lengths are the Euclidean `dist2` and areas the signed double area `sda` of
`DDG/Geometry/Plane.lean` (the metric of the type `ℝ × ℝ` is the sup metric and is never used).

* `integral_pow_mul_sub_pow` — the Beta integral `∫₀ᴸ tᵐ (L − t)ⁿ dt = Lᵐ⁺ⁿ⁺¹ m! n!/(m + n + 1)!`;
* `stdTri` — the standard triangle `{(s, t) | s, t ≥ 0, s + t ≤ 1}`; `integral_stdTri` — Fubini:
  `∫_{stdTri} g = ∫₀¹ ∫₀^{1−s} g(s, t) dt ds` for continuous `g`;
* `triMap p` — the affine map `(s, t) ↦ p₀ + s (p₁ − p₀) + t (p₂ − p₀)`, `triangle p` its image of
  the standard triangle, equal to the convex hull of the three corners (`triangle_eq_convexHull`),
  so it does not depend on their order (`triangle_perm`);
* `integral_triangle` — the change of variables `∫_T g = |D| ∫_{stdTri} g ∘ triMap`,
  `D = sda p₀ p₁ p₂` the Jacobian determinant (`det_triLin`), for every function `g`;
* `volume_real_triangle` — the area of the triangle is `|D|/2`, which is the Heron area of the
  Euclidean side lengths (`volume_real_triangle_eq_triArea`);
* `bary p i` — the barycentric (area) coordinates; on the triangle they are `(1 − s − t, s, t)`
  (`bary_triMap`), they sum to one and reproduce the point (`sum_bary`, `sum_bary_smul`);
* `integral_bary_slice` — the slice formula: for `(c, j, k)` a cyclic order of the corners and a
  continuous `G : ℝ → E`, `∫_T λ_c^{e_c} λ_j^q λ_k^r G(λ_c) = |D| q! r!/(q + r + 1)! ·
  ∫₀¹ s^{e_c} (1 − s)^{q+r+1} G(s) ds`;
* `integral_bary_monomial` — the area-coordinate rule
  `∫_T λ₀ᵃ λ₁ᵇ λ₂ᶜ = 2|T| a! b! c!/(a + b + c + 2)!`.

The triangle is any nondegenerate triangle of the plane (`sda p₀ p₁ p₂ ≠ 0`, either orientation).
An intrinsic triangle (three edge lengths) is placed by the isometric layout of `layout_isometric`;
every result below is a function of the edge lengths only, through `triArea` and `cornerCot`.
-/

noncomputable section

namespace DDG

namespace TriangleIntegral

open MeasureTheory Set Real
open scoped Nat

/-! ## The Beta integral -/

/-- **The Beta integral** with natural exponents, on `[0, L]`:
`∫₀ᴸ tᵐ (L − t)ⁿ dt = Lᵐ⁺ⁿ⁺¹ · m! n!/(m + n + 1)!`. By induction on `n`, from
`(L − t)ⁿ⁺¹ = L (L − t)ⁿ − t (L − t)ⁿ`. -/
theorem integral_pow_mul_sub_pow (L : ℝ) (m n : ℕ) :
    ∫ t in (0 : ℝ)..L, t ^ m * (L - t) ^ n =
      L ^ (m + n + 1) * ((m ! * n ! : ℕ) / ((m + n + 1)! : ℕ) : ℝ) := by
  induction n generalizing m with
  | zero =>
    simp only [pow_zero, mul_one, integral_pow, add_zero, Nat.factorial_zero]
    rw [Nat.factorial_succ]
    push_cast
    have : (m ! : ℝ) ≠ 0 := by positivity
    field_simp
    ring
  | succ n ih =>
    have e : ∀ t : ℝ, t ^ m * (L - t) ^ (n + 1) =
        L * (t ^ m * (L - t) ^ n) - t ^ (m + 1) * (L - t) ^ n := fun t => by ring
    simp_rw [e]
    rw [intervalIntegral.integral_sub, intervalIntegral.integral_const_mul, ih m, ih (m + 1)]
    · rw [show m + (n + 1) + 1 = m + n + 1 + 1 by ring, show m + 1 + n + 1 = m + n + 1 + 1 by ring,
        Nat.factorial_succ (m + n + 1), Nat.factorial_succ m, Nat.factorial_succ n]
      push_cast
      have h1 : ((m + n + 1)! : ℝ) ≠ 0 := by positivity
      have h2 : ((m + n + 1 : ℝ) + 1) ≠ 0 := by positivity
      field_simp
      ring
    · exact (by fun_prop : Continuous fun t : ℝ => L * (t ^ m * (L - t) ^ n)).intervalIntegrable _ _
    · exact (by fun_prop : Continuous fun t : ℝ => t ^ (m + 1) * (L - t) ^ n).intervalIntegrable _ _

/-! ## The standard triangle and Fubini -/

/-- The standard triangle `{(s, t) | 0 ≤ s, 0 ≤ t, s + t ≤ 1}`. -/
def stdTri : Set (ℝ × ℝ) := {q | 0 ≤ q.1 ∧ 0 ≤ q.2 ∧ q.1 + q.2 ≤ 1}

theorem isClosed_stdTri : IsClosed stdTri := by
  unfold stdTri
  refine (isClosed_le continuous_const continuous_fst).inter
    ((isClosed_le continuous_const continuous_snd).inter
      (isClosed_le (continuous_fst.add continuous_snd) continuous_const))

theorem measurableSet_stdTri : MeasurableSet stdTri := isClosed_stdTri.measurableSet

theorem isCompact_stdTri : IsCompact stdTri := by
  refine ((isCompact_Icc (a := (0 : ℝ)) (b := 1)).prod
    (isCompact_Icc (a := (0 : ℝ)) (b := 1))).of_isClosed_subset isClosed_stdTri ?_
  rintro ⟨s, t⟩ ⟨hs, ht, hst⟩
  exact ⟨⟨hs, by linarith⟩, ⟨ht, by linarith⟩⟩

/-- The slices of the standard triangle: `(s, t) ∈ stdTri ↔ s ∈ [0, 1] ∧ t ∈ [0, 1 − s]`. -/
theorem mem_stdTri_iff (s t : ℝ) : (s, t) ∈ stdTri ↔ s ∈ Icc (0 : ℝ) 1 ∧ t ∈ Icc 0 (1 - s) := by
  simp only [stdTri, mem_Icc]
  change (0 ≤ s ∧ 0 ≤ t ∧ s + t ≤ 1) ↔ _
  constructor
  · rintro ⟨hs, ht, hst⟩; exact ⟨⟨hs, by linarith⟩, ht, by linarith⟩
  · rintro ⟨⟨hs, _⟩, ht, hst⟩; exact ⟨hs, ht, by linarith⟩

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- **Fubini on the standard triangle**: for `g` continuous on the triangle,
`∫_{stdTri} g = ∫₀¹ (∫₀^{1−s} g(s, t) dt) ds`. -/
theorem integral_stdTri (g : ℝ × ℝ → E) (hg : ContinuousOn g stdTri) :
    ∫ q in stdTri, g q = ∫ s in (0 : ℝ)..1, ∫ t in (0 : ℝ)..(1 - s), g (s, t) := by
  have hint : Integrable (stdTri.indicator g) (volume.prod volume) := by
    rw [← Measure.volume_eq_prod]
    exact (hg.integrableOn_compact isCompact_stdTri).integrable_indicator measurableSet_stdTri
  rw [← integral_indicator measurableSet_stdTri, Measure.volume_eq_prod, integral_prod _ hint]
  have inner : ∀ s : ℝ, ∫ t, stdTri.indicator g (s, t) =
      (Icc (0 : ℝ) 1).indicator (fun s => ∫ t in (0 : ℝ)..(1 - s), g (s, t)) s := by
    intro s
    by_cases hs : s ∈ Icc (0 : ℝ) 1
    · rw [indicator_of_mem hs]
      have : (fun t => stdTri.indicator g (s, t)) =
          (Icc 0 (1 - s)).indicator (fun t => g (s, t)) := by
        funext t
        by_cases ht : t ∈ Icc 0 (1 - s)
        · rw [indicator_of_mem ((mem_stdTri_iff s t).2 ⟨hs, ht⟩), indicator_of_mem ht]
        · rw [indicator_of_notMem (fun h => ht ((mem_stdTri_iff s t).1 h).2),
            indicator_of_notMem ht]
      rw [this, integral_indicator measurableSet_Icc, integral_Icc_eq_integral_Ioc,
        ← intervalIntegral.integral_of_le (by linarith [hs.2])]
    · rw [indicator_of_notMem hs]
      have : (fun t => stdTri.indicator g (s, t)) = fun _ => 0 := by
        funext t
        exact indicator_of_notMem (fun h => hs ((mem_stdTri_iff s t).1 h).1) _
      rw [this, integral_zero]
  simp_rw [inner]
  rw [integral_indicator measurableSet_Icc, integral_Icc_eq_integral_Ioc,
    ← intervalIntegral.integral_of_le zero_le_one]

/-! ## The affine map onto a triangle and the change of variables -/

/-- The signed double area `D = (p₁ − p₀) × (p₂ − p₀)` of the triangle `p₀ p₁ p₂`. -/
abbrev D (p : Fin 3 → ℝ × ℝ) : ℝ := sda (p 0) (p 1) (p 2)

/-- The affine map of the standard triangle onto the triangle `p₀ p₁ p₂`:
`(s, t) ↦ p₀ + s (p₁ − p₀) + t (p₂ − p₀)`. -/
def triMap (p : Fin 3 → ℝ × ℝ) (q : ℝ × ℝ) : ℝ × ℝ := p 0 + q.1 • (p 1 - p 0) + q.2 • (p 2 - p 0)

/-- The linear part of `triMap`. -/
def triLin (p : Fin 3 → ℝ × ℝ) : ℝ × ℝ →L[ℝ] ℝ × ℝ :=
  (ContinuousLinearMap.fst ℝ ℝ ℝ).smulRight (p 1 - p 0) +
    (ContinuousLinearMap.snd ℝ ℝ ℝ).smulRight (p 2 - p 0)

/-- The triangle with corners `p₀, p₁, p₂`: the image of the standard triangle. -/
def triangle (p : Fin 3 → ℝ × ℝ) : Set (ℝ × ℝ) := triMap p '' stdTri

theorem triMap_eq (p : Fin 3 → ℝ × ℝ) (q : ℝ × ℝ) : triMap p q = p 0 + triLin p q := by
  simp only [triMap, triLin, add_apply, ContinuousLinearMap.smulRight_apply,
    ContinuousLinearMap.coe_fst', ContinuousLinearMap.coe_snd', add_assoc]

theorem hasFDerivAt_triMap (p : Fin 3 → ℝ × ℝ) (q : ℝ × ℝ) :
    HasFDerivAt (triMap p) (triLin p) q := by
  have h : triMap p = fun q => p 0 + triLin p q := funext (triMap_eq p)
  rw [h]
  exact (triLin p).hasFDerivAt.const_add (p 0)

/-- **The Jacobian determinant of the triangle map is the signed double area.** -/
theorem det_triLin (p : Fin 3 → ℝ × ℝ) : (triLin p).det = D p := by
  rw [ContinuousLinearMap.det, ← LinearMap.det_toMatrix (Module.Basis.finTwoProd ℝ),
    Matrix.det_fin_two]
  simp only [LinearMap.toMatrix_apply, Module.Basis.finTwoProd_zero, Module.Basis.finTwoProd_one,
    Module.Basis.coe_finTwoProd_repr, triLin, ContinuousLinearMap.toLinearMap_add,
    ContinuousLinearMap.coe_smulRight, LinearMap.add_apply, LinearMap.smulRight_apply,
    ContinuousLinearMap.coe_coe, ContinuousLinearMap.coe_fst', ContinuousLinearMap.coe_snd',
    one_smul, zero_smul, add_zero, zero_add, Matrix.cons_val_zero, Matrix.cons_val_one,
    Prod.fst_sub, Prod.snd_sub, sda]
  ring

theorem triMap_injective {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) : Function.Injective (triMap p) := by
  rintro ⟨s, t⟩ ⟨s', t'⟩ h
  have h1 := congrArg Prod.fst h
  have h2 := congrArg Prod.snd h
  simp only [triMap, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, Prod.fst_sub,
    Prod.snd_sub, smul_eq_mul] at h1 h2
  have ha : (s - s') * D p = 0 := by
    simp only [D, sda]
    linear_combination ((p 2).2 - (p 0).2) * h1 - ((p 2).1 - (p 0).1) * h2
  have hb : (t - t') * D p = 0 := by
    simp only [D, sda]
    linear_combination ((p 1).1 - (p 0).1) * h2 - ((p 1).2 - (p 0).2) * h1
  have e1 := (mul_eq_zero.1 ha).resolve_right hD
  have e2 := (mul_eq_zero.1 hb).resolve_right hD
  ext <;> simp only <;> linarith

/-- **Change of variables to the standard triangle.** For a nondegenerate triangle and every
function `g` (no integrability needed: both sides are `0` when `g` is not integrable),
`∫_T g = |D| ∫_{stdTri} g ∘ triMap`. -/
theorem integral_triangle {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) (g : ℝ × ℝ → E) :
    ∫ x in triangle p, g x = |D p| • ∫ q in stdTri, g (triMap p q) := by
  rw [triangle, integral_image_eq_integral_abs_det_fderiv_smul volume measurableSet_stdTri
    (fun q _ => (hasFDerivAt_triMap p q).hasFDerivWithinAt) (triMap_injective hD).injOn g,
    det_triLin, integral_smul]

/-- The standard triangle has area `½`. -/
theorem volume_real_stdTri : volume.real stdTri = 1 / 2 := by
  have h := integral_stdTri (fun _ => (1 : ℝ)) continuousOn_const
  rw [setIntegral_const, smul_eq_mul, mul_one] at h
  rw [h]
  simp only [intervalIntegral.integral_const, sub_zero, smul_eq_mul, mul_one]
  rw [intervalIntegral.integral_sub intervalIntegrable_const
    (continuous_id'.intervalIntegrable (μ := volume) 0 1)]
  simp only [intervalIntegral.integral_const, integral_id]
  norm_num

/-- **The area of a triangle is half the absolute signed double area.** -/
theorem volume_real_triangle {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) :
    volume.real (triangle p) = |D p| / 2 := by
  have h := integral_triangle hD (fun _ => (1 : ℝ))
  rw [setIntegral_const, setIntegral_const, volume_real_stdTri] at h
  simp only [smul_eq_mul, mul_one] at h
  rw [h]
  ring

/-- **The area of a triangle is the Heron area of its Euclidean side lengths**
`ℓ₁₂ = |p₁p₂|`, `ℓ₂₀ = |p₂p₀|`, `ℓ₀₁ = |p₀p₁|`. -/
theorem volume_real_triangle_eq_triArea {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) :
    volume.real (triangle p) = triArea (dist2 (p 1) (p 2)) (dist2 (p 2) (p 0)) (dist2 (p 0) (p 1)) := by
  rw [volume_real_triangle hD, D, abs_sda_eq_two_mul_triArea]
  ring

/-- A constant integrates to its value times the area. -/
theorem integral_triangle_const [CompleteSpace E] {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) (c : E) :
    ∫ _ in triangle p, c = (|D p| / 2) • c := by
  rw [setIntegral_const, volume_real_triangle hD]

/-! ## The triangle is the convex hull of its corners -/

/-- **The triangle is the convex hull of its three corners.** -/
theorem triangle_eq_convexHull (p : Fin 3 → ℝ × ℝ) :
    triangle p = convexHull ℝ {p 0, p 1, p 2} := by
  ext x
  rw [convexHull_insert (by simp : ({p 1, p 2} : Set (ℝ × ℝ)).Nonempty), convexHull_pair,
    mem_convexJoin]
  constructor
  · rintro ⟨⟨s, t⟩, ⟨hs, ht, hst⟩, rfl⟩
    simp only at hs ht hst
    rcases (add_nonneg hs ht).lt_or_eq with hpos | hzero
    · refine ⟨p 0, mem_singleton _, (s / (s + t)) • p 1 + (t / (s + t)) • p 2,
        ⟨s / (s + t), t / (s + t), by positivity, by positivity, by field_simp, rfl⟩,
        ⟨1 - (s + t), s + t, by linarith, hpos.le, by ring, ?_⟩⟩
      ext <;> simp only [triMap, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
        Prod.fst_sub, Prod.snd_sub, smul_eq_mul] <;> field_simp <;> ring
    · have hs0 : s = 0 := by linarith
      have ht0 : t = 0 := by linarith
      refine ⟨p 0, mem_singleton _, p 1, left_mem_segment _ _ _, ?_⟩
      have : triMap (p := p) (s, t) = p 0 := by simp [triMap, hs0, ht0]
      rw [this]
      exact left_mem_segment _ _ _
  · rintro ⟨a0, ha0, y, ⟨c, d, hc, hd, hcd, rfl⟩, ⟨a, b, ha, hb, hab, rfl⟩⟩
    rw [mem_singleton_iff] at ha0
    subst ha0
    refine ⟨(b * c, b * d), ⟨mul_nonneg hb hc, mul_nonneg hb hd, ?_⟩, ?_⟩
    · simp only
      nlinarith
    · ext
      · simp only [triMap, Prod.fst_add, Prod.smul_fst, Prod.fst_sub, smul_eq_mul]
        linear_combination -(p 0).1 * hab - (p 0).1 * b * hcd
      · simp only [triMap, Prod.snd_add, Prod.smul_snd, Prod.snd_sub, smul_eq_mul]
        linear_combination -(p 0).2 * hab - (p 0).2 * b * hcd

/-! ## Barycentric coordinates -/

/-- The barycentric (area) coordinate of `x` for the corner `i`: the signed double area of the
triangle with the corner `i` replaced by `x`, over `D`. -/
def bary (p : Fin 3 → ℝ × ℝ) (i : Fin 3) (x : ℝ × ℝ) : ℝ :=
  sda (Function.update p i x 0) (Function.update p i x 1) (Function.update p i x 2) / D p

section Bary

variable (p : Fin 3 → ℝ × ℝ) (x : ℝ × ℝ)

theorem bary_zero : bary p 0 x = sda x (p 1) (p 2) / D p := by
  simp [bary]

theorem bary_one : bary p 1 x = sda (p 0) x (p 2) / D p := by
  simp [bary]

theorem bary_two : bary p 2 x = sda (p 0) (p 1) x / D p := by
  simp [bary]

variable {p}

/-- On the triangle map the barycentric coordinates are `(1 − s − t, s, t)`. -/
theorem bary_triMap_zero (hD : D p ≠ 0) (q : ℝ × ℝ) : bary p 0 (triMap p q) = 1 - q.1 - q.2 := by
  rw [bary_zero, div_eq_iff hD]
  simp only [D, sda, triMap, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
    Prod.fst_sub, Prod.snd_sub, smul_eq_mul]
  ring

theorem bary_triMap_one (hD : D p ≠ 0) (q : ℝ × ℝ) : bary p 1 (triMap p q) = q.1 := by
  rw [bary_one, div_eq_iff hD]
  simp only [D, sda, triMap, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
    Prod.fst_sub, Prod.snd_sub, smul_eq_mul]
  ring

theorem bary_triMap_two (hD : D p ≠ 0) (q : ℝ × ℝ) : bary p 2 (triMap p q) = q.2 := by
  rw [bary_two, div_eq_iff hD]
  simp only [D, sda, triMap, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
    Prod.fst_sub, Prod.snd_sub, smul_eq_mul]
  ring

/-- The barycentric coordinates sum to one. -/
theorem sum_bary (hD : D p ≠ 0) : bary p 0 x + bary p 1 x + bary p 2 x = 1 := by
  rw [bary_zero, bary_one, bary_two, ← add_div, ← add_div, div_eq_one_iff_eq hD]
  simp only [D, sda]
  ring

/-- The barycentric coordinates reproduce the point: `Σ λᵢ(x) pᵢ = x`. -/
theorem sum_bary_smul (hD : D p ≠ 0) : bary p 0 x • p 0 + bary p 1 x • p 1 + bary p 2 x • p 2 = x := by
  rw [bary_zero, bary_one, bary_two]
  ext <;> simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul] <;>
    field_simp <;> simp only [D, sda] <;> ring

/-- `λᵢ(pⱼ) = δᵢⱼ`. -/
theorem bary_corner (hD : D p ≠ 0) (i j : Fin 3) : bary p i (p j) = if i = j then 1 else 0 := by
  fin_cases i <;> fin_cases j <;>
    simp only [Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk, Fin.isValue, bary_zero, bary_one,
      bary_two] <;> split_ifs with h <;>
    first
    | exact div_self hD
    | exact absurd h (by decide)
    | (rw [div_eq_zero_iff]; left; unfold sda; ring)

end Bary

/-! ## Cyclic relabelling of the corners -/

/-- The corners in the cyclic order `p₁, p₂, p₀`. -/
def rot (p : Fin 3 → ℝ × ℝ) : Fin 3 → ℝ × ℝ := fun i => p (i + 1)

@[simp] theorem rot_zero (p : Fin 3 → ℝ × ℝ) : rot p 0 = p 1 := rfl
@[simp] theorem rot_one (p : Fin 3 → ℝ × ℝ) : rot p 1 = p 2 := rfl
@[simp] theorem rot_two (p : Fin 3 → ℝ × ℝ) : rot p 2 = p 0 := rfl

theorem D_rot (p : Fin 3 → ℝ × ℝ) : D (rot p) = D p := by
  simp only [D, rot_zero, rot_one, rot_two, sda_rotate]

theorem triangle_rot (p : Fin 3 → ℝ × ℝ) : triangle (rot p) = triangle p := by
  rw [triangle_eq_convexHull, triangle_eq_convexHull, rot_zero, rot_one, rot_two]
  congr 1
  ext y
  simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
  tauto

theorem bary_rot_zero (p : Fin 3 → ℝ × ℝ) (x : ℝ × ℝ) : bary (rot p) 0 x = bary p 1 x := by
  rw [bary_zero, bary_one, D_rot, rot_one, rot_two, sda_rotate]

theorem bary_rot_one (p : Fin 3 → ℝ × ℝ) (x : ℝ × ℝ) : bary (rot p) 1 x = bary p 2 x := by
  rw [bary_one, bary_two, D_rot, rot_zero, rot_two, sda_rotate]

theorem bary_rot_two (p : Fin 3 → ℝ × ℝ) (x : ℝ × ℝ) : bary (rot p) 2 x = bary p 0 x := by
  rw [bary_two, bary_zero, D_rot, rot_zero, rot_one, sda_rotate]

/-! ## The slice formula and the area-coordinate rule -/

/-- **The slice formula at the corner `1`.** For a continuous `G : ℝ → E`,
`∫_T λ₀^{e₀} λ₁^{e₁} λ₂^{e₂} G(λ₁) = |D| e₂! e₀!/(e₂ + e₀ + 1)! · ∫₀¹ s^{e₁} (1 − s)^{e₂+e₀+1} G(s) ds`:
on the slice `λ₁ = s` the other two coordinates are `t` and `(1 − s) − t`, and the inner integral
is a Beta integral. -/
theorem integral_bary_slice_one [CompleteSpace E] {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0)
    (e : Fin 3 → ℕ) (G : ℝ → E) (hG : Continuous G) :
    ∫ x in triangle p, (bary p 0 x ^ e 0 * bary p 1 x ^ e 1 * bary p 2 x ^ e 2) • G (bary p 1 x) =
      (|D p| * (((e 2)! * (e 0)! : ℕ) / ((e 2 + e 0 + 1)! : ℕ) : ℝ)) •
        ∫ s in (0 : ℝ)..1, (s ^ e 1 * (1 - s) ^ (e 2 + e 0 + 1)) • G s := by
  rw [integral_triangle hD]
  simp_rw [bary_triMap_zero hD, bary_triMap_one hD, bary_triMap_two hD]
  rw [integral_stdTri _ (by fun_prop : Continuous fun q : ℝ × ℝ =>
    ((1 - q.1 - q.2) ^ e 0 * q.1 ^ e 1 * q.2 ^ e 2) • G q.1).continuousOn]
  set C : ℝ := (((e 2)! * (e 0)! : ℕ) / ((e 2 + e 0 + 1)! : ℕ) : ℝ)
  have inner : (fun s : ℝ => ∫ t in (0 : ℝ)..(1 - s),
      ((1 - (s, t).1 - (s, t).2) ^ e 0 * (s, t).1 ^ e 1 * (s, t).2 ^ e 2) • G (s, t).1) =
      fun s => C • ((s ^ e 1 * (1 - s) ^ (e 2 + e 0 + 1)) • G s) := by
    funext s
    simp only
    rw [intervalIntegral.integral_smul_const]
    have e1 : (fun t : ℝ => (1 - s - t) ^ e 0 * s ^ e 1 * t ^ e 2) =
        fun t => s ^ e 1 * (t ^ e 2 * ((1 - s) - t) ^ e 0) := by
      funext t; ring
    rw [e1, intervalIntegral.integral_const_mul, integral_pow_mul_sub_pow, smul_smul]
    congr 1
    ring
  rw [inner, intervalIntegral.integral_smul, smul_smul]

/-- **The slice formula at any corner** `c`, with `(c, c + 1, c + 2)` the cyclic order (the formula
of `monomial_integral` in `crates/ddg-intrinsic-apps/src/vector/knoeppel.rs`, for
`G(s) = e^{iΩs}`): `∫_T λ₀^{e₀} λ₁^{e₁} λ₂^{e₂} G(λ_c) =
|D| q! r!/(q + r + 1)! · ∫₀¹ s^{e_c} (1 − s)^{q+r+1} G(s) ds`, `q = e_{c+1}`, `r = e_{c+2}`. -/
theorem integral_bary_slice [CompleteSpace E] {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0)
    (e : Fin 3 → ℕ) (G : ℝ → E) (hG : Continuous G) (c : Fin 3) :
    ∫ x in triangle p, (bary p 0 x ^ e 0 * bary p 1 x ^ e 1 * bary p 2 x ^ e 2) • G (bary p c x) =
      (|D p| * (((e (c + 1))! * (e (c + 2))! : ℕ) / ((e (c + 1) + e (c + 2) + 1)! : ℕ) : ℝ)) •
        ∫ s in (0 : ℝ)..1, (s ^ e c * (1 - s) ^ (e (c + 1) + e (c + 2) + 1)) • G s := by
  fin_cases c
  · -- `c = 0`: the slice formula at the corner `1` of `rot (rot p)`.
    have h := integral_bary_slice_one (p := rot (rot p)) (by rwa [D_rot, D_rot])
      (fun i => e (i + 2)) G hG
    simp only [bary_rot_zero, bary_rot_one, bary_rot_two, D_rot, triangle_rot, Fin.reduceAdd,
      Fin.isValue] at h
    simp only [Fin.zero_eta, Fin.isValue, zero_add]
    rw [← h]
    congr 1
    funext x
    congr 1
    ring
  · exact integral_bary_slice_one hD e G hG
  · -- `c = 2`: the slice formula at the corner `1` of `rot p`.
    have h := integral_bary_slice_one (p := rot p) (by rwa [D_rot]) (fun i => e (i + 1)) G hG
    simp only [bary_rot_zero, bary_rot_one, bary_rot_two, D_rot, triangle_rot, Fin.reduceAdd,
      Fin.isValue] at h
    simp only [Fin.reduceFinMk, Fin.isValue, Fin.reduceAdd]
    rw [← h]
    congr 1
    funext x
    congr 1
    ring

/-- **The area-coordinate rule** (Zienkiewicz–Taylor; Ciarlet Sec. 4.1):
`∫_T λ₀ᵃ λ₁ᵇ λ₂ᶜ = |D| a! b! c!/(a + b + c + 2)! = 2|T| a! b! c!/(a + b + c + 2)!`. -/
theorem integral_bary_monomial {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) (a b c : ℕ) :
    ∫ x in triangle p, bary p 0 x ^ a * bary p 1 x ^ b * bary p 2 x ^ c =
      |D p| * ((a ! * b ! * c ! : ℕ) / ((a + b + c + 2)! : ℕ) : ℝ) := by
  have h := integral_bary_slice_one hD ![a, b, c] (fun _ => (1 : ℝ)) continuous_const
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, smul_eq_mul,
    mul_one, Matrix.head_cons, Matrix.tail_cons] at h
  rw [h]
  have hb := integral_pow_mul_sub_pow 1 b (c + a + 1)
  rw [one_pow, one_mul] at hb
  rw [hb, show b + (c + a + 1) + 1 = a + b + c + 2 by ring]
  push_cast
  have h1 : ((c + a + 1)! : ℝ) ≠ 0 := by positivity
  have h2 : ((a + b + c + 2)! : ℝ) ≠ 0 := by positivity
  field_simp

/-- The area-coordinate rule with the area of the triangle: `2|T| a! b! c!/(a + b + c + 2)!`. -/
theorem integral_bary_monomial_area {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) (a b c : ℕ) :
    ∫ x in triangle p, bary p 0 x ^ a * bary p 1 x ^ b * bary p 2 x ^ c =
      2 * volume.real (triangle p) * ((a ! * b ! * c ! : ℕ) / ((a + b + c + 2)! : ℕ) : ℝ) := by
  rw [integral_bary_monomial hD, volume_real_triangle hD]
  ring

end TriangleIntegral

end DDG
