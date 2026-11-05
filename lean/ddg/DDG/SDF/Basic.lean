import DDG.Geometry.Vec3

/-!
# Points as Tori: vectors, frames and gradients

"Points as Tori: Fast Pointwise Signed Distance for Point Clouds", Feng, Gkioulekas, Crane,
ACM TOG 45(4), Article 53, 2026 ("the paper"); contract `crates/ddg-points-as-tori/SPEC.md`
("SPEC", formula ids `pat-…`, §13).

Vectors are `RC.V3 = Fin 3 → ℝ` with `dotProduct` (`⬝ᵥ`), `crossProduct` (`⨯₃`) and the Euclidean
norm `RC.vnorm v = √(v·v)`. Gradients are Fréchet derivatives `HasFDerivAt f (dualOf ∇f) x`, where
`dualOf a` is the linear form `h ↦ a·h` (Mathlib's norm on `Fin 3 → ℝ` is the sup norm; the
derivative does not depend on the choice of norm).

A `Frame` is an orthonormal right-handed frame `(s, t, n)` (`s × t = n`); `comb α β γ` is the
vector `α s + β t + γ n`. Dot and cross products of combinations reduce to coordinates
(`comb_dot`, `comb_cross`).
-/

noncomputable section

namespace DDG

namespace PAT

open Matrix Real
open RC (V3 vnorm vnorm_nonneg vnorm_sq vnorm_mul_self vnorm_eq_zero_iff vnorm_pos_iff vnorm_smul
  dot_self_nonneg)

/-! ## Linear forms and gradients -/

/-- The linear form `h ↦ a·h` (the differential with gradient `a`). -/
def dualOf (a : V3) : V3 →L[ℝ] ℝ :=
  LinearMap.toContinuousLinearMap
    { toFun := fun h => a ⬝ᵥ h
      map_add' := fun x y => dotProduct_add a x y
      map_smul' := fun c x => by simp [dotProduct_smul] }

@[simp] theorem dualOf_apply (a h : V3) : dualOf a h = a ⬝ᵥ h := rfl

theorem dualOf_add (a b : V3) : dualOf (a + b) = dualOf a + dualOf b := by
  ext h; simp [add_dotProduct]

theorem dualOf_smul (c : ℝ) (a : V3) : dualOf (c • a) = c • dualOf a := by
  ext h; simp [smul_dotProduct]

/-- The linear map `h ↦ h × u`. -/
def crossRight (u : V3) : V3 →L[ℝ] V3 :=
  LinearMap.toContinuousLinearMap
    { toFun := fun h => h ⨯₃ u
      map_add' := fun x y => by simp [map_add, LinearMap.add_apply]
      map_smul' := fun c x => by simp [map_smul, LinearMap.smul_apply] }

@[simp] theorem crossRight_apply (u h : V3) : crossRight u h = h ⨯₃ u := rfl

theorem hasFDerivAt_dot_const (a x : V3) : HasFDerivAt (fun y => a ⬝ᵥ y) (dualOf a) x :=
  (dualOf a).hasFDerivAt

/-- `∇(y·y) = 2y`. -/
theorem hasFDerivAt_dot_self (x : V3) :
    HasFDerivAt (fun y : V3 => y ⬝ᵥ y) (dualOf ((2 : ℝ) • x)) x := by
  have hc : ∀ i : Fin 3, HasFDerivAt (fun y : V3 => y i)
      (ContinuousLinearMap.proj (R := ℝ) (φ := fun _ : Fin 3 => ℝ) i) x :=
    fun i => hasFDerivAt_apply (𝕜 := ℝ) (F' := fun _ : Fin 3 => ℝ) i x
  have := HasFDerivAt.sum (u := Finset.univ) (fun i _ => (hc i).mul (hc i))
  simp only [dotProduct]
  convert this using 1
  · funext y; simp
  · ext h
    simp [dotProduct, Fin.sum_univ_three]
    ring

/-- Chain rule for the Euclidean norm: `∇|f| = (Df)ᵀ f/|f|` where `f ≠ 0`. -/
theorem hasFDerivAt_vnorm_comp {f : V3 → V3} {f' : V3 →L[ℝ] V3} {x : V3} (hf : HasFDerivAt f f' x)
    (h0 : f x ≠ 0) :
    HasFDerivAt (fun y => vnorm (f y)) ((vnorm (f x))⁻¹ • (dualOf (f x)).comp f') x := by
  have hpos : 0 < vnorm (f x) := vnorm_pos_iff.mpr h0
  have hne : f x ⬝ᵥ f x ≠ 0 := by rw [← vnorm_mul_self]; positivity
  have h2 := (hasFDerivAt_dot_self (f x)).comp x hf
  have := h2.sqrt hne
  refine this.congr_fderiv ?_
  ext h
  simp only [_root_.smul_apply, ContinuousLinearMap.comp_apply, dualOf_apply,
    smul_dotProduct, smul_eq_mul, vnorm, Function.comp]
  field_simp

/-! ## Orthonormal frames -/

/-- An orthonormal right-handed frame `(s, t, n)` of `ℝ³`. -/
structure Frame where
  s : V3
  t : V3
  n : V3
  ss : s ⬝ᵥ s = 1
  tt : t ⬝ᵥ t = 1
  nn : n ⬝ᵥ n = 1
  st : s ⬝ᵥ t = 0
  sn : s ⬝ᵥ n = 0
  tn : t ⬝ᵥ n = 0
  cross : s ⨯₃ t = n

namespace Frame

variable (F : Frame)

theorem ts : F.t ⬝ᵥ F.s = 0 := by rw [dotProduct_comm]; exact F.st
theorem ns : F.n ⬝ᵥ F.s = 0 := by rw [dotProduct_comm]; exact F.sn
theorem nt : F.n ⬝ᵥ F.t = 0 := by rw [dotProduct_comm]; exact F.tn

theorem cross_tn : F.t ⨯₃ F.n = F.s := by
  rw [← F.cross, cross_cross_eq_smul_sub_smul', F.tt, F.st]; simp

theorem cross_ns : F.n ⨯₃ F.s = F.t := by
  rw [← F.cross]
  have := cross_cross_eq_smul_sub_smul F.s F.t F.s
  rw [this, F.ss, F.ts]; simp

theorem cross_ts : F.t ⨯₃ F.s = -F.n := by
  rw [← cross_anticomm, F.cross]
theorem cross_nt : F.n ⨯₃ F.t = -F.s := by
  rw [← cross_anticomm, F.cross_tn]
theorem cross_sn : F.s ⨯₃ F.n = -F.t := by
  rw [← cross_anticomm, F.cross_ns]

/-- The vector with frame coordinates `(α, β, γ)`. -/
def comb (α β γ : ℝ) : V3 := α • F.s + β • F.t + γ • F.n

/-- Dot products of combinations are coordinate dot products. -/
theorem comb_dot (α β γ α' β' γ' : ℝ) :
    F.comb α β γ ⬝ᵥ F.comb α' β' γ' = α * α' + β * β' + γ * γ' := by
  simp only [comb, add_dotProduct, dotProduct_add, smul_dotProduct, dotProduct_smul, F.ss, F.tt,
    F.nn, F.st, F.sn, F.tn, F.ts, F.ns, F.nt, smul_eq_mul]
  ring

/-- Cross products of combinations are coordinate cross products. -/
theorem comb_cross (α β γ α' β' γ' : ℝ) :
    F.comb α β γ ⨯₃ F.comb α' β' γ' =
      F.comb (β * γ' - γ * β') (γ * α' - α * γ') (α * β' - β * α') := by
  simp only [comb, map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply,
    cross_self, F.cross, F.cross_tn, F.cross_ns, F.cross_ts, F.cross_nt, F.cross_sn]
  ext i
  simp only [Pi.add_apply, Pi.smul_apply, Pi.neg_apply, Pi.zero_apply, smul_eq_mul]
  ring

theorem comb_add (α β γ α' β' γ' : ℝ) :
    F.comb α β γ + F.comb α' β' γ' = F.comb (α + α') (β + β') (γ + γ') := by
  simp only [comb, add_smul]; abel

theorem comb_smul (c α β γ : ℝ) : c • F.comb α β γ = F.comb (c * α) (c * β) (c * γ) := by
  simp only [comb, smul_add, smul_smul]

theorem comb_sub (α β γ α' β' γ' : ℝ) :
    F.comb α β γ - F.comb α' β' γ' = F.comb (α - α') (β - β') (γ - γ') := by
  simp only [comb, sub_smul]; abel

theorem vnorm_comb (α β γ : ℝ) : vnorm (F.comb α β γ) = √(α ^ 2 + β ^ 2 + γ ^ 2) := by
  rw [vnorm, comb_dot]; ring_nf

theorem comb_s : F.comb 1 0 0 = F.s := by simp [comb]
theorem comb_t : F.comb 0 1 0 = F.t := by simp [comb]
theorem comb_n : F.comb 0 0 1 = F.n := by simp [comb]

end Frame

end PAT

end DDG
