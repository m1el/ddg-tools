import DdgProofs.RepulsiveCurves.Kernel

/-!
# Repulsive Curves: the near-pair quadrature (WP16; SPEC §3.4, §17o; not in the paper)

The pair integral `∫_I ∫_J k(x, y, T_I) dx dy` that the paper's endpoint rule (Eq. 18)
approximates is evaluated accurately for close edge pairs and blended into the paper's term.
Machine-checked here:

* `rc-near-h-constant`: for a unit `T`, `|T × (d + uT)| = |T × d|` and `P_T (d + uT) = P_T d`;
  with `x = a + uT`, `|x − y|² = (u − T·(y − a))² + |P_T (y − a)|²`. So along the edge `I`
  the numerator of the kernel is the constant `h^α`, `h = |P_T (y − a)|`, and the integrand is
  `h^α ((u − u₀)² + h²)^{−β/2}`.
* `rc-near-inner`: the antiderivative after `u = h tan θ`:
  `d/du [h^{1−β} S(arctan(u/h))] = (u² + h²)^{−β/2}` for any `S` with `S' = cos^{β−2}`, hence
  `∫_a^b (u² + h²)^{−β/2} du = h^{1−β} (S(arctan(b/h)) − S(arctan(a/h)))`, and
  `S(θ) = ∫_0^θ cos^{β−2}` has that derivative.
* `rc-near-outer`: the sinh substitution of the graded outer quadrature,
  `∫_c^{c+L} f = ∫_0^{arsinh(L/δ)} f(c + δ sinh w) δ cosh w dw` (and the mirrored half).
* `rc-near-blend`: the quintic `W(τ) = τ³(10 − 15τ + 6τ²)`: `W(0) = 0`, `W(1) = 1`,
  `W' = 30τ²(1 − τ)²`, `W'' = 60τ(1 − τ)(1 − 2τ)`, both zero at `0` and `1` (C² blend), and the
  weight `w(r) = W((c₂ − r)/(c₂ − c₁))` with `w' = −30τ²(1−τ)²/(c₂ − c₁)`.
* `rc-near-differential`: the chain rule of the blended pair term
  `e = e_p + w(r)(e_n − e_p)` along any motion.
* `rc-near-divergence`: for two perpendicular segments of half-length 1 crossing at distance
  `d ≤ 1`, the pair integral is at least `4 · 3^{−β/2} d^{2+α−β}` — it diverges as `d → 0`
  because `2 + α − β < 0` (−1 for (3, 6), −½ for (2, 4.5)); the paper's endpoint rule stays
  bounded there.

The series that evaluates `S` and the Gauss–Legendre error of the graded rule are cited
(`CITED.md`: `rc-near-cos-series`, `rc-near-closest`).
-/

noncomputable section

namespace DDG

namespace RC

open Matrix Real

/-! ## `rc-near-h-constant` -/

/-- **`P_T (d + uT) = P_T d`** for a unit `T`.

[formula: rc-near-h-constant] -/
theorem proj_add_smul_self {T : V3} (hT : T ⬝ᵥ T = 1) (d : V3) (u : ℝ) :
    proj T (d + u • T) = proj T d := by
  unfold proj
  rw [dotProduct_add, dotProduct_smul, smul_eq_mul, hT, mul_one, add_smul]
  abel

/-- **`|T × (d + uT)| = |T × d|`**: the numerator of the kernel is constant along an edge.

[formula: rc-near-h-constant] -/
theorem cross_add_smul_self (T d : V3) (u : ℝ) :
    crossProduct T (d + u • T) = crossProduct T d := by
  rw [map_add, map_smul, cross_self, smul_zero, add_zero]

/-- With `x = a + uT` (unit `T`) and `r = y − a`:
`|x − y|² = (u − T·r)² + |P_T r|²`, and `|T × (x − y)| = |P_T r|`.

[formula: rc-near-h-constant] -/
theorem edge_point_distance {T : V3} (hT : T ⬝ᵥ T = 1) (a y : V3) (u : ℝ) :
    vnorm (a + u • T - y) ^ 2 = (u - T ⬝ᵥ (y - a)) ^ 2 + vnorm (proj T (y - a)) ^ 2 ∧
      vnorm (crossProduct T (a + u • T - y)) = vnorm (proj T (y - a)) := by
  constructor
  · rw [proj_vnorm_sq hT, vnorm_sq]
    have e : a + u • T - y = u • T - (y - a) := by abel
    rw [e]
    simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, smul_eq_mul, hT,
      dotProduct_comm (y - a) T]
    ring
  · have e : a + u • T - y = -(y - a) + u • T := by abel
    rw [e, cross_add_smul_self, map_neg, vnorm_neg]
    exact (cross_vnorm_sq_unit hT _).2

/-! ## `rc-near-inner` -/

/-- The integrand of the inner integral after the substitution: `cos(arctan x)^m`
`= (1 + x²)^{−m/2}`. -/
theorem cos_arctan_rpow (x m : ℝ) : Real.cos (Real.arctan x) ^ m = (1 + x ^ 2) ^ (-m / 2) := by
  have hA : 0 < 1 + x ^ 2 := by positivity
  rw [Real.cos_arctan, one_div, Real.inv_rpow (Real.sqrt_nonneg _), Real.sqrt_eq_rpow,
    ← Real.rpow_mul hA.le, ← Real.rpow_neg hA.le]
  congr 1
  ring

/-- **The inner antiderivative.** For `h > 0` and any `S` with `S' = cos^m`, `m = β − 2`:
`d/du [h^{1−β} S(arctan(u/h))] = (u² + h²)^{−β/2}` (the substitution `u = h tan θ`).

[formula: rc-near-inner] -/
theorem near_inner_hasDerivAt {m β h : ℝ} (hβ : β = m + 2) (hh : 0 < h) {S : ℝ → ℝ}
    (hS : ∀ θ, HasDerivAt S (Real.cos θ ^ m) θ) (u : ℝ) :
    HasDerivAt (fun u => h ^ (1 - β) * S (Real.arctan (u / h))) ((u ^ 2 + h ^ 2) ^ (-β / 2)) u := by
  have h1 : HasDerivAt (fun u => u / h) (1 / h) u := (hasDerivAt_id u).div_const h
  have h2 := (Real.hasDerivAt_arctan (u / h)).comp u h1
  have h3 := (hS (Real.arctan (u / h))).comp u h2
  have h4 := h3.const_mul (h ^ (1 - β))
  refine h4.congr_deriv ?_
  rw [cos_arctan_rpow]
  set A := 1 + (u / h) ^ 2 with hAdef
  have hA : 0 < A := by positivity
  have hsum : u ^ 2 + h ^ 2 = h ^ 2 * A := by
    rw [hAdef]; field_simp; ring
  rw [hsum, Real.mul_rpow (by positivity) hA.le]
  have hh2 : (h ^ 2) ^ (-β / 2) = h ^ (-β) := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hh.le]; congr 1; push_cast; ring
  have hAβ : A ^ (-β / 2) = A ^ (-m / 2) * A⁻¹ := by
    rw [hβ, show -(m + 2) / 2 = -m / 2 + (-1 : ℝ) by ring, Real.rpow_add hA, Real.rpow_neg_one]
  have hhβ : h ^ (-β) = h ^ (1 - β) * h⁻¹ := by
    rw [show -β = (1 - β) + (-1 : ℝ) by ring, Real.rpow_add hh, Real.rpow_neg_one]
  rw [hh2, hAβ, hhβ]
  simp only [one_div]
  ring

/-- `S(θ) = ∫_0^θ cos^m` has `S' = cos^m` (`m ≥ 0`), so it is an admissible `S` above.

[formula: rc-near-inner] -/
theorem cos_power_integral_hasDerivAt {m : ℝ} (hm : 0 ≤ m) (θ : ℝ) :
    HasDerivAt (fun θ => ∫ φ in (0 : ℝ)..θ, Real.cos φ ^ m) (Real.cos θ ^ m) θ := by
  have hc : Continuous fun φ => Real.cos φ ^ m :=
    Real.continuous_cos.rpow_const fun _ => Or.inr hm
  exact (hc.integral_hasStrictDerivAt 0 θ).hasDerivAt

/-- **The inner integral in closed form**:
`∫_a^b (u² + h²)^{−β/2} du = h^{1−β} (S(arctan(b/h)) − S(arctan(a/h)))`, `S = ∫_0^· cos^{β−2}`.
With `u ↦ u − u₀` and the factor `h^α` this is `F(y) = h^{1+α−β}(S(θ₁) − S(θ₀))`
(`edge_point_distance`).

[formula: rc-near-inner] -/
theorem near_inner_integral {m β h : ℝ} (hβ : β = m + 2) (hm : 0 ≤ m) (hh : 0 < h) (a b : ℝ) :
    ∫ u in a..b, (u ^ 2 + h ^ 2) ^ (-β / 2) =
      h ^ (1 - β) * ((∫ φ in (0 : ℝ)..Real.arctan (b / h), Real.cos φ ^ m) -
        ∫ φ in (0 : ℝ)..Real.arctan (a / h), Real.cos φ ^ m) := by
  have hd := near_inner_hasDerivAt hβ hh (cos_power_integral_hasDerivAt hm)
  have hcont : Continuous fun u : ℝ => (u ^ 2 + h ^ 2) ^ (-β / 2) :=
    (by fun_prop : Continuous fun u : ℝ => u ^ 2 + h ^ 2).rpow_const
      fun u => Or.inl (by positivity)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun u _ => hd u)
    (hcont.intervalIntegrable a b)]
  ring

/-! ## `rc-near-outer` -/

/-- **The sinh grading** toward the left end `c` of a half-panel of length `L ≥ 0` with scale
`δ > 0`: `∫_c^{c+L} f = ∫_0^{arsinh(L/δ)} f(c + δ sinh w) δ cosh w dw`.

[formula: rc-near-outer] -/
theorem sinh_grading_left {f : ℝ → ℝ} (hf : Continuous f) (c L δ : ℝ) (hδ : 0 < δ) :
    ∫ t in c..c + L, f t =
      ∫ w in (0 : ℝ)..Real.arsinh (L / δ), f (c + δ * Real.sinh w) * (δ * Real.cosh w) := by
  have hg : ∀ w ∈ Set.uIcc (0 : ℝ) (Real.arsinh (L / δ)),
      HasDerivAt (fun w => c + δ * Real.sinh w) (δ * Real.cosh w) w := fun w _ =>
    ((Real.hasDerivAt_sinh w).const_mul δ).const_add c
  have hg' : ContinuousOn (fun w => δ * Real.cosh w) (Set.uIcc 0 (Real.arsinh (L / δ))) :=
    (continuous_const.mul Real.continuous_cosh).continuousOn
  have := intervalIntegral.integral_comp_mul_deriv hg hg' hf
  simp only [Function.comp, Real.sinh_zero, mul_zero, add_zero, Real.sinh_arsinh] at this
  have hL : δ * (L / δ) = L := by field_simp
  rw [hL] at this
  exact this.symm

/-- The mirrored half (graded toward the right end `c`):
`∫_{c−L}^{c} f = ∫_0^{arsinh(L/δ)} f(c − δ sinh w) δ cosh w dw`.

[formula: rc-near-outer] -/
theorem sinh_grading_right {f : ℝ → ℝ} (hf : Continuous f) (c L δ : ℝ) (hδ : 0 < δ) :
    ∫ t in c - L..c, f t =
      ∫ w in (0 : ℝ)..Real.arsinh (L / δ), f (c - δ * Real.sinh w) * (δ * Real.cosh w) := by
  have hg : ∀ w ∈ Set.uIcc (0 : ℝ) (Real.arsinh (L / δ)),
      HasDerivAt (fun w => c - δ * Real.sinh w) (-(δ * Real.cosh w)) w := fun w _ =>
    ((Real.hasDerivAt_sinh w).const_mul δ).const_sub c
  have hg' : ContinuousOn (fun w => -(δ * Real.cosh w)) (Set.uIcc 0 (Real.arsinh (L / δ))) :=
    (continuous_const.mul Real.continuous_cosh).neg.continuousOn
  have := intervalIntegral.integral_comp_mul_deriv hg hg' hf
  simp only [Function.comp, Real.sinh_zero, mul_zero, sub_zero, Real.sinh_arsinh] at this
  have hL : δ * (L / δ) = L := by field_simp
  rw [hL] at this
  rw [intervalIntegral.integral_symm c (c - L), ← this, ← intervalIntegral.integral_neg]
  congr 1
  funext w
  ring

/-! ## `rc-near-blend` -/

/-- The quintic smoothstep `W(τ) = 10τ³ − 15τ⁴ + 6τ⁵ = τ³(10 − 15τ + 6τ²)`. -/
def smoother (τ : ℝ) : ℝ := 10 * τ ^ 3 - 15 * τ ^ 4 + 6 * τ ^ 5

/-- **`W(0) = 0`, `W(1) = 1`, `W = τ³(10 − 15τ + 6τ²)`.**

[formula: rc-near-blend] -/
theorem smoother_values :
    smoother 0 = 0 ∧ smoother 1 = 1 ∧ ∀ τ, smoother τ = τ ^ 3 * (10 - 15 * τ + 6 * τ ^ 2) := by
  refine ⟨by norm_num [smoother], by norm_num [smoother], fun τ => by unfold smoother; ring⟩

/-- **`W' = 30 τ² (1 − τ)²`** (zero at `0` and `1`).

[formula: rc-near-blend] -/
theorem smoother_hasDerivAt (τ : ℝ) : HasDerivAt smoother (30 * τ ^ 2 * (1 - τ) ^ 2) τ := by
  have h := (((hasDerivAt_pow 3 τ).const_mul 10).sub ((hasDerivAt_pow 4 τ).const_mul 15)).add
    ((hasDerivAt_pow 5 τ).const_mul 6)
  exact h.congr_deriv (by norm_num; ring)

/-- **`W'' = 60 τ (1 − τ)(1 − 2τ)`** (zero at `0` and `1`: the blend is C²).

[formula: rc-near-blend] -/
theorem smoother_deriv_hasDerivAt (τ : ℝ) :
    HasDerivAt (fun τ => 30 * τ ^ 2 * (1 - τ) ^ 2) (60 * τ * (1 - τ) * (1 - 2 * τ)) τ := by
  have h := ((hasDerivAt_pow 2 τ).const_mul 30).mul
    ((((hasDerivAt_id τ).const_sub 1)).pow 2)
  exact h.congr_deriv (by simp only [id]; norm_num; ring)

/-- **The blend weight** `w(r) = W((c₂ − r)/(c₂ − c₁))` (inside the band) has
`w'(r) = −30 τ² (1 − τ)² / (c₂ − c₁)`.

[formula: rc-near-blend] -/
theorem blend_weight_hasDerivAt {c₁ c₂ : ℝ} (r : ℝ) :
    HasDerivAt (fun r => smoother ((c₂ - r) / (c₂ - c₁)))
      (-(30 * ((c₂ - r) / (c₂ - c₁)) ^ 2 * (1 - (c₂ - r) / (c₂ - c₁)) ^ 2) / (c₂ - c₁)) r := by
  have hτ : HasDerivAt (fun r => (c₂ - r) / (c₂ - c₁)) (-1 / (c₂ - c₁)) r :=
    ((hasDerivAt_id r).const_sub c₂).div_const (c₂ - c₁) |>.congr_deriv (by simp [neg_div])
  refine ((smoother_hasDerivAt ((c₂ - r) / (c₂ - c₁))).comp r hτ).congr_deriv ?_
  ring

/-! ## `rc-near-differential` -/

/-- **The chain rule of the blended pair term** `e = e_p + w(r)(e_n − e_p)` along a motion `s`:
`e' = (1 − w) e_p' + w e_n' + w'(r) r' (e_n − e_p)`.

[formula: rc-near-differential] -/
theorem blended_pair_hasDerivAt {ep en r : ℝ → ℝ} {w : ℝ → ℝ} {ep' en' r' w' s : ℝ}
    (hp : HasDerivAt ep ep' s) (hn : HasDerivAt en en' s) (hr : HasDerivAt r r' s)
    (hw : HasDerivAt w w' (r s)) :
    HasDerivAt (fun s => ep s + w (r s) * (en s - ep s))
      ((1 - w (r s)) * ep' + w (r s) * en' + w' * r' * (en s - ep s)) s := by
  have h := hp.add ((hw.comp s hr).mul (hn.sub hp))
  refine h.congr_deriv ?_
  simp only [Function.comp, Pi.sub_apply]
  ring

/-! ## `rc-near-divergence` -/

/-- The pair integral of two perpendicular segments of half-length 1 crossing in their middles
at distance `d`: `I = [−1, 1] e_x`, `J = d e_z + [−1, 1] e_y`, `T_I = e_x`, so
`|T_I × (x − y)| = √(v² + d²)` and `|x − y|² = u² + v² + d²`. -/
def perpPair (α β d : ℝ) : ℝ :=
  ∫ v in (-1 : ℝ)..1, ∫ u in (-1 : ℝ)..1, (v ^ 2 + d ^ 2) ^ (α / 2) * (u ^ 2 + v ^ 2 + d ^ 2) ^ (-β / 2)

theorem perpPair_integrand_continuous {α β d : ℝ} (hd : 0 < d) :
    Continuous (Function.uncurry fun (v u : ℝ) =>
      (v ^ 2 + d ^ 2) ^ (α / 2) * (u ^ 2 + v ^ 2 + d ^ 2) ^ (-β / 2)) := by
  refine Continuous.mul ?_ ?_
  · exact (by fun_prop : Continuous fun p : ℝ × ℝ => p.1 ^ 2 + d ^ 2).rpow_const
      fun p => Or.inl (by positivity)
  · exact (by fun_prop : Continuous fun p : ℝ × ℝ => p.2 ^ 2 + p.1 ^ 2 + d ^ 2).rpow_const
      fun p => Or.inl (by positivity)

/-- **Divergence of the pair integral at a mid-edge crossing**: for `0 < d ≤ 1`, `α, β ≥ 0`,
`perpPair α β d ≥ 4 · 3^{−β/2} d^{2+α−β}` — unbounded as `d → 0` when `β > α + 2` (the exponent
is −1 for (3, 6), −½ for (2, 4.5)). The paper's endpoint rule (Eq. 18) of the same pair is bounded
as `d → 0` (all four endpoint pairs are at distance `≥ √2`).

[formula: rc-near-divergence] -/
theorem perpPair_lower_bound {α β d : ℝ} (hα : 0 ≤ α) (hβ : 0 ≤ β) (hd : 0 < d) (hd1 : d ≤ 1) :
    4 * (3 : ℝ) ^ (-β / 2) * d ^ (2 + α - β) ≤ perpPair α β d := by
  set f : ℝ → ℝ → ℝ := fun v u =>
    (v ^ 2 + d ^ 2) ^ (α / 2) * (u ^ 2 + v ^ 2 + d ^ 2) ^ (-β / 2) with hf
  have hfc := perpPair_integrand_continuous (α := α) (β := β) hd
  have hf0 : ∀ v u, 0 ≤ f v u := fun v u => by
    simp only [hf]; positivity
  set c₀ : ℝ := (d ^ 2) ^ (α / 2) * (3 * d ^ 2) ^ (-β / 2) with hc₀
  -- Pointwise bound on the square |u|, |v| ≤ d.
  have hbound : ∀ v u, |v| ≤ d → |u| ≤ d → c₀ ≤ f v u := by
    intro v u hv hu
    have hv2 : v ^ 2 ≤ d ^ 2 := sq_le_sq' (abs_le.mp hv).1 (abs_le.mp hv).2
    have hu2 : u ^ 2 ≤ d ^ 2 := sq_le_sq' (abs_le.mp hu).1 (abs_le.mp hu).2
    have hd2 : 0 < d ^ 2 := by positivity
    refine mul_le_mul ?_ ?_ (by positivity) (by positivity)
    · exact Real.rpow_le_rpow hd2.le (by nlinarith [sq_nonneg v]) (by linarith)
    · exact Real.rpow_le_rpow_of_nonpos (by positivity) (by nlinarith)
        (by linarith)
  -- The inner integral.
  set g : ℝ → ℝ := fun v => ∫ u in (-1 : ℝ)..1, f v u with hg
  have hgc : Continuous g :=
    intervalIntegral.continuous_parametric_intervalIntegral_of_continuous' hfc (-1) 1
  have hg0 : ∀ v, 0 ≤ g v := fun v =>
    intervalIntegral.integral_nonneg (by norm_num) fun u _ => hf0 v u
  have hgv : ∀ v, |v| ≤ d → 2 * d * c₀ ≤ g v := by
    intro v hv
    have hfv : Continuous (f v) := hfc.comp (Continuous.prodMk_right v)
    calc 2 * d * c₀ = ∫ _u in -d..d, c₀ := by
            rw [intervalIntegral.integral_const, smul_eq_mul]; ring
      _ ≤ ∫ u in -d..d, f v u :=
          intervalIntegral.integral_mono_on (by linarith) intervalIntegrable_const
            (hfv.intervalIntegrable _ _) fun u hu => hbound v u hv (abs_le.mpr hu)
      _ ≤ g v :=
          intervalIntegral.integral_mono_interval (by linarith) (by linarith) hd1
            (Filter.Eventually.of_forall fun u => hf0 v u) (hfv.intervalIntegrable _ _)
  have hmain : 2 * d * (2 * d * c₀) ≤ perpPair α β d := by
    calc 2 * d * (2 * d * c₀) = ∫ _v in -d..d, 2 * d * c₀ := by
            rw [intervalIntegral.integral_const, smul_eq_mul]; ring
      _ ≤ ∫ v in -d..d, g v :=
          intervalIntegral.integral_mono_on (by linarith) intervalIntegrable_const
            (hgc.intervalIntegrable _ _) fun v hv => hgv v (abs_le.mpr hv)
      _ ≤ perpPair α β d :=
          intervalIntegral.integral_mono_interval (by linarith) (by linarith) hd1
            (Filter.Eventually.of_forall hg0) (hgc.intervalIntegrable _ _)
  -- The constant: 4 d² (d²)^{α/2} (3d²)^{−β/2} = 4 · 3^{−β/2} d^{2+α−β}.
  have e1 : (d ^ 2) ^ (α / 2) = d ^ α := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hd.le]; congr 1; push_cast; ring
  have e2 : (3 * d ^ 2) ^ (-β / 2) = (3 : ℝ) ^ (-β / 2) * d ^ (-β) := by
    rw [Real.mul_rpow (by norm_num) (by positivity), ← Real.rpow_natCast,
      ← Real.rpow_mul hd.le]
    congr 2; push_cast; ring
  have e3 : d ^ (2 + α - β) = d ^ 2 * d ^ α * d ^ (-β) := by
    rw [show 2 + α - β = ((2 : ℕ) : ℝ) + α + -β by push_cast; ring, Real.rpow_add hd,
      Real.rpow_add hd, Real.rpow_natCast]
  have : 4 * (3 : ℝ) ^ (-β / 2) * d ^ (2 + α - β) = 2 * d * (2 * d * c₀) := by
    rw [hc₀, e1, e2, e3]; ring
  rw [this]
  exact hmain

end RC

end DDG
