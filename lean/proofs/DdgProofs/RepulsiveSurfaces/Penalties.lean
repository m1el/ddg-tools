import DdgProofs.RepulsiveSurfaces.Constraints
import DDG.Geometry.Vec3

/-!
# Repulsive Surfaces: penalty values and gradients (paper Sec. 7.3; SPEC §5.2, D21, D24)

Every penalty's gradient is the exact derivative of its printed value (D24), along any
differentiable motion `Γ`:

* `rs-penalties`:
  - area / volume deviation `(X/X₀ − 1)²`: derivative `2(X/X₀ − 1)/X₀ · X'` with `X' = ∇A·Γ'`
    (`rs-area-gradient`) or `∇V·Γ'` (`rs-volume-gradient`);
  - total area `A^p` (real `p`, `A > 0`) and signed total volume `V^k` (natural `k`);
  - squared error `Σ_v |f_v − f⁰_v|²`: gradient `2(f_v − f⁰_v)`;
  - boundary length `(L − Σ_{e∈∂M} ℓ_e)²`: derivative `−2(L − L_∂) Σ_e T_e·(Γ'_{e₂} − Γ'_{e₁})`;
  - boundary turning angle `θ_v²/ℓ_v` with `θ_v = arccos(T_in·T_out)` (for `0 < θ_v < π`) and
    `ℓ_v = (ℓ_in + ℓ_out)/2`, and `θ_v²` alone (`dual_length: false`);
  - implicit obstacle `Σ_v q_v d(f_v)^{−p}` (`d > 0`) and attractor `Σ_v q_v |d(f_v)|^p` (`p > 1`,
    also at `d = 0`), with the quadrature weights `q_v` moving too (D21: `q_v = a_v`).
* `rs-willmore` (Sec. 7.3.5): `E_W = fᵀ W M⁻¹ W f` restricted to the rows `I` (interior vertices)
  is `Σ_{v∈I} |y_v|²/m_v` with `y = W f`; along a motion where `W` and `m` vary too,
  `E_W' = Σ_{v∈I} (2 y_v·y_v'/m_v − m_v' |y_v|²/m_v²)`, `y_v' = Σ_w (W'_vw f_w + W_vw Γ'_w)`
  (`hasDerivAt_willmore`, the structure; the cotan and mass derivatives themselves are checked by
  finite differences, SPEC §13.2); with `W`, `M` frozen the gradient is `2 Wᵀ M_I⁻¹ W f`.
-/

noncomputable section

namespace DDG

namespace RS

open Matrix Real Finset
open RC (V3 vnorm vnorm_nonneg vnorm_pos_iff hasDerivAt_dot hasDerivAt_len hasDerivAt_tangent
  tangent len sgn)

/-! ## Scalar chain rules -/

section scalar

variable {X : ℝ → ℝ} {X' t : ℝ}

/-- **Deviation penalty** `(X/X₀ − 1)²`: derivative `2(X/X₀ − 1)/X₀ · X'`.

[formula: rs-penalties] -/
theorem hasDerivAt_deviation (hX : HasDerivAt X X' t) (X₀ : ℝ) :
    HasDerivAt (fun s => (X s / X₀ - 1) ^ 2) (2 * (X t / X₀ - 1) / X₀ * X') t := by
  have := ((hX.div_const X₀).sub_const 1).pow 2
  refine this.congr_deriv ?_
  push_cast; ring

/-- **Power penalty** `X^p` for real `p` and `X > 0`: derivative `p X^{p−1} X'`.

[formula: rs-penalties] -/
theorem hasDerivAt_power (hX : HasDerivAt X X' t) (hpos : 0 < X t) (p : ℝ) :
    HasDerivAt (fun s => X s ^ p) (p * X t ^ (p - 1) * X') t := by
  have := hX.rpow_const (p := p) (Or.inl hpos.ne')
  refine this.congr_deriv ?_
  ring

/-- **Signed power penalty** `X^k` for natural `k` (signed volume): derivative `k X^{k−1} X'`.

[formula: rs-penalties] -/
theorem hasDerivAt_npower (hX : HasDerivAt X X' t) (k : ℕ) :
    HasDerivAt (fun s => X s ^ k) (k * X t ^ (k - 1) * X') t := by
  have := hX.pow k
  refine this.congr_deriv ?_
  ring

end scalar

/-! ## Mesh penalties -/

section mesh

variable {V F : Type*} [Fintype V] [Fintype F] [DecidableEq V]
variable {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ} {fv : F → Fin 3 → V}

/-- **Area deviation** `(A/A₀ − 1)²` (Sec. 7.3.1): gradient `2(A/A₀ − 1)/A₀ ∇A`.

[formula: rs-penalties] -/
theorem hasDerivAt_areaDeviation (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (hnd : ∀ S, Nondeg (Γ t) fv S) (A₀ : ℝ) :
    HasDerivAt (fun s => (totalArea (Γ s) fv / A₀ - 1) ^ 2)
      (∑ v, ((2 * (totalArea (Γ t) fv / A₀ - 1) / A₀) • areaGrad (Γ t) fv v) ⬝ᵥ Γ' v) t := by
  refine (hasDerivAt_deviation (hasDerivAt_totalArea hΓ hnd) A₀).congr_deriv ?_
  simp only [smul_dotProduct, smul_eq_mul, Finset.mul_sum]

/-- **Volume deviation** `(V/V₀ − 1)²` (Sec. 7.3.1): gradient `2(V/V₀ − 1)/V₀ ∇V` (the exact
gradient; the reference's is 3×, D9).

[formula: rs-penalties] -/
theorem hasDerivAt_volumeDeviation (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) (V₀ : ℝ) :
    HasDerivAt (fun s => (volume (Γ s) fv / V₀ - 1) ^ 2)
      (∑ v, ((2 * (volume (Γ t) fv / V₀ - 1) / V₀) • volGrad (Γ t) fv v) ⬝ᵥ Γ' v) t := by
  refine (hasDerivAt_deviation (hasDerivAt_volume hΓ) V₀).congr_deriv ?_
  simp only [smul_dotProduct, smul_eq_mul, Finset.mul_sum]

/-- **Total area** `A^p` (`p = 1`: minimal surfaces, Fig. 9; `p = 2`: the reference's
`minimize area`): gradient `p A^{p−1} ∇A`.

[formula: rs-penalties] -/
theorem hasDerivAt_totalAreaPow (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (hnd : ∀ S, Nondeg (Γ t) fv S) (hA : 0 < totalArea (Γ t) fv) (p : ℝ) :
    HasDerivAt (fun s => totalArea (Γ s) fv ^ p)
      (∑ v, ((p * totalArea (Γ t) fv ^ (p - 1)) • areaGrad (Γ t) fv v) ⬝ᵥ Γ' v) t := by
  refine (hasDerivAt_power (hasDerivAt_totalArea hΓ hnd) hA p).congr_deriv ?_
  simp only [smul_dotProduct, smul_eq_mul, Finset.mul_sum]

/-- **Total signed volume** `V^k`: gradient `k V^{k−1} ∇V`.

[formula: rs-penalties] -/
theorem hasDerivAt_totalVolumePow (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) (k : ℕ) :
    HasDerivAt (fun s => volume (Γ s) fv ^ k)
      (∑ v, (((k : ℝ) * volume (Γ t) fv ^ (k - 1)) • volGrad (Γ t) fv v) ⬝ᵥ Γ' v) t := by
  refine (hasDerivAt_npower (hasDerivAt_volume hΓ) k).congr_deriv ?_
  simp only [smul_dotProduct, smul_eq_mul, Finset.mul_sum]

omit [DecidableEq V] in
/-- **Squared error** `Σ_v |f_v − f⁰_v|²` (Fig. 21–22): gradient `2(f_v − f⁰_v)`.

[formula: rs-penalties] -/
theorem hasDerivAt_squaredError (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) (f₀ : V → V3) :
    HasDerivAt (fun s => ∑ v, (Γ s v - f₀ v) ⬝ᵥ (Γ s v - f₀ v))
      (∑ v, ((2 : ℝ) • (Γ t v - f₀ v)) ⬝ᵥ Γ' v) t := by
  have h := HasDerivAt.fun_sum (u := Finset.univ) fun v _ =>
    hasDerivAt_dot ((hΓ v).sub_const (f₀ v)) ((hΓ v).sub_const (f₀ v))
  refine h.congr_deriv (Finset.sum_congr rfl fun v _ => ?_)
  rw [smul_dotProduct, dotProduct_comm (Γ' v)]; simp only [smul_eq_mul]; ring

end mesh

/-! ## Boundary penalties -/

section boundary

variable {V E : Type*} [Fintype E] {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
  {ep : E → Fin 2 → V}

/-- **Boundary length** `(L − Σ_{e∈∂M} ℓ_e)²` (Sec. 7.3.4) over the boundary edges `E` with
endpoints `ep`: derivative `−2(L − L_∂) Σ_e Σ_a s_a T_e·Γ'_{e_a}`.

[formula: rs-penalties] -/
theorem hasDerivAt_boundaryLength (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (hpos : ∀ e, 0 < len (Γ t) ep e) (L : ℝ) :
    HasDerivAt (fun s => (L - ∑ e, len (Γ s) ep e) ^ 2)
      (-2 * (L - ∑ e, len (Γ t) ep e) *
        ∑ e, ∑ a, sgn a * (tangent (Γ t) ep e ⬝ᵥ Γ' (ep e a))) t := by
  have h := HasDerivAt.fun_sum (u := Finset.univ) fun e _ => hasDerivAt_len (ep := ep) hΓ (hpos e)
  have := ((h.const_sub L).pow 2)
  refine this.congr_deriv ?_
  push_cast; ring

/-- **Turning-angle chain rule**: for `c(s) ∈ (−1, 1)` and `ℓ(s) > 0`,
`d/ds arccos(c)²/ℓ = (2 θ θ' ℓ − θ² ℓ')/ℓ²` with `θ = arccos c`, `θ' = −c'/√(1 − c²)`; and
`d/ds arccos(c)² = 2θθ'`.

[formula: rs-penalties] -/
theorem hasDerivAt_turning {c ℓ : ℝ → ℝ} {c' ℓ' : ℝ} (hc : HasDerivAt c c' t)
    (hℓ : HasDerivAt ℓ ℓ' t) (h1 : -1 < c t) (h2 : c t < 1) (hl : 0 < ℓ t) :
    HasDerivAt (fun s => arccos (c s) ^ 2 / ℓ s)
      ((2 * arccos (c t) * (-c' / √(1 - c t ^ 2)) * ℓ t - arccos (c t) ^ 2 * ℓ') / ℓ t ^ 2) t ∧
    HasDerivAt (fun s => arccos (c s) ^ 2) (2 * arccos (c t) * (-c' / √(1 - c t ^ 2))) t := by
  have hθ : HasDerivAt (fun s => arccos (c s)) (-c' / √(1 - c t ^ 2)) t := by
    have := (Real.hasDerivAt_arccos h1.ne' h2.ne).comp t hc
    refine this.congr_deriv ?_
    ring
  have hθ2 := hθ.pow 2
  constructor
  · have := hθ2.div hℓ hl.ne'
    refine this.congr_deriv ?_
    simp only [Pi.pow_apply]; push_cast; ring
  · refine hθ2.congr_deriv ?_
    push_cast; ring

omit [Fintype E] in
/-- **Boundary turning angle at a vertex** (Sec. 7.3.4): with the incoming boundary edge `I` and
the outgoing one `J` (positive lengths, `T_I·T_J ∈ (−1, 1)`), `θ = arccos(T_I·T_J)` and
`ℓ_v = (ℓ_I + ℓ_J)/2`, `d/ds θ²/ℓ_v` is given by `hasDerivAt_turning` with
`c' = T_I'·T_J + T_I·T_J'` (`T'` of `RC.hasDerivAt_tangent`) and `ℓ_v' = (ℓ_I' + ℓ_J')/2`.

[formula: rs-penalties] -/
theorem hasDerivAt_turningPenalty (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) {I J : E}
    (hI : 0 < len (Γ t) ep I) (hJ : 0 < len (Γ t) ep J)
    (h1 : -1 < tangent (Γ t) ep I ⬝ᵥ tangent (Γ t) ep J)
    (h2 : tangent (Γ t) ep I ⬝ᵥ tangent (Γ t) ep J < 1) :
    let TI' := (len (Γ t) ep I)⁻¹ • ∑ a, sgn a • RC.proj (tangent (Γ t) ep I) (Γ' (ep I a))
    let TJ' := (len (Γ t) ep J)⁻¹ • ∑ a, sgn a • RC.proj (tangent (Γ t) ep J) (Γ' (ep J a))
    let c' := TI' ⬝ᵥ tangent (Γ t) ep J + tangent (Γ t) ep I ⬝ᵥ TJ'
    let ℓ' := ((∑ a, sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a))) +
      ∑ a, sgn a * (tangent (Γ t) ep J ⬝ᵥ Γ' (ep J a))) / 2
    let c := tangent (Γ t) ep I ⬝ᵥ tangent (Γ t) ep J
    let ℓv := (len (Γ t) ep I + len (Γ t) ep J) / 2
    HasDerivAt (fun s => arccos (tangent (Γ s) ep I ⬝ᵥ tangent (Γ s) ep J) ^ 2 /
        ((len (Γ s) ep I + len (Γ s) ep J) / 2))
      ((2 * arccos c * (-c' / √(1 - c ^ 2)) * ℓv - arccos c ^ 2 * ℓ') / ℓv ^ 2) t := by
  intro TI' TJ' c' ℓ' c ℓv
  have hc := hasDerivAt_dot (hasDerivAt_tangent (ep := ep) (I := I) hΓ hI)
    (hasDerivAt_tangent (ep := ep) (I := J) hΓ hJ)
  have hl := ((hasDerivAt_len (ep := ep) (I := I) hΓ hI).add
    (hasDerivAt_len (ep := ep) (I := J) hΓ hJ)).div_const 2
  exact (hasDerivAt_turning hc hl h1 h2 (half_pos (add_pos hI hJ))).1

end boundary

/-! ## Implicit potentials (Sec. 7.3.3, D21) -/

section implicit

variable {V : Type*} [Fintype V] {t : ℝ}

/-- **Implicit obstacle** `Σ_v q_v d_v^{−p}` (`d_v = d(f_v) > 0`, `q_v` the quadrature weight,
moving with the mesh for `q_v = a_v`): derivative `Σ_v (q_v' d_v^{−p} − p q_v d_v^{−p−1} d_v')`.

[formula: rs-penalties] -/
theorem hasDerivAt_implicitObstacle {q d : ℝ → V → ℝ} {q' d' : V → ℝ}
    (hq : ∀ v, HasDerivAt (fun s => q s v) (q' v) t) (hd : ∀ v, HasDerivAt (fun s => d s v) (d' v) t)
    (hpos : ∀ v, 0 < d t v) (p : ℝ) :
    HasDerivAt (fun s => ∑ v, q s v * d s v ^ (-p))
      (∑ v, (q' v * d t v ^ (-p) + q t v * (-p * d t v ^ (-p - 1) * d' v))) t := by
  refine HasDerivAt.fun_sum fun v _ => ?_
  have := (hq v).mul ((hd v).rpow_const (p := -p) (Or.inl (hpos v).ne'))
  refine this.congr_deriv ?_
  ring

/-- **Implicit attractor** `Σ_v q_v |d_v|^p` for `p > 1` (also at `d_v = 0`): derivative
`Σ_v (q_v' |d_v|^p + p q_v |d_v|^{p−2} d_v d_v')`.

[formula: rs-penalties] -/
theorem hasDerivAt_implicitAttractor {q d : ℝ → V → ℝ} {q' d' : V → ℝ}
    (hq : ∀ v, HasDerivAt (fun s => q s v) (q' v) t) (hd : ∀ v, HasDerivAt (fun s => d s v) (d' v) t)
    {p : ℝ} (hp : 1 < p) :
    HasDerivAt (fun s => ∑ v, q s v * |d s v| ^ p)
      (∑ v, (q' v * |d t v| ^ p + q t v * (p * |d t v| ^ (p - 2) * d t v * d' v))) t := by
  refine HasDerivAt.fun_sum fun v _ => ?_
  have := (hq v).mul ((hasDerivAt_abs_rpow (d t v) hp).comp t (hd v))
  exact this

end implicit

/-! ## Willmore (Sec. 7.3.5) -/

section willmore

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- `y_v = (W f)_v = Σ_w W_vw f_w`. -/
def wApply (W : V → V → ℝ) (f : V → V3) (v : V) : V3 := ∑ w, W v w • f w

/-- The Willmore term `fᵀ W M_I⁻¹ W f = Σ_{v∈I} |(W f)_v|²/m_v` over the rows `I` (interior
vertices; boundary rows dropped as the reference).

[formula: rs-willmore] -/
def willmore (I : Finset V) (W : V → V → ℝ) (m : V → ℝ) (f : V → V3) : ℝ :=
  ∑ v ∈ I, (wApply W f v ⬝ᵥ wApply W f v) / m v

/-- **The Willmore term is `fᵀ W M⁻¹ W f`** on the interior rows (componentwise:
`Σ_c (W f_c)ᵀ diag(1_I/m) (W f_c)`).

[formula: rs-willmore] -/
theorem willmore_eq_matrix (I : Finset V) (W : Matrix V V ℝ) (m : V → ℝ) (f : V → V3) :
    willmore I W m f = ∑ c : Fin 3, (W *ᵥ fun v => f v c) ⬝ᵥ
      (diagonal fun v => if v ∈ I then (m v)⁻¹ else 0) *ᵥ (W *ᵥ fun v => f v c) := by
  have hD : ∀ x : V → ℝ, (diagonal fun v => if v ∈ I then (m v)⁻¹ else 0) *ᵥ x =
      fun v => (if v ∈ I then (m v)⁻¹ else 0) * x v := fun x => funext fun v => mulVec_diagonal _ _ v
  have hI : willmore I W m f =
      ∑ v, if v ∈ I then (wApply W f v ⬝ᵥ wApply W f v) / m v else 0 := by
    unfold willmore; rw [Finset.sum_ite_mem, Finset.univ_inter]
  have hy : ∀ v c, wApply W f v c = (W *ᵥ fun v => f v c) v := by
    intro v c; simp [wApply, mulVec, dotProduct, Finset.sum_apply]
  rw [hI]
  simp only [hD, dotProduct]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun v _ => ?_
  simp only [Fin.sum_univ_three, hy]
  split_ifs <;> ring

omit [DecidableEq V] in
/-- **The derivative of the Willmore term** (structure, SPEC §5.2): along a motion where the
positions, the stiffness `W` and the masses `m > 0` all move,
`E_W' = Σ_{v∈I} (2 y_v·y_v'/m_v − m_v'|y_v|²/m_v²)` with `y_v' = Σ_w (W'_vw f_w + W_vw Γ'_w)`.

[formula: rs-willmore] -/
theorem hasDerivAt_willmore (I : Finset V) {Wc : ℝ → V → V → ℝ} {W' : V → V → ℝ}
    {m : ℝ → V → ℝ} {m' : V → ℝ} {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
    (hW : ∀ v w, HasDerivAt (fun s => Wc s v w) (W' v w) t)
    (hm : ∀ v, HasDerivAt (fun s => m s v) (m' v) t) (hmpos : ∀ v, 0 < m t v)
    (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) :
    HasDerivAt (fun s => willmore I (Wc s) (m s) (Γ s))
      (∑ v ∈ I, (2 * (wApply (Wc t) (Γ t) v ⬝ᵥ
          ∑ w, (W' v w • Γ t w + Wc t v w • Γ' w)) / m t v -
        m' v * (wApply (Wc t) (Γ t) v ⬝ᵥ wApply (Wc t) (Γ t) v) / m t v ^ 2)) t := by
  refine HasDerivAt.fun_sum fun v _ => ?_
  have hy : HasDerivAt (fun s => wApply (Wc s) (Γ s) v)
      (∑ w, (W' v w • Γ t w + Wc t v w • Γ' w)) t :=
    HasDerivAt.fun_sum fun w _ => ((hW v w).smul (hΓ w)).congr_deriv (add_comm _ _)
  have := (hasDerivAt_dot hy hy).div (hm v) (hmpos v).ne'
  refine this.congr_deriv ?_
  rw [dotProduct_comm (∑ w, _) (wApply (Wc t) (Γ t) v)]
  field_simp
  ring

omit [DecidableEq V] in
/-- **Frozen `W` and `M`**: the gradient of `E_W` is `g_w = 2 Σ_{v∈I} W_vw y_v/m_v`, i.e.
`2 Wᵀ M_I⁻¹ W f` (`= 2 W M_I⁻¹ W f` for the symmetric cotan matrix).

[formula: rs-willmore] -/
theorem hasDerivAt_willmore_frozen (I : Finset V) (W : V → V → ℝ) {m : V → ℝ}
    (hmpos : ∀ v, 0 < m v) {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
    (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) :
    HasDerivAt (fun s => willmore I W m (Γ s))
      (∑ w, (∑ v ∈ I, (2 * W v w / m v) • wApply W (Γ t) v) ⬝ᵥ Γ' w) t := by
  have h := hasDerivAt_willmore (Wc := fun _ => W) (m := fun _ => m) (W' := fun _ _ => 0)
    (m' := fun _ => 0) I (fun v w => hasDerivAt_const _ _) (fun v => hasDerivAt_const _ _)
    hmpos hΓ
  refine h.congr_deriv ?_
  simp only [zero_smul, zero_add, zero_mul, zero_div, sub_zero, dotProduct_sum, sum_dotProduct,
    dotProduct_smul, smul_dotProduct, smul_eq_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun v _ => ?_
  rw [div_eq_mul_inv, Finset.mul_sum, Finset.sum_mul]
  refine Finset.sum_congr rfl fun w _ => ?_
  ring

end willmore

end RS

end DDG
