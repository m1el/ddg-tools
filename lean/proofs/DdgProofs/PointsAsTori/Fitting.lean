import DdgProofs.PointsAsTori.Basic

/-!
# Points as Tori: torus fitting from the six coefficients (Sec. 4.1, App. C; SPEC §7.1)

The Monge patch at `p` over the frame `(ŝ, t̂, n)` is
`Q*(s, t) = p + s ŝ + t t̂ + Q(s, t) n`, `Q = a00 + a10 s + a01 t + a11 st + a20 s² + a02 t²`
(Eq. 18). With `S = 1 + a10² + a01²`, `A = √S`:

```
I  = [[1 + a10², a10 a01], [a10 a01, 1 + a01²]],   II = (1/A) [[2a20, a11], [a11, 2a02]]    (App. C)
H  = (a02(1 + a10²) + a20(1 + a01²) − a11 a10 a01)/A³,   K = (4 a02 a20 − a11²)/A⁴        (Eq. 19)
v1 = a20 a01² − a02 a10² − a02 + a20,   v2 = a11(1 + a10²) − 2 a01 a10 a20
term = A³ H,   Δ = S(a11² − 4 a02 a20) + term² = A⁶ (H² − K) ≥ 0
κ± = (term ± √Δ)/A³ = H ± √(H² − K),   ev± = (v1 ± √Δ, v2)                 (N§6.1)
n* = (n − a10 ŝ − a01 t̂)/A                                                  (Eq. 21)
```

Proved here: the partial derivatives and the fundamental forms of `Q*` at `0`
(`fundamental_forms`), `H = tr(I⁻¹II)/2` and `K = det II/det I` (`mean_gauss`), the eigenvector
pair and Eq. 20 (`principal_pair`, `principal_pair_eq20`), the Rayleigh quotients
(`rayleigh`), the 3D directions (`principal_directions`), the normal (`surface_normal`), the radii
(`torus_radii`), the sign rule (`torus_sign`), the center and axis (`torus_frame`), and the scale
covariance of the whole fit (`fit_scale`, the output scaling of SPEC §4.2).

**Degenerate eigenvectors** (`principal_pair_degenerate`): `(v1 ± √Δ, v2)` can be the zero vector
for one sign away from umbilics (`v2 = 0`, `v1 ≠ 0`); the identities then hold trivially and the
Rayleigh quotient is undefined. The fit below uses the exact formulas without the 1e-10 clamps
of D9; the clamps are absolute thresholds and break exact scale covariance when active.
-/

noncomputable section

namespace DDG

namespace PAT

open Matrix Real
open RC (V3 vnorm vnorm_nonneg vnorm_sq vnorm_mul_self vnorm_eq_zero_iff vnorm_pos_iff vnorm_smul
  dot_self_nonneg)

/-- The six coefficients of the height function `Q` (Eq. 18; `a_{nm}` multiplies `sⁿ tᵐ`). -/
structure Coef where
  a00 : ℝ
  a10 : ℝ
  a01 : ℝ
  a11 : ℝ
  a20 : ℝ
  a02 : ℝ

namespace Coef

variable (a : Coef)

/-- The height function `Q(s, t)`. -/
def Q (s t : ℝ) : ℝ :=
  a.a00 + a.a10 * s + a.a01 * t + a.a11 * s * t + a.a20 * s ^ 2 + a.a02 * t ^ 2

def S : ℝ := 1 + a.a10 ^ 2 + a.a01 ^ 2
def A : ℝ := √a.S
def E : ℝ := 1 + a.a10 ^ 2
def F : ℝ := a.a10 * a.a01
def G : ℝ := 1 + a.a01 ^ 2
/-- Second fundamental form entries `e, f, g` (App. C). -/
def e : ℝ := 2 * a.a20 / a.A
def f : ℝ := a.a11 / a.A
def g : ℝ := 2 * a.a02 / a.A
/-- Mean curvature, Eq. 19. -/
def H : ℝ := (a.a02 * (1 + a.a10 ^ 2) + a.a20 * (1 + a.a01 ^ 2) - a.a11 * a.a10 * a.a01) / a.A ^ 3
/-- Gauss curvature, Eq. 19. -/
def K : ℝ := (4 * a.a02 * a.a20 - a.a11 ^ 2) / a.A ^ 4
def v1 : ℝ := a.a20 * a.a01 ^ 2 - a.a02 * a.a10 ^ 2 - a.a02 + a.a20
def v2 : ℝ := a.a11 * (a.a10 ^ 2 + 1) - 2 * a.a01 * a.a10 * a.a20
def term : ℝ := a.a02 * (a.a10 ^ 2 + 1) - a.a01 * a.a10 * a.a11 + a.a20 * (a.a01 ^ 2 + 1)
/-- The (unclamped) squared discriminant `Δ`; the code takes `disc = √max(Δ, 1e-10)`. -/
def Δ : ℝ := (1 + a.a01 ^ 2 + a.a10 ^ 2) * (a.a11 ^ 2 - 4 * a.a02 * a.a20) + a.term ^ 2
/-- `κ₊ = (term + √Δ)/A³`. -/
def kp : ℝ := (a.term + √a.Δ) / a.A ^ 3
/-- `κ₋ = (term − √Δ)/A³`. -/
def km : ℝ := (a.term - √a.Δ) / a.A ^ 3
/-- `ev₊ = (v1 + √Δ, v2)`. -/
def evp : Fin 2 → ℝ := ![a.v1 + √a.Δ, a.v2]
/-- `ev₋ = (v1 − √Δ, v2)`. -/
def evm : Fin 2 → ℝ := ![a.v1 - √a.Δ, a.v2]

theorem S_pos : 0 < a.S := by unfold S; positivity
theorem A_pos : 0 < a.A := Real.sqrt_pos.mpr a.S_pos
theorem A_sq : a.A ^ 2 = a.S := Real.sq_sqrt a.S_pos.le
theorem A_ne : a.A ≠ 0 := a.A_pos.ne'

/-- `Δ ≥ 0`: `E²Δ = (S a20 − F² a20 + EF a11 − E² a02)² + S (E a11 − 2F a20)²`. -/
theorem Δ_nonneg : 0 ≤ a.Δ := by
  have hE : 0 < a.E := by unfold E; positivity
  have key : a.E ^ 2 * a.Δ = (a.S * a.a20 - a.F ^ 2 * a.a20 + a.E * a.F * a.a11 - a.E ^ 2 * a.a02) ^ 2
      + a.S * (a.E * a.a11 - 2 * a.F * a.a20) ^ 2 := by
    unfold Δ term E F S; ring
  have h2 : 0 ≤ a.E ^ 2 * a.Δ := by rw [key]; have := a.S_pos; positivity
  by_contra h
  push Not at h
  have := mul_neg_of_pos_of_neg (pow_pos hE 2) h
  linarith

end Coef

/-! ## The Monge patch and its fundamental forms -/

/-- The Monge patch `Q*(s, t) = p + s ŝ + t t̂ + Q(s, t) n` (Eq. 18). -/
def patch (Fr : Frame) (p : V3) (a : Coef) (s t : ℝ) : V3 := p + Fr.comb s t (a.Q s t)

/-- The unit normal of Eq. 21, `n* = (n − a10 ŝ − a01 t̂)/A`. -/
def nstar (Fr : Frame) (a : Coef) : V3 := a.A⁻¹ • (Fr.n - a.a10 • Fr.s - a.a01 • Fr.t)

theorem nstar_comb (Fr : Frame) (a : Coef) :
    nstar Fr a = Fr.comb (-a.a10 / a.A) (-a.a01 / a.A) (1 / a.A) := by
  ext i; simp only [nstar, Frame.comb, Pi.smul_apply, Pi.sub_apply, Pi.add_apply, smul_eq_mul]
  field_simp; try ring

theorem hasDerivAt_comb (Fr : Frame) {α β γ : ℝ → ℝ} {α' β' γ' τ : ℝ} (hα : HasDerivAt α α' τ)
    (hβ : HasDerivAt β β' τ) (hγ : HasDerivAt γ γ' τ) :
    HasDerivAt (fun x => Fr.comb (α x) (β x) (γ x)) (Fr.comb α' β' γ') τ :=
  ((hα.smul_const Fr.s).add (hβ.smul_const Fr.t)).add (hγ.smul_const Fr.n)

theorem hasDerivAt_quad (c0 c1 c2 x : ℝ) :
    HasDerivAt (fun y => c0 + c1 * y + c2 * y ^ 2) (c1 + 2 * c2 * x) x := by
  have h1 : HasDerivAt (fun y => c1 * y) c1 x := by simpa using (hasDerivAt_id x).const_mul c1
  have h2 : HasDerivAt (fun y => c2 * y ^ 2) (c2 * (2 * x)) x := by
    simpa using (hasDerivAt_pow 2 x).const_mul c2
  have := (h1.const_add c0).add h2
  convert this using 1; ring

theorem hasDerivAt_lin (c0 c1 x : ℝ) : HasDerivAt (fun y => c0 + c1 * y) c1 x := by
  simpa using ((hasDerivAt_id x).const_mul c1).const_add c0

/-- **The first and second fundamental forms of the Monge patch** (App. C). The partial
derivatives of `Q*` are `Q*_s = ŝ + Q_s n`, `Q*_t = t̂ + Q_t n`, the second partials are
`Q*_ss = 2a20 n`, `Q*_st = a11 n`, `Q*_tt = 2a02 n`, and at `(0, 0)`:
`E = 1 + a10²`, `F = a10 a01`, `G = 1 + a01²`, `e = n*·Q*_ss = 2a20/A`, `f = a11/A`,
`g = 2a02/A`.

[formula: pat-fundamental-forms] -/
theorem fundamental_forms (Fr : Frame) (p : V3) (a : Coef) (s t : ℝ) :
    HasDerivAt (fun x => patch Fr p a x t) (Fr.comb 1 0 (a.a10 + a.a11 * t + 2 * a.a20 * s)) s ∧
    HasDerivAt (fun y => patch Fr p a s y) (Fr.comb 0 1 (a.a01 + a.a11 * s + 2 * a.a02 * t)) t ∧
    HasDerivAt (fun x => Fr.comb 1 0 (a.a10 + a.a11 * t + 2 * a.a20 * x)) (Fr.comb 0 0 (2 * a.a20)) s ∧
    HasDerivAt (fun y => Fr.comb 1 0 (a.a10 + a.a11 * y + 2 * a.a20 * s)) (Fr.comb 0 0 a.a11) t ∧
    HasDerivAt (fun y => Fr.comb 0 1 (a.a01 + a.a11 * s + 2 * a.a02 * y)) (Fr.comb 0 0 (2 * a.a02)) t ∧
    (Fr.comb 1 0 a.a10 ⬝ᵥ Fr.comb 1 0 a.a10 = a.E ∧ Fr.comb 1 0 a.a10 ⬝ᵥ Fr.comb 0 1 a.a01 = a.F ∧
      Fr.comb 0 1 a.a01 ⬝ᵥ Fr.comb 0 1 a.a01 = a.G) ∧
    (nstar Fr a ⬝ᵥ Fr.comb 0 0 (2 * a.a20) = a.e ∧ nstar Fr a ⬝ᵥ Fr.comb 0 0 a.a11 = a.f ∧
      nstar Fr a ⬝ᵥ Fr.comb 0 0 (2 * a.a02) = a.g) := by
  have hc0 : ∀ x : ℝ, HasDerivAt (fun _ : ℝ => (0 : ℝ)) 0 x := fun x => hasDerivAt_const x 0
  have hc1 : ∀ x : ℝ, HasDerivAt (fun _ : ℝ => (1 : ℝ)) 0 x := fun x => hasDerivAt_const x 1
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · have hQ : (fun x => a.Q x t) = fun x => (a.a00 + a.a01 * t + a.a02 * t ^ 2) +
        (a.a10 + a.a11 * t) * x + a.a20 * x ^ 2 := by funext x; unfold Coef.Q; ring
    have h := hasDerivAt_comb Fr (hasDerivAt_id s) (hasDerivAt_const s t)
      (hQ ▸ hasDerivAt_quad (a.a00 + a.a01 * t + a.a02 * t ^ 2) (a.a10 + a.a11 * t) a.a20 s)
    have := h.const_add p
    convert this using 2 <;> (try simp [patch]) <;> (try (congr 1; ring))
  · have hQ : (fun y => a.Q s y) = fun y => (a.a00 + a.a10 * s + a.a20 * s ^ 2) +
        (a.a01 + a.a11 * s) * y + a.a02 * y ^ 2 := by funext y; unfold Coef.Q; ring
    have h := hasDerivAt_comb Fr (hasDerivAt_const t s) (hasDerivAt_id t)
      (hQ ▸ hasDerivAt_quad (a.a00 + a.a10 * s + a.a20 * s ^ 2) (a.a01 + a.a11 * s) a.a02 t)
    have := h.const_add p
    convert this using 2 <;> (try simp [patch]) <;> (try (congr 1; ring))
  · exact hasDerivAt_comb Fr (hc1 s) (hc0 s)
      ((hasDerivAt_lin (a.a10 + a.a11 * t) (2 * a.a20) s).congr_of_eventuallyEq
        (Filter.Eventually.of_forall fun y => by ring))
  · exact hasDerivAt_comb Fr (hc1 t) (hc0 t)
      ((hasDerivAt_lin (a.a10 + 2 * a.a20 * s) a.a11 t).congr_of_eventuallyEq
        (Filter.Eventually.of_forall fun y => by ring))
  · exact hasDerivAt_comb Fr (hc0 t) (hc1 t)
      ((hasDerivAt_lin (a.a01 + a.a11 * s) (2 * a.a02) t).congr_of_eventuallyEq
        (Filter.Eventually.of_forall fun y => by ring))
  · simp only [Fr.comb_dot, Coef.E, Coef.F, Coef.G]; refine ⟨by ring, by ring, by ring⟩
  · simp only [nstar_comb, Fr.comb_dot, Coef.e, Coef.f, Coef.g]
    refine ⟨by ring, by ring, by ring⟩

/-- **Eq. 21 is the unit normal of `Q*` at `0`**: `Q*_s × Q*_t = A n*`, `|Q*_s × Q*_t| = A`, so
`n* = (Q*_s × Q*_t)/|Q*_s × Q*_t|`; `n*` is a unit vector orthogonal to `Q*_s` and `Q*_t`.

[formula: pat-surface-normal] -/
theorem surface_normal (Fr : Frame) (a : Coef) :
    Fr.comb 1 0 a.a10 ⨯₃ Fr.comb 0 1 a.a01 = a.A • nstar Fr a ∧
    vnorm (Fr.comb 1 0 a.a10 ⨯₃ Fr.comb 0 1 a.a01) = a.A ∧
    nstar Fr a = (vnorm (Fr.comb 1 0 a.a10 ⨯₃ Fr.comb 0 1 a.a01))⁻¹ •
      (Fr.comb 1 0 a.a10 ⨯₃ Fr.comb 0 1 a.a01) ∧
    nstar Fr a ⬝ᵥ nstar Fr a = 1 ∧ nstar Fr a ⬝ᵥ Fr.comb 1 0 a.a10 = 0 ∧
      nstar Fr a ⬝ᵥ Fr.comb 0 1 a.a01 = 0 := by
  have hA := a.A_ne
  have hA2 := a.A_sq
  have hcross : Fr.comb 1 0 a.a10 ⨯₃ Fr.comb 0 1 a.a01 = a.A • nstar Fr a := by
    rw [Fr.comb_cross, nstar_comb, Fr.comb_smul]
    congr 1 <;> field_simp <;> ring
  have hnorm : vnorm (Fr.comb 1 0 a.a10 ⨯₃ Fr.comb 0 1 a.a01) = a.A := by
    rw [Fr.comb_cross, Fr.vnorm_comb]
    unfold Coef.A Coef.S; congr 1; ring
  refine ⟨hcross, hnorm, ?_, ?_, ?_, ?_⟩
  · rw [hnorm, hcross, smul_smul, inv_mul_cancel₀ hA, one_smul]
  · rw [nstar_comb, Fr.comb_dot]; field_simp; rw [hA2]; unfold Coef.S; ring
  · rw [nstar_comb, Fr.comb_dot]; field_simp; try ring
  · rw [nstar_comb, Fr.comb_dot]; field_simp; try ring

/-! ## Mean and Gauss curvature -/

/-- The first fundamental form at `0`. -/
def Coef.I1 (a : Coef) : Matrix (Fin 2) (Fin 2) ℝ := !![a.E, a.F; a.F, a.G]
/-- The second fundamental form at `0`. -/
def Coef.II (a : Coef) : Matrix (Fin 2) (Fin 2) ℝ := !![a.e, a.f; a.f, a.g]

theorem Coef.det_I1 (a : Coef) : a.I1.det = a.S := by
  simp [Coef.I1, Matrix.det_fin_two_of, Coef.E, Coef.F, Coef.G, Coef.S]; ring

/-- **`H = tr(I⁻¹II)/2` and `K = det II/det I` are Eq. 19**, and `H² − K = Δ/A⁶`, so the two
principal curvatures `H ± √(H² − K)` are `κ± = (term ± √Δ)/A³`.

[formula: pat-mean-gauss] -/
theorem mean_gauss (a : Coef) :
    (a.I1⁻¹ * a.II).trace / 2 = a.H ∧ a.II.det / a.I1.det = a.K ∧
    a.H ^ 2 - a.K = a.Δ / a.A ^ 6 ∧
    a.H + √(a.H ^ 2 - a.K) = a.kp ∧ a.H - √(a.H ^ 2 - a.K) = a.km := by
  have hA := a.A_ne
  have hA2 := a.A_sq
  have hS : a.S ≠ 0 := a.S_pos.ne'
  have hdet : a.I1.det ≠ 0 := by rw [a.det_I1]; exact hS
  have hinv : a.I1⁻¹ = a.S⁻¹ • !![a.G, -a.F; -a.F, a.E] := by
    rw [Matrix.inv_def, a.det_I1, Ring.inverse_eq_inv', Coef.I1, Matrix.adjugate_fin_two_of]
  have hHK : a.H ^ 2 - a.K = a.Δ / a.A ^ 6 := by
    have h4 : a.A ^ 4 = a.S ^ 2 := by rw [← hA2]; ring
    have h6 : a.A ^ 6 = a.S ^ 3 := by rw [← hA2]; ring
    have h6' : (a.A ^ 3) ^ 2 = a.S ^ 3 := by rw [← hA2]; ring
    unfold Coef.H Coef.K
    rw [div_pow, h6', h4, h6]
    field_simp
    unfold Coef.Δ Coef.term Coef.S; ring
  have hsq : √(a.H ^ 2 - a.K) = √a.Δ / a.A ^ 3 := by
    rw [hHK, Real.sqrt_div' _ (by positivity), show a.A ^ 6 = (a.A ^ 3) ^ 2 by ring,
      Real.sqrt_sq (by have := a.A_pos; positivity)]
  have hH : a.H = a.term / a.A ^ 3 := by unfold Coef.H Coef.term; ring_nf
  refine ⟨?_, ?_, hHK, ?_, ?_⟩
  · rw [hinv, Matrix.smul_mul, Matrix.trace_smul, Coef.II]
    simp only [Matrix.trace_fin_two, Matrix.mul_apply, Fin.sum_univ_two, Matrix.of_apply,
      Matrix.cons_val', Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.empty_val',
      Matrix.cons_val_fin_one, Matrix.head_cons, smul_eq_mul]
    unfold Coef.H Coef.e Coef.f Coef.g Coef.E Coef.F Coef.G
    have h3 : a.A ^ 3 = a.S * a.A := by rw [← hA2]; ring
    rw [h3]
    field_simp
    ring
  · rw [a.det_I1, Coef.II, Matrix.det_fin_two_of]
    unfold Coef.K Coef.e Coef.f Coef.g
    have h4 : a.A ^ 4 = a.S ^ 2 := by rw [← hA2]; ring
    rw [h4]
    field_simp
    rw [hA2]; ring
  · rw [hsq, hH, Coef.kp]; ring
  · rw [hsq, hH, Coef.km]; ring

/-! ## Principal curvatures and directions -/

/-- The quadratic (bilinear) form of a symmetric 2×2 matrix `[[m00, m01], [m01, m11]]`. -/
def bil (m00 m01 m11 : ℝ) (w w' : Fin 2 → ℝ) : ℝ :=
  w 0 * (m00 * w' 0 + m01 * w' 1) + w 1 * (m01 * w' 0 + m11 * w' 1)

/-- `w` solves `(II − κ I) w = 0` (both rows). -/
def Coef.IsEig (a : Coef) (κ : ℝ) (w : Fin 2 → ℝ) : Prop :=
  (a.e - κ * a.E) * w 0 + (a.f - κ * a.F) * w 1 = 0 ∧
    (a.f - κ * a.F) * w 0 + (a.g - κ * a.G) * w 1 = 0

/-- **The eigenvector pair** (N§6.1): for either square root `d` of `Δ` (`d = ±√Δ`),
`(II − κ I)(v1 + d, v2) = 0` with `κ = (term + d)/A³`. In particular `ev±` are generalized
eigenvectors of `(II, I)` for `κ±` (`principal_pair`). The identity is polynomial modulo
`d² = Δ`: `A³·row₁ = −E (d² − Δ)`, `A³·row₂ = −F (d² − Δ)`.

[formula: pat-principal-pair] -/
theorem principal_pair_root (a : Coef) {d : ℝ} (hd : d ^ 2 = a.Δ) :
    a.IsEig ((a.term + d) / a.A ^ 3) ![a.v1 + d, a.v2] := by
  have hA := a.A_ne
  have hA2 := a.A_sq
  have k1 : a.S * (2 * a.a20 * (a.v1 + d) + a.a11 * a.v2) -
      (a.term + d) * (a.E * (a.v1 + d) + a.F * a.v2) = 0 := by
    unfold Coef.Δ Coef.term at hd; unfold Coef.S Coef.v1 Coef.v2 Coef.E Coef.F Coef.term
    linear_combination (-(a.a10 ^ 2 + 1)) * hd
  have k2 : a.S * (a.a11 * (a.v1 + d) + 2 * a.a02 * a.v2) -
      (a.term + d) * (a.F * (a.v1 + d) + a.G * a.v2) = 0 := by
    unfold Coef.Δ Coef.term at hd; unfold Coef.S Coef.v1 Coef.v2 Coef.F Coef.G Coef.term
    linear_combination (-(a.a01 * a.a10)) * hd
  constructor
  · simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Coef.e, Coef.f]
    have : (2 * a.a20 / a.A - (a.term + d) / a.A ^ 3 * a.E) * (a.v1 + d) +
        (a.a11 / a.A - (a.term + d) / a.A ^ 3 * a.F) * a.v2 =
        (a.A ^ 2 * (2 * a.a20 * (a.v1 + d) + a.a11 * a.v2) -
          (a.term + d) * (a.E * (a.v1 + d) + a.F * a.v2)) / a.A ^ 3 := by
      field_simp; try ring
    rw [this, hA2, k1, zero_div]
  · simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Coef.f, Coef.g]
    have : (a.a11 / a.A - (a.term + d) / a.A ^ 3 * a.F) * (a.v1 + d) +
        (2 * a.a02 / a.A - (a.term + d) / a.A ^ 3 * a.G) * a.v2 =
        (a.A ^ 2 * (a.a11 * (a.v1 + d) + 2 * a.a02 * a.v2) -
          (a.term + d) * (a.F * (a.v1 + d) + a.G * a.v2)) / a.A ^ 3 := by
      field_simp; try ring
    rw [this, hA2, k2, zero_div]

/-- `ev₊`, `ev₋` solve `(II − κ± I) ev± = 0`.

[formula: pat-principal-pair] -/
theorem principal_pair (a : Coef) : a.IsEig a.kp a.evp ∧ a.IsEig a.km a.evm := by
  have hs := Real.sq_sqrt a.Δ_nonneg
  refine ⟨principal_pair_root a hs, ?_⟩
  have := principal_pair_root a (d := -√a.Δ) (by rw [neg_sq, hs])
  simpa [Coef.km, Coef.evm, sub_eq_add_neg] using this

/-- **Eq. 20**: `w = (κ a10 a01 A − a11, 2a20 − κ(1 + a10²)A)` solves `(II − κ I) w = 0` for
`κ = κ±` (the paper's form; it is the solution of the first row, the code's `ev±` of the
second). The characteristic polynomial vanishes: `A⁶ det(II − κI) = S (d² − Δ)`.

[formula: pat-principal-pair] -/
theorem principal_pair_eq20 (a : Coef) {d : ℝ} (hd : d ^ 2 = a.Δ) :
    a.IsEig ((a.term + d) / a.A ^ 3)
      ![(a.term + d) / a.A ^ 3 * a.a10 * a.a01 * a.A - a.a11,
        2 * a.a20 - (a.term + d) / a.A ^ 3 * (1 + a.a10 ^ 2) * a.A] := by
  have hA := a.A_ne
  have hA2 := a.A_sq
  set κ := (a.term + d) / a.A ^ 3 with hκ
  have hchar : (a.e - κ * a.E) * (a.g - κ * a.G) - (a.f - κ * a.F) ^ 2 = 0 := by
    have k : (2 * a.a20 * a.S - (a.term + d) * a.E) * (2 * a.a02 * a.S - (a.term + d) * a.G) -
        (a.a11 * a.S - (a.term + d) * a.F) ^ 2 = 0 := by
      unfold Coef.Δ at hd; unfold Coef.S Coef.E Coef.F Coef.G Coef.term
      unfold Coef.term at hd
      linear_combination (1 + a.a01 ^ 2 + a.a10 ^ 2) * hd
    have : (a.e - κ * a.E) * (a.g - κ * a.G) - (a.f - κ * a.F) ^ 2 =
        ((2 * a.a20 * a.A ^ 2 - (a.term + d) * a.E) * (2 * a.a02 * a.A ^ 2 - (a.term + d) * a.G) -
        (a.a11 * a.A ^ 2 - (a.term + d) * a.F) ^ 2) / a.A ^ 6 := by
      simp only [hκ, Coef.e, Coef.f, Coef.g]; field_simp; try ring
    rw [this, hA2, k, zero_div]
  constructor
  · simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Coef.e, Coef.f, Coef.E, Coef.F]
    field_simp; try ring
  · simp only [Matrix.cons_val_zero, Matrix.cons_val_one]
    have h1 : κ * a.a10 * a.a01 * a.A - a.a11 = -a.A * (a.f - κ * a.F) := by
      simp only [Coef.f, Coef.F]; field_simp; try ring
    have h2 : 2 * a.a20 - κ * (1 + a.a10 ^ 2) * a.A = a.A * (a.e - κ * a.E) := by
      simp only [Coef.e, Coef.E]; field_simp; try ring
    rw [h1, h2]
    linear_combination a.A * hchar

/-- **The pair can degenerate off the umbilics**: for `a10 = a01 = a11 = 0` and `a20 < a02`,
`ev₊ = (v1 + √Δ, v2) = 0` although `κ₊ ≠ κ₋` (symmetrically `ev₋ = 0` for `a20 > a02`).
The reference clamps the zero component to `±1e-10`, which turns `ev₊` into the eigenvector of
`κ₋`; both Rayleigh quotients then equal `κ₋`. -/
theorem principal_pair_degenerate (a : Coef) (h10 : a.a10 = 0) (h01 : a.a01 = 0) (h11 : a.a11 = 0)
    (hlt : a.a20 < a.a02) : a.evp = 0 ∧ a.kp ≠ a.km := by
  have hΔ : a.Δ = (a.a02 - a.a20) ^ 2 := by
    unfold Coef.Δ Coef.term; rw [h10, h01, h11]; ring
  have hs : √a.Δ = a.a02 - a.a20 := by rw [hΔ, Real.sqrt_sq (by linarith)]
  have hA : a.A = 1 := by unfold Coef.A Coef.S; rw [h10, h01]; norm_num
  refine ⟨?_, ?_⟩
  · ext i; fin_cases i
    · simp [Coef.evp, hs, Coef.v1, h10, h01]
    · simp [Coef.evp, Coef.v2, h10, h01, h11]
  · simp only [Coef.kp, Coef.km, hs, hA]; intro h; norm_num at h; linarith

/-- `ev±` are nonzero whenever `v2 ≠ 0`. -/
theorem principal_pair_ne_zero (a : Coef) (h : a.v2 ≠ 0) : a.evp ≠ 0 ∧ a.evm ≠ 0 := by
  constructor <;> intro h0 <;> apply h <;> simpa [Coef.evp, Coef.evm] using congrFun h0 1

/-- **Rayleigh quotient** (N§6.1): if `(II − κ I) w = 0` then `wᵀ II w = κ wᵀ I w`; for
`wᵀ I w ≠ 0` (any `w ≠ 0`, `I` being positive definite, `bil_I_pos`) the quotient is `κ`.
The code's curvatures are these quotients of the normalized `ev±` (homogeneous of degree 0).

[formula: pat-curvature-rayleigh] -/
theorem rayleigh (a : Coef) {κ : ℝ} {w : Fin 2 → ℝ} (h : a.IsEig κ w)
    (hI : bil a.E a.F a.G w w ≠ 0) :
    bil a.e a.f a.g w w = κ * bil a.E a.F a.G w w ∧
      bil a.e a.f a.g w w / bil a.E a.F a.G w w = κ := by
  have key : bil a.e a.f a.g w w = κ * bil a.E a.F a.G w w := by
    obtain ⟨h1, h2⟩ := h
    unfold bil
    linear_combination w 0 * h1 + w 1 * h2
  exact ⟨key, by rw [key]; field_simp⟩

/-! ## The 3D principal directions -/

/-- The 3D direction `v(w) = w_x (ŝ + a10 n) + w_y (t̂ + a01 n)` (Eq. 20 and N§6.1). -/
def vdir (Fr : Frame) (a : Coef) (w : Fin 2 → ℝ) : V3 :=
  w 0 • (Fr.s + a.a10 • Fr.n) + w 1 • (Fr.t + a.a01 • Fr.n)

theorem vdir_comb (Fr : Frame) (a : Coef) (w : Fin 2 → ℝ) :
    vdir Fr a w = Fr.comb (w 0) (w 1) (a.a10 * w 0 + a.a01 * w 1) := by
  ext i; simp only [vdir, Frame.comb, Pi.smul_apply, Pi.add_apply, smul_eq_mul]; ring

/-- `v(w)·v(w') = wᵀ I w'`: `I` is the Gram matrix of `Q*_s, Q*_t`. -/
theorem vdir_dot (Fr : Frame) (a : Coef) (w w' : Fin 2 → ℝ) :
    vdir Fr a w ⬝ᵥ vdir Fr a w' = bil a.E a.F a.G w w' := by
  rw [vdir_comb, vdir_comb, Fr.comb_dot]; unfold bil Coef.E Coef.F Coef.G; ring

/-- `I` is positive definite. -/
theorem bil_I_pos (a : Coef) {w : Fin 2 → ℝ} (hw : w ≠ 0) : 0 < bil a.E a.F a.G w w := by
  have : bil a.E a.F a.G w w = w 0 ^ 2 + w 1 ^ 2 + (a.a10 * w 0 + a.a01 * w 1) ^ 2 := by
    unfold bil Coef.E Coef.F Coef.G; ring
  rw [this]
  have : w 0 ≠ 0 ∨ w 1 ≠ 0 := by
    by_contra h; push Not at h; apply hw; ext i; fin_cases i <;> simp [h.1, h.2]
  rcases this with h | h <;> positivity

/-- **The principal directions** (Eq. 20, App. C): `v(w)` is tangent to `Q*` at `0`
(`v·n* = 0`), and generalized eigenvectors of distinct eigenvalues are `I`-orthogonal, so
`v₊ ⊥ v₋` in `ℝ³` when `κ₊ ≠ κ₋`.

[formula: pat-principal-directions] -/
theorem principal_directions (Fr : Frame) (a : Coef) {κ κ' : ℝ} {w w' : Fin 2 → ℝ}
    (h : a.IsEig κ w) (h' : a.IsEig κ' w') (hne : κ ≠ κ') :
    vdir Fr a w ⬝ᵥ nstar Fr a = 0 ∧ bil a.E a.F a.G w w' = 0 ∧
      vdir Fr a w ⬝ᵥ vdir Fr a w' = 0 := by
  have hA := a.A_ne
  have horth : bil a.E a.F a.G w w' = 0 := by
    obtain ⟨h1, h2⟩ := h
    obtain ⟨h1', h2'⟩ := h'
    have : (κ - κ') * bil a.E a.F a.G w w' = 0 := by
      unfold bil
      linear_combination -(w' 0 * h1 + w' 1 * h2) + (w 0 * h1' + w 1 * h2')
    rcases mul_eq_zero.mp this with h0 | h0
    · exact absurd (sub_eq_zero.mp h0) hne
    · exact h0
  refine ⟨?_, horth, by rw [vdir_dot]; exact horth⟩
  rw [vdir_comb, nstar_comb, Fr.comb_dot]; field_simp; try ring

/-! ## Torus radii, sign, center and axis -/

/-- The sign rule, Eq. 22. -/
def torusSign (kmin kmax : ℝ) : ℝ :=
  if (kmax < 0 ∧ kmin < 0) ∨ (kmax < 0 ∧ 0 < kmin) then 1 else -1

/-- **Eq. 22 is "`κmax < 0` and `κmin ≠ 0`"**.

[formula: pat-torus-sign] -/
theorem torus_sign (kmin kmax : ℝ) :
    torusSign kmin kmax = if kmax < 0 ∧ kmin ≠ 0 then 1 else -1 := by
  unfold torusSign
  congr 1
  apply propext
  constructor
  · rintro (⟨h1, h2⟩ | ⟨h1, h2⟩)
    · exact ⟨h1, h2.ne⟩
    · exact ⟨h1, h2.ne'⟩
  · rintro ⟨h1, h2⟩
    rcases lt_or_gt_of_ne h2 with h | h
    · exact Or.inl ⟨h1, h⟩
    · exact Or.inr ⟨h1, h⟩

/-- **The radii** `r = 1/|κmax|`, `R = 1/|κmin| − sgn(κ₊κ₋) r` satisfy
`(R + sgn(κ₊κ₋) r)⁻¹ = |κmin|` (`sgn` is `Real.sign`, `0` at `0`, as `jnp.sign`); the same holds
with the D9 floor `|κ| ↦ max(|κ|, ε)`. For `|κmin| ≤ |κmax|` and `κmin ≠ 0`, `R ≥ 0`; for
`κ₊κ₋ > 0`, `R < r` (a spindle torus) iff `|κmax| < 2|κmin|`.

[formula: pat-torus-radii] -/
theorem torus_radii (kmin kmax σ ε : ℝ) :
    ((1 / |kmin| - Real.sign σ * (1 / |kmax|)) + Real.sign σ * (1 / |kmax|))⁻¹ = |kmin| ∧
    ((1 / max |kmin| ε - Real.sign σ * (1 / max |kmax| ε)) +
      Real.sign σ * (1 / max |kmax| ε))⁻¹ = max |kmin| ε ∧
    (kmin ≠ 0 → |kmin| ≤ |kmax| → 0 ≤ 1 / |kmin| - Real.sign σ * (1 / |kmax|)) ∧
    (0 < σ → kmin ≠ 0 → kmax ≠ 0 →
      (1 / |kmin| - Real.sign σ * (1 / |kmax|) < 1 / |kmax| ↔ |kmax| < 2 * |kmin|)) := by
  refine ⟨by simp, by simp, ?_, ?_⟩
  · intro h0 hle
    have hp : 0 < |kmin| := abs_pos.mpr h0
    have hq : 0 < |kmax| := hp.trans_le hle
    have h1 : 1 / |kmax| ≤ 1 / |kmin| := one_div_le_one_div_of_le hp hle
    rcases lt_trichotomy σ 0 with h | h | h
    · rw [Real.sign_of_neg h]; have : 0 < 1 / |kmax| := by positivity
      linarith [show 0 < 1 / |kmin| by positivity]
    · rw [h, Real.sign_zero]; simp only [zero_mul, sub_zero]; positivity
    · rw [Real.sign_of_pos h]; linarith
  · intro hσ h0 h1
    rw [Real.sign_of_pos hσ]
    have hp : 0 < |kmin| := abs_pos.mpr h0
    have hq : 0 < |kmax| := abs_pos.mpr h1
    rw [one_mul, sub_lt_iff_lt_add, ← two_mul, show 2 * (1 / |kmax|) = 2 / |kmax| by ring,
      div_lt_div_iff₀ hp hq]
    constructor <;> intro h <;> linarith

/-- **The torus axis** `u = (n* × v)/|n* × v|` for a unit `n*` and a nonzero tangent `v ⊥ n*`:
`|n* × v| = |v|`, `u` is a unit vector, `u ⊥ n*` and `u ⊥ v` (with the code's guard
`max(|n* × v|, 1e-10)` the same vector once `|v| ≥ 1e-10`). The center is
`c = p + a00 n − (sign/|κmin|) n*` (definition `fitCenter`).

[formula: pat-torus-frame] -/
theorem torus_frame {ns v : V3} (hn : ns ⬝ᵥ ns = 1) (hv : v ≠ 0) (hnv : ns ⬝ᵥ v = 0) :
    vnorm (ns ⨯₃ v) = vnorm v ∧
    ((vnorm (ns ⨯₃ v))⁻¹ • (ns ⨯₃ v)) ⬝ᵥ ((vnorm (ns ⨯₃ v))⁻¹ • (ns ⨯₃ v)) = 1 ∧
    ((vnorm (ns ⨯₃ v))⁻¹ • (ns ⨯₃ v)) ⬝ᵥ ns = 0 ∧
    ((vnorm (ns ⨯₃ v))⁻¹ • (ns ⨯₃ v)) ⬝ᵥ v = 0 ∧
    (1e-10 ≤ vnorm v → max (vnorm (ns ⨯₃ v)) 1e-10 = vnorm (ns ⨯₃ v)) := by
  have hnorm : vnorm (ns ⨯₃ v) = vnorm v := by
    unfold vnorm; rw [cross_dot_cross, hn, hnv]; ring_nf
  have hne : ns ⨯₃ v ≠ 0 := by
    intro h; apply hv; rw [← vnorm_eq_zero_iff, ← hnorm, h, RC.vnorm_zero]
  refine ⟨hnorm, RC.unit_dot_self hne, ?_, ?_, ?_⟩
  · rw [smul_dotProduct, dotProduct_comm, dot_self_cross, smul_zero]
  · rw [smul_dotProduct, dotProduct_comm, dot_cross_self, smul_zero]
  · intro h; rw [hnorm]; exact max_eq_left h

/-! ## The whole fit and its scale covariance -/

namespace Coef

variable (a : Coef)

def kmin : ℝ := if |a.kp| < |a.km| then a.kp else a.km
def kmax : ℝ := if |a.kp| < |a.km| then a.km else a.kp
def wmin : Fin 2 → ℝ := if |a.kp| < |a.km| then a.evp else a.evm
/-- Minor radius `r = 1/|κmax|`. -/
def minorR : ℝ := 1 / |a.kmax|
/-- Major radius `R = 1/|κmin| − sgn(κ₊κ₋) r`. -/
def majorR : ℝ := 1 / |a.kmin| - Real.sign (a.kp * a.km) * a.minorR
/-- The sign of the torus (Eq. 22). -/
def tsign : ℝ := torusSign a.kmin a.kmax

/-- The coefficients for the same network output on the cloud scaled by `c`: `a00 ↦ c a00`,
first order unchanged, second order `↦ ·/c` (SPEC §4.2 with `σ ↦ cσ`). -/
def scale (c : ℝ) : Coef := ⟨c * a.a00, a.a10, a.a01, a.a11 / c, a.a20 / c, a.a02 / c⟩

end Coef

/-- The torus center `c = p + a00 n − (sign/|κmin|) n*`. -/
def fitCenter (Fr : Frame) (p : V3) (a : Coef) : V3 :=
  p + a.a00 • Fr.n - (a.tsign / |a.kmin|) • nstar Fr a

/-- The torus axis `u = (n* × v_min)/|n* × v_min|`. -/
def fitAxis (Fr : Frame) (a : Coef) : V3 :=
  (vnorm (nstar Fr a ⨯₃ vdir Fr a a.wmin))⁻¹ • (nstar Fr a ⨯₃ vdir Fr a a.wmin)

/-- The torus SDF of Eq. 23–24 with the signed minor radius `r`
(`g(x) = sgn(r)(|d| − |r|)`, `d = (|(x − c) × u| − R, (x − c)·u)`). -/
def torusSDF (c u : V3) (R r : ℝ) (x : V3) : ℝ :=
  Real.sign r * (√((vnorm ((x - c) ⨯₃ u) - R) ^ 2 + ((x - c) ⬝ᵥ u) ^ 2) - |r|)

theorem sign_mul_pos' {x y : ℝ} (hy : 0 < y) : Real.sign (x * y) = Real.sign x := by
  rcases lt_trichotomy x 0 with h | h | h
  · rw [Real.sign_of_neg h, Real.sign_of_neg (mul_neg_of_neg_of_pos h hy)]
  · rw [h, zero_mul]
  · rw [Real.sign_of_pos h, Real.sign_of_pos (mul_pos h hy)]

section scaling

variable (a : Coef) {c : ℝ} (hc : 0 < c)
include hc

theorem scale_S : (a.scale c).S = a.S := rfl
theorem scale_A : (a.scale c).A = a.A := rfl
theorem scale_v1 : (a.scale c).v1 = a.v1 / c := by
  simp only [Coef.v1, Coef.scale]; field_simp
theorem scale_v2 : (a.scale c).v2 = a.v2 / c := by
  simp only [Coef.v2, Coef.scale]; field_simp
theorem scale_term : (a.scale c).term = a.term / c := by
  simp only [Coef.term, Coef.scale]; field_simp
theorem scale_sqrtΔ : √(a.scale c).Δ = √a.Δ / c := by
  have : (a.scale c).Δ = a.Δ / c ^ 2 := by
    have ht := scale_term a hc
    simp only [Coef.Δ, ht]; simp only [Coef.scale]; field_simp
  rw [this, Real.sqrt_div' _ (by positivity), Real.sqrt_sq hc.le]
theorem scale_kp : (a.scale c).kp = a.kp / c := by
  simp only [Coef.kp, scale_term a hc, scale_sqrtΔ a hc, scale_A a hc]; field_simp
theorem scale_km : (a.scale c).km = a.km / c := by
  simp only [Coef.km, scale_term a hc, scale_sqrtΔ a hc, scale_A a hc]; field_simp

theorem scale_lt : (|(a.scale c).kp| < |(a.scale c).km|) ↔ (|a.kp| < |a.km|) := by
  rw [scale_kp a hc, scale_km a hc, abs_div, abs_div, abs_of_pos hc]
  exact div_lt_div_iff_of_pos_right hc

theorem scale_kmin : (a.scale c).kmin = a.kmin / c := by
  unfold Coef.kmin; rw [if_congr (scale_lt a hc) rfl rfl]; split_ifs
  · exact scale_kp a hc
  · exact scale_km a hc
theorem scale_kmax : (a.scale c).kmax = a.kmax / c := by
  unfold Coef.kmax; rw [if_congr (scale_lt a hc) rfl rfl]; split_ifs
  · exact scale_km a hc
  · exact scale_kp a hc
theorem scale_wmin : (a.scale c).wmin = c⁻¹ • a.wmin := by
  unfold Coef.wmin; rw [if_congr (scale_lt a hc) rfl rfl]
  split_ifs <;> ext i <;> fin_cases i <;>
    simp [Coef.evp, Coef.evm, scale_v1 a hc, scale_v2 a hc, scale_sqrtΔ a hc] <;> ring
theorem scale_tsign : (a.scale c).tsign = a.tsign := by
  have hc' : 0 < c⁻¹ := inv_pos.mpr hc
  unfold Coef.tsign torusSign
  rw [scale_kmin a hc, scale_kmax a hc]
  simp only [div_eq_mul_inv, mul_neg_iff, mul_pos_iff, hc', not_lt.mpr hc'.le, and_true,
    and_false, or_false, false_or]
theorem scale_minor : (a.scale c).minorR = c * a.minorR := by
  unfold Coef.minorR; rw [scale_kmax a hc, abs_div, abs_of_pos hc]; field_simp
theorem scale_major : (a.scale c).majorR = c * a.majorR := by
  unfold Coef.majorR
  have hs : Real.sign ((a.scale c).kp * (a.scale c).km) = Real.sign (a.kp * a.km) := by
    rw [scale_kp a hc, scale_km a hc,
      show a.kp / c * (a.km / c) = a.kp * a.km * (c⁻¹ ^ 2) by field_simp]
    exact sign_mul_pos' (by positivity : (0:ℝ) < c⁻¹ ^ 2)
  rw [hs, scale_minor a hc, scale_kmin a hc, abs_div, abs_of_pos hc]
  rw [one_div_div, mul_sub, one_div]
  ring

theorem scale_center (Fr : Frame) (p : V3) :
    fitCenter Fr (c • p) (a.scale c) = c • fitCenter Fr p a := by
  unfold fitCenter
  rw [scale_tsign a hc, scale_kmin a hc, abs_div, abs_of_pos hc, div_div_eq_mul_div]
  have hn : nstar Fr (a.scale c) = nstar Fr a := rfl
  rw [hn]
  simp only [Coef.scale, smul_sub, smul_add, smul_smul]
  congr 2
  ring

theorem scale_axis (Fr : Frame) : fitAxis Fr (a.scale c) = fitAxis Fr a := by
  unfold fitAxis
  have hn : nstar Fr (a.scale c) = nstar Fr a := rfl
  have hv : vdir Fr (a.scale c) (a.scale c).wmin = c⁻¹ • vdir Fr a a.wmin := by
    rw [scale_wmin a hc]
    ext i; simp only [vdir, Coef.scale, Pi.smul_apply, Pi.add_apply, smul_eq_mul]; ring
  rw [hn, hv, map_smul, vnorm_smul, abs_of_pos (inv_pos.mpr hc), smul_smul, mul_inv, inv_inv]
  congr 1
  field_simp

theorem torusSDF_scale (ctr u x : V3) (R r : ℝ) :
    torusSDF (c • ctr) u (c * R) (c * r) (c • x) = c * torusSDF ctr u R r x := by
  unfold torusSDF
  have hx : (c • (x - ctr)) ⨯₃ u = c • ((x - ctr) ⨯₃ u) := by simp [map_smul]
  rw [← smul_sub, hx, vnorm_smul, abs_of_pos hc, smul_dotProduct, smul_eq_mul]
  rw [show (c * vnorm ((x - ctr) ⨯₃ u) - c * R) ^ 2 + (c * ((x - ctr) ⬝ᵥ u)) ^ 2 =
      c ^ 2 * ((vnorm ((x - ctr) ⨯₃ u) - R) ^ 2 + ((x - ctr) ⬝ᵥ u) ^ 2) by ring]
  rw [Real.sqrt_mul (by positivity), Real.sqrt_sq hc.le, abs_mul, abs_of_pos hc,
    mul_comm c r, sign_mul_pos' hc]
  ring

end scaling

/-- **Output scaling** (SPEC §4.2): with the network outputs fixed, scaling the cloud by
`c > 0` (about the origin; `σ ↦ cσ`, so `a ↦ a.scale c`) scales the torus center and both radii
by `c`, keeps the sign and the axis, and scales the SDF: `g'(c x) = c g(x)`. This holds for the
exact formulas; the absolute 1e-10 clamps of D9 break it only when they are active.

[formula: pat-output-scaling] -/
theorem fit_scale (Fr : Frame) (p : V3) (a : Coef) {c : ℝ} (hc : 0 < c) (x : V3) :
    fitCenter Fr (c • p) (a.scale c) = c • fitCenter Fr p a ∧
    fitAxis Fr (a.scale c) = fitAxis Fr a ∧
    (a.scale c).majorR = c * a.majorR ∧ (a.scale c).minorR = c * a.minorR ∧
    (a.scale c).tsign = a.tsign ∧
    torusSDF (fitCenter Fr (c • p) (a.scale c)) (fitAxis Fr (a.scale c)) (a.scale c).majorR
        ((a.scale c).tsign * (a.scale c).minorR) (c • x) =
      c * torusSDF (fitCenter Fr p a) (fitAxis Fr a) a.majorR (a.tsign * a.minorR) x := by
  refine ⟨scale_center a hc Fr p, scale_axis a hc Fr, scale_major a hc, scale_minor a hc,
    scale_tsign a hc, ?_⟩
  rw [scale_center a hc Fr p, scale_axis a hc Fr, scale_major a hc, scale_minor a hc,
    scale_tsign a hc, show a.tsign * (c * a.minorR) = c * (a.tsign * a.minorR) by ring,
    torusSDF_scale hc]

end PAT

end DDG
