import DdgProofs.Developability.Vec

/-!
# Developability: invariances of the energy and the identities of its gradient (SPEC §3)

`crates/ddg-developability/SPEC.md` §3 (first paragraph) and §10.4. All six energies depend on the
positions only through unit normals and angles, so they are invariant under translations, rotations
and uniform scalings of the mesh. For any differentiable function `E` of the positions
`f : Fin n → ℝ³` with gradient `g` (`dE(f)[δ] = Σ_v g_v · δ_v`) this gives

```
Σ_v g_v = 0,      Σ_v f_v × g_v = 0,      Σ_v f_v · g_v = 0.
```

(translation, rotation, scaling). Translations and scalings need only the derivative along a line;
rotations use the chain rule along the rotation `s ↦ R(s) f` (Fréchet derivative, `HasFDerivAt`).
The invariance is assumed only along the orbits through `f`, about the three coordinate axes.

Positions here are `Fin n → Fin 3 → ℝ` (a normed space, so that `HasFDerivAt` applies); the cross
and dot products are Mathlib's `crossProduct` and `dotProduct`.
-/

noncomputable section

namespace DdgProofs

namespace Dev

open Real Filter Topology

/-- Positions of `n` vertices. -/
abbrev Positions (n : ℕ) := Fin n → Fin 3 → ℝ

/-- Rotation by the angle `s` about the `x` axis. -/
def rotX (s : ℝ) (p : Fin 3 → ℝ) : Fin 3 → ℝ :=
  ![p 0, cos s * p 1 - sin s * p 2, sin s * p 1 + cos s * p 2]

/-- Rotation by the angle `s` about the `y` axis. -/
def rotY (s : ℝ) (p : Fin 3 → ℝ) : Fin 3 → ℝ :=
  ![cos s * p 0 + sin s * p 2, p 1, -(sin s * p 0) + cos s * p 2]

/-- Rotation by the angle `s` about the `z` axis. -/
def rotZ (s : ℝ) (p : Fin 3 → ℝ) : Fin 3 → ℝ :=
  ![cos s * p 0 - sin s * p 1, sin s * p 0 + cos s * p 1, p 2]

section

variable {n : ℕ} {E : Positions n → ℝ} {E' : Positions n →L[ℝ] ℝ} {f : Positions n}

/-- If `E` is constant near `0` along a curve `γ` through `f` (`γ 0 = f`) with velocity `γ'`,
then `dE(f)[γ'] = 0`. -/
theorem fderiv_apply_eq_zero_of_eventually_const (hE : HasFDerivAt E E' f) {γ : ℝ → Positions n}
    {γ' : Positions n} (hγ : HasDerivAt γ γ' 0) (h0 : γ 0 = f)
    (hc : ∀ᶠ s in 𝓝 (0 : ℝ), E (γ s) = E f) : E' γ' = 0 := by
  have h1 : HasDerivAt (E ∘ γ) (E' γ') 0 := hE.comp_hasDerivAt_of_eq 0 hγ h0.symm
  have h2 : HasDerivAt (E ∘ γ) 0 0 :=
    (hasDerivAt_const (0 : ℝ) (E f)).congr_of_eventuallyEq (by
      filter_upwards [hc] with s hs using hs)
  exact h1.unique h2

/-- Derivative of a curve of positions, vertex by vertex and coordinate by coordinate. -/
theorem hasDerivAt_positions {γ : ℝ → Positions n} {γ' : Positions n} {t : ℝ}
    (h : ∀ v c, HasDerivAt (fun s => γ s v c) (γ' v c) t) : HasDerivAt γ γ' t :=
  hasDerivAt_pi.mpr fun v => hasDerivAt_pi.mpr fun c => h v c

theorem rot_hasDerivAt (a b : ℝ) :
    HasDerivAt (fun s => cos s * a - sin s * b) (-b) 0 ∧
      HasDerivAt (fun s => sin s * a + cos s * b) a 0 := by
  constructor
  · have := ((hasDerivAt_cos 0).mul_const a).sub ((hasDerivAt_sin 0).mul_const b)
    exact this.congr_deriv (by simp)
  · have := ((hasDerivAt_sin 0).mul_const a).add ((hasDerivAt_cos 0).mul_const b)
    exact this.congr_deriv (by simp)

theorem rotY_hasDerivAt (a b : ℝ) :
    HasDerivAt (fun s => cos s * a + sin s * b) b 0 ∧
      HasDerivAt (fun s => -(sin s * a) + cos s * b) (-a) 0 := by
  constructor
  · have := ((hasDerivAt_cos 0).mul_const a).add ((hasDerivAt_sin 0).mul_const b)
    exact this.congr_deriv (by simp)
  · have := (((hasDerivAt_sin 0).mul_const a).neg).add ((hasDerivAt_cos 0).mul_const b)
    exact this.congr_deriv (by simp)

/-- **Translation invariance ⇒ the gradient sums to zero.** -/
theorem sum_grad_eq_zero (hE : HasFDerivAt E E' f) (g : Positions n)
    (hg : ∀ δ, E' δ = ∑ v, ∑ c, g v c * δ v c)
    (htr : ∀ (a : Fin 3 → ℝ) (s : ℝ), E (f + s • fun _ => a) = E f) :
    ∑ v, g v = 0 := by
  funext c
  have hγ : HasDerivAt (fun s : ℝ => f + s • fun (_ : Fin n) => Pi.single c (1 : ℝ))
      (fun _ => Pi.single c 1) 0 := by
    simpa using ((hasDerivAt_id (0 : ℝ)).smul_const
      (fun (_ : Fin n) => (Pi.single c (1 : ℝ) : Fin 3 → ℝ))).const_add f
  have := fderiv_apply_eq_zero_of_eventually_const hE hγ (by simp)
    (Eventually.of_forall fun s => htr _ s)
  rw [hg] at this
  simpa [Finset.sum_apply, Pi.single_apply] using this

/-- **Scale invariance ⇒ `Σ_v f_v · g_v = 0`.** -/
theorem sum_dot_grad_eq_zero (hE : HasFDerivAt E E' f) (g : Positions n)
    (hg : ∀ δ, E' δ = ∑ v, ∑ c, g v c * δ v c)
    (hsc : ∀ s : ℝ, 0 < s → E (s • f) = E f) :
    ∑ v, dotProduct (f v) (g v) = 0 := by
  have hγ : HasDerivAt (fun s : ℝ => (1 + s) • f) f 0 := by
    simpa using ((hasDerivAt_id (0 : ℝ)).const_add 1).smul_const f
  have hev : ∀ᶠ s in 𝓝 (0 : ℝ), E ((1 + s) • f) = E f := by
    have : ∀ᶠ s in 𝓝 (0 : ℝ), -1 < s := lt_mem_nhds (by norm_num)
    filter_upwards [this] with s hs using hsc _ (by linarith)
  have := fderiv_apply_eq_zero_of_eventually_const hE hγ (by simp) hev
  rw [hg] at this
  simpa [dotProduct, mul_comm] using this

/-- **Rotation invariance ⇒ `Σ_v f_v × g_v = 0`** (rotations about the three coordinate axes). -/
theorem sum_cross_grad_eq_zero (hE : HasFDerivAt E E' f) (g : Positions n)
    (hg : ∀ δ, E' δ = ∑ v, ∑ c, g v c * δ v c)
    (hx : ∀ s, E (fun v => rotX s (f v)) = E f) (hy : ∀ s, E (fun v => rotY s (f v)) = E f)
    (hz : ∀ s, E (fun v => rotZ s (f v)) = E f) :
    ∑ v, crossProduct (f v) (g v) = 0 := by
  have ev : ∀ p : Fin 3 → ℝ, ![p 0, p 1, p 2] = p := fun p => by
    funext c; fin_cases c <;> rfl
  -- about `x`
  have dX : E' (fun v => ![0, -(f v 2), f v 1]) = 0 := by
    refine fderiv_apply_eq_zero_of_eventually_const hE (γ := fun s v => rotX s (f v))
      (hasDerivAt_positions fun v c => ?_) ?_ (Eventually.of_forall hx)
    · fin_cases c
      · simpa [rotX] using hasDerivAt_const (0 : ℝ) (f v 0)
      · simpa [rotX] using (rot_hasDerivAt (f v 1) (f v 2)).1
      · simpa [rotX] using (rot_hasDerivAt (f v 1) (f v 2)).2
    · funext v; simp [rotX, ev]
  have dY : E' (fun v => ![f v 2, 0, -(f v 0)]) = 0 := by
    refine fderiv_apply_eq_zero_of_eventually_const hE (γ := fun s v => rotY s (f v))
      (hasDerivAt_positions fun v c => ?_) ?_ (Eventually.of_forall hy)
    · fin_cases c
      · simpa [rotY] using (rotY_hasDerivAt (f v 0) (f v 2)).1
      · simpa [rotY] using hasDerivAt_const (0 : ℝ) (f v 1)
      · simpa [rotY] using (rotY_hasDerivAt (f v 0) (f v 2)).2
    · funext v; simp [rotY, ev]
  have dZ : E' (fun v => ![-(f v 1), f v 0, 0]) = 0 := by
    refine fderiv_apply_eq_zero_of_eventually_const hE (γ := fun s v => rotZ s (f v))
      (hasDerivAt_positions fun v c => ?_) ?_ (Eventually.of_forall hz)
    · fin_cases c
      · simpa [rotZ] using (rot_hasDerivAt (f v 0) (f v 1)).1
      · simpa [rotZ] using (rot_hasDerivAt (f v 0) (f v 1)).2
      · simpa [rotZ] using hasDerivAt_const (0 : ℝ) (f v 2)
    · funext v; simp [rotZ, ev]
  rw [hg] at dX dY dZ
  simp only [Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] at dX dY dZ
  funext c
  rw [Finset.sum_apply]
  fin_cases c <;> simp only [cross_apply, Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons, Pi.zero_apply]
  · rw [← dX]
    refine Finset.sum_congr rfl fun v _ => by ring
  · rw [← dY]
    refine Finset.sum_congr rfl fun v _ => by ring
  · rw [← dZ]
    refine Finset.sum_congr rfl fun v _ => by ring

/-- **Invariances of the energy and its gradient** (SPEC §3, §10.4): let `E` be differentiable at
the positions `f` with gradient `g` (`dE(f)[δ] = Σ_v g_v · δ_v`), and invariant along the
translations, the rotations about the three coordinate axes and the uniform scalings of `f`.
Then `Σ_v g_v = 0`, `Σ_v f_v × g_v = 0` and `Σ_v f_v · g_v = 0`: the gradient flow does not
translate, rotate or scale the mesh to first order.

[formula: dev-gradient-invariance] -/
theorem gradient_invariance (hE : HasFDerivAt E E' f) (g : Positions n)
    (hg : ∀ δ, E' δ = ∑ v, ∑ c, g v c * δ v c)
    (htr : ∀ (a : Fin 3 → ℝ) (s : ℝ), E (f + s • fun _ => a) = E f)
    (hx : ∀ s, E (fun v => rotX s (f v)) = E f) (hy : ∀ s, E (fun v => rotY s (f v)) = E f)
    (hz : ∀ s, E (fun v => rotZ s (f v)) = E f) (hsc : ∀ s : ℝ, 0 < s → E (s • f) = E f) :
    ∑ v, g v = 0 ∧ ∑ v, crossProduct (f v) (g v) = 0 ∧ ∑ v, dotProduct (f v) (g v) = 0 :=
  ⟨sum_grad_eq_zero hE g hg htr, sum_cross_grad_eq_zero hE g hg hx hy hz,
    sum_dot_grad_eq_zero hE g hg hsc⟩

end

end Dev

end DdgProofs
