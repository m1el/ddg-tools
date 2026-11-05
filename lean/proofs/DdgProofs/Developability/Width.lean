import DdgProofs.Developability.Vec

/-!
# Developability: the branching-free intrinsic width (SPEC §3.7, paper App. B.5.1)

`crates/ddg-developability/SPEC.md` §3.7, decision D7. The points are the log-mapped normals
`Ñ_σ ∈ T_{n0}S²` of a star (`Ddg.Dev.logMap`); here they are arbitrary points `P 0 … P n`.

* `pairWidthSq_eq_sq_widthAlong`: `ψ(u) = max_{σ,τ} ⟨P_σ − P_τ, u⟩² = (max_σ ⟨P_σ,u⟩ − min_σ ⟨P_σ,u⟩)²`.
* `hasDerivAt3_altitudeDir`: the unit altitude `u = n0 × ŵ`, `w = x_b − x_a`, has the derivative
  `du = −(u·dw) ŵ/|w|` for an in-plane motion (`w ⊥ n0`, `dw ⊥ n0`, `n0` fixed and unit). The
  paper prints `∇_{x_a} u* = −(N_i × w) wᵀ/|w|³`, the negated transpose of this Jacobian (D7).
* `gramWidth_eq_sq` and `hasDerivAt_gramWidth`: the energy `Δ² = ((x_s − x_t)·u)²` of the active
  pair is, for points in the plane `⊥ n0`, the Gram expression
  `W = |d|² − (d·w)²/|w|²` (`d = x_s − x_t`), which depends only on inner products of the points;
  along **arbitrary** curves of the four points its derivative is `Σ z_q · x_q'` with the
  sensitivities of `Ddg.Dev.minWidthSensitivities`. This is why the out-of-plane parts of `dÑ` (the
  motion of `n0` itself) do not matter (SPEC §3.7).

That the minimum of `ψ` over unit `u ⊥ n0` is attained at a hull-edge normal (rotating calipers)
is cited (`dev-min-width-hull-edge`, CITED.md).
-/

noncomputable section

namespace DDG

namespace Dev

open scoped DDG.Oracle
open Real Ddg.Dev

/-! ## Finite maxima and minima -/

section maxmin

variable (f : ℕ → ℝ)

theorem width_maxTo1_zero : maxTo1 0 f = f 0 := rfl
theorem width_maxTo1_succ (n : ℕ) : maxTo1 (n + 1) f = max (maxTo1 n f) (f (n + 1)) := rfl
theorem width_minTo1_zero : minTo1 0 f = f 0 := rfl
theorem width_minTo1_succ (n : ℕ) : minTo1 (n + 1) f = min (minTo1 n f) (f (n + 1)) := rfl

theorem width_le_maxTo1 : ∀ (n i : ℕ), i ≤ n → f i ≤ maxTo1 n f
  | 0, i, hi => by rw [Nat.le_zero.mp hi, width_maxTo1_zero]
  | n + 1, i, hi => by
    rw [width_maxTo1_succ]
    rcases Nat.lt_or_eq_of_le hi with h | rfl
    · exact le_trans (width_le_maxTo1 n i (Nat.lt_succ_iff.mp h)) (le_max_left _ _)
    · exact le_max_right _ _

theorem width_maxTo1_le {b : ℝ} : ∀ (n : ℕ), (∀ i ≤ n, f i ≤ b) → maxTo1 n f ≤ b
  | 0, h => by rw [width_maxTo1_zero]; exact h 0 le_rfl
  | n + 1, h => by
    rw [width_maxTo1_succ]
    exact max_le (width_maxTo1_le n fun i hi => h i (Nat.le_succ_of_le hi)) (h _ le_rfl)

theorem width_exists_maxTo1 : ∀ (n : ℕ), ∃ i ≤ n, maxTo1 n f = f i
  | 0 => ⟨0, le_rfl, rfl⟩
  | n + 1 => by
    obtain ⟨i, hi, e⟩ := width_exists_maxTo1 n
    rw [width_maxTo1_succ]
    rcases le_total (maxTo1 n f) (f (n + 1)) with h | h
    · exact ⟨n + 1, le_rfl, max_eq_right h⟩
    · exact ⟨i, Nat.le_succ_of_le hi, by rw [max_eq_left h, e]⟩

theorem width_minTo1_le : ∀ (n i : ℕ), i ≤ n → minTo1 n f ≤ f i
  | 0, i, hi => by rw [Nat.le_zero.mp hi, width_minTo1_zero]
  | n + 1, i, hi => by
    rw [width_minTo1_succ]
    rcases Nat.lt_or_eq_of_le hi with h | rfl
    · exact le_trans (min_le_left _ _) (width_minTo1_le n i (Nat.lt_succ_iff.mp h))
    · exact min_le_right _ _

theorem width_exists_minTo1 : ∀ (n : ℕ), ∃ i ≤ n, minTo1 n f = f i
  | 0 => ⟨0, le_rfl, rfl⟩
  | n + 1 => by
    obtain ⟨i, hi, e⟩ := width_exists_minTo1 n
    rw [width_minTo1_succ]
    rcases le_total (minTo1 n f) (f (n + 1)) with h | h
    · exact ⟨i, Nat.le_succ_of_le hi, by rw [min_eq_left h, e]⟩
    · exact ⟨n + 1, le_rfl, min_eq_right h⟩

end maxmin

/-! ## The width along a direction -/

/-- **Width along `u` of the pairs = squared width** (SPEC §3.7): over the `n + 1` points
`P 0 … P n`, `max_{σ,τ} ⟨P_σ − P_τ, u⟩² = (max_σ ⟨P_σ, u⟩ − min_σ ⟨P_σ, u⟩)²`. So `ψ(u)` is the
squared width of the point set in the direction `u`, and minimizing it over unit `u` is the
minimum width of the convex hull.

[formula: dev-intrinsic-width] -/
theorem pairWidthSq_eq_sq_widthAlong (n : ℕ) (P : ℕ → V3) (u : V3) :
    pairWidthSq n P u = widthAlong n P u ^ 2 := by
  set p : ℕ → ℝ := fun m => dot (P m) u with hp
  have hd : ∀ i j, dot (sub (P i) (P j)) u = p i - p j := by
    intro i j; simp only [hp, dot_def, sub_x, sub_y, sub_z]; ring
  have hw : widthAlong n P u = maxTo1 n p - minTo1 n p := rfl
  have hps : pairWidthSq n P u =
      maxTo1 n (fun i => maxTo1 n (fun j => (p i - p j) * (p i - p j))) := by
    unfold pairWidthSq; simp only [hd]
  rw [hps, hw]
  apply le_antisymm
  · apply width_maxTo1_le
    intro i hi
    apply width_maxTo1_le
    intro j hj
    have a1 := width_le_maxTo1 p n i hi
    have a2 := width_minTo1_le p n i hi
    have a3 := width_le_maxTo1 p n j hj
    have a4 := width_minTo1_le p n j hj
    rw [← sq]
    exact sq_le_sq' (by linarith) (by linarith)
  · obtain ⟨i, hi, ei⟩ := width_exists_maxTo1 p n
    obtain ⟨j, hj, ej⟩ := width_exists_minTo1 p n
    rw [ei, ej, sq]
    exact le_trans (width_le_maxTo1 (fun j => (p i - p j) * (p i - p j)) n j hj)
      (width_le_maxTo1 (fun i => maxTo1 n (fun j => (p i - p j) * (p i - p j))) n i hi)

/-! ## An orthonormal frame `v, n0 × v, n0` -/

/-- Completeness of the orthonormal frame `(v, m × v, m)` (`m`, `v` unit and orthogonal):
`X = (v·X) v + ((m×v)·X)(m×v) + (m·X) m`. -/
theorem frame_expand (m v X : V3) (hm : norm2 m = 1) (hv : norm2 v = 1) (hmv : dot m v = 0) :
    X = add (add (smul (dot v X) v) (smul (dot (cross m v) X) (cross m v))) (smul (dot m X) m) := by
  rw [norm2_def] at hm hv
  rw [dot_def] at hmv
  ext <;> simp only [add_x, add_y, add_z, smul_x, smul_y, smul_z, cross_x, cross_y, cross_z,
      dot_def]
  · linear_combination (-X.x + (m.x * X.x + m.y * X.y + m.z * X.z) * m.x) * hv +
      (-(v.x * v.x + v.y * v.y + v.z * v.z) * X.x + (v.x * X.x + v.y * X.y + v.z * X.z) * v.x) * hm +
      ((m.x * v.x + m.y * v.y + m.z * v.z) * X.x - (m.x * X.x + m.y * X.y + m.z * X.z) * v.x -
        (v.x * X.x + v.y * X.y + v.z * X.z) * m.x) * hmv
  · linear_combination (-X.y + (m.x * X.x + m.y * X.y + m.z * X.z) * m.y) * hv +
      (-(v.x * v.x + v.y * v.y + v.z * v.z) * X.y + (v.x * X.x + v.y * X.y + v.z * X.z) * v.y) * hm +
      ((m.x * v.x + m.y * v.y + m.z * v.z) * X.y - (m.x * X.x + m.y * X.y + m.z * X.z) * v.y -
        (v.x * X.x + v.y * X.y + v.z * X.z) * m.y) * hmv
  · linear_combination (-X.z + (m.x * X.x + m.y * X.y + m.z * X.z) * m.z) * hv +
      (-(v.x * v.x + v.y * v.y + v.z * v.z) * X.z + (v.x * X.x + v.y * X.y + v.z * X.z) * v.z) * hm +
      ((m.x * v.x + m.y * v.y + m.z * v.z) * X.z - (m.x * X.x + m.y * X.y + m.z * X.z) * v.z -
        (v.x * X.x + v.y * X.y + v.z * X.z) * m.z) * hmv

/-- In the plane `⊥ m`: `X − (v·X) v = ((m×v)·X)(m×v)`. -/
theorem sub_proj_eq (m v X : V3) (hm : norm2 m = 1) (hv : norm2 v = 1) (hmv : dot m v = 0)
    (hmX : dot m X = 0) :
    sub X (smul (dot v X) v) = smul (dot (cross m v) X) (cross m v) := by
  have h := frame_expand m v X hm hv hmv
  rw [hmX] at h
  ext
  · have := congrArg Ddg.Point.x h; simp at this ⊢; linarith
  · have := congrArg Ddg.Point.y h; simp at this ⊢; linarith
  · have := congrArg Ddg.Point.z h; simp at this ⊢; linarith

/-- `m × (m × v) = −v` for unit `m ⊥ v`. -/
theorem cross_cross_self (m v : V3) (hm : norm2 m = 1) (hmv : dot m v = 0) :
    cross m (cross m v) = neg v := by
  rw [cross_cross, hmv]
  have : dot m m = 1 := hm
  rw [this]
  ext <;> simp

/-! ## Derivative of the altitude direction -/

/-- **Derivative of the unit altitude** (SPEC §3.7, corrected B.5.1 formula, D7): with `n0` fixed
and unit, `w(t) = x_b(t) − x_a(t)`, `w ≠ 0`, `w ⊥ n0` and `w' ⊥ n0` at `t`, the direction
`u = n0 × w/|w|` has the derivative `−(u·w') ŵ/|w|` (`ŵ = w/|w|`), i.e. the Jacobians
`∂u/∂x_a = ŵ uᵀ/|w|`, `∂u/∂x_b = −ŵ uᵀ/|w|`.

[formula: dev-width-altitude-gradient] -/
theorem hasDerivAt3_altitudeDir {xa xb : ℝ → V3} {xa' xb' : V3} {t : ℝ} (n0 : V3)
    (ha : HasDerivAt3 xa xa' t) (hb : HasDerivAt3 xb xb' t) (hn0 : norm2 n0 = 1)
    (hw0 : sub (xb t) (xa t) ≠ zero) (hw : dot n0 (sub (xb t) (xa t)) = 0)
    (hw' : dot n0 (sub xb' xa') = 0) :
    HasDerivAt3 (fun s => altitudeDir n0 (xa s) (xb s))
      (altitudeDerivative n0 (xa t) (xb t) (sub xb' xa')) t := by
  set w := sub (xb t) (xa t) with hwdef
  set w' := sub xb' xa'
  have hW : HasDerivAt3 (fun s => sub (xb s) (xa s)) w' t := hb.sub ha
  have hN := (HasDerivAt3.const n0 t).cross (hW.normalize hw0)
  unfold altitudeDir
  refine hN.of_eq ?_
  -- algebra
  set v := normalize w with hv
  have hvn : norm2 v = 1 := norm2_normalize w hw0
  have hmv : dot n0 v = 0 := by
    rw [hv, normalize_def]; simp only [dot_def, smul_x, smul_y, smul_z] at hw ⊢
    linear_combination (1 / norm w) * hw
  have hsp := sub_proj_eq n0 v w' hn0 hvn hmv hw'
  have hcc := cross_cross_self n0 v hn0 hmv
  unfold altitudeDerivative altitudeDir
  rw [← hv]
  -- the derivative is n0 × ((w' − (v·w') v)/|w|)
  have e1 : cross n0 (smul (1 / norm w) (sub w' (smul (dot v w') v))) =
      smul (1 / norm w) (cross n0 (sub w' (smul (dot v w') v))) := by
    ext <;> simp <;> ring
  have e0 : cross (zero : V3) (normalize w) = zero := by ext <;> simp
  rw [e0, e1, hsp]
  have e2 : cross n0 (smul (dot (cross n0 v) w') (cross n0 v)) =
      smul (dot (cross n0 v) w') (cross n0 (cross n0 v)) := by
    ext <;> simp <;> ring
  rw [e2, hcc]
  ext <;> simp [dot_def] <;> ring

/-! ## The energy of the active pair as a Gram expression, and its derivative -/

/-- `W = |d|² − (d·w)²/|w|²` with `d = x_s − x_t`, `w = x_b − x_a`. -/
def gramWidth (xs xt xa xb : V3) : ℝ :=
  norm2 (sub xs xt) - dot (sub xs xt) (sub xb xa) ^ 2 / norm2 (sub xb xa)

/-- **The energy of the active pair is the Gram expression** (SPEC §3.7): for `d = x_s − x_t`
and `w = x_b − x_a ≠ 0` in the plane `⊥ n0` (`n0` unit), `((x_s − x_t)·u)² = W` with
`u = altitudeDir n0 x_a x_b`.

[formula: dev-min-width-gradient] -/
theorem gramWidth_eq_sq (n0 xs xt xa xb : V3) (hn0 : norm2 n0 = 1) (hw0 : sub xb xa ≠ zero)
    (hd : dot n0 (sub xs xt) = 0) (hw : dot n0 (sub xb xa) = 0) :
    gramWidth xs xt xa xb = dot (sub xs xt) (altitudeDir n0 xa xb) ^ 2 := by
  set w := sub xb xa
  set d := sub xs xt
  set v := normalize w with hv
  have hvn : norm2 v = 1 := norm2_normalize w hw0
  have hmv : dot n0 v = 0 := by
    rw [hv, normalize_def]; simp only [dot_def, smul_x, smul_y, smul_z] at hw ⊢
    linear_combination (1 / norm w) * hw
  have hexp := frame_expand n0 v d hn0 hvn hmv
  rw [hd] at hexp
  -- |d|² = (v·d)² + (u·d)²
  have hd2 : norm2 d = dot v d ^ 2 + dot (cross n0 v) d ^ 2 := by
    conv_lhs => rw [norm2_def, hexp]
    simp only [add_x, add_y, add_z, smul_x, smul_y, smul_z]
    have hu : norm2 (cross n0 v) = 1 := by
      have := norm2_cross_add_sq_dot n0 v
      rw [hn0, hvn, hmv] at this; linarith
    have hvu : dot v (cross n0 v) = 0 := by rw [dot_comm]; exact dot_cross_self_right n0 v
    rw [norm2_def] at hvn hu
    rw [dot_def] at hvu
    linear_combination (dot v d) ^ 2 * hvn + (dot (cross n0 v) d) ^ 2 * hu +
      2 * dot v d * dot (cross n0 v) d * hvu
  -- (d·w)²/|w|² = (v·d)²
  have hw2 : 0 < norm2 w := (norm2_pos_iff w).mpr hw0
  have hvd : dot v d = dot d w / norm w := by
    rw [hv, normalize_def]; simp only [dot_def, smul_x, smul_y, smul_z]; ring
  unfold gramWidth altitudeDir
  rw [hd2, hvd, dot_comm d (cross n0 v), ← norm_sq w]
  have : 0 < norm w := (norm_pos_iff' w).mpr hw0
  field_simp
  ring

/-- The algebra behind `hasDerivAt_gramWidth`: the sensitivities applied to a motion equal the
derivative of the Gram expression, for `d, w ⊥ n0`. -/
theorem minWidthSensitivities_apply (n0 xs xt xa xb xs' xt' xa' xb' : V3) (hn0 : norm2 n0 = 1)
    (hw0 : sub xb xa ≠ zero) (hd : dot n0 (sub xs xt) = 0) (hw : dot n0 (sub xb xa) = 0) :
    let z := minWidthSensitivities n0 xs xt xa xb
    dot z.1 xs' + dot z.2.1 xt' + dot z.2.2.1 xa' + dot z.2.2.2 xb' =
      2 * dot (sub xs xt) (sub xs' xt') -
        (2 * dot (sub xs xt) (sub xb xa) * (dot (sub xs' xt') (sub xb xa) +
            dot (sub xs xt) (sub xb' xa')) * norm2 (sub xb xa) -
          dot (sub xs xt) (sub xb xa) ^ 2 * (2 * dot (sub xb xa) (sub xb' xa'))) /
          norm2 (sub xb xa) ^ 2 := by
  intro z
  have hvn : norm2 (normalize (sub xb xa)) = 1 := norm2_normalize _ hw0
  have hmv : dot n0 (normalize (sub xb xa)) = 0 := by
    rw [normalize_def]; simp only [dot_def, smul_x, smul_y, smul_z] at hw ⊢
    linear_combination (1 / norm (sub xb xa)) * hw
  have hsp := sub_proj_eq n0 (normalize (sub xb xa)) (sub xs xt) hn0 hvn hmv hd
  have hz : dot z.1 xs' + dot z.2.1 xt' + dot z.2.2.1 xa' + dot z.2.2.2 xb' =
      2 * dot (smul (dot (cross n0 (normalize (sub xb xa))) (sub xs xt))
          (cross n0 (normalize (sub xb xa)))) (sub xs' xt') -
        2 * (dot (sub xs xt) (normalize (sub xb xa)) / norm (sub xb xa)) *
          dot (smul (dot (cross n0 (normalize (sub xb xa))) (sub xs xt))
            (cross n0 (normalize (sub xb xa)))) (sub xb' xa') := by
    simp only [z, minWidthSensitivities, altitudeDir]
    simp only [dot_def, smul_x, smul_y, smul_z, sub_x, sub_y, sub_z]
    ring
  rw [hz, ← hsp, normalize_def]
  have hpos : 0 < norm (sub xb xa) := (norm_pos_iff' _).mpr hw0
  rw [← norm_sq (sub xb xa)]
  simp only [dot_def, smul_x, smul_y, smul_z, sub_x, sub_y, sub_z]
  field_simp
  ring

/-- **Gradient of the branching-free intrinsic width** (SPEC §3.7): along arbitrary
differentiable curves `x_s, x_t, x_a, x_b` (no constraint on the motion), if at `t` the vectors
`d = x_s − x_t` and `w = x_b − x_a ≠ 0` lie in the plane `⊥ n0` (`n0` unit), then the Gram
expression `W` has the derivative `z_s·x_s' + z_t·x_t' + z_a·x_a' + z_b·x_b'` with
`(z_s, z_t, z_a, z_b) = minWidthSensitivities n0 x_s x_t x_a x_b`
(`z_s = 2Δu`, `z_t = −2Δu`, `z_a = 2Δ(d·ŵ)/|w| u`, `z_b = −z_a`, `Δ = d·u`).

[formula: dev-min-width-gradient] -/
theorem hasDerivAt_gramWidth {xs xt xa xb : ℝ → V3} {xs' xt' xa' xb' : V3} {t : ℝ} (n0 : V3)
    (hs : HasDerivAt3 xs xs' t) (ht : HasDerivAt3 xt xt' t) (ha : HasDerivAt3 xa xa' t)
    (hb : HasDerivAt3 xb xb' t) (hn0 : norm2 n0 = 1) (hw0 : sub (xb t) (xa t) ≠ zero)
    (hd : dot n0 (sub (xs t) (xt t)) = 0) (hw : dot n0 (sub (xb t) (xa t)) = 0) :
    let z := minWidthSensitivities n0 (xs t) (xt t) (xa t) (xb t)
    HasDerivAt (fun s => gramWidth (xs s) (xt s) (xa s) (xb s))
      (dot z.1 xs' + dot z.2.1 xt' + dot z.2.2.1 xa' + dot z.2.2.2 xb') t := by
  intro z
  have hD := hs.sub ht
  have hW := hb.sub ha
  have hw2 : 0 < norm2 (sub (xb t) (xa t)) := (norm2_pos_iff _).mpr hw0
  have hder := hD.norm2.sub (((hD.dot hW).pow 2).div hW.norm2 hw2.ne')
  unfold gramWidth
  refine hder.congr_deriv ?_
  rw [minWidthSensitivities_apply n0 _ _ _ _ xs' xt' xa' xb' hn0 hw0 hd hw]
  simp only [Pi.pow_apply]
  norm_num

end Dev

end DDG
