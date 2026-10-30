import DdgProofs.RepulsiveCurves.NearPair

/-!
# Repulsive Curves: second-neighbour rules of the near-pair quadrature (WP16b; SPEC §17p)

Options of `NearPairOptions` that are off by default; not in the paper. For a pair `I, K, J`
joined by one edge `K`, the one-sided rule replaces the blend of `rc-near-blend` by
`e = E^p + w E^p ψ(E^n/E^p − 1)` and the gated rule multiplies the correction by
`g = W(|s* − s_K|/σ) W(|t* − t_K|/σ)`. Machine-checked here:

* `rc-near-onesided`: the polynomial piece `ψ(z) = 6z³ − 8z⁴ + 3z⁵` of the smooth positive part
  has `ψ(0) = ψ'(0) = ψ''(0) = 0`, `ψ(1) = ψ'(1) = 1`, `ψ''(1) = 0` (so `ψ = 0` below 0, the
  polynomial on `[0, 1]` and `z` above join C²), `0 ≤ ψ(z) ≤ z` and `ψ' ≥ 0` on `[0, 1]`; the
  one-sided term is bounded by `E^p ≤ e ≤ max(E^p, E^n)` for `0 ≤ w ≤ 1`, `E^p > 0`; and its
  chain rule along a motion.
* `rc-near-gate`: the chain rule of the gated term `w g (E^n − E^p)` with `g` the product of two
  quintics `W` (`smoother`, `rc-near-blend`) of the scaled closest-point offsets.
* `rc-near-foldgate` (WP16c, SPEC §17r): the fold-gated weight `W = 1 − (1 − w_c)(1 − w G)` is
  in `[0, 1]` (so the one-sided bounds hold with `W` in place of `w`) and its chain rule. The
  gradient of `cos φ` (the normalized dot product of the two edge vectors) is the standard one and
  is checked by central differences (`tests/wp16c_band.rs`).

Not machine-checked (argued in SPEC §17p): that the composition with the clamped closest-point
parameters stays C¹ (the parameters are only Lipschitz where a clamp switches, and `W' = 0` at the
clamped value 0).
-/

noncomputable section

namespace DdgProofs

namespace RC

/-! ## `rc-near-onesided` -/

/-- The polynomial piece `ψ(z) = 6z³ − 8z⁴ + 3z⁵` of the smooth positive part. -/
def posPartPoly (z : ℝ) : ℝ := 6 * z ^ 3 - 8 * z ^ 4 + 3 * z ^ 5

/-- The smooth positive part: `0` for `z ≤ 0`, `ψ(z)` on `(0, 1)`, `z` for `z ≥ 1`. -/
def smoothPos (z : ℝ) : ℝ := if z ≤ 0 then 0 else if 1 ≤ z then z else posPartPoly z

/-- **`ψ' = 18z² − 32z³ + 15z⁴`.**

[formula: rc-near-onesided] -/
theorem posPartPoly_hasDerivAt (z : ℝ) :
    HasDerivAt posPartPoly (18 * z ^ 2 - 32 * z ^ 3 + 15 * z ^ 4) z := by
  have h := (((hasDerivAt_pow 3 z).const_mul 6).sub ((hasDerivAt_pow 4 z).const_mul 8)).add
    ((hasDerivAt_pow 5 z).const_mul 3)
  exact h.congr_deriv (by norm_num; ring)

/-- **`ψ'' = 36z − 96z² + 60z³`.**

[formula: rc-near-onesided] -/
theorem posPartPoly_deriv_hasDerivAt (z : ℝ) :
    HasDerivAt (fun z => 18 * z ^ 2 - 32 * z ^ 3 + 15 * z ^ 4) (36 * z - 96 * z ^ 2 + 60 * z ^ 3)
      z := by
  have h := (((hasDerivAt_pow 2 z).const_mul 18).sub ((hasDerivAt_pow 3 z).const_mul 32)).add
    ((hasDerivAt_pow 4 z).const_mul 15)
  exact h.congr_deriv (by norm_num; ring)

/-- **The junctions are C²:** `ψ(0) = ψ'(0) = ψ''(0) = 0` (the zero branch) and
`ψ(1) = ψ'(1) = 1`, `ψ''(1) = 0` (the branch `z`).

[formula: rc-near-onesided] -/
theorem posPartPoly_junctions :
    posPartPoly 0 = 0 ∧ (18 * (0 : ℝ) ^ 2 - 32 * 0 ^ 3 + 15 * 0 ^ 4) = 0 ∧
      (36 * (0 : ℝ) - 96 * 0 ^ 2 + 60 * 0 ^ 3) = 0 ∧
    posPartPoly 1 = 1 ∧ (18 * (1 : ℝ) ^ 2 - 32 * 1 ^ 3 + 15 * 1 ^ 4) = 1 ∧
      (36 * (1 : ℝ) - 96 * 1 ^ 2 + 60 * 1 ^ 3) = 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> norm_num [posPartPoly]

/-- **`0 ≤ ψ(z) ≤ z` on `[0, 1]`:** `ψ = z³ (3(z − 4/3)² + 2/3)` and
`z − ψ = z (1 − z)³ (1 + 3z)`.

[formula: rc-near-onesided] -/
theorem posPartPoly_bounds {z : ℝ} (h0 : 0 ≤ z) (h1 : z ≤ 1) :
    0 ≤ posPartPoly z ∧ posPartPoly z ≤ z := by
  have e1 : posPartPoly z = z ^ 3 * (3 * (z - 4 / 3) ^ 2 + 2 / 3) := by unfold posPartPoly; ring
  have e2 : z - posPartPoly z = z * (1 - z) ^ 3 * (1 + 3 * z) := by unfold posPartPoly; ring
  refine ⟨?_, ?_⟩
  · rw [e1]; positivity
  · have : 0 ≤ z * (1 - z) ^ 3 * (1 + 3 * z) := by
      have : 0 ≤ 1 - z := by linarith
      positivity
    linarith

/-- **`ψ' ≥ 0`:** `ψ' = z² (15(z − 16/15)² + 14/15)`.

[formula: rc-near-onesided] -/
theorem posPartPoly_deriv_nonneg (z : ℝ) : 0 ≤ 18 * z ^ 2 - 32 * z ^ 3 + 15 * z ^ 4 := by
  have e : 18 * z ^ 2 - 32 * z ^ 3 + 15 * z ^ 4 = z ^ 2 * (15 * (z - 16 / 15) ^ 2 + 14 / 15) := by
    ring
  rw [e]; positivity

/-- **`0 ≤ smoothPos z ≤ max z 0`.**

[formula: rc-near-onesided] -/
theorem smoothPos_bounds (z : ℝ) : 0 ≤ smoothPos z ∧ smoothPos z ≤ max z 0 := by
  unfold smoothPos
  split_ifs with h0 h1
  · exact ⟨le_rfl, le_max_right _ _⟩
  · exact ⟨by linarith, le_max_left _ _⟩
  · have hb := posPartPoly_bounds (le_of_lt (not_le.mp h0)) (le_of_lt (not_le.mp h1))
    exact ⟨hb.1, hb.2.trans (le_max_left _ _)⟩

/-- **The one-sided term never lowers the paper's pair and never exceeds the larger of the two
integrals:** for `E^p > 0` and `0 ≤ w ≤ 1`,
`E^p ≤ E^p + w E^p ψ(E^n/E^p − 1) ≤ max(E^p, E^n)`.

[formula: rc-near-onesided] -/
theorem onesided_bounds {ep en w : ℝ} (hep : 0 < ep) (hw0 : 0 ≤ w) (hw1 : w ≤ 1) :
    ep ≤ ep + w * ep * smoothPos (en / ep - 1) ∧
      ep + w * ep * smoothPos (en / ep - 1) ≤ max ep en := by
  obtain ⟨h0, h1⟩ := smoothPos_bounds (en / ep - 1)
  refine ⟨by have := mul_nonneg (mul_nonneg hw0 hep.le) h0; linarith, ?_⟩
  have hz : ep * max (en / ep - 1) 0 = max (en - ep) 0 := by
    rw [mul_max_of_nonneg _ _ hep.le, mul_zero, mul_sub, mul_div_cancel₀ _ hep.ne', mul_one]
  have : w * ep * smoothPos (en / ep - 1) ≤ max (en - ep) 0 := by
    calc w * ep * smoothPos (en / ep - 1) ≤ 1 * ep * max (en / ep - 1) 0 := by
          apply mul_le_mul (mul_le_mul_of_nonneg_right hw1 hep.le) h1 h0 (by positivity)
      _ = max (en - ep) 0 := by rw [one_mul, hz]
  rcases le_total en ep with h | h
  · rw [max_eq_right (by linarith : en - ep ≤ 0)] at this
    rw [max_eq_left h]; linarith
  · rw [max_eq_left (by linarith : 0 ≤ en - ep)] at this
    rw [max_eq_right h]; linarith

/-- **The chain rule of the one-sided correction** `c = w(r) E^p ψ(E^n/E^p − 1)` along a motion
`s` (with `ψ` differentiable at the argument, derivative `ψ'`):
`c' = w' r' E^p ψ + w (ψ E^p' + ψ' (E^n' − (E^n/E^p) E^p'))`.

[formula: rc-near-onesided] -/
theorem onesided_hasDerivAt {ep en r : ℝ → ℝ} {w ψ : ℝ → ℝ} {ep' en' r' w' ψ' s : ℝ}
    (hp : HasDerivAt ep ep' s) (hn : HasDerivAt en en' s) (hr : HasDerivAt r r' s)
    (hw : HasDerivAt w w' (r s)) (hψ : HasDerivAt ψ ψ' (en s / ep s - 1)) (hep : ep s ≠ 0) :
    HasDerivAt (fun s => w (r s) * (ep s * ψ (en s / ep s - 1)))
      (w' * r' * (ep s * ψ (en s / ep s - 1)) +
        w (r s) * (ψ (en s / ep s - 1) * ep' +
          ψ' * (en' - en s / ep s * ep'))) s := by
  have hq : HasDerivAt (fun s => en s / ep s - 1) ((en' * ep s - en s * ep') / ep s ^ 2) s :=
    (hn.div hp hep).sub_const 1
  have h := (hw.comp s hr).mul (hp.mul (hψ.comp s hq))
  have key : ep s * ((en' * ep s - en s * ep') / ep s ^ 2) = en' - en s / ep s * ep' := by
    field_simp
  refine h.congr_deriv ?_
  simp only [Function.comp, Pi.mul_apply]
  generalize ψ (en s / ep s - 1) = q
  rw [show ep s * (ψ' * ((en' * ep s - en s * ep') / ep s ^ 2)) =
    ψ' * (ep s * ((en' * ep s - en s * ep') / ep s ^ 2)) by ring, key]
  ring

/-! ## `rc-near-gate` -/

/-- **The chain rule of the gated correction** `c = w(r) G(a/σ) G(b/σ) (E^n − E^p)`, `G` the
quintic `W`, along a motion `s` (`a`, `b` the offsets of the closest points from the ends of the
joining edge, differentiable at `s`).

[formula: rc-near-gate] -/
theorem gated_hasDerivAt {ep en r a b : ℝ → ℝ} {w : ℝ → ℝ} {ep' en' r' w' a' b' s σ : ℝ}
    (hp : HasDerivAt ep ep' s) (hn : HasDerivAt en en' s) (hr : HasDerivAt r r' s)
    (hw : HasDerivAt w w' (r s)) (ha : HasDerivAt a a' s) (hb : HasDerivAt b b' s) :
    HasDerivAt (fun s => w (r s) * (smoother (a s / σ) * smoother (b s / σ)) * (en s - ep s))
      (w' * r' * (smoother (a s / σ) * smoother (b s / σ)) * (en s - ep s) +
        w (r s) * (30 * (a s / σ) ^ 2 * (1 - a s / σ) ^ 2 * (a' / σ) * smoother (b s / σ) +
          smoother (a s / σ) * (30 * (b s / σ) ^ 2 * (1 - b s / σ) ^ 2 * (b' / σ))) *
          (en s - ep s) +
        w (r s) * (smoother (a s / σ) * smoother (b s / σ)) * (en' - ep')) s := by
  have hga := (smoother_hasDerivAt (a s / σ)).comp s (ha.div_const σ)
  have hgb := (smoother_hasDerivAt (b s / σ)).comp s (hb.div_const σ)
  have h := ((hw.comp s hr).mul (hga.mul hgb)).mul (hn.sub hp)
  refine h.congr_deriv ?_
  simp only [Function.comp, Pi.mul_apply, Pi.sub_apply]
  ring

/-- **The either gate** `g = 1 − (1 − G(a/σ))(1 − G(b/σ))` along a motion `s`:
`g' = (1 − G(b/σ)) G'(a/σ) a'/σ + (1 − G(a/σ)) G'(b/σ) b'/σ`.

[formula: rc-near-gate] -/
theorem either_gate_hasDerivAt {a b : ℝ → ℝ} {a' b' s σ : ℝ} (ha : HasDerivAt a a' s)
    (hb : HasDerivAt b b' s) :
    HasDerivAt (fun s => 1 - (1 - smoother (a s / σ)) * (1 - smoother (b s / σ)))
      ((1 - smoother (b s / σ)) * (30 * (a s / σ) ^ 2 * (1 - a s / σ) ^ 2 * (a' / σ)) +
        (1 - smoother (a s / σ)) * (30 * (b s / σ) ^ 2 * (1 - b s / σ) ^ 2 * (b' / σ))) s := by
  have hga := (smoother_hasDerivAt (a s / σ)).comp s (ha.div_const σ)
  have hgb := (smoother_hasDerivAt (b s / σ)).comp s (hb.div_const σ)
  have h := ((hga.const_sub 1).mul (hgb.const_sub 1)).const_sub 1
  refine h.congr_deriv ?_
  simp only [Function.comp]
  ring

/-- **The either gate is a weight that shuts only at a corner:** for `x, y ∈ [0, 1]`,
`0 ≤ 1 − (1 − x)(1 − y) ≤ 1`, and it is `0` iff `x = y = 0`.

[formula: rc-near-gate] -/
theorem either_gate_bounds {x y : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) (hy0 : 0 ≤ y) (hy1 : y ≤ 1) :
    0 ≤ 1 - (1 - x) * (1 - y) ∧ 1 - (1 - x) * (1 - y) ≤ 1 ∧
      (1 - (1 - x) * (1 - y) = 0 ↔ x = 0 ∧ y = 0) := by
  have h1 : 0 ≤ (1 - x) * (1 - y) := mul_nonneg (by linarith) (by linarith)
  have h2 : (1 - x) * (1 - y) ≤ 1 := by nlinarith
  refine ⟨by linarith, by linarith, ⟨fun h => ?_, fun ⟨hx, hy⟩ => by simp [hx, hy]⟩⟩
  have hxy : x + y = x * y := by linarith
  have hx : x = 0 := by nlinarith
  exact ⟨hx, by subst hx; linarith⟩

/-! ## `rc-near-foldgate` (WP16c, SPEC §17r; an option, off by default) -/

/-- **The fold-gated weight is a weight:** for `w_c, w, G ∈ [0, 1]`,
`0 ≤ W = 1 − (1 − w_c)(1 − w G) ≤ 1`, so with `onesided_bounds` the fold-gated term keeps
`E^p ≤ e ≤ max(E^p, E^n)`; and `W = w` wherever the contact weight vanishes and the gate is open.

[formula: rc-near-foldgate] -/
theorem foldgate_bounds {wc w g : ℝ} (hc0 : 0 ≤ wc) (hc1 : wc ≤ 1) (hw0 : 0 ≤ w) (hw1 : w ≤ 1)
    (hg0 : 0 ≤ g) (hg1 : g ≤ 1) :
    0 ≤ 1 - (1 - wc) * (1 - w * g) ∧ 1 - (1 - wc) * (1 - w * g) ≤ 1 ∧
      1 - (1 - 0) * (1 - w * 1) = w := by
  have h0 : 0 ≤ w * g := mul_nonneg hw0 hg0
  have h1 : w * g ≤ 1 := by nlinarith
  obtain ⟨a, b, _⟩ := either_gate_bounds hc0 hc1 h0 h1
  exact ⟨a, b, by ring⟩

/-- **The chain rule of the fold-gated weight** `W = 1 − (1 − w_c(r))(1 − w(r) G((cos φ − c₀)/δ))`
along a motion `s` (`G` the quintic `smoother`, `δ = cos_open − cos_shut`):
`W' = (1 − w G) w_c' r' + (1 − w_c)(w' r' G + w G' cos φ'/δ)`.

[formula: rc-near-foldgate] -/
theorem foldgate_hasDerivAt {r c : ℝ → ℝ} {wc w : ℝ → ℝ} {r' c' wc' w' s c₀ δ : ℝ}
    (hr : HasDerivAt r r' s) (hc : HasDerivAt c c' s) (hwc : HasDerivAt wc wc' (r s))
    (hw : HasDerivAt w w' (r s)) :
    HasDerivAt
      (fun s => 1 - (1 - wc (r s)) * (1 - w (r s) * smoother ((c s - c₀) / δ)))
      ((1 - w (r s) * smoother ((c s - c₀) / δ)) * (wc' * r') +
        (1 - wc (r s)) * (w' * r' * smoother ((c s - c₀) / δ) +
          w (r s) * (30 * ((c s - c₀) / δ) ^ 2 * (1 - (c s - c₀) / δ) ^ 2 * (c' / δ)))) s := by
  have hg := (smoother_hasDerivAt ((c s - c₀) / δ)).comp s ((hc.sub_const c₀).div_const δ)
  have h := (((hwc.comp s hr).const_sub 1).mul (((hw.comp s hr).mul hg).const_sub 1)).const_sub 1
  refine h.congr_deriv ?_
  simp only [Function.comp, Pi.mul_apply]
  generalize smoother ((c s - c₀) / δ) = q
  ring

end RC

end DdgProofs
