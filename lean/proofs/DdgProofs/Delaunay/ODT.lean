import DDG.Intrinsic.Move
import DDG.FEM.Tiling

/-!
# Optimal Delaunay triangulations (Chen–Holst): the vertex update and the interpolation error

`crates/ddg-intrinsic/SPEC.md` §4 (Optimal Delaunay triangulation), N20, N21; the paper
Sec. 4.3; L. Chen, M. Holst, *Efficient mesh optimization schemes based on Optimal Delaunay
Triangulations*, CMAME 200 (2011), Lemma 4.1, Lemma 4.2 (Eq. 4.3–4.6), Lemma 4.3, Eq. 4.13,
Sec. 4.2.1 of the preprint (barycentres at the boundary).

Plane points are pairs `ℝ × ℝ` with `sda` (doubled signed area), `sqDist2`, `dist2` of
`DDG/Geometry/Plane.lean` and `nsq`, `rot`, `toDual` of `DDG/Geometry/PlaneLemmas.lean`;
the integrals are those of `DDG/FEM/Elements.lean` and `Tiling.lean`.

* **`it-odt-update`** — `weightedBarycenter` is `A (a + b + c)/3` in the frame of the
  circumcentre vector of `DDG/Intrinsic/Move.lean`
  (`weightedBarycenter_eq`). For a vertex `x` with a closed counter-clockwise link
  `p 0, …, p (n−1), p n = p 0` and `|T_k| = sda x (p k) (p (k+1)) / 2`, the energy of Chen and
  Holst (Lemma 4.1 with `d = 2`, the star's part) `S(x) = Σ_k (|x|² + |p_k|² + |p_{k+1}|²) |T_k| / 3`
  has the derivative `(2/3) Σ_k (|T_k| x − |T_k| c_k)` (Eq. 4.5; `hasFDerivAt_odtEnergy`), with
  `|T_k| c_k` the polynomial `circumVec`, which divided by `|T_k|` is equidistant from the three
  corners (`circumVec_equidistant`) and in the frame of a triangle is
  `weightedCircumcenter` (`circumVec_layout`). `S` is exactly quadratic, with Hessian
  `(2/3) |ω| I` (Lemma 4.3; `odtEnergy_add`); hence the update of Eq. 4.13,
  `odtDisplacement` of `w = Σ_k |T_k| (c_k − x) = −(3/2) ∇S` by the area `|ω|` of the star,
  changes the energy by exactly `−(3 |∇S|²/(4|ω|)) α (2 − α)`: a descent step for every
  `0 < α < 2`, the minimum of `S` for `α = 1` (`odtEnergy_step`, `odtEnergy_step_le`).
* **`it-odt-interpolation-error`** — Chen–Holst Lemma 4.1 on one triangle, for `u = |x|²`:
  `∫_T |u_I − u| = (A/3) Σ|pᵢ|² − ∫_T |x|² = (A/12)(ℓ₀₁² + ℓ₁₂² + ℓ₂₀²)`
  (`integral_abs_odt_error`, `integral_odt_error_eq_sides`), and the star energy `odtEnergy` above is `Σ_k ∫_{T_k} u_I` (`odtEnergy_eq_sum_integral`).
* **`it-odt-energy`** — Chen–Holst Lemma 4.1 for a triangulation (`d = 2`, uniform density): for
  any `u_I` that is the linear interpolant of `u = |x|²` on every triangle,
  `E = ∫_Ω |u_I − u| = Σ_k (|T_k|/3) Σ_{i} |q_{k,i}|² − ∫_Ω |x|²`
  (`IsTiling.integral_abs_odt_error`), `= Σ_k (|T_k|/12)(ℓ₀₁² + ℓ₁₂² + ℓ₂₀²)` (a function of the
  edge lengths, `IsTiling.integral_abs_odt_error_eq_sides`), and with the vertices of a mesh
  `= (1/3) Σ_i |x_i|² |ω_i| − ∫_Ω |x|²`, `|ω_i|` the area of the star of `x_i`
  (`IsTiling.integral_abs_odt_error_vertices`). For a mesh whose triangles are a fan around `x`
  (counter-clockwise) plus fixed triangles, `E(x) = S(x) + C − ∫_Ω |y|²` with `S = odtEnergy` above and `C` independent of `x` (`odt_error_star`); so whenever the moved fan
  still tiles the same region with the others, the update of Eq. 4.13 lowers the interpolation
  error itself, by exactly what it lowers `S`, and `α = 1` minimizes `E` over all positions with
  that property (`odt_error_step`, `odt_error_step_le`).

Not covered here: what the minimizing triangulations are and the boundary rule as a remedy
(`CITED.md`, `it-odt-delaunay-optimality`). The identities are about a planar star with signed areas; that the star of an intrinsic vertex
unfolds into such a star is SPEC N19.

(Library step S2b: assembled from `DDG/Intrinsic/Move.lean` § `it-odt-update`,
`DDG/FEM/Elements.lean` § `it-odt-interpolation-error` and
`DDG/FEM/Tiling.lean` § `it-odt-energy`; the declarations keep their
namespaces `DDG.Intrinsic` and `DDG.TriangleIntegral` until the renames of step S4.)

Own definitions: `circumVec`, `odtEnergy`, `odtGradient`, `starArea`, `odtInterp`, `starAreaOf`,
`fan`.
-/

noncomputable section

namespace DDG

namespace Intrinsic

open Real Finset

/-! ## `it-odt-update`: the barycentre and the displacement -/

/-- The area `A` of the triangle `abc` times the vector from `a` to its barycentre, in the frame
of `a` with `b = (ℓ_ab, 0)` and `c = (c_x, h)` above: `(A (ℓ_ab + c_x)/3, A h/3)`,
`c_x = diamondX ℓ_ab ℓ_ac ℓ_bc`, `h = 2A/ℓ_ab`. -/
def weightedBarycenter (lab lac lbc A : ℝ) : ℝ × ℝ :=
  (A * (lab + diamondX lab lac lbc) / 3, A * (2 * A / lab) / 3)

/-- The update of the optimal Delaunay triangulation at a vertex (Chen and Holst 2011,
Eq. 4.13): the sum `w = Σ_T |T| (c_T − x)` over the triangles of the star, divided by the area
`total` of the star and scaled by the step `step`. -/
def odtDisplacement (w : ℝ × ℝ) (total step : ℝ) : ℝ × ℝ :=
  (step * w.1 / total, step * w.2 / total)

/-- **The weighted barycentre** is `A (a + b + c)/3` with `a = 0`, `b = (ℓ_ab, 0)` and
`c = (c_x, 2A/ℓ_ab)`, `c_x = diamondX ℓ_ab ℓ_ac ℓ_bc` (the frame of the circumcentre).

[formula: it-odt-update] -/
theorem weightedBarycenter_eq (lab lac lbc A : ℝ) :
    weightedBarycenter lab lac lbc A =
      A • ((1 / 3 : ℝ) • ((0, 0) + (lab, 0) + (diamondX lab lac lbc, 2 * A / lab))) := by
  simp only [weightedBarycenter, Prod.smul_mk, Prod.mk_add_mk,
    smul_eq_mul]
  ext <;> dsimp only <;> ring

theorem odtDisplacement_eq (w : ℝ × ℝ) (total step : ℝ) :
    odtDisplacement w total step = (step / total) • w := by
  simp only [odtDisplacement, Prod.smul_def, smul_eq_mul]
  ext <;> dsimp only <;> ring

/-! ## `it-odt-update`: the circumcentre of a triangle of the star, without division -/

/-- `|T| c_T` for the triangle `(x, b, c)`, `|T| = sda x b c / 2` and `c_T` its circumcentre, as a
polynomial: `¼ [|x|² rot(c − b) − |b|² rot(c − x) + |c|² rot(b − x)]`. -/
def circumVec (x b c : ℝ × ℝ) : ℝ × ℝ :=
  (1 / 4 : ℝ) • (nsq x • rot (c - b) - nsq b • rot (c - x) + nsq c • rot (b - x))

theorem circumVec_fst (x b c : ℝ × ℝ) :
    (circumVec x b c).1 = (nsq x * (-(c.2 - b.2)) - nsq b * (-(c.2 - x.2)) +
      nsq c * (-(b.2 - x.2))) / 4 := by
  simp only [circumVec, rot, Prod.smul_fst, Prod.fst_add, Prod.fst_sub, Prod.snd_sub,
    smul_eq_mul]
  ring

theorem circumVec_snd (x b c : ℝ × ℝ) :
    (circumVec x b c).2 = (nsq x * (c.1 - b.1) - nsq b * (c.1 - x.1) + nsq c * (b.1 - x.1)) / 4 := by
  simp only [circumVec, rot, Prod.smul_snd, Prod.snd_add, Prod.fst_sub, Prod.snd_sub,
    smul_eq_mul]
  ring

/-- **`circumVec / |T|` is the circumcentre**: for a triangle `(x, b, c)` with
`|T| = sda x b c / 2 ≠ 0`, the point `o = circumVec x b c / |T|` is equidistant from the three
corners.

[formula: it-odt-update] -/
theorem circumVec_equidistant {x b c : ℝ × ℝ} (hT : sda x b c ≠ 0) :
    let T := sda x b c / 2
    let o : ℝ × ℝ := ((circumVec x b c).1 / T, (circumVec x b c).2 / T)
    sqDist2 x o = sqDist2 b o ∧ sqDist2 x o = sqDist2 c o := by
  intro T o
  have hT2 : T ≠ 0 := by simp only [T]; exact div_ne_zero hT two_ne_zero
  have h1 : o.1 * T = (circumVec x b c).1 := by simp only [o]; field_simp
  have h2 : o.2 * T = (circumVec x b c).2 := by simp only [o]; field_simp
  rw [circumVec_fst] at h1
  rw [circumVec_snd] at h2
  have hTd : T * 2 = sda x b c := by simp only [T]; ring
  clear_value o T
  unfold sqDist2
  unfold nsq at h1 h2
  unfold sda at hTd
  constructor
  · apply mul_right_cancel₀ hT2
    linear_combination (-2 * (x.1 - b.1)) * h1 + (-2 * (x.2 - b.2)) * h2 +
      ((x.1 ^ 2 + x.2 ^ 2 - b.1 ^ 2 - b.2 ^ 2) / 2) * hTd
  · apply mul_right_cancel₀ hT2
    linear_combination (-2 * (x.1 - c.1)) * h1 + (-2 * (x.2 - c.2)) * h2 +
      ((x.1 ^ 2 + x.2 ^ 2 - c.1 ^ 2 - c.2 ^ 2) / 2) * hTd

/-- **In the frame of a triangle it is `weightedCircumcenter`**: with `x = 0`,
`b = (ℓ_ab, 0)`, `c = (c_x, h)`, `c_x = (ℓ_ab² + ℓ_ac² − ℓ_bc²)/2ℓ_ab`, `|c| = ℓ_ac` and
`A = ℓ_ab h / 2`.

[formula: it-odt-update] -/
theorem circumVec_layout {lab lac lbc cx h : ℝ} (hab : lab ≠ 0)
    (hcx : cx = (lab ^ 2 + lac ^ 2 - lbc ^ 2) / (2 * lab)) (hc : cx ^ 2 + h ^ 2 = lac ^ 2) :
    circumVec (0, 0) (lab, 0) (cx, h) = weightedCircumcenter lab lac lbc (lab * h / 2) := by
  have hcx' : cx * (2 * lab) = lab ^ 2 + lac ^ 2 - lbc ^ 2 := by rw [hcx]; field_simp
  rw [weightedCircumcenter_eq]
  ext
  · rw [circumVec_fst]; unfold nsq; dsimp only; ring
  · rw [circumVec_snd]; unfold nsq; dsimp only
    linear_combination (lab / 4) * hc - (lab / 8) * hcx'

/-! ## `it-odt-update`: the energy of the star and its derivative -/

/-- The energy of the star (Chen and Holst, Lemma 4.1 with `d = 2`, the part that depends on the
centre `x`): `S(x) = Σ_k (|x|² + |p_k|² + |p_{k+1}|²) |T_k| / 3`, `|T_k| = sda x p_k p_{k+1} / 2`. -/
def odtEnergy (p : ℕ → ℝ × ℝ) (n : ℕ) (x : ℝ × ℝ) : ℝ :=
  ∑ k ∈ range n, (nsq x + (nsq (p k) + nsq (p (k + 1)))) * (sda x (p k) (p (k + 1)) / 2) / 3

/-- The gradient of Chen and Holst's Eq. 4.5 (`d = 2`): `(2/3) Σ_k (|T_k| x − |T_k| c_k)`. -/
def odtGradient (p : ℕ → ℝ × ℝ) (n : ℕ) (x : ℝ × ℝ) : ℝ × ℝ :=
  (2 / 3 : ℝ) • ∑ k ∈ range n, ((sda x (p k) (p (k + 1)) / 2) • x - circumVec x (p k) (p (k + 1)))

/-- The area of the star (signed): `|ω| = Σ_k |T_k|`. -/
def starArea (p : ℕ → ℝ × ℝ) (n : ℕ) (x : ℝ × ℝ) : ℝ :=
  ∑ k ∈ range n, sda x (p k) (p (k + 1)) / 2

/-- One term: the derivative of `(|y|² + B) sda(y, b, c) / 6`. -/
def triGrad (B : ℝ) (b c x : ℝ × ℝ) : ℝ × ℝ :=
  ((2 * x.1 * sda x b c + (nsq x + B) * (b.2 - c.2)) / 6,
    (2 * x.2 * sda x b c + (nsq x + B) * (c.1 - b.1)) / 6)

theorem hasFDerivAt_triTerm (B : ℝ) (b c x : ℝ × ℝ) :
    HasFDerivAt (fun y : ℝ × ℝ => (nsq y + B) * (sda y b c / 2) / 3) (toDual (triGrad B b c x)) x := by
  have h1 : HasFDerivAt (fun y : ℝ × ℝ => y.1) (ContinuousLinearMap.fst ℝ ℝ ℝ) x := hasFDerivAt_fst
  have h2 : HasFDerivAt (fun y : ℝ × ℝ => y.2) (ContinuousLinearMap.snd ℝ ℝ ℝ) x := hasFDerivAt_snd
  have hn := ((h1.pow 2).add (h2.pow 2)).add_const B
  have hs := ((h1.const_sub b.1).mul (h2.const_sub c.2)).sub
    ((h2.const_sub b.2).mul (h1.const_sub c.1))
  have hf := (hn.mul hs).const_mul (1 / 6 : ℝ)
  convert hf using 1
  · funext y
    simp only [Pi.add_apply, Pi.mul_apply, Pi.sub_apply, nsq, sda]
    ring
  · refine ContinuousLinearMap.ext fun v => ?_
    simp only [triGrad, toDual_apply, nsq, sda, add_apply,
      smul_apply, ContinuousLinearMap.coe_fst', ContinuousLinearMap.coe_snd',
      neg_apply, sub_apply, smul_eq_mul, Pi.add_apply,
      Pi.mul_apply, Pi.sub_apply, nsmul_eq_mul]
    push_cast
    ring

/-- Per triangle (the identity behind Eq. 4.5): with `F(y) = ½ (|y|² rot(y − x) + |x|² rot(y))`,
`6 triGrad = 2 [2 (|T| x − |T| c_T)] + |x|² rot(c − b) + 2 (F(c) − F(b))`, componentwise. -/
def forceF (x y : ℝ × ℝ) : ℝ × ℝ :=
  ((nsq y * (-(y.2 - x.2)) + nsq x * (-y.2)) / 2, (nsq y * (y.1 - x.1) + nsq x * y.1) / 2)

theorem triGrad_split (b c x : ℝ × ℝ) :
    triGrad (nsq b + nsq c) b c x =
      ((2 / 3 : ℝ) * ((sda x b c / 2) * x.1 - (circumVec x b c).1) +
          (nsq x * (b.2 - c.2) / 2 + (forceF x c).1 - (forceF x b).1) / 3,
        (2 / 3 : ℝ) * ((sda x b c / 2) * x.2 - (circumVec x b c).2) +
          (nsq x * (c.1 - b.1) / 2 + (forceF x c).2 - (forceF x b).2) / 3) := by
  rw [circumVec_fst, circumVec_snd]
  unfold triGrad forceF sda nsq
  ext <;> dsimp only <;> ring

/-- The sum of the `rot (p (k+1) − p k)` around a closed link is zero, componentwise. -/
theorem sum_link_sides {p : ℕ → ℝ × ℝ} {n : ℕ} (hp : p n = p 0) :
    ∑ k ∈ range n, ((p k).2 - (p (k + 1)).2) = 0 ∧
      ∑ k ∈ range n, ((p (k + 1)).1 - (p k).1) = 0 := by
  constructor
  · rw [sum_range_sub' (fun k => (p k).2), hp, sub_self]
  · rw [sum_range_sub (fun k => (p k).1), hp, sub_self]

/-- **The sum of the terms is the gradient of Eq. 4.5** on a closed link: the forces `F` telescope
and `Σ rot(p_{k+1} − p_k) = 0`. -/
theorem sum_triGrad {p : ℕ → ℝ × ℝ} {n : ℕ} (hp : p n = p 0) (x : ℝ × ℝ) :
    ∑ k ∈ range n, triGrad (nsq (p k) + nsq (p (k + 1))) (p k) (p (k + 1)) x =
      odtGradient p n x := by
  obtain ⟨s2, s1⟩ := sum_link_sides hp
  have tF1 : ∑ k ∈ range n, ((forceF x (p (k + 1))).1 - (forceF x (p k)).1) = 0 := by
    rw [sum_range_sub (fun k => (forceF x (p k)).1), hp, sub_self]
  have tF2 : ∑ k ∈ range n, ((forceF x (p (k + 1))).2 - (forceF x (p k)).2) = 0 := by
    rw [sum_range_sub (fun k => (forceF x (p k)).2), hp, sub_self]
  simp only [triGrad_split]
  unfold odtGradient
  ext
  · simp only [Prod.fst_sum, Prod.smul_fst, Prod.fst_sub, smul_eq_mul, sum_add_distrib, ← sum_div,
      ← mul_sum]
    have e : ∑ k ∈ range n, (nsq x * ((p k).2 - (p (k + 1)).2) / 2 + (forceF x (p (k + 1))).1 -
        (forceF x (p k)).1) = nsq x / 2 * ∑ k ∈ range n, ((p k).2 - (p (k + 1)).2) +
          ∑ k ∈ range n, ((forceF x (p (k + 1))).1 - (forceF x (p k)).1) := by
      rw [mul_sum, ← sum_add_distrib]
      exact sum_congr rfl fun k _ => by ring
    rw [e, s2, tF1]
    simp
  · simp only [Prod.snd_sum, Prod.smul_snd, Prod.snd_sub, smul_eq_mul, sum_add_distrib, ← sum_div,
      ← mul_sum]
    have e : ∑ k ∈ range n, (nsq x * ((p (k + 1)).1 - (p k).1) / 2 + (forceF x (p (k + 1))).2 -
        (forceF x (p k)).2) = nsq x / 2 * ∑ k ∈ range n, ((p (k + 1)).1 - (p k).1) +
          ∑ k ∈ range n, ((forceF x (p (k + 1))).2 - (forceF x (p k)).2) := by
      rw [mul_sum, ← sum_add_distrib]
      exact sum_congr rfl fun k _ => by ring
    rw [e, s1, tF2]
    simp

/-- **The gradient of the energy of the star (Chen and Holst, Lemma 4.2, Eq. 4.5, `d = 2`)**: on a
closed link `p 0, …, p n = p 0`, the energy `S` has at every `x` the derivative
`v ↦ ⟨(2/3) Σ_k (|T_k| x − |T_k| c_k), v⟩`, `|T_k| c_k = circumVec x (p k) (p (k+1))`.

[formula: it-odt-update] -/
theorem hasFDerivAt_odtEnergy {p : ℕ → ℝ × ℝ} {n : ℕ} (hp : p n = p 0) (x : ℝ × ℝ) :
    HasFDerivAt (odtEnergy p n) (toDual (odtGradient p n x)) x := by
  have hsum := HasFDerivAt.sum (u := range n) (x := x)
    (fun k _ => hasFDerivAt_triTerm (nsq (p k) + nsq (p (k + 1))) (p k) (p (k + 1)) x)
  have e : (∑ k ∈ range n, fun y => (nsq y + (nsq (p k) + nsq (p (k + 1)))) *
      (sda y (p k) (p (k + 1)) / 2) / 3) = odtEnergy p n := by
    funext y
    rw [Finset.sum_apply]
    rfl
  have e2 : ∑ k ∈ range n, toDual (triGrad (nsq (p k) + nsq (p (k + 1))) (p k) (p (k + 1)) x) =
      toDual (odtGradient p n x) := by
    rw [← sum_triGrad hp x]
    refine ContinuousLinearMap.ext fun v => ?_
    simp only [_root_.sum_apply, toDual_apply, Prod.fst_sum,
      Prod.snd_sum, sum_mul, ← sum_add_distrib]
  rw [← e, ← e2]
  exact hsum

/-- The same as directional derivatives: `t ↦ S(x + t e)` has the derivative `⟨G, e⟩` at `0`.

[formula: it-odt-update] -/
theorem hasDerivAt_odtEnergy_line {p : ℕ → ℝ × ℝ} {n : ℕ} (hp : p n = p 0) (x e : ℝ × ℝ) :
    HasDerivAt (fun t : ℝ => odtEnergy p n (x + t • e))
      ((odtGradient p n x).1 * e.1 + (odtGradient p n x).2 * e.2) 0 := by
  have hl : HasDerivAt (fun t : ℝ => x + t • e) e 0 := by
    simpa using ((hasDerivAt_id (0 : ℝ)).smul_const e).const_add x
  have hS := hasFDerivAt_odtEnergy hp x
  have hS' : HasFDerivAt (odtEnergy p n) (toDual (odtGradient p n x)) (x + (0 : ℝ) • e) := by
    simpa using hS
  have := hS'.comp_hasDerivAt (0 : ℝ) hl
  rw [toDual_apply] at this
  exact this

/-! ## `it-odt-update`: the energy is quadratic; the step of Eq. 4.13 -/

/-- The area of the star does not depend on its centre (`sum_sda_closed_link`). -/
theorem starArea_indep {p : ℕ → ℝ × ℝ} {n : ℕ} (hp : p n = p 0) (x y : ℝ × ℝ) :
    starArea p n x = starArea p n y := by
  unfold starArea
  rw [← sum_div, ← sum_div, sum_sda_closed_link hp x y]

/-- Per triangle: the exact expansion of one term of `S` around `x`. -/
theorem triTerm_add (B : ℝ) (b c x v : ℝ × ℝ) :
    (nsq (x + v) + B) * (sda (x + v) b c / 2) / 3 =
      (nsq x + B) * (sda x b c / 2) / 3 + toDual (triGrad B b c x) v +
        sda x b c / 2 * nsq v / 3 +
        (2 * (x.1 * v.1 + x.2 * v.2) + nsq v) * ((b.2 - c.2) * v.1 + (c.1 - b.1) * v.2) / 6 := by
  rw [toDual_apply]
  unfold triGrad nsq sda
  simp only [Prod.fst_add, Prod.snd_add]
  ring

/-- **The energy of the star is exactly quadratic** (Chen and Holst, Lemma 4.3: the Hessian is
`(2/3)|ω|` times the identity): `S(x + v) = S(x) + ⟨∇S(x), v⟩ + |ω| |v|² / 3`.

[formula: it-odt-update] -/
theorem odtEnergy_add {p : ℕ → ℝ × ℝ} {n : ℕ} (hp : p n = p 0) (x v : ℝ × ℝ) :
    odtEnergy p n (x + v) =
      odtEnergy p n x + toDual (odtGradient p n x) v + starArea p n x * nsq v / 3 := by
  obtain ⟨s2, s1⟩ := sum_link_sides hp
  have hg : toDual (odtGradient p n x) v =
      ∑ k ∈ range n, toDual (triGrad (nsq (p k) + nsq (p (k + 1))) (p k) (p (k + 1)) x) v := by
    rw [← sum_triGrad hp x]
    simp only [toDual_apply, Prod.fst_sum, Prod.snd_sum, sum_mul, ← sum_add_distrib]
  have hc : ∑ k ∈ range n, (2 * (x.1 * v.1 + x.2 * v.2) + nsq v) *
      (((p k).2 - (p (k + 1)).2) * v.1 + ((p (k + 1)).1 - (p k).1) * v.2) / 6 = 0 := by
    have : ∑ k ∈ range n, (2 * (x.1 * v.1 + x.2 * v.2) + nsq v) *
        (((p k).2 - (p (k + 1)).2) * v.1 + ((p (k + 1)).1 - (p k).1) * v.2) / 6 =
        (2 * (x.1 * v.1 + x.2 * v.2) + nsq v) / 6 *
          ((∑ k ∈ range n, ((p k).2 - (p (k + 1)).2)) * v.1 +
            (∑ k ∈ range n, ((p (k + 1)).1 - (p k).1)) * v.2) := by
      rw [sum_mul, sum_mul, ← sum_add_distrib, mul_sum]
      exact sum_congr rfl fun k _ => by ring
    rw [this, s1, s2]
    ring
  unfold odtEnergy starArea
  rw [sum_congr rfl fun k _ => triTerm_add (nsq (p k) + nsq (p (k + 1))) (p k) (p (k + 1)) x v,
    sum_add_distrib, sum_add_distrib, sum_add_distrib, hc, add_zero, hg, sum_mul, sum_div]

/-- The displacement of Eq. 4.13 is `−(3α/(2|ω|))` times the gradient: with
`w = Σ_k (|T_k| c_k − |T_k| x) = −(3/2) ∇S`, `odtDisplacement w |ω| α = −(3α/(2|ω|)) ∇S`.

[formula: it-odt-update] -/
theorem odtDisplacement_eq_gradient (p : ℕ → ℝ × ℝ) (n : ℕ) (x : ℝ × ℝ) (step : ℝ) :
    let w := ∑ k ∈ range n, (circumVec x (p k) (p (k + 1)) - (sda x (p k) (p (k + 1)) / 2) • x)
    odtDisplacement w (starArea p n x) step =
      (-(3 * step / (2 * starArea p n x))) • odtGradient p n x := by
  intro w
  rw [odtDisplacement_eq]
  unfold odtGradient
  have : w = (-1 : ℝ) • ∑ k ∈ range n, ((sda x (p k) (p (k + 1)) / 2) • x -
      circumVec x (p k) (p (k + 1))) := by
    simp only [w, smul_sum, neg_one_smul, neg_sub]
  rw [this, smul_smul, smul_smul]
  congr 1
  by_cases h : starArea p n x = 0
  · simp [h]
  · field_simp

/-- **The step of Eq. 4.13 changes the energy by exactly `−(3|∇S|²/(4|ω|)) α (2 − α)`** (for a
star of nonzero area): moving `x` by `odtDisplacement w |ω| α`, `w = Σ_k |T_k| (c_k − x)`.

[formula: it-odt-update] -/
theorem odtEnergy_step {p : ℕ → ℝ × ℝ} {n : ℕ} (hp : p n = p 0) (x : ℝ × ℝ)
    (hω : starArea p n x ≠ 0) (step : ℝ) :
    let w := ∑ k ∈ range n, (circumVec x (p k) (p (k + 1)) - (sda x (p k) (p (k + 1)) / 2) • x)
    odtEnergy p n (x + odtDisplacement w (starArea p n x) step) =
      odtEnergy p n x -
        3 * nsq (odtGradient p n x) / (4 * starArea p n x) * (step * (2 - step)) := by
  intro w
  rw [odtDisplacement_eq_gradient, odtEnergy_add hp, toDual_apply]
  unfold nsq
  simp only [Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  field_simp
  ring

/-- **A descent step for `0 < α < 2`** on a star of positive area, strict unless `∇S = 0`; and
`α = 1` minimizes `S` over all centres.

[formula: it-odt-update] -/
theorem odtEnergy_step_le {p : ℕ → ℝ × ℝ} {n : ℕ} (hp : p n = p 0) (x : ℝ × ℝ)
    (hω : 0 < starArea p n x) {step : ℝ} (h0 : 0 ≤ step) (h2 : step ≤ 2) :
    let w := ∑ k ∈ range n, (circumVec x (p k) (p (k + 1)) - (sda x (p k) (p (k + 1)) / 2) • x)
    odtEnergy p n (x + odtDisplacement w (starArea p n x) step) ≤ odtEnergy p n x ∧
      ∀ y, odtEnergy p n (x + odtDisplacement w (starArea p n x) 1) ≤ odtEnergy p n y := by
  intro w
  have hn : 0 ≤ nsq (odtGradient p n x) := by unfold nsq; positivity
  refine ⟨?_, fun y => ?_⟩
  · rw [odtEnergy_step hp x hω.ne' step]
    have : 0 ≤ 3 * nsq (odtGradient p n x) / (4 * starArea p n x) * (step * (2 - step)) := by
      apply mul_nonneg
      · positivity
      · exact mul_nonneg h0 (by linarith)
    linarith
  · rw [odtEnergy_step hp x hω.ne' 1]
    have hy := odtEnergy_add hp x (y - x)
    rw [add_sub_cancel] at hy
    rw [hy, toDual_apply]
    set G := odtGradient p n x
    set ω := starArea p n x
    -- `⟨G, v⟩ + ω|v|²/3 ≥ −3|G|²/(4ω)`: complete the square
    have hsq : 0 ≤ nsq ((3 / (2 * ω)) • G + (y - x)) := by unfold nsq; positivity
    unfold nsq at hsq ⊢
    simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul,
      Prod.fst_sub, Prod.snd_sub] at hsq ⊢
    have key : -(3 * (G.1 ^ 2 + G.2 ^ 2) / (4 * ω) * (1 * (2 - 1))) ≤
        G.1 * (y.1 - x.1) + G.2 * (y.2 - x.2) + ω * ((y.1 - x.1) ^ 2 + (y.2 - x.2) ^ 2) / 3 := by
      have e : G.1 * (y.1 - x.1) + G.2 * (y.2 - x.2) + ω * ((y.1 - x.1) ^ 2 + (y.2 - x.2) ^ 2) / 3 +
          3 * (G.1 ^ 2 + G.2 ^ 2) / (4 * ω) * (1 * (2 - 1)) =
          ω / 3 * ((3 / (2 * ω) * G.1 + (y.1 - x.1)) ^ 2 + (3 / (2 * ω) * G.2 + (y.2 - x.2)) ^ 2) := by
        field_simp
        ring
      have : 0 ≤ ω / 3 * ((3 / (2 * ω) * G.1 + (y.1 - x.1)) ^ 2 +
          (3 / (2 * ω) * G.2 + (y.2 - x.2)) ^ 2) := by positivity
      linarith
    linarith

end Intrinsic

namespace TriangleIntegral

open MeasureTheory Set Real

section ChenHolst

variable {p : Fin 3 → ℝ × ℝ}

/-! ## `it-odt-interpolation-error`: Chen–Holst Lemma 4.1 on a triangle -/

open Intrinsic in
/-- The linear interpolant of `u(x) = |x|²`. -/
def odtInterp (p : Fin 3 → ℝ × ℝ) : ℝ × ℝ → ℝ := interp p fun i => nsq (p i)

open Intrinsic in
/-- `u_I − u = Σ_{i<j} λᵢ λⱼ |pᵢ − pⱼ|²` for `u = |x|²` (everywhere in the plane). -/
theorem odtInterp_sub (hD : D p ≠ 0) (x : ℝ × ℝ) :
    odtInterp p x - nsq x = bary p 0 x * bary p 1 x * sqDist2 (p 0) (p 1) +
      bary p 1 x * bary p 2 x * sqDist2 (p 1) (p 2) + bary p 2 x * bary p 0 x * sqDist2 (p 2) (p 0) := by
  have hsum := sum_bary x hD
  have hx := sum_bary_smul x hD
  have hx1 := congrArg Prod.fst hx
  have hx2 := congrArg Prod.snd hx
  simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul] at hx1 hx2
  simp only [odtInterp, interp, Fin.sum_univ_three, nsq, sqDist2]
  set l0 := bary p 0 x
  set l1 := bary p 1 x
  set l2 := bary p 2 x
  rw [← hx1, ← hx2]
  linear_combination (-(l0 * ((p 0).1 ^ 2 + (p 0).2 ^ 2) + l1 * ((p 1).1 ^ 2 + (p 1).2 ^ 2) +
    l2 * ((p 2).1 ^ 2 + (p 2).2 ^ 2)) +
    (l0 * (p 0).1 + l1 * (p 1).1 + l2 * (p 2).1) * (0 : ℝ)) * hsum

open Intrinsic in
/-- **The interpolant of the convex `|x|²` lies above it on the triangle.** -/
theorem nsq_le_odtInterp (hD : D p ≠ 0) {x : ℝ × ℝ} (hx : x ∈ triangle p) : nsq x ≤ odtInterp p x := by
  have h := odtInterp_sub hD x
  have := bary_nonneg_of_mem hD hx
  have h0 := this 0
  have h1 := this 1
  have h2 := this 2
  have := sqDist2_nonneg (p 0) (p 1)
  have := sqDist2_nonneg (p 1) (p 2)
  have := sqDist2_nonneg (p 2) (p 0)
  have : 0 ≤ odtInterp p x - nsq x := by rw [h]; positivity
  linarith

theorem continuous_odtInterp (p : Fin 3 → ℝ × ℝ) : Continuous (odtInterp p) :=
  continuous_interp p _

theorem continuous_nsq : Continuous Intrinsic.nsq := by
  unfold Intrinsic.nsq; fun_prop

open Intrinsic in
/-- **The interpolation error of `|x|²` on a triangle** in closed form:
`∫_T |u_I − u| = (A/12)(ℓ₁₂² + ℓ₂₀² + ℓ₀₁²)`.

[formula: it-odt-interpolation-error] -/
theorem integral_odt_error_eq_sides (hD : D p ≠ 0) :
    ∫ x in triangle p, |odtInterp p x - nsq x| =
      area p / 12 * (sideLen p 0 ^ 2 + sideLen p 1 ^ 2 + sideLen p 2 ^ 2) := by
  rw [setIntegral_congr_fun (measurableSet_triangle p) (fun x hx =>
    abs_of_nonneg (sub_nonneg.2 (nsq_le_odtInterp hD hx)))]
  simp_rw [odtInterp_sub hD]
  have hi : ∀ (i j : Fin 3) (c : ℝ), Integrable (fun x => bary p i x * bary p j x * c)
      (volume.restrict (triangle p)) := fun i j c =>
    integrableOn_triangle p (((continuous_bary p i).mul (continuous_bary p j)).mul continuous_const)
  rw [integral_add ?_ (hi 2 0 _), integral_add (hi 0 1 _) (hi 1 2 _),
    integral_mul_const, integral_mul_const, integral_mul_const, integral_bary_mul_bary hD,
    integral_bary_mul_bary hD, integral_bary_mul_bary hD, sideLen_zero, sideLen_one, sideLen_two,
    dist2_sq, dist2_sq, dist2_sq]
  · simp only [Fin.isValue, Fin.reduceEq, ite_false]
    rw [show sqDist2 (p 2) (p 0) = sqDist2 (p 0) (p 2) by unfold sqDist2; ring]
    ring
  · exact (hi 0 1 _).add (hi 1 2 _)

open Intrinsic in
/-- **Chen–Holst Lemma 4.1 on one triangle** (`d = 2`, uniform density): for `u(x) = |x|²` and its
linear interpolant `u_I`, `∫_T |u_I − u| = (1/3) Σᵢ |pᵢ|² |T| − ∫_T |x|²`.

[formula: it-odt-interpolation-error] -/
theorem integral_abs_odt_error (hD : D p ≠ 0) :
    ∫ x in triangle p, |odtInterp p x - nsq x| =
      (nsq (p 0) + nsq (p 1) + nsq (p 2)) * area p / 3 - ∫ x in triangle p, nsq x := by
  rw [setIntegral_congr_fun (measurableSet_triangle p) (fun x hx =>
    abs_of_nonneg (sub_nonneg.2 (nsq_le_odtInterp hD hx))),
    integral_sub (integrableOn_triangle p (continuous_odtInterp p))
      (integrableOn_triangle p continuous_nsq), odtInterp, integral_interp hD]
  ring

open Intrinsic in
/-- **The star energy `odtEnergy` is the integral of the interpolant**: with
counter-clockwise triangles `T_k = (x, p_k, p_{k+1})`, `odtEnergy p n x = Σ_k ∫_{T_k} u_I`; with
`integral_abs_odt_error` this is Lemma 4.1 for the star, `Σ_k ∫_{T_k} |u_I − u| + Σ_k ∫_{T_k} |y|²`.

[formula: it-odt-interpolation-error] -/
theorem odtEnergy_eq_sum_integral (q : ℕ → ℝ × ℝ) (n : ℕ) (x : ℝ × ℝ)
    (hpos : ∀ k < n, 0 < sda x (q k) (q (k + 1))) :
    odtEnergy q n x = ∑ k ∈ Finset.range n, ∫ y in triangle ![x, q k, q (k + 1)],
      odtInterp ![x, q k, q (k + 1)] y := by
  unfold odtEnergy
  refine Finset.sum_congr rfl fun k hk => ?_
  have h := hpos k (Finset.mem_range.1 hk)
  have hD : D ![x, q k, q (k + 1)] = sda x (q k) (q (k + 1)) := rfl
  rw [odtInterp, integral_interp (by rw [hD]; exact h.ne'),
    show area ![x, q k, q (k + 1)] = |D ![x, q k, q (k + 1)]| / 2 by
      rw [abs_D_eq_two_mul_area]; ring, hD, abs_of_pos h]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons]
  ring

open Intrinsic in
/-- **Lemma 4.1 for a star** (as a sum over its triangles): `S(x) = Σ_k ∫_{T_k} |u_I − u| +
Σ_k ∫_{T_k} |y|²`, so on a fixed tiling of the star region `S` and the interpolation error differ by
a constant. -/
theorem odtEnergy_eq_sum_error_add (q : ℕ → ℝ × ℝ) (n : ℕ) (x : ℝ × ℝ)
    (hpos : ∀ k < n, 0 < sda x (q k) (q (k + 1))) :
    odtEnergy q n x = ∑ k ∈ Finset.range n, ((∫ y in triangle ![x, q k, q (k + 1)],
      |odtInterp ![x, q k, q (k + 1)] y - nsq y|) + ∫ y in triangle ![x, q k, q (k + 1)], nsq y) := by
  rw [odtEnergy_eq_sum_integral q n x hpos]
  refine Finset.sum_congr rfl fun k hk => ?_
  have h := hpos k (Finset.mem_range.1 hk)
  have hD : D ![x, q k, q (k + 1)] ≠ 0 := h.ne'
  rw [integral_abs_odt_error hD, odtInterp, integral_interp hD]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons]
  ring

end ChenHolst

/-! ## `it-odt-energy`: Chen–Holst Lemma 4.1 over a triangulation -/

section Odt

open Intrinsic

variable {ι : Type*} {s : Finset ι} {q : ι → Fin 3 → ℝ × ℝ} {Ω : Set (ℝ × ℝ)}

theorem integrableOn_odt_error (p : Fin 3 → ℝ × ℝ) {g : ℝ × ℝ → ℝ}
    (hg : ∀ x ∈ triangle p, g x = odtInterp p x) :
    IntegrableOn (fun x => |g x - nsq x|) (triangle p) :=
  (integrableOn_triangle p ((continuous_odtInterp p).sub continuous_nsq).abs).congr_fun
    (fun x hx => by simp only [hg x hx, Pi.sub_apply]) (measurableSet_triangle _)

/-- The interpolation error over a tiling is the sum of the errors of the tiles, for any `u_I`
that is the linear interpolant of `|x|²` on every tile. -/
theorem IsTiling.integral_abs_odt_error_eq_sum (h : IsTiling s q Ω) {g : ℝ × ℝ → ℝ}
    (hg : ∀ k ∈ s, ∀ x ∈ triangle (q k), g x = odtInterp (q k) x) :
    ∫ x in Ω, |g x - nsq x| = ∑ k ∈ s, ∫ x in triangle (q k), |odtInterp (q k) x - nsq x| := by
  rw [h.integral_eq_sum fun k hk => integrableOn_odt_error (q k) (hg k hk)]
  refine Finset.sum_congr rfl fun k hk =>
    setIntegral_congr_fun (measurableSet_triangle _) fun x hx => ?_
  simp only [hg k hk x hx]

/-- **Chen–Holst Lemma 4.1 for a triangulation** (`d = 2`, uniform density): over a tiling of `Ω`,
`E = ∫_Ω |u_I − u| = Σ_k (|T_k|/3) Σ_i |q_{k,i}|² − ∫_Ω |x|²` for `u = |x|²`.

[formula: it-odt-energy] -/
theorem IsTiling.integral_abs_odt_error (h : IsTiling s q Ω) {g : ℝ × ℝ → ℝ}
    (hg : ∀ k ∈ s, ∀ x ∈ triangle (q k), g x = odtInterp (q k) x) :
    ∫ x in Ω, |g x - nsq x| =
      (∑ k ∈ s, (nsq (q k 0) + nsq (q k 1) + nsq (q k 2)) * area (q k) / 3) - ∫ x in Ω, nsq x := by
  rw [h.integral_abs_odt_error_eq_sum hg, h.integral_eq_sum_of_continuous continuous_nsq,
    ← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun k hk =>
    _root_.DDG.TriangleIntegral.integral_abs_odt_error (h.nondeg k hk)

/-- **The interpolation error is a function of the edge lengths**:
`E = Σ_k (|T_k|/12)(ℓ₀₁² + ℓ₁₂² + ℓ₂₀²)` over a tiling.

[formula: it-odt-energy] -/
theorem IsTiling.integral_abs_odt_error_eq_sides (h : IsTiling s q Ω) {g : ℝ × ℝ → ℝ}
    (hg : ∀ k ∈ s, ∀ x ∈ triangle (q k), g x = odtInterp (q k) x) :
    ∫ x in Ω, |g x - nsq x| = ∑ k ∈ s, area (q k) / 12 *
      (sideLen (q k) 0 ^ 2 + sideLen (q k) 1 ^ 2 + sideLen (q k) 2 ^ 2) := by
  rw [h.integral_abs_odt_error_eq_sum hg]
  exact Finset.sum_congr rfl fun k hk => integral_odt_error_eq_sides (h.nondeg k hk)

/-- The area of the star of the vertex `i`: `Σ_k Σ_c [v_{k,c} = i] A_k` (each corner of a tile at
`i` counts the tile once). -/
def starAreaOf {V : Type*} [DecidableEq V] (s : Finset ι) (v : ι → Fin 3 → V) (A : ι → ℝ)
    (i : V) : ℝ :=
  ∑ k ∈ s, ∑ c : Fin 3, if v k c = i then A k else 0

/-- Regrouping a sum over the corners of the tiles by vertex. -/
theorem sum_corners_eq_sum_vertices {V : Type*} [Fintype V] [DecidableEq V] (s : Finset ι)
    (v : ι → Fin 3 → V) (A : ι → ℝ) (F : V → ℝ) :
    ∑ k ∈ s, (F (v k 0) + F (v k 1) + F (v k 2)) * A k = ∑ i, F i * starAreaOf s v A i := by
  simp only [starAreaOf, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [Finset.sum_comm]
  simp only [mul_ite, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true,
    Fin.sum_univ_three]
  ring

/-- **Lemma 4.1 as Chen and Holst print it**: for the tiles of a mesh with vertices `X` and
corner map `v`, `E = (1/3) Σ_i |x_i|² |ω_i| − ∫_Ω |x|²`, `|ω_i|` the area of the star of `x_i`.

[formula: it-odt-energy] -/
theorem IsTiling.integral_abs_odt_error_vertices {V : Type*} [Fintype V] [DecidableEq V]
    (h : IsTiling s q Ω) (v : ι → Fin 3 → V) (X : V → ℝ × ℝ)
    (hq : ∀ k ∈ s, ∀ c, q k c = X (v k c)) {g : ℝ × ℝ → ℝ}
    (hg : ∀ k ∈ s, ∀ x ∈ triangle (q k), g x = odtInterp (q k) x) :
    ∫ x in Ω, |g x - nsq x| =
      (∑ i, nsq (X i) * starAreaOf s v (fun k => area (q k)) i) / 3 - ∫ x in Ω, nsq x := by
  rw [h.integral_abs_odt_error hg, ← sum_corners_eq_sum_vertices, Finset.sum_div]
  congr 1
  refine Finset.sum_congr rfl fun k hk => ?_
  rw [hq k hk 0, hq k hk 1, hq k hk 2]

/-- The fan triangle `(x, q_k, q_{k+1})` of a star. -/
abbrev fan (q : ℕ → ℝ × ℝ) (x : ℝ × ℝ) (k : ℕ) : Fin 3 → ℝ × ℝ := ![x, q k, q (k + 1)]

theorem fan_term (q : ℕ → ℝ × ℝ) (x : ℝ × ℝ) (k : ℕ) (hpos : 0 < sda x (q k) (q (k + 1))) :
    (nsq (fan q x k 0) + nsq (fan q x k 1) + nsq (fan q x k 2)) * area (fan q x k) / 3 =
      (nsq x + (nsq (q k) + nsq (q (k + 1)))) * (sda x (q k) (q (k + 1)) / 2) / 3 := by
  have hD : D (fan q x k) = sda x (q k) (q (k + 1)) := rfl
  rw [show area (fan q x k) = |D (fan q x k)| / 2 by rw [abs_D_eq_two_mul_area]; ring, hD,
    abs_of_pos hpos]
  simp only [fan, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons]
  ring

/-- **The star energy is the interpolation error up to a constant.** A mesh whose tiles are the
counter-clockwise fan `(x, q_k, q_{k+1})`, `k < n`, around `x` and fixed triangles `r_j`, `j ∈ t`,
tiling `Ω`: `E(x) = S(x) + Σ_j (|T_j|/3) Σ_i |r_{j,i}|² − ∫_Ω |y|²`, with `S = odtEnergy`
(above), the only term that depends on `x`.

[formula: it-odt-energy] -/
theorem odt_error_star {κ : Type*} {n : ℕ} {t : Finset κ} {r : κ → Fin 3 → ℝ × ℝ}
    {q : ℕ → ℝ × ℝ} {x : ℝ × ℝ} (hpos : ∀ k < n, 0 < sda x (q k) (q (k + 1)))
    (h : IsTiling ((Finset.range n).disjSum t) (Sum.elim (fan q x) r) Ω) {g : ℝ × ℝ → ℝ}
    (hg : ∀ k ∈ (Finset.range n).disjSum t, ∀ y ∈ triangle (Sum.elim (fan q x) r k),
      g y = odtInterp (Sum.elim (fan q x) r k) y) :
    ∫ y in Ω, |g y - nsq y| = odtEnergy q n x +
      (∑ j ∈ t, (nsq (r j 0) + nsq (r j 1) + nsq (r j 2)) * area (r j) / 3) - ∫ y in Ω, nsq y := by
  rw [h.integral_abs_odt_error hg, Finset.sum_disjSum]
  simp only [Sum.elim_inl, Sum.elim_inr]
  congr 2
  unfold odtEnergy
  exact Finset.sum_congr rfl fun k hk => fan_term q x k (hpos k (Finset.mem_range.1 hk))

/-- **Moving the centre changes the interpolation error by the change of the star energy**, for
two positions whose fans tile the same region with the fixed triangles.

[formula: it-odt-energy] -/
theorem odt_error_sub {κ : Type*} {n : ℕ} {t : Finset κ} {r : κ → Fin 3 → ℝ × ℝ}
    {q : ℕ → ℝ × ℝ} {x x' : ℝ × ℝ} (hpos : ∀ k < n, 0 < sda x (q k) (q (k + 1)))
    (hpos' : ∀ k < n, 0 < sda x' (q k) (q (k + 1)))
    (h : IsTiling ((Finset.range n).disjSum t) (Sum.elim (fan q x) r) Ω)
    (h' : IsTiling ((Finset.range n).disjSum t) (Sum.elim (fan q x') r) Ω) {g g' : ℝ × ℝ → ℝ}
    (hg : ∀ k ∈ (Finset.range n).disjSum t, ∀ y ∈ triangle (Sum.elim (fan q x) r k),
      g y = odtInterp (Sum.elim (fan q x) r k) y)
    (hg' : ∀ k ∈ (Finset.range n).disjSum t, ∀ y ∈ triangle (Sum.elim (fan q x') r k),
      g' y = odtInterp (Sum.elim (fan q x') r k) y) :
    (∫ y in Ω, |g' y - nsq y|) - ∫ y in Ω, |g y - nsq y| = odtEnergy q n x' - odtEnergy q n x := by
  rw [odt_error_star hpos' h' hg', odt_error_star hpos h hg]
  ring

/-- **The update of Eq. 4.13 lowers the interpolation error.** With a closed counter-clockwise
link, the step `x' = x + odtDisplacement(w, |ω|, α)` for `0 ≤ α ≤ 2`, if the moved fan still tiles
`Ω` with the fixed triangles, does not increase `E = ∫_Ω |u_I − u|`; for `α = 1`, `E(x') ≤ E(y)`
for every position `y` whose fan tiles `Ω` with them.

[formula: it-odt-energy] -/
theorem odt_error_step_le {κ : Type*} {n : ℕ} {t : Finset κ} {r : κ → Fin 3 → ℝ × ℝ}
    {q : ℕ → ℝ × ℝ} (hq : q n = q 0) {x x' : ℝ × ℝ} (hω : 0 < starArea q n x)
    {step : ℝ} (h0 : 0 ≤ step) (h2 : step ≤ 2)
    (hx' : x' = x + odtDisplacement
      (∑ k ∈ Finset.range n, (circumVec x (q k) (q (k + 1)) - (sda x (q k) (q (k + 1)) / 2) • x))
      (starArea q n x) step)
    (hpos : ∀ k < n, 0 < sda x (q k) (q (k + 1)))
    (hpos' : ∀ k < n, 0 < sda x' (q k) (q (k + 1)))
    (h : IsTiling ((Finset.range n).disjSum t) (Sum.elim (fan q x) r) Ω)
    (h' : IsTiling ((Finset.range n).disjSum t) (Sum.elim (fan q x') r) Ω) {g g' : ℝ × ℝ → ℝ}
    (hg : ∀ k ∈ (Finset.range n).disjSum t, ∀ y ∈ triangle (Sum.elim (fan q x) r k),
      g y = odtInterp (Sum.elim (fan q x) r k) y)
    (hg' : ∀ k ∈ (Finset.range n).disjSum t, ∀ y ∈ triangle (Sum.elim (fan q x') r k),
      g' y = odtInterp (Sum.elim (fan q x') r k) y) :
    (∫ y in Ω, |g' y - nsq y|) ≤ ∫ y in Ω, |g y - nsq y| ∧
      (step = 1 → ∀ (y : ℝ × ℝ) (gy : ℝ × ℝ → ℝ), (∀ k < n, 0 < sda y (q k) (q (k + 1))) →
        IsTiling ((Finset.range n).disjSum t) (Sum.elim (fan q y) r) Ω →
        (∀ k ∈ (Finset.range n).disjSum t, ∀ z ∈ triangle (Sum.elim (fan q y) r k),
          gy z = odtInterp (Sum.elim (fan q y) r k) z) →
        (∫ z in Ω, |g' z - nsq z|) ≤ ∫ z in Ω, |gy z - nsq z|) := by
  have hS := odtEnergy_step_le hq x hω h0 h2
  simp only [← hx'] at hS
  have e := odt_error_sub hpos hpos' h h' hg hg'
  refine ⟨by linarith [hS.1], fun h1 y gy hposy hy hgy => ?_⟩
  have hmin := (odtEnergy_step_le hq x hω zero_le_one one_le_two).2 y
  simp only [← h1, ← hx'] at hmin
  have ey := odt_error_sub hposy hpos' hy h' hgy hg'
  linarith

end Odt

end TriangleIntegral

end DDG
