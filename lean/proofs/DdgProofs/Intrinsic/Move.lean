import DdgProofs.Intrinsic.Diamond

/-!
# Intrinsic triangulations: moving a vertex, the optimal Delaunay triangulation

`crates/ddg-intrinsic/SPEC.md` §4 (Move a vertex, Optimal Delaunay triangulation), N19, N20, N21;
the paper Sec. 3.3.3, Alg. 10, Sec. 4.3; L. Chen, M. Holst, *Efficient mesh optimization schemes
based on Optimal Delaunay Triangulations*, CMAME 200 (2011), Lemma 4.2 (Eq. 4.3–4.6), Eq. 4.13,
Sec. 4.2.1 of the preprint (barycentres at the boundary).

Plane points are pairs `ℝ × ℝ` with `sda` (doubled signed area) and `sqDist2`, `dist2` of
`DdgProofs/Bff/Quality.lean`.

* **`it-move-length`** — the oracle's `movedLength l r (cos α) (sin α)` is the paper's law of
  cosines `√(r² + ℓ² − 2 r ℓ cos α)` (`movedLength_eq_law_of_cosines`) and the distance between
  the moved end `r (cos φ, sin φ)` and the far end `ℓ (cos ψ, sin ψ)` when `α = ψ − φ`
  (`movedLength_eq_dist2`).
* **`it-move-area`** — `movedDoubleArea (sda 0 j k) v (k − j)` is the doubled signed area of
  the triangle `(v, j, k)` (`movedDoubleArea_eq_sda`); around a closed link the doubled signed
  areas of the star add up to a number that does not depend on the centre
  (`sum_sda_closed_link`): moving the vertex keeps the area of its star.
* **`it-circumcenter-vector`** — `weightedCircumcenter ℓ_ab ℓ_ac ℓ_bc A / A` is equidistant
  from `a = 0`, `b = (ℓ_ab, 0)` and `c = (c_x, h)` (`weightedCircumcenter_equidistant`), also
  stated with `A` the Heron area of the three lengths and `c` the oracle's apex
  (`weightedCircumcenter_equidistant_heron`).
* **`it-odt-update`** — `weightedBarycenter` is `A (a + b + c)/3` in the same frame
  (`weightedBarycenter_eq`). For a vertex `x` with a closed counter-clockwise link
  `p 0, …, p (n−1), p n = p 0` and `|T_k| = sda x (p k) (p (k+1)) / 2`, the energy of Chen and
  Holst (Lemma 4.1 with `d = 2`, the star's part) `S(x) = Σ_k (|x|² + |p_k|² + |p_{k+1}|²) |T_k| / 3`
  has the derivative `(2/3) Σ_k (|T_k| x − |T_k| c_k)` (Eq. 4.5; `hasFDerivAt_odtEnergy`), with
  `|T_k| c_k` the polynomial `circumVec`, which divided by `|T_k|` is equidistant from the three
  corners (`circumVec_equidistant`) and in the frame of a triangle is the oracle's
  `weightedCircumcenter` (`circumVec_layout`). `S` is exactly quadratic, with Hessian
  `(2/3) |ω| I` (Lemma 4.3; `odtEnergy_add`); hence the update of Eq. 4.13, the oracle's
  `odtDisplacement` of `w = Σ_k |T_k| (c_k − x) = −(3/2) ∇S` by the area `|ω|` of the star,
  changes the energy by exactly `−(3 |∇S|²/(4|ω|)) α (2 − α)`: a descent step for every
  `0 < α < 2`, the minimum of `S` for `α = 1` (`odtEnergy_step`, `odtEnergy_step_le`).

Not covered here: that `S` is the interpolation error of `|x|²` up to a constant (Lemma 4.1;
proved over a tiling in `DdgProofs/TriangleIntegral/Tiling.lean`, `it-odt-energy`), what the
minimizing triangulations are and the boundary rule as a remedy (`CITED.md`,
`it-odt-delaunay-optimality`). The
identities are about a planar star with signed areas; that the star of an intrinsic vertex
unfolds into such a star is SPEC N19.

Own definitions: `rot`, `nsq`, `circumVec`, `odtEnergy`, `odtGradient`, `toDual`.
-/

noncomputable section

namespace DdgProofs

namespace Intrinsic

open scoped DdgProofs.Oracle
open Real Finset

/-! ## Plane helpers -/

/-- The quarter turn `rot u = (−u₂, u₁)`. -/
def rot (u : ℝ × ℝ) : ℝ × ℝ := (-u.2, u.1)

/-- The squared norm `|u|² = u₁² + u₂²`. -/
def nsq (u : ℝ × ℝ) : ℝ := u.1 ^ 2 + u.2 ^ 2

/-- The linear form `v ↦ G₁ v₁ + G₂ v₂` (the gradient `G` as a derivative). -/
def toDual (G : ℝ × ℝ) : (ℝ × ℝ) →L[ℝ] ℝ :=
  G.1 • ContinuousLinearMap.fst ℝ ℝ ℝ + G.2 • ContinuousLinearMap.snd ℝ ℝ ℝ

@[simp] theorem toDual_apply (G v : ℝ × ℝ) : toDual G v = G.1 * v.1 + G.2 * v.2 := by
  simp [toDual]

/-! ## `it-move-length` -/

theorem movedLength_eq (l r c s : ℝ) :
    Ddg.movedLength l r c s = √((l - r * c) ^ 2 + (r * s) ^ 2) := by
  unfold Ddg.movedLength
  simp only [Oracle.sqrt_eq]
  congr 1
  ring

/-- **The length of a moved edge is the paper's law of cosines** (Alg. 10):
`movedLength ℓ r (cos α) (sin α) = √(r² + ℓ² − 2 r ℓ cos α)`.

[formula: it-move-length] -/
theorem movedLength_eq_law_of_cosines (l r α : ℝ) :
    Ddg.movedLength l r (cos α) (sin α) = √(r ^ 2 + l ^ 2 - 2 * r * l * cos α) := by
  rw [movedLength_eq]
  congr 1
  linear_combination r ^ 2 * sin_sq_add_cos_sq α

/-- **… and the distance from the moved end to the far end**: the vertex moves from the origin
to `r (cos φ, sin φ)`, the far end of the edge is at `ℓ (cos ψ, sin ψ)`; with `α = ψ − φ` (the
angle from the motion to the edge) the new length is their distance.

[formula: it-move-length] -/
theorem movedLength_eq_dist2 (l r φ ψ : ℝ) :
    Ddg.movedLength l r (cos (ψ - φ)) (sin (ψ - φ)) =
      dist2 (r * cos φ, r * sin φ) (l * cos ψ, l * sin ψ) := by
  rw [movedLength_eq]
  unfold dist2
  congr 1
  dsimp only
  rw [cos_sub, sin_sub]
  linear_combination (r ^ 2 * (sin φ ^ 2 + cos φ ^ 2) - l ^ 2) * sin_sq_add_cos_sq ψ

/-- The length does not change for a zero motion. -/
theorem movedLength_zero {l : ℝ} (hl : 0 ≤ l) (c s : ℝ) : Ddg.movedLength l 0 c s = l := by
  rw [movedLength_eq]
  simp [sqrt_sq hl]

/-! ## `it-move-area` -/

theorem movedDoubleArea_eq (A2 : ℝ) (v d : ℝ × ℝ) :
    Ddg.movedDoubleArea A2 v d = A2 - (v.1 * d.2 - v.2 * d.1) := rfl

/-- **The doubled signed area after the move**: for the triangle `(0, j, k)` of the star, the
vertex at the origin moved to `v`, `sda v j k = movedDoubleArea (sda 0 j k) v (k − j)`.

[formula: it-move-area] -/
theorem movedDoubleArea_eq_sda (v j k : ℝ × ℝ) :
    sda v j k = Ddg.movedDoubleArea (sda (0, 0) j k) v (k - j) := by
  rw [movedDoubleArea_eq]
  unfold sda
  simp only [Prod.fst_sub, Prod.snd_sub]
  ring

/-- The doubled signed area relative to another centre: `sda v a b = sda w a b + (v − w) × (a − b)`. -/
theorem sda_change_centre (v w a b : ℝ × ℝ) :
    sda v a b = sda w a b + ((v.1 - w.1) * (a.2 - b.2) - (v.2 - w.2) * (a.1 - b.1)) := by
  unfold sda
  ring

/-- **The star keeps its area**: around a closed link (`p n = p 0`) the sum of the doubled signed
areas of the triangles `(v, p k, p (k+1))` does not depend on `v` (the differences telescope).

[formula: it-move-area] -/
theorem sum_sda_closed_link {p : ℕ → ℝ × ℝ} {n : ℕ} (hp : p n = p 0) (v w : ℝ × ℝ) :
    ∑ k ∈ range n, sda v (p k) (p (k + 1)) = ∑ k ∈ range n, sda w (p k) (p (k + 1)) := by
  rw [← sub_eq_zero, ← sum_sub_distrib]
  have h : ∀ k ∈ range n, sda v (p k) (p (k + 1)) - sda w (p k) (p (k + 1)) =
      (v.1 - w.1) * ((p k).2 - (p (k + 1)).2) - (v.2 - w.2) * ((p k).1 - (p (k + 1)).1) := by
    intro k _
    rw [sda_change_centre v w]
    ring
  rw [sum_congr rfl h, sum_sub_distrib, ← mul_sum, ← mul_sum,
    sum_range_sub' (fun k => (p k).2), sum_range_sub' (fun k => (p k).1), hp]
  ring

/-! ## `it-circumcenter-vector` -/

theorem weightedCircumcenter_eq (lab lac lbc A : ℝ) :
    Ddg.weightedCircumcenter lab lac lbc A =
      (A * lab / 2, lab * (lac * lac + lbc * lbc - lab * lab) / 8) := by
  have e8 : (@OfNat.ofNat ℝ 8 Ddg.instOfNatScalar) = (8 : ℝ) := by rw [Oracle.ofNat_eq]
  simp only [Ddg.weightedCircumcenter, Oracle.two_eq, e8]

/-- **The weighted circumcentre divided by the area is the circumcentre**: in the frame of `a`
with `b = (ℓ_ab, 0)` and `c = (c_x, h)`, `c_x = (ℓ_ab² + ℓ_ac² − ℓ_bc²)/2ℓ_ab`, `h = 2A/ℓ_ab`,
`A ≠ 0`, `ℓ_ab ≠ 0` and `|c| = ℓ_ac`, the point `o = weightedCircumcenter … / A` has the same
squared distance from `a = 0`, `b` and `c`.

[formula: it-circumcenter-vector] -/
theorem weightedCircumcenter_equidistant {lab lac lbc A cx h : ℝ} (hA : A ≠ 0) (hab : lab ≠ 0)
    (hcx : cx = (lab ^ 2 + lac ^ 2 - lbc ^ 2) / (2 * lab)) (hh : h = 2 * A / lab)
    (hc : cx ^ 2 + h ^ 2 = lac ^ 2) :
    let o : ℝ × ℝ := ((Ddg.weightedCircumcenter lab lac lbc A).1 / A,
      (Ddg.weightedCircumcenter lab lac lbc A).2 / A)
    sqDist2 (0, 0) o = sqDist2 (lab, 0) o ∧ sqDist2 (0, 0) o = sqDist2 (cx, h) o := by
  intro o
  have o1 : o.1 = lab / 2 := by
    simp only [o, weightedCircumcenter_eq]
    field_simp
  have o2 : o.2 * A = lab * (lac * lac + lbc * lbc - lab * lab) / 8 := by
    simp only [o, weightedCircumcenter_eq]
    field_simp
  have hcx' : cx * (2 * lab) = lab ^ 2 + lac ^ 2 - lbc ^ 2 := by
    rw [hcx]; field_simp
  have hh' : h * lab = 2 * A := by rw [hh]; field_simp
  unfold sqDist2
  dsimp only
  constructor
  · rw [o1]; ring
  · -- `|o − c|² − |o|² = |c|² − 2 o·c`, and `o·c = (ℓ_ab c_x + 2 o₂ h)/2 = ℓ_ac²/2`
    have key : (o.2 * A) * (h * lab) = (lab * (lac * lac + lbc * lbc - lab * lab) / 8) * (2 * A) := by
      rw [o2, hh']
    have key2 : o.2 * h = (lac * lac + lbc * lbc - lab * lab) / 4 := by
      have : (o.2 * h) * (A * lab) = ((lac * lac + lbc * lbc - lab * lab) / 4) * (A * lab) := by
        linear_combination key
      exact mul_right_cancel₀ (mul_ne_zero hA hab) this
    rw [o1]
    linear_combination (-1 : ℝ) * hc + hcx' / 2 + 2 * key2

/-- **The same with the Heron area**: for a triangle with the sides `ℓ_ab, ℓ_ac, ℓ_bc` (strict
triangle inequalities) and `A = triArea`, `c` the oracle's apex
`(diamondX ℓ_ab ℓ_ac ℓ_bc, diamondH ℓ_ab ℓ_ac ℓ_bc)`, the point
`weightedCircumcenter ℓ_ab ℓ_ac ℓ_bc A / A` is equidistant from `0`, `(ℓ_ab, 0)` and `c`.

[formula: it-circumcenter-vector] -/
theorem weightedCircumcenter_equidistant_heron {lab lac lbc : ℝ} (ht : IsTriangle lab lac lbc) :
    let A := triArea lab lac lbc
    let o : ℝ × ℝ := ((Ddg.weightedCircumcenter lab lac lbc A).1 / A,
      (Ddg.weightedCircumcenter lab lac lbc A).2 / A)
    sqDist2 (0, 0) o = sqDist2 (lab, 0) o ∧
      sqDist2 (0, 0) o = sqDist2 (Ddg.diamondX lab lac lbc, Ddg.diamondH lab lac lbc) o := by
  intro A o
  have hab : 0 < lab := by linarith [ht.lt_b, ht.lt_c]
  have hA : 0 < A := triArea_pos ht
  obtain ⟨-, -, -, hb, -, -, -⟩ := diamond_isometric ht ht
  have hc : Ddg.diamondX lab lac lbc ^ 2 + Ddg.diamondH lab lac lbc ^ 2 = lac ^ 2 := by
    have h0 := dist2_sq' (0, 0) (apexK lab lac lbc)
    rw [hb] at h0
    rw [h0]
    simp [apexK]
  exact weightedCircumcenter_equidistant hA.ne' hab.ne'
    (by rw [diamondX_eq]; ring) (by rw [diamondH_eq]) hc

/-! ## `it-odt-update`: the barycentre and the displacement -/

/-- **The weighted barycentre** is `A (a + b + c)/3` with `a = 0`, `b = (ℓ_ab, 0)` and
`c = (c_x, 2A/ℓ_ab)`, `c_x` the oracle's `diamondX` (the frame of the circumcentre).

[formula: it-odt-update] -/
theorem weightedBarycenter_eq (lab lac lbc A : ℝ) :
    Ddg.weightedBarycenter lab lac lbc A =
      A • ((1 / 3 : ℝ) • ((0, 0) + (lab, 0) + (Ddg.diamondX lab lac lbc, 2 * A / lab))) := by
  simp only [Ddg.weightedBarycenter, Oracle.two_eq, Oracle.three_eq, Prod.smul_mk, Prod.mk_add_mk,
    smul_eq_mul]
  ext <;> dsimp only <;> ring

theorem odtDisplacement_eq (w : ℝ × ℝ) (total step : ℝ) :
    Ddg.odtDisplacement w total step = (step / total) • w := by
  simp only [Ddg.odtDisplacement, Prod.smul_def, smul_eq_mul]
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

/-- **In the frame of a triangle it is the oracle's `weightedCircumcenter`**: with `x = 0`,
`b = (ℓ_ab, 0)`, `c = (c_x, h)`, `c_x = (ℓ_ab² + ℓ_ac² − ℓ_bc²)/2ℓ_ab`, `|c| = ℓ_ac` and
`A = ℓ_ab h / 2`.

[formula: it-odt-update] -/
theorem circumVec_layout {lab lac lbc cx h : ℝ} (hab : lab ≠ 0)
    (hcx : cx = (lab ^ 2 + lac ^ 2 - lbc ^ 2) / (2 * lab)) (hc : cx ^ 2 + h ^ 2 = lac ^ 2) :
    circumVec (0, 0) (lab, 0) (cx, h) = Ddg.weightedCircumcenter lab lac lbc (lab * h / 2) := by
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
    Ddg.odtDisplacement w (starArea p n x) step =
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
    odtEnergy p n (x + Ddg.odtDisplacement w (starArea p n x) step) =
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
    odtEnergy p n (x + Ddg.odtDisplacement w (starArea p n x) step) ≤ odtEnergy p n x ∧
      ∀ y, odtEnergy p n (x + Ddg.odtDisplacement w (starArea p n x) 1) ≤ odtEnergy p n y := by
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

end DdgProofs
