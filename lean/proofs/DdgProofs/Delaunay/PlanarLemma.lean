import DdgProofs.Delaunay.PlanarTermination

/-!
# Planar Delaunay triangulations, part 4: Delaunay's lemma

**Delaunay's lemma.** A triangulation of a finite point set `S` that is locally Delaunay at every
interior edge is Delaunay: no point of `S` lies strictly inside the circumcircle of any triangle
(`delaunay_lemma`). With the termination of Lawson's algorithm
(`DdgProofs/Delaunay/PlanarTermination.lean`): from every triangulation of `S` a sequence of
Delaunay flips reaches a Delaunay triangulation of `S` (`exists_flips_to_delaunay`); in
particular a Delaunay triangulation of `S` exists as soon as some triangulation of `S` does
(`exists_delaunay_of_triangulation`; the existence of a first triangulation is a hypothesis).

## The proof: walking towards the point

Let `t₀ ∈ T` have `p ∈ S` strictly inside its circumcircle. Pick `x` in the interior of `t₀` such
that the line `xp` contains no point of `S` other than `p` (`exists_generic_interior`: finitely
many lines through `p` have empty interior). Walk along `γ(τ) = x + τ (p − x)`, `τ ∈ [0, 1]`.
The function `circleLift t p = inCircle(t, p) / sda(t)` (the negative power of `p` with respect
to the circumcircle of `t`; positive inside) is taken as a level `φ₀ = circleLift t₀ p > 0`, and
`τ*` is the last parameter at which `γ` is in a triangle `t` with `circleLift t p ≥ φ₀`
(a maximum, the set being closed). Then:

* `τ* < 1`, since `p` is in no such triangle (it would be a corner, and corners have power `0`);
* `γ` leaves `t` at `τ*` through the relative interior of an edge `ab` (`exit_edge`; not through
  a corner, which would be a point of `S` on the line `xp`), and `p` is strictly beyond `ab`;
* just after `τ*`, `γ` is in a triangle `t' ≠ t` that contains `γ(τ*)` (cover, finitely many closed
  triangles); a triangle containing a relative interior point of the edge `ab` of `t` is `t` or
  the triangle on the other side of `ab` (**adjacency**, `eq_insert_of_mem_edge`: from disjoint
  interiors and the vertex condition), so `t' = {a, b, d}` with `d` across `ab`;
* the edge `ab` is locally Delaunay, so `circleLift t' p ≥ circleLift t p ≥ φ₀`
  (`circleLift_le_of_not_inCircumcircle`: across a Delaunay edge the lifted plane goes up) —
  contradicting the maximality of `τ*`.

This is the walk argument of Edelsbrunner's proof (power of `p` monotone along the walk), with
the genericity of the walk made explicit and the "next triangle" given by the adjacency lemma.

## References

* B. Delaunay, *Sur la sphère vide*, Bull. Acad. Sci. URSS (1934), 793–800 (the lemma, in any
  dimension).
* H. Edelsbrunner, *Geometry and Topology for Mesh Generation*, Cambridge University Press, 2001,
  Sec. 1.2 (Delaunay lemma; proof by walking towards the point).
* M. de Berg, O. Cheong, M. van Kreveld, M. Overmars, *Computational Geometry: Algorithms and
  Applications*, 3rd ed., Springer, 2008, Ch. 9, Theorem 9.8 (a triangulation is legal iff it is a
  Delaunay triangulation; proof by angle maximization).
* R. Sibson, *Locally equiangular triangulations*, The Computer Journal 21(3), 243–245, 1978.

The bibliographic data were not checked against the works here.
-/

noncomputable section

namespace DdgProofs

namespace Delaunay

namespace Planar

open MeasureTheory Set Filter Topology

/-- If `x` lies on the lines `ab` and `ca` of a nondegenerate triangle, it is `a`. -/
theorem eq_of_sda_eq_zero' {a b c x : ℝ × ℝ} (h : sda a b c ≠ 0) (h2 : sda a x c = 0)
    (h3 : sda a b x = 0) : x = a := by
  have h1 : sda x b c = sda a b c := by have := sda_four a b c x; linarith
  rw [sda_cramer h x, h1, h2, h3, div_self h]
  simp

/-! ## Points on a line, compared along it -/

/-- Two points `P ∈ K`, `Q ∈ K'` on the same side of `y ∈ K ∩ K'` on a line through `y` (`N P = N y
+ t_P n`, `N > 0`): one lies between `y` and the other, hence in both convex sets; if that forces
equality either way, `P = Q`. -/
theorem eq_of_same_side {K K' : Set (ℝ × ℝ)} (hK : Convex ℝ K) (hK' : Convex ℝ K')
    {y P Q n : ℝ × ℝ} {N tP tQ : ℝ} (hN : 0 < N) (hyK : y ∈ K) (hyK' : y ∈ K') (hP : P ∈ K)
    (hQ : Q ∈ K') (hP1 : N * P.1 = N * y.1 + tP * n.1) (hP2 : N * P.2 = N * y.2 + tP * n.2)
    (hQ1 : N * Q.1 = N * y.1 + tQ * n.1) (hQ2 : N * Q.2 = N * y.2 + tQ * n.2) (hs : 0 < tP * tQ)
    (h1 : Q ∈ K → Q = P) (h2 : P ∈ K' → P = Q) : P = Q := by
  have htP : tP ≠ 0 := by rintro rfl; simp at hs
  have htQ : tQ ≠ 0 := by rintro rfl; simp at hs
  have hθ : 0 < tQ / tP := by
    have : tQ / tP = (tP * tQ) / (tP * tP) := by field_simp
    rw [this]; exact div_pos hs (mul_self_pos.2 htP)
  rcases le_or_gt (tQ / tP) 1 with hle | hgt
  · refine (h1 ?_).symm
    have e : Q = (1 - tQ / tP) • y + (tQ / tP) • P := by
      ext
      · simp only [Prod.fst_add, Prod.smul_fst, smul_eq_mul]
        have eA : Q.1 = y.1 + tQ * n.1 / N := by field_simp; linarith
        have eB : P.1 = y.1 + tP * n.1 / N := by field_simp; linarith
        rw [eA, eB]; field_simp; ring
      · simp only [Prod.snd_add, Prod.smul_snd, smul_eq_mul]
        have eA : Q.2 = y.2 + tQ * n.2 / N := by field_simp; linarith
        have eB : P.2 = y.2 + tP * n.2 / N := by field_simp; linarith
        rw [eA, eB]; field_simp; ring
    rw [e]
    exact hK hyK hP (by linarith) hθ.le (by ring)
  · refine h2 ?_
    have hθ' : 0 < tP / tQ := by
      have : tP / tQ = (tP * tQ) / (tQ * tQ) := by field_simp
      rw [this]; exact div_pos hs (mul_self_pos.2 htQ)
    have hle' : tP / tQ ≤ 1 := by
      have : tP / tQ = (tQ / tP)⁻¹ := (inv_div _ _).symm
      rw [this]; exact (inv_lt_one_of_one_lt₀ hgt).le
    have e : P = (1 - tP / tQ) • y + (tP / tQ) • Q := by
      ext
      · simp only [Prod.fst_add, Prod.smul_fst, smul_eq_mul]
        have eA : P.1 = y.1 + tP * n.1 / N := by field_simp; linarith
        have eB : Q.1 = y.1 + tQ * n.1 / N := by field_simp; linarith
        rw [eA, eB]; field_simp; ring
      · simp only [Prod.snd_add, Prod.smul_snd, smul_eq_mul]
        have eA : P.2 = y.2 + tP * n.2 / N := by field_simp; linarith
        have eB : Q.2 = y.2 + tQ * n.2 / N := by field_simp; linarith
        rw [eA, eB]; field_simp; ring
    rw [e]
    exact hK' hyK' hQ (by linarith) hθ'.le (by ring)

/-- `|b − a|²`. -/
def lineN (a b : ℝ × ℝ) : ℝ := (b.1 - a.1) ^ 2 + (b.2 - a.2) ^ 2

/-- The coordinate of `z` along the line through `y` with direction `b − a` (times `|b − a|²`). -/
def lineT (a b y z : ℝ × ℝ) : ℝ := (z.1 - y.1) * (b.1 - a.1) + (z.2 - y.2) * (b.2 - a.2)

theorem line_param {a b y z : ℝ × ℝ}
    (hl : (b.1 - a.1) * (z.2 - y.2) - (b.2 - a.2) * (z.1 - y.1) = 0) :
    lineN a b * z.1 = lineN a b * y.1 + lineT a b y z * (b - a).1 ∧
      lineN a b * z.2 = lineN a b * y.2 + lineT a b y z * (b - a).2 := by
  simp only [lineN, lineT, Prod.fst_sub, Prod.snd_sub]
  constructor
  · linear_combination (-(b.2 - a.2)) * hl
  · linear_combination (b.1 - a.1) * hl

theorem lineN_pos {a b : ℝ × ℝ} (h : a ≠ b) : 0 < lineN a b := by
  unfold lineN
  by_contra hle
  push Not at hle
  have h1 : b.1 - a.1 = 0 := by nlinarith [sq_nonneg (b.1 - a.1), sq_nonneg (b.2 - a.2)]
  have h2 : b.2 - a.2 = 0 := by nlinarith [sq_nonneg (b.1 - a.1), sq_nonneg (b.2 - a.2)]
  exact h (Prod.ext (by linarith) (by linarith))

theorem neg_of_mul_sq {P d N Q : ℝ} (h : P * d ^ 2 = N * Q) (hN : 0 < N) (hQ : Q < 0) : P < 0 := by
  by_contra hP; push Not at hP; nlinarith [mul_nonneg hP (sq_nonneg d), mul_neg_of_pos_of_neg hN hQ]

theorem neg_of_pos_mul {K t M : ℝ} (hK : 0 < K) (h : K * t = M) (hM : M < 0) : t < 0 := by
  by_contra h'; push Not at h'; have := mul_nonneg hK.le h'; linarith

theorem pos_of_pos_mul {K t M : ℝ} (hK : 0 < K) (h : K * t = M) (hM : 0 < M) : 0 < t := by
  by_contra h'; push Not at h'; have := mul_nonpos_of_nonneg_of_nonpos hK.le h'; linarith

/-! ## Adjacency across an edge -/

set_option maxHeartbeats 1000000 in
/-- The ordered core of `eq_insert_of_mem_edge`: `y` is in the relative interior of the edge `ab` of
`{a,b,c}` and on the edge `v₁v₂` of another triangle `{v₁,v₂,v₃}` of `T`. -/
theorem eq_insert_of_mem_edge_aux {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {a b c y v1 v2 v3 : ℝ × ℝ} (hu : {a, b, c} ∈ T)
    (hy1 : 0 < sda a b c * sda y b c) (hy2 : 0 < sda a b c * sda a y c) (hy3 : sda a b y = 0)
    (hyS : y ∉ S) (hu' : {v1, v2, v3} ∈ T) (hne : ({v1, v2, v3} : Finset (ℝ × ℝ)) ≠ {a, b, c})
    (hg1 : 0 ≤ sda v1 v2 v3 * sda y v2 v3) (hg2 : 0 ≤ sda v1 v2 v3 * sda v1 y v3)
    (hg3 : sda v1 v2 y = 0) :
    ∃ d, ({v1, v2, v3} : Finset (ℝ × ℝ)) = {a, b, d} ∧ sda a b c * sda a b d < 0 := by
  have hD : sda a b c ≠ 0 := hT.sda_ne_zero hu
  have hD' : sda v1 v2 v3 ≠ 0 := hT.sda_ne_zero hu'
  have hab := sda_ne_zero_left hD
  have hv12 := sda_ne_zero_left hD'
  have hv1S : v1 ∈ S := hT.mem_S hu' (by simp)
  have hv2S : v2 ∈ S := hT.mem_S hu' (by simp)
  have haS : a ∈ S := hT.mem_S hu (by simp)
  have hbS : b ∈ S := hT.mem_S hu (by simp)
  have hyv1 : y ≠ v1 := fun h => hyS (h ▸ hv1S)
  have hyv2 : y ≠ v2 := fun h => hyS (h ▸ hv2S)
  have hg1' : 0 < sda v1 v2 v3 * sda y v2 v3 := lt_of_le_of_ne hg1 fun h0 =>
    hyv2 (eq_of_sda_eq_zero hD' ((mul_eq_zero.1 h0.symm).resolve_left hD') hg3)
  have hg2' : 0 < sda v1 v2 v3 * sda v1 y v3 := lt_of_le_of_ne hg2 fun h0 =>
    hyv1 (eq_of_sda_eq_zero' hD' ((mul_eq_zero.1 h0.symm).resolve_left hD') hg3)
  have hyu : y ∈ convexHull ℝ {a, b, c} := (mem_hull3_iff hD).2 ⟨hy1.le, hy2.le, by rw [hy3, mul_zero]⟩
  have hyu' : y ∈ convexHull ℝ {v1, v2, v3} :=
    (mem_hull3_iff hD').2 ⟨hg1, hg2, by rw [hg3, mul_zero]⟩
  have hdisj : Disjoint (interior (convexHull ℝ {v1, v2, v3})) (interior (convexHull ℝ {a, b, c})) := by
    have : Disjoint (interior (hull {v1, v2, v3})) (interior (hull {a, b, c})) :=
      hT.disjoint hu' hu hne
    rwa [hull_triple, hull_triple] at this
  -- no direction points into both triangles
  have nocommon : ∀ w : ℝ × ℝ, 0 < sda a b c * ((b.1 - a.1) * w.2 - (b.2 - a.2) * w.1) →
      0 < sda v1 v2 v3 * ((v2.1 - v1.1) * w.2 - (v2.2 - v1.2) * w.1) → False := by
    intro w hw1 hw2
    have hO : IsOpen {z : ℝ × ℝ | 0 < sda a b c * sda z b c ∧ 0 < sda a b c * sda a z c ∧
        0 < sda v1 v2 v3 * sda z v2 v3 ∧ 0 < sda v1 v2 v3 * sda v1 z v3} :=
      (isOpen_lt continuous_const (continuous_const.mul (continuous_sda₁ b c))).inter
        ((isOpen_lt continuous_const (continuous_const.mul (continuous_sda₂ a c))).inter
        ((isOpen_lt continuous_const (continuous_const.mul (continuous_sda₁ v2 v3))).inter
        (isOpen_lt continuous_const (continuous_const.mul (continuous_sda₂ v1 v3)))))
    obtain ⟨s, hs, hspos⟩ := ((exists_pos_near (hO.mem_nhds ⟨hy1, hy2, hg1', hg2'⟩) w).and
      self_mem_nhdsWithin).exists
    have hspos' : (0 : ℝ) < s := hspos
    have es1 : sda a b (y + s • w) = s * ((b.1 - a.1) * w.2 - (b.2 - a.2) * w.1) := by
      rw [sda_add₃, hy3]; simp only [Prod.smul_fst, Prod.smul_snd, smul_eq_mul]; ring
    have es2 : sda v1 v2 (y + s • w) = s * ((v2.1 - v1.1) * w.2 - (v2.2 - v1.2) * w.1) := by
      rw [sda_add₃, hg3]; simp only [Prod.smul_fst, Prod.smul_snd, smul_eq_mul]; ring
    have i1 := mem_interior_hull3 hD hs.1 hs.2.1 (by rw [es1, mul_left_comm]; exact mul_pos hspos' hw1)
    have i2 := mem_interior_hull3 hD' hs.2.2.1 hs.2.2.2
      (by rw [es2, mul_left_comm]; exact mul_pos hspos' hw2)
    exact Set.disjoint_left.1 hdisj i2 i1
  have hN := lineN_pos hab
  have hE := lineN_pos hv12
  unfold lineN at hN hE
  -- the edges are parallel
  have hX : (b.1 - a.1) * (v2.2 - v1.2) - (b.2 - a.2) * (v2.1 - v1.1) = 0 := by
    obtain ⟨X, hXdef⟩ : ∃ X, X = (b.1 - a.1) * (v2.2 - v1.2) - (b.2 - a.2) * (v2.1 - v1.1) :=
      ⟨_, rfl⟩
    rw [← hXdef]
    by_contra hX
    refine nocommon ((sda a b c * (v2.1 - v1.1) - sda v1 v2 v3 * (b.1 - a.1)) / X,
      (sda a b c * (v2.2 - v1.2) - sda v1 v2 v3 * (b.2 - a.2)) / X) ?_ ?_
    · have : (b.1 - a.1) * ((sda a b c * (v2.2 - v1.2) - sda v1 v2 v3 * (b.2 - a.2)) / X) -
          (b.2 - a.2) * ((sda a b c * (v2.1 - v1.1) - sda v1 v2 v3 * (b.1 - a.1)) / X) =
          sda a b c * X / X := by rw [hXdef]; ring
      simp only; rw [this, mul_div_assoc, div_self hX, mul_one]; exact mul_self_pos.2 hD
    · have : (v2.1 - v1.1) * ((sda a b c * (v2.2 - v1.2) - sda v1 v2 v3 * (b.2 - a.2)) / X) -
          (v2.2 - v1.2) * ((sda a b c * (v2.1 - v1.1) - sda v1 v2 v3 * (b.1 - a.1)) / X) =
          sda v1 v2 v3 * X / X := by rw [hXdef]; ring
      simp only; rw [this, mul_div_assoc, div_self hX, mul_one]; exact mul_self_pos.2 hD'
  -- ... and point in opposite directions relative to the two triangles
  have hdot : sda v1 v2 v3 * sda a b c *
      ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) < 0 := by
    have hle : sda v1 v2 v3 * sda a b c *
        ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) ≤ 0 := by
      by_contra hpos
      push Not at hpos
      refine nocommon (-(sda a b c * (b.2 - a.2)), sda a b c * (b.1 - a.1)) ?_ ?_
      · have : sda a b c * ((b.1 - a.1) * (sda a b c * (b.1 - a.1)) -
            (b.2 - a.2) * -(sda a b c * (b.2 - a.2))) =
            (sda a b c * sda a b c) * ((b.1 - a.1) ^ 2 + (b.2 - a.2) ^ 2) := by ring
        simp only; rw [this]; exact mul_pos (mul_self_pos.2 hD) hN
      · have : sda v1 v2 v3 * ((v2.1 - v1.1) * (sda a b c * (b.1 - a.1)) -
            (v2.2 - v1.2) * -(sda a b c * (b.2 - a.2))) = sda v1 v2 v3 * sda a b c *
            ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) := by ring
        simp only; rw [this]; exact hpos
    have hne0 : (v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2) ≠ 0 := by
      intro h0
      have : ((b.1 - a.1) * (v2.2 - v1.2) - (b.2 - a.2) * (v2.1 - v1.1)) ^ 2 +
          ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) ^ 2 =
          ((b.1 - a.1) ^ 2 + (b.2 - a.2) ^ 2) * ((v2.1 - v1.1) ^ 2 + (v2.2 - v1.2) ^ 2) := by ring
      rw [hX, h0] at this
      have := mul_pos hN hE
      linarith
    exact lt_of_le_of_ne hle (mul_ne_zero (mul_ne_zero hD' hD) hne0)
  -- `v₁` and `v₂` are on the line `ab`
  have hg3' := hg3
  unfold sda at hg3'
  have hq : (b.1 - a.1) * (y.2 - v1.2) - (b.2 - a.2) * (y.1 - v1.1) = 0 := by
    have k1 : ((b.1 - a.1) * (y.2 - v1.2) - (b.2 - a.2) * (y.1 - v1.1)) * (v2.1 - v1.1) = 0 := by
      linear_combination (y.1 - v1.1) * hX + (b.1 - a.1) * hg3'
    have k2 : ((b.1 - a.1) * (y.2 - v1.2) - (b.2 - a.2) * (y.1 - v1.1)) * (v2.2 - v1.2) = 0 := by
      linear_combination (y.2 - v1.2) * hX + (b.2 - a.2) * hg3'
    have : ((b.1 - a.1) * (y.2 - v1.2) - (b.2 - a.2) * (y.1 - v1.1)) *
        ((v2.1 - v1.1) ^ 2 + (v2.2 - v1.2) ^ 2) = 0 := by
      linear_combination (v2.1 - v1.1) * k1 + (v2.2 - v1.2) * k2
    exact (mul_eq_zero.1 this).resolve_right hE.ne'
  have hy3' := hy3
  unfold sda at hy3'
  have la : (b.1 - a.1) * (a.2 - y.2) - (b.2 - a.2) * (a.1 - y.1) = 0 := by linear_combination -hy3'
  have lb : (b.1 - a.1) * (b.2 - y.2) - (b.2 - a.2) * (b.1 - y.1) = 0 := by linear_combination -hy3'
  have l1 : (b.1 - a.1) * (v1.2 - y.2) - (b.2 - a.2) * (v1.1 - y.1) = 0 := by linear_combination -hq
  have l2 : (b.1 - a.1) * (v2.2 - y.2) - (b.2 - a.2) * (v2.1 - y.1) = 0 := by
    linear_combination hX - hq
  have hsv1 : sda a b v1 = 0 := by unfold sda; linear_combination hy3' - hq
  have hsv2 : sda a b v2 = 0 := by unfold sda; linear_combination hy3' + l2
  -- the third corner is across `ab`
  have hv3 : sda a b c * sda a b v3 < 0 := by
    have hid : sda a b v3 * ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) =
        ((b.1 - a.1) ^ 2 + (b.2 - a.2) ^ 2) * sda v1 v2 v3 := by
      unfold sda
      linear_combination ((b.2 - a.2) * (v3.2 - v1.2) + (b.1 - a.1) * (v3.1 - v1.1)) * hX +
        ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) * hy3' -
        ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) * hq
    have : (sda a b c * sda a b v3) * ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) ^ 2 =
        ((b.1 - a.1) ^ 2 + (b.2 - a.2) ^ 2) * (sda v1 v2 v3 * sda a b c *
          ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2))) := by
      linear_combination sda a b c * ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) * hid
    exact neg_of_mul_sq this hN hdot
  -- positions along the line
  have tA : sda a b c * lineT a b y a = -(sda a y c) * lineN a b := by
    unfold lineT lineN sda; linear_combination ((a.1 - c.1) * (b.1 - a.1) + (a.2 - c.2) * (b.2 - a.2)) * hy3'
  have tB : sda a b c * lineT a b y b = sda y b c * lineN a b := by
    unfold lineT lineN sda; linear_combination ((b.1 - c.1) * (b.1 - a.1) + (b.2 - c.2) * (b.2 - a.2)) * hy3'
  have tV1 : sda v1 v2 v3 * lineT a b y v1 = -(sda v1 y v3) *
      ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) := by
    unfold lineT sda; linear_combination ((v1.1 - v3.1) * (b.1 - a.1) + (v1.2 - v3.2) * (b.2 - a.2)) * hg3'
  have tV2 : sda v1 v2 v3 * lineT a b y v2 = sda y v2 v3 *
      ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) := by
    unfold lineT sda; linear_combination ((v2.1 - v3.1) * (b.1 - a.1) + (v2.2 - v3.2) * (b.2 - a.2)) * hg3'
  have hNl : 0 < lineN a b := lineN_pos hab
  have hDD := mul_self_pos.2 hD
  have hDD' := mul_self_pos.2 hD'
  have sA : lineT a b y a < 0 := by
    have : (sda a b c * sda a b c) * lineT a b y a = -(sda a b c * sda a y c) * lineN a b := by
      linear_combination sda a b c * tA
    exact neg_of_pos_mul hDD this (by rw [neg_mul]; exact neg_neg_of_pos (mul_pos hy2 hNl))
  have sB : 0 < lineT a b y b := by
    have : (sda a b c * sda a b c) * lineT a b y b = (sda a b c * sda y b c) * lineN a b := by
      linear_combination sda a b c * tB
    exact pos_of_pos_mul hDD this (mul_pos hy1 hNl)
  have pa := line_param la
  have pb := line_param lb
  have p1 := line_param l1
  have p2 := line_param l2
  have hK := convex_convexHull ℝ ({a, b, c} : Set (ℝ × ℝ))
  have hK' := convex_convexHull ℝ ({v1, v2, v3} : Set (ℝ × ℝ))
  have haK : a ∈ convexHull ℝ ({a, b, c} : Set (ℝ × ℝ)) := subset_convexHull ℝ _ (by simp)
  have hbK : b ∈ convexHull ℝ ({a, b, c} : Set (ℝ × ℝ)) := subset_convexHull ℝ _ (by simp)
  have hv1K : v1 ∈ convexHull ℝ ({v1, v2, v3} : Set (ℝ × ℝ)) := subset_convexHull ℝ _ (by simp)
  have hv2K : v2 ∈ convexHull ℝ ({v1, v2, v3} : Set (ℝ × ℝ)) := subset_convexHull ℝ _ (by simp)
  have vtx : ∀ z ∈ S, z ∈ convexHull ℝ ({a, b, c} : Set (ℝ × ℝ)) → z = a ∨ z = b ∨ z = c := by
    intro z hz hzK
    have := hT.vertex _ hu z hz (by rwa [hull_triple])
    simpa using this
  have vtx' : ∀ z ∈ S, z ∈ convexHull ℝ ({v1, v2, v3} : Set (ℝ × ℝ)) → z = v1 ∨ z = v2 ∨ z = v3 := by
    intro z hz hzK
    have := hT.vertex _ hu' z hz (by rwa [hull_triple])
    simpa using this
  have hc0 : ∀ z, sda a b z = 0 → z ≠ c := by rintro z hz rfl; exact hD hz
  have hv30 : ∀ z, sda a b z = 0 → z ≠ v3 := by
    rintro z hz rfl; rw [hz, mul_zero] at hv3; exact lt_irrefl _ hv3
  have hsa : sda a b a = 0 := sda_self₁₃ _ _
  have hsb : sda a b b = 0 := sda_self₂₃ _ _
  rcases lt_or_gt_of_ne (show (v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2) ≠ 0 by
      rintro h0; rw [h0, mul_zero] at hdot; exact lt_irrefl _ hdot) with hneg | hpos
  · -- `v₁` on `b`'s side, `v₂` on `a`'s side
    have s1 : 0 < lineT a b y v1 := by
      have : (sda v1 v2 v3 * sda v1 v2 v3) * lineT a b y v1 = -(sda v1 v2 v3 * sda v1 y v3) *
          ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) := by
        linear_combination sda v1 v2 v3 * tV1
      exact pos_of_pos_mul hDD' this (by rw [neg_mul]; exact neg_pos.2 (mul_neg_of_pos_of_neg hg2' hneg))
    have s2 : lineT a b y v2 < 0 := by
      have : (sda v1 v2 v3 * sda v1 v2 v3) * lineT a b y v2 = (sda v1 v2 v3 * sda y v2 v3) *
          ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) := by
        linear_combination sda v1 v2 v3 * tV2
      exact neg_of_pos_mul hDD' this (mul_neg_of_pos_of_neg hg1' hneg)
    have e2 : a = v2 := eq_of_same_side hK hK' hNl hyu hyu' haK hv2K pa.1 pa.2 p2.1 p2.2
      (mul_pos_of_neg_of_neg sA s2)
      (fun h => by
        rcases vtx v2 hv2S h with h' | h' | h'
        · exact h'
        · exfalso; have := congrArg (lineT a b y) h'; linarith
        · exact absurd h' (hc0 v2 hsv2))
      (fun h => by
        rcases vtx' a haS h with h' | h' | h'
        · exfalso; have := congrArg (lineT a b y) h'; linarith
        · exact h'
        · exact absurd h' (hv30 a hsa))
    have e1 : b = v1 := eq_of_same_side hK hK' hNl hyu hyu' hbK hv1K pb.1 pb.2 p1.1 p1.2
      (mul_pos sB s1)
      (fun h => by
        rcases vtx v1 hv1S h with h' | h' | h'
        · exfalso; have := congrArg (lineT a b y) h'; linarith
        · exact h'
        · exact absurd h' (hc0 v1 hsv1))
      (fun h => by
        rcases vtx' b hbS h with h' | h' | h'
        · exact h'
        · exfalso; have := congrArg (lineT a b y) h'; linarith
        · exact absurd h' (hv30 b hsb))
    refine ⟨v3, ?_, hv3⟩
    rw [← e1, ← e2, Finset.insert_comm]
  · -- `v₁` on `a`'s side, `v₂` on `b`'s side
    have s1 : lineT a b y v1 < 0 := by
      have : (sda v1 v2 v3 * sda v1 v2 v3) * lineT a b y v1 = -(sda v1 v2 v3 * sda v1 y v3) *
          ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) := by
        linear_combination sda v1 v2 v3 * tV1
      exact neg_of_pos_mul hDD' this (by rw [neg_mul]; exact neg_neg_of_pos (mul_pos hg2' hpos))
    have s2 : 0 < lineT a b y v2 := by
      have : (sda v1 v2 v3 * sda v1 v2 v3) * lineT a b y v2 = (sda v1 v2 v3 * sda y v2 v3) *
          ((v2.1 - v1.1) * (b.1 - a.1) + (v2.2 - v1.2) * (b.2 - a.2)) := by
        linear_combination sda v1 v2 v3 * tV2
      exact pos_of_pos_mul hDD' this (mul_pos hg1' hpos)
    have e1 : a = v1 := eq_of_same_side hK hK' hNl hyu hyu' haK hv1K pa.1 pa.2 p1.1 p1.2
      (mul_pos_of_neg_of_neg sA s1)
      (fun h => by
        rcases vtx v1 hv1S h with h' | h' | h'
        · exact h'
        · exfalso; have := congrArg (lineT a b y) h'; linarith
        · exact absurd h' (hc0 v1 hsv1))
      (fun h => by
        rcases vtx' a haS h with h' | h' | h'
        · exact h'
        · exfalso; have := congrArg (lineT a b y) h'; linarith
        · exact absurd h' (hv30 a hsa))
    have e2 : b = v2 := eq_of_same_side hK hK' hNl hyu hyu' hbK hv2K pb.1 pb.2 p2.1 p2.2
      (mul_pos sB s2)
      (fun h => by
        rcases vtx v2 hv2S h with h' | h' | h'
        · exfalso; have := congrArg (lineT a b y) h'; linarith
        · exact h'
        · exact absurd h' (hc0 v2 hsv2))
      (fun h => by
        rcases vtx' b hbS h with h' | h' | h'
        · exfalso; have := congrArg (lineT a b y) h'; linarith
        · exact h'
        · exact absurd h' (hv30 b hsb))
    exact ⟨v3, by rw [← e1, ← e2], hv3⟩

/-- **Adjacency.** Let `y` be a relative interior point of the edge `ab` of the triangle
`{a, b, c}` of a triangulation of `S`, not itself a point of `S`. Every other triangle of the
triangulation that contains `y` is `{a, b, d}` with `d` strictly on the other side of `ab`.
(From disjoint interiors: the other triangle cannot reach into `{a, b, c}` near `y`, so its edge
through `y` lies on the line `ab` and its third corner across; from the vertex condition: its
two corners on that line are exactly `a` and `b`.) -/
theorem eq_insert_of_mem_edge {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {a b c y : ℝ × ℝ} (hu : {a, b, c} ∈ T)
    (hy1 : 0 < sda a b c * sda y b c) (hy2 : 0 < sda a b c * sda a y c) (hy3 : sda a b y = 0)
    (hyS : y ∉ S) {u' : Finset (ℝ × ℝ)} (hu' : u' ∈ T) (hne : u' ≠ {a, b, c}) (hyu' : y ∈ hull u') :
    ∃ d, u' = {a, b, d} ∧ sda a b c * sda a b d < 0 := by
  have hD := hT.sda_ne_zero hu
  obtain ⟨v1, v2, v3, hD', rfl⟩ := hT.tri u' hu'
  rw [hull_triple] at hyu'
  obtain ⟨g1, g2, g3⟩ := (mem_hull3_iff hD').1 hyu'
  have hyu : y ∈ convexHull ℝ {a, b, c} := (mem_hull3_iff hD).2 ⟨hy1.le, hy2.le, by rw [hy3, mul_zero]⟩
  have hnotint : ¬ (0 < sda v1 v2 v3 * sda y v2 v3 ∧ 0 < sda v1 v2 v3 * sda v1 y v3 ∧
      0 < sda v1 v2 v3 * sda v1 v2 y) := by
    rintro ⟨h1, h2, h3⟩
    have hi := mem_interior_hull3 hD' h1 h2 h3
    have hcl := hull3_subset_closure_interior hD hyu
    rw [mem_closure_iff] at hcl
    obtain ⟨z, hz1, hz2⟩ := hcl _ isOpen_interior hi
    have : Disjoint (interior (hull {v1, v2, v3})) (interior (hull {a, b, c})) := hT.disjoint hu' hu hne
    rw [hull_triple, hull_triple] at this
    exact Set.disjoint_left.1 this hz1 hz2
  have zero : ∀ {t : ℝ}, 0 ≤ sda v1 v2 v3 * t → ¬ 0 < sda v1 v2 v3 * t → t = 0 := fun h h' =>
    (mul_eq_zero.1 (le_antisymm (not_lt.1 h') h)).resolve_left hD'
  by_cases h3 : 0 < sda v1 v2 v3 * sda v1 v2 y
  · by_cases h1 : 0 < sda v1 v2 v3 * sda y v2 v3
    · -- the second coordinate vanishes: rotate to `v₃ v₁ v₂`
      have h2 := zero g2 fun h2 => hnotint ⟨h1, h2, h3⟩
      have e : ({v1, v2, v3} : Finset (ℝ × ℝ)) = {v3, v1, v2} := by ext; simp; tauto
      have r : sda v3 v1 v2 = sda v1 v2 v3 := by unfold sda; ring
      obtain ⟨d, hd, hs⟩ := eq_insert_of_mem_edge_aux hT hu hy1 hy2 hy3 hyS (e ▸ hu') (e ▸ hne)
        (by rw [r, show sda y v1 v2 = sda v1 v2 y by unfold sda; ring]; exact g3)
        (by rw [r, show sda v3 y v2 = sda y v2 v3 by unfold sda; ring]; exact g1)
        (by rw [show sda v3 v1 y = sda v1 y v3 by unfold sda; ring]; exact h2)
      exact ⟨d, e.trans hd, hs⟩
    · -- the first coordinate vanishes: rotate to `v₂ v₃ v₁`
      have h1 := zero g1 h1
      have e : ({v1, v2, v3} : Finset (ℝ × ℝ)) = {v2, v3, v1} := by ext; simp; tauto
      have r : sda v2 v3 v1 = sda v1 v2 v3 := by unfold sda; ring
      obtain ⟨d, hd, hs⟩ := eq_insert_of_mem_edge_aux hT hu hy1 hy2 hy3 hyS (e ▸ hu') (e ▸ hne)
        (by rw [r, show sda y v3 v1 = sda v1 y v3 by unfold sda; ring]; exact g2)
        (by rw [r, show sda v2 y v1 = sda v1 v2 y by unfold sda; ring]; exact g3)
        (by rw [show sda v2 v3 y = sda y v2 v3 by unfold sda; ring]; exact h1)
      exact ⟨d, e.trans hd, hs⟩
  · exact eq_insert_of_mem_edge_aux hT hu hy1 hy2 hy3 hyS hu' hne g1 g2 (zero g3 h3)

/-! ## A generic starting point for the walk -/

theorem interior_biUnion_eq_empty {ι : Type*} (G : Finset ι) (L : ι → Set (ℝ × ℝ))
    (hc : ∀ i ∈ G, IsClosed (L i)) (he : ∀ i ∈ G, interior (L i) = ∅) :
    interior (⋃ i ∈ G, L i) = ∅ := by
  classical
  induction G using Finset.induction_on with
  | empty => simp
  | insert i G hi IH =>
    rw [Finset.set_biUnion_insert, Set.union_comm,
      interior_union_isClosed_of_interior_empty
        (isClosed_biUnion_finset fun j hj => hc j (Finset.mem_insert_of_mem hj))
        (he i (Finset.mem_insert_self _ _))]
    exact IH (fun j hj => hc j (Finset.mem_insert_of_mem hj)) (fun j hj => he j (Finset.mem_insert_of_mem hj))

/-- A line through `p` and another point has empty interior. -/
theorem interior_line_eq_empty {p s : ℝ × ℝ} (hps : s ≠ p) :
    interior {z : ℝ × ℝ | sda z p s = 0} = ∅ := by
  rw [Set.eq_empty_iff_forall_notMem]
  intro z hz
  obtain ⟨τ, hτ, hτpos⟩ := ((exists_pos_near (isOpen_interior.mem_nhds hz) (s.2 - p.2, p.1 - s.1)).and
    self_mem_nhdsWithin).exists
  have hτpos' : (0 : ℝ) < τ := hτpos
  have h0 : sda z p s = 0 := by have h := interior_subset hz; exact h
  have h1 : sda (z + τ • (s.2 - p.2, p.1 - s.1)) p s = 0 := by have h := interior_subset hτ; exact h
  have e : sda (z + τ • (s.2 - p.2, p.1 - s.1)) p s = sda z p s - τ * lineN s p := by
    unfold sda lineN; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
      smul_eq_mul]; ring
  rw [e, h0] at h1
  have := mul_pos hτpos' (lineN_pos hps)
  linarith

/-- **A generic point.** Inside a nondegenerate triangle there is a point `x` such that the line
`xp` meets none of the finitely many points of `F ∌ p`. -/
theorem exists_generic_interior {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) (p : ℝ × ℝ)
    (F : Finset (ℝ × ℝ)) (hF : p ∉ F) :
    ∃ x ∈ interior (convexHull ℝ {a, b, c}), ∀ s ∈ F, sda x p s ≠ 0 := by
  by_contra hno
  push Not at hno
  have hsub : interior (convexHull ℝ ({a, b, c} : Set (ℝ × ℝ))) ⊆
      interior (⋃ s ∈ F, {z : ℝ × ℝ | sda z p s = 0}) := by
    refine interior_maximal (fun x hx => ?_) isOpen_interior
    obtain ⟨s, hs, h0⟩ := hno x hx
    exact Set.mem_biUnion hs h0
  rw [interior_biUnion_eq_empty F _ (fun s _ => isClosed_eq (continuous_sda₁ p s) continuous_const)
    (fun s hs => interior_line_eq_empty fun h => hF (h ▸ hs))] at hsub
  obtain ⟨x, hx⟩ := interior_hull3_nonempty h
  exact hsub hx

/-! ## Leaving a triangle -/

theorem eventually_nonneg_along {f : ℝ × ℝ → ℝ} (hf : Continuous f)
    (hfa : ∀ u v : ℝ × ℝ, ∀ s t : ℝ, s + t = 1 → f (s • u + t • v) = s * f u + t * f v)
    {γ : ℝ → ℝ × ℝ} {x p : ℝ × ℝ} (hγ : ∀ τ, γ τ = x + τ • (p - x)) {τ' : ℝ} (hτ : τ' < 1)
    (hy : 0 ≤ f (γ τ')) (hp : f (γ τ') = 0 → 0 ≤ f p) : ∀ᶠ τ in 𝓝[>] τ', 0 ≤ f (γ τ) := by
  rcases hy.lt_or_eq with hpos | hzero
  · have hc : Continuous fun τ => f (γ τ) := by
      have : γ = fun τ => x + τ • (p - x) := funext hγ
      rw [this]; fun_prop
    exact ((hc.tendsto τ').eventually (lt_mem_nhds hpos)).filter_mono nhdsWithin_le_nhds |>.mono
      fun _ h => h.le
  · have hfp := hp hzero.symm
    filter_upwards [Ioo_mem_nhdsGT hτ] with τ hτ'
    have h1 : (0 : ℝ) < 1 - τ' := by linarith
    have e : γ τ = ((1 - τ) / (1 - τ')) • γ τ' + ((τ - τ') / (1 - τ')) • p := by
      rw [hγ, hγ]
      ext <;> simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, Prod.fst_sub,
        Prod.snd_sub, smul_eq_mul] <;> field_simp <;> ring
    rw [e, hfa _ _ _ _ (by field_simp; ring), ← hzero]
    have : 0 ≤ (τ - τ') / (1 - τ') := div_nonneg (by linarith [hτ'.1]) h1.le
    nlinarith

/-- **Leaving a triangle through an edge.** If the walk `γ` is in the triangle `a b c` at `τ'` and
outside it on `(τ', 1]`, then at `τ'` it is on one of the three edges, and `p = γ 1` is strictly
beyond that edge. -/
theorem exit_edge {a b c x p : ℝ × ℝ} (h : sda a b c ≠ 0) {γ : ℝ → ℝ × ℝ}
    (hγ : ∀ τ, γ τ = x + τ • (p - x)) {τ' : ℝ} (hτ : τ' < 1) (hy : γ τ' ∈ convexHull ℝ {a, b, c})
    (hout : ∀ τ ∈ Ioc τ' 1, γ τ ∉ convexHull ℝ {a, b, c}) :
    (sda a b (γ τ') = 0 ∧ sda a b c * sda a b p < 0) ∨
      (sda (γ τ') b c = 0 ∧ sda a b c * sda p b c < 0) ∨
      (sda a (γ τ') c = 0 ∧ sda a b c * sda a p c < 0) := by
  by_contra hno
  push Not at hno
  obtain ⟨n3, n1, n2⟩ := hno
  obtain ⟨y1, y2, y3⟩ := (mem_hull3_iff h).1 hy
  have zero : ∀ {t : ℝ}, sda a b c * t = 0 → t = 0 := fun h0 => (mul_eq_zero.1 h0).resolve_left h
  have e1 := eventually_nonneg_along (f := fun z => sda a b c * sda z b c)
    (continuous_const.mul (continuous_sda₁ b c))
    (fun u v s t hst => by rw [sda_combo₁ _ _ _ _ hst]; ring) hγ hτ y1
    (fun h0 => n1 (zero h0))
  have e2 := eventually_nonneg_along (f := fun z => sda a b c * sda a z c)
    (continuous_const.mul (continuous_sda₂ a c))
    (fun u v s t hst => by rw [sda_combo₂ _ _ _ _ hst]; ring) hγ hτ y2
    (fun h0 => n2 (zero h0))
  have e3 := eventually_nonneg_along (f := fun z => sda a b c * sda a b z)
    (continuous_const.mul (continuous_sda₃ a b))
    (fun u v s t hst => by rw [sda_combo₃ _ _ _ _ hst]; ring) hγ hτ y3
    (fun h0 => n3 (zero h0))
  obtain ⟨τ, ⟨h1, h2, h3⟩, hτ'⟩ := ((e1.and (e2.and e3)).and (Ioo_mem_nhdsGT hτ)).exists
  exact hout τ ⟨hτ'.1, hτ'.2.le⟩ ((mem_hull3_iff h).2 ⟨h1, h2, h3⟩)

/-! ## Delaunay's lemma -/

/-- **Delaunay's lemma.** A triangulation of a finite point set `S` that is locally Delaunay at
every interior edge is Delaunay: no point of `S` lies strictly inside the circumcircle of any of
its triangles.

[formula: delaunay-lemma] -/
theorem delaunay_lemma {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (hloc : IsLocallyDelaunay T) : IsDelaunay S T := by
  intro a0 b0 c0 ht0 p hp hin
  have hD0 := hT.sda_ne_zero ht0
  -- a point of `S` strictly inside a circumcircle is not in the triangle
  have notin : ∀ a b c, {a, b, c} ∈ T → 0 < circleLift a b c p → p ∉ convexHull ℝ {a, b, c} := by
    intro a b c hu hpos hpm
    have := hT.vertex _ hu p hp (by rwa [hull_triple])
    simp only [Finset.mem_insert, Finset.mem_singleton] at this
    rcases this with rfl | rfl | rfl <;> simp [circleLift, inCircle_self₁, inCircle_self₂,
      inCircle_self₃] at hpos
  have hφ0 : 0 < circleLift a0 b0 c0 p := (circleLift_pos_iff hD0 p).2 hin
  obtain ⟨x, hx, hgen⟩ := exists_generic_interior hD0 p (S.erase p) (Finset.notMem_erase p S)
  obtain ⟨γ, hγ⟩ : ∃ γ : ℝ → ℝ × ℝ, ∀ τ, γ τ = x + τ • (p - x) := ⟨_, fun _ => rfl⟩
  have hγf : γ = fun τ => x + τ • (p - x) := funext hγ
  have hγc : Continuous γ := by rw [hγf]; fun_prop
  have hxp : x ≠ p := fun h => notin a0 b0 c0 ht0 hφ0 (h ▸ interior_subset hx)
  have hxS : x ∈ convexHull ℝ (S : Set (ℝ × ℝ)) :=
    hT.hull_subset ht0 (by rw [hull_triple]; exact interior_subset hx)
  have hpS : p ∈ convexHull ℝ (S : Set (ℝ × ℝ)) := subset_convexHull ℝ _ hp
  have hγS : ∀ τ ∈ Icc (0 : ℝ) 1, γ τ ∈ convexHull ℝ (S : Set (ℝ × ℝ)) := by
    intro τ hτ
    have e : γ τ = (1 - τ) • x + τ • p := by
      rw [hγ, smul_sub, sub_smul, one_smul]; abel
    rw [e]
    exact convex_convexHull ℝ _ hxS hpS (by linarith [hτ.2]) hτ.1 (by ring)
  have hγp : ∀ τ, γ τ = p → τ = 1 := by
    intro τ h
    by_contra hne
    have h1 := congrArg Prod.fst h
    have h2 := congrArg Prod.snd h
    simp only [hγ, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, Prod.fst_sub,
      Prod.snd_sub, smul_eq_mul] at h1 h2
    have hne' : (1 - τ) ≠ 0 := sub_ne_zero.2 (Ne.symm hne)
    apply hxp
    ext
    · have : (1 - τ) * (x.1 - p.1) = 0 := by linarith
      linarith [(mul_eq_zero.1 this).resolve_left hne']
    · have : (1 - τ) * (x.2 - p.2) = 0 := by linarith
      linarith [(mul_eq_zero.1 this).resolve_left hne']
  have hγS' : ∀ τ < 1, γ τ ∉ S := by
    intro τ hτ hmem
    have hne : γ τ ≠ p := fun h => absurd (hγp τ h) hτ.ne
    refine hgen (γ τ) (Finset.mem_erase.2 ⟨hne, hmem⟩) ?_
    rw [hγ]; unfold sda; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
      Prod.fst_sub, Prod.snd_sub, smul_eq_mul]; ring
  -- the last parameter at which the walk is in a triangle at level `φ₀` or higher
  set φ0 := circleLift a0 b0 c0 p
  let Q : Set ((ℝ × ℝ) × (ℝ × ℝ) × (ℝ × ℝ)) :=
    {q | ({q.1, q.2.1, q.2.2} : Finset (ℝ × ℝ)) ∈ T ∧ φ0 ≤ circleLift q.1 q.2.1 q.2.2 p}
  have hQ : Q.Finite := by
    refine (Finset.finite_toSet (S ×ˢ (S ×ˢ S))).subset fun q hq => ?_
    simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe]
    exact ⟨hT.mem_S hq.1 (by simp), hT.mem_S hq.1 (by simp), hT.mem_S hq.1 (by simp)⟩
  let J : Set ℝ := Icc 0 1 ∩ ⋃ q ∈ Q, γ ⁻¹' convexHull ℝ {q.1, q.2.1, q.2.2}
  have hJc : IsClosed J :=
    isClosed_Icc.inter (hQ.isClosed_biUnion fun q _ => (isClosed_hull3 _ _ _).preimage hγc)
  have hJcpt : IsCompact J := isCompact_Icc.of_isClosed_subset hJc Set.inter_subset_left
  have memJ : ∀ τ a b c, τ ∈ Icc (0 : ℝ) 1 → {a, b, c} ∈ T → φ0 ≤ circleLift a b c p →
      γ τ ∈ convexHull ℝ {a, b, c} → τ ∈ J := fun τ a b c hτ hu hφ hm =>
    ⟨hτ, Set.mem_biUnion (x := (a, b, c)) ⟨hu, hφ⟩ hm⟩
  have h0J : (0 : ℝ) ∈ J := memJ 0 a0 b0 c0 ⟨le_rfl, zero_le_one⟩ ht0 le_rfl
    (by rw [hγ, zero_smul, add_zero]; exact interior_subset hx)
  have hτJ : sSup J ∈ J := hJcpt.sSup_mem ⟨0, h0J⟩
  have hmax : ∀ τ ∈ J, τ ≤ sSup J := fun τ hτ => le_csSup hJcpt.bddAbove hτ
  set τ' := sSup J
  obtain ⟨⟨hτ0, hτ1⟩, hτU⟩ := hτJ
  obtain ⟨⟨a, b, c⟩, ⟨hu, hφ⟩, hyu⟩ := Set.mem_iUnion₂.1 hτU
  simp only [Set.mem_preimage] at hyu
  simp only at hu hφ
  have hout : ∀ a b c, {a, b, c} ∈ T → φ0 ≤ circleLift a b c p →
      ∀ τ ∈ Ioc τ' 1, γ τ ∉ convexHull ℝ {a, b, c} := fun a b c hu hφ τ hτ hm =>
    not_le.2 hτ.1 (hmax τ (memJ τ a b c ⟨by linarith [hτ.1], hτ.2⟩ hu hφ hm))
  -- one step of the walk is impossible
  have step : ∀ a b c, {a, b, c} ∈ T → φ0 ≤ circleLift a b c p → γ τ' ∈ convexHull ℝ {a, b, c} →
      τ' < 1 → sda a b (γ τ') = 0 → sda a b c * sda a b p < 0 → False := by
    intro a b c hu hφ hyu hτ1 hy3 hpside
    have hD := hT.sda_ne_zero hu
    have hyS := hγS' τ' hτ1
    obtain ⟨hy1, hy2, -⟩ := (mem_hull3_iff hD).1 hyu
    have hy1' : 0 < sda a b c * sda (γ τ') b c := lt_of_le_of_ne hy1 fun h0 => hyS (by
      rw [eq_of_sda_eq_zero hD ((mul_eq_zero.1 h0.symm).resolve_left hD) hy3]
      exact hT.mem_S hu (by simp))
    have hy2' : 0 < sda a b c * sda a (γ τ') c := lt_of_le_of_ne hy2 fun h0 => hyS (by
      rw [eq_of_sda_eq_zero' hD ((mul_eq_zero.1 h0.symm).resolve_left hD) hy3]
      exact hT.mem_S hu (by simp))
    -- the next triangle along the walk
    obtain ⟨u', hu', hyu', τ'', hτ'', hu''⟩ :
        ∃ u' ∈ T, γ τ' ∈ hull u' ∧ ∃ τ'' ∈ Ioc τ' 1, γ τ'' ∈ hull u' := by
      have hcov : Ioc τ' 1 ⊆ ⋃ u ∈ (T : Set (Finset (ℝ × ℝ))), γ ⁻¹' hull u ∩ Ioc τ' 1 := by
        intro τ hτ
        have := hγS τ ⟨by linarith [hτ.1], hτ.2⟩
        rw [← hT.cover] at this
        obtain ⟨u, hu, hm⟩ := Set.mem_iUnion₂.1 this
        exact Set.mem_biUnion (Finset.mem_coe.2 hu) ⟨hm, hτ⟩
      have hcl : τ' ∈ closure (⋃ u ∈ (T : Set (Finset (ℝ × ℝ))), γ ⁻¹' hull u ∩ Ioc τ' 1) :=
        closure_mono hcov (by rw [closure_Ioc hτ1.ne]; exact ⟨le_rfl, hτ1.le⟩)
      rw [T.finite_toSet.closure_biUnion] at hcl
      obtain ⟨u, hu, hcu⟩ := Set.mem_iUnion₂.1 hcl
      have hne : (γ ⁻¹' hull u ∩ Ioc τ' 1).Nonempty := by
        by_contra hemp
        rw [Set.not_nonempty_iff_eq_empty] at hemp
        rw [hemp, closure_empty] at hcu
        exact hcu
      obtain ⟨τ'', hm, hτ''⟩ := hne
      have hcA : closure (γ ⁻¹' hull u ∩ Ioc τ' 1) ⊆ γ ⁻¹' hull u :=
        closure_minimal Set.inter_subset_left
          (((Set.toFinite _).isCompact_convexHull (𝕜 := ℝ)).isClosed.preimage hγc)
      exact ⟨u, hu, hcA hcu, τ'', hτ'', hm⟩
    have hne : u' ≠ {a, b, c} := by
      rintro rfl
      rw [hull_triple] at hu''
      exact hout a b c hu hφ τ'' hτ'' hu''
    obtain ⟨d, rfl, hd⟩ := eq_insert_of_mem_edge hT hu hy1' hy2' hy3 hyS hu' hne hyu'
    have hcd : c ≠ d := by
      rintro rfl
      exact absurd hd (not_lt.2 (mul_self_nonneg _))
    have hnot := hloc a b c d hu hu' hcd
    have hpd : 0 < sda a b d * sda a b p := by
      have := mul_pos_of_neg_of_neg hd hpside
      have e : sda a b c * sda a b d * (sda a b c * sda a b p) =
          (sda a b c * sda a b c) * (sda a b d * sda a b p) := by ring
      rw [e] at this
      exact pos_of_mul_pos_right this (mul_self_nonneg _)
    have hle := circleLift_le_of_not_inCircumcircle hd hnot hpd
    have hJ'' : τ'' ∈ J := memJ τ'' a b d ⟨by linarith [hτ''.1], hτ''.2⟩ hu' (hφ.trans hle)
      (by rwa [hull_triple] at hu'')
    exact absurd (hmax τ'' hJ'') (not_le.2 hτ''.1)
  -- the walk does not reach `p` inside a triangle at level `φ₀`
  have hτ1' : τ' < 1 := by
    refine lt_of_le_of_ne hτ1 fun h => notin a b c hu (hφ0.trans_le hφ) ?_
    have : γ τ' = p := by rw [hγ, h, one_smul]; abel
    rw [this] at hyu; exact hyu
  have hDabc := hT.sda_ne_zero hu
  rcases exit_edge hDabc hγ hτ1' hyu (hout a b c hu hφ) with ⟨h3, hs⟩ | ⟨h1, hs⟩ | ⟨h2, hs⟩
  · exact step a b c hu hφ hyu hτ1' h3 hs
  · refine step b c a (by convert hu using 1; ext; simp; tauto) (by rwa [circleLift_rotate])
      (by rw [← hull_rot a b c]; exact hyu) hτ1' (by rw [← h1]; unfold sda; ring) ?_
    have : sda b c a * sda b c p = sda a b c * sda p b c := by unfold sda; ring
    rwa [this]
  · refine step c a b (by convert hu using 1; ext; simp; tauto)
      (by rw [circleLift_rotate b c a, circleLift_rotate a b c]; exact hφ)
      (by rw [← hull_rot b c a, ← hull_rot a b c]; exact hyu)
      hτ1' (by rw [← h2]; unfold sda; ring) ?_
    have : sda c a b * sda c a p = sda a b c * sda a p c := by unfold sda; ring
    rwa [this]

/-- For a triangulation, Delaunay and locally Delaunay are the same. -/
theorem isDelaunay_iff_isLocallyDelaunay {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) : IsDelaunay S T ↔ IsLocallyDelaunay T :=
  ⟨IsDelaunay.isLocallyDelaunay hT, delaunay_lemma hT⟩

/-- **Delaunay = empty circumdisks**: `T` is Delaunay iff no point of `S` is closer to the
circumcenter of a triangle than its corners. -/
theorem isDelaunay_iff_empty_circumdisk {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) :
    IsDelaunay S T ↔ ∀ a b c, {a, b, c} ∈ T → ∀ p ∈ S,
      sqDist2 (circumcenter a b c) a ≤ sqDist2 (circumcenter a b c) p := by
  refine forall₄_congr fun a b c hu => forall₂_congr fun p _ => ?_
  rw [inCircumcircle_iff_lt (hT.sda_ne_zero hu), not_lt]

/-- **Lawson's algorithm computes a Delaunay triangulation**: from every triangulation of `S`, a
sequence of Delaunay flips reaches a Delaunay triangulation of `S`.

[formula: planar-flip-termination] -/
theorem exists_flips_to_delaunay {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) :
    ∃ T', Relation.ReflTransGen DelaunayFlip T T' ∧ IsTriangulation S T' ∧ IsDelaunay S T' := by
  obtain ⟨T', h, hT', hloc⟩ := exists_flips_to_locallyDelaunay hT
  exact ⟨T', h, hT', delaunay_lemma hT' hloc⟩

/-- **Existence of a Delaunay triangulation**, given a triangulation of `S` (a hypothesis here). -/
theorem exists_delaunay_of_triangulation {S : Finset (ℝ × ℝ)}
    (h : ∃ T, IsTriangulation S T) : ∃ T, IsTriangulation S T ∧ IsDelaunay S T := by
  obtain ⟨T, hT⟩ := h
  obtain ⟨T', -, hT', hD⟩ := exists_flips_to_delaunay hT
  exact ⟨T', hT', hD⟩

/-- A single nondegenerate triangle is a triangulation of its three corners (so the hypotheses
are satisfiable), and it is Delaunay. -/
theorem isTriangulation_triangle {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) :
    IsTriangulation {a, b, c} {{a, b, c}} := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro t ht; rw [Finset.mem_singleton.1 ht]; exact ⟨a, b, c, h, rfl⟩
  · intro t ht; rw [Finset.mem_singleton.1 ht]
  · simp [hull]
  · simp
  · intro t ht s hs _; rwa [Finset.mem_singleton.1 ht]

theorem isDelaunay_triangle {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) :
    IsDelaunay {a, b, c} {{a, b, c}} := by
  have hT := isTriangulation_triangle h
  refine delaunay_lemma hT fun i j k l h1 h2 hkl => ?_
  exact absurd ((Finset.mem_singleton.1 h1).trans (Finset.mem_singleton.1 h2).symm)
    (IsTriangulation.triple_ne hkl (hT.sda_ne_zero h2))

end Planar

end Delaunay

end DdgProofs
