import DDG.Intrinsic.Barycentric
import DDG.Intrinsic.VertexRemoval

/-!
# Intrinsic simplification: removing a flat vertex (paper Sec. 4.2, Sec. 7.1 "Vertex removal")

`crates/ddg-intrinsic-simplification/SPEC.md` §3.5, S13, S14, S16, §13; the paper Sec. 4.2 and
its insets; Gillespie, Sharp, Crane 2021, App. D, Thm. D.1.

Plane points are pairs of reals; `dist2`, `sda` (twice the signed area) are those of
`DDG/Geometry/Plane.lean`, `Intrinsic.barycentric` the barycentric coordinates of
`DDG/Intrinsic/Barycentric.lean` (`it-barycentric-coordinates`).
`polar r t = (r cos t, r sin t)`. The removed vertex `i` is laid out at the origin, its first
neighbour on the positive `x`-axis, the others at the angles accumulated from the corner angles
at `i` computed from the lengths (law of cosines).

* **`iem-degree-three-barycentric`** — a flat degree-3 star (faces `iab`, `ibc`, `ica`, all strict
  triangles, angle sum `2π` at `i`) lays out in the plane with the given lengths and the three
  faces counter-clockwise (`degreeThree_layout`); the merged face `abc` is counter-clockwise with
  the outer lengths, and the barycentric coordinates of `i` in it are the Heron areas of the faces
  `ibc`, `ica`, `iab` divided by their sum (`degreeThree_barycentric`), all positive. A point of
  one of the three faces is re-expressed in `abc` by linear substitution, with non-negative
  coordinates summing to one (`subst_barycentric`).
* **`iem-boundary-merge`** — a flat boundary vertex of degree 2 (faces `iab`, `ica`, angle sum `π`
  at `i`, `b` and `c` its boundary neighbours) lies on the segment `bc`, at
  `(ℓ_ic b + ℓ_ib c)/(ℓ_ib + ℓ_ic)`; `ℓ_bc = ℓ_ib + ℓ_ic`; the merged face `abc` is
  counter-clockwise; the barycentric coordinates of `i` in it are
  `(0, ℓ_ic/(ℓ_ib + ℓ_ic), ℓ_ib/(ℓ_ib + ℓ_ic))` (`boundaryMerge`).
* **`iem-flip-to-degree-three`** — interior, simplicial: `ic-vertex-removal` (Thm. D.1,
  `flatStar_exists_flippable`, a restatement). The **boundary analogue** (`FlatHalfStar`): a flat
  boundary vertex (angle sum `π`) of degree `d ≥ 3` has an interior spoke whose quadrilateral is
  **strictly** convex (`FlatHalfStar.exists_flippable`): every inner corner is below `π` (the
  angles at `i` sum to `π`), and the `d − 1` outer corners sum to less than `(d − 1)π`. Unlike the
  interior case (the cross star, `DDG.Intrinsic.cross_star_no_strict`) strict convexity, the flip
  condition of SPEC §3.2, always holds for one of them. In a simplicial half-star the flip of spoke
  `i j_k` adds the edge `j_{k−1} j_{k+1}`, not incident to `i`: the degree drops by one, so
  repeated flips reach degree 2 (the combinatorial half, as for Thm. D.1; Δ-complexes:
  `ic-vertex-removal-delta`, open).
-/

noncomputable section

namespace DDG

namespace IntrinsicSimplification

open Real Finset

/-! ## Plane layouts -/

/-- The point at distance `r` and angle `t` from the origin. -/
def polar (r t : ℝ) : ℝ × ℝ := (r * cos t, r * sin t)

theorem dist2_polar (r t s u : ℝ) :
    dist2 (polar r t) (polar s u) = √(r ^ 2 + s ^ 2 - 2 * r * s * cos (u - t)) := by
  unfold dist2 polar
  congr 1
  rw [cos_sub]
  have h1 := sin_sq_add_cos_sq t
  have h2 := sin_sq_add_cos_sq u
  linear_combination r ^ 2 * h1 + s ^ 2 * h2

theorem dist2_origin_polar {r : ℝ} (hr : 0 ≤ r) (t : ℝ) : dist2 (0, 0) (polar r t) = r := by
  unfold dist2 polar
  have : (r * cos t - 0) ^ 2 + (r * sin t - 0) ^ 2 = r ^ 2 := by
    linear_combination r ^ 2 * sin_sq_add_cos_sq t
  rw [this, sqrt_sq hr]

theorem sda_origin_polar (r t s u : ℝ) :
    sda (0, 0) (polar r t) (polar s u) = r * s * sin (u - t) := by
  unfold sda polar
  rw [sin_sub]
  ring

/-- The law of cosines from the lengths: `√(r² + s² − 2rs cos θ) = ℓ` for the angle `θ` between
the sides `r`, `s` of a triangle with the third side `ℓ`. -/
theorem law_of_cosines {l r s : ℝ} (h : IsTriangle l r s) :
    √(r ^ 2 + s ^ 2 - 2 * r * s * cos (cornerAngle l r s)) = l := by
  have hr := h.pos_b
  have hs := h.pos_c
  rw [cos_cornerAngle h]
  have : r ^ 2 + s ^ 2 - 2 * r * s * ((r ^ 2 + s ^ 2 - l ^ 2) / (2 * r * s)) = l ^ 2 := by
    field_simp; ring
  rw [this, sqrt_sq h.pos_a.le]

theorem sin_cornerAngle_pos {l r s : ℝ} (h : IsTriangle l r s) : 0 < sin (cornerAngle l r s) :=
  sin_pos_of_pos_of_lt_pi (cornerAngle_pos h) (cornerAngle_lt_pi h)

/-- `r s sin θ = 2·Area` for the angle `θ` between the sides `r`, `s`. -/
theorem mul_sin_cornerAngle {l r s : ℝ} (h : IsTriangle l r s) :
    r * s * sin (cornerAngle l r s) = 2 * triArea l r s := by
  have hr := h.pos_b
  have hs := h.pos_c
  rw [sin_cornerAngle h]
  field_simp

/-! ## `iem-degree-three-barycentric` -/

section DegreeThree

variable {ra rb rc ab bc ca : ℝ}

/-- The layout of a degree-3 star: `i` at the origin, `a` on the `x`-axis, `b` and `c` at the
accumulated angles. -/
def layoutB (ra rb ab : ℝ) : ℝ × ℝ := polar rb (cornerAngle ab ra rb)
/-- see `layoutB` -/
def layoutC (ra rb rc ab bc : ℝ) : ℝ × ℝ := polar rc (cornerAngle ab ra rb + cornerAngle bc rb rc)

/-- **The layout of a flat degree-3 star.** With the faces `iab`, `ibc`, `ica` strict triangles
(spokes `ra rb rc`, outer sides `ab bc ca`) and the angle sum `2π` at `i`, the points
`P = 0, A = (ra, 0), B, C` have the six given distances and the three faces are
counter-clockwise, with twice the Heron area as signed area.

[formula: iem-degree-three-barycentric] -/
theorem degreeThree_layout (h1 : IsTriangle ab ra rb) (h2 : IsTriangle bc rb rc)
    (h3 : IsTriangle ca rc ra)
    (hflat : cornerAngle ab ra rb + cornerAngle bc rb rc + cornerAngle ca rc ra = 2 * π) :
    let P : ℝ × ℝ := (0, 0)
    let A := polar ra 0
    let B := layoutB ra rb ab
    let C := layoutC ra rb rc ab bc
    dist2 P A = ra ∧ dist2 P B = rb ∧ dist2 P C = rc ∧
      dist2 A B = ab ∧ dist2 B C = bc ∧ dist2 C A = ca ∧
      sda P A B = 2 * triArea ab ra rb ∧ sda P B C = 2 * triArea bc rb rc ∧
      sda P C A = 2 * triArea ca rc ra := by
  intro P A B C
  have hθ3 : cornerAngle ab ra rb + cornerAngle bc rb rc = 2 * π - cornerAngle ca rc ra := by
    linarith
  refine ⟨dist2_origin_polar h1.pos_b.le 0, dist2_origin_polar h1.pos_c.le _,
    dist2_origin_polar h2.pos_c.le _, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · simp only [A, B, layoutB, dist2_polar, sub_zero]
    exact law_of_cosines h1
  · simp only [B, C, layoutB, layoutC, dist2_polar, add_sub_cancel_left]
    exact law_of_cosines h2
  · simp only [C, A, layoutC, dist2_polar, zero_sub, cos_neg, hθ3, cos_two_pi_sub]
    have := law_of_cosines h3
    rwa [show rc ^ 2 + ra ^ 2 - 2 * rc * ra * cos (cornerAngle ca rc ra) =
      ra ^ 2 + rc ^ 2 - 2 * rc * ra * cos (cornerAngle ca rc ra) by ring] at this ⊢
  · simp only [P, A, B, layoutB, sda_origin_polar, sub_zero]
    exact mul_sin_cornerAngle h1
  · simp only [P, B, C, layoutB, layoutC, sda_origin_polar, add_sub_cancel_left]
    exact mul_sin_cornerAngle h2
  · simp only [P, C, A, layoutC, sda_origin_polar, zero_sub, sin_neg, hθ3, sin_two_pi_sub,
      neg_neg]
    exact mul_sin_cornerAngle h3

/-- **`iem-degree-three-barycentric`.** In a flat degree-3 star the merged face `abc` (laid out as
above) is counter-clockwise, its signed area is the sum of the three, and the barycentric
coordinates of `i` in it are `(A_ibc, A_ica, A_iab)/(A_ibc + A_ica + A_iab)` (Heron areas of the
three faces from their lengths), all positive.

[formula: iem-degree-three-barycentric] -/
theorem degreeThree_barycentric (h1 : IsTriangle ab ra rb) (h2 : IsTriangle bc rb rc)
    (h3 : IsTriangle ca rc ra)
    (hflat : cornerAngle ab ra rb + cornerAngle bc rb rc + cornerAngle ca rc ra = 2 * π) :
    let A := polar ra 0
    let B := layoutB ra rb ab
    let C := layoutC ra rb rc ab bc
    let S := triArea bc rb rc + triArea ca rc ra + triArea ab ra rb
    sda A B C = 2 * S ∧ 0 < S ∧
      Intrinsic.barycentric A.1 A.2 B.1 B.2 C.1 C.2 0 0 =
        (triArea bc rb rc / S, triArea ca rc ra / S, triArea ab ra rb / S) ∧
      0 < triArea bc rb rc / S ∧ 0 < triArea ca rc ra / S ∧ 0 < triArea ab ra rb / S := by
  intro A B C S
  obtain ⟨-, -, -, -, -, -, sAB, sBC, sCA⟩ := degreeThree_layout h1 h2 h3 hflat
  have t1 := triArea_pos h1
  have t2 := triArea_pos h2
  have t3 := triArea_pos h3
  have hsplit := Intrinsic.sda_split A B C ((0 : ℝ), (0 : ℝ))
  have e1 : sda A ((0 : ℝ), (0 : ℝ)) C = sda ((0 : ℝ), (0 : ℝ)) C A := by unfold sda; ring
  have e2 : sda A B ((0 : ℝ), (0 : ℝ)) = sda ((0 : ℝ), (0 : ℝ)) A B := by unfold sda; ring
  have hABC : sda A B C = 2 * S := by
    rw [← hsplit, e1, e2, sBC, sCA, sAB]; ring
  have hS : 0 < S := by positivity
  refine ⟨hABC, hS, ?_, by positivity, by positivity, by positivity⟩
  have := Intrinsic.barycentric_eq A B C ((0 : ℝ), (0 : ℝ))
  simp only at this
  rw [this, e1, e2, sBC, sCA, sAB, hABC]
  simp only [Prod.mk.injEq]
  refine ⟨?_, ?_, ?_⟩ <;> field_simp

end DegreeThree

/-- **The points of a removed face are re-expressed by linear substitution**: if `q` has the
coordinates `(β_p, β_a, β_b)` in a face `p a b` and `p` has the coordinates `(λ_a, λ_b, λ_c)` in
`a b c`, then `q` has `(β_a + β_p λ_a, β_b + β_p λ_b, β_p λ_c)` in `a b c`: non-negative and summing
to one when both are.

[formula: iem-degree-three-barycentric] -/
theorem subst_barycentric {E : Type*} [AddCommGroup E] [Module ℝ E] {p a b c q : E}
    {βp βa βb la lb lc : ℝ} (hq : q = βp • p + βa • a + βb • b) (hp : p = la • a + lb • b + lc • c) :
    q = (βa + βp * la) • a + (βb + βp * lb) • b + (βp * lc) • c ∧
      ((βp + βa + βb = 1 → la + lb + lc = 1 →
        (βa + βp * la) + (βb + βp * lb) + βp * lc = 1) ∧
      (0 ≤ βp → 0 ≤ βa → 0 ≤ βb → 0 ≤ la → 0 ≤ lb → 0 ≤ lc →
        0 ≤ βa + βp * la ∧ 0 ≤ βb + βp * lb ∧ 0 ≤ βp * lc)) := by
  refine ⟨?_, fun h1 h2 => ?_, fun h1 h2 h3 h4 h5 h6 => ⟨by positivity, by positivity, by positivity⟩⟩
  · rw [hq, hp]; module
  · linear_combination h1 + βp * h2

/-! ## `iem-boundary-merge` -/

section BoundaryMerge

variable {ra rb rc ab ca : ℝ}

/-- **`iem-boundary-merge`.** A flat boundary vertex `i` of degree 2: faces `iab` and `ica`
(spoke `ra = ℓ_ia`, boundary spokes `rb = ℓ_ib`, `rc = ℓ_ic`, outer sides `ab`, `ca`), strict
triangles with the angles at `i` summing to `π`. Laid out with `i = 0`, `a` on the `x`-axis, `b`
counter-clockwise and `c` clockwise of it, `i` lies on the segment `bc` at
`(rc·b + rb·c)/(rb + rc)`, the new boundary edge has length `ℓ_bc = rb + rc`, the merged face `abc`
is counter-clockwise with the sides `ab`, `rb + rc`, `ca`, and the barycentric coordinates of `i`
in it are `(0, rc/(rb + rc), rb/(rb + rc))`.

[formula: iem-boundary-merge] -/
theorem boundaryMerge (h1 : IsTriangle ab ra rb) (h2 : IsTriangle ca rc ra)
    (hflat : cornerAngle ab ra rb + cornerAngle ca rc ra = π) :
    let A := polar ra 0
    let B := polar rb (cornerAngle ab ra rb)
    let C := polar rc (-cornerAngle ca rc ra)
    dist2 A B = ab ∧ dist2 C A = ca ∧ dist2 B C = rb + rc ∧
      ((rc / (rb + rc)) * B.1 + (rb / (rb + rc)) * C.1,
        (rc / (rb + rc)) * B.2 + (rb / (rb + rc)) * C.2) = ((0 : ℝ), (0 : ℝ)) ∧
      0 < sda A B C ∧
      Intrinsic.barycentric A.1 A.2 B.1 B.2 C.1 C.2 0 0 = (0, rc / (rb + rc), rb / (rb + rc)) := by
  intro A B C
  have hb := h1.pos_c
  have hc := h2.pos_b
  set θ1 := cornerAngle ab ra rb
  set θ2 := cornerAngle ca rc ra
  have hθ : -θ2 = θ1 - π := by linarith
  have hC : C = (-(rc * cos θ1), -(rc * sin θ1)) := by
    simp only [C, polar]
    rw [show -cornerAngle ca rc ra = θ1 - π by linarith, cos_sub_pi, sin_sub_pi]
    ext <;> simp
  have hBC : dist2 B C = rb + rc := by
    rw [hC]
    unfold dist2
    simp only [B, polar]
    have : (-(rc * cos θ1) - rb * cos θ1) ^ 2 + (-(rc * sin θ1) - rb * sin θ1) ^ 2 =
        (rb + rc) ^ 2 := by linear_combination (rb + rc) ^ 2 * sin_sq_add_cos_sq θ1
    rw [this, sqrt_sq (by positivity)]
  have hmid : ((rc / (rb + rc)) * B.1 + (rb / (rb + rc)) * C.1,
      (rc / (rb + rc)) * B.2 + (rb / (rb + rc)) * C.2) = ((0 : ℝ), (0 : ℝ)) := by
    rw [hC]
    simp only [B, polar, Prod.mk.injEq]
    constructor <;> field_simp <;> ring
  have sAB : sda (0, 0) A B = 2 * triArea ab ra rb := by
    simp only [A, B, sda_origin_polar, sub_zero]; exact mul_sin_cornerAngle h1
  have sCA : sda (0, 0) C A = 2 * triArea ca rc ra := by
    simp only [C, A, sda_origin_polar, zero_sub, neg_neg]; exact mul_sin_cornerAngle h2
  have sBC : sda (0, 0) B C = 0 := by
    simp only [B, C, sda_origin_polar]
    rw [show -θ2 - θ1 = -π by linarith, sin_neg, sin_pi]; ring
  have hsplit := Intrinsic.sda_split A B C ((0 : ℝ), (0 : ℝ))
  have e1 : sda A ((0 : ℝ), (0 : ℝ)) C = sda ((0 : ℝ), (0 : ℝ)) C A := by unfold sda; ring
  have e2 : sda A B ((0 : ℝ), (0 : ℝ)) = sda ((0 : ℝ), (0 : ℝ)) A B := by unfold sda; ring
  have hABC : sda A B C = 2 * triArea ab ra rb + 2 * triArea ca rc ra := by
    rw [← hsplit, e1, e2, sBC, sCA, sAB]; ring
  have hpos : 0 < sda A B C := by
    rw [hABC]; have := triArea_pos h1; have := triArea_pos h2; positivity
  refine ⟨?_, ?_, hBC, hmid, hpos, ?_⟩
  · simp only [A, B, dist2_polar, sub_zero]; exact law_of_cosines h1
  · simp only [C, A, dist2_polar, zero_sub, neg_neg]
    have := law_of_cosines h2
    rwa [show rc ^ 2 + ra ^ 2 - 2 * rc * ra * cos θ2 = ra ^ 2 + rc ^ 2 - 2 * rc * ra * cos θ2 by
      ring] at this ⊢
  · have hsum : (0 : ℝ) + rc / (rb + rc) + rb / (rb + rc) = 1 := by field_simp; ring
    have hx := congrArg Prod.fst hmid
    have hy := congrArg Prod.snd hmid
    simp only at hx hy
    exact Intrinsic.barycentric_unique (p := ((0 : ℝ), (0 : ℝ))) hpos.ne' hsum
      (by rw [zero_mul, zero_add]; exact hx) (by rw [zero_mul, zero_add]; exact hy)

end BoundaryMerge

/-! ## `iem-flip-to-degree-three` -/

/-- **Interior, simplicial: Thm. D.1** (`ic-vertex-removal`): a flat star of degree `≥ 4` has a
spoke whose quadrilateral is convex (corners at most `π`). Not strictly, in general
(`DDG.Intrinsic.cross_star_no_strict`).

[formula: iem-flip-to-degree-three] -/
theorem flatStar_exists_flippable {n : ℕ} (S : DDG.Intrinsic.FlatStar n) (hn : 3 ≤ n) :
    ∃ k, S.inner k ≤ π ∧ S.outer k ≤ π :=
  S.exists_flippable hn

/-- A flat boundary star of degree `n + 3`: triangles `t_k = (i, j_k, j_{k+1})`, `k < n + 3`, with
the angles `θ_k` at `i`, `α_k` at `j_k`, `β_k` at `j_{k+1}`; `j_0` and `j_{n+3}` are the boundary
neighbours; the angles at `i` sum to `π`. The interior spokes are `i j_{s+1}`, `s < n + 2`. -/
structure FlatHalfStar (n : ℕ) where
  /-- angle at `i` of `t_k` -/
  θ : ℕ → ℝ
  /-- angle at `j_k` of `t_k` -/
  α : ℕ → ℝ
  /-- angle at `j_{k+1}` of `t_k` -/
  β : ℕ → ℝ
  θ_pos : ∀ k < n + 3, 0 < θ k
  α_pos : ∀ k < n + 3, 0 < α k
  β_pos : ∀ k < n + 3, 0 < β k
  triangle : ∀ k < n + 3, θ k + α k + β k = π
  flat : ∑ k ∈ range (n + 3), θ k = π

namespace FlatHalfStar

variable {n : ℕ} (S : FlatHalfStar n)

/-- The angle at `i` of the quadrilateral of the interior spoke `i j_{s+1}`. -/
def inner (s : ℕ) : ℝ := S.θ s + S.θ (s + 1)

/-- The angle at `j_{s+1}` of the quadrilateral of the interior spoke `i j_{s+1}`. -/
def outer (s : ℕ) : ℝ := S.β s + S.α (s + 1)

theorem inner_lt (s : ℕ) (hs : s < n + 2) : S.inner s < π := by
  have hmem1 : s ∈ range (n + 3) := mem_range.mpr (by omega)
  have hmem2 : s + 1 ∈ (range (n + 3)).erase s := by
    simp only [mem_erase, mem_range]; omega
  have hne : (((range (n + 3)).erase s).erase (s + 1)).Nonempty := by
    by_cases h : s = 0
    · exact ⟨2, by simp only [mem_erase, mem_range]; omega⟩
    · exact ⟨0, by simp only [mem_erase, mem_range]; omega⟩
  have hrest : 0 < ∑ k ∈ ((range (n + 3)).erase s).erase (s + 1), S.θ k :=
    Finset.sum_pos (fun k hk => S.θ_pos k (by simp only [mem_erase, mem_range] at hk; omega)) hne
  have h := S.flat
  rw [← add_sum_erase _ _ hmem1, ← add_sum_erase _ _ hmem2] at h
  unfold inner
  linarith

theorem sum_outer : ∑ s ∈ range (n + 2), S.outer s =
    ((n + 2 : ℕ) : ℝ) * π - S.β (n + 2) - S.α 0 := by
  have htri : ∑ k ∈ range (n + 3), (S.α k + S.β k) = ∑ k ∈ range (n + 3), (π - S.θ k) :=
    Finset.sum_congr rfl (fun k hk => by have := S.triangle k (mem_range.mp hk); linarith)
  rw [Finset.sum_sub_distrib, S.flat, Finset.sum_const, card_range, nsmul_eq_mul,
    Finset.sum_add_distrib, Finset.sum_range_succ' S.α, Finset.sum_range_succ S.β] at htri
  unfold outer
  rw [Finset.sum_add_distrib]
  push_cast at htri ⊢
  linarith

/-- **`iem-flip-to-degree-three`, boundary analogue.** A flat boundary vertex of degree `n + 3 ≥ 3`
has an interior spoke `i j_{s+1}` whose quadrilateral is strictly convex: inner corner and outer
corner both below `π` (the other two corners are triangle angles).

[formula: iem-flip-to-degree-three] -/
theorem exists_flippable : ∃ s < n + 2, S.inner s < π ∧ S.outer s < π := by
  by_contra h
  push Not at h
  have hall : ∀ s ∈ range (n + 2), π ≤ S.outer s := fun s hs =>
    h s (mem_range.mp hs) (S.inner_lt s (mem_range.mp hs))
  have hge : ∑ s ∈ range (n + 2), π ≤ ∑ s ∈ range (n + 2), S.outer s := Finset.sum_le_sum hall
  rw [S.sum_outer, Finset.sum_const, card_range, nsmul_eq_mul] at hge
  have := S.β_pos (n + 2) (by omega)
  have := S.α_pos 0 (by omega)
  linarith

end FlatHalfStar

end IntrinsicSimplification

end DDG
