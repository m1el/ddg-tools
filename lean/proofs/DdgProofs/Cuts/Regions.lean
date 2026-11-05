import DdgProofs.Cuts.OnCut
import DdgProofs.Cuts.Crossing

/-!
# Cut flow: `N` regions

`crates/ddg-variational-cuts/SPEC.md` §8 (paper App. A), D20, D21. A state is one function `φ_k`
per region; the region of a vertex is the function that is smallest there (ties: the lowest
index). The theorems are about the oracle's `crossingParamRegions`, `triplePointCoord`,
`tripleRate`, `projectionShift`, `projectedValue`, `turningAngle`, `lengthForce` and, for the
linear systems, about `Cuts.fluxMatrix`. `IsRegion`, `IsTwoSmallest`, `regionRate` are this
package's own definitions of what SPEC says in words.

* **`cut-edge-crossing-regions`** — with `i = region(a) ≠ j = region(b)`: the differences
  `ψ_a = φ_i(a) − φ_j(a) ≤ 0 ≤ ψ_b` cannot both vanish (`region_differences`), `t` is where the
  interpolants of `φ_i` and `φ_j` are equal (`crossingParamRegions_interp`), `0 ≤ t ≤ 1`, with
  `t < 1` if `i < j` and `0 < t` if `j < i` (`crossingParamRegions_mem`); on a state `[φ, −φ]`
  the regions and the parameter are those of the two-region code (`isRegion_two_zero_iff`,
  `isRegion_two_one_iff`, `crossingParamRegions_two`).
* **`cut-triple-point`** — the barycentric coordinates `triplePointCoord` of the centroid of the
  three cut points are positive and sum to `1` (`triplePoint_coord_sum`, `triplePoint_coord_pos`,
  `triplePoint_eq_centroid`); the six triangles `(c_k, p_k, m)`, `(c_k, m, p_{k−1})` have
  positive determinants that sum to `1` (`tripleFace_dets`, `tripleFace_dets_pos`), so their areas
  `|det| · A_f` sum to `A_f`; each corner quad `(c_k, p_k, m, p_{k−1})` is convex: also the two
  triangles of its other diagonal are positively oriented, the centroid lying beyond the line
  `p_{k−1} p_k` (`tripleFace_convex`).
* **`cut-region-projection`** — after subtracting the mean of the two smallest entries the
  smallest is `−d`, the second `+d`, every other entry is `≥ d`, exactly one entry is negative
  when the two smallest differ, all differences, the regions and the crossings are kept, and the
  projection is idempotent (`projected_*`); on `[φ, −φ]` nothing changes
  (`projectionShift_two`).
* **`cut-region-extension`** — the rate at a triple point is the oracle's `smoothedSpeed` of the
  two neighbouring cut points: a convex combination, exact for a rate linear in arclength
  (`tripleRate_mem`, `tripleRate_linear`); with the substitution `φ̇_j := −φ̇_k` a state
  `[φ, −φ]` gets the rates `[σ, −σ]` of the two-region code and keeps its form under the explicit
  step, the length step and the symmetry projection (`regionRate_two`, `explicitStep_neg`,
  `lengthStep_neg`, `classMean_neg`, `harmonicExtension_neg`).
* **`cut-length-step-regions`** — the systems of §8.5 are the on-cut systems on the boundary of
  a region, so the M-matrix facts, uniqueness, `Σ ℓ λ⁰ = Σ Θ` on the boundary of every region,
  and the bounds hold (`regionChain_*`); the two sides of a curve get opposite forces and
  opposite rates (`curve_loads_opposite`, `onCut_rates_opposite`); the loads `γ − β` of the
  three regions at a triple point sum to zero when their angles sum to `2π`, and each vanishes
  at a balanced junction (`junction_loads_sum`, `junction_load_balanced`); the explicit step
  changes `φ_k − φ_j` at the rate `δ` of the interface and reduces to the two-region step
  (`regionStep_difference`, `regionStep_two`).
-/

noncomputable section

namespace DDG

namespace Cuts

open Matrix
open scoped DDG.Oracle
open Real

/-! ## The region of a vertex (SPEC §8.1) -/

section Region

variable {K : Type*} [LinearOrder K]

/-- `k` is the region of a vertex with the values `x_l = φ_l(vertex)`: `x_k` is the smallest,
and among equal values the index is the lowest. -/
def IsRegion (x : K → ℝ) (k : K) : Prop := (∀ l, x k ≤ x l) ∧ ∀ l, l < k → x k < x l

/-- The region of a vertex is unique. -/
theorem IsRegion.unique {x : K → ℝ} {k k' : K} (h : IsRegion x k) (h' : IsRegion x k') :
    k = k' := by
  rcases lt_trichotomy k k' with hlt | heq | hgt
  · exact absurd (h'.2 k hlt) (not_lt.mpr (h.1 k'))
  · exact heq
  · exact absurd (h.2 k' hgt) (not_lt.mpr (h'.1 k))

/-- Every vertex has a region. -/
theorem exists_isRegion [Fintype K] [Nonempty K] (x : K → ℝ) : ∃ k, IsRegion x k := by
  classical
  obtain ⟨k₀, -, hk₀⟩ := Finset.exists_min_image Finset.univ x Finset.univ_nonempty
  have hS : (Finset.univ.filter fun k => ∀ l, x k ≤ x l).Nonempty :=
    ⟨k₀, by simpa using fun l => hk₀ l (Finset.mem_univ l)⟩
  refine ⟨(Finset.univ.filter fun k => ∀ l, x k ≤ x l).min' hS, ?_, ?_⟩
  · have := Finset.min'_mem _ hS
    simpa using this
  · intro l hl
    have hmem := Finset.min'_mem _ hS
    have hmin : ∀ l', x ((Finset.univ.filter fun k => ∀ l, x k ≤ x l).min' hS) ≤ x l' := by
      simpa using hmem
    rcases (hmin l).lt_or_eq with h | h
    · exact h
    · exfalso
      have hl' : l ∈ Finset.univ.filter fun k => ∀ l, x k ≤ x l := by
        simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        intro l'
        rw [← h]
        exact hmin l'
      exact absurd (Finset.min'_le _ l hl') (not_le.mpr hl)

/-- Subtracting the same number from all values at a vertex does not change its region. -/
theorem isRegion_sub_iff (x : K → ℝ) (m : ℝ) (k : K) :
    IsRegion (fun l => x l - m) k ↔ IsRegion x k := by
  simp only [IsRegion, sub_le_sub_iff_right, sub_lt_sub_iff_right]

/-- **Two regions**: on a state `[φ, −φ]` a vertex is in region 0 iff `φ ≤ 0` … -/
theorem isRegion_two_zero_iff (φ : ℝ) : IsRegion ![φ, -φ] (0 : Fin 2) ↔ φ ≤ 0 := by
  constructor
  · intro h
    have h1 : φ ≤ -φ := h.1 1
    linarith
  · intro h
    refine ⟨fun l => ?_, fun l hl => absurd hl (Fin.not_lt_zero l)⟩
    fin_cases l
    · exact le_rfl
    · show φ ≤ -φ
      linarith

/-- … and in region 1 iff `φ > 0`: the rule of SPEC §3.2. -/
theorem isRegion_two_one_iff (φ : ℝ) : IsRegion ![φ, -φ] (1 : Fin 2) ↔ 0 < φ := by
  constructor
  · intro h
    have h1 : -φ < φ := h.2 0 (by decide)
    linarith
  · intro h
    refine ⟨fun l => ?_, fun l hl => ?_⟩
    · fin_cases l
      · show -φ ≤ φ
        linarith
      · exact le_rfl
    · fin_cases l
      · show -φ < φ
        linarith
      · exact absurd hl (by decide)

/-! ## Crossings between two regions (SPEC §8.2, paper Eq. 31) -/

/-- The oracle's `crossingParamRegions` is the level-crossing parameter of the two
differences. -/
theorem crossingParamRegions_eq (a b c d : ℝ) :
    Ddg.crossingParamRegions a b c d = levelCrossingParam (a - b) (c - d) := rfl

/-- The crossing parameter does not change when both values are negated. -/
theorem crossingParam_neg (a b : ℝ) : levelCrossingParam (-a) (-b) = levelCrossingParam a b := by
  rw [crossingParam_eq, crossingParam_eq]
  have : -a - -b = -(a - b) := by ring
  rw [this, neg_div_neg_eq]

/-- **The two differences.** If `i` is the region of `a` and `j ≠ i` the region of `b`, then
`ψ_a = φ_i(a) − φ_j(a) ≤ 0 ≤ ψ_b = φ_i(b) − φ_j(b)`, with `ψ_b > 0` if `i < j` and `ψ_a < 0` if
`j < i` (the tie rule): they cannot both vanish.

[formula: cut-edge-crossing-regions] -/
theorem region_differences {xa xb : K → ℝ} {i j : K} (hi : IsRegion xa i) (hj : IsRegion xb j)
    (hij : i ≠ j) :
    xa i - xa j ≤ 0 ∧ 0 ≤ xb i - xb j ∧ (i < j → 0 < xb i - xb j) ∧ (j < i → xa i - xa j < 0) ∧
      xa i - xa j ≠ xb i - xb j := by
  have h1 : xa i - xa j ≤ 0 := by linarith [hi.1 j]
  have h2 : 0 ≤ xb i - xb j := by linarith [hj.1 i]
  have h3 : i < j → 0 < xb i - xb j := fun h => by linarith [hj.2 i h]
  have h4 : j < i → xa i - xa j < 0 := fun h => by linarith [hi.2 j h]
  refine ⟨h1, h2, h3, h4, ?_⟩
  rcases lt_or_gt_of_ne hij with h | h
  · have := h3 h
    intro heq
    linarith
  · have := h4 h
    intro heq
    linarith

/-- **The cut point is where the interpolants of `φ_i` and `φ_j` are equal** (paper Eq. 31).

[formula: cut-edge-crossing-regions] -/
theorem crossingParamRegions_interp {xa xb : K → ℝ} {i j : K} (hi : IsRegion xa i)
    (hj : IsRegion xb j) (hij : i ≠ j) :
    (1 - Ddg.crossingParamRegions (xa i) (xa j) (xb i) (xb j)) * xa i +
        Ddg.crossingParamRegions (xa i) (xa j) (xb i) (xb j) * xb i =
      (1 - Ddg.crossingParamRegions (xa i) (xa j) (xb i) (xb j)) * xa j +
        Ddg.crossingParamRegions (xa i) (xa j) (xb i) (xb j) * xb j := by
  have hne := (region_differences hi hj hij).2.2.2.2
  have h := interp_crossingParam hne
  rw [crossingParamRegions_eq]
  linear_combination h

/-- **The parameter lies in `[0, 1]`**: in `[0, 1)` if `i < j` and in `(0, 1]` if `j < i`.

[formula: cut-edge-crossing-regions] -/
theorem crossingParamRegions_mem {xa xb : K → ℝ} {i j : K} (hi : IsRegion xa i)
    (hj : IsRegion xb j) (hij : i ≠ j) :
    0 ≤ Ddg.crossingParamRegions (xa i) (xa j) (xb i) (xb j) ∧
      Ddg.crossingParamRegions (xa i) (xa j) (xb i) (xb j) ≤ 1 ∧
      (i < j → Ddg.crossingParamRegions (xa i) (xa j) (xb i) (xb j) < 1) ∧
      (j < i → 0 < Ddg.crossingParamRegions (xa i) (xa j) (xb i) (xb j)) := by
  obtain ⟨h1, h2, h3, h4, -⟩ := region_differences hi hj hij
  rw [crossingParamRegions_eq]
  rcases lt_or_gt_of_ne hij with h | h
  · obtain ⟨ha, hb⟩ := crossingParam_mem_Ico h1 (h3 h)
    exact ⟨ha, hb.le, fun _ => hb, fun h' => absurd h (not_lt.mpr h'.le)⟩
  · have hm := crossingParam_mem_Ioc (neg_pos.mpr (h4 h)) (neg_nonpos.mpr h2)
    rw [crossingParam_neg] at hm
    exact ⟨hm.1.le, hm.2, fun h' => absurd h (not_lt.mpr h'.le), fun _ => hm.1⟩

end Region

/-- **On a state `[φ, −φ]` the crossing is the two-region crossing** `φ_a/(φ_a − φ_b)`, whichever
end of the edge is in region 0.

[formula: cut-edge-crossing-regions] -/
theorem crossingParamRegions_two (φa φb : ℝ) :
    Ddg.crossingParamRegions φa (-φa) φb (-φb) = Ddg.crossingParam φa φb ∧
      Ddg.crossingParamRegions (-φa) φa (-φb) φb = Ddg.crossingParam φa φb := by
  constructor
  · rw [crossingParamRegions_eq, Oracle.crossingParam_eq, crossingParam_eq, crossingParam_eq]
    have e1 : φa - -φa = 2 * φa := by ring
    have e2 : φa - -φa - (φb - -φb) = 2 * (φa - φb) := by ring
    rw [e2, e1, mul_div_mul_left _ _ (two_ne_zero)]
  · rw [crossingParamRegions_eq, Oracle.crossingParam_eq, crossingParam_eq, crossingParam_eq]
    have e1 : -φa - φa = -2 * φa := by ring
    have e2 : -φa - φa - (-φb - φb) = -2 * (φa - φb) := by ring
    rw [e2, e1, mul_div_mul_left _ _ (by norm_num : (-2 : ℝ) ≠ 0)]

/-! ## The triple point of a face with three regions (SPEC §8.2) -/

/-- The oracle's `triplePointCoord` at `ℝ`. -/
theorem triplePointCoord_eq (sNext sPrev : ℝ) :
    Ddg.triplePointCoord sNext sPrev = ((1 - sNext) + sPrev) / 3 := by
  simp only [Ddg.triplePointCoord, Oracle.one_eq, Oracle.three_eq]

/-- **The barycentric coordinates of the triple point sum to `1`.**

[formula: cut-triple-point] -/
theorem triplePoint_coord_sum (s0 s1 s2 : ℝ) :
    Ddg.triplePointCoord s0 s2 + Ddg.triplePointCoord s1 s0 + Ddg.triplePointCoord s2 s1 = 1 := by
  simp only [triplePointCoord_eq]
  ring

/-- **They are positive**, so the triple point lies inside the face (for cut points that are not
at the far end of their edges).

[formula: cut-triple-point] -/
theorem triplePoint_coord_pos {sNext sPrev : ℝ} (h1 : sNext < 1) (h2 : 0 < sPrev) :
    0 < Ddg.triplePointCoord sNext sPrev := by
  rw [triplePointCoord_eq]
  have : 0 < (1 - sNext) + sPrev := by linarith
  positivity

/-- **The triple point is the centroid of the three cut points**: with
`p_k = (1 − s_k) c_k + s_k c_{k+1}`, `(p_0 + p_1 + p_2)/3` has the barycentric coordinates
`triplePointCoord` with respect to the corners. -/
theorem triplePoint_eq_centroid {E : Type*} [AddCommGroup E] [Module ℝ E] (c0 c1 c2 : E)
    (s0 s1 s2 : ℝ) :
    (1 / 3 : ℝ) • (((1 - s0) • c0 + s0 • c1) + ((1 - s1) • c1 + s1 • c2) +
        ((1 - s2) • c2 + s2 • c0)) =
      Ddg.triplePointCoord s0 s2 • c0 + Ddg.triplePointCoord s1 s0 • c1 +
        Ddg.triplePointCoord s2 s1 • c2 := by
  simp only [triplePointCoord_eq]
  module

/-- **The six triangles of a face with three regions.** Barycentric rows: corners
`(1,0,0), (0,1,0), (0,0,1)`, cut points `p_0 = (1−s_0, s_0, 0)`, `p_1 = (0, 1−s_1, s_1)`,
`p_2 = (s_2, 0, 1−s_2)`, triple point `m = (m_0, m_1, m_2)`. The determinants of
`(c_k, p_k, m)` and `(c_k, m, p_{k−1})` are `s_k m_{k+2}` and `m_{k+1} (1 − s_{k−1})`, and the six
sum to `m_0 + m_1 + m_2 = 1`: with `cut-subtriangle-area` the areas sum to the area of the face.

[formula: cut-triple-point] -/
theorem tripleFace_dets (s0 s1 s2 : ℝ) {m0 m1 m2 : ℝ} (hm : m0 + m1 + m2 = 1) :
    Ddg.det3 1 0 0 (1 - s0) s0 0 m0 m1 m2 = s0 * m2 ∧
    Ddg.det3 1 0 0 m0 m1 m2 s2 0 (1 - s2) = m1 * (1 - s2) ∧
    Ddg.det3 0 1 0 0 (1 - s1) s1 m0 m1 m2 = s1 * m0 ∧
    Ddg.det3 0 1 0 m0 m1 m2 (1 - s0) s0 0 = m2 * (1 - s0) ∧
    Ddg.det3 0 0 1 s2 0 (1 - s2) m0 m1 m2 = s2 * m1 ∧
    Ddg.det3 0 0 1 m0 m1 m2 0 (1 - s1) s1 = m0 * (1 - s1) ∧
    s0 * m2 + m1 * (1 - s2) + s1 * m0 + m2 * (1 - s0) + s2 * m1 + m0 * (1 - s1) = 1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · unfold Ddg.det3; ring
  · unfold Ddg.det3; ring
  · unfold Ddg.det3; ring
  · unfold Ddg.det3; ring
  · unfold Ddg.det3; ring
  · unfold Ddg.det3; ring
  · linear_combination hm

/-- **The six triangles are positively oriented** for cut points strictly inside their edges and
the triple point at their centroid.

[formula: cut-triple-point] -/
theorem tripleFace_dets_pos {s0 s1 s2 : ℝ} (h0 : 0 < s0) (h0' : s0 < 1) (h1 : 0 < s1)
    (h1' : s1 < 1) (h2 : 0 < s2) (h2' : s2 < 1) :
    0 < s0 * Ddg.triplePointCoord s2 s1 ∧
    0 < Ddg.triplePointCoord s1 s0 * (1 - s2) ∧
    0 < s1 * Ddg.triplePointCoord s0 s2 ∧
    0 < Ddg.triplePointCoord s2 s1 * (1 - s0) ∧
    0 < s2 * Ddg.triplePointCoord s1 s0 ∧
    0 < Ddg.triplePointCoord s0 s2 * (1 - s1) := by
  have m0 := triplePoint_coord_pos h0' h2
  have m1 := triplePoint_coord_pos h1' h0
  have m2 := triplePoint_coord_pos h2' h1
  exact ⟨mul_pos h0 m2, mul_pos m1 (by linarith), mul_pos h1 m0, mul_pos m2 (by linarith),
    mul_pos h2 m1, mul_pos m0 (by linarith)⟩

/-- **The corner quads are convex.** For the quad `(c_k, p_k, m, p_{k−1})` the triangles of the
other diagonal, `(c_k, p_k, p_{k−1})` and `(p_k, m, p_{k−1})`, are positively oriented too: their
determinants are `s_k (1 − s_{k−1})` and `((1−s_0)(1−s_1)(1−s_2) + s_0 s_1 s_2)/3` (a third of
the triangle of the three cut points, whose determinant is `(1−s_0)(1−s_1)(1−s_2) + s_0 s_1 s_2`),
so the centroid lies beyond the line through the two cut points of the corner.

[formula: cut-triple-point] -/
theorem tripleFace_convex {s0 s1 s2 : ℝ} (h0 : 0 < s0) (h0' : s0 < 1) (h1 : 0 < s1)
    (h1' : s1 < 1) (h2 : 0 < s2) (h2' : s2 < 1) :
    Ddg.det3 (1 - s0) s0 0 0 (1 - s1) s1 s2 0 (1 - s2) =
        (1 - s0) * (1 - s1) * (1 - s2) + s0 * s1 * s2 ∧
    0 < (1 - s0) * (1 - s1) * (1 - s2) + s0 * s1 * s2 ∧
    Ddg.det3 1 0 0 (1 - s0) s0 0 s2 0 (1 - s2) = s0 * (1 - s2) ∧
    Ddg.det3 (1 - s0) s0 0 (Ddg.triplePointCoord s0 s2) (Ddg.triplePointCoord s1 s0)
        (Ddg.triplePointCoord s2 s1) s2 0 (1 - s2) =
      ((1 - s0) * (1 - s1) * (1 - s2) + s0 * s1 * s2) / 3 ∧
    Ddg.det3 0 1 0 0 (1 - s1) s1 (1 - s0) s0 0 = s1 * (1 - s0) ∧
    Ddg.det3 0 (1 - s1) s1 (Ddg.triplePointCoord s0 s2) (Ddg.triplePointCoord s1 s0)
        (Ddg.triplePointCoord s2 s1) (1 - s0) s0 0 =
      ((1 - s0) * (1 - s1) * (1 - s2) + s0 * s1 * s2) / 3 ∧
    Ddg.det3 0 0 1 s2 0 (1 - s2) 0 (1 - s1) s1 = s2 * (1 - s1) ∧
    Ddg.det3 s2 0 (1 - s2) (Ddg.triplePointCoord s0 s2) (Ddg.triplePointCoord s1 s0)
        (Ddg.triplePointCoord s2 s1) 0 (1 - s1) s1 =
      ((1 - s0) * (1 - s1) * (1 - s2) + s0 * s1 * s2) / 3 := by
  have hpos : 0 < (1 - s0) * (1 - s1) * (1 - s2) + s0 * s1 * s2 := by
    have a0 : 0 < 1 - s0 := by linarith
    have a1 : 0 < 1 - s1 := by linarith
    have a2 : 0 < 1 - s2 := by linarith
    positivity
  refine ⟨?_, hpos, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · unfold Ddg.det3; ring
  · unfold Ddg.det3; ring
  · simp only [triplePointCoord_eq]; unfold Ddg.det3; ring
  · unfold Ddg.det3; ring
  · simp only [triplePointCoord_eq]; unfold Ddg.det3; ring
  · unfold Ddg.det3; ring
  · simp only [triplePointCoord_eq]; unfold Ddg.det3; ring

/-! ## The projection of the redistancing (SPEC §8.4, Losasso et al. 2006) -/

section Projection

variable {K : Type*}

/-- `i0` and `i1` are the positions of the smallest and the second smallest of the values `x`
at a vertex. -/
structure IsTwoSmallest (x : K → ℝ) (i0 i1 : K) : Prop where
  ne : i0 ≠ i1
  le : x i0 ≤ x i1
  rest : ∀ k, k ≠ i0 → x i1 ≤ x k

/-- With at least two functions the two smallest exist. -/
theorem exists_isTwoSmallest [Fintype K] (h : ∃ a b : K, a ≠ b) (x : K → ℝ) :
    ∃ i0 i1, IsTwoSmallest x i0 i1 := by
  classical
  obtain ⟨a, b, hab⟩ := h
  obtain ⟨i0, -, h0⟩ := Finset.exists_min_image Finset.univ x ⟨a, Finset.mem_univ a⟩
  have hne : (Finset.univ.erase i0).Nonempty := by
    by_cases ha : a = i0
    · exact ⟨b, Finset.mem_erase.mpr ⟨fun hb => hab (ha.trans hb.symm), Finset.mem_univ b⟩⟩
    · exact ⟨a, Finset.mem_erase.mpr ⟨ha, Finset.mem_univ a⟩⟩
  obtain ⟨i1, hi1, h1⟩ := Finset.exists_min_image (Finset.univ.erase i0) x hne
  exact ⟨i0, i1, (Finset.ne_of_mem_erase hi1).symm, h0 i1 (Finset.mem_univ i1),
    fun k hk => h1 k (Finset.mem_erase.mpr ⟨hk, Finset.mem_univ k⟩)⟩

/-- The oracle's `projectionShift`, `projectedValue` at `ℝ`. -/
theorem projectionShift_eq (a b : ℝ) : Ddg.projectionShift a b = (a + b) / 2 := by
  simp only [Ddg.projectionShift, Oracle.two_eq]

theorem projectedValue_eq (x m : ℝ) : Ddg.projectedValue x m = x - m := rfl

/-- **After the projection the smallest entry is `−d` and the second `+d`**,
`d = (second − smallest)/2 ≥ 0`.

[formula: cut-region-projection] -/
theorem projected_two_smallest {x : K → ℝ} {i0 i1 : K} (h : IsTwoSmallest x i0 i1) :
    Ddg.projectedValue (x i0) (Ddg.projectionShift (x i0) (x i1)) = -((x i1 - x i0) / 2) ∧
      Ddg.projectedValue (x i1) (Ddg.projectionShift (x i0) (x i1)) = (x i1 - x i0) / 2 ∧
      0 ≤ (x i1 - x i0) / 2 := by
  rw [projectedValue_eq, projectedValue_eq, projectionShift_eq]
  refine ⟨by ring, by ring, ?_⟩
  have := h.le
  linarith

/-- **Every other entry is at least the second.**

[formula: cut-region-projection] -/
theorem projected_rest {x : K → ℝ} {i0 i1 : K} (h : IsTwoSmallest x i0 i1) {k : K}
    (hk : k ≠ i0) :
    (x i1 - x i0) / 2 ≤ Ddg.projectedValue (x k) (Ddg.projectionShift (x i0) (x i1)) := by
  rw [projectedValue_eq, projectionShift_eq]
  have := h.rest k hk
  linarith

/-- **Exactly one entry is negative when the two smallest differ.**

[formula: cut-region-projection] -/
theorem projected_neg_iff {x : K → ℝ} {i0 i1 : K} (h : IsTwoSmallest x i0 i1)
    (hlt : x i0 < x i1) (k : K) :
    Ddg.projectedValue (x k) (Ddg.projectionShift (x i0) (x i1)) < 0 ↔ k = i0 := by
  constructor
  · intro hneg
    by_contra hk
    have := projected_rest h hk
    linarith
  · rintro rfl
    rw [(projected_two_smallest h).1]
    linarith

/-- **All differences are kept**, at the vertex … -/
theorem projected_sub (x y m : ℝ) :
    Ddg.projectedValue x m - Ddg.projectedValue y m = x - y := by
  simp only [projectedValue_eq]
  ring

/-- … hence the region of the vertex is kept … -/
theorem isRegion_projected [LinearOrder K] (x : K → ℝ) (m : ℝ) (k : K) :
    IsRegion (fun l => Ddg.projectedValue (x l) m) k ↔ IsRegion x k :=
  isRegion_sub_iff x m k

/-- … **and no crossing moves**: the crossing parameter of an edge uses only the differences
`φ_i − φ_j` at its two ends, each of which is shifted as a whole.

[formula: cut-region-projection] -/
theorem crossingParamRegions_projected (a b c d ma mb : ℝ) :
    Ddg.crossingParamRegions (Ddg.projectedValue a ma) (Ddg.projectedValue b ma)
        (Ddg.projectedValue c mb) (Ddg.projectedValue d mb) =
      Ddg.crossingParamRegions a b c d := by
  simp only [crossingParamRegions_eq, projectedValue_eq]
  congr 1 <;> ring

/-- **The projection is idempotent**: the projected values have the same two smallest, and their
mean is `0`. -/
theorem projected_isTwoSmallest {x : K → ℝ} {i0 i1 : K} (h : IsTwoSmallest x i0 i1) (m : ℝ) :
    IsTwoSmallest (fun k => Ddg.projectedValue (x k) m) i0 i1 :=
  ⟨h.ne, by simp only [projectedValue_eq]; linarith [h.le],
    fun k hk => by simp only [projectedValue_eq]; linarith [h.rest k hk]⟩

theorem projectionShift_projected (a b : ℝ) :
    Ddg.projectionShift (Ddg.projectedValue a (Ddg.projectionShift a b))
      (Ddg.projectedValue b (Ddg.projectionShift a b)) = 0 := by
  simp only [projectionShift_eq, projectedValue_eq]
  ring

/-- **For two regions nothing changes**: on `[φ, −φ]` the mean of the two (smallest) entries is
`0`.

[formula: cut-region-projection] -/
theorem projectionShift_two (φ : ℝ) :
    Ddg.projectionShift φ (-φ) = 0 ∧ Ddg.projectionShift (-φ) φ = 0 ∧
      Ddg.projectedValue φ 0 = φ := by
  simp only [projectionShift_eq, projectedValue_eq]
  refine ⟨by ring, by ring, by ring⟩

end Projection

/-! ## Rates of the functions (SPEC §8.3, §8.4) -/

/-- The oracle's `tripleRate` is the arclength interpolation `smoothedSpeed` between the two
cut points next to the triple point on the boundary of the region. -/
theorem tripleRate_eq_smoothedSpeed (l l' s s' : ℝ) :
    Ddg.tripleRate l l' s s' = Ddg.smoothedSpeed l l' s s' := rfl

/-- **The rate at a triple point lies between the rates of the two cut points next to it.**

[formula: cut-region-extension] -/
theorem tripleRate_mem {l l' : ℝ} (hl : 0 ≤ l) (hl' : 0 ≤ l') (h : 0 < l + l') (s s' : ℝ) :
    min s s' ≤ Ddg.tripleRate l l' s s' ∧ Ddg.tripleRate l l' s s' ≤ max s s' := by
  rw [tripleRate_eq_smoothedSpeed]
  exact smoothedSpeed_mem hl hl' h s s'

/-- **It is exact for a rate that is linear in arclength** along the boundary of the region
through the triple point.

[formula: cut-region-extension] -/
theorem tripleRate_linear {l l' : ℝ} (h : l + l' ≠ 0) (α β : ℝ) :
    Ddg.tripleRate l l' (α - β * l) (α + β * l') = α := by
  rw [tripleRate_eq_smoothedSpeed]
  exact smoothedSpeed_linear h α β

/-- The rate at a triple point is odd in the rates (it is linear). -/
theorem tripleRate_neg (l l' s s' : ℝ) :
    Ddg.tripleRate l l' (-s) (-s') = -Ddg.tripleRate l l' s s' := by
  simp only [tripleRate_eq_smoothedSpeed, smoothedSpeed_eq]
  ring

/-- The rates of all functions at a vertex of the region `own`, from the extension `e` of the
rates of `φ_own` (SPEC §8.4): `φ̇_own = e` and `φ̇_j = −e` for every other `j`. -/
def regionRate {K : Type*} [DecidableEq K] (own : K) (e : ℝ) (j : K) : ℝ :=
  if j = own then e else -e

theorem regionRate_own {K : Type*} [DecidableEq K] (own : K) (e : ℝ) :
    regionRate own e own = e := by
  simp [regionRate]

theorem regionRate_other {K : Type*} [DecidableEq K] {own j : K} (h : j ≠ own) (e : ℝ) :
    regionRate own e j = -e := by
  simp [regionRate, h]

/-- The difference `φ_own − φ_j` changes at the rate `2e`; the functions of the other regions
move together (their differences are kept). -/
theorem regionRate_sub {K : Type*} [DecidableEq K] {own j j' : K} (h : j ≠ own) (h' : j' ≠ own)
    (e : ℝ) :
    regionRate own e own - regionRate own e j = 2 * e ∧
      regionRate own e j - regionRate own e j' = 0 := by
  rw [regionRate_own, regionRate_other h, regionRate_other h']
  constructor <;> ring

/-- **Two regions: the rates are `[σ, −σ]` at every vertex.** With `σ` the two-region extension
of `s`: at a vertex of region 0 the extension of the rates of `φ_0` is `e = σ`; at a vertex of
region 1 the boundary rates of `φ_1` are `−s`, their extension is `e = −σ`
(`harmonicExtension_neg`); in both cases `φ̇_0 = σ`, `φ̇_1 = −σ`.

[formula: cut-region-extension] -/
theorem regionRate_two (σ : ℝ) :
    regionRate (0 : Fin 2) σ = ![σ, -σ] ∧ regionRate (1 : Fin 2) (-σ) = ![σ, -σ] := by
  constructor <;> funext j <;> fin_cases j <;> simp [regionRate]

/-- The harmonic extension is odd in its data: if `σ` extends `s` then `−σ` extends `−s`. -/
theorem harmonicExtension_neg {V : Type*} [Fintype V] {L : Matrix V V ℝ} {I : V → Prop}
    {s σ : V → ℝ} (h : (∀ v, ¬I v → σ v = s v) ∧ ∀ v, I v → (L *ᵥ σ) v = 0) :
    (∀ v, ¬I v → (-σ) v = (-s) v) ∧ ∀ v, I v → (L *ᵥ (-σ)) v = 0 := by
  refine ⟨fun v hv => by simp [h.1 v hv], fun v hv => ?_⟩
  rw [Matrix.mulVec_neg]
  simp [h.2 v hv]

/-- **The explicit step keeps the form `[φ, −φ]`**: `(−φ) + τ(−σ) = −(φ + τσ)`.

[formula: cut-region-extension] -/
theorem explicitStep_neg (φ τ σ : ℝ) :
    Ddg.explicitStep (-φ) τ (-σ) = -Ddg.explicitStep φ τ σ := by
  simp only [explicitStep_eq]
  ring

/-- **The length step keeps the form `[φ, −φ]`**: if `φ′` solves `(W_V + τL_w) φ′ = W_V φ̃` then
`−φ′` solves it for `−φ̃` (and it is the only solution, `lengthStep_existsUnique`). -/
theorem lengthStep_neg {V : Type*} [Fintype V] [DecidableEq V] {A : V → ℝ} {L : Matrix V V ℝ}
    {τ : ℝ} {φt φ' : V → ℝ}
    (h : (Matrix.diagonal A + τ • L) *ᵥ φ' = Matrix.diagonal A *ᵥ φt) :
    (Matrix.diagonal A + τ • L) *ᵥ (-φ') = Matrix.diagonal A *ᵥ (-φt) := by
  rw [Matrix.mulVec_neg, Matrix.mulVec_neg, h]

/-- **The symmetry projection keeps the form `[φ, −φ]`**: the class mean is odd. -/
theorem classMean_neg {V κ : Type*} [Fintype V] [DecidableEq κ] (cls : V → κ) (φ : V → ℝ) :
    classMean cls (fun v => -φ v) = fun v => -classMean cls φ v := by
  funext v
  unfold classMean
  rw [Finset.sum_neg_distrib, neg_div]

/-! ## The default stages on networks (SPEC §8.5, D20) -/

/-- **A curve between two regions: opposite turning.** If the interior angles of the two
regions at a cut point sum to `2π`, their turning angles `π − β` are opposite. -/
theorem curve_loads_opposite {βi βj : ℝ} (h : βi + βj = 2 * π) :
    Ddg.turningAngle π βj = -Ddg.turningAngle π βi := by
  simp only [turningAngle_eq]
  linarith

/-- **The loads of the three regions at a triple point sum to zero** when their interior angles
there sum to `2π`: `Σ_X (2π/3 − β^X) = 0`.

[formula: cut-length-step-regions] -/
theorem junction_loads_sum {β0 β1 β2 : ℝ} (h : β0 + β1 + β2 = 2 * π) :
    Ddg.turningAngle (2 * π / 3) β0 + Ddg.turningAngle (2 * π / 3) β1 +
      Ddg.turningAngle (2 * π / 3) β2 = 0 := by
  simp only [turningAngle_eq]
  linarith

/-- **At a balanced junction every load vanishes**: three curves meeting at 120°.

[formula: cut-length-step-regions] -/
theorem junction_load_balanced : Ddg.turningAngle (2 * π / 3) (2 * π / 3) = 0 := by
  simp only [turningAngle_eq]
  ring

/-- A region with a sharp corner (`β < 2π/3`) has a positive load there (it retreats). -/
theorem junction_load_pos {β : ℝ} (h : β < 2 * π / 3) : 0 < Ddg.turningAngle (2 * π / 3) β := by
  simp only [turningAngle_eq]
  linarith

section Chain

variable {C : Type*} [Fintype C] [DecidableEq C]

/-- **The systems of §8.5 are those of §3.10 on the boundary nodes of a region**: with `C` the
boundary copies of all regions, positive `ℓ`, symmetric nonnegative couplings, both
`ℓ_v λ⁰_v + Σ κ_σ (λ⁰_v − λ⁰_u) = Θ_v` and
`ℓ_v ρ_v + τ Σ (w̄_σ/|σ|) (ρ_v − ρ_u) = ℓ_v (s_v + λ⁰_v)` have exactly one solution.

[formula: cut-length-step-regions] -/
theorem regionChain_existsUnique {ℓ : C → ℝ} (hℓ : ∀ v, 0 < ℓ v) {κ m : C → C → ℝ}
    (hκ : ∀ p q, κ p q = κ q p) (hκ0 : ∀ p q, 0 ≤ κ p q) (hm : ∀ p q, m p q = m q p)
    (hm0 : ∀ p q, 0 ≤ m p q) {τ : ℝ} (hτ : 0 ≤ τ) (Θ s : C → ℝ) :
    ∃! lam : C → ℝ, (fluxMatrix ℓ κ *ᵥ lam = fun v => Ddg.fluxRhs (Θ v)) ∧
      ∃! ρ : C → ℝ, fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ ρ =
        fun v => Ddg.rateRhs (ℓ v) (s v) (lam v) := by
  obtain ⟨lam, hlam, hu⟩ := onCut_force_existsUnique hℓ hκ hκ0 Θ
  refine ⟨lam, ⟨hlam, onCut_rate_existsUnique hℓ hm hm0 hτ s lam⟩, fun y hy => hu y hy.1⟩

/-- **The total force on the boundary of every region is kept**: for the set `X` of the boundary
nodes of one region (no coupling leaves it), `Σ_{v ∈ X} ℓ_v λ⁰_v = Σ_{v ∈ X} Θ_v`, and
`Σ_{v ∈ X} ℓ_v ρ_v = Σ_{v ∈ X} ℓ_v (s_v + λ⁰_v)`.

[formula: cut-length-step-regions] -/
theorem regionChain_sums {ℓ Θ s lam ρ : C → ℝ} {κ m : C → C → ℝ} {τ : ℝ}
    (hκ : ∀ p q, κ p q = κ q p) (hm : ∀ p q, m p q = m q p)
    (hlam : fluxMatrix ℓ κ *ᵥ lam = fun v => Ddg.fluxRhs (Θ v))
    (hρ : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ ρ = fun v => Ddg.rateRhs (ℓ v) (s v) (lam v))
    (X : Finset C) (hXκ : ∀ p ∈ X, ∀ q, q ∉ X → κ p q = 0)
    (hXm : ∀ p ∈ X, ∀ q, q ∉ X → m p q = 0) :
    ∑ v ∈ X, ℓ v * lam v = ∑ v ∈ X, Θ v ∧
      ∑ v ∈ X, ℓ v * ρ v = ∑ v ∈ X, ℓ v * (s v + lam v) :=
  ⟨fluxRecovery_sum_piece hκ hlam X hXκ, onCut_rate_sum_piece hm hρ X hXm⟩

/-- **The rate of a region lies between the smallest and the largest `s + λ⁰` on its
boundary**, and a region in equilibrium (`s + λ⁰ = 0` on its boundary) does not move.

[formula: cut-length-step-regions] -/
theorem regionChain_bounds {ℓ s lam ρ : C → ℝ} {m : C → C → ℝ} {τ : ℝ} (hℓ : ∀ v, 0 < ℓ v)
    (hm0 : ∀ p q, 0 ≤ m p q) (hτ : 0 ≤ τ)
    (hρ : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ ρ = fun v => Ddg.rateRhs (ℓ v) (s v) (lam v))
    (X : Finset C) (hXm : ∀ p ∈ X, ∀ q, q ∉ X → m p q = 0) {lo hi : ℝ}
    (hlo : ∀ v ∈ X, lo ≤ s v + lam v) (hhi : ∀ v ∈ X, s v + lam v ≤ hi) :
    ∀ v ∈ X, lo ≤ ρ v ∧ ρ v ≤ hi :=
  onCut_rate_bounds_piece hℓ hm0 hτ hρ X hXm hlo hhi

/-- **Opposite rates on the two sides of a curve.** If two chains have the same dual lengths and
couplings (the two sides of a closed curve between two regions) and opposite `s` and `λ⁰`
(`s^j = −s^i`, and `λ⁰^j = −λ⁰^i` by `fluxRecovery_neg` from the opposite forces), their rates
are opposite: `ρ^j = −ρ^i`.

[formula: cut-length-step-regions] -/
theorem onCut_rates_opposite {ℓ s lam r s' lam' r' : C → ℝ} {m : C → C → ℝ} {τ : ℝ}
    (hℓ : ∀ p, 0 < ℓ p) (hm : ∀ p q, m p q = m q p) (hm0 : ∀ p q, 0 ≤ m p q) (hτ : 0 ≤ τ)
    (h : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r = fun p => Ddg.rateRhs (ℓ p) (s p) (lam p))
    (h' : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r' =
      fun p => Ddg.rateRhs (ℓ p) (s' p) (lam' p))
    (hs : ∀ p, s' p = -s p) (hl : ∀ p, lam' p = -lam p) : r' = -r := by
  have e : (fun p => Ddg.rateRhs (ℓ p) (s' p) (lam' p)) =
      fun p => Ddg.fluxRhs (-(ℓ p * (s p + lam p))) := by
    funext p
    rw [rateRhs_eq, fluxRhs_eq, hs p, hl p]
    ring
  rw [e] at h'
  exact fluxRecovery_neg hℓ (fun p q => by rw [hm p q]) (fun p q => mul_nonneg hτ (hm0 p q)) h h'

end Chain

/-- **The explicit step of §8.5 moves an interface at its own rate.** At a vertex of region `k`
with crossings to region `j`: `φ_k += τ e`, `φ_j += τ (e − δ)`, so `φ_k − φ_j` changes by `τ δ`,
`δ` the weighted mean of `ρ(p_a) − ρ(p_b)` over the crossings.

[formula: cut-length-step-regions] -/
theorem regionStep_difference (φk φj τ e δ : ℝ) :
    Ddg.explicitStep φk τ e - Ddg.explicitStep φj τ (e - δ) = (φk - φj) + τ * δ := by
  simp only [explicitStep_eq]
  ring

/-- On a curve between two regions the two copies of a cut point have opposite rates, so the
weighted mean `δ` of `ρ(p_a) − ρ(p_b)` is twice the weighted mean `e` of `ρ(p_a)`. -/
theorem wmean_sub_neg {ι : Type*} (s : Finset ι) (c ρ : ι → ℝ) :
    wmean s c (fun p => ρ p - -ρ p) = 2 * wmean s c ρ := by
  rw [← wmean_smul]
  congr 1
  funext p
  ring

/-- **Two regions: the step of §8.5 is the two-region step.** With `δ = 2e` the crossed function
gets `−τ e`: a state `[φ, −φ]` goes to `[φ + τe, −(φ + τe)]`.

[formula: cut-length-step-regions] -/
theorem regionStep_two (φ τ e : ℝ) :
    Ddg.explicitStep (-φ) τ (e - 2 * e) = -Ddg.explicitStep φ τ e := by
  simp only [explicitStep_eq]
  ring

end Cuts

end DDG
