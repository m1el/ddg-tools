import DdgProofs.Intrinsic.Diamond

/-!
# Intrinsic triangulations: the Delaunay condition, the cotan weight, and what a flip does to the angles

`crates/ddg-intrinsic/SPEC.md` §4 (Delaunay flipping), §5 items 2 and 4, N6; the paper Sec. 4,
Eq. 4, Sec. 4.1.

Notation of `DdgProofs/Intrinsic/Diamond.lean`: an edge `ij` of length `a` with the triangles
`ijk` (`b = ℓ_ik`, `c = ℓ_jk`) and `jil` (`d = ℓ_il`, `e = ℓ_jl`); the angles opposite the edge
are `θ_k = cornerAngle a b c` and `θ_l = cornerAngle a d e`. After a flip the edge `kl` of length
`f` has the diamond `(f, c, e, b, d)` (apexes `j` and `i`).

* **`it-delaunay-condition`**, **`it-cotan-weight`** — for two angles in `(0, π)`:
  `cot θ_k + cot θ_l = sin(θ_k + θ_l)/(sin θ_k sin θ_l)`, so the cotan weight is `≥ 0`, `> 0`,
  `= 0` exactly when `θ_k + θ_l` is `≤ π`, `< π`, `= π` (`cotanWeight_nonneg_iff`, …); the same
  from the five lengths (`delaunay_iff_cotanWeight_nonneg`). An edge that violates the condition
  has a strictly convex diamond (`diamondConvex_of_not_delaunay`); two corners of one triangle
  never violate it (`delaunay_of_same_triangle`). A flip of a strictly convex diamond negates
  the excess `θ_k + θ_l − π` (`flip_delaunayExcess`), so the flip of a violating edge gives an
  edge that satisfies the condition strictly (`flip_of_not_delaunay`).
* **`it-angle-normalization`** — a flip keeps the angle sum at each of its four vertices: the
  new corner at `i` is the sum of the two old ones, the old corner at `k` is the sum of the two
  new ones, and likewise at `j` and `l` (`flip_angle_sums`).

Termination of the flipping algorithm, existence and uniqueness of the intrinsic Delaunay
triangulation and Rippa's theorem are cited (`CITED.md`).
-/

noncomputable section

namespace DdgProofs

namespace Intrinsic

open scoped DdgProofs.Oracle
open Real

/-! ## The oracle's definitions at `ℝ` -/

theorem delaunayExcess_eq (θk θl : ℝ) : Ddg.delaunayExcess θk θl = θk + θl - π := by
  simp only [Ddg.delaunayExcess, Oracle.pi_eq]

theorem delaunayHolds_iff (tol θk θl : ℝ) :
    Ddg.delaunayHolds tol θk θl = true ↔ Ddg.delaunayExcess θk θl ≤ tol := by
  simp only [Ddg.delaunayHolds, decide_eq_true_eq]

theorem cotanWeight_eq (x y : ℝ) : Ddg.cotanWeight x y = x + y := rfl

theorem cotanWeightFromLengths_eq (a b c d e : ℝ) :
    Ddg.cotanWeightFromLengths a b c d e = cornerCot a b c + cornerCot a d e := by
  simp only [Ddg.cotanWeightFromLengths, Ddg.cotanWeight, Oracle.cornerCot_eq]

theorem oppositeAngles_eq {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    Ddg.oppositeAngles a b c d e = (cornerAngle a b c, cornerAngle a d e) := by
  simp only [Ddg.oppositeAngles, Oracle.cornerAngle_eq hk, Oracle.cornerAngle_eq hl]

/-- The corner angle does not depend on the order of the two adjacent sides. -/
theorem cornerAngle_swap (a b c : ℝ) : cornerAngle a c b = cornerAngle a b c := by
  unfold cornerAngle
  congr 1
  ring

/-! ## Cotangents of two angles -/

/-- `cot x + cot y = sin(x + y)/(sin x sin y)`. -/
theorem cot_add_cot {x y : ℝ} (hx : sin x ≠ 0) (hy : sin y ≠ 0) :
    cot x + cot y = sin (x + y) / (sin x * sin y) := by
  rw [cot_eq_cos_div_sin, cot_eq_cos_div_sin, sin_add]
  field_simp
  ring

/-- For `0 < z < 2π`: `sin z ≥ 0 ⟺ z ≤ π`. -/
theorem sin_nonneg_iff_le_pi {z : ℝ} (h0 : 0 < z) (h2 : z < 2 * π) : 0 ≤ sin z ↔ z ≤ π := by
  constructor
  · intro h
    by_contra hcon
    have h1 := Real.sin_pos_of_pos_of_lt_pi (sub_pos.mpr (not_le.mp hcon))
      (by linarith : z - π < π)
    rw [Real.sin_sub_pi] at h1
    linarith
  · exact fun h => Real.sin_nonneg_of_nonneg_of_le_pi h0.le h

theorem div_nonneg_iff_of_pos {x y : ℝ} (hy : 0 < y) : 0 ≤ x / y ↔ 0 ≤ x := by
  constructor
  · intro h
    have := mul_nonneg h hy.le
    rwa [div_mul_cancel₀ _ hy.ne'] at this
  · exact fun h => div_nonneg h hy.le

/-- **The cotan weight has the sign of `π − (θ_k + θ_l)`**: for angles in `(0, π)` it is
`sin(θ_k + θ_l)/(sin θ_k sin θ_l)` with a positive denominator.

[formula: it-cotan-weight] -/
theorem cotanWeight_eq_sin_div {θk θl : ℝ} (hk0 : 0 < θk) (hk : θk < π) (hl0 : 0 < θl)
    (hl : θl < π) :
    Ddg.cotanWeight (cot θk) (cot θl) = sin (θk + θl) / (sin θk * sin θl) ∧
      0 < sin θk * sin θl := by
  have s1 := Real.sin_pos_of_pos_of_lt_pi hk0 hk
  have s2 := Real.sin_pos_of_pos_of_lt_pi hl0 hl
  exact ⟨by rw [cotanWeight_eq, cot_add_cot s1.ne' s2.ne'], mul_pos s1 s2⟩

/-- **The Delaunay condition is the nonnegativity of the cotan weight**: for two angles in
`(0, π)`, `cot θ_k + cot θ_l ≥ 0 ⟺ θ_k + θ_l ≤ π`.

[formula: it-cotan-weight] -/
theorem cotanWeight_nonneg_iff {θk θl : ℝ} (hk0 : 0 < θk) (hk : θk < π) (hl0 : 0 < θl)
    (hl : θl < π) :
    0 ≤ Ddg.cotanWeight (cot θk) (cot θl) ↔ Ddg.delaunayExcess θk θl ≤ 0 := by
  obtain ⟨e, hpos⟩ := cotanWeight_eq_sin_div hk0 hk hl0 hl
  rw [e, delaunayExcess_eq, div_nonneg_iff_of_pos hpos,
    sin_nonneg_iff_le_pi (by linarith) (by linarith)]
  constructor <;> intro h <;> linarith

/-- The strict form: `cot θ_k + cot θ_l > 0 ⟺ θ_k + θ_l < π`.

[formula: it-cotan-weight] -/
theorem cotanWeight_pos_iff {θk θl : ℝ} (hk0 : 0 < θk) (hk : θk < π) (hl0 : 0 < θl)
    (hl : θl < π) :
    0 < Ddg.cotanWeight (cot θk) (cot θl) ↔ Ddg.delaunayExcess θk θl < 0 := by
  obtain ⟨e, hpos⟩ := cotanWeight_eq_sin_div hk0 hk hl0 hl
  rw [e, delaunayExcess_eq, div_pos_iff_of_pos_right hpos,
    sin_pos_iff_lt_pi (by linarith) (by linarith)]
  constructor <;> intro h <;> linarith

/-- The weight vanishes exactly at a tie, `θ_k + θ_l = π` (the four points are cocircular). -/
theorem cotanWeight_eq_zero_iff {θk θl : ℝ} (hk0 : 0 < θk) (hk : θk < π) (hl0 : 0 < θl)
    (hl : θl < π) :
    Ddg.cotanWeight (cot θk) (cot θl) = 0 ↔ Ddg.delaunayExcess θk θl = 0 := by
  have h1 := cotanWeight_nonneg_iff hk0 hk hl0 hl
  have h2 := cotanWeight_pos_iff hk0 hk hl0 hl
  constructor
  · intro h
    have a1 := h1.mp h.ge
    have a2 : ¬ Ddg.delaunayExcess θk θl < 0 := fun hlt => (h2.mpr hlt).ne' h
    linarith [not_lt.mp a2]
  · intro h
    have a1 := h1.mpr h.le
    have a2 : ¬ 0 < Ddg.cotanWeight (cot θk) (cot θl) := fun hpos => (h2.mp hpos).ne h
    linarith [not_lt.mp a2]

/-! ## From the five lengths -/

/-- The oracle's weight from the lengths is the sum of the cotangents of the two opposite
angles. -/
theorem cotanWeightFromLengths_eq_cot {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    Ddg.cotanWeightFromLengths a b c d e =
      Ddg.cotanWeight (cot (cornerAngle a b c)) (cot (cornerAngle a d e)) := by
  rw [cotanWeightFromLengths_eq, cotanWeight_eq, cornerCot_eq_cot hk, cornerCot_eq_cot hl]

/-- **The Delaunay condition from the five lengths**: the oracle's test with tolerance `0` on
its `oppositeAngles` holds exactly when its `cotanWeightFromLengths` is nonnegative.

[formula: it-delaunay-condition] -/
theorem delaunay_iff_cotanWeight_nonneg {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    Ddg.delaunayHolds 0 (Ddg.oppositeAngles a b c d e).1 (Ddg.oppositeAngles a b c d e).2 =
        true ↔ 0 ≤ Ddg.cotanWeightFromLengths a b c d e := by
  rw [delaunayHolds_iff, oppositeAngles_eq hk hl, cotanWeightFromLengths_eq_cot hk hl,
    cotanWeight_nonneg_iff (cornerAngle_pos hk) (cornerAngle_lt_pi hk) (cornerAngle_pos hl)
      (cornerAngle_lt_pi hl)]

/-- The same as a statement about the angles: `θ_k + θ_l ≤ π ⟺` weight `≥ 0`.

[formula: it-delaunay-condition] -/
theorem angle_sum_le_pi_iff_cotanWeight_nonneg {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    cornerAngle a b c + cornerAngle a d e ≤ π ↔ 0 ≤ Ddg.cotanWeightFromLengths a b c d e := by
  rw [cotanWeightFromLengths_eq_cot hk hl,
    cotanWeight_nonneg_iff (cornerAngle_pos hk) (cornerAngle_lt_pi hk) (cornerAngle_pos hl)
      (cornerAngle_lt_pi hl), delaunayExcess_eq]
  constructor <;> intro h <;> linarith

/-- An edge violates the condition exactly when its cotan weight is negative. -/
theorem pi_lt_angle_sum_iff_cotanWeight_neg {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    π < cornerAngle a b c + cornerAngle a d e ↔ Ddg.cotanWeightFromLengths a b c d e < 0 := by
  rw [← not_le, angle_sum_le_pi_iff_cotanWeight_nonneg hk hl, not_le]

/-- **Two corners of one triangle sum to less than `π`**: an edge with the same triangle on both
sides (an endpoint of degree one) satisfies the condition strictly (the paper's remark in
Sec. 4.1). -/
theorem delaunay_of_same_triangle {a b c : ℝ} (h : IsTriangle a b c) :
    Ddg.delaunayExcess (cornerAngle a b c) (cornerAngle b c a) < 0 := by
  have h1 := cornerAngle_add_cornerAngle_add_cornerAngle h
  have h2 := cornerAngle_pos h.rotate.rotate
  rw [delaunayExcess_eq]
  linarith

/-- The angles of the quadrilateral at the two ends of the edge sum to `2π − (θ_k + θ_l)`. -/
theorem quad_angle_sum {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    (cornerAngle c a b + cornerAngle e a d) + (cornerAngle b c a + cornerAngle d e a) =
      2 * π - (cornerAngle a b c + cornerAngle a d e) := by
  have h1 := cornerAngle_add_cornerAngle_add_cornerAngle hk
  have h2 := cornerAngle_add_cornerAngle_add_cornerAngle hl
  linarith

/-- An edge whose opposite angles sum to `π` or more (a violation, or a tie) has a strictly
convex diamond. -/
theorem diamondConvex_of_pi_le {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) (h : π ≤ cornerAngle a b c + cornerAngle a d e) :
    DiamondConvex a b c d e := by
  have hq := quad_angle_sum hk hl
  have a1 := cornerAngle_pos hk.rotate.rotate
  have a3 := cornerAngle_pos hl.rotate.rotate
  have b1 := cornerAngle_pos hk.rotate
  have b3 := cornerAngle_pos hl.rotate
  rw [diamondConvex_iff_angles hk hl]
  constructor <;> linarith

/-- **An edge that violates the Delaunay condition can be flipped**: if `θ_k + θ_l > π` the
diamond is strictly convex (the paper's remark in Sec. 4.1).

[formula: it-delaunay-condition] -/
theorem diamondConvex_of_not_delaunay {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) (h : π < cornerAngle a b c + cornerAngle a d e) :
    DiamondConvex a b c d e :=
  diamondConvex_of_pi_le hk hl h.le

/-- The same from the weight: an edge with a negative cotan weight can be flipped.

[formula: it-cotan-weight] -/
theorem diamondConvex_of_cotanWeight_neg {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) (h : Ddg.cotanWeightFromLengths a b c d e < 0) :
    DiamondConvex a b c d e :=
  diamondConvex_of_not_delaunay hk hl ((pi_lt_angle_sum_iff_cotanWeight_neg hk hl).mpr h)

/-! ## What a flip does to the angles -/

/-- The new corner at `i` (in the triangle `ilk`) is the sum of the two old corners at `i`,
when that sum is at most `π`. -/
theorem flip_angle_i {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle c a b + cornerAngle e a d ≤ π) :
    cornerAngle (Ddg.flippedLengthFromLengths a b c d e) b d =
      cornerAngle c a b + cornerAngle e a d := by
  have hb := hk.pos_b
  have hd := hl.pos_b
  have a1 := cornerAngle_pos hk.rotate.rotate
  have a3 := cornerAngle_pos hl.rotate.rotate
  have hsq := flippedLengthFromLengths_sq hk hl
  have e1 : (b ^ 2 + d ^ 2 - Ddg.flippedLengthFromLengths a b c d e ^ 2) / (2 * b * d) =
      cos (cornerAngle c a b + cornerAngle e a d) := by
    rw [hsq]
    field_simp
    ring
  rw [show cornerAngle (Ddg.flippedLengthFromLengths a b c d e) b d =
      arccos ((b ^ 2 + d ^ 2 - Ddg.flippedLengthFromLengths a b c d e ^ 2) / (2 * b * d))
      from rfl, e1, Real.arccos_cos (by linarith) h]

/-- The new corner at `j` (in the triangle `jkl`) is the sum of the two old corners at `j`. -/
theorem flip_angle_j {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle b c a + cornerAngle d e a ≤ π) :
    cornerAngle (Ddg.flippedLengthFromLengths a b c d e) c e =
      cornerAngle b c a + cornerAngle d e a := by
  have hc := hk.pos_c
  have he := hl.pos_c
  have b1 := cornerAngle_pos hk.rotate
  have b3 := cornerAngle_pos hl.rotate
  have hsq := flippedLengthFromLengths_sq' hk hl
  have e1 : (c ^ 2 + e ^ 2 - Ddg.flippedLengthFromLengths a b c d e ^ 2) / (2 * c * e) =
      cos (cornerAngle b c a + cornerAngle d e a) := by
    rw [hsq]
    field_simp
    ring
  rw [show cornerAngle (Ddg.flippedLengthFromLengths a b c d e) c e =
      arccos ((c ^ 2 + e ^ 2 - Ddg.flippedLengthFromLengths a b c d e ^ 2) / (2 * c * e))
      from rfl, e1, Real.arccos_cos (by linarith) h]

/-- **The diamond of the flipped edge is strictly convex too** (it is the same quadrilateral):
with `f` the flipped length, `DiamondConvex f c e b d`. -/
theorem flip_convex {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) :
    DiamondConvex (Ddg.flippedLengthFromLengths a b c d e) c e b d := by
  have ha := hk.pos_a
  obtain ⟨hk0, hl0, -, e2, e3, e4, e5⟩ := diamond_isometric hk hl
  obtain ⟨c1, c2⟩ := diamondConvex_sda hconv
  have hJ : 0 < sda (apexK a b c) (apexL a d e) (a, 0) := by
    rw [sda_rotate]
    exact c2
  have hI : sda (apexK a b c) (apexL a d e) (0, 0) < 0 := by
    rw [sda_rotate, sda_swap_right]
    linarith
  obtain ⟨r1, r2⟩ := flipArea_realized hJ hI
  rw [dist2_apex, dist2_comm (apexK a b c) (a, 0), dist2_comm (apexL a d e) (a, 0),
    dist2_comm (apexK a b c) (0, 0), dist2_comm (apexL a d e) (0, 0), e2, e3, e4, e5] at r1 r2
  obtain ⟨s1, s2⟩ := diamond_sda a (Ddg.diamondX a b c) (Ddg.diamondH a b c)
    (Ddg.diamondX a d e) (Ddg.diamondH a d e)
  have t1 : sda (apexK a b c) (0, 0) (a, 0) = a * Ddg.diamondH a b c := by
    rw [sda_rotate, sda_rotate]
    exact s1
  have t2 : sda (apexL a d e) (a, 0) (0, 0) = a * Ddg.diamondH a d e := by
    rw [sda_rotate, sda_swap_right]
    unfold apexL
    rw [s2]
    ring
  unfold DiamondConvex
  rw [flipConvex_iff, r1, r2, t1, t2]
  exact ⟨mul_pos ha hk0, mul_pos ha hl0⟩

/-- **A flip keeps the angle sum at each of the four vertices of the diamond.** For a strictly
convex diamond with flipped length `f`: the corner of the new triangle `ilk` at `i` is the sum of
the two old corners at `i`; the corner of `jkl` at `j` is the sum of the two old corners at `j`;
the old corner at `k` is the sum of the corners of `jkl` and `ilk` at `k`; the old corner at `l`
is the sum of their corners at `l`. So the totals `Θ` that the angle normalization uses do not
change (SPEC §3, §5 item 2).

[formula: it-angle-normalization] -/
theorem flip_angle_sums {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) :
    cornerAngle (Ddg.flippedLengthFromLengths a b c d e) b d =
        cornerAngle c a b + cornerAngle e a d ∧
      cornerAngle (Ddg.flippedLengthFromLengths a b c d e) c e =
        cornerAngle b c a + cornerAngle d e a ∧
      cornerAngle e (Ddg.flippedLengthFromLengths a b c d e) c +
          cornerAngle d (Ddg.flippedLengthFromLengths a b c d e) b = cornerAngle a b c ∧
      cornerAngle c e (Ddg.flippedLengthFromLengths a b c d e) +
          cornerAngle b d (Ddg.flippedLengthFromLengths a b c d e) = cornerAngle a d e := by
  obtain ⟨h1, h2⟩ := (diamondConvex_iff_angles hk hl).mp hconv
  obtain ⟨t1, t2⟩ := flip_isTriangle hk hl hconv
  have hk' := t2.swap
  have hconv' := flip_convex hk hl hconv
  obtain ⟨g1, g2⟩ := (diamondConvex_iff_angles hk' t1).mp hconv'
  have hback := flip_back hk hl hconv
  have k1 := flip_angle_i hk' t1 g1.le
  have k2 := flip_angle_j hk' t1 g2.le
  rw [hback] at k1 k2
  refine ⟨flip_angle_i hk hl h1.le, flip_angle_j hk hl h2.le, ?_, ?_⟩
  · rw [← k1, cornerAngle_swap]
  · rw [← k2, cornerAngle_swap]

/-- **A flip negates the Delaunay excess**: for a strictly convex diamond, the opposite angles
of the new edge `kl` (at `j` and at `i`) sum to `2π − (θ_k + θ_l)`.

[formula: it-delaunay-condition] -/
theorem flip_delaunayExcess {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) :
    Ddg.delaunayExcess (cornerAngle (Ddg.flippedLengthFromLengths a b c d e) c e)
        (cornerAngle (Ddg.flippedLengthFromLengths a b c d e) b d) =
      -Ddg.delaunayExcess (cornerAngle a b c) (cornerAngle a d e) := by
  obtain ⟨h1, h2, -, -⟩ := flip_angle_sums hk hl hconv
  have hq := quad_angle_sum hk hl
  rw [delaunayExcess_eq, delaunayExcess_eq, h1, h2]
  linarith

/-- **After the flip of a violating edge the new edge satisfies the condition strictly**: if
`θ_k + θ_l > π`, the diamond is strictly convex, the two new triangles satisfy the strict
triangle inequalities, the oracle's excess of the new diamond `(f, c, e, b, d)` is negative and
its cotan weight positive.

[formula: it-delaunay-condition] -/
theorem flip_of_not_delaunay {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : π < cornerAngle a b c + cornerAngle a d e) :
    DiamondConvex a b c d e ∧
      Ddg.delaunayExcess
          (Ddg.oppositeAngles (Ddg.flippedLengthFromLengths a b c d e) c e b d).1
          (Ddg.oppositeAngles (Ddg.flippedLengthFromLengths a b c d e) c e b d).2 < 0 ∧
      0 < Ddg.cotanWeightFromLengths (Ddg.flippedLengthFromLengths a b c d e) c e b d := by
  have hconv := diamondConvex_of_not_delaunay hk hl h
  obtain ⟨t1, t2⟩ := flip_isTriangle hk hl hconv
  have hk' := t2.swap
  have hx := flip_delaunayExcess hk hl hconv
  have hneg : Ddg.delaunayExcess
      (cornerAngle (Ddg.flippedLengthFromLengths a b c d e) c e)
      (cornerAngle (Ddg.flippedLengthFromLengths a b c d e) b d) < 0 := by
    rw [hx, delaunayExcess_eq]
    linarith
  refine ⟨hconv, ?_, ?_⟩
  · rw [oppositeAngles_eq hk' t1]
    exact hneg
  · rw [cotanWeightFromLengths_eq_cot hk' t1]
    exact (cotanWeight_pos_iff (cornerAngle_pos hk') (cornerAngle_lt_pi hk') (cornerAngle_pos t1)
      (cornerAngle_lt_pi t1)).mpr hneg

/-- At a tie (`θ_k + θ_l = π`) both diagonals satisfy the condition with equality: the flip
gives the excess `0` again (SPEC N6). -/
theorem flip_of_tie {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle a b c + cornerAngle a d e = π) :
    DiamondConvex a b c d e ∧
      Ddg.delaunayExcess (cornerAngle (Ddg.flippedLengthFromLengths a b c d e) c e)
        (cornerAngle (Ddg.flippedLengthFromLengths a b c d e) b d) = 0 := by
  have hconv := diamondConvex_of_pi_le hk hl h.ge
  refine ⟨hconv, ?_⟩
  rw [flip_delaunayExcess hk hl hconv, delaunayExcess_eq]
  linarith

end Intrinsic

end DdgProofs
