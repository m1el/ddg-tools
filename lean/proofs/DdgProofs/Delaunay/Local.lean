import DdgProofs.Intrinsic.Delaunay

/-!
# Delaunay triangulations of piecewise flat surfaces, layer 1: one diamond

Everything in this file is about a single *diamond*: an edge `ij` of length `a` with the two
triangles `ijk` and `jil` on its sides, given intrinsically by five lengths

* `a = ℓ_ij`, `b = ℓ_ik`, `c = ℓ_jk` (the triangle `ijk`, `IsTriangle a b c`),
* `d = ℓ_il`, `e = ℓ_jl` (the triangle `jil`, `IsTriangle a d e`).

Corner angles are those of `DdgProofs/Triangle.lean` (`cornerAngle x y z` is the angle opposite
the side `x`): `θ_k = cornerAngle a b c`, `θ_l = cornerAngle a d e` are the angles opposite the
edge; at `i` the diamond has the corners `cornerAngle c a b` (in `ijk`) and `cornerAngle e a d`
(in `jil`), at `j` the corners `cornerAngle b c a` and `cornerAngle d e a`. A diamond is a flat
quadrilateral, so the statements hold verbatim for an interior edge of any triangulation of a
piecewise flat surface (also when vertices coincide on the surface, e.g. `k = l`): nothing here
uses an embedding of the surface.

After a flip the edge `kl` has length `f = flipLength a b c d e`; its diamond is
`(f, c, e, b, d)` (apex `j` with `ℓ_kj = c`, `ℓ_lj = e`; apex `i` with `ℓ_ki = b`,
`ℓ_li = d`); the new triangles are `ilk` with sides `(f, b, d)` and `jkl` with sides `(f, e, c)`.

## Main statements

* (a) **Equivalent local Delaunay criteria.** `IsLocallyDelaunay a b c d e` is
  `θ_k + θ_l ≤ π` (B–S Def. 8). It is equivalent to `cot θ_k + cot θ_l ≥ 0`
  (`isLocallyDelaunay_iff_cot_nonneg`; strict and tie forms `pi_lt_iff_cot_neg`,
  `eq_pi_iff_cot_eq_zero`), and, in the planar layout of the diamond, to "`l` does not lie strictly
  inside the circumcircle of `i, j, k`" (`isLocallyDelaunay_iff_not_inside_circumcircle`) and
  symmetrically for `k` and the circle through `i, j, l`
  (`isLocallyDelaunay_iff_not_inside_circumcircle'`). The underlying identity is the power of
  `l` with respect to the circumcircle of `ijk`: `2 A_jil (cot θ_k + cot θ_l)`
  (`power_apexBelow`, B–S Lemma 9 in lengths).
* (b) **The flip of a non-Delaunay edge.** If `θ_k + θ_l > π` (or `= π`) the diamond is strictly
  convex (`isConvexDiamond_of_pi_le`), the two new triangles are nondegenerate
  (`flipLength_isTriangle`), the new edge is locally Delaunay, strictly
  (`flip_isLocallyDelaunay`, `flip_opposite_angles`: the new opposite angles sum to
  `2π − (θ_k + θ_l)`), the new diamond is strictly convex (`flip_isConvexDiamond`), flipping back
  gives the old edge (`flipLength_flipLength`), the total area is kept (`flip_area`), and the
  corner angles add up at all four vertices (`flip_corner_angles`): at `i` and `j` the new corner
  is the sum of the two old ones, at `k` and `l` the old corner is the sum of the two new ones.
  So the cone angle (the sum of the corner angles) at every vertex of the surface is unchanged.
  `flip_of_not_isLocallyDelaunay` collects these for use by a triangulation-level theory.

The potential that decreases under a flip (c) is in `DdgProofs/Delaunay/Potential.lean`,
Rippa's theorem for one diamond (d) in `DdgProofs/Delaunay/Rippa.lean`, Ruppert's encroachment
criterion (e) in `DdgProofs/Delaunay/Encroachment.lean`.

## References

* A. I. Bobenko, B. A. Springborn, *A discrete Laplace–Beltrami operator for simplicial
  surfaces*, Discrete Comput. Geom. 38(4), 740–756, 2007, arXiv:math/0503219 ("B–S"; numbering
  of arXiv v3): Definition 8 and Lemma 9 (locally Delaunay ⟺ `θ_k + θ_l ≤ π` ⟺ the cotan
  weight is `≥ 0` ⟺ empty circumcircle of the unfolded diamond), Proposition 11 (a non-Delaunay
  edge can be flipped).
* C. L. Lawson, *Software for C¹ surface interpolation*, in: Mathematical Software III,
  Academic Press, 1977, 161–194 (the flip algorithm; the angle form of the local criterion).
* R. Sibson, *Locally equiangular triangulations*, The Computer Journal 21(3), 243–245, 1978
  (local and global Delaunay criteria in the plane).
* C. Indermitte, Th. M. Liebling, M. Troyanov, H. Clémençon, *Voronoi diagrams on piecewise flat
  surfaces and an application to biological growth*, Theoret. Comput. Sci. 263(1–2), 263–274,
  2001, Sec. 4 (legal edges and flips on piecewise flat surfaces).
* M. Fisher, B. Springborn, A. I. Bobenko, P. Schröder, *An algorithm for the construction of
  intrinsic Delaunay triangulations with applications to digital geometry processing*,
  Computing 81(2–3), 199–213, 2007 (the intrinsic flip algorithm on edge lengths).
* N. Sharp, Y. Soliman, K. Crane, *Navigating intrinsic triangulations*, ACM TOG 38(4), 2019,
  Sec. 3.3.1 (Fig. 9: the flipped length) and Sec. 4 (Eq. 4).

B–S, Indermitte et al. and the Signpost paper were read (`CITED.md`, "Navigating Intrinsic
Triangulations"); the bibliographic data of Lawson, Sibson and Fisher et al. were not checked
against the works here.

## Dependencies

Mathlib, `DdgProofs/Triangle.lean` (the corner angle, cotangent and area from three lengths)
and, for the proofs only, the planar layout of a diamond in `DdgProofs/Intrinsic/Diamond.lean`
and `DdgProofs/Intrinsic/Delaunay.lean` (which are stated about the oracle's `Ddg.*`
definitions). The statements of this file mention only `Triangle.lean` definitions and the plane
`ℝ × ℝ` with `dist2` (`DdgProofs/Bff/Quality.lean`); `flipLength_eq_flippedLengthFromLengths`
is the bridge to the oracle's flipped length.
-/

noncomputable section

namespace DDG

namespace Delaunay

open Real
open scoped DDG.Oracle

variable {a b c d e : ℝ}

/-! ## Definitions -/

/-- The edge `ij` of the diamond `(a, b, c, d, e)` is **locally Delaunay**: the two angles
opposite it sum to at most `π` (B–S Def. 8, Lemma 9; Signpost paper Eq. 4). -/
def IsLocallyDelaunay (a b c d e : ℝ) : Prop := cornerAngle a b c + cornerAngle a d e ≤ π

/-- The diamond is **strictly convex**: its angles at the two ends `i`, `j` of the edge are less
than `π`, so that the other diagonal `kl` lies inside it. -/
def IsConvexDiamond (a b c d e : ℝ) : Prop :=
  cornerAngle c a b + cornerAngle e a d < π ∧ cornerAngle b c a + cornerAngle d e a < π

/-- **The length of the flipped edge `kl`** from the five lengths:
`ℓ_kl² = b² + d² − [(a² + b² − c²)(a² + d² − e²) − 16 A_ijk A_jil] / 2a²`
(the law of cosines at `i` with the angle `cornerAngle c a b + cornerAngle e a d`, written in
lengths; `flipLength_sq_law_of_cosines`). -/
def flipLength (a b c d e : ℝ) : ℝ :=
  √(b ^ 2 + d ^ 2 - ((a ^ 2 + b ^ 2 - c ^ 2) * (a ^ 2 + d ^ 2 - e ^ 2) -
    16 * triArea a b c * triArea a d e) / (2 * a ^ 2))

/-- The apex `k` of the planar layout `i = (0, 0)`, `j = (a, 0)`, above the edge. -/
def apexAbove (a b c : ℝ) : ℝ × ℝ := ((a ^ 2 + b ^ 2 - c ^ 2) / (2 * a), 2 * triArea a b c / a)

/-- The apex `l` of the planar layout, below the edge. -/
def apexBelow (a d e : ℝ) : ℝ × ℝ :=
  ((a ^ 2 + d ^ 2 - e ^ 2) / (2 * a), -(2 * triArea a d e / a))

/-- The circumcentre of `i = (0, 0)`, `j = (a, 0)` and `k = apexAbove a b c`:
`(a/2, (a/2) cot θ_k)`. -/
def circumcenterAbove (a b c : ℝ) : ℝ × ℝ := (a / 2, a * cornerCot a b c / 2)

/-- The circumcentre of `i`, `j` and `l = apexBelow a d e`: `(a/2, −(a/2) cot θ_l)`. -/
def circumcenterBelow (a d e : ℝ) : ℝ × ℝ := (a / 2, -(a * cornerCot a d e / 2))

/-! ## Bridges to `DdgProofs/Intrinsic` -/

theorem apexAbove_eq (a b c : ℝ) : apexAbove a b c = Intrinsic.apexK a b c := by
  unfold apexAbove Intrinsic.apexK
  rw [Intrinsic.diamondX_eq, Intrinsic.diamondH_eq]
  congr 1
  ring

theorem apexBelow_eq (a d e : ℝ) : apexBelow a d e = Intrinsic.apexL a d e := by
  unfold apexBelow Intrinsic.apexL
  rw [Intrinsic.diamondX_eq, Intrinsic.diamondH_eq]
  congr 1
  ring

/-- `flipLength` is the oracle's (and the Rust code's) flipped length.

[formula: it-flip-length] -/
theorem flipLength_eq_flippedLengthFromLengths (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    flipLength a b c d e = Ddg.flippedLengthFromLengths a b c d e := by
  rw [flipLength, ← Intrinsic.flippedLengthFromLengths_sq_lengths hk hl, Real.sqrt_sq]
  rw [Intrinsic.flippedLengthFromLengths_eq, Intrinsic.flippedLength_eq]
  exact Real.sqrt_nonneg _

theorem isConvexDiamond_iff_diamondConvex (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    IsConvexDiamond a b c d e ↔ Intrinsic.DiamondConvex a b c d e :=
  (Intrinsic.diamondConvex_iff_angles hk hl).symm

theorem flipLength_nonneg (a b c d e : ℝ) : 0 ≤ flipLength a b c d e := Real.sqrt_nonneg _

/-- **The paper's form of the flipped length** (Signpost paper Fig. 9): the law of cosines at `i`,
`ℓ_kl² = ℓ_ik² + ℓ_il² − 2 ℓ_ik ℓ_il cos(θ_i^{jk} + θ_i^{lj})`, and the same at `j`. -/
theorem flipLength_sq_law_of_cosines (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    flipLength a b c d e ^ 2 =
        b ^ 2 + d ^ 2 - 2 * b * d * cos (cornerAngle c a b + cornerAngle e a d) ∧
      flipLength a b c d e ^ 2 =
        c ^ 2 + e ^ 2 - 2 * c * e * cos (cornerAngle b c a + cornerAngle d e a) := by
  rw [flipLength_eq_flippedLengthFromLengths hk hl]
  exact ⟨Intrinsic.flippedLengthFromLengths_sq hk hl, Intrinsic.flippedLengthFromLengths_sq' hk hl⟩

/-! ## (a) The local Delaunay criteria -/

/-- **Locally Delaunay ⟺ nonnegative cotan weight** (B–S Lemma 9): `θ_k + θ_l ≤ π ⟺
cot θ_k + cot θ_l ≥ 0`, with the cotangents computed from the lengths
(`cornerCot_eq_cot`).

[formula: it-delaunay-condition] -/
theorem isLocallyDelaunay_iff_cot_nonneg (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    IsLocallyDelaunay a b c d e ↔ 0 ≤ cornerCot a b c + cornerCot a d e := by
  rw [IsLocallyDelaunay, Intrinsic.angle_sum_le_pi_iff_cotanWeight_nonneg hk hl,
    Intrinsic.cotanWeightFromLengths_eq]

/-- The strict form: `θ_k + θ_l > π ⟺ cot θ_k + cot θ_l < 0`. -/
theorem pi_lt_iff_cot_neg (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    π < cornerAngle a b c + cornerAngle a d e ↔ cornerCot a b c + cornerCot a d e < 0 := by
  rw [Intrinsic.pi_lt_angle_sum_iff_cotanWeight_neg hk hl, Intrinsic.cotanWeightFromLengths_eq]

/-- `¬ IsLocallyDelaunay` is `θ_k + θ_l > π`. -/
theorem not_isLocallyDelaunay_iff :
    ¬ IsLocallyDelaunay a b c d e ↔ π < cornerAngle a b c + cornerAngle a d e := by
  rw [IsLocallyDelaunay, not_le]

/-- The tie: `θ_k + θ_l = π ⟺ cot θ_k + cot θ_l = 0` (the four points are cocircular). -/
theorem eq_pi_iff_cot_eq_zero (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    cornerAngle a b c + cornerAngle a d e = π ↔ cornerCot a b c + cornerCot a d e = 0 := by
  have h := Intrinsic.cotanWeight_eq_zero_iff (cornerAngle_pos hk) (cornerAngle_lt_pi hk)
    (cornerAngle_pos hl) (cornerAngle_lt_pi hl)
  rw [Intrinsic.cotanWeight_eq, ← cornerCot_eq_cot hk, ← cornerCot_eq_cot hl,
    Intrinsic.delaunayExcess_eq] at h
  rw [h]
  constructor <;> intro h' <;> linarith

/-- The criterion does not depend on which end of the edge is called `i`, nor on which triangle
is called `ijk`. -/
theorem isLocallyDelaunay_comm :
    IsLocallyDelaunay a d e b c ↔ IsLocallyDelaunay a b c d e := by
  unfold IsLocallyDelaunay
  rw [add_comm]

theorem isLocallyDelaunay_mirror :
    IsLocallyDelaunay a c b e d ↔ IsLocallyDelaunay a b c d e := by
  unfold IsLocallyDelaunay
  rw [Intrinsic.cornerAngle_swap a b c, Intrinsic.cornerAngle_swap a d e]

/-! ### The empty circumcircle in the planar layout -/

/-- **The planar layout of the diamond**: `i = (0, 0)`, `j = (a, 0)`, `k = apexAbove a b c`
above the edge and `l = apexBelow a d e` below it have the five given distances, and the
distance `kl` is `flipLength a b c d e`. -/
theorem layout_isometric (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    0 < (apexAbove a b c).2 ∧ (apexBelow a d e).2 < 0 ∧
      dist2 (0, 0) (apexAbove a b c) = b ∧ dist2 (a, 0) (apexAbove a b c) = c ∧
      dist2 (0, 0) (apexBelow a d e) = d ∧ dist2 (a, 0) (apexBelow a d e) = e ∧
      dist2 (apexAbove a b c) (apexBelow a d e) = flipLength a b c d e := by
  obtain ⟨h1, h2, -, h4, h5, h6, h7⟩ := Intrinsic.diamond_isometric hk hl
  rw [apexAbove_eq, apexBelow_eq, Intrinsic.dist2_apex, flipLength_eq_flippedLengthFromLengths hk hl]
  refine ⟨h1, ?_, h4, h5, h6, h7, rfl⟩
  simp only [Intrinsic.apexL]
  linarith

private theorem sq_of_dist2 {p q : ℝ × ℝ} {r : ℝ} (h : dist2 p q = r) :
    r ^ 2 = (q.1 - p.1) ^ 2 + (q.2 - p.2) ^ 2 := by
  rw [← h, Intrinsic.dist2_sq']

/-- The squared distances of the layout as polynomial identities in its coordinates. -/
private theorem layout_sq (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    b ^ 2 = (apexAbove a b c).1 ^ 2 + (apexAbove a b c).2 ^ 2 ∧
      c ^ 2 = ((apexAbove a b c).1 - a) ^ 2 + (apexAbove a b c).2 ^ 2 ∧
      d ^ 2 = (apexBelow a d e).1 ^ 2 + (apexBelow a d e).2 ^ 2 ∧
      e ^ 2 = ((apexBelow a d e).1 - a) ^ 2 + (apexBelow a d e).2 ^ 2 := by
  obtain ⟨-, -, h3, h4, h5, h6, -⟩ := layout_isometric hk hl
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [sq_of_dist2 h3]; ring
  · rw [sq_of_dist2 h4]; ring
  · rw [sq_of_dist2 h5]; ring
  · rw [sq_of_dist2 h6]; ring

/-- `circumcenterAbove a b c` is at the same distance from `i`, `j` and `k`. -/
theorem circumcenterAbove_equidistant (hk : IsTriangle a b c) :
    dist2 (circumcenterAbove a b c) (a, 0) = dist2 (circumcenterAbove a b c) (0, 0) ∧
      dist2 (circumcenterAbove a b c) (apexAbove a b c) =
        dist2 (circumcenterAbove a b c) (0, 0) := by
  obtain ⟨h1, h2, -, -⟩ := layout_sq hk hk
  have ha := hk.pos_a
  have hT := triArea_pos hk
  constructor
  · unfold dist2 circumcenterAbove
    congr 1
    ring
  · unfold dist2
    congr 1
    have hc : a * cornerCot a b c * (2 * triArea a b c / a) = (b ^ 2 + c ^ 2 - a ^ 2) / 2 := by
      unfold cornerCot
      field_simp
      ring
    unfold circumcenterAbove
    simp only [apexAbove] at h1 h2 ⊢
    nlinarith [hc]

/-- `circumcenterBelow a d e` is at the same distance from `i`, `j` and `l`. -/
theorem circumcenterBelow_equidistant (hl : IsTriangle a d e) :
    dist2 (circumcenterBelow a d e) (a, 0) = dist2 (circumcenterBelow a d e) (0, 0) ∧
      dist2 (circumcenterBelow a d e) (apexBelow a d e) =
        dist2 (circumcenterBelow a d e) (0, 0) := by
  obtain ⟨-, -, h1, h2⟩ := layout_sq hl hl
  have ha := hl.pos_a
  have hT := triArea_pos hl
  constructor
  · unfold dist2 circumcenterBelow
    congr 1
    ring
  · unfold dist2
    congr 1
    have hc : a * cornerCot a d e * (2 * triArea a d e / a) = (d ^ 2 + e ^ 2 - a ^ 2) / 2 := by
      unfold cornerCot
      field_simp
      ring
    unfold circumcenterBelow
    simp only [apexBelow] at h1 h2 ⊢
    nlinarith [hc]

/-- **The power of `l` with respect to the circumcircle of `i, j, k`** is
`|O l|² − R² = 2 A_jil (cot θ_k + cot θ_l)`: negative exactly when the edge is not locally
Delaunay. -/
theorem power_apexBelow (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    sqDist2 (circumcenterAbove a b c) (apexBelow a d e) -
        sqDist2 (circumcenterAbove a b c) (0, 0) =
      2 * triArea a d e * (cornerCot a b c + cornerCot a d e) := by
  obtain ⟨-, -, h1, h2⟩ := layout_sq hk hl
  have ha := hk.pos_a
  have hT := triArea_pos hl
  have hc : a * cornerCot a d e * (2 * triArea a d e / a) = (d ^ 2 + e ^ 2 - a ^ 2) / 2 := by
    unfold cornerCot
    field_simp
    ring
  unfold sqDist2 circumcenterAbove
  simp only [apexBelow] at h1 h2 ⊢
  have e1 : 2 * triArea a d e * (cornerCot a b c + cornerCot a d e) =
      a * cornerCot a b c * (2 * triArea a d e / a) +
        a * cornerCot a d e * (2 * triArea a d e / a) := by
    field_simp
  rw [e1, hc]
  nlinarith

/-- The power of `k` with respect to the circumcircle of `i, j, l`:
`2 A_ijk (cot θ_k + cot θ_l)`. -/
theorem power_apexAbove (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    sqDist2 (circumcenterBelow a d e) (apexAbove a b c) -
        sqDist2 (circumcenterBelow a d e) (0, 0) =
      2 * triArea a b c * (cornerCot a b c + cornerCot a d e) := by
  obtain ⟨h1, h2, -, -⟩ := layout_sq hk hl
  have ha := hk.pos_a
  have hT := triArea_pos hk
  have hc : a * cornerCot a b c * (2 * triArea a b c / a) = (b ^ 2 + c ^ 2 - a ^ 2) / 2 := by
    unfold cornerCot
    field_simp
    ring
  unfold sqDist2 circumcenterBelow
  simp only [apexAbove] at h1 h2 ⊢
  have e1 : 2 * triArea a b c * (cornerCot a b c + cornerCot a d e) =
      a * cornerCot a b c * (2 * triArea a b c / a) +
        a * cornerCot a d e * (2 * triArea a b c / a) := by
    field_simp
  rw [e1, hc]
  nlinarith

private theorem dist2_lt_iff (p q r : ℝ × ℝ) :
    dist2 p q < dist2 p r ↔ sqDist2 p q - sqDist2 p r < 0 := by
  rw [dist2, dist2, Real.sqrt_lt_sqrt_iff (by positivity), sub_neg]
  rfl

private theorem dist2_eq_iff (p q r : ℝ × ℝ) :
    dist2 p q = dist2 p r ↔ sqDist2 p q - sqDist2 p r = 0 := by
  rw [dist2, dist2, Real.sqrt_inj (by positivity) (by positivity), sub_eq_zero]
  rfl

/-- **Not locally Delaunay ⟺ `l` lies strictly inside the circumcircle of `i, j, k`** in the
planar layout of the diamond (B–S Lemma 9; Sibson 1978 in the plane).

[formula: it-delaunay-condition] -/
theorem inside_circumcircle_iff (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    dist2 (circumcenterAbove a b c) (apexBelow a d e) < dist2 (circumcenterAbove a b c) (0, 0) ↔
      ¬ IsLocallyDelaunay a b c d e := by
  rw [dist2_lt_iff, power_apexBelow hk hl, not_isLocallyDelaunay_iff, pi_lt_iff_cot_neg hk hl]
  have hT := triArea_pos hl
  constructor
  · intro h
    by_contra h'
    push Not at h'
    nlinarith
  · intro h
    nlinarith

/-- **Locally Delaunay ⟺ `l` is not strictly inside the circumcircle of `i, j, k`.** -/
theorem isLocallyDelaunay_iff_not_inside_circumcircle (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    IsLocallyDelaunay a b c d e ↔
      dist2 (circumcenterAbove a b c) (0, 0) ≤ dist2 (circumcenterAbove a b c) (apexBelow a d e) := by
  rw [← not_lt, inside_circumcircle_iff hk hl, not_not]

/-- **A tie ⟺ `l` lies on the circumcircle of `i, j, k`.** -/
theorem on_circumcircle_iff (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    dist2 (circumcenterAbove a b c) (apexBelow a d e) = dist2 (circumcenterAbove a b c) (0, 0) ↔
      cornerAngle a b c + cornerAngle a d e = π := by
  rw [dist2_eq_iff, power_apexBelow hk hl, eq_pi_iff_cot_eq_zero hk hl]
  have hT := (triArea_pos hl).ne'
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h | h
    · exfalso; exact hT (by linarith)
    · exact h
  · intro h
    rw [h, mul_zero]

/-- The symmetric criterion: not locally Delaunay ⟺ `k` lies strictly inside the circumcircle
of `i, j, l`. -/
theorem inside_circumcircle_iff' (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    dist2 (circumcenterBelow a d e) (apexAbove a b c) < dist2 (circumcenterBelow a d e) (0, 0) ↔
      ¬ IsLocallyDelaunay a b c d e := by
  rw [dist2_lt_iff, power_apexAbove hk hl, not_isLocallyDelaunay_iff, pi_lt_iff_cot_neg hk hl]
  have hT := triArea_pos hk
  constructor
  · intro h
    by_contra h'
    push Not at h'
    nlinarith
  · intro h
    nlinarith

theorem isLocallyDelaunay_iff_not_inside_circumcircle' (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    IsLocallyDelaunay a b c d e ↔
      dist2 (circumcenterBelow a d e) (0, 0) ≤ dist2 (circumcenterBelow a d e) (apexAbove a b c) := by
  rw [← not_lt, inside_circumcircle_iff' hk hl, not_not]

/-! ## (b) The flip -/

/-- **An edge with `θ_k + θ_l ≥ π` has a strictly convex diamond** (B–S Prop. 11 for one edge):
the angles of the quadrilateral at `i` and `j` sum to `2π − (θ_k + θ_l) ≤ π` together, and each
is positive. -/
theorem isConvexDiamond_of_pi_le (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : π ≤ cornerAngle a b c + cornerAngle a d e) : IsConvexDiamond a b c d e :=
  (isConvexDiamond_iff_diamondConvex hk hl).mpr (Intrinsic.diamondConvex_of_pi_le hk hl h)

/-- **A non-Delaunay edge has a strictly convex diamond.** -/
theorem isConvexDiamond_of_not_isLocallyDelaunay (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : ¬ IsLocallyDelaunay a b c d e) : IsConvexDiamond a b c d e :=
  isConvexDiamond_of_pi_le hk hl (not_isLocallyDelaunay_iff.mp h).le

/-- **The two triangles a flip creates are nondegenerate**: `ilk` with the sides `(f, b, d)` and
`jkl` with the sides `(f, e, c)`, `f = flipLength a b c d e`. -/
theorem flipLength_isTriangle (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    IsTriangle (flipLength a b c d e) b d ∧ IsTriangle (flipLength a b c d e) e c := by
  rw [flipLength_eq_flippedLengthFromLengths hk hl]
  exact Intrinsic.flip_isTriangle hk hl ((isConvexDiamond_iff_diamondConvex hk hl).mp hconv)

/-- **Flipping back gives the original edge**: in the diamond `(f, c, e, b, d)` of the new edge,
the flipped length is `a`. (The diamond of the edge after two flips is `(a, e, d, c, b)`, the
original one read from `j` with `l` first; see `isLocallyDelaunay_comm`,
`isLocallyDelaunay_mirror`.) -/
theorem flipLength_flipLength (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    flipLength (flipLength a b c d e) c e b d = a := by
  obtain ⟨t1, t2⟩ := flipLength_isTriangle hk hl hconv
  rw [flipLength_eq_flippedLengthFromLengths t2.swap t1, flipLength_eq_flippedLengthFromLengths hk hl]
  exact Intrinsic.flip_back hk hl ((isConvexDiamond_iff_diamondConvex hk hl).mp hconv)

/-- **The diamond of the new edge is strictly convex** (it is the same quadrilateral). -/
theorem flip_isConvexDiamond (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    IsConvexDiamond (flipLength a b c d e) c e b d := by
  obtain ⟨t1, t2⟩ := flipLength_isTriangle hk hl hconv
  rw [isConvexDiamond_iff_diamondConvex t2.swap t1, flipLength_eq_flippedLengthFromLengths hk hl]
  exact Intrinsic.flip_convex hk hl ((isConvexDiamond_iff_diamondConvex hk hl).mp hconv)

/-- **A flip keeps the total area**: `A_ilk + A_jkl = A_ijk + A_jil`. -/
theorem flip_area (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    triArea (flipLength a b c d e) b d + triArea (flipLength a b c d e) e c =
      triArea a b c + triArea a d e := by
  rw [flipLength_eq_flippedLengthFromLengths hk hl]
  exact Intrinsic.flip_area hk hl ((isConvexDiamond_iff_diamondConvex hk hl).mp hconv)

/-- **A flip adds the corner angles at the four vertices of the diamond.** With
`f = flipLength a b c d e`, for a strictly convex diamond:
* at `i`: the corner of `ilk` (opposite `f`) is the sum of the old corners of `ijk` and `jil`;
* at `j`: the corner of `jkl` (opposite `f`) is the sum of the old corners;
* at `k`: the old corner `θ_k` is the sum of the corners of `jkl` (opposite `e`) and of `ilk`
  (opposite `d`);
* at `l`: the old corner `θ_l` is the sum of the corners of `jkl` (opposite `c`) and of `ilk`
  (opposite `b`).

Hence the sum of the corner angles at each of `i, j, k, l` (its cone angle, if all its corners
are counted) is the same before and after the flip.

[formula: it-angle-normalization] -/
theorem flip_corner_angles (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    cornerAngle (flipLength a b c d e) b d = cornerAngle c a b + cornerAngle e a d ∧
      cornerAngle (flipLength a b c d e) c e = cornerAngle b c a + cornerAngle d e a ∧
      cornerAngle e (flipLength a b c d e) c + cornerAngle d (flipLength a b c d e) b =
        cornerAngle a b c ∧
      cornerAngle c e (flipLength a b c d e) + cornerAngle b d (flipLength a b c d e) =
        cornerAngle a d e := by
  rw [flipLength_eq_flippedLengthFromLengths hk hl]
  exact Intrinsic.flip_angle_sums hk hl ((isConvexDiamond_iff_diamondConvex hk hl).mp hconv)

/-- **The opposite angles of the new edge sum to `2π − (θ_k + θ_l)`.** -/
theorem flip_opposite_angles (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    cornerAngle (flipLength a b c d e) c e + cornerAngle (flipLength a b c d e) b d =
      2 * π - (cornerAngle a b c + cornerAngle a d e) := by
  obtain ⟨h1, h2, -, -⟩ := flip_corner_angles hk hl hconv
  have hq := Intrinsic.quad_angle_sum hk hl
  rw [h1, h2]
  linarith

/-- **After the flip of a non-Delaunay edge the new edge is locally Delaunay, strictly**: its
opposite angles sum to less than `π`. -/
theorem flip_isLocallyDelaunay (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : ¬ IsLocallyDelaunay a b c d e) :
    cornerAngle (flipLength a b c d e) c e + cornerAngle (flipLength a b c d e) b d < π := by
  rw [flip_opposite_angles hk hl (isConvexDiamond_of_not_isLocallyDelaunay hk hl h)]
  have := not_isLocallyDelaunay_iff.mp h
  linarith

/-- At a tie both diagonals are locally Delaunay, with equality. -/
theorem flip_of_tie (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle a b c + cornerAngle a d e = π) :
    IsConvexDiamond a b c d e ∧
      cornerAngle (flipLength a b c d e) c e + cornerAngle (flipLength a b c d e) b d = π := by
  have hconv := isConvexDiamond_of_pi_le hk hl h.ge
  refine ⟨hconv, ?_⟩
  rw [flip_opposite_angles hk hl hconv, h]
  ring

/-- **The flip of a non-Delaunay edge, all at once** (the statement a triangulation-level theory
uses). For a diamond `(a, b, c, d, e)` of two nondegenerate triangles with `θ_k + θ_l > π`, and
`f = flipLength a b c d e`:
1. the diamond is strictly convex;
2. the new triangles `(f, b, d)` (`ilk`) and `(f, e, c)` (`jkl`) are nondegenerate;
3. the new edge is locally Delaunay, strictly;
4. flipping the new edge gives back the length `a`;
5. the corner angles at `i, j, k, l` add up as in `flip_corner_angles`;
6. the total area is unchanged.

[formula: it-delaunay-condition] -/
theorem flip_of_not_isLocallyDelaunay (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : ¬ IsLocallyDelaunay a b c d e) :
    IsConvexDiamond a b c d e ∧
      IsTriangle (flipLength a b c d e) b d ∧ IsTriangle (flipLength a b c d e) e c ∧
      cornerAngle (flipLength a b c d e) c e + cornerAngle (flipLength a b c d e) b d < π ∧
      flipLength (flipLength a b c d e) c e b d = a ∧
      (cornerAngle (flipLength a b c d e) b d = cornerAngle c a b + cornerAngle e a d ∧
        cornerAngle (flipLength a b c d e) c e = cornerAngle b c a + cornerAngle d e a ∧
        cornerAngle e (flipLength a b c d e) c + cornerAngle d (flipLength a b c d e) b =
          cornerAngle a b c ∧
        cornerAngle c e (flipLength a b c d e) + cornerAngle b d (flipLength a b c d e) =
          cornerAngle a d e) ∧
      triArea (flipLength a b c d e) b d + triArea (flipLength a b c d e) e c =
        triArea a b c + triArea a d e := by
  have hconv := isConvexDiamond_of_not_isLocallyDelaunay hk hl h
  obtain ⟨t1, t2⟩ := flipLength_isTriangle hk hl hconv
  exact ⟨hconv, t1, t2, flip_isLocallyDelaunay hk hl h, flipLength_flipLength hk hl hconv,
    flip_corner_angles hk hl hconv, flip_area hk hl hconv⟩

/-! ## Coordinates of a strictly convex diamond

The identities of `Potential.lean` and `Rippa.lean` are polynomial identities in the coordinates
of the layout: `k = (X, H)`, `l = (Y, −G)`. -/

/-- **Every length and area of the diamond and of its flip in the coordinates of the layout**:
`k = (X, H)`, `l = (Y, −G)` with `H, G > 0`; the doubled areas of the new triangles are
`X G + Y H > 0` (`ilk`) and `(a − X) G + (a − Y) H > 0` (`jkl`). -/
theorem layout_coords (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    ∃ X H Y G : ℝ, 0 < H ∧ 0 < G ∧ 0 < X * G + Y * H ∧ 0 < (a - X) * G + (a - Y) * H ∧
      b ^ 2 = X ^ 2 + H ^ 2 ∧ c ^ 2 = (a - X) ^ 2 + H ^ 2 ∧ d ^ 2 = Y ^ 2 + G ^ 2 ∧
      e ^ 2 = (a - Y) ^ 2 + G ^ 2 ∧ flipLength a b c d e ^ 2 = (X - Y) ^ 2 + (H + G) ^ 2 ∧
      2 * triArea a b c = a * H ∧ 2 * triArea a d e = a * G ∧
      2 * triArea (flipLength a b c d e) b d = X * G + Y * H ∧
      2 * triArea (flipLength a b c d e) e c = (a - X) * G + (a - Y) * H := by
  have ha := hk.pos_a
  have hconv' := (isConvexDiamond_iff_diamondConvex hk hl).mp hconv
  obtain ⟨h1, h2, -, e2, e3, e4, e5⟩ := Intrinsic.diamond_isometric hk hl
  obtain ⟨c1, c2⟩ := Intrinsic.diamondConvex_sda hconv'
  have t1 := abs_sda_eq_two_mul_triArea (0, 0) (Intrinsic.apexL a d e) (Intrinsic.apexK a b c)
  have t2 := abs_sda_eq_two_mul_triArea (a, 0) (Intrinsic.apexK a b c) (Intrinsic.apexL a d e)
  rw [abs_of_pos c1, Intrinsic.dist2_comm (Intrinsic.apexL a d e) (Intrinsic.apexK a b c),
    Intrinsic.dist2_apex, Intrinsic.dist2_comm (Intrinsic.apexK a b c) (0, 0), e2, e4] at t1
  rw [abs_of_pos c2, Intrinsic.dist2_apex, Intrinsic.dist2_comm (Intrinsic.apexL a d e) (a, 0),
    e5, e3] at t2
  rw [← flipLength_eq_flippedLengthFromLengths hk hl] at t1 t2
  have f1 := Intrinsic.dist2_apex a b c d e
  rw [← flipLength_eq_flippedLengthFromLengths hk hl] at f1
  refine ⟨Ddg.diamondX a b c, Ddg.diamondH a b c, Ddg.diamondX a d e, Ddg.diamondH a d e,
    h1, h2, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · have := c1; unfold sda Intrinsic.apexK Intrinsic.apexL at this; simp only at this; linarith
  · have := c2; unfold sda Intrinsic.apexK Intrinsic.apexL at this; simp only at this; linarith
  · rw [sq_of_dist2 e2]; simp only [Intrinsic.apexK]; ring
  · rw [sq_of_dist2 e3]; simp only [Intrinsic.apexK]; ring
  · rw [sq_of_dist2 e4]; simp only [Intrinsic.apexL]; ring
  · rw [sq_of_dist2 e5]; simp only [Intrinsic.apexL]; ring
  · rw [sq_of_dist2 f1]; simp only [Intrinsic.apexK, Intrinsic.apexL]; ring
  · rw [Intrinsic.diamondH_eq]; field_simp
  · rw [Intrinsic.diamondH_eq]; field_simp
  · rw [← t1]; unfold sda Intrinsic.apexK Intrinsic.apexL; simp only; ring
  · rw [← t2]; unfold sda Intrinsic.apexK Intrinsic.apexL; simp only; ring

end Delaunay

end DDG
