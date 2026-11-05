import DdgProofs.Intrinsic.FlipAngles

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
and `DdgProofs/Intrinsic/FlipAngles.lean`. The statements of the `Delaunay` part of this file mention only `Triangle.lean` definitions and
the plane `ℝ × ℝ` with `dist2` (`DdgProofs/Bff/Quality.lean`);
`flipLength_eq_flippedLengthFromLengths` relates `flipLength` to the flipped length
`Intrinsic.flippedLengthFromLengths` of the Rust code. The section on the intrinsic criterion
(namespace `DDG.Intrinsic`, moved here from `DdgProofs/Intrinsic/Delaunay.lean`) defines the
closed forms of the Rust code on `ℝ` (`delaunayExcess`, `delaunayHolds`, `cotanWeight`,
`cotanWeightFromLengths`, `oppositeAngles`); the oracle's copies are related to them by
one-way bridges in the app layer (`Oracle/Intrinsic.lean`).
-/

noncomputable section

namespace DDG

/-! ## The intrinsic Delaunay criterion of one diamond

(Moved here from `DdgProofs/Intrinsic/Delaunay.lean` in library step S2b; the names keep the
`DDG.Intrinsic` namespace until step S4.)

* **`it-delaunay-condition`**, **`it-cotan-weight`** — `cotanWeight` (the Rust function with
  this tag) is the sum `cot θ_k + cot θ_l`: twice the cotan weight
  `w = ½ (cot θ_k + cot θ_l)` of the edge (SPEC N6, N17), and minus twice the entry of the cotan
  Laplacian, `L_ij = −w` (`offDiag_eq_neg_half_cotanWeight`, `cotanWeight_eq_neg_two_mul_offDiag`).
  For two angles in `(0, π)`: `cot θ_k + cot θ_l = sin(θ_k + θ_l)/(sin θ_k sin θ_l)`, so the sum,
  and with it `w`, is `≥ 0`, `> 0`, `= 0` exactly when `θ_k + θ_l` is `≤ π`, `< π`, `= π`
  (`cotanWeight_nonneg_iff`, `cotanWeight_pos_iff`, `cotanWeight_eq_zero_iff`,
  `half_cotanWeight_sign`); the same from the five lengths
  (`delaunay_iff_cotanWeight_nonneg`). An edge that violates the condition
  has a strictly convex diamond (`diamondConvex_of_not_delaunay`); two corners of one triangle
  never violate it (`delaunay_of_same_triangle`). A flip of a strictly convex diamond negates
  the excess `θ_k + θ_l − π` (`flip_delaunayExcess`), so the flip of a violating edge gives an
  edge that satisfies the condition strictly (`flip_of_not_delaunay`).

Termination of the flipping algorithm, existence and uniqueness of the intrinsic Delaunay
triangulation and Rippa's theorem are cited (`CITED.md`).
-/

namespace Intrinsic

open Real

/-! ### The formulas -/

/-- Excess of the opposite angles over `π`: `θ_k + θ_l − π`; the edge is locally Delaunay iff it
is `≤ 0`. -/
def delaunayExcess (thetaK thetaL : ℝ) : ℝ := thetaK + thetaL - π

/-- The Delaunay condition with a tolerance: `θ_k + θ_l − π ≤ tol` (the test of the Rust code, a
`Bool`). -/
def delaunayHolds (tol thetaK thetaL : ℝ) : Bool := decide (delaunayExcess thetaK thetaL ≤ tol)

/-- The sum `cot θ_k + cot θ_l` of the cotangents of the two angles opposite an edge: twice the
cotan weight `w = ½ (cot θ_k + cot θ_l)` of SPEC N6, N17, and minus twice the off-diagonal entry
of the (positive semidefinite) cotan Laplacian. -/
def cotanWeight (cotK cotL : ℝ) : ℝ := cotK + cotL

/-- The cotan weight from the five lengths of the diamond. -/
def cotanWeightFromLengths (a b c d e : ℝ) : ℝ := cotanWeight (cornerCot a b c) (cornerCot a d e)

/-- The opposite angles from the five lengths: `(θ_k, θ_l)`. -/
def oppositeAngles (a b c d e : ℝ) : ℝ × ℝ := (cornerAngle a b c, cornerAngle a d e)

theorem delaunayExcess_eq (θk θl : ℝ) : delaunayExcess θk θl = θk + θl - π := rfl

theorem delaunayHolds_iff (tol θk θl : ℝ) :
    delaunayHolds tol θk θl = true ↔ delaunayExcess θk θl ≤ tol := by
  simp only [delaunayHolds, decide_eq_true_eq]

theorem cotanWeight_eq (x y : ℝ) : cotanWeight x y = x + y := rfl

theorem cotanWeightFromLengths_eq (a b c d e : ℝ) :
    cotanWeightFromLengths a b c d e = cornerCot a b c + cornerCot a d e := rfl

theorem oppositeAngles_eq (a b c d e : ℝ) :
    oppositeAngles a b c d e = (cornerAngle a b c, cornerAngle a d e) := rfl

/-- The order in which the Rust code evaluates the excess, `(θ_k − π) + θ_l`, is the same
number. -/
theorem delaunayExcess_eq' (θk θl : ℝ) : delaunayExcess θk θl = (θk - π) + θl := by
  rw [delaunayExcess_eq]
  ring

/-! ### Cotangents of two angles -/

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

/-- **The sum of the two cotangents has the sign of `π − (θ_k + θ_l)`**: for angles in `(0, π)`
`cotanWeight` (`cot θ_k + cot θ_l`, twice the weight `w` of SPEC N6) is
`sin(θ_k + θ_l)/(sin θ_k sin θ_l)` with a positive denominator.

[formula: it-cotan-weight] -/
theorem cotanWeight_eq_sin_div {θk θl : ℝ} (hk0 : 0 < θk) (hk : θk < π) (hl0 : 0 < θl)
    (hl : θl < π) :
    cotanWeight (cot θk) (cot θl) = sin (θk + θl) / (sin θk * sin θl) ∧
      0 < sin θk * sin θl := by
  have s1 := Real.sin_pos_of_pos_of_lt_pi hk0 hk
  have s2 := Real.sin_pos_of_pos_of_lt_pi hl0 hl
  exact ⟨by rw [cotanWeight_eq, cot_add_cot s1.ne' s2.ne'], mul_pos s1 s2⟩

/-- **The Delaunay condition is the nonnegativity of the cotan weight**: for two angles in
`(0, π)`, `cot θ_k + cot θ_l ≥ 0 ⟺ θ_k + θ_l ≤ π`. (`cotanWeight` is `2w`; for the
weight `w` itself see `half_cotanWeight_sign`.)

[formula: it-cotan-weight] -/
theorem cotanWeight_nonneg_iff {θk θl : ℝ} (hk0 : 0 < θk) (hk : θk < π) (hl0 : 0 < θl)
    (hl : θl < π) :
    0 ≤ cotanWeight (cot θk) (cot θl) ↔ delaunayExcess θk θl ≤ 0 := by
  obtain ⟨e, hpos⟩ := cotanWeight_eq_sin_div hk0 hk hl0 hl
  rw [e, delaunayExcess_eq, div_nonneg_iff_of_pos hpos,
    sin_nonneg_iff_le_pi (by linarith) (by linarith)]
  constructor <;> intro h <;> linarith

/-- The strict form: `cot θ_k + cot θ_l > 0 ⟺ θ_k + θ_l < π`.

[formula: it-cotan-weight] -/
theorem cotanWeight_pos_iff {θk θl : ℝ} (hk0 : 0 < θk) (hk : θk < π) (hl0 : 0 < θl)
    (hl : θl < π) :
    0 < cotanWeight (cot θk) (cot θl) ↔ delaunayExcess θk θl < 0 := by
  obtain ⟨e, hpos⟩ := cotanWeight_eq_sin_div hk0 hk hl0 hl
  rw [e, delaunayExcess_eq, div_pos_iff_of_pos_right hpos,
    sin_pos_iff_lt_pi (by linarith) (by linarith)]
  constructor <;> intro h <;> linarith

/-- The weight vanishes exactly at a tie, `θ_k + θ_l = π` (the four points are cocircular). -/
theorem cotanWeight_eq_zero_iff {θk θl : ℝ} (hk0 : 0 < θk) (hk : θk < π) (hl0 : 0 < θl)
    (hl : θl < π) :
    cotanWeight (cot θk) (cot θl) = 0 ↔ delaunayExcess θk θl = 0 := by
  have h1 := cotanWeight_nonneg_iff hk0 hk hl0 hl
  have h2 := cotanWeight_pos_iff hk0 hk hl0 hl
  constructor
  · intro h
    have a1 := h1.mp h.ge
    have a2 : ¬ delaunayExcess θk θl < 0 := fun hlt => (h2.mpr hlt).ne' h
    linarith [not_lt.mp a2]
  · intro h
    have a1 := h1.mpr h.le
    have a2 : ¬ 0 < cotanWeight (cot θk) (cot θl) := fun hpos => (h2.mp hpos).ne h
    linarith [not_lt.mp a2]

/-- **The weight `w = ½ (cot θ_k + cot θ_l)` of SPEC N6 is half `cotanWeight` and
has the same sign**: `w ≥ 0 ⟺ θ_k + θ_l ≤ π`, `w > 0 ⟺ θ_k + θ_l < π`,
`w = 0 ⟺ θ_k + θ_l = π`, for two angles in `(0, π)`.

[formula: it-cotan-weight] -/
theorem half_cotanWeight_sign {θk θl : ℝ} (hk0 : 0 < θk) (hk : θk < π) (hl0 : 0 < θl)
    (hl : θl < π) :
    (0 ≤ cotanWeight (cot θk) (cot θl) / 2 ↔ delaunayExcess θk θl ≤ 0) ∧
      (0 < cotanWeight (cot θk) (cot θl) / 2 ↔ delaunayExcess θk θl < 0) ∧
      (cotanWeight (cot θk) (cot θl) / 2 = 0 ↔ delaunayExcess θk θl = 0) := by
  rw [← cotanWeight_nonneg_iff hk0 hk hl0 hl, ← cotanWeight_pos_iff hk0 hk hl0 hl,
    ← cotanWeight_eq_zero_iff hk0 hk hl0 hl]
  refine ⟨?_, ?_, ?_⟩ <;> constructor <;> intro h <;> linarith

/-! ### From the five lengths -/

/-- The weight from the lengths is the sum of the cotangents of the two opposite
angles. -/
theorem cotanWeightFromLengths_eq_cot {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    cotanWeightFromLengths a b c d e =
      cotanWeight (cot (cornerAngle a b c)) (cot (cornerAngle a d e)) := by
  rw [cotanWeightFromLengths_eq, cotanWeight_eq, cornerCot_eq_cot hk, cornerCot_eq_cot hl]

/-- **The Delaunay condition from the five lengths**: the test `delaunayHolds` with tolerance
`0` on `oppositeAngles` holds exactly when `cotanWeightFromLengths` is nonnegative.

[formula: it-delaunay-condition] -/
theorem delaunay_iff_cotanWeight_nonneg {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    delaunayHolds 0 (oppositeAngles a b c d e).1 (oppositeAngles a b c d e).2 =
        true ↔ 0 ≤ cotanWeightFromLengths a b c d e := by
  rw [delaunayHolds_iff, oppositeAngles_eq, cotanWeightFromLengths_eq_cot hk hl,
    cotanWeight_nonneg_iff (cornerAngle_pos hk) (cornerAngle_lt_pi hk) (cornerAngle_pos hl)
      (cornerAngle_lt_pi hl)]

/-- The same as a statement about the angles: `θ_k + θ_l ≤ π ⟺` weight `≥ 0`.

[formula: it-delaunay-condition] -/
theorem angle_sum_le_pi_iff_cotanWeight_nonneg {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    cornerAngle a b c + cornerAngle a d e ≤ π ↔ 0 ≤ cotanWeightFromLengths a b c d e := by
  rw [cotanWeightFromLengths_eq_cot hk hl,
    cotanWeight_nonneg_iff (cornerAngle_pos hk) (cornerAngle_lt_pi hk) (cornerAngle_pos hl)
      (cornerAngle_lt_pi hl), delaunayExcess_eq]
  constructor <;> intro h <;> linarith

/-- An edge violates the condition exactly when its cotan weight is negative. -/
theorem pi_lt_angle_sum_iff_cotanWeight_neg {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    π < cornerAngle a b c + cornerAngle a d e ↔ cotanWeightFromLengths a b c d e < 0 := by
  rw [← not_le, angle_sum_le_pi_iff_cotanWeight_nonneg hk hl, not_le]

/-- **Two corners of one triangle sum to less than `π`**: an edge with the same triangle on both
sides (an endpoint of degree one) satisfies the condition strictly (the paper's remark in
Sec. 4.1). -/
theorem delaunay_of_same_triangle {a b c : ℝ} (h : IsTriangle a b c) :
    delaunayExcess (cornerAngle a b c) (cornerAngle b c a) < 0 := by
  have h1 := cornerAngle_add_cornerAngle_add_cornerAngle h
  have h2 := cornerAngle_pos h.rotate.rotate
  rw [delaunayExcess_eq]
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
    (hl : IsTriangle a d e) (h : cotanWeightFromLengths a b c d e < 0) :
    DiamondConvex a b c d e :=
  diamondConvex_of_not_delaunay hk hl ((pi_lt_angle_sum_iff_cotanWeight_neg hk hl).mpr h)

/-! ### The flip of a violating edge -/

/-- **A flip negates the Delaunay excess**: for a strictly convex diamond, the opposite angles
of the new edge `kl` (at `j` and at `i`) sum to `2π − (θ_k + θ_l)`.

[formula: it-delaunay-condition] -/
theorem flip_delaunayExcess {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) :
    delaunayExcess (cornerAngle (flippedLengthFromLengths a b c d e) c e)
        (cornerAngle (flippedLengthFromLengths a b c d e) b d) =
      -delaunayExcess (cornerAngle a b c) (cornerAngle a d e) := by
  obtain ⟨h1, h2, -, -⟩ := flip_angle_sums hk hl hconv
  have hq := quad_angle_sum hk hl
  rw [delaunayExcess_eq, delaunayExcess_eq, h1, h2]
  linarith

/-- **After the flip of a violating edge the new edge satisfies the condition strictly**: if
`θ_k + θ_l > π`, the diamond is strictly convex, the two new triangles satisfy the strict
triangle inequalities, the excess of the new diamond `(f, c, e, b, d)` is negative and
its cotan weight positive.

[formula: it-delaunay-condition] -/
theorem flip_of_not_delaunay {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : π < cornerAngle a b c + cornerAngle a d e) :
    DiamondConvex a b c d e ∧
      delaunayExcess
          (oppositeAngles (flippedLengthFromLengths a b c d e) c e b d).1
          (oppositeAngles (flippedLengthFromLengths a b c d e) c e b d).2 < 0 ∧
      0 < cotanWeightFromLengths (flippedLengthFromLengths a b c d e) c e b d := by
  have hconv := diamondConvex_of_not_delaunay hk hl h
  obtain ⟨t1, t2⟩ := flip_isTriangle hk hl hconv
  have hk' := t2.swap
  have hx := flip_delaunayExcess hk hl hconv
  have hneg : delaunayExcess
      (cornerAngle (flippedLengthFromLengths a b c d e) c e)
      (cornerAngle (flippedLengthFromLengths a b c d e) b d) < 0 := by
    rw [hx, delaunayExcess_eq]
    linarith
  refine ⟨hconv, ?_, ?_⟩
  · rw [oppositeAngles_eq]
    exact hneg
  · rw [cotanWeightFromLengths_eq_cot hk' t1]
    exact (cotanWeight_pos_iff (cornerAngle_pos hk') (cornerAngle_lt_pi hk') (cornerAngle_pos t1)
      (cornerAngle_lt_pi t1)).mpr hneg

/-- At a tie (`θ_k + θ_l = π`) both diagonals satisfy the condition with equality: the flip
gives the excess `0` again (SPEC N6). -/
theorem flip_of_tie {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle a b c + cornerAngle a d e = π) :
    DiamondConvex a b c d e ∧
      delaunayExcess (cornerAngle (flippedLengthFromLengths a b c d e) c e)
        (cornerAngle (flippedLengthFromLengths a b c d e) b d) = 0 := by
  have hconv := diamondConvex_of_pi_le hk hl h.ge
  refine ⟨hconv, ?_⟩
  rw [flip_delaunayExcess hk hl hconv, delaunayExcess_eq]
  linarith

/-! ### The entry of the cotan Laplacian -/

section Laplacian

variable {V F : Type*} [DecidableEq V] [Fintype F]

/-- **The entry of the cotan Laplacian is minus the cotan weight `w`, i.e. minus half
`cotanWeight`** (`cotan-laplacian`, `DdgProofs/Cotan.lean`; `L` positive semidefinite):
if only the two faces `f₁ ≠ f₂` of a mesh have the edge `{i, j}`, then
`L_ij = −½ (cot θ_k + cot θ_l) = −w`, with the two cotangents those of the corners opposite the
edge.

[formula: it-cotan-weight] -/
theorem offDiag_eq_neg_half_cotanWeight (M : TriMesh V F) {i j : V} {f₁ f₂ : F} (hne : f₁ ≠ f₂)
    (h : ∀ f, f ≠ f₁ → f ≠ f₂ → M.cotOpp f i j = 0) :
    M.offDiag i j = -(cotanWeight (M.cotOpp f₁ i j) (M.cotOpp f₂ i j)) / 2 := by
  unfold TriMesh.offDiag
  rw [cotanWeight_eq, Finset.sum_eq_add f₁ f₂ hne (fun f _ hf => h f hf.1 hf.2)
    (fun hf => absurd (Finset.mem_univ _) hf) (fun hf => absurd (Finset.mem_univ _) hf)]

/-- **The function tagged `it-cotan-weight` returns `2w = −2 L_ij`**: the sum of the two
cotangents is twice the weight of SPEC N6 and minus twice the entry of the cotan Laplacian.

[formula: it-cotan-weight] -/
theorem cotanWeight_eq_neg_two_mul_offDiag (M : TriMesh V F) {i j : V} {f₁ f₂ : F}
    (hne : f₁ ≠ f₂) (h : ∀ f, f ≠ f₁ → f ≠ f₂ → M.cotOpp f i j = 0) :
    cotanWeight (M.cotOpp f₁ i j) (M.cotOpp f₂ i j) = -2 * M.offDiag i j := by
  rw [offDiag_eq_neg_half_cotanWeight M hne h]
  ring

/-- **The off-diagonal entry is nonpositive exactly when the cotan weight of the edge is
nonnegative**, that is (by `cotanWeight_nonneg_iff`) when the edge satisfies the Delaunay
condition: the Laplacian of a Delaunay triangulation has no positive off-diagonal entry
(SPEC §5 item 4, N6). -/
theorem offDiag_nonpos_iff (M : TriMesh V F) {i j : V} {f₁ f₂ : F} (hne : f₁ ≠ f₂)
    (h : ∀ f, f ≠ f₁ → f ≠ f₂ → M.cotOpp f i j = 0) :
    M.offDiag i j ≤ 0 ↔ 0 ≤ cotanWeight (M.cotOpp f₁ i j) (M.cotOpp f₂ i j) := by
  rw [offDiag_eq_neg_half_cotanWeight M hne h]
  constructor <;> intro h' <;> linarith

end Laplacian

end Intrinsic

namespace Delaunay

open Real

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

/-- `flipLength` is the flipped length of the Rust code.

[formula: it-flip-length] -/
theorem flipLength_eq_flippedLengthFromLengths (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    flipLength a b c d e = Intrinsic.flippedLengthFromLengths a b c d e := by
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
  refine ⟨Intrinsic.diamondX a b c, Intrinsic.diamondH a b c, Intrinsic.diamondX a d e, Intrinsic.diamondH a d e,
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
