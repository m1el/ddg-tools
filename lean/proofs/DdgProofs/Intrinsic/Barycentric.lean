import DdgProofs.Intrinsic.Diamond
import DdgProofs.Cuts.Crossing

/-!
# Intrinsic triangulations: barycentric coordinates, distances to the corners, the circumcentre

`crates/ddg-intrinsic/SPEC.md` §4 (Insert, Point query, Delaunay refinement), N10, N14; the paper
Sec. 3.3.2, Sec. 4.2, Alg. 6, 7, 8.

A triangle with the sides `ℓ₀₁, ℓ₁₂, ℓ₂₀` is laid out as `v₀ = (0, 0)`, `v₁ = (ℓ₀₁, 0)`,
`v₂ = apexK ℓ₀₁ ℓ₂₀ ℓ₁₂` (`DdgProofs/Intrinsic/Diamond.lean`).

* **`it-barycentric-coordinates`** — the oracle's `barycentric` of a plane triangle of nonzero
  area: the three ratios of signed areas sum to one, reproduce the point, and are the only
  numbers that do (`barycentric_sum`, `barycentric_reproduces`, `barycentric_unique`); the form
  of the Rust code (normalized by the sum of the three areas) is the same
  (`barycentric_eq_normalized`). `pointFromDistances` returns the point from its distances to
  the three corners (`pointFromDistances_eq_self`).
* **`it-point-distances`** — the oracle's `cornerDistances` are the plane distances of the
  point `b₁ v₁ + b₂ v₂` to the corners (`cornerDistances_eq`), and for `b₀ + b₁ + b₂ = 1` they
  are the oracle's `baryLength` of the three sides and the coordinate differences: the formula
  of the paper's Alg. 6, by the theorem of `cut-barycentric-length`
  (`cornerDistances_eq_baryLength`). On a side (`b₂ = 0`) the distances to its two ends are
  `ℓ₀₁ b₁` and `ℓ₀₁ b₀`, the oracle's `edgeSplitLengths`, and they sum to `ℓ₀₁`
  (`edgeSplitLengths_eq_cornerDistances`, `edgeSplitLengths_sum`). With `b_i` the weight of
  `i`, the distance from `i` is `ℓ_ij b_j`: the paper's Alg. 8 prints `ℓ_ip ← ℓ_ij b_i`, which
  contradicts its Alg. 6 (`alg6_distance_on_edge`, `alg8_as_printed_differs`).
* **`it-circumcenter`** — the oracle's `circumcenterBary`: the weights sum to `16 A²`, the
  coordinates to one (`circumWeight_sum`, `circumcenterBary_sum`); the point has the distance
  `a₀a₁a₂/(4A)` to all three corners (`circumcenter_equidistant`); the coordinate of a corner
  is negative exactly when the angle there is obtuse, zero exactly when it is right, and at
  most one is negative (`circumcenterBary_neg_iff`, `circumcenterBary_eq_zero_iff`,
  `circumcenterBary_pos_of_neg`).

Own definition: `layoutPoint` (the point `b₁ v₁ + b₂ v₂` of the layout).
-/

noncomputable section

namespace DdgProofs

namespace Intrinsic

open scoped DdgProofs.Oracle
open Real

/-! ## Barycentric coordinates in a plane triangle -/

theorem signedArea2_eq_sda (a b c : ℝ × ℝ) :
    Ddg.signedArea2 a.1 a.2 b.1 b.2 c.1 c.2 = sda a b c := by
  unfold Ddg.signedArea2 sda
  ring

theorem barycentric_eq (a b c p : ℝ × ℝ) :
    Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2 =
      (sda p b c / sda a b c, sda a p c / sda a b c, sda a b p / sda a b c) := by
  simp only [Ddg.barycentric, signedArea2_eq_sda]

/-- The three triangles that a point makes with the sides add up to the whole (signed). -/
theorem sda_split (a b c p : ℝ × ℝ) : sda p b c + sda a p c + sda a b p = sda a b c := by
  unfold sda
  ring

/-- **The barycentric coordinates sum to one.**

[formula: it-barycentric-coordinates] -/
theorem barycentric_sum {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) (p : ℝ × ℝ) :
    (Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2).1 +
        (Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2).2.1 +
        (Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2).2.2 = 1 := by
  rw [barycentric_eq]
  dsimp only
  rw [← add_div, ← add_div, sda_split, div_self h]

/-- **They reproduce the point**: `b₀ a + b₁ b + b₂ c = p`.

[formula: it-barycentric-coordinates] -/
theorem barycentric_reproduces {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) (p : ℝ × ℝ) :
    ((Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2).1 * a.1 +
        (Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2).2.1 * b.1 +
        (Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2).2.2 * c.1,
      (Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2).1 * a.2 +
        (Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2).2.1 * b.2 +
        (Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2).2.2 * c.2) = p := by
  have k1 : sda p b c * a.1 + sda a p c * b.1 + sda a b p * c.1 = p.1 * sda a b c := by
    unfold sda
    ring
  have k2 : sda p b c * a.2 + sda a p c * b.2 + sda a b p * c.2 = p.2 * sda a b c := by
    unfold sda
    ring
  rw [barycentric_eq]
  ext
  · dsimp only
    rw [div_mul_eq_mul_div, div_mul_eq_mul_div, div_mul_eq_mul_div, ← add_div, ← add_div, k1,
      mul_div_assoc, div_self h, mul_one]
  · dsimp only
    rw [div_mul_eq_mul_div, div_mul_eq_mul_div, div_mul_eq_mul_div, ← add_div, ← add_div, k2,
      mul_div_assoc, div_self h, mul_one]

/-- **They are the only such numbers**: if `u + v + w = 1` and `u a + v b + w c = p` then
`(u, v, w)` are the oracle's coordinates. -/
theorem barycentric_unique {a b c p : ℝ × ℝ} (h : sda a b c ≠ 0) {u v w : ℝ}
    (hs : u + v + w = 1) (hx : u * a.1 + v * b.1 + w * c.1 = p.1)
    (hy : u * a.2 + v * b.2 + w * c.2 = p.2) :
    Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2 = (u, v, w) := by
  have hw : w = 1 - u - v := by linarith
  have e0 : sda p b c = u * sda a b c := by
    unfold sda
    rw [← hx, ← hy, hw]
    ring
  have e1 : sda a p c = v * sda a b c := by
    unfold sda
    rw [← hx, ← hy, hw]
    ring
  have e2 : sda a b p = w * sda a b c := by
    unfold sda
    rw [← hx, ← hy, hw]
    ring
  rw [barycentric_eq, e0, e1, e2, mul_div_assoc, mul_div_assoc, mul_div_assoc, div_self h,
    mul_one, mul_one, mul_one]

/-- The form of the Rust code: the three areas `[p, b, c]`, `[p, c, a]`, `[p, a, b]` divided
by their sum. In exact arithmetic it is the oracle's. -/
theorem barycentric_eq_normalized (a b c p : ℝ × ℝ) :
    (sda p b c / (sda p b c + sda p c a + sda p a b),
        sda p c a / (sda p b c + sda p c a + sda p a b),
        sda p a b / (sda p b c + sda p c a + sda p a b)) =
      Ddg.barycentric a.1 a.2 b.1 b.2 c.1 c.2 p.1 p.2 := by
  have e1 : sda p c a = sda a p c := by
    unfold sda
    ring
  have e2 : sda p a b = sda a b p := by
    unfold sda
    ring
  rw [barycentric_eq, e1, e2, sda_split]

/-- **A point from its distances to the three corners**: in the layout of the triangle
`(ℓ_ab, ℓ_bc, ℓ_ca)` the oracle's `pointFromDistances` of the distances of a point `p` is `p`.

[formula: it-barycentric-coordinates] -/
theorem pointFromDistances_eq_self {lab lbc lca : ℝ} (h : IsTriangle lab lbc lca) (p : ℝ × ℝ) :
    Ddg.pointFromDistances lab lbc lca (dist2 (0, 0) p) (dist2 (lab, 0) p)
      (dist2 (apexK lab lca lbc) p) = p := by
  have hl := h.pos_a
  obtain ⟨hcy, -, -, -, -, -, -⟩ := diamond_isometric h.swap h.swap
  have hra : dist2 (0, 0) p ^ 2 = p.1 ^ 2 + p.2 ^ 2 := by
    rw [dist2_sq']
    dsimp only
    ring
  have hrb : dist2 (lab, 0) p ^ 2 = (p.1 - lab) ^ 2 + p.2 ^ 2 := by
    rw [dist2_sq']
    dsimp only
    ring
  have hrc : dist2 (apexK lab lca lbc) p ^ 2 =
      (p.1 - Ddg.diamondX lab lca lbc) ^ 2 + (p.2 - Ddg.diamondH lab lca lbc) ^ 2 := by
    rw [dist2_sq']
    rfl
  have hx : (lab * lab + dist2 (0, 0) p * dist2 (0, 0) p -
      dist2 (lab, 0) p * dist2 (lab, 0) p) / (2 * lab) = p.1 := by
    rw [div_eq_iff (by positivity)]
    linear_combination hra - hrb
  have hy : (Ddg.diamondX lab lca lbc * Ddg.diamondX lab lca lbc +
      Ddg.diamondH lab lca lbc * Ddg.diamondH lab lca lbc +
      dist2 (0, 0) p * dist2 (0, 0) p -
      dist2 (apexK lab lca lbc) p * dist2 (apexK lab lca lbc) p -
      2 * Ddg.diamondX lab lca lbc * ((lab * lab + dist2 (0, 0) p * dist2 (0, 0) p -
        dist2 (lab, 0) p * dist2 (lab, 0) p) / (2 * lab))) /
      (2 * Ddg.diamondH lab lca lbc) = p.2 := by
    rw [hx, div_eq_iff (by positivity)]
    linear_combination hra - hrc
  simp only [Ddg.pointFromDistances, Oracle.two_eq]
  exact Prod.ext hx hy

/-! ## Distances of a point to the corners -/

/-- The point with the barycentric coordinates `(·, b₁, b₂)` in the layout: `b₁ v₁ + b₂ v₂`. -/
def layoutPoint (l01 l12 l20 b1 b2 : ℝ) : ℝ × ℝ :=
  (b1 * l01 + b2 * Ddg.diamondX l01 l20 l12, b2 * Ddg.diamondH l01 l20 l12)

/-- **The oracle's `cornerDistances` are the plane distances of the point to the three corners
of the layout.** (They do not depend on `b₀`.)

[formula: it-point-distances] -/
theorem cornerDistances_eq (l01 l12 l20 b0 b1 b2 : ℝ) :
    Ddg.cornerDistances l01 l12 l20 b0 b1 b2 =
      (dist2 (0, 0) (layoutPoint l01 l12 l20 b1 b2),
        dist2 (l01, 0) (layoutPoint l01 l12 l20 b1 b2),
        dist2 (apexK l01 l20 l12) (layoutPoint l01 l12 l20 b1 b2)) := by
  unfold dist2 layoutPoint apexK
  simp only [Ddg.cornerDistances, Ddg.planeDistance, Oracle.sqrt_eq, Oracle.zero_eq]
  refine Prod.ext ?_ (Prod.ext ?_ ?_) <;>
  · dsimp only
    congr 1
    ring

/-- The plane distance as a distance of complex numbers. -/
theorem dist2_eq_dist_toC (p q : ℝ × ℝ) : dist2 p q = dist (toC p) (toC q) := by
  rw [Complex.dist_eq_re_im]
  unfold dist2 toC
  congr 1
  dsimp only
  ring

/-- **The barycentric length formula in the plane**, from the theorem of
`cut-barycentric-length` (`baryLength_eq_dist`): for plane points `v₀, v₁, v₂` and coordinates
`x`, `y` with sum one, the distance of `Σ y_k v_k` and `Σ x_k v_k` is the oracle's `baryLength`
of the three distances and `x − y`. -/
theorem dist2_bary_eq_baryLength (v0 v1 v2 : ℝ × ℝ) {x0 x1 x2 y0 y1 y2 : ℝ}
    (hx : x0 + x1 + x2 = 1) (hy : y0 + y1 + y2 = 1) :
    dist2 (y0 * v0.1 + y1 * v1.1 + y2 * v2.1, y0 * v0.2 + y1 * v1.2 + y2 * v2.2)
        (x0 * v0.1 + x1 * v1.1 + x2 * v2.1, x0 * v0.2 + x1 * v1.2 + x2 * v2.2) =
      Ddg.baryLength (dist2 v0 v1) (dist2 v1 v2) (dist2 v2 v0) (x0 - y0) (x1 - y1) (x2 - y2) := by
  have key := Cuts.baryLength_eq_dist (toC v0) (toC v1) (toC v2) hx hy
  rw [dist2_eq_dist_toC v0 v1, dist2_eq_dist_toC v1 v2, dist2_eq_dist_toC v2 v0, key,
    dist2_eq_dist_toC, dist_comm]
  congr 1 <;>
  · apply Complex.ext <;> simp [toC]

/-- **The distances to the corners are the formula of the paper's Alg. 6**: for a triangle and
`b₀ + b₁ + b₂ = 1`, the oracle's `cornerDistances` are its `baryLength` (of
`cut-barycentric-length`) for the differences `b − e_k`.

[formula: it-point-distances] -/
theorem cornerDistances_eq_baryLength {l01 l12 l20 : ℝ} (h : IsTriangle l01 l12 l20)
    {b0 b1 b2 : ℝ} (hb : b0 + b1 + b2 = 1) :
    Ddg.cornerDistances l01 l12 l20 b0 b1 b2 =
      (Ddg.baryLength l01 l12 l20 (b0 - 1) b1 b2, Ddg.baryLength l01 l12 l20 b0 (b1 - 1) b2,
        Ddg.baryLength l01 l12 l20 b0 b1 (b2 - 1)) := by
  obtain ⟨-, -, e1, e2, e3, -, -⟩ := diamond_isometric h.swap h.swap
  have e2' : dist2 (apexK l01 l20 l12) (0, 0) = l20 := by rw [dist2_comm]; exact e2
  have k0 := dist2_bary_eq_baryLength (0, 0) (l01, 0) (apexK l01 l20 l12) hb
    (by norm_num : (1 : ℝ) + 0 + 0 = 1)
  have k1 := dist2_bary_eq_baryLength (0, 0) (l01, 0) (apexK l01 l20 l12) hb
    (by norm_num : (0 : ℝ) + 1 + 0 = 1)
  have k2 := dist2_bary_eq_baryLength (0, 0) (l01, 0) (apexK l01 l20 l12) hb
    (by norm_num : (0 : ℝ) + 0 + 1 = 1)
  rw [e1, e3, e2'] at k0 k1 k2
  have eP : (b0 * ((0 : ℝ), (0 : ℝ)).1 + b1 * (l01, (0 : ℝ)).1 + b2 * (apexK l01 l20 l12).1,
      b0 * ((0 : ℝ), (0 : ℝ)).2 + b1 * (l01, (0 : ℝ)).2 + b2 * (apexK l01 l20 l12).2) =
      layoutPoint l01 l12 l20 b1 b2 := by
    unfold layoutPoint apexK
    ext <;> simp
  rw [eP] at k0 k1 k2
  have q0 : (1 * ((0 : ℝ), (0 : ℝ)).1 + 0 * (l01, (0 : ℝ)).1 + 0 * (apexK l01 l20 l12).1,
      1 * ((0 : ℝ), (0 : ℝ)).2 + 0 * (l01, (0 : ℝ)).2 + 0 * (apexK l01 l20 l12).2) =
      ((0 : ℝ), (0 : ℝ)) := by
    ext <;> simp
  have q1 : (0 * ((0 : ℝ), (0 : ℝ)).1 + 1 * (l01, (0 : ℝ)).1 + 0 * (apexK l01 l20 l12).1,
      0 * ((0 : ℝ), (0 : ℝ)).2 + 1 * (l01, (0 : ℝ)).2 + 0 * (apexK l01 l20 l12).2) =
      (l01, (0 : ℝ)) := by
    ext <;> simp
  have q2 : (0 * ((0 : ℝ), (0 : ℝ)).1 + 0 * (l01, (0 : ℝ)).1 + 1 * (apexK l01 l20 l12).1,
      0 * ((0 : ℝ), (0 : ℝ)).2 + 0 * (l01, (0 : ℝ)).2 + 1 * (apexK l01 l20 l12).2) =
      apexK l01 l20 l12 := by
    ext <;> simp
  rw [q0, sub_zero, sub_zero] at k0
  rw [q1, sub_zero, sub_zero] at k1
  rw [q2, sub_zero, sub_zero] at k2
  rw [cornerDistances_eq, k0, k1, k2]

/-! ## A point on a side: the lengths of an edge split -/

theorem edgeSplitLengths_eq (lij bi bj : ℝ) :
    Ddg.edgeSplitLengths lij bi bj = (lij * bj, lij * bi) := rfl

/-- **The two lengths of an edge split sum to the length of the edge** (`b_i + b_j = 1`).

[formula: it-point-distances] -/
theorem edgeSplitLengths_sum {bi bj : ℝ} (h : bi + bj = 1) (lij : ℝ) :
    (Ddg.edgeSplitLengths lij bi bj).1 + (Ddg.edgeSplitLengths lij bi bj).2 = lij := by
  rw [edgeSplitLengths_eq]
  dsimp only
  rw [← mul_add, add_comm, h, mul_one]

/-- **The paper's Alg. 6 on an edge**: the distance from the corner `i` to the point with the
barycentric weights `(b_i, b_j, 0)` is `ℓ_ij b_j` (`b_i + b_j = 1`, `ℓ_ij, b_j ≥ 0`): the weight
of the *other* end times the length. -/
theorem alg6_distance_on_edge {lij : ℝ} (hl : 0 ≤ lij) {bi bj : ℝ} (hb : bi + bj = 1)
    (hj : 0 ≤ bj) (ljk lki : ℝ) :
    Ddg.baryLength lij ljk lki (bi - 1) bj 0 = lij * bj := by
  have e : Ddg.baryLength2 lij ljk lki (bi - 1) bj 0 = (lij * bj) ^ 2 := by
    rw [Cuts.baryLength2_eq, show bi - 1 = -bj by linarith]
    ring
  rw [Cuts.baryLength_eq, e, max_eq_right (sq_nonneg _), sqrt_sq (mul_nonneg hl hj)]

/-- … and from the corner `j` it is `ℓ_ij b_i`. -/
theorem alg6_distance_on_edge' {lij : ℝ} (hl : 0 ≤ lij) {bi bj : ℝ} (hb : bi + bj = 1)
    (hi : 0 ≤ bi) (ljk lki : ℝ) :
    Ddg.baryLength lij ljk lki bi (bj - 1) 0 = lij * bi := by
  have e : Ddg.baryLength2 lij ljk lki bi (bj - 1) 0 = (lij * bi) ^ 2 := by
    rw [Cuts.baryLength2_eq, show bj - 1 = -bi by linarith]
    ring
  rw [Cuts.baryLength_eq, e, max_eq_right (sq_nonneg _), sqrt_sq (mul_nonneg hl hi)]

/-- **The oracle's `edgeSplitLengths` are the distances of the point on the side to its two
ends**: `cornerDistances` with `b₂ = 0` gives `ℓ₀₁ b₁` from corner `0` and `ℓ₀₁ b₀` from corner
`1`.

[formula: it-point-distances] -/
theorem edgeSplitLengths_eq_cornerDistances {l01 l12 l20 : ℝ} (h : IsTriangle l01 l12 l20)
    {b0 b1 : ℝ} (hb : b0 + b1 = 1) (h0 : 0 ≤ b0) (h1 : 0 ≤ b1) :
    (Ddg.cornerDistances l01 l12 l20 b0 b1 0).1 = (Ddg.edgeSplitLengths l01 b0 b1).1 ∧
      (Ddg.cornerDistances l01 l12 l20 b0 b1 0).2.1 = (Ddg.edgeSplitLengths l01 b0 b1).2 := by
  rw [cornerDistances_eq_baryLength h (by linarith : b0 + b1 + 0 = 1), edgeSplitLengths_eq]
  exact ⟨alg6_distance_on_edge h.pos_a.le hb h1 _ _, alg6_distance_on_edge' h.pos_a.le hb h0 _ _⟩

/-- **The paper's Alg. 8 as printed contradicts its Alg. 6**: for the point with the weights
`(b_i, b_j) = (3/4, 1/4)` on an edge of length `1`, Alg. 6 gives the distance `1/4 = ℓ_ij b_j`
from `i`; Alg. 8, line 4, prints `ℓ_ip ← ℓ_ij b_i = 3/4`. -/
theorem alg8_as_printed_differs (ljk lki : ℝ) :
    Ddg.baryLength 1 ljk lki ((3 / 4 : ℝ) - 1) (1 / 4) 0 = 1 / 4 ∧ (1 : ℝ) * (3 / 4) ≠ 1 / 4 := by
  refine ⟨?_, by norm_num⟩
  rw [alg6_distance_on_edge zero_le_one (by norm_num) (by norm_num)]
  norm_num

/-! ## The circumcentre

`a₀, a₁, a₂` are the sides opposite the corners `0, 1, 2`; the angle at corner `0` is
`cornerAngle a₀ a₁ a₂`. The layout has `ℓ₀₁ = a₂`, `ℓ₁₂ = a₀`, `ℓ₂₀ = a₁`. -/

theorem circumWeight_eq (ak ak1 ak2 : ℝ) :
    Ddg.circumWeight ak ak1 ak2 = ak * ak * (ak1 * ak1 + ak2 * ak2 - ak * ak) := rfl

/-- The sum of the three weights of the circumcentre. -/
def circumSum (a0 a1 a2 : ℝ) : ℝ :=
  Ddg.circumWeight a0 a1 a2 + Ddg.circumWeight a1 a2 a0 + Ddg.circumWeight a2 a0 a1

theorem circumcenterBary_eq (a0 a1 a2 : ℝ) :
    Ddg.circumcenterBary a0 a1 a2 =
      (Ddg.circumWeight a0 a1 a2 / circumSum a0 a1 a2,
        Ddg.circumWeight a1 a2 a0 / circumSum a0 a1 a2,
        Ddg.circumWeight a2 a0 a1 / circumSum a0 a1 a2) := rfl

/-- **The weights sum to `16 A²`** (so the normalization divides by a positive number). -/
theorem circumWeight_sum {a0 a1 a2 : ℝ} (h : IsTriangle a0 a1 a2) :
    circumSum a0 a1 a2 = 16 * triArea a0 a1 a2 ^ 2 := by
  rw [sixteen_mul_triArea_sq h]
  unfold circumSum
  simp only [circumWeight_eq]
  ring

theorem circumSum_pos {a0 a1 a2 : ℝ} (h : IsTriangle a0 a1 a2) : 0 < circumSum a0 a1 a2 := by
  rw [circumWeight_sum h]
  have := triArea_pos h
  positivity

/-- **The barycentric coordinates of the circumcentre sum to one.**

[formula: it-circumcenter] -/
theorem circumcenterBary_sum {a0 a1 a2 : ℝ} (h : IsTriangle a0 a1 a2) :
    (Ddg.circumcenterBary a0 a1 a2).1 + (Ddg.circumcenterBary a0 a1 a2).2.1 +
      (Ddg.circumcenterBary a0 a1 a2).2.2 = 1 := by
  rw [circumcenterBary_eq]
  dsimp only
  rw [← add_div, ← add_div]
  exact div_self (circumSum_pos h).ne'

/-- For `0 < x < π`: `cos x < 0 ⟺ x > π/2`. -/
theorem cos_neg_iff {x : ℝ} (h0 : 0 < x) (hπ : x < π) : cos x < 0 ↔ π / 2 < x := by
  constructor
  · intro h
    by_contra hcon
    have := Real.cos_nonneg_of_neg_pi_div_two_le_of_le (by linarith [pi_pos]) (not_lt.mp hcon)
    linarith
  · exact fun h => Real.cos_neg_of_pi_div_two_lt_of_lt h (by linarith [pi_pos])

/-- The weight of a corner is negative exactly when its angle is obtuse. -/
theorem circumWeight_neg_iff {a0 a1 a2 : ℝ} (h : IsTriangle a0 a1 a2) :
    Ddg.circumWeight a0 a1 a2 < 0 ↔ π / 2 < cornerAngle a0 a1 a2 := by
  have ha := h.pos_a
  have hb := h.pos_b
  have hc := h.pos_c
  have hpos : 0 < a0 * a0 := mul_pos ha ha
  have e : a1 * a1 + a2 * a2 - a0 * a0 = a1 ^ 2 + a2 ^ 2 - a0 ^ 2 := by ring
  rw [← cos_neg_iff (cornerAngle_pos h) (cornerAngle_lt_pi h), cos_cornerAngle h,
    circumWeight_eq, e, div_lt_iff₀ (by positivity), zero_mul]
  constructor
  · intro h'
    by_contra hcon
    exact absurd (mul_nonneg hpos.le (not_lt.mp hcon)) (not_le.mpr h')
  · exact fun h' => mul_neg_of_pos_of_neg hpos h'

/-- The weight of a corner vanishes exactly when its angle is a right angle. -/
theorem circumWeight_eq_zero_iff {a0 a1 a2 : ℝ} (h : IsTriangle a0 a1 a2) :
    Ddg.circumWeight a0 a1 a2 = 0 ↔ cornerAngle a0 a1 a2 = π / 2 := by
  have ha := h.pos_a
  have hb := h.pos_b
  have hc := h.pos_c
  have e : a1 * a1 + a2 * a2 - a0 * a0 = a1 ^ 2 + a2 ^ 2 - a0 ^ 2 := by ring
  have hcos := cos_cornerAngle h
  have h0 := cornerAngle_pos h
  have hπ := cornerAngle_lt_pi h
  rw [circumWeight_eq, e]
  constructor
  · intro hw
    have h1 : a1 ^ 2 + a2 ^ 2 - a0 ^ 2 = 0 := by
      rcases mul_eq_zero.mp hw with h' | h'
      · exact absurd h' (mul_pos ha ha).ne'
      · exact h'
    have h2 : cos (cornerAngle a0 a1 a2) = cos (π / 2) := by
      rw [hcos, h1, zero_div, Real.cos_pi_div_two]
    exact Real.injOn_cos ⟨h0.le, hπ.le⟩ ⟨by positivity, by linarith [pi_pos]⟩ h2
  · intro hθ
    rw [hθ, Real.cos_pi_div_two] at hcos
    have h1 : a1 ^ 2 + a2 ^ 2 - a0 ^ 2 = 0 := by
      have := hcos.symm
      rw [div_eq_zero_iff] at this
      rcases this with h' | h'
      · exact h'
      · exact absurd h' (by positivity)
    rw [h1, mul_zero]

/-- **A coordinate of the circumcentre is negative exactly when the angle at that corner is
obtuse** (the circumcentre then lies beyond the opposite side; the tracing rule of the
refinement relies on this).

[formula: it-circumcenter] -/
theorem circumcenterBary_neg_iff {a0 a1 a2 : ℝ} (h : IsTriangle a0 a1 a2) :
    ((Ddg.circumcenterBary a0 a1 a2).1 < 0 ↔ π / 2 < cornerAngle a0 a1 a2) ∧
      ((Ddg.circumcenterBary a0 a1 a2).2.1 < 0 ↔ π / 2 < cornerAngle a1 a2 a0) ∧
      ((Ddg.circumcenterBary a0 a1 a2).2.2 < 0 ↔ π / 2 < cornerAngle a2 a0 a1) := by
  have hs := circumSum_pos h
  rw [circumcenterBary_eq]
  dsimp only
  rw [div_lt_iff₀ hs, div_lt_iff₀ hs, div_lt_iff₀ hs, zero_mul]
  exact ⟨circumWeight_neg_iff h, circumWeight_neg_iff h.rotate,
    circumWeight_neg_iff h.rotate.rotate⟩

/-- … and zero exactly when that angle is a right angle (the circumcentre is on the opposite
side). -/
theorem circumcenterBary_eq_zero_iff {a0 a1 a2 : ℝ} (h : IsTriangle a0 a1 a2) :
    ((Ddg.circumcenterBary a0 a1 a2).1 = 0 ↔ cornerAngle a0 a1 a2 = π / 2) ∧
      ((Ddg.circumcenterBary a0 a1 a2).2.1 = 0 ↔ cornerAngle a1 a2 a0 = π / 2) ∧
      ((Ddg.circumcenterBary a0 a1 a2).2.2 = 0 ↔ cornerAngle a2 a0 a1 = π / 2) := by
  have hs := (circumSum_pos h).ne'
  rw [circumcenterBary_eq]
  dsimp only
  rw [div_eq_zero_iff, div_eq_zero_iff, div_eq_zero_iff]
  refine ⟨?_, ?_, ?_⟩
  · rw [← circumWeight_eq_zero_iff h]
    exact ⟨fun h' => h'.resolve_right hs, Or.inl⟩
  · rw [← circumWeight_eq_zero_iff h.rotate]
    exact ⟨fun h' => h'.resolve_right hs, Or.inl⟩
  · rw [← circumWeight_eq_zero_iff h.rotate.rotate]
    exact ⟨fun h' => h'.resolve_right hs, Or.inl⟩

/-- **At most one coordinate is negative**: if the coordinate of corner `0` is negative the
other two are positive (a triangle has at most one obtuse angle). -/
theorem circumcenterBary_pos_of_neg {a0 a1 a2 : ℝ} (h : IsTriangle a0 a1 a2)
    (hneg : (Ddg.circumcenterBary a0 a1 a2).1 < 0) :
    0 < (Ddg.circumcenterBary a0 a1 a2).2.1 ∧ 0 < (Ddg.circumcenterBary a0 a1 a2).2.2 := by
  have ha := h.pos_a
  have hb := h.pos_b
  have hc := h.pos_c
  have hs := circumSum_pos h
  rw [circumcenterBary_eq] at hneg ⊢
  dsimp only at hneg ⊢
  rw [div_lt_iff₀ hs, zero_mul, circumWeight_eq] at hneg
  have hlt : a1 * a1 + a2 * a2 - a0 * a0 < 0 := by
    by_contra hcon
    exact absurd (mul_nonneg (mul_pos ha ha).le (not_lt.mp hcon)) (not_le.mpr hneg)
  refine ⟨div_pos ?_ hs, div_pos ?_ hs⟩
  · rw [circumWeight_eq]
    exact mul_pos (mul_pos hb hb) (by nlinarith [mul_pos hc hc])
  · rw [circumWeight_eq]
    exact mul_pos (mul_pos hc hc) (by nlinarith [mul_pos hb hb])

/-- Algebra behind `circumcenter_equidistant`: clearing the denominator `S` of the three
coordinates `β_k = w_k/S`. -/
theorem circum_aux0 {S P x y z w0 w1 w2 β0 β1 β2 : ℝ} (hS : S ≠ 0) (h0 : β0 * S = w0)
    (h1 : β1 * S = w1) (h2 : β2 * S = w2)
    (poly : -((w0 - S) * w1 * x + w1 * w2 * y + w2 * (w0 - S) * z) = P * S) :
    -((β0 - 1) * β1 * x + β1 * β2 * y + β2 * (β0 - 1) * z) = P / S := by
  rw [eq_div_iff hS]
  apply mul_right_cancel₀ hS
  rw [← poly, ← h0, ← h1, ← h2]
  ring

theorem circum_aux1 {S P x y z w0 w1 w2 β0 β1 β2 : ℝ} (hS : S ≠ 0) (h0 : β0 * S = w0)
    (h1 : β1 * S = w1) (h2 : β2 * S = w2)
    (poly : -(w0 * (w1 - S) * x + (w1 - S) * w2 * y + w2 * w0 * z) = P * S) :
    -(β0 * (β1 - 1) * x + (β1 - 1) * β2 * y + β2 * β0 * z) = P / S := by
  rw [eq_div_iff hS]
  apply mul_right_cancel₀ hS
  rw [← poly, ← h0, ← h1, ← h2]
  ring

theorem circum_aux2 {S P x y z w0 w1 w2 β0 β1 β2 : ℝ} (hS : S ≠ 0) (h0 : β0 * S = w0)
    (h1 : β1 * S = w1) (h2 : β2 * S = w2)
    (poly : -(w0 * w1 * x + w1 * (w2 - S) * y + (w2 - S) * w0 * z) = P * S) :
    -(β0 * β1 * x + β1 * (β2 - 1) * y + (β2 - 1) * β0 * z) = P / S := by
  rw [eq_div_iff hS]
  apply mul_right_cancel₀ hS
  rw [← poly, ← h0, ← h1, ← h2]
  ring

/-- **The point with these coordinates is equidistant from the three corners**: in the layout
of the triangle (`ℓ₀₁ = a₂`, `ℓ₁₂ = a₀`, `ℓ₂₀ = a₁`) the oracle's `cornerDistances` of the
oracle's `circumcenterBary` are all the circumradius `a₀a₁a₂/(4A)`.

[formula: it-circumcenter] -/
theorem circumcenter_equidistant {a0 a1 a2 : ℝ} (h : IsTriangle a0 a1 a2) :
    Ddg.cornerDistances a2 a0 a1 (Ddg.circumcenterBary a0 a1 a2).1
        (Ddg.circumcenterBary a0 a1 a2).2.1 (Ddg.circumcenterBary a0 a1 a2).2.2 =
      (a0 * a1 * a2 / (4 * triArea a0 a1 a2), a0 * a1 * a2 / (4 * triArea a0 a1 a2),
        a0 * a1 * a2 / (4 * triArea a0 a1 a2)) := by
  have ha := h.pos_a
  have hb := h.pos_b
  have hc := h.pos_c
  have hT := triArea_pos h
  have hs := circumSum_pos h
  have hsE := circumWeight_sum h
  have hR : 0 ≤ a0 * a1 * a2 / (4 * triArea a0 a1 a2) := by positivity
  have hR2 : (a0 * a1 * a2 / (4 * triArea a0 a1 a2)) ^ 2 =
      (a0 * a1 * a2) ^ 2 / circumSum a0 a1 a2 := by
    rw [hsE, div_pow]
    congr 1
    ring
  have root : ∀ x : ℝ, x = (a0 * a1 * a2) ^ 2 / circumSum a0 a1 a2 →
      √(max 0 x) = a0 * a1 * a2 / (4 * triArea a0 a1 a2) := by
    intro x hx
    rw [hx, ← hR2, max_eq_right (sq_nonneg _), sqrt_sq hR]
  rw [cornerDistances_eq_baryLength h.rotate.rotate (circumcenterBary_sum h),
    Cuts.baryLength_eq, Cuts.baryLength_eq, Cuts.baryLength_eq]
  refine Prod.ext ?_ (Prod.ext ?_ ?_)
  · apply root
    rw [Cuts.baryLength2_eq, circumcenterBary_eq]
    dsimp only
    refine circum_aux0 hs.ne' (div_mul_cancel₀ _ hs.ne') (div_mul_cancel₀ _ hs.ne')
      (div_mul_cancel₀ _ hs.ne') ?_
    unfold circumSum
    simp only [circumWeight_eq]
    ring
  · apply root
    rw [Cuts.baryLength2_eq, circumcenterBary_eq]
    dsimp only
    refine circum_aux1 hs.ne' (div_mul_cancel₀ _ hs.ne') (div_mul_cancel₀ _ hs.ne')
      (div_mul_cancel₀ _ hs.ne') ?_
    unfold circumSum
    simp only [circumWeight_eq]
    ring
  · apply root
    rw [Cuts.baryLength2_eq, circumcenterBary_eq]
    dsimp only
    refine circum_aux2 hs.ne' (div_mul_cancel₀ _ hs.ne') (div_mul_cancel₀ _ hs.ne')
      (div_mul_cancel₀ _ hs.ne') ?_
    unfold circumSum
    simp only [circumWeight_eq]
    ring

end Intrinsic

end DdgProofs
