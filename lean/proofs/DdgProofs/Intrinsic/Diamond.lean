import DdgProofs.Oracle.Scalar
import DdgProofs.Bff.Quality

/-!
# Intrinsic triangulations: unfolding a triangle and a diamond, the flipped edge

`crates/ddg-intrinsic/SPEC.md` §4 (Flip), N5; "Navigating Intrinsic Triangulations" (Sharp,
Soliman, Crane 2019) Sec. 3.3.1, Fig. 9, Alg. 5.

The diamond of an edge `ij` of length `a`: the triangles `ijk` and `jil`, with `b = ℓ_ik`,
`c = ℓ_jk`, `d = ℓ_il`, `e = ℓ_jl`. The *standard layout* is `i = (0, 0)`, `j = (a, 0)`,
`k = (x_k, h_k)`, `l = (x_l, −h_l)` with `x = diamondX`, `h = diamondH`. Plane points are pairs
`ℝ × ℝ` with the Euclidean distance `dist2` and the doubled signed area `sda` of
`DdgProofs/Bff/Quality.lean`.

All tagged theorems are about the oracle's definitions (`Ddg.diamondX`, `diamondH`,
`flippedLength`, `flippedLengthFromLengths`, `flipArea1`, `flipArea2`, `flipConvex`,
`unfoldApex`) at `ℝ`. This package's own: `rigid` (a rotation followed by a translation).
-/

noncomputable section

namespace DdgProofs

namespace Intrinsic

open scoped DdgProofs.Oracle
open Real

/-! ## Plane helpers -/

theorem dist2_comm (p q : ℝ × ℝ) : dist2 p q = dist2 q p := by
  unfold dist2
  congr 1
  ring

theorem dist2_nonneg (p q : ℝ × ℝ) : 0 ≤ dist2 p q := sqrt_nonneg _

/-- `dist2 p q ^ 2` written out. -/
theorem dist2_sq' (p q : ℝ × ℝ) : dist2 p q ^ 2 = (q.1 - p.1) ^ 2 + (q.2 - p.2) ^ 2 :=
  dist2_sq p q

/-- A triangle of nonzero signed area has two different first corners. -/
theorem dist2_pos_of_sda_ne {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) : 0 < dist2 a b := by
  rcases (dist2_nonneg a b).lt_or_eq with hpos | hzero
  · exact hpos
  · exfalso
    have h2 := dist2_sq' a b
    rw [← hzero] at h2
    have h1 : (b.1 - a.1) ^ 2 = 0 := by nlinarith [sq_nonneg (b.1 - a.1), sq_nonneg (b.2 - a.2)]
    have h3 : (b.2 - a.2) ^ 2 = 0 := by nlinarith [sq_nonneg (b.1 - a.1), sq_nonneg (b.2 - a.2)]
    have e1 : b.1 - a.1 = 0 := by simpa using h1
    have e2 : b.2 - a.2 = 0 := by simpa using h3
    apply h
    unfold sda
    rw [e1, e2]
    ring

/-- Nonnegative numbers whose Heron product is positive satisfy the strict triangle
inequalities. -/
theorem isTriangle_of_heron_pos {x y z : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) (hz : 0 ≤ z)
    (hP : 0 < (x + y + z) * (y - x + z) * (x - y + z) * (x - z + y)) : IsTriangle x y z := by
  have hS : 0 ≤ x + y + z := by linarith
  refine ⟨?_, ?_, ?_⟩
  · by_contra hcon
    have h1 : y - x + z ≤ 0 := by linarith
    have h2 : 0 ≤ x - y + z := by linarith
    have h3 : 0 ≤ x - z + y := by linarith
    have : (x + y + z) * (y - x + z) * (x - y + z) * (x - z + y) =
        (y - x + z) * ((x + y + z) * (x - y + z) * (x - z + y)) := by ring
    rw [this] at hP
    have := mul_nonpos_of_nonpos_of_nonneg h1 (mul_nonneg (mul_nonneg hS h2) h3)
    linarith
  · by_contra hcon
    have h1 : x - y + z ≤ 0 := by linarith
    have h2 : 0 ≤ y - x + z := by linarith
    have h3 : 0 ≤ x - z + y := by linarith
    have : (x + y + z) * (y - x + z) * (x - y + z) * (x - z + y) =
        (x - y + z) * ((x + y + z) * (y - x + z) * (x - z + y)) := by ring
    rw [this] at hP
    have := mul_nonpos_of_nonpos_of_nonneg h1 (mul_nonneg (mul_nonneg hS h2) h3)
    linarith
  · by_contra hcon
    have h1 : x - z + y ≤ 0 := by linarith
    have h2 : 0 ≤ y - x + z := by linarith
    have h3 : 0 ≤ x - y + z := by linarith
    have : (x + y + z) * (y - x + z) * (x - y + z) * (x - z + y) =
        (x - z + y) * ((x + y + z) * (y - x + z) * (x - y + z)) := by ring
    rw [this] at hP
    have := mul_nonpos_of_nonpos_of_nonneg h1 (mul_nonneg (mul_nonneg hS h2) h3)
    linarith

/-- **The sides of a plane triangle of nonzero area satisfy the strict triangle
inequalities.** -/
theorem isTriangle_of_sda_ne {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) :
    IsTriangle (dist2 b c) (dist2 c a) (dist2 a b) := by
  apply isTriangle_of_heron_pos (dist2_nonneg _ _) (dist2_nonneg _ _) (dist2_nonneg _ _)
  have key : (dist2 b c + dist2 c a + dist2 a b) * (dist2 c a - dist2 b c + dist2 a b) *
      (dist2 b c - dist2 c a + dist2 a b) * (dist2 b c - dist2 a b + dist2 c a) =
      (2 * sda a b c) ^ 2 := by
    calc _ = 2 * dist2 b c ^ 2 * dist2 c a ^ 2 + 2 * dist2 c a ^ 2 * dist2 a b ^ 2 +
          2 * dist2 a b ^ 2 * dist2 b c ^ 2 - (dist2 b c ^ 2) ^ 2 - (dist2 c a ^ 2) ^ 2 -
          (dist2 a b ^ 2) ^ 2 := by ring
      _ = _ := by
        rw [dist2_sq, dist2_sq, dist2_sq]
        linear_combination (-1 : ℝ) * sda_sq_mul_four a b c
  rw [key]
  positivity

/-- The Heron area of the sides of a plane triangle, in the order of the oracle's diamond
(`a = |ij|`, `b = |ik|`, `c = |jk|`), is half the absolute signed double area. -/
theorem two_mul_triArea_dist2 (i j k : ℝ × ℝ) :
    2 * triArea (dist2 i j) (dist2 i k) (dist2 j k) = |sda i j k| := by
  rw [abs_sda_eq_two_mul_triArea i j k, dist2_comm k i]
  unfold triArea
  congr 3
  ring

/-- The triangle inequalities do not depend on the order of the last two sides. -/
theorem _root_.DdgProofs.IsTriangle.swap {a b c : ℝ} (h : IsTriangle a b c) : IsTriangle a c b :=
  ⟨by linarith [h.lt_a], by linarith [h.lt_c], by linarith [h.lt_b]⟩

/-! ## The oracle's definitions at `ℝ` -/

theorem diamondX_eq (a b c : ℝ) : Ddg.diamondX a b c = (a * a + (b - c) * (b + c)) / (2 * a) := by
  simp only [Ddg.diamondX, Oracle.two_eq]

theorem diamondH_eq (a b c : ℝ) : Ddg.diamondH a b c = 2 * triArea a b c / a := by
  simp only [Ddg.diamondH, Oracle.two_eq, Oracle.triArea_eq]

/-- The oracle's diamond is the reference layout of `DdgProofs/Bff/Quality.lean` (the triangle
`i, j, k` with `ℓ₀₁ = a`, `ℓ₁₂ = c`, `ℓ₂₀ = b`). -/
theorem diamondX_eq_layoutX (a b c : ℝ) : Ddg.diamondX a b c = layoutX a c b := by
  rw [diamondX_eq, layoutX]

theorem diamondH_eq_layoutY (a b c : ℝ) : Ddg.diamondH a b c = layoutY a (triArea a c b) := by
  rw [diamondH_eq, layoutY, triArea_swap]

theorem flippedLength_eq (xk hk xl hl : ℝ) :
    Ddg.flippedLength xk hk xl hl = √((xk - xl) * (xk - xl) + (hk + hl) * (hk + hl)) := by
  simp only [Ddg.flippedLength, Oracle.sqrt_eq]

/-- **The flipped length is the distance of the two apexes** `k = (x_k, h_k)` and
`l = (x_l, −h_l)`.

[formula: it-flip-length] -/
theorem flippedLength_eq_dist2 (xk hk xl hl : ℝ) :
    Ddg.flippedLength xk hk xl hl = dist2 (xk, hk) (xl, -hl) := by
  rw [flippedLength_eq]
  unfold dist2
  congr 1
  ring

theorem flippedLengthFromLengths_eq (a b c d e : ℝ) :
    Ddg.flippedLengthFromLengths a b c d e =
      Ddg.flippedLength (Ddg.diamondX a b c) (Ddg.diamondH a b c) (Ddg.diamondX a d e)
        (Ddg.diamondH a d e) := rfl

theorem flipArea1_eq (xk hk xl hl : ℝ) : Ddg.flipArea1 xk hk xl hl = xk * hl + xl * hk := rfl

theorem flipArea2_eq (a xk hk xl hl : ℝ) :
    Ddg.flipArea2 a xk hk xl hl = (a - xk) * hl + (a - xl) * hk := rfl

/-- `flipArea1` is the signed double area of the triangle `i l k` of the standard layout. -/
theorem flipArea1_eq_sda (xk hk xl hl : ℝ) :
    Ddg.flipArea1 xk hk xl hl = sda (0, 0) (xl, -hl) (xk, hk) := by
  rw [flipArea1_eq]
  unfold sda
  ring

/-- `flipArea2` is the signed double area of the triangle `j k l` of the standard layout. -/
theorem flipArea2_eq_sda (a xk hk xl hl : ℝ) :
    Ddg.flipArea2 a xk hk xl hl = sda (a, 0) (xk, hk) (xl, -hl) := by
  rw [flipArea2_eq]
  unfold sda
  ring

/-- The oracle's `flipConvex`: both new triangles have positive signed area. -/
theorem flipConvex_iff (a xk hk xl hl : ℝ) :
    Ddg.flipConvex a xk hk xl hl = true ↔
      0 < Ddg.flipArea1 xk hk xl hl ∧ 0 < Ddg.flipArea2 a xk hk xl hl := by
  simp only [Ddg.flipConvex, Bool.and_eq_true, decide_eq_true_eq, Oracle.zero_eq]

theorem unfoldApex_eq (ax ay dx dy lab lbc lca : ℝ) :
    Ddg.unfoldApex ax ay dx dy lab lbc lca =
      (ax + (lab * lab + lca * lca - lbc * lbc) / (2 * (lab * lab)) * dx -
          2 * triArea lab lbc lca / (lab * lab) * dy,
        ay + (lab * lab + lca * lca - lbc * lbc) / (2 * (lab * lab)) * dy +
          2 * triArea lab lbc lca / (lab * lab) * dx) := by
  simp only [Ddg.unfoldApex, Oracle.two_eq, Oracle.triArea_eq]

/-! ## The standard layout of a diamond -/

/-- The apex `k = (x_k, h_k)` of the standard layout, above the edge. -/
def apexK (a b c : ℝ) : ℝ × ℝ := (Ddg.diamondX a b c, Ddg.diamondH a b c)

/-- The apex `l = (x_l, −h_l)` of the standard layout, below the edge. -/
def apexL (a d e : ℝ) : ℝ × ℝ := (Ddg.diamondX a d e, -Ddg.diamondH a d e)

/-- The diamond with the five lengths is strictly convex in the oracle's sense: `flipConvex` of
its standard layout. -/
def DiamondConvex (a b c d e : ℝ) : Prop :=
  Ddg.flipConvex a (Ddg.diamondX a b c) (Ddg.diamondH a b c) (Ddg.diamondX a d e)
    (Ddg.diamondH a d e) = true

/-- **The standard layout of a diamond has the five given lengths, and the two apexes lie on
opposite sides of the edge**: for triangles `(a, b, c)` and `(a, d, e)`, with `i = (0, 0)`,
`j = (a, 0)`, `k = (x_k, h_k)`, `l = (x_l, −h_l)`: `h_k, h_l > 0`, `|ij| = a`, `|ik| = b`,
`|jk| = c`, `|il| = d`, `|jl| = e`.

[formula: it-unfold-diamond] -/
theorem diamond_isometric {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    0 < Ddg.diamondH a b c ∧ 0 < Ddg.diamondH a d e ∧
      dist2 (0, 0) (a, 0) = a ∧
      dist2 (0, 0) (apexK a b c) = b ∧ dist2 (a, 0) (apexK a b c) = c ∧
      dist2 (0, 0) (apexL a d e) = d ∧ dist2 (a, 0) (apexL a d e) = e := by
  obtain ⟨k1, k2, k3, k4⟩ := layout_isometric hk.swap rfl
  obtain ⟨l1, -, l3, l4⟩ := layout_isometric hl.swap rfl
  unfold apexK apexL
  rw [diamondX_eq_layoutX, diamondH_eq_layoutY, diamondX_eq_layoutX, diamondH_eq_layoutY]
  refine ⟨k1, l1, k2, k3, k4, Eq.trans ?_ l3, Eq.trans ?_ l4⟩
  · unfold dist2
    congr 1
    dsimp only
    ring
  · unfold dist2
    congr 1
    dsimp only
    ring

/-- Signed double areas of the two triangles of the standard layout: `a h_k` and `−a h_l`. -/
theorem diamond_sda (a xk hk xl hl : ℝ) :
    sda (0, 0) (a, 0) (xk, hk) = a * hk ∧ sda (0, 0) (a, 0) (xl, -hl) = -(a * hl) := by
  constructor <;>
  · unfold sda
    dsimp only
    ring

/-- The standard values of an arbitrary layout: the abscissa of the apex is the projection of
`k − i` on the direction of the edge. -/
theorem diamondX_realized {i j : ℝ × ℝ} (hij : 0 < dist2 i j) (k : ℝ × ℝ) :
    Ddg.diamondX (dist2 i j) (dist2 i k) (dist2 j k) =
      ((j.1 - i.1) * (k.1 - i.1) + (j.2 - i.2) * (k.2 - i.2)) / dist2 i j := by
  have ha := dist2_sq' i j
  have hb := dist2_sq' i k
  have hc := dist2_sq' j k
  rw [diamondX_eq, div_eq_div_iff (by positivity) hij.ne']
  linear_combination dist2 i j * (ha + hb - hc)

/-- … and the height is the absolute signed double area over the length of the edge. -/
theorem diamondH_realized (i j k : ℝ × ℝ) :
    Ddg.diamondH (dist2 i j) (dist2 i k) (dist2 j k) = |sda i j k| / dist2 i j := by
  rw [diamondH_eq, two_mul_triArea_dist2]

/-! ## Uniqueness up to a rigid motion -/

/-- The rotation with cosine `cs` and sine `sn`, followed by the translation by `t`. -/
def rigid (cs sn : ℝ) (t p : ℝ × ℝ) : ℝ × ℝ :=
  (t.1 + cs * p.1 - sn * p.2, t.2 + sn * p.1 + cs * p.2)

/-- The origin of the frame goes to the translation. -/
theorem rigid_zero (cs sn : ℝ) (t : ℝ × ℝ) : rigid cs sn t (0, 0) = t := by
  unfold rigid
  ext <;> simp

/-- With the direction of `b − a` the far end of the edge goes to `b`. -/
theorem rigid_edge (A B : ℝ × ℝ) {lab : ℝ} (hl : lab ≠ 0) :
    rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (lab, 0) = B := by
  unfold rigid
  ext <;>
  · dsimp only
    field_simp
    ring

/-- A rigid motion keeps distances. -/
theorem dist2_rigid {cs sn : ℝ} (h : cs ^ 2 + sn ^ 2 = 1) (t p q : ℝ × ℝ) :
    dist2 (rigid cs sn t p) (rigid cs sn t q) = dist2 p q := by
  unfold dist2 rigid
  congr 1
  dsimp only
  linear_combination ((q.1 - p.1) ^ 2 + (q.2 - p.2) ^ 2) * h

/-- A rigid motion keeps signed areas (it keeps the orientation). -/
theorem sda_rigid {cs sn : ℝ} (h : cs ^ 2 + sn ^ 2 = 1) (t p q r : ℝ × ℝ) :
    sda (rigid cs sn t p) (rigid cs sn t q) (rigid cs sn t r) = sda p q r := by
  unfold sda rigid
  dsimp only
  linear_combination ((q.1 - p.1) * (r.2 - p.2) - (q.2 - p.2) * (r.1 - p.1)) * h

/-- **The layout of a diamond is unique up to a rigid motion.** Four plane points `i, j, k, l`
with `k` to the left and `l` to the right of the directed line `i → j` are the standard layout
of their five distances `|ij|, |ik|, |jk|, |il|, |jl|`, rotated and translated.

[formula: it-unfold-diamond] -/
theorem diamond_unique {i j k l : ℝ × ℝ} (hk : 0 < sda i j k) (hl : sda i j l < 0) :
    ∃ cs sn : ℝ, cs ^ 2 + sn ^ 2 = 1 ∧
      i = rigid cs sn i (0, 0) ∧ j = rigid cs sn i (dist2 i j, 0) ∧
      k = rigid cs sn i (apexK (dist2 i j) (dist2 i k) (dist2 j k)) ∧
      l = rigid cs sn i (apexL (dist2 i j) (dist2 i l) (dist2 j l)) := by
  have ha : 0 < dist2 i j := dist2_pos_of_sda_ne hk.ne'
  have ha2 := dist2_sq' i j
  have hne : dist2 i j ≠ 0 := ha.ne'
  refine ⟨(j.1 - i.1) / dist2 i j, (j.2 - i.2) / dist2 i j, ?_, ?_, ?_, ?_, ?_⟩
  · field_simp
    linarith
  · unfold rigid
    ext <;> simp
  · unfold rigid
    ext <;>
    · dsimp only
      field_simp
      ring
  · unfold rigid apexK
    rw [diamondX_realized ha k, diamondH_realized, abs_of_pos hk]
    unfold sda
    ext
    · dsimp only
      field_simp
      linear_combination (k.1 - i.1) * ha2
    · dsimp only
      field_simp
      linear_combination (k.2 - i.2) * ha2
  · unfold rigid apexL
    rw [diamondX_realized ha l, diamondH_realized, abs_of_neg hl]
    unfold sda
    ext
    · dsimp only
      field_simp
      linear_combination (l.1 - i.1) * ha2
    · dsimp only
      field_simp
      linear_combination (l.2 - i.2) * ha2

/-- **The length formula gives the distance of the apexes in every layout**: for plane points
with `k` left and `l` right of `i → j`,
`flippedLengthFromLengths |ij| |ik| |jk| |il| |jl| = |kl|`.

[formula: it-flip-length] -/
theorem flippedLengthFromLengths_realized {i j k l : ℝ × ℝ} (hk : 0 < sda i j k)
    (hl : sda i j l < 0) :
    Ddg.flippedLengthFromLengths (dist2 i j) (dist2 i k) (dist2 j k) (dist2 i l) (dist2 j l) =
      dist2 k l := by
  obtain ⟨cs, sn, h1, -, -, ek, el⟩ := diamond_unique hk hl
  have e : dist2 k l = dist2 (rigid cs sn i (apexK (dist2 i j) (dist2 i k) (dist2 j k)))
      (rigid cs sn i (apexL (dist2 i j) (dist2 i l) (dist2 j l))) := congrArg₂ dist2 ek el
  rw [e, dist2_rigid h1, flippedLengthFromLengths_eq, flippedLength_eq_dist2]
  rfl

/-- **The two signed areas of `flipConvex` are those of the triangles `i l k` and `j k l` in
every layout.** -/
theorem flipArea_realized {i j k l : ℝ × ℝ} (hk : 0 < sda i j k) (hl : sda i j l < 0) :
    Ddg.flipArea1 (Ddg.diamondX (dist2 i j) (dist2 i k) (dist2 j k))
        (Ddg.diamondH (dist2 i j) (dist2 i k) (dist2 j k))
        (Ddg.diamondX (dist2 i j) (dist2 i l) (dist2 j l))
        (Ddg.diamondH (dist2 i j) (dist2 i l) (dist2 j l)) = sda i l k ∧
      Ddg.flipArea2 (dist2 i j) (Ddg.diamondX (dist2 i j) (dist2 i k) (dist2 j k))
        (Ddg.diamondH (dist2 i j) (dist2 i k) (dist2 j k))
        (Ddg.diamondX (dist2 i j) (dist2 i l) (dist2 j l))
        (Ddg.diamondH (dist2 i j) (dist2 i l) (dist2 j l)) = sda j k l := by
  obtain ⟨cs, sn, h1, ei, ej, ek, el⟩ := diamond_unique hk hl
  have e1 : sda i l k = sda (rigid cs sn i (0, 0))
      (rigid cs sn i (apexL (dist2 i j) (dist2 i l) (dist2 j l)))
      (rigid cs sn i (apexK (dist2 i j) (dist2 i k) (dist2 j k))) :=
    congr (congr (congrArg sda ei) el) ek
  have e2 : sda j k l = sda (rigid cs sn i (dist2 i j, 0))
      (rigid cs sn i (apexK (dist2 i j) (dist2 i k) (dist2 j k)))
      (rigid cs sn i (apexL (dist2 i j) (dist2 i l) (dist2 j l))) :=
    congr (congr (congrArg sda ej) ek) el
  rw [e1, e2, sda_rigid h1, sda_rigid h1, flipArea1_eq_sda, flipArea2_eq_sda]
  exact ⟨rfl, rfl⟩

/-! ## Unfolding one triangle -/

/-- The third vertex of a counter-clockwise plane triangle is determined by the first two and
the three lengths: it is the apex of the standard layout, rotated into the direction of the
edge. -/
theorem apex_unique {i j k : ℝ × ℝ} (hk : 0 < sda i j k) :
    k = rigid ((j.1 - i.1) / dist2 i j) ((j.2 - i.2) / dist2 i j) i
      (apexK (dist2 i j) (dist2 i k) (dist2 j k)) := by
  have ha : 0 < dist2 i j := dist2_pos_of_sda_ne hk.ne'
  have ha2 := dist2_sq' i j
  have hne : dist2 i j ≠ 0 := ha.ne'
  unfold rigid apexK
  rw [diamondX_realized ha k, diamondH_realized, abs_of_pos hk]
  unfold sda
  ext
  · dsimp only
    field_simp
    linear_combination (k.1 - i.1) * ha2
  · dsimp only
    field_simp
    linear_combination (k.2 - i.2) * ha2

/-- The oracle's `unfoldApex` is the apex of the standard layout, rotated by the direction
`d/ℓ_ab` and translated to `a`. -/
theorem unfoldApex_eq_rigid (A : ℝ × ℝ) (dx dy : ℝ) {lab : ℝ} (hl : lab ≠ 0) (lbc lca : ℝ) :
    Ddg.unfoldApex A.1 A.2 dx dy lab lbc lca =
      rigid (dx / lab) (dy / lab) A (apexK lab lca lbc) := by
  rw [unfoldApex_eq]
  unfold rigid apexK
  rw [diamondX_eq, diamondH_eq, triArea_swap]
  ext
  · dsimp only
    field_simp
    ring
  · dsimp only
    field_simp
    ring

/-- In the frame of the edge the oracle's `unfoldApex` is the apex of the diamond:
`unfoldApex (0,0) (a,0) a c b = (x_k, h_k)`. This is the Rust `triangle_layout`. -/
theorem unfoldApex_layout {a : ℝ} (ha : a ≠ 0) (b c : ℝ) :
    Ddg.unfoldApex 0 0 a 0 a c b = apexK a b c := by
  have := unfoldApex_eq_rigid (0, 0) a 0 ha c b
  dsimp only at this
  rw [this]
  unfold rigid
  ext <;>
  · dsimp only
    field_simp
    ring

/-- **The unfolded triangle has the three given lengths and lies to the left of `a → b`.** For
plane points `a, b` with `|ab| = ℓ_ab` and a triangle `(ℓ_ab, ℓ_bc, ℓ_ca)`, the point
`c = unfoldApex a (b − a) ℓ_ab ℓ_bc ℓ_ca` has `|ac| = ℓ_ca`, `|bc| = ℓ_bc`, and the signed double
area of `a, b, c` is `2A > 0`.

[formula: it-unfold-apex] -/
theorem unfoldApex_isometric {A B : ℝ × ℝ} {lab lbc lca : ℝ} (hAB : dist2 A B = lab)
    (h : IsTriangle lab lbc lca) :
    dist2 A (Ddg.unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) lab lbc lca) = lca ∧
      dist2 B (Ddg.unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) lab lbc lca) = lbc ∧
      sda A B (Ddg.unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) lab lbc lca) =
        2 * triArea lab lbc lca ∧
      0 < triArea lab lbc lca := by
  have hl : lab ≠ 0 := h.pos_a.ne'
  have hd := dist2_sq' A B
  rw [hAB] at hd
  have hcs : ((B.1 - A.1) / lab) ^ 2 + ((B.2 - A.2) / lab) ^ 2 = 1 := by
    field_simp
    linarith
  obtain ⟨-, -, -, k3, k4, -, -⟩ := diamond_isometric h.swap h.swap
  have hs := (diamond_sda lab (Ddg.diamondX lab lca lbc) (Ddg.diamondH lab lca lbc) 0 0).1
  rw [unfoldApex_eq_rigid A _ _ hl]
  refine ⟨?_, ?_, ?_, triArea_pos h⟩
  · have e : dist2 A (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) =
        dist2 (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (0, 0))
          (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) := by
      rw [rigid_zero]
    rw [e, dist2_rigid hcs]
    exact k3
  · have e : dist2 B (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) =
        dist2 (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (lab, 0))
          (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) := by
      rw [rigid_edge A B hl]
    rw [e, dist2_rigid hcs]
    exact k4
  · have e : sda A B (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) =
        sda (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (0, 0))
          (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (lab, 0))
          (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) := by
      rw [rigid_zero, rigid_edge A B hl]
    rw [e, sda_rigid hcs]
    unfold apexK
    rw [hs, diamondH_eq, triArea_swap]
    field_simp

/-- **Every counter-clockwise plane triangle is the one `unfoldApex` builds from its lengths**:
if `c` lies to the left of `a → b` then
`c = unfoldApex a (b − a) |ab| |bc| |ca|`. So the unfolded apex is unique.

[formula: it-unfold-apex] -/
theorem unfoldApex_unique {A B C : ℝ × ℝ} (h : 0 < sda A B C) :
    Ddg.unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) (dist2 A B) (dist2 B C) (dist2 C A) = C := by
  have ha : dist2 A B ≠ 0 := (dist2_pos_of_sda_ne h.ne').ne'
  rw [unfoldApex_eq_rigid A _ _ ha, dist2_comm C A]
  exact (apex_unique h).symm

end Intrinsic

end DdgProofs
