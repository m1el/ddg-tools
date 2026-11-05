import DDG.Geometry.PlaneLemmas

/-!
# Intrinsic triangulations: unfolding a triangle and a diamond, the flipped edge

`crates/ddg-intrinsic/SPEC.md` §4 (Flip), N5; "Navigating Intrinsic Triangulations" (Sharp,
Soliman, Crane 2019) Sec. 3.3.1, Fig. 9, Alg. 5.

The diamond of an edge `ij` of length `a`: the triangles `ijk` and `jil`, with `b = ℓ_ik`,
`c = ℓ_jk`, `d = ℓ_il`, `e = ℓ_jl`. The *standard layout* is `i = (0, 0)`, `j = (a, 0)`,
`k = (x_k, h_k)`, `l = (x_l, −h_l)` with `x = diamondX`, `h = diamondH`. Plane points are pairs
`ℝ × ℝ` with the Euclidean distance `dist2` and the doubled signed area `sda` of
`DDG/Geometry/Plane.lean`, with the plane lemmas of
`DDG/Geometry/PlaneLemmas.lean`.

Definitions (the closed forms of the Rust code, on `ℝ`): `diamondX`, `diamondH`,
`flippedLength`, `flippedLengthFromLengths`, `flipArea1`, `flipArea2`, `flipConvex`,
`unfoldApex`, and `rigid` (a rotation followed by a translation). The oracle's scalar-generic
copies are related to them by one-way bridges in the app layer (`Oracle/Intrinsic.lean`).
-/

noncomputable section

namespace DDG

namespace Intrinsic

open Real

/-! ## The formulas -/

/-- Abscissa of the apex of the triangle with base `a` (from the left end of the base) and the
other sides `b` (to the left end) and `c` (to the right end): `x = (a² + (b − c)(b + c))/2a`. -/
def diamondX (a b c : ℝ) : ℝ := (a * a + (b - c) * (b + c)) / (2 * a)

/-- Height of the apex over the base: `h = 2A/a`. -/
def diamondH (a b c : ℝ) : ℝ := 2 * triArea a b c / a

/-- Length of the flipped edge: `ℓ_kl = √((x_k − x_l)² + (h_k + h_l)²)`. -/
def flippedLength (xk hk xl hl : ℝ) : ℝ :=
  √((xk - xl) * (xk - xl) + (hk + hl) * (hk + hl))

/-- The flipped length from the five lengths `a, b, c, d, e` of the diamond. -/
def flippedLengthFromLengths (a b c d e : ℝ) : ℝ :=
  flippedLength (diamondX a b c) (diamondH a b c) (diamondX a d e) (diamondH a d e)

/-- Twice the signed area of the triangle `i l k` that a flip creates: `x_k h_l + x_l h_k`. -/
def flipArea1 (xk hk xl hl : ℝ) : ℝ := xk * hl + xl * hk

/-- Twice the signed area of the triangle `j k l` that a flip creates:
`(a − x_k) h_l + (a − x_l) h_k`. -/
def flipArea2 (a xk hk xl hl : ℝ) : ℝ := (a - xk) * hl + (a - xl) * hk

/-- The two triangles form a strictly convex quadrilateral: both new triangles have positive
area (the test of the Rust code, a `Bool`). -/
def flipConvex (a xk hk xl hl : ℝ) : Bool :=
  decide (0 < flipArea1 xk hk xl hl) && decide (0 < flipArea2 a xk hk xl hl)

/-- Third vertex of the triangle `a b c` laid out to the left of the directed edge `a → b`,
`d = b − a`: `c = a + α d + β d⊥`, `α = (ℓ_ab² + ℓ_ca² − ℓ_bc²)/2ℓ_ab²`, `β = 2A/ℓ_ab²`,
`d⊥ = (−d_y, d_x)`. -/
def unfoldApex (ax ay dx dy lab lbc lca : ℝ) : ℝ × ℝ :=
  (ax + (lab * lab + lca * lca - lbc * lbc) / (2 * (lab * lab)) * dx -
      2 * triArea lab lbc lca / (lab * lab) * dy,
    ay + (lab * lab + lca * lca - lbc * lbc) / (2 * (lab * lab)) * dy +
      2 * triArea lab lbc lca / (lab * lab) * dx)

theorem diamondX_eq (a b c : ℝ) : diamondX a b c = (a * a + (b - c) * (b + c)) / (2 * a) := rfl

theorem diamondH_eq (a b c : ℝ) : diamondH a b c = 2 * triArea a b c / a := rfl

/-- The diamond is the reference layout of `DDG/Geometry/Plane.lean` (the triangle
`i, j, k` with `ℓ₀₁ = a`, `ℓ₁₂ = c`, `ℓ₂₀ = b`). -/
theorem diamondX_eq_layoutX (a b c : ℝ) : diamondX a b c = layoutX a c b := by
  rw [diamondX_eq, layoutX]

theorem diamondH_eq_layoutY (a b c : ℝ) : diamondH a b c = layoutY a (triArea a c b) := by
  rw [diamondH_eq, layoutY, triArea_swap]

theorem flippedLength_eq (xk hk xl hl : ℝ) :
    flippedLength xk hk xl hl = √((xk - xl) * (xk - xl) + (hk + hl) * (hk + hl)) := rfl

/-- **The flipped length is the distance of the two apexes** `k = (x_k, h_k)` and
`l = (x_l, −h_l)`.

[formula: it-flip-length] -/
theorem flippedLength_eq_dist2 (xk hk xl hl : ℝ) :
    flippedLength xk hk xl hl = dist2 (xk, hk) (xl, -hl) := by
  rw [flippedLength_eq]
  unfold dist2
  congr 1
  ring

theorem flippedLengthFromLengths_eq (a b c d e : ℝ) :
    flippedLengthFromLengths a b c d e =
      flippedLength (diamondX a b c) (diamondH a b c) (diamondX a d e) (diamondH a d e) := rfl

theorem flipArea1_eq (xk hk xl hl : ℝ) : flipArea1 xk hk xl hl = xk * hl + xl * hk := rfl

theorem flipArea2_eq (a xk hk xl hl : ℝ) :
    flipArea2 a xk hk xl hl = (a - xk) * hl + (a - xl) * hk := rfl

/-- `flipArea1` is the signed double area of the triangle `i l k` of the standard layout. -/
theorem flipArea1_eq_sda (xk hk xl hl : ℝ) :
    flipArea1 xk hk xl hl = sda (0, 0) (xl, -hl) (xk, hk) := by
  rw [flipArea1_eq]
  unfold sda
  ring

/-- `flipArea2` is the signed double area of the triangle `j k l` of the standard layout. -/
theorem flipArea2_eq_sda (a xk hk xl hl : ℝ) :
    flipArea2 a xk hk xl hl = sda (a, 0) (xk, hk) (xl, -hl) := by
  rw [flipArea2_eq]
  unfold sda
  ring

/-- `flipConvex`: both new triangles have positive signed area. -/
theorem flipConvex_iff (a xk hk xl hl : ℝ) :
    flipConvex a xk hk xl hl = true ↔
      0 < flipArea1 xk hk xl hl ∧ 0 < flipArea2 a xk hk xl hl := by
  simp only [flipConvex, Bool.and_eq_true, decide_eq_true_eq]

theorem unfoldApex_eq (ax ay dx dy lab lbc lca : ℝ) :
    unfoldApex ax ay dx dy lab lbc lca =
      (ax + (lab * lab + lca * lca - lbc * lbc) / (2 * (lab * lab)) * dx -
          2 * triArea lab lbc lca / (lab * lab) * dy,
        ay + (lab * lab + lca * lca - lbc * lbc) / (2 * (lab * lab)) * dy +
          2 * triArea lab lbc lca / (lab * lab) * dx) := rfl

/-! ## The standard layout of a diamond -/

/-- The apex `k = (x_k, h_k)` of the standard layout, above the edge. -/
def apexK (a b c : ℝ) : ℝ × ℝ := (diamondX a b c, diamondH a b c)

/-- The apex `l = (x_l, −h_l)` of the standard layout, below the edge. -/
def apexL (a d e : ℝ) : ℝ × ℝ := (diamondX a d e, -diamondH a d e)

/-- The diamond with the five lengths is strictly convex: `flipConvex` of
its standard layout. -/
def DiamondConvex (a b c d e : ℝ) : Prop :=
  flipConvex a (diamondX a b c) (diamondH a b c) (diamondX a d e)
    (diamondH a d e) = true

/-- **The standard layout of a diamond has the five given lengths, and the two apexes lie on
opposite sides of the edge**: for triangles `(a, b, c)` and `(a, d, e)`, with `i = (0, 0)`,
`j = (a, 0)`, `k = (x_k, h_k)`, `l = (x_l, −h_l)`: `h_k, h_l > 0`, `|ij| = a`, `|ik| = b`,
`|jk| = c`, `|il| = d`, `|jl| = e`.

[formula: it-unfold-diamond] -/
theorem diamond_isometric {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    0 < diamondH a b c ∧ 0 < diamondH a d e ∧
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
    diamondX (dist2 i j) (dist2 i k) (dist2 j k) =
      ((j.1 - i.1) * (k.1 - i.1) + (j.2 - i.2) * (k.2 - i.2)) / dist2 i j := by
  have ha := dist2_sq' i j
  have hb := dist2_sq' i k
  have hc := dist2_sq' j k
  rw [diamondX_eq, div_eq_div_iff (by positivity) hij.ne']
  linear_combination dist2 i j * (ha + hb - hc)

/-- … and the height is the absolute signed double area over the length of the edge. -/
theorem diamondH_realized (i j k : ℝ × ℝ) :
    diamondH (dist2 i j) (dist2 i k) (dist2 j k) = |sda i j k| / dist2 i j := by
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

/-- Without `cs² + sn² = 1` the map is a similarity: distances are scaled by `√(cs² + sn²)`. -/
theorem dist2_rigid_sq (cs sn : ℝ) (t p q : ℝ × ℝ) :
    dist2 (rigid cs sn t p) (rigid cs sn t q) ^ 2 = (cs ^ 2 + sn ^ 2) * dist2 p q ^ 2 := by
  rw [dist2_sq', dist2_sq']
  unfold rigid
  dsimp only
  ring

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
    flippedLengthFromLengths (dist2 i j) (dist2 i k) (dist2 j k) (dist2 i l) (dist2 j l) =
      dist2 k l := by
  obtain ⟨cs, sn, h1, -, -, ek, el⟩ := diamond_unique hk hl
  have e : dist2 k l = dist2 (rigid cs sn i (apexK (dist2 i j) (dist2 i k) (dist2 j k)))
      (rigid cs sn i (apexL (dist2 i j) (dist2 i l) (dist2 j l))) := congrArg₂ dist2 ek el
  rw [e, dist2_rigid h1, flippedLengthFromLengths_eq, flippedLength_eq_dist2]
  rfl

/-- **The two signed areas of `flipConvex` are those of the triangles `i l k` and `j k l` in
every layout.** -/
theorem flipArea_realized {i j k l : ℝ × ℝ} (hk : 0 < sda i j k) (hl : sda i j l < 0) :
    flipArea1 (diamondX (dist2 i j) (dist2 i k) (dist2 j k))
        (diamondH (dist2 i j) (dist2 i k) (dist2 j k))
        (diamondX (dist2 i j) (dist2 i l) (dist2 j l))
        (diamondH (dist2 i j) (dist2 i l) (dist2 j l)) = sda i l k ∧
      flipArea2 (dist2 i j) (diamondX (dist2 i j) (dist2 i k) (dist2 j k))
        (diamondH (dist2 i j) (dist2 i k) (dist2 j k))
        (diamondX (dist2 i j) (dist2 i l) (dist2 j l))
        (diamondH (dist2 i j) (dist2 i l) (dist2 j l)) = sda j k l := by
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

/-- `unfoldApex` is the apex of the standard layout, rotated by the direction
`d/ℓ_ab` and translated to `a`. -/
theorem unfoldApex_eq_rigid (A : ℝ × ℝ) (dx dy : ℝ) {lab : ℝ} (hl : lab ≠ 0) (lbc lca : ℝ) :
    unfoldApex A.1 A.2 dx dy lab lbc lca =
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

/-- In the frame of the edge `unfoldApex` is the apex of the diamond:
`unfoldApex (0,0) (a,0) a c b = (x_k, h_k)`. This is the Rust `triangle_layout`. -/
theorem unfoldApex_layout {a : ℝ} (ha : a ≠ 0) (b c : ℝ) :
    unfoldApex 0 0 a 0 a c b = apexK a b c := by
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
    dist2 A (unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) lab lbc lca) = lca ∧
      dist2 B (unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) lab lbc lca) = lbc ∧
      sda A B (unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) lab lbc lca) =
        2 * triArea lab lbc lca ∧
      0 < triArea lab lbc lca := by
  have hl : lab ≠ 0 := h.pos_a.ne'
  have hd := dist2_sq' A B
  rw [hAB] at hd
  have hcs : ((B.1 - A.1) / lab) ^ 2 + ((B.2 - A.2) / lab) ^ 2 = 1 := by
    field_simp
    linarith
  obtain ⟨-, -, -, k3, k4, -, -⟩ := diamond_isometric h.swap h.swap
  have hs := (diamond_sda lab (diamondX lab lca lbc) (diamondH lab lca lbc) 0 0).1
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

/-- **If the placed edge does not have exactly the length `ℓ_ab`, the unfolded triangle is the
given one scaled by `|b − a|/ℓ_ab`**: the angles are those of the three given lengths, the size
follows the placed edge.

[formula: it-unfold-apex] -/
theorem unfoldApex_similar (A B : ℝ × ℝ) {lab lbc lca : ℝ} (h : IsTriangle lab lbc lca) :
    dist2 A (unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) lab lbc lca) ^ 2 =
        (dist2 A B / lab) ^ 2 * lca ^ 2 ∧
      dist2 B (unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) lab lbc lca) ^ 2 =
        (dist2 A B / lab) ^ 2 * lbc ^ 2 := by
  have hl : lab ≠ 0 := h.pos_a.ne'
  have hd := dist2_sq' A B
  have hcs : ((B.1 - A.1) / lab) ^ 2 + ((B.2 - A.2) / lab) ^ 2 = (dist2 A B / lab) ^ 2 := by
    rw [div_pow, div_pow, div_pow, hd]
    ring
  obtain ⟨-, -, -, k3, k4, -, -⟩ := diamond_isometric h.swap h.swap
  rw [unfoldApex_eq_rigid A _ _ hl]
  constructor
  · have e : dist2 A (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) =
        dist2 (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (0, 0))
          (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) := by
      rw [rigid_zero]
    rw [e, dist2_rigid_sq, hcs, k3]
  · have e : dist2 B (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) =
        dist2 (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (lab, 0))
          (rigid ((B.1 - A.1) / lab) ((B.2 - A.2) / lab) A (apexK lab lca lbc)) := by
      rw [rigid_edge A B hl]
    rw [e, dist2_rigid_sq, hcs, k4]

/-- **Every counter-clockwise plane triangle is the one `unfoldApex` builds from its lengths**:
if `c` lies to the left of `a → b` then
`c = unfoldApex a (b − a) |ab| |bc| |ca|`. So the unfolded apex is unique.

[formula: it-unfold-apex] -/
theorem unfoldApex_unique {A B C : ℝ × ℝ} (h : 0 < sda A B C) :
    unfoldApex A.1 A.2 (B.1 - A.1) (B.2 - A.2) (dist2 A B) (dist2 B C) (dist2 C A) = C := by
  have ha : dist2 A B ≠ 0 := (dist2_pos_of_sda_ne h.ne').ne'
  rw [unfoldApex_eq_rigid A _ _ ha, dist2_comm C A]
  exact (apex_unique h).symm

/-! ## The diamond in terms of the corner angles

In the triangle `ijk` with `a = |ij|`, `b = |ik|`, `c = |jk|` the angle at `i` is
`cornerAngle c a b`, at `j` it is `cornerAngle b c a`, at `k` it is `cornerAngle a b c`; in
`jil` with `d = |il|`, `e = |jl|`: `cornerAngle e a d` at `i`, `cornerAngle d e a` at `j`,
`cornerAngle a d e` at `l`. -/

theorem diamondX_eq_cos {a b c : ℝ} (h : IsTriangle a b c) :
    diamondX a b c = b * cos (cornerAngle c a b) := by
  have ha := h.pos_a
  have hb := h.pos_b
  rw [cos_cornerAngle h.rotate.rotate, diamondX_eq]
  field_simp
  ring

theorem diamondH_eq_sin {a b c : ℝ} (h : IsTriangle a b c) :
    diamondH a b c = b * sin (cornerAngle c a b) := by
  have ha := h.pos_a
  have hb := h.pos_b
  rw [sin_cornerAngle h.rotate.rotate, diamondH_eq, ← triArea_rotate c a b]
  field_simp

theorem diamondX_eq_cos' {a b c : ℝ} (h : IsTriangle a b c) :
    a - diamondX a b c = c * cos (cornerAngle b c a) := by
  have ha := h.pos_a
  have hc := h.pos_c
  rw [cos_cornerAngle h.rotate, diamondX_eq]
  field_simp
  ring

theorem diamondH_eq_sin' {a b c : ℝ} (h : IsTriangle a b c) :
    diamondH a b c = c * sin (cornerAngle b c a) := by
  have ha := h.pos_a
  have hc := h.pos_c
  rw [sin_cornerAngle h.rotate, diamondH_eq, triArea_rotate a b c]
  field_simp

/-- **The first area of `flipConvex` is `ℓ_ik ℓ_il sin` of the angle of the quadrilateral at
`i`.** -/
theorem flipArea1_eq_sin {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    flipArea1 (diamondX a b c) (diamondH a b c) (diamondX a d e)
        (diamondH a d e) = b * d * sin (cornerAngle c a b + cornerAngle e a d) := by
  rw [flipArea1_eq, diamondX_eq_cos hk, diamondH_eq_sin hk, diamondX_eq_cos hl,
    diamondH_eq_sin hl, sin_add]
  ring

/-- **The second is `ℓ_jk ℓ_jl sin` of the angle of the quadrilateral at `j`.** -/
theorem flipArea2_eq_sin {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    flipArea2 a (diamondX a b c) (diamondH a b c) (diamondX a d e)
        (diamondH a d e) = c * e * sin (cornerAngle b c a + cornerAngle d e a) := by
  rw [flipArea2_eq, diamondX_eq_cos' hk, diamondX_eq_cos' hl, diamondH_eq_sin' hk,
    diamondH_eq_sin' hl, sin_add]
  ring

/-- For `0 < x < 2π`: `sin x > 0 ⟺ x < π`. -/
theorem sin_pos_iff_lt_pi {x : ℝ} (h0 : 0 < x) (h2 : x < 2 * π) : 0 < sin x ↔ x < π := by
  constructor
  · intro h
    by_contra hcon
    have h1 := Real.sin_nonneg_of_nonneg_of_le_pi (sub_nonneg.mpr (not_lt.mp hcon))
      (by linarith : x - π ≤ π)
    rw [Real.sin_sub_pi] at h1
    linarith
  · exact fun h => Real.sin_pos_of_pos_of_lt_pi h0 h

/-- **Strict convexity in angles**: the diamond is strictly convex exactly when the angles of
the quadrilateral at both ends of the edge are less than `π`.

[formula: it-flip-convexity] -/
theorem diamondConvex_iff_angles {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    DiamondConvex a b c d e ↔
      cornerAngle c a b + cornerAngle e a d < π ∧ cornerAngle b c a + cornerAngle d e a < π := by
  have hb := hk.pos_b
  have hc := hk.pos_c
  have hd := hl.pos_b
  have he := hl.pos_c
  have a1 := cornerAngle_pos hk.rotate.rotate
  have a2 := cornerAngle_lt_pi hk.rotate.rotate
  have a3 := cornerAngle_pos hl.rotate.rotate
  have a4 := cornerAngle_lt_pi hl.rotate.rotate
  have b1 := cornerAngle_pos hk.rotate
  have b2 := cornerAngle_lt_pi hk.rotate
  have b3 := cornerAngle_pos hl.rotate
  have b4 := cornerAngle_lt_pi hl.rotate
  unfold DiamondConvex
  rw [flipConvex_iff, flipArea1_eq_sin hk hl, flipArea2_eq_sin hk hl,
    mul_pos_iff_of_pos_left (mul_pos hb hd), mul_pos_iff_of_pos_left (mul_pos hc he),
    sin_pos_iff_lt_pi (by linarith) (by linarith), sin_pos_iff_lt_pi (by linarith) (by linarith)]

/-! ## The length of the flipped edge -/

/-- The distance of the two apexes of the standard layout is the flipped length. -/
theorem dist2_apex (a b c d e : ℝ) :
    dist2 (apexK a b c) (apexL a d e) = flippedLengthFromLengths a b c d e := by
  rw [flippedLengthFromLengths_eq, flippedLength_eq_dist2]
  rfl

/-- **The paper's formula for the flipped edge** (Fig. 9, Alg. 5): the law of cosines at `i`
with the angle `θ_i^{jk} + θ_i^{lj}` between `il` and `ik`,
`ℓ_kl² = ℓ_ik² + ℓ_il² − 2 ℓ_ik ℓ_il cos θ_i^{lk}`.

[formula: it-flip-length] -/
theorem flippedLengthFromLengths_sq {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    flippedLengthFromLengths a b c d e ^ 2 =
      b ^ 2 + d ^ 2 - 2 * b * d * cos (cornerAngle c a b + cornerAngle e a d) := by
  rw [flippedLengthFromLengths_eq, flippedLength_eq, sq_sqrt (add_nonneg (mul_self_nonneg _) (mul_self_nonneg _)),
    diamondX_eq_cos hk, diamondH_eq_sin hk, diamondX_eq_cos hl, diamondH_eq_sin hl, cos_add]
  have h1 := sin_sq_add_cos_sq (cornerAngle c a b)
  have h2 := sin_sq_add_cos_sq (cornerAngle e a d)
  linear_combination b ^ 2 * h1 + d ^ 2 * h2

/-- The same at `j`: `ℓ_kl² = ℓ_jk² + ℓ_jl² − 2 ℓ_jk ℓ_jl cos(θ_j^{ki} + θ_j^{il})`. -/
theorem flippedLengthFromLengths_sq' {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    flippedLengthFromLengths a b c d e ^ 2 =
      c ^ 2 + e ^ 2 - 2 * c * e * cos (cornerAngle b c a + cornerAngle d e a) := by
  have e1 := diamondX_eq_cos' hk
  have e2 := diamondX_eq_cos' hl
  have e3 : diamondX a b c - diamondX a d e =
      e * cos (cornerAngle d e a) - c * cos (cornerAngle b c a) := by linarith
  rw [flippedLengthFromLengths_eq, flippedLength_eq, sq_sqrt (add_nonneg (mul_self_nonneg _) (mul_self_nonneg _)), e3,
    diamondH_eq_sin' hk, diamondH_eq_sin' hl, cos_add]
  have h1 := sin_sq_add_cos_sq (cornerAngle b c a)
  have h2 := sin_sq_add_cos_sq (cornerAngle d e a)
  linear_combination c ^ 2 * h1 + e ^ 2 * h2

/-- **The flipped length in the five lengths**:
`ℓ_kl² = b² + d² − [(a² + b² − c²)(a² + d² − e²) − 16 A_ijk A_jil] / 2a²`.

[formula: it-flip-length] -/
theorem flippedLengthFromLengths_sq_lengths {a b c d e : ℝ} (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    flippedLengthFromLengths a b c d e ^ 2 =
      b ^ 2 + d ^ 2 - ((a ^ 2 + b ^ 2 - c ^ 2) * (a ^ 2 + d ^ 2 - e ^ 2) -
        16 * triArea a b c * triArea a d e) / (2 * a ^ 2) := by
  have ha := hk.pos_a
  have h2a : 2 * a ^ 2 ≠ 0 := by positivity
  have hA := sixteen_mul_triArea_sq hk
  have hB := sixteen_mul_triArea_sq hl
  have hx : diamondX a b c * (2 * a) = a ^ 2 + b ^ 2 - c ^ 2 := by
    rw [diamondX_eq, div_mul_cancel₀ _ (by positivity)]
    ring
  have hxl : diamondX a d e * (2 * a) = a ^ 2 + d ^ 2 - e ^ 2 := by
    rw [diamondX_eq, div_mul_cancel₀ _ (by positivity)]
    ring
  have hh : diamondH a b c * a = 2 * triArea a b c := by
    rw [diamondH_eq, div_mul_cancel₀ _ ha.ne']
  have hhl : diamondH a d e * a = 2 * triArea a d e := by
    rw [diamondH_eq, div_mul_cancel₀ _ ha.ne']
  rw [flippedLengthFromLengths_eq, flippedLength_eq, sq_sqrt (add_nonneg (mul_self_nonneg _) (mul_self_nonneg _))]
  have key : ((diamondX a b c - diamondX a d e) *
        (diamondX a b c - diamondX a d e) +
      (diamondH a b c + diamondH a d e) * (diamondH a b c + diamondH a d e)) *
        (4 * a ^ 2) =
      (b ^ 2 + d ^ 2) * (4 * a ^ 2) - 2 * ((a ^ 2 + b ^ 2 - c ^ 2) * (a ^ 2 + d ^ 2 - e ^ 2) -
        16 * triArea a b c * triArea a d e) := by
    linear_combination
      (diamondX a b c * (2 * a) - diamondX a d e * (2 * a) + (a ^ 2 + b ^ 2 - c ^ 2) -
          (a ^ 2 + d ^ 2 - e ^ 2)) * hx -
        (diamondX a b c * (2 * a) - diamondX a d e * (2 * a) + (a ^ 2 + b ^ 2 - c ^ 2) -
          (a ^ 2 + d ^ 2 - e ^ 2)) * hxl +
        4 * (diamondH a b c * a + diamondH a d e * a + 2 * triArea a b c +
          2 * triArea a d e) * hh +
        4 * (diamondH a b c * a + diamondH a d e * a + 2 * triArea a b c +
          2 * triArea a d e) * hhl + hA + hB
  have e2 : ((a ^ 2 + b ^ 2 - c ^ 2) * (a ^ 2 + d ^ 2 - e ^ 2) -
      16 * triArea a b c * triArea a d e) / (2 * a ^ 2) =
      b ^ 2 + d ^ 2 - ((diamondX a b c - diamondX a d e) *
        (diamondX a b c - diamondX a d e) +
      (diamondH a b c + diamondH a d e) * (diamondH a b c + diamondH a d e)) := by
    rw [div_eq_iff h2a]
    linear_combination (1 / 2 : ℝ) * key
  linarith

theorem diamondX_swap {a : ℝ} (ha : a ≠ 0) (b c : ℝ) :
    diamondX a c b = a - diamondX a b c := by
  rw [diamondX_eq, diamondX_eq]
  field_simp
  ring

theorem diamondH_swap (a b c : ℝ) : diamondH a c b = diamondH a b c := by
  rw [diamondH_eq, diamondH_eq, triArea_swap]

/-- The flipped length does not change when the diamond is mirrored (`i ↔ j`). -/
theorem flippedLengthFromLengths_mirror {a : ℝ} (ha : a ≠ 0) (b c d e : ℝ) :
    flippedLengthFromLengths a c b e d = flippedLengthFromLengths a b c d e := by
  rw [flippedLengthFromLengths_eq, flippedLengthFromLengths_eq, flippedLength_eq,
    flippedLength_eq, diamondX_swap ha b c, diamondX_swap ha d e, diamondH_swap a b c,
    diamondH_swap a d e]
  congr 1
  ring

/-- … nor when the two triangles are exchanged (`k ↔ l`). -/
theorem flippedLengthFromLengths_comm (a b c d e : ℝ) :
    flippedLengthFromLengths a d e b c = flippedLengthFromLengths a b c d e := by
  rw [flippedLengthFromLengths_eq, flippedLengthFromLengths_eq, flippedLength_eq,
    flippedLength_eq]
  congr 1
  ring

/-! ## The flip of a strictly convex diamond -/

/-- The two conditions of `flipConvex` on the points of the standard layout. -/
theorem diamondConvex_sda {a b c d e : ℝ} (h : DiamondConvex a b c d e) :
    0 < sda (0, 0) (apexL a d e) (apexK a b c) ∧ 0 < sda (a, 0) (apexK a b c) (apexL a d e) := by
  obtain ⟨c1, c2⟩ := (flipConvex_iff _ _ _ _ _).mp h
  rw [flipArea1_eq_sda] at c1
  rw [flipArea2_eq_sda] at c2
  exact ⟨c1, c2⟩

/-- **Flipping back gives the original length.** In the diamond of the new edge `kl` (length
`f`, apexes `j` and `i`, sides `ℓ_kj = c`, `ℓ_lj = e`, `ℓ_ki = b`, `ℓ_li = d`) the flipped
length is `a`.

[formula: it-flip-length] -/
theorem flip_back {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) :
    flippedLengthFromLengths (flippedLengthFromLengths a b c d e) c e b d = a := by
  obtain ⟨-, -, e1, e2, e3, e4, e5⟩ := diamond_isometric hk hl
  obtain ⟨c1, c2⟩ := diamondConvex_sda hconv
  have hJ : 0 < sda (apexK a b c) (apexL a d e) (a, 0) := by
    rw [sda_rotate]
    exact c2
  have hI : sda (apexK a b c) (apexL a d e) (0, 0) < 0 := by
    rw [sda_rotate, sda_swap_right]
    linarith
  have h := flippedLengthFromLengths_realized hJ hI
  rw [dist2_apex, dist2_comm (apexK a b c) (a, 0), dist2_comm (apexL a d e) (a, 0),
    dist2_comm (apexK a b c) (0, 0), dist2_comm (apexL a d e) (0, 0),
    dist2_comm (a, 0) (0, 0), e1, e2, e3, e4, e5] at h
  exact h

/-- **The two triangles a flip creates satisfy the strict triangle inequalities**: `ilk` with
the sides `f, b, d` and `jkl` with the sides `f, e, c` (`f` the flipped length).

[formula: it-flip-length] -/
theorem flip_isTriangle {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) :
    IsTriangle (flippedLengthFromLengths a b c d e) b d ∧
      IsTriangle (flippedLengthFromLengths a b c d e) e c := by
  obtain ⟨-, -, -, e2, e3, e4, e5⟩ := diamond_isometric hk hl
  obtain ⟨c1, c2⟩ := diamondConvex_sda hconv
  have t1 := isTriangle_of_sda_ne c1.ne'
  have t2 := isTriangle_of_sda_ne c2.ne'
  rw [dist2_comm (apexL a d e) (apexK a b c), dist2_apex, dist2_comm (apexK a b c) (0, 0), e2,
    e4] at t1
  rw [dist2_apex, dist2_comm (apexL a d e) (a, 0), e5, e3] at t2
  exact ⟨t1, t2⟩

/-- **A flip keeps the total area**: `A_ilk + A_jkl = A_ijk + A_jil`.

[formula: it-flip-length] -/
theorem flip_area {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) :
    triArea (flippedLengthFromLengths a b c d e) b d +
        triArea (flippedLengthFromLengths a b c d e) e c =
      triArea a b c + triArea a d e := by
  have ha := hk.pos_a
  obtain ⟨-, -, -, e2, e3, e4, e5⟩ := diamond_isometric hk hl
  obtain ⟨c1, c2⟩ := diamondConvex_sda hconv
  have t1 := abs_sda_eq_two_mul_triArea (0, 0) (apexL a d e) (apexK a b c)
  have t2 := abs_sda_eq_two_mul_triArea (a, 0) (apexK a b c) (apexL a d e)
  rw [abs_of_pos c1, dist2_comm (apexL a d e) (apexK a b c), dist2_apex,
    dist2_comm (apexK a b c) (0, 0), e2, e4] at t1
  rw [abs_of_pos c2, dist2_apex, dist2_comm (apexL a d e) (a, 0), e5, e3] at t2
  have s1 : sda (0, 0) (apexL a d e) (apexK a b c) + sda (a, 0) (apexK a b c) (apexL a d e) =
      a * (diamondH a b c + diamondH a d e) := by
    unfold sda apexK apexL
    dsimp only
    ring
  have hh : a * diamondH a b c = 2 * triArea a b c := by
    rw [diamondH_eq, mul_div_cancel₀ _ ha.ne']
  have hhl : a * diamondH a d e = 2 * triArea a d e := by
    rw [diamondH_eq, mul_div_cancel₀ _ ha.ne']
  linarith

/-! ## Strict convexity -/

/-- **`flipConvex` holds exactly when the two diagonals of the unfolded quadrilateral cross in
their interiors**: for `a, h_k, h_l > 0`, the segments `ij` (from `(0,0)` to `(a,0)`) and `kl`
(from `(x_k, h_k)` to `(x_l, −h_l)`) meet in a point `(s a, 0)` with `0 < s < 1` that is
`k + t (l − k)` with `0 < t < 1`.

[formula: it-flip-convexity] -/
theorem flipConvex_iff_diagonals_cross {a xk hk xl hl : ℝ} (ha : 0 < a) (hhk : 0 < hk)
    (hhl : 0 < hl) :
    flipConvex a xk hk xl hl = true ↔
      ∃ s t : ℝ, 0 < s ∧ s < 1 ∧ 0 < t ∧ t < 1 ∧
        (s * a, (0 : ℝ)) = (xk + t * (xl - xk), hk + t * (-hl - hk)) := by
  rw [flipConvex_iff, flipArea1_eq, flipArea2_eq]
  have hs : 0 < hk + hl := by linarith
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨(xk * hl + xl * hk) / (a * (hk + hl)), hk / (hk + hl), by positivity, ?_,
      by positivity, ?_, ?_⟩
    · rw [div_lt_one (by positivity)]
      nlinarith
    · rw [div_lt_one hs]
      linarith
    · ext
      · dsimp only
        field_simp
        ring
      · dsimp only
        field_simp
        ring
  · rintro ⟨s, t, hs0, hs1, ht0, ht1, heq⟩
    obtain ⟨e1, e2⟩ := Prod.mk.inj heq
    have ht : t * (hk + hl) = hk := by linarith
    have f1 : xk * hl + xl * hk = (hk + hl) * (s * a) := by
      linear_combination (-(xl - xk)) * ht - (hk + hl) * e1
    have f2 : (a - xk) * hl + (a - xl) * hk = (hk + hl) * ((1 - s) * a) := by
      linear_combination (-1 : ℝ) * f1
    rw [f1, f2]
    exact ⟨by positivity, mul_pos hs (mul_pos (by linarith) ha)⟩

/-- **`flipConvex` holds exactly when the quadrilateral `i, l, j, k` is strictly convex**: all
four turns of the counter-clockwise quadrilateral are positive (`a, h_k, h_l > 0`).

[formula: it-flip-convexity] -/
theorem flipConvex_iff_turns {a xk hk xl hl : ℝ} (ha : 0 < a) (hhk : 0 < hk) (hhl : 0 < hl) :
    flipConvex a xk hk xl hl = true ↔
      0 < sda (0, 0) (xl, -hl) (a, 0) ∧ 0 < sda (xl, -hl) (a, 0) (xk, hk) ∧
        0 < sda (a, 0) (xk, hk) (0, 0) ∧ 0 < sda (xk, hk) (0, 0) (xl, -hl) := by
  have e1 : sda (0, 0) (xl, -hl) (a, 0) = a * hl := by
    unfold sda
    dsimp only
    ring
  have e2 : sda (xl, -hl) (a, 0) (xk, hk) = (a - xk) * hl + (a - xl) * hk := by
    unfold sda
    dsimp only
    ring
  have e3 : sda (a, 0) (xk, hk) (0, 0) = a * hk := by
    unfold sda
    dsimp only
    ring
  have e4 : sda (xk, hk) (0, 0) (xl, -hl) = xk * hl + xl * hk := by
    unfold sda
    dsimp only
    ring
  rw [flipConvex_iff, flipArea1_eq, flipArea2_eq, e1, e2, e3, e4]
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨mul_pos ha hhl, h2, mul_pos ha hhk, h1⟩
  · rintro ⟨-, h2, -, h1⟩
    exact ⟨h1, h2⟩

end Intrinsic

end DDG
