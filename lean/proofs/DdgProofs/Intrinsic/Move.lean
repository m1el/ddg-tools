import DdgProofs.Intrinsic.Diamond

/-!
# Intrinsic triangulations: moving a vertex, the circumcentre vector

`crates/ddg-intrinsic/SPEC.md` §4 (Move a vertex, Optimal Delaunay triangulation), N19, N20, N21;
the paper Sec. 3.3.3, Alg. 10, Sec. 4.3; L. Chen, M. Holst, *Efficient mesh optimization schemes
based on Optimal Delaunay Triangulations*, CMAME 200 (2011), Lemma 4.2 (Eq. 4.3–4.6), Eq. 4.13,
Sec. 4.2.1 of the preprint (barycentres at the boundary).

Plane points are pairs `ℝ × ℝ` with `sda` (doubled signed area) and `sqDist2`, `dist2` of
`DdgProofs/Bff/Quality.lean`.

* **`it-move-length`** — `movedLength l r (cos α) (sin α)` is the paper's law of
  cosines `√(r² + ℓ² − 2 r ℓ cos α)` (`movedLength_eq_law_of_cosines`) and the distance between
  the moved end `r (cos φ, sin φ)` and the far end `ℓ (cos ψ, sin ψ)` when `α = ψ − φ`
  (`movedLength_eq_dist2`).
* **`it-move-area`** — `movedDoubleArea (sda 0 j k) v (k − j)` is the doubled signed area of
  the triangle `(v, j, k)` (`movedDoubleArea_eq_sda`); around a closed link the doubled signed
  areas of the star add up to a number that does not depend on the centre
  (`sum_sda_closed_link`): moving the vertex keeps the area of its star.
* **`it-circumcenter-vector`** — `weightedCircumcenter ℓ_ab ℓ_ac ℓ_bc A / A` is equidistant
  from `a = 0`, `b = (ℓ_ab, 0)` and `c = (c_x, h)` (`weightedCircumcenter_equidistant`), also
  stated with `A` the Heron area of the three lengths and `c` the apex `apexK`
  (`weightedCircumcenter_equidistant_heron`).

The optimal Delaunay triangulation update (`it-odt-update`: the barycentre, the circumcentre
vector `circumVec`, the star energy and the step of Eq. 4.13) is in `DdgProofs/Delaunay/ODT.lean`.

Definitions (the closed forms of the Rust code, on `ℝ`): `movedLength`, `movedDoubleArea`,
`weightedCircumcenter` (the plane helpers `rot`, `nsq`, `toDual` are in
`DdgProofs/Geometry/PlaneLemmas.lean`).
-/

noncomputable section

namespace DDG

namespace Intrinsic

open Real Finset

/-! ## `it-move-length` -/

/-- Length of the edge `ij` after its end `i` moves by a vector of length `r` that makes the
angle `α` with the edge, given by `cosA = cos α` and `sinA = sin α`:
`ℓ' = √((ℓ − r cos α)² + (r sin α)²)`. -/
def movedLength (l r cosA sinA : ℝ) : ℝ :=
  √((l - r * cosA) * (l - r * cosA) + (r * sinA) * (r * sinA))

theorem movedLength_eq (l r c s : ℝ) :
    movedLength l r c s = √((l - r * c) ^ 2 + (r * s) ^ 2) := by
  unfold movedLength
  congr 1
  ring

/-- **The length of a moved edge is the paper's law of cosines** (Alg. 10):
`movedLength ℓ r (cos α) (sin α) = √(r² + ℓ² − 2 r ℓ cos α)`.

[formula: it-move-length] -/
theorem movedLength_eq_law_of_cosines (l r α : ℝ) :
    movedLength l r (cos α) (sin α) = √(r ^ 2 + l ^ 2 - 2 * r * l * cos α) := by
  rw [movedLength_eq]
  congr 1
  linear_combination r ^ 2 * sin_sq_add_cos_sq α

/-- **… and the distance from the moved end to the far end**: the vertex moves from the origin
to `r (cos φ, sin φ)`, the far end of the edge is at `ℓ (cos ψ, sin ψ)`; with `α = ψ − φ` (the
angle from the motion to the edge) the new length is their distance.

[formula: it-move-length] -/
theorem movedLength_eq_dist2 (l r φ ψ : ℝ) :
    movedLength l r (cos (ψ - φ)) (sin (ψ - φ)) =
      dist2 (r * cos φ, r * sin φ) (l * cos ψ, l * sin ψ) := by
  rw [movedLength_eq]
  unfold dist2
  congr 1
  dsimp only
  rw [cos_sub, sin_sub]
  linear_combination (r ^ 2 * (sin φ ^ 2 + cos φ ^ 2) - l ^ 2) * sin_sq_add_cos_sq ψ

/-- The length does not change for a zero motion. -/
theorem movedLength_zero {l : ℝ} (hl : 0 ≤ l) (c s : ℝ) : movedLength l 0 c s = l := by
  rw [movedLength_eq]
  simp [sqrt_sq hl]

/-! ## `it-move-area` -/

/-- Twice the signed area of a triangle of the star of a vertex after the vertex (at the origin)
moves by `v`: `2A − v × d`, `A2 = 2A` twice the area before and `d = k − j` the side of the
link. -/
def movedDoubleArea (A2 : ℝ) (v d : ℝ × ℝ) : ℝ := A2 - (v.1 * d.2 - v.2 * d.1)

theorem movedDoubleArea_eq (A2 : ℝ) (v d : ℝ × ℝ) :
    movedDoubleArea A2 v d = A2 - (v.1 * d.2 - v.2 * d.1) := rfl

/-- **The doubled signed area after the move**: for the triangle `(0, j, k)` of the star, the
vertex at the origin moved to `v`, `sda v j k = movedDoubleArea (sda 0 j k) v (k − j)`.

[formula: it-move-area] -/
theorem movedDoubleArea_eq_sda (v j k : ℝ × ℝ) :
    sda v j k = movedDoubleArea (sda (0, 0) j k) v (k - j) := by
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

/-- The area `A` of the triangle `abc` times the vector from `a` to its circumcentre, in the
frame of `a` with `b = (ℓ_ab, 0)` and `c` above: `(A ℓ_ab/2, ℓ_ab (ℓ_ac² + ℓ_bc² − ℓ_ab²)/8)`. -/
def weightedCircumcenter (lab lac lbc A : ℝ) : ℝ × ℝ :=
  (A * lab / 2, lab * (lac * lac + lbc * lbc - lab * lab) / 8)

theorem weightedCircumcenter_eq (lab lac lbc A : ℝ) :
    weightedCircumcenter lab lac lbc A =
      (A * lab / 2, lab * (lac * lac + lbc * lbc - lab * lab) / 8) := by
  rfl

/-- **The weighted circumcentre divided by the area is the circumcentre**: in the frame of `a`
with `b = (ℓ_ab, 0)` and `c = (c_x, h)`, `c_x = (ℓ_ab² + ℓ_ac² − ℓ_bc²)/2ℓ_ab`, `h = 2A/ℓ_ab`,
`A ≠ 0`, `ℓ_ab ≠ 0` and `|c| = ℓ_ac`, the point `o = weightedCircumcenter … / A` has the same
squared distance from `a = 0`, `b` and `c`.

[formula: it-circumcenter-vector] -/
theorem weightedCircumcenter_equidistant {lab lac lbc A cx h : ℝ} (hA : A ≠ 0) (hab : lab ≠ 0)
    (hcx : cx = (lab ^ 2 + lac ^ 2 - lbc ^ 2) / (2 * lab)) (hh : h = 2 * A / lab)
    (hc : cx ^ 2 + h ^ 2 = lac ^ 2) :
    let o : ℝ × ℝ := ((weightedCircumcenter lab lac lbc A).1 / A,
      (weightedCircumcenter lab lac lbc A).2 / A)
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
triangle inequalities) and `A = triArea`, `c` the apex
`(diamondX ℓ_ab ℓ_ac ℓ_bc, diamondH ℓ_ab ℓ_ac ℓ_bc)`, the point
`weightedCircumcenter ℓ_ab ℓ_ac ℓ_bc A / A` is equidistant from `0`, `(ℓ_ab, 0)` and `c`.

[formula: it-circumcenter-vector] -/
theorem weightedCircumcenter_equidistant_heron {lab lac lbc : ℝ} (ht : IsTriangle lab lac lbc) :
    let A := triArea lab lac lbc
    let o : ℝ × ℝ := ((weightedCircumcenter lab lac lbc A).1 / A,
      (weightedCircumcenter lab lac lbc A).2 / A)
    sqDist2 (0, 0) o = sqDist2 (lab, 0) o ∧
      sqDist2 (0, 0) o = sqDist2 (diamondX lab lac lbc, diamondH lab lac lbc) o := by
  intro A o
  have hab : 0 < lab := by linarith [ht.lt_b, ht.lt_c]
  have hA : 0 < A := triArea_pos ht
  obtain ⟨-, -, -, hb, -, -, -⟩ := diamond_isometric ht ht
  have hc : diamondX lab lac lbc ^ 2 + diamondH lab lac lbc ^ 2 = lac ^ 2 := by
    have h0 := dist2_sq' (0, 0) (apexK lab lac lbc)
    rw [hb] at h0
    rw [h0]
    simp [apexK]
  exact weightedCircumcenter_equidistant hA.ne' hab.ne'
    (by rw [diamondX_eq]; ring) (by rw [diamondH_eq]) hc

end Intrinsic

end DDG
