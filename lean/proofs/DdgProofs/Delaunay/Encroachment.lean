import DDG.Geometry.Triangle

/-!
# Delaunay refinement: encroachment and the circumradius

Ruppert's algorithm splits a boundary segment when a vertex *encroaches* on it, i.e. lies inside
its diametral circle (the circle with the segment as diameter). In a Delaunay triangulation the
test is made on the angle opposite the segment: the vertex is strictly inside the diametral circle
exactly when the angle it subtends is obtuse, on it exactly when the angle is right (Thales).

## Main statements

* `inner_vsub_vsub_eq_dist_midpoint_sq` — for points `a, b, p` of a Euclidean space,
  `⟪a − p, b − p⟫ = |p m|² − (|ab|/2)²` with `m` the midpoint of `ab`.
* `dist_midpoint_lt_iff` — **`p` lies strictly inside the diametral circle of `ab` iff
  `∠ a p b > π/2`**; `dist_midpoint_eq_iff` (on the circle iff the angle is right, the strict
  converse of Thales' theorem, `EuclideanGeometry.Sphere.thales_theorem` in Mathlib),
  `dist_midpoint_gt_iff` (outside iff acute). No hypothesis on `p`: if `p = a` or `p = b` the
  angle is `π/2` (Mathlib's convention) and `p` lies on the circle.
* `pi_div_two_lt_cornerAngle_iff` — the same in lengths: the corner opposite a side `s` of a
  triangle with the other sides `x, y` is obtuse iff `x² + y² < s²`, and
  `encroaches_iff_cornerAngle` — for Euclidean points, `p` is strictly inside the diametral
  circle of `ab` iff the corner angle at `p` of the triangle `a p b`, computed from its three
  lengths, exceeds `π/2`. This is the test the Rust refinement evaluates (on the intrinsic
  triangle, whose layout in the plane is a Euclidean triangle with the same lengths).
* Circumradius (Chew, Ruppert, Shewchuk's radius-edge ratio): `circumradius a b c = abc/4A`
  equals `a/(2 sin θ_a)` (`circumradius_eq_div_sin`, the law of sines), a side is shorter than the
  circumradius iff the opposite angle has `sin θ < 1/2` (`lt_circumradius_iff`), so an angle below
  `π/6` is opposite a side shorter than the circumradius (`lt_circumradius_of_lt_pi_div_six`).

The termination of Ruppert's and Chew's algorithms is not proved here (`CITED.md`,
`it-ruppert-encroachment`, `it-chew-second-algorithm`).

## References

* J. Ruppert, *A Delaunay refinement algorithm for quality 2-dimensional mesh generation*,
  J. Algorithms 18(3), 548–585, 1995.
* J. R. Shewchuk, *Delaunay refinement algorithms for triangular mesh generation*, Comput. Geom.
  22(1–3), 21–74, 2002, Sec. 3.1 (encroachment) and Sec. 3 (radius-edge ratio); *Delaunay
  Refinement Mesh Generation*, Ph.D. thesis, CMU-CS-97-137, 1997, Sec. 5.3.2 (encroachment
  diagnosed by an obtuse angle opposite the segment).
* L. P. Chew, *Guaranteed-quality mesh generation for curved surfaces*, SCG '93, 274–280, 1993.

Dependencies: Mathlib and `DDG/Geometry/Triangle.lean` only.
-/

noncomputable section

namespace DDG

namespace Delaunay

open Real EuclideanGeometry

section Euclidean

variable {V P : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [MetricSpace P]
  [NormedAddTorsor V P]

/-- `⟪a − p, b − p⟫ = |p m|² − (|ab|/2)²`, `m` the midpoint of `ab` (Apollonius). -/
theorem inner_vsub_vsub_eq_dist_midpoint_sq (a b p : P) :
    inner ℝ (a -ᵥ p) (b -ᵥ p) = dist p (midpoint ℝ a b) ^ 2 - (dist a b / 2) ^ 2 := by
  have h1 := dist_sq_add_dist_sq_eq_two_mul_dist_midpoint_sq_add_half_dist_sq p a b
  have h2 := norm_sub_sq_real (a -ᵥ p) (b -ᵥ p)
  rw [vsub_sub_vsub_cancel_right, ← dist_eq_norm_vsub V, ← dist_eq_norm_vsub V,
    ← dist_eq_norm_vsub V] at h2
  rw [dist_comm p a, dist_comm p b] at h1
  linarith

/-- **Encroachment**: `p` lies strictly inside the circle with diameter `ab` iff the angle
`∠ a p b` is obtuse (Shewchuk 1997, Sec. 5.3.2). -/
theorem dist_midpoint_lt_iff (a b p : P) :
    dist p (midpoint ℝ a b) < dist a b / 2 ↔ π / 2 < ∠ a p b := by
  rw [EuclideanGeometry.angle, ← InnerProductGeometry.inner_neg_iff_pi_div_two_lt_angle,
    inner_vsub_vsub_eq_dist_midpoint_sq, sub_neg]
  have h1 := dist_nonneg (x := p) (y := midpoint ℝ a b)
  have h2 : 0 ≤ dist a b / 2 := by have := dist_nonneg (x := a) (y := b); linarith
  constructor
  · intro h; nlinarith
  · intro h; nlinarith

/-- **Thales' theorem and its converse**: `p` lies on the circle with diameter `ab` iff
`∠ a p b = π/2`. -/
theorem dist_midpoint_eq_iff (a b p : P) :
    dist p (midpoint ℝ a b) = dist a b / 2 ↔ ∠ a p b = π / 2 := by
  rw [EuclideanGeometry.angle, ← InnerProductGeometry.inner_eq_zero_iff_angle_eq_pi_div_two,
    inner_vsub_vsub_eq_dist_midpoint_sq, sub_eq_zero]
  have h1 := dist_nonneg (x := p) (y := midpoint ℝ a b)
  have h2 : 0 ≤ dist a b / 2 := by have := dist_nonneg (x := a) (y := b); linarith
  constructor
  · intro h; rw [h]
  · intro h
    have := (sq_eq_sq₀ h1 h2).mp h
    exact this

/-- `p` lies strictly outside the circle with diameter `ab` iff `∠ a p b` is acute. -/
theorem dist_midpoint_gt_iff (a b p : P) :
    dist a b / 2 < dist p (midpoint ℝ a b) ↔ ∠ a p b < π / 2 := by
  have h1 := dist_midpoint_lt_iff a b p
  have h2 := dist_midpoint_eq_iff a b p
  constructor
  · intro h
    rcases lt_trichotomy (∠ a p b) (π / 2) with h' | h' | h'
    · exact h'
    · exact absurd (h2.mpr h') h.ne'
    · exact absurd (h1.mpr h') (not_lt.mpr h.le)
  · intro h
    rcases lt_trichotomy (dist a b / 2) (dist p (midpoint ℝ a b)) with h' | h' | h'
    · exact h'
    · exact absurd (h2.mp h'.symm) h.ne
    · exact absurd (h1.mp h') (not_lt.mpr h.le)

end Euclidean

/-! ## In lengths -/

variable {s x y : ℝ}

/-- **The corner opposite `s` is obtuse iff `x² + y² < s²`.**

[formula: it-encroachment-angle] -/
theorem pi_div_two_lt_cornerAngle_iff (h : IsTriangle s x y) :
    π / 2 < cornerAngle s x y ↔ x ^ 2 + y ^ 2 < s ^ 2 := by
  have hx := h.pos_b
  have hy := h.pos_c
  rw [cornerAngle, Real.pi_div_two_lt_arccos, div_neg_iff]
  constructor
  · rintro (⟨-, h2⟩ | ⟨h1, -⟩)
    · exact absurd h2 (not_lt.mpr (by positivity))
    · linarith
  · intro h'
    right
    exact ⟨by linarith, by positivity⟩

/-- **The test of the refinement in lengths**: for points `a, b, p` of a Euclidean space with
`p ≠ a, b`, `p` lies strictly inside the diametral circle of `ab` iff the corner angle at `p`,
computed from the three distances, exceeds `π/2`.

[formula: it-encroachment-angle] -/
theorem encroaches_iff_cornerAngle {V P : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    [MetricSpace P] [NormedAddTorsor V P] {a b p : P} (ha : a ≠ p) (hb : b ≠ p) :
    dist p (midpoint ℝ a b) < dist a b / 2 ↔
      π / 2 < cornerAngle (dist a b) (dist a p) (dist b p) := by
  rw [dist_midpoint_lt_iff, angle_eq_cornerAngle ha hb]

/-! ## The circumradius -/

/-- The circumradius of a triangle from its lengths: `R = abc/(4A)`. -/
def circumradius (a b c : ℝ) : ℝ := a * b * c / (4 * triArea a b c)

variable {a b c : ℝ}

theorem circumradius_rotate (a b c : ℝ) : circumradius b c a = circumradius a b c := by
  unfold circumradius
  rw [triArea_rotate]
  ring

theorem circumradius_pos (h : IsTriangle a b c) : 0 < circumradius a b c := by
  have := h.pos_a
  have := h.pos_b
  have := h.pos_c
  have := triArea_pos h
  unfold circumradius
  positivity

/-- **The law of sines**: `R = a/(2 sin θ_a)`, `θ_a` the angle opposite `a`. -/
theorem circumradius_eq_div_sin (h : IsTriangle a b c) :
    circumradius a b c = a / (2 * sin (cornerAngle a b c)) := by
  have := h.pos_b
  have := h.pos_c
  have := (triArea_pos h).ne'
  rw [sin_cornerAngle h, circumradius]
  field_simp
  ring

/-- A side is shorter than the circumradius iff the sine of the opposite angle is below `1/2`. -/
theorem lt_circumradius_iff (h : IsTriangle a b c) :
    a < circumradius a b c ↔ 2 * sin (cornerAngle a b c) < 1 := by
  have ha := h.pos_a
  have hs : 0 < 2 * sin (cornerAngle a b c) :=
    mul_pos two_pos (sin_pos_of_pos_of_lt_pi (cornerAngle_pos h) (cornerAngle_lt_pi h))
  rw [circumradius_eq_div_sin h, lt_div_iff₀ hs]
  constructor
  · intro h'; nlinarith
  · intro h'; nlinarith

/-- **An angle below `π/6` is opposite a side shorter than the circumradius** (the radius-edge
ratio of such a triangle exceeds `1`; Chew 1993, Shewchuk 2002 Sec. 3). -/
theorem lt_circumradius_of_lt_pi_div_six (h : IsTriangle a b c)
    (hθ : cornerAngle a b c < π / 6) : a < circumradius a b c := by
  rw [lt_circumradius_iff h]
  have h0 := cornerAngle_pos h
  have : sin (cornerAngle a b c) < sin (π / 6) :=
    sin_lt_sin_of_lt_of_le_pi_div_two (by linarith [pi_pos]) (by linarith [pi_pos]) hθ
  rw [sin_pi_div_six] at this
  linarith

end Delaunay

end DDG
