import DDG.NormalCoords.Frames.Angle
import DDG.Intrinsic.Strip
import DDG.Intrinsic.Signpost

/-!
# Input-aligned frames at an input vertex (`ic-input-aligned-frame`, part 2 of 4)

`crates/ddg-integer-coords/SPEC.md` §4.12 (I26); the Rust `IntegerTriangulation::frame_directions`,
`along_direction` (case `a ∈ V⁰`) and `segment_direction` (case `ChordEndPlace::Corner`),
`crates/ddg-integer-coords/src/frames.rs`; the frame of IT-SPEC N2, N9 (`InputSurface::new`,
`crates/ddg-intrinsic/src/input.rs`).

**The signposts' frame at an input vertex** (`Star`). Around an input vertex `a` the input corners
`θ₀, θ₁, …` follow each other counter-clockwise from the reference halfedge (signpost `0`; at a
boundary vertex the boundary halfedge that leaves `a`); `Θ = Σθ`, `s = angleScale = 2π/Θ`
(`π/Θ` on the boundary). A tangent direction at `a` is a point of the cone of angle `Θ`: its
polar angle `ψ ∈ [0, Θ]` from the reference halfedge is `θ₀ + … + θ_{j−1} + α` for a direction
inside input corner `j` at the angle `α ∈ [0, θ_j]` from the input halfedge `j` that starts the
corner, measured in the layout of the corner's input face (the input surface is flat inside a
face). The signposts store `s ψ`, reduced to `[0, 2π)` inside (`frameAngle`); the input
halfedge `j` itself has `inSign j = s (θ₀ + … + θ_{j−1})` (`InputSurface::new`: `running * scale`),
which is the fold of Alg. 4's update from the reference (`inSign_eq_fold`, by
`signpost_fold_no_wrap` and `signpost_fold_boundary` of `it-signpost-update`). The polar angle
does not depend on how a direction on the seam between two corners is written
(`pre_succ`: end of corner `j` = start of corner `j + 1`; `frameAngle_full_turn`: the end of the
last corner is the reference again).

**Indices.** `j` is the position in the signposts' fan (from the reference halfedge); the
integers number the input halfedges by their local index (SPEC §2.2), another counter-clockwise
numbering. The Rust code reads `signpost(out)` and the corner angle of the halfedge `out` itself,
so only the halfedge matters: the local index `wedgeIndex r` names the halfedge `out`
(`Frames/Wedge.lean`), and `j` is its position in the fan; the corner of `out` is the one
counter-clockwise after it in both numberings.

**What the Rust code computes** for `h = a → b` (`rustCorner`, `rustAlong`):
* `n ≥ 0`: `φ = inSign(x) + s · clamp(angle_from(x, d), 0, θ_x)`, then `wrap` inside and
  `min(·, π)` on the boundary, with `x` the input halfedge of the wedge the roundabout names
  (`Frames/Wedge.lean`: it is the input direction just before `h`, so `h` leaves into corner `j`
  of `x`) and `d` the direction of the edge's straight segment at `a` in the layout of that
  corner's input face;
* `n = −1`: `inSign(g_r)`, or `π` for the arriving boundary edge (`r = deg₀ − 1`).

**Results.**
* `frame_input_vertex` (**`ic-input-aligned-frame`**, `n ≥ 0`): if `d` points into the corner
  (`d = λx + μy`, `λ, μ ≥ 0`, `d ≠ 0`, with `x`, `y` the two sides of the corner in its
  counter-clockwise layout and `θ_j = angleFrom x y`, the corner angle of the lengths,
  `cornerAngle_eq_angleFrom`), then `α = angle_from(x, d)` is the angle that turns `x` into `d`,
  `0 ≤ α ≤ θ_j` (so the clamp and the `min(·, π)` do nothing), and the Rust value is
  `frameAngle (θ₀ + … + θ_{j−1} + α)`: the signposts' direction of `d`, in `[0, 2π)` inside and
  `[0, π]` on the boundary.
* `cone_of_hit`, `cone_of_bary`: the hypothesis on `d` is what the strip gives. If the edge
  crosses an input edge (`n ≥ 1`), its segment from the corner `p` crosses the side opposite `p`
  (tail `b = p + x` below, head `c = p + y` above, `Chord::hit`) at a positive segment parameter
  (`ic-strip-intersection`, `Strip.hit_point`, `Strip.hit_mem`); if it crosses none (`n = 0`), its
  head lies in the corner's face (barycentric coordinates `≥ 0`).
* `rustAlong_eq` (`n = −1`): the Rust value is `frameAngle` of the polar angle of the input
  halfedge `r` (`π` for the arriving boundary edge, whose polar angle is `Θ`).

What the hypotheses of `frame_input_vertex` take from outside: that the intrinsic edge leaves `a`
along the first piece of the straight segment of its unfolded strip (`it-straightest-geodesic`,
cited: a geodesic through no vertex is the straight line of the unfolding of the faces it
crosses), and that the strip is the one the integers decide (`ic-normal-isotopy`, cited; the
crossings' input edges by `ic-extract-intrinsic-edge`).
-/

noncomputable section

namespace DDG.IntegerCoords

namespace Frames

open Real

/-- The star of an input vertex as the signposts see it: the input corner angles counter-clockwise
from the reference halfedge, and whether the vertex is on the boundary. -/
structure Star where
  θs : List ℝ
  bdry : Bool

namespace Star

variable (S : Star)

/-- The angle sum `Θ`. -/
def Θ : ℝ := S.θs.sum

/-- `s = 2π/Θ` (`π/Θ` on the boundary), the oracle's (and the Rust) `angle_scale`. -/
def scale : ℝ := Intrinsic.angleScale S.bdry S.Θ

/-- Input corner `j`. -/
def θ (j : ℕ) : ℝ := S.θs.getD j 0

/-- The corners before input halfedge `j`: its polar angle. -/
def pre (j : ℕ) : ℝ := (S.θs.take j).sum

/-- `InputSurface::signpost` of the `j`-th outgoing input halfedge: `running * scale`. -/
def inSign (j : ℕ) : ℝ := S.pre j * S.scale

/-- **The signposts' frame**: the stored direction of the tangent direction of polar angle `ψ`,
`s ψ`, reduced to `[0, 2π)` at an interior vertex. -/
def frameAngle (ψ : ℝ) : ℝ := if S.bdry then S.scale * ψ else wrap (S.scale * ψ)

/-- `segment_direction`, case `Corner`: corner `j`, raw angle `α = angle_from(x, d)`. -/
def rustCorner (j : ℕ) (α : ℝ) : ℝ :=
  if S.bdry then min (S.inSign j + S.scale * max 0 (min α (S.θ j))) π
  else wrap (S.inSign j + S.scale * max 0 (min α (S.θ j)))

/-- `along_direction`, case `a ∈ V⁰`: the input halfedge `r` (`π` for the arriving boundary edge,
`r = deg₀ − 1`, which is the number of corners at a boundary vertex). -/
def rustAlong (r : ℕ) : ℝ := if S.bdry ∧ r = S.θs.length then π else S.inSign r

/-- The corners are positive (a nondegenerate input mesh) and there is one. -/
structure Ok : Prop where
  pos : ∀ θ ∈ S.θs, 0 < θ
  ne : S.θs ≠ []

variable {S}

lemma Ok.Θ_pos (h : S.Ok) : 0 < S.Θ := List.sum_pos _ h.pos h.ne

lemma Ok.scale_pos (h : S.Ok) : 0 < S.scale := Intrinsic.angleScale_pos _ h.Θ_pos

lemma pre_nonneg (h : S.Ok) (j : ℕ) : 0 ≤ S.pre j :=
  List.sum_nonneg fun x hx => (h.pos x (List.mem_of_mem_take hx)).le

/-- **The seam between two corners**: the end of corner `j` is the start of corner `j + 1`. -/
theorem pre_succ {j : ℕ} (hj : j < S.θs.length) : S.pre (j + 1) = S.pre j + S.θ j := by
  simp only [pre, θ, List.getD_eq_getElem _ _ hj]
  exact List.sum_take_succ _ _ hj

lemma pre_le (h : S.Ok) (j : ℕ) : S.pre j ≤ S.Θ := by
  have e := List.sum_take_add_sum_drop S.θs j
  have : 0 ≤ (S.θs.drop j).sum :=
    List.sum_nonneg fun x hx => (h.pos x (List.mem_of_mem_drop hx)).le
  simp only [pre, Θ]; linarith

lemma pre_lt (h : S.Ok) {j : ℕ} (hj : j < S.θs.length) : S.pre j < S.Θ := by
  have e := List.sum_take_add_sum_drop S.θs j
  have : 0 < (S.θs.drop j).sum :=
    List.sum_pos _ (fun x hx => h.pos x (List.mem_of_mem_drop hx))
      (by rw [Ne, List.drop_eq_nil_iff]; omega)
  simp only [pre, Θ]; linarith

lemma pre_length : S.pre S.θs.length = S.Θ := by simp [pre, Θ]

lemma scale_mul_Θ (h : S.Ok) : S.scale * S.Θ = if S.bdry then π else 2 * π := by
  have := h.Θ_pos.ne'
  cases hb : S.bdry
  · simp only [scale, hb, Intrinsic.angleScale_interior]; field_simp; simp
  · simp only [scale, hb, Intrinsic.angleScale_boundary]; field_simp; simp

/-- **A full turn is the reference again**: at an interior vertex the polar angle `Θ` (the end of
the last corner) has the frame angle `0` of the reference halfedge. -/
theorem frameAngle_full_turn (h : S.Ok) (hb : S.bdry = false) :
    S.frameAngle S.Θ = 0 ∧ S.frameAngle 0 = 0 := by
  have e := scale_mul_Θ h
  simp only [hb] at e
  simp only [frameAngle, hb, e, mul_zero]
  exact ⟨wrap_two_pi, wrap_eq_self le_rfl Real.two_pi_pos⟩

/-- **The input signposts are Alg. 4's fold** from the reference halfedge (`φ = 0`) over the
corners before it (`it-signpost-update`), for every outgoing input halfedge `j < deg`.
[formula: ic-input-aligned-frame] -/
theorem inSign_eq_fold (h : S.Ok) {j : ℕ} (hj : j < S.θs.length) :
    S.inSign j = (S.θs.take j).foldl (fun φ θ => Intrinsic.signpostUpdate S.bdry φ θ S.Θ) 0 := by
  cases hb : S.bdry
  · have hlt := pre_lt h hj
    have hΘ := h.Θ_pos
    have hsm : 2 * π / S.Θ * S.pre j < 2 * π := by
      rw [div_mul_eq_mul_div, div_lt_iff₀ hΘ]
      nlinarith [Real.pi_pos]
    rw [Intrinsic.signpost_fold_no_wrap hΘ _ 0 le_rfl
      (fun x hx => (h.pos x (List.mem_of_mem_take hx)).le) (by simpa [pre] using hsm)]
    simp only [inSign, scale, hb, Intrinsic.angleScale_interior, pre]; ring
  · rw [Intrinsic.signpost_fold_boundary]
    simp only [inSign, scale, hb, Intrinsic.angleScale_boundary, pre]; ring

/-- **`n = −1` at an input vertex** (`along_direction`): the Rust value is the signposts' frame
angle of the input halfedge `r` (`r < deg₀` inside; on the boundary `r ≤` the number of corners,
`r = deg₀ − 1` being the arriving boundary edge, polar angle `Θ`, direction `π`).
[formula: ic-input-aligned-frame] -/
theorem rustAlong_eq (h : S.Ok) {r : ℕ} (hr : r < S.θs.length ∨ (S.bdry = true ∧ r = S.θs.length)) :
    S.rustAlong r = S.frameAngle (S.pre r) := by
  have hs := h.scale_pos
  cases hb : S.bdry
  · have hr' : r < S.θs.length := by simpa [hb] using hr
    have hlt := pre_lt h hr'
    have e := scale_mul_Θ h
    simp only [hb, Bool.false_eq_true, ite_false] at e
    have h0 := pre_nonneg h r
    simp only [rustAlong, frameAngle, hb, Bool.false_eq_true, false_and, ite_false, inSign]
    rw [wrap_eq_self (by positivity) (by nlinarith)]; ring
  · by_cases hl : r = S.θs.length
    · subst hl
      have e := scale_mul_Θ h
      simp only [hb, ite_true] at e
      simp only [rustAlong, frameAngle, hb, and_self, ite_true, pre_length, e]
    · simp only [rustAlong, frameAngle, hb, true_and, hl, ite_false, ite_true, inSign]; ring

/-- **The arithmetic of the corner case**: for `0 ≤ α ≤ θ_j` in a corner `j` of the star, the Rust
value (`inSign j + s · clamp(α, 0, θ_j)`, wrapped inside, `min(·, π)` on the boundary) is the
signposts' frame angle of the polar angle `θ₀ + … + θ_{j−1} + α`; it lies in `[0, 2π)` inside
and in `[0, π]` on the boundary. [formula: ic-input-aligned-frame] -/
theorem rustCorner_eq (h : S.Ok) {j : ℕ} (hj : j < S.θs.length) {α : ℝ} (h0 : 0 ≤ α)
    (h1 : α ≤ S.θ j) :
    S.rustCorner j α = S.frameAngle (S.pre j + α) ∧
      0 ≤ S.rustCorner j α ∧
      (if S.bdry then S.rustCorner j α ≤ π else S.rustCorner j α < 2 * π) := by
  have hc : max 0 (min α (S.θ j)) = α := by rw [min_eq_left h1, max_eq_right h0]
  have hs := h.scale_pos
  have hle : S.pre j + α ≤ S.Θ := by
    have := pre_le h (j + 1); rw [pre_succ hj] at this; linarith
  have hp := pre_nonneg h j
  have e := scale_mul_Θ h
  have hsle : S.scale * (S.pre j + α) ≤ S.scale * S.Θ := mul_le_mul_of_nonneg_left hle hs.le
  have hsum : S.inSign j + S.scale * α = S.scale * (S.pre j + α) := by
    simp only [inSign]; ring
  cases hb : S.bdry
  · simp only [rustCorner, frameAngle, hb, Bool.false_eq_true, ite_false, hc, hsum]
    exact ⟨trivial, (wrap_mem _).1, (wrap_mem _).2⟩
  · simp only [hb, ite_true] at e
    have : S.scale * (S.pre j + α) ≤ π := by linarith
    simp only [rustCorner, frameAngle, hb, ite_true, hc, hsum, min_eq_left this]
    exact ⟨trivial, by positivity, this⟩

/-! ## The geometry of the corner -/

/-- **The corner angle of the lengths is `angleFrom`**: in a counter-clockwise layout of a
triangle with the corner at `p`, sides `x = b − p`, `y = c − p` of lengths `ℓ_x`, `ℓ_y` and
opposite side `‖y − x‖ = ℓ_o`, the corner angle (`cornerAngle`, law of cosines, which the Rust
`corner_angle` evaluates in Kahan's form, `corner-angle`) is `angleFrom x y`.
[formula: ic-input-aligned-frame] -/
theorem cornerAngle_eq_angleFrom {x y : ℝ × ℝ} (hxy : 0 < cross x y) {lo lx ly : ℝ}
    (hx : ‖toC x‖ = lx) (hy : ‖toC y‖ = ly) (ho : ‖toC (y - x)‖ = lo) :
    cornerAngle lo lx ly = angleFrom x y := by
  have hx0 : x ≠ 0 := by rintro rfl; simp [cross] at hxy
  have hy0 : y ≠ 0 := by rintro rfl; simp [cross] at hxy
  rw [angleFrom_eq_arccos hx0 hy0 hxy.le, cornerAngle, hx, hy]
  congr 1
  have nx : lx ^ 2 = x.1 ^ 2 + x.2 ^ 2 := by
    rw [← hx, Complex.sq_norm, Complex.normSq_apply]; simp [toC]; ring
  have ny : ly ^ 2 = y.1 ^ 2 + y.2 ^ 2 := by
    rw [← hy, Complex.sq_norm, Complex.normSq_apply]; simp [toC]; ring
  have no : lo ^ 2 = (y.1 - x.1) ^ 2 + (y.2 - x.2) ^ 2 := by
    rw [← ho, Complex.sq_norm, Complex.normSq_apply]; simp [toC]; ring
  have hlx : lx ≠ 0 := by rw [← hx]; exact norm_ne_zero_iff.mpr (toC_ne_zero hx0)
  have hly : ly ≠ 0 := by rw [← hy]; exact norm_ne_zero_iff.mpr (toC_ne_zero hy0)
  rw [nx, ny, no]
  field_simp
  simp only [dot]; ring

/-- **A segment that crosses the side opposite the corner** (`Chord::hit` on the first exit, tail
`b` strictly below the segment's line, head `c` strictly above, at a positive segment parameter;
`ic-strip-intersection`) **leaves the corner into it**: `d = λ(b − p) + μ(c − p)`, `λ, μ ≥ 0`.
At the head of the edge the same holds for the reversed segment and the last crossing
(`Strip.hit_reverse`). [formula: ic-input-aligned-frame] -/
theorem cone_of_hit {p d b c : ℝ × ℝ} {len : ℝ} (hlen : 0 < len) (hd : len ^ 2 = d.1 ^ 2 + d.2 ^ 2)
    (ha : (Strip.frame p d len b).2 < 0) (hb : 0 < (Strip.frame p d len c).2)
    (hs : 0 < Strip.hitS len (Strip.frame p d len b) (Strip.frame p d len c)) :
    ∃ l m : ℝ, 0 ≤ l ∧ 0 ≤ m ∧
      d = (l * (b - p).1 + m * (c - p).1, l * (b - p).2 + m * (c - p).2) := by
  have hne : (Strip.frame p d len b).2 ≠ (Strip.frame p d len c).2 := by intro e; linarith
  obtain ⟨e1, e2⟩ := Strip.hit_point hlen hd hne
  obtain ⟨t0, t1⟩ := Strip.hit_mem ha hb
  set t := Strip.hitT (Strip.frame p d len b) (Strip.frame p d len c)
  set s := Strip.hitS len (Strip.frame p d len b) (Strip.frame p d len c)
  refine ⟨(1 - t) / s, t / s, div_nonneg (by linarith) hs.le, div_nonneg t0.le hs.le, ?_⟩
  simp only [Prod.fst_sub, Prod.snd_sub]
  ext
  · field_simp; linarith
  · field_simp; linarith

/-- **A segment that ends in the corner's face** (no crossing, `n = 0`): the head
`q = p + d = β₀ p + β₁ b + β₂ c` with `β ≥ 0`, `Σβ = 1`, so `d = β₁(b − p) + β₂(c − p)`. -/
theorem cone_of_bary {p d b c : ℝ × ℝ} {β₀ β₁ β₂ : ℝ} (hsum : β₀ + β₁ + β₂ = 1)
    (hq : (p.1 + d.1, p.2 + d.2) = (β₀ * p.1 + β₁ * b.1 + β₂ * c.1, β₀ * p.2 + β₁ * b.2 + β₂ * c.2)) :
    d = (β₁ * (b - p).1 + β₂ * (c - p).1, β₁ * (b - p).2 + β₂ * (c - p).2) := by
  have h1 := congrArg Prod.fst hq
  have h2 := congrArg Prod.snd hq
  simp only at h1 h2
  simp only [Prod.fst_sub, Prod.snd_sub]
  have : β₀ = 1 - β₁ - β₂ := by linarith
  subst this
  ext
  · simp only; linarith
  · simp only; linarith

/-- **`ic-input-aligned-frame` at an input vertex, `n ≥ 0`.** Corner `j` of the star is laid out
counter-clockwise with sides `x` (the input halfedge `j`, which the roundabout names:
`Frames/Wedge.lean`) and `y`, its angle `θ_j = angleFrom x y` (the corner angle of the lengths,
`cornerAngle_eq_angleFrom`), and the edge's segment leaves the corner into it (`d = λx + μy`,
`λ, μ ≥ 0`, `d ≠ 0`: `cone_of_hit`, `cone_of_bary`). Then `α = angle_from(x, d)` is the angle
that turns `x` into `d` (the only one in `(−π, π]`, `angleFrom_unique`), it lies in `[0, θ_j]`
(the Rust clamp does nothing, SPEC I5: zero clamps in exact arithmetic), and the Rust value
`segment_direction` is the signposts' frame angle of the direction `d`, polar angle
`θ₀ + … + θ_{j−1} + α`: in `[0, 2π)` inside, in `[0, π]` on the boundary (where the `min(·, π)`
does nothing). [formula: ic-input-aligned-frame] -/
theorem frame_input_vertex (h : S.Ok) {j : ℕ} (hj : j < S.θs.length) {x y d : ℝ × ℝ}
    (hxy : 0 < cross x y) (hθ : S.θ j = angleFrom x y) {l m : ℝ} (hl : 0 ≤ l) (hm : 0 ≤ m)
    (hd : d = (l * x.1 + m * y.1, l * x.2 + m * y.2)) (hne : d ≠ 0) :
    (‖toC x‖ : ℂ) * toC d = ‖toC d‖ * Complex.exp (angleFrom x d * Complex.I) * toC x ∧
      0 ≤ angleFrom x d ∧ angleFrom x d ≤ S.θ j ∧
      S.rustCorner j (angleFrom x d) = S.frameAngle (S.pre j + angleFrom x d) ∧
      0 ≤ S.rustCorner j (angleFrom x d) ∧
      (if S.bdry then S.rustCorner j (angleFrom x d) ≤ π
        else S.rustCorner j (angleFrom x d) < 2 * π) := by
  obtain ⟨a0, a1⟩ := angleFrom_le_of_cone hxy hl hm hd hne
  rw [← hθ] at a1
  exact ⟨angleFrom_spec x d, a0, a1, rustCorner_eq h hj a0 a1⟩

end Star

end Frames

end DDG.IntegerCoords
