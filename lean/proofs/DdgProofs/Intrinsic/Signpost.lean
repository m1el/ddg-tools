import DdgProofs.Intrinsic.FlipAngles

/-!
# Intrinsic triangulations: angle normalization and the signpost update

`crates/ddg-intrinsic/SPEC.md` §3 (angle normalization), §4 (Init, Flip), N2; the paper Sec. 2.3
(Eq. 1, 2), Sec. 3.2.1 (Eq. 3), Alg. 1, 2, 4.

* **`it-angle-normalization`** — `angleScale` is `2π/Θ` at an interior vertex and
  `π/Θ` at a boundary vertex (`angleScale_interior`, `angleScale_boundary`); the normalized
  corner angles around a vertex sum to `2π`, respectively `π`
  (`sum_normalized_interior`, `sum_normalized_boundary`); a vertex with `Θ = 2π` (inserted in
  the interior) or `Θ = π` (inserted on the boundary) has the scale `1`; the scale is undone by
  `Θ/2π` (Alg. 2, `angleScale_undo`). That a flip keeps `Θ` at its four vertices is
  `flip_angle_sums` in `DdgProofs/Intrinsic/FlipAngles.lean`.
* **`it-signpost-update`** — `signpostUpdate` at an interior vertex is
  `φ_ij + (2π/Θ) θ` reduced once by `2π` (`signpostUpdate_interior`): for a signpost in
  `[0, 2π)` and a corner `0 ≤ θ ≤ Θ` the result lies in `[0, 2π)` and is the same direction
  (`signpostUpdate_mem`, `signpostUpdate_direction`). Updates compose: two corners in a row
  give the update by their sum (`signpostUpdate_add`), a list of corners the update by its sum
  (`signpost_fold`), and all corners around the vertex return to the start
  (`signpost_around_interior`). From the reference halfedge (`φ = 0`) no reduction happens
  before the turn is complete: the signposts are the cumulative sums of Eq. 2
  (`signpost_fold_no_wrap`). At a boundary vertex nothing is reduced, and the corners from the
  first boundary halfedge to the last add up to the direction `π` (`signpost_fold_boundary`,
  `signpost_boundary_end`). After a flip, the update through the new edge agrees with the old
  one (`flip_signpost_at_k`, `flip_signpost_at_l`).

Definitions (the closed forms of the Rust code, on `ℝ`): `angleScale`, `reduceFullTurn`,
`signpostUpdate`. The fold over a list of corners is `List.foldl` of the update (the loop of
Alg. 4).
-/

noncomputable section

namespace DDG

namespace Intrinsic

open Real

/-! ## The formulas -/

/-- Scale that makes the angles around a vertex of total angle `Θ` add up to `2π` (interior) or
`π` (boundary). -/
def angleScale (isBoundary : Bool) (total : ℝ) : ℝ :=
  if isBoundary then π / total else 2 * π / total

/-- One subtraction of a full turn: for `x ∈ [0, 4π)` the representative in `[0, 2π)`. -/
def reduceFullTurn (x : ℝ) : ℝ := if 2 * π ≤ x then x - 2 * π else x

/-- Signpost of the next halfedge around a vertex: `φ_ik = φ_ij + s θ_i^{jk}`, `s` the scale of
the vertex, reduced to `[0, 2π)` at an interior vertex (not reduced at a boundary vertex). -/
def signpostUpdate (isBoundary : Bool) (phi theta total : ℝ) : ℝ :=
  if isBoundary then phi + angleScale isBoundary total * theta
  else reduceFullTurn (phi + angleScale isBoundary total * theta)

theorem angleScale_interior (Θ : ℝ) : angleScale false Θ = 2 * π / Θ := by
  simp [angleScale]

theorem angleScale_boundary (Θ : ℝ) : angleScale true Θ = π / Θ := by
  simp [angleScale]

theorem reduceFullTurn_eq (x : ℝ) :
    reduceFullTurn x = if 2 * π ≤ x then x - 2 * π else x := by
  rfl

theorem signpostUpdate_interior (φ θ Θ : ℝ) :
    signpostUpdate false φ θ Θ = reduceFullTurn (φ + 2 * π / Θ * θ) := by
  simp [signpostUpdate, angleScale_interior]

theorem signpostUpdate_boundary (φ θ Θ : ℝ) :
    signpostUpdate true φ θ Θ = φ + π / Θ * θ := by
  simp [signpostUpdate, angleScale_boundary]

/-! ## Angle normalization -/

/-- **The normalized angles around an interior vertex sum to `2π`**: with `Θ` the sum of the
corner angles `θ_n` (`Θ ≠ 0`), `Σ_n (2π/Θ) θ_n = 2π`.

[formula: it-angle-normalization] -/
theorem sum_normalized_interior {ι : Type*} (s : Finset ι) (θ : ι → ℝ)
    (hΘ : ∑ n ∈ s, θ n ≠ 0) :
    ∑ n ∈ s, angleScale false (∑ m ∈ s, θ m) * θ n = 2 * π := by
  rw [← Finset.mul_sum, angleScale_interior, div_mul_cancel₀ _ hΘ]

/-- **At a boundary vertex they sum to `π`** (SPEC N2).

[formula: it-angle-normalization] -/
theorem sum_normalized_boundary {ι : Type*} (s : Finset ι) (θ : ι → ℝ)
    (hΘ : ∑ n ∈ s, θ n ≠ 0) :
    ∑ n ∈ s, angleScale true (∑ m ∈ s, θ m) * θ n = π := by
  rw [← Finset.mul_sum, angleScale_boundary, div_mul_cancel₀ _ hΘ]

/-- A vertex inserted in the interior (`Θ = 2π`) and one inserted on the boundary (`Θ = π`)
have the scale `1`. -/
theorem angleScale_flat : angleScale false (2 * π) = 1 ∧ angleScale true π = 1 := by
  rw [angleScale_interior, angleScale_boundary]
  exact ⟨div_self (by positivity), div_self pi_ne_zero⟩

/-- The scale is positive for a positive angle sum. -/
theorem angleScale_pos (b : Bool) {Θ : ℝ} (hΘ : 0 < Θ) : 0 < angleScale b Θ := by
  cases b
  · rw [angleScale_interior]
    positivity
  · rw [angleScale_boundary]
    positivity

/-- A normalized corner angle lies in `[0, 2π]` (interior), for `0 ≤ θ ≤ Θ`. -/
theorem normalized_mem {Θ θ : ℝ} (hΘ : 0 < Θ) (h0 : 0 ≤ θ) (h1 : θ ≤ Θ) :
    0 ≤ 2 * π / Θ * θ ∧ 2 * π / Θ * θ ≤ 2 * π := by
  refine ⟨by positivity, ?_⟩
  rw [div_mul_eq_mul_div, div_le_iff₀ hΘ]
  nlinarith [pi_pos]

/-- **The scale is undone by `Θ/2π`** (the paper's Alg. 2, line 7): `Θ/2π · (2π/Θ · θ) = θ`.

[formula: it-angle-normalization] -/
theorem angleScale_undo {Θ : ℝ} (hΘ : Θ ≠ 0) (θ : ℝ) :
    Θ / (2 * π) * (angleScale false Θ * θ) = θ := by
  rw [angleScale_interior]
  have := pi_ne_zero
  field_simp

/-! ## The signpost update at an interior vertex -/

/-- One reduction brings `[0, 4π)` into `[0, 2π)`, by subtracting `0` or `2π`. -/
theorem reduceFullTurn_mem {x : ℝ} (h0 : 0 ≤ x) (h4 : x < 4 * π) :
    0 ≤ reduceFullTurn x ∧ reduceFullTurn x < 2 * π ∧
      (reduceFullTurn x = x ∨ reduceFullTurn x = x - 2 * π) := by
  rw [reduceFullTurn_eq]
  split_ifs with h
  · exact ⟨by linarith, by linarith, Or.inr rfl⟩
  · exact ⟨h0, not_le.mp h, Or.inl rfl⟩

/-- **The updated signpost lies in `[0, 2π)`** for a signpost in `[0, 2π)` and a corner
`0 ≤ θ ≤ Θ`, and differs from `φ_ij + (2π/Θ) θ` by `0` or `2π`.

[formula: it-signpost-update] -/
theorem signpostUpdate_mem {φ θ Θ : ℝ} (hφ0 : 0 ≤ φ) (hφ : φ < 2 * π) (hΘ : 0 < Θ)
    (h0 : 0 ≤ θ) (h1 : θ ≤ Θ) :
    0 ≤ signpostUpdate false φ θ Θ ∧ signpostUpdate false φ θ Θ < 2 * π ∧
      (signpostUpdate false φ θ Θ = φ + 2 * π / Θ * θ ∨
        signpostUpdate false φ θ Θ = φ + 2 * π / Θ * θ - 2 * π) := by
  obtain ⟨n0, n1⟩ := normalized_mem hΘ h0 h1
  rw [signpostUpdate_interior]
  exact reduceFullTurn_mem (by linarith) (by linarith)

/-- **The update gives the direction of the next halfedge** (the paper's Eq. 3): as a unit
vector, the updated signpost is `φ_ij + (2π/Θ) θ_i^{jk}`.

[formula: it-signpost-update] -/
theorem signpostUpdate_direction {φ θ Θ : ℝ} (hφ0 : 0 ≤ φ) (hφ : φ < 2 * π) (hΘ : 0 < Θ)
    (h0 : 0 ≤ θ) (h1 : θ ≤ Θ) :
    cos (signpostUpdate false φ θ Θ) = cos (φ + 2 * π / Θ * θ) ∧
      sin (signpostUpdate false φ θ Θ) = sin (φ + 2 * π / Θ * θ) := by
  obtain ⟨-, -, h | h⟩ := signpostUpdate_mem hφ0 hφ hΘ h0 h1
  · rw [h]
    exact ⟨rfl, rfl⟩
  · rw [h]
    exact ⟨Real.cos_sub_two_pi _, Real.sin_sub_two_pi _⟩

/-- **Two updates in a row are the update by the sum of the two corners** (`θ₁ + θ₂ ≤ Θ`,
`φ < 2π`).

[formula: it-signpost-update] -/
theorem signpostUpdate_add {φ θ₁ θ₂ Θ : ℝ} (hφ : φ < 2 * π) (hΘ : 0 < Θ)
    (h1 : 0 ≤ θ₁) (h2 : 0 ≤ θ₂) (hs : θ₁ + θ₂ ≤ Θ) :
    signpostUpdate false (signpostUpdate false φ θ₁ Θ) θ₂ Θ =
      signpostUpdate false φ (θ₁ + θ₂) Θ := by
  have hx1 : 0 ≤ 2 * π / Θ * θ₁ := by positivity
  have hx2 : 0 ≤ 2 * π / Θ * θ₂ := by positivity
  have hx := (normalized_mem hΘ (by linarith : 0 ≤ θ₁ + θ₂) hs).2
  have hsum : 2 * π / Θ * (θ₁ + θ₂) = 2 * π / Θ * θ₁ + 2 * π / Θ * θ₂ := by ring
  rw [hsum] at hx
  simp only [signpostUpdate_interior, reduceFullTurn_eq, hsum]
  split_ifs <;> linarith

/-- An update by the corner `0` changes nothing. -/
theorem signpostUpdate_zero {φ : ℝ} (hφ : φ < 2 * π) (Θ : ℝ) :
    signpostUpdate false φ 0 Θ = φ := by
  rw [signpostUpdate_interior, reduceFullTurn_eq, mul_zero, add_zero, ite_eq_right (not_le.mpr hφ)]

/-- **A full turn returns to the start**: the update by the whole angle `Θ` is the identity on
`[0, 2π)`.

[formula: it-signpost-update] -/
theorem signpostUpdate_full_turn {φ Θ : ℝ} (hφ0 : 0 ≤ φ) (hΘ : Θ ≠ 0) :
    signpostUpdate false φ Θ Θ = φ := by
  have e : 2 * π / Θ * Θ = 2 * π := div_mul_cancel₀ _ hΘ
  rw [signpostUpdate_interior, reduceFullTurn_eq, e, ite_eq_left (by linarith)]
  ring

/-- **Updating over a list of corners is the update by their sum** (corners `≥ 0`, sum `≤ Θ`,
start in `[0, 2π)`). -/
theorem signpost_fold {Θ : ℝ} (hΘ : 0 < Θ) :
    ∀ (θs : List ℝ) (φ : ℝ), 0 ≤ φ → φ < 2 * π → (∀ θ ∈ θs, 0 ≤ θ) → θs.sum ≤ Θ →
      θs.foldl (fun φ θ => signpostUpdate false φ θ Θ) φ =
        signpostUpdate false φ θs.sum Θ
  | [], φ, _, hφ, _, _ => by
    rw [List.foldl_nil, List.sum_nil, signpostUpdate_zero hφ]
  | θ :: rest, φ, hφ0, hφ, hpos, hsum => by
    have hθ : 0 ≤ θ := hpos θ (List.mem_cons_self ..)
    have hrest : ∀ x ∈ rest, 0 ≤ x := fun x hx => hpos x (List.mem_cons_of_mem _ hx)
    have hrs : 0 ≤ rest.sum := List.sum_nonneg hrest
    rw [List.sum_cons] at hsum
    obtain ⟨m0, m1, -⟩ := signpostUpdate_mem hφ0 hφ hΘ hθ (by linarith)
    rw [List.foldl_cons, signpost_fold hΘ rest _ m0 m1 hrest (by linarith), List.sum_cons,
      signpostUpdate_add hφ hΘ hθ hrs hsum]

/-- **Going once around an interior vertex returns to the start**: updating the signpost
`φ ∈ [0, 2π)` over all corners at the vertex (each `≥ 0`, their sum `Θ`) gives `φ` again.

[formula: it-signpost-update] -/
theorem signpost_around_interior {Θ : ℝ} (hΘ : 0 < Θ) {θs : List ℝ} (hpos : ∀ θ ∈ θs, 0 ≤ θ)
    (hsum : θs.sum = Θ) {φ : ℝ} (hφ0 : 0 ≤ φ) (hφ : φ < 2 * π) :
    θs.foldl (fun φ θ => signpostUpdate false φ θ Θ) φ = φ := by
  rw [signpost_fold hΘ θs φ hφ0 hφ hpos hsum.le, hsum, signpostUpdate_full_turn hφ0 hΘ.ne']

/-- **Before the turn is complete nothing is reduced** (the paper's Eq. 2): starting from
`φ ≥ 0`, if `φ + (2π/Θ) Σθ < 2π` the signpost after the corners `θs` is `φ + (2π/Θ) Σθ`; from
the reference halfedge (`φ = 0`) the signposts are the cumulative sums of the normalized
angles. -/
theorem signpost_fold_no_wrap {Θ : ℝ} (hΘ : 0 < Θ) :
    ∀ (θs : List ℝ) (φ : ℝ), 0 ≤ φ → (∀ θ ∈ θs, 0 ≤ θ) → φ + 2 * π / Θ * θs.sum < 2 * π →
      θs.foldl (fun φ θ => signpostUpdate false φ θ Θ) φ = φ + 2 * π / Θ * θs.sum
  | [], φ, _, _, _ => by
    rw [List.foldl_nil, List.sum_nil, mul_zero, add_zero]
  | θ :: rest, φ, hφ0, hpos, hlt => by
    have hθ : 0 ≤ θ := hpos θ (List.mem_cons_self ..)
    have hrest : ∀ x ∈ rest, 0 ≤ x := fun x hx => hpos x (List.mem_cons_of_mem _ hx)
    have hrs : 0 ≤ 2 * π / Θ * rest.sum := mul_nonneg (by positivity) (List.sum_nonneg hrest)
    have hx : 0 ≤ 2 * π / Θ * θ := by positivity
    rw [List.sum_cons, mul_add] at hlt
    have e : signpostUpdate false φ θ Θ = φ + 2 * π / Θ * θ := by
      rw [signpostUpdate_interior, reduceFullTurn_eq, ite_eq_right (not_le.mpr (by linarith))]
    rw [List.foldl_cons, e, signpost_fold_no_wrap hΘ rest _ (by linarith) hrest (by linarith),
      List.sum_cons]
    ring

/-! ## The signpost update at a boundary vertex -/

/-- **At a boundary vertex nothing is reduced**: after the corners `θs` the signpost is
`φ + (π/Θ) Σθ`. -/
theorem signpost_fold_boundary (Θ : ℝ) :
    ∀ (θs : List ℝ) (φ : ℝ),
      θs.foldl (fun φ θ => signpostUpdate true φ θ Θ) φ = φ + π / Θ * θs.sum
  | [], φ => by rw [List.foldl_nil, List.sum_nil, mul_zero, add_zero]
  | θ :: rest, φ => by
    rw [List.foldl_cons, signpostUpdate_boundary, signpost_fold_boundary Θ rest, List.sum_cons]
    ring

/-- **The boundary edge that arrives has the direction `π`** (SPEC N2): from the boundary
halfedge that leaves (`φ = 0`) over all corners at the vertex (sum `Θ ≠ 0`).

[formula: it-signpost-update] -/
theorem signpost_boundary_end {Θ : ℝ} (hΘ : Θ ≠ 0) {θs : List ℝ} (hsum : θs.sum = Θ) :
    θs.foldl (fun φ θ => signpostUpdate true φ θ Θ) 0 = π := by
  rw [signpost_fold_boundary, hsum, zero_add, div_mul_cancel₀ _ hΘ]

/-! ## A flip and the signposts (the paper's Alg. 5) -/

/-- **After a flip the signposts at `k` stay consistent**: updating across the new edge `kl`
(by the corner of `jkl` at `k`, then by that of `ilk`) gives the same signpost as the old
update by the corner of `ijk` at `k`. So the signposts of the halfedges that the flip does not
touch need no change. -/
theorem flip_signpost_at_k {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) {φ Θ : ℝ} (hφ : φ < 2 * π) (hΘ : 0 < Θ)
    (hle : cornerAngle a b c ≤ Θ) :
    signpostUpdate false
        (signpostUpdate false φ (cornerAngle e (flippedLengthFromLengths a b c d e) c) Θ)
        (cornerAngle d (flippedLengthFromLengths a b c d e) b) Θ =
      signpostUpdate false φ (cornerAngle a b c) Θ := by
  obtain ⟨-, -, h3, -⟩ := flip_angle_sums hk hl hconv
  obtain ⟨t1, t2⟩ := flip_isTriangle hk hl hconv
  have p1 : 0 ≤ cornerAngle e (flippedLengthFromLengths a b c d e) c := by
    rw [cornerAngle_swap e c]
    exact (cornerAngle_pos t2.rotate).le
  have p2 : 0 ≤ cornerAngle d (flippedLengthFromLengths a b c d e) b :=
    (cornerAngle_pos t1.rotate.rotate).le
  rw [signpostUpdate_add hφ hΘ p1 p2 (by linarith), h3]

/-- The same at `l`: through the new edge `lk` (corners of `jkl` and `ilk` at `l`). -/
theorem flip_signpost_at_l {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) {φ Θ : ℝ} (hφ : φ < 2 * π) (hΘ : 0 < Θ)
    (hle : cornerAngle a d e ≤ Θ) :
    signpostUpdate false
        (signpostUpdate false φ (cornerAngle c e (flippedLengthFromLengths a b c d e)) Θ)
        (cornerAngle b d (flippedLengthFromLengths a b c d e)) Θ =
      signpostUpdate false φ (cornerAngle a d e) Θ := by
  obtain ⟨-, -, -, h4⟩ := flip_angle_sums hk hl hconv
  obtain ⟨t1, t2⟩ := flip_isTriangle hk hl hconv
  have p1 : 0 ≤ cornerAngle c e (flippedLengthFromLengths a b c d e) := by
    rw [cornerAngle_swap c _ e]
    exact (cornerAngle_pos t2.rotate.rotate).le
  have p2 : 0 ≤ cornerAngle b d (flippedLengthFromLengths a b c d e) :=
    (cornerAngle_pos t1.rotate).le
  rw [signpostUpdate_add hφ hΘ p1 p2 (by linarith), h4]

end Intrinsic

end DDG
