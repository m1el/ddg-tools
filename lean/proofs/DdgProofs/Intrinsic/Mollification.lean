import DdgProofs.Triangle

/-!
# Intrinsic mollification of edge lengths

Transcription of `crates/ddg-bff/src/mollify.rs` into definitions on real numbers.

## `intrinsic-mollification` — `mollify_edge_lengths`

* `margin_iff_mollifyDeficit_le` — a face satisfies the inequalities with margin at `ε` iff
  its deficit is `≤ ε`;
* `mollifyNeeded_le_iff` — `needed ≤ ε` iff `0 ≤ ε` and the face satisfies them at `ε`;
* `mollifyEpsilonLoop_le_iff` — the loop returns the least `ε ≥ 0` at which all faces do;
* `mollifyEpsilonLoop_eq`, `mollifyEpsilonLoop_ofFn` — it is `ε* = max(0, max deficit)`;
* `mollifyEpsilon_le_iff`, `mollifyEpsilon_nonneg`, `mollifyEpsilon_feasible`,
  `mollifyEpsilon_le_of_feasible` — `ε* ≥ 0`, `ε*` is admissible, and every admissible
  `ε ≥ 0` is `≥ ε*`;
* `mollifyEpsilon_eq_zero_iff` — `ε* = 0` iff all faces already satisfy the inequalities;
* `isTriangle_of_margin`, `isTriangle_mollified` — for `δ > 0` the mollified sides of every
  face satisfy the strict triangle inequalities (nothing else is needed).

Nothing is claimed about floating-point rounding.

Nothing is claimed about floating-point rounding.
-/

noncomputable section

namespace DDG

open Real

/-! ## Intrinsic mollification -/

/-- The triangle inequalities of a face with sides `a, b, c` hold with margin `δ` after every
side is lengthened by `ε`: `(x+ε) + (y+ε) ≥ (z+ε) + δ` for each side `z`. -/
def MarginHolds (a b c δ ε : ℝ) : Prop :=
  (b + ε) + (c + ε) ≥ (a + ε) + δ ∧ (c + ε) + (a + ε) ≥ (b + ε) + δ ∧
    (a + ε) + (b + ε) ≥ (c + ε) + δ

/-- By how much a face misses the margin: the largest of `ℓ_k − ℓ_{k+1} − ℓ_{k+2} + δ`. -/
def mollifyDeficit (a b c δ : ℝ) : ℝ :=
  max (a - b - c + δ) (max (b - c - a + δ) (c - a - b + δ))

/-- `needed` as the loop over the three corners computes it: starting from `0`, the maximum
with `ℓ_k − ℓ_{k+1} − ℓ_{k+2} + δ` for `k = 0, 1, 2`. -/
def mollifyNeeded (a b c δ : ℝ) : ℝ :=
  max (max (max 0 (a - b - c + δ)) (b - c - a + δ)) (c - a - b + δ)

/-- One pass of the loop over the faces: `if needed > 0 { if needed > ε { ε = needed } }`. -/
def mollifyStep (δ ε : ℝ) (l : ℝ × ℝ × ℝ) : ℝ :=
  if 0 < mollifyNeeded l.1 l.2.1 l.2.2 δ then
    if ε < mollifyNeeded l.1 l.2.1 l.2.2 δ then mollifyNeeded l.1 l.2.1 l.2.2 δ else ε
  else ε

/-- `epsilon` as the loop over the faces computes it, starting from `0`. A face is the triple
of its side lengths. -/
def mollifyEpsilonLoop (δ : ℝ) (faces : List (ℝ × ℝ × ℝ)) : ℝ := faces.foldl (mollifyStep δ) 0

/-- **The margin inequalities of a face hold at `ε` iff its deficit is at most `ε`** (for every
real `ε`).

[formula: intrinsic-mollification] -/
theorem margin_iff_mollifyDeficit_le (a b c δ ε : ℝ) :
    MarginHolds a b c δ ε ↔ mollifyDeficit a b c δ ≤ ε := by
  unfold MarginHolds mollifyDeficit
  rw [max_le_iff, max_le_iff]
  constructor
  · rintro ⟨h1, h2, h3⟩
    exact ⟨by linarith, by linarith, by linarith⟩
  · rintro ⟨h1, h2, h3⟩
    exact ⟨by linarith, by linarith, by linarith⟩

/-- The loop over the corners starts from `0`: `needed = max(0, deficit)`. -/
theorem mollifyNeeded_eq (a b c δ : ℝ) :
    mollifyNeeded a b c δ = max 0 (mollifyDeficit a b c δ) := by
  unfold mollifyNeeded mollifyDeficit
  rw [max_assoc, max_assoc]

/-- **`needed ≤ ε` iff `ε ≥ 0` and the margin inequalities of the face hold at `ε`.**

[formula: intrinsic-mollification] -/
theorem mollifyNeeded_le_iff (a b c δ ε : ℝ) :
    mollifyNeeded a b c δ ≤ ε ↔ 0 ≤ ε ∧ MarginHolds a b c δ ε := by
  rw [mollifyNeeded_eq, max_le_iff, margin_iff_mollifyDeficit_le]

theorem mollifyNeeded_nonneg (a b c δ : ℝ) : 0 ≤ mollifyNeeded a b c δ := by
  rw [mollifyNeeded_eq]
  exact le_max_left _ _

private theorem mollifyStep_nonneg {δ ε₀ : ℝ} (h0 : 0 ≤ ε₀) (l : ℝ × ℝ × ℝ) :
    0 ≤ mollifyStep δ ε₀ l := by
  unfold mollifyStep
  split_ifs with h1 h2
  · exact h1.le
  · exact h0
  · exact h0

private theorem mollifyStep_le_iff {δ ε₀ : ℝ} (h0 : 0 ≤ ε₀) (l : ℝ × ℝ × ℝ) (ε : ℝ) :
    mollifyStep δ ε₀ l ≤ ε ↔ ε₀ ≤ ε ∧ mollifyNeeded l.1 l.2.1 l.2.2 δ ≤ ε := by
  have hn := mollifyNeeded_nonneg l.1 l.2.1 l.2.2 δ
  unfold mollifyStep
  split_ifs with h1 h2
  · exact ⟨fun h => ⟨by linarith, h⟩, fun h => h.2⟩
  · exact ⟨fun h => ⟨h, by linarith⟩, fun h => h.1⟩
  · exact ⟨fun h => ⟨h, by linarith⟩, fun h => h.1⟩

private theorem foldl_mollifyStep_le_iff (δ : ℝ) (faces : List (ℝ × ℝ × ℝ)) :
    ∀ ε₀ : ℝ, 0 ≤ ε₀ → ∀ ε : ℝ, faces.foldl (mollifyStep δ) ε₀ ≤ ε ↔
      ε₀ ≤ ε ∧ ∀ l ∈ faces, mollifyNeeded l.1 l.2.1 l.2.2 δ ≤ ε := by
  induction faces with
  | nil =>
    intro ε₀ _ ε
    simp
  | cons l ls ih =>
    intro ε₀ h0 ε
    rw [List.foldl_cons, ih _ (mollifyStep_nonneg h0 l), mollifyStep_le_iff h0,
      List.forall_mem_cons, and_assoc]

/-- **The loop computes the least admissible `ε`**: `epsilon ≤ ε` iff `ε ≥ 0` and the margin
inequalities of all faces hold at `ε`. So `epsilon` is itself `≥ 0` and admissible (take
`ε = epsilon`), and every admissible `ε ≥ 0` is at least `epsilon`.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilonLoop_le_iff (δ : ℝ) (faces : List (ℝ × ℝ × ℝ)) (ε : ℝ) :
    mollifyEpsilonLoop δ faces ≤ ε ↔
      0 ≤ ε ∧ ∀ l ∈ faces, MarginHolds l.1 l.2.1 l.2.2 δ ε := by
  unfold mollifyEpsilonLoop
  rw [foldl_mollifyStep_le_iff δ faces 0 le_rfl]
  constructor
  · rintro ⟨h0, h⟩
    exact ⟨h0, fun l hl => ((mollifyNeeded_le_iff _ _ _ _ _).mp (h l hl)).2⟩
  · rintro ⟨h0, h⟩
    exact ⟨h0, fun l hl => (mollifyNeeded_le_iff _ _ _ _ _).mpr ⟨h0, h l hl⟩⟩

section Family

variable {F : Type*} [Fintype F] [Nonempty F]

/-- `ε* = max(0, max over the faces of the deficit)` for a finite nonempty family of faces
`l f = (a, b, c)`: the closed form of what `mollify_edge_lengths` adds to every edge. -/
def mollifyEpsilon (l : F → ℝ × ℝ × ℝ) (δ : ℝ) : ℝ :=
  max 0 (Finset.univ.sup' Finset.univ_nonempty fun f =>
    mollifyDeficit (l f).1 (l f).2.1 (l f).2.2 δ)

/-- **`ε* ≤ ε` iff `ε ≥ 0` and all faces satisfy the inequalities with margin at `ε`**: `ε*` is
the least such `ε`.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_le_iff (l : F → ℝ × ℝ × ℝ) (δ ε : ℝ) :
    mollifyEpsilon l δ ≤ ε ↔ 0 ≤ ε ∧ ∀ f, MarginHolds (l f).1 (l f).2.1 (l f).2.2 δ ε := by
  unfold mollifyEpsilon
  rw [max_le_iff, Finset.sup'_le_iff]
  simp only [Finset.mem_univ, true_implies, margin_iff_mollifyDeficit_le]

/-- **The loop computes `ε*`**, in whatever order (and however often) it visits the faces.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilonLoop_eq (l : F → ℝ × ℝ × ℝ) (δ : ℝ) {faces : List (ℝ × ℝ × ℝ)}
    (hfaces : ∀ x, x ∈ faces ↔ ∃ f, l f = x) :
    mollifyEpsilonLoop δ faces = mollifyEpsilon l δ := by
  refine eq_of_forall_ge_iff fun ε => ?_
  rw [mollifyEpsilonLoop_le_iff, mollifyEpsilon_le_iff]
  constructor
  · rintro ⟨h0, h⟩
    exact ⟨h0, fun f => h _ ((hfaces _).mpr ⟨f, rfl⟩)⟩
  · rintro ⟨h0, h⟩
    refine ⟨h0, fun x hx => ?_⟩
    obtain ⟨f, rfl⟩ := (hfaces x).mp hx
    exact h f

/-- The loop over the faces `0, …, n` in this order computes `ε*`. (For no faces at all the
loop returns `0`.)

[formula: intrinsic-mollification] -/
theorem mollifyEpsilonLoop_ofFn {n : ℕ} (l : Fin (n + 1) → ℝ × ℝ × ℝ) (δ : ℝ) :
    mollifyEpsilonLoop δ (List.ofFn l) = mollifyEpsilon l δ :=
  mollifyEpsilonLoop_eq l δ fun _ => List.mem_ofFn

/-- **`ε* ≥ 0`.**

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_nonneg (l : F → ℝ × ℝ × ℝ) (δ : ℝ) : 0 ≤ mollifyEpsilon l δ :=
  ((mollifyEpsilon_le_iff l δ _).mp le_rfl).1

/-- **`ε*` is admissible**: after lengthening every side by `ε*`, every face satisfies its three
triangle inequalities with margin `δ`.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_feasible (l : F → ℝ × ℝ × ℝ) (δ : ℝ) (f : F) :
    MarginHolds (l f).1 (l f).2.1 (l f).2.2 δ (mollifyEpsilon l δ) :=
  ((mollifyEpsilon_le_iff l δ _).mp le_rfl).2 f

/-- **`ε*` is the least admissible amount**: every `ε ≥ 0` for which all faces satisfy the
inequalities with margin is at least `ε*`.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_le_of_feasible (l : F → ℝ × ℝ × ℝ) (δ : ℝ) {ε : ℝ} (hε : 0 ≤ ε)
    (h : ∀ f, MarginHolds (l f).1 (l f).2.1 (l f).2.2 δ ε) : mollifyEpsilon l δ ≤ ε :=
  (mollifyEpsilon_le_iff l δ ε).mpr ⟨hε, h⟩

/-- The lengths are left unchanged (`ε* = 0`) exactly when all faces already satisfy the
inequalities with margin.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_eq_zero_iff (l : F → ℝ × ℝ × ℝ) (δ : ℝ) :
    mollifyEpsilon l δ = 0 ↔ ∀ f, MarginHolds (l f).1 (l f).2.1 (l f).2.2 δ 0 := by
  constructor
  · intro h
    exact ((mollifyEpsilon_le_iff l δ 0).mp h.le).2
  · intro h
    exact le_antisymm ((mollifyEpsilon_le_iff l δ 0).mpr ⟨le_rfl, h⟩) (mollifyEpsilon_nonneg l δ)

end Family

/-- **A positive margin gives a nondegenerate triangle.** If the margin inequalities hold at
`ε` with `δ > 0`, the lengthened sides satisfy the strict triangle inequalities, and each of
them is at least `δ`. Nothing else is needed: no sign condition on `a, b, c` or `ε`.

[formula: intrinsic-mollification] -/
theorem isTriangle_of_margin {a b c δ ε : ℝ} (hδ : 0 < δ) (h : MarginHolds a b c δ ε) :
    IsTriangle (a + ε) (b + ε) (c + ε) ∧ δ ≤ a + ε ∧ δ ≤ b + ε ∧ δ ≤ c + ε := by
  obtain ⟨h1, h2, h3⟩ := h
  exact ⟨⟨by linarith, by linarith, by linarith⟩, by linarith, by linarith, by linarith⟩

/-- **The mollified metric is a metric of nondegenerate triangles** whenever `δ > 0`: for every
face, the sides lengthened by `ε*` satisfy `IsTriangle` (so they are positive and the face has
positive area). For `δ = 0` only the weak inequalities hold, and a degenerate face stays
degenerate (see the example below).

[formula: intrinsic-mollification] -/
theorem isTriangle_mollified {F : Type*} [Fintype F] [Nonempty F] (l : F → ℝ × ℝ × ℝ) {δ : ℝ}
    (hδ : 0 < δ) (f : F) :
    IsTriangle ((l f).1 + mollifyEpsilon l δ) ((l f).2.1 + mollifyEpsilon l δ)
        ((l f).2.2 + mollifyEpsilon l δ) ∧
      0 < triArea ((l f).1 + mollifyEpsilon l δ) ((l f).2.1 + mollifyEpsilon l δ)
        ((l f).2.2 + mollifyEpsilon l δ) := by
  have h := (isTriangle_of_margin hδ (mollifyEpsilon_feasible l δ f)).1
  exact ⟨h, triArea_pos h⟩

/-- `δ > 0` cannot be dropped: the degenerate face `(1, 1, 2)` needs no lengthening for
`δ = 0` and is not a triangle. -/
example : mollifyNeeded 1 1 2 0 = 0 ∧ ¬IsTriangle ((1 : ℝ) + 0) (1 + 0) (2 + 0) := by
  refine ⟨by norm_num [mollifyNeeded], fun h => ?_⟩
  have := h.lt_c
  norm_num at this

end DDG
