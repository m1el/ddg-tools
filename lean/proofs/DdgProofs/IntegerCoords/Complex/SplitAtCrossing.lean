import DdgProofs.IntegerCoords.SplitAtCrossing
import DdgProofs.IntegerCoords.Complex.Merge

/-!
# Integer coordinates: the faces of a split at a crossing, and the merge back (`ic-split-at-crossing`)

`crates/ddg-integer-coords/SPEC.md` §4.5 ("At a crossing", I28), §4.6 (the merge); the Rust
`formulas::split_side_at_crossing`, `split::split_integers`, `merge::merge_edge_normal`.

A valid face `abc` (side counts `n`, admissible corner counts `c`) is split at crossing `j` of side
`0` (`a→b`) by the new vertex `m` on the arc `A` through that crossing (`SplitAtCrossing.lean`). In
the numbering of `Complex/Merge.lean` (`X2 = a`, `X0 = b`, `X1 = c`) the two new faces are

* `A = (m, b, c)`: side `0` = `m→b` (the half with `n₀ − 1 − j` crossings), `1` = `b→c` (`n₁`),
  `2` = `c→m` (the spoke, `spoke n c j` = the arcs crossing `mc`, `chordCount`);
* `B = (m, c, a)`: side `0` = `m→c` (the spoke), `1` = `c→a` (`n₂`), `2` = `a→m` (the half with `j`
  crossings);

and `mS = 1` (the spoke runs along an input edge, `n_mc = −1`) exactly when `A` comes from `c`. The
crossing `j` itself belongs to neither half: `A` ends at `m`.

* `split_side`: both new faces are valid (explicit admissible corner counts `cA`, `cB`, so Eq. 1 is
  even: `Merge.Side`), the side is straight (`Merge.Straight`, the condition of
  `check_straight_arcs`), there is exactly **one curve end at `m`** on this side (`Merge.ends = 1`:
  `A`), and merging the two faces back (`Merge.mergeN`, `Merge.mergeC`) gives the original face, side
  counts and corner counts: the split is the inverse of the merge, the `+1` case.
* `emanatingBefore_eq_newFace`: the Rust `emanating_before` (the arcs from `c` ending on `ab` before
  `m`) is Eq. 2 of the new face `B` at `c`, the arcs leaving `c` between `c→a` and `c→m`: the
  roundabout `r_cm = r_ca + n⁻_ca + emanating_before` is Eq. 4's form (`ic-roundabout-update`).
* `split_at_crossing_merge` (both sides of the edge, the second face `b a d` with crossing
  `n₀ − 1 − j`): the halves glue (`a→m` and `m→a` have `j`, `m→b` and `b→m` have `n₀ − 1 − j`),
  and `Merge.merge_normal` applies with `τ = 1`, `on_edge`: the Rust `merge_edge_normal(n₀ − 1 − j,
  j, true)` returns `n₀`, and both merged faces are the original ones.
-/

set_option linter.unusedSimpArgs false

namespace DDG.IntegerCoords

open Ddg.IntegerCoords (cornerCoord faceValid nPlus)

namespace SplitAtCrossing

variable (n c : Fin 3 → ℕ) (j : ℕ)

/-- The crossings of the spoke `mc`: the arcs crossing the chord from `m` to `c`. -/
def spoke : ℕ := chordCount n c (j + 1) 2

/-- Side counts of `A = (m, b, c)`. -/
def nA : Fin 3 → ℕ := ![n 0 - 1 - j, n 1, spoke n c j]

/-- Side counts of `B = (m, c, a)`. -/
def nB : Fin 3 → ℕ := ![spoke n c j, n 2, j]

/-- `n⁻` of the spoke: `1` when `A` comes from `c` (the spoke runs along its input edge). -/
def mS : ℕ := if canon n c 0 j = .corner 2 then 1 else 0

/-- Corner counts of `A` (corners `m, b, c`). If `A` cuts `a`: the arcs around `a` beyond `A` cut
`m`, those around `b` cut `b`, and those around `c` and from `a` cut `c`. If `A` comes from `c`:
only the arcs around `b`. If `A` cuts `b`: the arcs around `b` beyond `A` (`n₀ − 1 − j`) cut `b`,
and those around `b` inside `A`, around `c` and from `a` cut `c`. -/
def cA : Fin 3 → ℕ :=
  if j < c 0 then ![c 0 - 1 - j, c 1, c 2 + em n c 0]
  else if j < c 0 + em n c 2 then ![0, c 1, 0]
  else ![0, n 0 - 1 - j, c 2 + em n c 0 + (j - c 0 - em n c 2)]

/-- Corner counts of `B` (corners `m, c, a`), mirrored. -/
def cB : Fin 3 → ℕ :=
  if j < c 0 then ![0, c 0 - 1 - j + c 2 + em n c 1, j]
  else if j < c 0 + em n c 2 then ![0, 0, c 0]
  else ![j - c 0 - em n c 2, c 2, c 0]

/-- The facts about an admissible face that the case analyses use. -/
lemma facts {n c : Fin 3 → ℕ} (hc : Admissible n c) :
    c 0 + c 1 ≤ n 0 ∧ c 1 + c 2 ≤ n 1 ∧ c 2 + c 0 ≤ n 2 ∧
    em n c 0 + c 1 + c 2 = n 1 ∧ em n c 1 + c 2 + c 0 = n 2 ∧ em n c 2 + c 0 + c 1 = n 0 ∧
    (0 < em n c 0 → c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0) ∧
    (0 < em n c 1 → c 1 = 0 ∧ em n c 2 = 0 ∧ em n c 0 = 0) ∧
    (0 < em n c 2 → c 2 = 0 ∧ em n c 0 = 0 ∧ em n c 1 = 0) := by
  obtain ⟨a0, a1, a2, b0, b1, b2, x0, x1, x2⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := hc.em_add
  exact ⟨a0, a1, a2, m0, m1, m2, fun h => ⟨b0 h, x0 h⟩, fun h => ⟨b1 h, x1 h⟩,
    fun h => ⟨b2 h, (x2 h).1, (x2 h).2⟩⟩

lemma em_vec (x y : Fin 3 → ℕ) :
    em x y 0 = x 1 - y 1 - y 2 ∧ em x y 1 = x 2 - y 2 - y 0 ∧ em x y 2 = x 0 - y 0 - y 1 :=
  em_lits x y

/-- Admissibility of a face given by literal vectors, as nine linear facts. -/
lemma admissible_vec {x0 x1 x2 y0 y1 y2 : ℕ}
    (h : y0 + y1 ≤ x0 ∧ y1 + y2 ≤ x1 ∧ y2 + y0 ≤ x2 ∧
      (0 < x1 - y1 - y2 → y0 = 0) ∧ (0 < x2 - y2 - y0 → y1 = 0) ∧ (0 < x0 - y0 - y1 → y2 = 0) ∧
      (0 < x1 - y1 - y2 → x2 - y2 - y0 = 0 ∧ x0 - y0 - y1 = 0) ∧
      (0 < x2 - y2 - y0 → x0 - y0 - y1 = 0 ∧ x1 - y1 - y2 = 0) ∧
      (0 < x0 - y0 - y1 → x1 - y1 - y2 = 0 ∧ x2 - y2 - y0 = 0)) :
    Admissible ![x0, x1, x2] ![y0, y1, y2] := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9⟩ := h
  refine ⟨fun k => ?_, fun k => ?_, fun k => ?_⟩ <;> rcases fin3_cases k with rfl | rfl | rfl <;>
    simp only [em, l1, l2, l3, l4, l5, l6, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] <;>
    assumption

variable {n c j}

/-- **The two new faces are admissible** with the corner counts `cA`, `cB`.
[formula: ic-split-at-crossing] -/
theorem admissible_halves (hc : Admissible n c) (hj : j < n 0) :
    Admissible (nA n c j) (cA n c j) ∧ Admissible (nB n c j) (cB n c j) := by
  have hT := chordCount_closed hc hj
  obtain ⟨a0, a1, a2, m0, m1, m2, e0, e1, e2⟩ := facts hc
  unfold nA nB cA cB spoke
  rw [hT]
  generalize em n c 0 = f0 at *
  generalize em n c 1 = f1 at *
  generalize em n c 2 = f2 at *
  generalize c 0 = x0 at *
  generalize c 1 = x1 at *
  generalize c 2 = x2 at *
  generalize n 0 = n0 at *
  generalize n 1 = n1 at *
  generalize n 2 = n2 at *
  rcases Glue.em_cases e0 e1 e2 with ⟨z1, z2, z3⟩ | ⟨z1, z2, z3⟩ | ⟨z1, z2, z3⟩ | ⟨z1, z2, z3⟩ <;>
    subst z1 z2 z3 <;> split_ifs <;>
    exact ⟨admissible_vec (by omega), admissible_vec (by omega)⟩

/-- A face with admissible corner counts is valid (Eq. 1 even). -/
lemma faceValid_of_admissible {m d : Fin 3 → ℕ} (h : Admissible m d) :
    faceValid (m 0) (m 1) (m 2) = true :=
  (FaceModel.exists_valid_iff m).1 ⟨⟨m, canon m d⟩, rfl, canon_valid h⟩

/-- The corner counts of Eq. 1 are the admissible ones. -/
lemma cc_eq {m d : Fin 3 → ℕ} (h : Admissible m d) : cc m = d := by
  funext k
  simp only [cc, ← h.c_eq k, Int.toNat_natCast]

/-- The case analysis of `split_side`: the emanating corner of the face, the block of the
crossing, then linear arithmetic on the literal vectors. -/
macro "sac_cases" : tactic => `(tactic| (
  simp only [em, fin3_lits.1, fin3_lits.2.1, fin3_lits.2.2.1, fin3_lits.2.2.2.1,
    fin3_lits.2.2.2.2.1, fin3_lits.2.2.2.2.2, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] <;>
  rcases Glue.em_cases ‹0 < _ → _ ∧ _ ∧ _› ‹0 < _ → _ ∧ _ ∧ _› ‹0 < _ → _ ∧ _ ∧ _› with
    ⟨z1, z2, z3⟩ | ⟨z1, z2, z3⟩ | ⟨z1, z2, z3⟩ | ⟨z1, z2, z3⟩ <;>
  split_ifs <;>
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons] <;> omega))

/-- **`ic-split-at-crossing`: one side of the split.** For a valid face (admissible `c`) split at
crossing `j < n₀` of side `0`: the faces `A = (m, b, c)`, `B = (m, c, a)` with the halves `n₀ − 1 − j`
and `j` (the crossing `j` counted by neither), the spoke `spoke n c j` (`n_mc`, or `0` with
`mS = 1` when the spoke runs along the input edge of the arc through `m`) and the old sides `n₁`,
`n₂` unchanged, are valid faces glued along the spoke (`Merge.Side`), straight at `m`
(`Merge.Straight`), with exactly one curve end at `m` (`Merge.ends = 1`), and merging them back
gives the original face: side counts `n` (the merged edge `n₀ − 1 − j + j + 1`) and corner counts
`c`. [formula: ic-split-at-crossing] -/
theorem split_side (hc : Admissible n c) (hj : j < n 0) :
    Merge.Side (nA n c j) (nB n c j) (mS n c j) ∧ Merge.Straight (nA n c j) (nB n c j) ∧
      Merge.ends (nA n c j) (nB n c j) (mS n c j) = 1 ∧
      Merge.mergeN (nA n c j) (nB n c j) (mS n c j) = n ∧
      Merge.mergeC (nA n c j) (nB n c j) = c := by
  obtain ⟨hA, hB⟩ := admissible_halves hc hj
  have eA := cc_eq hA
  have eB := cc_eq hB
  have hfan := canon_side0_eq_corner hc j
  have hT := chordCount_closed hc hj
  obtain ⟨a0, a1, a2, m0, m1, m2, e0, e1, e2⟩ := facts hc
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  refine ⟨⟨faceValid_of_admissible hA, faceValid_of_admissible hB, rfl, ?_, ?_⟩, ?_, ?_, ?_, ?_⟩
  · unfold mS; split_ifs <;> omega
  · intro h
    unfold mS at h
    split_ifs at h with h'
    · exact chordCount_fan hc hj h'
  all_goals
    try funext k
    simp only [Merge.Straight, Merge.ends, Merge.mergeN, Merge.mergeC, Merge.ea, Merge.eb, eA, eB]
    try simp only [mS, hfan]
    unfold nA nB cA cB spoke
    rw [hT]
    generalize em n c 0 = f0 at *
    generalize em n c 1 = f1 at *
    generalize em n c 2 = f2 at *
    generalize c = c0 at *
    generalize n = n0 at *
  · sac_cases
  · sac_cases
  · rcases fin3_cases k with rfl | rfl | rfl <;> sac_cases
  · rcases fin3_cases k with rfl | rfl | rfl <;> sac_cases

/-- **The Rust `emanating_before` is Eq. 2 of the new face `B = (m, c, a)` at `c`**: the arcs from
`c` that end on `ab` before `m` are the arcs leaving `c` in the face between `c→a` and `c→m`, so the
roundabout `r_cm = r_ca + n⁻_ca + emanating_before` is Eq. 4's form (`ic-roundabout-update`).
[formula: ic-split-at-crossing] -/
theorem emanatingBefore_eq_newFace (hc : Admissible n c) (hj : j < n 0) :
    ((Finset.range (em n c 2)).filter (fun s => c 0 + s < j)).card =
      em (nB n c j) (cc (nB n c j)) 1 := by
  rw [emanatingBefore_eq j, cc_eq (admissible_halves hc hj).2]
  obtain ⟨a0, a1, a2, m0, m1, m2, e0, e1, e2⟩ := facts hc
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [em, l1, l2, l3, l4, l5, l6] at m0 m1 m2 e0 e1 e2 ⊢
  unfold nB cB
  simp only [em, l1, l2, l3, l4, l5, l6]
  split_ifs <;>
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      Matrix.head_cons, Matrix.tail_cons] <;> omega

/-- **Both sides of the edge, and the merge back.** The edge `ab` (side `0` of the face `abc`, counts
`n`) is side `0` of the face `b a d` from `b` (counts `n'`, `n'₀ = n₀`); crossing `j` from `a` is
crossing `j' = n₀ − 1 − j` from `b` (`ic-reverse-crossing`). Split on both sides: the halves glue
(`a→m` has `j` crossings on both sides, `m→b` has `n₀ − 1 − j`), each side is straight with one curve
end at `m`, and the Rust merge of the halves `n₁ = n₀ − 1 − j` (`m→b`), `n₂ = j` (`a→m`) with
`on_edge` (`Merge.merge_normal`, `τ = 1`) returns `n₀`, the merged faces being the original ones.
[formula: ic-split-at-crossing] -/
theorem split_at_crossing_merge {n' c' : Fin 3 → ℕ} (hc : Admissible n c) (hc' : Admissible n' c')
    (hn : n' 0 = n 0) (hj : j < n 0) :
    nB n c j 2 = nA n' c' (n 0 - 1 - j) 0 ∧ nA n c j 0 = nB n' c' (n 0 - 1 - j) 2 ∧
      Merge.mergeNormal ((n 0 - 1 - j : ℕ) : ℤ) (j : ℤ) true = some (n 0 : ℤ) ∧
      Merge.mergeN (nA n c j) (nB n c j) (mS n c j) = n ∧
      Merge.mergeN (nA n' c' (n 0 - 1 - j)) (nB n' c' (n 0 - 1 - j)) (mS n' c' (n 0 - 1 - j)) = n' ∧
      Merge.mergeC (nA n c j) (nB n c j) = c ∧
      Merge.mergeC (nA n' c' (n 0 - 1 - j)) (nB n' c' (n 0 - 1 - j)) = c' := by
  have hj' : n 0 - 1 - j < n' 0 := by omega
  obtain ⟨sU, tU, eU, nU, cU⟩ := split_side hc hj
  obtain ⟨sL, tL, eL, nL, cL⟩ := split_side hc' hj'
  have g1 : nB n c j 2 = nA n' c' (n 0 - 1 - j) 0 := by
    simp only [nA, nB, Matrix.cons_val_zero, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons]
    omega
  have g2 : nA n c j 0 = nB n' c' (n 0 - 1 - j) 2 := by
    simp only [nA, nB, Matrix.cons_val_zero, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons]
  refine ⟨g1, g2, ?_, nU, nL, cU, cL⟩
  have np : ∀ m : ℕ, nPlus (m : ℤ) = m := fun m => by simp only [nPlus]; omega
  obtain ⟨-, -, x, hx, hxn, -, -⟩ := Merge.merge_normal sU tU sL tL (n1 := ((n 0 - 1 - j : ℕ) : ℤ))
    (n2 := (j : ℤ)) (by omega) (by omega)
    (by rw [np]; simp only [nA, Matrix.cons_val_zero])
    (by rw [np, ← g2]; simp only [nA, Matrix.cons_val_zero])
    (by rw [np]; simp only [nB, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons])
    (by rw [np, ← g1]; simp only [nB, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons])
    eU eL (onEdge := true) (by omega) (by omega) (fun _ => by simp) le_rfl
  rw [nU] at hxn
  simp only [nPlus] at hxn
  rw [hx]
  congr 1
  omega

end SplitAtCrossing

end DDG.IntegerCoords
