import DDG.NormalCoords.Complex.Quad

/-!
# Integer coordinates, tier C: merging the two halves of a split edge (`ic-merge-edge-normal`)

`crates/ddg-integer-coords/SPEC.md` §4.6 ("The integers of the merge", I27); the Rust
`merge::merge_edge_normal`, `check_straight_arcs` (`crates/ddg-integer-coords/src/merge.rs`).

An inserted vertex `m` has two edges on one geodesic through it, the *halves* `m X0` and `X2 m`.
On one side of them lie the faces `A = (m, X0, X1)` and `B = (m, X1, X2)`, glued along the spoke
`m X1`; merging the halves replaces them by the face `M = (X2, X0, X1)` whose side `X2→X0` is the
union `X2→m` ∪ `m→X0`. In the face-model numbering:

* `A`: corners `m, X0, X1`; side `0` = `m→X0` (the first half), `1` = `X0→X1`, `2` = `X1→m` (the
  spoke);
* `B`: corners `m, X1, X2`; side `0` = `m→X1` (the spoke, the other way), `1` = `X1→X2`, `2` =
  `X2→m` (the second half);
* `M`: corners `X2, X0, X1`; side `0` = `X2→X0` (the merged edge), `1` = `X0→X1`, `2` = `X1→X2`.

`mS = n⁻` of the spoke (`1`: the spoke runs along an input edge). The **ends at `m` on this side**
are the arcs leaving corner `m` in `A` (to `X0X1`), the spoke if it runs along an input edge, and
the arcs leaving corner `m` in `B`: `ends = e_m(A) + mS + e_m(B)`. Merged, each of them crosses the
merged edge at the point `m`.

The **straightness condition** of the Rust merge (`check_straight_arcs`): no arc runs through the
spoke next to `m` from the half `m X0` or the corner `X0` to the half `m X2` or the corner `X2`:
`min(c_m(A) + e_X0(A), c_m(B) + e_X2(B)) = 0`.

* `MergeSide.admissible`: then the merged face, with side counts
  `(n⁺_{X2 m} + n⁺_{m X0} + ends, n_{X0 X1}, n_{X1 X2})`, is valid (its corner data, `mergeC`, is
  admissible), and its arcs are those of `A` and `B` joined through the spoke, the ends at `m`
  becoming crossings of the merged edge: so the merged edge has `n⁺₁ + n⁺₂ + ends` crossings.
* `merge_normal`: with the same number of ends on both sides (`τ`; SPEC §2.4 item 4: an inserted
  vertex on an input edge `E` has exactly two curve ends, one on each side by the second check of
  the merge, none inside a face), the Rust `merge_edge_normal n₁ n₂ on_edge` returns `n` with
  `n⁺ = n⁺₁ + n⁺₂ + τ` (`τ = 1` iff `m` lies on an input edge the halves do not run along), and
  `n = −1` exactly when both halves run along an input edge. [formula: ic-merge-edge-normal]
-/

namespace DDG.IntegerCoords

open DDG.IntegerCoords.Formulas (faceValid nPlus nMinus)

namespace Merge

variable (nA nB : Fin 3 → ℕ) (mS : ℕ)

/-- Arcs leaving the corners of `A` and of `B`. -/
def ea (k : Fin 3) : ℕ := em nA (cc nA) k
def eb (k : Fin 3) : ℕ := em nB (cc nB) k

/-- The curve ends at `m` on this side. -/
def ends : ℕ := ea nA 0 + mS + eb nB 0

/-- Side counts of the merged face `(X2, X0, X1)`. -/
def mergeN : Fin 3 → ℕ := ![nB 2 + nA 0 + ends nA nB mS, nA 1, nB 1]

/-- Corner counts of the merged face: at `X2` the arcs cutting `X2` in `B`, the arcs from the half
`m X0` through the spoke (cutting `m` in `A`) and the ends at `m` towards `X1X2`; at `X0` the same
mirrored; at `X1` the arcs around `X1` crossing the spoke. -/
def mergeC : Fin 3 → ℕ :=
  ![cc nB 2 + cc nA 0 + eb nB 0, cc nA 1 + cc nB 0 + ea nA 0, min (cc nA 2) (cc nB 1)]

/-- Two valid faces glued along the spoke, with `mS = n⁻` of the spoke. -/
structure Side : Prop where
  vA : faceValid (nA 0) (nA 1) (nA 2) = true
  vB : faceValid (nB 0) (nB 1) (nB 2) = true
  hS : nA 2 = nB 0
  hm : mS ≤ 1
  hmS : mS = 1 → nA 2 = 0

/-- The straightness condition of `check_straight_arcs` on this side. -/
def Straight : Prop := min (cc nA 0 + ea nA 1) (cc nB 0 + eb nB 2) = 0

variable {nA nB mS}

lemma Side.facts (h : Side nA nB mS) :
    nA 0 = cc nA 0 + cc nA 1 + ea nA 2 ∧ nA 1 = cc nA 1 + cc nA 2 + ea nA 0 ∧
    nA 2 = cc nA 2 + cc nA 0 + ea nA 1 ∧
    (0 < ea nA 0 → cc nA 0 = 0 ∧ ea nA 1 = 0 ∧ ea nA 2 = 0) ∧
    (0 < ea nA 1 → cc nA 1 = 0 ∧ ea nA 2 = 0 ∧ ea nA 0 = 0) ∧
    (0 < ea nA 2 → cc nA 2 = 0 ∧ ea nA 0 = 0 ∧ ea nA 1 = 0) ∧
    nB 0 = cc nB 0 + cc nB 1 + eb nB 2 ∧ nB 1 = cc nB 1 + cc nB 2 + eb nB 0 ∧
    nB 2 = cc nB 2 + cc nB 0 + eb nB 1 ∧
    (0 < eb nB 0 → cc nB 0 = 0 ∧ eb nB 1 = 0 ∧ eb nB 2 = 0) ∧
    (0 < eb nB 1 → cc nB 1 = 0 ∧ eb nB 2 = 0 ∧ eb nB 0 = 0) ∧
    (0 < eb nB 2 → cc nB 2 = 0 ∧ eb nB 0 = 0 ∧ eb nB 1 = 0) := by
  have hA := admissible_cc h.vA
  have hB := admissible_cc h.vB
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  have a0 := hA.n_eq 0; have a1 := hA.n_eq 1; have a2 := hA.n_eq 2
  have b0 := hB.n_eq 0; have b1 := hB.n_eq 1; have b2 := hB.n_eq 2
  have p0 := hA.2.1 0; have p1 := hA.2.1 1; have p2 := hA.2.1 2
  have q0 := hA.2.2 0; have q1 := hA.2.2 1; have q2 := hA.2.2 2
  have r0 := hB.2.1 0; have r1 := hB.2.1 1; have r2 := hB.2.1 2
  have s0 := hB.2.2 0; have s1 := hB.2.2 1; have s2 := hB.2.2 2
  simp only [l1, l2, l3, l4, l5, l6] at *
  simp only [ea, eb]
  exact ⟨a0, a1, a2, fun h => ⟨p0 h, q0 h⟩, fun h => ⟨p1 h, q1 h⟩, fun h => ⟨p2 h, q2 h⟩,
    b0, b1, b2, fun h => ⟨r0 h, s0 h⟩, fun h => ⟨r1 h, s1 h⟩, fun h => ⟨r2 h, s2 h⟩⟩

set_option maxHeartbeats 1000000 in
/-- **The merged face is valid** (one side of the merge): under the straightness condition, the
corner data `mergeC` is admissible for the side counts `mergeN`; in particular the merged edge
carries `n⁺_{X2 m} + n⁺_{m X0} + ends` crossings (the ends at `m` on this side become crossings).
[formula: ic-merge-edge-normal] -/
theorem Side.admissible (h : Side nA nB mS) (hs : Straight nA nB) :
    Admissible (mergeN nA nB mS) (mergeC nA nB) := by
  obtain ⟨fA0, fA1, fA2, eA0, eA1, eA2, fB0, fB1, fB2, eB0, eB1, eB2⟩ := h.facts
  have hS := h.hS
  have hm := h.hm
  have hmS := h.hmS
  unfold Straight at hs
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  refine ⟨fun k => ?_, fun k => ?_, fun k => ?_⟩ <;>
    rcases fin3_cases k with rfl | rfl | rfl <;>
    simp only [em, mergeN, mergeC, ends, l1, l2, l3, l4, l5, l6, Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] <;>
    generalize cc nA 0 = x0 at * <;> generalize cc nA 1 = x1 at * <;>
    generalize cc nA 2 = x2 at * <;> generalize cc nB 0 = y0 at * <;>
    generalize cc nB 1 = y1 at * <;> generalize cc nB 2 = y2 at * <;>
    generalize ea nA 0 = a0 at * <;> generalize ea nA 1 = a1 at * <;>
    generalize ea nA 2 = a2 at * <;> generalize eb nB 0 = b0 at * <;>
    generalize eb nB 1 = b1 at * <;> generalize eb nB 2 = b2 at * <;>
    rcases Glue.em_cases eA0 eA1 eA2 with ⟨z1, z2, z3⟩ | ⟨z1, z2, z3⟩ | ⟨z1, z2, z3⟩ |
      ⟨z1, z2, z3⟩ <;>
    rcases Glue.em_cases eB0 eB1 eB2 with ⟨w1, w2, w3⟩ | ⟨w1, w2, w3⟩ | ⟨w1, w2, w3⟩ |
      ⟨w1, w2, w3⟩ <;>
    subst z1 z2 z3 w1 w2 w3 <;> omega

/-- The Rust `merge_edge_normal` (without the `i64` overflow): `−1` if both halves run along an
input edge (and `m` lies on it), `n₁ + n₂ + [m on an input edge]` if both cross, else refused. -/
def mergeNormal (n1 n2 : ℤ) (onEdge : Bool) : Option ℤ :=
  if n1 = -1 ∧ n2 = -1 then (if onEdge then some (-1) else none)
  else if 0 ≤ n1 ∧ 0 ≤ n2 then some (n1 + n2 + if onEdge then 1 else 0) else none

/-- **`ic-merge-edge-normal`: the normal coordinate of two merged halves.** The halves `m X0` (`n₁`)
and `X2 m` (`n₂`) of an inserted vertex `m`, with the faces `A, B` on one side and `A', B'` on the
other (the same numbering, the halves reversed: `n⁺₁` is side `0` of `A` and side `2` of `B'`,
`n⁺₂` side `2` of `B` and side `0` of `A'`), both sides straight (`check_straight_arcs`), and the
realization at `m` (SPEC §2.4 item 4, the second check of the merge): the same number `τ ≤ 1` of
curve ends at `m` on both sides, `τ = 1` exactly when `m` lies on an input edge the halves do not
run along, and both halves along an input edge only when `m` lies on it (then no other end). The
halves agree in running along an input edge or not (else `HalvesDisagree`). Then the Rust
`merge_edge_normal n₁ n₂ on_edge` returns `n` with `n⁺` = the crossings of the merged edge on
either side (`n⁺₁ + n⁺₂ + τ`, the merged faces being valid, `Side.admissible`), and `n = −1`
exactly when both halves run along an input edge (the merged edge then runs along it).
[formula: ic-merge-edge-normal] -/
theorem merge_normal {nA nB nA' nB' : Fin 3 → ℕ} {mS mS' : ℕ}
    (hU : Side nA nB mS) (hsU : Straight nA nB) (hL : Side nA' nB' mS') (hsL : Straight nA' nB')
    {n1 n2 : ℤ} (g1 : -1 ≤ n1) (g2 : -1 ≤ n2)
    (h1 : nPlus n1 = nA 0) (h1' : nPlus n1 = nB' 2) (h2 : nPlus n2 = nB 2) (h2' : nPlus n2 = nA' 0)
    {τ : ℕ} (hτ : ends nA nB mS = τ) (hτ' : ends nA' nB' mS' = τ) {onEdge : Bool}
    (agree : n1 = -1 ↔ n2 = -1)
    (along : n1 = -1 → onEdge = true ∧ τ = 0) (cross : 0 ≤ n1 → (τ = 1 ↔ onEdge = true))
    (hτ1 : τ ≤ 1) :
    Admissible (mergeN nA nB mS) (mergeC nA nB) ∧ Admissible (mergeN nA' nB' mS') (mergeC nA' nB') ∧
      ∃ n, mergeNormal n1 n2 onEdge = some n ∧ (nPlus n : ℤ) = mergeN nA nB mS 0 ∧
        (nPlus n : ℤ) = mergeN nA' nB' mS' 0 ∧ (n = -1 ↔ n1 = -1) := by
  refine ⟨hU.admissible hsU, hL.admissible hsL, ?_⟩
  simp only [mergeN, Matrix.cons_val_zero]
  rcases eq_or_lt_of_le g1 with e1 | e1
  · -- both halves along an input edge
    have e2 := agree.1 e1.symm
    obtain ⟨ho, ht⟩ := along e1.symm
    refine ⟨-1, by simp [mergeNormal, ← e1, e2, ho], ?_, ?_, by simp [← e1]⟩
    · rw [← e1] at h1; rw [e2] at h2
      simp only [nPlus] at h1 h2 ⊢; push_cast; rw [hτ]; omega
    · rw [← e1] at h1'; rw [e2] at h2'
      simp only [nPlus] at h1' h2' ⊢; push_cast; rw [hτ']; omega
  · have n1ne : n1 ≠ -1 := by omega
    have e2 : 0 ≤ n2 := by
      rcases eq_or_lt_of_le g2 with e | e
      · exact absurd (agree.2 e.symm) n1ne
      · omega
    have hc := cross (by omega)
    have t : ((if onEdge then 1 else 0) : ℤ) = τ := by
      cases onEdge <;> simp at hc ⊢ <;> omega
    have a1 : n1 = (nA 0 : ℤ) := by simp only [nPlus] at h1; omega
    have a2 : n2 = (nB 2 : ℤ) := by simp only [nPlus] at h2; omega
    have a1' : n1 = (nB' 2 : ℤ) := by simp only [nPlus] at h1'; omega
    have a2' : n2 = (nA' 0 : ℤ) := by simp only [nPlus] at h2'; omega
    refine ⟨n1 + n2 + if onEdge then 1 else 0, ?_, ?_, ?_, ?_⟩
    · simp only [mergeNormal]; rw [if_neg (by omega), if_pos ⟨by omega, e2⟩]
    · rw [t]; simp only [nPlus]; push_cast; rw [hτ]; omega
    · rw [t]; simp only [nPlus]; push_cast; rw [hτ']; omega
    · constructor
      · intro h; split_ifs at h <;> omega
      · intro h; exact absurd h n1ne

end Merge

end DDG.IntegerCoords
