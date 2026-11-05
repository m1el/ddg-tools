import DDG.NormalCoords.SplitEdge
import DDG.NormalCoords.Complex.Basic

/-!
# Integer coordinates: splitting an edge exactly at a crossing (`ic-split-at-crossing`, tier B)

`crates/ddg-integer-coords/SPEC.md` §4.5 ("At a crossing", I28); the Rust
`formulas::split_side_at_crossing` (`crates/ddg-integer-coords/src/formulas.rs`), used by
`split_integers` (`src/split.rs`).

A face `abc` (corners `0, 1, 2`, side `0` = `a→b`, the face model of `FaceModel.lean`) is split on
side `0` **at** its crossing `j` (`j < n₀`, counted from `a`): the new vertex `m` is the boundary
point `side 0 j`, of rank `j + 1`, and lies on the arc `A` of the face model through that crossing
(`canon n c 0 j` is its other end). The new straight edge `mc` crosses exactly the arcs whose ends
interleave with `{m, c}` (`ic-chord-interleave`); `A` shares the end `m` with it and the arcs from
`c` share the end `c`, so neither crosses it. `chordCount n c r k` counts the arcs of the face model
whose ends interleave with the chord from rank `r` to corner `k`.

* `canon_side0`: the arc through `m` cuts `a` (`j < c_a`), comes from `c`
  (`c_a ≤ j < c_a + e_c`) or cuts `b` (`j ≥ c_a + e_c`).
* `chordCount_eq_sepCount`: the arcs crossing `mc` are those separating the gap next to `m` on the
  side where `A` does not separate it from `c` (gap `j + 1` if `A` cuts `a`, else gap `j`) from
  `c`; so (`chordCount_closed`, from WP9's `splitSide_eq`) their number is the Rust
  `split_side(n, 0, j + 1)` if `A` cuts `a`, `0` if `A` comes from `c`, and
  `split_side(n, 0, j)` if `A` cuts `b`.
* `splitSideAtCrossing_eq` (**`ic-split-at-crossing`**, side `s` of a face with any normal
  coordinates `n`): the Lean transcription `splitSideAtCrossingR` of the Rust function returns
  `(−1 if A comes from c, else the number of arcs crossing mc; the arcs from c ending on the split
  side before m)`; in the `−1` case no arc crosses `mc` (it runs along the input edge of `A`).
  `splitSideAtCrossingR_eq_none`: `none` exactly for an invalid face or `j` out of range.

The two new faces and the merge back are in `Complex/SplitAtCrossing.lean`.
-/

namespace DDG.NormalCoords

open DDG.NormalCoords.Formulas (cornerNumerator cornerCoord emanating nPlus faceValid)

/-- Decidability of interleaving (named, to stay out of the way of other instances). -/
instance decInterleaveSAC (a b x y : ℕ) : Decidable (Interleave a b x y) := by
  unfold Interleave; infer_instance

/-! ## The arcs crossing a chord from a boundary point to a corner -/

/-- The number of arcs of the face model with side counts `n` and corner counts `c` (enumerated by
the families `cornerArc`, `fanArc`, each arc once: `FaceModel.Valid.arcs`) that **cross the chord**
from the boundary point of rank `r` to corner `k`: their ends interleave with `{r, rank of k}`. An
arc with an end at `r` or at corner `k` does not cross it. [formula: ic-split-at-crossing] -/
def chordCount (n c : Fin 3 → ℕ) (r : ℕ) (k : Fin 3) : ℕ :=
  ∑ i : Fin 3, ((Finset.range (c i)).filter
      (fun t => Interleave r (cornerRank n k) (cornerArc n i t).1 (cornerArc n i t).2)).card +
    ∑ κ : Fin 3, ((Finset.range (em n c κ)).filter
      (fun s => Interleave r (cornerRank n k) (fanArc n c κ s).1 (fanArc n c κ s).2)).card

/-- **Crossing a chord from a boundary point = separating the gap next to it.** For an arc
`(x, y)` without an end at `r`, crossing the chord `{r, q}` is separating either gap next to `r`
(`g = r` or `g = r − 1`) from `q`. -/
lemma interleave_iff_separates {r q x y g : ℕ} (hx : x ≠ r) (hy : y ≠ r)
    (hg : g = r ∨ g + 1 = r) : Interleave r q x y ↔ Separates (x, y) g q := by
  simp only [Interleave, Separates, GapIn, Btw]
  omega

/-- The arc through crossing `j` of side `0`: it cuts corner `a` (`j < c₀`, to side `2`), comes
from corner `c` (`c₀ ≤ j < c₀ + e₂`) or cuts corner `b` (to side `1`). -/
lemma canon_side0 {n c : Fin 3 → ℕ} (hc : Admissible n c) (j : ℕ) :
    canon n c 0 j = if j < c 0 then .side 2 (n 2 - 1 - j)
      else if j < c 0 + em n c 2 then .corner 2 else .side 1 (n 0 - 1 - j) := by
  have e := hc.em_add
  obtain ⟨l1, l2, -⟩ := fin3_lits
  unfold canon
  simp only [l1, l2]
  by_cases h1 : j < c 0
  · simp only [h1, ↓reduceIte]
  · simp only [h1, ↓reduceIte]
    by_cases h2 : j < c 0 + em n c 2
    · rw [ite_eq_left (by omega), ite_eq_left h2]
    · rw [ite_eq_right (by omega), ite_eq_right h2]

/-- The arc through `m` comes from `c` exactly in the middle block. -/
lemma canon_side0_eq_corner {n c : Fin 3 → ℕ} (hc : Admissible n c) (j : ℕ) :
    canon n c 0 j = .corner 2 ↔ c 0 ≤ j ∧ j < c 0 + em n c 2 := by
  rw [canon_side0 hc]
  split_ifs with h1 h2
  · exact ⟨fun h => (by cases h), fun h => (by omega)⟩
  · exact ⟨fun _ => ⟨by omega, h2⟩, fun _ => rfl⟩
  · exact ⟨fun h => (by cases h), fun h => (by omega)⟩

/-- The gap next to `m` on the side where the arc `A` through `m` does not separate it from `c`. -/
def crossGap (c : Fin 3 → ℕ) (j : ℕ) : ℕ := if j < c 0 then j + 1 else j

lemma crossGap_cases (c : Fin 3 → ℕ) (j : ℕ) :
    (j < c 0 ∧ crossGap c j = j + 1) ∨ (c 0 ≤ j ∧ crossGap c j = j) := by
  unfold crossGap; split_ifs <;> omega

/-- An arc around a corner with an end at `m` (that is, `A`) does not separate `crossGap` from
`c`. -/
lemma not_sep_cornerArc {n c : Fin 3 → ℕ} (hc : Admissible n c) {j : ℕ} (hj : j < n 0)
    (i : Fin 3) {t : ℕ} (ht : t < c i)
    (he : (cornerArc n i t).1 = j + 1 ∨ (cornerArc n i t).2 = j + 1) :
    ¬ Separates (cornerArc n i t) (crossGap c j) (cornerRank n 2) := by
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := hc.em_add
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  rcases crossGap_cases c j with ⟨hg1, hg⟩ | ⟨hg1, hg⟩ <;> rw [hg] <;>
  rcases fin3_cases i with rfl | rfl | rfl <;>
    simp only [Separates, cornerArc, Pt.rank, GapIn, Btw, l2, l4, l6,
      cornerRank_zero, cornerRank_one, cornerRank_two, Fin.isValue] at ht he ⊢ <;> omega

/-- An arc from a corner with an end at `m` (that is, `A`, from `c`) does not separate `crossGap`
from `c`. -/
lemma not_sep_fanArc {n c : Fin 3 → ℕ} (hc : Admissible n c) {j : ℕ} (hj : j < n 0)
    (κ : Fin 3) {s : ℕ} (hs : s < em n c κ)
    (he : (fanArc n c κ s).1 = j + 1 ∨ (fanArc n c κ s).2 = j + 1) :
    ¬ Separates (fanArc n c κ s) (crossGap c j) (cornerRank n 2) := by
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := hc.em_add
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  rcases crossGap_cases c j with ⟨hg1, hg⟩ | ⟨hg1, hg⟩ <;> rw [hg] <;>
  rcases fin3_cases κ with rfl | rfl | rfl <;>
    simp only [Separates, fanArc, Pt.rank, GapIn, Btw, l1, l3, l5,
      cornerRank_zero, cornerRank_one, cornerRank_two, Fin.isValue] at hs he ⊢ <;> omega

/-- Per arc: crossing `mc` is separating `crossGap` from `c`. -/
lemma interleave_iff_sep_gap {x y j g q : ℕ} (hg : g = j + 1 ∨ g = j)
    (hA : (x = j + 1 ∨ y = j + 1) → ¬ Separates (x, y) g q) :
    Interleave (j + 1) q x y ↔ Separates (x, y) g q := by
  by_cases h : x = j + 1 ∨ y = j + 1
  · have := hA h
    simp only [this, iff_false]
    unfold Interleave; omega
  · exact interleave_iff_separates (by omega) (by omega) (by omega)

/-- **The arcs crossing the new edge `mc` are the arcs separating `crossGap` from `c`.**
[formula: ic-split-at-crossing] -/
theorem chordCount_eq_sepCount {n c : Fin 3 → ℕ} (hc : Admissible n c) {j : ℕ} (hj : j < n 0) :
    chordCount n c (j + 1) 2 = sepCount n c (crossGap c j) 2 := by
  have hg : crossGap c j = j + 1 ∨ crossGap c j = j := by
    rcases crossGap_cases c j with ⟨-, h⟩ | ⟨-, h⟩ <;> omega
  unfold chordCount sepCount
  congr 1
  · refine Finset.sum_congr rfl (fun i _ => ?_)
    congr 1
    refine Finset.filter_congr (fun t ht => ?_)
    exact interleave_iff_sep_gap hg
      (not_sep_cornerArc hc hj i (Finset.mem_range.mp ht))
  · refine Finset.sum_congr rfl (fun κ _ => ?_)
    congr 1
    refine Finset.filter_congr (fun s hs => ?_)
    exact interleave_iff_sep_gap hg
      (not_sep_fanArc hc hj κ (Finset.mem_range.mp hs))

/-- **The number of arcs crossing `mc`, in closed form**: with `c`, `e` the corner and emanating
counts, `c₂ + e₀ + e₁ + (c₀ − j − 1)` if `A` cuts `a`, `0` if `A` comes from `c`, and
`c₂ + e₀ + e₁ + (j − c₀ − e₂)` if `A` cuts `b`; and these are WP9's `split_side` values (the right
side of `splitSide_eq`) at `before = j + 1`, (`0` in the middle block), `before = j`.
[formula: ic-split-at-crossing] -/
theorem chordCount_closed {n c : Fin 3 → ℕ} (hc : Admissible n c) {j : ℕ} (hj : j < n 0) :
    chordCount n c (j + 1) 2 =
      if j < c 0 then c 2 + em n c 0 + em n c 1 + (c 0 - j - 1)
      else if j < c 0 + em n c 2 then 0
      else c 2 + em n c 0 + em n c 1 + (j - c 0 - em n c 2) := by
  have hb : crossGap c j ≤ n 0 := by unfold crossGap; split_ifs <;> omega
  have h := splitSide_eq hc hb
  rw [← chordCount_eq_sepCount hc hj] at h
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  have k0 := hc.coords 0
  have k1 := hc.coords 1
  have k2 := hc.coords 2
  simp only [l1, l2, l3, l4, l5, l6] at k0 k1 k2
  rw [k0.1, k1.1, k2.1, k0.2, k1.2] at h
  obtain ⟨a0, -, -, -, -, b2, -, -, x2⟩ := hc.lits
  obtain ⟨-, -, m2⟩ := hc.em_add
  rcases crossGap_cases c j with ⟨hg1, hg⟩ | ⟨hg1, hg⟩ <;> rw [hg] at h
  · rw [ite_eq_left hg1]; omega
  · rw [ite_eq_right (by omega)]
    by_cases h2 : j < c 0 + em n c 2
    · rw [ite_eq_left h2]
      have := b2 (by omega)
      have := x2 (by omega)
      omega
    · rw [ite_eq_right h2]; omega

/-- In the middle block (`A` comes from `c`) no arc crosses `mc`: the new edge runs along the input
edge of `A`. [formula: ic-split-at-crossing] -/
theorem chordCount_fan {n c : Fin 3 → ℕ} (hc : Admissible n c) {j : ℕ} (hj : j < n 0)
    (hA : canon n c 0 j = .corner 2) : chordCount n c (j + 1) 2 = 0 := by
  rw [chordCount_closed hc hj]
  obtain ⟨h1, h2⟩ := (canon_side0_eq_corner hc j).1 hA
  rw [ite_eq_right (by omega), ite_eq_left h2]

/-! ## The arcs from `c` ending on the split side before `m` -/

/-- The arcs from `c` that end on `ab` before `m` (`A` excluded: it ends at `m`): the Rust
`emanating_before`, `clamp(before − c_a, 0, e_c)` with `before` = `j + 1` or `j`.
[formula: ic-split-at-crossing] -/
theorem emanatingBefore_at {n c : Fin 3 → ℕ} (j : ℕ) :
    (((Finset.range (em n c 2)).filter (fun s => c 0 + s < j)).card : ℤ) =
      max 0 (min ((crossGap c j : ℤ) - c 0) (em n c 2)) := by
  rw [emanatingBefore_eq j]
  unfold crossGap
  split_ifs <;> omega

/-! ## The Rust function -/

/-- The Rust `FaceCoordinates::of(n).c[k]` (Eq. 1). -/
def fcC (n : Fin 3 → ℤ) (k : Fin 3) : ℤ := cornerCoord (n (k + 1)) (n k) (n (k + 2))

/-- The Rust `FaceCoordinates::of(n).e[k]` (`face_emanating`, Eq. 2). -/
def fcE (n : Fin 3 → ℤ) (k : Fin 3) : ℤ := emanating (n (k + 1)) (n k) (n (k + 2))

/-- The Rust `formulas::split_side(n, s, before)` without the `i64` overflow: `(to_opposite,
emanating_before)`, `none` for an invalid face (`FaceCoordinates::of` fails exactly when a
numerator of Eq. 1 is odd, `faceValid`; it is never negative) or `before ∉ [0, n⁺_s]`. Rust's
`clamp(0, e_c)` is `max 0 (min · e_c)` since `e_c ≥ 0`. [formula: ic-split-crossed-edge] -/
def splitSideR (n : Fin 3 → ℤ) (s : Fin 3) (before : ℤ) : Option (ℤ × ℤ) :=
  if faceValid (n 0) (n 1) (n 2) = true then
    if 0 ≤ before ∧ before ≤ nPlus (n s) then
      some (fcC n (s + 2) + fcE n s + fcE n (s + 1) + max (fcC n s - before) 0 +
          max (before - (nPlus (n s) - fcC n (s + 1))) 0,
        max 0 (min (before - fcC n s) (fcE n (s + 2))))
    else none
  else none

/-- The Rust `formulas::split_side_at_crossing(n, s, j)` without the `i64` overflow:
`(to_opposite, emanating_before)`. [formula: ic-split-at-crossing] -/
def splitSideAtCrossingR (n : Fin 3 → ℤ) (s : Fin 3) (j : ℤ) : Option (ℤ × ℤ) :=
  if faceValid (n 0) (n 1) (n 2) = true then
    if 0 ≤ j ∧ j < nPlus (n s) then
      match splitSideR n s (if j < fcC n s then j + 1 else j) with
      | some x => some (if fcC n s ≤ j ∧ j < fcC n s + fcE n (s + 2) then -1 else x.1, x.2)
      | none => none
    else none
  else none

lemma faceValid_cycle (a b c : ℤ) : faceValid b c a = faceValid a b c := by
  simp only [faceValid]
  cases cornerNumerator b a c % 2 == 0 <;> cases cornerNumerator c b a % 2 == 0 <;>
    cases cornerNumerator a c b % 2 == 0 <;> rfl

lemma fin3_rot_lits : (∀ s : Fin 3, 0 + s = s) ∧ (∀ k s : Fin 3, k + 1 + s = k + s + 1) ∧
    (∀ k s : Fin 3, k + 2 + s = k + s + 2) ∧ (∀ s : Fin 3, 1 + s = s + 1) ∧
    (∀ s : Fin 3, 2 + s = s + 2) := by decide

lemma faceValid_rotSide (n : Fin 3 → ℤ) (s : Fin 3) :
    faceValid (n (0 + s)) (n (1 + s)) (n (2 + s)) = faceValid (n 0) (n 1) (n 2) := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  have z : (0 : Fin 3) + 0 = 0 ∧ (1 : Fin 3) + 0 = 1 ∧ (2 : Fin 3) + 0 = 2 := by decide
  rcases fin3_cases s with rfl | rfl | rfl
  · rw [z.1, z.2.1, z.2.2]
  · rw [l1, l3, l5, faceValid_cycle]
  · rw [l2, l4, l6, faceValid_cycle (n 1), faceValid_cycle]

/-- The Rust function on side `s` is the same function on side `0` of the face numbered from
corner `s`. -/
theorem splitSideAtCrossingR_rot (n : Fin 3 → ℤ) (s : Fin 3) (j : ℤ) :
    splitSideAtCrossingR n s j = splitSideAtCrossingR (fun k => n (k + s)) 0 j := by
  obtain ⟨r0, r1, r2, r3, r4⟩ := fin3_rot_lits
  have hc : ∀ k, fcC (fun k => n (k + s)) k = fcC n (k + s) := fun k => by
    simp only [fcC, r1, r2]
  have he : ∀ k, fcE (fun k => n (k + s)) k = fcE n (k + s) := fun k => by
    simp only [fcE, r1, r2]
  simp only [splitSideAtCrossingR, splitSideR, hc, he]
  rw [faceValid_rotSide]
  simp only [r0, r3, r4]

lemma cornerCoord_nPlus (a b c : ℤ) :
    cornerCoord a b c = cornerCoord (nPlus a) (nPlus b) (nPlus c) := by
  have h : ∀ x : ℤ, nPlus (nPlus x) = nPlus x := fun x => by simp only [nPlus]; omega
  simp only [cornerCoord, cornerNumerator, emanating, h]

lemma emanating_nPlus (a b c : ℤ) :
    emanating a b c = emanating (nPlus a) (nPlus b) (nPlus c) := by
  have h : ∀ x : ℤ, nPlus (nPlus x) = nPlus x := fun x => by simp only [nPlus]; omega
  simp only [emanating, h]

lemma nPlus_eq_toNat (x : ℤ) : nPlus x = (x.toNat : ℤ) := by
  simp only [nPlus]; omega

/-- **`ic-split-at-crossing`: the Rust `split_side_at_crossing` on side `0`.** For a valid face
with normal coordinates `n` (face model: side counts `N = n⁺`, corner counts `c` = Eq. 1) and a
crossing `j < n⁺₀` of side `0`, the Rust function returns `(n_mc, emanating_before)` with
`n_mc = −1` if the arc `A` through crossing `j` comes from corner `2` (then no arc crosses `mc`,
`chordCount_fan`: `mc` runs along the input edge of `A`), else `n_mc` = the number of arcs of the
face model crossing the chord from `m` (rank `j + 1`) to corner `2` (`A` not counted: it ends at
`m`); `emanating_before` is the number of arcs from corner `2` ending on side `0` before `m` (`A`
not counted). [formula: ic-split-at-crossing] -/
theorem splitSideAtCrossingR_side0 {n : Fin 3 → ℤ} (hv : faceValid (n 0) (n 1) (n 2) = true)
    {j : ℕ} (hj : (j : ℤ) < nPlus (n 0)) :
    splitSideAtCrossingR n 0 j = some
      (if canon (fun k => (n k).toNat) (cc fun k => (n k).toNat) 0 j = .corner 2 then -1
        else (chordCount (fun k => (n k).toNat) (cc fun k => (n k).toNat) (j + 1) 2 : ℤ),
       (((Finset.range (em (fun k => (n k).toNat) (cc fun k => (n k).toNat) 2)).filter
          (fun s => cc (fun k => (n k).toNat) 0 + s < j)).card : ℤ)) := by
  set N : Fin 3 → ℕ := fun k => (n k).toNat with hN
  have hNk : ∀ k, nPlus (n k) = (N k : ℤ) := fun k => nPlus_eq_toNat _
  have hvN : faceValid (N 0) (N 1) (N 2) = true := by
    rw [← hNk, ← hNk, ← hNk, ← faceValid_nPlus]; exact hv
  have hc := admissible_cc hvN
  set c := cc N with hcdef
  have hC : ∀ k, fcC n k = c k := fun k => by
    rw [fcC, cornerCoord_nPlus, hNk, hNk, hNk, (hc.coords k).1]
  have hE : ∀ k, fcE n k = em N c k := fun k => by
    rw [fcE, emanating_nPlus, hNk, hNk, hNk, (hc.coords k).2]
  have hj' : j < N 0 := by rw [hNk] at hj; exact_mod_cast hj
  have hcl := chordCount_closed hc hj'
  have heb := emanatingBefore_at (n := N) (c := c) j
  have hfan := canon_side0_eq_corner hc j
  obtain ⟨a0, -⟩ := hc.lits
  obtain ⟨-, -, m2⟩ := hc.em_add
  obtain ⟨l1, l2, -⟩ := fin3_lits
  have k0 := hC 0; have k1 := hC 1; have k2 := hC 2
  have g0 := hE 0; have g1 := hE 1; have g2 := hE 2
  simp only [splitSideAtCrossingR, splitSideR, hv, ↓reduceIte, l1, l2, hNk 0, k0, k1, k2,
    g0, g1, g2]
  rw [ite_eq_left ⟨by omega, by omega⟩]
  by_cases h1 : j < c 0
  · rw [crossGap, ite_eq_left h1] at heb
    rw [ite_eq_left (show (j : ℤ) < c 0 by omega), ite_eq_left ⟨by omega, by omega⟩]
    simp only [Option.some.injEq, Prod.mk.injEq]
    refine ⟨?_, ?_⟩
    · rw [ite_eq_right (by omega), ite_eq_right (fun h => by have := (hfan.1 h).1; omega), hcl,
        ite_eq_left h1]
      push_cast
      omega
    · rw [heb]; push_cast; omega
  · rw [crossGap, ite_eq_right h1] at heb
    rw [ite_eq_right (show ¬ (j : ℤ) < c 0 by omega), ite_eq_left ⟨by omega, by omega⟩]
    simp only [Option.some.injEq, Prod.mk.injEq]
    refine ⟨?_, by rw [heb]⟩
    by_cases h2 : j < c 0 + em N c 2
    · rw [ite_eq_left ⟨by omega, by omega⟩, ite_eq_left (hfan.2 ⟨by omega, h2⟩)]
    · rw [ite_eq_right (by omega), ite_eq_right (fun h => h2 (hfan.1 h).2), hcl, ite_eq_right h1,
        ite_eq_right h2]
      push_cast
      omega

/-- **`ic-split-at-crossing`: the Rust `split_side_at_crossing(n, s, j)`.** For a valid face with
normal coordinates `n`, a side `s` and a crossing `0 ≤ j < n⁺_s` of it, the Rust function returns,
in the face model of the face numbered from corner `s` (side counts `N k = n⁺_{k+s}`, corner
counts `c` = Eq. 1; `a, b, c` = corners `s, s + 1, s + 2`): `n_mc = −1` if the arc through
crossing `j` comes from `c` (`canon N c 0 j = corner 2`; then no arc crosses `mc`,
`chordCount_fan`), else the number of arcs crossing the chord from `m` to `c`
(`chordCount N c (j + 1) 2`, which is `split_side` at `before = j + 1` if the arc cuts `a`, at
`before = j` if it cuts `b`: `chordCount_eq_sepCount`, `chordCount_closed`); and
`emanating_before` = the arcs from `c` ending on `ab` before `m`. [formula: ic-split-at-crossing] -/
theorem splitSideAtCrossingR_eq {n : Fin 3 → ℤ} (s : Fin 3)
    (hv : faceValid (n 0) (n 1) (n 2) = true) {j : ℕ} (hj : (j : ℤ) < nPlus (n s)) :
    splitSideAtCrossingR n s j = some
      (if canon (fun k => (n (k + s)).toNat) (cc fun k => (n (k + s)).toNat) 0 j = .corner 2 then -1
        else (chordCount (fun k => (n (k + s)).toNat) (cc fun k => (n (k + s)).toNat) (j + 1) 2 : ℤ),
       (((Finset.range (em (fun k => (n (k + s)).toNat) (cc fun k => (n (k + s)).toNat) 2)).filter
          (fun t => cc (fun k => (n (k + s)).toNat) 0 + t < j)).card : ℤ)) := by
  obtain ⟨r0, -⟩ := fin3_rot_lits
  rw [splitSideAtCrossingR_rot]
  refine splitSideAtCrossingR_side0 (n := fun k => n (k + s)) ?_ ?_
  · rw [faceValid_rotSide]; exact hv
  · simp only [r0]; exact hj

/-- **`none` exactly for an invalid face or a crossing out of range** (`0 ≤ j < n⁺_s`).
[formula: ic-split-at-crossing] -/
theorem splitSideAtCrossingR_eq_none (n : Fin 3 → ℤ) (s : Fin 3) (j : ℤ) :
    splitSideAtCrossingR n s j = none ↔
      faceValid (n 0) (n 1) (n 2) = false ∨ j < 0 ∨ nPlus (n s) ≤ j := by
  by_cases hv : faceValid (n 0) (n 1) (n 2) = true
  · by_cases hj : 0 ≤ j ∧ j < nPlus (n s)
    · obtain ⟨k, rfl⟩ := Int.eq_ofNat_of_zero_le hj.1
      rw [splitSideAtCrossingR_eq s hv hj.2]
      simp only [reduceCtorEq, false_iff, not_or, Bool.not_eq_false]
      exact ⟨hv, by omega, by omega⟩
    · simp only [splitSideAtCrossingR, hv, ↓reduceIte, hj, true_iff]
      omega
  · simp only [Bool.not_eq_true] at hv
    simp only [splitSideAtCrossingR, hv, Bool.false_eq_true, ↓reduceIte, true_or]

end DDG.NormalCoords
