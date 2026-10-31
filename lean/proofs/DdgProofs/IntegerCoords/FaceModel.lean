import DdgProofs.IntegerCoords.Formulas
import Mathlib.Data.Fintype.Fin
import Mathlib.Data.Finset.Card
import Mathlib.Tactic.FinCases

/-!
# Integer coordinates: the face model (tier B)

`crates/ddg-integer-coords/SPEC.md` §9.1–9.2; the paper Sec. 2.3, Fig. 6, 7, 11, Supp. Alg. 1, 9.
The design and how tier C builds on it: `DESIGN.md` in this directory.

A **face model** is one triangle of `T¹` with the pieces of the curves of `T⁰` inside it. Corners
`0, 1, 2`; side `k` runs from corner `k` to corner `k + 1` (the halfedge `3f + k`), and carries
`n k` crossings numbered `0, 1, …` from corner `k` (the paper's combinatorial crossings `(ij→, p)`).
The boundary of the triangle is the cycle

    corner 0, side 0 (0 … n₀−1), corner 1, side 1 (0 … n₁−1), corner 2, side 2 (0 … n₂−1)

and `Pt.rank` numbers it `0, 1, …, n₀ + n₁ + n₂ + 2` in this order. Every crossing is the end of
exactly one arc; `partner k p` is the other end of the arc through crossing `p` of side `k` (a
crossing or a corner). An arc between two corners would run along a side (an edge with `n = −1`)
and is not part of a face model.

Validity (`FaceModel.Valid`) is purely combinatorial:

* *normal*: an arc from side `k` ends on one of the two other sides or at the opposite corner
  `k + 2` (a straight segment cannot meet the same side twice, nor join a corner to a side through
  it);
* *symmetric*: the partner of a partner is the crossing itself;
* *disjoint*: the endpoints of two arcs do not interleave on the boundary cycle (`Interleave`; arcs
  that share a corner do not interleave). For straight arcs this is the same as not crossing
  (`ic-chord-interleave`, `DdgProofs/IntegerCoords/Chords.lean`).

Results.

* **The structure theorem** (`FaceModel.Valid.eq_canon`): every valid model is the canonical one
  (`canon n c`, the oracle's `Ddg.IntegerCoords.partner`) for corner counts `c` that are admissible
  (`Admissible`): on side `k` the first `c k` crossings cut corner `k` and reach side `k + 2` in
  reverse order, the last `c (k+1)` cut corner `k + 1` and reach side `k + 1` in reverse order, the
  ones between come from corner `k + 2`; at most one corner has arcs to its opposite side, and then
  no arc cuts that corner. Conversely every admissible `canon n c` is valid (`canon_valid`).
* **`ic-corner-coordinates`, `ic-emanating`**: in a valid model the numbers of arcs that cut
  corner `k` and that leave corner `k` (`cutCount`, `emCount`) are exactly Eq. 1 and Eq. 2 of the
  side counts (`Valid.cornerNumerator_eq`, `Valid.cornerCoord_eq`, `Valid.emanating_eq`).
* **`ic-face-decoding`** (Fig. 11): the side counts determine the model (`Valid.partner_eq`), and a
  model with given side counts exists if and only if the face is valid in the sense of Eq. 1, i.e.
  the oracle's `faceValid` (`exists_valid_iff`); the oracle's `partner` is the model
  (`oracle_partner_eq`).
* **`ic-cevians-cross`** (combinatorial half): arcs from two different corners to their opposite
  sides interleave, so at most one corner emanates (`Valid.emanating_unique`), and an emanating
  corner is cut by no arc (`Valid.cut_eq_zero_of_emanating`).
* **`ic-trace-step`**: Alg. 1's three cases follow the unique arc through the crossing
  (`Valid.traceStep_correct`).
-/

namespace DdgProofs

namespace IntegerCoords

open Ddg.IntegerCoords (cornerNumerator cornerCoord emanating nPlus)

/-! ## Boundary points and interleaving -/

/-- A point of the boundary of a face: a corner, or crossing `p` of side `k`. -/
inductive Pt where
  | corner (k : Fin 3)
  | side (k : Fin 3) (p : ℕ)
  deriving DecidableEq

/-- Position of corner `k` on the boundary cycle. -/
def cornerRank (n : Fin 3 → ℕ) (k : Fin 3) : ℕ :=
  if k = 0 then 0 else if k = 1 then n 0 + 1 else n 0 + n 1 + 2

/-- Position of a boundary point on the boundary cycle. -/
def Pt.rank (n : Fin 3 → ℕ) : Pt → ℕ
  | .corner k => cornerRank n k
  | .side k p => cornerRank n k + 1 + p

/-- `x` lies strictly between `a` and `b` (in either order). -/
def Btw (a b x : ℕ) : Prop := (a < x ∧ x < b) ∨ (b < x ∧ x < a)

/-- The chords `{a, b}` and `{x, y}` of a cycle (positions on it) **interleave**: none of `x, y` is
an end of the first, and exactly one of them lies between `a` and `b`. Interleaving does not depend
on where the cycle is cut open; chords with a common end never interleave. -/
def Interleave (a b x y : ℕ) : Prop :=
  x ≠ a ∧ x ≠ b ∧ y ≠ a ∧ y ≠ b ∧ ((Btw a b x ∧ ¬ Btw a b y) ∨ (¬ Btw a b x ∧ Btw a b y))

/-! ## Face models -/

/-- A face model: the side counts and the other end of the arc through every crossing (the values
at `p ≥ n k` are irrelevant). -/
structure FaceModel where
  /-- crossings per side (the `n⁺` of the three edges) -/
  n : Fin 3 → ℕ
  /-- the other end of the arc through crossing `p` of side `k` -/
  partner : Fin 3 → ℕ → Pt

namespace FaceModel

/-- Validity: normal, symmetric, pairwise disjoint arcs. [formula: ic-face-decoding] -/
structure Valid (M : FaceModel) : Prop where
  normal : ∀ k p, p < M.n k →
    (∃ q, q < M.n (k + 2) ∧ M.partner k p = .side (k + 2) q) ∨
      (∃ q, q < M.n (k + 1) ∧ M.partner k p = .side (k + 1) q) ∨ M.partner k p = .corner (k + 2)
  symm : ∀ k p k' q, p < M.n k → M.partner k p = .side k' q → M.partner k' q = .side k p
  disjoint : ∀ k p k' p', p < M.n k → p' < M.n k' →
    ¬ Interleave ((Pt.side k p).rank M.n) ((M.partner k p).rank M.n)
      ((Pt.side k' p').rank M.n) ((M.partner k' p').rank M.n)

end FaceModel

/-! ## The canonical model -/

/-- Arcs from corner `k` to the opposite side `k + 1`, for corner counts `c`. -/
def em (n c : Fin 3 → ℕ) (k : Fin 3) : ℕ := n (k + 1) - c (k + 1) - c (k + 2)

/-- Corner counts `c` fit the side counts `n`: the two corner blocks of every side fit in it, at
most one corner emanates, and no arc cuts an emanating corner. -/
def Admissible (n c : Fin 3 → ℕ) : Prop :=
  (∀ k, c k + c (k + 1) ≤ n k) ∧ (∀ k, 0 < em n c k → c k = 0) ∧
    (∀ k, 0 < em n c k → em n c (k + 1) = 0 ∧ em n c (k + 2) = 0)

/-- The canonical partner (Fig. 11, Alg. 9): on side `k`, crossings `p < c k` cut corner `k`,
crossings `p ≥ n k − c (k+1)` cut corner `k + 1`, the others come from corner `k + 2`.
[formula: ic-face-decoding] -/
def canon (n c : Fin 3 → ℕ) (k : Fin 3) (p : ℕ) : Pt :=
  if p < c k then .side (k + 2) (n (k + 2) - 1 - p)
  else if p + c (k + 1) < n k then .corner (k + 2)
  else .side (k + 1) (n k - 1 - p)

/-! ### Rank arithmetic, one corner at a time -/

@[simp] lemma fin3_add_one_add_one (k : Fin 3) : k + 1 + 1 = k + 2 := by
  fin_cases k <;> rfl

@[simp] lemma fin3_add_two_add_one (k : Fin 3) : k + 2 + 1 = k := by
  fin_cases k <;> rfl

@[simp] lemma fin3_add_one_add_two (k : Fin 3) : k + 1 + 2 = k := by
  fin_cases k <;> rfl

@[simp] lemma fin3_add_two_add_two (k : Fin 3) : k + 2 + 2 = k + 1 := by
  fin_cases k <;> rfl

lemma fin3_add_one_ne (k : Fin 3) : k + 1 ≠ k := by fin_cases k <;> decide
lemma fin3_add_two_ne (k : Fin 3) : k + 2 ≠ k := by fin_cases k <;> decide
lemma fin3_add_two_ne_add_one (k : Fin 3) : k + 2 ≠ k + 1 := by fin_cases k <;> decide

/-- Simp set that turns ranks of concrete corners and sides into arithmetic. -/
lemma rank_corner (n : Fin 3 → ℕ) (k : Fin 3) : (Pt.corner k).rank n = cornerRank n k := rfl
lemma rank_side (n : Fin 3 → ℕ) (k : Fin 3) (p : ℕ) :
    (Pt.side k p).rank n = cornerRank n k + 1 + p := rfl

@[simp] lemma cornerRank_zero (n : Fin 3 → ℕ) : cornerRank n 0 = 0 := rfl
@[simp] lemma cornerRank_one (n : Fin 3 → ℕ) : cornerRank n 1 = n 0 + 1 := rfl
@[simp] lemma cornerRank_two (n : Fin 3 → ℕ) : cornerRank n 2 = n 0 + n 1 + 2 := rfl


/-! ### Monotonicity of the arc types along a side -/

/-- Close a goal about the ranks of boundary points of the sides around a corner variable `k`, by
cases on `k` and linear arithmetic. -/
macro "rank_omega" k:ident : tactic =>
  `(tactic| (fin_cases $k:ident <;> simp [cornerRank, Interleave, Btw, Pt.rank] at * <;> omega))

namespace FaceModel

variable {M : FaceModel}

/-- The partner of a crossing is on another side, with a position in range. -/
lemma Valid.side_lt (hM : M.Valid) {k j : Fin 3} {p q : ℕ} (hp : p < M.n k)
    (h : M.partner k p = .side j q) : q < M.n j ∧ (j = k + 2 ∨ j = k + 1) := by
  rcases hM.normal k p hp with ⟨q', hq', h'⟩ | ⟨q', hq', h'⟩ | h'
  · rw [h] at h'; cases h'; exact ⟨hq', Or.inl rfl⟩
  · rw [h] at h'; cases h'; exact ⟨hq', Or.inr rfl⟩
  · rw [h] at h'; cases h'

/-- Two crossings with the same partner crossing are equal. -/
lemma Valid.inj (hM : M.Valid) {k k' j : Fin 3} {p p' q : ℕ} (hp : p < M.n k) (hp' : p' < M.n k')
    (h : M.partner k p = .side j q) (h' : M.partner k' p' = .side j q) : k = k' ∧ p = p' := by
  have e1 := hM.symm k p j q hp h
  have e2 := hM.symm k' p' j q hp' h'
  rw [e1] at e2; cases e2; exact ⟨rfl, rfl⟩

/-- The interleaving condition for two crossings whose partners are known. -/
lemma Valid.disj (hM : M.Valid) {k k' : Fin 3} {p p' : ℕ} {x y : Pt} (hp : p < M.n k)
    (hp' : p' < M.n k') (hx : M.partner k p = x) (hy : M.partner k' p' = y) :
    ¬ Interleave ((Pt.side k p).rank M.n) (x.rank M.n) ((Pt.side k' p').rank M.n) (y.rank M.n) := by
  have := hM.disjoint k p k' p' hp hp'
  rwa [hx, hy] at this

end FaceModel

end IntegerCoords

end DdgProofs

namespace DdgProofs.IntegerCoords.FaceModel

variable {M : FaceModel}

/-- If crossing `p'` of side `k` cuts corner `k` (partner `q'` on side `k + 2`), so does every
earlier crossing `p < p'`, with a later partner `q > q'`. -/
lemma Valid.L_mono (hM : M.Valid) {k : Fin 3} {p p' q' : ℕ} (hpp : p < p') (hp' : p' < M.n k)
    (hq' : M.partner k p' = .side (k + 2) q') :
    ∃ q, q < M.n (k + 2) ∧ q' < q ∧ M.partner k p = .side (k + 2) q := by
  have hp : p < M.n k := by omega
  have hq'lt := (hM.side_lt hp' hq').1
  rcases hM.normal k p hp with ⟨q, hq, h⟩ | ⟨q, hq, h⟩ | h
  · refine ⟨q, hq, ?_, h⟩
    have hne : q ≠ q' := by
      rintro rfl; have := (hM.inj hp hp' h hq').2; omega
    have hd := hM.disj hp hp' h hq'
    rank_omega k
  · exfalso
    have hd := hM.disj hp hp' h hq'
    rank_omega k
  · exfalso
    have hd := hM.disj hp hp' h hq'
    rank_omega k

/-- If crossing `p` of side `k` cuts corner `k + 1` (partner `q` on side `k + 1`), so does every
later crossing `p' > p`, with an earlier partner `q' < q`. -/
lemma Valid.R_mono (hM : M.Valid) {k : Fin 3} {p p' q : ℕ} (hpp : p < p') (hp' : p' < M.n k)
    (hq : M.partner k p = .side (k + 1) q) :
    ∃ q', q' < M.n (k + 1) ∧ q' < q ∧ M.partner k p' = .side (k + 1) q' := by
  have hp : p < M.n k := by omega
  have hqlt := (hM.side_lt hp hq).1
  rcases hM.normal k p' hp' with ⟨q', hq', h⟩ | ⟨q', hq', h⟩ | h
  · exfalso
    have hd := hM.disj hp hp' hq h
    rank_omega k
  · refine ⟨q', hq', ?_, h⟩
    have hne : q' ≠ q := by
      rintro rfl; have := (hM.inj hp hp' hq h).2; omega
    have hd := hM.disj hp hp' hq h
    rank_omega k
  · exfalso
    have hd := hM.disj hp hp' hq h
    rank_omega k

/-- Crossing `p` of side `k` cuts corner `k` (its partner is on side `k + 2`). -/
def IsL (M : FaceModel) (k : Fin 3) (p : ℕ) : Prop := ∃ q, M.partner k p = .side (k + 2) q

/-- Crossing `p` of side `k` cuts corner `k + 1` (its partner is on side `k + 1`). -/
def IsR (M : FaceModel) (k : Fin 3) (p : ℕ) : Prop := ∃ q, M.partner k p = .side (k + 1) q

end DdgProofs.IntegerCoords.FaceModel

namespace DdgProofs.IntegerCoords

/-! ### Initial segments and strictly decreasing maps of natural numbers -/

/-- A predicate on `[0, n)` that is closed downwards holds exactly on an initial segment. -/
lemma exists_prefix (P : ℕ → Prop) (n : ℕ) (hP : ∀ p p', p < p' → p' < n → P p' → P p) :
    ∃ a ≤ n, ∀ p < n, (P p ↔ p < a) := by
  classical
  induction n with
  | zero => exact ⟨0, le_rfl, fun p hp => absurd hp (Nat.not_lt_zero _)⟩
  | succ n ih =>
    obtain ⟨a, ha, hiff⟩ := ih (fun p p' h h' => hP p p' h (by omega))
    by_cases hlt : a < n
    · refine ⟨a, by omega, fun p hp => ?_⟩
      rcases Nat.lt_succ_iff_lt_or_eq.mp hp with hp | rfl
      · exact hiff p hp
      · constructor
        · intro hPn
          have h1 := hP a p hlt (by omega) hPn
          have h2 := (hiff a hlt).1 h1
          omega
        · intro h; omega
    · have han : a = n := by omega
      subst han
      by_cases hPn : P a
      · refine ⟨a + 1, le_rfl, fun p hp => ?_⟩
        rcases Nat.lt_succ_iff_lt_or_eq.mp hp with hp | rfl
        · rw [hiff p hp]; omega
        · simp [hPn]
      · refine ⟨a, by omega, fun p hp => ?_⟩
        rcases Nat.lt_succ_iff_lt_or_eq.mp hp with hp | rfl
        · exact hiff p hp
        · simp [hPn]

/-- A predicate on `[0, n)` that is closed upwards holds exactly on a final segment. -/
lemma exists_suffix (P : ℕ → Prop) (n : ℕ) (hP : ∀ p p', p < p' → p' < n → P p → P p') :
    ∃ b ≤ n, ∀ p < n, (P p ↔ n - b ≤ p) := by
  obtain ⟨a, ha, hiff⟩ := exists_prefix (fun p => P (n - 1 - p)) n
    (fun p p' h h' hp' => hP (n - 1 - p') (n - 1 - p) (by omega) (by omega) hp')
  refine ⟨a, ha, fun p hp => ?_⟩
  have := hiff (n - 1 - p) (by omega)
  have e : n - 1 - (n - 1 - p) = p := by omega
  rw [e] at this
  rw [this]; omega

/-- A strictly decreasing map of `[0, a)` into `[0, N)` satisfies `f p + p < N`. -/
lemma anti_add_lt {f : ℕ → ℕ} {a N : ℕ} (hanti : ∀ p p', p < p' → p' < a → f p' < f p)
    (hlt : ∀ p < a, f p < N) : ∀ p < a, f p + p < N := by
  intro p hp
  induction p with
  | zero => simpa using hlt 0 hp
  | succ m ih =>
    have h1 := ih (by omega)
    have h2 := hanti m (m + 1) (by omega) hp
    omega

/-- A strictly decreasing map of `[0, a)` into `[L, ∞)` satisfies `f p ≥ L + (a − 1 − p)`. -/
lemma anti_ge {f : ℕ → ℕ} {a L : ℕ} (hanti : ∀ p p', p < p' → p' < a → f p' < f p)
    (hge : ∀ p < a, L ≤ f p) : ∀ p < a, L + (a - 1 - p) ≤ f p := by
  intro p hp
  induction h : a - 1 - p generalizing p with
  | zero => have := hge p hp; omega
  | succ m ih =>
    have h1 := ih (p + 1) (by omega) (by omega)
    have h2 := hanti p (p + 1) (by omega) (by omega)
    omega

end DdgProofs.IntegerCoords

namespace DdgProofs.IntegerCoords.FaceModel

variable {M : FaceModel}

/-- The position of a boundary point on its side (`0` for a corner). -/
def _root_.DdgProofs.IntegerCoords.Pt.pos : Pt → ℕ
  | .corner _ => 0
  | .side _ p => p

/-- **Blocks.** On every side the crossings that cut corner `k` form an initial segment of length
`a`, those that cut corner `k + 1` a final segment of length `b`, and `a + b ≤ n k`. -/
lemma Valid.blocks (hM : M.Valid) (k : Fin 3) :
    ∃ a b, a + b ≤ M.n k ∧ (∀ p < M.n k, M.IsL k p ↔ p < a) ∧
      (∀ p < M.n k, M.IsR k p ↔ M.n k - b ≤ p) := by
  obtain ⟨a, ha, hL⟩ := exists_prefix (M.IsL k) (M.n k) (fun p p' h h' ⟨q', hq'⟩ => by
    obtain ⟨q, -, -, hq⟩ := hM.L_mono h h' hq'
    exact ⟨q, hq⟩)
  obtain ⟨b, hb, hR⟩ := exists_suffix (M.IsR k) (M.n k) (fun p p' h h' ⟨q, hq⟩ => by
    obtain ⟨q', -, -, hq'⟩ := hM.R_mono h h' hq
    exact ⟨q', hq'⟩)
  refine ⟨a, b, ?_, hL, hR⟩
  by_contra hlt
  have hp : M.n k - b < M.n k := by omega
  obtain ⟨q, hq⟩ := (hL _ hp).2 (by omega)
  obtain ⟨q', hq'⟩ := (hR _ hp).2 le_rfl
  rw [hq] at hq'
  exact fin3_add_two_ne_add_one k (Pt.side.inj hq').1

/-- **Matching of the corner blocks.** The crossings `p < a` of side `k` that cut corner `k` are
joined to the last `a` crossings of side `k + 2` in reverse order, and these are exactly the
crossings of side `k + 2` that cut corner `k`. -/
lemma Valid.match_blocks (hM : M.Valid) {k : Fin 3} {a b : ℕ}
    (hL : ∀ p < M.n k, M.IsL k p ↔ p < a) (hab : a ≤ M.n k)
    (hR : ∀ p < M.n (k + 2), M.IsR (k + 2) p ↔ M.n (k + 2) - b ≤ p) (hb : b ≤ M.n (k + 2)) :
    a = b ∧ ∀ p < a, M.partner k p = .side (k + 2) (M.n (k + 2) - 1 - p) := by
  -- the map f : [0, a) → side k + 2
  have hf0 : ∀ p < a, ∃ q, M.partner k p = .side (k + 2) q := fun p hp => (hL p (by omega)).2 hp
  choose! f hf using hf0
  have hflt : ∀ p < a, f p < M.n (k + 2) := fun p hp => (hM.side_lt (by omega) (hf p hp)).1
  have hfge : ∀ p < a, M.n (k + 2) - b ≤ f p := by
    intro p hp
    have hs := hM.symm k p (k + 2) (f p) (by omega) (hf p hp)
    have hr : M.IsR (k + 2) (f p) := ⟨p, by rw [hs, fin3_add_two_add_one]⟩
    exact (hR (f p) (hflt p hp)).1 hr
  have hfanti : ∀ p p', p < p' → p' < a → f p' < f p := by
    intro p p' h h'
    obtain ⟨q, -, hlt, hq⟩ := hM.L_mono h (by omega) (hf p' h')
    rw [hf p (by omega)] at hq
    rw [(Pt.side.inj hq).2]; exact hlt
  -- the map g : the final block of side k + 2 → side k
  have hg0 : ∀ s < b, ∃ q, M.partner (k + 2) (M.n (k + 2) - b + s) = .side k q := by
    intro s hs
    obtain ⟨q, hq⟩ := (hR (M.n (k + 2) - b + s) (by omega)).2 (by omega)
    rw [fin3_add_two_add_one] at hq
    exact ⟨q, hq⟩
  choose! g hg using hg0
  have hglt : ∀ s < b, g s < a := by
    intro s hs
    have hlt := (hM.side_lt (by omega) (hg s hs)).1
    have hs' := hM.symm (k + 2) (M.n (k + 2) - b + s) k (g s) (by omega) (hg s hs)
    exact (hL (g s) hlt).1 ⟨_, hs'⟩
  have hganti : ∀ s s', s < s' → s' < b → g s' < g s := by
    intro s s' h h'
    have hgs : M.partner (k + 2) (M.n (k + 2) - b + s) = .side (k + 2 + 1) (g s) := by
      rw [fin3_add_two_add_one]; exact hg s (by omega)
    obtain ⟨q', -, hlt, hq'⟩ :=
      hM.R_mono (k := k + 2) (show M.n (k + 2) - b + s < M.n (k + 2) - b + s' by omega)
        (by omega) hgs
    rw [fin3_add_two_add_one, hg s' h'] at hq'
    rw [(Pt.side.inj hq').2]; exact hlt
  -- sizes
  have hba : b ≤ a := by
    rcases Nat.eq_zero_or_pos b with h0 | hpos
    · omega
    · have := anti_add_lt hganti hglt (b - 1) (by omega); omega
  have hab' : a ≤ b := by
    rcases Nat.eq_zero_or_pos a with h0 | hpos
    · omega
    · have := anti_add_lt (f := fun p => f p - (M.n (k + 2) - b)) (a := a) (N := b)
        (fun p p' h h' => by
          have := hfanti p p' h h'; have := hfge p (by omega); have := hfge p' h'; omega)
        (fun p hp => by have := hflt p hp; have := hfge p hp; omega) (a - 1) (by omega)
      omega
  refine ⟨le_antisymm hab' hba, fun p hp => ?_⟩
  have h1 := anti_add_lt hfanti hflt p hp
  have h2 := anti_ge hfanti hfge p hp
  have : f p = M.n (k + 2) - 1 - p := by omega
  rw [hf p hp, this]

/-- **The structure theorem.** A valid face model is the canonical model of its corner counts,
and these are admissible. [formula: ic-face-decoding] -/
theorem Valid.eq_canon (hM : M.Valid) :
    ∃ c, Admissible M.n c ∧ ∀ k p, p < M.n k → M.partner k p = canon M.n c k p := by
  choose a b hab hL hR using hM.blocks
  have hmatch : ∀ k, a k = b (k + 2) ∧
      ∀ p < a k, M.partner k p = .side (k + 2) (M.n (k + 2) - 1 - p) := fun k =>
    hM.match_blocks (hL k) (by have := hab k; omega) (hR (k + 2)) (by have := hab (k + 2); omega)
  have hb : ∀ k, b k = a (k + 1) := by
    intro k; have := (hmatch (k + 1)).1; rw [fin3_add_one_add_two] at this; exact this.symm
  -- the middle block comes from the opposite corner
  have hmid : ∀ k p, a k ≤ p → p + a (k + 1) < M.n k → M.partner k p = .corner (k + 2) := by
    intro k p h1 h2
    have hp : p < M.n k := by omega
    rcases hM.normal k p hp with ⟨q, -, h⟩ | ⟨q, -, h⟩ | h
    · have := (hL k p hp).1 ⟨q, h⟩; omega
    · have := (hR k p hp).1 ⟨q, h⟩; rw [hb] at this; omega
    · exact h
  have hadm1 : ∀ k, a k + a (k + 1) ≤ M.n k := fun k => by rw [← hb]; exact hab k
  -- the final block
  have hfin : ∀ k p, p < M.n k → M.n k ≤ p + a (k + 1) →
      M.partner k p = .side (k + 1) (M.n k - 1 - p) := by
    intro k p hp h
    obtain ⟨q, hq⟩ := (hR k p hp).2 (by rw [hb]; omega)
    have hqlt := (hM.side_lt hp hq).1
    have hql : M.IsL (k + 1) q := ⟨p, by rw [hM.symm k p (k + 1) q hp hq, fin3_add_one_add_two]⟩
    have hqa := (hL (k + 1) q hqlt).1 hql
    have e1 := (hmatch (k + 1)).2 q hqa
    rw [fin3_add_one_add_two, hM.symm k p (k + 1) q hp hq] at e1
    have e2 := (Pt.side.inj e1).2
    have := hadm1 k
    rw [hq]; congr 1; omega
  -- an emanating corner is cut by no arc
  have hadm2 : ∀ k, 0 < em M.n a k → a k = 0 := by
    intro k he
    unfold em at he
    by_contra hak
    have h1 : M.partner (k + 1) (a (k + 1)) = .corner k := by
      rw [hmid (k + 1) (a (k + 1)) le_rfl (by rw [fin3_add_one_add_one]; omega),
        fin3_add_one_add_two]
    have h2 := (hmatch k).2 0 (by omega)
    have hd := hM.disj (k := k + 1) (k' := k) (by omega)
      (show 0 < M.n k by have := hadm1 k; omega) h1 h2
    have hn2 : 0 < M.n (k + 2) := by have := hadm1 (k + 2); rw [fin3_add_two_add_one] at this; omega
    have hnk1 : a (k + 1) < M.n (k + 1) := by omega
    clear h1 h2 hmid hfin hmatch hL hR hab hb hadm1
    rank_omega k
  -- at most one corner emanates
  have hadm3 : ∀ k, 0 < em M.n a k → em M.n a (k + 1) = 0 ∧ em M.n a (k + 2) = 0 := by
    intro k he
    unfold em at he ⊢
    rw [fin3_add_one_add_one, fin3_add_one_add_two, fin3_add_two_add_one, fin3_add_two_add_two] at *
    have h1 : M.partner (k + 1) (a (k + 1)) = .corner k := by
      rw [hmid (k + 1) (a (k + 1)) le_rfl (by rw [fin3_add_one_add_one]; omega),
        fin3_add_one_add_two]
    constructor
    · by_contra hne
      have h2 : M.partner (k + 2) (a (k + 2)) = .corner (k + 1) := by
        rw [hmid (k + 2) (a (k + 2)) le_rfl (by rw [fin3_add_two_add_one]; omega),
          fin3_add_two_add_two]
      have hd := hM.disj (k := k + 1) (k' := k + 2) (by omega) (by omega) h1 h2
      clear h1 h2 hmid hfin hmatch hL hR hab hb hadm1 hadm2
      rank_omega k
    · by_contra hne
      have h2 : M.partner k (a k) = .corner (k + 2) :=
        hmid k (a k) le_rfl (by omega)
      have hd := hM.disj (k := k + 1) (k' := k) (by omega) (by omega) h1 h2
      clear h1 h2 hmid hfin hmatch hL hR hab hb hadm1 hadm2
      rank_omega k
  refine ⟨a, ⟨hadm1, hadm2, hadm3⟩, fun k p hp => ?_⟩
  unfold canon
  split_ifs with h1 h2
  · exact (hmatch k).2 p h1
  · exact hmid k p (by omega) h2
  · exact hfin k p hp (by omega)

end DdgProofs.IntegerCoords.FaceModel

namespace DdgProofs.IntegerCoords

lemma canon_normal {n c : Fin 3 → ℕ} (h : Admissible n c) (k : Fin 3) (p : ℕ) (hp : p < n k) :
    (∃ q, q < n (k + 2) ∧ canon n c k p = .side (k + 2) q) ∨
      (∃ q, q < n (k + 1) ∧ canon n c k p = .side (k + 1) q) ∨ canon n c k p = .corner (k + 2) := by
  have a1 := h.1 (k + 2)
  have a2 := h.1 (k + 1)
  rw [fin3_add_two_add_one] at a1
  rw [fin3_add_one_add_one] at a2
  unfold canon
  split_ifs with e1 e2
  · left; exact ⟨_, by omega, rfl⟩
  · right; right; rfl
  · right; left; exact ⟨_, by omega, rfl⟩

lemma canon_symm {n c : Fin 3 → ℕ} (h : Admissible n c) (k : Fin 3) (p : ℕ) (k' : Fin 3) (q : ℕ)
    (hp : p < n k) (hq : canon n c k p = .side k' q) : canon n c k' q = .side k p := by
  have a0 := h.1 k
  have a1 := h.1 (k + 2)
  have a2 := h.1 (k + 1)
  rw [fin3_add_two_add_one] at a1
  rw [fin3_add_one_add_one] at a2
  unfold canon at hq
  by_cases e1 : p < c k
  · simp only [e1, ↓reduceIte] at hq
    obtain ⟨rfl, rfl⟩ := Pt.side.inj hq
    unfold canon
    rw [fin3_add_two_add_one]
    split_ifs <;> first | omega | (congr 1; omega)
  · by_cases e2 : p + c (k + 1) < n k
    · simp only [e1, e2, ↓reduceIte] at hq; cases hq
    · simp only [e1, e2, ↓reduceIte] at hq
      obtain ⟨rfl, rfl⟩ := Pt.side.inj hq
      unfold canon
      rw [fin3_add_one_add_two]
      split_ifs <;> first | omega | (congr 1; omega)

lemma fin3_cases (k : Fin 3) : k = 0 ∨ k = 1 ∨ k = 2 := by fin_cases k <;> simp

/-- Fin 3 additions of literals. -/
lemma fin3_lits : (0 : Fin 3) + 1 = 1 ∧ (0 : Fin 3) + 2 = 2 ∧ (1 : Fin 3) + 1 = 2 ∧
    (1 : Fin 3) + 2 = 0 ∧ (2 : Fin 3) + 1 = 0 ∧ (2 : Fin 3) + 2 = 1 := by decide

lemma rank_canon (n c : Fin 3 → ℕ) (k : Fin 3) (p : ℕ) :
    (canon n c k p).rank n = if p < c k then cornerRank n (k + 2) + 1 + (n (k + 2) - 1 - p)
      else if p + c (k + 1) < n k then cornerRank n (k + 2)
      else cornerRank n (k + 1) + 1 + (n k - 1 - p) := by
  unfold canon; split_ifs <;> rfl

lemma canon_disjoint {n c : Fin 3 → ℕ} (h : Admissible n c) (k : Fin 3) (p : ℕ) (k' : Fin 3)
    (p' : ℕ) (hp : p < n k) (hp' : p' < n k') :
    ¬ Interleave ((Pt.side k p).rank n) ((canon n c k p).rank n) ((Pt.side k' p').rank n)
      ((canon n c k' p').rank n) := by
  obtain ⟨h1, h2, h3⟩ := h
  have a0 := h1 0; have a1 := h1 1; have a2 := h1 2
  have b0 := h2 0; have b1 := h2 1; have b2 := h2 2
  have d0 := h3 0; have d1 := h3 1; have d2 := h3 2
  simp only [em] at b0 b1 b2 d0 d1 d2
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [l1, l2, l3, l4, l5, l6] at a0 a1 a2 b0 b1 b2 d0 d1 d2
  clear h1 h2 h3
  rw [rank_canon, rank_canon]
  rcases fin3_cases k with rfl | rfl | rfl <;> rcases fin3_cases k' with rfl | rfl | rfl <;>
    simp only [Pt.rank, l1, l2, l3, l4, l5, l6, cornerRank_zero, cornerRank_one, cornerRank_two,
      Interleave, Btw] <;> split_ifs <;> omega

/-- **Every admissible canonical model is valid** (the existence half of `ic-face-decoding`).
[formula: ic-face-decoding] -/
theorem canon_valid {n c : Fin 3 → ℕ} (h : Admissible n c) : (⟨n, canon n c⟩ : FaceModel).Valid :=
  ⟨canon_normal h, canon_symm h, canon_disjoint h⟩

end DdgProofs.IntegerCoords
