import DdgProofs.IntegerCoords.FaceModel
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.IntervalCases

/-!
# Integer coordinates: face split (tier B, SPEC I15)

`crates/ddg-integer-coords/SPEC.md` §4.4, decision I15, §9.2; the paper Sec. 3.3, Supp. Alg. 5, 6.

A vertex `m` inserted in a face is joined to the three corners by new edges. For straight arcs, the
segment from `m` to corner `j` crosses exactly the arcs that separate `m` from `j` (an arc that ends
at `j` is met only at `j`). The new normal coordinate `n_{mj}` is therefore the number of arcs of the
face model that separate the region of `m` from corner `j`. Everything here is combinatorial:

* a **gap** `g < N` (`N = n₀ + n₁ + n₂ + 3` boundary points) is the piece of the boundary between the
  points of rank `g` and `g + 1` (cyclically); every region of the arc arrangement meets the boundary
  in at least one gap, so gaps represent the regions;
* an arc `{a, b}` **separates** gap `g` from the boundary point `r` (`Separates`) when `r` is not an
  end of it and exactly one of `g`, `r` lies between `a` and `b`;
* the arcs of a valid face model are the families `cornerArc k t` (`t < c_k`: crossing `t` of side
  `k` to crossing `n_{k+2} − 1 − t` of side `k + 2`) and `fanArc κ s` (`s < e_κ`: corner `κ` to
  crossing `c_{κ+1} + s` of side `κ + 1`), each arc exactly once (`FaceModel.Valid.arcs`, from the
  structure theorem);
* `sepCount n c g j` is the number of arcs that separate gap `g` from corner `j`.

Results.

* **`ic-face-split-normal`** — for every gap `g` of a valid face, the oracle's `faceSplit?` on the
  region of `g` (`gapRegion`, in the oracle's labels `corner k t` / `middle s`) returns the numbers
  of arcs that separate `g` from the three corners (`faceSplit?_gapRegion`; per corner
  `splitCount_gapRegion`, which also shows that the region is one of the face's). For every region
  label the three new faces `(k, k + 1, m)` are valid in the sense of Eq. 1 (`splitCount_faceValid`):
  the counts of SPEC I15 are valid by construction, whichever region floating point picks.
* Not proved here: that every label of the oracle's `regions` list is the region of some gap (checked
  by the oracle's `ic-selftest`), and the edge split as a face split at a
  gap on a side without crossings (Eq. 6 is proved arithmetically, `splitEdgeNormal_eq`).
-/

set_option linter.unusedSimpArgs false
set_option linter.unusedTactic false
set_option linter.unreachableTactic false
set_option linter.unnecessarySeqFocus false

namespace DdgProofs

namespace IntegerCoords

open Ddg.IntegerCoords (cornerNumerator cornerCoord emanating nPlus faceValid)

/-! ## Gaps, separation, the arc families -/

/-- Number of boundary points (and of gaps) of a face. -/
def total (n : Fin 3 → ℕ) : ℕ := n 0 + n 1 + n 2 + 3

/-- Gap `g` (between the boundary points of rank `g` and `g + 1`) lies between the ends `a`, `b` of
a chord. -/
def GapIn (a b g : ℕ) : Prop := (a ≤ g ∧ g < b) ∨ (b ≤ g ∧ g < a)

instance (a b g : ℕ) : Decidable (GapIn a b g) := by unfold GapIn; infer_instance

instance (a b x : ℕ) : Decidable (Btw a b x) := by unfold Btw; infer_instance

/-- The chord with ends `ab = (a, b)` separates gap `g` from the boundary point `r`. -/
def Separates (ab : ℕ × ℕ) (g r : ℕ) : Prop :=
  r ≠ ab.1 ∧ r ≠ ab.2 ∧ (GapIn ab.1 ab.2 g ↔ ¬ Btw ab.1 ab.2 r)

instance (ab : ℕ × ℕ) (g r : ℕ) : Decidable (Separates ab g r) := by
  unfold Separates; infer_instance

/-- The `t`-th arc around corner `k` (from the corner), as the ranks of its ends. -/
def cornerArc (n : Fin 3 → ℕ) (k : Fin 3) (t : ℕ) : ℕ × ℕ :=
  ((Pt.side k t).rank n, (Pt.side (k + 2) (n (k + 2) - 1 - t)).rank n)

/-- The `s`-th arc from corner `κ` to the opposite side (from corner `κ + 1`), as ranks. -/
def fanArc (n c : Fin 3 → ℕ) (κ : Fin 3) (s : ℕ) : ℕ × ℕ :=
  ((Pt.corner κ).rank n, (Pt.side (κ + 1) (c (κ + 1) + s)).rank n)

/-- The number of arcs of the face model with side counts `n` and corner counts `c` that separate
gap `g` from corner `j`. [formula: ic-face-split-normal] -/
def sepCount (n c : Fin 3 → ℕ) (g : ℕ) (j : Fin 3) : ℕ :=
  ∑ k : Fin 3, ((Finset.range (c k)).filter
      (fun t => Separates (cornerArc n k t) g (cornerRank n j))).card +
    ∑ κ : Fin 3, ((Finset.range (em n c κ)).filter
      (fun s => Separates (fanArc n c κ s) g (cornerRank n j))).card

/-- **The arcs of a valid face model are the two families, each arc once.** Crossing `t < c_k` of
side `k` is joined to crossing `n_{k+2} − 1 − t` of side `k + 2`; crossing `c_{κ+1} + s` (`s < e_κ`)
of side `κ + 1` is joined to corner `κ`; and every crossing of side `k` is one of: the first end of
an arc around corner `k` (`p < c_k`), a fan end (`p = c_k + s`, `s < e_{k+2}`), or the second end of
an arc around corner `k + 1` (`p = n_k − 1 − t`, `t < c_{k+1}`). [formula: ic-face-split-normal] -/
theorem FaceModel.Valid.arcs {M : FaceModel} (hM : M.Valid) :
    ∃ c, Admissible M.n c ∧
      (∀ k t, t < c k → M.partner k t = .side (k + 2) (M.n (k + 2) - 1 - t)) ∧
      (∀ κ s, s < em M.n c κ → M.partner (κ + 1) (c (κ + 1) + s) = .corner κ) ∧
      (∀ k p, p < M.n k → p < c k ∨ (∃ s, s < em M.n c (k + 2) ∧ p = c k + s) ∨
        (∃ t, t < c (k + 1) ∧ p = M.n k - 1 - t)) := by
  obtain ⟨c, hc, e⟩ := hM.eq_canon
  refine ⟨c, hc, fun k t ht => ?_, fun κ s hs => ?_, fun k p hp => ?_⟩
  · have := hc.le k
    rw [e k t (by omega)]; unfold canon; simp only [ht, ↓reduceIte]
  · have h1 := hc.le (κ + 1)
    unfold em at hs
    rw [fin3_add_one_add_one] at h1
    rw [e (κ + 1) _ (by omega)]; unfold canon
    rw [fin3_add_one_add_one, fin3_add_one_add_two]
    split_ifs <;> first | rfl | omega
  · have h1 := hc.le k
    by_cases h : p < c k
    · exact Or.inl h
    · by_cases h' : p + c (k + 1) < M.n k
      · refine Or.inr (Or.inl ⟨p - c k, ?_, by omega⟩)
        unfold em; rw [fin3_add_two_add_one, fin3_add_two_add_two]; omega
      · exact Or.inr (Or.inr ⟨M.n k - 1 - p, by omega, by omega⟩)

/-! ## Levels of a gap -/

/-- The side of the boundary that gap `g` lies on (its slot is `g − cornerRank n k`). -/
def gapSide (n : Fin 3 → ℕ) (g : ℕ) : Fin 3 :=
  if g ≤ n 0 then 0 else if g ≤ n 0 + n 1 + 1 then 1 else 2

lemma gapSide_cases (n : Fin 3 → ℕ) (g : ℕ) :
    (gapSide n g = 0 ∧ g ≤ n 0) ∨ (gapSide n g = 1 ∧ n 0 + 1 ≤ g ∧ g ≤ n 0 + n 1 + 1) ∨
      (gapSide n g = 2 ∧ n 0 + n 1 + 2 ≤ g) := by
  unfold gapSide; split_ifs <;> simp only [true_and, Fin.reduceEq, false_and, or_false, false_or] <;> omega

/-- The number of arcs around corner `k` that separate gap `g` from corner `k`. -/
def tau (n c : Fin 3 → ℕ) (k : Fin 3) (g : ℕ) : ℕ :=
  if gapSide n g = k then min (g - cornerRank n k) (c k)
  else if gapSide n g = k + 2 then min (c k) (n (k + 2) - (g - cornerRank n (k + 2)))
  else c k

/-- The number of arcs from corner `κ` that separate gap `g` from corner `κ + 1`. -/
def phi (n c : Fin 3 → ℕ) (κ : Fin 3) (g : ℕ) : ℕ :=
  if gapSide n g = κ + 1 then min (em n c κ) (g - cornerRank n (κ + 1) - c (κ + 1))
  else if gapSide n g = κ + 2 then em n c κ
  else 0

lemma tau_le (n c : Fin 3 → ℕ) (k : Fin 3) (g : ℕ) : tau n c k g ≤ c k := by
  unfold tau; split_ifs <;> omega

lemma phi_le (n c : Fin 3 → ℕ) (κ : Fin 3) (g : ℕ) : phi n c κ g ≤ em n c κ := by
  unfold phi; split_ifs <;> omega

/-- Admissibility at the three literal corners. -/
lemma Admissible.lits {n c : Fin 3 → ℕ} (h : Admissible n c) :
    c 0 + c 1 ≤ n 0 ∧ c 1 + c 2 ≤ n 1 ∧ c 2 + c 0 ≤ n 2 ∧
    (0 < em n c 0 → c 0 = 0) ∧ (0 < em n c 1 → c 1 = 0) ∧ (0 < em n c 2 → c 2 = 0) ∧
    (0 < em n c 0 → em n c 1 = 0 ∧ em n c 2 = 0) ∧ (0 < em n c 1 → em n c 2 = 0 ∧ em n c 0 = 0) ∧
    (0 < em n c 2 → em n c 0 = 0 ∧ em n c 1 = 0) := by
  obtain ⟨h1, h2, h3⟩ := h
  have d1 := h3 1
  have d2 := h3 2
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  rw [l3, l4] at d1
  rw [l5, l6] at d2
  exact ⟨h1 0, h1 1, by have := h1 2; rwa [l5] at this, h2 0, h2 1, h2 2, h3 0, d1, d2⟩

/-- The em values as explicit differences. -/
lemma em_lits (n c : Fin 3 → ℕ) :
    em n c 0 = n 1 - c 1 - c 2 ∧ em n c 1 = n 2 - c 2 - c 0 ∧ em n c 2 = n 0 - c 0 - c 1 := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  refine ⟨?_, ?_, ?_⟩ <;> simp only [em, l1, l2, l3, l4, l5, l6]

/-- The emanation cases of an admissible face with the oracle's decoding, as equalities (no
implications left for `omega` to split). -/
lemma Admissible.split_cases {n c : Fin 3 → ℕ} (hc : Admissible n c) {d : Ddg.IntegerCoords.Decoded}
    (hk : (em n c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0 ∧ d.kappa = 0 ∧ d.e = 0) ∨
        (0 < em n c 0 ∧ d.kappa = 0 ∧ d.e = em n c 0) ∨
        (0 < em n c 1 ∧ d.kappa = 1 ∧ d.e = em n c 1) ∨
        (0 < em n c 2 ∧ d.kappa = 2 ∧ d.e = em n c 2)) :
    (em n c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0 ∧ d.kappa = 0 ∧ d.e = 0) ∨
      (0 < em n c 0 ∧ c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0 ∧ d.kappa = 0 ∧ d.e = em n c 0) ∨
      (0 < em n c 1 ∧ c 1 = 0 ∧ em n c 2 = 0 ∧ em n c 0 = 0 ∧ d.kappa = 1 ∧ d.e = em n c 1) ∨
      (0 < em n c 2 ∧ c 2 = 0 ∧ em n c 0 = 0 ∧ em n c 1 = 0 ∧ d.kappa = 2 ∧ d.e = em n c 2) := by
  obtain ⟨-, -, -, b0, b1, b2, x0, x1, x2⟩ := hc.lits
  rcases hk with h | ⟨e, h1, h2⟩ | ⟨e, h1, h2⟩ | ⟨e, h1, h2⟩
  · exact Or.inl h
  · exact Or.inr (Or.inl ⟨e, b0 e, (x0 e).1, (x0 e).2, h1, h2⟩)
  · exact Or.inr (Or.inr (Or.inl ⟨e, b1 e, (x1 e).1, (x1 e).2, h1, h2⟩))
  · exact Or.inr (Or.inr (Or.inr ⟨e, b2 e, (x2 e).1, (x2 e).2, h1, h2⟩))

/-- The em values, additively (no truncated subtraction). -/
lemma Admissible.em_add {n c : Fin 3 → ℕ} (hc : Admissible n c) :
    em n c 0 + c 1 + c 2 = n 1 ∧ em n c 1 + c 2 + c 0 = n 2 ∧ em n c 2 + c 0 + c 1 = n 0 := by
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := em_lits n c
  omega

/-- **Separation by an arc around a corner.** For `t < c_k`, the arc `cornerArc k t` separates
gap `g` from corner `k` exactly for the `τ_k(g)` innermost arcs (`t < τ`), and from the two other
corners exactly for the others. -/
lemma sep_cornerArc {n c : Fin 3 → ℕ} (hc : Admissible n c) (k j : Fin 3) {g : ℕ}
    (hg : g < total n) {t : ℕ} (ht : t < c k) :
    Separates (cornerArc n k t) g (cornerRank n j) ↔
      (if j = k then t < tau n c k g else tau n c k g ≤ t) := by
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  unfold total at hg
  rcases gapSide_cases n g with ⟨hs, hgs⟩ | ⟨hs, hgs⟩ | ⟨hs, hgs⟩ <;>
  rcases fin3_cases k with rfl | rfl | rfl <;> rcases fin3_cases j with rfl | rfl | rfl <;>
    simp only [Separates, cornerArc, tau, hs, Pt.rank, GapIn, Btw, l2, l4, l6,
      cornerRank_zero, cornerRank_one, cornerRank_two, Fin.isValue, reduceIte,
      Fin.reduceEq] at ht ⊢ <;> omega

/-- **Separation by an arc from the emanating corner.** For `s < e_κ`, `fanArc κ s` separates gap
`g` from corner `κ + 1` exactly for `s < φ(g)`, from corner `κ + 2` exactly for the others, and
from corner `κ` (where it ends) never. -/
lemma sep_fanArc {n c : Fin 3 → ℕ} (hc : Admissible n c) (κ j : Fin 3) {g : ℕ}
    (_hg : g < total n) {s : ℕ} (hs : s < em n c κ) :
    Separates (fanArc n c κ s) g (cornerRank n j) ↔
      (if j = κ + 1 then s < phi n c κ g else if j = κ + 2 then phi n c κ g ≤ s else False) := by
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := em_lits n c
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  rcases gapSide_cases n g with ⟨hq, hgs⟩ | ⟨hq, hgs⟩ | ⟨hq, hgs⟩ <;>
  rcases fin3_cases κ with rfl | rfl | rfl <;> rcases fin3_cases j with rfl | rfl | rfl <;>
    simp only [Separates, fanArc, phi, hq, Pt.rank, GapIn, Btw, l1, l2, l3, l4, l5, l6,
      cornerRank_zero, cornerRank_one, cornerRank_two, Fin.isValue, reduceIte,
      Fin.reduceEq, m0, m1, m2] at hs ⊢ <;> first | omega | simp

lemma card_filter_lt (c τ : ℕ) (P : ℕ → Prop) [DecidablePred P] (hτ : τ ≤ c)
    (h : ∀ t < c, P t ↔ t < τ) : ((Finset.range c).filter P).card = τ := by
  rw [Finset.filter_congr (q := fun t => t < τ) (fun t ht => h t (Finset.mem_range.mp ht))]
  have : (Finset.range c).filter (fun t => t < τ) = Finset.range τ := by
    ext t; simp only [Finset.mem_filter, Finset.mem_range]; omega
  rw [this, Finset.card_range]

lemma card_filter_ge (c τ : ℕ) (P : ℕ → Prop) [DecidablePred P]
    (h : ∀ t < c, P t ↔ τ ≤ t) : ((Finset.range c).filter P).card = c - τ := by
  rw [Finset.filter_congr (q := fun t => τ ≤ t) (fun t ht => h t (Finset.mem_range.mp ht))]
  have : (Finset.range c).filter (fun t => τ ≤ t) = Finset.Ico τ c := by
    ext t; simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_Ico]; omega
  rw [this, Nat.card_Ico]

lemma card_filter_false (c : ℕ) (P : ℕ → Prop) [DecidablePred P]
    (h : ∀ t < c, ¬ P t) : ((Finset.range c).filter P).card = 0 := by
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  exact fun t ht => h t (Finset.mem_range.mp ht)

/-- The number of arcs around corner `k` that separate gap `g` from corner `j`. -/
lemma card_cornerArc {n c : Fin 3 → ℕ} (hc : Admissible n c) (k j : Fin 3) {g : ℕ}
    (hg : g < total n) :
    ((Finset.range (c k)).filter (fun t => Separates (cornerArc n k t) g (cornerRank n j))).card =
      if j = k then tau n c k g else c k - tau n c k g := by
  split_ifs with hjk
  · exact card_filter_lt _ _ _ (tau_le n c k g) (fun t ht => by
      rw [sep_cornerArc hc k j hg ht]; simp only [hjk, ↓reduceIte])
  · exact card_filter_ge _ _ _ (fun t ht => by rw [sep_cornerArc hc k j hg ht]; simp only [hjk, ↓reduceIte])

/-- The number of arcs from corner `κ` that separate gap `g` from corner `j`. -/
lemma card_fanArc {n c : Fin 3 → ℕ} (hc : Admissible n c) (κ j : Fin 3) {g : ℕ}
    (hg : g < total n) :
    ((Finset.range (em n c κ)).filter (fun s => Separates (fanArc n c κ s) g (cornerRank n j))).card =
      if j = κ + 1 then phi n c κ g else if j = κ + 2 then em n c κ - phi n c κ g else 0 := by
  split_ifs with h1 h2
  · exact card_filter_lt _ _ _ (phi_le n c κ g) (fun s hs => by
      rw [sep_fanArc hc κ j hg hs]; simp only [h1, ↓reduceIte])
  · exact card_filter_ge _ _ _ (fun s hs => by
      rw [sep_fanArc hc κ j hg hs, h2]; simp only [fin3_add_two_ne_add_one, ↓reduceIte])
  · exact card_filter_false _ _ (fun s hs => by
      rw [sep_fanArc hc κ j hg hs]; simp only [h1, h2, ↓reduceIte, not_false_eq_true])

/-- **The separation counts in closed form.** [formula: ic-face-split-normal] -/
theorem sepCount_eq {n c : Fin 3 → ℕ} (hc : Admissible n c) {g : ℕ} (hg : g < total n)
    (j : Fin 3) :
    sepCount n c g j =
      ∑ k : Fin 3, (if j = k then tau n c k g else c k - tau n c k g) +
        ∑ κ : Fin 3, (if j = κ + 1 then phi n c κ g
          else if j = κ + 2 then em n c κ - phi n c κ g else 0) := by
  unfold sepCount
  congr 1
  · exact Finset.sum_congr rfl (fun k _ => card_cornerArc hc k j hg)
  · exact Finset.sum_congr rfl (fun κ _ => card_fanArc hc κ j hg)

/-- **The region of a gap**, in the oracle's labels (`Ddg.IntegerCoords.Region`): on side `k` at
slot `q` (between crossings `q − 1` and `q`), the corner region `(k, q)` if `q < c_k`, the corner
region `(k + 1, n_k − q)` if `q > n_k − c_{k+1}`, else a middle region: the fan region `q − c_k` if
side `k` receives the fan (`κ = k + 2`), the last fan region `e` if the fan's corner is `k + 1`, the
first (`0`) if it is `k`, and the central region `0` if nothing emanates.
[formula: ic-face-split-normal] -/
def gapRegion (n c : Fin 3 → ℕ) (g : ℕ) : Ddg.IntegerCoords.Region :=
  let k := gapSide n g
  let q := g - cornerRank n k
  if q < c k then .corner k.val q
  else if n k < q + c (k + 1) then .corner (k + 1).val (n k - q)
  else .middle (if 0 < em n c (k + 2) then q - c k
    else if 0 < em n c (k + 1) then em n c (k + 1) else 0)

/-! ## The oracle's decoded face -/

lemma oracle_side (n : Fin 3 → ℕ) :
    Ddg.IntegerCoords.side (n 0) (n 1) (n 2) 0 = n 0 ∧ Ddg.IntegerCoords.side (n 0) (n 1) (n 2) 1 = n 1 ∧
      Ddg.IntegerCoords.side (n 0) (n 1) (n 2) 2 = n 2 := by
  simp [Ddg.IntegerCoords.side]

/-- The oracle's per-corner Eq. 1 and 2 on an admissible face are `c` and `em`. -/
lemma oracle_corner {n c : Fin 3 → ℕ} (hc : Admissible n c) :
    (Ddg.IntegerCoords.faceCorner (n 0) (n 1) (n 2) 0 = c 0 ∧
      Ddg.IntegerCoords.faceCorner (n 0) (n 1) (n 2) 1 = c 1 ∧
      Ddg.IntegerCoords.faceCorner (n 0) (n 1) (n 2) 2 = c 2) ∧
    (Ddg.IntegerCoords.faceEmanating (n 0) (n 1) (n 2) 0 = em n c 0 ∧
      Ddg.IntegerCoords.faceEmanating (n 0) (n 1) (n 2) 1 = em n c 1 ∧
      Ddg.IntegerCoords.faceEmanating (n 0) (n 1) (n 2) 2 = em n c 2) := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  have h0 := hc.coords 0
  have h1 := hc.coords 1
  have h2 := hc.coords 2
  simp only [l1, l2, l3, l4, l5, l6] at h0 h1 h2
  simp only [Ddg.IntegerCoords.faceCorner, Ddg.IntegerCoords.faceEmanating, Ddg.IntegerCoords.side,
    Ddg.IntegerCoords.next, Ddg.IntegerCoords.prev]
  norm_num
  exact ⟨⟨h0.1, h1.1, h2.1⟩, ⟨h0.2, h1.2, h2.2⟩⟩

/-- **The oracle's decoding of an admissible face**: `n⁺`, `c` per corner, and the emanating corner
with its count (`κ = 0`, `e = 0` if none). -/
lemma oracle_decode {n c : Fin 3 → ℕ} (hc : Admissible n c) :
    ∃ d, Ddg.IntegerCoords.decode? (n 0) (n 1) (n 2) = some d ∧
      d.n 0 = n 0 ∧ d.n 1 = n 1 ∧ d.n 2 = n 2 ∧ d.c 0 = c 0 ∧ d.c 1 = c 1 ∧ d.c 2 = c 2 ∧
      ((em n c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0 ∧ d.kappa = 0 ∧ d.e = 0) ∨
        (0 < em n c 0 ∧ d.kappa = 0 ∧ d.e = em n c 0) ∨
        (0 < em n c 1 ∧ d.kappa = 1 ∧ d.e = em n c 1) ∨
        (0 < em n c 2 ∧ d.kappa = 2 ∧ d.e = em n c 2)) := by
  have hv := (FaceModel.exists_valid_iff n).1 ⟨_, rfl, canon_valid hc⟩
  obtain ⟨⟨k0, k1, k2⟩, ⟨f0, f1, f2⟩⟩ := oracle_corner hc
  obtain ⟨s0, s1, s2⟩ := oracle_side n
  have np : ∀ m : ℕ, nPlus (m : ℤ) = m := fun m => by simp only [nPlus]; omega
  unfold Ddg.IntegerCoords.decode?
  simp only [hv, Bool.not_true, Bool.false_eq_true, ↓reduceIte, List.map_cons, List.map_nil,
    List.find?_cons, f0, f1, f2, Int.toNat_natCast]
  by_cases e0 : 0 < em n c 0
  · simp only [e0, decide_true]
    refine ⟨_, rfl, ?_⟩
    simp only [s0, s1, s2, np, k0, k1, k2, Int.toNat_natCast, true_and, and_true, true_or,
      or_true, and_false, false_and, or_false, false_or, OfNat.ofNat_ne_zero, zero_ne_one] <;>
      omega
  · by_cases e1 : 0 < em n c 1
    · simp only [e0, e1, decide_true, decide_false]
      refine ⟨_, rfl, ?_⟩
      simp only [s0, s1, s2, np, k0, k1, k2, Int.toNat_natCast, true_and, and_true, true_or,
      or_true, and_false, false_and, or_false, false_or, OfNat.ofNat_ne_zero, zero_ne_one] <;>
      omega
    · by_cases e2 : 0 < em n c 2
      · simp only [e0, e1, e2, decide_true, decide_false]
        refine ⟨_, rfl, ?_⟩
        simp only [s0, s1, s2, np, k0, k1, k2, Int.toNat_natCast, true_and, and_true, true_or,
      or_true, and_false, false_and, or_false, false_or, OfNat.ofNat_ne_zero, zero_ne_one] <;>
        omega
      · simp only [e0, e1, e2, decide_false, List.find?_nil]
        refine ⟨_, rfl, ?_⟩
        simp only [s0, s1, s2, np, k0, k1, k2, Int.toNat_natCast, true_and, and_true, true_or,
      or_true, and_false, false_and, or_false, false_or, OfNat.ofNat_ne_zero, zero_ne_one] <;>
        omega

/-! ## The face split -/

/-- The levels and the region of a gap on side 0. -/
lemma values_side0 {n c : Fin 3 → ℕ} {g : ℕ} (h : gapSide n g = 0) :
    tau n c 0 g = min g (c 0) ∧ tau n c 1 g = min (c 1) (n 0 - g) ∧ tau n c 2 g = c 2 ∧
    phi n c 0 g = 0 ∧ phi n c 1 g = em n c 1 ∧ phi n c 2 g = min (em n c 2) (g - c 0) ∧
    gapRegion n c g = if g < c 0 then .corner 0 g else if n 0 < g + c 1 then .corner 1 (n 0 - g)
      else .middle (if 0 < em n c 2 then g - c 0 else if 0 < em n c 1 then em n c 1 else 0) := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [tau, phi, gapRegion, h, l1, l2, l3, l4, l5, l6, cornerRank_zero, Fin.isValue,
    Fin.reduceEq, reduceIte, Nat.sub_zero, Fin.val_zero, Fin.val_one, true_and]

lemma values_side1 {n c : Fin 3 → ℕ} {g : ℕ} (h : gapSide n g = 1) :
    tau n c 0 g = c 0 ∧ tau n c 1 g = min (g - (n 0 + 1)) (c 1) ∧
    tau n c 2 g = min (c 2) (n 1 - (g - (n 0 + 1))) ∧
    phi n c 0 g = min (em n c 0) (g - (n 0 + 1) - c 1) ∧ phi n c 1 g = 0 ∧ phi n c 2 g = em n c 2 ∧
    gapRegion n c g = if g - (n 0 + 1) < c 1 then .corner 1 (g - (n 0 + 1))
      else if n 1 < g - (n 0 + 1) + c 2 then .corner 2 (n 1 - (g - (n 0 + 1)))
      else .middle (if 0 < em n c 0 then g - (n 0 + 1) - c 1
        else if 0 < em n c 2 then em n c 2 else 0) := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [tau, phi, gapRegion, h, l1, l2, l3, l4, l5, l6, cornerRank_one, Fin.isValue,
    Fin.reduceEq, reduceIte, Fin.val_one, Fin.val_two, true_and]

lemma values_side2 {n c : Fin 3 → ℕ} {g : ℕ} (h : gapSide n g = 2) :
    tau n c 0 g = min (c 0) (n 2 - (g - (n 0 + n 1 + 2))) ∧ tau n c 1 g = c 1 ∧
    tau n c 2 g = min (g - (n 0 + n 1 + 2)) (c 2) ∧
    phi n c 0 g = em n c 0 ∧ phi n c 1 g = min (em n c 1) (g - (n 0 + n 1 + 2) - c 2) ∧
    phi n c 2 g = 0 ∧
    gapRegion n c g = if g - (n 0 + n 1 + 2) < c 2 then .corner 2 (g - (n 0 + n 1 + 2))
      else if n 2 < g - (n 0 + n 1 + 2) + c 0 then .corner 0 (n 2 - (g - (n 0 + n 1 + 2)))
      else .middle (if 0 < em n c 1 then g - (n 0 + n 1 + 2) - c 2
        else if 0 < em n c 0 then em n c 0 else 0) := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [tau, phi, gapRegion, h, l1, l2, l3, l4, l5, l6, cornerRank_two, Fin.isValue,
    Fin.reduceEq, reduceIte, Fin.val_two, Fin.val_zero, true_and]

set_option maxHeartbeats 2000000 in
lemma splitCount_gapRegion_side0 {n c : Fin 3 → ℕ} (hc : Admissible n c) {g : ℕ} (hg : g < total n)
    {d : Ddg.IntegerCoords.Decoded} (dc0 : d.c 0 = c 0) (dc1 : d.c 1 = c 1) (dc2 : d.c 2 = c 2)
    (hk : (em n c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0 ∧ d.kappa = 0 ∧ d.e = 0) ∨
        (0 < em n c 0 ∧ d.kappa = 0 ∧ d.e = em n c 0) ∨
        (0 < em n c 1 ∧ d.kappa = 1 ∧ d.e = em n c 1) ∨
        (0 < em n c 2 ∧ d.kappa = 2 ∧ d.e = em n c 2)) (j : Fin 3)
    (hs : gapSide n g = 0) (hgs : g ≤ n 0) :
    Ddg.IntegerCoords.splitCount d (gapRegion n c g) j.val = sepCount n c g j ∧
      (gapRegion n c g).valid d = true := by
  rw [sepCount_eq hc hg]
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := hc.em_add
  have hk := hc.split_cases hk
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  unfold total at hg
  simp only [Fin.sum_univ_three]
  obtain ⟨t0, t1, t2, f0, f1, f2, hr⟩ := values_side0 (c := c) hs
  rw [t0, t1, t2, f0, f1, f2, hr]
  clear t0 t1 t2 f0 f1 f2 hr
  rcases hk with ⟨e0, e1, e2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ |
    ⟨e0, z0, z1, z2, hκ, he⟩ <;>
  rcases fin3_cases j with rfl | rfl | rfl <;>
  simp only [l1, l2, l3, l4, l5, l6, Fin.isValue, Fin.reduceEq, reduceIte, Fin.val_zero,
    Fin.val_one, Fin.val_two] <;>
  split_ifs <;>
  simp only [Ddg.IntegerCoords.splitCount, Ddg.IntegerCoords.Region.valid,
    Ddg.IntegerCoords.next, Ddg.IntegerCoords.prev, hκ, he, dc0, dc1, dc2, beq_iff_eq,
    Bool.and_eq_true, decide_eq_true_eq, Bool.or_eq_true] <;>
  (try split_ifs) <;> first | omega | contradiction

set_option maxHeartbeats 2000000 in
lemma splitCount_gapRegion_side1 {n c : Fin 3 → ℕ} (hc : Admissible n c) {g : ℕ} (hg : g < total n)
    {d : Ddg.IntegerCoords.Decoded} (dc0 : d.c 0 = c 0) (dc1 : d.c 1 = c 1) (dc2 : d.c 2 = c 2)
    (hk : (em n c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0 ∧ d.kappa = 0 ∧ d.e = 0) ∨
        (0 < em n c 0 ∧ d.kappa = 0 ∧ d.e = em n c 0) ∨
        (0 < em n c 1 ∧ d.kappa = 1 ∧ d.e = em n c 1) ∨
        (0 < em n c 2 ∧ d.kappa = 2 ∧ d.e = em n c 2)) (j : Fin 3)
    (hs : gapSide n g = 1) (hgs : n 0 + 1 ≤ g ∧ g ≤ n 0 + n 1 + 1) :
    Ddg.IntegerCoords.splitCount d (gapRegion n c g) j.val = sepCount n c g j ∧
      (gapRegion n c g).valid d = true := by
  rw [sepCount_eq hc hg]
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := hc.em_add
  have hk := hc.split_cases hk
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  unfold total at hg
  simp only [Fin.sum_univ_three]
  obtain ⟨t0, t1, t2, f0, f1, f2, hr⟩ := values_side1 (c := c) hs
  rw [t0, t1, t2, f0, f1, f2, hr]
  clear t0 t1 t2 f0 f1 f2 hr
  rcases hk with ⟨e0, e1, e2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ |
    ⟨e0, z0, z1, z2, hκ, he⟩ <;>
  rcases fin3_cases j with rfl | rfl | rfl <;>
  simp only [l1, l2, l3, l4, l5, l6, Fin.isValue, Fin.reduceEq, reduceIte, Fin.val_zero,
    Fin.val_one, Fin.val_two] <;>
  split_ifs <;>
  simp only [Ddg.IntegerCoords.splitCount, Ddg.IntegerCoords.Region.valid,
    Ddg.IntegerCoords.next, Ddg.IntegerCoords.prev, hκ, he, dc0, dc1, dc2, beq_iff_eq,
    Bool.and_eq_true, decide_eq_true_eq, Bool.or_eq_true] <;>
  (try split_ifs) <;> first | omega | contradiction

set_option maxHeartbeats 2000000 in
lemma splitCount_gapRegion_side2 {n c : Fin 3 → ℕ} (hc : Admissible n c) {g : ℕ} (hg : g < total n)
    {d : Ddg.IntegerCoords.Decoded} (dc0 : d.c 0 = c 0) (dc1 : d.c 1 = c 1) (dc2 : d.c 2 = c 2)
    (hk : (em n c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0 ∧ d.kappa = 0 ∧ d.e = 0) ∨
        (0 < em n c 0 ∧ d.kappa = 0 ∧ d.e = em n c 0) ∨
        (0 < em n c 1 ∧ d.kappa = 1 ∧ d.e = em n c 1) ∨
        (0 < em n c 2 ∧ d.kappa = 2 ∧ d.e = em n c 2)) (j : Fin 3)
    (hs : gapSide n g = 2) (hgs : n 0 + n 1 + 2 ≤ g) :
    Ddg.IntegerCoords.splitCount d (gapRegion n c g) j.val = sepCount n c g j ∧
      (gapRegion n c g).valid d = true := by
  rw [sepCount_eq hc hg]
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := hc.em_add
  have hk := hc.split_cases hk
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  unfold total at hg
  simp only [Fin.sum_univ_three]
  obtain ⟨t0, t1, t2, f0, f1, f2, hr⟩ := values_side2 (c := c) hs
  rw [t0, t1, t2, f0, f1, f2, hr]
  clear t0 t1 t2 f0 f1 f2 hr
  rcases hk with ⟨e0, e1, e2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ |
    ⟨e0, z0, z1, z2, hκ, he⟩ <;>
  rcases fin3_cases j with rfl | rfl | rfl <;>
  simp only [l1, l2, l3, l4, l5, l6, Fin.isValue, Fin.reduceEq, reduceIte, Fin.val_zero,
    Fin.val_one, Fin.val_two] <;>
  split_ifs <;>
  simp only [Ddg.IntegerCoords.splitCount, Ddg.IntegerCoords.Region.valid,
    Ddg.IntegerCoords.next, Ddg.IntegerCoords.prev, hκ, he, dc0, dc1, dc2, beq_iff_eq,
    Bool.and_eq_true, decide_eq_true_eq, Bool.or_eq_true] <;>
  (try split_ifs) <;> first | omega | contradiction

/-- The oracle's counts for the region of a gap are the separation counts, and the region is one of
the face's. -/
theorem splitCount_gapRegion {n c : Fin 3 → ℕ} (hc : Admissible n c) {g : ℕ} (hg : g < total n)
    {d : Ddg.IntegerCoords.Decoded} (dc0 : d.c 0 = c 0) (dc1 : d.c 1 = c 1) (dc2 : d.c 2 = c 2)
    (hk : (em n c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0 ∧ d.kappa = 0 ∧ d.e = 0) ∨
        (0 < em n c 0 ∧ d.kappa = 0 ∧ d.e = em n c 0) ∨
        (0 < em n c 1 ∧ d.kappa = 1 ∧ d.e = em n c 1) ∨
        (0 < em n c 2 ∧ d.kappa = 2 ∧ d.e = em n c 2)) (j : Fin 3) :
    Ddg.IntegerCoords.splitCount d (gapRegion n c g) j.val = sepCount n c g j ∧
      (gapRegion n c g).valid d = true := by
  rcases gapSide_cases n g with ⟨hs, hgs⟩ | ⟨hs, hgs⟩ | ⟨hs, hgs⟩
  · exact splitCount_gapRegion_side0 hc hg dc0 dc1 dc2 hk j hs hgs
  · exact splitCount_gapRegion_side1 hc hg dc0 dc1 dc2 hk j hs hgs
  · exact splitCount_gapRegion_side2 hc hg dc0 dc1 dc2 hk j hs hgs

/-- **`ic-face-split-normal`: the oracle's face split is the count of separating arcs.** For a valid
face (admissible corner counts `c`, e.g. those of a valid face model by the structure theorem) and any
gap `g` of its boundary, the oracle's `faceSplit?` on the region of `g` returns the numbers of arcs
that separate `g` from the corners `0, 1, 2` — the normal coordinates of the edges from a vertex
inserted in that region to the corners, joined by arcs disjoint from the others (straight segments).
[formula: ic-face-split-normal] -/
theorem faceSplit?_gapRegion {n c : Fin 3 → ℕ} (hc : Admissible n c) {g : ℕ} (hg : g < total n) :
    Ddg.IntegerCoords.faceSplit? (n 0) (n 1) (n 2) (gapRegion n c g) =
      some (sepCount n c g 0, sepCount n c g 1, sepCount n c g 2) := by
  obtain ⟨d, hd, -, -, -, dc0, dc1, dc2, hk⟩ := oracle_decode hc
  have h0 := splitCount_gapRegion hc hg dc0 dc1 dc2 hk 0
  have h1 := splitCount_gapRegion hc hg dc0 dc1 dc2 hk 1
  have h2 := splitCount_gapRegion hc hg dc0 dc1 dc2 hk 2
  simp only [Fin.val_zero, Fin.val_one, Fin.val_two] at h0 h1 h2
  unfold Ddg.IntegerCoords.faceSplit?
  rw [hd]
  simp only [h0.2, ↓reduceIte, h0.1, h1.1, h2.1]

set_option maxHeartbeats 4000000 in
/-- **The three new faces are valid, whichever region is chosen** (SPEC I15). For every region `R`
of a valid face, the new face `(k, k + 1, m)` with the sides `n_k` (old), `n_{m,k+1}`, `n_{mk}`
(the oracle's `splitCount`) is valid in the sense of Eq. 1 (`faceValid`). [formula: ic-face-split-normal] -/
theorem splitCount_faceValid {n c : Fin 3 → ℕ} (hc : Admissible n c)
    {d : Ddg.IntegerCoords.Decoded} (dc0 : d.c 0 = c 0) (dc1 : d.c 1 = c 1) (dc2 : d.c 2 = c 2)
    (hk : (em n c 0 = 0 ∧ em n c 1 = 0 ∧ em n c 2 = 0 ∧ d.kappa = 0 ∧ d.e = 0) ∨
        (0 < em n c 0 ∧ d.kappa = 0 ∧ d.e = em n c 0) ∨
        (0 < em n c 1 ∧ d.kappa = 1 ∧ d.e = em n c 1) ∨
        (0 < em n c 2 ∧ d.kappa = 2 ∧ d.e = em n c 2))
    (R : Ddg.IntegerCoords.Region) (hR : R.valid d = true) (k : Fin 3) :
    faceValid (n k) (Ddg.IntegerCoords.splitCount d R (k + 1).val)
      (Ddg.IntegerCoords.splitCount d R k.val) = true := by
  rw [faceValid_iff]
  have np : ∀ m : ℕ, nPlus (m : ℤ) = m := fun m => by simp only [nPlus]; omega
  simp only [emanating, np]
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := hc.em_add
  have hk := hc.split_cases hk
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  rcases R with ⟨k0, t⟩ | s
  · simp only [Ddg.IntegerCoords.Region.valid, Bool.and_eq_true, decide_eq_true_eq] at hR
    obtain ⟨hk0, ht⟩ := hR
    interval_cases k0 <;>
    rcases hk with ⟨e0, e1, e2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ |
      ⟨e0, z0, z1, z2, hκ, he⟩ <;>
    rcases fin3_cases k with rfl | rfl | rfl <;>
    simp only [l1, l2, l3, l4, l5, l6, Fin.isValue, Fin.val_zero, Fin.val_one, Fin.val_two,
      Ddg.IntegerCoords.splitCount, Ddg.IntegerCoords.next, Ddg.IntegerCoords.prev, hκ, he,
      dc0, dc1, dc2, beq_iff_eq, Bool.and_eq_true, Bool.or_eq_true] at ht ⊢ <;>
    (try split_ifs) <;> (try contradiction) <;> omega
  · simp only [Ddg.IntegerCoords.Region.valid, decide_eq_true_eq] at hR
    rcases hk with ⟨e0, e1, e2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ | ⟨e0, z0, z1, z2, hκ, he⟩ |
      ⟨e0, z0, z1, z2, hκ, he⟩ <;>
    rcases fin3_cases k with rfl | rfl | rfl <;>
    simp only [l1, l2, l3, l4, l5, l6, Fin.isValue, Fin.val_zero, Fin.val_one, Fin.val_two,
      Ddg.IntegerCoords.splitCount, Ddg.IntegerCoords.next, Ddg.IntegerCoords.prev, hκ, he,
      dc0, dc1, dc2, beq_iff_eq, Bool.and_eq_true, Bool.or_eq_true] at hR ⊢ <;>
    (try split_ifs) <;> (try contradiction) <;> omega

end IntegerCoords

end DdgProofs
