import DdgProofs.IntegerCoords.FaceModel
import Mathlib.Algebra.BigOperators.Fin

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

* **`ic-face-split-normal`** — for every gap `g` of a valid face, its region in the oracle's labels
  (`gapRegion`) is one of the face's regions (`gapRegion_valid`), and the oracle's `splitCount` of
  that region is the number of arcs separating `g` from each corner (`sepCount_eq_splitCount`); every
  region label of the oracle is the region of some gap (`exists_gap_of_valid`). The three new faces
  `(m, k, k + 1)` are valid whichever region is chosen (`splitCount_faceValid`): the counts of SPEC
  I15 are valid by construction. With `faceSplit?` (`faceSplit?_eq`) this is the oracle's whole face
  split.
* **`ic-split-edge-normal`** in the same terms: a vertex on a side without crossings (an edge with
  `n = −1`) is separated from the opposite corner by every arc, and Eq. 6 is their number
  (`sepCount_split_edge`).
-/

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

/-- Gap `g` lies on side `k` (between corner `k` and corner `k + 1`, both included). -/
def OnSide (n : Fin 3 → ℕ) (k : Fin 3) (g : ℕ) : Prop :=
  cornerRank n k ≤ g ∧ g ≤ cornerRank n k + n k

instance (n : Fin 3 → ℕ) (k : Fin 3) (g : ℕ) : Decidable (OnSide n k g) := by
  unfold OnSide; infer_instance

/-- The number of arcs around corner `k` that separate gap `g` from corner `k`. -/
def tau (n c : Fin 3 → ℕ) (k : Fin 3) (g : ℕ) : ℕ :=
  if OnSide n k g then min (g - cornerRank n k) (c k)
  else if OnSide n (k + 2) g then min (c k) (n (k + 2) - (g - cornerRank n (k + 2)))
  else c k

/-- The number of arcs from corner `κ` that separate gap `g` from corner `κ + 1`. -/
def phi (n c : Fin 3 → ℕ) (κ : Fin 3) (g : ℕ) : ℕ :=
  if OnSide n (κ + 1) g then min (em n c κ) (g - cornerRank n (κ + 1) - c (κ + 1))
  else if OnSide n (κ + 2) g then em n c κ
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
  rcases fin3_cases k with rfl | rfl | rfl <;> rcases fin3_cases j with rfl | rfl | rfl <;>
    simp only [Separates, cornerArc, tau, OnSide, Pt.rank, GapIn, Btw, l2, l4, l6,
      cornerRank_zero, cornerRank_one, cornerRank_two, Fin.isValue, reduceIte,
      Fin.reduceEq] at ht ⊢ <;> split_ifs <;> omega

/-- **Separation by an arc from the emanating corner.** For `s < e_κ`, `fanArc κ s` separates gap
`g` from corner `κ + 1` exactly for `s < φ(g)`, from corner `κ + 2` exactly for the others, and
from corner `κ` (where it ends) never. -/
lemma sep_fanArc {n c : Fin 3 → ℕ} (hc : Admissible n c) (κ j : Fin 3) {g : ℕ}
    (hg : g < total n) {s : ℕ} (hs : s < em n c κ) :
    Separates (fanArc n c κ s) g (cornerRank n j) ↔
      (if j = κ + 1 then s < phi n c κ g else if j = κ + 2 then phi n c κ g ≤ s else False) := by
  obtain ⟨a0, a1, a2, -⟩ := hc.lits
  obtain ⟨m0, m1, m2⟩ := em_lits n c
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  unfold total at hg
  rcases fin3_cases κ with rfl | rfl | rfl <;> rcases fin3_cases j with rfl | rfl | rfl <;>
    simp only [Separates, fanArc, phi, OnSide, Pt.rank, GapIn, Btw, l1, l2, l3, l4, l5, l6,
      cornerRank_zero, cornerRank_one, cornerRank_two, Fin.isValue, reduceIte,
      Fin.reduceEq, m0, m1, m2] at hs ⊢ <;> first | (split_ifs <;> omega) | simp

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

end IntegerCoords

end DdgProofs
