import DdgProofs.IntegerCoords.Complex.Basic

/-!
# Integer coordinates, tier C: `TraceFrom` terminates (`ic-trace-terminates`)

`crates/ddg-integer-coords/SPEC.md` §3.1, §9.2; Supp. Alg. 1–2; the Rust `trace_from`,
`trace_with`, `extract_curve` (`crates/ddg-integer-coords/src/trace.rs`).

`iter z m` is the `m`-th crossing of the walk from `z` (`none` once the curve has ended at a
corner). The walk is the Rust `trace_with` loop: by `Valid.next_eq_oracle` each step is the
oracle's `traceStep` on the integers.

* `Valid.trace_terminates`: from any crossing, the walk either **ends at a corner** after `m`
  crossings, or **returns to its start** after `m` crossings (a closed curve); in both cases the
  `m` crossings visited are pairwise different *and* no crossing is visited in both directions
  (the Rust tracer's test, SPEC I18: a repeated `(halfedge, index)` pair in either direction), so
  `2m ≤ Σ_sides n⁺`, i.e. `m ≤ Σ_edges n⁺` (each interior edge is two sides; boundary edges
  have `n = −1`). The first repetition, if any, is the start crossing itself, in the same
  direction: exactly the dichotomy that SPEC §3.1 relies on ("a walk either reaches a vertex or
  returns to its start, which is detected") to avoid the global bound.
* `Valid.trace_from_end`: from the first crossing of a curve leaving a vertex (`IsEnd`, every
  arc of `input_directions`/`curve_ends`) the walk always ends at a corner.
* `Valid.trace_terminates_of_noClosed`: in a state without closed curves (SPEC §2.4 item 4: every
  crossing belongs to the curve of an input edge, which ends at vertices; what `certify` checks)
  the walk from **every** crossing ends at a corner. Closed normal curves exist in complex models
  in general (a torus), so this hypothesis cannot be dropped; it is not a property of the faces.

The proofs use only that the step map is injective and commutes with reversal
(`Valid.next_inj`, `Valid.next_rev`), which come from the face pairing and the gluing being
involutions (`DESIGN.md` item 2).
-/

namespace DDG.IntegerCoords

namespace Cx

variable {F : Type} (C : Cx F)

/-- The `m`-th crossing of the walk from `z` (`none` after the curve has ended at a corner). -/
def iter (z : Xing F) : ℕ → Option (Xing F)
  | 0 => some z
  | m + 1 => (iter z m).bind fun y => (C.next y).getLeft?

/-- All crossings, as a finite set. -/
def xings [Fintype F] : Finset (Xing F) :=
  (Finset.univ.sigma fun s : F × Fin 3 => Finset.range (C.np s)).map
    (Equiv.sigmaEquivProd (F × Fin 3) ℕ).toEmbedding

/-- No closed curves: no crossing comes back to itself (SPEC §2.4 item 4: every crossing lies on
the curve of an input edge between two vertices). -/
def NoClosed : Prop := ∀ z, C.IsX z → ∀ m, 0 < m → C.iter z m ≠ some z

/-- The walk from `z` visits `m ≥ 1` crossings, pairwise different and never one of them in both
directions; hence `2m ≤ Σ n⁺` over the sides. -/
structure Visits [Fintype F] (z : Xing F) (m : ℕ) : Prop where
  pos : 0 < m
  defined : ∀ i < m, ∃ y, C.iter z i = some y ∧ C.IsX y
  distinct : ∀ i j y y', i < m → j < m → C.iter z i = some y → C.iter z j = some y' →
    (y = y' → i = j) ∧ y ≠ C.rev y'
  bound : 2 * m ≤ ∑ s, C.np s

variable {C}

lemma card_xings [Fintype F] : C.xings.card = ∑ s, C.np s := by
  simp [xings, Finset.card_sigma]

lemma mem_xings [Fintype F] {z : Xing F} : z ∈ C.xings ↔ C.IsX z := by
  constructor
  · intro h
    simp only [xings, Finset.mem_map, Finset.mem_sigma, Finset.mem_univ, Finset.mem_range,
      true_and] at h
    obtain ⟨⟨s, p⟩, hp, rfl⟩ := h
    exact hp
  · intro h
    simp only [xings, Finset.mem_map, Finset.mem_sigma, Finset.mem_univ, Finset.mem_range,
      true_and]
    exact ⟨⟨z.1, z.2⟩, h, rfl⟩

lemma iter_succ {z : Xing F} {m : ℕ} :
    C.iter z (m + 1) = (C.iter z m).bind fun y => (C.next y).getLeft? := rfl

lemma iter_succ_of {z y y' : Xing F} {m : ℕ} (h : C.iter z m = some y) (hn : C.next y = .inl y') :
    C.iter z (m + 1) = some y' := by
  rw [iter_succ, h]; simp [hn]

lemma iter_succ_some {z y' : Xing F} {m : ℕ} (h : C.iter z (m + 1) = some y') :
    ∃ y, C.iter z m = some y ∧ C.next y = .inl y' := by
  rw [iter_succ] at h
  rcases e : C.iter z m with _ | y
  · rw [e] at h; cases h
  · rw [e] at h
    refine ⟨y, rfl, ?_⟩
    simp only [Option.bind_some] at h
    rcases e' : C.next y with a | b
    · rw [e'] at h; simp at h; rw [h]
    · rw [e'] at h; simp at h

lemma iter_succ_none {z y : Xing F} {m : ℕ} (h : C.iter z m = some y)
    (h' : C.iter z (m + 1) = none) : ∃ c, C.next y = .inr c := by
  rw [iter_succ, h] at h'
  simp only [Option.bind_some] at h'
  rcases e : C.next y with a | b
  · rw [e] at h'; simp at h'
  · exact ⟨b, rfl⟩

lemma iter_none_succ {z : Xing F} {m : ℕ} (h : C.iter z m = none) : C.iter z (m + 1) = none := by
  rw [iter_succ, h]; rfl

lemma iter_none_mono {z : Xing F} {m m' : ℕ} (h : C.iter z m = none) (hm : m ≤ m') :
    C.iter z m' = none := by
  induction hm with
  | refl => exact h
  | step _ ih => exact iter_none_succ ih

lemma Valid.iter_isX (hC : C.Valid) {z y : Xing F} (hz : C.IsX z) {m : ℕ}
    (h : C.iter z m = some y) : C.IsX y := by
  induction m generalizing y with
  | zero => simp only [iter, Option.some.injEq] at h; rw [← h]; exact hz
  | succ m ih =>
    obtain ⟨x, hx, hn⟩ := iter_succ_some h
    exact hC.next_isX (ih hx) hn

/-- Walking back: equal crossings at steps `i + 1` and `j + 1` come from equal crossings. -/
lemma Valid.iter_back (hC : C.Valid) {z a : Xing F} (hz : C.IsX z) {i j : ℕ}
    (hi : C.iter z (i + 1) = some a) (hj : C.iter z (j + 1) = some a) : C.iter z i = C.iter z j := by
  obtain ⟨y, hy, hny⟩ := iter_succ_some hi
  obtain ⟨y', hy', hny'⟩ := iter_succ_some hj
  rw [hy, hy', hC.next_inj (hC.iter_isX hz hy) (hC.iter_isX hz hy') hny hny']

/-- **The first repetition is the start.** If the walk visits a crossing twice, at steps
`i < j`, then it is back at `z` after `j − i` steps. -/
lemma Valid.iter_periodic (hC : C.Valid) {z a : Xing F} (hz : C.IsX z) :
    ∀ {i j : ℕ}, i < j → C.iter z i = some a → C.iter z j = some a → C.iter z (j - i) = some z := by
  intro i
  induction i generalizing a with
  | zero =>
    intro j _ hi hj
    simp only [iter, Option.some.injEq] at hi
    rw [Nat.sub_zero, hj, hi]
  | succ i ih =>
    intro j hij hi hj
    obtain ⟨j', rfl⟩ : ∃ j', j = j' + 1 := ⟨j - 1, by omega⟩
    have hb := hC.iter_back hz hi hj
    obtain ⟨y, hy, -⟩ := iter_succ_some hi
    rw [hy] at hb
    have := ih (a := y) (by omega) hy hb.symm
    rwa [show j' + 1 - (i + 1) = j' - i by omega]

/-- **No crossing in both directions** (the palindrome argument): the walk never meets the
reversal of a crossing it has visited. -/
lemma Valid.iter_ne_rev (hC : C.Valid) {z : Xing F} (hz : C.IsX z) :
    ∀ d i (a : Xing F), C.iter z i = some a → C.iter z (i + d + 1) ≠ some (C.rev a) := by
  intro d
  induction d using Nat.strong_induction_on with
  | _ d ih =>
    intro i a ha hj
    have hax := hC.iter_isX hz ha
    rcases d with _ | _ | d
    · obtain ⟨y, hy, hn⟩ := iter_succ_some hj
      rw [ha] at hy
      cases hy
      exact hC.next_ne_rev hax hn
    · obtain ⟨c, hc, hnc⟩ := iter_succ_some hj
      obtain ⟨y, hy, hny⟩ := iter_succ_some hc
      rw [ha] at hy
      cases hy
      have hcx := hC.next_isX hax hny
      have e := hC.next_rev hcx hnc
      rw [hC.rev_rev hax] at e
      have : c = C.rev c := Sum.inl.inj (hny.symm.trans e)
      exact hC.rev_ne hcx this.symm
    · -- the steps after `a` and before `rev a`
      obtain ⟨c, hc, hnc⟩ := iter_succ_some hj
      have hb : ∃ b, C.iter z (i + 1) = some b ∧ C.next a = .inl b := by
        have h1 : C.iter z (i + 1) ≠ none := by
          intro h
          have := iter_none_mono h (show i + 1 ≤ i + (d + 2) + 1 by omega)
          rw [this] at hj; cases hj
        obtain ⟨b, hb⟩ := Option.ne_none_iff_exists'.mp h1
        obtain ⟨y, hy, hny⟩ := iter_succ_some hb
        rw [ha] at hy; cases hy
        exact ⟨b, hb, hny⟩
      obtain ⟨b, hb, hnb⟩ := hb
      have hcx := hC.iter_isX hz hc
      have e := hC.next_rev hcx hnc
      rw [hC.rev_rev hax] at e
      have hbc : b = C.rev c := Sum.inl.inj (hnb.symm.trans e)
      have hc' : C.iter z (i + 1 + d + 1) = some (C.rev b) := by
        rw [hbc, hC.rev_rev hcx, ← hc]; congr 1; omega
      exact ih d (by omega) (i + 1) b hb hc'

/-- **`ic-trace-terminates`: `TraceFrom` ends at a corner or comes back to its start, without
repeating a crossing.** From any crossing `z` of a valid complex model, the walk of Supp. Alg. 1
(the Rust `trace_with`; each step is `traceStep`, `Valid.next_eq_oracle`) visits `m ≥ 1`
crossings, pairwise different and none of them in both directions, so `2m ≤ Σ_sides n⁺`
(`m ≤ Σ_edges n⁺`), and then either ends at a corner `c` (the `m`-th crossing steps to a corner)
or is back at `z`, in the same direction (a closed curve: the only way the Rust tracer can see a
crossing twice). [formula: ic-trace-terminates] -/
theorem Valid.trace_terminates [Fintype F] [DecidableEq F] (hC : C.Valid) {z : Xing F}
    (hz : C.IsX z) :
    (∃ m y c, C.Visits z m ∧ C.iter z (m - 1) = some y ∧ C.next y = .inr c) ∨
      (∃ m, C.Visits z m ∧ C.iter z m = some z) := by
  classical
  -- some step ends the walk or returns to `z`
  have hex : ∃ m, 0 < m ∧ (C.iter z m = none ∨ C.iter z m = some z) := by
    by_cases hall : ∀ i ≤ C.xings.card, C.iter z i ≠ none
    · have hmaps : Set.MapsTo (fun i => (C.iter z i).getD z) ↑(Finset.range (C.xings.card + 1))
          ↑C.xings := by
        intro i hi
        simp only [Finset.coe_range, Set.mem_Iio] at hi
        obtain ⟨y, hy⟩ := Option.ne_none_iff_exists'.mp (hall i (by omega))
        simp only [Finset.mem_coe, hy, Option.getD_some]
        exact mem_xings.2 (hC.iter_isX hz hy)
      obtain ⟨i, hi, j, hj, hij, he⟩ :=
        Finset.exists_ne_map_eq_of_card_lt_of_maps_to (by simp) hmaps
      simp only [Finset.mem_range] at hi hj
      obtain ⟨y, hy⟩ := Option.ne_none_iff_exists'.mp (hall i (by omega))
      obtain ⟨y', hy'⟩ := Option.ne_none_iff_exists'.mp (hall j (by omega))
      simp only [hy, hy', Option.getD_some] at he
      subst he
      rcases Nat.lt_or_gt_of_ne hij with h | h
      · exact ⟨j - i, by omega, Or.inr (hC.iter_periodic hz h hy hy')⟩
      · exact ⟨i - j, by omega, Or.inr (hC.iter_periodic hz h hy' hy)⟩
    · simp only [not_forall, not_not] at hall
      obtain ⟨i, -, hi⟩ := hall
      refine ⟨i, ?_, Or.inl hi⟩
      rcases Nat.eq_zero_or_pos i with rfl | h
      · simp [iter] at hi
      · exact h
  let m := Nat.find hex
  have hm : 0 < m ∧ (C.iter z m = none ∨ C.iter z m = some z) := Nat.find_spec hex
  have hmin : ∀ k, 0 < k → k < m → C.iter z k ≠ none ∧ C.iter z k ≠ some z := by
    intro k hk hkm
    have := Nat.find_min hex hkm
    simp only [not_and, not_or] at this
    exact this hk
  -- the visited crossings
  have hdef : ∀ i < m, ∃ y, C.iter z i = some y ∧ C.IsX y := by
    intro i hi
    rcases Nat.eq_zero_or_pos i with rfl | hpos
    · exact ⟨z, rfl, hz⟩
    · obtain ⟨y, hy⟩ := Option.ne_none_iff_exists'.mp (hmin i hpos hi).1
      exact ⟨y, hy, hC.iter_isX hz hy⟩
  have hdist : ∀ i j y y', i < m → j < m → C.iter z i = some y → C.iter z j = some y' →
      (y = y' → i = j) ∧ y ≠ C.rev y' := by
    intro i j y y' hi hj hy hy'
    refine ⟨fun e => ?_, fun e => ?_⟩
    · subst e
      by_contra hne
      rcases Nat.lt_or_gt_of_ne hne with h | h
      · exact (hmin (j - i) (by omega) (by omega)).2 (hC.iter_periodic hz h hy hy')
      · exact (hmin (i - j) (by omega) (by omega)).2 (hC.iter_periodic hz h hy' hy)
    · have hy'x := hC.iter_isX hz hy'
      rcases lt_trichotomy i j with h | rfl | h
      · have := hC.iter_ne_rev hz (j - i - 1) i y hy
        rw [show i + (j - i - 1) + 1 = j by omega, hy', e, hC.rev_rev hy'x] at this
        exact this rfl
      · rw [hy] at hy'; cases hy'
        exact hC.rev_ne hy'x e.symm
      · have := hC.iter_ne_rev hz (i - j - 1) j y' hy'
        rw [show j + (i - j - 1) + 1 = i by omega, hy, e] at this
        exact this rfl
  have hbound : 2 * m ≤ ∑ s, C.np s := by
    rw [← card_xings]
    let g : ℕ → Xing F := fun i => (C.iter z i).getD z
    let S := (Finset.range m).image g
    have hg : ∀ i < m, C.iter z i = some (g i) := by
      intro i hi
      obtain ⟨y, hy, -⟩ := hdef i hi
      simp only [g, hy, Option.getD_some]
    have hS : S.card = m := by
      rw [Finset.card_image_of_injOn, Finset.card_range]
      intro i hi j hj e
      simp only [Finset.coe_range, Set.mem_Iio] at hi hj
      exact (hdist i j (g i) (g j) hi hj (hg i hi) (hg j hj)).1 e
    have hSx : ∀ y ∈ S, C.IsX y := by
      intro y hy
      obtain ⟨i, hi, rfl⟩ := Finset.mem_image.mp hy
      exact hC.iter_isX hz (hg i (Finset.mem_range.mp hi))
    have hR : (S.image C.rev).card = m := by
      rw [Finset.card_image_of_injOn, hS]
      intro a ha b hb e
      have := congrArg C.rev e
      rwa [hC.rev_rev (hSx a ha), hC.rev_rev (hSx b hb)] at this
    have hdisj : Disjoint S (S.image C.rev) := by
      rw [Finset.disjoint_left]
      intro a ha hb'
      obtain ⟨b, hb, e⟩ := Finset.mem_image.mp hb'
      obtain ⟨i, hi, rfl⟩ := Finset.mem_image.mp ha
      obtain ⟨j, hj, rfl⟩ := Finset.mem_image.mp hb
      simp only [Finset.mem_range] at hi hj
      exact (hdist i j (g i) (g j) hi hj (hg i hi) (hg j hj)).2 e.symm
    have hsub : S ∪ S.image C.rev ⊆ C.xings := by
      intro y hy
      rcases Finset.mem_union.mp hy with h | h
      · exact mem_xings.2 (hSx y h)
      · obtain ⟨b, hb, rfl⟩ := Finset.mem_image.mp h
        exact mem_xings.2 (hC.rev_isX (hSx b hb))
    have := Finset.card_le_card hsub
    rw [Finset.card_union_of_disjoint hdisj, hS, hR] at this
    omega
  have hvis : C.Visits z m := ⟨hm.1, hdef, hdist, hbound⟩
  rcases hm.2 with h | h
  · left
    obtain ⟨y, hy, -⟩ := hdef (m - 1) (by omega)
    obtain ⟨c, hc⟩ := iter_succ_none hy (by rwa [show m - 1 + 1 = m by omega])
    exact ⟨m, y, c, hvis, hy, hc⟩
  · exact Or.inr ⟨m, hvis, h⟩

/-- A walk that comes back to its start was entered through a crossing. -/
lemma Valid.not_periodic_of_end (hC : C.Valid) {z : Xing F} (hz : C.IsX z) (he : C.IsEnd z)
    {m : ℕ} (hm : 0 < m) : C.iter z m ≠ some z := by
  intro h
  obtain ⟨y, hy, hn⟩ := iter_succ_some (show C.iter z (m - 1 + 1) = some z by
    rwa [show m - 1 + 1 = m by omega])
  exact hC.not_next_eq_end (hC.iter_isX hz hy) he hn

/-- **From the first crossing of a curve leaving a vertex, `TraceFrom` ends at a corner**, after
at most `Σ_edges n⁺` crossings, none visited twice or in both directions.
[formula: ic-trace-terminates] -/
theorem Valid.trace_from_end [Fintype F] [DecidableEq F] (hC : C.Valid) {z : Xing F}
    (hz : C.IsX z) (he : C.IsEnd z) :
    ∃ m y c, C.Visits z m ∧ C.iter z (m - 1) = some y ∧ C.next y = .inr c := by
  rcases hC.trace_terminates hz with h | ⟨m, hv, h⟩
  · exact h
  · exact absurd h (hC.not_periodic_of_end hz he hv.pos)

/-- **Without closed curves (SPEC §2.4 item 4), `TraceFrom` from every crossing ends at a
corner**, after at most `Σ_edges n⁺` crossings, none visited twice or in both directions.
[formula: ic-trace-terminates] -/
theorem Valid.trace_terminates_of_noClosed [Fintype F] [DecidableEq F] (hC : C.Valid)
    (hnc : C.NoClosed) {z : Xing F} (hz : C.IsX z) :
    ∃ m y c, C.Visits z m ∧ C.iter z (m - 1) = some y ∧ C.next y = .inr c := by
  rcases hC.trace_terminates hz with h | ⟨m, hv, h⟩
  · exact h
  · exact absurd h (hnc z hz m hv.pos)

end Cx

end DDG.IntegerCoords
