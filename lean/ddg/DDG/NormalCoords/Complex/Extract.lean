import DDG.NormalCoords.Complex.Trace
import DDG.NormalCoords.Complex.Labels

/-!
# Integer coordinates, tier C: the input edges crossed by an intrinsic edge
(`ic-extract-intrinsic-edge`, the combinatorial part)

`crates/ddg-integer-coords/SPEC.md` §3.1 (`ExtractIntrinsicEdge`); the Rust
`extract_intrinsic_edge` (`crates/ddg-integer-coords/src/geometry.rs`), which identifies the input
edge through each crossing `p` of the edge (Alg. 1 backwards to the end of the curve, then Alg. 3),
stopping early at a crossing of the edge already identified.

* `Valid.end_of_crossing`: in a state without closed curves (SPEC §2.4 item 4), every crossing `z`
  lies on the curve of exactly one end `e` (an arc leaving a corner): the walk from `e` passes
  through `z` after `k` steps, for exactly one `(e, k)` (`Valid.end_unique`), and ends at a corner
  (`trace_from_end`). The end is found by walking backwards from `z` (Alg. 2, the backward half of
  `extract_curve`): the last crossing `y` before the corner gives `e = rev y`.
* `Valid.end_is_arc`: an end is the `q`-th arc leaving its corner `c` (`e = (side c + 1,
  (emList c)[q])`), so Alg. 3 gives it the label `lab c (n⁻_c + q)` (`Valid.identify_arc`).

So, for an edge with `n⁺` crossings, `ExtractIntrinsicEdge` returns, for `p = 0, …, n⁺ − 1` (in
the order of the crossings along the halfedge from its tail), the input edge of the unique curve
through crossing `p`; two crossings on the same curve get the same input edge (the Rust shortcut is
exact). The positions along the edge are the strip geometry (`ic-strip-intersection`).
-/

namespace DDG.NormalCoords

namespace Cx

variable {F : Type} {C : Cx F}

/-- The last crossing before a corner, reversed, is an end: the arc of its face through it comes
from that corner. -/
lemma isEnd_rev_of_next_inr {y : Xing F} {c : F × Fin 3} (h : C.next y = .inr c) :
    C.IsEnd (C.rev y) := by
  obtain ⟨e, -⟩ := next_eq_inr h
  exact ⟨c.2, e⟩

/-- Walking backwards: if the walk from `z` reaches `w` after `k` steps, the walk from `rev w`
reaches `rev z` after `k` steps. -/
lemma Valid.iter_rev (hC : C.Valid) {z : Xing F} (hz : C.IsX z) :
    ∀ {k : ℕ} {w : Xing F}, C.iter z k = some w → C.iter (C.rev w) k = some (C.rev z) := by
  intro k
  induction k with
  | zero =>
    intro w h
    simp only [iter, Option.some.injEq] at h
    rw [← h]; rfl
  | succ k ih =>
    intro w h
    obtain ⟨y, hy, hn⟩ := iter_succ_some h
    have hyx := hC.iter_isX hz hy
    have e := hC.next_rev hyx hn
    -- the walk from `rev w` first goes to `rev y`, then `k` more steps to `rev z`
    have ih' := ih hy
    have hwx : C.IsX w := hC.next_isX hyx hn
    -- iter (rev w) (k+1) = iter (rev y) k, by the first step
    have step : ∀ (j : ℕ) (u v : Xing F), C.next u = .inl v →
        C.iter u (j + 1) = C.iter v j := by
      intro j
      induction j with
      | zero => intro u v huv; simp [iter, huv]
      | succ j ihj =>
        intro u v huv
        rw [iter_succ, ihj u v huv, ← iter_succ]
    rw [step k _ _ e, ih']

/-- Two ends whose walks meet at the same crossing are the same end, after the same number of
steps. -/
lemma Valid.end_unique (hC : C.Valid) {e e' z : Xing F} (he : C.IsX e) (he' : C.IsX e')
    (hee : C.IsEnd e) (hee' : C.IsEnd e') :
    ∀ {k k' : ℕ}, C.iter e k = some z → C.iter e' k' = some z → e = e' ∧ k = k' := by
  intro k
  induction k generalizing z with
  | zero =>
    intro k' h h'
    simp only [iter, Option.some.injEq] at h
    subst h
    rcases k' with _ | k'
    · simp only [iter, Option.some.injEq] at h'; exact ⟨h'.symm, rfl⟩
    · obtain ⟨y, hy, hn⟩ := iter_succ_some h'
      exact absurd hn (hC.not_next_eq_end (hC.iter_isX he' hy) hee)
  | succ k ih =>
    intro k' h h'
    obtain ⟨y, hy, hn⟩ := iter_succ_some h
    rcases k' with _ | k'
    · simp only [iter, Option.some.injEq] at h'
      subst h'
      exact absurd hn (hC.not_next_eq_end (hC.iter_isX he hy) hee')
    · obtain ⟨y', hy', hn'⟩ := iter_succ_some h'
      have := hC.next_inj (hC.iter_isX he hy) (hC.iter_isX he' hy') hn hn'
      subst this
      obtain ⟨h1, h2⟩ := ih hy hy'
      exact ⟨h1, by omega⟩

/-- **`ic-extract-intrinsic-edge` (combinatorial part): every crossing lies on the curve of exactly
one end.** In a valid state without closed curves, for every crossing `z` there is an end `e` (an
arc leaving a corner) and a `k` with the walk from `e` at `z` after `k` steps, and the walk from
`e` ends at a corner; `(e, k)` is unique, and `e` is the reversal of the last crossing of the walk
backwards from `z` (Alg. 2). [formula: ic-extract-intrinsic-edge] -/
theorem Valid.end_of_crossing [Fintype F] [DecidableEq F] (hC : C.Valid) (hnc : C.NoClosed)
    {z : Xing F} (hz : C.IsX z) :
    ∃ e k, C.IsX e ∧ C.IsEnd e ∧ C.iter e k = some z ∧
      (∃ m y c, C.Visits e m ∧ C.iter e (m - 1) = some y ∧ C.next y = .inr c) ∧
      ∀ e' k', C.IsX e' → C.IsEnd e' → C.iter e' k' = some z → e' = e ∧ k' = k := by
  have hrz := hC.rev_isX hz
  obtain ⟨m, y, c, hv, hy, hn⟩ := hC.trace_terminates_of_noClosed hnc hrz
  have hyx := hC.iter_isX hrz hy
  have hee := isEnd_rev_of_next_inr hn
  have hex := hC.rev_isX hyx
  have hit := hC.iter_rev hrz hy
  rw [hC.rev_rev hz] at hit
  refine ⟨C.rev y, m - 1, hex, hee, hit, hC.trace_from_end hex hee, ?_⟩
  intro e' k' he' hee' h'
  exact hC.end_unique he' hex hee' hee h' hit

/-- **An end is an arc leaving its corner**: in a valid state, an end `e` on side `(f, k)` comes
from corner `c = (f, k + 2)` and is the `q`-th arc of `emList c`, so Alg. 3 identifies it by
`identifyIndex` = `lab c (n⁻_c + q)` (`Valid.identify_arc`). [formula: ic-extract-intrinsic-edge] -/
theorem Valid.end_is_arc (hC : C.Valid) {e : Xing F} (he : C.IsX e) (hee : C.IsEnd e) :
    ∃ q, ∃ hq : q < C.emc (e.1.1, e.1.2 + 2),
      (C.emList (e.1.1, e.1.2 + 2))[q]'hq = e.2 ∧ e.1.2 = e.1.2 + 2 + 1 := by
  obtain ⟨c, hc⟩ := hee
  have hp : e.2 < (C.M e.1.1).n e.1.2 := by rw [M_n]; exact he
  -- the corner of an arc to side `k` is `k + 2`
  have hck : c = e.1.2 + 2 := by
    rcases (hC.M_valid e.1.1).normal e.1.2 e.2 hp with ⟨q, -, h⟩ | ⟨q, -, h⟩ | h
    · rw [hc] at h; cases h
    · rw [hc] at h; cases h
    · rw [hc] at h; exact Pt.corner.inj h
  subst hck
  have hmem : e.2 ∈ C.emList (e.1.1, e.1.2 + 2) := by
    simp only [emList, List.mem_filter, List.mem_range, fin3_add_two_add_one, decide_eq_true_eq]
    exact ⟨he, hc⟩
  obtain ⟨q, hq, hqe⟩ := List.getElem_of_mem hmem
  exact ⟨q, hq, hqe, (fin3_add_two_add_one _).symm⟩

end Cx

end DDG.NormalCoords
