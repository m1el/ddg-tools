import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Data.Finset.Card

/-!
# Integer coordinates: removing an inserted vertex by flips (tier B′, Thm. D.1)

`crates/ddg-integer-coords/SPEC.md` §4.6, I8, §9.2; the paper Sec. 3.5, Supp. App. D, Thm. D.1,
Fig. 23.

**The star.** An inserted vertex `i` of degree `d` with neighbours `j_0, …, j_{d−1}` in
counter-clockwise order (indices mod `d`) and triangles `t_k = (i, j_k, j_{k+1})`. Triangle `t_k` has
the angle `θ_k` at `i`, `α_k` at `j_k`, `β_k` at `j_{k+1}`; all positive, `θ_k + α_k + β_k = π`, and
`Σ θ_k = 2π` (`i` is flat: an inserted vertex has no curvature). Flipping the edge `i j_k` (shared by
`t_{k−1}` and `t_k`) needs the quadrilateral `i, j_{k−1}, j_k, j_{k+1}` to be convex: its angles are
the *inner* corner `θ_{k−1} + θ_k` at `i`, the *outer* corner `β_{k−1} + α_k` at `j_k`, and `α_{k−1}`,
`β_k` (each `< π`, being triangle angles). In a simplicial complex the flip removes the edge `i j_k`
and adds `j_{k−1} j_{k+1}`, which is not incident to `i`: the degree of `i` drops by one, and every
vertex keeps degree at least two (the neighbours of `i` are distinct). So everything Thm. D.1 needs is
the angle statement, which is what is proved here, exactly as the paper's proof goes.

* **`ic-vertex-removal`** (`exists_flippable`) — for `d ≥ 4` there is an edge `i j_k` whose inner
  and outer corners are both at most `π`. Steps: the outer corners sum to `(d − 2)π`, so at least three
  are at most `π` (`three_le_card_outer`); two inner corners at spokes that are not adjacent cover
  disjoint angles at `i`, so at most two inner corners exceed `π` (`card_inner_le_two`).
* The paper's "can be flipped" means a convex quadrilateral with angles **at most** `π`. With *strict*
  convexity (the flip condition of IT-SPEC N5, `TooFlat`/convexity tolerance) the statement is false:
  the degree-4 star with right angles at `i` and half-right angles elsewhere has every inner corner
  equal to `π` (`cross_star_no_strict`), as the paper itself remarks ("the equality case … degree four
  cross configuration"). If no inner and no outer corner equals `π`, a strictly convex flip exists
  (`exists_flippable_strict`).
* **The Δ-complex case** (`ic-vertex-removal-delta`, open; `CITED.md`). `FlatStar` is the cyclic
  sequence of the `d` triangle corners at `i`; it does not assume that the neighbours `j_k` or the
  triangles `t_k` are distinct (a triangle with two corners at `i` appears twice, with its three
  angles permuted). The proof above uses only the corner cycle, so its conclusion holds in a
  Δ-complex as well: some spoke `i j_k` has a convex quadrilateral. What is not proved, and what the
  paper leaves as a conjecture, is the combinatorial half: in a Δ-complex that spoke can be a loop
  (both ends at `i`), `j_{k−1}` or `j_{k+1}` can be `i` itself (the new edge is then a loop at `i`),
  and `j_k` can have degree one, so the flip may be impossible or may not lower the degree. A proof
  needs a progress measure beyond the degree, as Sharp and Crane 2020a (FlipOut), App. A, Thm. A.5–A.8,
  use for a wedge of angle `< π` (the first twice-incident edge moves towards the start); for the full
  `2π` star the inner corners are not automatically below `π`, which is where that argument stops.
-/

namespace DdgProofs

namespace IntegerCoords

open Finset Real

variable {n : ℕ}

/-- A flat star of degree `n + 1` around an inserted vertex: the three angles of each triangle. -/
structure FlatStar (n : ℕ) where
  /-- angle at the centre `i` of triangle `t_k = (i, j_k, j_{k+1})` -/
  θ : Fin (n + 1) → ℝ
  /-- angle at `j_k` -/
  α : Fin (n + 1) → ℝ
  /-- angle at `j_{k+1}` -/
  β : Fin (n + 1) → ℝ
  θ_pos : ∀ k, 0 < θ k
  α_pos : ∀ k, 0 < α k
  β_pos : ∀ k, 0 < β k
  triangle : ∀ k, θ k + α k + β k = π
  flat : ∑ k, θ k = 2 * π

namespace FlatStar

variable (S : FlatStar n)

/-- The angle at `i` of the quadrilateral of the edge `i j_k`. -/
def inner (k : Fin (n + 1)) : ℝ := S.θ (k - 1) + S.θ k

/-- The angle at `j_k` of the quadrilateral of the edge `i j_k`. -/
def outer (k : Fin (n + 1)) : ℝ := S.β (k - 1) + S.α k

/-- The outer corners sum to `(d − 2)π`. -/
lemma sum_outer : ∑ k, S.outer k = ((n + 1 : ℕ) - 2 : ℝ) * π := by
  have h1 : ∑ k, S.β (k - 1) = ∑ k, S.β k :=
    Fintype.sum_equiv (Equiv.subRight 1) _ _ (fun _ => rfl)
  have h2 : ∑ k : Fin (n + 1), (S.θ k + S.α k + S.β k) = ∑ _k : Fin (n + 1), π :=
    Finset.sum_congr rfl (fun k _ => S.triangle k)
  simp only [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul] at h2
  simp only [outer, Finset.sum_add_distrib, h1]
  rw [S.flat] at h2
  push_cast at h2 ⊢
  linarith

/-- **At least three outer corners are at most `π`** (for `d ≥ 3`). -/
lemma three_le_card_outer (hn : 2 ≤ n) :
    3 ≤ (univ.filter (fun k => S.outer k ≤ π)).card := by
  by_contra hlt
  push Not at hlt
  set O := univ.filter (fun k => S.outer k ≤ π) with hO
  have hc : (univ.filter (fun k => ¬ S.outer k ≤ π)).card = n + 1 - O.card := by
    have := Finset.card_filter_add_card_filter_not (s := (univ : Finset (Fin (n + 1))))
      (fun k => S.outer k ≤ π)
    rw [Finset.card_univ, Fintype.card_fin] at this
    rw [hO]; omega
  have hpos : ∀ k, 0 < S.outer k := fun k => add_pos (S.β_pos _) (S.α_pos _)
  -- the sum over the complement exceeds `(#complement) π`
  have hbig : ((n + 1 - O.card : ℕ) : ℝ) * π <
      ∑ k ∈ univ.filter (fun k => ¬ S.outer k ≤ π), S.outer k := by
    rw [← hc]
    have hne : (univ.filter (fun k => ¬ S.outer k ≤ π)).Nonempty := by
      rw [← Finset.card_pos, hc]; omega
    have := Finset.sum_lt_sum_of_nonempty hne (f := fun _ => π) (g := S.outer)
      (fun k hk => by simpa using (Finset.mem_filter.mp hk).2)
    simpa [Finset.sum_const, nsmul_eq_mul] using this
  have hle : ∑ k ∈ univ.filter (fun k => ¬ S.outer k ≤ π), S.outer k ≤ ∑ k, S.outer k :=
    Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
      (fun k _ _ => (hpos k).le)
  rw [S.sum_outer] at hle
  have hcast : ((n + 1 - O.card : ℕ) : ℝ) ≥ ((n + 1 : ℕ) - 2 : ℝ) := by
    have : n + 1 - O.card ≥ n + 1 - 2 := by omega
    have h2 : ((n + 1 - 2 : ℕ) : ℝ) = ((n + 1 : ℕ) - 2 : ℝ) := by
      rw [Nat.cast_sub (by omega)]; push_cast; ring
    rw [← h2]; exact_mod_cast this
  nlinarith [Real.pi_pos]

/-- Two spokes `a`, `b` are adjacent around `i`. -/
def Adj (a b : Fin (n + 1)) : Prop := b = a + 1 ∨ a = b + 1

lemma val_add_one' (a : Fin (n + 1)) :
    (a + 1).val = a.val + 1 ∨ (a.val = n ∧ (a + 1).val = 0) := by
  rw [Fin.val_add_one]
  split_ifs with h
  · right; exact ⟨by rw [h]; simp, rfl⟩
  · left; rfl

/-- Among three distinct spokes of a vertex of degree at least four, two are not adjacent. -/
lemma not_all_adj (hn : 3 ≤ n) {a b c : Fin (n + 1)} (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (h1 : Adj a b) (h2 : Adj a c) (h3 : Adj b c) : False := by
  have va := a.isLt; have vb := b.isLt; have vc := c.isLt
  have dab : a.val ≠ b.val := fun h => hab (Fin.ext h)
  have dac : a.val ≠ c.val := fun h => hac (Fin.ext h)
  have dbc : b.val ≠ c.val := fun h => hbc (Fin.ext h)
  have ea := val_add_one' a; have eb := val_add_one' b; have ec := val_add_one' c
  unfold Adj at h1 h2 h3
  rcases h1 with h1 | h1 <;> rcases h2 with h2 | h2 <;> rcases h3 with h3 | h3 <;>
    have g1 := congrArg Fin.val h1 <;> have g2 := congrArg Fin.val h2 <;>
    have g3 := congrArg Fin.val h3 <;> omega

/-- The inner corners of two non-adjacent distinct spokes cover disjoint angles at `i`. -/
lemma inner_add_inner_le (hn : 1 ≤ n) (a b : Fin (n + 1)) (hab : a ≠ b) (h1 : ¬ Adj a b) :
    S.inner a + S.inner b ≤ 2 * π := by
  unfold Adj at h1
  push Not at h1
  have d1 : a - 1 ≠ a := by
    intro h; have := congrArg (· + 1) h; simp at this; omega
  have d2 : b - 1 ≠ b := by
    intro h; have := congrArg (· + 1) h; simp at this; omega
  have d3 : a - 1 ≠ b := fun h => h1.2 (by rw [← h]; simp)
  have d4 : a - 1 ≠ b - 1 := fun h => hab (by simpa using h)
  have d5 : a ≠ b - 1 := fun h => h1.1 (by rw [h]; simp)
  have hsub : ({a - 1, a, b - 1, b} : Finset (Fin (n + 1))) ⊆ univ := subset_univ _
  have hle := Finset.sum_le_sum_of_subset_of_nonneg hsub
    (f := S.θ) (fun k _ _ => (S.θ_pos k).le)
  rw [S.flat] at hle
  rw [Finset.sum_insert (by simp [d1, d3, d4]), Finset.sum_insert (by simp [hab, d5]),
    Finset.sum_insert (by simp [d2]), Finset.sum_singleton] at hle
  unfold inner
  linarith

/-- **At most two inner corners exceed `π`** (for `d ≥ 4`). -/
lemma card_inner_le_two (hn : 3 ≤ n) : (univ.filter (fun k => π < S.inner k)).card ≤ 2 := by
  by_contra hlt
  push Not at hlt
  obtain ⟨a, ha, b, hb, c, hc, hab, hac, hbc⟩ := Finset.two_lt_card.mp hlt
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha hb hc
  by_cases h1 : Adj a b
  · by_cases h2 : Adj a c
    · by_cases h3 : Adj b c
      · exact not_all_adj hn hab hac hbc h1 h2 h3
      · have := S.inner_add_inner_le (by omega) b c hbc h3; linarith
    · have := S.inner_add_inner_le (by omega) a c hac h2; linarith
  · have := S.inner_add_inner_le (by omega) a b hab h1; linarith

/-- **`ic-vertex-removal` (Thm. D.1).** A flat star (angle sum `2π`) of degree `d ≥ 4` in a simplicial
complex has an edge `i j_k` whose quadrilateral is convex: its inner corner `θ_{k−1} + θ_k` and its
outer corner `β_{k−1} + α_k` are at most `π` (the other two corners are triangle angles, below `π`).
Flipping it lowers the degree of `i` by one (see the module docstring for the combinatorics).
[formula: ic-vertex-removal] -/
theorem exists_flippable (hn : 3 ≤ n) : ∃ k, S.inner k ≤ π ∧ S.outer k ≤ π := by
  have hO := S.three_le_card_outer (by omega)
  have hI := S.card_inner_le_two hn
  obtain ⟨k, hk, hk'⟩ := Finset.exists_mem_notMem_of_card_lt_card (by omega :
    (univ.filter (fun k => π < S.inner k)).card < (univ.filter (fun k => S.outer k ≤ π)).card)
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, not_lt] at hk hk'
  exact ⟨k, hk', hk⟩

/-- Without corners of exactly `π`, the quadrilateral is strictly convex.
[formula: ic-vertex-removal] -/
theorem exists_flippable_strict (hn : 3 ≤ n) (hi : ∀ k, S.inner k ≠ π) (ho : ∀ k, S.outer k ≠ π) :
    ∃ k, S.inner k < π ∧ S.outer k < π := by
  obtain ⟨k, h1, h2⟩ := S.exists_flippable hn
  exact ⟨k, lt_of_le_of_ne h1 (hi k), lt_of_le_of_ne h2 (ho k)⟩

end FlatStar

/-- The degree-4 cross: right angles at the centre, half-right angles at the neighbours. -/
noncomputable def crossStar : FlatStar 3 where
  θ := fun _ => π / 2
  α := fun _ => π / 4
  β := fun _ => π / 4
  θ_pos := fun _ => by positivity
  α_pos := fun _ => by positivity
  β_pos := fun _ => by positivity
  triangle := fun _ => by ring
  flat := by simp [Finset.sum_const]; ring

/-- **The strict version of Thm. D.1 is false**: in the degree-4 cross every inner corner is exactly
`π`, so no edge at the centre has a strictly convex quadrilateral (the paper's "equality case").
[formula: ic-vertex-removal] -/
theorem cross_star_no_strict : ¬ ∃ k, crossStar.inner k < π ∧ crossStar.outer k < π := by
  rintro ⟨k, h, -⟩
  simp only [FlatStar.inner, crossStar] at h
  linarith

end IntegerCoords

end DdgProofs
