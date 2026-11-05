import Mathlib.Data.Real.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Tactic.Positivity
import Mathlib.Data.Finset.Max
import Mathlib.Order.Interval.Finset.Fin
import Mathlib.Tactic.Linarith

/-!
# Delaunay layer 4: the counting lemma for chords through a fan

A finite sequence of chords `i = 0, …, n − 1`, each with a *pivot* `v i` (a point) and a
*length* `lam i ≥ 0`. A maximal run of consecutive chords with the same pivot is a *block*; two
consecutive chords with different pivots form a *switch*. If

* every run of consecutive chords with the same pivot has at most `M` chords (`hrun`), and
* every switch costs at least `c₀ > 0` of length, `c₀ ≤ lam i + lam (i+1)` (`hsw`),

then `n ≤ M · (1 + 2 Λ / c₀)` with `Λ = Σ lam i` (`fan_count`): there are at most `2Λ/c₀`
switches, since each chord belongs to at most two switches, hence at most `1 + 2Λ/c₀` blocks of
at most `M` chords each.

This is the combinatorial half of the bound on the number of edges of `s` crossed by an edge
of a flipped triangulation (`Fan.lean`); the geometric half supplies `M` (a straight segment
through a fan of corners at one vertex subtends less than `π`, so a block has fewer than
`π / α_min` chords) and `c₀` (`sin α_min · ℓ_min`).
-/

namespace DDG.Delaunay.Surface

open Finset

/-- The switches: the indices `i` with `i + 1 < n` and `v i ≠ v (i+1)`. -/
noncomputable def switches {α : Type*} {n : ℕ} (v : Fin n → α) : Finset (Fin n) :=
  open scoped Classical in univ.filter fun i => ∃ h : i.1 + 1 < n, v i ≠ v ⟨i.1 + 1, h⟩

/-- Membership in `switches`. -/
theorem mem_switches {α : Type*} {n : ℕ} (v : Fin n → α) (i : Fin n) :
    i ∈ switches v ↔ ∃ h : i.1 + 1 < n, v i ≠ v ⟨i.1 + 1, h⟩ := by
  classical
  unfold switches
  rw [mem_filter]
  simp only [mem_univ, true_and]

/-- At most `2Λ/c₀` switches. -/
theorem card_switches_le {α : Type*} {n : ℕ} (v : Fin n → α) (lam : Fin n → ℝ) {c₀ : ℝ}
    (hc : 0 < c₀) (hlam : ∀ i, 0 ≤ lam i)
    (hsw : ∀ i : Fin n, ∀ h : i.1 + 1 < n, v i ≠ v ⟨i.1 + 1, h⟩ → c₀ ≤ lam i + lam ⟨i.1 + 1, h⟩) :
    ((switches v).card : ℝ) ≤ 2 * (∑ i, lam i) / c₀ := by
  classical
  set S := switches v
  let succ' : Fin n → Fin n := fun i => if h : i.1 + 1 < n then ⟨i.1 + 1, h⟩ else i
  have hsucc : ∀ i ∈ S, ∃ h : i.1 + 1 < n, succ' i = ⟨i.1 + 1, h⟩ ∧ v i ≠ v ⟨i.1 + 1, h⟩ := by
    intro i hi
    obtain ⟨h, hne⟩ := (mem_switches v i).1 hi
    exact ⟨h, by simp [succ', h], hne⟩
  have hinj : Set.InjOn succ' S := by
    intro a ha b hb hab
    obtain ⟨h1, e1, -⟩ := hsucc a ha
    obtain ⟨h2, e2, -⟩ := hsucc b hb
    rw [e1, e2] at hab
    apply Fin.ext
    have := congrArg Fin.val hab
    simp only at this
    omega
  have h1 : ∑ i ∈ S, lam i ≤ ∑ i, lam i :=
    sum_le_sum_of_subset_of_nonneg (subset_univ _) (fun i _ _ => hlam i)
  have h2 : ∑ i ∈ S, lam (succ' i) ≤ ∑ i, lam i := by
    rw [← sum_image hinj]
    exact sum_le_sum_of_subset_of_nonneg (subset_univ _) (fun i _ _ => hlam i)
  have h3 : ∑ i ∈ S, c₀ ≤ ∑ i ∈ S, (lam i + lam (succ' i)) := by
    apply sum_le_sum
    intro i hi
    obtain ⟨h, e, hne⟩ := hsucc i hi
    rw [e]
    exact hsw i h hne
  rw [sum_const, nsmul_eq_mul, sum_add_distrib] at h3
  rw [le_div_iff₀ hc]
  linarith

/-- At most `M · (1 + #switches)` chords. -/
theorem card_le_mul_switches {α : Type*} {n M : ℕ} (v : Fin n → α)
    (hrun : ∀ i j : Fin n, i ≤ j → (∀ m, i ≤ m → m ≤ j → v m = v i) → j.1 + 1 - i.1 ≤ M) :
    n ≤ M * (1 + (switches v).card) := by
  classical
  set S := switches v
  -- the block index of chord `i`: the number of switches before it
  let b : Fin n → ℕ := fun i => (S.filter (· < i)).card
  have hb_le : ∀ i, b i ≤ S.card := fun i => card_le_card (filter_subset _ _)
  have hb_lt : ∀ i j m : Fin n, i ≤ m → m < j → m ∈ S → b i < b j := by
    intro i j m him hmj hmS
    apply card_lt_card
    rw [ssubset_iff_of_subset]
    · exact ⟨m, mem_filter.2 ⟨hmS, hmj⟩, fun h => absurd (mem_filter.1 h).2 (not_lt.2 him)⟩
    · intro x hx
      rw [mem_filter] at hx ⊢
      exact ⟨hx.1, lt_of_lt_of_le hx.2 (le_of_lt (lt_of_le_of_lt him hmj))⟩
  have hfib : ∀ r ∈ univ.image b, (univ.filter fun x => b x = r).card ≤ M := by
    intro r hr
    set F := univ.filter fun x => b x = r
    have hne : F.Nonempty := by
      obtain ⟨x, -, hx⟩ := mem_image.1 hr
      exact ⟨x, mem_filter.2 ⟨mem_univ _, hx⟩⟩
    set i := F.min' hne
    set j := F.max' hne
    have hiF : b i = r := (mem_filter.1 (F.min'_mem hne)).2
    have hjF : b j = r := (mem_filter.1 (F.max'_mem hne)).2
    have hij : i ≤ j := F.min'_le _ (F.max'_mem hne)
    have hconst : ∀ k : ℕ, ∀ m : Fin n, m.1 = i.1 + k → m ≤ j → v m = v i := by
      intro k
      induction k with
      | zero =>
        intro m hm _
        rw [show m = i from Fin.ext (by simpa using hm)]
      | succ k ih =>
        intro m hm hmj
        have hk : i.1 + k < n := by omega
        have h' : (⟨i.1 + k, hk⟩ : Fin n).1 + 1 < n := by simp only; omega
        have e : (⟨(⟨i.1 + k, hk⟩ : Fin n).1 + 1, h'⟩ : Fin n) = m := Fin.ext (by simp only; omega)
        have hprev := ih ⟨i.1 + k, hk⟩ rfl (by rw [Fin.le_def] at hmj ⊢; simp only; omega)
        by_contra hne'
        have hmS : (⟨i.1 + k, hk⟩ : Fin n) ∈ S :=
          (mem_switches v _).2 ⟨h', by rw [e, hprev]; exact fun h => hne' h.symm⟩
        have := hb_lt i j ⟨i.1 + k, hk⟩ (by rw [Fin.le_def]; simp only; omega)
          (by rw [Fin.lt_def]; rw [Fin.le_def] at hmj; simp only; omega) hmS
        omega
    have hsub : F ⊆ Icc i j := by
      intro x hx
      rw [mem_Icc]
      exact ⟨F.min'_le x hx, F.le_max' x hx⟩
    calc F.card ≤ (Icc i j).card := card_le_card hsub
      _ = j.1 + 1 - i.1 := Fin.card_Icc i j
      _ ≤ M := hrun i j hij (fun m him hmj =>
          hconst (m.1 - i.1) m (by rw [Fin.le_def] at him; omega) hmj)
  have h1 := card_le_mul_card_image univ M hfib
  have h2 : (univ.image b).card ≤ S.card + 1 := by
    calc (univ.image b).card ≤ (range (S.card + 1)).card :=
          card_le_card (fun x hx => by
            obtain ⟨y, -, rfl⟩ := mem_image.1 hx
            exact mem_range.2 (Nat.lt_succ_of_le (hb_le y)))
      _ = S.card + 1 := card_range _
  rw [card_univ, Fintype.card_fin] at h1
  calc n ≤ M * (univ.image b).card := h1
    _ ≤ M * (S.card + 1) := Nat.mul_le_mul_left _ h2
    _ = M * (1 + S.card) := by rw [Nat.add_comm]

/-- **The counting lemma.** -/
theorem fan_count {α : Type*} {n M : ℕ} (v : Fin n → α) (lam : Fin n → ℝ) {c₀ : ℝ}
    (hc : 0 < c₀) (hlam : ∀ i, 0 ≤ lam i)
    (hrun : ∀ i j : Fin n, i ≤ j → (∀ m, i ≤ m → m ≤ j → v m = v i) → j.1 + 1 - i.1 ≤ M)
    (hsw : ∀ i : Fin n, ∀ h : i.1 + 1 < n, v i ≠ v ⟨i.1 + 1, h⟩ → c₀ ≤ lam i + lam ⟨i.1 + 1, h⟩) :
    (n : ℝ) ≤ M * (1 + 2 * (∑ i, lam i) / c₀) := by
  have h1 := card_le_mul_switches v hrun
  have h2 := card_switches_le v lam hc hlam hsw
  have h1' : (n : ℝ) ≤ M * (1 + (switches v).card) := by exact_mod_cast h1
  have hM : (0 : ℝ) ≤ M := by positivity
  nlinarith

end DDG.Delaunay.Surface
