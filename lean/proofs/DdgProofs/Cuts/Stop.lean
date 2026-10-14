import Mathlib

/-!
# Cut flow: the energy criterion of the stop rule and the lowest state of a run

`crates/ddg-variational-cuts/SPEC.md` §4 ("Stopping"), D27, D29. The ids are `cut-stop-energy`
(the Rust `energy_stop`) and `cut-run-lowest` (`run_outcome`). There is no oracle definition
(the oracle replays single steps); the definitions below are **this package's own**, written
after the Rust functions.

With `E_k` the energy of the state after `k` steps, `W` the window and `ε` the tolerance:

```
E_{k−W} − E_k >  ε |E_k|        the run goes on
|E_{k−W} − E_k| ≤ ε |E_k|       EnergyStalled
E_k − E_{k−W} >  ε |E_k|        EnergyRose: the flow is set back to the state of lowest energy
```

* **The three cases are exclusive and exhaustive** for `ε ≥ 0` (`energy_cases`); for `ε < 0`
  the first and the third can hold together (`energy_cases_overlap_of_neg`; the Rust code does
  not evaluate the criterion for `ε ≤ 0`). `energyVerdict` is the chain of comparisons of the
  Rust code; it returns the case of the table (`energyVerdict_goesOn_iff`,
  `energyVerdict_stalled_iff`, `energyVerdict_rose_iff`). A run stops by the criterion exactly
  when the energy has fallen by at most `ε |E_k|` over the window, as before D27
  (`energyVerdict_stops_iff`). With `ε = 0` the table would stop every run whose energy has not
  strictly fallen (`energyVerdict_zero_stops_iff`).
* **The lowest state.** `lowestState E n` scans the states `0, …, n` and replaces the state it
  holds only by one of strictly lower energy (the Rust loop). It returns the earliest state of
  least energy (`lowestState_spec`), and that state is unique (`lowestState_unique`).
  `runOutcome` is `run_outcome`: it compares the lowest of the states before the steps with the
  state after the last step. Without a vanished cut the result is `lowestState` over all states
  of the run (`runOutcome_eq_lowestState`, `runOutcome_no_steps`); with `EnergyRose` the state
  returned is that one: its energy is at most that of every state of the run, and less than that
  of every earlier one (`runOutcome_rose`); with another reason the run stays in its last state
  (`runOutcome_not_rose`); if the cut has vanished the last state does not count
  (`runOutcome_vanished`). After `EnergyRose` the lowest energy is strictly below `E_k`
  (`lowestState_lt_of_rose`).

All of this is exact real arithmetic. A NaN energy fails both strict comparisons of the Rust
code and is reported as `EnergyStalled`; that is outside these statements.
-/

noncomputable section

namespace DdgProofs

namespace Cuts

/-! ## The energy criterion -/

/-- The three outcomes of the energy criterion of the stop rule. -/
inductive EnergyVerdict
  | goesOn
  | stalled
  | rose
  deriving DecidableEq

/-- The energy criterion as the Rust code evaluates it (`energy_stop`), for the tolerance `ε`,
the energy `Eb = E_{k−W}` at the start of the window and the energy `Ek = E_k` at its end:
first "fell by more than `ε |E_k|`", then "rose by more than that", else "stood still". -/
def energyVerdict (ε Eb Ek : ℝ) : EnergyVerdict :=
  if ε * |Ek| < Eb - Ek then .goesOn else if ε * |Ek| < Ek - Eb then .rose else .stalled

/-- **The three cases of the energy criterion are exclusive and exhaustive** (`ε ≥ 0`): exactly
one of `E_{k−W} − E_k > ε|E_k|`, `|E_{k−W} − E_k| ≤ ε|E_k|`, `E_k − E_{k−W} > ε|E_k|` holds.

[formula: cut-stop-energy] -/
theorem energy_cases {ε : ℝ} (hε : 0 ≤ ε) (Eb Ek : ℝ) :
    (ε * |Ek| < Eb - Ek ∧ ¬ |Eb - Ek| ≤ ε * |Ek| ∧ ¬ ε * |Ek| < Ek - Eb) ∨
      (¬ ε * |Ek| < Eb - Ek ∧ |Eb - Ek| ≤ ε * |Ek| ∧ ¬ ε * |Ek| < Ek - Eb) ∨
      (¬ ε * |Ek| < Eb - Ek ∧ ¬ |Eb - Ek| ≤ ε * |Ek| ∧ ε * |Ek| < Ek - Eb) := by
  have ht : 0 ≤ ε * |Ek| := mul_nonneg hε (abs_nonneg _)
  by_cases hA : ε * |Ek| < Eb - Ek
  · left
    refine ⟨hA, fun hB => ?_, fun hC => ?_⟩
    · have := (abs_le.mp hB).2
      linarith
    · linarith
  · by_cases hC : ε * |Ek| < Ek - Eb
    · right; right
      refine ⟨hA, fun hB => ?_, hC⟩
      have := (abs_le.mp hB).1
      linarith
    · right; left
      exact ⟨hA, abs_le.mpr ⟨by linarith [not_lt.mp hC], not_lt.mp hA⟩, hC⟩

/-- The three cases are exhaustive for every `ε`; the second excludes the other two for every
`ε`. -/
theorem energy_cases_exhaustive (ε Eb Ek : ℝ) :
    ε * |Ek| < Eb - Ek ∨ |Eb - Ek| ≤ ε * |Ek| ∨ ε * |Ek| < Ek - Eb := by
  by_cases hA : ε * |Ek| < Eb - Ek
  · exact Or.inl hA
  · by_cases hC : ε * |Ek| < Ek - Eb
    · exact Or.inr (Or.inr hC)
    · exact Or.inr (Or.inl (abs_le.mpr ⟨by linarith [not_lt.mp hC], not_lt.mp hA⟩))

/-- For a negative tolerance the first and the third case can hold together (`ε = −1`,
`E_{k−W} = E_k = 1`): the hypothesis `ε ≥ 0` of `energy_cases` is needed. -/
theorem energy_cases_overlap_of_neg :
    ∃ ε Eb Ek : ℝ, ε * |Ek| < Eb - Ek ∧ ε * |Ek| < Ek - Eb :=
  ⟨-1, 1, 1, by norm_num, by norm_num⟩

/-- The run goes on exactly in the first case. -/
theorem energyVerdict_goesOn_iff (ε Eb Ek : ℝ) :
    energyVerdict ε Eb Ek = .goesOn ↔ ε * |Ek| < Eb - Ek := by
  unfold energyVerdict
  split_ifs with h1 h2
  · exact ⟨fun _ => h1, fun _ => rfl⟩
  · exact ⟨fun h => absurd h (by decide), fun h => absurd h h1⟩
  · exact ⟨fun h => absurd h (by decide), fun h => absurd h h1⟩

/-- `EnergyStalled` exactly in the second case, `|E_{k−W} − E_k| ≤ ε|E_k|`. -/
theorem energyVerdict_stalled_iff (ε Eb Ek : ℝ) :
    energyVerdict ε Eb Ek = .stalled ↔ |Eb - Ek| ≤ ε * |Ek| := by
  unfold energyVerdict
  split_ifs with h1 h2
  · refine ⟨fun h => absurd h (by decide), fun h => ?_⟩
    have := (abs_le.mp h).2
    exact absurd this (not_le.mpr h1)
  · refine ⟨fun h => absurd h (by decide), fun h => ?_⟩
    have := (abs_le.mp h).1
    linarith
  · exact ⟨fun _ => abs_le.mpr ⟨by linarith [not_lt.mp h2], not_lt.mp h1⟩, fun _ => rfl⟩

/-- `EnergyRose` exactly in the third case, `E_k − E_{k−W} > ε|E_k|` (`ε ≥ 0`). -/
theorem energyVerdict_rose_iff {ε : ℝ} (hε : 0 ≤ ε) (Eb Ek : ℝ) :
    energyVerdict ε Eb Ek = .rose ↔ ε * |Ek| < Ek - Eb := by
  have ht : 0 ≤ ε * |Ek| := mul_nonneg hε (abs_nonneg _)
  unfold energyVerdict
  split_ifs with h1 h2
  · exact ⟨fun h => absurd h (by decide), fun h => by linarith⟩
  · exact ⟨fun _ => h2, fun _ => rfl⟩
  · exact ⟨fun h => absurd h (by decide), fun h => absurd h h2⟩

/-- **When a run stops is unchanged by D27**: the criterion ends the run (stalled or rose)
exactly when the energy has fallen by at most `ε|E_k|` over the window. -/
theorem energyVerdict_stops_iff (ε Eb Ek : ℝ) :
    energyVerdict ε Eb Ek ≠ .goesOn ↔ Eb - Ek ≤ ε * |Ek| := by
  rw [Ne, energyVerdict_goesOn_iff, not_lt]

/-- With the tolerance `0` the table stops every run whose energy has not strictly fallen over
the window. (The Rust code does not evaluate the criterion for `ε ≤ 0`.) -/
theorem energyVerdict_zero_stops_iff (Eb Ek : ℝ) : energyVerdict 0 Eb Ek ≠ .goesOn ↔ Eb ≤ Ek := by
  rw [energyVerdict_stops_iff, zero_mul, sub_nonpos]

/-! ## The state of lowest energy -/

/-- The lowest energy among the states `0, …, n` with the energies `E`, and its state: the scan
of the Rust code, which replaces the state it holds only by one of strictly lower energy. -/
def lowestState (E : ℕ → ℝ) : ℕ → ℝ × ℕ
  | 0 => (E 0, 0)
  | n + 1 => if E (n + 1) < (lowestState E n).1 then (E (n + 1), n + 1) else lowestState E n

theorem lowestState_zero (E : ℕ → ℝ) : lowestState E 0 = (E 0, 0) := rfl

theorem lowestState_succ (E : ℕ → ℝ) (n : ℕ) :
    lowestState E (n + 1) =
      if E (n + 1) < (lowestState E n).1 then (E (n + 1), n + 1) else lowestState E n := rfl

/-- **`lowestState` is the earliest state of least energy**: its state `j ≤ n` has the energy
it reports, `E_j ≤ E_i` for all `i ≤ n`, and `E_j < E_i` for all `i < j`. -/
theorem lowestState_spec (E : ℕ → ℝ) (n : ℕ) :
    (lowestState E n).2 ≤ n ∧ (lowestState E n).1 = E (lowestState E n).2 ∧
      (∀ i, i ≤ n → (lowestState E n).1 ≤ E i) ∧
      ∀ i, i < (lowestState E n).2 → (lowestState E n).1 < E i := by
  induction n with
  | zero =>
    refine ⟨le_rfl, rfl, fun i hi => ?_, fun i hi => absurd hi (Nat.not_lt_zero i)⟩
    rw [Nat.le_zero.mp hi]
    exact le_rfl
  | succ n ih =>
    obtain ⟨h1, h2, h3, h4⟩ := ih
    by_cases hlt : E (n + 1) < (lowestState E n).1
    · have e : lowestState E (n + 1) = (E (n + 1), n + 1) := by
        rw [lowestState_succ, ite_eq_left hlt]
      rw [e]
      refine ⟨le_rfl, rfl, fun i hi => ?_, fun i hi => ?_⟩
      · rcases Nat.lt_or_ge i (n + 1) with hi' | hi'
        · exact hlt.le.trans (h3 i (Nat.lt_succ_iff.mp hi'))
        · rw [le_antisymm hi hi']
      · exact lt_of_lt_of_le hlt (h3 i (Nat.lt_succ_iff.mp hi))
    · have e : lowestState E (n + 1) = lowestState E n := by
        rw [lowestState_succ, ite_eq_right hlt]
      rw [e]
      refine ⟨h1.trans (Nat.le_succ n), h2, fun i hi => ?_, h4⟩
      rcases Nat.lt_or_ge i (n + 1) with hi' | hi'
      · exact h3 i (Nat.lt_succ_iff.mp hi')
      · rw [le_antisymm hi hi']
        exact not_lt.mp hlt

/-- The earliest state of least energy is unique: `lowestState` returns it. -/
theorem lowestState_unique {E : ℕ → ℝ} {n j : ℕ} (hj : j ≤ n) (hmin : ∀ i, i ≤ n → E j ≤ E i)
    (hfirst : ∀ i, i < j → E j < E i) : (lowestState E n).2 = j := by
  obtain ⟨h1, h2, h3, h4⟩ := lowestState_spec E n
  rcases Nat.lt_trichotomy (lowestState E n).2 j with h | h | h
  · have a := hfirst _ h
    have b := h3 j hj
    rw [h2] at b
    linarith
  · exact h
  · have a := h4 j h
    have b := hmin _ h1
    rw [h2] at a
    linarith

/-- `run_outcome` of the Rust code. `lowest`: the lowest energy among the states before the
steps and its state (none if no step was taken); `last`: the energy after the last of `n` steps;
`rose`, `vanished`: the stop reason is `EnergyRose`, `CutVanished`. Returns the lowest energy of
the run, its state, and the state the run returns. -/
def runOutcome (rose vanished : Bool) (lowest : Option (ℝ × ℕ)) (last : ℝ) (n : ℕ) :
    ℝ × ℕ × ℕ :=
  let best : ℝ × ℕ :=
    match lowest with
    | some (e, s) => if e ≤ last ∨ vanished = true then (e, s) else (last, n)
    | none => (last, n)
  (best.1, best.2, if rose = true then best.2 else n)

/-- A run without a step: the start state is the lowest and the one returned. -/
theorem runOutcome_no_steps (rose : Bool) (E0 : ℝ) :
    runOutcome rose false none E0 0 = (E0, 0, 0) := by
  cases rose <;> rfl

/-- **Without a vanished cut the lowest state of a run is `lowestState` over all its states**
(the state after the last step included; ties go to the earlier state), and the state returned
is that one for `EnergyRose` and the last one otherwise. -/
theorem runOutcome_eq_lowestState (rose : Bool) (E : ℕ → ℝ) (n : ℕ) :
    runOutcome rose false (some (lowestState E n)) (E (n + 1)) (n + 1) =
      ((lowestState E (n + 1)).1, (lowestState E (n + 1)).2,
        if rose = true then (lowestState E (n + 1)).2 else n + 1) := by
  by_cases hlt : E (n + 1) < (lowestState E n).1
  · have e : lowestState E (n + 1) = (E (n + 1), n + 1) := by
      rw [lowestState_succ, ite_eq_left hlt]
    have hn : ¬ ((lowestState E n).1 ≤ E (n + 1) ∨ false = true) := by
      rintro (h | h)
      · exact absurd h (not_le.mpr hlt)
      · exact absurd h (by decide)
    simp only [runOutcome, e, ite_eq_right hn]
  · have e : lowestState E (n + 1) = lowestState E n := by
      rw [lowestState_succ, ite_eq_right hlt]
    have hp : (lowestState E n).1 ≤ E (n + 1) ∨ false = true := Or.inl (not_lt.mp hlt)
    simp only [runOutcome, e, ite_eq_left hp]

/-- **After `EnergyRose` the state returned has the lowest energy of the run** (every state of
such a run has a cut): with `j` the state returned after `n + 1` steps, `j ≤ n + 1`, the lowest
energy reported is `E_j`, `E_j ≤ E_i` for every state `i ≤ n + 1`, and `E_j < E_i` for every
earlier state `i < j`.

[formula: cut-run-lowest] -/
theorem runOutcome_rose (E : ℕ → ℝ) (n : ℕ) :
    (runOutcome true false (some (lowestState E n)) (E (n + 1)) (n + 1)).2.2 =
        (runOutcome true false (some (lowestState E n)) (E (n + 1)) (n + 1)).2.1 ∧
      (runOutcome true false (some (lowestState E n)) (E (n + 1)) (n + 1)).2.2 ≤ n + 1 ∧
      (runOutcome true false (some (lowestState E n)) (E (n + 1)) (n + 1)).1 =
        E (runOutcome true false (some (lowestState E n)) (E (n + 1)) (n + 1)).2.2 ∧
      (∀ i, i ≤ n + 1 →
        E (runOutcome true false (some (lowestState E n)) (E (n + 1)) (n + 1)).2.2 ≤ E i) ∧
      ∀ i, i < (runOutcome true false (some (lowestState E n)) (E (n + 1)) (n + 1)).2.2 →
        E (runOutcome true false (some (lowestState E n)) (E (n + 1)) (n + 1)).2.2 < E i := by
  obtain ⟨h1, h2, h3, h4⟩ := lowestState_spec E (n + 1)
  rw [runOutcome_eq_lowestState]
  simp only [ite_true]
  exact ⟨trivial, h1, h2, fun i hi => h2 ▸ h3 i hi, fun i hi => h2 ▸ h4 i hi⟩

/-- With a stop reason other than `EnergyRose` the run stays in its last state. -/
theorem runOutcome_not_rose (vanished : Bool) (lowest : Option (ℝ × ℕ)) (last : ℝ) (n : ℕ) :
    (runOutcome false vanished lowest last n).2.2 = n := by
  simp [runOutcome]

/-- **If the cut has vanished the last state does not count**: the lowest state reported is the
lowest among the states before the steps, which have a cut, whatever `last` is; the run stays
in its last state. -/
theorem runOutcome_vanished (E : ℕ → ℝ) (n : ℕ) (last : ℝ) :
    runOutcome false true (some (lowestState E n)) last (n + 1) =
      ((lowestState E n).1, (lowestState E n).2, n + 1) := by
  simp [runOutcome]

/-- **After `EnergyRose` the lowest energy of the run is strictly below `E_k`**, the energy at
the end of the window (`ε ≥ 0`, `k ≤ n`): the run is set back, or its last state is lower. -/
theorem lowestState_lt_of_rose {ε : ℝ} (hε : 0 ≤ ε) {E : ℕ → ℝ} {k n : ℕ} (W : ℕ) (hk : k ≤ n)
    (h : energyVerdict ε (E (k - W)) (E k) = .rose) : (lowestState E n).1 < E k := by
  have h1 := (energyVerdict_rose_iff hε _ _).mp h
  have ht : 0 ≤ ε * |E k| := mul_nonneg hε (abs_nonneg _)
  have h2 := (lowestState_spec E n).2.2.1 (k - W) ((Nat.sub_le k W).trans hk)
  linarith

end Cuts

end DdgProofs
