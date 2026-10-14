import DdgProofs.Cuts.VertexLimit

/-!
# Cut flow: the stop rule of a run (energy criterion, lowest state, motion criterion)

`crates/ddg-variational-cuts/SPEC.md` §4 ("Stopping"), D27, D29, D31. The ids are
`cut-stop-energy` (the Rust `energy_stop`), `cut-run-lowest` (`run_outcome`) and
`cut-stop-motion` (`motion_criterion`). There is no oracle definition (the oracle replays single
steps); the definitions below are **this package's own**, written after the Rust functions.

* **`cut-stop-motion`** (D31) — `motionCriterion tol moved predicted` is
  `tol > 0 ∧ moved ≤ tol ∧ predicted ≤ tol`. It implies the criterion before D31
  (`moved ≤ tol`), never holds for `tol ≤ 0`, and is monotone in the predicted motion
  (`motionCriterion_spec`, `motionCriterion_mono`). With the length step on the cut the
  predicted motion is `τ` times the largest absolute weighted-mean rate at a cut-adjacent vertex
  (`predictedMotionOnCut`); for `τ ≥ 0` and weights `≥ 0` with a positive sum at every such
  vertex it is at most `τ · R` for every bound `R` of the rates, in particular at most
  `τ · max_p |r_p|` (`predictedMotionOnCut_le`, `predictedMotionOnCut_le_max`); it is `0` when
  all rates are `0`, for instance at a rest state `s + λ⁰ = 0`
  (`predictedMotionOnCut_zero`, `predictedMotionOnCut_rest`), and then the criterion is the old
  one, so a step that moved nothing converges for a positive tolerance
  (`motionCriterion_zero_rates`, `motionCriterion_rest`). With the length step on the mesh the
  predicted motion is the largest `|φ′ − φ|` over the cut-adjacent vertices before the local
  limit (`predictedMotionOnMesh`); the local limit does not exceed it
  (`limited_le_predictedMotionOnMesh`).

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

/-! ## The motion criterion (SPEC §4, D31) -/

section Motion

open scoped DdgProofs.Oracle

/-- The motion criterion of the stop rule as the Rust code evaluates it (`motion_criterion`):
the tolerance is positive, the step moved the cut by at most the tolerance, and the predicted
motion is at most the tolerance too. -/
def motionCriterion (tol moved predicted : ℝ) : Prop := 0 < tol ∧ moved ≤ tol ∧ predicted ≤ tol

/-- **The motion criterion**: it holds iff the tolerance is positive and both the cut motion
and the predicted motion are at most it; it implies the criterion before D31 (the cut motion is
at most the tolerance); and with a tolerance of zero or less a run never converges by motion,
also not after a step that moved nothing.

[formula: cut-stop-motion] -/
theorem motionCriterion_spec (tol moved predicted : ℝ) :
    (motionCriterion tol moved predicted ↔ 0 < tol ∧ moved ≤ tol ∧ predicted ≤ tol) ∧
      (motionCriterion tol moved predicted → moved ≤ tol) ∧
      (tol ≤ 0 → ¬ motionCriterion tol moved predicted) :=
  ⟨Iff.rfl, fun h => h.2.1, fun h0 h => absurd h.1 (not_lt.mpr h0)⟩

/-- A smaller predicted motion (and a smaller cut motion) keeps the criterion: the first version
of the prediction in D31, `τ · max_p |r_p|`, is the stricter one (`predictedMotionOnCut_le_max`). -/
theorem motionCriterion_mono {tol moved moved' predicted predicted' : ℝ}
    (h : motionCriterion tol moved predicted) (hm : moved' ≤ moved)
    (hp : predicted' ≤ predicted) : motionCriterion tol moved' predicted' :=
  ⟨h.1, hm.trans h.2.1, hp.trans h.2.2⟩

/-- The largest absolute value of `x` over a finite set, `0` over the empty set (the fold of the
Rust code). -/
def maxAbs {ι : Type*} (A : Finset ι) (x : ι → ℝ) : ℝ := A.fold max 0 fun i => |x i|

theorem maxAbs_le_iff {ι : Type*} {A : Finset ι} {x : ι → ℝ} {B : ℝ} :
    maxAbs A x ≤ B ↔ 0 ≤ B ∧ ∀ i ∈ A, |x i| ≤ B :=
  Finset.fold_max_le B

theorem maxAbs_nonneg {ι : Type*} (A : Finset ι) (x : ι → ℝ) : 0 ≤ maxAbs A x :=
  (maxAbs_le_iff.mp le_rfl).1

theorem abs_le_maxAbs {ι : Type*} {A : Finset ι} (x : ι → ℝ) {i : ι} (hi : i ∈ A) :
    |x i| ≤ maxAbs A x :=
  (maxAbs_le_iff.mp le_rfl).2 i hi

theorem maxAbs_zero {ι : Type*} (A : Finset ι) : maxAbs A (fun _ => (0 : ℝ)) = 0 :=
  le_antisymm (maxAbs_le_iff.mpr ⟨le_rfl, fun _ _ => by rw [abs_zero]⟩) (maxAbs_nonneg _ _)

variable {V P : Type*} [Fintype P]

/-- **The predicted motion with the length step on the cut**: `τ` times the largest absolute
rate that a cut-adjacent vertex `i ∈ A` gets as the weighted mean of the rates of the cut points
on its edges, `Σ_p c_ip r_ip / Σ_p c_ip` with `c_ip = ℓ_p β_ip` (the first line of the local
extension of SPEC §3.8). `r i p` is the rate of the cut point `p` as the vertex `i` sees it (on
a network the rate of the copy of `p` on the side of `i`; with two regions it does not depend on
`i`). -/
def predictedMotionOnCut (τ : ℝ) (A : Finset V) (c r : V → P → ℝ) : ℝ :=
  τ * maxAbs A fun i => wmean Finset.univ (c i) (r i)

/-- **On the cut the predicted motion is at most `τ` times any bound of the rates** (`τ ≥ 0`;
at every cut-adjacent vertex the weights are `≥ 0` and not all zero): each vertex gets a convex
combination of the rates of its cut points.

[formula: cut-stop-motion] -/
theorem predictedMotionOnCut_le {τ : ℝ} (hτ : 0 ≤ τ) {A : Finset V} {c r : V → P → ℝ}
    (hc : ∀ i ∈ A, ∀ p, 0 ≤ c i p) (hpos : ∀ i ∈ A, 0 < ∑ p, c i p) {R : ℝ} (hR0 : 0 ≤ R)
    (hR : ∀ i ∈ A, ∀ p, c i p ≠ 0 → |r i p| ≤ R) :
    predictedMotionOnCut τ A c r ≤ τ * R := by
  unfold predictedMotionOnCut
  refine mul_le_mul_of_nonneg_left (maxAbs_le_iff.mpr ⟨hR0, fun i hi => ?_⟩) hτ
  have h := wmean_mem (s := Finset.univ) (c := c i) (x := r i) (fun p _ => hc i hi p)
    (hpos i hi) (lo := -R) (hi := R) (fun p _ hp => (abs_le.mp (hR i hi p hp)).1)
    (fun p _ hp => (abs_le.mp (hR i hi p hp)).2)
  exact abs_le.mpr h

/-- **… in particular at most `τ · max_p |r_p|`** over the cut points (two regions: one rate
per cut point).

[formula: cut-stop-motion] -/
theorem predictedMotionOnCut_le_max {τ : ℝ} (hτ : 0 ≤ τ) {A : Finset V} {c : V → P → ℝ}
    (hc : ∀ i ∈ A, ∀ p, 0 ≤ c i p) (hpos : ∀ i ∈ A, 0 < ∑ p, c i p) (r : P → ℝ) :
    predictedMotionOnCut τ A c (fun _ => r) ≤ τ * maxAbs Finset.univ r :=
  predictedMotionOnCut_le hτ hc hpos (maxAbs_nonneg _ _)
    fun _ _ p _ => abs_le_maxAbs r (Finset.mem_univ p)

/-- The predicted motion is nonnegative for `τ ≥ 0`. -/
theorem predictedMotionOnCut_nonneg {τ : ℝ} (hτ : 0 ≤ τ) (A : Finset V) (c r : V → P → ℝ) :
    0 ≤ predictedMotionOnCut τ A c r :=
  mul_nonneg hτ (maxAbs_nonneg _ _)

/-- **With all rates zero the predicted motion is zero**, whatever the weights and `τ`. -/
theorem predictedMotionOnCut_zero (τ : ℝ) (A : Finset V) (c : V → P → ℝ) :
    predictedMotionOnCut τ A c (fun _ _ => 0) = 0 := by
  unfold predictedMotionOnCut
  have e : (fun i => wmean Finset.univ (c i) fun _ => (0 : ℝ)) = fun _ => (0 : ℝ) := by
    funext i
    unfold wmean
    simp
  rw [e, maxAbs_zero, mul_zero]

/-- **With all rates zero the criterion is the one before D31**; in particular a step that
moved nothing converges for every positive tolerance.

[formula: cut-stop-motion] -/
theorem motionCriterion_zero_rates (tol moved τ : ℝ) (A : Finset V) (c : V → P → ℝ) :
    (motionCriterion tol moved (predictedMotionOnCut τ A c fun _ _ => 0) ↔
        0 < tol ∧ moved ≤ tol) ∧
      (0 < tol → motionCriterion tol 0 (predictedMotionOnCut τ A c fun _ _ => 0)) := by
  rw [predictedMotionOnCut_zero]
  exact ⟨⟨fun h => ⟨h.1, h.2.1⟩, fun h => ⟨h.1, h.2, h.1.le⟩⟩, fun h => ⟨h, h.le, h.le⟩⟩

/-- **The predicted motion with the length step on the mesh**: the largest `|φ′(i) − φ(i)|` over
the vertices `A` that are cut-adjacent before the step, `φ′` after the explicit and the length
step and before the local limit. -/
def predictedMotionOnMesh (A : Finset V) (φ φ' : V → ℝ) : ℝ := maxAbs A fun i => φ' i - φ i

/-- If the two steps change no cut-adjacent vertex the predicted motion is zero. -/
theorem predictedMotionOnMesh_zero {A : Finset V} {φ φ' : V → ℝ} (h : ∀ i ∈ A, φ' i = φ i) :
    predictedMotionOnMesh A φ φ' = 0 :=
  le_antisymm (maxAbs_le_iff.mpr ⟨le_rfl, fun i hi => by rw [h i hi, sub_self, abs_zero]⟩)
    (maxAbs_nonneg _ _)

/-- **The local limit at the vertices does not exceed the predicted motion**: with factors in
`[0, 1]`, the largest change after the limit over the cut-adjacent vertices is at most the
largest change before it. -/
theorem limited_le_predictedMotionOnMesh {A : Finset V} (φ φ' : V → ℝ) {f : V → ℝ}
    (hf0 : ∀ i ∈ A, 0 ≤ f i) (hf1 : ∀ i ∈ A, f i ≤ 1) :
    maxAbs A (fun i => Ddg.limitedValue (φ i) (φ' i) (f i) - φ i) ≤
      predictedMotionOnMesh A φ φ' :=
  maxAbs_le_iff.mpr ⟨maxAbs_nonneg _ _, fun i hi =>
    (abs_limitedValue_sub_le (hf0 i hi) (hf1 i hi) (φ i) (φ' i)).trans
      (abs_le_maxAbs (fun i => φ' i - φ i) hi)⟩

end Motion

section MotionRest

open Matrix

variable {V C : Type*} [Fintype C] [DecidableEq C]

/-- **At a rest state the predicted motion is zero, whatever the step**: if `s + λ⁰ = 0` at
every cut point, the rate of the system on the cut (SPEC §3.10) is zero (`onCut_rest`), and so
is the predicted motion. -/
theorem predictedMotionOnCut_rest {ℓ s lam r : C → ℝ} {m : C → C → ℝ} {τ : ℝ}
    (hℓ : ∀ p, 0 < ℓ p) (hm0 : ∀ p q, 0 ≤ m p q) (hτ : 0 ≤ τ)
    (h : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r = fun p => Ddg.fluxRhs (ℓ p * (s p + lam p)))
    (h0 : ∀ p, s p + lam p = 0) (A : Finset V) (c : V → C → ℝ) :
    predictedMotionOnCut τ A c (fun _ => r) = 0 := by
  have hr : r = fun _ => 0 := funext fun p => onCut_rest hℓ hm0 hτ h h0 p
  rw [hr]
  exact predictedMotionOnCut_zero τ A c

/-- **A run at a rest state converges by motion**: with `s + λ⁰ = 0` at every cut point, a step
that moved nothing satisfies the criterion for every positive tolerance.

[formula: cut-stop-motion] -/
theorem motionCriterion_rest {ℓ s lam r : C → ℝ} {m : C → C → ℝ} {τ tol : ℝ}
    (hℓ : ∀ p, 0 < ℓ p) (hm0 : ∀ p q, 0 ≤ m p q) (hτ : 0 ≤ τ)
    (h : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r = fun p => Ddg.fluxRhs (ℓ p * (s p + lam p)))
    (h0 : ∀ p, s p + lam p = 0) (htol : 0 < tol) (A : Finset V) (c : V → C → ℝ) :
    motionCriterion tol 0 (predictedMotionOnCut τ A c fun _ => r) := by
  rw [predictedMotionOnCut_rest hℓ hm0 hτ h h0]
  exact ⟨htol, htol.le, htol.le⟩

end MotionRest

end Cuts

end DdgProofs
