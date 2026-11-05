import DdgProofs.Developability.Vec
import DdgProofs.Optimization.Wolfe

/-!
# Developability: the weak Wolfe line search

`crates/ddg-developability/SPEC.md` §4.2 (L-BFGS: a pair is stored only if `sᵀy > 0`) and §4.3
(Lewis–Overton Algorithm 2.6; A. S. Lewis, M. L. Overton, *Nonsmooth optimization via
quasi-Newton methods*, Math. Program. A 141, 2013). The paper (Sec. 4.3) names the method.

* `weakWolfe_curvature_of_curvatureHolds` — the BFGS curvature condition `sᵀy > 0` for a step
  that satisfies the oracle's `curvatureHolds` (the statement on real vectors,
  `Dev.weakWolfe_curvature`, is in `DdgProofs.Optimization.Wolfe`);
* `lineSearch_spec` — about the oracle's `Ddg.Dev.lineSearch` (the loop of Algorithm 2.6 with the
  limits of SPEC §4.3): an `accepted t` outcome satisfies `S(t)` and `C(t)`; a `capped t` outcome
  satisfies `S(t)`; a `weak t` outcome satisfies `S(t)` and `t > 0`.

Termination and existence of weak Wolfe steps (L&O Thm. 2.5, 2.7) are cited
(`dev-lewis-overton-termination`), not proved; the loop here is bounded by its fuel.
-/

noncomputable section

namespace DDG

namespace Dev

open scoped DDG.Oracle

/-! ## The curvature condition -/

/-- **Weak Wolfe implies the curvature condition**, with the weak Wolfe condition as the oracle's
`curvatureHolds c2 (gᵀd) (g₊ᵀd)` (bridge to `Dev.weakWolfe_curvature`). -/
theorem weakWolfe_curvature_of_curvatureHolds {n : ℕ} (g gp d : Fin n → ℝ) {t c2 : ℝ} (ht : 0 < t)
    (hdesc : g ⬝ᵥ d < 0) (hc2 : c2 < 1)
    (hC : Ddg.Dev.curvatureHolds c2 (g ⬝ᵥ d) (gp ⬝ᵥ d) = true) :
    t * (c2 - 1) * (g ⬝ᵥ d) < (t • d) ⬝ᵥ (gp - g) ∧ 0 < (t • d) ⬝ᵥ (gp - g) := by
  simp only [Ddg.Dev.curvatureHolds, decide_eq_true_eq] at hC
  exact weakWolfe_curvature g gp d ht hdesc hc2 hC


/-! ## The line search -/

section lineSearch

open Ddg.Dev

variable (eval : ℝ → Option (ℝ × ℝ)) (c1 c2 s : ℝ)

/-- `S(t)`: the energy is defined at `x + t d` and `h(t) < c1 s t` (Armijo, strict). -/
def ArmijoAt (t : ℝ) : Prop := ∃ h gd, eval t = some (h, gd) ∧ h < c1 * s * t

/-- `C(t)`: the energy is defined at `x + t d` and `g(x + t d)ᵀd > c2 s` (weak Wolfe). -/
def CurvatureAt (t : ℝ) : Prop := ∃ h gd, eval t = some (h, gd) ∧ c2 * s < gd

/-- What each outcome guarantees. -/
def OutcomeOk : LineSearchOutcome ℝ → Prop
  | .accepted t => ArmijoAt eval c1 s t ∧ CurvatureAt eval c2 s t
  | .capped t => ArmijoAt eval c1 s t
  | .weak t => ArmijoAt eval c1 s t ∧ 0 < t
  | .failed => True

theorem outcomeOk_weak {lo : ℝ} (hpos : 0 < lo) (hlo : lo = 0 ∨ ArmijoAt eval c1 s lo) :
    OutcomeOk eval c1 c2 s (.weak lo) := by
  refine ⟨?_, hpos⟩
  rcases hlo with h0 | h0
  · exact absurd h0 hpos.ne'
  · exact h0

theorem outcomeOk_failed : OutcomeOk eval c1 c2 s .failed := trivial

/-- The invariant "`lo = 0` or `S(lo)`" gives a correct outcome of every run of the loop. -/
theorem lineSearchLoop_ok (tmax : ℝ) (maxD maxB : ℕ) :
    ∀ (fuel : ℕ) (lo : ℝ) (hi : Option ℝ) (t : ℝ) (nd nb : ℕ) (trials : List ℝ),
      (lo = 0 ∨ ArmijoAt eval c1 s lo) →
      OutcomeOk eval c1 c2 s
        (lineSearchLoop eval c1 c2 s tmax maxD maxB fuel lo hi t nd nb trials).1 := by
  intro fuel
  induction fuel with
  | zero =>
    intro lo hi t nd nb trials hlo
    simp only [lineSearchLoop]
    split_ifs with h
    · exact outcomeOk_weak eval c1 c2 s (by simpa using h) hlo
    · exact outcomeOk_failed eval c1 c2 s
  | succ fuel ih =>
    intro lo hi t nd nb trials hlo
    simp only [lineSearchLoop]
    -- the evaluation at `t`
    rcases he : eval t with _ | ⟨h, gd⟩
    · -- undefined energy: `S(t)` fails, `β ← t`, bisection or fallback
      simp only [Bool.false_and, Bool.false_eq_true, ite_false]
      split_ifs with hb hp
      · exact ih _ _ _ _ _ _ hlo
      · exact outcomeOk_weak eval c1 c2 s (by simpa using hp) hlo
      · exact outcomeOk_failed eval c1 c2 s
    · simp only
      by_cases hS : armijoHolds c1 s t h = true
      · have hSt : ArmijoAt eval c1 s t :=
          ⟨h, gd, he, by simpa [armijoHolds] using hS⟩
        by_cases hC : curvatureHolds c2 s gd = true
        · simp only [hS, hC, Bool.and_self, ite_true]
          exact ⟨hSt, h, gd, he, by simpa [curvatureHolds] using hC⟩
        · have hC' : curvatureHolds c2 s gd = false := by simpa using hC
          simp only [hS, hC', Bool.and_false, Bool.false_eq_true, ite_false, ite_true]
          cases hi with
          | none =>
            simp only
            split_ifs with h1 h2 hp
            · exact hSt
            · exact ih _ _ _ _ _ _ (Or.inr hSt)
            · exact outcomeOk_weak eval c1 c2 s (by simpa using hp) (Or.inr hSt)
            · exact outcomeOk_failed eval c1 c2 s
          | some b =>
            simp only
            split_ifs with h1 hp
            · exact ih _ _ _ _ _ _ (Or.inr hSt)
            · exact outcomeOk_weak eval c1 c2 s (by simpa using hp) (Or.inr hSt)
            · exact outcomeOk_failed eval c1 c2 s
      · have hS' : armijoHolds c1 s t h = false := by simpa using hS
        simp only [hS', Bool.false_and, Bool.false_eq_true, ite_false]
        split_ifs with hb hp
        · exact ih _ _ _ _ _ _ hlo
        · exact outcomeOk_weak eval c1 c2 s (by simpa using hp) hlo
        · exact outcomeOk_failed eval c1 c2 s

/-- **The Lewis–Overton line search returns what it claims** (SPEC §4.3): for the oracle's
`lineSearch eval c1 c2 s t_max maxDoublings maxBisections`, an outcome `accepted t` satisfies
`S(t)` (`h(t) < c1 s t`) and `C(t)` (`g(x + t d)ᵀd > c2 s`), `capped t` satisfies `S(t)`, and
`weak t` satisfies `S(t)` with `t > 0`. (A trial where the energy is undefined counts as a failure
of `S`.) Nothing is claimed about `failed`, nor about termination within the limits.

[formula: dev-lewis-overton-line-search] -/
theorem lineSearch_spec (tmax : ℝ) (maxD maxB : ℕ) :
    OutcomeOk eval c1 c2 s (lineSearch eval c1 c2 s tmax maxD maxB).1 :=
  lineSearchLoop_ok eval c1 c2 s tmax maxD maxB _ _ _ _ _ _ _ (Or.inl (by simp))

/-- The accepted step of the line search satisfies both conditions (corollary of
`lineSearch_spec`). -/
theorem lineSearch_accepted {tmax t : ℝ} {maxD maxB : ℕ}
    (h : (lineSearch eval c1 c2 s tmax maxD maxB).1 = .accepted t) :
    ArmijoAt eval c1 s t ∧ CurvatureAt eval c2 s t := by
  have := lineSearch_spec eval c1 c2 s tmax maxD maxB
  rw [h] at this
  exact this

end lineSearch

end Dev

end DDG
