import DdgOracle.Scalar

/-!
# Finite sums

Two structurally recursive sums, chosen so that they unfold by `rfl`/`simp` and match
`Finset.sum_range_succ` and `List.sum_cons` in Mathlib.
-/

namespace Ddg

variable {α : Type} [Scalar α]

/-- `sumTo n f = f 0 + f 1 + ⋯ + f (n - 1)`, associated to the left. -/
def sumTo : Nat → (Nat → α) → α
  | 0, _ => 0
  | n + 1, f => sumTo n f + f n

/-- Sum of a list, `x₀ + (x₁ + (⋯ + 0))`. -/
def sumList : List α → α
  | [] => 0
  | x :: xs => x + sumList xs

end Ddg
