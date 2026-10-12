/-!
# The scalar interface

Every formula of the oracle is written once, against the small interface `Ddg.Scalar`, and is
then used at two types:

* `Float`, for execution (the instance is at the end of this file);
* the real numbers, in `lean/proofs`, where theorems are proved about the very same definitions.

The class carries operations only, no laws: `Float` satisfies none of the field axioms, and the
proofs package gets its laws from the instance it supplies for `ℝ`.

The instances that make numerals and `if a < b then … else …` work for a generic scalar are
`scoped`, so that importing this package does not add global instances to a downstream project.
They are active inside `namespace Ddg` and after `open Ddg`.
-/

namespace Ddg

/-- Arithmetic, order, `sqrt`, `atan` and `π`: everything the formulas of the oracle need.

Field by field, the intended instance at the real numbers is the obvious one
(`ofNat n = (n : ℝ)`, `sqrt = Real.sqrt`, `atan = Real.arctan`, `pi = Real.pi`, `max`/`min` the
lattice operations, the decidability instances classical). -/
class Scalar (α : Type) extends Add α, Sub α, Mul α, Div α, Neg α, LT α, LE α, Max α, Min α where
  /-- Embedding of the natural numbers; numerals such as `(2 : α)` elaborate to `ofNat 2`. -/
  ofNat : Nat → α
  /-- Square root. -/
  sqrt : α → α
  /-- Inverse tangent, with values in `(-π/2, π/2)`. -/
  atan : α → α
  /-- The constant `π`. -/
  pi : α
  /-- Decidability of `<`, used by the active sets. -/
  decLt : (a b : α) → Decidable (a < b)
  /-- Decidability of `≤`. -/
  decLe : (a b : α) → Decidable (a ≤ b)

/-- Numerals for a generic scalar: `(n : α) = Scalar.ofNat n`. -/
scoped instance (priority := low) instOfNatScalar {α : Type} [Scalar α] {n : Nat} : OfNat α n :=
  ⟨Scalar.ofNat n⟩

/-- `a < b` is decidable for a generic scalar. -/
scoped instance (priority := low) instDecidableLtScalar {α : Type} [Scalar α] (a b : α) :
    Decidable (a < b) :=
  Scalar.decLt a b

/-- `a ≤ b` is decidable for a generic scalar. -/
scoped instance (priority := low) instDecidableLeScalar {α : Type} [Scalar α] (a b : α) :
    Decidable (a ≤ b) :=
  Scalar.decLe a b

/-- Absolute value, `|x| = max x (-x)`. -/
def Scalar.abs {α : Type} [Scalar α] (x : α) : α := max x (-x)

/-- The IEEE-754 double nearest to `π`. -/
def floatPi : Float := Float.ofBits 0x400921FB54442D18

/-- `Float.ofNat`, computed by the machine conversion whenever `n < 2⁶⁴` (both round to
nearest). `Float.ofNat` itself goes through arbitrary-precision arithmetic on every call, which
dominates the run time when the numerals `0`, `1`, `2` of the generic formulas are evaluated
millions of times. -/
def floatOfNat (n : Nat) : Float :=
  let k := n.toUInt64
  if k.toNat = n then k.toFloat else Float.ofNat n

/-- `0.0`. A `Float` literal inside a loop is evaluated anew, and slowly, in every iteration;
hot loops use this constant instead. -/
def fzero : Float := 0.0

/-- The executable instance. `max a b = if a ≤ b then b else a` and
`min a b = if a ≤ b then a else b`; `atan`, `sqrt` are the C library functions. -/
instance : Scalar Float where
  max := fun a b => if a ≤ b then b else a
  min := fun a b => if a ≤ b then a else b
  ofNat := floatOfNat
  sqrt := Float.sqrt
  atan := Float.atan
  pi := floatPi
  decLt := fun a b => Float.decLt a b
  decLe := fun a b => Float.decLe a b

end Ddg
