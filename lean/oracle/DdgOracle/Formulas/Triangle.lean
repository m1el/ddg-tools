import DdgOracle.Scalar

/-!
# Euclidean triangles from edge lengths

All triangle quantities are functions of the three edge lengths only (the algorithm is intrinsic,
paper Sec. 5.1.1). Throughout, `a b c` are the side lengths and the corner in question is the one
**opposite the first argument `a`**, so the corner sits between the sides `b` and `c`.
-/

namespace Ddg

variable {α : Type} [Scalar α]

/-- A point of `ℝ³`. -/
structure Point (α : Type) where
  /-- first coordinate -/
  x : α
  /-- second coordinate -/
  y : α
  /-- third coordinate -/
  z : α

/-- Euclidean distance `‖p − q‖`. -/
def dist (p q : Point α) : α :=
  Scalar.sqrt ((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z))

/-- Interior angle opposite the side `a` of a triangle with side lengths `a b c`.

This is the half-angle form of the law of cosines,
`tan²(θ/2) = (a−b+c)(a+b−c) / ((a+b+c)(b+c−a))`, which is equivalent to
`cos θ = (b² + c² − a²) / (2bc)` for a nondegenerate triangle and numerically stable for small
angles (unlike `arccos`). -/
def cornerAngle (a b c : α) : α :=
  2 * Scalar.atan (Scalar.sqrt (((a - b + c) * (a - c + b)) / ((a + b + c) * (b - a + c))))

/-- Area of the triangle with side lengths `a b c` (Heron's formula). -/
def triArea (a b c : α) : α :=
  Scalar.sqrt ((a + b + c) * (b - a + c) * (a - b + c) * (a - c + b)) / 4

/-- Cotangent of the interior angle opposite the side `a`:
`cot θ = cos θ / sin θ = (b² + c² − a²) / (4 · area)`. -/
def cornerCot (a b c : α) : α :=
  (b * b + c * c - a * a) / (4 * triArea a b c)

end Ddg
