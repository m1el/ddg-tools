import DdgOracle.Dense

/-!
# Dense symmetric matrices for the cut flow

A square matrix in a `FloatArray`, triangle-wise assembly of cotan stiffness matrices and the
solution of a principal subsystem by LU.
-/

namespace Ddg.Cuts

/-- Row-major square matrix. -/
structure Mat where
  /-- dimension -/
  n : Nat
  /-- entries -/
  d : FloatArray
  deriving Inhabited

namespace Mat

/-- The zero matrix. -/
def zero (n : Nat) : Mat := ⟨n, zeros (n * n)⟩

/-- Entry `(i, j)`. -/
def get (m : Mat) (i j : Nat) : Float := m.d.get! (i * m.n + j)

/-- Add `x` to the entry `(i, j)`. -/
def addAt (m : Mat) (i j : Nat) (x : Float) : Mat :=
  ⟨m.n, m.d.set! (i * m.n + j) (m.d.get! (i * m.n + j) + x)⟩

/-- Add `w` times the cotan stiffness matrix of one triangle with nodes `n0 n1 n2` and corner
cotangents `c0 c1 c2` (corner `k` at node `nk`): the off-diagonal entry of the edge opposite
corner `k` is `−½ w c_k`, rows sum to zero. -/
def addTri (m : Mat) (n0 n1 n2 : Nat) (c0 c1 c2 w : Float) : Mat :=
  let h := w / 2
  let m := m.addAt n1 n2 (-(h * c0))
  let m := m.addAt n2 n1 (-(h * c0))
  let m := m.addAt n2 n0 (-(h * c1))
  let m := m.addAt n0 n2 (-(h * c1))
  let m := m.addAt n0 n1 (-(h * c2))
  let m := m.addAt n1 n0 (-(h * c2))
  let m := m.addAt n0 n0 (h * (c1 + c2))
  let m := m.addAt n1 n1 (h * (c0 + c2))
  m.addAt n2 n2 (h * (c0 + c1))

/-- `M x`. -/
def mulVec (m : Mat) (x : FloatArray) : FloatArray := Id.run do
  let mut out := zeros m.n
  for i in [0:m.n] do
    let mut s := fzero
    for j in [0:m.n] do
      s := s + m.get i j * x.get! j
    out := out.set! i s
  return out

/-- Solve the principal subsystem `M[idx, idx] y = rhs` (`rhs` indexed like `idx`).
`none` if it is singular. -/
def solveSub (m : Mat) (idx : Array Nat) (rhs : FloatArray) : Option FloatArray := Id.run do
  let k := idx.size
  if k == 0 then return some (zeros 0)
  let mut a := zeros (k * k)
  for r in [0:k] do
    for c in [0:k] do
      a := a.set! (r * k + c) (m.get idx[r]! idx[c]!)
  match luFactor k a with
  | none => return none
  | some lu => return some (lu.solve rhs)

end Mat
end Ddg.Cuts
