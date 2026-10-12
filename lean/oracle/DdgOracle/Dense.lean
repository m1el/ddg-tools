import DdgOracle.Scalar

/-!
# Dense linear algebra

Row-major dense matrices in a `FloatArray`, LU factorization with partial pivoting, and an
independent Gauss–Jordan solver with complete pivoting used only by the brute-force cross-check.
-/

namespace Ddg

/-- A `FloatArray` of `n` zeros. -/
def zeros (n : Nat) : FloatArray := Id.run do
  let mut a := FloatArray.emptyWithCapacity n
  for _ in [0:n] do
    a := a.push fzero
  return a

/-- `FloatArray` with entries `f 0, …, f (n-1)`. -/
def FloatArray.ofFn' (n : Nat) (f : Nat → Float) : FloatArray := Id.run do
  let mut a := FloatArray.emptyWithCapacity n
  for i in [0:n] do
    a := a.push (f i)
  return a

/-- Largest absolute value of the entries. -/
def maxAbs (a : FloatArray) : Float := Id.run do
  let mut m := fzero
  for i in [0:a.size] do
    let v := (a.get! i).abs
    if v > m || v.isNaN then m := v
  return m

/-- LU factorization `P A = L U` of an `n × n` matrix. -/
structure LU where
  /-- dimension -/
  n : Nat
  /-- `L` (unit lower, below the diagonal) and `U` (upper), packed row-major -/
  a : FloatArray
  /-- row `i` of `P A` is row `perm[i]` of `A` -/
  perm : Array Nat

/-- LU factorization with partial pivoting. `none` if a pivot is zero or not finite. -/
def luFactor (n : Nat) (a0 : FloatArray) : Option LU := Id.run do
  let mut a := a0
  let mut perm : Array Nat := Array.range n
  for k in [0:n] do
    let mut p := k
    let mut best := (a.get! (k * n + k)).abs
    for i in [k+1:n] do
      let v := (a.get! (i * n + k)).abs
      if v > best then
        best := v
        p := i
    if !(best > fzero && best.isFinite) then return none
    if p != k then
      for j in [0:n] do
        let x := a.get! (k * n + j)
        a := a.set! (k * n + j) (a.get! (p * n + j))
        a := a.set! (p * n + j) x
      let pk := perm[k]!
      perm := (perm.set! k perm[p]!).set! p pk
    let piv := a.get! (k * n + k)
    let rk := k * n
    for i in [k+1:n] do
      let ri := i * n
      let f := a.get! (ri + k) / piv
      if f != fzero then
        a := a.set! (ri + k) f
        for j in [k+1:n] do
          a := a.set! (ri + j) (a.get! (ri + j) - f * a.get! (rk + j))
  return some { n, a, perm }

/-- Solve `A x = b` from the factorization. -/
def LU.solve (lu : LU) (b : FloatArray) : FloatArray := Id.run do
  let n := lu.n
  let mut x := FloatArray.emptyWithCapacity n
  for i in [0:n] do
    x := x.push (b.get! lu.perm[i]!)
  for i in [0:n] do
    let ri := i * n
    let mut s := x.get! i
    for j in [0:i] do
      s := s - lu.a.get! (ri + j) * x.get! j
    x := x.set! i s
  for ii in [0:n] do
    let i := n - 1 - ii
    let ri := i * n
    let mut s := x.get! i
    for j in [i+1:n] do
      s := s - lu.a.get! (ri + j) * x.get! j
    x := x.set! i (s / lu.a.get! (ri + i))
  return x

/-- Solve `A x = b` by Gauss–Jordan elimination with complete pivoting. `none` if a pivot is
smaller than `tol` times the largest entry of `A`. Shares no code with `luFactor`. -/
def gaussJordanSolve (n : Nat) (a0 b0 : FloatArray) (tol : Float := 1e-11) : Option FloatArray :=
  Id.run do
  let mut a := a0
  let mut b := b0
  let scale := maxAbs a0
  if !(scale > fzero && scale.isFinite) then return none
  -- `col[k]` is the column eliminated at step `k`, `row[k]` the pivot row
  let mut rowUsed : Array Bool := Array.replicate n false
  let mut colUsed : Array Bool := Array.replicate n false
  let mut pivRow : Array Nat := Array.replicate n 0
  let mut pivCol : Array Nat := Array.replicate n 0
  for k in [0:n] do
    let mut best := fzero
    let mut pr := n
    let mut pc := n
    for i in [0:n] do
      if !rowUsed[i]! then
        for j in [0:n] do
          if !colUsed[j]! then
            let v := (a.get! (i * n + j)).abs
            if v > best then
              best := v
              pr := i
              pc := j
    if !(best > tol * scale) || pr == n then return none
    rowUsed := rowUsed.set! pr true
    colUsed := colUsed.set! pc true
    pivRow := pivRow.set! k pr
    pivCol := pivCol.set! k pc
    let piv := a.get! (pr * n + pc)
    for i in [0:n] do
      if i != pr then
        let f := a.get! (i * n + pc) / piv
        if f != fzero then
          for j in [0:n] do
            a := a.set! (i * n + j) (a.get! (i * n + j) - f * a.get! (pr * n + j))
          b := b.set! i (b.get! i - f * b.get! pr)
  let mut x := zeros n
  for k in [0:n] do
    let pr := pivRow[k]!
    let pc := pivCol[k]!
    x := x.set! pc (b.get! pr / a.get! (pr * n + pc))
  return some x

end Ddg
