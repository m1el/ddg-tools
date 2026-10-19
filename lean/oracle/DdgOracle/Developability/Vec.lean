import DdgOracle.Hex

/-!
# Small vectors and a Jacobi eigensolver for the developability oracle

Plain `Float` code, written for this file only (the scalar-generic formulas of
`DdgOracle/Formulas/Developability.lean` are deliberately not used, so that the two
implementations are independent). The eigensolver is the cyclic Jacobi method of Golub and Van Loan,
*Matrix Computations*, Sec. 8.5, on a dense symmetric matrix of size 2 or 3.
-/

namespace Ddg.Developability

/-- A vector of `R³`. -/
structure V3 where
  /-- first coordinate -/
  x : Float
  /-- second coordinate -/
  y : Float
  /-- third coordinate -/
  z : Float
  deriving Inhabited

namespace V3

/-- zero vector -/
def zero : V3 := ⟨0, 0, 0⟩
/-- sum -/
def add (a b : V3) : V3 := ⟨a.x + b.x, a.y + b.y, a.z + b.z⟩
/-- difference -/
def sub (a b : V3) : V3 := ⟨a.x - b.x, a.y - b.y, a.z - b.z⟩
/-- negation -/
def neg (a : V3) : V3 := ⟨-a.x, -a.y, -a.z⟩
/-- scalar multiple -/
def smul (s : Float) (a : V3) : V3 := ⟨s * a.x, s * a.y, s * a.z⟩
/-- dot product -/
def dot (a b : V3) : Float := a.x * b.x + a.y * b.y + a.z * b.z
/-- cross product -/
def cross (a b : V3) : V3 :=
  ⟨a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x⟩
/-- squared norm -/
def norm2 (a : V3) : Float := dot a a
/-- Euclidean norm -/
def norm (a : V3) : Float := Float.sqrt (norm2 a)
/-- the vector divided by its norm -/
def unit (a : V3) : V3 := smul (1 / norm a) a
/-- `true` iff all coordinates are finite -/
def finite (a : V3) : Bool := a.x.isFinite && a.y.isFinite && a.z.isFinite
/-- largest absolute coordinate -/
def maxAbs (a : V3) : Float := max a.x.abs (max a.y.abs a.z.abs)

instance : Add V3 := ⟨add⟩
instance : Sub V3 := ⟨sub⟩
instance : Neg V3 := ⟨neg⟩
instance : HMul Float V3 V3 := ⟨smul⟩

/-- Coordinate `i ∈ {0,1,2}`. -/
def get (a : V3) (i : Nat) : Float := if i == 0 then a.x else if i == 1 then a.y else a.z

end V3

/-- Sum of an array of vectors, in index order. -/
def sumV (xs : Array V3) : V3 := xs.foldl (· + ·) V3.zero

/-- `x` with the sign of `s` (`+1` for `s = 0`). -/
def signOf (s : Float) : Float := if s < 0 then -1 else 1

/-- Eigen-decomposition of the symmetric `n × n` matrix `a` (row-major, `n ≤ 3` in practice) by
cyclic Jacobi sweeps. Returns the eigenvalues in ascending order and the matching unit eigenvectors
`vecs[k]` (as arrays of length `n`). -/
def symEig (n : Nat) (a0 : Array Float) : Array Float × Array (Array Float) := Id.run do
  let mut a := a0
  let mut v : Array Float := Array.ofFn (n := n * n) fun i => if i.val / n == i.val % n then 1.0 else 0.0
  for _ in [0:60] do
    let mut rotated := false
    for p in [0:n] do
      for q in [p + 1:n] do
        let apq := a[p * n + q]!
        let scale := Float.sqrt (a[p * n + p]!.abs * a[q * n + q]!.abs)
        if apq != 0 && apq.abs > 1e-300 && apq.abs > 1e-18 * scale then
          rotated := true
          let tau := (a[q * n + q]! - a[p * n + p]!) / (2 * apq)
          let t := signOf tau / (tau.abs + Float.sqrt (1 + tau * tau))
          let c := 1 / Float.sqrt (1 + t * t)
          let s := t * c
          for k in [0:n] do
            let akp := a[k * n + p]!
            let akq := a[k * n + q]!
            a := (a.set! (k * n + p) (c * akp - s * akq)).set! (k * n + q) (s * akp + c * akq)
          for k in [0:n] do
            let apk := a[p * n + k]!
            let aqk := a[q * n + k]!
            a := (a.set! (p * n + k) (c * apk - s * aqk)).set! (q * n + k) (s * apk + c * aqk)
          a := (a.set! (p * n + q) 0).set! (q * n + p) 0
          for k in [0:n] do
            let vkp := v[k * n + p]!
            let vkq := v[k * n + q]!
            v := (v.set! (k * n + p) (c * vkp - s * vkq)).set! (k * n + q) (s * vkp + c * vkq)
    if !rotated then break
  let idx := (Array.range n).qsort fun i j => a[i * n + i]! < a[j * n + j]!
  let vals := idx.map fun i => a[i * n + i]!
  let vecs := idx.map fun j => Array.ofFn (n := n) fun k => v[k.val * n + j]!
  return (vals, vecs)

end Ddg.Developability
