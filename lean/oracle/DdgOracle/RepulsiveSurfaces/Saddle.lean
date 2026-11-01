import DdgOracle.RepulsiveSurfaces.Constraints

/-!
# The nested saddle and the preconditioner (SPEC §6.1, §6.2, §6.3)

* `Ap`: the dense matrix `[[A₃, C_Qᵀ], [C_Q, 0]]` (positional rows inside), factored by LU.
* `Nested`: the Schur rows outside, `x = y + Z S⁻¹ C_N y`, `S = −C_N Z`, `Z = A′⁻¹ C_Nᵀ`.
* The correction is `h = −Z S⁻¹ Φ_N`, so that `C_N h = Φ_N` (SPEC §6.2 prints `h = A′⁻¹ C_Nᵀ S⁻¹ Φ`,
  which gives `C_N h = −Φ`; see NOTES-WP11.md).
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

/-- Dense row of a sparse Jacobian row. -/
def denseRow (n3 : Nat) (r : CRow) : FloatArray := Id.run do
  let mut x := zeros n3
  for (j, v) in r.jac do x := x.set! j (x.get! j + v)
  return x

/-- `[[A⊗I₃, C_Qᵀ],[C_Q, 0]]`, its LU and the size data. -/
structure Ap where
  /-- `3 |V|` -/
  n3 : Nat
  /-- positional rows -/
  kq : Nat
  /-- positional Jacobian, `kq × n3` -/
  cq : Mat
  /-- the factorization -/
  lu : LU
  /-- smallest over largest pivot of the equilibrated matrix -/
  pivot : Float

/-- Assemble and factor `A′`. `SingularSaddle` when the equilibrated pivot ratio is below `1e-13`. -/
def buildAp (nv : Nat) (a : Mat) (q : Array CRow) : Except String Ap := do
  let n3 := 3 * nv
  let kq := q.size
  let n := n3 + kq
  let mut m := zeros (n * n)
  for v in [0:nv] do
    for w in [0:nv] do
      let x := a.get! (v * nv + w)
      for d in [0:3] do
        m := m.set! ((3 * v + d) * n + 3 * w + d) x
  let mut cq := zeros (kq * n3)
  for r in [0:kq] do
    let row := denseRow n3 q[r]!
    for j in [0:n3] do
      cq := cq.set! (r * n3 + j) (row.get! j)
      m := m.set! ((n3 + r) * n + j) (row.get! j)
      m := m.set! (j * n + n3 + r) (row.get! j)
  let some lu := luFactor n m | throw "SingularSaddle: the saddle matrix is singular"
  -- equilibrate: primal unknowns by 1/sqrt(A_jj), constraint rows to unit norm
  let mut dd := Array.replicate n 1.0
  for j in [0:n3] do
    let x := m.get! (j * n + j)
    dd := dd.set! j (if x > 0.0 then 1.0 / Float.sqrt x else 1.0)
  for r in [0:kq] do
    let mut s := 0.0
    for j in [0:n3] do
      let x := m.get! ((n3 + r) * n + j) * dd[j]!
      s := s + x * x
    dd := dd.set! (n3 + r) (if s > 0.0 then 1.0 / Float.sqrt s else 1.0)
  let mut me := zeros (n * n)
  for i in [0:n] do
    for j in [0:n] do
      me := me.set! (i * n + j) (dd[i]! * m.get! (i * n + j) * dd[j]!)
  let some lue := luFactor n me | throw "SingularSaddle: the equilibrated saddle matrix is singular"
  let mut lo := 1.0e300
  let mut hi := 0.0
  for i in [0:n] do
    let x := (lue.a.get! (i * n + i)).abs
    if x < lo then lo := x
    if x > hi then hi := x
  let ratio := if n == 0 then 1.0 else lo / hi
  if !(ratio > 1.0e-13) then throw s!"SingularSaddle: pivot ratio {ratio} (redundant rows or a singular metric on the kernel of C)"
  return { n3, kq, cq, lu, pivot := ratio }

/-- Solve `A′ [y; λ] = [b; 0]`; returns `(y, λ)`. -/
def Ap.solve (ap : Ap) (b : FloatArray) : FloatArray × FloatArray :=
  let rhs := FloatArray.ofFn' (ap.n3 + ap.kq) fun i => if i < ap.n3 then b.get! i else 0.0
  let x := ap.lu.solve rhs
  (FloatArray.ofFn' ap.n3 fun i => x.get! i, FloatArray.ofFn' ap.kq fun i => x.get! (ap.n3 + i))

/-- Solve `A′ x = rhs` for a full right-hand side. -/
def Ap.solveFull (ap : Ap) (rhs : FloatArray) : FloatArray := ap.lu.solve rhs

/-- Dot product. -/
def dotF (a b : FloatArray) : Float := Id.run do
  let mut s := 0.0
  for i in [0:a.size] do s := s + a.get! i * b.get! i
  return s

/-- The nested saddle: `A′` plus the Schur rows `C_N`. -/
structure Nested where
  /-- the inner system -/
  ap : Ap
  /-- number of Schur rows -/
  kn : Nat
  /-- `C_N`, `kn × n3` -/
  cn : Mat
  /-- `Z = A′⁻¹ C_Nᵀ`, position block, one array per row -/
  z : Array FloatArray
  /-- `S = −C_N Z` -/
  s : Mat
  /-- LU of `S` (none when `kn = 0`) -/
  sLu : Option LU

/-- Build the nested saddle from the rows `(kind, positional, row)`. -/
def buildNested (nv : Nat) (a : Mat) (rows : Array (String × Bool × CRow)) : Except String Nested := do
  let q := (rows.filter (·.2.1)).map (·.2.2)
  let nrows := (rows.filter (fun r => !r.2.1)).map (·.2.2)
  let ap ← buildAp nv a q
  let kn := nrows.size
  let n3 := ap.n3
  let mut cn := zeros (kn * n3)
  let mut z : Array FloatArray := #[]
  for r in [0:kn] do
    let row := denseRow n3 nrows[r]!
    for j in [0:n3] do cn := cn.set! (r * n3 + j) (row.get! j)
    z := z.push (ap.solve row).1
  let mut s := zeros (kn * kn)
  for r in [0:kn] do
    for p in [0:kn] do
      let row := FloatArray.ofFn' n3 fun j => cn.get! (r * n3 + j)
      s := s.set! (r * kn + p) (-(dotF row z[p]!))
  let sLu ← if kn == 0 then pure none else
    match luFactor kn s with
    | some l => pure (some l)
    | none => throw "SingularSaddle: the Schur matrix S is singular"
  return { ap, kn, cn, z, s, sLu }

/-- `C_N v`. -/
def Nested.applyCn (nst : Nested) (v : FloatArray) : FloatArray :=
  FloatArray.ofFn' nst.kn fun r => Id.run do
    let mut s := 0.0
    for j in [0:nst.ap.n3] do s := s + nst.cn.get! (r * nst.ap.n3 + j) * v.get! j
    return s

/-- Combination `Σ_q c_q Z_q`. -/
private def Nested.combine (nst : Nested) (c : FloatArray) : FloatArray := Id.run do
  let mut x := zeros nst.ap.n3
  for q in [0:nst.kn] do
    for j in [0:nst.ap.n3] do
      x := x.set! j (x.get! j + c.get! q * nst.z[q]!.get! j)
  return x

/-- The direction `x = y + Z S⁻¹ C_N y`, `y = A′⁻¹ g` (SPEC §6.1). -/
def Nested.direction (nst : Nested) (g : FloatArray) : FloatArray :=
  let y := (nst.ap.solve g).1
  match nst.sLu with
  | none => y
  | some l =>
    let c := l.solve (nst.applyCn y)
    let zc := nst.combine c
    FloatArray.ofFn' nst.ap.n3 fun j => y.get! j + zc.get! j

/-- The Newton correction `h = −Z S⁻¹ Φ_N`, with `C_N h = Φ_N` (SPEC §6.2); the step is `f ← f − h`. -/
def Nested.correction (nst : Nested) (phi : FloatArray) : FloatArray :=
  match nst.sLu with
  | none => zeros nst.ap.n3
  | some l =>
    let c := l.solve phi
    let zc := nst.combine c
    FloatArray.ofFn' nst.ap.n3 fun j => -(zc.get! j)

/-- Positions as a flat vector `3 v + axis`. -/
def flatten (a : Array V3) : FloatArray := Id.run do
  let mut x := FloatArray.emptyWithCapacity (3 * a.size)
  for p in a do
    x := x.push p.x
    x := x.push p.y
    x := x.push p.z
  return x

/-- Inverse of `flatten`. -/
def unflatten (x : FloatArray) : Array V3 :=
  (Array.range (x.size / 3)).map fun v => ⟨x.get! (3 * v), x.get! (3 * v + 1), x.get! (3 * v + 2)⟩

/-! ## Preconditioner (SPEC §6.3) -/

/-- Result of the `precond` task. -/
structure Precond where
  /-- `ε = ε_rel tr W / tr M` -/
  eps : Float
  /-- `L_{2−s}` (`|V| × |V|`) -/
  lfr : Mat
  /-- the factored `L′` -/
  ap : Ap
  /-- the result `P⁻¹ r` -/
  z : FloatArray

/-- `ε = ε_rel · tr W / tr M` and `L′ = [[(W + εM) ⊗ I₃, C_Qᵀ],[C_Q, 0]]`. -/
def buildPrecondMatrix (nv : Nat) (w : Mat) (m : Array Float) (epsRel : Float) : Float × Mat := Id.run do
  let mut trw := 0.0
  let mut trm := 0.0
  for i in [0:nv] do
    trw := trw + w.get! (i * nv + i)
    trm := trm + m[i]!
  let eps := epsRel * trw / trm
  let mut a := w
  for i in [0:nv] do
    a := a.set! (i * nv + i) (a.get! (i * nv + i) + eps * m[i]!)
  return (eps, a)

/-- `P⁻¹ r = L′⁻¹ blockdiag(L_{2−s} ⊗ I₃, 0) L′⁻¹ r`; the multiplier part is zeroed between the solves. -/
def applyPrecond (nv : Nat) (ap : Ap) (lfr : Mat) (r : FloatArray) : FloatArray :=
  let p1 := (ap.solve r).1
  let v := FloatArray.ofFn' (3 * nv) fun idx => Id.run do
    let i := idx / 3
    let d := idx % 3
    let mut s := 0.0
    for j in [0:nv] do s := s + lfr.get! (i * nv + j) * p1.get! (3 * j + d)
    return s
  (ap.solve v).1

end Ddg.RepulsiveSurfaces
