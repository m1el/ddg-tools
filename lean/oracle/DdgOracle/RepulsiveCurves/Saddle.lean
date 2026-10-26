import DdgOracle.RepulsiveCurves.Constraints

/-!
# One projected Sobolev gradient and one chord projection (SPEC §6.1, §6.2)

`[formula: rc-saddle-gradient]`, `[formula: rc-descent-direction]`, `[formula: rc-chord-projection]`.
The saddle matrix `[Ā Cᵀ; C 0]` is assembled densely, `Ā = A ⊗ I` on the active coordinates, and
factored once by the LU of `Dense.lean`.
-/

namespace Ddg.RepulsiveCurves

/-- The factored saddle system at a base curve. -/
structure Saddle where
  /-- number of unknowns `d|V|` -/
  nd : Nat
  /-- number of constraint rows -/
  k : Nat
  /-- active coordinates -/
  coords : Array Nat
  /-- the Jacobian, `k × nd` -/
  c : FloatArray
  /-- LU of the saddle matrix -/
  lu : LU

/-- Assemble and factor `[Ā Cᵀ; C 0]`. -/
def Saddle.build (nv : Nat) (coords : Array Nat) (a : FloatArray) (rows : Array Row) : Except String Saddle := do
  let nc := coords.size
  let nd := nv * nc
  let k := rows.size
  let c := denseJacobian nv coords rows
  let n := nd + k
  let mut m := zeros (n * n)
  for v in [0:nv] do
    for w in [0:nv] do
      let x := a.get! (v * nv + w)
      for s in [0:nc] do
        m := m.set! ((v * nc + s) * n + (w * nc + s)) x
  for r in [0:k] do
    for j in [0:nd] do
      let x := c.get! (r * nd + j)
      m := m.set! ((nd + r) * n + j) x
      m := m.set! (j * n + nd + r) x
  match luFactor n m with
  | some lu =>
    -- a pivot below 1e-14 of the largest means a redundant constraint row (as in the Rust solver)
    let mut pmax : Float := 0.0
    let mut pmin : Float := 1e300
    for i in [0:n] do
      let p := (lu.a.get! (i * n + i)).abs
      if p > pmax then pmax := p
      if p < pmin then pmin := p
    if !(pmin > 1e-14 * pmax) then throw "SingularSaddle: a pivot of the saddle matrix is below 1e-14 of the largest"
    return { nd, k, coords, c, lu }
  | none => throw "SingularSaddle: the saddle matrix is singular"

/-- Solution `(x, λ)` of `[Ā Cᵀ; C 0] (x, λ) = (a, b)`. -/
def Saddle.solve (s : Saddle) (a b : FloatArray) : FloatArray × FloatArray := Id.run do
  let mut rhs := zeros (s.nd + s.k)
  for i in [0:s.nd] do rhs := rhs.set! i (a.get! i)
  for i in [0:s.k] do rhs := rhs.set! (s.nd + i) (b.get! i)
  let sol := s.lu.solve rhs
  let mut x := zeros s.nd
  let mut lam := zeros s.k
  for i in [0:s.nd] do x := x.set! i (sol.get! i)
  for i in [0:s.k] do lam := lam.set! i (sol.get! (s.nd + i))
  return (x, lam)

/-- Flatten a vertex field to the interleaved unknowns over `coords`. -/
def gather (coords : Array Nat) (f : Array V3) : FloatArray := Id.run do
  let mut x := FloatArray.emptyWithCapacity (f.size * coords.size)
  for p in f do
    for c in coords do x := x.push (p.get c)
  return x

/-- Back to a vertex field (0 along the held axis). -/
def scatter (coords : Array Nat) (x : FloatArray) : Array V3 := Id.run do
  let nc := coords.size
  let nv := x.size / nc
  let mut f : Array V3 := #[]
  for v in [0:nv] do
    let mut p := V3.zero
    for s in [0:nc] do p := p.set coords[s]! (x.get! (v * nc + s))
    f := f.push p
  return f

/-- Result of the projected gradient. -/
structure GradResult where
  /-- `g̃` per vertex -/
  g : Array V3
  /-- multipliers -/
  lam : FloatArray
  /-- `‖g̃‖₂` -/
  gnorm : Float
  /-- `dF · g̃` -/
  dfg : Float

/-- `[Ā Cᵀ; C 0] (g̃, λ) = (dFᵀ, 0)` (SPEC §6.1). -/
def projectedGradient (s : Saddle) (dE : Array V3) : GradResult :=
  let a := gather s.coords dE
  let (x, lam) := s.solve a (zeros s.k)
  let gnorm := Float.sqrt (x.foldl (fun acc v => acc + v * v) 0.0)
  let dfg := Id.run do
    let mut t : Float := 0.0
    for i in [0:s.nd] do t := t + a.get! i * x.get! i
    return t
  { g := scatter s.coords x, lam, gnorm, dfg }

/-- How a projection ended (SPEC §6.2). -/
inductive ProjStatus where
  /-- all rows within tolerance -/
  | converged
  /-- iteration limit -/
  | maxIterations
  /-- the largest relative violation grew by more than 10× -/
  | diverged
  /-- a tangent site ended with `s · X ≤ 0` -/
  | tangentFlipped (cons : Nat)
  deriving Inhabited

/-- Name in the text format. -/
def ProjStatus.name : ProjStatus → String
  | .converged => "converged"
  | .maxIterations => "max-iterations"
  | .diverged => "diverged"
  | .tangentFlipped c => s!"tangent-flipped {c}"

/-- Result of the chord projection. -/
structure ProjResult where
  /-- final iterate -/
  pos : Array V3
  /-- how it ended -/
  status : ProjStatus
  /-- chord iterations (solves) -/
  iterations : Nat
  /-- `max_r |Φ_r|` at the end -/
  maxViolation : Float
  /-- `max_r |Φ_r|/tol_r` per iterate, starting with the trial point -/
  history : Array Float
  deriving Inhabited

/-- The chord projection of SPEC §6.2 from `trial`, with `Ā`, `C` frozen in `s`. `tau` is the
length tolerance `τ_len` (already multiplied by `ℓ̄₀` where it is relative); `lbar` is `ℓ̄₀`. The row
tolerances are `τ_len` (length), `τ_len/ℓ̄₀` (tangent), `τ_len |∇f(γ_v)|` (surface, at the
current iterate) (§5.1). -/
def chordProject (net : Net) (rc : Array RCons) (s : Saddle) (trial : Array V3) (tau lbar : Float)
    (maxit : Nat) : ProjResult := Id.run do
  let mut pos := trial
  let mut hist : Array Float := #[]
  let mut it := 0
  let zero := zeros s.nd
  let mut result : Option ProjResult := none
  while result.isNone do
    let rows := evalRows net pos s.coords rc
    let mut rel : Float := 0.0
    let mut maxv : Float := 0.0
    for r in rows do
      let tol := if r.unit == 0 then tau else if r.unit == 1 then tau / lbar else tau * r.gscale
      let q := r.phi.abs / tol
      if q.isNaN || rel.isNaN then rel := Float.ofBits 0x7ff8000000000000
      else if q > rel then rel := q
      if r.phi.abs > maxv then maxv := r.phi.abs
    hist := hist.push rel
    let fin := fun (st : ProjStatus) => some { pos, status := st, iterations := it, maxViolation := maxv, history := hist : ProjResult }
    if rel ≤ 1.0 then
      -- orientation of the tangent sites: s · X > 0
      let mut flipped : Option Nat := none
      for ci in [0:rc.size] do
        match rc[ci]! with
        | .tangent site x _ _ =>
          match tangentOf net pos site with
          | .ok (t, _) => if flipped.isNone && !(t.dot x > 0.0) then flipped := some ci
          | .error _ => pure ()
        | _ => pure ()
      result := match flipped with
        | some ci => fin (.tangentFlipped ci)
        | none => fin .converged
    else if !rel.isFinite || (it > 0 && rel > 10.0 * hist[0]!) then
      result := fin .diverged
    else if it ≥ maxit then
      result := fin .maxIterations
    else
      let mut b := zeros s.k
      for i in [0:s.k] do b := b.set! i (-(rows[i]!.phi))
      let (x, _) := s.solve zero b
      let dx := scatter s.coords x
      pos := (Array.range pos.size).map fun v => pos[v]! + dx[v]!
      it := it + 1
  return result.get!

end Ddg.RepulsiveCurves
