import DdgOracle.Formulas
import DdgOracle.Operators
import DdgOracle.Dense
import DdgOracle.Hex

/-!
# The two-stage semismooth Newton solver (paper Fig. 12), dense

Stage I solves the regularized problems for `γ = 1, 10⁻¹, …, 10⁻¹⁰`, Stage II the exact problem
(`γ = 0`), each started from the previous solution. All of them are solved by the same routine,
`newtonStage`, as zeros of `Problem.Fimp γ` in the unknowns `(u, φ, µ)`; residuals and Newton
matrices are the scalar-generic definitions of `DdgOracle.Formulas.Optimality`, evaluated at
`Float`.

What is added to the paper, and why, is documented where it happens:

* `normalizeSite`: forbidden and prescribed cones as degenerate bounds;
* `nextState` (in `Formulas.Threshold`), the single-change mode of `newtonStage` and the
  regularized fallback in `newton`: safeguards against cycling and non-uniqueness;
* `balanceShift`, `newtonMatrix` (bordering) and `inclusionShift`: closed meshes;
* `snap` and the refinement steps: accuracy of the final iterate.

The solver never judges its own result: `DdgOracle.Certificate` does.
-/

namespace Ddg

/-- Parameters of the canonical problem, per vertex of the mesh (values at boundary vertices are
read but have no effect). -/
structure Params where
  /-- `λ ≥ 0` -/
  lambda : Float
  /-- distortion weights `wE_i ≥ 0` -/
  wE : FloatArray
  /-- cone weights `wR_i ∈ [0, +∞]` -/
  wR : FloatArray
  /-- lower bounds `lo_i ≤ 0` (possibly `−∞`) -/
  lo : FloatArray
  /-- upper bounds `hi_i ≥ 0` (possibly `+∞`) -/
  hi : FloatArray
  /-- prescribed cones: `(vertex, angle)` -/
  cones : Array (Nat × Float)
  deriving Inhabited

/-- Threshold and bounds `(t, lo, hi)` of one vertex after resolving the special cases:
a prescribed cone `c` is `(0, c, c)`, a forbidden cone (`wR = +∞`) is `(0, 0, 0)`, and
otherwise `t = λ·wR`. Both special cases have no cone cost and a constant proximal map. -/
def normalizeSite (lambda wR lo hi : Float) (prescribed : Option Float) : Float × Float × Float :=
  match prescribed with
  | some c => (0.0, c, c)
  | none => if wR == posInf then (0.0, 0.0, 0.0) else (lambda * wR, lo, hi)

/-- Check the parameters against the mesh. -/
def validateParams (ops : Operators) (p : Params) : Except String Unit := do
  let nv := ops.nv
  if !(p.lambda.isFinite && p.lambda ≥ 0) then throw "lambda must be finite and nonnegative"
  if p.wE.size != nv || p.wR.size != nv || p.lo.size != nv || p.hi.size != nv then
    throw s!"expected {nv} per-vertex parameter lines"
  for v in [0:nv] do
    let wE := p.wE.get! v
    let wR := p.wR.get! v
    let lo := p.lo.get! v
    let hi := p.hi.get! v
    if !(wE.isFinite && wE ≥ 0) then throw s!"vertex {v}: wE must be finite and nonnegative"
    if !(wR ≥ 0) then throw s!"vertex {v}: wR must be in [0, +inf]"
    if !(lo ≤ 0) then throw s!"vertex {v}: lo must be in [-inf, 0]"
    if !(hi ≥ 0) then throw s!"vertex {v}: hi must be in [0, +inf]"
  let mut seen : Array Bool := Array.replicate nv false
  for (v, c) in p.cones do
    if v ≥ nv then throw s!"prescribed cone at vertex {v}: index outside 0..{nv - 1}"
    if ops.topo.isBoundary[v]! then throw s!"prescribed cone at vertex {v}: boundary vertices carry no cone"
    if seen[v]! then throw s!"prescribed cone at vertex {v}: given twice"
    if !c.isFinite then throw s!"prescribed cone at vertex {v}: angle must be finite"
    seen := seen.set! v true

/-- The canonical problem on the interior vertices, in arrays. -/
structure Setup where
  /-- number of interior vertices -/
  n : Nat
  /-- no boundary: `L` is singular and Gauss–Bonnet constrains `Σ µ` -/
  closed : Bool
  /-- vertex index of each interior index -/
  interior : Array Nat
  /-- `L_II`, row-major `n × n` -/
  L : FloatArray
  /-- lumped masses -/
  mass : FloatArray
  /-- distortion weights -/
  wE : FloatArray
  /-- angle defects -/
  omega : FloatArray
  /-- thresholds `t_i` (see `normalizeSite`) -/
  t : FloatArray
  /-- lower bounds -/
  lo : FloatArray
  /-- upper bounds -/
  hi : FloatArray

namespace Setup

/-- The same data as a scalar-generic `Problem`. -/
def problem (s : Setup) : Problem Float where
  n := s.n
  L := fun i j => s.L.get! (i * s.n + j)
  mass := fun i => s.mass.get! i
  wE := fun i => s.wE.get! i
  omega := fun i => s.omega.get! i
  t := fun i => s.t.get! i
  lo := fun i => s.lo.get! i
  hi := fun i => s.hi.get! i

/-- `Σ Ω_i`; on a closed mesh this is `2πχ` and the value `Σ µ_i` must take. -/
def target (s : Setup) : Float := sumTo s.n fun i => s.omega.get! i

end Setup

/-- Restrict the operators and parameters to the interior vertices. -/
def buildSetup (ops : Operators) (p : Params) : Setup := Id.run do
  let nv := ops.nv
  let mut interior : Array Nat := #[]
  for v in [0:nv] do
    if !ops.topo.isBoundary[v]! then interior := interior.push v
  let n := interior.size
  let mut prescribed : Array (Option Float) := Array.replicate nv none
  for (v, c) in p.cones do
    prescribed := prescribed.set! v (some c)
  let mut L := FloatArray.emptyWithCapacity (n * n)
  for i in [0:n] do
    for j in [0:n] do
      L := L.push (ops.L interior[i]! interior[j]!)
  let mut t := FloatArray.emptyWithCapacity n
  let mut lo := FloatArray.emptyWithCapacity n
  let mut hi := FloatArray.emptyWithCapacity n
  for i in [0:n] do
    let v := interior[i]!
    let (ti, li, hi') := normalizeSite p.lambda (p.wR.get! v) (p.lo.get! v) (p.hi.get! v) prescribed[v]!
    t := t.push ti
    lo := lo.push li
    hi := hi.push hi'
  return {
    n, closed := ops.topo.closed, interior, L
    mass := FloatArray.ofFn' n fun i => ops.mass.get! interior[i]!
    wE := FloatArray.ofFn' n fun i => p.wE.get! interior[i]!
    omega := FloatArray.ofFn' n fun i => ops.curvature.get! interior[i]!
    t, lo, hi }

/-! ## Helpers -/

/-- A `FloatArray` as a function. -/
def vecFn (a : FloatArray) : Nat → Float := fun i => a.get! i

/-- Row-major `N × N` matrix with entries `entry r c`. -/
def assemble (N : Nat) (entry : Nat → Nat → Float) : FloatArray := Id.run do
  let mut a := FloatArray.emptyWithCapacity (N * N)
  for r in [0:N] do
    for c in [0:N] do
      a := a.push (entry r c)
  return a

/-- `[f 0, …, f (n-1)]`. -/
def arrayOf {β : Type} (n : Nat) (f : Nat → β) : Array β := Id.run do
  let mut a := Array.mkEmpty n
  for i in [0:n] do
    a := a.push (f i)
  return a

/-- No vertex is in `A⁺ ∪ A⁻`. -/
def noneActive (states : Array ConeState) : Bool := states.all fun s => !s.isActive

/-- `(|L| |v|)_i`, the magnitude of the terms of `(L v)_i`. -/
def absApply (s : Setup) (v : FloatArray) (i : Nat) : Float := Id.run do
  let mut acc := fzero
  for j in [0:s.n] do
    acc := acc + (s.L.get! (i * s.n + j)).abs * (v.get! j).abs
  return acc

/-- Relative size at which a residual counts as rounding noise. -/
def residualEps : Float := 1e-13

/-- Tolerance on `|Σ µ − 2πχ|` when deciding whether fixed cone angles satisfy Gauss–Bonnet. -/
def gaussBonnetTol : Float := 1e-9

/-! ## Closed meshes: the constant of `φ`

On a closed mesh `L` annihilates constants. Summing the state equation gives Gauss–Bonnet,
`Σ µ = Σ Ω = 2πχ`; summing the adjoint equation gives `Σ wE m u = 0`, which fixes the constant of
`u`; and `φ` is determined by the adjoint equation only up to a constant `c`, the multiplier of
the Gauss–Bonnet constraint. In the Newton system `c` is pinned by the rows `φ_i = ±t_i` of the
active vertices, so the system is regular exactly when the active set is nonempty.

When the active set is empty all `µ_i` sit on flat pieces of `P` and `c` cannot be found by
linearization. It is then found from the constraint it belongs to: `balanceShift` returns the
constant of smallest magnitude for which `Σ_i P_i(x_i + c) = 2πχ`. This is the exact proximal
step of `g + indicator{Σ µ = 2πχ}`, and for Stage I the exact line minimization of the dual
objective along the null direction of its Hessian.
-/

/-- Smallest-magnitude `c` with `Σ_i clamp(k_i · P_{t_i}(x_i + c)) = target`, or `none` if no
such `c` exists (the bounds are incompatible with Gauss–Bonnet). The left-hand side is a
nondecreasing piecewise linear function of `c`; it is evaluated at all of its breakpoints and
inverted exactly on the segment that brackets `target`. -/
def balanceShift (s : Setup) (k x : Nat → Float) (target : Float) : Option Float := Id.run do
  let n := s.n
  let t := vecFn s.t
  let lo := vecFn s.lo
  let hi := vecFn s.hi
  let h := fun (c : Float) =>
    sumTo n fun i => regularizedMeasure (k i) (t i) (lo i) (hi i) (x i + c)
  let mut breaks : Array Float := #[]
  let mut slopeLeft := 0.0
  let mut slopeRight := 0.0
  for i in [0:n] do
    if lo i < hi i then
      breaks := (breaks.push (-(t i) - x i)).push (t i - x i)
      if (hi i).isFinite then breaks := breaks.push (t i + hi i / k i - x i)
      else slopeRight := slopeRight + k i
      if (lo i).isFinite then breaks := breaks.push (-(t i) + lo i / k i - x i)
      else slopeLeft := slopeLeft + k i
  if breaks.size == 0 then
    return if (h 0.0 - target).abs ≤ gaussBonnetTol then some 0.0 else none
  let bs := breaks.qsort (· < ·)
  let hv := bs.map h
  let last := bs.size - 1
  let mut tgt := target
  if slopeLeft == 0.0 && tgt < hv[0]! then
    if hv[0]! - tgt ≤ gaussBonnetTol then tgt := hv[0]! else return none
  if slopeRight == 0.0 && tgt > hv[last]! then
    if tgt - hv[last]! ≤ gaussBonnetTol then tgt := hv[last]! else return none
  let mut c1 := 0.0
  let mut c2 := 0.0
  if tgt < hv[0]! then
    c1 := bs[0]! - (hv[0]! - tgt) / slopeLeft
    c2 := c1
  else if tgt > hv[last]! then
    c1 := bs[last]! + (tgt - hv[last]!) / slopeRight
    c2 := c1
  else
    -- leftmost root: first breakpoint with value ≥ target
    let mut k1 := 0
    for j in [0:bs.size] do
      if hv[j]! ≥ tgt then
        k1 := j
        break
    if k1 == 0 then c1 := if slopeLeft > 0.0 then bs[0]! else negInf
    else
      c1 := bs[k1-1]! + (tgt - hv[k1-1]!) / (hv[k1]! - hv[k1-1]!) * (bs[k1]! - bs[k1-1]!)
    -- rightmost root: last breakpoint with value ≤ target
    let mut k2 := last
    for jj in [0:bs.size] do
      let j := last - jj
      if hv[j]! ≤ tgt then
        k2 := j
        break
    if k2 == last then c2 := if slopeRight > 0.0 then bs[last]! else posInf
    else
      c2 := bs[k2]! + (tgt - hv[k2]!) / (hv[k2+1]! - hv[k2]!) * (bs[k2+1]! - bs[k2]!)
  return some (if 0.0 < c1 then c1 else if c2 < 0.0 then c2 else 0.0)

/-- With all `µ_i` on flat pieces (`µ_i ∈ {0, lo_i, hi_i}`), the constant `c` to add to `φ` so
that `φ_i + c ∈ ∂g_i(µ_i)` for all `i`: the admissible constants form an interval; return its
point nearest to `0`, or, if the interval is empty, its (reversed) midpoint, which minimizes the
largest violation. This is how complementarity fixes the constant of `φ` when no vertex is
active. -/
def inclusionShift (s : Setup) (phi mu : FloatArray) : Float := Id.run do
  let mut cLo := negInf
  let mut cHi := posInf
  for i in [0:s.n] do
    let t := s.t.get! i
    let m := mu.get! i
    let p := phi.get! i
    if s.lo.get! i < m then
      let a := (if 0 < m then t else -t) - p
      if a > cLo then cLo := a
    if m < s.hi.get! i then
      let b := (if m < 0 then -t else t) - p
      if b < cHi then cHi := b
  if cLo ≤ cHi then
    return if 0.0 < cLo then cLo else if cHi < 0.0 then cHi else 0.0
  else
    return (cLo + cHi) / 2

/-- Position of the natural index `x` of a `3n` vector (blocks `u, φ, µ`, resp. rows
`state, adjoint, cone`) in the ordering used for the LU factorization (`µ, u, φ`, resp.
`cone, state, adjoint`), and vice versa: `x ↦ 2n + x` for `x < n`, `x ↦ x − n` otherwise, read
from the factorization order to the natural one.

Putting the cone block first lets plain LU with partial pivoting eliminate it at no cost (each
of its columns has at most two nonzeros), so the factorization costs as much as that of a
`2n × 2n` matrix although the Newton system is assembled in full. -/
def blockPerm (n x : Nat) : Nat := if x < n then 2 * n + x else x - n

/-- The `3n × 3n` Newton matrix in factorization order, with one extra column and row if
`bordered`. The border removes the null space of the empty-active-set system on a closed mesh:
the new column puts a slack `α·m_i` into the state rows (so the system is solvable even if
Gauss–Bonnet holds only up to rounding; `α` is the violation), the new row requires
`Σ m_i q_i = 0` (the update of `φ` has zero mean). -/
def newtonMatrix (s : Setup) (entry : Nat → Nat → Float) (bordered : Bool) : FloatArray := Id.run do
  let n := s.n
  let N := 3 * n
  let M := if bordered then N + 1 else N
  let mut a := FloatArray.emptyWithCapacity (M * M)
  for r in [0:N] do
    let ro := blockPerm n r
    for c in [0:N] do
      a := a.push (entry ro (blockPerm n c))
    if bordered then a := a.push (if ro < n then s.mass.get! ro else fzero)
  if bordered then
    for c in [0:N] do
      let co := blockPerm n c
      a := a.push (if n ≤ co && co < 2 * n then s.mass.get! (co - n) else fzero)
    a := a.push fzero
  return a

/-! ## Newton iterations -/

/-- Why a Newton iteration stopped. -/
inductive Stop where
  /-- all residuals are at rounding level -/
  | residual
  /-- the active sets did not change (paper Sec. 5.3.3) -/
  | states
  /-- iteration cap reached -/
  | cap
  /-- the Newton matrix is singular -/
  | singular
  /-- closed mesh: the bounds are incompatible with Gauss–Bonnet -/
  | infeasible
  /-- the exact Newton system is singular (the minimizer is not unique); the result is the
  minimizer of a regularized problem with `γ ≤ 10⁻¹¹` -/
  | regularized
  deriving Repr, BEq, Inhabited

/-- Name of a stop reason in the text format. -/
def Stop.name : Stop → String
  | .residual => "residual"
  | .states => "active_sets"
  | .cap => "cap"
  | .singular => "singular"
  | .infeasible => "infeasible"
  | .regularized => "regularized"

/-- Iteration caps and the regularization sequence. -/
structure Options where
  /-- Stage I regularization parameters, in order -/
  gammas : Array Float := #[1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
  /-- Newton steps per regularized problem -/
  maxIterOne : Nat := 60
  /-- Newton steps for the exact problem -/
  maxIterTwo : Nat := 100
  /-- steps of the single-change fallback, per problem -/
  maxIterSingle : Nat := 400
  /-- regularization parameters tried when the exact Newton system is singular -/
  gammasSingular : Array Float := #[1e-11, 1e-12, 1e-13]
  /-- extra Newton steps of Stage II on unchanged active sets (iterative refinement, reusing
  the LU factorization) -/
  refineLimit : Nat := 2

/-- After a Newton step with the given states, the components that the cone rows of the Newton
system determine in closed form are assigned exactly instead of through the LU solve:
`µ_i ∈ {0, lo_i, hi_i}` off the active sets and, for the exact problem, `φ_i = ±t_i` on them.
The cone angles the solver reports are therefore exactly zero where there is no cone. -/
def snap (s : Setup) (exact : Bool) (states : Array ConeState) (phi mu : FloatArray) :
    FloatArray × FloatArray := Id.run do
  let mut phi := phi
  let mut mu := mu
  for i in [0:s.n] do
    match states[i]! with
    | .zero => mu := mu.set! i fzero
    | .lower => mu := mu.set! i (s.lo.get! i)
    | .upper => mu := mu.set! i (s.hi.get! i)
    | .pos => if exact then phi := phi.set! i (s.t.get! i)
    | .neg => if exact then phi := phi.set! i (-(s.t.get! i))
  return (phi, mu)

/-- Result of one Newton solve. -/
structure StageResult where
  /-- log scale factor -/
  u : FloatArray
  /-- adjoint -/
  phi : FloatArray
  /-- cone angles -/
  mu : FloatArray
  /-- the pieces used in the last Newton step -/
  sigma : Option (Array ConeState)
  /-- number of Newton steps -/
  iters : Nat
  /-- why the iteration stopped -/
  stop : Stop

/-- Semismooth Newton for `Fimp γ (u, φ, µ) = 0`; `γ = 0` is the exact system `F = 0`
(`ComputeOptimalMeasure` of Fig. 12), `γ > 0` a regularized problem
(`ComputeRegularizedMeasure`). `sigma0` are the pieces of the last Newton step that produced the
starting point, if any. With `single` only one vertex changes its piece per step (the one whose
piece is most wrong): slow, but it breaks the cycles the simultaneous update can run into.

Each step: determine the pieces `σ` (`nextState`), solve the linear system
`DFpiece γ σ · y = −Fpiece γ σ (u, φ, µ)`, add `y`. Stop when the residual of `Fimp γ` is at
rounding level, or when the pieces did not change (then, for `refineLimit` more steps, the same
factorization is reused for iterative refinement). -/
def newtonStage (s : Setup) (gamma : Float) (maxIter refineLimit : Nat)
    (u0 phi0 mu0 : FloatArray) (sigma0 : Option (Array ConeState)) (single : Bool := false) :
    StageResult := Id.run do
  let p := s.problem
  let n := s.n
  let exact := gamma == 0.0
  let scale := p.coneScale gamma
  let mut u := u0
  let mut phi := phi0
  let mut mu := mu0
  let mut sigmaPrev := sigma0
  let mut stepped := false
  let mut cache : Option LU := none
  let mut refined := 0
  let mut iters := 0
  for _ in [0:maxIter + 1] do
    -- the pieces the current point lies on, and the safeguarded choice
    let targetOf := fun (phi mu : FloatArray) =>
      arrayOf n fun i =>
        (p.impState gamma (vecFn phi) (vecFn mu) i).canon (s.t.get! i) (s.lo.get! i) (s.hi.get! i)
    let chooseFrom := fun (prev : Option (Array ConeState)) (target : Array ConeState)
        (phi mu : FloatArray) => Id.run do
      let all := arrayOf n fun i =>
        nextState (s.t.get! i) (s.lo.get! i) (s.hi.get! i) (prev.map fun a => a[i]!) target[i]!
      match single, prev with
      | true, some old =>
        -- single-change mode: only the vertex whose piece is most wrong moves
        let mut worst := fzero
        let mut which := n
        for i in [0:n] do
          if all[i]! != old[i]! then
            let x := phi.get! i + mu.get! i
            let k := scale i
            let t := s.t.get! i
            let lo := s.lo.get! i
            let hi := s.hi.get! i
            let d := (old[i]!.value k t lo hi x - regularizedMeasure k t lo hi x).abs
            if which == n || d > worst then
              worst := d
              which := i
        return if which == n then all else old.set! which all[which]!
      | _, _ => return all
    let mut target := targetOf phi mu
    let mut sigma := chooseFrom sigmaPrev target phi mu
    let mut shifted := false
    let mut bordered := false
    if s.closed && noneActive sigma then
      -- all µ_i are fixed by σ: they must satisfy Gauss–Bonnet
      let fixedSum := fun (sigma : Array ConeState) (phi mu : FloatArray) => sumTo n fun i =>
        sigma[i]!.value (scale i) (s.t.get! i) (s.lo.get! i) (s.hi.get! i) (phi.get! i + mu.get! i)
      if (fixedSum sigma phi mu - s.target).abs > gaussBonnetTol then
        match balanceShift s scale (fun i => phi.get! i + mu.get! i) s.target with
        | none => return { u, phi, mu, sigma := sigmaPrev, iters, stop := .infeasible }
        | some c =>
          phi := FloatArray.ofFn' n fun i => phi.get! i + c
          target := targetOf phi mu
          sigma := chooseFrom sigmaPrev target phi mu
          shifted := true
          if noneActive sigma && (fixedSum sigma phi mu - s.target).abs > gaussBonnetTol then
            -- the safeguard kept every vertex on a fixed piece: take the balanced pieces
            sigma := target
      bordered := noneActive sigma
    let residual :=
      if exact then FloatArray.ofFn' (3 * n) (p.F (vecFn u) (vecFn phi) (vecFn mu))
      else FloatArray.ofFn' (3 * n) (p.Fimp gamma (vecFn u) (vecFn phi) (vecFn mu))
    -- residual at rounding level?
    let mut small := true
    for i in [0:n] do
      let s1 := absApply s u i + (mu.get! i).abs + (s.omega.get! i).abs
      let s2 := absApply s phi i + (s.wE.get! i * s.mass.get! i * u.get! i).abs
      let s3 := (mu.get! i).abs + scale i * ((phi.get! i).abs + (mu.get! i).abs + s.t.get! i)
      if !((residual.get! i).abs ≤ residualEps * s1) then small := false
      if !((residual.get! (n + i)).abs ≤ residualEps * s2) then small := false
      if !((residual.get! (2 * n + i)).abs ≤ residualEps * s3) then small := false
    if small then return { u, phi, mu, sigma := sigmaPrev, iters, stop := .residual }
    if stepped && !shifted && sigmaPrev == some sigma then
      if refined ≥ refineLimit then return { u, phi, mu, sigma := sigmaPrev, iters, stop := .states }
      refined := refined + 1
    else
      refined := 0
      cache := none
    if iters ≥ maxIter then return { u, phi, mu, sigma := sigmaPrev, iters, stop := .cap }
    -- on the pieces the point lies on, the affine map is the residual itself
    let onPieces := sigma == target
    let rhsResidual :=
      if onPieces then residual
      else FloatArray.ofFn' (3 * n)
        (p.Fpiece gamma (fun i => sigma[i]!) (vecFn u) (vecFn phi) (vecFn mu))
    let lu? : Option LU := match cache with
      | some lu => some lu
      | none =>
        let entry :=
          if exact && onPieces then p.DF (vecFn phi) (vecFn mu)
          else p.DFpiece gamma (fun i => sigma[i]!)
        luFactor (if bordered then 3 * n + 1 else 3 * n) (newtonMatrix s entry bordered)
    match lu? with
    | none => return { u, phi, mu, sigma := sigmaPrev, iters, stop := .singular }
    | some lu =>
      cache := some lu
      let rhs := FloatArray.ofFn' lu.n fun r =>
        if r < 3 * n then -(rhsResidual.get! (blockPerm n r)) else fzero
      let y := lu.solve rhs
      -- natural index `x` sits at position `x + n` (u, φ) or `x − 2n` (µ)
      u := FloatArray.ofFn' n fun i => u.get! i + y.get! (n + i)
      phi := FloatArray.ofFn' n fun i => phi.get! i + y.get! (2 * n + i)
      mu := FloatArray.ofFn' n fun i => mu.get! i + y.get! i
      (phi, mu) := snap s exact sigma phi mu
      if bordered && exact then
        let c := inclusionShift s phi mu
        phi := FloatArray.ofFn' n fun i => phi.get! i + c
      sigmaPrev := some sigma
      stepped := true
      iters := iters + 1
  return { u, phi, mu, sigma := sigmaPrev, iters, stop := .cap }

/-- Raw result of the two stages, on the interior vertices. -/
structure NewtonResult where
  /-- log scale factor -/
  u : FloatArray
  /-- adjoint -/
  phi : FloatArray
  /-- cone angles -/
  mu : FloatArray
  /-- Newton steps of each regularized problem -/
  itersOne : Array Nat
  /-- Newton steps of Stage II -/
  itersTwo : Nat
  /-- how Stage II stopped -/
  stop : Stop

/-- The paper's two-stage method (Fig. 12): the regularized problems for the `γ` of
`opts.gammas`, each started from the solution of the previous one, then the exact problem
(`γ = 0`) started from the last regularized solution `(u, φ, µ)`. -/
def newton (s : Setup) (opts : Options := {}) : NewtonResult := Id.run do
  let n := s.n
  let mut u := zeros n
  let mut phi := zeros n
  let mut mu := zeros n
  let mut sigma : Option (Array ConeState) := none
  let mut itersOne : Array Nat := #[]
  for gamma in opts.gammas do
    let mut r := newtonStage s gamma opts.maxIterOne 0 u phi mu sigma
    if r.stop == .cap then
      let r' := newtonStage s gamma opts.maxIterSingle 0 u phi mu sigma true
      r := { r' with iters := r.iters + r'.iters }
    itersOne := itersOne.push r.iters
    if r.stop == .infeasible then
      return { u, phi, mu, itersOne, itersTwo := 0, stop := .infeasible }
    -- a regularized problem that did not converge still provides a starting point
    if r.stop != .singular then
      u := r.u
      phi := r.phi
      mu := r.mu
      sigma := r.sigma
  let mut r := newtonStage s 0.0 opts.maxIterTwo opts.refineLimit u phi mu sigma
  if r.stop == .cap then
    let r' := newtonStage s 0.0 opts.maxIterSingle opts.refineLimit u phi mu sigma true
    r := { r' with iters := r.iters + r'.iters }
  if r.stop == .singular then
    -- The exact system has no unique solution: the minimizer is not unique (this needs
    -- `wE_i = 0` somewhere). The regularized problems are strictly convex in `µ`; their
    -- minimizers converge to the minimizer of the exact problem with the smallest `Σ µ_i²/m_i`.
    let mut found := false
    for gamma in opts.gammasSingular do
      let rr := newtonStage s gamma opts.maxIterOne 0 u phi mu sigma
      if rr.stop == .residual || rr.stop == .states then
        u := rr.u
        phi := rr.phi
        mu := rr.mu
        sigma := rr.sigma
        found := true
    if found then
      return { u, phi, mu, itersOne, itersTwo := r.iters, stop := .regularized }
  -- Enforce the cone equation `µ = P(φ + µ)` exactly on the returned point: on the pieces it
  -- lies on, `µ_i ∈ {0, lo_i, hi_i}` or `φ_i = ±t_i`. This changes values by rounding errors
  -- only (the iteration stops when the residual is at rounding level) and removes them.
  let pieces := arrayOf n fun i =>
    (s.problem.state (vecFn r.phi) (vecFn r.mu) i).canon (s.t.get! i) (s.lo.get! i) (s.hi.get! i)
  let (phi2, mu2) :=
    if r.stop == .residual || r.stop == .states then snap s true pieces r.phi r.mu
    else (r.phi, r.mu)
  return { u := r.u, phi := phi2, mu := mu2, itersOne, itersTwo := r.iters, stop := r.stop }

end Ddg
