import DdgOracle.Solver

/-!
# Optimality certificate

Given any candidate `(u, φ, µ)` on the vertices of a mesh, check the optimality conditions of the
canonical problem directly. Nothing here depends on how the candidate was computed:

* the Laplacian is applied face by face from the corner cotangents, not through the dense matrix
  the Newton solver uses;
* the inclusion is checked from the original per-vertex parameters, not from the solver's
  active sets.

Since the problem is convex, a candidate with zero residuals is a global minimizer; the duality
gap is a second, independent witness of the same fact.
-/

namespace Ddg

/-- Residuals of the optimality conditions and the quantities derived from a candidate. -/
structure Certificate where
  /-- `½ Σ wE_i m_i u_i² + λ Σ wR_i |µ_i|` -/
  objective : Float
  /-- `½ Σ wE_i m_i u_i²` -/
  distortion : Float
  /-- `λ Σ wR_i |µ_i|` over the vertices that are neither forbidden nor prescribed -/
  coneCost : Float
  /-- `max_i |(L u)_i − µ_i + Ω_i|` over interior vertices -/
  resPrimal : Float
  /-- `max_i |(L φ)_i + wE_i m_i u_i|` over interior vertices -/
  resAdjoint : Float
  /-- `max_i dist(φ_i, ∂g_i(µ_i))`, with `µ_i` snapped to `0`, `lo_i`, `hi_i` within `snapTol` -/
  resInclusion : Float
  /-- `max_i dist(µ_i, [lo_i, hi_i])`; for forbidden and prescribed cones `|µ_i − value|` -/
  resBounds : Float
  /-- primal objective minus dual objective -/
  gap : Float
  /-- `Σ µ_i` over interior vertices -/
  sumMu : Float
  /-- `Σ Ω_i` over interior vertices (`2πχ` on a closed mesh) -/
  sumOmega : Float
  /-- closed mesh: `φ + c` satisfies the optimality conditions for `c` in
  `[shiftLo, shiftHi]`; both are `0` on a mesh with boundary -/
  shiftLo : Float
  /-- see `shiftLo` -/
  shiftHi : Float
  deriving Inhabited

/-- Threshold on the residuals and on the relative duality gap for the status `converged`. -/
def certificateTol : Float := 1e-9

/-- A cone angle within `snapTol` of `0` or of one of its bounds is treated as equal to it when
the inclusion is checked and when the interval `[shiftLo, shiftHi]` is computed. Without it the
verdict at a vertex with `µ_i = ±10⁻¹⁷` would depend on the sign of a rounding error. -/
def snapTol : Float := 1e-9

namespace Certificate

/-- All optimality conditions hold to `certificateTol`. -/
def ok (c : Certificate) : Bool :=
  c.resPrimal ≤ certificateTol && c.resAdjoint ≤ certificateTol &&
  c.resInclusion ≤ certificateTol && c.resBounds ≤ certificateTol &&
  c.gap.abs ≤ certificateTol * (if c.objective.abs > 1.0 then c.objective.abs else 1.0)

/-- The first condition that fails, for error messages. -/
def failure (c : Certificate) : String :=
  if !(c.resPrimal ≤ certificateTol) then "state equation residual too large"
  else if !(c.resAdjoint ≤ certificateTol) then "adjoint equation residual too large"
  else if !(c.resBounds ≤ certificateTol) then "cone angle bounds violated"
  else if !(c.resInclusion ≤ certificateTol) then "subdifferential inclusion violated"
  else "duality gap too large"

end Certificate

/-- `L v` for a vector on all vertices, accumulated face by face:
`(L v)_i = ½ Σ_{faces ijk} (cot θ_k (v_i − v_j) + cot θ_j (v_i − v_k))`. -/
def applyLaplacian (ops : Operators) (v : FloatArray) : FloatArray := Id.run do
  let mut out := zeros ops.nv
  for fi in [0:ops.topo.faces.size] do
    let f := ops.topo.faces[fi]!
    let g := ops.faceGeom[fi]!
    let wa := g.cotA / 2
    let wb := g.cotB / 2
    let wc := g.cotC / 2
    let va := v.get! f.a
    let vb := v.get! f.b
    let vc := v.get! f.c
    out := out.set! f.a (out.get! f.a + (wc * (va - vb) + wb * (va - vc)))
    out := out.set! f.b (out.get! f.b + (wc * (vb - va) + wa * (vb - vc)))
    out := out.set! f.c (out.get! f.c + (wb * (vc - va) + wa * (vc - vb)))
  return out

/-- Check the optimality conditions of the canonical problem for `(u, φ, µ)` given on all
vertices (zero at boundary vertices). -/
def certify (ops : Operators) (p : Params) (u phi mu : FloatArray) : Certificate := Id.run do
  let nv := ops.nv
  let closed := ops.topo.closed
  let mut prescribed : Array (Option Float) := Array.replicate nv none
  for (v, c) in p.cones do
    prescribed := prescribed.set! v (some c)
  let lu := applyLaplacian ops u
  let lphi := applyLaplacian ops phi
  let mut distortion := 0.0
  let mut coneCost := 0.0
  let mut dual := 0.0
  let mut resPrimal := 0.0
  let mut resAdjoint := 0.0
  let mut resInclusion := 0.0
  let mut resBounds := 0.0
  let mut sumMu := 0.0
  let mut sumOmega := 0.0
  let mut shiftLo := negInf
  let mut shiftHi := posInf
  let worse := fun (acc x : Float) => if x > acc || x.isNaN then x else acc
  for v in [0:nv] do
    if ops.topo.isBoundary[v]! then continue
    let w := p.wE.get! v * ops.mass.get! v
    let ui := u.get! v
    let fi := phi.get! v
    let mi := mu.get! v
    let om := ops.curvature.get! v
    distortion := distortion + w * (ui * ui) / 2
    dual := dual + om * fi
    sumMu := sumMu + mi
    sumOmega := sumOmega + om
    resPrimal := worse resPrimal (lu.get! v - mi + om).abs
    resAdjoint := worse resAdjoint (lphi.get! v + w * ui).abs
    match prescribed[v]! with
    | some c =>
      -- g = indicator of {c}: g*(φ) = c φ, no condition on φ
      resBounds := worse resBounds (mi - c).abs
      dual := dual - c * fi
    | none =>
      if p.wR.get! v == posInf then
        -- g = indicator of {0}: g* = 0, no condition on φ
        resBounds := worse resBounds mi.abs
      else
        let t := p.lambda * p.wR.get! v
        let lo := p.lo.get! v
        let hi := p.hi.get! v
        coneCost := coneCost + t * mi.abs
        resBounds := worse resBounds (boundViolation lo hi mi)
        resInclusion := worse resInclusion (subdiffViolationTol snapTol t lo hi mi fi)
        -- g*(φ) = sup over µ in [lo, hi] of (φ µ − t |µ|): attained at 0, lo or hi. An infinite
        -- bound contributes nothing here; exceeding ±t there is reported by `resInclusion`.
        let mut conj := 0.0
        if hi.isFinite && hi * (fi - t) > conj then conj := hi * (fi - t)
        if lo.isFinite && lo * (fi + t) > conj then conj := lo * (fi + t)
        dual := dual - conj
        -- interval of constants that may be added to φ
        let atLo := mi ≤ lo + snapTol
        let atHi := mi ≥ hi - snapTol
        if !atLo then shiftLo := worse shiftLo ((if mi > snapTol then t else -t) - fi)
        if !atHi then
          let b := (if mi < -snapTol then -t else t) - fi
          if b < shiftHi then shiftHi := b
  let objective := distortion + coneCost
  dual := dual - distortion
  if !closed then
    shiftLo := 0.0
    shiftHi := 0.0
  return {
    objective, distortion, coneCost, resPrimal, resAdjoint, resInclusion, resBounds
    gap := objective - dual, sumMu, sumOmega, shiftLo, shiftHi }

/-! ## The complete solve -/

/-- Outcome of a solve. -/
inductive Status where
  /-- the certificate holds -/
  | converged
  /-- no certified solution was found; the message says why -/
  | notConverged (reason : String)
  /-- closed mesh whose bounds are incompatible with Gauss–Bonnet -/
  | infeasible (reason : String)
  deriving Repr, Inhabited

/-- Result of a solve, on all vertices. -/
structure SolveResult where
  /-- outcome -/
  status : Status
  /-- number of interior vertices -/
  nInterior : Nat
  /-- mesh without boundary -/
  closed : Bool
  /-- log scale factor (0 at boundary vertices) -/
  u : FloatArray
  /-- adjoint (0 at boundary vertices) -/
  phi : FloatArray
  /-- cone angles (0 at boundary vertices) -/
  mu : FloatArray
  /-- Newton steps of each regularized problem -/
  itersOne : Array Nat
  /-- Newton steps of Stage II -/
  itersTwo : Nat
  /-- how Stage II stopped -/
  stop : Stop
  /-- the optimality certificate of `(u, φ, µ)` -/
  cert : Certificate
  deriving Inhabited

/-- Scatter a vector on the interior vertices to all vertices. -/
def scatter (s : Setup) (nv : Nat) (x : FloatArray) : FloatArray := Id.run do
  let mut out := zeros nv
  for i in [0:s.n] do
    out := out.set! s.interior[i]! (x.get! i)
  return out

/-- Solve the canonical problem by the two-stage Newton method and certify the result.
`Except.error` reports invalid parameters. -/
def solve (ops : Operators) (p : Params) (opts : Options := {}) : Except String SolveResult := do
  validateParams ops p
  let s := buildSetup ops p
  let nv := ops.nv
  if s.closed then
    let w := sumTo s.n fun i => s.wE.get! i * s.mass.get! i
    if !(w > 0) then
      throw "closed mesh with all wE = 0: the additive constant of u is undetermined"
    let sumLo := sumTo s.n fun i => s.lo.get! i
    let sumHi := sumTo s.n fun i => s.hi.get! i
    if sumLo > s.target + gaussBonnetTol || sumHi < s.target - gaussBonnetTol then
      let z := zeros nv
      return {
        status := .infeasible "the cone angle bounds are incompatible with Gauss-Bonnet (sum of mu must equal 2 pi chi)"
        nInterior := s.n, closed := true, u := z, phi := z, mu := z
        itersOne := #[], itersTwo := 0, stop := .infeasible, cert := certify ops p z z z }
  let r := newton s opts
  let u := scatter s nv r.u
  let phi := scatter s nv r.phi
  let mu := scatter s nv r.mu
  let cert := certify ops p u phi mu
  let status : Status :=
    if r.stop == .infeasible then
      .infeasible "the cone angle bounds are incompatible with Gauss-Bonnet (sum of mu must equal 2 pi chi)"
    else if cert.ok then .converged
    else
      let how := match r.stop with
        | .cap => "iteration cap reached"
        | .singular => "singular Newton system"
        | _ => "Newton iteration stopped"
      .notConverged s!"{how}; {cert.failure}"
  return {
    status, nInterior := s.n, closed := s.closed, u, phi, mu
    itersOne := r.itersOne, itersTwo := r.itersTwo, stop := r.stop, cert }

end Ddg
