import DdgOracle.Solver

/-!
# Brute-force cross-check for tiny problems

Enumerate, for every interior vertex, which equation holds there
(`µ_i = 0`, `µ_i = lo_i`, `µ_i = hi_i`, `φ_i = t_i` or `φ_i = −t_i`), solve the resulting linear
optimality system, and keep the assignments whose solution satisfies every optimality condition
of the canonical problem. By convexity all of them are the same minimizer.

This path shares no code with the Newton solver beyond the problem data: it assembles its own
system and solves it with Gauss–Jordan elimination with complete pivoting.
-/

namespace Ddg

/-- The equation imposed at one vertex. -/
inductive Choice where
  /-- `µ_i = value` -/
  | muFixed (value : Float)
  /-- `φ_i = value` -/
  | phiFixed (value : Float)
  deriving Inhabited

/-- The equations that can hold at a vertex with threshold `t` and bounds `[lo, hi]`. -/
def choicesOf (t lo hi : Float) : Array Choice := Id.run do
  if lo == hi then return #[.muFixed lo]
  let mut cs : Array Choice := #[]
  if lo ≤ 0 && 0 ≤ hi && t > 0 then cs := cs.push (.muFixed 0.0)
  if lo.isFinite && lo != 0 then cs := cs.push (.muFixed lo)
  if hi.isFinite && hi != 0 then cs := cs.push (.muFixed hi)
  cs := cs.push (.phiFixed t)
  if t > 0 then cs := cs.push (.phiFixed (-t))
  return cs

/-- A minimizer found by enumeration. -/
structure BruteResult where
  /-- log scale factor -/
  u : FloatArray
  /-- adjoint -/
  phi : FloatArray
  /-- cone angles -/
  mu : FloatArray
  /-- objective value -/
  objective : Float
  /-- number of assignments whose solution satisfies all optimality conditions -/
  numValid : Nat
  /-- largest difference in `u` or `µ` between any two valid assignments -/
  spread : Float

/-- Does `(φ, µ)` satisfy bounds and inclusion at every vertex, up to `tol`? -/
def bruteValid (s : Setup) (phi mu : FloatArray) (tol : Float) : Bool := Id.run do
  for i in [0:s.n] do
    let t := s.t.get! i
    let lo := s.lo.get! i
    let hi := s.hi.get! i
    let m := mu.get! i
    let p := phi.get! i
    if !(lo - tol ≤ m && m ≤ hi + tol) then return false
    let atLo := m ≤ lo + tol
    let atHi := m ≥ hi - tol
    if !atLo && !((if m > tol then t else -t) - tol ≤ p) then return false
    if !atHi && !(p ≤ (if m < -tol then -t else t) + tol) then return false
  return true

/-- Enumerate all assignments. `none` if no assignment is valid or there are more than `limit`
of them. -/
def bruteForce (s : Setup) (limit : Nat := 400000) (tol : Float := 1e-8) : Option BruteResult :=
  Id.run do
  let n := s.n
  let choices := arrayOf n fun i => choicesOf (s.t.get! i) (s.lo.get! i) (s.hi.get! i)
  let mut total := 1
  for c in choices do
    total := total * c.size
    if total > limit then return none
  let w := fun i => s.wE.get! i * s.mass.get! i
  let mut best : Option BruteResult := none
  for code in [0:total] do
    -- decode the assignment (mixed radix)
    let mut pick : Array Choice := Array.mkEmpty n
    let mut rest := code
    let mut anyPhi := false
    for i in [0:n] do
      let c := choices[i]!
      let ch := c[rest % c.size]!
      rest := rest / c.size
      pick := pick.push ch
      match ch with
      | .phiFixed _ => anyPhi := true
      | .muFixed _ => pure ()
    -- On a closed mesh with every µ fixed the system is singular: one extra unknown α (slack of
    -- Gauss–Bonnet, must vanish) and one extra equation (Σ m φ = 0) make it regular.
    let extra := if s.closed && !anyPhi then 1 else 0
    let N := 3 * n + extra
    let mut a := zeros (N * N)
    let mut b := zeros N
    for i in [0:n] do
      for j in [0:n] do
        let l := s.L.get! (i * n + j)
        a := a.set! (i * N + j) l
        a := a.set! ((n + i) * N + (n + j)) l
      a := a.set! (i * N + (2 * n + i)) (-1.0)
      b := b.set! i (-(s.omega.get! i))
      a := a.set! ((n + i) * N + i) (w i)
      match pick[i]! with
      | .muFixed v =>
        a := a.set! ((2 * n + i) * N + (2 * n + i)) 1.0
        b := b.set! (2 * n + i) v
      | .phiFixed v =>
        a := a.set! ((2 * n + i) * N + (n + i)) 1.0
        b := b.set! (2 * n + i) v
      if extra == 1 then
        a := a.set! (i * N + 3 * n) (s.mass.get! i)
        a := a.set! (3 * n * N + (n + i)) (s.mass.get! i)
    match gaussJordanSolve N a b with
    | none => continue
    | some x =>
      if extra == 1 && !((x.get! (3 * n)).abs ≤ tol) then continue
      let u := FloatArray.ofFn' n fun i => x.get! i
      let mut phi := FloatArray.ofFn' n fun i => x.get! (n + i)
      let mu := FloatArray.ofFn' n fun i => x.get! (2 * n + i)
      if extra == 1 then
        -- the constant of φ is free: move it into the admissible interval if there is one
        let mut cLo := negInf
        let mut cHi := posInf
        for i in [0:n] do
          let t := s.t.get! i
          let m := mu.get! i
          if !(m ≤ s.lo.get! i + tol) then
            let v := (if m > tol then t else -t) - phi.get! i
            if v > cLo then cLo := v
          if !(m ≥ s.hi.get! i - tol) then
            let v := (if m < -tol then -t else t) - phi.get! i
            if v < cHi then cHi := v
        let c := if 0.0 < cLo then cLo else if cHi < 0.0 then cHi else 0.0
        phi := FloatArray.ofFn' n fun i => phi.get! i + c
      if !bruteValid s phi mu tol then continue
      let objective := s.problem.objective (vecFn u) (vecFn mu)
      match best with
      | none => best := some { u, phi, mu, objective, numValid := 1, spread := 0.0 }
      | some r =>
        let mut d := r.spread
        for i in [0:n] do
          let du := (u.get! i - r.u.get! i).abs
          let dm := (mu.get! i - r.mu.get! i).abs
          if du > d then d := du
          if dm > d then d := dm
        best := some { r with numValid := r.numValid + 1, spread := d }
  return best

end Ddg
