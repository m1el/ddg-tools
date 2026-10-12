import DdgOracle.Format
import DdgOracle.BruteForce

/-!
# `selftest`

Checks that do not need any external reference:

* the formulas against closed-form values, and the proximal identities on exact grids;
* Gauss–Bonnet and the structure of the Laplacian on every kind of generated mesh;
* the Newton solver against the brute-force enumeration, against a closed form, against random
  feasible perturbations (no perturbation may lower the objective), and against the paper's
  qualitative claims;
* invariance under scaling and relabeling;
* the text format, including recovery from malformed input.
-/

namespace Ddg

/-- Collects the outcome of the checks. -/
structure Tally where
  /-- number of checks that passed -/
  passed : Nat := 0
  /-- number of checks that failed -/
  failed : Nat := 0

/-- Record and print one check. -/
def check (tally : IO.Ref Tally) (name : String) (ok : Bool) (detail : String := "") : IO Unit := do
  if ok then
    tally.modify fun t => { t with passed := t.passed + 1 }
    IO.println s!"ok    {name}{if detail.isEmpty then "" else "  [" ++ detail ++ "]"}"
  else
    tally.modify fun t => { t with failed := t.failed + 1 }
    IO.println s!"FAIL  {name}{if detail.isEmpty then "" else "  [" ++ detail ++ "]"}"

/-- Short decimal rendering for diagnostics (never used for data exchange). -/
def sci (x : Float) : String :=
  if x == 0.0 then "0"
  else if x.isNaN then "NaN"
  else if x.isInf then (if x > 0 then "inf" else "-inf")
  else
    let e := (Float.log10 x.abs).floor
    let m := x / Float.exp (e * Float.log 10.0)
    let digits := (Float.round (m * 100.0)) / 100.0
    s!"{digits}e{e.toInt64}".replace "000000e" "e" |>.replace "0000e" "e"

/-- Largest absolute difference of two vectors. -/
def maxDiff (a b : FloatArray) : Float := Id.run do
  let mut d := 0.0
  for i in [0:a.size] do
    let x := (a.get! i - b.get! i).abs
    if x > d || x.isNaN then d := x
  return d

/-- Solve a case; `none` if the mesh or the parameters are invalid. -/
def solveCase (c : Case) : Option (Operators × SolveResult) :=
  match c.mesh.operators with
  | .error _ => none
  | .ok ops =>
    match solve ops c.params with
    | .error _ => none
    | .ok r => some (ops, r)

/-- `true` iff the status is `converged`. -/
def SolveResult.converged (r : SolveResult) : Bool :=
  match r.status with
  | .converged => true
  | _ => false

/-- Restriction of a vector on all vertices to the interior vertices. -/
def gather (s : Setup) (x : FloatArray) : FloatArray :=
  FloatArray.ofFn' s.n fun i => x.get! s.interior[i]!

/-- Objective value of the feasible point with cone angles `mu` (on the interior vertices):
solve the state equation for `u` by Gauss–Jordan elimination (on a closed mesh together with
the optimal constant, `Σ wE m u = 0`) and evaluate. Independent of the Newton solver. -/
def objectiveOf (s : Setup) (mu : FloatArray) : Option Float := Id.run do
  let n := s.n
  if n == 0 then return some 0.0
  let N := if s.closed then n + 1 else n
  let mut a := zeros (N * N)
  let mut b := zeros N
  for i in [0:n] do
    for j in [0:n] do
      a := a.set! (i * N + j) (s.L.get! (i * n + j))
    b := b.set! i (mu.get! i - s.omega.get! i)
    if s.closed then
      let w := s.wE.get! i * s.mass.get! i
      a := a.set! (i * N + n) w
      a := a.set! (n * N + i) w
  match gaussJordanSolve N a b 1e-14 with
  | none => return none
  | some x => return some (s.problem.objective (vecFn x) (vecFn mu))

/-- Try `trials` random feasible perturbations of the cone angles `mu` and return the largest
decrease of the objective found (`0` if none lowers it). -/
def bestPerturbation (s : Setup) (mu : FloatArray) (trials seed : Nat) : Float := Id.run do
  let n := s.n
  let some base := objectiveOf s mu | return posInf
  let mut rng : Rng := ⟨(seed * 31 + 5).toUInt64⟩
  let mut worst := 0.0
  let clampTo := fun (i : Nat) (x : Float) =>
    let lo := s.lo.get! i
    let hi := s.hi.get! i
    if x < lo then lo else if x > hi then hi else x
  for trial in [0:trials] do
    let eps := #[1e-1, 1e-2, 1e-4, 1e-6][trial % 4]!
    let mut m := mu
    if s.closed then
      -- move curvature between pairs of vertices: Σ µ stays fixed
      for _ in [0:1 + trial % 3] do
        let (i, r1) := rng.below n
        let (j, r2) := r1.below n
        let (d, r3) := r2.range (-eps) eps
        rng := r3
        if i != j then
          let di := clampTo i (m.get! i + d) - m.get! i
          let dj := -(clampTo j (m.get! j - di) - m.get! j)
          let d' := if di.abs < dj.abs then di else dj
          let mi := m.get! i + d'
          let mj := m.get! j - d'
          if mi == clampTo i mi && mj == clampTo j mj then
            m := (m.set! i mi).set! j mj
    else
      for i in [0:n] do
        let (d, r1) := rng.range (-eps) eps
        let (pick, r2) := r1.below 3
        rng := r2
        if pick == 0 || trial % 5 == 0 then m := m.set! i (clampTo i (m.get! i + d))
    match objectiveOf s m with
    | none => pure ()
    | some v => if base - v > worst then worst := base - v
  return worst

/-! ## The checks -/

private def checkFormulas (tally : IO.Ref Tally) : IO Unit := do
  -- hex round trip
  let samples : Array Float := #[0.0, -0.0, 1.0, -1.0, floatPi, 1e-310, 1.7976931348623157e308,
    posInf, negInf, 0.1, 1.0 / 3.0, -2.5e-17]
  let mut ok := true
  for x in samples do
    match hexToFloat? (floatToHex x) with
    | some y => if y.toBits != x.toBits then ok := false
    | none => ok := false
  check tally "hex: round trip is bit exact"
    (ok && floatToHex 1.0 == "3ff0000000000000" && floatToHex posInf == "7ff0000000000000"
      && floatToHex negInf == "fff0000000000000" && floatToHex (-0.0) == "8000000000000000"
      && hexToFloat? "3FF0000000000000" == some 1.0 && (hexToFloat? "3ff").isNone
      && (hexToFloat? "3ff000000000000g").isNone)
  -- triangles
  let third := floatPi / 3
  check tally "triangle: equilateral angle, cotangent, area"
    ((cornerAngle (1.0 : Float) 1.0 1.0 - third).abs ≤ 1e-15
      && (cornerCot (1.0 : Float) 1.0 1.0 - 1 / Float.sqrt 3).abs ≤ 1e-15
      && (triArea (1.0 : Float) 1.0 1.0 - Float.sqrt 3 / 4).abs ≤ 1e-15)
  check tally "triangle: 3-4-5 right triangle"
    ((cornerAngle (5.0 : Float) 3.0 4.0 - floatPi / 2).abs ≤ 1e-15
      && (cornerCot (5.0 : Float) 3.0 4.0).abs ≤ 1e-15
      && (cornerCot (3.0 : Float) 4.0 5.0 - 4.0 / 3.0).abs ≤ 1e-15
      && triArea (3.0 : Float) 4.0 5.0 == 6.0)
  let s := cornerAngle (2.0 : Float) 3.0 4.0 + cornerAngle (3.0 : Float) 4.0 2.0
    + cornerAngle (4.0 : Float) 2.0 3.0
  check tally "triangle: angles of a scalene triangle sum to pi" ((s - floatPi).abs ≤ 1e-15)
  check tally "triangle: small angle is accurate (needle with apex angle 1e-7)"
    ((cornerAngle (1e-7 : Float) 1.0 1.0 - 2 * Float.asin 0.5e-7).abs ≤ 1e-15)
  -- proximal identities on a grid of exactly representable numbers
  let grid : Array Float := #[-3.0, -2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
  let boxes : Array (Float × Float) := #[(negInf, posInf), (0.0, posInf), (negInf, 0.0),
    (-1.0, 1.0), (0.0, 1.0), (-0.5, 2.0), (0.0, 0.0), (1.0, 1.0), (-0.5, -0.5)]
  let mut fixedPoint := true
  let mut pieces := true
  let mut count := 0
  for t in #[0.0, 0.25, 0.5, 1.0] do
    for (lo, hi) in boxes do
      for mu in grid do
        for phi in grid do
          count := count + 1
          let isFixed := mu == boundedThreshold t lo hi (phi + mu)
          let inSub := subdiffViolation t lo hi mu phi == 0.0 && boundViolation lo hi mu == 0.0
          if isFixed != inSub then fixedPoint := false
      for x in grid do
        for k in #[1.0, 0.5, 4.0] do
          let st := coneState k t lo hi x
          if regularizedMeasure k t lo hi x != st.value k t lo hi x then pieces := false
          if (st.canon t lo hi).value k t lo hi x != st.value k t lo hi x then pieces := false
  check tally "prox: mu = P(phi + mu) iff phi in subdifferential (all sign and bound cases)"
    fixedPoint s!"{count} points"
  check tally "prox: clamp(k P_t(x)) equals the value of its piece, also after canonicalization" pieces
  -- regularized: µ = clamp((m/γ) P(φ)) iff µ = clamp(c P(φ + µ)), c = m/(m+γ)
  let mut implicit := true
  for t in #[0.0, 0.5] do
    for (lo, hi) in #[(negInf, posInf), (0.0, 1.0), (-1.0, 2.0)] do
      for phi in grid do
        -- m = 1, γ = 1: k = 1, c = 1/2
        let mu := regularizedMeasure 1.0 t lo hi phi
        if mu != regularizedMeasure 0.5 t lo hi (phi + mu) then implicit := false
  check tally "prox: eliminated and implicit regularized cone conditions agree" implicit

private def meshZoo : Array (String × RawMesh) := #[
  ("fan", fan 6 0.4),
  ("grid", gridDisk 6 5 fun x y => 0.3 * Float.sin (3 * x) * Float.cos (2 * y)),
  ("disk", (polar 3 9 0.0 1.0 fun x y => 0.5 * (1 - x * x - y * y)).jitter 0.02 3),
  ("annulus", (polar 3 10 0.4 1.0 fun x y => 0.2 * x * y).jitter 0.02 4),
  ("tetrahedron", tetrahedron),
  ("octahedron", bipyramid 4 1.0),
  ("cube", cube),
  ("sphere", (uvSphere 5 7 fun p => ⟨1.3 * p.x, p.y, 0.8 * p.z⟩).jitter 0.03 5),
  ("torus", (torus 7 5 1.0 0.4 fun a b => 0.2 * Float.sin (2 * a) * Float.cos b).jitter 0.01 6)]

private def checkOperators (tally : IO.Ref Tally) : IO Unit := do
  for (name, m) in meshZoo do
    match m.operators with
    | .error e => check tally s!"operators[{name}]: mesh accepted" false e
    | .ok ops =>
      let nv := ops.nv
      let total := sumTo nv fun v => ops.curvature.get! v
      let gb := (total - 2 * floatPi * Float.ofInt ops.topo.chi).abs
      check tally s!"operators[{name}]: Gauss-Bonnet, sum of curvature = 2 pi chi" (gb ≤ 1e-12)
        s!"chi={ops.topo.chi} error={sci gb}"
      let msum := sumTo nv fun v => ops.mass.get! v
      let mut sym := 0.0
      let mut rows := 0.0
      for i in [0:nv] do
        let mut r := 0.0
        for j in [0:nv] do
          r := r + ops.L i j
          let d := (ops.L i j - ops.L j i).abs
          if d > sym then sym := d
        if r.abs > rows then rows := r.abs
      -- dense matrix against the face-by-face application, and x'Lx >= 0
      let mut rng : Rng := ⟨17⟩
      let mut applyErr := 0.0
      let mut minQuad := posInf
      for _ in [0:8] do
        let mut x := FloatArray.emptyWithCapacity nv
        for _ in [0:nv] do
          let (v, r) := rng.range (-1.0) 1.0
          rng := r
          x := x.push v
        let y := applyLaplacian ops x
        let mut quad := 0.0
        for i in [0:nv] do
          let mut acc := 0.0
          for j in [0:nv] do
            acc := acc + ops.L i j * x.get! j
          if (acc - y.get! i).abs > applyErr then applyErr := (acc - y.get! i).abs
          quad := quad + x.get! i * acc
        if quad < minQuad then minQuad := quad
      check tally s!"operators[{name}]: L symmetric, zero row sums, PSD, matches face-wise apply; masses sum to 1"
        (sym == 0.0 && rows ≤ 1e-12 && applyErr ≤ 1e-12 && minQuad ≥ 0.0 && (msum - 1).abs ≤ 1e-14)
        s!"rows={sci rows} apply={sci applyErr} min x'Lx={sci minQuad}"
  -- known values
  if let .ok ops := tetrahedron.operators then
    check tally "operators: regular tetrahedron has defect pi and mass 1/4 at every vertex"
      ((List.range 4).all fun v =>
        (ops.curvature.get! v - floatPi).abs ≤ 1e-14 && (ops.mass.get! v - 0.25).abs ≤ 1e-15)
  if let .ok ops := cube.operators then
    check tally "operators: cube has defect pi/2 at every vertex"
      ((List.range 8).all fun v => (ops.curvature.get! v - floatPi / 2).abs ≤ 1e-14)
  if let .ok ops := (gridDisk 5 5 fun _ _ => 0.0).operators then
    let mut flat := true
    let mut turning := 0.0
    for v in [0:ops.nv] do
      if ops.topo.isBoundary[v]! then turning := turning + ops.curvature.get! v
      else if (ops.curvature.get! v).abs > 1e-14 then flat := false
    check tally "operators: flat square has zero defect inside and total boundary turning 2 pi"
      (flat && (turning - 2 * floatPi).abs ≤ 1e-13)
  -- rejected meshes
  let rejects : Array (String × Nat × Array Face × String) := #[
    ("edge in three faces", 5, #[⟨0, 1, 2⟩, ⟨1, 0, 3⟩, ⟨0, 1, 4⟩], "directed edge (0, 1)"),
    ("inconsistent orientation", 4, #[⟨0, 1, 2⟩, ⟨1, 3, 2⟩, ⟨0, 1, 3⟩], "directed edge (0, 1)"),
    ("index out of range", 3, #[⟨0, 1, 7⟩], "outside 0..2"),
    ("repeated vertex", 3, #[⟨0, 1, 1⟩], "repeated vertex"),
    ("unused vertex", 4, #[⟨0, 1, 2⟩], "vertex 3 is not used"),
    ("two components", 6, #[⟨0, 1, 2⟩, ⟨3, 4, 5⟩], "not connected"),
    ("bow tie vertex", 5, #[⟨0, 1, 2⟩, ⟨0, 3, 4⟩], "vertex 0")]
  for (what, nv, faces, needle) in rejects do
    match buildTopology nv faces with
    | .ok _ => check tally s!"operators: rejects {what}" false "accepted"
    | .error e => check tally s!"operators: rejects {what}" ((e.splitOn needle).length > 1) e
  let flatTri : RawMesh := { pos := #[⟨0, 0, 0⟩, ⟨1, 0, 0⟩, ⟨2, 0, 0⟩, ⟨0, 1, 0⟩],
                             faces := #[⟨0, 1, 2⟩, ⟨0, 2, 3⟩] }
  match flatTri.operators with
  | .ok _ => check tally "operators: rejects a degenerate triangle" false "accepted"
  | .error e => check tally "operators: rejects a degenerate triangle" true e

/-- Compare the Newton result of a case with the brute-force enumeration. -/
private def checkBrute (tally : IO.Ref Tally) (c : Case) : IO Unit := do
  match c.mesh.operators with
  | .error e => check tally s!"brute[{c.name}]" false e
  | .ok ops =>
    match solve ops c.params with
    | .error e => check tally s!"brute[{c.name}]" false e
    | .ok r =>
      let s := buildSetup ops c.params
      match bruteForce s with
      | none => check tally s!"brute[{c.name}]" false "enumeration found no valid assignment"
      | some b =>
        let du := maxDiff (gather s r.u) b.u
        let dm := maxDiff (gather s r.mu) b.mu
        let dobj := (r.cert.objective - b.objective).abs
        check tally s!"brute[{c.name}]: Newton = enumeration"
          (r.converged && du ≤ 1e-8 && dm ≤ 1e-8 && dobj ≤ 1e-9 && b.spread ≤ 1e-7)
          s!"n={s.n} valid assignments={b.numValid} |du|={sci du} |dmu|={sci dm} newton steps={r.itersOne.foldl (· + ·) 0}+{r.itersTwo}"

private def bruteCases : Array Case := Id.run do
  let bumpy := fun (x y : Float) => 0.5 * Float.exp (-6 * ((x - 0.4) * (x - 0.4) + (y - 0.6) * (y - 0.6)))
  let disk4 := (gridDisk 4 4 bumpy).jitter 0.03 11          -- 4 interior
  let disk6 := (gridDisk 5 4 bumpy).jitter 0.03 12          -- 6 interior
  let disk7 := (polar 2 6 0.0 1.0 fun x y => 0.6 * (1 - x * x - y * y) + 0.2 * x).jitter 0.03 13  -- 7 interior
  let ann := (polar 2 6 0.4 1.0 fun x y => 0.4 * x * y + 0.3 * x * x).jitter 0.03 14     -- 6 interior
  let tet := tetrahedron.jitter 0.2 15
  let oct := (bipyramid 4 0.8).jitter 0.1 16
  let bip5 := (bipyramid 5 0.5).jitter 0.05 17               -- 7 vertices
  let cub := cube.jitter 0.15 18
  let u := fun (m : RawMesh) (l : Float) => uniformParams m.pos.size l
  let half := floatPi / 2
  let mut cs : Array Case := #[]
  -- uniform, several λ, every topology
  for l in #[0.002, 0.02, 0.2, 2.0] do
    cs := cs.push { name := s!"disk7/lambda={sci l}", mesh := disk7, params := u disk7 l }
    cs := cs.push { name := s!"annulus/lambda={sci l}", mesh := ann, params := u ann l }
    cs := cs.push { name := s!"octahedron/lambda={sci l}", mesh := oct, params := u oct l }
    cs := cs.push { name := s!"cube/lambda={sci l}", mesh := cub, params := u cub l }
  -- weights
  cs := cs.push { name := "disk7/weighted", mesh := disk7, params := (u disk7 0.01).withRandomWeights 4.0 1 }
  cs := cs.push { name := "cube/weighted", mesh := cub, params := (u cub 0.05).withRandomWeights 4.0 2 }
  cs := cs.push { name := "bipyramid5/weighted", mesh := bip5, params := (u bip5 0.3).withRandomWeights 3.0 3 }
  -- free sites
  cs := cs.push { name := "disk7/free", mesh := disk7, params := (u disk7 0.05).withFree #[0] }
  cs := cs.push { name := "cube/free", mesh := cub, params := (u cub 0.3).withFree #[2] }
  cs := cs.push { name := "tetrahedron/free", mesh := tet, params := (u tet 0.5).withFree #[1] }
  cs := cs.push { name := "octahedron/two-free", mesh := oct, params := (u oct 0.1).withFree #[0, 3] }
  -- forbidden regions
  cs := cs.push { name := "disk7/forbidden", mesh := disk7, params := (u disk7 0.002).withForbidden fun v => v % 2 == 0 }
  cs := cs.push { name := "cube/forbidden", mesh := cub, params := (u cub 0.2).withForbidden fun v => v < 4 }
  -- bounds
  cs := cs.push { name := "disk6/positive", mesh := disk6, params := (u disk6 0.002).withBounds 0.0 posInf }
  cs := cs.push { name := "disk6/negative", mesh := disk6, params := (u disk6 0.002).withBounds negInf 0.0 }
  cs := cs.push { name := "annulus/positive", mesh := ann, params := (u ann 0.002).withBounds 0.0 posInf }
  cs := cs.push { name := "disk4/half-pi", mesh := disk4, params := (u disk4 0.001).withBounds (-half) half }
  cs := cs.push { name := "disk6/tight", mesh := disk6, params := (u disk6 0.0005).withBounds (-0.05) 0.08 }
  -- Σ µ = 4π with six angles ≤ 2π/3: every vertex sits at its bound, no vertex is active
  cs := cs.push { name := "octahedron/all-at-bound", mesh := oct, params := (u oct 0.2).withBounds 0.0 (2 * floatPi / 3) }
  cs := cs.push { name := "octahedron/positive", mesh := oct, params := (u oct 0.2).withBounds 0.0 posInf }
  cs := cs.push { name := "octahedron/tight", mesh := oct, params := (u oct 0.05).withBounds (-1.0) 2.2 }
  cs := cs.push { name := "tetrahedron/capped", mesh := tet, params := (u tet 0.1).withBounds 0.0 3.3 }
  -- prescribed cones
  cs := cs.push { name := "disk7/prescribed", mesh := disk7, params := (u disk7 0.05).withCones #[(0, 0.7), (3, -0.4)] }
  cs := cs.push { name := "cube/prescribed", mesh := cub, params := (u cub 0.2).withCones #[(1, 2.0), (6, -0.5)] }
  -- everything at once
  cs := cs.push { name := "disk6/mixed", mesh := disk6,
                  params := (((u disk6 0.003).withRandomWeights 2.0 5).withBounds (-0.3) 0.4).withCones #[(6, 0.2)] }
  cs := cs.push { name := "octahedron/mixed", mesh := oct,
                  params := ((((u oct 0.1).withRandomWeights 2.0 6).withBounds (-1.0) 2.6).withFree #[5]).withCones #[(2, 1.0)] }
  return cs

private def checkSolver (tally : IO.Ref Tally) : IO Unit := do
  -- closed form: one interior vertex
  let fanMesh := fan 7 0.6
  if let .ok ops := fanMesh.operators then
    let a := ops.L 0 0
    let w := ops.mass.get! 0
    let om := ops.curvature.get! 0
    let mut worst := 0.0
    for l in #[0.0, 1e-4, 1e-3, 1e-2, 0.1, 1.0] do
      for (lo, hi) in #[(negInf, posInf), (0.0, 0.3), (-1.0, 0.0), (0.0, posInf)] do
        let p := (uniformParams ops.nv l).withBounds lo hi
        if let .ok r := solve ops p then
          let expected := boundedThreshold (l * a * a / w) lo hi om
          let d := (r.mu.get! 0 - expected).abs + ((r.u.get! 0) - (expected - om) / a).abs
          if d > worst || !r.converged then worst := if r.converged then d else posInf
    check tally "solver: one interior vertex matches mu = clamp(P_{lambda L^2/m}(Omega))" (worst ≤ 1e-12)
      s!"defect={sci om} worst error={sci worst}"
  -- brute force
  for c in bruteCases do
    checkBrute tally c
  -- random feasible perturbations never lower the objective
  let half := floatPi / 2
  let medium : Array Case := #[
    { name := "grid 7x6", mesh := (gridDisk 7 6 fun x y => 0.4 * Float.sin (4 * x) * Float.sin (3 * y)).jitter 0.01 21,
      params := uniformParams 42 0.003 },
    { name := "annulus weighted", mesh := (polar 3 10 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.2 * y * y).jitter 0.01 22,
      params := (uniformParams 40 0.002).withRandomWeights 3.0 7 },
    { name := "sphere free site", mesh := (uvSphere 5 7 fun p => ⟨1.4 * p.x, p.y, 0.7 * p.z⟩).jitter 0.03 23,
      params := (uniformParams 30 0.05).withFree #[0] },
    { name := "sphere half-pi bounds", mesh := (uvSphere 5 7 fun p => ⟨1.4 * p.x, p.y, 0.7 * p.z⟩).jitter 0.03 23,
      params := (uniformParams 30 0.05).withBounds (-half) half },
    { name := "torus", mesh := (torus 8 5 1.0 0.45 fun a b => 0.25 * Float.sin (2 * a) * Float.cos b).jitter 0.01 24,
      params := uniformParams 40 0.004 },
    { name := "torus positive only", mesh := (torus 8 5 1.0 0.45 fun a b => 0.25 * Float.sin (2 * a) * Float.cos b).jitter 0.01 24,
      params := (uniformParams 40 0.004).withBounds 0.0 posInf }]
  for c in medium do
    match solveCase c with
    | none => check tally s!"perturbation[{c.name}]" false "invalid case"
    | some (ops, r) =>
      let s := buildSetup ops c.params
      let mu := gather s r.mu
      let recomputed := (objectiveOf s mu).getD posInf
      let gain := bestPerturbation s mu 400 31
      let cones := (List.range ops.nv).filter (fun v => r.mu.get! v != 0.0) |>.length
      check tally s!"perturbation[{c.name}]: no feasible perturbation lowers the objective"
        (r.converged && gain ≤ 1e-12 && (recomputed - r.cert.objective).abs ≤ 1e-10)
        s!"nv={ops.nv} cones={cones} objective={sci r.cert.objective} best decrease={sci gain}"
      -- the paper's residual map (Eq. 20) vanishes at the result
      let pr := s.problem
      let res := FloatArray.ofFn' (3 * s.n) (pr.F (vecFn (gather s r.u)) (vecFn (gather s r.phi)) (vecFn mu))
      check tally s!"residual[{c.name}]: F(u, phi, mu) = 0 (Eq. 20)" (maxAbs res ≤ 1e-11)
        s!"max |F|={sci (maxAbs res)}"

private def checkClaims (tally : IO.Ref Tally) : IO Unit := do
  -- large λ on a disk: no cones, u is the cone-free flattening
  let disk := (polar 4 10 0.0 1.0 fun x y => 0.6 * (1 - x * x - y * y)).jitter 0.01 41
  if let some (ops, r) := solveCase { name := "", mesh := disk, params := uniformParams disk.pos.size 10.0 } then
    let noCones := (List.range ops.nv).all fun v => r.mu.get! v == 0.0
    let lu := applyLaplacian ops r.u
    let mut res := 0.0
    for v in [0:ops.nv] do
      if !ops.topo.isBoundary[v]! then
        let d := (lu.get! v + ops.curvature.get! v).abs
        if d > res then res := d
    check tally "claim: large lambda on a disk gives no cones (u solves L u = -Omega)"
      (r.converged && noCones && res ≤ 1e-11) s!"lambda=10 nv={ops.nv}"
  if let some (ops, r) := solveCase { name := "", mesh := disk, params := uniformParams disk.pos.size 0.01 } then
    let cones := (List.range ops.nv).filter (fun v => r.mu.get! v != 0.0) |>.length
    check tally "claim: a smaller lambda on the same disk gives a few cones"
      (r.converged && cones > 0 && cones ≤ ops.nv / 4) s!"lambda=0.01 cones={cones} of {ops.nv}"
  -- closed convex surface without free site: every vertex is a cone, µ = Ω, u = 0
  let ell := uvSphere 6 8 fun p => ⟨1.5 * p.x, 1.1 * p.y, 0.8 * p.z⟩
  let nv := ell.pos.size
  match solveCase { name := "", mesh := ell, params := uniformParams nv 0.05 },
        solveCase { name := "", mesh := ell, params := (uniformParams nv 0.05).withFree #[0] } with
  | some (ops, r), some (_, rf) =>
    let mut dmu := 0.0
    let mut umax := 0.0
    let mut allPos := true
    for v in [0:nv] do
      let d := (r.mu.get! v - ops.curvature.get! v).abs
      if d > dmu then dmu := d
      if (r.u.get! v).abs > umax then umax := (r.u.get! v).abs
      if !(r.mu.get! v > 0) then allPos := false
    check tally "claim: closed convex surface, no free site: a cone at every vertex, mu = Omega, u = 0"
      (r.converged && allPos && dmu ≤ 1e-10 && umax ≤ 1e-10
        && (r.cert.objective - 0.05 * 4 * floatPi).abs ≤ 1e-10)
      s!"nv={nv} objective={sci r.cert.objective} = lambda 4 pi"
    let cones := (List.range nv).filter (fun v => v != 0 && rf.mu.get! v != 0.0) |>.length
    let mut paid := 0.0
    for v in [1:nv] do
      paid := paid + (rf.mu.get! v).abs
    check tally "claim: one free site on the same surface removes most cones"
      (rf.converged && cones < nv / 2 && paid < 4 * floatPi - 1.0 && rf.cert.objective < r.cert.objective
        && (rf.cert.sumMu - 4 * floatPi).abs ≤ 1e-10)
      s!"cones besides the free site={cones} of {nv - 1}, their total |angle|={sci paid}, angle at free site={sci (rf.mu.get! 0)}"
  | _, _ => check tally "claim: closed convex surface" false "invalid case"
  -- cube with angles in [0, π/2]: eight cones of π/2, no vertex strictly inside its bounds
  if let some (_, r) := solveCase { name := "", mesh := cube, params := (uniformParams 8 0.3).withBounds 0.0 (floatPi / 2) } then
    let allHalf := (List.range 8).all fun v => r.mu.get! v == floatPi / 2
    check tally "closed: cube with bounds [0, pi/2] has eight cones of pi/2; constant of phi is an interval"
      (r.converged && allHalf && r.cert.shiftHi - r.cert.shiftLo > 1.0)
      s!"phi_shift=[{sci r.cert.shiftLo}, {sci r.cert.shiftHi}]"
  -- infeasible bounds on a closed mesh
  let infeasible := fun (p : Params) (m : RawMesh) =>
    match solveCase { name := "", mesh := m, params := p } with
    | some (_, r) => (match r.status with | .infeasible _ => true | _ => false)
    | none => false
  check tally "closed: bounds incompatible with Gauss-Bonnet are reported as infeasible"
    (infeasible ((uniformParams 8 0.3).withBounds 0.0 0.1) cube
      && infeasible ((uniformParams 8 0.3).withForbidden fun _ => true) cube
      && infeasible ((uniformParams 8 0.3).withBounds negInf 0.0) cube
      && infeasible ((uniformParams 8 0.3).withCones ((Array.range 8).map fun v => (v, 1.0))) cube)
  if let some (_, r) := solveCase { name := "", mesh := cube, params := (uniformParams 8 0.3).withCones ((Array.range 8).map fun v => (v, if v == 0 then 2 * floatPi else 2 * floatPi / 7)) } then
    check tally "closed: all cones prescribed with total 4 pi is solved (no free cone angle left)"
      (r.converged && r.mu.get! 0 == 2 * floatPi) s!"objective={sci r.cert.objective}"
  -- flat torus: nothing to do
  let flatTorus : RawMesh := Id.run do
    -- a flat torus has no embedding; a regular torus mesh with forbidden cones everywhere is the
    -- cone-free flattening instead
    return torus 6 5 1.0 0.4 fun _ _ => 0.0
  if let some (ops, r) := solveCase { name := "", mesh := flatTorus, params := (uniformParams 30 0.3).withForbidden fun _ => true } then
    let noCones := (List.range ops.nv).all fun v => r.mu.get! v == 0.0
    check tally "closed: torus with cones forbidden everywhere is feasible (chi = 0) and has no cones"
      (r.converged && noCones) s!"objective={sci r.cert.objective}"
  -- prescribed cones are met exactly
  let sph := (uvSphere 5 7 id).jitter 0.03 43
  if let some (_, r) := solveCase { name := "", mesh := sph, params := (uniformParams 30 0.2).withCones #[(3, floatPi), (20, -0.5)] } then
    check tally "prescribed cones are attained exactly and Gauss-Bonnet holds"
      (r.converged && r.mu.get! 3 == floatPi && r.mu.get! 20 == -0.5 && (r.cert.sumMu - 4 * floatPi).abs ≤ 1e-10)

private def checkInvariance (tally : IO.Ref Tally) : IO Unit := do
  let cases : Array Case := #[
    { name := "disk", mesh := (polar 3 9 0.0 1.0 fun x y => 0.6 * (1 - x * x - y * y) + 0.2 * x * y).jitter 0.02 51,
      params := (uniformParams 28 0.004).withRandomWeights 2.0 8 },
    { name := "sphere", mesh := (uvSphere 5 7 fun p => ⟨1.4 * p.x, p.y, 0.7 * p.z⟩).jitter 0.03 52,
      params := ((uniformParams 30 0.05).withRandomWeights 2.0 9).withFree #[4] }]
  for c in cases do
    let nv := c.mesh.pos.size
    match solveCase c with
    | none => check tally s!"invariance[{c.name}]" false "invalid case"
    | some (_, r) =>
      -- uniform scaling
      let scaled := { c with mesh := c.mesh.mapPos fun p => ⟨7.5 * p.x, 7.5 * p.y, 7.5 * p.z⟩ }
      match solveCase scaled with
      | none => check tally s!"invariance[{c.name}]: scaling" false "invalid case"
      | some (_, rs) =>
        let d := maxDiff r.mu rs.mu + maxDiff r.u rs.u + maxDiff r.phi rs.phi
        check tally s!"invariance[{c.name}]: uniform scaling by 7.5 leaves u, phi, mu unchanged"
          (r.converged && rs.converged && d ≤ 1e-9) s!"difference={sci d}"
      -- relabeling
      let perm := randomPerm nv 77
      let permute := fun (x : FloatArray) => Id.run do
        let mut y := zeros nv
        for v in [0:nv] do
          y := y.set! perm[v]! (x.get! v)
        return y
      let p := c.params
      let relabeled : Case := {
        name := c.name, mesh := c.mesh.relabel perm
        params := { lambda := p.lambda, wE := permute p.wE, wR := permute p.wR, lo := permute p.lo,
                    hi := permute p.hi, cones := p.cones.map fun (v, a) => (perm[v]!, a) } }
      match solveCase relabeled with
      | none => check tally s!"invariance[{c.name}]: relabeling" false "invalid case"
      | some (_, rp) =>
        let d := maxDiff (permute r.mu) rp.mu + maxDiff (permute r.u) rp.u + maxDiff (permute r.phi) rp.phi
        check tally s!"invariance[{c.name}]: relabeling the vertices permutes u, phi, mu"
          (rp.converged && d ≤ 1e-9) s!"difference={sci d}"

/-- Run text through the streaming front end. -/
def runText (text : String) (fixed : Option Task) : IO String := do
  let out ← IO.mkRef ""
  let input ← IO.mkRef text.toUTF8
  let pos ← IO.mkRef 0
  let inStream : IO.FS.Stream := {
    flush := pure ()
    read := fun _ => pure ByteArray.empty
    write := fun _ => pure ()
    getLine := do
      let bytes ← input.get
      let start ← pos.get
      if start ≥ bytes.size then return ""
      let mut stop := start
      while stop < bytes.size && bytes[stop]! != 10 do
        stop := stop + 1
      let stop' := if stop < bytes.size then stop + 1 else stop
      pos.set stop'
      return String.fromUTF8! (bytes.extract start stop')
    putStr := fun _ => pure ()
    isTty := pure false }
  let outStream : IO.FS.Stream := {
    flush := pure ()
    read := fun _ => pure ByteArray.empty
    write := fun _ => pure ()
    getLine := pure ""
    putStr := fun s => out.modify (· ++ s)
    isTty := pure false }
  let _ ← runStream inStream outStream fixed
  out.get

private def checkFormat (tally : IO.Ref Tally) : IO Unit := do
  let mesh := (gridDisk 4 4 fun x y => 0.3 * x * y).jitter 0.02 61
  let c : Case := { name := "roundtrip", mesh,
                    params := ((uniformParams 16 0.01).withFree #[5]).withCones #[(6, 0.25)] }
  let text := printCase c
  let out ← runText text (some .solve)
  match solveCase c with
  | none => check tally "format: round trip" false "invalid case"
  | some (ops, r) =>
    check tally "format: a printed case parses back and solves to the identical output"
      (out == printSolve c.name ops r && r.converged)
    let ops' ← runText (printCase c none false) (some .operators)
    check tally "format: operators output for a case without parameter section"
      (ops' == printOperators c.name ops)
  -- batch mode with task lines, comments and blank lines
  let batch := "# two tasks for one mesh\n\n" ++ printCase c (some .operators) ++ "\n" ++ printCase { c with name := "second" } (some .solve)
  let out ← runText batch none
  check tally "format: batch mode runs the task named by each case"
    ((out.splitOn "task operators\nstatus ok").length == 2 && (out.splitOn "case second\ntask solve\nstatus converged").length == 2)
  -- malformed input: one error block per bad case, later cases unaffected
  let good := printCase { c with name := "good" }
  let bad := "case broken\nmesh 3 1\nv 3ff0000000000000 0 0\nend\n"
  let noEnd := "case unterminated\nmesh 3 1\n"
  let out ← runText (bad ++ good ++ "garbage line\nend\n" ++ noEnd ++ good) (some .solve)
  let blocks := (out.splitOn "\nend\n").filter (· != "")
  check tally "format: malformed cases yield error blocks and do not disturb the following cases"
    (blocks.length == 5
      && (out.splitOn "case broken\ntask solve\nstatus error").length == 2
      && (out.splitOn "case ?\ntask solve\nstatus error").length == 2
      && (out.splitOn "case unterminated\ntask solve\nstatus error missing 'end'").length == 2
      && (out.splitOn "case good\ntask solve\nstatus converged").length == 3)
  let invalid := "case nonmanifold\nmesh 5 3\n" ++ "".intercalate ((List.range 5).map fun i => s!"v {floatToHex (Float.ofNat i)} {floatToHex (Float.ofNat (i * i))} {floatToHex 0.0}\n") ++ "f 0 1 2\nf 1 0 3\nf 0 1 4\nend\n"
  let out ← runText invalid (some .operators)
  check tally "format: an invalid mesh is an error with a precise message"
    ((out.splitOn "status error directed edge (0, 1) occurs in faces 0 and 2").length == 2)

/-- Random sweep over meshes of every topology and parameters with every feature. Problems
with `wE = 0` somewhere are left out here: their minimizer need not be unique. -/
private def checkSweep (tally : IO.Ref Tally) (count : Nat) : IO Unit := do
  let mut converged := 0
  let mut infeasible := 0
  let mut failed : Array String := #[]
  let mut worstGain := 0.0
  let mut steps := 0
  let mut used := 0
  let mut seed := 0
  while used < count do
    seed := seed + 1
    let (kind, m) := randomMesh (12 + seed % 30) seed
    let p := randomParams m seed
    if (List.range p.wE.size).any (fun v => p.wE.get! v == 0.0) then continue
    used := used + 1
    match solveCase { name := "", mesh := m, params := p } with
    | none => failed := failed.push s!"{kind}#{seed}: invalid"
    | some (ops, r) =>
      match r.status with
      | .converged =>
        converged := converged + 1
        steps := steps + r.itersTwo + r.itersOne.foldl (· + ·) 0
        if seed % 4 == 0 then
          let s := buildSetup ops p
          let g := bestPerturbation s (gather s r.mu) 60 seed
          if g > worstGain then worstGain := g
      | .infeasible _ => infeasible := infeasible + 1
      | .notConverged why => failed := failed.push s!"{kind}#{seed}: {why}"
  check tally s!"sweep: {count} random cases (all topologies, weights, free sites, forbidden regions, bounds, prescribed cones)"
    (failed.size == 0 && worstGain ≤ 1e-12)
    s!"converged={converged} infeasible={infeasible} failed={failed.size} mean Newton steps={steps / (if converged == 0 then 1 else converged)} best perturbation decrease={sci worstGain}"
  for f in failed do
    IO.println s!"      {f}"

/-- Run all checks; returns the number of failures. -/
def selftest (sweep : Nat := 300) : IO Nat := do
  let tally ← IO.mkRef ({} : Tally)
  checkFormulas tally
  checkOperators tally
  checkSolver tally
  checkClaims tally
  checkInvariance tally
  checkFormat tally
  checkSweep tally sweep
  let t ← tally.get
  IO.println s!"selftest: {t.passed} passed, {t.failed} failed"
  return t.failed

end Ddg
