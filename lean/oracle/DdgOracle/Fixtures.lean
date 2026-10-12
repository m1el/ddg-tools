import DdgOracle.Selftest

/-!
# Golden fixtures

A fixed list of cases covering every feature of the canonical problem. `fixtures write DIR`
writes, for every case `NAME`,

* `NAME.in` — the case in the input format,
* `NAME.operators.out` — the output of `ddg-oracle operators NAME.in`,
* `NAME.solve.out` — the output of `ddg-oracle solve NAME.in`,

and `MANIFEST.tsv` with one line `NAME <tab> STATUS <tab> DESCRIPTION` per case.
`fixtures check DIR` recomputes everything and compares it with the files.
-/

namespace Ddg

/-- A fixture: a case and what it is for. -/
structure Fixture where
  /-- the case; its name is the file stem -/
  case : Case
  /-- one-line description for the manifest -/
  what : String

private def fx (name what : String) (mesh : RawMesh) (params : Params) : Fixture :=
  { case := { name, mesh, params }, what }

/-- Curvature-adapted weights of the paper's Sec. 8.1: `wE = r`, `wR = 1/r` with
`r = 1 / (|Ω|/m + ε)`, here rescaled so that the median of `r` is `1`. -/
def curvatureWeights (mesh : RawMesh) (lambda eps : Float) : Params :=
  match mesh.operators with
  | .error _ => uniformParams mesh.pos.size lambda
  | .ok ops =>
    let nv := ops.nv
    let raw := arrayOf nv fun v => 1 / ((ops.curvature.get! v).abs / ops.mass.get! v + eps)
    let median := (raw.qsort (· < ·))[nv / 2]!
    { uniformParams nv lambda with
      wE := FloatArray.ofFn' nv fun v => raw[v]! / median
      wR := FloatArray.ofFn' nv fun v => median / raw[v]! }

/-- The fixture list. -/
def fixtures : Array Fixture := Id.run do
  let half := floatPi / 2
  -- meshes
  let fan7 := fan 7 0.6
  let grid := (gridDisk 7 6 fun x y => 0.35 * Float.sin (4 * x) * Float.sin (3 * y)
    + 0.3 * Float.exp (-12 * ((x - 0.3) * (x - 0.3) + (y - 0.6) * (y - 0.6)))).jitter 0.01 101
  let disk := (polar 4 10 0.0 1.0 fun x y => 0.6 * (1 - x * x - y * y) + 0.15 * x * y).jitter 0.01 102
  let ann := (polar 4 12 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.25 * y * y).jitter 0.01 103
  let ell := (uvSphere 6 8 fun p => ⟨1.4 * p.x, 1.1 * p.y, 0.8 * p.z⟩).jitter 0.02 104
  -- strongly non-convex closed surface: total positive curvature 23, negative -10.5
  let knob := (uvSphere 8 10 fun p =>
    let r := 1 + 0.45 * Float.sin (3 * p.x) * Float.sin (3 * p.y) + 0.35 * Float.cos (4 * p.z)
    ⟨r * p.x, r * p.y, r * p.z⟩).jitter 0.004 105
  let cigar := (uvSphere 5 6 fun p => ⟨0.5 * p.x, 0.5 * p.y, 2.5 * p.z⟩).jitter 0.004 109
  let tor := (torus 9 6 1.0 0.42 fun a b => 0.22 * Float.sin (2 * a) * Float.cos b).jitter 0.01 106
  let bigDisk := (gridDisk 12 11 fun x y => 0.3 * Float.sin (5 * x) * Float.cos (4 * y) + 0.2 * x * x).jitter 0.004 107
  let bigSphere := (uvSphere 10 14 fun p =>
    let r := 1 + 0.25 * Float.cos (4 * p.z) + 0.15 * p.x * p.y
    ⟨1.2 * r * p.x, r * p.y, 0.9 * r * p.z⟩).jitter 0.005 108
  let u := fun (m : RawMesh) (l : Float) => uniformParams m.pos.size l
  let inner := fun (m : RawMesh) (k : Nat) => (interiorVertices m)[k]!
  let mut fs : Array Fixture := #[]
  -- uniform weights, several λ, every topology
  fs := fs.push (fx "fan-uniform" "disk with one interior vertex, lambda = 0.02: closed form mu = P(Omega)" fan7 (u fan7 0.02))
  fs := fs.push (fx "grid-lambda-0.5" "disk (grid), uniform, lambda = 0.5 (paper default): no cones" grid (u grid 0.5))
  fs := fs.push (fx "grid-lambda-0.01" "disk (grid), uniform, lambda = 0.01" grid (u grid 0.01))
  fs := fs.push (fx "grid-lambda-0.002" "disk (grid), uniform, lambda = 0.002: several cones of both signs" grid (u grid 0.002))
  fs := fs.push (fx "grid-lambda-0" "disk (grid), lambda = 0: mu = Omega everywhere, u = 0" grid (u grid 0.0))
  fs := fs.push (fx "disk-lambda-0.01" "disk (polar), uniform, lambda = 0.01" disk (u disk 0.01))
  fs := fs.push (fx "disk-lambda-0.1" "disk (polar), uniform, lambda = 0.1: no cones" disk (u disk 0.1))
  fs := fs.push (fx "annulus-lambda-0.005" "annulus, uniform, lambda = 0.005" ann (u ann 0.005))
  fs := fs.push (fx "annulus-lambda-0.0005" "annulus, uniform, lambda = 0.0005" ann (u ann 0.0005))
  fs := fs.push (fx "ellipsoid-lambda-0.5" "closed convex (ellipsoid), uniform: a cone at every vertex, mu = Omega, u = 0" ell (u ell 0.5))
  fs := fs.push (fx "knob-lambda-0.5" "closed non-convex, lambda = 0.5: positive cones only (the same for every lambda down to about 0.001)" knob (u knob 0.5))
  fs := fs.push (fx "knob-lambda-0.0003" "closed non-convex, lambda = 0.0003: cones of both signs" knob (u knob 0.0003))
  fs := fs.push (fx "knob-lambda-0.0001" "closed non-convex, lambda = 0.0001: many cones of both signs" knob (u knob 0.0001))
  fs := fs.push (fx "torus-lambda-0.5" "torus, uniform, lambda = 0.5: no cones, constant of phi is an interval" tor (u tor 0.5))
  fs := fs.push (fx "torus-lambda-0.01" "torus, uniform, lambda = 0.01: cones of both signs summing to 0" tor (u tor 0.01))
  fs := fs.push (fx "torus-lambda-0.001" "torus, uniform, lambda = 0.001" tor (u tor 0.001))
  fs := fs.push (fx "tetrahedron-uniform" "regular tetrahedron: mu = pi at every vertex" tetrahedron (u tetrahedron 0.5))
  fs := fs.push (fx "cube-uniform" "cube: mu = pi/2 at every vertex" cube (u cube 0.5))
  -- weights
  fs := fs.push (fx "grid-weighted" "disk, random wE and wR in [1/3, 3]" grid ((u grid 0.004).withRandomWeights 3.0 201))
  fs := fs.push (fx "annulus-weighted" "annulus, random wE and wR in [1/3, 3]" ann ((u ann 0.002).withRandomWeights 3.0 202))
  fs := fs.push (fx "knob-weighted" "closed non-convex, random wE and wR in [1/4, 4]" knob ((u knob 0.05).withRandomWeights 4.0 203))
  fs := fs.push (fx "torus-weighted" "torus, random wE and wR in [1/3, 3]" tor ((u tor 0.004).withRandomWeights 3.0 204))
  fs := fs.push (fx "knob-curvature-weights" "closed non-convex, curvature-adapted weights of Sec. 8.1 (wE = r, wR = 1/r)" knob (curvatureWeights knob 0.05 1.0))
  fs := fs.push (fx "grid-curvature-weights" "disk, curvature-adapted weights of Sec. 8.1 (wE = r, wR = 1/r)" grid (curvatureWeights grid 0.004 1.0))
  -- free sites
  fs := fs.push (fx "ellipsoid-free-site" "closed convex with one free site (wR = 0): far fewer cones" ell ((u ell 0.05).withFree #[0]))
  fs := fs.push (fx "knob-free-site-lambda-0.5" "closed non-convex, one free site, lambda = 0.5: the free site takes all of 4 pi" knob ((u knob 0.5).withFree #[0]))
  fs := fs.push (fx "knob-free-site-lambda-0.05" "closed non-convex, one free site, lambda = 0.05" knob ((u knob 0.05).withFree #[0]))
  fs := fs.push (fx "knob-two-free-sites" "closed non-convex, two free sites" knob ((u knob 0.05).withFree #[0, 35]))
  fs := fs.push (fx "disk-free-site" "disk with a free site at the centre" disk ((u disk 0.05).withFree #[0]))
  fs := fs.push (fx "torus-free-site" "torus with one free site" tor ((u tor 0.01).withFree #[7]))
  -- forbidden regions
  fs := fs.push (fx "grid-forbidden-half" "disk, cones forbidden (wR = +inf) where x < 0.5" grid ((u grid 0.002).withForbidden fun v => grid.pos[v]!.x < 0.5))
  fs := fs.push (fx "knob-forbidden-north" "closed non-convex, cones forbidden on the upper half" knob ((u knob 0.05).withForbidden fun v => knob.pos[v]!.z > 0.0))
  fs := fs.push (fx "torus-forbidden-all" "torus, cones forbidden everywhere: the cone-free flattening" tor ((u tor 0.01).withForbidden fun _ => true))
  fs := fs.push (fx "grid-visibility" "disk, Sec. 8.1 visibility: wE = 0 where x < 0.5, cones forbidden where x >= 0.5 (u, mu need not be unique where wE = 0)"
    grid (((u grid 0.002).withForbidden fun v => grid.pos[v]!.x ≥ 0.5).withZeroDistortion fun v => grid.pos[v]!.x < 0.5))
  -- bounds
  fs := fs.push (fx "grid-positive-only" "disk, mu >= 0" grid ((u grid 0.002).withBounds 0.0 posInf))
  fs := fs.push (fx "knob-positive-only" "closed non-convex, lambda = 0.0001 and mu >= 0 (Sec. 8.2): compare knob-lambda-0.0001" knob ((u knob 0.0001).withBounds 0.0 posInf))
  fs := fs.push (fx "torus-positive-only" "torus, mu >= 0: forces mu = 0 (chi = 0)" tor ((u tor 0.01).withBounds 0.0 posInf))
  fs := fs.push (fx "grid-half-pi" "disk, -pi/2 <= mu <= pi/2, small lambda" grid ((u grid 0.0005).withBounds (-half) half))
  fs := fs.push (fx "knob-half-pi" "closed non-convex, lambda = 0.0001 and -pi/2 <= mu <= pi/2 (Sec. 8.2): some cones at pi/2" knob ((u knob 0.0001).withBounds (-half) half))
  fs := fs.push (fx "knob-free-sites-half-pi" "closed non-convex, two free sites whose angles are capped at pi/2" knob (((u knob 0.05).withFree #[0, 35]).withBounds (-half) half))
  fs := fs.push (fx "cigar-zero-to-half-pi" "closed convex (elongated), 0 <= mu <= pi/2 (Fig. 23): the two tips are capped" cigar ((u cigar 0.05).withBounds 0.0 half))
  fs := fs.push (fx "cube-zero-to-half-pi" "cube, 0 <= mu <= pi/2: eight cones of pi/2, no vertex strictly inside its bounds, constant of phi is an interval" cube ((u cube 0.3).withBounds 0.0 half))
  fs := fs.push (fx "disk-tight-bounds" "disk, -0.05 <= mu <= 0.08: many cones at their bounds" disk ((u disk 0.001).withBounds (-0.05) 0.08))
  -- prescribed cones
  fs := fs.push (fx "disk-prescribed" "disk, two prescribed cones" disk ((u disk 0.02).withCones #[(inner disk 3, 0.6), (inner disk 17, -0.3)]))
  fs := fs.push (fx "knob-prescribed" "closed non-convex, three prescribed cones (pi, pi/2, -pi/2)" knob ((u knob 0.2).withCones #[(0, floatPi), (20, half), (41, -half)]))
  fs := fs.push (fx "cube-all-prescribed" "cube, every cone prescribed (total 4 pi): nothing left to optimize" cube
    ((u cube 0.3).withCones ((Array.range 8).map fun v => (v, if v == 0 then 2 * floatPi else 2 * floatPi / 7))))
  -- everything at once
  fs := fs.push (fx "grid-mixed" "disk: weights, a free site, a forbidden strip, bounds and a prescribed cone" grid
    ((((((u grid 0.003).withRandomWeights 2.0 205).withBounds (-0.4) 0.5).withFree #[inner grid 4]).withForbidden
      (fun v => grid.pos[v]!.y > 0.75)).withCones #[(inner grid 10, 0.2)]))
  fs := fs.push (fx "knob-mixed" "closed: weights, a free site, a forbidden cap, bounds and a prescribed cone" knob
    ((((((u knob 0.05).withRandomWeights 2.0 206).withBounds (-half) half).withFree #[25]).withForbidden
      (fun v => knob.pos[v]!.z > 0.8)).withCones #[(40, 1.0)]))
  -- larger meshes
  fs := fs.push (fx "big-disk" "disk with 132 vertices, lambda = 0.003" bigDisk (u bigDisk 0.003))
  fs := fs.push (fx "big-sphere" "closed with 128 vertices, lambda = 0.1" bigSphere (u bigSphere 0.1))
  fs := fs.push (fx "big-sphere-weighted-bounded" "closed with 128 vertices, weights and -pi/2 <= mu <= pi/2" bigSphere
    (((u bigSphere 0.03).withRandomWeights 2.0 207).withBounds (-half) half))
  -- hard cases
  fs := fs.push (fx "ellipsoid-zero-distortion-cap" "closed convex, wE = 0 where z > 0.3: the minimizer is not unique (u is free where wE = 0), the exact Newton system is singular; compare the objective only" ell
    ((u ell 0.05).withZeroDistortion fun v => ell.pos[v]!.z > 0.3))
  let (_, hard) := randomMesh 30 175
  fs := fs.push (fx "bipyramid-active-set-cycle" "closed bipyramid, random weights, bounds -pi/2..pi/2, forbidden vertices, one prescribed cone: the simultaneous active-set update cycles on this input"
    hard (randomParams hard 175))
  -- outcomes other than `converged`
  fs := fs.push (fx "cube-infeasible-bounds" "cube, 0 <= mu <= 0.1: Gauss-Bonnet cannot hold, status infeasible" cube ((u cube 0.3).withBounds 0.0 0.1))
  fs := fs.push (fx "error-nonmanifold" "three triangles on one edge: status error"
    { pos := #[⟨0, 0, 0⟩, ⟨1, 0, 0⟩, ⟨0, 1, 0⟩, ⟨0, -1, 0⟩, ⟨0, 0, 1⟩], faces := #[⟨0, 1, 2⟩, ⟨1, 0, 3⟩, ⟨0, 1, 4⟩] } (uniformParams 5 0.5))
  fs := fs.push (fx "error-degenerate-triangle" "a triangle with collinear vertices: status error"
    { pos := #[⟨0, 0, 0⟩, ⟨1, 0, 0⟩, ⟨2, 0, 0⟩, ⟨0, 1, 0⟩], faces := #[⟨0, 1, 2⟩, ⟨0, 2, 3⟩] } (uniformParams 4 0.5))
  fs := fs.push (fx "error-negative-weight" "wE < 0 at a vertex: operators succeed, solve reports status error"
    fan7 { u fan7 0.5 with wE := (u fan7 0.5).wE.set! 0 (-1.0) })
  return fs

/-- The three texts of a fixture: input, operators output, solve output. -/
def Fixture.render (f : Fixture) : IO (String × String × String) := do
  let input := printCase f.case
  return (input, ← runText input (some .operators), ← runText input (some .solve))

/-- Status word of an output block. -/
def statusOf (output : String) : String :=
  match (output.splitOn "\n").find? (·.startsWith "status ") with
  | some l => (tokens l)[1]!
  | none => "?"

/-- Write all fixtures to `dir`. -/
def writeFixtures (dir : System.FilePath) : IO Unit := do
  IO.FS.createDirAll dir
  let mut manifest := ""
  for f in fixtures do
    let (input, ops, sol) ← f.render
    IO.FS.writeFile (dir / s!"{f.case.name}.in") input
    IO.FS.writeFile (dir / s!"{f.case.name}.operators.out") ops
    IO.FS.writeFile (dir / s!"{f.case.name}.solve.out") sol
    manifest := manifest ++ s!"{f.case.name}\t{statusOf ops}/{statusOf sol}\t{f.what}\n"
    IO.println s!"{f.case.name}: operators {statusOf ops}, solve {statusOf sol}"
  IO.FS.writeFile (dir / "MANIFEST.tsv") manifest
  IO.println s!"wrote {fixtures.size} fixtures to {dir}"

/-- Compare two output blocks: identical line structure and tags, integers equal, floats equal
up to `tol` (absolute, and relative for large values). The `iterations` line is not compared.
Returns a description of the first difference. -/
def compareOutputs (expected actual : String) (tol : Float := 1e-9) : Option String := Id.run do
  let le := (expected.splitOn "\n").toArray
  let la := (actual.splitOn "\n").toArray
  if le.size != la.size then return some s!"{le.size} lines expected, {la.size} produced"
  for i in [0:le.size] do
    let te := tokens le[i]!
    let ta := tokens la[i]!
    if te.size == 0 && ta.size == 0 then continue
    if te.size != ta.size then return some s!"line {i + 1}: '{la[i]!}' instead of '{le[i]!}'"
    if te[0]! == "iterations" && ta[0]! == "iterations" then continue
    if te[0]! == "status" then
      if te[1]! != ta[1]! then return some s!"line {i + 1}: '{la[i]!}' instead of '{le[i]!}'"
      continue
    for j in [0:te.size] do
      match hexToFloat? te[j]!, hexToFloat? ta[j]! with
      | some x, some y =>
        let scale := if x.abs > 1.0 then x.abs else 1.0
        if !(x == y || (x - y).abs ≤ tol * scale) then
          return some s!"line {i + 1}, token {j + 1}: {sci y} instead of {sci x} ('{le[i]!}')"
      | _, _ =>
        if te[j]! != ta[j]! then return some s!"line {i + 1}: '{la[i]!}' instead of '{le[i]!}'"
  return none

/-- Recompute all fixtures and compare them with the files in `dir`. Returns the number of
differences. -/
def checkFixtures (dir : System.FilePath) : IO Nat := do
  let mut bad := 0
  for f in fixtures do
    let (input, ops, sol) ← f.render
    let name := f.case.name
    let mut problems : Array String := #[]
    for (suffix, text, exact) in [(".in", input, true), (".operators.out", ops, false), (".solve.out", sol, false)] do
      let path := dir / s!"{name}{suffix}"
      if !(← path.pathExists) then
        problems := problems.push s!"{name}{suffix} is missing"
      else
        let stored ← IO.FS.readFile path
        if exact then
          if stored != text then problems := problems.push s!"{name}{suffix} differs from the generated input"
        else if let some d := compareOutputs stored text then
          problems := problems.push s!"{name}{suffix}: {d}"
    if problems.size == 0 then IO.println s!"ok    {name}"
    else
      bad := bad + 1
      for p in problems do
        IO.println s!"FAIL  {p}"
  IO.println s!"fixtures: {fixtures.size - bad} of {fixtures.size} match"
  return bad

end Ddg
