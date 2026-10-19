import DdgOracle.Developability.Selftest

/-!
# Golden fixtures of the developability oracle

`fixtures/dev-NAME.in|out` and `fixtures/dev-MANIFEST.tsv`. Meshes are small; each energy task prints
all six energies.
-/

namespace Ddg.Developability

open Ddg

private def run1 {α : Type} (g : StateM Rng α) (seed : Nat) : α := (g.run ⟨seed.toUInt64⟩).1

private def mk (name task : String) (mesh : DMesh) : DevCase := { name, task, mesh, energies := EnergyKind.all }

/-- The fixtures: case and description. -/
def fixtures : Array (DevCase × String) := Id.run do
  let mut fs : Array (DevCase × String) := #[]
  let add := fun (fs : Array (DevCase × String)) (c : DevCase) (what : String) => fs.push (c, what)
  -- regular cones (SPEC 10.3)
  for (k, h) in [(4, 0.6), (5, 0.6), (6, 2.0), (8, 0.3), (3, 1.0)] do
    fs := add fs (mk s!"dev-cone-{k}" "energy" (regularCone k h)) s!"regular cone, valence {k}, height {h}: closed forms of SPEC 10.3 (ties: the margins are zero)"
  -- hinges and flat stars
  fs := add fs (mk "dev-hinge-exact" "energy" (run1 (hingeStar 2 3 0.8 0.0) 1)) "exact hinge star (two planar arcs): all energies vanish"
  fs := add fs (mk "dev-hinge-noisy" "energy" (run1 (hingeStar 3 2 0.8 1e-6) 2)) "hinge star with noise 1e-6"
  fs := add fs (mk "dev-flat-noisy" "energy" (run1 (interiorStar 6 0.0 1e-4) 3)) "nearly flat star (heights 1e-4)"
  -- random stars
  fs := add fs (mk "dev-star-3" "energy" (run1 (interiorStar 3 0.5 0.4) 4)) "random interior star of valence 3 (a triple point: not counted by default)"
  fs := add fs (mk "dev-star-6" "energy" (run1 (interiorStar 6 0.7 0.4) 5)) "random interior star of valence 6"
  fs := add fs (mk "dev-star-10" "energy" (run1 (interiorStar 10 0.4 0.5) 6)) "random interior star of valence 10"
  fs := add fs (mk "dev-star-folded" "energy" (run1 (interiorStar 7 0.2 1.8) 7)) "strongly folded star of valence 7 (not embedded)"
  fs := add fs (mk "dev-fan-4" "energy" (run1 (boundaryFan 4 4.0 0.8 0.4) 8)) "boundary fan of 4 faces"
  fs := add fs (mk "dev-fan-8" "energy" (run1 (boundaryFan 8 5.5 0.5 0.5) 9)) "boundary fan of 8 faces"
  fs := add fs { mk "dev-fan-8-noboundary" "energy" (run1 (boundaryFan 8 5.5 0.5 0.5) 9) with includeBoundary := false } "the same fan with boundary vertices excluded from the domain"
  fs := add fs { mk "dev-star-3-all" "energy" (run1 (interiorStar 3 0.5 0.4) 4) with omit3 := false } "valence-3 star with the omission switched off"
  -- small meshes
  let ico := ofRaw ((Ddg.Cuts.icosphere 0).jitter 0.08 11)
  fs := add fs (mk "dev-icosphere" "energy" ico) "jittered icosahedron (closed, valence 5)"
  fs := add fs (mk "dev-cube" "energy" (ofRaw (cube.jitter 0.1 12))) "jittered triangulated cube (valence-3 corners are omitted)"
  fs := add fs { mk "dev-cube-all" "energy" (ofRaw (cube.jitter 0.1 12)) with omit3 := false } "the same cube with the valence-3 omission off"
  fs := add fs (mk "dev-tetra" "energy" (ofRaw (tetrahedron.jitter 0.05 13))) "tetrahedron: every vertex has valence 3, so the domain is empty"
  fs := add fs (mk "dev-octahedron" "energy" (ofRaw ((bipyramid 4 1.0).jitter 0.1 14))) "jittered octahedron (valence 4)"
  let crumpled := ofRaw ((gridDisk 5 5 fun x y => 0.3 * Float.sin (5 * x) * y).jitter 0.25 15)
  fs := add fs (mk "dev-crumpled" "energy" crumpled) "crumpled 5x5 grid (open, 9 interior vertices)"
  fs := add fs { mk "dev-crumpled-noboundary" "energy" crumpled with includeBoundary := false } "the same grid, boundary vertices excluded"
  fs := add fs (mk "dev-annulus" "energy" (ofRaw ((polar 2 6 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.2 * y * y).jitter 0.05 16))) "annulus with two boundary loops"
  fs := add fs (mk "dev-torus" "energy" (ofRaw ((torus 4 3 2.0 0.7 fun a b => 0.1 * Float.sin (a + 2 * b)).jitter 0.05 17))) "small torus"
  -- undefined intrinsic energy: a double cover has a zero vertex normal
  let cover : DMesh := { pos := #[⟨0, 0, 0⟩, ⟨1, 0, 0⟩, ⟨0, 1, 0⟩], faces := #[(0, 1, 2), (0, 2, 1)] }
  fs := add fs (mk "dev-double-cover" "energy" cover) "two coincident faces with opposite orientation: zero vertex normal, the intrinsic energies are undefined"
  fs := add fs (mk "dev-degenerate" "energy" { crumpled with pos := crumpled.pos.set! 6 crumpled.pos[7]! }) "two coincident vertices: degenerate faces"
  fs := add fs (mk "dev-nonmanifold" "energy" { pos := crumpled.pos, faces := crumpled.faces.push (0, 1, 2) }) "a face added twice: the directed edge is used twice"
  -- one step
  for e in EnergyKind.all do
    fs := add fs { mk s!"dev-step-{e.name}" "step" crumpled with energies := #[e] } s!"first L-BFGS step ({e.name}) on the crumpled grid, default options"
  let big : StepOpts := { initialStep := 3.0, maxDisplacement := 2.0, c1 := 1e-4, c2 := 0.9 }
  fs := add fs { mk "dev-step-bisect" "step" ico with energies := #[.covariance], stepOpts := big } "covariance step with a large initial step: bisection"
  let small : StepOpts := { initialStep := 0.001, maxDisplacement := 0.2 }
  fs := add fs { mk "dev-step-double" "step" ico with energies := #[.combinatorial], stepOpts := small } "combinatorial step with a tiny initial step: doubling up to the cap"
  fs := add fs { mk "dev-step-pinned" "step" crumpled with energies := #[.intrinsic], stepOpts := { pinned := #[0, 1, 2, 6, 12] } } "intrinsic step with pinned vertices"
  fs := add fs { mk "dev-step-open" "step" (ofRaw ((polar 2 6 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.2 * y * y).jitter 0.05 16)) with energies := #[.covarianceMax], includeBoundary := false } "lambda_max step on an annulus, boundary excluded"
  -- subdivision
  fs := add fs (mk "dev-subdivide-tetra" "subdivide" (ofRaw (tetrahedron.jitter 0.05 13))) "4-1 subdivision of a tetrahedron"
  fs := add fs (mk "dev-subdivide-fan" "subdivide" (run1 (boundaryFan 5 4.0 0.8 0.4) 8)) "4-1 subdivision of a boundary fan (boundary edges are split too)"
  fs := add fs (mk "dev-subdivide-crumpled" "subdivide" crumpled) "4-1 subdivision of the crumpled grid"
  return fs

private def statusOf (output : String) : String :=
  match (output.splitOn "\n").find? (·.startsWith "status ") with
  | some l => (tokens l)[1]!
  | none => "?"

/-- Compare two outputs token by token: floats to `tol` (absolute and relative), NaN equal to NaN. -/
def compareDev (expected actual : String) (tol : Float := 1e-9) : Option String := Id.run do
  let le := (expected.splitOn "\n").toArray
  let la := (actual.splitOn "\n").toArray
  if le.size != la.size then return some s!"{le.size} lines expected, {la.size} produced"
  for i in [0:le.size] do
    let te := tokens le[i]!
    let ta := tokens la[i]!
    if te.size != ta.size then return some s!"line {i + 1}: '{la[i]!}' instead of '{le[i]!}'"
    for j in [0:te.size] do
      match hexToFloat? te[j]!, hexToFloat? ta[j]! with
      | some x, some y =>
        let scale := if x.abs > 1.0 then x.abs else 1.0
        if !(x == y || (x.isNaN && y.isNaN) || (x - y).abs ≤ tol * scale) then
          return some s!"line {i + 1}, token {j + 1}: {sci y} instead of {sci x} ('{le[i]!}')"
      | _, _ =>
        if te[j]! != ta[j]! then return some s!"line {i + 1}: '{la[i]!}' instead of '{le[i]!}'"
  return none

/-- Write the fixtures. -/
def writeFixtures (dir : System.FilePath) : IO Unit := do
  IO.FS.createDirAll dir
  let mut manifest := ""
  for (c, what) in fixtures do
    let output := runCase c
    IO.FS.writeFile (dir / s!"{c.name}.in") (printCase c)
    IO.FS.writeFile (dir / s!"{c.name}.out") output
    manifest := manifest ++ s!"{c.name}\t{c.task}\t{statusOf output}\t{what}\n"
  IO.FS.writeFile (dir / "MANIFEST.tsv") manifest
  IO.println s!"wrote {fixtures.size} fixtures to {dir}"

/-- Check the fixtures; returns the number of differences. -/
def checkFixtures (dir : System.FilePath) : IO Nat := do
  let mut bad := 0
  for (c, _) in fixtures do
    let name := c.name
    let mut problems : Array String := #[]
    for (suffix, text, exact) in [(".in", printCase c, true), (".out", runCase c, false)] do
      let path := dir / s!"{name}{suffix}"
      if !(← path.pathExists) then
        problems := problems.push s!"{name}{suffix} is missing"
      else
        let stored ← IO.FS.readFile path
        if exact then
          if stored != text then problems := problems.push s!"{name}{suffix} differs from the generated input"
        else if let some d := compareDev stored text then
          problems := problems.push s!"{name}{suffix}: {d}"
    if problems.size == 0 then IO.println s!"ok    {name}"
    else
      bad := bad + 1
      for p in problems do IO.println s!"FAIL  {p}"
  IO.println s!"developability fixtures: {fixtures.size - bad} of {fixtures.size} match"
  return bad

end Ddg.Developability
