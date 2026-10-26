import DdgOracle.RepulsiveCurves.Selftest

/-!
# Golden fixtures of the `rc` oracle

`fixtures/rc/NAME.in|out` and `fixtures/rc/MANIFEST.tsv`: circle, trefoil, open curve with tangent
sites, tripod junction, Hopf link, planar curves, sphere constraints, two exponent pairs, both
quadratures, projection endings, and invalid inputs.
-/

namespace Ddg.RepulsiveCurves

private def run1 {α : Type} (g : StateM Rng α) (seed : Nat) : α := (g.run ⟨seed.toUInt64⟩).1

private def mk (name task : String) (net : Net) : RCCase := { name, task, net }

/-- The current tangent of a site plus an offset: a feasible direction `X` that is not equal to `s`. -/
private def dirOf (net : Net) (site : Site) (off : V3) : V3 :=
  match tangentOf net net.pos site with
  | .ok (t, _) => t + off
  | .error _ => off

private def ex36 : Expo := ⟨3.0, 6.0⟩
private def ex245 : Expo := ⟨2.0, 4.5⟩

/-- The fixtures: case and description. -/
def fixtures : Array (RCCase × String) := Id.run do
  let mut fs : Array (RCCase × String) := #[]
  let add := fun (fs : Array (RCCase × String)) (c : RCCase) (what : String) => fs.push (c, what)
  let circ := circle 24 1.0
  let tre80 := trefoil 80 0.25
  let tre40 := trefoil 40 0.25
  let helix := run1 (helix 16 0.04) 31
  let tri := run1 (tripod 4 0.03) 32
  let hop := run1 (hopf 12 0.05) 33
  let pl := (run1 (wobbleLoop 18 0.3 3 0.0 0.03) 34)
  let planar : Net := { pl with pos := pl.pos.map fun p => ⟨p.x, p.y, 0.0⟩ }
  let planarY : Net := { pl with pos := pl.pos.map fun p => ⟨p.x, 0.5, p.y⟩ }
  let bumpy := run1 (wobbleLoop 14 0.3 5 0.3 0.05) 35
  -- energy
  fs := add fs (mk "rc-circle-energy" "energy" circ) "regular 24-gon, (3,6): radial equal-norm differential"
  fs := add fs { mk "rc-circle-energy-245" "energy" circ with ex := ex245 } "regular 24-gon, (2,4.5)"
  fs := add fs (mk "rc-trefoil80-energy" "energy" tre80) "trefoil with 80 vertices, (3,6)"
  fs := add fs { mk "rc-trefoil80-energy-245" "energy" tre80 with ex := ex245 } "trefoil with 80 vertices, (2,4.5)"
  fs := add fs { mk "rc-helix-energy" "energy" helix with ex := ex245 } "open helical arc, (2,4.5)"
  fs := add fs (mk "rc-tripod-energy" "energy" tri) "3-junction graph, (3,6)"
  fs := add fs { mk "rc-hopf-energy" "energy" hop with ex := ⟨2.5, 5.5⟩ } "two linked loops, (2.5,5.5)"
  fs := add fs { mk "rc-planar-energy" "energy" planar with planar := some 2 } "planar loop in z = 0: no differential along z"
  -- gram
  fs := add fs (mk "rc-circle-gram" "gram" circ) "24-gon: B, B0, A, K^B, K0, both quadratures, (3,6)"
  fs := add fs { mk "rc-trefoil40-gram" "gram" tre40 with ex := ex245 } "trefoil with 40 vertices, (2,4.5), both quadratures"
  fs := add fs { mk "rc-trefoil80-gram" "gram" tre80 with allMatrices := false } "trefoil with 80 vertices: A only, both quadratures"
  fs := add fs (mk "rc-helix-gram" "gram" helix) "open helix, (3,6)"
  fs := add fs { mk "rc-tripod-gram" "gram" tri with ex := ex245 } "tripod, (2,4.5)"
  fs := add fs (mk "rc-hopf-gram" "gram" hop) "Hopf link, (3,6): B0 couples the two components"
  fs := add fs { mk "rc-planar-gram" "gram" planar with planar := some 2 } "planar loop, (3,6)"
  -- saddle: one case per constraint type, plus combinations
  fs := add fs { mk "rc-circle-saddle-lengths" "saddle" circ with cons := #[.barycenter none, .edgeLengths none] } "24-gon, barycenter + all edge lengths: the projected gradient vanishes (analytic)"
  fs := add fs { mk "rc-circle-saddle-point" "saddle" circ with cons := #[.point 0 none, .totalLength none] } "24-gon, point 0 + total length"
  fs := add fs { mk "rc-trefoil-saddle" "saddle" tre80 with cons := #[.barycenter none, .totalLength none], quad := .trapezoid } "trefoil (80), barycenter + total length, trapezoid"
  fs := add fs { mk "rc-trefoil-saddle-245" "saddle" tre40 with ex := ex245, cons := #[.barycenter none, .edgeLengths none] } "trefoil (40), (2,4.5), barycenter + all edge lengths"
  fs := add fs { mk "rc-helix-saddle-tangent" "saddle" helix with cons := #[.barycenter none, .totalLength none, .tangent (.vertex 0) none, .tangent (.vertex 15) (some (dirOf helix (.vertex 15) ⟨0.2, -0.3, 0.1⟩)), .tangent (.edge 7) none] } "open curve: tangents at both ends (degree-1 vertices) and at an edge"
  fs := add fs { mk "rc-helix-saddle-tangent-v" "saddle" helix with ex := ex245, quad := .trapezoid, cons := #[.barycenter none, .edgeLengths (some #[2, 3, 9]), .tangent (.vertex 5) (some (dirOf helix (.vertex 5) ⟨-0.1, 0.2, 0.1⟩))] } "degree-2 vertex tangent (sum of unit edge vectors), subset of edge lengths"
  fs := add fs { mk "rc-tripod-saddle" "saddle" tri with cons := #[.barycenter none, .edgeLengths none] } "tripod, barycenter + all edge lengths"
  fs := add fs { mk "rc-tripod-saddle-tangent" "saddle" tri with ex := ex245, cons := #[.barycenter none, .totalLength none, .tangent (.vertex 4) none, .tangent (.vertex 8) none, .tangent (.vertex 12) none] } "tripod, tangents at the three free ends"
  fs := add fs { mk "rc-hopf-saddle" "saddle" hop with ex := ⟨2.5, 5.5⟩, cons := #[.barycenter none, .totalLength none] } "Hopf link, barycenter + total length"
  fs := add fs { mk "rc-hopf-saddle-points" "saddle" hop with cons := #[.point 0 none, .point 13 none, .edgeLengths (some #[0, 1, 2, 3, 12, 13])] } "Hopf link, two point constraints + some edge lengths"
  let (c0, r0) := centroidRadius bumpy
  fs := add fs { mk "rc-bumpy-saddle-sphere" "saddle" bumpy with cons := #[.barycenter none, .sphere c0 r0 (some #[0, 3, 6, 9, 12])] } "wobbly loop, barycenter + implicit sphere |x - c|^2 - R^2 on five vertices"
  fs := add fs { mk "rc-bumpy-saddle-sphere-all" "saddle" bumpy with ex := ex245, cons := #[.sphere c0 r0 none, .totalLength none, .barycenter none] } "sphere on all vertices, total length and barycenter (k = 18 rows), (2,4.5)"
  fs := add fs { mk "rc-planar-saddle" "saddle" planar with planar := some 2, cons := #[.barycenter none, .edgeLengths none] } "planar loop (z held), barycenter (2 rows) + edge lengths"
  fs := add fs { mk "rc-planar-saddle-y" "saddle" planarY with planar := some 1, ex := ex245, quad := .trapezoid, cons := #[.point 0 none, .totalLength none, .sphere ⟨0.0, 0.5, 0.0⟩ 1.0 (some #[2, 5])] } "loop in the plane y = 0.5 (y held), point + total length + sphere"
  -- projection
  let big := fun (n : Net) => (1.0 : Float) * 0.0 + (n.nv.toFloat * 0.0 + 0.3)
  fs := add fs { mk "rc-circle-project-point" "project" circ with cons := #[.point 0 none, .totalLength none], tau := 0.3, tauRel := true } "24-gon: step of 0.3 edge lengths, projection onto point + total length"
  fs := add fs { mk "rc-trefoil-project" "project" tre80 with cons := #[.barycenter none, .totalLength none], tau := big tre80, tauRel := true } "trefoil (80): chord projection, default tolerance"
  fs := add fs { mk "rc-trefoil-project-lengths" "project" tre40 with ex := ex245, quad := .trapezoid, cons := #[.barycenter none, .edgeLengths none], tau := 0.1, tauRel := true, tolVal := 1e-8, maxit := 6 } "trefoil (40): all edge lengths, tolerance 1e-8 edge lengths, up to 6 iterations"
  fs := add fs { mk "rc-helix-project-tangent" "project" helix with cons := #[.barycenter none, .totalLength none, .tangent (.vertex 0) none, .tangent (.vertex 15) none, .tangent (.edge 7) none], tau := 0.4, tauRel := true, tolVal := 1e-6, maxit := 5 } "open curve: tangent rows are dimensionless (tolerance tau/lbar)"
  fs := add fs { mk "rc-helix-project-maxit" "project" helix with cons := #[.barycenter none, .edgeLengths none], tau := 1.0, tauRel := true, tolVal := 1e-12, maxit := 1 } "iteration limit 1 with a tight tolerance: max-iterations"
  fs := add fs { mk "rc-helix-project-abs" "project" helix with ex := ex245, cons := #[.barycenter none, .totalLength none], tau := 0.2, tauRel := true, tolMean := false, tolVal := 1e-5 } "absolute tolerance"
  fs := add fs { mk "rc-tripod-project" "project" tri with cons := #[.barycenter none, .edgeLengths none, .tangent (.vertex 4) none], tau := 0.08, tauRel := true, tolVal := 1e-7, maxit := 4 } "tripod: lengths + a tangent"
  fs := add fs { mk "rc-hopf-project" "project" hop with ex := ⟨2.5, 5.5⟩, cons := #[.barycenter none, .totalLength none], tau := 0.3, tauRel := true } "Hopf link"
  fs := add fs { mk "rc-bumpy-project-sphere" "project" bumpy with cons := #[.barycenter none, .sphere c0 r0 (some #[0, 3, 6, 9, 12])], tau := 0.1, tauRel := true, tolVal := 1e-6, maxit := 5 } "surface rows: tolerance scaled by |grad f|"
  fs := add fs { mk "rc-planar-project" "project" planar with planar := some 2, cons := #[.barycenter none, .edgeLengths none], tau := 0.3, tauRel := true, tolVal := 1e-7, maxit := 5 } "planar projection"
  fs := add fs { mk "rc-circle-project-diverge" "project" circ with cons := #[.point 0 none, .edgeLengths (some #[0, 1, 2, 3, 4, 5, 6, 7])], tau := 40.0, tauRel := true, maxit := 3 } "very large step (40 edge lengths): the iteration limit, whatever the ending"
  -- errors (status error KIND ...)
  fs := add fs { mk "rc-error-exponents" "energy" circ with ex := ⟨3.0, 4.5⟩ } "beta < alpha + 2: InvalidExponents"
  fs := add fs { mk "rc-error-exponents-upper" "energy" circ with ex := ⟨2.0, 5.0⟩ } "beta >= 2 alpha + 1: InvalidExponents"
  fs := add fs (mk "rc-error-selfloop" "energy" { circ with edges := circ.edges.set! 3 (3, 3) }) "edge (3,3): SelfLoop"
  fs := add fs (mk "rc-error-duplicate" "energy" { circ with edges := circ.edges.push (1, 0) }) "edge (1,0) twice: DuplicateEdge"
  fs := add fs (mk "rc-error-isolated" "energy" { circ with pos := circ.pos.push ⟨5.0, 5.0, 5.0⟩ }) "a vertex without edges: IsolatedVertex"
  fs := add fs (mk "rc-error-zero-edge" "energy" { circ with pos := circ.pos.set! 4 circ.pos[3]! }) "coincident neighbours: ZeroLengthEdge"
  fs := add fs (mk "rc-error-coincident" "energy" { circ with pos := circ.pos.set! 10 circ.pos[3]! }) "coincident non-adjacent vertices: CoincidentVertices"
  fs := add fs (mk "rc-error-no-pairs" "energy" (circle 3 1.0)) "a triangle has no disjoint edges: NoInteractingPairs"
  fs := add fs (mk "rc-error-out-of-range" "energy" { circ with edges := circ.edges.set! 3 (3, 99) }) "edge to vertex 99: VertexOutOfRange"
  fs := add fs { mk "rc-error-not-planar" "energy" circ with planar := some 2, net := { circ with pos := circ.pos.set! 3 ⟨circ.pos[3]!.x, circ.pos[3]!.y, 0.5⟩ } } "planar axis but one vertex off the plane: NotPlanar"
  fs := add fs { mk "rc-error-singular" "saddle" circ with cons := #[.point 0 none, .point 0 none] } "the same point twice: the saddle matrix is singular"
  fs := add fs { mk "rc-error-redundant" "saddle" bumpy with cons := #[.point 1 none, .sphere c0 r0 (some #[1, 4])] } "a surface row at a point-constrained vertex is redundant (D9): singular (the surface value at vertex 1 is not zero, but the row is dependent)"
  fs := add fs { mk "rc-error-tangent-degree" "saddle" tri with cons := #[.barycenter none, .tangent (.vertex 0) none] } "tangent at the degree-3 junction: InvalidTangentSite"
  fs := add fs { mk "rc-error-tangent-flipped" "saddle" helix with cons := #[.barycenter none, .tangent (.edge 3) (some ((-1.0 : Float) * (helix.pos[4]! - helix.pos[3]!)))] } "tangent opposite to its direction: InfeasibleConstraint"
  return fs

private def statusOf (output : String) : String :=
  match (output.splitOn "\n").find? (·.startsWith "status ") with
  | some l => (tokens l)[1]!
  | none => "?"

/-- Compare two outputs token by token: floats to `tol` (relative to `max(|x|, 1)`), NaN equal to NaN. -/
def compareRC (expected actual : String) (tol : Float := 1e-9) : Option String := Id.run do
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
          return some s!"line {i + 1}, token {j + 1}: differs ('{le[i]!}')"
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
    let st := statusOf output
    let detail := if st == "error" then ((output.splitOn "\n").find? (·.startsWith "status ")).getD "" else
      match (output.splitOn "\n").find? (·.startsWith "result ") with | some l => l | none => ""
    manifest := manifest ++ s!"{c.name}\t{c.task}\t{st}\t{what}\t{detail}\n"
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
        else if let some d := compareRC stored text then
          problems := problems.push s!"{name}{suffix}: {d}"
    if problems.size == 0 then IO.println s!"ok    {name}"
    else
      bad := bad + 1
      for p in problems do IO.println s!"FAIL  {p}"
  IO.println s!"rc fixtures: {fixtures.size - bad} of {fixtures.size} match"
  return bad

end Ddg.RepulsiveCurves
