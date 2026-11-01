import DdgOracle.RepulsiveSurfaces.Selftest

/-!
# Golden fixtures of the `rs` oracle

`fixtures/rs/NAME.in|out` and `fixtures/rs/MANIFEST.tsv`: every task on perturbed polyhedra, tori,
a patch with boundary, a Möbius band and two components, several exponent pairs, and invalid inputs.
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

private def mk (name task : String) (mesh : Mesh) : RSCase := { name, task, mesh }

private def ex612 : Ex := ⟨6.0, 12.0⟩
private def ex48 : Ex := ⟨4.0, 8.0⟩
private def ex36 : Ex := ⟨3.0, 6.0⟩
private def ex25 : Ex := ⟨2.0, 5.0⟩
private def ex375 : Ex := ⟨3.0, 7.5⟩

private def constU (n : Nat) (seed : Nat) : Array Float := Id.run do
  let mut rng : Rng := ⟨seed.toUInt64 * 31 + 5⟩
  let mut a : Array Float := #[]
  for _ in [0:n] do
    let (x, r) := rng.range (-1.0) 1.0
    a := a.push x
    rng := r
  return a

private def randV3 (n seed : Nat) : Array V3 := Id.run do
  let u := constU (3 * n) seed
  return (Array.range n).map fun i => ⟨u[3 * i]!, u[3 * i + 1]!, u[3 * i + 2]!⟩

private def aOf (m : Mesh) : Float := match m.validate with
  | .ok t => totalArea t (faceGeo t.faces m.pos) none
  | .error _ => 0.0

private def vOf (m : Mesh) : Float := match m.validate with
  | .ok t => totalVolume t m.pos none
  | .error _ => 0.0

/-- The fixtures: case and description. -/
def fixtures : Array (RSCase × String) := Id.run do
  let mut fs : Array (RSCase × String) := #[]
  let add := fun (fs : Array (RSCase × String)) (c : RSCase) (what : String) => fs.push (c, what)
  let tet := tetra.jitter 0.12 101
  let oct := octa.jitter 0.1 102
  let ico := icosa.jitter 0.12 103
  let tor := (torus 5 3 2.0 0.8).jitter 0.05 104
  let pat := (grid 3 2 0.3).jitter 0.05 105
  let mob := (mobius 8).jitter 0.03 106
  let two := (disjoint (Mesh.scale tetra 0.6) octa ⟨4.0, 0.3, 0.2⟩).jitter 0.05 107
  let flat := grid 3 3 0.0
  -- energy
  fs := add fs { mk "rs-tetra-energy" "energy" tet with ex := ex612 } "perturbed tetrahedron, (6,12)"
  fs := add fs { mk "rs-octa-energy-48" "energy" oct with ex := ex48 } "perturbed octahedron, (4,8) subcritical"
  fs := add fs { mk "rs-icosa-energy" "energy" ico with ex := ex612 } "perturbed icosahedron, (6,12)"
  fs := add fs { mk "rs-icosa-energy-36" "energy" ico with ex := ex36 } "perturbed icosahedron, (3,6)"
  fs := add fs { mk "rs-torus-energy" "energy" tor with ex := ex612 } "perturbed 5x3 torus, (6,12)"
  fs := add fs { mk "rs-patch-energy-25" "energy" pat with ex := ex25 } "open patch with boundary, (2,5)"
  fs := add fs { mk "rs-mobius-energy" "energy" mob with ex := ex612 } "Moebius band (non-orientable), (6,12)"
  fs := add fs { mk "rs-two-energy-375" "energy" two with ex := ex375 } "two components (tetrahedron and octahedron), (3,7.5)"
  fs := add fs { mk "rs-flat-energy" "energy" flat with ex := ex612 } "flat 3x3 patch: N.(X-Y) = 0, energy and gradient exactly 0"
  -- df
  fs := add fs { mk "rs-tetra-df" "df" tet with u := #[tet.pos.map (fun p => 0.3 * p.x - 0.5 * p.y + 0.2 * p.z), constU 4 1] } "tetrahedron: a linear function (tangential part of c) and a random one"
  fs := add fs { mk "rs-icosa-df" "df" ico with u := #[constU 12 2, Array.replicate 12 1.0] } "icosahedron: random u, and u = 1 (D_f 1 = 0)"
  fs := add fs { mk "rs-patch-df" "df" pat with u := #[constU pat.nv 3] } "patch with boundary: random u"
  fs := add fs { mk "rs-mobius-df" "df" mob with u := #[constU mob.nv 4] } "Moebius band: random u (D_f is orientation free)"
  -- metric
  fs := add fs { mk "rs-tetra-metric" "metric" tet with ex := ex612, tau := 0.6 } "tetrahedron: L_0.6, B, B0, A, W, M, (6,12)"
  fs := add fs { mk "rs-octa-metric-nolow" "metric" oct with ex := ex48, tau := 0.25, mo := { clo := 0.0 } } "octahedron: low_order off (A = B), (4,8), tau 0.25"
  fs := add fs { mk "rs-icosa-metric-willmore" "metric" ico with ex := ex612, tau := 1.2, mo := { wtp := 2.0, clo := 1.0, ww := 0.3 } } "icosahedron: weights 2, 1, 0.3 (adds W M^-1 W), tau 1.2"
  fs := add fs { mk "rs-torus-metric" "metric" tor with ex := ex36, tau := 0.9 } "torus 5x3, (3,6)"
  fs := add fs { mk "rs-patch-metric" "metric" pat with ex := ex25, tau := 0.5 } "open patch, (2,5): boundary vertices"
  fs := add fs { mk "rs-mobius-metric" "metric" mob with ex := ex612 } "Moebius band, (6,12)"
  fs := add fs { mk "rs-two-metric-375" "metric" two with ex := ex375, tau := 0.4 } "two components: block structure through the kernel"
  -- constraints
  fs := add fs { mk "rs-two-cons-bary" "constraints" two with cons := #[.barycenter none none] } "per-component barycenters of two components (6 rows), target = current"
  fs := add fs { mk "rs-icosa-cons-pins" "constraints" ico with cons := #[.pins #[2, 7] (some #[ico.pos[2]! + ⟨0.1, 0.0, -0.05⟩, ico.pos[7]!]), .barycenter (some 0) (some ⟨0.05, 0.0, 0.0⟩)] } "two pins with targets and a barycenter with target"
  fs := add fs { mk "rs-icosa-cons-area-volume" "constraints" ico with cons := #[.area none (some (aOf ico * 1.05)), .volume none (some (vOf ico * 0.9))] } "total area x1.05 and volume x0.9: values and Jacobians"
  fs := add fs { mk "rs-two-cons-comp" "constraints" two with cons := #[.area (some 1) none, .volume (some 0) (some 1.0)] } "area of component 1, volume of component 0 (target 1)"
  fs := add fs { mk "rs-patch-cons-open" "constraints" pat with cons := #[.area none none, .volume none none, .barycenter none none] } "open patch: area and the signed cone volume over the origin"
  fs := add fs { mk "rs-mobius-cons-area" "constraints" mob with cons := #[.area none (some (aOf mob * 1.1)), .barycenter none none] } "Moebius band: area and barycenter (volume would be refused)"
  -- saddle
  fs := add fs { mk "rs-icosa-saddle" "saddle" ico with ex := ex612, cons := #[.barycenter none none, .area none (some (aOf ico * 1.03)), .volume none (some (vOf ico * 0.97))] } "icosahedron: barycenter in A', area and volume as Schur rows"
  fs := add fs { mk "rs-tetra-saddle-pins" "saddle" tet with ex := ex36, cons := #[.pins #[0] none, .pins #[1, 2] (some #[tet.pos[1]! + ⟨0.0, 0.05, 0.0⟩, tet.pos[2]!]), .volume none (some (vOf tet * 1.05))] } "tetrahedron: three pins and a volume row, (3,6)"
  fs := add fs { mk "rs-torus-saddle" "saddle" tor with ex := ex612, mo := { wtp := 1.0, clo := 0.0, ww := 0.0 }, cons := #[.barycenter none none, .volume none (some (vOf tor * 1.04))] } "torus: low_order off, barycenter and volume"
  fs := add fs { mk "rs-patch-saddle" "saddle" pat with ex := ex25, cons := #[.pins #[0, 5] none, .area none (some (aOf pat * 1.02))] } "open patch: two pinned vertices and an area row, (2,5)"
  fs := add fs { mk "rs-two-saddle" "saddle" two with ex := ex375, cons := #[.barycenter none none, .area none (some (aOf two * 1.02))] } "two components: per-component barycenters, one total-area row"
  fs := add fs { mk "rs-mobius-saddle" "saddle" mob with ex := ex612, cons := #[.barycenter none none, .area none (some (aOf mob * 0.95))] } "Moebius band: barycenter and area"
  fs := add fs { mk "rs-octa-saddle-none" "saddle" oct with ex := ex48, cons := #[.barycenter none none] } "octahedron: no Schur rows (k = 0), (4,8)"
  -- precond
  fs := add fs { mk "rs-icosa-precond" "precond" ico with ex := ex612, cons := #[.barycenter none none], r := randV3 12 7 } "icosahedron: eps_rel 1e-10, barycenter in L'"
  fs := add fs { mk "rs-patch-precond" "precond" pat with ex := ex25, epsRel := 1e-6, cons := #[.pins #[0] none, .barycenter none none], r := randV3 pat.nv 8 } "open patch: pin and barycenter, eps_rel 1e-6, (2,5)"
  fs := add fs { mk "rs-two-precond" "precond" two with ex := ex375, cons := #[.barycenter none none], r := randV3 two.nv 9 } "two components: per-component barycenter rows"
  -- step
  fs := add fs { mk "rs-icosa-step" "step" ico with ex := ex612, cons := #[.barycenter none none, .area none none, .volume none none] } "icosahedron: area and volume held, start: farthest vertex moves one mean edge"
  fs := add fs { mk "rs-torus-step-ref" "step" tor with ex := ex612, so := { start := .reference, newtonFixed := true, newton := 1 }, cons := #[.barycenter none none, .volume none (some (vOf tor * 1.02))] } "torus: reference start 2/|x|, one Newton step, volume target x1.02"
  fs := add fs { mk "rs-patch-step" "step" pat with ex := ex25, so := { start := .maxdisp 0.5 }, cons := #[.pins #[0, 3] none, .area none (some (aOf pat * 1.01))] } "open patch: pins and area target x1.01"
  fs := add fs { mk "rs-mobius-step" "step" mob with ex := ex612, cons := #[.barycenter none none, .area none none] } "Moebius band: area held"
  fs := add fs { mk "rs-two-step" "step" two with ex := ex375, so := { armijo := 0.3 }, cons := #[.barycenter none none] } "two components: Armijo c = 0.3"
  fs := add fs { mk "rs-icosa-step-stuck" "step" ico with ex := ex612, so := { start := .abs 1.0e6, maxHalv := 3 }, cons := #[.barycenter none none] } "huge start, 3 halvings allowed: stuck, positions unchanged"
  fs := add fs { mk "rs-tetra-step-flatmix" "step" tet with ex := ex48, so := { start := .abs 0.05 }, cons := #[.pins #[0] none, .area none none] } "tetrahedron: fixed start 0.05, one pin, area held"
  -- remesh-ops
  let opCase := fun (nm kind : String) (m : Mesh) (a b : Nat) => { mk nm "remesh-ops" m with op := some (kind, a, b) }
  fs := add fs (opCase "rs-icosa-split" "split" ico 0 1) "split an interior edge of the icosahedron: new vertex 12, faces 20..21 appended"
  fs := add fs (opCase "rs-patch-split-boundary" "split" pat 0 1) "split a boundary edge of the patch (one face)"
  fs := add fs (opCase "rs-icosa-collapse" "collapse" ico 0 1) "collapse an edge of the icosahedron to its midpoint (vertex 1 removed)"
  fs := add fs { opCase "rs-icosa-collapse-pinned" "collapse" ico 3 4 with pinned := #[4] } "collapse with vertex 4 pinned: survivor 3 sits at vertex 4"
  fs := add fs (opCase "rs-patch-collapse-boundary" "collapse" pat 5 9) "collapse an interior edge from a boundary vertex: the boundary vertex keeps its position"
  fs := add fs { opCase "rs-patch-collapse-both-pinned" "collapse" pat 5 6 with pinned := #[5, 6] } "both ends pinned: BothPinned"
  fs := add fs (opCase "rs-tetra-collapse-link" "collapse" tet 0 1) "collapse on a tetrahedron: the result would repeat a face (NonManifoldResult)"
  fs := add fs (opCase "rs-strip-collapse-link" "collapse" ((grid 1 3 0.3).jitter 0.05 108) 2 3) "rung of a one-quad-wide strip joins two boundary vertices: LinkCondition"
  fs := add fs (opCase "rs-icosa-flip" "flip" ico 0 1) "flip an interior edge of the icosahedron (no legality checks)"
  fs := add fs { opCase "rs-icosa-flip-delaunay" "flip" ico 0 1 with delaunay := true } "same flip with the Delaunay check: result depends on the opposite angles"
  fs := add fs (opCase "rs-patch-flip-boundary" "flip" pat 0 1) "boundary edge: BoundaryEdge"
  fs := add fs (opCase "rs-mobius-flip" "flip" mob 0 1) "Moebius band seam edge: its two faces traverse it the same way, InconsistentOrientation"
  fs := add fs (opCase "rs-mobius-flip-ok" "flip" mob 0 3) "Moebius band diagonal edge: an ordinary flip"
  fs := add fs (opCase "rs-icosa-split-pinned-note" "split" ico 4 9) "split (pins play no role in a split)"
  fs := add fs { opCase "rs-torus-flip-foldover" "flip" tor 0 1 with foldover := true } "torus flip with the foldover check"
  fs := add fs (opCase "rs-error-no-edge" "split" ico 0 3) "vertices 0 and 3 of the icosahedron share no edge: NoSuchEdge"
  -- errors
  fs := add fs (mk "rs-error-degenerate" "energy" { pos := #[⟨0.0, 0.0, 0.0⟩, ⟨1.0, 0.0, 0.0⟩, ⟨2.0, 0.0, 0.0⟩, ⟨0.0, 1.0, 0.0⟩], faces := #[⟨0, 1, 2⟩, ⟨0, 2, 3⟩] }) "a face with collinear vertices: DegenerateFace"
  fs := add fs (mk "rs-error-nonmanifold" "energy" { pos := #[⟨0.0, 0.0, 0.0⟩, ⟨1.0, 0.0, 0.0⟩, ⟨0.0, 1.0, 0.0⟩, ⟨0.0, 0.0, 1.0⟩, ⟨0.0, -1.0, 0.0⟩], faces := #[⟨0, 1, 2⟩, ⟨0, 1, 3⟩, ⟨0, 1, 4⟩] }) "three faces on one edge: NonManifoldEdge"
  fs := add fs (mk "rs-error-unreferenced" "energy" { pos := tet.pos.push ⟨5.0, 5.0, 5.0⟩, faces := tet.faces }) "a vertex without faces: UnreferencedVertex"
  fs := add fs (mk "rs-error-range" "energy" { tet with faces := tet.faces.set! 2 ⟨0, 2, 9⟩ }) "face index 9: VertexOutOfRange"
  fs := add fs { mk "rs-error-exponents" "energy" tet with ex := ⟨6.0, 8.0⟩ } "beta = alpha + 2 is not allowed: InvalidExponents"
  fs := add fs { mk "rs-error-same-barycenter" "energy" { pos := #[⟨-1.0, -1.0, 0.0⟩, ⟨1.0, -1.0, 0.0⟩, ⟨0.0, 2.0, 0.0⟩, ⟨1.0, 1.0, 0.5⟩, ⟨-1.0, 1.0, 0.5⟩, ⟨0.0, -2.0, -1.0⟩], faces := #[⟨0, 1, 2⟩, ⟨3, 4, 5⟩] } with ex := ex612 } "two faces with the same barycenter: NonFiniteEnergy"
  fs := add fs { mk "rs-error-volume-mobius" "constraints" mob with cons := #[.volume none none] } "volume on a non-orientable band: NonOrientable"
  fs := add fs { mk "rs-error-translation" "saddle" ico with cons := #[.area none (some (aOf ico * 1.03))] } "no barycenter or pin: TranslationNotFixed"
  fs := add fs { mk "rs-error-singular" "saddle" ico with cons := #[.pins #[0] none, .pins #[0] none] } "the same pin twice: SingularSaddle"
  fs := add fs { mk "rs-error-pin-range" "constraints" ico with cons := #[.pins #[99] none] } "pin at vertex 99: ConstraintOutOfRange"
  fs := add fs { mk "rs-error-component" "constraints" tet with cons := #[.area (some 3) none] } "component 3 of a connected mesh: ConstraintOutOfRange"
  return fs

/-- Does the printed input parse back to a case with the same output? -/
def roundTrips (c : RSCase) : Bool :=
  let lines := ((printCase c).splitOn "\n").toArray.filterMap fun l =>
    let t := tokens l
    if t.isEmpty || t[0]! == "case" || t[0]! == "end" then none else some t
  match parseBody c.name lines with
  | .error _ => false
  | .ok c2 => runCase c2 == runCase c && printCase c2 == printCase c

private def statusOf (output : String) : String :=
  match (output.splitOn "\n").find? (·.startsWith "status ") with
  | some l => (tokens l)[1]!
  | none => "?"

/-- Compare two outputs token by token: floats to `tol` (relative to `max(|x|, 1)`), NaN equal to NaN. -/
def compareRS (expected actual : String) (tol : Float := 1e-9) : Option String := Id.run do
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
        else if let some d := compareRS stored text then
          problems := problems.push s!"{name}{suffix}: {d}"
    if problems.size == 0 then IO.println s!"ok    {name}"
    else
      bad := bad + 1
      for p in problems do IO.println s!"FAIL  {p}"
  for (c, _) in fixtures do
    if !roundTrips c then
      bad := bad + 1
      IO.println s!"FAIL  {c.name}: the printed input does not parse back to the same case"
  IO.println s!"rs fixtures: {fixtures.size} cases, {bad} problems"
  return bad

end Ddg.RepulsiveSurfaces
