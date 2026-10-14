import DdgOracle.Fixtures
import DdgOracle.Intrinsic.Selftest

/-!
# Golden fixtures of the intrinsic operations

`fixtures/intrinsic/NAME.in|out` and `MANIFEST.tsv`.
-/

namespace Ddg.Intrinsic

open Ddg Ddg.Cuts

/-- The fixtures: case and description. -/
def fixtures : Array (ICase × String) := Id.run do
  let warp := fun (m : RawMesh) => m.mapPos fun p =>
    ⟨p.x * (1 + 0.05 * Float.sin (5 * p.y + p.z)), p.y * (1 + 0.04 * Float.cos (4 * p.x - 2 * p.z)), p.z⟩
  let ico1 := warp (icosphere 1)
  let ico2n := (icosphere 2).jitter 0.05 3
  let grid := sheared 9 8
  let annulus := (polar 4 14 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.2 * y * y).jitter 0.03 7
  let tube := (twistedTorus 12 3 1.5 0.03).jitter 0.001 3
  let twisted := (twistedTorus 36 3 2.0 0.06).jitter 0.002 17
  let mut fs : Array (ICase × String) := #[]
  let add := fun (fs : Array (ICase × String)) (c : ICase) (what : String) => fs.push (c, what)
  -- flip
  let cf1 := mkICase "it-flip-ico1" "flip" ico1
  fs := add fs { cf1 with flip := firstEdgeWith cf1 none, again := true } "flip one edge of an icosphere (level 1), then flip the new edge back"
  let cf2 := mkICase "it-flip-grid" "flip" grid
  fs := add fs { cf2 with flip := firstEdgeWith cf2 none, again := true } "flip an interior edge of a sheared grid with boundary, then again"
  fs := add fs { cf2 with name := "it-flip-refused-boundary", flip := firstEdgeWith cf2 (some .boundary) } "a boundary edge: refused"
  fs := add fs { cf2 with name := "it-flip-refused-not-convex", flip := firstEdgeWith cf2 (some .notConvex) } "an edge whose two triangles form a non-convex quadrilateral: refused"
  let cf3 := mkICase "it-flip-tube" "flip" tube
  fs := add fs { cf3 with flip := firstEdgeWith cf3 none, again := true } "flip an edge of the thin tube torus"
  -- Delaunay
  fs := add fs (mkICase "it-delaunay-grid" "delaunay" grid) "Delaunay flipping of a sheared grid with boundary"
  fs := add fs (mkICase "it-delaunay-ico2" "delaunay" ico2n) "Delaunay flipping of a noisy icosphere (level 2)"
  fs := add fs (mkICase "it-delaunay-annulus" "delaunay" annulus) "Delaunay flipping of an annulus (two boundary loops)"
  fs := add fs (mkICase "it-delaunay-tube" "delaunay" tube) "Delaunay flipping of a thin tube torus: the result has edges from a vertex to itself and double edges"
  fs := add fs (mkICase "it-delaunay-twisted" "delaunay" twisted) "Delaunay flipping of a long twisted torus"
  fs := add fs { mkICase "it-delaunay-grid-cotan" "delaunay" grid with crit := .cotan } "the reference's criterion (cotan weight < -1e-6) on the sheared grid"
  -- trace
  let tr2 := mkICase "it-trace-ico2" "trace" ico2n
  fs := add fs { tr2 with traces := #[(0, 0.3, 0.35), (7, 1.7, 0.6), (40, 4.0, 0.25), (100, 5.5, 0.9), (150, 2.9, 0.05)] } "traces from vertices of a noisy icosphere (several crossings)"
  let trg := mkICase "it-trace-grid" "trace" grid
  fs := add fs { trg with traces := #[(10, 0.4, 0.5), (20, 1.0, 3.0), (31, 0.2, 0.4), (40, 3.0, 0.6)] } "traces on a sheared grid with boundary: ending inside, at the boundary"
  -- trace of edges
  fs := add fs { mkICase "it-edgetrace-ico2" "edgetrace" ico2n with predelaunay := true } "after Delaunay flipping, every edge traced over the input mesh for its length"
  fs := add fs { mkICase "it-edgetrace-tube" "edgetrace" tube with predelaunay := true } "the same on the thin tube torus (self-edges, double edges)"
  let ce := mkICase "it-edgetrace-grid" "edgetrace" grid
  fs := add fs { ce with preflips := (match triOf ce with
      | none => #[]
      | some t => Id.run do
          let mut out : Array (Nat × Nat) := #[]
          let mut tt := t
          for h in [0:t.nh] do
            if out.size < 3 && t.eid[h]! == h then
              match flipHalfedge tt h with
              | .ok (t', _) =>
                out := out.push (tt.tail[h]!, tt.head h)
                tt := t'
              | .error _ => pure ()
          return out) } "three flips by hand, then every edge traced over the input mesh"
  -- circumcenter
  let mild := ((gridDisk 10 9 fun x y => 0.1 * Float.sin (4 * x) * y).mapPos fun p => ⟨p.x + 0.6 * p.y, p.y * 0.8, p.z⟩).jitter 0.04 9
  fs := add fs { mkICase "it-circumcenter-grid" "circumcenter" mild with circAll := true } "circumcenters of all faces of a sheared grid (many outside the face), traced from the barycentre"
  fs := add fs { mkICase "it-circumcenter-ico2" "circumcenter" ico2n with circAll := true } "circumcenters of all faces of a noisy icosphere (some outside their face, traced across edges)"
  -- moving an inserted vertex
  return fs ++ moveFixtures

/-- Write the fixtures. -/
def writeFixtures (dir : System.FilePath) : IO Unit := do
  IO.FS.createDirAll dir
  let mut manifest := ""
  for (c, what) in fixtures do
    let input := printCase c
    let output := runCase c
    IO.FS.writeFile (dir / s!"{c.name}.in") input
    IO.FS.writeFile (dir / s!"{c.name}.out") output
    manifest := manifest ++ s!"{c.name}\t{statusOf output}\t{what}\n"
    IO.println s!"{c.name}: {statusOf output}"
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
        else if let some d := compareOutputs stored text then
          problems := problems.push s!"{name}{suffix}: {d}"
    if problems.size == 0 then IO.println s!"ok    {name}"
    else
      bad := bad + 1
      for p in problems do IO.println s!"FAIL  {p}"
  IO.println s!"intrinsic fixtures: {fixtures.size - bad} of {fixtures.size} match"
  return bad

end Ddg.Intrinsic
