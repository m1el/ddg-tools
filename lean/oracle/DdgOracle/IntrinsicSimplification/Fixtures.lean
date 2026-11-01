import DdgOracle.IntrinsicSimplification.Selftest

/-!
# Golden fixtures of the `iem` oracle

`fixtures/iem/NAME.in|out` and `MANIFEST.tsv` (name, status of the run, what the case is). Checking
recomputes both files and compares them byte for byte (the floats are bit patterns; the run is
deterministic on one machine, `exp`/`sin`/`cos`/`atan`/`log` being the C library's).
-/

namespace Ddg.IEM

private def statusOf (output : String) : String :=
  match ((output.splitOn "\n").find? (·.startsWith "flatten ")) <|>
      ((output.splitOn "\n").find? (·.startsWith "solve ")) with
  | some l => (toks l)[1]!
  | none => match (output.splitOn "\n").find? (·.startsWith "status ") with
    | some l => (toks l)[1]!
    | none => "?"

/-- Write the fixtures. -/
def writeFixtures (dir : System.FilePath) : IO Unit := do
  IO.FS.createDirAll dir
  let mut manifest := "name\tstatus\twhat\n"
  for (c, what) in allCases do
    let output := runCase c
    IO.FS.writeFile (dir / s!"{c.name}.in") (printCase c)
    IO.FS.writeFile (dir / s!"{c.name}.out") output
    manifest := manifest ++ s!"{c.name}\t{statusOf output}\t{what}\n"
    IO.println s!"{c.name}: {statusOf output}"
  IO.FS.writeFile (dir / "MANIFEST.tsv") manifest
  IO.println s!"wrote {allCases.size} fixtures to {dir}"

/-- Check the fixtures (exact equality); returns the number of differences. -/
def checkFixtures (dir : System.FilePath) : IO Nat := do
  let mut bad := 0
  for (c, _) in allCases do
    let mut problems : Array String := #[]
    for (suffix, text) in [(".in", printCase c), (".out", runCase c)] do
      let path := dir / s!"{c.name}{suffix}"
      if !(← path.pathExists) then problems := problems.push s!"{c.name}{suffix} is missing"
      else if (← IO.FS.readFile path) != text then
        problems := problems.push s!"{c.name}{suffix} differs from the recomputed text"
    if problems.size == 0 then IO.println s!"ok    {c.name}"
    else
      bad := bad + 1
      for p in problems do IO.println s!"FAIL  {p}"
  IO.println s!"iem fixtures: {allCases.size - bad} of {allCases.size} match"
  return bad

end Ddg.IEM
