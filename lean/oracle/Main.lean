import DdgOracle
open Ddg

def usage : String :=
"ddg-oracle: dense Float reference implementation of
  'Optimal Cone Singularities for Conformal Flattening' (Soliman, Slepcev, Crane 2018)

usage:
  ddg-oracle operators [IN [OUT]]   curvatures, masses and cotan Laplacian of every case
  ddg-oracle solve     [IN [OUT]]   solve the canonical problem for every case
  ddg-oracle batch     [IN [OUT]]   like the above; each case names its task in a 'task' line
  ddg-oracle selftest  [N [SIZE]]   internal consistency checks, with a sweep over N random cases
                                    (default 300) of up to about SIZE vertices (default 41)
  ddg-oracle fixtures write DIR     regenerate the golden fixtures
  ddg-oracle fixtures check DIR     recompute the fixtures and compare with the files (tolerance 1e-9)
  ddg-oracle bench                  throughput on meshes of about 50 and 200 vertices
  ddg-oracle cuts      [IN [OUT]]   one step of the cut flow of 'Variational Surface Cutting' (Sharp, Crane 2018)
  ddg-oracle cuts-selftest          internal checks of the cut flow
  ddg-oracle cuts-fixtures write DIR | check DIR   golden fixtures of the cut flow (fixtures/cuts)
  ddg-oracle cuts-default-fixtures write DIR | check DIR   fixtures of the default stages (fixtures/cuts-default)
  ddg-oracle regions   [IN [OUT]]   one step of the cut flow for N regions (SPEC section 8)
  ddg-oracle regions-selftest       internal checks of the N-region flow
  ddg-oracle regions-fixtures write DIR | check DIR   golden fixtures (fixtures/regions)
  ddg-oracle cuts-robust-fixtures write DIR | check DIR   fixtures of the robust statistics (fixtures/cuts-robust)
  ddg-oracle regions-robust-fixtures write DIR | check DIR   same for N regions (fixtures/regions-robust)
  ddg-oracle intrinsic [IN [OUT]]   single operations of 'Navigating Intrinsic Triangulations'
                                    (flip, delaunay, trace, edgetrace, circumcenter, move)
  ddg-oracle intrinsic-selftest     internal checks of the intrinsic operations
  ddg-oracle intrinsic-fixtures write DIR | check DIR   golden fixtures (fixtures/intrinsic)
  ddg-oracle dev [IN [OUT]]         energies, one L-BFGS step and 4-1 subdivision of 'Developability of Triangle Meshes'
                                    (Stein, Grinspun, Crane 2018); tasks energy, step, subdivide
  ddg-oracle dev-gen SEED COUNT [FIRST]   write COUNT seeded random cases for `dev` to standard output
  ddg-oracle dev-selftest           internal checks of the developability oracle
  ddg-oracle dev-fixtures write DIR | check DIR   golden fixtures (fixtures/dev-*)
  ddg-oracle rc [IN [OUT]]          projected Sobolev gradient and chord projection of 'Repulsive Curves'
                                    (Yu, Schumacher, Crane 2021); tasks energy, gram, saddle, project
  ddg-oracle rc-gen SEED COUNT [FIRST]   write COUNT seeded random cases for `rc` to standard output
  ddg-oracle rc-selftest [N]        internal checks of the repulsive-curves oracle (N seeded networks, default 120)
  ddg-oracle rc-fixtures write DIR | check DIR   golden fixtures (fixtures/rc)
  ddg-oracle rs [IN [OUT]]          dense tangent-point energy, metric, constraints, saddle, preconditioner, one step and
                                    remeshing operations of 'Repulsive Surfaces' (Yu, Brakensiek, Schumacher, Crane 2021);
                                    tasks energy, df, metric, constraints, saddle, precond, step, remesh-ops
  ddg-oracle rs-gen SEED COUNT [FIRST]   write COUNT seeded random cases for `rs` to standard output
  ddg-oracle rs-selftest [N]        internal checks of the repulsive-surfaces oracle (N seeded meshes, default 120)
  ddg-oracle rs-fixtures write DIR | check DIR   golden fixtures (fixtures/rs)
  ddg-oracle ic [IN [OUT]]          integer formulas of 'Integer Coordinates for Intrinsic Geometry Processing'
                                    (Gillespie, Sharp, Crane 2021); tasks corner, flip, splitedge, trace,
                                    tracecurve, identify, subdivide, regions, split
  ddg-oracle ic-selftest            internal checks of the integer-coordinates oracle
  ddg-oracle ic-fixtures write DIR | check DIR   golden fixtures (fixtures/ic, exact equality)
  ddg-oracle iem [IN [OUT]]         one star of 'Surface Simplification using Intrinsic Error Metrics'
                                    (Liu et al. 2023); tasks step, flatten, weights, removal
  ddg-oracle iem-selftest           internal checks of the intrinsic-simplification oracle
  ddg-oracle iem-fixtures write DIR | check DIR   golden fixtures (fixtures/iem)
  ddg-oracle pat [IN [OUT]]         features, torus fitting, per-torus SDF, blend and SHC/SSPD baselines of 'Points as Tori'
                                    (Feng, Crane 2025); tasks features, fit, sdf, blend, baselines
  ddg-oracle pat-selftest [N]       internal checks of the points-as-tori oracle (N seeded cases, default 60)
  ddg-oracle pat-fixtures write DIR | check DIR   golden fixtures (fixtures/pat)
  ddg-oracle version

IN and OUT default to standard input and output ('-'). All three processing commands read any
number of cases from one stream and flush after every case, so one process can serve a fuzzer.
The text format is specified in README.md."

def openInput (path : String) : IO IO.FS.Stream := do
  if path == "-" then IO.getStdin
  else return IO.FS.Stream.ofHandle (← IO.FS.Handle.mk path .read)

def openOutput (path : String) : IO IO.FS.Stream := do
  if path == "-" then IO.getStdout
  else return IO.FS.Stream.ofHandle (← IO.FS.Handle.mk path .write)

def process (fixed : Option Task) (rest : List String) : IO UInt32 := do
  if rest.length > 2 then
    IO.eprintln usage
    return 2
  let input ← openInput (rest.getD 0 "-")
  let output ← openOutput (rest.getD 1 "-")
  let _ ← runStream input output fixed
  output.flush
  return 0

def main (args : List String) : IO UInt32 := do
  match args with
  | "operators" :: rest => process (some .operators) rest
  | "solve" :: rest => process (some .solve) rest
  | "batch" :: rest => process none rest
  | ["selftest"] => return if (← selftest) == 0 then 0 else 1
  | ["selftest", n] =>
    match n.toNat? with
    | some k => return if (← selftest k) == 0 then 0 else 1
    | none => IO.eprintln usage; return 2
  | ["selftest", n, size] =>
    match n.toNat?, size.toNat? with
    | some k, some sz => return if (← selftest k sz) == 0 then 0 else 1
    | _, _ => IO.eprintln usage; return 2
  | ["fixtures", "write", dir] => writeFixtures dir; return 0
  | ["fixtures", "check", dir] => return if (← checkFixtures dir) == 0 then 0 else 1
  | "cuts" :: rest =>
    if rest.length > 2 then
      IO.eprintln usage
      return 2
    let input ← openInput (rest.getD 0 "-")
    let output ← openOutput (rest.getD 1 "-")
    let _ ← Cuts.runStream input output
    output.flush
    return 0
  | ["cuts-selftest"] => return if (← Cuts.selftest) + (← Cuts.selftestRobust) == 0 then 0 else 1
  | ["cuts-fixtures", "write", dir] => Cuts.writeFixtures dir; return 0
  | ["cuts-fixtures", "check", dir] => return if (← Cuts.checkFixtures dir) == 0 then 0 else 1
  | ["cuts-default-fixtures", "write", dir] => Cuts.writeList Cuts.defaultFixtures dir; return 0
  | ["cuts-default-fixtures", "check", dir] =>
    return if (← Cuts.checkList Cuts.defaultFixtures dir) == 0 then 0 else 1
  | "regions" :: rest =>
    if rest.length > 2 then
      IO.eprintln usage
      return 2
    let input ← openInput (rest.getD 0 "-")
    let output ← openOutput (rest.getD 1 "-")
    let _ ← Regions.runStream input output
    output.flush
    return 0
  | ["regions-fixtures", "write", dir] => Regions.writeFixtures dir; return 0
  | ["regions-fixtures", "check", dir] => return if (← Regions.checkFixtures dir) == 0 then 0 else 1
  | ["regions-selftest"] => return if (← Regions.selftest) + (← Regions.selftestRobust) == 0 then 0 else 1
  | ["cuts-robust-fixtures", "write", dir] => Cuts.writeList Cuts.robustFixtures dir; return 0
  | ["cuts-robust-fixtures", "check", dir] =>
    return if (← Cuts.checkList Cuts.robustFixtures dir) == 0 then 0 else 1
  | ["regions-robust-fixtures", "write", dir] => Regions.writeRobust dir; return 0
  | ["regions-robust-fixtures", "check", dir] => return if (← Regions.checkRobust dir) == 0 then 0 else 1
  | "intrinsic" :: rest =>
    if rest.length > 2 then
      IO.eprintln usage
      return 2
    let input ← openInput (rest.getD 0 "-")
    let output ← openOutput (rest.getD 1 "-")
    let _ ← Intrinsic.runStream input output
    output.flush
    return 0
  | ["intrinsic-selftest"] => return if (← Intrinsic.selftest) == 0 then 0 else 1
  | ["intrinsic-fixtures", "write", dir] => Intrinsic.writeFixtures dir; return 0
  | ["intrinsic-fixtures", "check", dir] => return if (← Intrinsic.checkFixtures dir) == 0 then 0 else 1
  | "dev" :: rest =>
    if rest.length > 2 then
      IO.eprintln usage
      return 2
    let input ← openInput (rest.getD 0 "-")
    let output ← openOutput (rest.getD 1 "-")
    let _ ← Developability.runStream input output
    output.flush
    return 0
  | "dev-gen" :: seed :: count :: rest =>
    match seed.toNat?, count.toNat?, (rest.getD 0 "0").toNat? with
    | some s, some n, some first =>
      let out ← IO.getStdout
      for i in [first:first + n] do
        out.putStr (Developability.printCase (Developability.genCase s i))
      out.flush
      return 0
    | _, _, _ => IO.eprintln usage; return 2
  | ["dev-selftest"] => return if (← Developability.selftest) == 0 then 0 else 1
  | ["dev-fixtures", "write", dir] => Developability.writeFixtures dir; return 0
  | ["dev-fixtures", "check", dir] => return if (← Developability.checkFixtures dir) == 0 then 0 else 1
  | "rc" :: rest =>
    if rest.length > 2 then
      IO.eprintln usage
      return 2
    let input ← openInput (rest.getD 0 "-")
    let output ← openOutput (rest.getD 1 "-")
    let _ ← RepulsiveCurves.runStream input output
    output.flush
    return 0
  | "rc-gen" :: seed :: count :: rest =>
    match seed.toNat?, count.toNat?, (rest.getD 0 "0").toNat? with
    | some s, some n, some first =>
      let out ← IO.getStdout
      for i in [first:first + n] do
        out.putStr (RepulsiveCurves.printCase (RepulsiveCurves.genCase s i))
      out.flush
      return 0
    | _, _, _ => IO.eprintln usage; return 2
  | ["rc-selftest"] =>
    let mut bad ← RepulsiveCurves.selftest
    for i in [0:300] do
      if !RepulsiveCurves.roundTrips (RepulsiveCurves.genCase 5 i) then
        IO.println s!"FAIL  generated case {i} does not round-trip"
        bad := bad + 1
    return if bad == 0 then 0 else 1
  | ["rc-selftest", n] =>
    match n.toNat? with
    | some k => return if (← RepulsiveCurves.selftest k) == 0 then 0 else 1
    | none => IO.eprintln usage; return 2
  | ["rc-fixtures", "write", dir] => RepulsiveCurves.writeFixtures dir; return 0
  | ["rc-fixtures", "check", dir] => return if (← RepulsiveCurves.checkFixtures dir) == 0 then 0 else 1
  | "rs" :: rest =>
    if rest.length > 2 then
      IO.eprintln usage
      return 2
    let input ← openInput (rest.getD 0 "-")
    let output ← openOutput (rest.getD 1 "-")
    let _ ← RepulsiveSurfaces.runStream input output
    output.flush
    return 0
  | "rs-gen" :: seed :: count :: rest =>
    match seed.toNat?, count.toNat?, (rest.getD 0 "0").toNat? with
    | some s, some n, some first =>
      let out ← IO.getStdout
      for i in [first:first + n] do
        out.putStr (RepulsiveSurfaces.printCase (RepulsiveSurfaces.genCase s i))
      out.flush
      return 0
    | _, _, _ => IO.eprintln usage; return 2
  | ["rs-selftest"] => return if (← RepulsiveSurfaces.selftest) == 0 then 0 else 1
  | ["rs-selftest", n] =>
    match n.toNat? with
    | some k => return if (← RepulsiveSurfaces.selftest k) == 0 then 0 else 1
    | none => IO.eprintln usage; return 2
  | ["rs-fixtures", "write", dir] => RepulsiveSurfaces.writeFixtures dir; return 0
  | ["rs-fixtures", "check", dir] => return if (← RepulsiveSurfaces.checkFixtures dir) == 0 then 0 else 1
  | "pat" :: rest =>
    if rest.length > 2 then
      IO.eprintln usage
      return 2
    let input ← openInput (rest.getD 0 "-")
    let output ← openOutput (rest.getD 1 "-")
    let _ ← PointsAsTori.runStream input output
    output.flush
    return 0
  | ["pat-selftest"] => return if (← PointsAsTori.selftest) == 0 then 0 else 1
  | ["pat-selftest", n] =>
    match n.toNat? with
    | some k => return if (← PointsAsTori.selftest k) == 0 then 0 else 1
    | none => IO.eprintln usage; return 2
  | ["pat-fixtures", "write", dir] => PointsAsTori.writeFixtures dir; return 0
  | ["pat-fixtures", "check", dir] => return if (← PointsAsTori.checkFixtures dir) == 0 then 0 else 1
  | "ic" :: rest =>
    if rest.length > 2 then
      IO.eprintln usage
      return 2
    let input ← openInput (rest.getD 0 "-")
    let output ← openOutput (rest.getD 1 "-")
    let _ ← IntegerCoords.runStream input output
    output.flush
    return 0
  | ["ic-selftest"] => return if (← IntegerCoords.selftest) == 0 then 0 else 1
  | ["ic-fixtures", "write", dir] => IntegerCoords.writeFixtures dir; return 0
  | ["ic-fixtures", "check", dir] => return if (← IntegerCoords.checkFixtures dir) == 0 then 0 else 1
  | "iem" :: rest =>
    if rest.length > 2 then
      IO.eprintln usage
      return 2
    let input ← openInput (rest.getD 0 "-")
    let output ← openOutput (rest.getD 1 "-")
    let _ ← IEM.runStream input output
    output.flush
    return 0
  | ["iem-selftest"] => return if (← IEM.selftest) == 0 then 0 else 1
  | ["iem-fixtures", "write", dir] => IEM.writeFixtures dir; return 0
  | ["iem-fixtures", "check", dir] => return if (← IEM.checkFixtures dir) == 0 then 0 else 1
  | ["bench"] => bench; return 0
  | ["version"] => IO.println "ddg-oracle 0.1.0 (text format 1)"; return 0
  | _ => IO.eprintln usage; return 2
