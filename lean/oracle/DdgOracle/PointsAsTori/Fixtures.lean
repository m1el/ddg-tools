import DdgOracle.PointsAsTori.Selftest

/-!
# Golden fixtures of the `pat` oracle

`fixtures/pat/NAME.in|out` and `fixtures/pat/MANIFEST.tsv`.
-/

namespace Ddg.PointsAsTori

private def sphereCoef : Array Float := #[0.0, 0.0, 0.0, 0.0, -0.5, -0.5]

/-- Random cloud: points near the ellipsoid `(1.0, 0.8, 0.6)`, normals roughly outward. -/
private def ellipsoidCloud (n : Nat) : G (Array V3 × Array V3) := do
  let mut ps := #[]
  let mut ns := #[]
  for _ in [0:n] do
    let u ← unitVec
    let j ← uni 0.9 1.1
    ps := ps.push ⟨u.x * j, u.y * 0.8 * j, u.z * 0.6 * j⟩
    let g : V3 := ⟨u.x, u.y / 0.8, u.z / 0.6⟩
    ns := ns.push (g.sdiv g.norm)
  return (ps, ns)

private def run {α : Type} [Inhabited α] (seed : Nat) (g : G α) : α := (g.run (mkRng seed)).1

private def randCoefs (n : Nat) : G (Array (Array Float)) := do
  let mut a := #[]
  for _ in [0:n] do a := a.push (← randCoef)
  return a

private def randAreas (n : Nat) : G (Array Float) := do
  let mut a := #[]
  for _ in [0:n] do a := a.push (← uni 0.05 0.12)
  return a

private def normalize (v : V3) : V3 := v.sdiv v.norm

private def withCoefs (c : PatCase) (cs : Array (Array Float)) : PatCase :=
  { c with coefs := (Array.range cs.size).map fun i => (i, cs[i]!) }

private def spreadQueries (extra : Array V3) (n : Nat) (lo hi : Float) (seed : Nat) : Array V3 :=
  extra ++ run seed (randQueries n lo hi)

/-- The fixtures: case and description. -/
def fixtures : Array (PatCase × String) := Id.run do
  let mut fs : Array (PatCase × String) := #[]
  -- ---- features
  let (fp, fn) := run 11 (ellipsoidCloud 12)
  fs := fs.push ({ cloudCase "pat-features-even" "features" fp fn with center := 0, nbr := some #[0, 3, 1, 7, 2, 5], y := some #[-0.3, 0.1, -0.2, 0.7, -0.5, 0.4] }, "6 explicit neighbors (even median), network output y to coefficients")
  let (fp2, fn2) := run 12 (ellipsoidCloud 40)
  let fn2 := fn2.set! 0 (normalize ⟨0.95, 0.1, -0.2⟩)
  fs := fs.push ({ cloudCase "pat-features-frame-branch" "features" fp2 fn2 with center := 0, nbhdK := some 9, y := some #[0.2, -0.1, 0.3, 0.1, -0.6, -0.4] }, "normal with |c_x| > 0.9 (second branch of the frame rule), auto neighborhood of 9")
  let (fp3, fn3) := fibonacci 100 1.0 V3.zero
  fs := fs.push ({ cloudCase "pat-features-k64" "features" fp3 fn3 with center := 3, nbhdK := some 64, y := some #[-1.0, 0.0, 0.0, 0.0, -0.5, -0.5] }, "100-point sphere, 64 tokens as in the network")
  let dupP := #[(⟨0.0, 0.0, 0.0⟩ : V3), ⟨0.0, 0.0, 0.0⟩, ⟨1.0, 0.0, 0.0⟩, ⟨0.0, 2.0, 0.0⟩, ⟨0.0, 0.0, 3.0⟩, ⟨1.0, 1.0, 1.0⟩, ⟨-1.0, 0.5, 0.0⟩]
  let dupN := dupP.map fun _ => (⟨0.0, 0.0, 1.0⟩ : V3)
  fs := fs.push ({ cloudCase "pat-features-duplicate" "features" dupP dupN with center := 1, nbhdK := some 5 },
    "exact duplicate of the center: self first, then the duplicate (distance 0) by index")
  fs := fs.push ({ cloudCase "pat-features-fewer-than-k" "features" (dupP.extract 2 7) (dupN.extract 2 7) with center := 0, nbhdK := some 64 },
    "5 points, k = 64: the neighborhood has 5 tokens")
  fs := fs.push ({ cloudCase "pat-error-features-sigma" "features" #[(⟨1.0, 1.0, 1.0⟩ : V3), ⟨1.0, 1.0, 1.0⟩, ⟨1.0, 1.0, 1.0⟩] #[ez3, ez3, ez3] with center := 0, nbhdK := some 3 },
    "all points coincide: DegenerateNeighborhood")
  fs := fs.push ({ cloudCase "pat-error-nonunit-normal" "features" #[(⟨0.0, 0.0, 0.0⟩ : V3), ⟨1.0, 0.0, 0.0⟩] #[ez3, (⟨0.0, 0.0, 2.0⟩ : V3)] with center := 0, nbhdK := some 2 },
    "normal of length 2: NonUnitNormal")
  -- ---- fit
  let fitPts : Array V3 := #[⟨0.3, -0.2, 0.5⟩, ⟨-0.4, 0.1, 0.2⟩, ⟨0.0, 0.0, 0.0⟩, ⟨1.0, 2.0, 3.0⟩, ⟨-1.0, 0.5, -0.3⟩, ⟨0.2, 0.2, -0.7⟩]
  let fitNrm : Array V3 := #[ez3, normalize ⟨0.95, 0.2, -0.1⟩, normalize ⟨-0.97, 0.1, 0.1⟩, normalize ⟨1.0, 1.0, 1.0⟩, ey3, normalize ⟨0.0, -0.3, 0.9⟩]
  let fitBase := cloudCase "x" "fit" fitPts fitNrm
  let six := fun (a : Array Float) => (Array.range 6).map fun _ => a
  fs := fs.push (withCoefs { fitBase with name := "pat-fit-sphere" } (six sphereCoef), "tangent-sphere initialization a02 = a20 = -0.5 at six points (both frame branches): R = 0, r = 1, center p - n")
  fs := fs.push (withCoefs { fitBase with name := "pat-fit-sphere-concave" } #[#[0.0, 0.0, 0.0, 0.0, 0.5, 0.5], #[0.02, 0.0, 0.0, 0.0, 0.5, 0.5], #[0.0, 0.0, 0.0, 0.0, -0.25, -0.25], #[-0.01, 0.0, 0.0, 0.0, 1.5, 1.5], #[0.0, 0.0, 0.0, 0.0, 0.1, 0.1], #[0.0, 0.0, 0.0, 0.0, -2.0, -2.0]],
    "umbilic points, positive curvature (sign -1) and other radii")
  fs := fs.push (withCoefs { fitBase with name := "pat-fit-saddle" } #[#[0.0, 0.0, 0.0, 0.0, -0.3, 0.4], #[0.0, -0.1, 0.2, 0.1, -0.3, 0.4], #[0.01, 0.15, -0.2, 0.3, 0.5, -0.8], #[0.0, 0.0, 0.0, 0.6, 0.0, 0.0], #[0.0, 0.05, 0.05, -0.2, 0.7, -0.2], #[-0.02, 0.2, 0.1, 0.0, -0.9, 0.35]],
    "saddles (kappa+ kappa- < 0), tilted tangent planes, pure a11")
  fs := fs.push (withCoefs { fitBase with name := "pat-fit-ellipsoid-like" } #[#[0.0, 0.0, 0.0, 0.0, -0.5, -0.9], #[0.0, 0.1, -0.1, 0.05, -0.4, -0.7], #[0.0, 0.0, 0.0, 0.0, 0.8, 0.3], #[0.01, 0.0, 0.0, 0.2, -0.6, -0.4], #[0.0, -0.2, 0.25, 0.0, 0.3, 0.9], #[0.0, 0.0, 0.0, 0.0, -0.15, -1.2]],
    "elliptic points: kappa+ kappa- > 0, spindle (|kmin| > |kmax|/2) and ring cases")
  fs := fs.push (withCoefs { fitBase with name := "pat-fit-cylinder" } #[#[0.0, 0.0, 0.0, 0.0, 0.0, -0.5], #[0.0, 0.0, 0.0, 0.0, -0.25, 0.0], #[0.0, 0.0, 0.0, 0.0, 0.0, 0.5], #[0.0, 0.0, 0.3, 0.0, 0.0, -0.5], #[0.0, 0.0, 0.0, 0.0, 1.0e-12, -0.5], #[0.02, 0.0, 0.0, 0.0, 0.0, 0.0]],
    "cylinder-like (one curvature exactly 0 or 1e-12): curvature floor, huge major radius; a plane with offset")
  fs := fs.push (withCoefs { fitBase with name := "pat-fit-nearplanar" } #[#[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], #[0.0, 0.0, 0.0, 1.0e-13, 1.0e-12, -1.0e-12], #[0.0, 0.0, 0.0, 0.0, 3.0e-11, 3.0e-11], #[0.0, 0.0, 0.1, 0.0, 0.0, 0.0], #[0.0, 0.0, 0.0, 0.0, -1.0e-11, 2.0e-11], #[0.0, 1.0e-9, -1.0e-9, 2.0e-9, 5.0e-9, -4.0e-9]],
    "planar and near-planar coefficients: disc clamp, eigenvector clamps, curvature floor all active")
  fs := fs.push ({ withCoefs { fitBase with name := "pat-fit-nearplanar-eps" } #[#[0.0, 0.0, 0.0, 1.0e-7, 1.0e-6, -1.0e-6], #[0.0, 0.0, 0.0, 0.0, 3.0e-5, 3.0e-5], #[0.0, 0.0, 0.0, 0.0, 0.0, -2.0e-5]] with fitOpts := { eps := 1.0e-6, floor := 1.0e-4 } },
    "custom clamp options eps = 1e-6, floor = 1e-4")
  let (rp, rn) := run 21 (ellipsoidCloud 48)
  fs := fs.push (withCoefs (cloudCase "pat-fit-random" "fit" rp rn) (run 22 (randCoefs 48)),
    "48 random coefficient vectors at the points of an ellipsoid cloud")
  fs := fs.push ({ cloudCase "pat-error-fit-nonunit" "fit" #[(⟨0.0, 0.0, 0.0⟩ : V3)] #[(⟨0.0, 0.0, 1.5⟩ : V3)] with coefs := #[(0, sphereCoef)] }, "normal of length 1.5: NonUnitNormal")
  -- ---- sdf
  let q1 : Array V3 := #[⟨1.0, 0.0, 0.0⟩, ⟨0.0, 1.0, 0.0⟩, ⟨-0.6, 0.8, 0.0⟩, ⟨0.0, 0.0, 0.5⟩, ⟨0.0, 0.0, -2.0⟩, ⟨0.0, 0.0, 0.0⟩, ⟨1.25, 0.0, 0.0⟩, ⟨0.75, 0.0, 0.0⟩, ⟨1.0, 0.0, 0.25⟩, ⟨2.0, 0.0, 0.0⟩, ⟨0.3, 0.4, 0.1⟩, ⟨-3.0, 1.0, 2.0⟩]
  let sdfCase := fun (name : String) (ts : Array Torus) (qs : Array V3) => { cloudCase name "sdf" #[] #[] with tori := ts, queries := qs }
  fs := fs.push (sdfCase "pat-sdf-torus" #[{ c := V3.zero, u := ez3, R := 1.0, r := 0.25 }] q1,
    "z-axis torus R = 1, r = 0.25: queries on the core circle (1,0,0), (0,1,0), (-0.6,0.8,0) [kind 2], on the axis (0,0,0.5), (0,0,-2), (0,0,0) [kind 1], inside and outside the tube")
  fs := fs.push (sdfCase "pat-sdf-inverted" #[{ c := ⟨0.5, -0.5, 0.25⟩, u := ez3, R := 1.0, r := -0.25 }, { c := V3.zero, u := normalize ⟨1.0, 1.0, 0.0⟩, R := 0.7, r := -0.3 }] (q1 ++ #[⟨1.5, -0.5, 0.25⟩, ⟨0.5, -0.5, 0.75⟩]),
    "negative minor radius (sign -1): the SDF is negated; a tilted second torus; core circle at (1.5,-0.5,0.25) and axis point (0.5,-0.5,0.75)")
  fs := fs.push (sdfCase "pat-sdf-spindle-sphere" #[{ c := V3.zero, u := ez3, R := 0.2, r := 0.5 }, { c := ⟨0.1, 0.2, 0.3⟩, u := ey3, R := 0.0, r := 1.0 }, { c := ⟨0.1, 0.2, 0.3⟩, u := ex3, R := 0.0, r := -2.0 }]
      (q1 ++ #[⟨0.1, 0.2, 0.3⟩, ⟨0.1, 1.2, 0.3⟩, ⟨0.1, 0.2, 1.3⟩, ⟨0.4, 0.6, 0.3⟩]),
    "R < |r| (spindle) and R = 0 (sphere of radius |r|; the center is the only core point, kind 2; points on the axis kind 1)")
  let randTori := run 31 (do
    let mut ts := #[]
    for _ in [0:6] do
      let c ← vecIn (-1.0) 1.0
      let u ← unitVec
      let R ← uni 0.0 2.0
      let r0 ← uni 0.1 1.2
      let sg ← below 2
      ts := ts.push ({ c, u, R, r := if sg == 0 then r0 else -r0 } : Torus)
    return ts)
  fs := fs.push (sdfCase "pat-sdf-random" randTori (run 32 (randQueries 40 (-2.5) 2.5)), "six random tori with skewed axes, 40 random queries each")
  -- ---- blend on the unit sphere with tangent-sphere coefficients: the SDF is |x| - 1
  let (sp, sn) := fibonacci 50 1.0 V3.zero
  let sphereQ := spreadQueries #[sp[0]!, sp[7]!, V3.zero, ⟨1.0, 0.0, 0.0⟩, ⟨0.0, 0.0, 3.0⟩, sp[10]! + ⟨0.01, 0.0, 0.0⟩] 200 (-2.0) 2.0 41
  let sphereBase := withCoefs (cloudCase "x" "blend" sp sn) ((Array.range 50).map fun _ => sphereCoef)
  fs := fs.push ({ sphereBase with name := "pat-blend-sphere-default", queries := sphereQ },
    "50-point unit sphere, tangent-sphere coefficients, defaults (knn 32, PerQuery C = 64); queries include a data point (singR), the center (all tori on their core, singG) and near-data-point queries")
  fs := fs.push ({ sphereBase with name := "pat-blend-sphere-all", queries := sphereQ, opts := { nbhd := .all } },
    "same cloud, all points (Eq. 26 exactly)")
  fs := fs.push ({ sphereBase with name := "pat-blend-sphere-constant", queries := sphereQ, opts := { nbhd := .knn 16, lam := .constant 20.0 } },
    "constant lambda = 20, knn 16 (the exact gradient equals the frozen one)")
  let sphereAreas := run 42 (randAreas 50)
  let sphereAreaBase := { sphereBase with areas := sphereAreas }
  fs := fs.push ({ sphereAreaBase with name := "pat-blend-sphere-paper", queries := sphereQ, opts := { nbhd := .radius none 32, lam := .paper 64.0 1000.0 20, areas := true } },
    "paper heuristic: lambda = 1e3/D, R_eval = 2C/lambda (D is printed in the paper line), kmean = 20, area weights; empty radius sets fall back to 32-NN")
  fs := fs.push ({ sphereBase with name := "pat-blend-sphere-radius", queries := sphereQ, opts := { nbhd := .radius (some 0.9) 5 } },
    "fixed radius 0.9, fallback 5, PerQuery")
  let (ep, en) := run 51 (ellipsoidCloud 50)
  let ec := run 52 (randCoefs 50)
  let ellBase := withCoefs (cloudCase "x" "blend" ep en) ec
  let ellQ := spreadQueries #[ep[3]!, ep[20]!, ep[3]! + ⟨0.0, 0.02, 0.0⟩, ⟨0.0, 0.0, 0.0⟩] 200 (-1.5) 1.5 53
  fs := fs.push ({ ellBase with name := "pat-blend-random-default", queries := ellQ },
    "50 points of an ellipsoid cloud with random coefficients, defaults; queries include data points and the origin")
  fs := fs.push ({ ellBase with name := "pat-blend-random-knn8-iso", queries := ellQ, opts := { nbhd := .knn 8, lam := .perQuery 16.0, iso := 0.1 } },
    "knn 8, C = 16, isovalue 0.1 (subtracted at the end)")
  fs := fs.push ({ ellBase with name := "pat-blend-random-areas", queries := ellQ, areas := sphereAreas, opts := { areas := true, nbhd := .knn 20, lam := .constant 8.0 } },
    "area weights, knn 20, constant lambda 8")
  -- explicit tori and singular sets
  let t3 : Array Torus := #[{ c := V3.zero, u := ez3, R := 1.0, r := 0.25 }, { c := ⟨0.5, 0.0, 0.0⟩, u := ez3, R := 0.5, r := -0.2 }, { c := ⟨0.0, 0.5, 0.0⟩, u := ex3, R := 0.4, r := 0.3 }]
  fs := fs.push ({ cloudCase "pat-blend-singular" "blend" #[(⟨1.0, 0.0, 0.25⟩ : V3), ⟨1.0, 0.0, 0.0⟩, ⟨0.5, 0.0, 0.5⟩] #[ex3, ex3, ez3] with tori := t3, opts := { nbhd := .all, lam := .constant 2.0 }, queries := #[⟨1.0, 0.0, 0.0⟩, ⟨0.0, 0.0, 0.5⟩, ⟨0.5, 0.5, 0.0⟩, ⟨0.0, 0.5, 0.4⟩, ⟨1.0, 0.0, 0.25⟩, ⟨0.3, 0.1, 0.2⟩] },
    "three explicit tori, all points: queries on a core circle (singG), on an axis, at a data point (singR)")
  fs := fs.push ({ cloudCase "pat-blend-duplicate" "blend" #[(⟨0.0, 0.0, 0.0⟩ : V3), ⟨0.0, 0.0, 0.0⟩, ⟨1.0, 0.0, 0.0⟩, ⟨0.0, 1.0, 0.0⟩] #[ez3, ez3, ez3, ez3] with tori := #[{ c := ⟨0.0, 0.0, -1.0⟩, u := ez3, R := 0.0, r := 1.0 }, { c := ⟨0.0, 0.0, -2.0⟩, u := ez3, R := 0.0, r := 2.0 }, { c := ⟨1.0, 0.0, -1.0⟩, u := ez3, R := 0.0, r := 1.0 }, { c := ⟨0.0, 1.0, -1.0⟩, u := ez3, R := 0.0, r := 1.0 }], opts := { nbhd := .knn 3, lam := .perQuery 8.0 }, queries := #[⟨0.0, 0.0, 0.0⟩, ⟨0.1, 0.1, 0.1⟩, ⟨0.5, 0.5, 0.2⟩] },
    "two coincident data points (ties by index), knn 3: the order of the sums")
  fs := fs.push ({ cloudCase "pat-error-blend-single" "blend" #[(⟨0.0, 0.0, 0.0⟩ : V3)] #[ez3] with tori := #[{ c := ⟨0.0, 0.0, -1.0⟩, u := ez3, R := 0.0, r := 1.0 }], queries := #[V3.zero] },
    "one point, query at the point: sigma = 0, DegenerateQuery")
  fs := fs.push ({ cloudCase "pat-error-blend-count" "blend" #[(⟨0.0, 0.0, 0.0⟩ : V3), ⟨1.0, 0.0, 0.0⟩] #[ez3, ez3] with coefs := #[(0, sphereCoef)], queries := #[V3.zero] },
    "one coefficient line for two points: ValueCount")
  -- ---- baselines
  let areaSphere := sp.map fun _ => 4.0 * 3.141592653589793 / 50.0
  let baseBase := { cloudCase "x" "baselines" sp sn with areas := areaSphere }
  fs := fs.push ({ baseBase with name := "pat-baselines-sphere-default", queries := sphereQ },
    "SHC and SSPD on the unit sphere with areas 4 pi / 50, defaults; queries include a data point (skipped SHC term)")
  fs := fs.push ({ baseBase with name := "pat-baselines-sphere-all", queries := sphereQ, opts := { nbhd := .all, lam := .constant 6.0 } },
    "all points, constant lambda 6")
  fs := fs.push ({ baseBase with name := "pat-baselines-random-paper", areas := sphereAreas, pts := ep, nrm := en, queries := ellQ.extract 2 ellQ.size, opts := { nbhd := .radius none 20, lam := .paper 64.0 1000.0 20 } },
    "ellipsoid cloud, random areas, paper rule with radius auto")
  return fs
where
  ez3 : V3 := ⟨0.0, 0.0, 1.0⟩
  ey3 : V3 := ⟨0.0, 1.0, 0.0⟩
  ex3 : V3 := ⟨1.0, 0.0, 0.0⟩

/-- Does the printed input parse back to a case with the same output? -/
def roundTrips (c : PatCase) : Bool :=
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
def compareOut (expected actual : String) (tol : Float := 1e-9) : Option String := Id.run do
  let le := (expected.splitOn "\n").toArray
  let la := (actual.splitOn "\n").toArray
  if le.size != la.size then return some s!"{le.size} lines expected, {la.size} produced"
  for i in [0:le.size] do
    let te := tokens le[i]!
    let ta := tokens la[i]!
    if te.size != ta.size then return some s!"line {i + 1}: '{la[i]!}' instead of '{le[i]!}'"
    for j in [0:te.size] do
      if te[j]!.length == 16 && ta[j]!.length == 16 then
        match hexToFloat? te[j]!, hexToFloat? ta[j]! with
        | some x, some y =>
          let scale := if x.abs > 1.0 then x.abs else 1.0
          if !(x == y || (x.isNaN && y.isNaN) || (x - y).abs ≤ tol * scale) then
            return some s!"line {i + 1}, token {j + 1}: differs ('{le[i]!}')"
        | _, _ => if te[j]! != ta[j]! then return some s!"line {i + 1}: '{la[i]!}' instead of '{le[i]!}'"
      else if te[j]! != ta[j]! then return some s!"line {i + 1}: '{la[i]!}' instead of '{le[i]!}'"
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
    let detail := if st == "error" then ((output.splitOn "\n").find? (·.startsWith "status ")).getD "" else ""
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
        else if let some d := compareOut stored text then
          problems := problems.push s!"{name}{suffix}: {d}"
    if problems.size == 0 then IO.println s!"ok    {name}"
    else
      bad := bad + 1
      for p in problems do IO.println s!"FAIL  {p}"
  for (c, _) in fixtures do
    if !roundTrips c then
      bad := bad + 1
      IO.println s!"FAIL  {c.name}: the printed input does not parse back to the same case"
  IO.println s!"pat fixtures: {fixtures.size} cases, {bad} problems"
  return bad

end Ddg.PointsAsTori
