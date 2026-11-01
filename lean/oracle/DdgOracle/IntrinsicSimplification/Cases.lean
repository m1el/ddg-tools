import DdgOracle.IntrinsicSimplification.Format

/-!
# Cases of the `iem` oracle

Stars built from spokes and angles at the removed vertex (the outer sides by the law of cosines, so
every face is a strict triangle before flattening), regular cones with a closed-form flattening, a
Δ-complex star with a self-edge, failure cases, and seeded random stars (SplitMix64, `Ddg.Rng`).
-/

namespace Ddg.IEM

open Ddg

/-- Tokens of a line. -/
def toks (s : String) : Array String := ((s.splitOn " ").filter (· ≠ "")).toArray

private def h (x : Float) : String := floatToHex x

/-- A `face` line. -/
def faceLine (v0 v1 v2 : Nat) (l0 l1 l2 : Float) : Array String :=
  toks s!"face {v0} {v1} {v2} {h l0} {h l1} {h l2}"

/-- A star with spokes `r` and angles `θ` at the removed vertex: faces `(0, j_k, j_{k+1})` with
`j_k = k + 1`; closed (`interior`, the last face returns to `j_0`) or open (boundary, `r` has one
more entry than `θ`). Side 0 is the outer side `√(r_k² + r_{k+1}² − 2 r_k r_{k+1} cos θ_k)`. -/
def fanFaces (r θ : Array Float) (interior : Bool) : Array SFace := Id.run do
  let n := θ.size
  let mut out : Array SFace := #[]
  for k in [0:n] do
    let k1 := if interior then (k + 1) % n else k + 1
    let a := r[k]!
    let b := r[k1]!
    let s := Float.sqrt (a * a + b * b - 2.0 * a * b * Float.cos θ[k]!)
    -- corner 0 = i, corner 1 = j_k, corner 2 = j_{k+1}; side 1 joins corners 2, 0 (spoke b),
    -- side 2 joins corners 0, 1 (spoke a)
    out := out.push { v := #[0, k + 1, k1 + 1], l := #[s, b, a] }
  return out

/-- The face lines of a star. -/
def faceLines (fs : Array SFace) : Array (Array String) :=
  fs.map fun f => faceLine f.v[0]! f.v[1]! f.v[2]! f.l[0]! f.l[1]! f.l[2]!

/-- The target line. -/
def targetLine (interior : Bool) : Array String :=
  toks (if interior then "target interior" else "target boundary")

/-- Random positive angles summing to `total`, each below `0.9 π`. -/
def randomAngles (rng : Rng) (n : Nat) (total : Float) : Array Float × Rng := Id.run do
  let mut r := rng
  let mut w : Array Float := #[]
  for _ in [0:n] do
    let (x, r') := r.range 0.5 1.5
    w := w.push x
    r := r'
  let s := w.foldl (· + ·) 0.0
  let a := w.map fun x => total * x / s
  let a := a.map fun x => if x < 0.9 * floatPi then x else 0.9 * floatPi
  return (a, r)

/-- Random spokes in `[0.6, 1.4]`. -/
def randomSpokes (rng : Rng) (n : Nat) : Array Float × Rng := Id.run do
  let mut r := rng
  let mut out : Array Float := #[]
  for _ in [0:n] do
    let (x, r') := r.range 0.6 1.4
    out := out.push x
    r := r'
  return (out, r)

/-- A random star: degree `deg`, total angle `total` at the vertex. -/
def randomStar (seed : Nat) (deg : Nat) (interior : Bool) (total : Float) : Array SFace :=
  let rng : Rng := ⟨seed.toUInt64 * 0x9E3779B97F4A7C15 + 17⟩
  let (θ, rng) := randomAngles rng deg total
  let (r, _) := randomSpokes rng (if interior then deg else deg + 1)
  fanFaces r θ interior

/-- The regular cone of degree `n`: spokes `r`, outer sides `s`. Its flattening is
`u = 2 log(s / (2 r sin(π/n)))` (all angles at the apex `2π/n`). -/
def coneStar (n : Nat) (r s : Float) : Array SFace :=
  (List.range n).toArray.map fun k =>
    { v := #[0, k + 1, (k + 1) % n + 1], l := #[s, r, r] }

/-- A Δ-complex star: a face with two corners at `i` (a self-edge `ii`, side 2) and two ordinary
faces. -/
def selfEdgeStar : Array SFace :=
  #[{ v := #[0, 0, 1], l := #[1.1, 0.9, 1.3] },
    { v := #[0, 1, 2], l := #[0.8, 1.0, 1.2] },
    { v := #[0, 2, 1], l := #[1.05, 0.95, 1.15] }]

/-- `removal` lines for a star: neighbours with the angular coordinates of the flattened star (the
accumulated angles at `i`, scaled to `2π` or `π`) and random frames at the neighbours; two
curvature channels and an area-like channel with random masses and error vectors. -/
def removalLines (fs : Array SFace) (interior : Bool) (seed : Nat) (memory : Bool) :
    Array (Array String) := Id.run do
  let target := if interior then 2.0 * floatPi else floatPi
  let fl := flatten fs target
  let mut lines := #[targetLine interior] ++ faceLines fs
  -- neighbours 1..n in label order; φ_ij from the flattened angles at `i`
  let n := fs.foldl (fun m f => max m (max f.v[1]! f.v[2]!)) 0
  let mut rng : Rng := ⟨seed.toUInt64 + 99⟩
  let mut acc := 0.0
  let mut phis : Array Float := Array.replicate (n + 1) 0.0
  let mut lens : Array Float := Array.replicate (n + 1) 1.0
  for f in fs do
    let L := f.sides fl.u
    -- the spoke to corner 1 is side 2, the spoke to corner 2 is side 1
    if f.v[0]! == 0 && f.v[1]! != 0 then
      phis := phis.set! f.v[1]! acc
      lens := lens.set! f.v[1]! f.l[2]!
    acc := acc + angleAt L 0
    if f.v[0]! == 0 && f.v[2]! != 0 then
      lens := lens.set! f.v[2]! f.l[1]!
  if !interior then phis := phis.set! (n) floatPi
  for j in [1:n + 1] do
    let (pj, r') := rng.range 0.0 (2.0 * floatPi)
    rng := r'
    lines := lines.push (toks s!"nbr {j} {h phis[j]!} {h pj} {h lens[j]!}")
  let Ki := target - thetaI fs 0.0
  let chans : Array (Float × Float) :=
    #[(1.0, if Ki > 0.0 then Ki else 0.0), (1.0, if Ki < 0.0 then -Ki else 0.0), (0.25, 0.1)]
  for (w, mi) in chans do
    let (tx, r1) := rng.range (-0.2) 0.2
    let (ty, r2) := r1.range (-0.2) 0.2
    rng := r2
    lines := lines.push (toks s!"channel {h w} {h mi} {h tx} {h ty}")
    for j in [1:n + 1] do
      let (m, r1) := rng.range 0.0 0.3
      let (x, r2) := r1.range (-0.2) 0.2
      let (y, r3) := r2.range (-0.2) 0.2
      rng := r3
      lines := lines.push (toks s!"mass {j} {h m} {h x} {h y}")
  if !memory then lines := lines.push (toks "memory 0")
  return lines

/-- A case. -/
def mkCase (name task : String) (lines : Array (Array String)) : ICase := { name, task, lines }

/-- All cases of the fixtures, with what they are. -/
def allCases : Array (ICase × String) := Id.run do
  let mut out : Array (ICase × String) := #[]
  -- regular cones: closed form u = 2 log(s/(2 r sin(π/n)))
  for (n, s) in [(3, 1.2), (4, 1.0), (5, 1.0), (6, 0.9), (8, 1.0), (6, 1.0)] do
    let fs := coneStar n 1.0 s
    let what := s!"regular cone, degree {n}, spokes 1, outer sides {s}: flattening u = 2 log(s/(2 sin(π/{n})))"
    out := out.push (mkCase s!"iem-flatten-cone-{n}-{s}" "flatten"
      (#[targetLine true] ++ faceLines fs), what)
  -- random interior and boundary stars
  for seed in [0:12] do
    let deg := 3 + seed % 6
    let total := 2.0 * floatPi * (0.8 + 0.4 * Float.ofNat (seed % 5) / 4.0)
    let fs := randomStar seed deg true total
    out := out.push (mkCase s!"iem-flatten-interior-{seed}" "flatten"
      (#[targetLine true] ++ faceLines fs),
      s!"random interior star, degree {deg}, angle sum {total}")
  for seed in [0:8] do
    let deg := 2 + seed % 4
    let total := floatPi * (0.6 + 0.8 * Float.ofNat (seed % 5) / 4.0)
    let fs := randomStar (100 + seed) deg false total
    out := out.push (mkCase s!"iem-flatten-boundary-{seed}" "flatten"
      (#[targetLine false] ++ faceLines fs),
      s!"random boundary star, degree {deg}, angle sum {total}")
  -- Δ-complex star with a self-edge
  out := out.push (mkCase "iem-flatten-self-edge" "flatten"
    (#[targetLine true] ++ faceLines selfEdgeStar),
    "a face with two corners at i (self-edge ii) and two ordinary faces")
  -- failures
  out := out.push (mkCase "iem-flatten-ear" "flatten"
    (#[targetLine false] ++ faceLines (randomStar 7 1 false 2.0)),
    "an ear (boundary vertex of degree 1): the root makes the face degenerate (status invalid)")
  out := out.push (mkCase "iem-flatten-degree-two" "flatten"
    (#[targetLine true] ++ faceLines (randomStar 8 2 true (1.6 * floatPi))),
    "an interior vertex of degree 2: the root makes both faces degenerate (status invalid)")
  out := out.push (mkCase "iem-flatten-self-face" "flatten"
    (#[targetLine true, faceLine 0 0 0 1.0 1.1 0.9]),
    "a self-face: its angle sum is pi for every u, target 2 pi (status noroot)")
  -- Newton steps
  for (u, nm) in [(0.0, "0"), (0.3, "0.3"), (-0.4, "-0.4"), (2.5, "2.5")] do
    out := out.push (mkCase s!"iem-step-interior-{nm}" "step"
      (#[targetLine true] ++ faceLines (randomStar 3 6 true (2.3 * floatPi)) ++ #[toks s!"u {h u}"]),
      s!"one Newton step of a random interior star (degree 6) at u = {nm}")
  out := out.push (mkCase "iem-step-self-edge" "step"
    (#[targetLine true] ++ faceLines selfEdgeStar ++ #[toks s!"u {h 0.1}"]),
    "Newton step of the self-edge star: hess (the derivative) differs from the per-corner sum of Eq. 4 (paperhess)")
  out := out.push (mkCase "iem-step-boundary" "step"
    (#[targetLine false] ++ faceLines (randomStar 104 4 false (0.8 * floatPi)) ++
      #[toks s!"u {h (-0.2)}"]),
    "Newton step of a random boundary star")
  -- weights
  out := out.push (mkCase "iem-weights-basic" "weights"
    #[toks s!"ki {h 0.3}", toks s!"dk 1 {h 0.1}", toks s!"dk 2 {h (-0.05)}", toks s!"dk 3 {h 0.25}"],
    "Eq. 9: |dK_j| / sum |dK_l|")
  out := out.push (mkCase "iem-weights-flat" "weights"
    #[toks s!"ki {h 1e-13}", toks s!"dk 1 {h 1e-14}", toks s!"dk 2 {h 0.0}", toks s!"dk 4 {h 0.0}"],
    "a flat vertex (|K_i| <= 1e-12): equal shares")
  -- removals
  for seed in [0:6] do
    let deg := 3 + seed % 4
    let total := 2.0 * floatPi * (0.85 + 0.1 * Float.ofNat seed)
    let fs := randomStar (200 + seed) deg true total
    out := out.push (mkCase s!"iem-removal-interior-{seed}" "removal"
      (removalLines fs true seed true),
      s!"one removal of a random interior star (degree {deg}), three channels, memory")
  for seed in [0:3] do
    let fs := randomStar (300 + seed) (2 + seed) false (floatPi * (0.7 + 0.25 * Float.ofNat seed))
    out := out.push (mkCase s!"iem-removal-boundary-{seed}" "removal"
      (removalLines fs false (10 + seed) true),
      s!"one removal of a random boundary star (degree {2 + seed}), three channels, memory")
  out := out.push (mkCase "iem-removal-memoryless" "removal"
    (removalLines (randomStar 201 4 true (1.7 * floatPi)) true 1 false),
    "memory 0: cost equals the memoryless form sum alpha m_i l_ij")
  out := out.push (mkCase "iem-removal-flat" "removal"
    (removalLines (coneStar 6 1.0 1.0) true 2 true),
    "an already flat vertex (regular hexagon): equal weights 1/6")
  out := out.push (mkCase "iem-removal-ear" "removal"
    (removalLines (randomStar 7 1 false 2.0) false 3 true),
    "an ear: the flattening is invalid, cost +inf")
  return out

end Ddg.IEM
