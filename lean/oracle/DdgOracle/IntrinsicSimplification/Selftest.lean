import DdgOracle.IntrinsicSimplification.Cases

/-!
# Self-test of the `iem` oracle

Float checks, independent of the Lean proofs (which prove the general statements):

* the regular cones flatten to the closed form `u = 2 log(s/(2 r sin(π/n)))`;
* the Newton denominator is minus the derivative of `Θ_i` (central differences), also with a
  self-edge, where the per-corner sum printed in Eq. 4 is not;
* `Θ_i(u)` is non-increasing on a grid of `u` (extended angles);
* flattening conserves curvature (`Σ_j ΔK_j = K_i`), the weights are non-negative and sum to one,
  the cost without memory equals `Σ α_j m_i ℓ̃_ij`;
* Eq. 10's projective map: the weights `e^{−u}` (Springborn et al. 2008, App. D) send the
  circumcircle of the face to the circumcircle of the flattened face, the printed `e^{+u}` does not;
* every fixture case gives a block without error, and reading back its printed input gives the
  same output.
-/

namespace Ddg.IEM

open Ddg

private def fail (msg : String) : IO Nat := do
  IO.println s!"FAIL  {msg}"
  return 1

/-- Plane layout of a triangle with sides `l0` (opposite corner 0) … : corner 0 at the origin,
corner 1 on the x axis. -/
def layout (l0 l1 l2 : Float) : Array (Float × Float) :=
  let x := (l2 * l2 + l1 * l1 - l0 * l0) / (2.0 * l2)
  #[(0.0, 0.0), (l2, 0.0), (x, Float.sqrt (l1 * l1 - x * x))]

/-- Circumcentre and radius. -/
def circum (p : Array (Float × Float)) : (Float × Float) × Float :=
  let (ax, ay) := p[0]!
  let (bx, by_) := p[1]!
  let (cx, cy) := p[2]!
  let d := 2.0 * (ax * (by_ - cy) + bx * (cy - ay) + cx * (ay - by_))
  let a2 := ax * ax + ay * ay
  let b2 := bx * bx + by_ * by_
  let c2 := cx * cx + cy * cy
  let ux := (a2 * (by_ - cy) + b2 * (cy - ay) + c2 * (ay - by_)) / d
  let uy := (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
  ((ux, uy), Float.sqrt ((ax - ux) * (ax - ux) + (ay - uy) * (ay - uy)))

/-- Barycentric coordinates of `q` in `p`. -/
def bary (p : Array (Float × Float)) (q : Float × Float) : Float × Float × Float :=
  let s := fun (a b c : Float × Float) => (b.1 - a.1) * (c.2 - a.2) - (b.2 - a.2) * (c.1 - a.1)
  let t := s p[0]! p[1]! p[2]!
  (s q p[1]! p[2]! / t, s p[0]! q p[2]! / t, s p[0]! p[1]! q / t)

/-- Largest deviation from the new circumcircle of the images of points of the old one, under the
projective map with the weight `w0` at corner 0 (`u` applied at corner 0). -/
def circleError (l0 l1 l2 u w0 : Float) : Float := Id.run do
  let p := layout l0 l1 l2
  let e := Float.exp (u / 2.0)
  let q := layout l0 (l1 * e) (l2 * e)
  let ((cx, cy), r) := circum p
  let ((dx, dy), r') := circum q
  let mut worst := 0.0
  for k in [0:12] do
    let t := 2.0 * floatPi * Float.ofNat k / 12.0
    let pt := (cx + r * Float.cos t, cy + r * Float.sin t)
    let (b0, b1, b2) := bary p pt
    let (c0, c1, c2) := iemProjInterp w0 1.0 1.0 b0 b1 b2
    let x := c0 * q[0]!.1 + c1 * q[1]!.1 + c2 * q[2]!.1
    let y := c0 * q[0]!.2 + c1 * q[1]!.2 + c2 * q[2]!.2
    let err := Float.abs (Float.sqrt ((x - dx) * (x - dx) + (y - dy) * (y - dy)) - r')
    if err > worst then worst := err
  return worst

/-- Run the self-test; returns the number of failures. -/
def selftest : IO Nat := do
  let mut bad := 0
  -- cones
  for (n, s) in [(3, 1.2), (4, 1.0), (5, 1.0), (6, 0.9), (8, 1.0), (6, 1.0), (7, 0.7)] do
    let r := flatten (coneStar n 1.0 s) (2.0 * floatPi)
    let exact := 2.0 * Float.log (s / (2.0 * Float.sin (floatPi / Float.ofNat n)))
    if r.status != .ok || Float.abs (r.u - exact) > 1e-12 then
      bad := bad + (← fail s!"cone {n} {s}: u = {r.u}, closed form {exact}, status {r.status.toStr}")
  -- derivative by central differences
  let stars : Array (String × Array SFace) := #[
    ("interior", randomStar 3 6 true (2.3 * floatPi)), ("boundary", randomStar 104 4 false (0.8 * floatPi)),
    ("cone", coneStar 5 1.0 0.8), ("self-edge", selfEdgeStar)]
  let mut paperGap := 0.0
  for (nm, fs) in stars do
    for u in [-0.3, -0.1, 0.0, 0.1, 0.2] do
      if allValid fs (u - 1e-4) && allValid fs (u + 1e-4) then
        let hh := 1e-6
        let fd := (thetaI fs (u + hh) - thetaI fs (u - hh)) / (2.0 * hh)
        let H := hessian fs u
        if Float.abs (fd + H) > 1e-6 * (1.0 + Float.abs H) then
          bad := bad + (← fail s!"{nm} u = {u}: dTheta/du = {fd}, -hessian = {-H}")
        let gap := Float.abs (fd + paperHessian fs u)
        if nm == "self-edge" then
          if gap > paperGap then paperGap := gap
        else if gap > 1e-6 * (1.0 + Float.abs H) then
          bad := bad + (← fail s!"{nm} u = {u}: the paper's denominator differs on a simplicial star")
  if paperGap < 1e-3 then
    bad := bad + (← fail s!"self-edge star: the per-corner denominator of Eq. 4 should differ from the derivative (gap {paperGap})")
  else IO.println s!"ok    self-edge star: |dTheta/du + paper denominator| up to {paperGap} (the corrected denominator matches to 1e-6)"
  -- monotone on a grid
  for (nm, fs) in stars do
    let mut prev := thetaI fs (-8.0)
    for k in [1:321] do
      let u := -8.0 + 0.05 * Float.ofNat k
      let th := thetaI fs u
      if th > prev + 1e-12 then
        bad := bad + (← fail s!"{nm}: Theta increases at u = {u} ({prev} -> {th})")
      prev := th
  -- conservation, weights, memoryless cost on the removal fixtures
  for (c, _) in allCases do
    if c.task == "removal" then
      match parseBody c with
      | .error e => bad := bad + (← fail s!"{c.name}: {e}")
      | .ok b =>
        let r := removal b.faces b.target b.nbrs b.chans b.memory
        if r.flat.status == .ok then
          let sdk := r.nbr.foldl (fun s x => s + x.1) 0.0
          if Float.abs (sdk - r.Ki) > 1e-10 then
            bad := bad + (← fail s!"{c.name}: sum dK = {sdk}, K_i = {r.Ki}")
          let sa := r.nbr.foldl (fun s x => s + x.2.1) 0.0
          if Float.abs (sa - 1.0) > 1e-12 || r.nbr.any (fun x => x.2.1 < 0.0) then
            bad := bad + (← fail s!"{c.name}: weights sum to {sa}")
          if !b.memory && Float.abs (r.cost - r.memoryless) > 1e-12 * (1.0 + r.cost) then
            bad := bad + (← fail s!"{c.name}: cost {r.cost}, memoryless form {r.memoryless}")
  -- projective interpolation
  let springborn := circleError 1.0 0.8 0.9 0.7 (Float.exp (-0.7))
  let printed := circleError 1.0 0.8 0.9 0.7 (Float.exp 0.7)
  if springborn > 1e-12 || printed < 1e-3 then
    bad := bad + (← fail s!"circumcircle: weights e^-u err {springborn}, e^+u err {printed}")
  else IO.println s!"ok    Eq. 10: weights e^-u keep the circumcircle ({springborn}), the printed e^+u do not ({printed})"
  -- fixtures run and round-trip
  for (c, _) in allCases do
    let out := runCase c
    if (out.splitOn "status error").length > 1 then
      bad := bad + (← fail s!"{c.name}: error block")
    let text := printCase c
    let lines := (text.splitOn "\n").filter (· ≠ "") |>.map toks |>.toArray
    match parseCase c.name (lines.extract 1 (lines.size - 1)) with
    | .error e => bad := bad + (← fail s!"{c.name}: does not parse back: {e}")
    | .ok c' => if runCase c' != out then bad := bad + (← fail s!"{c.name}: does not round-trip")
  IO.println s!"iem selftest: {bad} failure(s), {allCases.size} cases"
  return bad

end Ddg.IEM
