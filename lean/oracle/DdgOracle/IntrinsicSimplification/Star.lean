import DdgOracle.Formulas
import DdgOracle.Hex
import DdgOracle.Generators

/-!
# Intrinsic simplification: one star, in `Float` (`crates/ddg-intrinsic-simplification/SPEC.md`)

"Surface Simplification using Intrinsic Error Metrics" (Liu, Gillespie, Chislett, Sharp, Jacobson,
Crane 2023), written from the paper and SPEC §3.4, §4, not from the authors' code. Every formula is
the scalar-generic one of `DdgOracle/Formulas/IntrinsicSimplification.lean` (proved about at `ℝ`
in `lean/proofs/DdgProofs/IntrinsicSimplification/`); this file evaluates them on a star, solves
the flattening and assembles one removal.

**The star.** Faces with three corner labels (`0` is the vertex `i` being removed, `≥ 1` its
neighbours; labels may repeat: a face with two corners `0` has a self-edge `ii`) and three sides
(side `k` opposite corner `k`). Flattening with scale `u` multiplies side `k` by `e^{u m/2}`, `m` the
number of its ends labelled `0` (Eq. 3 with `u_j = 0`).
-/

namespace Ddg.IEM

open Ddg

/-- A face of the star. -/
structure SFace where
  /-- corner labels (`0` = the vertex being removed) -/
  v : Array Nat
  /-- sides, side `k` opposite corner `k` -/
  l : Array Float
  deriving Inhabited

/-- The ends of side `k` that are the vertex `0`. -/
def SFace.mult (f : SFace) (k : Nat) : Float :=
  (if f.v[(k + 1) % 3]! == 0 then 1.0 else 0.0) + (if f.v[(k + 2) % 3]! == 0 then 1.0 else 0.0)

/-- The sides after flattening with scale `u` (Eq. 3). -/
def SFace.sides (f : SFace) (u : Float) : Array Float :=
  #[0, 1, 2].map fun k => f.l[k]! * Float.exp (u * f.mult k / 2.0)

/-- Strict triangle inequalities. -/
def isValid (L : Array Float) : Bool :=
  L[0]! < L[1]! + L[2]! && L[1]! < L[2]! + L[0]! && L[2]! < L[0]! + L[1]!

/-- The extended angle at corner `k` (`iemExtAngle`, CETM Sec. 3.1). -/
def angleAt (L : Array Float) (k : Nat) : Float :=
  iemExtAngle L[k]! L[(k + 1) % 3]! L[(k + 2) % 3]!

/-- The cotangent at corner `k`; `0` for a face that is not a strict triangle (SPEC §3.4). -/
def cotAt (L : Array Float) (k : Nat) : Float :=
  if isValid L then cornerCot L[k]! L[(k + 1) % 3]! L[(k + 2) % 3]! else 0.0

/-- `Θ` of the label `v` over the star faces at scale `u`. -/
def angleSumOf (faces : Array SFace) (v : Nat) (u : Float) : Float := Id.run do
  let mut s := 0.0
  for f in faces do
    let L := f.sides u
    for k in [0:3] do
      if f.v[k]! == v then s := s + angleAt L k
  return s

/-- `Θ_i(u)`. -/
def thetaI (faces : Array SFace) (u : Float) : Float := angleSumOf faces 0 u

/-- The denominator of Eq. 4, `½ Σ_f hessTerm f`: per corner at `i`, the cotangents opposite its
sides that are not self-edges (`iem-flatten-newton`; the theorem `hasDerivAt_angleSum`). -/
def hessian (faces : Array SFace) (u : Float) : Float := Id.run do
  let mut s := 0.0
  for f in faces do
    let L := f.sides u
    for k in [0:3] do
      if f.v[k]! == 0 then
        let k1 := (k + 1) % 3
        let k2 := (k + 2) % 3
        if f.v[k1]! != 0 then s := s + cotAt L k2
        if f.v[k2]! != 0 then s := s + cotAt L k1
  return 0.5 * s

/-- The denominator as printed in Eq. 4 (per corner at `i`, both other cotangents), for comparison:
equal to `hessian` unless the star has a self-edge at `i`. -/
def paperHessian (faces : Array SFace) (u : Float) : Float := Id.run do
  let mut s := 0.0
  for f in faces do
    let L := f.sides u
    for k in [0:3] do
      if f.v[k]! == 0 then s := s + cotAt L ((k + 1) % 3) + cotAt L ((k + 2) % 3)
  return 0.5 * s

/-- Outcome of the flattening. -/
inductive FlatStatus
  /-- converged, every face a strict triangle at the root -/
  | ok
  /-- converged, but some face is not a strict triangle at the root (Springborn's flip fallback is
  needed, SPEC §3.4) -/
  | invalid
  /-- no bracket: the target is not in the range of `Θ_i` on the scales tried -/
  | noroot
  /-- the iteration cap was reached -/
  | maxiter
  deriving BEq, Inhabited

/-- The status as a token. -/
def FlatStatus.toStr : FlatStatus → String
  | .ok => "ok"
  | .invalid => "invalid"
  | .noroot => "noroot"
  | .maxiter => "maxiter"

/-- Result of the flattening. -/
structure FlatResult where
  /-- outcome -/
  status : FlatStatus
  /-- Newton/bisection iterations -/
  iters : Nat
  /-- the scale -/
  u : Float
  /-- `Θ_i(u)` -/
  theta : Float
  deriving Inhabited

/-- Every face valid at `u`. -/
def allValid (faces : Array SFace) (u : Float) : Bool := faces.all fun f => isValid (f.sides u)

/-- **The flattening (SPEC §3.4, S9)**: a bracket `Θ(lo) ≥ Θ̃ ≥ Θ(hi)` grown from `u = 0` by doubling
steps (`|u| ≤ 256`), then Eq. 4's Newton step when it stays strictly inside the bracket and the
denominator is positive, bisection otherwise; stop at `|Θ − Θ̃| ≤ 1e-12·Θ̃`, at most 100 steps.
`iem-flatten-monotone` (the theorems `exists_root_of_bracket`, `angleSum_antitone`) is what makes the
bracket contain a root. -/
def flatten (faces : Array SFace) (target : Float) : FlatResult := Id.run do
  let tol := 1e-12 * target
  let th0 := thetaI faces 0.0
  let mut lo := 0.0
  let mut hi := 0.0
  let mut found := true
  if th0 ≥ target then
    -- grow `hi` until Θ(hi) ≤ target
    let mut step := 1.0
    found := false
    for _ in [0:9] do
      if !found then
        hi := step
        if thetaI faces hi ≤ target then found := true
        else lo := hi; step := 2.0 * step
  else
    let mut step := 1.0
    found := false
    for _ in [0:9] do
      if !found then
        lo := -step
        if thetaI faces lo ≥ target then found := true
        else hi := lo; step := 2.0 * step
  if !found then return { status := .noroot, iters := 0, u := 0.0, theta := th0 }
  let mut u := if th0 ≥ target then lo else hi
  let mut th := thetaI faces u
  let mut iters := 0
  let mut done : Bool := Float.abs (th - target) ≤ tol
  while !done && iters < 100 do
    iters := iters + 1
    let h := hessian faces u
    let un := iemNewtonStep u target th h
    let uNext := if h > 0.0 && lo < un && un < hi then un else 0.5 * (lo + hi)
    u := uNext
    th := thetaI faces u
    if th ≥ target then lo := u else hi := u
    if Float.abs (th - target) ≤ tol || hi - lo ≤ 1e-15 * (1.0 + Float.abs u) then done := true
  let st := if !done then FlatStatus.maxiter
    else if allValid faces u then FlatStatus.ok else FlatStatus.invalid
  return { status := st, iters, u, theta := th }

/-! ## One removal: weights, masses, error vectors, cost (Eq. 5, 8, 9, 6) -/

/-- A neighbour: its label, the angular coordinates `φ_ij` (at `i`) and `φ_ji` (at `j`) of the edge
used for the transport, after flattening (SPEC §4.3), and that edge's length before flattening. -/
structure Nbr where
  /-- label -/
  j : Nat
  /-- `φ_ij` -/
  phiIJ : Float
  /-- `φ_ji` -/
  phiJI : Float
  /-- `ℓ_ij` before flattening (the spoke; scaled by `e^{u/2}`) -/
  len : Float
  deriving Inhabited

/-- A channel: weight, mass and error vector at `i`, and per neighbour (in the order of the
neighbours) mass and error vector. -/
structure Channel where
  /-- weight `w^p` -/
  w : Float
  /-- `m_i` -/
  mi : Float
  /-- `t_i` -/
  ti : Float × Float
  /-- `(m_j, t_j)` per neighbour -/
  nb : Array (Float × (Float × Float))
  deriving Inhabited

/-- Result of one removal. -/
structure RemovalResult where
  /-- the flattening -/
  flat : FlatResult
  /-- `K_i = Θ̃_i − Θ_i` -/
  Ki : Float
  /-- per neighbour: `ΔK_j = K̃_j − K_j`, `α_j`, `ℓ̃_ij` -/
  nbr : Array (Float × Float × Float)
  /-- per channel, per neighbour: `m̃_j`, `t̃_j` -/
  chan : Array (Array (Float × (Float × Float)))
  /-- `C_i = Σ_j Σ_p w^p m̃_j ‖t̃_j‖` (Eq. 6) -/
  cost : Float
  /-- `Σ_p w^p Σ_j α_j m_i ℓ̃_ij` (Eq. 6, second form) -/
  memoryless : Float
  deriving Inhabited

/-- The flat-vertex threshold of S11. -/
def flatTol : Float := 1e-12

/-- **One removal** (SPEC §4.2–§4.4): flatten, the curvature changes of the neighbours
(`iem-flatten-conserves-curvature`), Eq. 9 (`iem-redistribution-weights`), Eq. 5 and Eq. 8 per channel
(`R_ij = e^{i(φ_ji + π − φ_ij)}`, `e_ji = ℓ̃_ij e^{iφ_ji}`), the cost of Eq. 6. With `memory = false`
the error vectors are taken as `0`. -/
def removal (faces : Array SFace) (target : Float) (nbrs : Array Nbr) (chans : Array Channel)
    (memory : Bool) : RemovalResult := Id.run do
  let fl := flatten faces target
  let Ki := target - thetaI faces 0.0
  let dK := nbrs.map fun n => -(angleSumOf faces n.j fl.u - angleSumOf faces n.j 0.0)
  let sumAbs := dK.foldl (fun s x => s + Float.abs x) 0.0
  let count := Float.ofNat nbrs.size
  let alpha := dK.map fun d => iemWeight d sumAbs Ki flatTol count
  let scale := Float.exp (fl.u / 2.0)
  let lens := nbrs.map fun n => n.len * scale
  let nbrOut := (List.range nbrs.size).toArray.map fun k => (dK[k]!, alpha[k]!, lens[k]!)
  let mut chanOut : Array (Array (Float × (Float × Float))) := #[]
  let mut cost := 0.0
  let mut memoryless := 0.0
  for c in chans do
    let ti := if memory then c.ti else (0.0, 0.0)
    let mut row : Array (Float × (Float × Float)) := #[]
    for k in [0:nbrs.size] do
      let n := nbrs[k]!
      let (mj, tj0) := c.nb[k]!
      let tj := if memory then tj0 else (0.0, 0.0)
      let rho := n.phiJI + floatPi - n.phiIJ
      let R := (Float.cos rho, Float.sin rho)
      let e := (lens[k]! * Float.cos n.phiJI, lens[k]! * Float.sin n.phiJI)
      let m' := iemMassUpdate mj alpha[k]! c.mi
      let t' := iemVectorUpdate mj alpha[k]! c.mi R ti e tj
      row := row.push (m', t')
      cost := cost + c.w * iemCostTerm m' t'
      memoryless := memoryless + c.w * alpha[k]! * c.mi * lens[k]!
    chanOut := chanOut.push row
  return { flat := fl, Ki, nbr := nbrOut, chan := chanOut, cost, memoryless }

end Ddg.IEM
