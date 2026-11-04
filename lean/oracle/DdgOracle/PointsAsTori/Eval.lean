import DdgOracle.PointsAsTori.Basic

/-!
# Points as Tori: neighborhoods, the blend (Eq. 25-26, App. D) and the SHC / SSPD baselines

SPEC §7.4 and §9. Sums run over the neighbor set in ascending `(ρ², index)` order.
-/

namespace Ddg.PointsAsTori

/-- Evaluation neighborhood (SPEC §7.4). `radius none` means `R_eval = 2C/λ` of the paper rule. -/
inductive Nbhd where
  | knn (k : Nat)
  | all
  | radius (r : Option Float) (fallback : Nat)

/-- The λ rule (SPEC §7.4). -/
inductive LamRule where
  | perQuery (c : Float)
  | constant (l : Float)
  | paper (c factor : Float) (kMean : Nat)

/-- Options of the evaluation. -/
structure EvalOpts where
  nbhd : Nbhd := .knn 32
  lam : LamRule := .perQuery 64.0
  areas : Bool := false
  iso : Float := 0.0

/-- `(ρ², index)` of every point, sorted ascending (ties by index). -/
def sortedDist2 (pts : Array V3) (x : V3) : Array (Float × Nat) :=
  let a := (Array.range pts.size).map fun j => ((pts[j]! - x).norm2, j)
  a.qsort fun u v => u.1 < v.1 || (u.1 == v.1 && u.2 < v.2)

/-- `D`: mean over the points of the mean distance to the `kMean` nearest points, self included. -/
def meanNeighborDistance (pts : Array V3) (kMean : Nat) : Float := Id.run do
  let n := pts.size
  let k := if kMean < n then kMean else n
  let mut tot := 0.0
  for i in [0:n] do
    let s := sortedDist2 pts pts[i]!
    let mut acc := 0.0
    for m in [0:k] do acc := acc + Float.sqrt s[m]!.1
    tot := tot + acc / k.toFloat
  return tot / n.toFloat

/-- A neighbor set with its shift and λ. -/
structure NSet where
  idx : Array Nat
  rho : Array Float
  sigma : Float
  lambda : Float
  /-- the exponent shift: `σ` for `PerQuery`, the nearest distance for the constant rules (no overflow) -/
  shift : Float

/-- Paper rule data: `(D, λ, R_eval)`. -/
def paperData (pts : Array V3) (c factor : Float) (kMean : Nat) : Float × Float × Float :=
  let D := meanNeighborDistance pts kMean
  let l := factor / D
  (D, l, 2.0 * c / l)

/-- The neighbor set, `σ_x` and `λ_x` of a query. -/
def neighborSet (pts : Array V3) (o : EvalOpts) (x : V3) : Except String NSet := do
  let s := sortedDist2 pts x
  let take := fun (k : Nat) => s.extract 0 (if k < s.size then k else s.size)
  let sel ← match o.nbhd with
    | .knn k => pure (take k)
    | .all => pure s
    | .radius r fb =>
      let rr ← match r with
        | some v => pure v
        | none => match o.lam with
          | .paper c f km => let (_, _, re) := paperData pts c f km; pure re
          | _ => throw "InvalidOptions: radius auto needs the paper lambda rule"
      let inside := s.filter fun u => Float.sqrt u.1 ≤ rr
      pure (if inside.isEmpty then take fb else inside)
  if sel.isEmpty then throw "DegenerateQuery: empty neighbor set"
  let idx := sel.map (·.2)
  let rho := sel.map fun u => Float.sqrt u.1
  let mut mx := 0.0
  for r in rho do if r > mx then mx := r
  let sigma := mx / 2.0
  match o.lam with
  | .perQuery _ => if !(sigma > 0.0) then throw "DegenerateQuery: shift sigma is 0 (lambda = C / sigma)"
  | _ => pure ()
  let lambda := match o.lam with
    | .perQuery c => c / sigma
    | .constant l => l
    | .paper _ f km => let (_, l, _) := paperData pts 0.0 f km; l
  let shift := match o.lam with
    | .perQuery _ => sigma
    | _ => rho[0]!
  return { idx, rho, sigma, lambda, shift }

/-- Result of the blend at one query. -/
structure BlendOut where
  phi : Float
  grad : V3
  gradExact : V3
  lap : Float
  nset : NSet
  singG : Nat
  singR : Nat

/-- The blend `φ`, the frozen gradient, the exact gradient (`PerQuery` only differs) and the
Laplacian with `λ`, `σ` frozen. `rh = 0` and the `2λ/ρ` term is dropped at `ρ = 0`. -/
def blend (pts : Array V3) (tori : Array Torus) (areas : Array Float) (o : EvalOpts) (x : V3) :
    Except String BlendOut := do
  let ns ← neighborSet pts o x
  let lam := ns.lambda
  let sig := ns.shift
  let mut W := 0.0
  let mut N := 0.0
  let mut singG := 0
  let mut singR := 0
  let mut ws : Array Float := #[]
  let mut gs : Array SdfOut := #[]
  let mut rhs : Array V3 := #[]
  for m in [0:ns.idx.size] do
    let i := ns.idx[m]!
    let rho := ns.rho[m]!
    let sd := tori[i]!.sdf x
    let w0 := Float.exp (-lam * (rho - sig))
    let w := if o.areas then areas[i]! * w0 else w0
    let rh := if rho == 0.0 then V3.zero else (x - pts[i]!).sdiv rho
    if rho == 0.0 then singR := singR + 1
    if sd.kind != 0 then singG := singG + 1
    W := W + w
    N := N + w * sd.g
    ws := ws.push w
    gs := gs.push sd
    rhs := rhs.push rh
  let phi0 := N / W
  let mut A := V3.zero
  let mut Bs := V3.zero
  let mut extra := 0.0
  let mut dN := 0.0
  let mut dW := 0.0
  for m in [0:ns.idx.size] do
    let w := ws[m]!
    let sd := gs[m]!
    let rh := rhs[m]!
    let rho := ns.rho[m]!
    A := A + (sd.grad - rh.smul (lam * sd.g)).smul w
    Bs := Bs + rh.smul w
    extra := extra + w * (sig - rho) * (sd.g - phi0)
    let lapRho := if rho == 0.0 then 0.0 else 2.0 * lam / rho
    dN := dN + w * ((lam * lam - lapRho) * sd.g - 2.0 * lam * rh.dot sd.grad + sd.lap)
    dW := dW + w * (lam * lam - lapRho)
  A := A.sdiv W
  let B := (Bs.smul lam).sdiv W
  let grad := A + B.smul phi0
  let last := ns.idx.size - 1
  let gradExact := match o.lam with
    | .perQuery c =>
      let dl := rhs[last]!.smul (-(c / (sig * sig)) * 0.5)
      grad + dl.smul (extra / W)
    | _ => grad
  let lap := (dN - phi0 * dW) / W + 2.0 * B.dot grad
  return { phi := phi0 - o.iso, grad, gradExact, lap, nset := ns, singG, singR }

/-- SHC (the shift plays the role of σ in the exponent and in the logarithm) (Eq. 28) and SSPD (Eq. 29) with their gradients (λ, σ frozen). -/
structure BaseOut where
  shc : Float
  shcGrad : V3
  sspd : Float
  sspdGrad : V3
  nset : NSet
  skipped : Nat

private def pi : Float := 3.141592653589793

/-- Terms with `ρ = 0` are left out of the SHC sum (and counted); in SSPD their `rh` is 0. -/
def baselines (pts nrm : Array V3) (areas : Array Float) (o : EvalOpts) (x : V3) :
    Except String BaseOut := do
  let ns ← neighborSet pts o x
  let lam := ns.lambda
  let sig := ns.shift
  let mut w := 0.0
  let mut gw := V3.zero
  let mut sn := 0.0
  let mut sW := 0.0
  let mut gn := V3.zero
  let mut gW := V3.zero
  let mut skipped := 0
  for m in [0:ns.idx.size] do
    let i := ns.idx[m]!
    let rho := ns.rho[m]!
    let d := x - pts[i]!
    let n := nrm[i]!
    let A := areas[i]!
    let e := Float.exp (-lam * (rho - sig))
    let dn := d.dot n
    -- SSPD
    let rh := if rho == 0.0 then V3.zero else d.sdiv rho
    sn := sn + A * dn * e
    sW := sW + A * e
    gn := gn + (n - rh.smul (lam * dn)).smul (A * e)
    gW := gW - rh.smul (lam * A * e)
    -- SHC (the shift plays the role of σ in the exponent and in the logarithm)
    if rho == 0.0 then skipped := skipped + 1
    else
      let rho3 := rho * rho * rho
      let rho5 := rho3 * rho * rho
      w := w + A * (lam * rho + 1.0) * dn / (2.0 * pi * rho3) * e
      let c1 := A * e / (2.0 * pi)
      gw := gw + (n.smul ((lam * rho + 1.0) / rho3) - d.smul (dn * (lam * lam * rho * rho + 3.0 * lam * rho + 3.0) / rho5)).smul c1
  if w == 0.0 then throw "DegenerateQuery: SHC sum is exactly 0"
  let shc := sgn0 w * (-(Float.log w.abs) / lam + sig)
  let shcGrad := gw.smul (-1.0 / (lam * w.abs))
  let sspd := sn / sW
  let sspdGrad := (gn - gW.smul sspd).sdiv sW
  return { shc, shcGrad, sspd, sspdGrad, nset := ns, skipped }

end Ddg.PointsAsTori
