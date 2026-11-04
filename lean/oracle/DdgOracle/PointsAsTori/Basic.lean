import DdgOracle.Hex

/-!
# Points as Tori: features, torus fitting, per-torus SDF

Plain `Float`, written from `crates/ddg-points-as-tori/SPEC.md` §4.1, §4.2, §7.1, §7.2 (formulas
`pat-median-scale`, `pat-local-frame`, `pat-features`, `pat-output-scaling`, `pat-fundamental-forms`
... `pat-torus-frame`, `pat-torus-sdf`, `pat-torus-sdf-gradient`, `pat-torus-sdf-laplacian`) and
`NOTES-REFERENCE.md` §6.1 (clamps). It shares no code with the Rust crate.
-/

namespace Ddg.PointsAsTori

/-- A 3-vector. -/
structure V3 where
  x : Float
  y : Float
  z : Float
  deriving Inhabited

namespace V3
def zero : V3 := ⟨0.0, 0.0, 0.0⟩
def add (a b : V3) : V3 := ⟨a.x + b.x, a.y + b.y, a.z + b.z⟩
def sub (a b : V3) : V3 := ⟨a.x - b.x, a.y - b.y, a.z - b.z⟩
def smul (s : Float) (a : V3) : V3 := ⟨s * a.x, s * a.y, s * a.z⟩
/-- Componentwise division by a scalar (not multiplication by the reciprocal). -/
def sdiv (a : V3) (s : Float) : V3 := ⟨a.x / s, a.y / s, a.z / s⟩
def dot (a b : V3) : Float := a.x * b.x + a.y * b.y + a.z * b.z
def cross (a b : V3) : V3 :=
  ⟨a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x⟩
def norm2 (a : V3) : Float := dot a a
def norm (a : V3) : Float := Float.sqrt (dot a a)
def isFinite (a : V3) : Bool := a.x.isFinite && a.y.isFinite && a.z.isFinite
end V3

instance : Add V3 := ⟨V3.add⟩
instance : Sub V3 := ⟨V3.sub⟩

/-- `+1` for `x ≥ 0`, else `−1` (the C++ `sign`; NaN gives `−1`). -/
def sgnP (x : Float) : Float := if x ≥ 0.0 then 1.0 else -1.0

/-- `+1`, `−1`, `0` (`jnp.sign`). -/
def sgn0 (x : Float) : Float := if x > 0.0 then 1.0 else if x < 0.0 then -1.0 else 0.0

/-- `pat-local-frame`: `(ŝ, t̂)` from the unit normal `c`. -/
def frame (c : V3) : V3 × V3 :=
  let u : V3 := if c.x.abs > 0.9 then ⟨0.0, 1.0, 0.0⟩ else ⟨1.0, 0.0, 0.0⟩
  let cu := c.cross u
  let s := cu.sdiv cu.norm
  (s, c.cross s)

/-- `n / max(‖n‖, 1e-10)`. -/
def unitGuard (n : V3) : V3 :=
  let l := n.norm
  n.sdiv (if l > 1.0e-10 then l else 1.0e-10)

/-- Median of a nonempty list: the middle element, or the mean of the two middle ones (`pat-median-scale`). -/
def median (a : Array Float) : Float :=
  let s := a.qsort (· < ·)
  let k := s.size
  if k % 2 == 1 then s[k / 2]! else (s[k / 2 - 1]! + s[k / 2]!) / 2.0

/-- Features of one neighborhood. -/
structure Features where
  sigma : Float
  s : V3
  t : V3
  c : V3
  rows : Array (Array Float)

/-- `pat-features`: neighborhood `idx` (first entry is the center `p_i`), points and normals. -/
def features (pts nrm : Array V3) (idx : Array Nat) : Features :=
  let p0 := pts[idx[0]!]!
  let n0 := unitGuard nrm[idx[0]!]!
  let (s, t) := frame n0
  let sigma := median (idx.map fun j => (pts[j]! - p0).norm)
  let rows := idx.map fun j =>
    let d := (pts[j]! - p0).sdiv sigma
    let nj := unitGuard nrm[j]!
    #[d.dot s, d.dot t, d.dot n0, nj.dot s, nj.dot t, nj.dot n0]
  { sigma, s, t, c := n0, rows }

/-- `pat-output-scaling`: network output `y` to coefficients `(a00, a01, a10, a11, a02, a20)`. -/
def outputScaling (sigma : Float) (y : Array Float) : Array Float :=
  #[sigma * y[0]!, y[1]!, y[2]!, y[3]! / sigma, y[4]! / sigma, y[5]! / sigma]

/-- Options of the fit (SPEC §7.1). -/
structure FitOpts where
  eps : Float := 1.0e-10
  floor : Float := 1.0e-10

/-- A torus: center, unit axis, major radius `R`, signed minor radius `r`. -/
structure Torus where
  c : V3
  u : V3
  R : Float
  r : Float
  deriving Inhabited

/-- Result of fitting one torus, with the diagnostics of §7.1. -/
structure FitResult where
  torus : Torus
  sign : Float
  A : Float
  H : Float
  K : Float
  kPaperPlus : Float
  kPaperMinus : Float
  kPlus : Float
  kMinus : Float
  kMin : Float
  kMax : Float
  vMin : V3
  vMax : V3
  nStar : V3
  clampDisc : Bool
  clampEig : Bool
  clampDenom : Bool
  clampCurv : Bool
  umbilic : Bool
  finite : Bool
  deriving Inhabited

/-- `x` with its sign kept and its magnitude raised to at least `e` (`sgn⁺(x) · max(|x|, e)`). -/
private def clampSign (x e : Float) : Float := sgnP x * (if x.abs > e then x.abs else e)

/-- Fit the torus at `p` with unit normal `n` from the coefficients
`a = (a00, a01, a10, a11, a02, a20)`. JAX clamp rule (SPEC D9). -/
def fitTorus (o : FitOpts) (p n : V3) (a : Array Float) : FitResult :=
  let a00 := a[0]!
  let a01 := a[1]!
  let a10 := a[2]!
  let a11 := a[3]!
  let a02 := a[4]!
  let a20 := a[5]!
  let (sh, th) := frame n
  let A2 := 1.0 + a01 * a01 + a10 * a10
  let A := Float.sqrt A2
  let v1 := a20 * a01 * a01 - a02 * a10 * a10 - a02 + a20
  let v2 := a11 * (a10 * a10 + 1.0) - 2.0 * a01 * a10 * a20
  let term := a02 * (a10 * a10 + 1.0) - a01 * a10 * a11 + a20 * (a01 * a01 + 1.0)
  let discArg := A2 * (a11 * a11 - 4.0 * a02 * a20) + term * term
  let clampDisc := discArg < o.eps
  let disc := Float.sqrt (if clampDisc then o.eps else discArg)
  let e0raw := v1 + disc
  let clampE0 := e0raw.abs < o.eps
  let e0 := clampSign e0raw o.eps
  let e1raw := v2
  let clampE1 := e1raw.abs < o.eps
  let e1 := clampSign e1raw o.eps
  -- ev₊ = (e0, v2), ev₋ = (v1 − disc, e1)
  let epx0 := e0
  let epy0 := v2
  let enx0 := v1 - disc
  let eny0 := e1
  let lp := Float.sqrt (epx0 * epx0 + epy0 * epy0)
  let ln := Float.sqrt (enx0 * enx0 + eny0 * eny0)
  let epx := epx0 / lp
  let epy := epy0 / lp
  let enx := enx0 / ln
  let eny := eny0 / ln
  let dS := sh + n.smul a10
  let dT := th + n.smul a01
  let vPos := dS.smul epx + dT.smul epy
  let vNeg := dS.smul enx + dT.smul eny
  let nStar := (n - sh.smul a10 - th.smul a01).sdiv A
  let i11 := 1.0 + a10 * a10
  let i12 := a10 * a01
  let i22 := 1.0 + a01 * a01
  let ii11 := 2.0 * a20 / A
  let ii12 := a11 / A
  let ii22 := 2.0 * a02 / A
  let rayleigh := fun (ex ey : Float) =>
    let numer := ex * (ii11 * ex + ii12 * ey) + ey * (ii12 * ex + ii22 * ey)
    let denom := ex * (i11 * ex + i12 * ey) + ey * (i12 * ex + i22 * ey)
    (numer, denom)
  let (np, dp) := rayleigh epx epy
  let (nn, dn) := rayleigh enx eny
  let clampDenom := dp.abs < o.eps || dn.abs < o.eps
  let kPlus := np / clampSign dp o.eps
  let kMinus := nn / clampSign dn o.eps
  let usePos := kPlus.abs < kMinus.abs
  let kMin := if usePos then kPlus else kMinus
  let kMax := if usePos then kMinus else kPlus
  let vMin := if usePos then vPos else vNeg
  let vMax := if usePos then vNeg else vPos
  let aMax := if kMax.abs > o.floor then kMax.abs else o.floor
  let aMin := if kMin.abs > o.floor then kMin.abs else o.floor
  let clampCurv := kMax.abs < o.floor || kMin.abs < o.floor
  let r := 1.0 / aMax
  let R := 1.0 / aMin - sgn0 (kMax * kMin) * r
  let s := if (kMax < 0.0 && kMin < 0.0) || (kMax < 0.0 && kMin > 0.0) then 1.0 else -1.0
  let center := p + n.smul a00 - nStar.smul (s / aMin)
  let uRaw := nStar.cross vMin
  let un := uRaw.norm
  let u := uRaw.sdiv (if un > o.eps then un else o.eps)
  -- the paper's Eq. 19, independent of the eigenvectors
  let detI := A2
  let H := ((i22 * ii11 - 2.0 * i12 * ii12 + i11 * ii22) / detI) / 2.0
  let K := (ii11 * ii22 - ii12 * ii12) / detI
  let h2 := H * H - K
  let sq := Float.sqrt (if h2 > 0.0 then h2 else 0.0)
  let t : Torus := { c := center, u, R, r := s * r }
  let vminHat := vMin.sdiv vMin.norm
  let vmaxHat := vMax.sdiv vMax.norm
  let fin := center.isFinite && u.isFinite && R.isFinite && r.isFinite
  { torus := t, sign := s, A, H, K, kPaperPlus := H + sq, kPaperMinus := H - sq, kPlus, kMinus, kMin, kMax
    vMin := vminHat, vMax := vmaxHat, nStar
    clampDisc, clampEig := clampE0 || clampE1, clampDenom, clampCurv, umbilic := clampDisc, finite := fin }

/-- Output of the per-torus SDF at a point. `kind`: 0 regular, 1 on the axis, 2 on the core circle. -/
structure SdfOut where
  g : Float
  grad : V3
  lap : Float
  kind : Nat
  deriving Inhabited

/-- `pat-torus-sdf`, `pat-torus-sdf-gradient`, `pat-torus-sdf-laplacian` (SPEC §7.2).
On the core circle (`‖d‖ = 0`, checked first) the gradient and the Laplacian are 0; on the axis
(`v = 0`) the gradient is `sgn(r) · z/‖d‖ · u` and the Laplacian is 0. -/
def Torus.sdf (T : Torus) (x : V3) : SdfOut :=
  let xc := x - T.c
  let v := (xc.cross T.u).norm
  let z := xc.dot T.u
  let d1 := v - T.R
  let dn := Float.sqrt (d1 * d1 + z * z)
  let sg := sgnP T.r
  let g := sg * (dn - T.r.abs)
  if dn == 0.0 then { g, grad := V3.zero, lap := 0.0, kind := 2 }
  else if v == 0.0 then { g, grad := T.u.smul (sg * (z / dn)), lap := 0.0, kind := 1 }
  else
    let G1 := (xc - T.u.smul z).sdiv v
    let w1 := d1 / dn
    let w2 := z / dn
    let grad := (G1.smul w1 + T.u.smul w2).smul sg
    let lap := sg * (1.0 / dn + d1 / (dn * v))
    { g, grad, lap, kind := 0 }

end Ddg.PointsAsTori
