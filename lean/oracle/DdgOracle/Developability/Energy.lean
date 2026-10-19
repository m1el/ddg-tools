import DdgOracle.Developability.Mesh

/-!
# The six developability energies and their (sub)gradients

Written from `crates/ddg-developability/SPEC.md` sections 1-3 only. Per vertex: the energy `e_v`, a
margin of the choice made (SPEC 3.8), the sensitivities `y_c = ∂e_v/∂N_σ`, `s_c = ∂e_v/∂θ_c` per
corner of the star and `z_v = ∂e_v/∂S_v`; the gradient is assembled from them face by face (SPEC 2).
-/

namespace Ddg.Developability

/-- The six energies. -/
inductive EnergyKind
  | combinatorial | combinatorialMax | covariance | covarianceMax | intrinsic | intrinsicMax
  deriving BEq, Inhabited

namespace EnergyKind

/-- Name in the text format. -/
def name : EnergyKind → String
  | .combinatorial => "combinatorial"
  | .combinatorialMax => "combinatorial-max"
  | .covariance => "covariance"
  | .covarianceMax => "covariance-max"
  | .intrinsic => "intrinsic"
  | .intrinsicMax => "intrinsic-max"

/-- Parse a name. -/
def ofName? (s : String) : Option EnergyKind :=
  [EnergyKind.combinatorial, .combinatorialMax, .covariance, .covarianceMax, .intrinsic, .intrinsicMax].find?
    (·.name == s)

/-- All six. -/
def all : Array EnergyKind :=
  #[.combinatorial, .combinatorialMax, .covariance, .covarianceMax, .intrinsic, .intrinsicMax]

end EnergyKind

/-- Options of an evaluation. -/
structure Opts where
  /-- the energy -/
  energy : EnergyKind := .covariance
  /-- count boundary vertices -/
  includeBoundary : Bool := true
  /-- leave out interior vertices with exactly three faces -/
  omitValenceThree : Bool := true
  deriving Inhabited

/-- NaN -/
def nan : Float := 0.0 / 0.0

/-- Result for one star. -/
structure StarResult where
  /-- `e_v` -/
  value : Float
  /-- margin of the choice -/
  margin : Float
  /-- the choice made, as integers (partition, candidate, pair) -/
  choice : Array Nat
  /-- `∂e_v/∂N_σ` per face of the star -/
  y : Array V3
  /-- `∂e_v/∂θ` per face of the star -/
  s : Array Float
  /-- `∂e_v/∂S_v` -/
  z : V3
  /-- `some reason` if the energy is undefined here -/
  undef : Option String := none
  deriving Inhabited

private def inf : Float := posInf

/-! ### Combinatorial energies -/

/-- The bipartitions of a star with `k` faces as pairs of index arrays, in the order of SPEC 3.2. -/
def bipartitions (k : Nat) (closed : Bool) : Array (Array Nat × Array Nat) := Id.run do
  let mut out := #[]
  if closed then
    for a in [0:k] do
      for b in [a + 1:k] do
        let f1 := (Array.range (b - a)).map (· + a + 1)
        let f2 := (Array.range k).filter fun i => !(a < i && i ≤ b)
        out := out.push (f1, f2)
  else
    for a in [0:k - 1] do
      out := out.push ((Array.range (a + 1)), (Array.range (k - 1 - a)).map (· + a + 1))
  return out

private def meanOf (ns : Array V3) (idx : Array Nat) : V3 :=
  V3.smul (1 / idx.size.toFloat) (idx.foldl (fun acc i => acc + ns[i]!) V3.zero)

private def partWidth (ns : Array V3) (idx : Array Nat) : Float :=
  let m := meanOf ns idx
  (idx.foldl (fun acc i => acc + (ns[i]! - m).norm2) 0) / idx.size.toFloat

private def partMax (ns : Array V3) (idx : Array Nat) : Float × Nat × Nat := Id.run do
  let mut best := (0.0, idx[0]!, idx[0]!)
  for p in [0:idx.size] do
    for q in [p + 1:idx.size] do
      let d := (ns[idx[p]!]! - ns[idx[q]!]!).norm2
      if d > best.1 then best := (d, idx[p]!, idx[q]!)
  return best

/-- Smallest value with its (first) index, and the smallest value among the others. -/
private def argminMargin (vals : Array Float) : Nat × Float × Float := Id.run do
  let mut bi := 0
  for i in [1:vals.size] do
    if vals[i]! < vals[bi]! then bi := i
  let mut second := inf
  for i in [0:vals.size] do
    if i != bi && vals[i]! < second then second := vals[i]!
  return (bi, vals[bi]!, second - vals[bi]!)

private def zeroResult (k : Nat) : StarResult :=
  { value := 0, margin := inf, choice := #[], y := Array.replicate k V3.zero,
    s := Array.replicate k 0, z := V3.zero }

/-- `E_P` of a star. -/
def combinatorialStar (ns : Array V3) (closed : Bool) : StarResult := Id.run do
  let k := ns.size
  let parts := bipartitions k closed
  if parts.size == 0 then return zeroResult k
  let vals := parts.map fun (f1, f2) => partWidth ns f1 + partWidth ns f2
  let (bi, best, margin) := argminMargin vals
  let (f1, f2) := parts[bi]!
  let mut y := Array.replicate k V3.zero
  for part in [f1, f2] do
    let m := meanOf ns part
    for i in part do
      y := y.set! i (V3.smul (2 / part.size.toFloat) (ns[i]! - m))
  return { value := best, margin, choice := #[bi], y, s := Array.replicate k 0, z := V3.zero }

/-- `π_max` of a star. -/
def combinatorialMaxStar (ns : Array V3) (closed : Bool) : StarResult := Id.run do
  let k := ns.size
  let parts := bipartitions k closed
  if parts.size == 0 then return zeroResult k
  let widths := parts.map fun (f1, f2) =>
    let w1 := partMax ns f1
    let w2 := partMax ns f2
    if w2.1 > w1.1 then w2 else w1
  let vals := widths.map (·.1)
  let (bi, best, margin) := argminMargin vals
  let (_, p, q) := widths[bi]!
  -- the active pair must also be unambiguous: second largest pair value in the chosen partition
  let (f1, f2) := parts[bi]!
  let mut other := 0.0
  let mut seen := false
  for part in [f1, f2] do
    for a in [0:part.size] do
      for b in [a + 1:part.size] do
        let pa := part[a]!
        let pb := part[b]!
        if !(pa == p && pb == q) then
          let d := (ns[pa]! - ns[pb]!).norm2
          if !seen || d > other then
            other := d
            seen := true
  let inner := if seen then best - other else inf
  let mut y := Array.replicate k V3.zero
  if p != q then
    let d := ns[p]! - ns[q]!
    y := (y.set! p (V3.smul 2 d)).set! q (V3.smul (-2) d)
  return { value := best, margin := min margin inner, choice := #[bi, p, q], y,
           s := Array.replicate k 0, z := V3.zero }

/-! ### Covariance energies -/

/-- `E_λ` of a star. -/
def covarianceStar (ns : Array V3) (th : Array Float) : StarResult := Id.run do
  let k := ns.size
  let mut a := Array.replicate 9 0.0
  for i in [0:k] do
    let n := ns[i]!
    let c := #[n.x, n.y, n.z]
    for r in [0:3] do
      for s in [0:3] do
        a := a.modify (r * 3 + s) (· + th[i]! * c[r]! * c[s]!)
  let (vals, vecs) := symEig 3 a
  let x : V3 := ⟨vecs[0]![0]!, vecs[0]![1]!, vecs[0]![2]!⟩
  let tr := th.foldl (· + ·) 0
  let margin := if tr > 0 then (vals[1]! - vals[0]!) / tr else inf
  let mut y := Array.replicate k V3.zero
  let mut s := Array.replicate k 0.0
  for i in [0:k] do
    let d := V3.dot x ns[i]!
    y := y.set! i (V3.smul (2 * th[i]! * d) x)
    s := s.set! i (d * d)
  return { value := vals[0]!, margin, choice := #[], y, s, z := V3.zero }

/-- `λ_max` of a star. -/
def covarianceMaxStar (ns : Array V3) : StarResult := Id.run do
  let k := ns.size
  -- candidates: (index, sites as (face, sign) triple)
  let mut cands : Array (Array (Nat × Float)) := #[]
  for i in [0:k] do
    for j in [i + 1:k] do
      for l in [j + 1:k] do
        for (s2, s3) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] do
          cands := cands.push #[(i, 1.0), (j, s2), (l, s3)]
  for i in [0:k] do
    for j in [i + 1:k] do
      cands := cands.push #[(i, 1.0), (i, -1.0), (j, 1.0)]
  let mut best : Option (Float × Float × Nat) := none
  let mut second := inf
  for ci in [0:cands.size] do
    let cd := cands[ci]!
    let a := V3.smul cd[0]!.2 ns[cd[0]!.1]!
    let b := V3.smul cd[1]!.2 ns[cd[1]!.1]!
    let c := V3.smul cd[2]!.2 ns[cd[2]!.1]!
    let n := V3.cross (b - a) (c - a)
    let nn := n.norm
    if nn.isFinite && nn != 0 then
      let v := V3.smul (1 / nn) n
      let mut phi := 0.0
      for i in [0:k] do
        let d := V3.dot v ns[i]!
        phi := max phi (d * d)
      let da := V3.dot v a
      let slack := phi - da * da
      match best with
      | none => best := some (phi, slack, ci)
      | some (bp, bs, _) =>
        if phi < bp || (phi == bp && slack < bs) then
          second := min second bp
          best := some (phi, slack, ci)
        else
          second := min second phi
  match best with
  | none => return zeroResult k
  | some (phi, _, ci) =>
    let cd := cands[ci]!
    let a := V3.smul cd[0]!.2 ns[cd[0]!.1]!
    let b := V3.smul cd[1]!.2 ns[cd[1]!.1]!
    let c := V3.smul cd[2]!.2 ns[cd[2]!.1]!
    let n := V3.cross (b - a) (c - a)
    let nn := n.norm
    let v := V3.smul (1 / nn) n
    let g := V3.dot a (V3.cross b c) / nn
    let edges := #[c - b, a - c, b - a]
    let mut y := Array.replicate k V3.zero
    for q in [0:3] do
      let eq := edges[q]!
      let ind := if q == 0 then v else V3.zero
      let dq := V3.smul (2 * g) (ind + V3.smul (V3.dot (V3.cross eq v) a / nn) v)
      let (face, sg) := cd[q]!
      y := y.modify face (· + V3.smul sg dq)
    return { value := phi, margin := second - phi, choice := #[ci], y, s := Array.replicate k 0, z := V3.zero }

/-! ### Intrinsic energies -/

/-- `G(φ) = φ / sin φ` and `H(φ) = (sin φ - φ cos φ) / sin³ φ`, with the series below `1e-3`. -/
def logGH (phi sinPhi cosPhi : Float) : Float × Float :=
  if phi < 1e-3 then (1 + phi * phi / 6 + 7 * phi ^ 4 / 360, 1.0 / 3 + 2 * phi * phi / 15)
  else (phi / sinPhi, (sinPhi - phi * cosPhi) / (sinPhi * sinPhi * sinPhi))

/-- Orthonormal tangent basis of a unit vector. -/
def tangentBasis (n0 : V3) : V3 × V3 :=
  let ax : V3 := if n0.x.abs ≤ n0.y.abs && n0.x.abs ≤ n0.z.abs then ⟨1, 0, 0⟩
    else if n0.y.abs ≤ n0.z.abs then ⟨0, 1, 0⟩ else ⟨0, 0, 1⟩
  let t1 := V3.unit (V3.cross n0 ax)
  (t1, V3.cross n0 t1)

/-- Common part of the intrinsic energies: the log-mapped normals; `Except` the reason it is
undefined. -/
def logNormals (ns : Array V3) (sv : V3) :
    Except String (V3 × Array V3 × Array (Float × Float × Float × Float)) := do
  let sn := sv.norm
  if !(sn.isFinite) || sn == 0 then throw "zero vertex normal"
  let n0 := V3.smul (1 / sn) sv
  let mut nt : Array V3 := #[]
  let mut gh : Array (Float × Float × Float × Float) := #[]
  for N in ns do
    let c := V3.dot n0 N
    let sinPhi := (V3.cross n0 N).norm
    let phi := Float.atan2 sinPhi c
    if 3.141592653589793 - phi < 1e-6 then throw "normal opposite to the vertex normal"
    let (g, h) := logGH phi sinPhi c
    nt := nt.push (V3.smul g (N - V3.smul c n0))
    gh := gh.push (g, h, c, phi)
  return (n0, nt, gh)

/-- Back-propagate sensitivities `zs` w.r.t. the log-mapped normals to `y` (w.r.t. `N`) and `z_v`. -/
def logAdjoint (ns : Array V3) (n0 : V3) (sn : Float) (gh : Array (Float × Float × Float × Float))
    (zs : Array V3) : Array V3 × V3 := Id.run do
  let mut y : Array V3 := #[]
  let mut yn0 := V3.zero
  for i in [0:ns.size] do
    let (g, h, c, _) := gh[i]!
    let z := zs[i]!
    let N := ns[i]!
    let kk := g * V3.dot z n0 + h * V3.dot z (N - V3.smul c n0)
    y := y.push (V3.smul g z - V3.smul kk n0)
    yn0 := yn0 + (V3.smul (-g * c) z - V3.smul kk N)
  let zv := V3.smul (1 / sn) (yn0 - V3.smul (V3.dot yn0 n0) n0)
  return (y, zv)

/-- Intrinsic energy `λ_min` of the log-mapped covariance. -/
def intrinsicStar (ns : Array V3) (th : Array Float) (sv : V3) : StarResult := Id.run do
  let k := ns.size
  match logNormals ns sv with
  | .error r => return { zeroResult k with value := nan, margin := nan, undef := some r }
  | .ok (n0, nt, gh) =>
    let (t1, t2) := tangentBasis n0
    let mut m := #[0.0, 0.0, 0.0, 0.0]
    for i in [0:k] do
      let p := V3.dot nt[i]! t1
      let q := V3.dot nt[i]! t2
      m := #[m[0]! + th[i]! * p * p, m[1]! + th[i]! * p * q, m[2]! + th[i]! * q * p, m[3]! + th[i]! * q * q]
    let (vals, vecs) := symEig 2 m
    let x := V3.smul vecs[0]![0]! t1 + V3.smul vecs[0]![1]! t2
    let tr := th.foldl (· + ·) 0
    let margin := if tr > 0 then (vals[1]! - vals[0]!) / tr else inf
    let mut zs : Array V3 := #[]
    let mut s : Array Float := #[]
    for i in [0:k] do
      let d := V3.dot x nt[i]!
      zs := zs.push (V3.smul (2 * th[i]! * d) x)
      s := s.push (d * d)
    let (y, zv) := logAdjoint ns n0 sv.norm gh zs
    return { value := vals[0]!, margin, choice := #[], y, s, z := zv }

/-- Branching-free intrinsic width. -/
def intrinsicMaxStar (ns : Array V3) (sv : V3) : StarResult := Id.run do
  let k := ns.size
  match logNormals ns sv with
  | .error r => return { zeroResult k with value := nan, margin := nan, undef := some r }
  | .ok (n0, nt, gh) =>
    let mut best : Option (Float × Nat × Nat) := none
    let mut second := inf
    for a in [0:k] do
      for b in [a + 1:k] do
        let w := nt[b]! - nt[a]!
        let wn := w.norm
        if wn != 0 && wn.isFinite then
          let u := V3.cross n0 (V3.smul (1 / wn) w)
          let mut mx := -inf
          let mut mn := inf
          for i in [0:k] do
            let d := V3.dot nt[i]! u
            mx := max mx d
            mn := min mn d
          let psi := (mx - mn) * (mx - mn)
          match best with
          | none => best := some (psi, a, b)
          | some (bp, _, _) =>
            if psi < bp then
              second := min second bp
              best := some (psi, a, b)
            else second := min second psi
    match best with
    | none => return { zeroResult k with }
    | some (psi, a, b) =>
      let w := nt[b]! - nt[a]!
      let wn := w.norm
      let wh := V3.smul (1 / wn) w
      let u := V3.cross n0 wh
      let mut s := 0
      let mut t := 0
      for i in [0:k] do
        if V3.dot nt[i]! u > V3.dot nt[s]! u then s := i
        if V3.dot nt[i]! u < V3.dot nt[t]! u then t := i
      let d := nt[s]! - nt[t]!
      let delta := V3.dot d u
      let mut zs := Array.replicate k V3.zero
      zs := zs.modify s (· + V3.smul (2 * delta) u)
      zs := zs.modify t (· - V3.smul (2 * delta) u)
      let ca := V3.smul (2 * delta * V3.dot d wh / wn) u
      zs := zs.modify a (· + ca)
      zs := zs.modify b (· - ca)
      let (y, zv) := logAdjoint ns n0 sv.norm gh zs
      return { value := psi, margin := second - psi, choice := #[a, b], y,
               s := Array.replicate k 0, z := zv }

/-! ### Evaluation of a mesh -/

/-- Why an evaluation failed. -/
inductive EvalError
  | degenerate (face : Nat)
  | undefinedAt (vertex : Nat) (reason : String)
  | invalid (msg : String)
  | nearDegenerate (face : Nat)

/-- Text of the error. -/
def EvalError.text : EvalError → String
  | .degenerate f => s!"degenerate-face {f}"
  | .undefinedAt v r => s!"undefined {v} {r.replace " " "-"}"
  | .invalid m => s!"invalid {m}"
  | .nearDegenerate f => s!"near-degenerate-face {f}"

/-- Result of an evaluation. -/
structure Eval where
  /-- `E` -/
  total : Float
  /-- `e_v` -/
  e : Array Float
  /-- margins -/
  margin : Array Float
  /-- in the domain -/
  counted : Array Bool
  /-- valences -/
  valence : Array Nat
  /-- boundary flag -/
  boundary : Array Bool
  /-- choice per vertex -/
  choice : Array (Array Nat)
  /-- gradient -/
  grad : Array V3

/-- Evaluate energy and subgradient. -/
def evaluate (m : DMesh) (topo : Topo) (o : Opts) : Except EvalError Eval := do
  let nv := topo.nv
  let mut fd : Array FaceData := #[]
  for fi in [0:m.faces.size] do
    match faceData m fi with
    | some d => fd := fd.push d
    | none => throw (.degenerate fi)
  let mut results : Array StarResult := #[]
  let mut counted : Array Bool := #[]
  for v in [0:nv] do
    let st := topo.stars[v]!
    let k := st.faces.size
    let ns := st.faces.map fun f => fd[f]!.N
    let th := Array.ofFn (n := k) fun i => fd[st.faces[i]!]!.theta[st.slots[i]!]!
    let sv := V3.smul 0.5 (st.faces.foldl (fun acc f => acc + fd[f]!.n) V3.zero)
    let r := match o.energy with
      | .combinatorial => combinatorialStar ns st.closed
      | .combinatorialMax => combinatorialMaxStar ns st.closed
      | .covariance => covarianceStar ns th
      | .covarianceMax => covarianceMaxStar ns
      | .intrinsic => intrinsicStar ns th sv
      | .intrinsicMax => intrinsicMaxStar ns sv
    -- SPEC D33: a boundary star with at most two faces contributes exactly nothing (value 0, no
    -- sensitivities, nothing to tie), except where the intrinsic energies are undefined (D12).
    let r := if !st.closed && k ≤ 2 && r.undef.isNone then zeroResult k else r
    let inDomain := (o.includeBoundary || st.closed) && !(o.omitValenceThree && st.closed && k == 3)
    if inDomain then
      if let some reason := r.undef then throw (.undefinedAt v reason)
    results := results.push r
    counted := counted.push inDomain
  -- gradient: corner sensitivities, then per face
  let nf := m.faces.size
  let mut yc : Array V3 := Array.replicate (3 * nf) V3.zero
  let mut sc : Array Float := Array.replicate (3 * nf) 0
  let mut zv : Array V3 := Array.replicate nv V3.zero
  let mut total := 0.0
  for v in [0:nv] do
    if counted[v]! then
      let st := topo.stars[v]!
      let r := results[v]!
      total := total + r.value
      zv := zv.set! v r.z
      for i in [0:st.faces.size] do
        yc := yc.set! (3 * st.faces[i]! + st.slots[i]!) r.y[i]!
        sc := sc.set! (3 * st.faces[i]! + st.slots[i]!) r.s[i]!
  let mut grad : Array V3 := Array.replicate nv V3.zero
  for fi in [0:nf] do
    let f := m.faces[fi]!
    let vs := #[f.1, f.2.1, f.2.2]
    let d := fd[fi]!
    let ysum := yc[3 * fi]! + yc[3 * fi + 1]! + yc[3 * fi + 2]!
    let zsum := zv[vs[0]!]! + zv[vs[1]!]! + zv[vs[2]!]!
    let w := V3.smul (1 / d.nn) (ysum - V3.smul (V3.dot ysum d.N) d.N) + V3.smul 0.5 zsum
    for p in [0:3] do
      let e := m.pos[vs[(p + 2) % 3]!]! - m.pos[vs[(p + 1) % 3]!]!
      grad := grad.modify vs[p]! (· + V3.cross w e)
    for c in [0:3] do
      let i := vs[c]!
      let j := vs[(c + 1) % 3]!
      let kk := vs[(c + 2) % 3]!
      let pi := m.pos[i]!
      let pj := m.pos[j]!
      let pk := m.pos[kk]!
      let sg := sc[3 * fi + c]!
      let gj := V3.smul (sg / (pi - pj).norm2) (V3.cross d.N (pi - pj))
      let gk := V3.smul (sg / (pk - pi).norm2) (V3.cross d.N (pk - pi))
      grad := grad.modify j (· + gj)
      grad := grad.modify kk (· + gk)
      grad := grad.modify i (· - (gj + gk))
  return { total, e := results.map (·.value), margin := results.map (·.margin), counted,
           valence := topo.stars.map (·.faces.size), boundary := topo.stars.map (!·.closed),
           choice := results.map (·.choice), grad }

end Ddg.Developability
