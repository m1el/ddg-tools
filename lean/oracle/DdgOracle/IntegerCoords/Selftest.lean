import DdgOracle.IntegerCoords.Cases

/-!
# Self-test of the `ic` oracle

Cheap runtime checks of the integer formulas, independent of the Lean proofs (which prove the
general statements): over small exhaustive ranges,

* Eq. 3 applied to the new edge returns the old count (the flip involution), and its doubled value
  is even;
* the decoded pairing of every valid face (`partner`) is an involution on the crossings, normal
  (other side or opposite corner), and its arcs do not interleave on the boundary cycle;
  `decode?` fails exactly on invalid faces;
* the face-split counts of every region are the numbers of arcs separating a point of that region
  from each corner, counted by brute force on the boundary cycle, and the three new faces are valid;
* Alg. 1's step agrees with the pairing.
-/

namespace Ddg.IntegerCoords

/-- Rank of a boundary point on the cycle corner 0, side 0, corner 1, side 1, corner 2, side 2. -/
def cornerRank (n : Nat → Nat) (k : Nat) : Nat :=
  if k == 0 then 0 else if k == 1 then n 0 + 1 else n 0 + n 1 + 2

/-- Rank of a boundary point. -/
def rank (n : Nat → Nat) : Pt → Nat
  | .corner k => cornerRank n k
  | .side k p => cornerRank n k + 1 + p

/-- `x` strictly between `a` and `b`. -/
def btw (a b x : Nat) : Bool := (a < x && x < b) || (b < x && x < a)

/-- The chords `{a, b}` and `{x, y}` interleave. -/
def interleave (a b x y : Nat) : Bool :=
  x != a && x != b && y != a && y != b && (btw a b x != btw a b y)

/-- The arcs of a decoded face as rank pairs, each once. -/
def arcs (d : Decoded) : List (Nat × Nat) := Id.run do
  let mut out := []
  for k in [0:3] do
    for p in [0:d.n k] do
      let a := rank d.n (.side k p)
      let b := rank d.n (partner d.n d.c k p)
      match partner d.n d.c k p with
      | .corner _ => out := (a, b) :: out
      | .side _ _ => if a < b then out := (a, b) :: out
  return out

/-- A gap (between ranks `g` and `g + 1`) inside a region. -/
def regionGap (d : Decoded) : Region → Nat
  | .corner k t => rank d.n (.corner k) + t
  | .middle s =>
    if d.e > 0 then rank d.n (.corner (next d.kappa)) + d.c (next d.kappa) + s
    else d.c 0

/-- Arcs separating gap `g` from corner `j`, by brute force. -/
def separating (d : Decoded) (g j : Nat) : Nat :=
  let r := rank d.n (.corner j)
  (arcs d).filter (fun (a, b) =>
    let lo := min a b
    let hi := max a b
    r != a && r != b && ((lo ≤ g && g < hi) != (lo < r && r < hi))) |>.length

/-- Run the self-test; returns the number of failures. -/
def selftest : IO Nat := do
  let mut bad := 0
  let mut checked := 0
  -- flip involution
  let rng := (List.range 7).map fun (i : Nat) => (i : Int) - 1
  for a in rng do for b in rng do for c in rng do for d in rng do for e in rng do
    if faceValid a b c && faceValid a d e then
      checked := checked + 1
      let x := flipNormal a b c d e
      if flipDoubled a b c d e % 2 != 0 then
        bad := bad + 1; IO.println s!"FAIL  flipDoubled odd for {a} {b} {c} {d} {e}"
      if flipNormal x c d e b != a then
        bad := bad + 1; IO.println s!"FAIL  flip back of {a} {b} {c} {d} {e} gives {flipNormal x c d e b}"
      if !(faceValid c d x && faceValid e b x) then
        bad := bad + 1; IO.println s!"FAIL  new faces of {a} {b} {c} {d} {e} not valid"
  IO.println s!"flip involution: {checked} valid pairs"
  -- faces
  let mut nfaces := 0
  for (n0, n1, n2) in triples (-1) 8 do
    match decode? n0 n1 n2 with
    | none =>
      if faceValid n0 n1 n2 then bad := bad + 1; IO.println s!"FAIL  decode? none on valid {n0} {n1} {n2}"
    | some d =>
      nfaces := nfaces + 1
      if !faceValid n0 n1 n2 then bad := bad + 1; IO.println s!"FAIL  decode? on invalid {n0} {n1} {n2}"
      let pre := s!"face {n0} {n1} {n2}"
      -- pairing
      for k in [0:3] do
        for p in [0:d.n k] do
          match partner d.n d.c k p with
          | .side k' q =>
            if k' == k || q ≥ d.n k' then bad := bad + 1; IO.println s!"FAIL  {pre}: arc {k} {p} not normal"
            else if partner d.n d.c k' q != .side k p then
              bad := bad + 1; IO.println s!"FAIL  {pre}: partner not an involution at {k} {p}"
          | .corner k' =>
            if k' != prev k then bad := bad + 1; IO.println s!"FAIL  {pre}: arc {k} {p} to a wrong corner"
      let as := arcs d
      for (a, b) in as do
        for (x, y) in as do
          if interleave a b x y then
            bad := bad + 1; IO.println s!"FAIL  {pre}: arcs ({a},{b}) and ({x},{y}) interleave"
      -- regions
      for r in regions d do
        match faceSplit? n0 n1 n2 r with
        | none => bad := bad + 1; IO.println s!"FAIL  {pre}: no counts for {r.toStr}"
        | some (m0, m1, m2) =>
          let ms : Nat → Nat := fun j => if j == 0 then m0 else if j == 1 then m1 else m2
          let g := regionGap d r
          for j in [0:3] do
            if separating d g j != ms j then
              bad := bad + 1
              IO.println s!"FAIL  {pre} {r.toStr}: count to corner {j} is {ms j}, brute force {separating d g j}"
          for k in [0:3] do
            if !faceValid (d.n k) (ms (next k)) (ms k) then
              bad := bad + 1; IO.println s!"FAIL  {pre} {r.toStr}: new face on side {k} not valid"
      -- tracing
      for s in [0:3] do
        for p in [0:d.n s] do
          let want := partner d.n d.c s (d.n s - 1 - p)
          let got := match traceStep (d.n s) (d.n (next s)) (d.n (prev s)) p with
            | .toIK q => Pt.side (next s) q.toNat
            | .toKJ q => Pt.side (prev s) q.toNat
            | .vertex => Pt.corner (prev s)
          if got != want then
            bad := bad + 1; IO.println s!"FAIL  {pre}: trace step through side {s} at {p}"
  IO.println s!"faces: {nfaces} valid faces with n in [-1, 8]"
  -- the fixtures run without an error block
  for (c, _) in allCases do
    if (runCase c).startsWith s!"case {c.name}\ntask {c.task}\nstatus error" then
      bad := bad + 1; IO.println s!"FAIL  case {c.name} gives an error block"
  IO.println (if bad == 0 then "ic selftest: all checks passed" else s!"ic selftest: {bad} failures")
  return bad

end Ddg.IntegerCoords
