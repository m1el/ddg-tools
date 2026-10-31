import DdgOracle.IntegerCoords.Format

/-!
# Cases of the `ic` oracle: hand-checked cases from the paper's figures, and generated ones

The figure cases were read off the figures of Gillespie, Sharp, Crane 2021 (rendered at 200 dpi);
`Fixtures.lean` lists, for each, what was checked by hand. The generated ones are exhaustive over
small counts or drawn from a fixed linear congruential generator (deterministic).
-/

namespace Ddg.IntegerCoords

/-- A case from its name, task and body lines (given as strings, split at spaces). -/
def mkCase (name task : String) (lines : List String) : ICCase :=
  { name, task, lines := (lines.map fun l => (l.splitOn " ").filter (· ≠ "") |>.toArray).toArray }

/-! ## Fig. 5 (right): a complete intrinsic triangulation with integer data

Input `T⁰`: the rectangle `[0, 2] × [0, 1]` with the vertices `TL 0, TM 1, TR 2, BL 3, BM 4, BR 5`,
its boundary, the middle vertical `TM–BM` and the diagonals `BL–TM`, `BM–TR`. `T¹`: an inserted
vertex `m = 6` inside the left square and `r = 7` on the right side (an edge split, so `TR–r` and
`r–BR` are `−1`); the labels of the figure are the normal coordinates. Faces (corners in
counter-clockwise order, halfedge `3f + k` from corner `k` to `k + 1`):

    F0 (m, BM, BR)  n = 1, −1, 3       F4 (BL, BM, m)  n = −1, 1, 0
    F1 (m, BR, r)   n = 3, −1, 3       F5 (TL, BL, m)  n = −1, 0, 0
    F2 (m, r, TM)   n = 3, 1, 0        F6 (TL, m, TM)  n = 0, 0, −1
    F3 (TM, r, TR)  n = 1, −1, −1
-/

/-- The faces of Fig. 5 as `face` lines. -/
def fig5Faces : List String :=
  ["face 1 -1 3", "face 3 -1 3", "face 3 1 0", "face 1 -1 -1", "face -1 1 0", "face -1 0 0",
    "face 0 0 -1"]

/-- The complex of Fig. 5 as `tracecurve` lines. -/
def fig5Complex : List String :=
  ["faces 7",
   "f 0 6 4 5", "f 1 6 5 7", "f 2 6 7 1", "f 3 1 7 2", "f 4 3 4 6", "f 5 0 3 6", "f 6 0 6 1",
   "tw 0 13 -1 3", "tw 1 2 -1 6", "tw 2 5 9 19", "tw 3 7 -1 -1", "tw 4 -1 0 16", "tw 5 -1 14 18",
   "tw 6 17 8 -1",
   "nf 0 1 -1 3", "nf 1 3 -1 3", "nf 2 3 1 0", "nf 3 1 -1 -1", "nf 4 -1 1 0", "nf 5 -1 0 0",
   "nf 6 0 0 -1"]

/-- A Δ-complex built by hand: face 0 `(i, j, j)` whose sides `i→j` (halfedge 0) and `j→i`
(halfedge 2) are glued to each other (a cone around `i`) and whose third side is the self-edge
`j→j` (halfedge 1), glued to the self-edge of face 1 `(j, x, j)` (halfedge 5); halfedges 3, 4 are
boundary. Vertices `i = 0, j = 1, x = 2`. -/
def deltaComplex (a b : Int) (n3 n4 : Int) : List String :=
  ["faces 2", "f 0 0 1 1", "f 1 1 2 1", "tw 0 2 5 0", "tw 1 -1 -1 1",
   s!"nf 0 {a} {b} {a}", s!"nf 1 {n3} {n4} {b}"]

/-- The cases taken from the paper's figures. -/
def figureCases : List (ICCase × String) := [
  (mkCase "ic-fig5-corner" "corner" fig5Faces,
    "Fig. 5 right, the seven faces of T1 with the printed labels: F0 (m,BM,BR) cuts corner m once (diagonal BL-TM) and BM emanates 2 (TM-BM, BM-TR) to BR-m; F2 (m,r,TM): TM emanates 2, one arc cuts corner r; F3: TR emanates 1; F4: BL emanates 1; F1: three arcs cut corner m"),
  (mkCase "ic-fig5-subdivide" "subdivide" fig5Faces,
    "Fig. 5 right, the pairing of the crossings in each face (Alg. 9) for the faces above; checked by hand against the drawing"),
  (mkCase "ic-fig5-tracecurve" "tracecurve" (fig5Complex ++
      ["start 3 0", "start 3 1", "start 3 2", "start 13 0", "start 7 0", "start 2 0", "start 6 0"]),
    "Fig. 5 right, whole complex: from m-BR (halfedge 3) the crossings 0, 1, 2 (from m) end at BL, BM, BM (BL-TM, TM-BM, BM-TR traced towards their lower ends); BL-TM traced from m-BM towards TM ends at TM; BM-TR from r-TM ends at TR; from BR-m (halfedge 2) crossing 0 (BM-TR) ends at TR; from m-r crossing 0 ends at BL"),
  (mkCase "ic-fig5-flip" "flip" ["quad 3 1 0 3 -1", "again"],
    "Fig. 5 right, flip m-r (faces m,r,TM and r,m,BR): the new edge TM-BR crosses BM-TR once, n = 1; flipping back gives 3"),
  (mkCase "ic-fig5-identify" "identify" ["walk 0 -1 3 1 0 0", "walk 0 -1 3 1 1 0"],
    "Fig. 5 right, Alg. 3 at the boundary vertex BM (T0 halfedges counter-clockwise from the boundary: BM->BR 0, BM->TR 1, BM->TM 2, BM->BL 3; r of BM->BR is 0, n = -1): the curves from BM crossing BR->m at 0 and 1 (from BR) are BM->TR (index 1) and BM->TM (index 2)"),
  (mkCase "ic-fig7-trace" "trace" ["cross 7 2 4 4", "cross 7 2 3 2", "cross 2 4 1 1"],
    "Fig. 7 (= Fig. 20), face jik read off the drawing: case 1 n_ij=7, n_ik=2, n_kj=4 (c_j=4, e_k=1, c_i=2), the teal crossing p=4 goes left to kj at 1; case 2 n_ij=7, n_ik=2, n_kj=3 (e_k=2), p=2 ends at k; case 3 n_ij=2, n_ik=4, n_kj=1 (c_i=2, c_k=1, e_j=1), p=1 goes right to ik at 1"),
  (mkCase "ic-fig9-flip" "flip" ["quad 3 0 0 3 -1", "roundabout 2 8 4 6", "again",
      "quad 3 0 0 3 -1", "roundabout 6 8 5 6"],
    "Fig. 9 / inset of Eq. 4: three curves from k cross ij and end on il (e^il_k = 3), lj runs along an input edge (n_lj = -1, n-_lj = 1), n-_ki = 0, e^jk_l = 0: n_kl = 0, r_kl = r_ki + 3, r_lk = r_lj + 1 (second quad: the same with the reduction modulo deg0)"),
  (mkCase "ic-fig10-split" "split" ["face 3 1 2 corner 0 1", "face -1 -1 -1 middle 0"],
    "Fig. 10: top, face ijk (sides i->j 3, j->k 1, k->i 2; c_i=2, c_j=1) split at the star between the two arcs around i: n_mi=1, n_mj=2, n_mk=1 as printed; bottom, a face whose side ij runs along a curve: all three counts 0"),
  (mkCase "ic-fig10-flip" "flip" ["quad 3 2 1 1 4", "again", "quad -1 0 0 -1 -1", "again"],
    "Fig. 10: top, after the split flip ij (faces i,j,m with n 3,2,1 and j,i,l with 3,1,4): n_lm = 2 as printed; bottom, ij with n = -1 flipped: n_lm = 1 (the new edge crosses the curve along ij, the failure the figure shows)"),
  (mkCase "ic-fig10-splitedge" "splitedge" ["edge -1 -1", "edge 2 -1", "edge 0 3"],
    "Fig. 10 bottom / Eq. 6: splitting the shared edge explicitly, n_mk = max(n_ki, n_jk, 0): 0 when both other sides are input edges"),
  (mkCase "ic-fig11-subdivide" "subdivide" ["face 4 3 5", "face 3 8 2"],
    "Fig. 11: case 1 (corners k,i,j; c_k=3, c_i=1, c_j=2, sides k->i 4, i->j 3, j->k 5) and case 2 (c_k=0, e_k=3, c_i=3, c_j=2, sides 3, 8, 2)"),
  (mkCase "ic-fig11-regions" "regions" ["face 4 3 5", "face 3 8 2"],
    "Fig. 11: the regions of the two faces: 3+1+2 corner regions and the central one; 3+2 corner regions and the fan of 4"),
  (mkCase "ic-fig22-split" "split" ["face 4 3 5 corner 0 2", "face 2 5 1 middle 1"],
    "Fig. 22: case 1 (same face as Fig. 11 case 1), m beyond 2 of the 3 arcs around k (nu_k=2, sigma_k=1): n_mk=2, n_mi=2, n_mj=3 as printed; case 2 (c_i=2, c_j=1, e_k=2, sides k->i 2, i->j 5, j->k 1), m in the fan after one emanating arc (epsilon_k=1): n_mk=0, n_mi=3, n_mj=2 as printed"),
  (mkCase "ic-delta-tracecurve" "tracecurve" (deltaComplex 1 0 (-1) (-1) ++ ["start 0 0"] ++
      ["bound 50"]),
    "Delta-complex by hand (face (i,j,j) with its sides i->j and j->i glued, self-edge j->j): one arc cuts the cone around i, a closed curve: the trace reports the repeated crossing (SPEC I18)"),
  (mkCase "ic-delta-tracecurve-2" "tracecurve" (deltaComplex 2 2 (-1) (-1) ++
      ["start 5 0", "start 5 1", "start 1 0", "start 1 1", "start 0 0", "start 0 1"]),
    "Delta-complex by hand, counts 2 on the glued sides and on the self-edge: x emanates two curves through the self-edge; traced through the cone")
]

/-! ## Generated cases -/

/-- A deterministic generator (64-bit LCG, Knuth's constants). -/
def lcg (s : UInt64) : UInt64 := s * 6364136223846793005 + 1442695040888963407

/-- `count` values in `[lo, hi]` from the seed. -/
def lcgInts (seed : UInt64) (count : Nat) (lo hi : Int) : Array Int := Id.run do
  let mut s := lcg (seed + 1)
  let mut out : Array Int := #[]
  let span := (hi - lo + 1).toNat
  for _ in [0:count] do
    s := lcg s
    out := out.push (lo + (((s >>> 33).toNat % span : Nat) : Int))
  return out

/-- All triples with entries in `[lo, hi]`. -/
def triples (lo hi : Int) : List (Int × Int × Int) := Id.run do
  let r := (List.range (hi - lo + 1).toNat).map fun (i : Nat) => lo + (i : Int)
  let mut out := []
  for a in r do
    for b in r do
      for c in r do
        out := (a, b, c) :: out
  return out.reverse

/-- Generated cases. -/
def generatedCases : List (ICCase × String) := Id.run do
  let all := triples (-1) 6
  let faceLines := all.map fun (a, b, c) => s!"face {a} {b} {c}"
  let small := (triples 0 4).filter fun (a, b, c) => faceValid a b c
  let traceLines := small.flatMap fun (a, b, c) =>
    (List.range a.toNat).map fun p => s!"cross {a} {b} {c} {p}"
  -- random valid quads
  let xs := lcgInts 2021 6000 (-1) 8
  let mut quads : List String := []
  let mut count := 0
  for i in [0:1000] do
    let q := (xs[5 * i]!, xs[5 * i + 1]!, xs[5 * i + 2]!, xs[5 * i + 3]!, xs[5 * i + 4]!)
    let (a, b, c, d, e) := q
    if count < 300 && faceValid a b c && faceValid a d e then
      let r := lcgInts (3000 + i.toUInt64) 4 0 9
      quads := quads ++ [s!"quad {a} {b} {c} {d} {e}",
        s!"roundabout {r[0]!} {(r[1]! + 1).toNat} {r[2]!} {(r[3]! + 1).toNat}", "again"]
      count := count + 1
  -- Only arguments whose face `ija` is valid (Eq. 1 numerator even): Alg. 3 is defined on valid
  -- states only, and the Rust side refuses an odd numerator (IC SPEC I14).
  let ws := lcgInts 77 3600 (-1) 9
  let walks := ((List.range 600).filterMap fun i =>
    let rr := ws[6 * i]! + 1
    let nai := ws[6 * i + 1]!
    let nij := ws[6 * i + 2]! + 1
    let nja := ws[6 * i + 3]!
    let p := ws[6 * i + 4]! % (nij + 1)
    let deg := (ws[6 * i + 5]! + 1).toNat
    if faceValid nij nja nai then some s!"walk {rr} {nai} {nij} {nja} {p} {deg}" else none).take 200
  let edges := ((triples (-1) 5).filter fun (_, _, c) => c == 0).map fun (a, b, _) => s!"edge {a} {b}"
  return [
    (mkCase "ic-corner-all" "corner" faceLines, "every face with n in [-1, 6]: validity, numerators, Eq. 1 (odd when invalid), Eq. 2"),
    (mkCase "ic-subdivide-all" "subdivide" faceLines, "every face with n in [-1, 6]: the decoded corner data and the pairing of Alg. 9"),
    (mkCase "ic-regions-all" "regions" faceLines, "every face with n in [-1, 6]: its regions"),
    (mkCase "ic-split-all" "split" faceLines, "every face with n in [-1, 6]: the face-split counts of every region (I15)"),
    (mkCase "ic-trace-all" "trace" traceLines, "every valid face with n in [0, 4] entered through side 0 at every crossing (Alg. 1)"),
    (mkCase "ic-flip-random" "flip" quads, "300 random valid quadrilaterals (n in [-1, 8], seed 2021) with Eq. 4 and the flip back"),
    (mkCase "ic-identify-random" "identify" walks, "200 random arguments of Alg. 3's index on valid faces (seed 77), with and without the reduction"),
    (mkCase "ic-splitedge-all" "splitedge" edges, "Eq. 6 for n in [-1, 5]")
  ]

/-- All fixture cases, in order. -/
def allCases : List (ICCase × String) := figureCases ++ generatedCases

end Ddg.IntegerCoords
