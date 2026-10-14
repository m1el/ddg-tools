import DdgOracle.Formulas.Cuts

/-!
# Scalar-generic closed forms of intrinsic triangulations (`crates/ddg-intrinsic/SPEC.md`)

"Navigating Intrinsic Triangulations", Sharp, Soliman, Crane 2019. Every definition carries the id of
the formula it belongs to as `[formula: it-…]`. They are plain closed forms; the executable
(`DdgOracle.Intrinsic.*`) uses them for every number it prints.

Conventions of the diamond. An edge `ij` of length `a` has the triangles `ijk` and `jil`; the sides are
`b = ℓ_ik`, `c = ℓ_jk` of the first and `d = ℓ_il`, `e = ℓ_jl` of the second. Unfolded into the plane
with `i = (0, 0)` and `j = (a, 0)`, the apex `k` is at `(x_k, h_k)` above and `l` at `(x_l, −h_l)`
below the edge.
-/

namespace Ddg

variable {α : Type} [Scalar α]

/-! ## Angle normalization and signpost update -/

/-- Scale that makes the angles around a vertex of total angle `Θ` add up to `2π` (interior) or `π`
(boundary). [formula: it-angle-normalization] -/
def angleScale (isBoundary : Bool) (total : α) : α :=
  if isBoundary then Scalar.pi / total else 2 * Scalar.pi / total

/-- One subtraction of a full turn: for `x ∈ [0, 4π)` the representative in `[0, 2π)`.
[formula: it-signpost-update] -/
def reduceFullTurn (x : α) : α := if 2 * Scalar.pi ≤ x then x - 2 * Scalar.pi else x

/-- Signpost of the next halfedge around a vertex: `φ_ik = φ_ij + s θ_i^{jk}`, `s` the scale of the
vertex, reduced to `[0, 2π)` at an interior vertex (not reduced at a boundary vertex, where the angles
stay in `[0, π]`). [formula: it-signpost-update] -/
def signpostUpdate (isBoundary : Bool) (phi theta total : α) : α :=
  let v := phi + angleScale isBoundary total * theta
  if isBoundary then v else reduceFullTurn v

/-! ## The diamond: unfolding and flipping -/

/-- Abscissa of the apex of the triangle with base `a` (from the left end of the base) and the other
sides `b` (to the left end) and `c` (to the right end): `x = (a² + (b − c)(b + c))/2a`.
[formula: it-unfold-diamond] -/
def diamondX (a b c : α) : α := (a * a + (b - c) * (b + c)) / (2 * a)

/-- Height of the apex over the base: `h = 2A/a`. [formula: it-unfold-diamond] -/
def diamondH (a b c : α) : α := 2 * triArea a b c / a

/-- Length of the flipped edge: `ℓ_kl² = (x_k − x_l)² + (h_k + h_l)²`. [formula: it-flip-length] -/
def flippedLength (xk hk xl hl : α) : α :=
  Scalar.sqrt ((xk - xl) * (xk - xl) + (hk + hl) * (hk + hl))

/-- The same from the five lengths `a, b, c, d, e` of the diamond. [formula: it-flip-length] -/
def flippedLengthFromLengths (a b c d e : α) : α :=
  flippedLength (diamondX a b c) (diamondH a b c) (diamondX a d e) (diamondH a d e)

/-- Twice the signed area of the triangle `i l k` that a flip creates: `x_k h_l + x_l h_k`.
[formula: it-flip-convexity] -/
def flipArea1 (xk hk xl hl : α) : α := xk * hl + xl * hk

/-- Twice the signed area of the triangle `j k l` that a flip creates:
`(a − x_k) h_l + (a − x_l) h_k`. [formula: it-flip-convexity] -/
def flipArea2 (a xk hk xl hl : α) : α := (a - xk) * hl + (a - xl) * hk

/-- The two triangles form a strictly convex quadrilateral: both new triangles have positive area.
[formula: it-flip-convexity] -/
def flipConvex (a xk hk xl hl : α) : Bool :=
  decide (0 < flipArea1 xk hk xl hl) && decide (0 < flipArea2 a xk hk xl hl)

/-! ## Delaunay condition and cotan weight -/

/-- Excess of the opposite angles over `π`: `θ_k + θ_l − π`; the edge is locally Delaunay iff it is
`≤ 0`. [formula: it-delaunay-condition] -/
def delaunayExcess (thetaK thetaL : α) : α := thetaK + thetaL - Scalar.pi

/-- The Delaunay condition with a tolerance: `θ_k + θ_l − π ≤ tol`. [formula: it-delaunay-condition] -/
def delaunayHolds (tol thetaK thetaL : α) : Bool := decide (delaunayExcess thetaK thetaL ≤ tol)

/-- The sum `cot θ_k + cot θ_l` of the cotangents of the two angles opposite an edge. With the
workspace's positive semidefinite cotan Laplacian this is *minus twice the off-diagonal entry*
`L_ij = −½ (cot θ_k + cot θ_l)`, and it is *twice the cotan weight* `w = ½ (cot θ_k + cot θ_l)`
that `crates/ddg-intrinsic/SPEC.md` (N6, N17) defines; the weight itself is `cotanWeight cotK cotL / 2`.
It has the sign of `π − (θ_k + θ_l)`. [formula: it-cotan-weight] -/
def cotanWeight (cotK cotL : α) : α := cotK + cotL

/-- The cotan weight from the five lengths of the diamond. [formula: it-cotan-weight] -/
def cotanWeightFromLengths (a b c d e : α) : α :=
  cotanWeight (cornerCot a b c) (cornerCot a d e)

/-- The opposite angles from the five lengths: `(θ_k, θ_l)`. [formula: it-delaunay-condition] -/
def oppositeAngles (a b c d e : α) : α × α := (cornerAngle a b c, cornerAngle a d e)

/-! ## Unfolding a triangle -/

/-- Third vertex of the triangle `a b c` laid out to the left of the directed edge `a → b`, `d = b − a`:
`c = a + α d + β d⊥`, `α = (ℓ_ab² + ℓ_ca² − ℓ_bc²)/2ℓ_ab²`, `β = 2A/ℓ_ab²`, `d⊥ = (−d_y, d_x)`.
[formula: it-unfold-apex] -/
def unfoldApex (ax ay dx dy lab lbc lca : α) : α × α :=
  let alpha := (lab * lab + lca * lca - lbc * lbc) / (2 * (lab * lab))
  let beta := 2 * triArea lab lbc lca / (lab * lab)
  (ax + alpha * dx - beta * dy, ay + alpha * dy + beta * dx)

/-! ## Tracing: ray against the edges of a triangle -/

/-- Parameter at which the segment from `(x_t, y_t)` to `(x_h, y_h)` crosses the x-axis (the ray of a
tracing step): `t = y_tail/(y_tail − y_head)`. [formula: it-ray-edge-intersection] -/
def axisCrossingT (yTail yHead : α) : α := yTail / (yTail - yHead)

/-- Abscissa of that crossing: `x = x_tail + t (x_head − x_tail)`. [formula: it-ray-edge-intersection] -/
def axisCrossingX (xTail xHead t : α) : α := xTail + t * (xHead - xTail)

/-- The path leaves through the edge from the new apex to the vertex on the other side of the axis:
`true` iff the apex `y > 0` (then the exit edge ends at the vertex with `y < 0`).
[formula: it-unfold-apex] -/
def apexAboveAxis (yApex : α) : Bool := decide (0 < yApex)

/-! ## Barycentric coordinates and distances -/

/-- Twice the signed area of the triangle `a b c`. [formula: it-barycentric-coordinates] -/
def signedArea2 (ax ay bx bY cx cy : α) : α := (bx - ax) * (cy - ay) - (cx - ax) * (bY - ay)

/-- Barycentric coordinates of `p` in the triangle `a b c` of a layout: ratios of signed areas.
[formula: it-barycentric-coordinates] -/
def barycentric (ax ay bx bY cx cy px py : α) : α × α × α :=
  let area := signedArea2 ax ay bx bY cx cy
  (signedArea2 px py bx bY cx cy / area, signedArea2 ax ay px py cx cy / area,
   signedArea2 ax ay bx bY px py / area)

/-- Position of a point from its distances `r_a, r_b, r_c` to the three corners of the triangle laid
out with `a = (0,0)`, `b = (ℓ_ab, 0)` and `c` above: `x = (ℓ_ab² + r_a² − r_b²)/2ℓ_ab`,
`y = (c_x² + c_y² + r_a² − r_c² − 2 c_x x)/(2 c_y)`. [formula: it-barycentric-coordinates] -/
def pointFromDistances (lab lbc lca ra rb rc : α) : α × α :=
  let cx := diamondX lab lca lbc
  let cy := diamondH lab lca lbc
  let x := (lab * lab + ra * ra - rb * rb) / (2 * lab)
  (x, (cx * cx + cy * cy + ra * ra - rc * rc - 2 * cx * x) / (2 * cy))

/-- Euclidean distance of two points of a layout. [formula: it-point-distances] -/
def planeDistance (px py qx qy : α) : α :=
  Scalar.sqrt ((px - qx) * (px - qx) + (py - qy) * (py - qy))

/-- Distances `d_k = |Σ_j b_j v_j − v_k|` of the point with barycentric coordinates `(b_0, b_1, b_2)` to
the corners of the triangle `(ℓ_01, ℓ_12, ℓ_20)` laid out with `v_0 = (0,0)`, `v_1 = (ℓ_01, 0)`, `v_2`
above. [formula: it-point-distances] -/
def cornerDistances (l01 l12 l20 b0 b1 b2 : α) : α × α × α :=
  let cx := diamondX l01 l20 l12
  let cy := diamondH l01 l20 l12
  let px := b1 * l01 + b2 * cx
  let py := b2 * cy
  (planeDistance px py 0 0, planeDistance px py l01 0, planeDistance px py cx cy)

/-- Lengths of the two halves of an edge `ij` split at the point with barycentric coordinates
`(b_i, b_j)`: `ℓ_ip = ℓ_ij b_j`, `ℓ_jp = ℓ_ij b_i`; the lengths to the opposite corners are
`cornerDistances` with `b_k = 0`. [formula: it-point-distances] -/
def edgeSplitLengths (lij bi bj : α) : α × α := (lij * bj, lij * bi)

/-! ## Circumcenter -/

/-- Unnormalized barycentric weight of the circumcenter at the corner opposite the side `a_k`, the
other sides being `a_{k+1}, a_{k+2}`: `a_k² (a_{k+1}² + a_{k+2}² − a_k²)`. [formula: it-circumcenter] -/
def circumWeight (ak ak1 ak2 : α) : α := ak * ak * (ak1 * ak1 + ak2 * ak2 - ak * ak)

/-- Barycentric coordinates of the circumcenter of the triangle whose sides opposite the corners
`0, 1, 2` are `a_0, a_1, a_2`. [formula: it-circumcenter] -/
def circumcenterBary (a0 a1 a2 : α) : α × α × α :=
  let w0 := circumWeight a0 a1 a2
  let w1 := circumWeight a1 a2 a0
  let w2 := circumWeight a2 a0 a1
  let s := w0 + w1 + w2
  (w0 / s, w1 / s, w2 / s)

end Ddg
