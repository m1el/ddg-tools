import DdgOracle.Formulas.Sum
import DdgOracle.Formulas.Triangle

/-!
# Scalar-generic formulas of "Developability of Triangle Meshes" (`crates/ddg-developability/SPEC.md`)

Stein, Grinspun, Crane, ACM TOG 37(4) 2018. Every definition carries the id of the formula it
belongs to as `[formula: dev-…]`; the theorems about them (at `ℝ`) are in
`lean/proofs/DdgProofs/Developability/`. They are written from the SPEC (which is written from the
paper's text only; clean room: the authors' code is never consulted).

Conventions (SPEC §1). Vectors of `ℝ³` are `Ddg.Point`. A face is a CCW triple `(f_i, f_j, f_k)`;
`n = (f_j − f_i) × (f_k − f_i)` is twice its vector area, `N = n/|n|` its unit normal, the corner
angle at `f_i` is `atan2(|a × b|, a · b)` with `a = f_j − f_i`, `b = f_k − f_i`, and the edge
opposite a corner is taken CCW: `e_i = f_k − f_j`, `e_j = f_i − f_k`, `e_k = f_j − f_i`.

Index sets are given as a count `n` and a function `Nat → _` (the entries `0 … n − 1`); finite sums
use `sumTo`/`vsumTo` (associated to the left), finite maxima `maxTo1` (of the entries `0 … n`).

Everything is exact closed forms, except two numerical devices that are marked as such and are not
covered by any theorem: the Taylor series of the log-map coefficients near `φ = 0`
(`logGSeries`, `logHSeries`) and the overflow guard of the Jacobi rotation (SPEC §3.10), which the
Rust code applies for `|θ| > 1e150` and which is not reproduced here.
-/

namespace Ddg

namespace Dev

variable {α : Type} [Scalar α]

/-! ## Vectors of `ℝ³` -/

/-- `p + q`. -/
def add (p q : Point α) : Point α := ⟨p.x + q.x, p.y + q.y, p.z + q.z⟩

/-- `p − q`. -/
def sub (p q : Point α) : Point α := ⟨p.x - q.x, p.y - q.y, p.z - q.z⟩

/-- `−p`. -/
def neg (p : Point α) : Point α := ⟨-p.x, -p.y, -p.z⟩

/-- `s p`. -/
def smul (s : α) (p : Point α) : Point α := ⟨s * p.x, s * p.y, s * p.z⟩

/-- The zero vector. -/
def zero : Point α := ⟨0, 0, 0⟩

/-- Dot product. -/
def dot (p q : Point α) : α := p.x * q.x + p.y * q.y + p.z * q.z

/-- Cross product. -/
def cross (p q : Point α) : Point α :=
  ⟨p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x⟩

/-- Squared length `|p|²`. -/
def norm2 (p : Point α) : α := dot p p

/-- Length `|p|`. -/
def norm (p : Point α) : α := Scalar.sqrt (norm2 p)

/-- `p/|p|` (undefined, i.e. a division by zero, for `p = 0`). -/
def normalize (p : Point α) : Point α := smul (1 / norm p) p

/-- `det(a, b, c) = (a × b) · c`. -/
def det3 (a b c : Point α) : α := dot (cross a b) c

/-- `f 0 + f 1 + ⋯ + f (n − 1)` for vectors, componentwise `sumTo`. -/
def vsumTo (n : Nat) (f : Nat → Point α) : Point α :=
  ⟨sumTo n (fun i => (f i).x), sumTo n (fun i => (f i).y), sumTo n (fun i => (f i).z)⟩

/-- `max (f 0) … (f n)`: the maximum of the `n + 1` values `f 0, …, f n`. -/
def maxTo1 : Nat → (Nat → α) → α
  | 0, f => f 0
  | n + 1, f => max (maxTo1 n f) (f (n + 1))

/-- `min (f 0) … (f n)`. -/
def minTo1 : Nat → (Nat → α) → α
  | 0, f => f 0
  | n + 1, f => min (minTo1 n f) (f (n + 1))

/-- The two-argument arc tangent with values in `(−π, π]` (`atan2(0, 0) = 0`, as in C):
the angle of the point `(x, y)`. For `y ≥ 0` it lies in `[0, π]`. [formula: dev-corner-angle] -/
def atan2 (y x : α) : α :=
  if 0 < x then Scalar.atan (y / x)
  else if x < 0 then
    (if y < 0 then Scalar.atan (y / x) - Scalar.pi else Scalar.atan (y / x) + Scalar.pi)
  else if 0 < y then Scalar.pi / 2
  else if y < 0 then -(Scalar.pi / 2)
  else 0

/-! ## Per face (SPEC §1, §2) -/

/-- `n = (f_j − f_i) × (f_k − f_i)`, twice the vector area of the face `(f_i, f_j, f_k)`.
[formula: dev-face-frame] -/
def faceNormal (fi fj fk : Point α) : Point α := cross (sub fj fi) (sub fk fi)

/-- Area `A = |n|/2`. [formula: dev-face-frame] -/
def faceArea (fi fj fk : Point α) : α := norm (faceNormal fi fj fk) / 2

/-- Unit normal `N = n/|n|`. [formula: dev-face-frame] -/
def faceUnitNormal (fi fj fk : Point α) : Point α := normalize (faceNormal fi fj fk)

/-- Corner angle at `f_i`: `atan2(|a × b|, a · b)`, `a = f_j − f_i`, `b = f_k − f_i`.
[formula: dev-corner-angle] -/
def cornerAngleAt (fi fj fk : Point α) : α :=
  atan2 (norm (cross (sub fj fi) (sub fk fi))) (dot (sub fj fi) (sub fk fi))

/-- `∇_{f_i} A = ½ N × (f_k − f_j)` (Eq. 7). [formula: dev-area-gradient] -/
def areaGradient (fi fj fk : Point α) : Point α :=
  smul (1 / 2) (cross (faceUnitNormal fi fj fk) (sub fk fj))

/-- Derivative of the unit normal for a motion `δ` of `f_i` (Eq. 8 with `1/(2A)`, D2):
`dN = (N · δ) (e_i × N)/(2A)`, `e_i = f_k − f_j`. [formula: dev-normal-gradient] -/
def normalDerivative (fi fj fk δ : Point α) : Point α :=
  smul (dot (faceUnitNormal fi fj fk) δ / (2 * faceArea fi fj fk))
    (cross (sub fk fj) (faceUnitNormal fi fj fk))

/-- `∇_{f_j} θ_i = N × (f_i − f_j)/|f_i − f_j|²` (Eq. 9 with the square, D3).
[formula: dev-angle-gradient] -/
def angleGradientJ (fi fj fk : Point α) : Point α :=
  smul (1 / norm2 (sub fi fj)) (cross (faceUnitNormal fi fj fk) (sub fi fj))

/-- `∇_{f_k} θ_i = N × (f_k − f_i)/|f_k − f_i|²`. [formula: dev-angle-gradient] -/
def angleGradientK (fi fj fk : Point α) : Point α :=
  smul (1 / norm2 (sub fk fi)) (cross (faceUnitNormal fi fj fk) (sub fk fi))

/-- `∇_{f_i} θ_i = −(∇_{f_j} θ_i + ∇_{f_k} θ_i)`. [formula: dev-angle-gradient] -/
def angleGradientI (fi fj fk : Point α) : Point α :=
  neg (add (angleGradientJ fi fj fk) (angleGradientK fi fj fk))

/-- Sensitivity with respect to `n` (SPEC §2): `w = (y − (y·N) N)/|n| + ζ`, from the sensitivity
`y` of the unit normal and `ζ` of `n` itself (`ζ = ½ Σ_u z_u` for the vertex normals).
[formula: dev-backprop-face] -/
def backpropW (fi fj fk y zeta : Point α) : Point α :=
  let N := faceUnitNormal fi fj fk
  add (smul (1 / norm (faceNormal fi fj fk)) (sub y (smul (dot y N) N))) zeta

/-- The gradient contributions `(g_i, g_j, g_k)` of one face (SPEC §2): normal term
`w × e_p`, then the angle terms of the corners `i, j, k` with sensitivities `s_i, s_j, s_k`.
[formula: dev-backprop-face] -/
def backpropFace (fi fj fk y zeta : Point α) (si sj sk : α) : Point α × Point α × Point α :=
  let w := backpropW fi fj fk y zeta
  -- normal term
  let gi := cross w (sub fk fj)
  let gj := cross w (sub fi fk)
  let gk := cross w (sub fj fi)
  -- corner i: angle at f_i of (f_i, f_j, f_k)
  let gi := add gi (smul si (angleGradientI fi fj fk))
  let gj := add gj (smul si (angleGradientJ fi fj fk))
  let gk := add gk (smul si (angleGradientK fi fj fk))
  -- corner j: angle at f_j of (f_j, f_k, f_i)
  let gj := add gj (smul sj (angleGradientI fj fk fi))
  let gk := add gk (smul sj (angleGradientJ fj fk fi))
  let gi := add gi (smul sj (angleGradientK fj fk fi))
  -- corner k: angle at f_k of (f_k, f_i, f_j)
  let gk := add gk (smul sk (angleGradientI fk fi fj))
  let gi := add gi (smul sk (angleGradientJ fk fi fj))
  let gj := add gj (smul sk (angleGradientK fk fi fj))
  (gi, gj, gk)

/-! ## Vertex normal (SPEC §1) -/

/-- Vector area of a star given by its faces: `S = ½ Σ n_σ`. [formula: dev-vertex-normal] -/
def vectorArea (n : Nat) (face : Nat → Point α × Point α × Point α) : Point α :=
  smul (1 / 2) (vsumTo n (fun m => faceNormal (face m).1 (face m).2.1 (face m).2.2))

/-- Vector area of the star of an interior vertex `v` with the cyclic ring `w_0 … w_{k−1}`
(faces `(v, w_m, w_{m+1 mod k})`). [formula: dev-vertex-normal] -/
def ringVectorArea (v : Point α) (k : Nat) (ring : Nat → Point α) : Point α :=
  vectorArea k (fun m => (v, ring m, ring ((m + 1) % k)))

/-- Vertex normal `N_v = S_v/|S_v|`. [formula: dev-vertex-normal] -/
def vertexNormal (n : Nat) (face : Nat → Point α × Point α × Point α) : Point α :=
  normalize (vectorArea n face)

/-! ## Combinatorial energies (SPEC §3.2, §3.3) -/

/-- The edge-connected bipartitions of an interior star of valence `k`, as pairs `(a, b)`,
`0 ≤ a < b ≤ k − 1`, in lexicographic order. [formula: dev-star-bipartitions] -/
def interiorBipartitions (k : Nat) : List (Nat × Nat) :=
  (List.range k).flatMap fun a => ((List.range k).filter (fun b => a < b)).map fun b => (a, b)

/-- The part `F_1 = {σ_{a+1}, …, σ_b}` of the bipartition `(a, b)`. [formula: dev-star-bipartitions] -/
def arcPart (a b : Nat) : List Nat := (List.range (b - a)).map fun i => a + 1 + i

/-- The other part `F_2 = {σ_{b+1}, …, σ_a}` (indices mod `k`). [formula: dev-star-bipartitions] -/
def arcComplement (k a b : Nat) : List Nat :=
  (List.range (k - (b - a))).map fun i => (b + 1 + i) % k

/-- The bipartitions of a boundary star with `k` faces: split after face `a`, `0 ≤ a ≤ k − 2`,
`F_1 = {σ_0 … σ_a}`. [formula: dev-star-bipartitions] -/
def boundarySplits (k : Nat) : List Nat := List.range (k - 1)

/-- Mean `N̄ = (1/n) Σ N_σ` of a part with `n` faces. [formula: dev-combinatorial-width] -/
def partMean (n : Nat) (N : Nat → Point α) : Point α := smul (1 / Scalar.ofNat n) (vsumTo n N)

/-- Width of a part, first form of Eq. 1: `(1/n) Σ |N_σ − N̄|²` (the form the code evaluates).
[formula: dev-combinatorial-width] -/
def partWidth (n : Nat) (N : Nat → Point α) : α :=
  sumTo n (fun i => norm2 (sub (N i) (partMean n N))) / Scalar.ofNat n

/-- Second form of Eq. 1, over **unordered** pairs (D1): `(1/n²) Σ_{σ1<σ2} |N_σ1 − N_σ2|²`.
[formula: dev-combinatorial-width] -/
def partWidthPairs (n : Nat) (N : Nat → Point α) : α :=
  sumTo n (fun i => sumTo i (fun j => norm2 (sub (N i) (N j)))) / (Scalar.ofNat n * Scalar.ofNat n)

/-- Third form of Eq. 1 for unit normals: `1 − |N̄|²`. [formula: dev-combinatorial-width] -/
def partWidthMean (n : Nat) (N : Nat → Point α) : α := 1 - norm2 (partMean n N)

/-- Sensitivity of the width with respect to `N_σ`: `(2/n)(N_σ − N̄)`.
[formula: dev-combinatorial-gradient] -/
def partWidthSensitivity (n : Nat) (N : Nat → Point α) (i : Nat) : Point α :=
  smul (2 / Scalar.ofNat n) (sub (N i) (partMean n N))

/-- Width of a part for `π_max` (Eq. 5): `max_{σ1,σ2} |N_σ1 − N_σ2|²` over the `n + 1` normals
`N 0 … N n` (a part with one face has width `0`). [formula: dev-combinatorial-max] -/
def partMaxWidth (n : Nat) (N : Nat → Point α) : α :=
  maxTo1 n (fun i => maxTo1 i (fun j => norm2 (sub (N i) (N j))))

/-- Sensitivities of the active pair `(σ1, σ2)` of `π_max`: `y_{σ1} = 2(N_σ1 − N_σ2)`,
`y_{σ2} = −y_{σ1}`. [formula: dev-combinatorial-max] -/
def pairSensitivity (N1 N2 : Point α) : Point α := smul 2 (sub N1 N2)

/-! ## Covariance energy (SPEC §3.4) -/

/-- A symmetric 3x3 matrix by its six entries. -/
structure Sym3 (α : Type) where
  /-- entry (0,0) -/
  xx : α
  /-- entry (1,1) -/
  yy : α
  /-- entry (2,2) -/
  zz : α
  /-- entry (0,1) = (1,0) -/
  xy : α
  /-- entry (0,2) = (2,0) -/
  xz : α
  /-- entry (1,2) = (2,1) -/
  yz : α

namespace Sym3

/-- `A v`. -/
def apply (A : Sym3 α) (v : Point α) : Point α :=
  ⟨A.xx * v.x + A.xy * v.y + A.xz * v.z,
   A.xy * v.x + A.yy * v.y + A.yz * v.z,
   A.xz * v.x + A.yz * v.y + A.zz * v.z⟩

/-- The quadratic form `uᵀ A u`. -/
def quad (A : Sym3 α) (u : Point α) : α := dot u (apply A u)

/-- `A + B`. -/
def add (A B : Sym3 α) : Sym3 α :=
  ⟨A.xx + B.xx, A.yy + B.yy, A.zz + B.zz, A.xy + B.xy, A.xz + B.xz, A.yz + B.yz⟩

/-- The zero matrix. -/
def zero : Sym3 α := ⟨0, 0, 0, 0, 0, 0⟩

/-- `θ N Nᵀ`. -/
def outer (θ : α) (N : Point α) : Sym3 α :=
  ⟨θ * (N.x * N.x), θ * (N.y * N.y), θ * (N.z * N.z), θ * (N.x * N.y), θ * (N.x * N.z),
   θ * (N.y * N.z)⟩

/-- `Σ_{m < n} A m`. -/
def sumTo : Nat → (Nat → Sym3 α) → Sym3 α
  | 0, _ => zero
  | n + 1, A => add (sumTo n A) (A n)

end Sym3

/-- `A_v = Σ_σ θ_σ N_σ N_σᵀ` (Eq. 4) over the `n` corners of a star.
[formula: dev-covariance-matrix] -/
def covMatrix (n : Nat) (θ : Nat → α) (N : Nat → Point α) : Sym3 α :=
  Sym3.sumTo n (fun m => Sym3.outer (θ m) (N m))

/-- The objective of Eq. 3: `Σ_σ θ_σ ⟨u, N_σ⟩²`. [formula: dev-covariance-matrix] -/
def covObjective (n : Nat) (θ : Nat → α) (N : Nat → Point α) (u : Point α) : α :=
  sumTo n (fun m => θ m * (dot u (N m) * dot u (N m)))

/-- Sensitivity of `λ` with respect to `N_σ` (Eq. 10): `y = 2 θ (x·N) x`.
[formula: dev-covariance-gradient] -/
def covSensitivityN (θ : α) (x N : Point α) : Point α := smul (2 * θ * dot x N) x

/-- Sensitivity of `λ` with respect to `θ_σ` (Eq. 10): `s = (x·N)²`.
[formula: dev-covariance-gradient] -/
def covSensitivityTheta (x N : Point α) : α := dot x N * dot x N

/-! ## Max covariance energy (SPEC §3.5) -/

/-- `φ(u) = max_σ ⟨u, N_σ⟩²` over the `n + 1` normals `N 0 … N n`.
[formula: dev-max-covariance-voronoi] -/
def maxCovObjective (n : Nat) (N : Nat → Point α) (u : Point α) : α :=
  maxTo1 n (fun m => dot u (N m) * dot u (N m))

/-- `n = (b − a) × (c − a)` for the sites `a, b, c`. [formula: dev-spherical-circumcenter] -/
def siteNormal (a b c : Point α) : Point α := cross (sub b a) (sub c a)

/-- The spherical circumcenter `v = n/|n|` of three unit vectors (D5).
[formula: dev-spherical-circumcenter] -/
def circumcenter (a b c : Point α) : Point α := normalize (siteNormal a b c)

/-- `g = det(a, b, c)/|n| = ⟨v, a⟩ = ⟨v, b⟩ = ⟨v, c⟩`. [formula: dev-spherical-circumcenter] -/
def circumG (a b c : Point α) : α := det3 a b c / norm (siteNormal a b c)

/-- `∂(g²)/∂q` for the site `q` with opposite edge `e_q` (SPEC §3.5):
`2g [[q = a] v + ((e_q × v)·a/|n|) v]`; `isA` says whether `q = a`.
[formula: dev-max-covariance-gradient] -/
def circumSensitivity (a b c eq : Point α) (isA : Bool) : Point α :=
  let v := circumcenter a b c
  let g := dot v a
  let coef := dot (cross eq v) a / norm (siteNormal a b c)
  smul (2 * g) (add (if isA then v else zero) (smul coef v))

/-- The three sensitivities `(∂/∂a, ∂/∂b, ∂/∂c)` of `g(a, b, c)²`, with `e_a = c − b`,
`e_b = a − c`, `e_c = b − a`. A site `q = s_q N_σ` passes `s_q ∂/∂q` to `N_σ` (D5).
[formula: dev-max-covariance-gradient] -/
def maxCovSensitivities (a b c : Point α) : Point α × Point α × Point α :=
  (circumSensitivity a b c (sub c b) true, circumSensitivity a b c (sub a c) false,
   circumSensitivity a b c (sub b a) false)

/-! ## Intrinsic width (SPEC §3.6) -/

/-- `φ = atan2(|n0 × N|, ⟨n0, N⟩)`, the angle from `n0` to `N`. [formula: dev-log-map] -/
def logAngle (n0 N : Point α) : α := atan2 (norm (cross n0 N)) (dot n0 N)

/-- `G = φ / sin φ`, with `sin φ = |n0 × N|` (exact form; undefined at `φ = 0`).
[formula: dev-log-map] -/
def logG (n0 N : Point α) : α := logAngle n0 N / norm (cross n0 N)

/-- `H = (sin φ − φ cos φ)/sin³ φ` (exact form). [formula: dev-log-map-gradient] -/
def logH (n0 N : Point α) : α :=
  let s := norm (cross n0 N)
  (s - logAngle n0 N * dot n0 N) / (s * s * s)

/-- Taylor polynomial `1 + φ²/6 + 7φ⁴/360` of `G` (numerical device for `φ < 1e-3`; not covered by
a theorem). -/
def logGSeries (φ : α) : α := 1 + φ * φ / 6 + 7 * (φ * φ * φ * φ) / 360

/-- Taylor polynomial `1/3 + 2φ²/15` of `H` (numerical device; not covered by a theorem). -/
def logHSeries (φ : α) : α := 1 / 3 + 2 * (φ * φ) / 15

/-- `Ñ = log_{n0}(N) = (φ/sin φ)(N − ⟨n0, N⟩ n0)`. [formula: dev-log-map] -/
def logMap (n0 N : Point α) : Point α := smul (logG n0 N) (sub N (smul (dot n0 N) n0))

/-- The exponential map of the unit sphere at `p`: `exp_p(v) = cos|v| p + (sin|v|/|v|) v`, written
with given `cos|v|` and `sin|v|/|v|` (the scalar interface has no `sin`, `cos`).
[formula: dev-log-map] -/
def expMapWith (cosv sincv : α) (p v : Point α) : Point α := add (smul cosv p) (smul sincv v)

/-- Adjoint of the log map for a sensitivity `z` of `Ñ` (SPEC §3.6, corrected Eq. 11, D6):
`K = G (z·n0) + H z·(N − c n0)`, `y_N = G z − K n0`, `y_{n0} = −G c z − K N`.
Returns `(y_N, y_{n0})`. [formula: dev-log-map-gradient] -/
def logMapAdjointWith (G H : α) (n0 N z : Point α) : Point α × Point α :=
  let c := dot n0 N
  let K := G * dot z n0 + H * dot z (sub N (smul c n0))
  (sub (smul G z) (smul K n0), sub (neg (smul (G * c) z)) (smul K N))

/-- The adjoint with the exact `G`, `H`. [formula: dev-log-map-gradient] -/
def logMapAdjoint (n0 N z : Point α) : Point α × Point α :=
  logMapAdjointWith (logG n0 N) (logH n0 N) n0 N z

/-- Sensitivity of a vertex normal `n0 = S/|S|` passed to the vector area:
`z_v = (y − (y·n0) n0)/|S|`. [formula: dev-log-map-gradient] -/
def vertexNormalSensitivity (S y : Point α) : Point α :=
  let n0 := normalize S
  smul (1 / norm S) (sub y (smul (dot y n0) n0))

/-! ## Branching-free intrinsic width (SPEC §3.7) -/

/-- `max_σ ⟨P_σ, u⟩ − min_σ ⟨P_σ, u⟩` over the `n + 1` points `P 0 … P n`.
[formula: dev-intrinsic-width] -/
def widthAlong (n : Nat) (P : Nat → Point α) (u : Point α) : α :=
  maxTo1 n (fun m => dot (P m) u) - minTo1 n (fun m => dot (P m) u)

/-- `ψ(u) = max_{σ,τ} ⟨P_σ − P_τ, u⟩²` over ordered pairs of the `n + 1` points.
[formula: dev-intrinsic-width] -/
def pairWidthSq (n : Nat) (P : Nat → Point α) (u : Point α) : α :=
  maxTo1 n (fun i => maxTo1 n (fun j => dot (sub (P i) (P j)) u * dot (sub (P i) (P j)) u))

/-- The unit altitude direction `u = n0 × ŵ`, `w = x_b − x_a`. [formula: dev-width-altitude-gradient] -/
def altitudeDir (n0 xa xb : Point α) : Point α := cross n0 (normalize (sub xb xa))

/-- Derivative of the altitude direction in the plane: `du = −(u·dw) ŵ/|w|`.
[formula: dev-width-altitude-gradient] -/
def altitudeDerivative (n0 xa xb dw : Point α) : Point α :=
  let w := sub xb xa
  smul (-(dot (altitudeDir n0 xa xb) dw) / norm w) (normalize w)

/-- Sensitivities `(z_s, z_t, z_a, z_b)` of `Δ²`, `Δ = (x_s − x_t)·u`, `u = n0 × ŵ`,
`w = x_b − x_a` (SPEC §3.7): `z_s = 2Δu`, `z_t = −2Δu`, `z_a = 2Δ (d·ŵ)/|w| u`, `z_b = −z_a`.
[formula: dev-min-width-gradient] -/
def minWidthSensitivities (n0 xs xt xa xb : Point α) : Point α × Point α × Point α × Point α :=
  let w := sub xb xa
  let u := altitudeDir n0 xa xb
  let d := sub xs xt
  let Δ := dot d u
  let k := 2 * Δ * dot d (normalize w) / norm w
  (smul (2 * Δ) u, smul (-(2 * Δ)) u, smul k u, smul (-k) u)

/-! ## Jacobi rotation (SPEC §3.10) -/

/-- Sign `±1` (`+1` at `0`). -/
def sgn (x : α) : α := if 0 ≤ x then 1 else -1

/-- Rutishauser's `t = tan φ`: the root of smaller modulus of `t² + 2θt − 1 = 0`,
`θ = (a_qq − a_pp)/(2 a_pq)`, `t = sgn θ/(|θ| + √(θ² + 1))`. [formula: dev-jacobi-rotation] -/
def jacobiT (app aqq apq : α) : α :=
  let θ := (aqq - app) / (2 * apq)
  sgn θ / (Scalar.abs θ + Scalar.sqrt (θ * θ + 1))

/-- `c = 1/√(t² + 1)` (and `s = t c`). [formula: dev-jacobi-rotation] -/
def jacobiC (t : α) : α := 1 / Scalar.sqrt (t * t + 1)

/-- The entries of a symmetric 3x3 matrix in the order of a rotation plane `(p, q)` and the third
index `r`. -/
structure PlaneEntries (α : Type) where
  /-- `a_pp` -/
  app : α
  /-- `a_qq` -/
  aqq : α
  /-- `a_rr` -/
  arr : α
  /-- `a_pq` -/
  apq : α
  /-- `a_rp` -/
  arp : α
  /-- `a_rq` -/
  arq : α

/-- One Jacobi rotation in the plane `(p, q)` as `eig::sym_eig3` applies it:
`a_pp ← a_pp − t a_pq`, `a_qq ← a_qq + t a_pq`, `a_pq ← 0`, `a_rp ← c a_rp − s a_rq`,
`a_rq ← s a_rp + c a_rq`. [formula: dev-jacobi-rotation] -/
def jacobiStep (m : PlaneEntries α) : PlaneEntries α :=
  let t := jacobiT m.app m.aqq m.apq
  let c := jacobiC t
  let s := t * c
  ⟨m.app - t * m.apq, m.aqq + t * m.apq, m.arr, 0, c * m.arp - s * m.arq, s * m.arp + c * m.arq⟩

/-! ## Line search (SPEC §4.2, §4.3) -/

/-- Armijo condition `S(t)`: `h(t) < c1 s t` (strict). [formula: dev-lewis-overton-line-search] -/
def armijoHolds (c1 s t h : α) : Bool := decide (h < c1 * s * t)

/-- Weak Wolfe condition `C(t)`: `g(x + t d)ᵀ d > c2 s`. [formula: dev-lewis-overton-line-search] -/
def curvatureHolds (c2 s gd : α) : Bool := decide (c2 * s < gd)

/-- Outcome of the weak Wolfe line search. -/
inductive LineSearchOutcome (α : Type) where
  /-- `t` satisfies `S` and `C` -/
  | accepted (t : α)
  /-- the step cap was reached; `t = α ≥ t_max` satisfies `S` -/
  | capped (t : α)
  /-- a limit was exceeded; `t = α > 0` satisfies `S` -/
  | weak (t : α)
  /-- a limit was exceeded with `α = 0` -/
  | failed

/-- The loop of Lewis–Overton Algorithm 2.6 (SPEC §4.3). `eval t` is `none` when the energy is
undefined at `x + t d` and `some (h(t), g(x + t d)ᵀ d)` otherwise. State: the bracket `lo`
(`α`), `hi` (`β`, `none` = `+∞`), the trial `t`, the numbers of doublings and bisections done;
`fuel` bounds the number of trials. Returns the outcome and the trials in order.
[formula: dev-lewis-overton-line-search] -/
def lineSearchLoop (eval : α → Option (α × α)) (c1 c2 s tmax : α) (maxDoublings maxBisections : Nat) :
    Nat → α → Option α → α → Nat → Nat → List α → LineSearchOutcome α × List α
  | 0, lo, _, _, _, _, trials =>
    (if 0 < lo then .weak lo else .failed, trials)
  | fuel + 1, lo, hi, t, nd, nb, trials =>
    let trials := trials ++ [t]
    let (sOk, cOk) := match eval t with
      | none => (false, false)
      | some (h, gd) => (armijoHolds c1 s t h, curvatureHolds c2 s gd)
    if sOk && cOk then (.accepted t, trials)
    else
      let lo' := if sOk then t else lo
      let hi' := if sOk then hi else some t
      match hi' with
      | some b =>
        if nb < maxBisections then
          lineSearchLoop eval c1 c2 s tmax maxDoublings maxBisections fuel lo' hi' ((lo' + b) / 2) nd (nb + 1) trials
        else (if 0 < lo' then .weak lo' else .failed, trials)
      | none =>
        if tmax ≤ lo' then (.capped lo', trials)
        else if nd < maxDoublings then
          lineSearchLoop eval c1 c2 s tmax maxDoublings maxBisections fuel lo' none (min (2 * lo') tmax) (nd + 1) nb trials
        else (if 0 < lo' then .weak lo' else .failed, trials)

/-- Lewis–Overton weak Wolfe line search with `α = 0`, `β = +∞`, `t = min(1, t_max)`, at most
`maxDoublings` doublings and `maxBisections` bisections. [formula: dev-lewis-overton-line-search] -/
def lineSearch (eval : α → Option (α × α)) (c1 c2 s tmax : α) (maxDoublings maxBisections : Nat) :
    LineSearchOutcome α × List α :=
  lineSearchLoop eval c1 c2 s tmax maxDoublings maxBisections (maxDoublings + maxBisections + 1)
    0 none (min 1 tmax) 0 0 []

/-! ## Subdivision and remeshing (SPEC §5.2, §5.3) -/

/-- Midpoint `(a + b)/2` of an edge (4-1 subdivision, no smoothing). [formula: dev-midpoint-subdivision] -/
def midpoint (a b : Point α) : Point α := smul (1 / 2) (add a b)

/-- Number of corner angles below `θ_s`. [formula: dev-small-angle-classification] -/
def smallAngleCount (θs a0 a1 a2 : α) : Nat :=
  (if a0 < θs then 1 else 0) + (if a1 < θs then 1 else 0) + (if a2 < θs then 1 else 0)

end Dev

end Ddg
