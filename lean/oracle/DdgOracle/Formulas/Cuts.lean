import DdgOracle.Formulas.Mesh

/-!
# Scalar-generic closed forms of the cut flow (`crates/ddg-variational-cuts/SPEC.md`)

Variational Surface Cutting, Sharp & Crane 2018. Every definition carries the id of the formula
of SPEC §10 that it belongs to as `[formula: <id>]`. Linear-algebra statements (the Yamabe system,
the harmonic extension, the length step as a matrix equation) are not definitions here: the
executable assembles the matrices from the entries defined below and solves them densely.

The Hopf–Lax update (`cut-redistance`) needs one square root and the order; both are in
`Ddg.Scalar`, so it is written generically like everything else.
-/

namespace Ddg

variable {α : Type} [Scalar α]

/-! ## Face curvature, `cut-face-curvature` (SPEC §3.1) -/

/-- Corner angle rescaled so that the angles around an interior vertex sum to `2π`:
`θ̃ = θ · 2π / Θ`. [formula: cut-face-curvature] -/
def normalizedAngle (theta total : α) : α := theta * (2 * Scalar.pi) / total

/-- `θ̃` at a boundary vertex of `M` is `θ` itself, at an interior vertex `normalizedAngle`.
[formula: cut-face-curvature] -/
def cornerNormalized (isBoundary : Bool) (theta total : α) : α :=
  if isBoundary then theta else normalizedAngle theta total

/-- Integrated Gaussian curvature of a face from its three normalized angles:
`Ω_f = θ̃_i + θ̃_j + θ̃_k − π`. [formula: cut-face-curvature] -/
def faceCurvature (t1 t2 t3 : α) : α := t1 + t2 + t3 - Scalar.pi

/-- Curvature density `k_f = Ω_f / A_f`. [formula: cut-face-curvature] -/
def faceCurvatureDensity (omega area : α) : α := omega / area

/-- `Ω_f` of the face `f` of a mesh with edge lengths `len`; `angleSum v = Θ_v` and `isB v` says
whether `v ∈ ∂M`. The sum over all faces is `Σ_{v interior} (2π − Θ_v)`
[formula: cut-face-curvature-sum]. -/
def Face.omega (len : Nat → Nat → α) (angleSum : Nat → α) (isB : Nat → Bool) (f : Face) : α :=
  faceCurvature (cornerNormalized (isB f.a) (f.angleA len) (angleSum f.a))
    (cornerNormalized (isB f.b) (f.angleB len) (angleSum f.b))
    (cornerNormalized (isB f.c) (f.angleC len) (angleSum f.c))

/-! ## Crossings, `cut-edge-crossing` (SPEC §3.2) -/

/-- Parameter `t = φ_a / (φ_a − φ_b)` of the zero of the linear interpolant on `[a, b]`.
[formula: cut-edge-crossing] -/
def crossingParam (phiA phiB : α) : α := phiA / (phiA - phiB)

/-- `t` moved into `[δ, 1 − δ]`. [formula: cut-edge-crossing] -/
def clampCrossing (delta t : α) : α := max delta (min (1 - delta) t)

/-- `t` was moved by `clampCrossing`. [formula: cut-edge-crossing] -/
def crossingMoved (delta t : α) : Bool := decide (t < delta) || decide (1 - delta < t)

/-! ## Lengths and areas inside a face (SPEC §3.3) -/

/-- Squared distance of two points of the face with sides `lab lbc lca` whose barycentric
difference is `(da, db, dc)` (`da + db + dc = 0`):
`−(da·db·ℓab² + db·dc·ℓbc² + dc·da·ℓca²)`. [formula: cut-barycentric-length] -/
def baryLength2 (lab lbc lca da db dc : α) : α :=
  -(da * db * (lab * lab) + db * dc * (lbc * lbc) + dc * da * (lca * lca))

/-- The distance; the argument of the root is floored at `0` against rounding.
[formula: cut-barycentric-length] -/
def baryLength (lab lbc lca da db dc : α) : α :=
  Scalar.sqrt (max 0 (baryLength2 lab lbc lca da db dc))

/-- Determinant of three barycentric coordinate rows. [formula: cut-subtriangle-area] -/
def det3 (x0 x1 x2 y0 y1 y2 z0 z1 z2 : α) : α :=
  x0 * (y1 * z2 - y2 * z1) - x1 * (y0 * z2 - y2 * z0) + x2 * (y0 * z1 - y1 * z0)

/-- Area of a sub-triangle: `|det| · A_f`. [formula: cut-subtriangle-area] -/
def subTriangleArea (det faceArea : α) : α := Scalar.abs det * faceArea

/-! ## Weights of the weighted Laplacian (SPEC §3.6, §3.10) -/

/-- Mean of a quantity over the three vertices of a triangle: `ā_t = (a_i + a_j + a_k)/3`.
[formula: cut-weighted-laplacian] -/
def triangleMean (x y z : α) : α := (x + y + z) / 3

/-- Entry `(i, j)`, `i ≠ j`, of `Σ_f w̄_f · (cotan Laplacian of the single face f)`: the faces
carry their own weight `wbar`. With `wbar = 1` it is `laplacianOffDiag`.
[formula: cut-weighted-laplacian] [formula: cut-length-step] -/
def weightedOffDiag (faces : List Face) (len : Nat → Nat → α) (wbar : Face → α) (i j : Nat) : α :=
  -(sumList (faces.map fun f => wbar f * f.cotOpposite len i j)) / 2

/-- Diagonal entry of that matrix on `nv` vertices (rows sum to zero).
[formula: cut-weighted-laplacian] -/
def weightedDiag (nv : Nat) (offDiag : Nat → Nat → α) (i : Nat) : α :=
  laplacianDiag nv offDiag i

/-- Entry of the matrix `W_V + τ L_w` of the length step, `A` the vertex area, `lw` the entry of
`L_w`. [formula: cut-length-step] -/
def lengthStepEntry (isDiag : Bool) (tau areaV lw : α) : α :=
  (if isDiag then areaV else 0) + tau * lw

/-- Entry of the matrix `W_V D⁻¹ + τ L` of the paper's penalty rule, `d = α_L + η_i`.
[formula: cut-penalty-step-paper] -/
def penaltyStepEntry (isDiag : Bool) (tau areaV d l : α) : α :=
  (if isDiag then areaV / d else 0) + tau * l

/-! ## Gradient densities (SPEC §3.5, §3.6) -/

/-- `∂u/∂n ≈ v_i / ℓ_i`. [formula: cut-normal-derivative] -/
def normalDerivative (v ell : α) : α := v / ell

/-- Lumped boundary length `ℓ_i = ½ Σ (boundary edges at i)` from the two (or one) lengths.
[formula: cut-normal-derivative] -/
def dualLength (edgeLengths : List α) : α := sumList edgeLengths / 2

/-- `G = (v/ℓ)²`. [formula: cut-gradient-dirichlet] -/
def gradDirichlet (v ell : α) : α := (v / ell) * (v / ell)

/-- `G = −a (v/ℓ)² − (v/ℓ)(q/ℓ)`. [formula: cut-gradient-scaled] -/
def gradScaled (a v q ell : α) : α := -(a * ((v / ell) * (v / ell))) - (v / ell) * (q / ell)

/-- `G = b(0) − (v/ℓ)(q/ℓ)`, `b0 = b(0)`. [formula: cut-gradient-pointwise] -/
def gradPointwise (b0 v q ell : α) : α := b0 - (v / ell) * (q / ell)

/-- Right-hand side of the adjoint at a boundary vertex of the scaled energy:
`r_i = −2 ((L_a u)_i − a_i v_i)`. [formula: cut-adjoint-scaled] -/
def adjointRhsScaledBoundary (lau a v : α) : α := -2 * (lau - a * v)

/-- Right-hand side of the adjoint at an interior vertex of the scaled energy:
`r_i = −2 (L_a u)_i`. [formula: cut-adjoint-scaled] -/
def adjointRhsScaledInterior (lau : α) : α := -2 * lau

/-- Right-hand side of the pointwise adjoint: `r_i = −A_i^S b′(u_i)`.
[formula: cut-adjoint-pointwise] -/
def adjointRhsPointwise (areaS bPrime : α) : α := -(areaS * bPrime)

/-- Integrated normal derivative of the adjoint at a boundary vertex: `q_i = (L_BI p_I)_i − r_i`.
[formula: cut-adjoint-scaled] -/
def adjointFlux (lbi r : α) : α := lbi - r

/-- `v_i = (L_BI u_I)_i + g_i`. [formula: cut-normal-derivative] -/
def normalFlux (lbi g : α) : α := lbi + g

/-! ## Curvature-based weight (SPEC §3.6) -/

/-- `K_i = (2π − Θ_i)/A_i` at an interior vertex. [formula: cut-feature-size] -/
def vertexCurvatureDensity (theta areaV : α) : α := (2 * Scalar.pi - theta) / areaV

/-- `r⁰ = 1 / (√|K| + 1/r_ε)`. [formula: cut-feature-size] -/
def featureSizeRaw (rEps K : α) : α := 1 / (Scalar.sqrt (Scalar.abs K) + 1 / rEps)

/-- `a = r³`. [formula: cut-feature-size] -/
def featureWeight (r : α) : α := r * r * r

/-! ## Speed, smoothing, step (SPEC §3.7 – §3.10, §4) -/

/-- `s_c = G_0 − G_1`. [formula: cut-speed] -/
def cutSpeed (g0 g1 : α) : α := g0 - g1

/-- Region total `G_r = w_r G + β_r`. [formula: cut-speed] -/
def regionDensity (w G beta : α) : α := w * G + beta

/-- Aspect ratio `(longest edge)² √3 / (4 area)`. [formula: cut-sliver-smoothing] -/
def aspectRatio (longest area : α) : α := longest * longest * Scalar.sqrt 3 / (4 * area)

/-- Arclength-weighted mean of the speeds of the two neighbours along the cut:
`(ℓ₂ s₁ + ℓ₁ s₂)/(ℓ₁ + ℓ₂)`. [formula: cut-sliver-smoothing] -/
def smoothedSpeed (l1 l2 s1 s2 : α) : α := (l2 * s1 + l1 * s2) / (l1 + l2)

/-- `φ̃ = φ + τ σ`. [formula: cut-explicit-step] -/
def explicitStep (phi tau sigma : α) : α := phi + tau * sigma

/-- `τ = min(τ_max, c · h / max|s|)`. [formula: cut-step-cap] -/
def stepCap (tauMax maxMotion h maxSpeed : α) : α := min tauMax (maxMotion * h / maxSpeed)

/-- Symmetry projection: the mean of `φ′` over a class. [formula: cut-symmetry-projection] -/
def classMean (values : List α) (count : α) : α := sumList values / count

/-! ## Local Hopf–Lax update (SPEC §3.12) -/

/-- Edge term `T_j + ℓ_ij`. [formula: cut-redistance] -/
def hopfLaxEdge (Tj lij : α) : α := Tj + lij

/-- `min_{λ∈[0,1]} λ T_j + (1−λ) T_k + |x_i − (λ x_j + (1−λ) x_k)|` for the triangle `i j k` with
sides `lij lik ljk`, in closed form from the intrinsic geometry.

With `d = x_k − x_i`, `e = x_j − x_k`, the function is `T_k + λ ΔT + |d + λ e|`, `ΔT = T_j − T_k`,
which is convex. Put `a = |e|`, `λ₀ = −d·e/a²` and `h² = |d|² − (d·e)²/a²` (the squared distance of
`x_i` from the line `jk`); `d·e = (lij² − lik² − ljk²)/2`. If `|ΔT| < a` the unconstrained minimum
is at `λ* = λ₀ − ΔT h / (a √(a² − ΔT²))`, and the minimum over `[0,1]` is the value at `λ*`
clamped to `[0,1]`; otherwise the minimum is at an end point of `[0,1]`.
[formula: cut-redistance] -/
def hopfLaxTriangle (lij lik ljk Tj Tk : α) : α :=
  let dd := lik * lik
  let de := (lij * lij - lik * lik - ljk * ljk) / 2
  let a2 := ljk * ljk
  let lam0 := -de / a2
  let h2 := max 0 (dd - de * de / a2)
  let h := Scalar.sqrt h2
  let a := ljk
  let dT := Tj - Tk
  let value := fun (lam : α) =>
    Tk + lam * dT + Scalar.sqrt (max 0 (dd + 2 * lam * de + lam * lam * a2))
  if Scalar.abs dT < a then
    let lam := lam0 - dT * h / (a * Scalar.sqrt (a2 - dT * dT))
    value (max 0 (min 1 lam))
  else
    min (value 0) (value 1)

/-! ## Additions for SPEC revision f1d56b5 (flux recovery, cancellation-free forms) -/

/-- Cancellation-free form of `Ω_f`: `Σ θ_i (2π − Θ_i)/Θ_i` over the interior corners `i` of the
face. Equal to `faceCurvature` of the normalized angles up to rounding, because the corner angles
of a Euclidean triangle sum to `π`. [formula: cut-face-curvature] -/
def Face.omegaCancellationFree (len : Nat → Nat → α) (angleSum : Nat → α) (isB : Nat → Bool)
    (f : Face) : α :=
  let term := fun (isBd : Bool) (theta total : α) =>
    if isBd then (0 : α) else theta * (2 * Scalar.pi - total) / total
  term (isB f.a) (f.angleA len) (angleSum f.a) + term (isB f.b) (f.angleB len) (angleSum f.b)
    + term (isB f.c) (f.angleC len) (angleSum f.c)

/-- Inner product of two barycentric differences `d`, `e` within a face with sides `lab lbc lca`
(the polarization of `baryLength2`):
`−½ [(d_a e_b + d_b e_a) ℓab² + (d_b e_c + d_c e_b) ℓbc² + (d_c e_a + d_a e_c) ℓca²]`.
[formula: cut-barycentric-length] -/
def baryInner (lab lbc lca da db dc ea eb ec : α) : α :=
  -(((da * eb + db * ea) * (lab * lab) + (db * ec + dc * eb) * (lbc * lbc)
      + (dc * ea + da * ec) * (lca * lca)) / 2)

/-- Cotangent at the corner `x` of a sub-triangle `xyz` from `⟨y − x, z − x⟩` and its area:
`⟨y − x, z − x⟩ / (2 · area)`. [formula: cut-barycentric-length] -/
def cotFromInner (inner area : α) : α := inner / (2 * area)

/-- `h_f² = (ℓab² + ℓbc² + ℓca²)/3`. [formula: cut-flux-recovery] -/
def faceMeanSquareEdge (lab lbc lca : α) : α := (lab * lab + lbc * lbc + lca * lca) / 3

/-- Coupling of a cut segment: `κ_σ = c² h_f² / |pq|`. [formula: cut-flux-recovery] -/
def fluxCoupling (c h2 segLen : α) : α := c * c * h2 / segLen

/-- Diagonal entry of the flux-recovery matrix at the cut vertex `p`:
`ℓ_p + Σ_{σ ∋ p} κ_σ`. [formula: cut-flux-recovery] -/
def fluxDiagEntry (ell kappaSum : α) : α := ell + kappaSum

/-- Off-diagonal entry `−κ_σ` of the flux-recovery matrix for the segment `σ = pq`.
[formula: cut-flux-recovery] -/
def fluxOffEntry (kappa : α) : α := -kappa

/-- Right-hand side of the flux-recovery system: the integrated flux (`v_p`, or `q_p` for the
adjoint). [formula: cut-flux-recovery] -/
def fluxRhs (flux : α) : α := flux

/-- Gradient densities written with the recovered derivatives `d`, `e`.
[formula: cut-gradient-dirichlet] -/
def gradDirichletD (d : α) : α := d * d

/-- [formula: cut-gradient-scaled] -/
def gradScaledD (a d e : α) : α := -(a * (d * d)) - d * e

/-- [formula: cut-gradient-pointwise] -/
def gradPointwiseD (b0 d e : α) : α := b0 - d * e

/-! ## Additions for SPEC revision da35966 (default stages, N regions) -/

/-- Corner angle from its cotangent: `θ = π/2 − arctan(cot θ)` (valid for `θ ∈ (0, π)`).
[formula: cut-length-step-on-cut] -/
def angleFromCot (cot : α) : α := Scalar.pi / 2 - Scalar.atan cot

/-- Turning angle of the cut at a cut point: `π − β`, `β` the interior angle of the lower region;
at a cut point on `∂M` the reflected cut turns by `π/2 − β`. In general `γ − β` with `γ = π`
(cut point), `π/2` (cut point on `∂M`), `2π/3` (triple point of §8.5).
[formula: cut-length-step-on-cut] [formula: cut-length-step-regions] -/
def turningAngle (gamma beta : α) : α := gamma - beta

/-- Derivative out of a region of a linear `w` across a segment `pq`, from the triangle of the
region at the segment with third vertex `z`: `h` its height over `pq`, `f` the foot at the
fraction `lam` of `pq` from `p`, `w_f = (1−lam) w_p + lam w_q`; the derivative is
`(w_f − w_z)/h`. [formula: cut-length-step-on-cut] -/
def normalDerivativeOfLinear (wp wq wz lam h : α) : α := ((1 - lam) * wp + lam * wq - wz) / h

/-- Fraction of `pq` at which the foot of the height from the third vertex lies, from the
cotangents of the triangle at `p` and `q`: `cot_p / (cot_p + cot_q)`.
[formula: cut-length-step-on-cut] -/
def footFraction (cotP cotQ : α) : α := cotP / (cotP + cotQ)

/-- Height of the third vertex over the segment: `2 · area / |pq|`.
[formula: cut-length-step-on-cut] -/
def heightOver (area segLen : α) : α := 2 * area / segLen

/-- The integrated length force at a cut point (or boundary node):
`Θ = w (γ − β) + Σ_σ ½ |σ| (∂w/∂n)_σ`; here for the two (or one) segments, given as the list of
`(|σ|, ∂w/∂n)`. [formula: cut-length-step-on-cut] [formula: cut-length-step-regions] -/
def lengthForce (w gamma beta : α) (segs : List (α × α)) : α :=
  w * turningAngle gamma beta + sumList (segs.map fun s => s.1 * s.2 / 2)

/-- Mean `w̄_σ = (w_p + w_q)/2` of a segment. [formula: cut-length-step-on-cut] -/
def segmentMeanWeight (wp wq : α) : α := (wp + wq) / 2

/-- Coupling of the rate system on a segment: `τ w̄_σ / |σ|`.
[formula: cut-length-step-on-cut] -/
def rateCoupling (tau wbar segLen : α) : α := tau * wbar / segLen

/-- Right-hand side of the rate system: `ℓ (s + λ⁰)`. [formula: cut-length-step-on-cut] -/
def rateRhs (ell s lam : α) : α := ell * (s + lam)

/-- `|∇φ|` in a crossed face whose lone vertex `a` is cut off by `pq`:
`|φ_a| · |pq| / (2 · area(a, p, q))`. [formula: cut-gradient-norm] -/
def gradientNormLone (phiA segLen area : α) : α := Scalar.abs phiA * segLen / (2 * area)

/-- `ḡ = Σ|pq| / Σ(|pq|/|∇φ|)`. [formula: cut-gradient-norm] -/
def meanGradient (sumLen sumLenOverGrad : α) : α := sumLen / sumLenOverGrad

/-- Normalized weighted mean `Σ w v / Σ w` of the pairs `(w_k, v_k)` of the two lists (the longer list
is cut to the length of the shorter one, so the denominator is the sum of the weights that are used):
the first line of the local extension (`w = ℓ_p β_ip`, `v = r_p`). [formula: cut-extension-local] -/
def weightedMean (weights values : List α) : α :=
  let pairs := weights.zip values
  sumList (pairs.map fun x => x.1 * x.2) / sumList (pairs.map fun x => x.1)

/-- Second line of the local extension: `Σ_j (σ_j/ℓ_ij) / Σ_j (1/ℓ_ij)` over the neighbours one ring
closer, given as `(σ_j, ℓ_ij)`. [formula: cut-extension-local] -/
def ringMean (nbrs : List (α × α)) : α :=
  sumList (nbrs.map fun x => x.1 / x.2) / sumList (nbrs.map fun x => 1 / x.2)

/-- Quadratic part of the gradient density `Q = G − β − w b(0)`. [formula: cut-roughening-rate] -/
def roughQ (G beta w b0 : α) : α := G - beta - w * b0

/-- Roughening rate `ω = 2 (Q_0 + Q_1)`. [formula: cut-roughening-rate] -/
def roughRate (q0 q1 : α) : α := 2 * (q0 + q1)

/-- Length weight that makes ripples longer than `ℓ` grow: `α_L = ω ℓ/(2π)`.
[formula: cut-roughening-rate] -/
def lengthWeightForRipples (omega ell : α) : α := omega * ell / (2 * Scalar.pi)

/-- Adaptive step: `min(τ_max, time_scale/ω, max_motion · h / max|ρ|)` (each of the last two only
if its denominator is positive). [formula: cut-step-cap] -/
def adaptiveStep (tauMax timeScale omega maxMotion h maxRho : α) : α :=
  let a := if 0 < omega then min tauMax (timeScale / omega) else tauMax
  if 0 < maxRho then min a (maxMotion * h / maxRho) else a

/-- Crossing parameter between regions `i = region a` and `j = region b` (paper Eq. 31):
`ψ_a/(ψ_a − ψ_b)` with `ψ_a = φ_i(a) − φ_j(a)`, `ψ_b = φ_i(b) − φ_j(b)`.
[formula: cut-edge-crossing-regions] -/
def crossingParamRegions (phiIa phiJa phiIb phiJb : α) : α :=
  (phiIa - phiJa) / ((phiIa - phiJa) - (phiIb - phiJb))

/-- Barycentric coordinate of the triple point at the corner `c_k` of a face whose cut points on
the edges leaving and arriving at `c_k` sit at the parameters `s_k` (from `c_k`) and `s_prev`
(from `c_{k−1}`, so the point is at `1 − s_prev` from `c_k`): `((1 − s_k) + s_prev)/3`.
[formula: cut-triple-point] -/
def triplePointCoord (sNext sPrev : α) : α := ((1 - sNext) + sPrev) / 3

/-- Rate at the copy of a triple point `m` in a region from the two cut points `p`, `p′` of the face
on its boundary: `(|mp′| s_p + |mp| s_p′)/(|mp| + |mp′|)`. [formula: cut-region-extension] -/
def tripleRate (lenMp lenMpPrime sp spPrime : α) : α :=
  (lenMpPrime * sp + lenMp * spPrime) / (lenMp + lenMpPrime)

/-- Mean of the two smallest values at a vertex (Losasso et al.): `m`.
[formula: cut-region-projection] -/
def projectionShift (smallest second : α) : α := (smallest + second) / 2

/-- Projected value `φ_k(i) = φ_k′(i) − m`. [formula: cut-region-projection] -/
def projectedValue (phi m : α) : α := phi - m

/-! ## Additions for SPEC revision 36f50cc (robust step statistics) -/

/-- Insertion of `(x, ℓ)` into a list sorted ascending by `x`, after all entries with an equal key
(so that equal values stay in their original order). [formula: cut-weighted-quantile] -/
def insertAscending (e : α × α) : List (α × α) → List (α × α)
  | [] => [e]
  | f :: fs => if e.1 < f.1 then e :: f :: fs else f :: insertAscending e fs

/-- Stable ascending sort by the first component. [formula: cut-weighted-quantile] -/
def sortAscending : List (α × α) → List (α × α)
  | [] => []
  | e :: es => insertAscending e (sortAscending es)

/-- Sum from the left, `((0 + x₀) + x₁) + ⋯`, the order of SPEC §1's running sums.
[formula: cut-weighted-quantile] -/
def sumLeft (xs : List α) : α := xs.foldl (fun a x => a + x) 0

/-- The first value of a list sorted ascending whose running weight reaches `target`; the last
value if none does; `0` for an empty list. [formula: cut-weighted-quantile] -/
def firstReaching (target : α) : α → List (α × α) → α
  | _, [] => 0
  | acc, [e] => if target ≤ acc + e.2 then e.1 else e.1
  | acc, e :: es => if target ≤ acc + e.2 then e.1 else firstReaching target (acc + e.2) es

/-- Weighted quantile of SPEC §1: the values `x_c` in ascending order (equal values in the order of the
list), the first `x_c` at which the weights summed so far (from the left, in that order) reach
`q · Σ ℓ`, `Σ ℓ` summed in the same order; `0` for an empty list. With `q = 1` the result is the
largest value (a trailing zero weight excepted). [formula: cut-weighted-quantile] -/
def weightedQuantile (q : α) (xs ls : List α) : α :=
  let sorted := sortAscending (xs.zip ls)
  firstReaching (q * sumLeft (sorted.map fun e => e.2)) 0 sorted

/-- Speed of a cut point in crossed edges per unit time, `|ρ|/h_c`. [formula: cut-step-local] -/
def speedInEdges (rho hc : α) : α := Scalar.abs rho / hc

/-- Robust step: `min(τ_max, time_scale/ω, max_motion/m̂)` with every term whose denominator is not
positive skipped. [formula: cut-step-local] -/
def robustStep (tauMax timeScale omega maxMotion mhat : α) : α :=
  let a := if 0 < omega then min tauMax (timeScale / omega) else tauMax
  if 0 < mhat then min a (maxMotion / mhat) else a

/-- Speed that the local limit of a cut point is measured with: `|r_c|` for the capped rule,
`max(|r_c|, ½ (ω_c − ω)₊)` for the adaptive one. [formula: cut-step-local] -/
def localSpeed (r omegaC omega : α) : α := max (Scalar.abs r) (max 0 (omegaC - omega) / 2)

/-- Factor of the local limit, `min(1, max_motion · h_c / (τ v_c))`; `1` where `v_c` or `τ` is not
positive. [formula: cut-step-local] -/
def localLimitFactor (maxMotion hc tau v : α) : α :=
  if 0 < tau * v then min 1 (maxMotion * hc / (tau * v)) else 1

end Ddg
