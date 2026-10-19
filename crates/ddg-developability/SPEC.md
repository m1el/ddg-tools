# Developability of Triangle Meshes — discrete specification

The contract between the Rust crate (`crates/ddg-developability`), the Lean oracle
(`lean/oracle`), the Lean proofs (`lean/proofs`), the Python binding (`crates/ddg-py`) and the
Blender add-on for

> Oded Stein, Eitan Grinspun, Keenan Crane. *Developability of Triangle Meshes.* ACM
> Transactions on Graphics 37(4), Article 77, 2018. <https://doi.org/10.1145/3197517.3201303>
> (`../DiscreteDevelopable/paper.pdf`; "the paper", equation numbers are the paper's)

**Clean room.** This file is written from the text of the paper only (formulas transcribed from
the page images, not from a text extraction). The authors' code is GPL-2.0 with an LGPL 3x3
eigensolver; it is quarantined outside this tree and must never be opened, searched, fetched or
paraphrased by anyone working on this crate. Kopp's paper (2008) may be read and cited, his
published C code may not be copied. Papers the text cites may be downloaded and cited; this
file uses Lewis & Overton 2013 (Algorithm 2.6 and Sec. 5.1 of the preprint
`cs.nyu.edu/~overton/papers/pdffiles/bfgs_inexactLS.pdf`), Nocedal & Wright 2006 (Ch. 3, 7),
Kopp 2008 (arXiv physics/0610206, abstract only), Erickson & Har-Peled 2004 and Sawhney & Crane
2017 (through `ddg-bff`).

Every formula below was re-derived; where the printed paper differs, the difference has an id
`D<n>` in §8 with the printed version, ours, and the evidence. If an implementation or a test
contradicts this file, stop and report it; do not quietly follow either side. Numeric checks
of every claim marked *(checked)* were run with a throw-away NumPy script (finite differences,
dense search), not part of the repository.

Conventions of `CLAUDE.md` apply: `f64` everywhere, stored indices `u32`, manifold oriented
meshes only (`ddg_core::TriMesh`), bad input gives an error enum and never a panic, no silent
clamping, normalization or rescaling (in particular, positions are never rescaled; tolerances
are expressed in units of the mean edge length, §4.1), every solve returns a report, no global
state (thread count is an option), formula tags `[formula: dev-…]` (§9).

## 1. Discrete setting

Mesh `M = (V, E, F)` with positions `f_v ∈ R³`. A face is a CCW triple `(i, j, k)` as stored in
`TriMesh::faces()`.

**Per face** (`σ = (i,j,k)`):

```
n_σ = (f_j − f_i) × (f_k − f_i)          (= f_i×f_j + f_j×f_k + f_k×f_i, cyclic)
A_σ = |n_σ| / 2                           area
N_σ = n_σ / |n_σ|                         unit normal          [formula: dev-face-frame]
e_i^σ = f_k − f_j,  e_j^σ = f_i − f_k,  e_k^σ = f_j − f_i     edge opposite each corner, CCW
θ_i^σ = atan2(|a × b|, a · b),  a = f_j − f_i, b = f_k − f_i  corner angle   [formula: dev-corner-angle]
```

The angle is computed with `atan2` from positions (accurate for small angles; D13), not from
edge lengths. A face with `n_σ = 0` (or non-finite) is degenerate: the normal is undefined.

**Vertex star.** `St(v)` = the faces incident to `v` in the CCW order of
`TriMesh::outgoing_halfedges(v)` (faces `σ_0 … σ_{k−1}`, `k = vertex_face_count(v)`). For an
interior vertex the faces form a **cycle** (`σ_m` and `σ_{m+1 mod k}` share the edge from `v`
to a neighbour); for a boundary vertex they form a **path** starting at the boundary halfedge
leaving `v` (`k − 1` shared edges). The *valence* of an interior vertex is `k`. Corner `c` of
the star is the corner of `σ_m` at `v`; corner ids are `3 f + slot` as in `ddg_core`.

**Vertex normal** (area-weighted, paper App. B.5):

```
S_v = ½ Σ_{σ ∈ St(v)} n_σ  (vector area),    N_v = S_v / |S_v|           [formula: dev-vertex-normal]
```

undefined if `S_v = 0`. For an interior vertex `S_v` does not depend on `f_v` (the terms in `f_v`
telescope around the cycle) — a check for the gradient code.

**Hinge (Def. 3.1).** `St(v)` is a hinge if it is embedded and its faces split into two
edge-connected regions of constant normal. *Theorem 3.3:* an embedded star is a hinge iff all
`N_σ`, `σ ∈ St(v)`, lie in a common plane through the origin. *Props. A.1–A.3:* a non-flat
interior hinge has two antiparallel edges; a discrete developable immersion without flat
vertices is discrete ruled; a valence-3 hinge is flat. These are cited, not used by code
(`dev-hinge-coplanar`, `dev-hinge-antiparallel`, `dev-hinge-ruled`, `dev-valence-three-flat`).
All six energies vanish exactly when the normals of the star are coplanar (for `E_P`: when a
two-arc split with constant normal per arc exists), so they also vanish on the non-embedded
configurations of Fig. 9 (Spike limit, DoubleCover); embeddedness is not checked (D28).

## 2. Derivatives of the basic quantities, and the adjoint form

Directional derivatives for a motion `δ` of one corner `p ∈ σ`, `e = e_p^σ` (all *checked*):

```
d n_σ   = e × δ                                                              (n is bilinear)
∇_{f_p} A_σ = ½ N_σ × e                                        (Eq. 7)     [formula: dev-area-gradient]
d N_σ   = (N_σ·δ) (e × N_σ) / (2 A_σ),  Jacobian  ∂N_σ/∂f_p = (e × N_σ) N_σᵀ / (2 A_σ)
                                                               (Eq. 8, D2) [formula: dev-normal-gradient]
∇_{f_j} θ_i^σ = N_σ × (f_i − f_j) / |f_i − f_j|²
∇_{f_k} θ_i^σ = N_σ × (f_k − f_i) / |f_k − f_i|²
∇_{f_i} θ_i^σ = −(∇_{f_j} θ_i^σ + ∇_{f_k} θ_i^σ)               (Eq. 9, D3) [formula: dev-angle-gradient]
```

Derivation of Eq. 8: `dN = (I − N Nᵀ) dn / |n|`; split `δ = (N·δ)N + t` with `t ⊥ N`; `e × t ∥ N`
is removed by the projector, `e × N` is tangent, so `dN = (N·δ)(e × N)/|n|`. Eq. 9: moving
`f_j` rotates `f_j − f_i` about `N`; the angle grows for clockwise motion `−N × (f_j − f_i)`, at
rate `1/|f_j − f_i|` per unit length, i.e. gradient `N × (f_i − f_j)/|f_i − f_j|²`.

**Adjoint (back-propagation) form.** Every energy is a sum of per-vertex terms `e_v` that depend
on the positions only through the `N_σ`, `θ_c` and (intrinsic energies) `S_v` of the star. The
evaluation computes, per star, the *sensitivities* `y_c = ∂e_v/∂N_σ ∈ R³` and `s_c = ∂e_v/∂θ_c`
per corner `c` of the star and `z_v = ∂e_v/∂S_v ∈ R³`, and then per face

```
y_σ = Σ_{c ∈ σ} y_c                     (the three corners of σ belong to three different stars)
w_σ = (y_σ − (y_σ·N_σ) N_σ) / |n_σ|  +  ½ Σ_{u ∈ σ} z_u         (sensitivity w.r.t. n_σ)
g_p += w_σ × e_p^σ                      for each corner p of σ   (since w·(e×δ) = δ·(w×e))
g_j += s_c N_σ × (f_i − f_j)/|f_i − f_j|²,  g_k += s_c N_σ × (f_k − f_i)/|f_k − f_i|²,
g_i −= (both)                            for each corner c = i of σ                [formula: dev-backprop-face]
```

(Eq. 7 is the special case `w = N/2`.) §3.9 fixes the order of these sums.

## 3. The energies

`E(f) = Σ_{v ∈ D} e_v(f)`. All six per-vertex energies are dimensionless, invariant under rigid
motions and uniform scaling of `f` (they depend only on unit normals and angles), hence
`Σ_v ∇_v E = 0`, `Σ_v f_v × ∇_v E = 0`, `Σ_v f_v · ∇_v E = 0` [formula: dev-gradient-invariance]
(checked by tests; the flow does not translate, rotate or scale the mesh to first order).

### 3.1 Domain `D` (D8)

* Every vertex is in `D`, interior and boundary (paper: "sum over all vertices", `E^λ := Σ_{i∈V}`),
  with option `include_boundary` (default `true`).
* **Valence 3** (Sec. 4.1.4 "Triple points", Fig. 16): interior vertices with exactly 3 faces are
  left out, option `omit_valence_three` (default `true`). A boundary vertex with ≤ 2 faces has
  `e_v = 0` for all six energies anyway (two normals are always coplanar, and `N_v` lies on the
  arc between them), so the rule concerns interior vertices only.
* Per-vertex energies are reported for every vertex; vertices outside `D` report the value of
  the formula with a flag `counted = false` (the cut, §6, uses them) and contribute nothing to `E`
  or its gradient. Where the formula is undefined at an uncounted vertex (§3.6) its entry is
  `NaN`; that is not an error. At a counted vertex it is (`DevError::UndefinedEnergy`).

### 3.2 Combinatorial energy `E_P` (Eq. 1–2)

**Bipartitions** `P_v` [formula: dev-star-bipartitions]: pairs `{F_1, F_2}` of complementary,
non-empty, edge-connected sets of faces of `St(v)` (D9):

* interior, valence `k`: complementary arcs of the cycle. Enumerated as `(a, b)`,
  `0 ≤ a < b ≤ k − 1`, `F_1 = {σ_{a+1}, …, σ_b}`, `F_2` = the rest; `|P_v| = k(k−1)/2`
  (3 for k = 3, 15 for k = 6);
* boundary, `k` faces: split of the path after face `a`, `0 ≤ a ≤ k − 2`, `F_1 = {σ_0 … σ_a}`;
  `|P_v| = k − 1`. A boundary star with one face has no bipartition and `e_v = 0`.

An empty part is not allowed: it would only add `1 − |N̄_{St}|²`, which is zero only for flat
stars, where a two-arc split is zero as well.

**Width of a part** (`n_p = |F_p|`, `N̄_p = (1/n_p) Σ_{σ∈F_p} N_σ`):

```
π(P) = Σ_{p=1,2} (1/n_p) Σ_{σ∈F_p} |N_σ − N̄_p|²
     = Σ_{p=1,2} (1/n_p²) Σ_{σ1<σ2 ∈ F_p} |N_σ1 − N_σ2|²      (unordered pairs; D1)
     = Σ_{p=1,2} (1 − |N̄_p|²)                                   [formula: dev-combinatorial-width]
e_v = E_P,v = min_{P ∈ P_v} π(P)                                 (Eq. 2)
```

The code evaluates the first form (no cancellation when the part is nearly flat). With prefix
sums of `N_σ` around the star every arc mean costs O(1); the variance of an arc needs its own
O(n_p) pass, so a star costs O(k³) worst case (k ≈ 6: negligible).

**Gradient** [formula: dev-combinatorial-gradient]: for the minimizing `P`,
`dπ = (2/n_p) Σ_σ (N_σ − N̄_p)·dN_σ − (2/n_p) dN̄_p·Σ_σ(N_σ − N̄_p)` and the last sum is 0, so

```
y_c = (2/n_p) (N_σ − N̄_p)    for the corner c of v in σ ∈ F_p;      s_c = 0.
```

The paper's B.2 form `∇|N_σ1 − N_σ2|² = 2(∇N_σ1 − ∇N_σ2)ᵀ(N_σ1 − N_σ2)` summed over unordered pairs
with `1/n_p²` gives the same (*checked*).

### 3.3 Combinatorial max energy `π_max` (Eq. 5)

```
π_max(P) = max_{p} max_{σ1, σ2 ∈ F_p} |N_σ1 − N_σ2|²,      e_v = min_{P ∈ P_v} π_max(P)
                                                                 [formula: dev-combinatorial-max]
```

(Eq. 5 writes `(N_2 − N_1)²` for the squared norm.) Evaluation: `D_{στ} = |N_σ − N_τ|²` for all
pairs of the star; arc maxima by extension of arcs, O(k³). Subgradient: the active pair
`(σ1, σ2)` of the minimizing partition: `y_{σ1} = 2(N_σ1 − N_σ2)`, `y_{σ2} = −2(N_σ1 − N_σ2)`,
all other `y = 0`, `s = 0`. A part with one face has width 0.

### 3.4 Covariance energy `E_λ` (Eq. 3–4)

```
A_v = Σ_{σ∈St(v)} θ_v^σ N_σ N_σᵀ                 (3x3, symmetric PSD)   [formula: dev-covariance-matrix]
e_v = λ_v = min_{|u|=1} Σ θ ⟨u, N_σ⟩² = λ_min(A_v)          (Eq. 3 = Eq. 4; Rayleigh)
```

Eigen-decomposition by `eig::sym_eig3` (cyclic Jacobi, already in the crate; D29). `x` = its unit
eigenvector of the smallest eigenvalue. The computed `λ` may be negative by rounding
(`≈ −ε‖A‖`); it is **not** clamped (D14).

**Gradient** (Eq. 10) [formula: dev-covariance-gradient]: for a simple smallest eigenvalue,
`dλ = xᵀ dA x` (`∇_A λ = x xᵀ`, [formula: dev-simple-eigenvalue-derivative], cited; the
supergradient inequality `λ_min(B) ≤ λ_min(A) + xᵀ(B − A)x` holds always and is proved,
[formula: dev-eigenvalue-supergradient]); with `dA = Σ dθ N Nᵀ + θ(dN Nᵀ + N dNᵀ)`:

```
y_c = 2 θ_v^σ (x·N_σ) x,      s_c = (x·N_σ)²
```

which with §2 is Eq. 10 (with the corrected Eq. 8–9). Angle weights make `A_v` invariant under
retessellation of a planar region: splitting a corner into two corners with the same normal
leaves `A_v` unchanged [formula: dev-angle-weight-retessellation].

### 3.5 Max covariance energy `λ_max` (Eq. 6, App. B.4)

```
φ(u) = max_{σ∈St(v)} ⟨u, N_σ⟩²,      e_v = λ_max,v = min_{|u|=1} φ(u)
```

**Why candidates suffice** [formula: dev-max-covariance-voronoi] (proof sketch, cited as
reasoning; *checked* against dense search on 200 random stars, never worse). Let
`S = {±N_σ}` (the sites). `⟨u,N⟩² = cos² ∠(u, ±N)`, so `φ(u) = cos² d(u)` with `d(u)` the
spherical distance from `u` to the nearest site; minimizing `φ` is the largest empty circle
problem. Inside the Voronoi cell of a site `s`, `d = ∠(u, s)` has no interior maximum unless
`u ⊥ s`, which inside the cell forces every normal `⊥ u` (all coplanar, `φ = 0`). On a Voronoi
edge (bisector of `s, t`, parametrized `u = cos α m + sin α ŵ`, `m = (s+t)/|s+t|`,
`ŵ = s×t/|s×t|`), `|⟨u,s⟩| = |cos α| ⟨m,s⟩` is monotone in `|α|`, so `d` is maximal at an
endpoint (a Voronoi vertex) or at `u = ±ŵ`, which again forces all normals `⊥ ŵ`. Hence the
minimum of `φ` is attained at a Voronoi vertex, i.e. at a point equidistant from (at least) three
sites, or — when all normals are parallel — `φ ≡ 0` on a great circle. **Midpoints of pairs are
never needed** (they are minima of `d` along edges). And since every unit vector `w` gives
`φ(w) ≥ min φ`, the minimum over any candidate set that contains the Voronoi vertices is exact:
extra candidates (non-Voronoi triples) are harmless.

**Circumcenter** [formula: dev-spherical-circumcenter]: for sites `a, b, c`,
`n = (b − a) × (c − a)`, `v = n/|n|` satisfies `⟨v,a⟩ = ⟨v,b⟩ = ⟨v,c⟩ = det(a,b,c)/|n| =: g`
(since `n·a = (b×c)·a`); `|n| = 2 A_abc` (Euclidean area of the triangle `abc`). The paper's
"spherical centroid" is this circumcenter (D5). The antipode `−v` is the other equidistant
point; `φ(v) = φ(−v)`.

**Candidates** (D11), with `n = |St(v)|`, faces indexed `0 … n−1` in star order:

1. for every face triple `σ1 < σ2 < σ3` and signs `(s2, s3) ∈ {+,−}²` (order `++, +−, −+, −−`):
   sites `(a, b, c) = (N_σ1, s2 N_σ2, s3 N_σ3)` (flipping all three signs gives the same `±v`);
2. for every pair `σ1 < σ2`: sites `(N_σ1, −N_σ1, N_σ2)`, i.e. `v ∥ N_σ1 × N_σ2` (the pole; it is
   the minimizer only when all normals are coplanar, then `φ = 0`).

`4·C(n,3) + C(n,2)` candidates (95 for n = 6). A candidate with `n = 0` exactly (or non-finite) is
skipped; nothing else is skipped (a nearly degenerate triple still yields a unit vector and an
upper bound). If no candidate remains, all normals are parallel or antiparallel and `e_v = 0`
with zero sensitivities. Otherwise for each candidate `φ(v)` (max over all `n` faces) and
`slack = φ(v) − ⟨v,a⟩²` (`≥ 0`; zero when `a,b,c` are among the nearest sites); the chosen
candidate minimizes the key `(φ(v), slack, enumeration index)` lexicographically.
`e_v = φ(v)` of the chosen candidate.

**Gradient** [formula: dev-max-covariance-gradient] (*checked*, relative FD error 6e-10). At a
generic minimizer `v` is the circumcenter of the active triple and `e_v = g(a,b,c)²` with
`g = det(a,b,c)/|n|` smooth. With `M = a`, `v = n/|n|` (this orientation), and the edges of the
triangle `abc` opposite each site in CCW order, `e_a = c − b`, `e_b = a − c`, `e_c = b − a`:
`dv = (v·dq)(e_q × v)/|n|` (the normal derivative of §2 for the "triangle" `abc`), so

```
∂e_v/∂q = 2g [ [q = a] v + ((e_q × v)·a / |n|) v ]       for q ∈ {a, b, c}
y_c(σ_q) += s_q ∂e_v/∂q        (q = s_q N_{σ_q}; for candidates of kind 2, a and b are the same face)
s_c = 0
```

The paper's formula omits the site sign `s_q` (D5); it is wrong whenever an active site is an
antipode `−N_σ` (relative error ≈ 2 in the check).

### 3.6 Intrinsic width (App. B.5)

With `n0 = N_v` (§1) and, per face of the star, `c = ⟨n0, N_σ⟩`, `φ_σ = atan2(|n0 × N_σ|, c)`:

```
Ñ_σ = log_{n0}(N_σ) = (φ_σ / sin φ_σ) (N_σ − c n0) ∈ T_{n0}S²     (= φ_σ v_σ of the paper)
                                                                  [formula: dev-log-map]
Ã_v = Σ_σ θ_v^σ Ñ_σ Ñ_σᵀ restricted to T_{n0}S²;  e_v = λ_min(Ã_v)
```

Implementation: an orthonormal tangent basis `t1, t2` of `n0` (any deterministic choice; the
energy and gradient do not depend on it), the 2x2 matrix `Tᵀ Ã T`, `eig::sym_eig2`, and the 3D
eigenvector `x = y1 t1 + y2 t2 ⊥ n0`. Because every `Ñ ⊥ n0`, `Ã x = λ x` in 3D and the
derivative of the basis drops out (`yᵀ(dTᵀÃT + TᵀÃdT)y = 2λ xᵀ dT y = 0` by `TᵀT = I`):
`dλ = xᵀ dÃ x` with the 3D `dÃ` [formula: dev-tangent-eigen-derivative]. Hence, as in §3.4,

```
∂e_v/∂Ñ_σ = 2 θ (x·Ñ_σ) x =: z_σ,      s_c = (x·Ñ_σ)².
```

**Log map adjoint** [formula: dev-log-map-gradient] (*checked*, including `φ` down to 1e-6).
With `G(φ) = φ / sin φ` and `H(φ) = G′(φ)/sin φ = (sin φ − φ cos φ)/sin³ φ` (both smooth at 0:
`G = 1 + φ²/6 + 7φ⁴/360`, `H = 1/3 + 2φ²/15 + O(φ⁴)`; use the series for `φ < 1e-3`):

```
dÑ = G (dN − c dn0 − dc n0) − H dc (N − c n0),   dc = n0·dN + N·dn0
⟹ for a sensitivity z:   K = G (z·n0) + H z·(N − c n0)
   ∂/∂N_σ :   y_c  = G z − K n0
   ∂/∂n0  :   y_n0 += −G c z − K N_σ          (summed over the faces of the star)
   n0 = S_v/|S_v| ⟹  z_v = (y_n0 − (y_n0·n0) n0) / |S_v|        (feeds w_σ in §2)
```

This is the corrected Eq. 11 (D6): in the basis `ν = n0 ×̂ N_σ`, `μ_i = ν ×̂ n0` (= the log
direction `v_σ`), `μ_f = ν ×̂ N_σ`, it reads `∂Ñ/∂N = μ_i μ_fᵀ + (φ/sin φ) ν νᵀ` (as printed) and
`∂Ñ/∂n0 = −(μ_i μ_iᵀ + φ n0 μ_iᵀ + (φ/tan φ) ν νᵀ)` (printed: `μ_i μ_fᵀ` in the first term,
which is wrong by the factor `cos φ` on that component, FD error 1.7). The paper's symbol
`v_i^jk` in the gradient paragraph denotes the binormal `ν`, not the log direction (notation).
The smooth form above has no `1/sin φ` at `φ = 0`.

**Undefined** (D12): `S_v = 0`; or a face with `π − φ_σ < 1e-6` (normal nearly opposite the
vertex normal; the log map is singular). Then the energy is undefined: `evaluate` returns
`DevError::UndefinedEnergy { vertex, reason }`, a line search treats the trial point as
`+∞` (§4.3).

### 3.7 Branching-free intrinsic width (App. B.5.1)

```
ψ(u) = max_{σ,τ ∈ St(v)} ⟨Ñ_σ − Ñ_τ, u⟩² = (max_σ ⟨Ñ_σ,u⟩ − min_σ ⟨Ñ_σ,u⟩)²,
e_v = min_{u ∈ T_{n0}, |u|=1} ψ(u) = (minimum width of conv{Ñ_σ})²     [formula: dev-intrinsic-width]
```

No angle weights. The minimum width of a planar convex polygon is attained in a direction
normal to one of its edges (rotating calipers) [formula: dev-min-width-hull-edge], cited. So
the candidates are `u_ab = n0 × (Ñ_b − Ñ_a)/|Ñ_b − Ñ_a|` for pairs `a < b` with `Ñ_a ≠ Ñ_b`
(exact). The paper enumerates `x_a = ±Ñ_a, x_b = ±Ñ_b`: a superset (every unit `u` is an upper
bound), same minimum (*checked*), 4× the work — we use the pairs (D7). The chosen candidate
minimizes `(ψ(u_ab), a, b)`; within it the maximizing pair is `s = argmax ⟨Ñ,u⟩`,
`t = argmin ⟨Ñ,u⟩` (lowest index on ties). No candidate (all `Ñ` equal, e.g. a flat star) gives
`e_v = 0`. Zero iff the `Ñ` are collinear; since `Σ A_σ sin φ_σ v_σ = 0` (definition of `n0`) a
line of `Ñ`s cannot avoid the origin, so this happens iff the normals are coplanar.

**Gradient** [formula: dev-min-width-gradient] (*checked*). With `w = Ñ_b − Ñ_a`,
`ŵ = w/|w|`, `u = n0 × ŵ`, `d = Ñ_s − Ñ_t`, `Δ = d·u` (`e_v = Δ²`): `du = −(u·dw) ŵ / |w|`
in the plane [formula: dev-width-altitude-gradient], i.e. `∂u/∂x_a = ŵ uᵀ/|w|`,
`∂u/∂x_b = −ŵ uᵀ/|w|` (Jacobian layout). Hence

```
z_s += 2Δ u,   z_t −= 2Δ u,   z_a += 2Δ (d·ŵ)/|w| u,   z_b −= 2Δ (d·ŵ)/|w| u,      s_c = 0
```

and then the log-map adjoint of §3.6. The paper prints `∇_{x_a} u* = −(N_i × w) wᵀ/|w|³`, which
as a Jacobian maps into `span(u)` (impossible for a unit vector); it is the negated transpose of
ours, i.e. right only in gradient layout with the opposite orientation `u = w × N_i/|w|` (D7).
The out-of-plane parts of `dÑ` do not matter (the energy depends on the Gram matrix of the
coplanar `Ñ`s).

### 3.8 Nonsmoothness and ties (D10)

All energies except `E_λ` and the intrinsic width are min/max of smooth pieces; `E_λ` and the
intrinsic width are nonsmooth at repeated smallest eigenvalues (e.g. every regular cone, §10.3).
The code returns the gradient of the active piece chosen by the deterministic rules above
(lowest key), which is an element of the Clarke generalized gradient; Lewis & Overton's
practice ("any tie-breaking rule … generally produces a subgradient", Sec. 5.1). The
`Evaluation` reports per vertex the *margin* of the choice — second-best minus best key value
(partition, candidate or pair), or the eigenvalue gap `λ_2 − λ_1` divided by `tr(A)`; `+∞` when
there is no alternative — and counts *near ties*, margins `≤ 1e-12 · max(1, e_v)`. Margins are
for reporting and tests only; they never change a result.

### 3.9 Evaluation, parallelism, determinism

1. Faces (parallel): `n, A, N, θ` for every face. A degenerate face is an error naming it
   (`DevError::DegenerateFace { face }`).
2. Stars (parallel): `e_v`, `y_c`, `s_c` per corner of the star (corner storage is per corner id,
   each corner belongs to one star), `z_v`.
3. Faces (parallel): `w_σ` and the three per-corner gradient contributions `q_c` (§2), in the
   fixed order: normal term, then angle terms of corners 0, 1, 2.
4. Vertices (parallel): `g_v = Σ q_c` over the corners at `v` in `outgoing_halfedges` order; pinned
   vertices get `0` (§4.8). `E = Σ_{v∈D} e_v` summed sequentially in vertex order.

Results are bitwise identical for any thread count (`threads: Option<usize>`, `None` = rayon's
default pool size, run in a local `rayon::ThreadPool`, no global pool configuration). Cost per
vertex: `E_λ` one `sym_eig3` (≈ 260 ns measured) + O(k); `λ_max` O(k⁴) ≈ 600 dot products at
k = 6; `E_P` O(k³).

### 3.10 Eigensolvers

`eig::sym_eig3` / `eig::sym_eig2` (cyclic Jacobi with Rutishauser's rotation, exact power-of-two
pre-scaling, ascending eigenvalues, right-handed orthonormal eigenvectors), already in the crate
[formula: dev-jacobi-rotation]. Kopp's paper compares Jacobi, QL, Cuppen, analytic (Cardano) and
hybrid methods and finds Jacobi the most accurate; the paper's remark that an accurate solver is
needed for `E_λ` is met by Jacobi (D29). Eigenvectors are meaningful only up to the gap
(Davis–Kahan); ties are handled by §3.8.

## 4. Optimization

Variables: all vertex positions (`3|V|` reals) minus pinned vertices. Energy and gradient by
§3.9. Three methods: L-BFGS with the Lewis–Overton weak Wolfe line search (the paper's best,
Sec. 4.3; default), gradient descent with a line search, fixed step.

### 4.1 Units (D17)

`ℓ̄` = mean edge length of the mesh at the start of each `optimize` call (reported). Positions
are never rescaled. A step is a displacement field `Δf = t d`; "size" of a direction is
`‖d‖_∞ = max_v |d_v|` (Euclidean per vertex). Dimensionless options:

* `initial_step` (default `0.1`): whenever the direction is steepest descent (first iteration,
  after a memory reset, gradient descent), `d = −(initial_step · ℓ̄ / ‖g‖_∞) g`, so the unit trial
  `t = 1` moves the fastest vertex by `0.1 ℓ̄`.
* `max_displacement` (default `0.5`): trial steps are capped at `t_max = max_displacement · ℓ̄ /
  ‖d‖_∞` (D19; not in the paper: the energy is scale invariant and bounded, so the doubling phase
  of a line search could otherwise accept an absurd step).
* `gradient_tolerance` (default `1e-8`): stop when `ℓ̄ · ‖g‖_∞ ≤ gradient_tolerance`
  (`ℓ̄ ‖g‖` is dimensionless).
* `fixed step` `σ` (§4.5) in units of `ℓ̄²`.

### 4.2 L-BFGS (Nocedal & Wright Alg. 7.4, 7.5) [formula: dev-lbfgs-two-loop]

Memory `m` (option `memory`, default 10). Pairs `s_k = x_{k+1} − x_k`, `y_k = g_{k+1} − g_k` (over
the free coordinates). Direction `d = −H g` by the two-loop recursion with
`H_0 = (s_kᵀy_k / y_kᵀy_k) I` for the newest pair. A pair is stored only if `s_kᵀ y_k > 0`
(guaranteed by a weak Wolfe step, [formula: dev-weak-wolfe-curvature]:
`y_kᵀ s_k = t (g_{k+1} − g_k)ᵀ d > t (c2 − 1) g_kᵀ d > 0`; can fail at a capped step or by
rounding, then the pair is skipped and counted). If `gᵀd ≥ 0` (not a descent direction) the memory
is reset and steepest descent (§4.1) used; counted. Memory is reset after every topology change
(subdivision, remeshing) and at the start of each round.

### 4.3 Lewis–Overton weak Wolfe line search (L&O Algorithm 2.6) [formula: dev-lewis-overton-line-search]

`h(t) = E(x + t d) − E(x)`, `s = gᵀd < 0`. A trial `t` is accepted when

```
S(t):  h(t) < c1 s t                         (Armijo; strict)
C(t):  g(x + t d)ᵀ d > c2 s                  (weak Wolfe; differentiability not checked, L&O Sec. 5.1)
```

Defaults `c1 = 0`, `c2 = 0.5` (L&O Sec. 5.1: "we set c1 = 0 … c2 = 1/2"); options, with
`0 ≤ c1 < c2 < 1` enforced (error otherwise). Algorithm:

```
α ← 0, β ← +∞, t ← min(1, t_max)
loop:
    evaluate E, g at x + t d       (undefined energy ⇒ S(t) fails; counted in undefined_trials)
    if S(t) fails:          β ← t
    else if C(t) fails:     α ← t
    else:                   accept t, stop
    if β < ∞:               t ← (α + β)/2            (bisection, at most max_bisections = 30)
    else if α ≥ t_max:      accept α, stop           (capped step, counted)
    else:                   t ← min(2α, t_max)        (doubling, at most max_doublings = 30)
limits exceeded: if α > 0, take α (it satisfies S; counted as a weak step, no BFGS pair if
sᵀy ≤ 0); otherwise the line search failed.
```

Termination for semi-algebraic `E` and the existence of weak Wolfe steps are L&O Thm. 2.5/2.7
(cited, [formula: dev-lewis-overton-termination]). On a failed line search L-BFGS resets its memory
and retries once from steepest descent; a second failure ends the run with status
`LineSearchFailed` (L&O terminate at "breakdown", Sec. 5.1).

### 4.4 Gradient descent

`d = −(initial_step ℓ̄ / ‖g‖_∞) g` each iteration, line search one of:

* `WeakWolfe` (§4.3) — the paper's "Armijo–Wolfe line search";
* `Backtracking { c1 = 1e-4, shrink = 0.5, max_shrinks = 50 }` (N&W Alg. 3.1, Armijo only,
  [formula: dev-armijo-backtracking], cited): start at `t = min(2 t_prev, t_max)` (`t_prev` = last
  accepted step in units of the current `d`, initially 1), shrink until
  `h(t) ≤ c1 s t`; failure ends the run (`LineSearchFailed`).

### 4.5 Fixed step

`x ← x − σ ℓ̄² g` (option `step = σ`, default `1e-3`), no line search, no decrease guarantee;
energy increases are counted in the report.

### 4.6 Stopping (D27)

Checked after every accepted step, first match wins; statuses:
`GradientTolerance` (§4.1), `EnergyStalled` (`E_{k} − E_{k+1} ≤ energy_tolerance · max(E_k, 1e-300)`
for `stall_iterations` consecutive iterations; defaults `1e-10`, `10`), `MaxIterations` (default
1000 per round), `LineSearchFailed`, `Cancelled` (the progress callback returned `Break`),
`NonDescent` (two consecutive non-descent steepest-descent directions — only possible through
rounding). The paper subdivides "once the norm of the energy gradient is below a given threshold"
(Sec. 4.2): that is `GradientTolerance` or `EnergyStalled`.

### 4.7 Report

`OptimizeReport { status, iterations, evaluations (energy+gradient), line_search_trials,
undefined_trials, capped_steps, weak_steps, skipped_pairs, memory_resets, energy_initial,
energy_final, energy_history: Vec<f64> (one per accepted iterate, including the start),
gradient_norm_history (ℓ̄‖g‖_∞), mean_edge_length (ℓ̄), max_vertex_displacement (total, in ℓ̄),
energy_increases (fixed step), near_ties (last evaluation), seconds }`.

### 4.8 Pinned vertices (D18)

Option `pinned: Vec<u32>` (not in the paper; for the Blender selection): their gradient entries are
zero (exact projection for equality constraints); out-of-range ids are an error.

### 4.9 Admissibility (D31)

`OptimizeOptions::min_face_angle` (default `1e-4` rad = 0.006°, `0` = off; `[0, π/3)`): a point is
*admissible* iff every face has its smallest corner angle `≥ min_face_angle`. The mesh
objective checks it before the energy (O(|F|), `sin θ_min = |n| / (product of the two longer
edges)`, no trigonometry): an inadmissible point is treated exactly like an undefined energy
(§3.6, D12): `+∞` inside a line search (§4.3), counted in `OptimizeReport::inadmissible_trials`
(`undefined_trials` stays for D12), and `DevError::NearDegenerateFace { face, min_angle, limit }`
at the starting point (an input below the floor is reported, never repaired). The starting point
is checked with a relative slack of `1e-6` on `sin θ` (a subdivision or a remeshing of an accepted
point may change the last digits); trial points are checked exactly. The generic [`minimize`]
ignores the option. *Why:* every energy of §3 is scale free and bounded, so squeezing a face to
zero area costs (almost) nothing; the weak Wolfe line search accepts such steps, the energy
decreases a little, and the face collapses until the cross product underflows (D31).

## 5. Flow driver, subdivision, remeshing

### 5.1 Flow (Sec. 4.2, 5.2)

```
flow(mesh, positions, FlowOptions { rounds: Vec<Round>, optimize, remesh: Option<RemeshOptions>,
     pinned, threads }, progress) -> FlowResult
Round { energy: EnergyOptions, max_iterations, subdivide_before: bool }
```

For each round: (subdivide if asked; pinned ids follow §5.2) → (remesh if enabled) → optimize, with
the remeshing check of §5.3 after every accepted iteration. The paper uses `E_P` in a coarse phase
and `E_λ` in the fine phase for Figs. 23–24; rounds allow that. `FlowOptions::uniform(energy,
subdivisions, iterations)` builds `subdivisions + 1` equal rounds. Result: final `TriMesh`,
positions, `vertex_source: Vec<u32>` (input vertex id or `NONE` for created vertices),
`face_parent: Vec<u32>` (for each final face the input face it descends from: the four children of a
subdivided face have its parent; a collapse removes two faces and the others keep theirs; a flip
rewrites the two faces of the quad in place, the new face `(a, q, p)` keeps the parent of the
flagged cap triangle `f` and `(q, b, p)` that of its neighbour `g`), and
`FlowReport { rounds: [ { optimize: OptimizeReport, remesh: Vec<RemeshReport>, n_vertices, n_faces,
energy_after_subdivision } ], cancelled, stopped: Option<FlowStop>, min_angle, bad_faces, seconds }`.
`min_angle` is the smallest corner angle of the result (`≥ min_face_angle`, §4.9) and `bad_faces`
the number of faces below the remeshing's small angle that remain (the ones remeshing could not
fix; the flow says so rather than hiding them). **Stop (D32):** if the energy of a round is
undefined at its starting state (after the subdivision, §3.6: the midpoint of a flat fold), the
flow does not run that round and the later ones, returns the state of the previous round, and sets
`stopped = Some(FlowStop { round, after_subdivision, vertex, reason })`; for round 0 it is an error.
A remeshing whose result is undefined or has a face below the floor is computed but not applied
(`RemeshReport::rejected = Some(reason)`).

### 5.2 4-1 subdivision [formula: dev-midpoint-subdivision]

Positions only, no smoothing (the paper's inset: regular 4-1 split; D20). New vertex of edge `e`
(edge ids of `TriMesh::edges()`): id `|V| + e`, position `(f_a + f_b)/2`, also on boundary edges.
Face `t = (a, b, c)` with midpoints `m_ab, m_bc, m_ca` becomes faces `4t … 4t+3`:
`(a, m_ab, m_ca)`, `(b, m_bc, m_ab)`, `(c, m_ca, m_bc)`, `(m_ab, m_bc, m_ca)` (orientation kept).
Old vertex ids are kept. Pinned: old pinned stay pinned; a midpoint is pinned iff both ends are.
**Energy is preserved** (to rounding) for all six energies: an old star keeps its normals, angles
and the direction of its vector area; a new star (interior: 6 faces, boundary: 3) consists of two
planar arcs, a hinge with energy 0. Valences: old unchanged, interior midpoints 6.
Returns `Subdivided { mesh, positions, edge_midpoint: Vec<u32>, face_parent: Vec<u32> }`.

### 5.3 Small-angle remeshing (Sec. 4.3, Fig. 18; D21)

Options `RemeshOptions { small_angle (default 1° = π/180), collapse_boundary_edges (true),
max_passes (10), check_interval (1 iteration) }`, `0 < small_angle < π/3` enforced. A triangle is
*bad* if it has a corner angle `< small_angle`; with that bound it has one or two such angles
[formula: dev-small-angle-classification]:

* **two small angles** (a cap; the third angle is `> π − 2·small_angle`): flip the edge opposite
  the large angle (its longest edge). Legal iff the edge is interior, the two opposite vertices are
  not already adjacent, each endpoint keeps ≥ 3 faces (interior) / ≥ 1 face (boundary), both new
  faces are non-degenerate and the smallest corner angle of the two faces increases.
  [formula: dev-flip-manifold] (cited)
* **one small angle** (a needle): collapse the edge opposite it (its shortest edge) `(p, q)` into
  one vertex. Legal iff the link condition holds (`link(p) ∩ link(q) = link(pq)`,
  [formula: dev-link-condition], cited), the component keeps > 4 vertices, every surviving face
  of the merged star stays non-degenerate with its normal turning by less than 90°, every
  affected interior vertex keeps ≥ 3 faces, and not both endpoints pinned. Position: the midpoint;
  the boundary endpoint if exactly one endpoint is on the boundary; the pinned endpoint if one is
  pinned. An interior edge with both endpoints on the boundary is refused (it would pinch the
  surface). A boundary edge is collapsed (to its midpoint) only if `collapse_boundary_edges`.
  The surviving vertex is the lower id; the other is removed. A collapse is also refused if a
  face of the merged star would get a smallest angle below `min_face_angle` (default `1e-4`,
  `0` = off, `< small_angle`) and below its present smallest angle (D31).

A pass collects the bad triangles, sorts them by `(smallest angle, face id)` and handles each that
is still present and still bad; passes repeat until nothing changes or `max_passes`. Vertices are
then compacted (`ddg_core::topology::compact_vertices`). The result is always a valid manifold
`TriMesh` (rebuilt and validated; a failure is a `DevError::Internal`, never a panic).
`RemeshReport { flips, collapses, refused_flips, refused_collapses, unresolved_bad_triangles,
min_angle_before, min_angle_after, energy_before, energy_after }` and `vertex_to_old: Vec<u32>`.
Inside a flow the check runs every `check_interval` accepted iterations, costs O(|F|), and only a
remesh that changed something resets the L-BFGS memory.

## 6. Cut and flatten (Sec. 4.4; D22, D23)

`cut_and_flatten(mesh, positions, &CutOptions) -> Result<DevelopableFlattening, DevError>`.

1. **Terminals.** Per-vertex energies `e_v` of `CutOptions.energy` (default `Covariance`)
   evaluated with `omit_valence_three = false` (a valence-3 triple point has positive energy unless
   flat and must be cut through), plus `extra_terminals`. Terminals: `e_v > epsilon` (default
   `1e-3`; Fig. 20 uses 0.01 … 0.001) **or** (D30, our addition) an interior vertex with
   `|Ω_v| > defect_tolerance` (`Ω_v = 2π − Σθ`, [formula: angle-defect]; default
   `Some(0.1)` rad, `None` = the paper's rule). A cut through a vertex is what lets the flattening
   absorb its angle defect isometrically; the energies are blind to a defect on a star whose
   normals are coplanar but which is not embedded (D28, D30). `energy` stays `Covariance` by default
   (D22). `star_embeddedness` (pub, also stored in the result) reports for every vertex the smallest
   cosine between a face normal of its star and `N_v`: `> 0` means the star is a graph over its
   tangent plane (embedded), `<= 0` or `NaN` flags a spike / needle / double cover / fold.
2. **Edge weights** [formula: dev-cut-edge-weight] (cited: design choice of Sec. 4.4): with `x_v` the
   unit eigenvector of the smallest eigenvalue of `A_v` (§3.4, for every vertex, whatever the
   energy) and `ê` the unit edge direction,
   `w_e = c_len · ℓ_e/ℓ̄ + c_align · (1 − ½(|ê·x_a| + |ê·x_b|))`, defaults `c_len = 1`,
   `c_align = 0` (open question 2). Non-negative by construction.
3. **Cut graph**: `ddg_bff::plan_cut(&CutPlanInput { mesh, edge_lengths: &w, cone_vertices:
   terminals, may_cut, must_cut })` — Mehlhorn's Steiner-tree approximation through the terminals
   plus tree–cotree handle loops (Erickson–Whittlesey greedy loops) and cuts between boundary
   loops. This replaces the paper's Erickson & Har-Peled (2004) cut (both are approximations of
   the NP-hard shortest cut; for genus 0 both reduce to a tree through the terminals, the paper's
   "minimal spanning tree"). A closed genus-0 component with fewer than 2 terminals gets the
   highest-energy vertices of that component added (ties: lower id) until it has 2 (`plan_cut`
   refuses otherwise); reported. `may_cut` (polygon edges for polygon input, §7.2) and `must_cut`
   (user seams) pass through.
4. **Cut**: `ddg_bff::cut::cut_mesh(mesh, &plan.cut_edges)`; charts = components
   (`TriMesh::split_components`); each must be a disk (guaranteed by `plan_cut`).
5. **Flatten** each chart with BFF and zero boundary scale factors ("isometry along the
   boundary"): `DiskBff::new(&chart, &Geometry::from_positions(..), DiskOptions)` then
   `flatten_with_scale_factors(&zeros, extension)`, `extension` default `Extension::Harmonic`
   (reproduces the boundary curve, i.e. boundary lengths, exactly; `Conjugate` is BFF's automatic
   mode). Charts are not packed (the binding packs, §7.2).
6. **Distortion**: per face `u_f = ½ ln(A_uv / A_f)` (`ddg_bff::quality::log_scale_factor`,
   [formula: face-log-scale-factor], existing), per vertex the area-weighted mean over its faces.

Output `DevelopableFlattening { cut_edges: Vec<u32> (edge ids, ascending), cut_vertex_pairs,
corner_uvs: Vec<[f64;2]> (3|F|), face_chart: Vec<u32>, face_log_scale: Vec<f64>,
vertex_log_scale: Vec<f64>, terminals: Vec<u32>, added_terminals: Vec<u32>, plan: CutPlan,
charts: Vec<ChartReport> (BFF's: flips, overlap, distortion), vertex_energy: Vec<f64>,
vertex_defect: Vec<f64> (Ω_v), star_embeddedness: Vec<f64> }`.
**Out of scope** (D23): gluing small "cracks" and re-flattening, automatic avoidance of global
overlap (the paper does both by hand); overlap is reported per chart.

## 7. Public API

### 7.1 Rust (`ddg_developability`)

```rust
pub enum Energy { Combinatorial, CombinatorialMax, Covariance, CovarianceMax, Intrinsic, IntrinsicMax }
pub struct EnergyOptions { pub energy: Energy, pub include_boundary: bool, pub omit_valence_three: bool,
                           pub threads: Option<usize> }                      // Default: Covariance, true, true, None
pub struct Evaluation { pub energy: f64, pub vertex_energy: Vec<f64>, pub counted: Vec<bool>,
                        pub gradient: Vec<[f64; 3]>, pub vertex_margin: Vec<f64>, pub near_ties: usize }
pub fn evaluate(mesh: &TriMesh, positions: &[[f64; 3]], opts: &EnergyOptions) -> Result<Evaluation, DevError>;
pub fn evaluate_energy(mesh, positions, opts) -> Result<(f64, Vec<f64>), DevError>;   // no gradient

pub enum Method { Lbfgs { memory: usize }, GradientDescent, FixedStep { step: f64 } }
pub enum LineSearch { WeakWolfe { c1: f64, c2: f64, max_doublings: u32, max_bisections: u32 },
                      Backtracking { c1: f64, shrink: f64, max_shrinks: u32 } }
pub struct OptimizeOptions { pub method: Method, pub line_search: LineSearch, pub max_iterations: usize,
    pub gradient_tolerance: f64, pub energy_tolerance: f64, pub stall_iterations: usize,
    pub initial_step: f64, pub max_displacement: f64, pub pinned: Vec<u32>,
    pub min_face_angle: f64 }                                                // 1e-4, §4.9
pub enum Control { Continue, Break }
pub struct Progress<'a> { pub round: usize, pub iteration: usize, pub energy: f64, pub gradient_norm: f64,
                          pub positions: &'a [[f64; 3]] }
pub fn optimize(mesh: &TriMesh, positions: &[[f64; 3]], energy: &EnergyOptions, opts: &OptimizeOptions,
                progress: Option<&mut dyn FnMut(&Progress) -> Control>) -> Result<Optimized, DevError>;
    // Optimized { positions: Vec<[f64;3]>, report: OptimizeReport }
pub fn subdivide(mesh: &TriMesh, positions: &[[f64; 3]]) -> Result<Subdivided, DevError>;
pub fn remesh(mesh: &TriMesh, positions: &[[f64; 3]], opts: &RemeshOptions, pinned: &[u32]) -> Result<Remeshed, DevError>;
pub fn flow(mesh: &TriMesh, positions: &[[f64; 3]], opts: &FlowOptions,
            progress: Option<&mut dyn FnMut(&Progress) -> Control>) -> Result<FlowResult, DevError>;
pub fn cut_and_flatten(mesh: &TriMesh, positions: &[[f64; 3]], opts: &CutOptions) -> Result<DevelopableFlattening, DevError>;
pub mod eig;   // sym_eig3, sym_eig2 (exists)

pub enum DevError {
    Mesh(MeshError), Geometry(GeometryError), Bff(BffError), Eig(EigError),
    LengthMismatch { what: &'static str, expected: usize, got: usize },
    NonFinitePosition { vertex: u32 }, DegenerateFace { face: u32 },
    NearDegenerateFace { face: u32, min_angle: f64, limit: f64 },            // §4.9
    UndefinedEnergy { vertex: u32, reason: &'static str },
    InvalidOption { name: &'static str, value: f64, requirement: &'static str },
    VertexOutOfRange { what: &'static str, vertex: u32, n_vertices: usize },
    Internal(String),
}
```

Nothing mutates its input; every function returns new arrays. Option structs implement `Default`
with the values of this file and `validate()`.

### 7.2 Python (`crates/ddg-py`, module `ddg_tools`) and Blender (D24)

Input as for the unwrap: `positions`, `corner_vertices`, `polygon_offsets` and Blender's
`loop_triangles` (`MeshInput` conventions of `ddg-uv`: the caller's triangulation is used, never one
chosen by us).

* `ddg_tools.developability(...) -> {vertex_energy, energy, counted}` (for a vertex-colour preview).
* `ddg_tools.make_developable(..., options) -> {positions, faces | None, vertex_source, report}`:
  runs `flow` on the triangle mesh of the loop triangles. `faces` is `None` iff the topology did not
  change (no subdivision, no remesh operation applied); then only positions change.
* `ddg_tools.developable_unwrap(..., options) -> {corner_uvs (per polygon corner), seam_edges (vertex
  pairs), polygon_chart, corner_log_scale, report}`: `cut_and_flatten` with `may_cut` = polygon edges
  (diagonals of the triangulation are never cut, so all triangles of a polygon share one chart and
  the UV of a polygon corner is well defined), then packing with `ddg_py::pack`.

**Rule for polygons and topology (D24, to confirm, open question 1).** The library never edits the
caller's mesh. Blender operator "Make Developable" (undoable; options energy, method, iterations,
subdivision rounds, remeshing on/off, small angle, pin selected vertices):
(a) topology unchanged → writes positions into the existing mesh, polygons untouched; the
developability holds for the triangulation used, so the operator reports the largest
non-planarity of non-triangle polygons; (b) topology changed (subdivision or an applied remesh
operation), or the option "Triangulate" is on → replaces the mesh data by the returned triangle
mesh (material index per face from the parent polygon; other attributes and UVs dropped, with a
warning). Operator "Developable Seams + UV": marks the developable cut (`cut_and_flatten`'s cut
edges, restricted to polygon edges) as seams, existing seams kept, and writes the UV map with the
`ddg-uv` pipeline of the Conformal Unwrap operator with "Cut Marked Seams" on (cones, extra cuts,
packing and the overlap report come from it; D23); it does not move vertices. The paper's plain
per-chart flatten stays available in Rust and Python.

## 8. Differences from the paper as printed, and decisions

* **D1** Eq. 1, second line `Σ_p (1/n_p²) Σ_{σ1,σ2∈F_p} |N_σ1 − N_σ2|²` equals the first line only if the
  pairs are unordered; over ordered pairs it is twice the first line
  (`Σ_{ordered} = 2n Σ|N − N̄|²`; *checked*: 0.1346 vs 0.2691). Ours: the first line, pairs unordered.
* **D2** Eq. 8 prints `∇_{f_i} N = (1/𝒜)((f_k − f_j) × N) Nᵀ` with `𝒜` the area (Eq. 7 uses the
  same `𝒜` with `½`); correct is `1/(2𝒜)` (*checked*, FD error 1.2 vs 6e-11). The B.4 formula
  uses `2A_M`, consistent with ours.
* **D3** Eq. 9 prints `N × (f_i − f_j)/|f_i − f_j|`; correct is `/|f_i − f_j|²` (gradient of an
  angle has units 1/length; *checked*).
* **D4** App. B.2 says ties give "a subgradient of the piecewise smooth energy `E^λ`"; meant is `E^P`.
* **D5** App. B.4: the subgradient omits the sign of antipodal sites (`−N_σ`): ours multiplies the
  term of site `q` by `s_q` and fixes the orientation `v = n̂`, `e_a = c − b, e_b = a − c, e_c = b − a`
  (*checked*: paper's version relative error ≈ 2 on stars whose active triple has an antipode).
  "Spherical centroid" is the spherical circumcenter; "`A_M` and `A_p` are the areas of triangles with
  normals `M` and `p`" means `A_M`, `A_σ`. Also "the triple of triangles whose normals define `v`":
  the sites, with sign.
* **D6** Eq. 11, second bracket: `μ_i μ_fᵀ` must be `μ_i μ_iᵀ` (*checked*); the symbol `v_i^jk`
  there is the binormal `n0 ×̂ N`, not the log direction defined just above. We use the
  equivalent smooth adjoint of §3.6, defined at `φ = 0`.
* **D7** App. B.5.1: the altitude gradient as printed is the negated transpose of the Jacobian of
  `u = n0 × w/|w|` (*checked*); the `±` enumeration of pairs is a superset of the needed hull-edge
  normals (same minimum, 4× cost); we enumerate pairs only.
* **D8** Domain: boundary vertices included (paper sums over all `V`; option); valence-3 = interior
  vertex with 3 faces, omitted by default (Sec. 4.1.4; option).
* **D9** Bipartitions: both parts non-empty; enumeration order and counts of §3.2.
* **D10** Ties and repeated eigenvalues: deterministic active piece (§3.8); near ties reported.
* **D11** `λ_max` candidates: face triples with 4 sign patterns + poles of pairs; only exactly
  degenerate triples skipped (paper: "omit redundant sites"); no pair midpoints (proof §3.5);
  all-parallel stars give 0.
* **D12** Intrinsic energies undefined for `S_v = 0` or `π − φ < 1e-6`: error, `+∞` in line search.
* **D13** Angles by `atan2` from positions, normals from the cross product, vertex normal from the
  vector area (all consistent with the derivatives of §2).
* **D14** Rounding can make `λ` slightly negative; not clamped.
* **D15** L-BFGS memory default 10 (paper silent; L&O study full BFGS).
* **D16** Line search defaults `c1 = 0`, `c2 = 0.5` and the doubling/bisection limits 30/30 (L&O
  Sec. 5.1 use "a (large) limit"); retry once from steepest descent after a failure.
* **D17** Units: step sizes and tolerances in units of the mean edge length `ℓ̄` (positions not
  rescaled); initial trial moves the fastest vertex by `0.1 ℓ̄`.
* **D18** Pinned vertices (extension).
* **D19** Step cap `max_displacement = 0.5 ℓ̄` (extension: guards the doubling phase).
* **D20** Subdivision: midpoint 4-1, no smoothing, boundary edges split too; energy preserved.
* **D21** Remeshing thresholds and legality rules of §5.3 (paper gives neither); default 1°;
  checked every iteration; memory reset only after a change.
* **D22** Cut: terminals from energies without the valence-3 omission; `plan_cut` instead of
  Erickson–Har-Peled; terminals added on closed genus-0 components with < 2; weight formula of §6
  (the paper names the ingredients but not the combination or the sign). The energy that selects
  the terminals stays `Covariance` by default (not the intrinsic width, see D30 for why).
* **D23** Flatten: per chart BFF with `u = 0` on the boundary, harmonic extension; crack gluing and
  overlap avoidance out of scope for `cut_and_flatten` (the paper's method). Its charts can overlap
  themselves globally (spot after the default flow: one chart, 48 % self-overlap, no flipped
  triangle). The user decided (2026-10-19) that the UV map of the bindings and Blender comes from
  the `ddg-uv` atlas pipeline instead, with the developable cut as marked seams (§7.2).
* **D24** Polygons and topology rule of §7.2 (confirmed by the user, §11).
* **D25** Evaluation is deterministic and bitwise independent of the thread count (§3.9).
* **D26** Gradient descent: either line search; fixed step in units of `ℓ̄²`.
* **D27** Stopping rules and statuses of §4.6 (L&O run "until breakdown" or a limit).
* **D28** Embeddedness of stars is not checked (Theorem 3.3 assumes it); energies vanish on
  non-embedded coplanar stars too.
* **D29** Eigensolver: cyclic Jacobi (`eig.rs`), not Kopp's hybrid.
* **D30** Cut terminals: besides `e_v > ε`, every interior vertex with `|Ω_v| > 0.1` rad is a
  terminal (`CutOptions::defect_tolerance`, `None` restores the paper's rule).
  *Paper:* terminals are the vertices with energy above ε (Sec. 4.4, Figs. 20–22 use the intrinsic
  energy, which the paper says remedies spikes, Sec. 5.2). *Ours:* with `E_λ` the flow leaves
  stars that are not embedded (D28: spike, needle, double cover, Fig. 9) whose normals are
  coplanar, hence energy about 0, but whose angle defect is large; they never become terminals and
  BFF with zero boundary scale must absorb their defect as area distortion. *Measured*
  (`NOTES-flow.md`, "Cut study", `examples/cut_study.sh`): spot control mesh + 2 subdivisions,
  `E_λ`, 300 it/round, ε = 1e-3: non-terminals carry Σ|Ω| = 138.2, of which 134.9 on
  non-embedded stars (86 of the 87 vertices with |Ω| > 0.1); 99 of 5112 UV triangles flipped
  (1.9 %), mean |u| 0.90. Over-run (3 x 3000 it): 222.2 of 222.2 on non-embedded stars, no
  terminal above ε, 195 flipped (5.0 %), mean |u| 6.1. With `defect_tolerance = 0.1`: 0 flipped,
  mean |u| 4.0e-3 (1.3e-3 over-run), cut length 4.1 % of the total edge length (paper rule: 2.1 %,
  and 0.1 % when over-run). Over the three control meshes and both budgets, 0 flipped in all 6
  cases (paper rule: up to 5.7 %), at the price of 1.1 (bob) to 2 (spot) to 5 (blub, where the paper rule
  found 2 terminals) times the cut length of the paper rule on the 300-iteration runs; on the dense meshes (23k to 57k faces) the criterion changes nothing
  (0 flips either way, cut length +0.01 %). The intrinsic energy for terminals alone
  (option (a)) is not enough (up to 0.5 % flipped for ε ≤ 1e-2 on E_λ-flowed control meshes, 5.5 % at ε = 1e-2 on intrinsic-flowed ones) and
  cuts more (4.2 % against 3.4 % on E_λ-flowed control meshes) than the defect criterion, and it
  fails with an `DevError::UndefinedEnergy` error when a star has an exactly antipodal normal (seen on spiky
  over-run meshes); the defect criterion is energy independent, cheap and cannot fail. 0.1 rad
  against 0.02: 0.02 lowers mean |u| by 2 to 4 times for +20 to 25 % cut length.
* **D31** Admissibility floor on the smallest corner angle (§4.9), default `1e-4` rad; also a collapse
  rule of the remeshing. *Paper:* silent (it removes faces below a small angle by remeshing, but a
  face that cannot be fixed, because the flip or collapse is illegal, is left to the optimizer).
  *Found (NOTES-flow.md, "Robustness study"):* `E_P` and the intrinsic energies squeeze faces to
  angles below `3e-9`: bare `optimize` without remeshing does it in 40 of 48 small test runs
  (perturbed icospheres, tori, bumpy grids, 6 seeds, `E_P` and intrinsic; 800 iterations), and the flow with remeshing leaves
  such faces on the blub, spot and bob control meshes (blub `E_P`, 300 and 3000 it/round; spot 3000)
  because the unfixable faces (collapse refused by the normal-turn or link condition, flip refused)
  are counted once as "known unfixable" and then the optimizer is free to squeeze them
  (`Geometry::from_positions` of the cut then rejects them: the triangle inequality fails by
  rounding). *Decision:* wall at `1e-4` as an undefined energy. *Measured on 24 of those runs
  with floors 1e-5, 1e-4, 1e-3:* the final energy equals the unconstrained one to within 1 % in
  most cases and is lower in several (the squeezed faces are worth nothing), the lowest angle is
  exactly the floor.
* **D32** Flows never end in the energy being undefined at a round start, and do not walk into it.
  *Found:* "undefined at vertex 881" in the blub intrinsic 3000 it/round flow is raised at the
  start of round 1 by the energy evaluation of the subdivided mesh. The state at the end of round 0
  is a flat pancake: pairs of big faces of nearly opposite normals along an edge (cosine of the
  normals `-0.9999995`, areas 0.45 and 0.41); the midpoint of such an edge has a vector area
  nearly zero and a face normal within `1e-6` of the antipode of the vertex normal (D12). The
  flow with `E_λ` reaches the same pancakes (energy `1e-11`, fold gap `4e-10`). *Rejected:* a
  second wall on the fold angle (`min_fold_angle` 1e-3, 1e-4, 1e-5, 1e-6; the same measurements):
  a hard wall in L-BFGS stops the whole run when one edge reaches it (`LineSearchFailed` at energy
  0.9 instead of 5e-8 for spot `E_λ`; 1.4 against 2.7e-11 for blub), and with 1e-6 the undefined
  states remained on spot intrinsic-max. *Decision:* the flow evaluates the energy of the next
  round on the candidate state before committing; if it is undefined, it stops before that round
  with `FlowReport::stopped` (round 0: error). A remeshing result that is undefined or below the
  floor is not applied (`rejected`). Nothing is repaired.
* **D33** A boundary star with at most two faces contributes exactly nothing (value 0, no
  sensitivities) for every energy (the intrinsic ones are still undefined where D12 says so) (§3.1 says its energy vanishes identically;
  the eigenvector of the rank-deficient matrix is arbitrary and, for nearly coplanar normals, fed
  noise of the size of the angle into the gradient; found by the oracle fuzz).

## 9. Formula ids

Rust functions carry `[formula: <id>]`; each id is a Lean theorem in
`lean/proofs/DdgProofs/Developability/*.lean` ("proved") or a heading in `lean/proofs/CITED.md`
("cited"). `scripts/check-formulas.sh` fails on a Rust id with neither. Planned Rust locations are
modules of this crate.

| id | statement | Rust (planned) | Lean / CITED |
|---|---|---|---|
| `dev-face-frame` | `N = n/|n|`, `A = |n|/2`, `n` cyclic form | `geometry::face_frame` | proved, `Basic.lean` |
| `dev-corner-angle` | `atan2(|a×b|, a·b)` is the angle of `a, b` | `geometry::corner_angle` | proved, `Basic.lean` |
| `dev-area-gradient` | Eq. 7 | `geometry::backprop_face` | proved, `Basic.lean` |
| `dev-normal-gradient` | Eq. 8 with `1/(2A)` | `geometry::backprop_face` | proved, `Basic.lean` |
| `dev-angle-gradient` | Eq. 9 with `|e|²` | `geometry::backprop_face` | proved, `Basic.lean` |
| `dev-backprop-face` | adjoint form of §2 equals the chain rule | `geometry::backprop_face` | proved, `Basic.lean` |
| `dev-vertex-normal` | vector area; independent of `f_v` at interior `v` | `geometry::vertex_normal` | proved, `Basic.lean` |
| `dev-gradient-invariance` | invariance ⇒ `Σg = 0`, `Σ f×g = 0`, `Σ f·g = 0` | `evaluate` (doc) | proved, `Invariance.lean` |
| `dev-hinge-coplanar` | Theorem 3.3 | — (docs) | cited (paper's proof; embeddedness not in Mathlib) |
| `dev-hinge-antiparallel` | Prop. A.1 | — | cited |
| `dev-hinge-ruled` | Prop. A.2 | — | cited |
| `dev-valence-three-flat` | Prop. A.3 (justifies the omission) | `energy::domain` | cited |
| `dev-star-bipartitions` | enumeration is exactly the edge-connected bipartitions; counts `k(k−1)/2`, `k−1` | `energy::bipartitions` | proved, `Combinatorial.lean` |
| `dev-combinatorial-width` | Eq. 1, three forms equal (unordered pairs) | `energy::combinatorial` | proved, `Combinatorial.lean` |
| `dev-combinatorial-gradient` | `∂π/∂N_σ = (2/n_p)(N_σ − N̄_p)` | `energy::combinatorial` | proved, `Combinatorial.lean` |
| `dev-combinatorial-max` | Eq. 5 and gradient of the active pair | `energy::combinatorial_max` | proved, `Combinatorial.lean` |
| `dev-covariance-matrix` | Eq. 3 = smallest eigenvalue of Eq. 4 (Rayleigh) | `energy::covariance` | proved, `Covariance.lean` |
| `dev-eigenvalue-supergradient` | `λ_min(B) ≤ λ_min(A) + xᵀ(B−A)x` | `energy::covariance` | proved, `Covariance.lean` |
| `dev-simple-eigenvalue-derivative` | `∇_A λ = x xᵀ` at a simple eigenvalue | `energy::covariance` | cited (Magnus 1985; Kato 1995 II §5; Mathlib has no eigenvalue perturbation) |
| `dev-covariance-gradient` | Eq. 10 = chain rule of `xᵀ dA x` | `energy::covariance` | proved, `Covariance.lean` |
| `dev-angle-weight-retessellation` | splitting a corner with equal normal keeps `A_v` | `energy::covariance` (doc) | proved, `Covariance.lean` |
| `dev-max-covariance-voronoi` | `min φ` attained at a candidate of §3.5 | `energy::covariance_max` | cited (proof sketch §3.5; spherical Voronoi not in Mathlib) |
| `dev-spherical-circumcenter` | `v·a = v·b = v·c = det(a,b,c)/|n|` | `energy::covariance_max` | proved, `MaxCovariance.lean` |
| `dev-max-covariance-gradient` | §3.5 sensitivities = derivative of `g²` | `energy::covariance_max` | proved, `MaxCovariance.lean` |
| `dev-log-map` | `exp_{n0}(Ñ) = N`, `|Ñ| = φ`, `Ñ ⊥ n0` | `energy::log_map` | proved, `Intrinsic.lean` |
| `dev-log-map-gradient` | smooth adjoint of §3.6 (= corrected Eq. 11) | `energy::log_map_adjoint` | proved, `Intrinsic.lean` |
| `dev-tangent-eigen-derivative` | basis terms vanish: `dλ = xᵀ dÃ x` | `energy::intrinsic` | proved, `Intrinsic.lean` |
| `dev-intrinsic-width` | `max_{σ,τ}⟨Ñ_σ−Ñ_τ,u⟩² = (max − min)²` | `energy::intrinsic_max` | proved, `Width.lean` |
| `dev-min-width-hull-edge` | min width ⊥ a hull edge | `energy::intrinsic_max` | cited (Houle & Toussaint 1988; Toussaint 1983) |
| `dev-width-altitude-gradient` | `du = −(u·dw) ŵ/|w|` | `energy::intrinsic_max` | proved, `Width.lean` |
| `dev-min-width-gradient` | sensitivities of §3.7 | `energy::intrinsic_max` | proved, `Width.lean` |
| `dev-jacobi-rotation` | Rutishauser rotation zeroes `a_pq`, keeps trace and Frobenius norm; one rotation diagonalizes 2x2 | `eig::sym_eig3`, `eig::sym_eig2` | proved, `Eig.lean` |
| `dev-lbfgs-two-loop` | two-loop recursion = `H_k g`, `H_k` SPD if all `sᵀy > 0` | `optimize::lbfgs_direction` | cited (N&W Alg. 7.4, Sec. 7.2) |
| `dev-weak-wolfe-curvature` | weak Wolfe ⇒ `sᵀy > 0` | `optimize::weak_wolfe` | proved, `Optimizer.lean` |
| `dev-lewis-overton-line-search` | the accepted `t` satisfies `S` and `C` | `optimize::weak_wolfe` | proved, `Optimizer.lean` |
| `dev-lewis-overton-termination` | existence of weak Wolfe steps, termination (L&O Thm. 2.5, 2.7) | `optimize::weak_wolfe` (doc) | cited |
| `dev-armijo-backtracking` | finite termination for a descent direction of a C¹ function | `optimize::backtracking` | cited (N&W Lemma 3.1 / Alg. 3.1) |
| `dev-midpoint-subdivision` | children coplanar with parent, old angles kept, new stars hinges | `subdivide` | proved, `Subdivision.lean` |
| `dev-small-angle-classification` | `θ_s < π/3` ⇒ ≤ 2 small angles; two small ⇒ the third `> π − 2θ_s` and is opposite the longest edge | `remesh::classify` | proved, `Remesh.lean` |
| `dev-flip-manifold` | flip of an interior edge with non-adjacent opposite vertices keeps a manifold | `remesh::flip` | cited (Edelsbrunner 2001, Ch. 2) |
| `dev-link-condition` | collapse keeps a 2-manifold iff the link condition holds | `remesh::collapse` | cited (Dey, Edelsbrunner, Guha, Nekhayev 1999) |
| `dev-cut-edge-weight` | §6 step 2 | `cut::edge_weights` | cited (paper Sec. 4.4, design choice) |
| `dev-regular-cone-values` | closed forms of §10.3 for `E_λ`, intrinsic, intrinsic max | — (tests) | proved, `RegularCone.lean` (if trig sums prove heavy: oracle only, no Rust tag) |

Reused, existing: `face-log-scale-factor` (ddg-bff), `tree-cotree-disk`, `cut-leaf-pruning`
(through `plan_cut`). 43 `dev-` ids: 30 proved, 13 cited.

**Oracle** (`lean/oracle`, no Mathlib, `DdgOracle/Developability.lean`; commands
`dev [IN [OUT]]`, `dev-selftest`, `dev-fixtures write|check DIR`, text format v1 with hex floats):
task `energy` (mesh, positions, energy options → per-vertex energies, total, gradient, the chosen
partition / candidate / pair / eigen-tie flag per vertex), task `step` (one first L-BFGS
iteration: scaled steepest descent + Algorithm 2.6 with the options → trial `t` sequence, accepted
`t`, new positions), task `subdivide`. Its own eigensolver (cyclic Jacobi written from Golub & Van
Loan). Fixtures: §10.3 cones, hinges, boundary fans, random stars, the crumpled mesh, a cube with
valence-3 corners.

## 10. Verification

### 10.1 Layers (CLAUDE.md, adapted)

There is no certifiable minimizer. Layers: (1) finite differences of every energy against its
gradient; (2) exact values (zeros, closed forms); (3) invariances; (4) optimizer behaviour;
(5) the Lean oracle (differential, golden fixtures); (6) the Lean proofs; (7) end-to-end runs.

### 10.2 Finite differences

For each energy, random stars (interior valence 3–9, boundary fans), the crumpled mesh and spot's
control mesh: central differences `h = 1e-6 ℓ̄` per coordinate of every vertex of a sample, relative
error `≤ 1e-6` (observed 1e-10 to 2e-8 in the check), **only** where the active piece has a margin:
perturb only vertices of stars whose `vertex_margin` exceeds `1e-6 · max(1, e_v)` (and of no
neighbouring star with a smaller margin), and verify that the same pieces are active at `±h`.
Separate FD tests of the building blocks: Eq. 7–9, `backprop_face` with random sensitivities, log
map adjoint including `φ ∈ {1e-6, 1e-3, 1, 3}`, circumcenter derivative, altitude derivative.

### 10.3 Exact values

* Hinge stars (two half-planes, random fan) and flat stars: all six energies `≤ 1e-28` (observed
  `1e-32`). Two-face boundary stars: 0.
* **Regular cone** of valence `k ≥ 3` (apex `(0,0,h)`, ring `(cos 2πm/k, sin 2πm/k, 0)`), with
  `sin²α = h²/(h² + cos²(π/k))`, `θ = 2 arcsin(sin(π/k)/√(1+h²))` (apex angle), `α = atan2(h, cos π/k)`
  (*all checked to 1e-12 for k = 4…8, h ∈ {0.1, 0.6, 2}*):
  - `E_λ = θ k min(sin²α / 2, cos²α)` (the smaller eigenvalue is double when `sin²α/2 < cos²α`);
  - `E_P = sin²α · min_{1≤m≤k−1} [2 − ρ_m² − ρ_{k−m}²]`, `ρ_m = sin(mπ/k)/(m sin(π/k))`;
  - `λ_max = min(cos²α, sin²α cos²(π/k))` for even `k` (odd `k`: oracle values only);
  - intrinsic `= θ k α² / 2` (double eigenvalue);
  - intrinsic max `= (2α cos(π/k))²` (even `k`), `(α (1 + cos(π/k)))²` (odd `k`).
* Valence-3 interior vertex: `counted = false` and contributes 0 to `E` by default.

### 10.4 Invariances

Rigid motions and uniform scaling (factors 1e-3, 1e3): energies equal to `1e-12` relative,
gradients rotate / scale by `1/s`; `Σ g = 0`, `Σ f×g = 0`, `Σ f·g = 0` to `1e-10 ‖f‖‖g‖`.
Retessellation (split an edge and its faces at an interior point, coplanar): `E_λ` and `λ_max`
unchanged (not `E_P`; not the intrinsic energies at the edge's ends, whose vertex normals change).
4-1 subdivision: all six energies unchanged to `1e-12` relative. Determinism: 1, 2, 7, 64 threads
give bitwise equal results.

### 10.5 Optimizer

Energy history non-increasing for L-BFGS and GD (strict decrease of accepted steps with
`c1 = 0`); every accepted weak Wolfe step satisfies `S` and `C`; no stored pair with `sᵀy ≤ 0`;
BFGS with this line search on `f(x) = |x|` behaves as L&O Sec. 4 describe (iterates alternate
in sign, `|x|` decreases); L-BFGS on a smooth quadratic converges superlinearly; a run with
`max_iterations = 0` returns the input unchanged; pinned vertices never move; cancellation via the
callback returns the last accepted state.

### 10.6 Topology

Subdivision: Euler characteristic, boundary loop count and orientation preserved, `|F| × 4`.
Remeshing: output always passes `TriMesh::new`; a constructed cap triangle is flipped, a needle
collapsed; refusals for each legality rule are hit by a dedicated test; min angle never decreases
by a flip.

### 10.7 End-to-end (meshes in `data/models/`; unzip `crumpleddevelopable.zip` there)

* `CrumpledDevelopable.obj` (512 faces; the paper's Fig. 4 / 7 input, every interior vertex
  flattenable): `E_λ` flow, L-BFGS, 500 iterations: energy falls by ≥ 90 %; report mean
  `|Ω|`, max vertex displacement in `ℓ̄`, dihedral-angle histogram (expect fewer large creases).
* Spot control mesh (`spot_control_mesh.obj`, fan-triangulated, 372 faces) with 2 subdivision
  rounds, `E_P` then `E_λ` (as Figs. 23–24), remeshing on: completes, energy decreases in every
  round, then `cut_and_flatten`: every chart a disk, every terminal on the cut, `0` flipped UV
  triangles (achieved with the default `defect_tolerance = 0.1`, D30; asserted `== 0`; the paper's
  rule alone, `None`, gives 1.9 % on this run and 5.0 % when over-run, asserted as the evidence
  for D30), mean `|u_f|` ≤ 0.02 asserted (3.97e-3 measured), cut length ≤ 6 % of the edge length
  (4.1 % measured).
* `bob_controlmesh.obj` (genus 1): flow + cut: the cut contains the handle loops (2 loop edges
  counted by `CutPlan::n_loop_edges`) and every chart is a disk.
* A jittered cube triangulated so its 8 corners have valence 3 (Fig. 16): with the omission the
  corners keep `|Ω| > π/4`; without it they round off (energy at the corners falls).
* All six energies × {L-BFGS, GD+WeakWolfe, GD+Backtracking, FixedStep} on spot's control mesh for
  50 iterations: no error, monotone where guaranteed.
* Robustness (`tests/robust.rs`, D31 to D33): fast tests on perturbed icospheres (the optimizer
  squeezes faces without the floor and not with it; a start below the floor is an error; flows
  return `Geometry`-valid meshes; two-face boundary stars contribute exactly nothing; `face_parent`),
  ignored: the cut-study reproducers (blub `E_P` 300 it, blub intrinsic 3000 it) and the over-runs
  (3 x 3000 it) of all six energies on spot / bob / blub control meshes complete without error and
  without a face below the floor.
* Independent runs are parallel jobs (one per mesh), one log each in `runs_tmp/`.

### 10.8 Benchmarks (`cargo bench -p ddg-developability`, single-thread numbers run alone)

* `energy`: each of the six energies, value+gradient, on `spot_triangulated` (5856 faces),
  `bob_tri` (10688), `blub_triangulated` (14208) and spot subdivided twice (93696 faces); report
  ns per vertex, 1 thread and all threads. Targets (1 thread): `E_λ`, intrinsic `≤ 2 µs/vertex`,
  `E_P`, `π_max` `≤ 2 µs`, `λ_max`, intrinsic max `≤ 5 µs`.
* `flow`: spot control + 2 subdivisions with `E_λ` (≈ Fig. 1: 1k triangles, two subdivisions),
  crumpled, bob control; report wall time, iterations, evaluations, final energy. The paper's
  floor: "seconds to a few minutes" single-threaded on a 2018 laptop; our target is ≤ 60 s
  single-threaded for the spot run at the default iteration budget, and record speed-up with
  threads.

## 11. Questions answered by the user (2026-10-19)

1. **Polygons and topology (D24):** confirmed as written in §7.2.
2. **Cut edge weights (D22):** default weight 0 (length only); the along-ruling term stays an option.
3. **Boundary vertices (D8):** included by default; excluding them stays an option.
4. **Safeguards (D19, D21):** keep the step cap `0.5 ℓ̄` and the 1° threshold; tune by
   measurement and record the numbers if runs show stalls or poor triangles.
5. **UV map (D23):** the bindings and Blender take the developable cut as seams and unwrap with the
   `ddg-uv` atlas pipeline (reuse, not the paper's plain BFF flatten, whose charts overlap).
