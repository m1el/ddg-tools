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

**Blocked vertices.** The vertices of every face whose smallest angle is below `2 · min_face_angle`
are *blocked*: their entries of the gradient and of the search direction are zero. The blocked
set is computed at the starting point and again at every accepted point, and stays fixed during
a line search (trial points do not change it). Without blocking the line search stops at the
floor (spot control mesh, intrinsic, 300 it/round: 164 of 900 iterations at E 16.7; with
blocking the full budget, E 5.05; D31).

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
* `ddg_tools.make_developable(..., options) -> {positions, faces | None, face_polygon | None,
  vertex_source, report}`: runs `flow` on the triangle mesh of the loop triangles, all rounds in one
  library call. `faces` is `None` iff the topology did not change (no subdivision, no remesh
  operation applied); then only positions change. `face_polygon[t]` is the polygon of the caller's
  mesh that `FlowResult::face_parent[t]` (an input triangle) belongs to: exact through any number of
  rounds, flips and collapses. The report carries the library's: `stopped` (D32), `bad_faces`,
  `min_angle`, `inadmissible_trials`, and per round `flips`, `collapses`, `remesh_rejected`.
  `min_face_angle` (§4.9) is an option; it is also the remeshing's collapse floor and must be below
  `small_angle`. An input face below it is a `MeshError` naming the polygon (with remeshing on the
  flow repairs the input first, and the message names a face of the current mesh).
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
  (§12.3, `examples/cut_study.sh`): spot control mesh + 2 subdivisions,
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
  *Found (§12.4):* `E_P` and the intrinsic energies squeeze faces to
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

### 10.8 Benchmarks (`cargo bench -p ddg-developability`, single-thread numbers run alone; measured values in §12.6, §12.7)

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

## 12. Studies and implementation notes (measured data)

Moved here from the progress notes of the implementation (`NOTES-*.md`, deleted). The logs
quoted as `runs_tmp/dev-*` are not committed (the directory is gitignored); the scripts that
regenerate them are `examples/study.sh`, `cut_study.sh`/`cut_study.py`, `robust_study.sh`/
`robust_summary.py`/`robust_compare.py` and the benches. The benchmark numbers of §12.6 and §12.7
(energy, eig and flow benches) were taken alone on a quiet machine: AMD Ryzen Threadripper 3970X (32 cores, 64 threads), 2026-10-19, nothing else running (load average 0.1 to 1.5 before and during; the multithread parts load the machine themselves). The seconds
of the studies (§12.1 to §12.5) were taken in parallel runs under load (60 to 150) and remain upper
bounds.

### 12.1 How the flow, the cut and the tests behave

* `flow`: per round subdivide (pinned ids by `Subdivided::pinned_after`, new vertices get
  `vertex_source = NONE`), remesh once, then `optimize` in segments. The progress callback checks
  every `check_interval` accepted iterations for triangles with a corner angle below the small
  angle; if there are more than the count that was left after the last remesh (those are
  unfixable by the legality rules), it breaks, remeshes, and restarts the optimizer with the
  remaining budget of the round (fresh L-BFGS memory). A segment that ends with `LineSearchFailed`
  or `NonDescent` also triggers one remesh and a restart if the remesh changed something.
  At the end of the budget a pending check is still remeshed, so the returned mesh has no fixable
  bad triangle. `RoundReport.optimize` is the merge of `segments` (histories concatenated, a
  restart's first entry is the energy right after the remesh). `RemeshReport.energy_before/after`
  are filled in.
* `FlowOptions::threads = Some(n)` runs the flow in a local pool; a user progress callback (not
  `Send`) stays on the calling thread and is served over channels (positions are copied per call).
  Results are bitwise independent of the thread count (tested with 1 and 3).
* `Progress.iteration` is cumulative within the round over restarts; `optimize.pinned` must be
  empty (error), use `FlowOptions::pinned`.
* `cut_and_flatten` follows §6. Closed genus-0 components with fewer than 2 terminals get their
  highest-energy vertices added up front; a retry loop on `BffError::SphereNeedsTwoCones` covers
  closed pieces created by `must_cut`. Faces with zero plane area get `face_log_scale = NaN`
  (excluded from the vertex mean and from `mean_abs_log_scale`); flipped faces use |A_uv| and are
  counted in `charts[i].distortion.flipped_triangles`. Extra output beyond the SPEC struct:
  `mean_abs_log_scale` (area-weighted mean of |u_f|), helpers `n_charts`, `flipped_triangles`,
  `is_injective`. `CutOptions` has `length_weight` (c_len), `alignment_weight` (c_align, default 0),
  `extension`, `overlap_samples_per_triangle`, `threads` in addition to the SPEC fields.
* No new `DevError` variant was needed.

**Tests.**

`cargo test -p ddg-developability --test flow --test cut` (fast, < 2 s in release).
Long ones: `cargo test --release -p ddg-developability --test flow --test cut -- --ignored --nocapture`
(about 15 s in release). Exact/analytic checks: a flat disk, a tube and an annulus with seams unroll
isometrically (mean |u| < 1e-6, 0 flips); `c_align` makes the tube cut follow a ruling; torus cut has
2 loop edges; closed spheres get 2 added terminals; noisy cylinder and cone sheets: E_λ falls from
about 2 to 5e-7 in 600 iterations; icosphere: top 20 % of vertices hold 78 % of the energy after
400 iterations (22 % before); subdivision preserves the round-to-round energy to 1e-9 for all six
energies; remeshing in a flow removes an injected needle and keeps `vertex_source` consistent.
Models (release, `--ignored`): crumpled E_λ 61.9 -> 9.1e-4 in 500 iterations (-99.99 %), dihedral
angles above 30 degrees 354 -> 35 of 800 edges, mean |angle defect| 0 -> 0.025 (the flow does not
preserve isometry; it only makes normals coplanar per star).

### 12.2 Findings before D30

* **(Resolved by D30, see §12.3.)** **Spec expectation not met on spot.** SPEC §10.7 expects at most 0.1 % flipped UV triangles after
  flow + cut on spot. Measured (control mesh, 2 subdivisions, E_λ, ε = 1e-3): 300 iterations per
  round: 99 of 5112 (1.9 %), mean |u| 0.90; 1000 per round: 86 of 4568 (1.9 %), mean |u| 3.4.
  Lower ε helps (300 iterations: ε = 1e-4: 25 flips, 1e-5: 6 flips = 0.12 %, mean |u| 0.085, with
  1086 terminals). Cause (measured with Σ|Ω|): the flow leaves a lot of canceling positive and
  negative angle defect at vertices whose energy is below ε (spot, 300 it: Σ|Ω| 80 at terminals,
  138 at the other interior vertices). E_λ is quadratic in the normal spread while the defect is
  not, so a small energy does not mean a small defect. BFF with zero boundary scale must absorb
  that as area distortion. The tests assert <= 5 % for spot/bob/blub (blub 15 %) and print the
  numbers.
* **(Resolved by D30.)** **Over-running the flow gives spikes.** Spot, 3 rounds, 3000 iterations each: E_λ = 1e-8, no
  terminal above ε, and the cut/flatten has 195 flipped triangles and mean |u| 6. The minimizers
  include non-embedded stars (spikes, double covers; SPEC D28, paper Fig. 9). The same shows on
  blub (E < 1e-3 everywhere after 2 x 500 iterations but Σ|Ω| = 99). Stopping early matters; an
  embeddedness check would be a feature request (out of scope).
* Remeshing is active all the time: 10 to 60 collapses per round on spot/blub/bob control meshes
  (needles keep forming), flips are rare. Face counts drop 10 to 15 % after the first round on
  control meshes. No stalls and no line-search failures were seen in the study with the 1 degree
  threshold and the 0.5 step cap; a tuning was not needed.
* Crumpled developable with E_λ: `cut_and_flatten` finds no terminals after the flow (E_λ ~ 6e-4),
  one chart, 0 flips, mean |u| 0.01.

### 12.3 Cut study (D30): flipped UV triangles after flow + cut

Question: why does the paper's terminal rule (`E_λ > ε`) leave 1.9 % flipped UV triangles on spot
(SPEC §10.7 expected 0.1 %), 5 % when the flow is over-run, and what fixes it at what cut length.

**Hypothesis (non-embedded stars carry the defect): confirmed.** `star_embeddedness` (smallest
cosine between a face normal and the vertex normal; <= 0 means the star is not a graph over its
tangent plane) on spot control + 2 subdivisions, `E_λ`, ε = 1e-3, paper rule:

| run | terminals | Σ\|Ω\| at terminals | Σ\|Ω\| at non-terminals | of which on non-embedded stars | non-terminals with \|Ω\| > 0.1 (non-embedded) | flipped | mean \|u\| |
|---|---|---|---|---|---|---|---|
| 3 x 300 it | 78 | 80.4 | 138.2 | 134.9 | 87 (86) | 99 (1.9 %) | 0.90 |
| 3 x 3000 it | 0 (+2 added) | 9.4 | 222.2 | 222.2 | 77 (77) | 195 (5.0 %) | 6.08 |

Non-embedded stars are 14 % of the interior vertices (360 of 2558), yet hold 98 % of the defect of
the non-terminals; only 68 % of the flipped triangles touch one (the distortion spreads through
the chart), so the flips are not located exactly at the spikes, but the defect that causes them is.

**Study.** `examples/cut_study.sh` (flows once per mesh x flow energy x budget, written to OBJ, then
`develop --cut-only` for 3 terminal energies x ε in {1e-2, 1e-3, 1e-4} x `defect_tolerance` in
{off, 0.1, 0.02}: 1188 cut jobs, `xargs -P 56`, 3 min wall), `examples/cut_study.py` aggregates
(`runs_tmp/dev-cut-study/{summary.md,all.csv}`). Flows: spot, bob, blub control meshes and
CrumpledDevelopable (control: 2 subdivisions, 300 and 3000 it/round; crumpled 500 and 3000), seven
dense meshes (1 subdivision, 300 it/round), flow energies `E_λ`, intrinsic, `E_P`. 4 of 45 flows are
unusable for the cut or the flow failed (see "Problems found").

Control meshes, flow `E_λ`, ε = 1e-3. Cells: terminals / flipped UV triangles / mean |u| / cut length
as % of the total edge length.

| mesh, it/round | paper: E_λ, off | intrinsic, off (a) | E_λ + tol 0.1 (b, default) | intrinsic + tol 0.1 (c) | E_λ + tol 0.02 |
|---|---|---|---|---|---|
| spot 300 | 78 / 99 (1.94 %) / 0.90 / 2.06 | 309 / 6 (0.12 %) / 4.2e-2 / 4.60 | 165 / 0 / 4.0e-3 / 4.06 | 329 / 0 / 2.4e-3 / 4.76 | 301 / 0 / 1.0e-3 / 4.93 |
| spot 3000 | 2 / 195 (5.0 %) / 6.08 / 0.13 | n/a (energy undefined) | 79 / 0 / 1.3e-3 / 4.05 | n/a | 89 / 0 / 6.2e-5 / 4.13 |
| bob 300 | 255 / 0 / 4.9e-2 / 2.18 | 339 / 0 / 1.1e-2 / 2.77 | 276 / 0 / 1.7e-2 / 2.37 | 345 / 0 / 1.0e-2 / 2.80 | 352 / 0 / 6.1e-3 / 2.70 |
| bob 3000 | 2 / 145 (1.75 %) / 0.59 / 0.96 | 115 / 1 / 0.14 / 2.18 | 90 / 0 / 4.6e-3 / 2.10 | 175 / 0 / 3.3e-3 / 2.49 | 168 / 0 / 9.3e-4 / 2.56 |
| blub 300 | 2 / 140 (4.88 %) / 2.49 / 0.76 | 165 / 8 (0.28 %) / 0.10 / 5.10 | 72 / 0 / 6.2e-3 / 3.73 | 194 / 0 / 5.1e-3 / 5.54 | 129 / 0 / 8.7e-4 / 4.79 |
| blub 3000 | 2 / 158 (5.67 %) / 7.93 / 0.13 | n/a | 59 / 0 / 1.8e-3 / 4.05 | n/a | 64 / 0 / 6.3e-6 / 4.28 |

Summary over all usable flows (worst case over meshes; "cut" is the mean over meshes):

* control meshes, flow `E_λ`, 300 it: paper rule ε = 1e-2 / 1e-3 / 1e-4: worst 4.9 / 4.9 / 4.9 % flipped,
  cut 0.8 / 1.7 / 4.0 %. Tolerance 0.1: 0 flipped for every ε, worst mean |u| 3.3e-2 / 1.7e-2 / 5.8e-3,
  cut 3.1 / 3.4 / 5.3 %. Tolerance 0.02: 0 flipped, worst mean |u| 1.1e-2 / 6e-3 / 8.6e-4,
  cut 4.0 / 4.1 / 5.8 %.
* control meshes, flow `E_λ`, 3000 it: paper rule 5.7 % flipped, mean |u| up to 7.9, cut 0.4 %;
  tolerance 0.1: 0 flipped, mean |u| <= 4.6e-3, cut 3.4 %.
* control meshes, flow intrinsic or `E_P`: the paper rule is worse for coarser ε: worst flipped
  0.1 to 0.6 % at ε = 1e-4, 0.4 to 3.4 % at 1e-3, 0.5 to 5.5 % at 1e-2 (worst with `E_λ` terminals
  on the intrinsic flow); tolerance 0.1: 0 flipped in every case.
* dense meshes (seven, 23k to 57k faces, three flow energies): 0 flipped in every variant;
  tolerance 0.1 changes neither mean |u| nor the cut length (3.88 % -> 3.89 % at ε = 1e-3); mean |u|
  is set by ε alone (medians 9e-5 at ε = 1e-4, 8e-3 at 1e-3, 0.2 at 1e-2).
* CrumpledDevelopable: 0 flipped everywhere; the tolerance adds one terminal (1 cut edge, 0.17 %
  of the length) and lowers mean |u| from 1e-2 to 1.4e-4 (500 it) or 1e-6 (3000 it).
* The terminal energy matters little once the tolerance is set (`E_λ`, intrinsic and `E_P` give the
  same 0 flips and mean |u| within a factor 2); intrinsic costs 20 to 30 % more cut length.

**Decision (D30):** `CutOptions::energy` stays `Covariance`, new `defect_tolerance` defaults to
`Some(0.1)` rad (`None` = the paper's rule). Reason: it removes all flips on all control meshes and
over-runs, it is free on dense meshes, it cannot error (the intrinsic energy is undefined when a
face normal is antipodal to the vertex normal, which happens exactly on the spiky stars), and the
cut stays at 2.1 to 4.1 % of the edge length on control meshes (the paper rule is 0.1 to 2 % there
only because it leaves the defect uncut). Tolerance 0.02 is the choice for lowest distortion
(mean |u| below 1e-3 on spot and blub at +20 % cut length).

Problems found in other modules: fixed, see §12.4 (D31, D32).

### 12.4 Robustness study (D31, D32, D33)

Two failures of the cut study, both found in the flow itself. Logs: `runs_tmp/dev-robust/`
(`repro/`, `sweep/`, `sweep2/`, `over/`, `study/`); scripts `examples/robust_study.sh`,
`robust_summary.py`, `robust_compare.py`.

**1. Degenerate faces with `E_P` (and the intrinsic energies).** Mechanism, measured:

* The energies are scale free: squeezing a face to zero area costs (almost) nothing. Bare
  `optimize` (no remeshing) on small meshes (perturbed icospheres level 1 and 2, tori, bumpy grids,
  6 seeds, `E_P` and intrinsic, 800 iterations): 40 of 48 runs end with a face of angle
  1e-9 to 1e-12; the weak Wolfe search accepts every such step because the energy goes down a bit.
* In the flow the remeshing does not stop it: the blub control mesh with `E_P`, round 0, iteration
  90: one needle appears, its collapse is refused ("normal turn of the merged star"), the face is
  counted as `known_unfixable`, and from then on the flow only re-checks for *more* bad faces than
  that. The unfixable face shrinks from 1.2e-2 rad through 1e-4 (iteration 120), 1e-7 (iteration
  157) to 1e-9 (iteration 175) while the energy changes in the 7th digit. Refusal reasons over the
  three failing flows (29 + 1, 19 + 16 + 40 + 5, 4 + 3 + 8 + 18 + 8 + 8 refused operations):
  collapse "normal turn of the merged star" and "link condition", flip "normal turn", "degenerate
  new face", "smallest angle does not increase", "opposite vertices adjacent": all legal refusals,
  not bugs, but they leave faces of angles (0, 0, 180) and (0, 90, 90) in place.
* The final faces have angles like (2e-8, 180, 1.8e-8) degrees-scale; `Geometry::from_positions`
  rejects them in the cut ("violates the triangle inequality", Heron with rounding).

*Fix (D31):* admissibility floor `min_face_angle = 1e-4` rad in the optimizer (an inadmissible
trial is `+∞` like an undefined energy; the starting point is an error), the same floor as a collapse
rule of the remeshing, and the vertices of faces within a factor 2 of the floor are blocked
(gradient and direction zeroed, like pinned) so that one squeezed face does not stop the whole
optimization. The first version of the floor (no blocking) made the line search fail at the wall:
spot control mesh, intrinsic, 300 it/round, ended after 164 iterations at E = 16.7 (budget 900;
E = 6.2 before the floor), and the next rounds made 0 iterations; with the blocking the same flow
reaches E = 5.05 in 900 iterations. The default flow (E_λ, 2 rounds x 500 it) on spot, bob and blub
(triangulated and control meshes) runs the whole budget (`runs_tmp/dev-robust/default/`). Cost of
the floor on the energy: on 24 small runs the final energy with floors 1e-5, 1e-4, 1e-3 is within
1 % of the unconstrained run in most cases (see SPEC D31).

**2. Intrinsic energy undefined (blub, 3000 it/round).** The error is raised at the start of
round 1: `evaluate_energy` of the *subdivided* mesh. Mechanism: the end state of round 0 (93
vertices, E = 0.79) is a flat pancake, many edges with two big faces of nearly opposite normals
(normal cosines -1.000000000, -0.999999486, -0.999998921 ..., areas 0.45 / 0.41); the midpoint of
such an edge has two (planar-arc) half-stars with nearly opposite normals, vector area nearly
zero, so D12 applies there although it did not at the old vertices. `E_λ` flows reach the same
pancakes (E = 1e-11, smallest fold gap 4e-10, spot / blub / bob covariance 3000). So the line search
treats undefined trials as `+∞` correctly; the undefined state appears only in the subdivision
(an operation outside the line search). Candidates measured:

* A second wall on the fold angle (`min_fold_angle` = 1e-3, 1e-4, 1e-5, 1e-6 rad). It does not
  work: a hard wall in L-BFGS ends the run when one edge reaches it (`LineSearchFailed`; blub `E_λ`
  with 1e-3: E = 1.4 after 0.05 s, unconstrained 2.7e-11; spot `E_λ` 1e-5: 0.021, unconstrained
  5.5e-8), and even at 1e-6 spot intrinsic-max still hit the undefined state (sweeps in
  `runs_tmp/dev-robust/sweep/`, `sweep2/`). Removed again.
* Decision (D32): evaluate the energy of the next round on the candidate state before committing. If
  it is undefined, the flow does not run that round, returns the state of the previous round and
  reports `FlowReport::stopped = Some(FlowStop { round, after_subdivision, vertex, reason })` (an
  error for round 0). A remeshing result with an undefined energy or a face below the floor is not
  applied (`RemeshReport::rejected`). Nothing is repaired.
* In the over-run set below this happens in 2 of 18 flows (spot intrinsic: before round 2; spot
  intrinsic-max: before round 1); the other 16 complete all rounds. Blub intrinsic 3000 now
  completes (E = 1.07, was an error).

**3. Found by the oracle fuzz (D33):** a boundary star with at most two faces contributed noise to
the `E_λ` gradient (the eigenvector of a rank-deficient matrix is arbitrary); it now contributes
exactly nothing.

#### Over-runs: 3 rounds of 3000 iterations, spot / bob / blub control meshes, six energies

All 18 flows complete without an error; none has a face below 1e-4 rad (`min angle`), the result of
each builds a `Geometry`. "faces < 1 deg" are the thin faces left by the remeshing (reported as
`FlowReport::bad_faces`). Before (cut study binary): blub `E_P`, spot `E_P`, bob `E_P`/intrinsic left
degenerate faces (angles 3e-10 to 2e-11 rad), blub intrinsic and spot intrinsic-max failed with
errors. `min fold gap` (pi minus the largest angle between adjacent normals) is a measurement, not a
constraint. Seconds are single-thread under a machine load of 60 to 150.

| run | faces | status | E_final | s | inadmissible | min angle | faces < 1 deg | min fold gap |
|---|---|---|---|---|---|---|---|---|
| blub_ctl-combinatorial-max | 3324 | max/max/max | 2.270763e0 | 20.905 | 4 | 2.13e-02 | 0 | 1.12e-03 |
| blub_ctl-combinatorial | 2144 | stalled/stalled/stalled | 1.129796e0 | 0.923 | 54 | 3.30e-02 | 0 | 3.41e-08 |
| blub_ctl-covariance-max | 3012 | max/max/max | 2.084945e-2 | 43.724 | 14 | 1.10e-04 | 9 | 4.17e-08 |
| blub_ctl-covariance | 2832 | max/lsf/max | 9.052823e-9 | 11.960 | 0 | 1.44e-04 | 9 | 9.51e-10 |
| blub_ctl-intrinsic-max | 3192 | max/max/max | 1.532855e0 | 31.878 | 16 | 1.15e-04 | 3 | 1.98e-04 |
| blub_ctl-intrinsic | 2860 | stalled/lsf/lsf | 1.070741e0 | 0.571 | 2 | 3.03e-02 | 0 | 3.13e-03 |
| bob_ctl-combinatorial-max | 10192 | max/max/max | 2.529696e0 | 86.508 | 0 | 2.00e-02 | 0 | 5.43e-01 |
| bob_ctl-combinatorial | 8596 | max/max/max | 4.286885e-1 | 86.455 | 2 | 1.90e-04 | 5 | 1.79e-01 |
| bob_ctl-covariance-max | 10462 | max/max/max | 4.346150e-1 | 123.544 | 0 | 1.93e-02 | 0 | 1.85e-05 |
| bob_ctl-covariance | 8366 | max/max/max | 2.417610e-2 | 71.141 | 3 | 1.11e-04 | 1 | 1.14e-06 |
| bob_ctl-intrinsic-max | 10152 | max/max/max | 2.002256e0 | 83.406 | 0 | 2.57e-02 | 0 | 1.11e-02 |
| bob_ctl-intrinsic | 9040 | max/lsf/stalled | 2.444759e0 | 49.355 | 13 | 1.18e-02 | 2 | 3.26e-05 |
| spot_ctl-combinatorial-max | 5442 | max/max/max | 1.008967e1 | 38.432 | 8 | 1.20e-04 | 2 | 2.18e-04 |
| spot_ctl-combinatorial | 3952 | max/max/max | 1.037054e0 | 35.798 | 47 | 1.03e-04 | 19 | 2.52e-05 |
| spot_ctl-covariance-max | 4790 | max/max/max | 2.640672e-1 | 61.758 | 0 | 1.78e-02 | 0 | 2.66e-08 |
| spot_ctl-covariance | 3808 | max/max/max | 7.550502e-4 | 21.302 | 20 | 1.05e-04 | 11 | 7.10e-09 |
| spot_ctl-intrinsic-max | 324 | lsf | 1.315386e1 | 1.728 | 0 | 4.80e-02 | 0 | 1.64e-03 |
| spot_ctl-intrinsic | 1260 | max/lsf | 3.322428e0 | 7.378 | 19 | 1.95e-04 | 2 | 4.56e-07 |

#### Before / after on the flows of the cut study (45 flows)

Same flows as `examples/cut_study.sh` (control meshes 300 and 3000 it/round, dense meshes 300
it/round), then the cut with covariance terminals, eps 1e-3, defect tolerance 0.1 (D30). Old =
`runs_tmp/dev-cut-study`, new = `runs_tmp/dev-robust/study`. "flipped uv / status": flipped UV
triangles after the cut, then the optimizer status per round (lsf = line search failed). All 45
flows are listed; the dense-mesh flows (300 it/round) are identical in energy, iterations and faces,
the crumpled ones differ in the 4th digit or by 2 faces (the floor changes a few accepted steps).
Seconds differ because of the machine load (about 60 to 150 during the new runs, less during the
old), not because of the code (the check is O(|F|) per evaluation, measured below 1 % of an
evaluation). The cut succeeds on all 45 (3 failed before; the 4th failure was the flow itself) and
flips 0 UV triangles everywhere.

| flow | E_final old | E_final new | iterations old | new | s old | s new | faces old | new | flipped uv / status old | flipped uv / status new | new min angle | new bad faces |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CrumpledDevelopable-combinatorial-3000 | 4.404057e-6 | 4.404057e-6 | 3000 | 3000 | 1.787 | 1.835 | 507 | 507 | 0 / max | 0 / max | 5.884e-2 | 0 |
| CrumpledDevelopable-combinatorial-500 | 1.478583e-3 | 1.478583e-3 | 500 | 500 | 0.322 | 0.322 | 510 | 510 | 0 / max | 0 / max | 4.100e-2 | 0 |
| CrumpledDevelopable-covariance-3000 | 9.034545e-8 | 1.189913e-7 | 3000 | 3000 | 1.262 | 1.368 | 497 | 497 | 0 / max | 0 / max | 4.890e-2 | 0 |
| CrumpledDevelopable-covariance-500 | 5.724026e-4 | 8.912045e-4 | 500 | 500 | 0.253 | 0.221 | 499 | 501 | 0 / max | 0 / max | 1.954e-2 | 0 |
| CrumpledDevelopable-intrinsic-3000 | 2.723743e-7 | 2.326685e-7 | 3000 | 3000 | 1.342 | 1.351 | 506 | 508 | 0 / max | 0 / max | 5.655e-2 | 0 |
| CrumpledDevelopable-intrinsic-500 | 6.472441e-5 | 7.021902e-5 | 500 | 500 | 0.239 | 0.194 | 508 | 508 | 0 / max | 0 / max | 5.745e-2 | 0 |
| blub_ctl-combinatorial-300 | 1.733067e0 | 1.678362e0 | 807 | 900 | 2.992 | 2.022 | 3330 | 3292 | cut fails / max/max/stalled | 0 / max/max/max | 8.917e-3 | 2 |
| blub_ctl-combinatorial-3000 | 3.733598e0 | 1.129796e0 | 602 | 2795 | 0.389 | 0.792 | 3204 | 2144 | cut fails / lsf/lsf/lsf | 0 / stalled/stalled/stalled | 3.296e-2 | 0 |
| blub_ctl-covariance-300 | 1.319488e-2 | 1.319488e-2 | 900 | 900 | 1.231 | 0.999 | 2868 | 2868 | 0 / max/max/max | 0 / max/max/max | 2.020e-2 | 0 |
| blub_ctl-covariance-3000 | 2.719293e-11 | 9.052823e-9 | 6056 | 7252 | 2.491 | 7.702 | 2788 | 2832 | 0 / max/max/lsf | 0 / max/lsf/max | 1.444e-4 | 9 |
| blub_ctl-intrinsic-300 | 2.142529e0 | 2.142529e0 | 900 | 900 | 1.811 | 1.979 | 3040 | 3040 | 0 / max/max/max | 0 / max/max/max | 1.767e-2 | 0 |
| blub_ctl-intrinsic-3000 | ERR error: flow failed: the energy is undefi | 1.070741e0 |  | 1360 |  | 0.383 |  | 2860 |  | 0 / stalled/lsf/lsf | 3.030e-2 | 0 |
| blub_quadrangulated-combinatorial-300 | 4.146214e0 | 4.146214e0 | 600 | 600 | 29.408 | 30.676 | 56806 | 56806 | 0 / max/max | 0 / max/max | 2.348e-2 | 0 |
| blub_quadrangulated-covariance-300 | 5.186006e0 | 5.186006e0 | 600 | 600 | 27.035 | 27.720 | 56598 | 56598 | 0 / max/max | 0 / max/max | 1.906e-2 | 0 |
| blub_quadrangulated-intrinsic-300 | 6.198339e0 | 6.198339e0 | 600 | 600 | 33.605 | 33.360 | 56604 | 56604 | 0 / max/max | 0 / max/max | 1.734e-2 | 4 |
| blub_triangulated-combinatorial-300 | 3.724887e0 | 3.724887e0 | 600 | 600 | 30.003 | 29.304 | 56822 | 56822 | 0 / max/max | 0 / max/max | 2.206e-2 | 0 |
| blub_triangulated-covariance-300 | 5.632107e0 | 5.632107e0 | 600 | 600 | 33.106 | 33.190 | 56614 | 56614 | 0 / max/max | 0 / max/max | 1.801e-2 | 0 |
| blub_triangulated-intrinsic-300 | 6.474141e0 | 6.474141e0 | 600 | 600 | 35.852 | 36.584 | 56634 | 56634 | 0 / max/max | 0 / max/max | 1.759e-2 | 0 |
| bob_ctl-combinatorial-300 | 1.229317e0 | 1.229317e0 | 900 | 900 | 5.960 | 5.657 | 10422 | 10422 | 0 / max/max/max | 0 / max/max/max | 1.813e-2 | 0 |
| bob_ctl-combinatorial-3000 | 9.014284e-1 | 4.286885e-1 | 2362 | 9000 | 6.799 | 45.902 | 9278 | 8596 | 0 / lsf/stalled/stalled | 0 / max/max/max | 1.899e-4 | 5 |
| bob_ctl-covariance-300 | 9.265637e-1 | 9.265637e-1 | 900 | 900 | 4.505 | 4.268 | 10116 | 10116 | 0 / max/max/max | 0 / max/max/max | 1.833e-2 | 0 |
| bob_ctl-covariance-3000 | 1.194679e-2 | 2.417610e-2 | 9000 | 9000 | 30.704 | 32.823 | 8292 | 8366 | 0 / max/max/max | 0 / max/max/max | 1.108e-4 | 1 |
| bob_ctl-intrinsic-300 | 1.639766e0 | 1.639766e0 | 900 | 900 | 6.247 | 5.448 | 10254 | 10254 | 0 / max/max/max | 0 / max/max/max | 2.454e-2 | 0 |
| bob_ctl-intrinsic-3000 | 4.165154e0 | 2.444759e0 | 2078 | 6576 | 6.214 | 30.871 | 9126 | 9040 | 0 / lsf/lsf/lsf | 0 / max/lsf/stalled | 1.185e-2 | 2 |
| bob_isotropic-combinatorial-300 | 2.863124e0 | 2.863124e0 | 600 | 600 | 14.504 | 14.575 | 24638 | 24638 | 0 / max/max | 0 / max/max | 2.761e-2 | 0 |
| bob_isotropic-covariance-300 | 4.627136e0 | 4.627136e0 | 600 | 600 | 14.970 | 14.220 | 24592 | 24592 | 0 / max/max | 0 / max/max | 1.890e-2 | 0 |
| bob_isotropic-intrinsic-300 | 4.626109e0 | 4.626109e0 | 600 | 600 | 15.563 | 14.990 | 24610 | 24610 | 0 / max/max | 0 / max/max | 1.895e-2 | 0 |
| bob_quad-combinatorial-300 | 3.491903e0 | 3.491903e0 | 600 | 600 | 23.304 | 23.207 | 42734 | 42734 | 0 / max/max | 0 / max/max | 2.393e-2 | 0 |
| bob_quad-covariance-300 | 7.352449e0 | 7.352449e0 | 600 | 600 | 22.278 | 22.826 | 42660 | 42660 | 0 / max/max | 0 / max/max | 2.074e-2 | 0 |
| bob_quad-intrinsic-300 | 8.249370e0 | 8.249370e0 | 600 | 600 | 21.949 | 24.263 | 42630 | 42630 | 0 / max/max | 0 / max/max | 2.951e-2 | 0 |
| bob_tri-combinatorial-300 | 3.076487e0 | 3.076487e0 | 600 | 600 | 23.108 | 23.044 | 42750 | 42750 | 0 / max/max | 0 / max/max | 2.283e-2 | 0 |
| bob_tri-covariance-300 | 7.024482e0 | 7.024482e0 | 600 | 600 | 22.363 | 23.569 | 42662 | 42662 | 0 / max/max | 0 / max/max | 3.613e-2 | 0 |
| bob_tri-intrinsic-300 | 7.126420e0 | 7.126420e0 | 600 | 600 | 23.987 | 25.133 | 42656 | 42656 | 0 / max/max | 0 / max/max | 2.476e-2 | 0 |
| spot_ctl-combinatorial-300 | 1.530048e0 | 1.530048e0 | 900 | 900 | 3.072 | 3.235 | 5512 | 5512 | 0 / max/max/max | 0 / max/max/max | 2.640e-2 | 0 |
| spot_ctl-combinatorial-3000 | 2.485319e0 | 1.037054e0 | 975 | 9000 | 1.213 | 23.911 | 5216 | 3952 | cut fails / lsf/lsf/lsf | 0 / max/max/max | 1.027e-4 | 19 |
| spot_ctl-covariance-300 | 3.154606e-1 | 3.154606e-1 | 900 | 900 | 2.819 | 2.431 | 5112 | 5112 | 0 / max/max/max | 0 / max/max/max | 1.805e-2 | 0 |
| spot_ctl-covariance-3000 | 1.092982e-8 | 7.550502e-4 | 9000 | 9000 | 12.045 | 14.817 | 3888 | 3808 | 0 / max/max/max | 0 / max/max/max | 1.047e-4 | 11 |
| spot_ctl-intrinsic-300 | 6.202401e0 | 5.054394e0 | 900 | 900 | 3.373 | 2.347 | 5584 | 5448 | 0 / max/max/max | 0 / max/max/max | 1.024e-4 | 7 |
| spot_ctl-intrinsic-3000 | 4.892556e0 | 3.322428e0 | 6313 | 5359 | 23.363 | 5.157 | 5484 | 1260 | 0 / stalled/max/max | 0 / max/lsf | 1.951e-4 | 2 |
| spot_quadrangulated-combinatorial-300 | 8.276170e0 | 8.276170e0 | 600 | 600 | 13.451 | 12.570 | 23392 | 23392 | 0 / max/max | 0 / max/max | 2.351e-2 | 0 |
| spot_quadrangulated-covariance-300 | 9.867035e0 | 9.867035e0 | 600 | 600 | 14.241 | 14.118 | 23044 | 23044 | 0 / max/max | 0 / max/max | 1.869e-2 | 0 |
| spot_quadrangulated-intrinsic-300 | 1.351481e1 | 1.351481e1 | 600 | 600 | 16.480 | 15.103 | 23008 | 23008 | 0 / max/max | 0 / max/max | 1.799e-2 | 0 |
| spot_triangulated-combinatorial-300 | 8.062792e0 | 8.062792e0 | 600 | 600 | 13.707 | 12.577 | 23388 | 23388 | 0 / max/max | 0 / max/max | 2.784e-2 | 0 |
| spot_triangulated-covariance-300 | 1.143797e1 | 1.143797e1 | 600 | 600 | 12.757 | 12.467 | 23168 | 23168 | 0 / max/max | 0 / max/max | 2.085e-2 | 0 |
| spot_triangulated-intrinsic-300 | 1.664337e1 | 1.664337e1 | 600 | 600 | 17.093 | 16.104 | 23278 | 23278 | 0 / max/max | 0 / max/max | 2.856e-2 | 0 |

Reading: where the floor is not reached nothing changes (identical energies). Where it is reached
the flows run further (bob `E_P` 3000: 9000 iterations, E 0.43 instead of 0.90 after 2362; spot `E_P`
3000: E 1.04 instead of 2.49 after 975), because the old runs ended in `LineSearchFailed` on
squeezed faces. `E_λ` over-runs of the control meshes end at a higher energy than before (spot
7.6e-4 instead of 1e-8, bob 2.4e-2 instead of 1.2e-2) because they no longer collapse faces to
reach the flat pancake; blub `E_λ` 3000 ends at 9e-9 instead of 2.7e-11. These are the
over-runs; with 300 it/round the energies are identical.

### 12.5 Study: every model, six energies, two budgets

`sh crates/ddg-developability/examples/study.sh 60 "100 400"` after
`cargo build --release -p ddg-developability --example develop`. Every OBJ of `data/models`
(11 files; polygons fan-triangulated) x 6 energies x iteration budgets 100 and 400 per round, two
rounds (coarse + one subdivision), 1 thread per job, 60 jobs in parallel, ε = 1e-3, default options
(remeshing on). 132 jobs in 3 min 7 s wall (machine load average about 60 from other agents; the
times below are single-thread seconds under that load). Logs: `runs_tmp/dev-flow-study/job-*.log`.
Result: all 132 jobs finished; no error, no non-finite energy or position, no line-search failure,
no `NonDescent`; 1 stall (`EnergyStalled` in round 0 of spot_control_mesh Intrinsic, budget 400);
undefined trials: 1 in total; cut and flatten succeeded in all 132.
Columns: `flips/collapses` are remesh operations over both rounds; `cut` is terminals(+added)/cut
edges/charts; `flipped uv` and `mean |u|` are from the cut at ε = 1e-3 on the final mesh. E0 is
the energy of the input (energy-specific). Energies of different kinds are not comparable.

Budget 400 iterations per round:

| mesh | energy | faces in->out | E0 | E_final | s | status | flips/collapses | cut: term(+add)/edges/charts | flipped uv | mean abs u |
|---|---|---|---|---|---|---|---|---|---|---|
| CrumpledDevelopable | Combinatorial | 512->2040 | 33.9 | 0.00173 | 4.116 | max_iter/max_iter | 0/1 | 0(+0)/0/1 | 0 | 0.000831 |
| CrumpledDevelopable | CombinatorialMax | 512->2044 | 126 | 0.0224 | 1.774 | max_iter/max_iter | 0/1 | 1(+0)/6/1 | 0 | 0.00115 |
| CrumpledDevelopable | Covariance | 512->2000 | 61.9 | 0.000639 | 3.534 | max_iter/max_iter | 1/8 | 0(+0)/0/1 | 0 | 0.011 |
| CrumpledDevelopable | CovarianceMax | 512->2024 | 18.3 | 0.00209 | 3.505 | max_iter/max_iter | 0/3 | 0(+0)/0/1 | 0 | 0.00149 |
| CrumpledDevelopable | Intrinsic | 512->2032 | 79.6 | 6.68e-05 | 1.305 | max_iter/max_iter | 0/2 | 0(+0)/0/1 | 0 | 2e-05 |
| CrumpledDevelopable | IntrinsicMax | 512->2040 | 91.3 | 0.00518 | 1.676 | max_iter/max_iter | 0/1 | 0(+0)/0/1 | 0 | 0.000933 |
| blub_control_mesh | Combinatorial | 220->812 | 23.7 | 2.11 | 1.526 | max_iter/max_iter | 6/25 | 125(+0)/152/1 | 57 | 1.45 |
| blub_control_mesh | CombinatorialMax | 220->852 | 98.2 | 4.62 | 1.637 | max_iter/max_iter | 5/5 | 280(+0)/301/1 | 3 | 0.0873 |
| blub_control_mesh | Covariance | 220->700 | 37.6 | 0.00983 | 0.988 | max_iter/max_iter | 7/33 | 2(+0)/16/1 | 36 | 2.52 |
| blub_control_mesh | CovarianceMax | 220->788 | 10.3 | 0.0641 | 2.324 | max_iter/max_iter | 0/13 | 16(+0)/41/1 | 57 | 2.13 |
| blub_control_mesh | Intrinsic | 220->746 | 89 | 1.48 | 1.002 | max_iter/max_iter | 2/19 | 137(+0)/146/1 | 35 | 1.22 |
| blub_control_mesh | IntrinsicMax | 220->832 | 76 | 3.5 | 2.103 | max_iter/max_iter | 0/6 | 227(+0)/250/1 | 19 | 0.113 |
| blub_quadrangulated | Combinatorial | 14208->56798 | 19 | 3.88 | 167.064 | max_iter/max_iter | 0/5 | 663(+0)/1690/1 | 0 | 0.0378 |
| blub_quadrangulated | CombinatorialMax | 14208->56828 | 61.8 | 15.8 | 140.425 | max_iter/max_iter | 0/2 | 3150(+0)/5354/1 | 0 | 0.00633 |
| blub_quadrangulated | Covariance | 14208->56542 | 45.7 | 4.79 | 169.838 | max_iter/max_iter | 0/49 | 1041(+0)/2104/1 | 0 | 0.0117 |
| blub_quadrangulated | CovarianceMax | 14208->56824 | 10.5 | 1.91 | 165.817 | max_iter/max_iter | 0/1 | 214(+0)/746/1 | 0 | 0.252 |
| blub_quadrangulated | Intrinsic | 14208->56536 | 49.1 | 5.91 | 179.217 | max_iter/max_iter | 2/49 | 1279(+0)/2512/1 | 0 | 0.00919 |
| blub_quadrangulated | IntrinsicMax | 14208->56832 | 44.4 | 8.22 | 105.901 | max_iter/max_iter | 0/0 | 1722(+0)/3110/1 | 0 | 0.00812 |
| blub_triangulated | Combinatorial | 14208->56812 | 19.6 | 3.56 | 152.299 | max_iter/max_iter | 0/4 | 587(+0)/1495/1 | 0 | 0.0374 |
| blub_triangulated | CombinatorialMax | 14208->56832 | 62.2 | 17.1 | 161.897 | max_iter/max_iter | 0/0 | 2958(+0)/4912/1 | 0 | 0.00788 |
| blub_triangulated | Covariance | 14208->56534 | 45.6 | 5.13 | 181.483 | max_iter/max_iter | 0/59 | 1168(+0)/2406/1 | 0 | 0.00964 |
| blub_triangulated | CovarianceMax | 14208->56824 | 10.7 | 2.04 | 167.955 | max_iter/max_iter | 0/1 | 259(+0)/1006/1 | 0 | 0.0919 |
| blub_triangulated | Intrinsic | 14208->56580 | 48.4 | 6.81 | 186.133 | max_iter/max_iter | 1/45 | 1375(+0)/2647/1 | 0 | 0.00867 |
| blub_triangulated | IntrinsicMax | 14208->56824 | 44.4 | 8.48 | 142.241 | max_iter/max_iter | 0/1 | 1865(+0)/3271/1 | 0 | 0.00878 |
| bob_controlmesh | Combinatorial | 668->2572 | 16.3 | 1.51 | 6.612 | max_iter/max_iter | 0/17 | 327(+0)/348/1 | 0 | 0.00112 |
| bob_controlmesh | CombinatorialMax | 668->2656 | 58 | 10.3 | 6.478 | max_iter/max_iter | 0/2 | 753(+0)/785/1 | 0 | 0.000141 |
| bob_controlmesh | Covariance | 668->2490 | 47.2 | 1.39 | 5.863 | max_iter/max_iter | 0/31 | 255(+0)/306/1 | 1 | 0.029 |
| bob_controlmesh | CovarianceMax | 668->2640 | 11.1 | 1.13 | 7.905 | max_iter/max_iter | 0/4 | 386(+0)/440/1 | 1 | 0.0143 |
| bob_controlmesh | Intrinsic | 668->2546 | 67.7 | 2.79 | 6.529 | max_iter/max_iter | 0/21 | 349(+0)/383/1 | 0 | 0.00122 |
| bob_controlmesh | IntrinsicMax | 668->2606 | 58.5 | 8.15 | 8.086 | max_iter/max_iter | 1/9 | 670(+0)/696/1 | 0 | 0.00984 |
| bob_isotropic | Combinatorial | 6174->24628 | 17.4 | 2.79 | 90.926 | max_iter/max_iter | 0/10 | 579(+0)/1136/1 | 0 | 0.0515 |
| bob_isotropic | CombinatorialMax | 6174->24696 | 56.8 | 11.5 | 74.044 | max_iter/max_iter | 0/0 | 2374(+0)/3028/1 | 0 | 0.00185 |
| bob_isotropic | Covariance | 6174->24554 | 43.8 | 4.3 | 105.770 | max_iter/max_iter | 2/26 | 1096(+0)/1558/1 | 0 | 0.00532 |
| bob_isotropic | CovarianceMax | 6174->24696 | 11.9 | 2.23 | 112.171 | max_iter/max_iter | 0/0 | 448(+0)/960/1 | 0 | 0.0488 |
| bob_isotropic | Intrinsic | 6174->24590 | 44.4 | 4.48 | 106.043 | max_iter/max_iter | 4/23 | 1163(+0)/1588/1 | 0 | 0.00482 |
| bob_isotropic | IntrinsicMax | 6174->24672 | 48.5 | 9.07 | 80.680 | max_iter/max_iter | 0/3 | 2103(+0)/2729/1 | 0 | 0.00202 |
| bob_quad | Combinatorial | 10688->42736 | 13.3 | 3.35 | 135.340 | max_iter/max_iter | 0/2 | 420(+0)/856/1 | 0 | 0.162 |
| bob_quad | CombinatorialMax | 10688->42752 | 38.2 | 13.3 | 91.447 | max_iter/max_iter | 0/0 | 3017(+0)/3895/1 | 0 | 0.0172 |
| bob_quad | Covariance | 10688->42640 | 42.1 | 6.96 | 142.175 | max_iter/max_iter | 3/23 | 2784(+0)/3217/1 | 0 | 0.00943 |
| bob_quad | CovarianceMax | 10688->42752 | 8.58 | 2.53 | 132.234 | max_iter/max_iter | 0/0 | 342(+0)/842/1 | 0 | 0.163 |
| bob_quad | Intrinsic | 10688->42584 | 42.9 | 7.58 | 144.470 | max_iter/max_iter | 2/36 | 3034(+0)/3434/1 | 0 | 0.00913 |
| bob_quad | IntrinsicMax | 10688->42752 | 34.5 | 10.3 | 131.256 | max_iter/max_iter | 0/0 | 2504(+0)/3180/1 | 0 | 0.0185 |
| bob_tri | Combinatorial | 10688->42738 | 13.6 | 3.01 | 136.723 | max_iter/max_iter | 0/4 | 400(+0)/838/1 | 0 | 0.159 |
| bob_tri | CombinatorialMax | 10688->42752 | 38 | 12.7 | 125.328 | max_iter/max_iter | 0/0 | 2784(+0)/3606/1 | 0 | 0.0159 |
| bob_tri | Covariance | 10688->42636 | 42.9 | 6.47 | 134.003 | max_iter/max_iter | 1/22 | 1765(+0)/2510/1 | 0 | 0.0103 |
| bob_tri | CovarianceMax | 10688->42752 | 8.81 | 2.32 | 131.083 | max_iter/max_iter | 0/0 | 291(+0)/775/1 | 0 | 0.158 |
| bob_tri | Intrinsic | 10688->42638 | 43.8 | 6.59 | 110.589 | max_iter/max_iter | 1/24 | 1701(+0)/2528/1 | 0 | 0.0102 |
| bob_tri | IntrinsicMax | 10688->42752 | 35.5 | 9.52 | 130.692 | max_iter/max_iter | 0/0 | 2195(+0)/2987/1 | 0 | 0.0163 |
| spot_control_mesh | Combinatorial | 372->1352 | 40.2 | 2.09 | 1.211 | max_iter/max_iter | 0/23 | 196(+0)/234/1 | 15 | 0.206 |
| spot_control_mesh | CombinatorialMax | 372->1424 | 170 | 17.5 | 1.069 | max_iter/max_iter | 0/8 | 535(+0)/544/1 | 0 | 0.00581 |
| spot_control_mesh | Covariance | 372->1236 | 99 | 0.294 | 2.772 | max_iter/max_iter | 3/42 | 70(+0)/123/1 | 55 | 0.871 |
| spot_control_mesh | CovarianceMax | 372->1382 | 24.3 | 1.02 | 5.375 | max_iter/max_iter | 1/14 | 146(+0)/177/1 | 12 | 0.336 |
| spot_control_mesh | Intrinsic | 372->1394 | 172 | 9.4 | 2.406 | stalled/max_iter | 4/17 | 437(+0)/451/1 | 0 | 0.000628 |
| spot_control_mesh | IntrinsicMax | 372->1438 | 168 | 21.8 | 3.795 | max_iter/max_iter | 0/10 | 533(+0)/549/1 | 0 | 0.000436 |
| spot_quadrangulated | Combinatorial | 5856->23362 | 42 | 8.22 | 71.721 | max_iter/max_iter | 0/13 | 1249(+0)/1866/1 | 0 | 0.00374 |
| spot_quadrangulated | CombinatorialMax | 5856->23416 | 126 | 31.6 | 77.321 | max_iter/max_iter | 0/4 | 3065(+0)/3839/1 | 0 | 0.000839 |
| spot_quadrangulated | Covariance | 5856->23000 | 111 | 8.78 | 80.551 | max_iter/max_iter | 7/83 | 1499(+0)/2009/1 | 2 | 0.00263 |
| spot_quadrangulated | CovarianceMax | 5856->23388 | 25 | 4.93 | 83.495 | max_iter/max_iter | 0/6 | 914(+0)/1398/1 | 0 | 0.0136 |
| spot_quadrangulated | Intrinsic | 5856->22902 | 115 | 13.1 | 76.713 | max_iter/max_iter | 7/99 | 1862(+0)/2457/1 | 0 | 0.00128 |
| spot_quadrangulated | IntrinsicMax | 5856->23334 | 104 | 23.2 | 67.407 | max_iter/max_iter | 0/15 | 2627(+0)/3200/1 | 0 | 0.00087 |
| spot_triangulated | Combinatorial | 5856->23360 | 39.1 | 7.64 | 60.386 | max_iter/max_iter | 1/11 | 1383(+0)/1860/1 | 0 | 0.00794 |
| spot_triangulated | CombinatorialMax | 5856->23406 | 123 | 30 | 55.407 | max_iter/max_iter | 0/3 | 3449(+0)/3961/1 | 0 | 0.000835 |
| spot_triangulated | Covariance | 5856->23086 | 112 | 10.6 | 74.507 | max_iter/max_iter | 5/64 | 1609(+0)/2078/1 | 0 | 0.00275 |
| spot_triangulated | CovarianceMax | 5856->23378 | 25.5 | 4.8 | 74.532 | max_iter/max_iter | 0/8 | 879(+0)/1363/1 | 0 | 0.0158 |
| spot_triangulated | Intrinsic | 5856->23256 | 118 | 15.6 | 79.991 | max_iter/max_iter | 1/27 | 2014(+0)/2596/1 | 0 | 0.00176 |
| spot_triangulated | IntrinsicMax | 5856->23314 | 105 | 24.8 | 57.766 | max_iter/max_iter | 0/16 | 2926(+0)/3500/1 | 0 | 0.000963 |

Reading the table: the dense meshes (5.8k to 14k faces, 23k to 57k after the subdivision) run in
60 to 190 s single-threaded for 800 iterations; the control meshes in 1 to 8 s. The cut
distortion is small (mean |u| ~ 1e-3 to 1e-2, 0 flips) on dense meshes, large on the coarse control
meshes (spot_control_mesh Covariance 0.87, blub_control_mesh 1 to 2.5 with 20 to 60 flips of ~800
faces): coarse meshes with few vertices cannot hide the residual curvature. Budget 100 per round
gives the same picture (`runs_tmp/dev-flow-study/table100.md`, not reproduced; e.g. blub_control_mesh
Covariance 85 flips, spot_control_mesh Covariance 97 flips, all dense meshes 0 flips).

### 12.6 Benchmarks of the flow

Measured alone on a quiet machine (AMD Ryzen Threadripper 3970X (32 cores, 64 threads), 2026-10-19, nothing else running (load average 0.1 to 1.5 before and during; the multithread parts load the machine themselves); `cargo bench -p ddg-developability --bench flow`,
log `runs_tmp/bench-quiet/flow.log`). "64 threads" is a local pool of 64 threads on the 32 physical
cores. Remeshing on, default options.

| job | faces | iterations | evaluations | 1 thread | 64 threads |
|---|---|---|---|---|---|
| spot_triangulated, E_λ, 50 it | 5852 | 50 | 62 | 0.27 s (5.4 ms/it, 4.4 ms/eval) | 0.28 s |
| spot_triangulated, intrinsic | 5850 | 50 | 67 | 0.32 s (6.3 ms/it) | 0.30 s |
| bob_tri, E_λ | 10688 | 50 | 58 | 0.45 s (9.0 ms/it, 7.7 ms/eval) | 0.36 s |
| bob_tri, intrinsic | 10688 | 50 | 57 | 0.48 s (9.6 ms/it) | 0.36 s |
| blub_triangulated, E_λ | 14208 | 50 | 61 | 0.62 s (12.5 ms/it, 10.2 ms/eval) | 0.46 s |
| blub_triangulated, intrinsic | 14208 | 50 | 53 | 0.60 s (12.0 ms/it) | 0.43 s |
| spot control + 2 subdivisions, E_λ, 3 x 300 it | 5112 | 900 | 1427 | 1.78 s (2.0 ms/it) | 2.80 s |
| crumpled, E_λ, 300 it | 503 | 300 | 432 | 0.13 s | 0.43 s |
| bob control + 2 subdivisions, E_λ, 3 x 300 it | 10116 | 900 | 1397 | 3.63 s (4.0 ms/it) | 4.21 s |

Per evaluation (value + gradient + line search overhead) about 0.7 us per face single-threaded
(7.7 ms for 10.7k faces). The paper: "seconds to a few minutes" single-threaded on a 3.4 GHz i7
laptop for its examples; our spot run (the SPEC target: <= 60 s single-threaded) takes 1.8 s for 900
iterations, i.e. within the paper's "seconds" and 30 times below the target. The long Blender runs
(§ of `blender/README.md`, 2 subdivisions, 1500 iterations on 92k to 226k faces, 1 thread) take 68
to 169 s. Threads do not help the flow at these sizes (5k to 14k faces: 0.8 to 1.4 times at best,
slower for the small ones: the pool overhead of every call); only the energy alone scales (§12.7).

### 12.7 Energies and optimizer: decisions that go beyond the text, observations

Energies (SPEC §1 to §3):

* Margins: `π_max` uses min(partition gap, active-pair gap within the chosen partition); the
  intrinsic width divides the eigen-gap by `Σθ` (= tr A_v), not by tr Ã (scale-free noise when
  nearly flat); boundary stars with at most two faces get margin `+∞` (energy identically 0, D33).
* `H(φ)` is evaluated by its series (three terms) below `1e-2` instead of `1e-3` (the closed form
  loses about ε/φ²).
* `threads: None` runs on the current rayon pool (the caller's or rayon's default); `Some(n)`
  builds a local pool unless already inside an `n`-thread pool.
* Pipeline (§3.9): faces, then stars (fixed 256-vertex chunks), then the per-face backpropagation,
  then the per-vertex gather in star order: bitwise identical for 1, 2, 7 and 64 threads (tested).
* `dev-regular-cone-values` is not tagged in Rust (the tests compare closed forms; Lean decides).
* The star layout (`Stars`) is rebuilt on every call; an evaluator that keeps it across optimizer
  iterations would save part of each call (unmeasured, not in the API).
* Two-face boundary stars (D33): `star_energy` zeroes value, margin-tie, `y`, `s`, `z` after
  computing them (so the intrinsic energies still report D12 undefined states). Before, `E_λ` took
  an arbitrary eigenvector of the rank-deficient matrix and fed sensitivities of the size of the
  angle between the normals into the gradient (found by the oracle fuzz; `tests/robust.rs`).
* Energy bench (`cargo bench -p ddg-developability --bench energy`, value + gradient; quiet
  machine, AMD Ryzen Threadripper 3970X (32 cores, 64 threads), 2026-10-19, nothing else running (load average 0.1 to 1.5 before and during; the multithread parts load the machine themselves); logs `runs_tmp/bench-quiet/energy-{4,16,64}.log`), single thread, ns per vertex
  (spot / bob / blub / spot subdivided twice, 93696 faces): `E_P` 1547 / 1180 / 1206 / 1608,
  `π_max` 1014 / 992 / 959 / 1478, `E_λ` 756 / 742 / 741 / 1200, `λ_max` 3274 / 2897 / 2970 / 3860,
  intrinsic 880 / 851 / 833 / 1297, intrinsic max 1210 / 1167 / 1168 / 1634. All meet the §10.8
  targets (<= 2 µs, <= 5 µs for the max energies). Speed-up over one thread (the three runs
  repeat the single-thread numbers within 5 %): spot-sub2 `E_λ` 2.2 (4 threads), 3.2 (16), 5.5
  (64); `λ_max` 3.5 / 10.4 / 8.1; the other four energies 2.2 to 2.8 / 3.2 to 4.3 / 3.8 to 7.0.
  The 5.8k to 14k-face meshes: `E_λ` 1.3 to 2.9 (4), 1.6 to 2.6 (16), 0.9 to 1.4 (64); `λ_max`
  3.1 to 3.6 / 4.3 to 6.7 / 3.0 to 5.1. The scaling is modest (a call is 2 to 80 ms; the pool is
  built per measurement and the gather is serial per vertex star). `eig` bench: `sym_eig3` 233 ns,
  `sym_eig2` 37 ns per call.

Optimizer, subdivision, remeshing (SPEC §4, §5):

* `src/optimize/`: the generic `minimize(objective, x0, Problem{length_scale, group}, opts,
  callback)` (`core.rs`), options, report and statuses (`options.rs`), the mesh entry point
  `optimize` (`mesh.rs`). L-BFGS (two-loop recursion, H0 from the newest pair), Lewis–Overton
  Algorithm 2.6, gradient descent with weak Wolfe or backtracking, fixed step; units of the mean
  edge length, step cap, pinned vertices.
* Objective trait: `EvalError::Undefined` (+∞ in a line search, an error at the start) against
  `Fatal`. The mesh objective maps `UndefinedEnergy`, `DegenerateFace`, `NonFinitePosition` at a
  trial point to undefined. Fixed step: an undefined trial point ends the run with
  `LineSearchFailed`. L-BFGS with `Backtracking` is allowed (a pair is skipped if `sᵀy <= 0`).
* A line-search failure with non-empty memory: reset the memory and retry once from steepest
  descent; a failure from steepest descent ends the run. A non-descent L-BFGS direction resets the
  memory; a non-descent steepest direction ends the run (`NonDescent`) at once.
* Remeshing beyond the text: flips are also refused when the normal of a new face turns by `>= 90`
  degrees against an old face of the quad (reflex quad); the adjacency rule is tested before the
  valence rule (the valence rules are implied by adjacency / the link condition at interior
  vertices and are kept as safety nets); a pinned vertex survives a collapse (SPEC: the lower id)
  and `Remeshed::pinned` gives the new pinned ids; `RemeshOptions::validate_each_operation` is for
  tests; `energy_before/after` are `Option<f64>`, filled in by the flow driver; refusal reasons
  are `&'static str` and unit-tested one by one (`cargo test --lib remesh`).
* BFGS on Nesterov's second Chebyshev–Rosenbrock function (abs inside) stalls at non-minimizing
  Clarke stationary points (known: Gürbüzbalaban–Overton); on the first (squared) variant, `n = 3`,
  L-BFGS reaches `f < 1e-4` in 2000 iterations from 20 of 20 random starts.
* The admissibility floor (D31) keeps the L-BFGS memory across a change of the blocked set; the
  curvature test `sᵀy > 0` skips the bad pairs (`skipped_pairs`). The first version of the floor
  without blocking failed in the line search at the wall (§12.4).

### 12.8 Open items

* The intrinsic energies with `omit_valence_three = false` error on undefined vertices (the cut
  stage propagates the error); not seen on any data model.
* Out of scope (D23): gluing cracks and avoiding global overlap after the cut (the paper does both
  by hand). The end-to-end numbers of the Blender operators are in `blender/README.md`
  ("Measured").
