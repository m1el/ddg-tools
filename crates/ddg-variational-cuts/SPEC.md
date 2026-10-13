# Variational Surface Cutting — discrete specification

The contract between the Rust crate (`crates/ddg-variational-cuts`), the Lean oracle
(`lean/oracle`) and the Lean proofs (`lean/proofs`) for

> Nicholas Sharp, Keenan Crane. *Variational Surface Cutting.* ACM Transactions on Graphics 37(4),
> 2018. <https://doi.org/10.1145/3197517.3201356> (`../VariationalCuts/paper.pdf`)

**Clean room.** This is written from the text of the paper and the author's two derivation notes
(`../VariationalCuts/cea_tutorial.pdf`, `cea_derivations.pdf`) only. The authors' code is GPL and
must never be opened, searched, fetched or paraphrased by anyone working on this crate. Papers the
text cites may be downloaded and cited.

Where the paper is silent or inconsistent, the answer below is derived from the energy in §1 and
the difference is listed in §9 with an id `D<n>`. If an implementation or a test contradicts this
file, stop and report it; do not quietly follow either side.

Workspace conventions of `CLAUDE.md` apply: the cotan Laplacian `L` is **positive semidefinite**,
`u` is the **true** log conformal scale factor (flattened length ≈ `e^u ×` surface length),
`f64`, `u32` stored indices, everything intrinsic (edge lengths), nothing silently clamped.

## 1. Smooth problem and signs

A surface `M` (with or without boundary) is cut by a curve network `γ` into regions `R`. In every
region `u` solves the Yamabe equation with zero boundary values (paper Eq. 3, `Δ` PSD):

```
Δu = −K  in R,      u = 0  on ∂R   (the cut and the boundary of M)
```

Energy (paper Eq. 5, 17, 19, 30, and the area term of Fig. 14):

```
E(γ) = Σ_R [ w_R · J_R(u)  +  β_R · Area(R) ]  +  ∫_γ (α_L + η) ds
J = ∫ |∇u|²          Dirichlet            (Eq. 2)
J = ∫ a |∇u|²        scaled Dirichlet     (Eq. 17),  a > 0 a fixed function on M
J = ∫ b(u)           pointwise            (Sec. 3.5.2);  Hencky is b(u) = u²  (Eq. 19)
```

`γ` is counted once in the length term (D3). `w_R ≥ 0` (default 1) and `β_R` (default 0) are per
region; `η ≥ 0` is a function on `M` (default 0).

**Gradient density.** If the boundary of one region `R` moves with outward normal speed `σ`, then
`d/dt [w_R J_R + β_R Area(R)] = ∫_∂R G_R σ ds` with (paper Table 1, `n` the outward normal of `R`)

```
Dirichlet          G = (∂u/∂n)²
scaled Dirichlet   G = −a (∂u/∂n)² − (∂u/∂n)(∂p/∂n),     Δp = −2 Δ_a u,   p = 0 on ∂R,   Δ_a := −∇·(a∇)
pointwise          G = b(u) − (∂u/∂n)(∂p/∂n),            Δp = −b′(u),     p = 0 on ∂R
G_R = w_R · G + β_R
```

Checks: with `a ≡ 1`, `p = −2u` and the scaled density equals the Dirichlet one. On a flat disk of
radius `ρ` with constant source `K`: `u = K (r² − ρ²)/4`, `∫|∇u|² = π K² ρ⁴ / 8`, whose derivative
in `ρ` is `(Kρ/2)² · 2πρ`; `∫u² = π K² ρ⁶ / 48`, whose derivative is `π K² ρ⁵ / 8 = −u′(ρ) p′(ρ) · 2πρ`
with `p′(ρ) = −K ρ³ / 8`.

**Flow in level-set form.** Two regions are described by one function `φ`, a signed distance
(`|∇φ| = 1`), with region `0 = {φ ≤ 0}` and region `1 = {φ > 0}`. `n = ∇φ` is the outward normal of
region 0. Moving the cut with speed `σ` along `n` changes the energy at the rate

```
dE = ∫_γ [ G_0 − G_1 + ∂(α_L+η)/∂n + (α_L+η) κ ] σ ds,        κ = div n = −Δφ
```

Steepest descent is `σ = −[…]`, and `φ_t = −σ`, so on `γ`

```
φ_t = (G_0 − G_1)  −  Δ_w φ,       w := α_L + η,     Δ_w = −∇·(w∇)                        (★)
```

In the paper the first term is extended harmonically away from `γ` and integrated explicitly, the
second is integrated implicitly on `M`, and then `φ` is made a signed distance again. Here both
terms are evaluated on `γ` by default, the second implicitly, and extended together (§3.8–§3.10,
D14–D16); the stages of the paper stay selectable. A region whose own boundary derivative is large
**shrinks** (D2).

**Roughening.** Linearizing (★) about a straight cut between two half planes with constant normal
derivatives shows what the two terms do to a ripple `δ(x) = ε cos(kx)` of the cut: where a region
bulges out, its `∂u/∂n` falls by the factor `1 − kδ` (its boundary values are carried by a
harmonic function that decays like `e^{−k·distance}`), on the other side it rises by `1 + kδ`, so

```
δ_t = ( ω k − w k² ) δ,          ω := 2 (Q_0 + Q_1),   Q_r := G_r − β_r − w_r b(0)
```

with `Q_r` the part of the gradient density that is quadratic in the normal derivatives
(`w_r (∂u/∂n)²` for the Dirichlet energy). The distortion term amplifies every ripple, the shorter
the faster; the length term damps it. Ripples longer than `2π w/ω` grow, the fastest at
`k = ω/(2w)`. With `w = 0` the flow is ill-posed (the paper's Sec. 3.2.1: the distortion alone is
driven to zero by longer and longer cuts). `ω` has the unit of a speed, is reported by every step,
and sets the scale of `α_L` and of the time step (§4, D18). On the sphere of unit area
(`R = 1/√(4π)`), for the Dirichlet energy, the ripple `ε cos(mθ)` of the equator grows at the rate
`(m − 1) [ 4/R³ − α_L (m + 1)/R² ]` exactly (the DtN map of a hemisphere is `m/R`, the base
solution contributes `u″ = K`, and the curvature of the ripple is `(m² − 1) ε/R²`); the equator is
stable iff `α_L > 4/(3R) = 4.73`. For `N` regions (§8) every region `k` has its own `φ_k`,
negative inside, with `φ_{k,t} = G_k − G_j − Δ_w φ_k` on the part of `∂R_k` shared with region `j`.

## 2. Discrete setting

`M = (V, E, F)`: a `ddg_core::TriMesh` with a `ddg_core::Geometry` (edge lengths). Edge `e` is
`[a, b]` with `a < b`. `∂M` vertices are `B_M`.

**Normalization.** The flow runs on the mesh scaled to unit total area: every length is multiplied
by `s = 1/√(total area)`. All parameters (`α_L`, `η`, `β`, `τ`, `r_ε`, `c`, a curvature override)
refer to the scaled surface. Level-set values cross the API in the caller's length unit and are
multiplied by `s` on the way in and divided by it on the way out. Reports state which unit each
number is in. Formulas below are in the scaled unit.

`A_f` face area, `A_i = ⅓ Σ_{f∋i} A_f` (lumped vertex area, `W_V = diag(A_i)`), `θ_i^f` corner angle,
`Θ_i = Σ_f θ_i^f`.

## 3. Stages of one evaluation (two regions)

Given `φ ∈ R^V`.

### 3.1 Face curvature `[formula: cut-face-curvature]`

```
θ̃_i^f = θ_i^f · 2π/Θ_i   (i interior),      θ̃_i^f = θ_i^f   (i ∈ B_M)
Ω_f   = θ̃_i^f + θ̃_j^f + θ̃_k^f − π                 integrated Gaussian curvature of face f = ijk
k_f   = Ω_f / A_f                                   density
```

Positive on a sphere (D1). Boundary vertices carry no Gaussian curvature (D9a).
`Σ_f Ω_f = Σ_{i∉B_M} (2π − Θ_i)` `[formula: cut-face-curvature-sum]`; on a closed mesh this is `2πχ`.
The identity needs every vertex to lie in a face (`Θ_i > 0`), which `TriMesh::new` guarantees.
A curvature override replaces `k_f` by a caller-supplied density (`Uniform(c)` or one value per face).

Because the corner angles of a Euclidean triangle sum to `π`, the same number is
`Ω_f = Σ θ_i^f (2π − Θ_i)/Θ_i` over the interior corners `i` of `f`. Implementations may evaluate
this form: it has no cancellation where the surface is flat, and differs from the first by rounding
only (D13e).

### 3.2 Regions and crossings `[formula: cut-edge-crossing]`

`region(i) = 0` if `φ_i ≤ 0`, else `1`. An edge `[a, b]` with `region(a) ≠ region(b)` is *crossed*;
its cut vertex sits at parameter `t` from `a` (point `(1−t)·a + t·b`):

```
t = φ_a / (φ_a − φ_b),   then moved into [δ, 1−δ]
```

`δ` is the option `min_crossing_offset` (default `1e-3`, `0 ≤ δ < ½`); the number of crossings that
were moved (`t < δ` or `t > 1 − δ`) is reported (D9b). With `δ = 0` nothing is moved, and a crossing
on a vertex is an error, because its triangle has no area. Vertices of crossed edges are the
*cut-adjacent* vertices `Adj`. A face has zero or
two crossed edges; in the second case it carries one *cut segment* joining its two cut vertices.
The cut `γ` is the union of the segments: closed polylines, or polylines ending on `∂M`.

### 3.3 Cut mesh

Every crossed face `abc` with `a` alone on its side, cut vertices `p ∈ ab`, `q ∈ ac`, is replaced by
the triangle `(a, p, q)` and the quad `(p, b, c, q)`; the quad is split by the shorter of its
diagonals `p–c`, `b–q` (tie: `p–c`; in floating point `p–c` iff `|pc|² ≤ |bq|²`). The triangles of
a crossed face are listed in the order `(a, p, q)`, then `(p, b, c)`, `(p, c, q)` or `(p, b, q)`,
`(b, c, q)`. Orientation is inherited. Sub-triangles get the region of their original vertices, the
density `k_f` and (for scaled Dirichlet) linearly interpolated `a` of their parent (`a` at a cut
vertex is `(1−t) a_x + t a_y` with the moved `t`). On symmetric meshes the two diagonals can be
equal up to rounding; the number of faces where they agree to a relative `1e-9` is reported, and
two implementations may differ there at the level of the discretization error (D13c).

Lengths inside a face come from barycentric coordinates `[formula: cut-barycentric-length]`: for two
points with barycentric difference `d = (d_a, d_b, d_c)`, `Σ d = 0`,

```
|x − y|² = −( d_a d_b ℓ_ab² + d_b d_c ℓ_bc² + d_c d_a ℓ_ca² )
```

and areas from `[formula: cut-subtriangle-area]`: `area = |det(barycentric coordinates)| · A_f`.

The angles of a sub-triangle are computed in the plane of its face, not from its three edge
lengths: with the inner product of two barycentric differences `d`, `e` (the polarization of the
formula above, same id)

```
⟨d, e⟩ = −½ [ (d_a e_b + d_b e_a) ℓ_ab² + (d_b e_c + d_c e_b) ℓ_bc² + (d_c e_a + d_a e_c) ℓ_ca² ]
```

the cotangent at the corner `x` of a sub-triangle `xyz` is `⟨y − x, z − x⟩ / (2 · area)`. Nothing
is subtracted that could cancel, and a sub-triangle is valid as soon as its determinant is positive
(D13a). Uncut faces keep the cotangents and areas of `M`.

The *submesh* `S_r` of region `r` consists of its (sub-)triangles. Its vertices are the original
vertices of the region and the cut vertices (each cut vertex belongs to both submeshes; the cut
mesh that is handed out has one copy of it per region, §7). `B_S` =
cut vertices ∪ (`B_M` ∩ `S`); `I_S` = the rest. A submesh may have several connected components
and any topology.

### 3.4 Scale factor `[formula: cut-yamabe-system]`

With `L` the cotan Laplacian of `S` and `g_i = Σ_{t∋i, t⊂S} k_t A_t / 3` (the paper's `W_F K`):

```
u_B = 0,        L_II u_I = −g_I
```

A component of `S` with no boundary vertex has no solution; this happens only for a closed
component of `M` that the cut does not touch. It is reported and left out of everything (D9c):
there `u = p = 0` (its rows of the system are those of the identity), `σ = 0`, it contributes
nothing to `J`, and its area is left out of `Area(S_r)` in the energy of §6. The report lists these
components and their area per region.

Remark for implementers: `I_S` over both regions is exactly the interior vertices of `M`, and the
pattern of the block-diagonal matrix is a subset of the pattern of `M`'s Laplacian, so one symbolic
factorization serves the whole flow, and one numeric factorization per step serves §3.4, §3.6, §3.8.

### 3.5 Normal derivative `[formula: cut-normal-derivative]`, `[formula: cut-flux-recovery]`

For `i ∈ B_S`, `v_i := ∮ ∂u/∂n · (hat function of i) ds` with `n` **outward** from `S`:

```
v_i = (L_BI u_I)_i + g_i
ℓ_i = ½ Σ (lengths of the boundary edges of S at i)
```

`v` is determined by `L u = −g + v` alone: it is the residual of the equation at the rows that
were not solved. (The paper's `v` in Eq. 26, 28 has the opposite sign, D6.) At a cut vertex inside `M` both regions
have the same `ℓ`; at a cut vertex on `∂M` the two boundary edges are a cut segment and a piece of
the boundary edge of `M` (D9d).

**Flux recovery.** The paper takes `(∂u/∂n)_i ≈ v_i/ℓ_i` at every cut vertex on its own. That is
not robust: `v_i` is tested with a hat function that reaches one mesh width into the region, so its
error is of the order of the mesh width whatever `ℓ_i` is, and where two or more consecutive
crossings lie close to one mesh vertex `ℓ_i` is as small as their distance. The quotient is then
wrong by a factor of up to `1/δ` (D12). Sums of `v` over a stretch of the cut are accurate, so the
normal derivatives `d_p ≈ (∂u/∂n)_p` at the cut vertices `C` of `S` are computed together. For every
cut segment `σ = pq` of `S`, lying in the face `f = abc` of `M`:

```
κ_σ = c² · h_f² / |pq|,       h_f² = (ℓ_ab² + ℓ_bc² + ℓ_ca²)/3   (scaled unit),      c = 0.5

ℓ_p d_p + Σ_{σ = pq at p} κ_σ (d_p − d_q) = v_p          for every p ∈ C
```

This is one linear system per region with the matrix `diag(ℓ) + Σ_σ κ_σ [[1,−1],[−1,1]]`:
symmetric, positive definite, an M-matrix; it decouples over the connected pieces of the cut. It is
one backward Euler step of the 1D heat flow on the cut, with lumped mass `ℓ` and the local time
`(c h_f)²`, applied to `v/ℓ`. Consequences: `Σ_p ℓ_p d_p = Σ_p v_p` on every piece of the cut (the
total flux is kept), `d` stays between the smallest and the largest `v_p/ℓ_p` of the piece, and cut
vertices much closer to each other than `c h_f` get one common value. `c` is the option
`FluxRecovery::Smoothed { length }`; `FluxRecovery::Nodal` is the paper's `d_p = v_p/ℓ_p`, which is
also `c = 0`.

The same map applied to `q` (§3.6) gives `e_p ≈ (∂p/∂n)_p`. At vertices of `B_M` that are not cut
vertices no derivative is needed; `v_i/ℓ_i` is reported there.

### 3.6 Energy, adjoint, gradient density

```
Dirichlet   J = uᵀ L u = −u_I · g_I                                   [formula: cut-dirichlet-energy]
            G_i = d_i²                                                 [formula: cut-gradient-dirichlet]

scaled      J = uᵀ L_a u,   L_a = Σ_t ā_t · (cotan Laplacian of the single triangle t),
                            ā_t = (a_i + a_j + a_k)/3            [formula: cut-weighted-laplacian]
            L_II p_I = r_I,   r_i = −2 (L_a u)_i            (i ∈ I)   [formula: cut-adjoint-scaled]
            q_i = (L_BI p_I)_i − r_i,   r_i = −2 [ (L_a u)_i − a_i v_i ]   (i ∈ B)
            G_i = −a_i d_i² − d_i e_i                                  [formula: cut-gradient-scaled]

pointwise   J = Σ_{i∈S} A_i^S b(u_i),   A_i^S = lumped vertex area within S
            L_II p_I = r_I,   r_i = −A_i^S b′(u_i)                     [formula: cut-adjoint-pointwise]
            q_i = (L_BI p_I)_i − r_i     (i ∈ B, where u_i = 0)
            G_i = b(0) − d_i e_i                                       [formula: cut-gradient-pointwise]
```

`d_i`, `e_i` are the normal derivatives of `u` and `p` at the cut vertex `i` recovered from `v` and
`q` as in §3.5. `G` is needed at cut vertices only. The discrete adjoint `p` is the exact
derivative of the discrete objective `J` with respect to the source `g`, at a fixed mesh and a
fixed set of interior vertices: `∂J/∂g_i = p_i` for `i ∈ I` (with `u_I = −L_II⁻¹ g_I`). It is
not a derivative with respect to the cut, which changes the mesh and the set `I`; the step from
`p` to `G` is the smooth formula of §1.

`L_a` is the exact stiffness matrix of `∫ a ∇f·∇g` for piecewise linear `a, f, g`, hence positive
semidefinite for `a ≥ 0` on every mesh (D11). `q` is the integrated outward normal derivative of
`p`, defined like `v` (the extra `W_V` of the paper's Sec. 5.6 is not there, D6). With `a ≡ 1`: `p = −2u`, `q = −2v`, and the scaled density is
the Dirichlet one.

Curvature-based weight `[formula: cut-feature-size]` (paper Sec. 3.5.1, 5.6), computed once on `M`:

```
K_i  = (2π − Θ_i)/A_i  (interior),  0 on B_M
r⁰_i = 1 / ( √|K_i| + 1/r_ε )                (D5)
(W_V + c L⁺) r = W_V r⁰,   c = 0.02
a_i  = r_i³
```

`L⁺` is the Laplacian of the edge weights `max(w_ij, 0)`, `w_ij = ½ Σ cot` over the one or two
faces at the edge `ij` (D17). `W_V + c L⁺` is then an M-matrix on every mesh, so
`min r⁰ ≤ r_i ≤ max r⁰`, `a_i > 0`, and `Σ A_i r_i = Σ A_i r⁰_i`. On a mesh without negative edge
weights `L⁺ = L`. A caller-supplied `a` is checked for `a_i > 0`.

Region totals: `G_r,i = w_r · G_i + β_r`. A region with `w_r = 0` needs no solve; implementations
may solve it anyway and report its `u` (nothing else changes, D13d).

### 3.7 Speed on the cut `[formula: cut-speed]`

At every cut vertex `c`: `s_c = G_0,c − G_1,c`.

Sliver smoothing `[formula: cut-sliver-smoothing]` (paper Sec. 6.2, D7) is optional and off by
default (`passes = 0`): with the flux recovery of §3.5 it changes nothing measurable, and without
it, it does not repair what it was meant to repair (D12). A cut vertex is a *sliver vertex* if an
incident cut-mesh triangle (of either region) has aspect ratio `> ρ₀` (default 5, any `ρ₀ > 0`),
where aspect ratio is `(longest edge)² · √3 / (4 · area)` (1 for an equilateral triangle). For
`passes` rounds, every sliver vertex simultaneously takes the arclength-weighted mean of its two
neighbours along the cut, `s_c ← (ℓ₂ s₁ + ℓ₁ s₂)/(ℓ₁ + ℓ₂)`, with `ℓ_k` the length of the segment
to neighbour `k`. A cut vertex on `∂M` (one neighbour) is left alone. The procedure of the paper is
`FluxRecovery::Nodal` with `passes = 1` (`EvalOptions::paper()`).

### 3.8 Extension `[formula: cut-extension-local]`, `[formula: cut-harmonic-extension]`

A rate `r_p` at the cut points (§3.9, §3.10) becomes a rate `σ_i` at every vertex of `M`.

**From the cut (default, `Extension::Local`, D15).** With `ℓ_p` the dual length of cut point `p`
on the cut (half the lengths of its segments) and `β_ip` the barycentric weight of the cut-adjacent
vertex `i` at `p` (`1 − t` or `t` on the edge of `p`):

```
σ_i = Σ_p ℓ_p β_ip r_p / Σ_p ℓ_p β_ip              over the crossed edges at i
σ_i = Σ_j (σ_j / ℓ_ij) / Σ_j (1 / ℓ_ij)             over the neighbours j of i in its region that
                                                    are one ring closer to the cut
```

The ring of a vertex is the number of uncrossed edges on the shortest path to a cut-adjacent vertex
of its region; a vertex that no such path reaches gets `σ = 0`. No linear system. `σ` lies between
the smallest and the largest `r`. The first line is the transpose of the interpolation
`Σ_i β_ip σ_i` that the explicit step applies at the cut points, normalized, so
`Σ_p ℓ_p Σ_i β_ip σ_i = Σ_p ℓ_p r_p`: the mean rate of the cut is the mean of `r`.

**Harmonic (paper Sec. 5.3.1, `Extension::Harmonic`).** In each submesh, with `r` as Dirichlet data
on the cut vertices `C` and all original vertices `O` of the submesh unknown (so `∂M` gets the
natural, zero-flux condition, D9e):

```
L_OO σ_O = −L_OC r_C
```

A component without cut vertices gets `σ = 0`. On a closed mesh `L_OO` is the matrix of §3.4 and
its factorization is reused.

### 3.9 Explicit step `[formula: cut-explicit-step]`, `[formula: cut-gradient-norm]`

`φ̃_i = φ_i + τ σ_i`, with `σ` the extension of `ḡ · r`. `ḡ` is 1 in the paper
(`SpeedScaling::Unit`: `φ_t = s` on the cut, exact while `|∇φ| = 1` there). By default
(`SpeedScaling::MeanGradient`, D16) `ḡ` is the mean of `|∇φ|` over the connected piece of the cut
that the cut point lies on: in a crossed face whose lone vertex `a` is cut off by the segment `pq`,
`|∇φ| = |φ_a| · |pq| / (2 · area(a, p, q))`, and `ḡ = Σ |pq| / Σ (|pq| / |∇φ|)` over the piece.

### 3.10 Length and penalty step `[formula: cut-length-step-on-cut]`, `[formula: cut-length-step]`

`w_i = α_L + η_i ≥ 0` at the vertices, interpolated linearly to the cut points. `w ≡ 0` skips the
step (then `r = s`).

**On the cut (default, `LengthStep::OnCut`, D14).** The integrated length force at cut point `p`:

```
Θ_p = w_p · θ_p + Σ_{σ = pq at p} ½ |pq| (∂w/∂n)_σ
```

* `θ_p = π − β_p` is the turning angle of the cut at `p`, with `β_p` the interior angle of the
  lower region at `p` (the sum of the corner angles of its cut-mesh triangles at its copy of `p`).
  At a cut point on `∂M` the cut is reflected in the boundary (the zero-flux condition):
  `Θ_p = w_p (π/2 − β_p)`. The angle itself is used, not `2 sin(θ_p/2)`: angles add up along the
  cut to its change of direction, so the force on a run of cut points crowded around a mesh vertex
  is the turning of the cut across the run, and a closed cut on a flat surface loses area at the
  exact rate `2π w`.
* `(∂w/∂n)_σ` is the derivative of the linear `w` of the face of `σ` across the segment, out of
  the lower region: with `a` the lone vertex of the face, `h_a` its height over `pq` and `f` the
  foot of that height (at the fraction `cot_p/(cot_p + cot_q)` of `pq` in the triangle `(a, p, q)`),
  `(w_f − w_a)/h_a`, negated if `a` is in the higher region. `PenaltyRule::PaperDiagonal` drops
  this term (D8).

The force per unit length `λ⁰` is recovered exactly like the normal derivative (§3.5, same `c`),
and the rate `r` takes the length force at the end of the step and `s` at its start (moving `p`
out by `δ` relative to its neighbours turns the cut by `Σ (δ_p − δ_q)/|pq|` more):

```
ℓ_p λ⁰_p + Σ_{σ = pq at p} ((c h_σ)² / |pq|) (λ⁰_p − λ⁰_q) = Θ_p
ℓ_p r_p  + τ Σ_{σ = pq at p} (w̄_σ / |pq|) (r_p − r_q)       = ℓ_p (s_p + λ⁰_p),     w̄_σ = ½ (w_p + w_q)
```

Both matrices are symmetric positive definite M-matrices with row sums `ℓ_p`:
`Σ ℓ_p λ⁰_p = Σ Θ_p`, `Σ ℓ_p r_p = Σ ℓ_p (s_p + λ⁰_p)`, `r` lies between the smallest and the
largest `s + λ⁰`, and a cut with `s + λ⁰ = 0` does not move, whatever `τ`. In the notation of §1 a
ripple is multiplied per step by `1 + τ (ωk − wk²)/(1 + τ w k²)`, which lies in `(0, 1)` for every
`k > ω/w` and every `τ`.

**On the mesh (paper Fig. 9 step 6, `LengthStep::OnMesh`).** After the explicit step with `r = s`,
with `L_w = Σ_f w̄_f · (cotan Laplacian of the single face f)`, `w̄_f` the mean of `w` over the
three vertices of `f` (the weighted Laplacian of §3.6 on `M`):

```
(W_V + τ L_w) φ′ = W_V φ̃
```

This is backward Euler for the second term of (★); the matrix is symmetric positive definite, and
is factorized once and reused while `τ` and `w` do not change. With `η = 0` it is the paper's step
(D4). The paper's own penalty rule (Sec. 7.1), `PenaltyRule::PaperDiagonal`
`[formula: cut-penalty-step-paper]`, is

```
(W_V D⁻¹ + τ L) φ′ = W_V D⁻¹ φ̃,   D = diag(α_L + η_i) > 0
```

which drops the `∂η/∂n` term of the gradient (D8). The default rule is the gradient rule.

### 3.11 Symmetry `[formula: cut-symmetry-projection]`

If vertex classes are given (paper Sec. 7), every vertex takes the arithmetic mean of `φ′` over its
class, after §3.10 and again after §3.12.

Without classes a symmetric state stays symmetric to rounding, unless the symmetry reverses the
orientation of the faces and the diagonals of a cut quad tie (§3.3, D13c): the tie is broken by
the roles of the corners, the two halves are triangulated differently, and the symmetry holds to
the discretization error only (measured `1.2·10⁻⁵` after 10 steps on icosphere(3) for two latitude
circles about the axis of the mesh; `4·10⁻¹⁶` about a generic axis). With classes it is exact.

### 3.12 Redistancing `[formula: cut-redistance]`

Regions and cut-adjacent vertices `Adj` are recomputed from `φ′`. Signs stay. Vertices in `Adj` keep
their values (paper Sec. 6.4). For the others, `φ_i = sign_i · T_i` where `T` is computed by fast
marching on `M`, restricted to each sign class:

```
seeds:      T_i = |φ′_i| for i ∈ Adj; all seeds are accepted from the start and never change
tentative:  T_i = min of
              min_{λ∈[0,1]} [ λ T_j + (1−λ) T_k + | x_i − (λ x_j + (1−λ) x_k) | ]
                  over triangles (i,j,k) with j, k accepted and of i's sign class
              T_j + ℓ_ij    over accepted neighbours j of i's sign class
accept:     repeatedly the not yet accepted vertex with the smallest tentative value
            (ties: the lower vertex index); its value is final from then on
```

The inner minimum has a closed form (a convex function of `λ` in the plane of the triangle; the
edge terms are `λ = 0, 1`). A vertex that no front reaches (its sign class has no cut-adjacent
vertex in its component) keeps `φ′_i`, and is counted in the report.

On a mesh without obtuse angles this is the unique solution of the scheme
`T_i = min over all incident triangles of the same class`. On obtuse triangles the single pass and
that fixed point differ, and iterating to the fixed point is not affordable (measured in
`ddg_core::distance`: thousands of re-openings per vertex at an angle of 160°), so the single pass
is the definition. It is the paper's fast marching (Kimmel, Sethian 1998) with the local update
written in Hopf–Lax form (Bornemann, Rasch 2006) and without unfolding of obtuse triangles (D10).

## 4. One flow step

```
                                                      default                 paper (selectable)
1  regions, crossings, cut mesh                       §3.2 §3.3               same
2  u, v, d (and p, q, e) in each region               §3.4–3.6                d = v/ℓ
3  s on the cut                                       §3.7                    + sliver smoothing
4  length force Θ, λ⁰, the step τ, the rate r         §3.10 on the cut        r = s
5  σ on V: extension of ḡ·r                           §3.8 from the cut, §3.9 ḡ of the piece
                                                                              harmonic, ḡ = 1
6  φ̃ = φ + τσ                                        §3.9                    same
7  —                                                                          (W_V + τ L_w) φ′ = W_V φ̃
8  symmetry                                           §3.11                   same
9  redistance, then symmetry                          §3.12, §3.11            same
```

`FlowOptions::new` is the left column, `FlowOptions::paper` the right one; each of the choices
(`LengthStep`, `Extension`, `SpeedScaling`, `FluxRecovery`, `PenaltyRule`) is an option of its own.

**Step size** `[formula: cut-step-cap]`. `ρ` is the speed of the cut that is known before the
step: `s + λ⁰` with the length step on the cut (it bounds `|r|`), `s` with the length step on the
mesh. `h` is the shortest crossed edge, `ω` the roughening rate of the state (§1).

```
Fixed(τ)                                    τ                                        (the paper's)
Capped { tau_max, max_motion }              min(tau_max, max_motion · h / max|ρ|)    (paper Sec. 7.4)
Adaptive { time_scale, max_motion, tau_max} min(tau_max, time_scale / ω, max_motion · h / max|ρ|)
```

The paper gives no values (D9f). The default is `Adaptive { 0.05, 0.25, 1 }`. Why:

* The distortion term is explicit, and a cut reacts to a displacement at rates of the order `ω`
  per unit of length of the unit-area surface (sphere: the equator returns at the rate
  `ω/R = 3.5 ω`). `time_scale = 0.05` keeps `τ ×` that rate at `0.18`. A fixed `tau_max` cannot do
  this: `ω` at the same cut of spot is `3.4·10³` for the Dirichlet energy, `1.5·10²` for Hencky,
  `17` for the scaled Dirichlet energy with `r_ε = 1`.
* The error is first order in `τ`. Circle of §6 item 4: relative error of the radius `+3.9·10⁻²`,
  `+1.1·10⁻²`, `+2.9·10⁻³`, `+9·10⁻⁴` for 6, 24, 96, 384 steps, independent of the mesh; it is the
  error of the explicit Euler method for the radius, `½ α_L τ ln(ρ₀/ρ)/ρ²`. Latitude circle of §6
  item 5 with `α_L = 8`: error of `ψ` `+2.1·10⁻³` and `+5.3·10⁻⁴` for 50 and 200 steps. Equilibria
  do not depend on `τ` (§3.10).
* `max_motion = 0.25`: in the bump experiment (README, Fig. 5; Hencky energy, three length
  weights, `Capped { 10³, max_motion }`) none of the three runs came to rest in 3000 steps with
  `0.5`, one did with `0.25`, all three with `0.1`; the energy at the end was the same to three
  digits. A step never moves the cut by more than `max_motion · h`.

**Length weight.** No value of `α_L` is right for every energy: `E` is not homogeneous in it
(D9f). What can be fixed is the shortest ripple that the distortion term may amplify:
`α_L = ω ℓ/(2π)` `[formula: cut-roughening-rate]` (`length_weight_for_ripples`) makes ripples
longer than `ℓ` grow and shorter ones decay (§1). Verified on the sphere, where the rates are known
exactly (§6 item 5). `ℓ` must be several mesh edges, or the flow amplifies the noise of the mesh at
the rate `ω k`; with `α_L = 0` it does so at every scale, and no step size helps. The example
program defaults to `ℓ = 0.25` with `ω` of the initial cut; a cut that starts far from the
curvature it will reach has a small `ω` there, and the weight must then be set from a state that
matters (README, Fig. 5).

**Feature size cap.** `r_ε = 1` (the size of the unit-area surface) is the default of the example
program: on spot, the ratio of the largest to the median `|s|` over a planar cut (the unevenness
that the rescaling is meant to remove; Dirichlet: 35 and 218 for two planes) is 34, 32, 29, 26, 21,
16 and 212, 205, 199, 211, 211, 141 for `r_ε = 0.01, 0.03, 0.1, 0.3, 1, 10`: a small `r_ε` turns
the rescaling off, as the formula says (D5), a large one makes `a = r_ε³` in flat regions outweigh
everything else. The Hencky energy is more even than any of them there (23 and 14).

**Stopping.** A run stops at `max_steps`; or when the cut motion of a step (largest
`|φ_new − φ_old|` over vertices cut-adjacent before or after the step) is at most
`motion_tolerance`; or when the energy has fallen by at most `energy_tolerance · |E|` over the last
`energy_window` steps; or when the cut has vanished (one region is empty). The reason is in the
report. Defaults: 1000 steps, `10⁻⁵`, `10⁻⁴` over 25 steps (D9f). The second criterion is not in
the paper. It is needed: on an irregular mesh `s + λ⁰` has a noise floor of the size of the two
terms that balance (the turning angles of a polyline through a triangle mesh, and `s` near
vertices), the cut keeps moving by a few percent of an edge per step, and the motion never falls
below a small tolerance while the energy stands still (spot, scaled Dirichlet: energy constant to
`3·10⁻⁵` from step 150 on, cut motion `1.3·10⁻³` per step for 1350 further steps). On smooth
meshes the first criterion ends the run (sphere, bump).

Every step reports the energy terms of §6 for its start state, `τ`, the largest `|s|` and `|λ⁰|`,
`ω`, the cut motion, `|∇φ|` in the crossed faces (minimum, mean, maximum), the number of cut
points, moved crossings, closed components without a cut, vertices the redistancing did not reach,
the residuals of the solves, and the time of each stage.

Between steps the caller may change `α_L`, `η`, the energy, `w_r`, `β_r`, the curvature, the
symmetry classes and every option (the schedules of Sec. 7.4). The state is `φ`; a flow can be
resumed from a stored `φ`.

## 5. Steering: what the caller controls

| Paper | Control |
|---|---|
| Eq. 5 | `α_L` length weight, may change between steps |
| Sec. 3.2, 3.5 | energy: Dirichlet, scaled Dirichlet (`a` curvature-based or supplied), Hencky, pointwise `b, b′` |
| Sec. 7.1 | `η` per vertex (the library takes numbers; a helper computes the visibility penalty `η_i = ν_i · n_i·E_i` of the paper from positions and a camera point) |
| Sec. 7, Fig. 2, 3, 12 | symmetry classes |
| Fig. 14 | per-region distortion weight `w_r` and area weight `β_r` |
| Sec. 7.4, Fig. 17, 18 | curvature override: uniform or per face |
| Sec. 6, 7.4 | step size rule, stop rule; where the length step acts, how the rate is extended and scaled (§3.8–§3.10) |
| Fig. 9 input | initial `φ⁰`: any values, as they are or redistanced; helpers (`init`) for a geodesic circle around a vertex, two-valued vertex labels, two-valued face labels (a vertex takes the label that covers more of its angle; a tie gives label 0). Clusters of normals are labels the caller computes; with more than two they belong to §8 |
| (workspace) | fixed seams: cut the mesh along them first, they are then part of `∂M` |

## 6. Discrete energy and verification targets

```
E(φ) = Σ_r [ w_r J_r + β_r Area(S_r) ] + Σ_segments ½ (w_p + w_q) |pq|
```

(`w` interpolated linearly to the cut vertices with the moved `t`; every segment once; `Area(S_r)`
without the components of §3.4 that have no solution). The flow is a discretization of the gradient
flow of the smooth energy, not the exact gradient flow of `E(φ)`; tests therefore use tolerances
that shrink under refinement.

1. `Σ_f Ω_f` identity (§3.1) to rounding; areas of the cut mesh sum to 1; Euler characteristics.
2. Unit sphere scaled to unit area (`R = 1/√(4π)`), cut along the latitude circle of polar angle
   `ψ₀`: in the cap, `u = 2 ln( cos(ψ₀/2) / cos(ψ/2) )`, `∂u/∂n = tan(ψ₀/2)/R`,
   `J = −8π ln cos(ψ₀/2) − 4π sin²(ψ₀/2)`. Convergence under icosphere refinement.
3. Finite differences: for a smooth `ψ`, `E(φ+εψ) − E(φ−εψ) ≈ −2ε ∫_γ [s + (length terms)] ψ ds`.
   This is the test of every sign in §1. For a `φ` that is not a signed distance, `−ε ψ ds` is
   replaced by the area `dA_c` that region 0 gains at cut vertex `c`, which is exact on the cut mesh
   (in a face of area `A_f` whose lone vertex has the crossings at the fractions `s_p`, `s_q` of
   its two edges, the corner triangle has the area `s_p s_q A_f`): `dE ≈ Σ_c s_c dA_c` for the
   distortion and area terms. Measured (`tests/finite_differences.rs`): the two sides agree in sign
   for every energy, for each region alone and with `w_r`, `β_r`, and their difference, relative
   to `Σ_c |s_c dA_c|`, falls by 3 to 4 per refinement (sphere, levels 2–6: `5·10⁻²` to
   `3·10⁻⁴`; bumpy disk with the cut ending on `∂M`, grids 16–128: `2·10⁻¹` to `4·10⁻³`).
4. Flat disk, circular cut of radius `ρ₀`: `u = 0`, `s = 0`, and the circle shrinks as
   `ρ(t) = √(ρ₀² − 2 α_L t)`. Measured (`examples/flow_study.rs circle`: unit square, cells
   jittered by 0.3, `ρ₀ = 0.3`, `α_L = 1`, `t = 0.03`, exact radius `0.17321`), relative error of
   the mean radius for 6, 24, 96, 384 steps:

   ```
   grid   default stages                               stages of the paper (recovered flux)
    32    +3.8e-2  +9.6e-3  +2.0e-3  −4e-5             −16 %  −17 %  −17 %  −20 %
    64    +3.9e-2  +1.1e-2  +2.9e-3  +8.6e-4           −12 %  −11 %  −11 %  −15 %
   128    +3.9e-2  +1.1e-2  +3.0e-3  +9.0e-4           −9.4 % −7.5 % −6.4 % −9.1 %
   256    +4.0e-2  +1.1e-2  +2.9e-3                    −8.3 % −5.8 % −4.0 %
   ```

   The default stages are first order in `τ` and their error does not depend on the mesh; the
   cut stays round to `5·10⁻³` (grid 32) – `5·10⁻⁴` (grid 256) of its radius. The stages of the
   paper do not converge to the law in `τ` and converge slowly in `h` (D14).
5. (a) Sphere of unit area, latitude circle of polar angle `ψ`: `G_0 − G_1 = −4 cos ψ/(R² sin² ψ)`
   and the curvature of the circle is `cot ψ/R`, so
   `ψ′ = cos ψ/(R² sin ψ) · (4/(R sin ψ) − α_L)`; for `α_L = 0`,
   `F(ψ(t)) = F(ψ₀) + 4t/R³` with `F(ψ) = ln(sec ψ + tan ψ) − sin ψ`. Equilibria: the equator, and
   for `α_L R > 4` the latitude `sin ψ* = 4/(α_L R)`. Measured (`flow_study latitude-ode`), error
   of the mean `ψ` on icosphere(3), (4), (5), with 50 and 200 steps:

   ```
   α_L  ψ₀    t       exact ψ   default stages                                   paper stages (level 5)
    8   0.8   4e-3    1.11512   +2.4e-3 +8e-4 | +2.2e-3 +5e-4 | +2.1e-3 +5.3e-4   −2.2e-3  −7.5e-3
    8   0.8   1.2e-2  1.35078   +4.8e-3 +9e-4 | +4.8e-3 +1.2e-3 | +4.7e-3 +1.1e-3 −5e-4    −6.0e-3
   20   0.5   4e-3    0.73163   +5.8e-3 +4.3e-3 | +7e-4 −7e-4 | +1.7e-3 +3.7e-4   −4.5e-2  −3.0e-2
   20   1.1   2e-2    0.79701   +2.2e-2 +2.3e-2 | +4.4e-3 +5.0e-3 | −2.7e-4 +2.8e-4 −8.5e-2 −8.5e-2
    0   0.8   1e-3    0.98303   +1.3e-4 −5.9e-4 | +2.4e-4 −3.8e-4 | (+0.07, +0.11: noise, see below)
   ```

   With `α_L = 20` the flow from `ψ = 0.5` and from `ψ = 1.1` reaches `0.7843` and `0.7930`
   (icosphere(5), `ψ* = 0.78805`); the stages of the paper stop at `0.7393` and `0.7409` there and
   at `0.6535` and `0.7273` on icosphere(3). With `α_L = 0` the flow is ill-posed (§1): the test
   is possible on a coarse mesh over a short time only (icosphere(3), `t = 10⁻³`: `10⁻⁴`–`6·10⁻⁴`);
   on icosphere(5) the noise of the mesh has grown to a spread of `0.1`–`0.2` in `ψ` by then with
   the default stages (`3·10⁻³` with the stages of the paper, which damp short ripples more).

   (b) Ripples of the equator, `ε cos(mθ)`: growth rate `(m − 1)[4/R³ − α_L (m + 1)/R²]` (§1).
   Measured (`flow_study modes`, `ε = 0.02`, 160 steps over the time `0.7/|rate|`), default stages,
   icosphere(4) and (5), against the exact rate:

   ```
   α_L     m = 2                  m = 3                  m = 4                   m = 6
   0       +172  n/a   (+178)     +343  +372  (+356)     +520  +529  (+535)      +839  +871   (+891)
   2       +112  +101  (+103)     +156  +150  (+155)     +160  +155  (+158)      +41   +10    (+11)
   4       +53   +19   (+27)      −10   −49   (−46)      −153  −207  (−219)      −730  −833   (−868)
   4.726   +34   −8    (0)        −63   −113  (−119)     −273  −341  (−356)      −1016 −1144  (−1188)
   6       +3    −52   (−48)      −157  −228  (−247)     −486  −577  (−596)      −1520 −1690  (−1748)
   10      −78   −174  (−199)     −485  −612  (−649)     −1157 −1317 (−1350)     −3103 −3405  (−3507)
   ```

   (`n/a`: with `α_L = 0` the noise of icosphere(5) outgrows the mode `m = 2` over this time.)
   The damping by the length term is 89 % of the exact one on icosphere(4) and 98 % on (5)
   (`m = 4`): second order in the mesh width, about `1.1 (kh)²` low. With the stages of the paper
   the equator turns stable at `α_L ≈ 4.0` instead of `4.73` on icosphere(5) (rate of `m = 2` at
   `α_L = 4.726`: `−28`). A great circle in general position with `α_L = 6` stays within
   `10⁻³` of itself over 20 steps of `2·10⁻⁴` on icosphere(4) (`|s| ≤ 0.27`, 2 % of `G_0`).

   (c) The energy decreases in every step for small steps (`τ ≤ 0.02/ω`, motion `≤ 0.1 h`), in
   each of the four energies, on a sphere and on a disk whose boundary the cut ends on. With the
   default step it rises in up to a tenth of the steps, by at most `6·10⁻⁴` of the total drop
   (sphere, README Fig. 1). (d) Symmetric input stays symmetric (§3.11). (e) One step does not
   depend on the numbering of vertices and faces (`10⁻¹⁶`) and scales with the mesh as stated in
   §2 (`10⁻¹⁵`), with the default stages and with the paper's. (f) A straight cut across a linear
   penalty `η = A x` on a flat square moves down the penalty at the speed `A` with the gradient
   rule (`x = 0.5000` after `t = 0.05` from `0.6` with `A = 2`; `0.5033` with the stages of the
   paper) and stays where it is with the paper's rule (D8). (g) One-sided distortion with an area
   reward (paper Fig. 14): a cap on the sphere gains
   `tan²(ψ/2)/R² + β₀ + α_L cot ψ/R` per unit of area; with `α_L = 2` and `β₀ = −8.638` the stable
   equilibrium is `ψ* = 1.2`, and caps started at `1.0` and `1.4` end at `1.165`, `1.222`
   (icosphere(3)), `1.202`, `1.205` (4), `1.195`, `1.201` (5).
6. One step agrees with the Lean oracle on small meshes: the oracle implements the stages of the
   paper with the flux recovery of §3.5, and `tests/oracle.rs` replays its 26 fixtures with
   `LengthStep::OnMesh`, `Extension::Harmonic`, `SpeedScaling::Unit`: evaluation (70907 numbers)
   to `10⁻¹²`, `d, e, G, s` (6218) to `3.4·10⁻¹²`, every stage of the step (24444) to
   `1.6·10⁻¹³`. The default stages of §3.8–§3.10 are not in the oracle.
7. (Milestone 2) `N = 2` reproduces the two-region code.

## 7. Output

* the state `φ` (caller's unit) and the region of every vertex;
* the cut as points (`edge, t`) and segments (`point, point, face`), with helpers that return 3D
  positions when the caller has them, per point or as ordered polylines: every connected piece
  is closed or ends on `∂M` at both ends; open pieces come first, in the order of their
  lower-numbered end, and start there; closed pieces follow in the order of their lowest-numbered
  point, start there and go on to the lower-numbered of its two neighbours;
* the cut mesh as a `TriMesh` + `Geometry` with maps to parent faces, original vertices or cut
  points, and the region of every face. This is what a flattening consumes. It is the cut
  *surface*: a cut point has one vertex per region that meets there (two for two regions), so the
  regions are separate pieces and every vertex has one region; original vertices keep their index.
  Its geometry is in the caller's length unit;
* per evaluation, on the vertices of that cut mesh: `u`, `v`, `ℓ`, the recovered `∂u/∂n` (`d`),
  `p`, `q`, `∂p/∂n` (`e`), `G_r`; per cut point `s` before and after sliver smoothing; the terms of
  `E`; the counters (moved crossings, sliver vertices, near-ties of the quad diagonal, closed
  components without a cut and their area) and the residuals of the solves. All in the scaled
  unit.

## 8. N regions (milestone 2)

Paper App. A. Everything in §3–§4 is the case `N = 2`; this section says what changes.

### 8.1 State

`φ_k(i)` for `k < N`, `N ≥ 2`: every `φ_k` is a signed distance to the boundary of region `R_k`,
negative inside. `region(i) = argmin_k φ_k(i)` (ties: the lowest `k`). A region may be empty or
have several components (the paper finds `N = 10` enough for any network, by colouring). A
two-region state `φ` is the state `φ_0 = φ`, `φ_1 = −φ`; on such a state every stage below must
give what §3–§4 give, for the default and for the paper's stages (§6 item 7).

### 8.2 Crossings, triple points, cut mesh `[formula: cut-edge-crossing-regions]`, `[formula: cut-triple-point]`

An edge `[a, b]` with `i = region(a) ≠ j = region(b)` carries one cut point, where the
interpolants of `φ^i` and `φ^j` are equal (paper Eq. 31), then moved into `[δ, 1−δ]`:

```
t = (φ_a^i − φ_a^j) / (φ_a^i − φ_a^j + φ_b^j − φ_b^i)
```

A third function that is smaller in between is ignored: one cut point per edge (App. A.2).

A face with two regions is cut as in §3.3. A face `abc` with three regions has three cut points
`p_ab`, `p_bc`, `p_ca` and a *triple point* `m` at their centroid; the cut consists of the three
segments from the cut points to `m`. The part of corner `a` is the quad `(a, p_ab, m, p_ca)`; it is
convex (`m` lies beyond the line `p_ab p_ca` seen from `a`, because `p_bc` does) and is split by
`a–m`. In the cut mesh a cut point has one copy per adjacent region (two), a triple point three.
Every segment knows its two regions.

### 8.3 Region solves and speed

§3.4–§3.6 hold per region. A triple point is an ordinary boundary node of each of its three
regions, with two boundary segments there; the flux recovery of §3.5 runs along the whole boundary
of a region, through triple points.

At a cut point `p` between regions `i < j`: `s_p = G_i,p − G_j,p`. It is the rate of `φ_i` there;
`φ_j` has `−s_p`. Triple points carry no rate of their own: they follow their three cut points.

### 8.4 The paper's stages (App. A.3)

* **Extension.** In region `k`, `e` is the harmonic extension (§3.8) of the rate of `φ_k` on the
  boundary of `R_k`: `±s_p` at cut points, and at a triple point the linear interpolation, in
  arclength, between its two neighbouring cut points on that boundary (not stated in the paper).
  For a vertex `i` of region `k`: `φ̇_k(i) = e_i` and `φ̇_j(i) = −e_i` for every `j ≠ k`. This is the
  paper's substitution `φ̇_j := −φ̇_k` on the boundary of `R_k`; where `R_k` borders `R_j` it is the
  true rate of `φ_j`.
* **Length step.** `(W_V + τ L_w) φ_k′ = W_V φ̃_k` for every `k`: one factorization, `N` right-hand
  sides.
* **Symmetry.** The class mean of §3.11 for every `φ_k`. The regions must be invariant under the
  symmetry; a symmetry that permutes regions is not supported.
* **Redistancing** `[formula: cut-region-projection]`. Regions are recomputed. `Adj_k` are the
  vertices of the edges with exactly one end in `R_k`, `Adj = ∪ Adj_k`.
  1. At every `i ∈ Adj` (Losasso et al. 2006, Sec. 3.1): with `m` the mean of the two smallest
     `φ_k′(i)`, `φ_k(i) = φ_k′(i) − m` for all `k`. The smallest is then `−d`, the second `+d`: one
     negative entry, and the two zero sets coincide.
  2. For every `k`: fast marching as in §3.12 with the seeds `Adj_k`, values `|φ_k|`, and the two
     classes "in `R_k`" and "not in `R_k`"; `φ_k = −T` in `R_k`, `+T` elsewhere. A vertex of `Adj`
     that is not in `Adj_k` is not a seed for `k` and is recomputed.
  For `N = 2` step 1 changes nothing and step 2 is §3.12.

### 8.5 The default stages

§3.8–§3.10 (length force on the cut, extension from the cut, mean-gradient scaling) are to be
carried over to networks. The rate at a cut point between `i < j` drives the difference
`φ_i − φ_j` at the two ends of its edge; the functions of the other regions at a cut-adjacent
vertex follow the substitution of §8.4, and redistancing is that of §8.4. What has to be derived is
the length force near a triple point: it is minus the first variation of `∫_γ w ds` over the
network, written so that (a) away from triple points it is §3.10, (b) three straight curves
meeting at 120° exert no force, and any other junction is driven toward that, (c) the total is
conserved by the recovery as in §3.10. The derivation goes here when it is done, with a D-entry
(the paper diffuses every `φ_k` on the mesh and has no junction rule).

### 8.6 Tests

* `N = 2` against the two-region code, both stage sets, to rounding.
* Pure length flow (`w_k = 0`): a honeycomb on a flat torus given by edge lengths alone is
  stationary, a perturbed one returns, and the length never increases; a triple junction relaxes to
  120°.
* With distortion: the tetrahedral partition of the sphere (4 regions, 4 triple points) and the
  cubical one (6 regions, 8 triple points) are stationary by symmetry; the energy decreases along
  the flow from generic `N`-region starts for every energy.
* Invariance under renumbering the vertices and under permuting the regions (up to the tie rules,
  which are to be stated); scaling as in §2.

### 8.7 Data and output

`w_k`, `β_k` per region; `α_L`, `η`, the energy and the curvature are shared. Helpers: a state from
`N`-valued vertex labels and from face labels; a clustering of normals as a start (the paper's
"greedy clustering of vertex normals", Sec. 7, is not specified further: say what is done). Output:
cut points, triple points (face and barycentric coordinates), segments with their two regions, the
curves between junctions as polylines, and the cut mesh with the region of every face.

## 9. Differences from the paper as printed

* **D1** Eq. 24 defines face curvature as `π − Σθ̃`, the negative of Eq. 23's convention (it sums to
  `−2πχ`); Eq. 28 uses it with the sign that makes the two cancel. Here `Ω_f = Σθ̃ − π` and the
  right-hand side is `−g`.
* **D2** Eq. 14–15, Fig. 8 and Sec. 6.2 write the speed as `σ_D⁺ − σ_D⁻` while Sec. 5.4 puts `M⁺`
  where `φ ≥ 0` and Sec. 4.1 sets `∇φ = n`. Read literally that is ascent. Here the region where
  `φ ≤ 0` enters with `+`. Verified by test 3 of §6.
* **D3** Eq. 4 has a factor ½ that Eq. 14, 22, 29 do not carry. Here `E_L` is the length of `γ`.
* **D4** Eq. 29 has `(I + τα_L L)` on the left; Fig. 9 step 6 has `(W_V + τα_L L)`. The second is
  dimensionally right and is used.
* **D5** Sec. 3.5.1 has `r = 1/(√|K| + r_ε⁻¹)`, Sec. 5.6 has `1/√(|K_i| + r_ε⁻¹)`. The first is
  dimensionally consistent and is used. `c = 0.02` is applied on the unit-area surface.
* **D6** Sec. 5.6 writes the adjoint source as `−2 W_V (L_a u + [0; y])`; `L_a u` is already
  integrated, so there is no `W_V`. The paper's `v` (Eq. 26) is minus the outward derivative.
* **D7** Sec. 6.2: "aspect ratio less than ρ₀ := 5" is read as greater than, with the aspect ratio of
  §3.7; "à la Sec. 5.6" is read as Sec. 5.3.1. The number of smoothing passes is not given.
* **D8** Sec. 7.1 states that the shape derivative of `∫_γ η ds` "is just ηκ"; by the paper's Eq. 8
  it is `∂η/∂n + ηκ`. The default step uses the full derivative (`L_w`); the paper's rule is
  selectable.
* **D9** Not stated in the paper: (a) curvature at boundary vertices; (b) a zero of `φ` at a vertex
  and crossings close to a vertex; (c) a closed component without a cut; (d) the dual length where
  the cut meets `∂M`; (e) the boundary condition of the extension on `∂M`; (f) values of `τ`, `α_L`,
  `r_ε`, a convergence test, the diagonal of the quad in a cut face, the units of `c`.
* **D11** Sec. 5.3 builds the inhomogeneous Laplacian by multiplying the cotangent of edge `ij` by
  `(a_i + a_j)/2`. Here every triangle's stiffness matrix is multiplied by the mean of `a` over its
  three vertices: this is the exact finite element matrix for piecewise linear `a`, and it stays
  positive semidefinite on obtuse triangles, which the edge rule does not guarantee. The two agree
  when `a` is constant.
* **D10** Redistancing is fast marching as in the paper, but the local update is the Hopf–Lax
  minimum over the opposite edge (Bornemann, Rasch 2006) and obtuse triangles are not unfolded
  (Kimmel, Sethian 1998 unfold them); both are first-order approximations of the same distance.
  Ties are broken by vertex index so that two implementations agree.

* **D12** Sec. 6.2 evaluates the gradient density per cut vertex from `v_i/ℓ_i` ("division by the
  dual boundary length is needed to convert the integrated values of `v_i` to pointwise values")
  and repairs artifacts by 1D Laplacian smoothing of the speed at the vertices of sliver triangles.
  When two or more consecutive crossings lie close to one mesh vertex this fails: `ℓ_i` of the cut
  vertices in between is as small as their distance, the error of `v_i` is not, and their
  neighbours are sliver vertices too, so the smoothing averages wrong values. Measured on cuts with
  a known solution (`examples/estimator_study.rs`; sphere caps, and a flat annulus with a source
  that varies along the cut): with `δ = 10⁻³`, wherever crossings crowd around a vertex, the
  paper's procedure has errors of `∂u/∂n` up to `6·10⁴ %` and of `s` up to `5·10⁷ %` at every
  refinement level, with one or three smoothing passes as without. §3.5 instead smooths the
  integrated fluxes along the cut, implicitly and with the mass `ℓ`, before dividing:
  `(diag(ℓ) + S) d = v`, and likewise for `q`. With `c = 0.5` the largest error of `∂u/∂n` is
  `0.2–5 %` and of `s` `0.1–8 %` on meshes of 8–40 thousand vertices (`1–22 %` and `3–33 %` on
  meshes of about 600 vertices), also when `φ = 0` exactly at vertices, and it falls under
  refinement; in the finite-difference test of §6 the error falls with the square of the mesh
  width, where the paper's procedure stalls near `10⁻²`. Alternatives that were measured and
  rejected (largest error of `∂u/∂n` over the same cases): replacing sliver values by
  interpolation between non-sliver neighbours (`30–120 %` where most cut vertices are slivers, no
  better than the paper's where all are), the gradient of `u` in the triangle on each segment
  (bounded, `20–130 %`), pooling runs of crossings within a fraction `θ` of one vertex
  (`θ = 0.25`: `1–90 %`), and couplings `∝ |σ|⁰` (not robust) or `|σ|⁻²` (as good, not simpler).
  The price of smoothing is a bias where `∂u/∂n` varies along the cut over a few edges: with 5
  edges per period the largest error is `15 %` (nodal, where it is usable: `8 %`). `c` is an
  option; sliver smoothing stays available and is off by default. The paper's procedure is
  `FluxRecovery::Nodal` with one pass.
* **D13** Not stated in the paper, and not in D9: (a) how the angles of the sub-triangles are
  computed: here from the barycentric inner product and the determinant area (§3.3), never from
  the three edge lengths, so a thin sub-triangle cannot violate the triangle inequality by
  rounding; (b) a closed component without a cut is left out of the area term of the energy too,
  and `u = p = 0` there is obtained by identity rows in the one system matrix; (c) the quad
  diagonal can tie on symmetric meshes: the rule is `p–c` iff `|pc|² ≤ |bq|²` in floating point,
  and near-ties are counted; (d) a region with `w_r = 0` is solved and reported like any other;
  (e) `Ω_f` is evaluated in the cancellation-free form of §3.1; (f) `min_crossing_offset = 0` is
  allowed and makes a crossing on a vertex an error instead of moving it.

* **D14** Fig. 9 step 6 integrates the length term on the mesh, `(W_V + τ α_L L) φ′ = W_V φ̃`, after
  redistancing that keeps the cut-adjacent values (Sec. 6.4). The curve-shortening law of §6
  item 4 then fails: the radius is 4–20 % too small, the error does not go to zero with `τ` and
  falls slowly with the mesh width (table in §6). Cause: the discrete Laplacian differentiates `φ`
  twice across the boundary between the cut-adjacent vertices, which keep their values, and their
  neighbours, which carry first-order fast-marching distances. The two differ by `O(κ h²)`, and
  `1/h²` times that is an `O(κ)` error of the curvature speed, with a sign. Only `τ ≫ h²` averages
  it out. Replacing the far field by exact distances to the cut polyline leaves `−9 %` (grid 32);
  redistancing every vertex from the polyline drifts with the number of steps
  (`examples/redistance_study.rs`). On the sphere the same step stops short of the balanced
  latitude by 6 % and shifts the stability threshold of the equator from 4.73 to about 4.0 (§6
  item 5). §3.10 therefore takes the length terms where the energy is defined, on the cut: the
  turning angle of the cut polyline is the exact first variation of its length. Measured and
  rejected on the way: the quotient of the exact gradients of length and swept area per crossing
  (ill-conditioned where a mesh edge is tangent to the cut: blow-up); `2 sin(θ/2)` instead of `θ`
  (not additive over crossings crowded at a vertex); no recovery of `Θ/ℓ` (kinks below the mesh
  width give `λ ∼ ε/ℓ²` and lock); the implicit system for `λ` alone, `r = s + λ`, which ignores
  that `s` moves the cut too (conditionally stable: `τ ω π/h ≲ 1`); taking the length force on
  the cut of `φ̃` after the explicit step (stable, but its equilibria depend on `τ`). The paper's
  step is `LengthStep::OnMesh`.
* **D15** Sec. 5.3.1 extends the speed harmonically through the cut mesh. Its triangles are
  arbitrarily thin, and through them the value at a cut-adjacent vertex is governed by cotangent
  weights that have nothing to do with which cut points are near it: measured, a vertex whose own
  crossings have the rates `+4`, `+16`, `+6` receives `−4` through a quad diagonal to a crossing
  with `−8`. The two ends of a crossed edge are then pushed apart, `|∇φ|` at the cut grows to
  10–100, and the cut zigzags; this happened with the length term on the cut, whose rate varies
  along the cut more than `s` does (with the paper's step on the mesh the implicit solve hides
  it). Clamping negative cotangents does not help (the weights stay unrelated to distance). An
  inverse-distance mean over the crossings of a vertex is stable but not conservative: crossings
  with a tiny dual length dominate, and the circle of §6 ran 0.3–2.3 % fast. The weights of §3.8
  are the transpose of the interpolation back to the cut, which keeps the mean rate. The paper's
  extension is `Extension::Harmonic`.
* **D16** Eq. 21 sets `φ_t` equal to the normal speed, which needs `|∇φ| = 1` at the cut. The
  redistancing keeps the cut-adjacent values, so `|∇φ|` in crossed faces is whatever the steps
  made it; vertices that just became cut-adjacent carry first-order distances, and the mean is
  0.97–0.99 on the meshes of §6 (the circle then runs about 2 % fast on grid 64 and 0.5 % on
  grid 128).
  §3.9 multiplies the rate by the mean of `|∇φ|` over the piece of the cut. Multiplying by the
  local `|∇φ|` was tried first and is unstable for large steps: a vertex whose value is too large
  gets a larger rate, which the cut leaves behind (`|∇φ|` grew from 1.05 to 9 in 15 steps on
  icosphere(5), `α_L = 20`, `τ = 8·10⁻⁵`). The paper's rule is `SpeedScaling::Unit`.
* **D17** Sec. 5.6 smooths the feature size with `(W_V + c L) r = W_V r⁰`. With the cotan Laplacian
  this is not an M-matrix on obtuse triangles, and `r` can leave the range of `r⁰`, even become
  negative, so that `a = r³` is not a weight. §3.6 uses the Laplacian of the edge weights
  `max(w_ij, 0)`; `tests/identities.rs` checks `min r⁰ ≤ r ≤ max r⁰` on a sheared mesh for
  `c` from 0.02 to 100.
* **D18** Sec. 3.2.1 says that distortion alone is minimized by longer and longer cuts and adds the
  length term. The linearization of §1 makes this quantitative, and has consequences the paper
  does not state: (a) with `α_L = 0` the flow has no well-posed limit, so the test "a latitude
  circle follows its ODE without length term" exists only on coarse meshes over short times, and
  the tests of §6 use `α_L > 0` wherever a time of order one is integrated; (b) `α_L` must exceed
  `ω h/π` times a small factor or the flow amplifies the noise of the mesh; (c) the scale of
  `α_L` and of `τ` is `ω`, which differs by orders of magnitude between the energies and, for the
  scaled Dirichlet energy, with `r_ε³`; hence `StepSize::Adaptive` and
  `length_weight_for_ripples` (§4); (d) the equator of the sphere of Fig. 1 is stable iff
  `α_L > 4/(3R)`, the tennis-ball seam is the mode `m = 2`. The rates were measured (§6 item 5b).
* **D19** Not stated in the paper, and not in D9: (a) `Capped` and `Adaptive` limit the motion by
  `s + λ⁰`, distortion and length together, when the length step is on the cut; (b) a run also
  stops when the energy stands still (§4); (c) the helpers for initial cuts (§5): distances of
  the geodesic circle are fast-marching distances and a vertex at the radius is inside; vertex
  labels put every crossing on an edge whose ends differ, at the midpoint on congruent
  triangles; face labels go to vertices by the larger angle, a tie to label 0; a function given
  by the caller is divided by the mean `|∇φ|` over its cut before it is redistanced, which moves
  no crossing; (d) the visibility penalty is `ν_i · max(n_i · e_i, 0)` with the area-weighted
  vertex normal and no occlusion test; (e) the order of the polylines of §7; (f) without classes,
  symmetry is kept to rounding only if no diagonal ties (§3.11).

Add to this list whenever something else turns up; never drop an entry.

## 10. Formula ids

Rust functions carry `[formula: <id>]`; the Lean theorem or the `CITED.md` heading carries the same
id. Discrete, to be proved in `lean/proofs`: `cut-face-curvature`, `cut-face-curvature-sum`,
`cut-edge-crossing`, `cut-barycentric-length`, `cut-subtriangle-area`, `cut-yamabe-system`,
`cut-normal-derivative`, `cut-flux-recovery`, `cut-dirichlet-energy`, `cut-weighted-laplacian`, `cut-gradient-dirichlet`, `cut-adjoint-scaled`,
`cut-gradient-scaled`, `cut-adjoint-pointwise`, `cut-gradient-pointwise`, `cut-feature-size`,
`cut-speed`, `cut-sliver-smoothing`, `cut-harmonic-extension`, `cut-explicit-step`,
`cut-length-step`, `cut-penalty-step-paper`, `cut-symmetry-projection`, `cut-redistance`,
`cut-step-cap`, and for the default stages `cut-extension-local` (convex combination, mean rate
kept), `cut-gradient-norm`, `cut-length-step-on-cut` (both matrices are M-matrices with row sums
`ℓ`; totals kept; `r` within the range of `s + λ⁰`; rest state independent of `τ`),
`cut-roughening-rate`. Analytic instances, provable with Mathlib calculus: `cut-disk-radial-solution`,
`cut-sphere-cap-solution`, `cut-curve-shortening-circle`, and the growth rate of a ripple of the
equator of §1 (`cut-roughening-rate`). Smooth theory, cited:
`cut-yamabe-equation` (Eq. 3), `cut-shape-derivative` (Eq. 6–8), `cut-cea-method` (Eq. 9–12,
Table 1), `cut-level-set-evolution` (Eq. 20–22).
