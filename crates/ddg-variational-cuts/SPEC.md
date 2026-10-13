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

The first term is extended harmonically away from `γ` and integrated explicitly; the second is
integrated implicitly; then `φ` is made a signed distance again. A region whose own boundary
derivative is large **shrinks** (D2). For `N` regions (§8) every region `k` has its own `φ_k`,
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
A curvature override replaces `k_f` by a caller-supplied density (`Uniform(c)` or one value per face).

### 3.2 Regions and crossings `[formula: cut-edge-crossing]`

`region(i) = 0` if `φ_i ≤ 0`, else `1`. An edge `[a, b]` with `region(a) ≠ region(b)` is *crossed*;
its cut vertex sits at parameter `t` from `a` (point `(1−t)·a + t·b`):

```
t = φ_a / (φ_a − φ_b),   then moved into [δ, 1−δ]
```

`δ` is the option `min_crossing_offset` (default `1e-3`); the number of crossings that were moved is
reported (D9b). Vertices of crossed edges are the *cut-adjacent* vertices `Adj`. A face has zero or
two crossed edges; in the second case it carries one *cut segment* joining its two cut vertices.
The cut `γ` is the union of the segments: closed polylines, or polylines ending on `∂M`.

### 3.3 Cut mesh

Every crossed face `abc` with `a` alone on its side, cut vertices `p ∈ ab`, `q ∈ ac`, is replaced by
the triangle `(a, p, q)` and the quad `(p, b, c, q)`; the quad is split by the shorter of its
diagonals `p–c`, `b–q` (tie: `p–c`). Orientation is inherited. Sub-triangles get the region of their
original vertices, the density `k_f` and (for scaled Dirichlet) linearly interpolated `a` of their
parent.

Lengths inside a face come from barycentric coordinates `[formula: cut-barycentric-length]`: for two
points with barycentric difference `d = (d_a, d_b, d_c)`, `Σ d = 0`,

```
|x − y|² = −( d_a d_b ℓ_ab² + d_b d_c ℓ_bc² + d_c d_a ℓ_ca² )
```

and areas from `[formula: cut-subtriangle-area]`: `area = |det(barycentric coordinates)| · A_f`.

The *submesh* `S_r` of region `r` consists of its (sub-)triangles. Its vertices are the original
vertices of the region and the cut vertices (each cut vertex belongs to both submeshes). `B_S` =
cut vertices ∪ (`B_M` ∩ `S`); `I_S` = the rest. A submesh may have several connected components
and any topology.

### 3.4 Scale factor `[formula: cut-yamabe-system]`

With `L` the cotan Laplacian of `S` and `g_i = Σ_{t∋i, t⊂S} k_t A_t / 3` (the paper's `W_F K`):

```
u_B = 0,        L_II u_I = −g_I
```

A component of `S` with no boundary vertex has no solution; this happens only for a closed
component of `M` that the cut does not touch. It is reported and left out of everything (D9c).

Remark for implementers: `I_S` over both regions is exactly the interior vertices of `M`, and the
pattern of the block-diagonal matrix is a subset of the pattern of `M`'s Laplacian, so one symbolic
factorization serves the whole flow, and one numeric factorization per step serves §3.4, §3.6, §3.8.

### 3.5 Integrated normal derivative `[formula: cut-normal-derivative]`

For `i ∈ B_S`, `v_i := ∮ ∂u/∂n · (hat function of i) ds` with `n` **outward** from `S`:

```
v_i = (L_BI u_I)_i + g_i
ℓ_i = ½ Σ (lengths of the boundary edges of S at i)
(∂u/∂n)_i ≈ v_i / ℓ_i
```

(The paper's `v` in Eq. 26, 28 has the opposite sign, D6.) At a cut vertex inside `M` both regions
have the same `ℓ`; at a cut vertex on `∂M` the two boundary edges are a cut segment and a piece of
the boundary edge of `M` (D9d).

### 3.6 Energy, adjoint, gradient density

```
Dirichlet   J = uᵀ L u = −u_I · g_I                                   [formula: cut-dirichlet-energy]
            G_i = (v_i/ℓ_i)²                                           [formula: cut-gradient-dirichlet]

scaled      J = uᵀ L_a u,   L_a = Σ_t ā_t · (cotan Laplacian of the single triangle t),
                            ā_t = (a_i + a_j + a_k)/3            [formula: cut-weighted-laplacian]
            L_II p_I = r_I,   r_i = −2 (L_a u)_i            (i ∈ I)   [formula: cut-adjoint-scaled]
            q_i = (L_BI p_I)_i − r_i,   r_i = −2 [ (L_a u)_i − a_i v_i ]   (i ∈ B)
            G_i = −a_i (v_i/ℓ_i)² − (v_i/ℓ_i)(q_i/ℓ_i)                [formula: cut-gradient-scaled]

pointwise   J = Σ_{i∈S} A_i^S b(u_i),   A_i^S = lumped vertex area within S
            L_II p_I = r_I,   r_i = −A_i^S b′(u_i)                     [formula: cut-adjoint-pointwise]
            q_i = (L_BI p_I)_i − r_i     (i ∈ B, where u_i = 0)
            G_i = b(0) − (v_i/ℓ_i)(q_i/ℓ_i)                           [formula: cut-gradient-pointwise]
```

`L_a` is the exact stiffness matrix of `∫ a ∇f·∇g` for piecewise linear `a, f, g`, hence positive
semidefinite for `a ≥ 0` on every mesh (D11). `q` is the integrated outward normal derivative of
`p`, defined like `v` (the extra `W_V` of the paper's Sec. 5.6 is not there, D6). With `a ≡ 1`: `p = −2u`, `q = −2v`, and the scaled density is
the Dirichlet one.

Curvature-based weight `[formula: cut-feature-size]` (paper Sec. 3.5.1, 5.6), computed once on `M`:

```
K_i  = (2π − Θ_i)/A_i  (interior),  0 on B_M
r⁰_i = 1 / ( √|K_i| + 1/r_ε )                (D5)
(W_V + c L) r = W_V r⁰,   c = 0.02
a_i  = r_i³
```

Region totals: `G_r,i = w_r · G_i + β_r`. A region with `w_r = 0` needs no solve.

### 3.7 Speed on the cut `[formula: cut-speed]`

At every cut vertex `c`: `s_c = G_0,c − G_1,c`.

Sliver smoothing `[formula: cut-sliver-smoothing]` (paper Sec. 6.2, D7): a cut vertex is a *sliver
vertex* if an incident cut-mesh triangle has aspect ratio `> ρ₀` (default 5), where aspect ratio is
`(longest edge)² · √3 / (4 · area)` (1 for an equilateral triangle). For `passes` rounds (default 1),
every sliver vertex simultaneously takes the arclength-weighted mean of its two neighbours along
the cut, `s_c ← (ℓ₂ s₁ + ℓ₁ s₂)/(ℓ₁ + ℓ₂)`, with `ℓ_k` the length of the segment to neighbour `k`.
A cut vertex on `∂M` (one neighbour) is left alone.

### 3.8 Harmonic extension `[formula: cut-harmonic-extension]`

In each submesh, with `s` as Dirichlet data on the cut vertices `C` and all original vertices `O`
of the submesh unknown (so `∂M` gets the natural, zero-flux condition, D9e):

```
L_OO σ_O = −L_OC s_C
```

A component without cut vertices gets `σ = 0`. On a closed mesh `L_OO` is the matrix of §3.4.
Every original vertex of `M` now has a value `σ_i`.

### 3.9 Explicit step `[formula: cut-explicit-step]`

`φ̃_i = φ_i + τ σ_i`.

### 3.10 Length and penalty step `[formula: cut-length-step]`

On `M`, with `w_i = α_L + η_i ≥ 0` and `L_w = Σ_f w̄_f · (cotan Laplacian of the single face f)`,
`w̄_f` the mean of `w` over the three vertices of `f` (the weighted Laplacian of §3.6 on `M`):

```
(W_V + τ L_w) φ′ = W_V φ̃
```

This is backward Euler for the second term of (★); the matrix is symmetric positive definite.
With `η = 0` it is the paper's step 6 of Fig. 9 (D4). The paper's own penalty rule (Sec. 7.1),
selectable as `PenaltyRule::PaperDiagonal` `[formula: cut-penalty-step-paper]`, is

```
(W_V D⁻¹ + τ L) φ′ = W_V D⁻¹ φ̃,   D = diag(α_L + η_i) > 0
```

which drops the `∂η/∂n` term of the gradient (D8). The default is the gradient rule above.
`w ≡ 0` skips the step.

### 3.11 Symmetry `[formula: cut-symmetry-projection]`

If vertex classes are given (paper Sec. 7), every vertex takes the arithmetic mean of `φ′` over its
class, after §3.10 and again after §3.12.

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
1  regions, crossings, cut mesh                       §3.2 §3.3
2  u, v (and p, q) in each region                     §3.4–3.6
3  s on the cut, sliver smoothing                     §3.7
4  σ on V by harmonic extension                       §3.8
5  φ̃ = φ + τσ                                        §3.9
6  (W_V + τ L_w) φ′ = W_V φ̃                          §3.10, then §3.11
7  redistance                                         §3.12, then §3.11
```

**Step size.** `Fixed(τ)` is the paper's. `Capped { tau_max, max_motion }` uses
`τ = min(tau_max, max_motion · h / max_c |s_c|)` with `h` the shortest crossed edge
`[formula: cut-step-cap]` (the paper's remark in Sec. 7.4). The paper gives no values (D9f).

**Stopping.** A run stops at `max_steps`, or when the cut motion of a step (largest `|φ_new − φ_old|`
over vertices cut-adjacent before or after the step) falls below a tolerance, or when the cut has
vanished (one region is empty). The reason is in the report. Every step reports the energy terms of
§6, `τ`, the cut motion and the counters named above.

Between steps the caller may change `α_L`, `η`, the energy, `w_r`, `β_r` and the step size (the
schedules of Sec. 7.4). The state is `φ`; a flow can be resumed from a stored `φ`.

## 5. Steering: what the caller controls

| Paper | Control |
|---|---|
| Eq. 5 | `α_L` length weight, may change between steps |
| Sec. 3.2, 3.5 | energy: Dirichlet, scaled Dirichlet (`a` curvature-based or supplied), Hencky, pointwise `b, b′` |
| Sec. 7.1 | `η` per vertex (the library takes numbers; a helper computes the visibility penalty `η_i = ν_i · n_i·E_i` of the paper from positions and a camera point) |
| Sec. 7, Fig. 2, 3, 12 | symmetry classes |
| Fig. 14 | per-region distortion weight `w_r` and area weight `β_r` |
| Sec. 7.4, Fig. 17, 18 | curvature override: uniform or per face |
| Sec. 6, 7.4 | step size rule, stop rule |
| Fig. 9 input | initial `φ⁰`: any values; helpers for a geodesic circle, vertex labels, face labels, normal clusters |
| (workspace) | fixed seams: cut the mesh along them first, they are then part of `∂M` |

## 6. Discrete energy and verification targets

```
E(φ) = Σ_r [ w_r J_r + β_r Area(S_r) ] + Σ_segments ½ (w_p + w_q) |pq|
```

(`w` interpolated linearly to the cut vertices). The flow is a discretization of the gradient flow
of the smooth energy, not the exact gradient flow of `E(φ)`; tests therefore use tolerances that
shrink under refinement.

1. `Σ_f Ω_f` identity (§3.1) to rounding; areas of the cut mesh sum to 1; Euler characteristics.
2. Unit sphere scaled to unit area (`R = 1/√(4π)`), cut along the latitude circle of polar angle
   `ψ₀`: in the cap, `u = 2 ln( cos(ψ₀/2) / cos(ψ/2) )`, `∂u/∂n = tan(ψ₀/2)/R`,
   `J = −8π ln cos(ψ₀/2) − 4π sin²(ψ₀/2)`. Convergence under icosphere refinement.
3. Finite differences: for a smooth `ψ`, `E(φ+εψ) − E(φ−εψ) ≈ −2ε ∫_γ [s + (length terms)] ψ ds`.
   This is the test of every sign in §1.
4. Flat disk, circular cut of radius `ρ₀`: `u = 0`, `s = 0`, and the circle shrinks as
   `ρ(t) = √(ρ₀² − 2 α_L t)`.
5. The energy decreases along the flow for small `τ`; the equator of a sphere is stationary;
   symmetric input stays symmetric; results do not depend on vertex numbering or rigid motion and
   scale as stated in §2.
6. One step agrees with the Lean oracle on small meshes.
7. (Milestone 2) `N = 2` reproduces the two-region code.

## 7. Output

* the state `φ` (caller's unit) and the region of every vertex;
* the cut as points (`edge, t`) and segments (`point, point, face`), with a helper that returns 3D
  positions when the caller has them;
* the cut mesh as a `TriMesh` + `Geometry` with maps to parent faces, original vertices or cut
  points, and the region of every face. This is what a flattening consumes.

## 8. N regions (milestone 2, outline; details are fixed when milestone 1 is done)

State `φ_k(i)`, `k < N`; `region(i) = argmin_k φ_k(i)` (ties: lowest `k`). Crossing of an edge
`[a, b]` with regions `i ≠ j` (paper Eq. 31):
`t = (φ_a^i − φ_a^j) / (φ_a^i − φ_a^j + φ_b^j − φ_b^i)`. A face with three regions gets a triple point
at the centroid of its three cut vertices. Speeds per region boundary as in §1. Inside region `k`,
`φ̇_k` is the harmonic extension of its boundary speed and `φ̇_j = −φ̇_k` for every `j ≠ k`
(App. A.3). §3.10 and §3.12 act on every `φ_k`; cut-adjacent vertices are projected as in
Losasso et al. 2006, Sec. 3.1. With `N = 2`, `φ_0 = −φ_1 = φ` must reproduce §3–4.

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

Add to this list whenever something else turns up; never drop an entry.

## 10. Formula ids

Rust functions carry `[formula: <id>]`; the Lean theorem or the `CITED.md` heading carries the same
id. Discrete, to be proved in `lean/proofs`: `cut-face-curvature`, `cut-face-curvature-sum`,
`cut-edge-crossing`, `cut-barycentric-length`, `cut-subtriangle-area`, `cut-yamabe-system`,
`cut-normal-derivative`, `cut-dirichlet-energy`, `cut-weighted-laplacian`, `cut-gradient-dirichlet`, `cut-adjoint-scaled`,
`cut-gradient-scaled`, `cut-adjoint-pointwise`, `cut-gradient-pointwise`, `cut-feature-size`,
`cut-speed`, `cut-sliver-smoothing`, `cut-harmonic-extension`, `cut-explicit-step`,
`cut-length-step`, `cut-penalty-step-paper`, `cut-symmetry-projection`, `cut-redistance`,
`cut-step-cap`. Analytic instances, provable with Mathlib calculus: `cut-disk-radial-solution`,
`cut-sphere-cap-solution`, `cut-curve-shortening-circle`. Smooth theory, cited:
`cut-yamabe-equation` (Eq. 3), `cut-shape-derivative` (Eq. 6–8), `cut-cea-method` (Eq. 9–12,
Table 1), `cut-level-set-evolution` (Eq. 20–22).
