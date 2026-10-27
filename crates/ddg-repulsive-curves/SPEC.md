# Repulsive Curves — contract

The contract between the Rust crate (`crates/ddg-repulsive-curves`), the Lean oracle
(`lean/oracle`), the Lean proofs (`lean/proofs`), the Python binding (`crates/ddg-py`) and the
Blender add-on for

> Chris Yu, Henrik Schumacher, Keenan Crane. *Repulsive Curves.* ACM Transactions on Graphics
> 40(2), Article 10, 2021. <https://doi.org/10.1145/3439429>
> (`../RepulsiveCurves/RepulsiveCurves.pdf`, "the paper"; equation numbers are the paper's;
> `RepulsiveCurves_Supplemental.pdf`, "Supp.")

**Reference.** `../RepulsiveCurves/repulsive-curves` (MIT, Copyright (c) 2019 Christopher Yu;
"the reference", cited as `ref: <file>`). It may be read and ported; ported code keeps the notice
in `LICENSES/repulsive-curves-MIT.txt`. Its submodules (`deps/geometry-central`, `deps/polyscope`,
`deps/libgmultigrid`) were empty when this file was written and were not read; before code from
one of them is ported, its licence goes into `LICENSES/` too. The reference differs from the
paper in several places, some of them bugs (§11). This file says in each case what we do and why.

Every formula below was re-derived. Claims marked *(checked)* were verified with a throw-away
NumPy script (identities to 1e-15, finite differences to 1e-8, eigenvalues; not in the
repository). If an implementation or a test contradicts this file, stop and report it; do not
quietly follow either side.

Conventions of `CLAUDE.md` apply: `f64` everywhere, stored indices `u32`, bad input gives an
error enum naming the element and never a panic, nothing is normalized, clamped or repaired
silently, every solve returns a report, the thread count is an option (explicit pool, never the
global `rayon` pool), formula tags `[formula: rc-…]` (§12).

## 1. The problem

A curve network is `γ : M → R³` where `M` is a union of intervals and circles, possibly glued at
junctions. The generalized tangent-point energy (Sec. 3.2, Eq. 3) is

```
ℰ^α_β(γ) = ∬_{M²} k^α_β(γ(x), γ(y), T(x)) dx_γ dy_γ,     k^α_β(p, q, T) = |T × (p − q)|^α / |p − q|^β
```

with `T` the unit tangent. It is finite on embedded `C¹` curves exactly when `γ ∈ W^{s,α}` for
`α > 1`, `β ∈ [α + 2, 2α + 1)`, `s = (β − 1)/α ∈ (1, 2)` (Blatt 2013; paper Lemma A.1). Its
differential has order `2s` (A.3). `σ := s − 1 ∈ (0, 1)`. For `β − α > 2` the energy is not
scale invariant; it decreases when the curve grows (§3.3), so a length constraint or penalty is
normally needed.

We minimize `F = Ê + Σ penalties` over the vertex positions subject to constraints `Φ(γ) = 0`,
by projected gradient descent in the discrete fractional Sobolev–Slobodeckij metric
`A = B + B⁰` (Sec. 4–5), with the line searches of Sec. 5.4, optionally accelerated (Sec. 6,
App. B). For the evaluation of Sec. 7 the same machinery also runs L-BFGS and nonlinear
conjugate gradients on a soft-constrained problem, and projected gradients in the L², H¹ and H²
metrics (Supp. 2–3).

Out of scope (user decision): multi-agent path planning (Sec. 8.5) and streamline visualization
(Sec. 8.6). The field-alignment penalty stays in the API.

## 2. Discrete setting and data model

### 2.1 Curve network

`G = (V, E)`, positions `γ_v ∈ R³`, edges `I = (i₁, i₂)` with the orientation the caller gives.
Nothing depends on the orientation: `k` is even in `T`, and `D_I` (§4.1) changes sign together
on both sides of every product. Allowed: several components, closed loops, open curves
(endpoints of degree 1), junctions (degree ≥ 3).

```
ℓ_I = |γ_{i₂} − γ_{i₁}|,   T_I = (γ_{i₂} − γ_{i₁}) / ℓ_I,   x_I = (γ_{i₁} + γ_{i₂}) / 2        [formula: rc-edge-quantities]
m_v = ½ Σ_{I ∋ v} ℓ_I                    (dual length; the reference divides by deg v, D23)   [formula: rc-dual-length]
𝒫 = {(I, J) ∈ E × E : I ∩ J = ∅}        ordered pairs of edges without a common vertex
```

`ℓ̄₀` is the mean edge length when a session starts. It is the length unit of the tolerances
that are relative (§10).

### 2.2 Dimension (D18)

Positions are `[f64; 3]`; the crate is not generic in the dimension. Planar networks are handled
by `Planar { axis }`: one coordinate (X, Y or Z) is held at its input value. In planar mode the
input must have exactly the same value of that coordinate at every vertex
(`RcError::NotPlanar` otherwise), the metric drops that coordinate block (`Ā` acts on two
coordinates per vertex instead of three, so the solve is smaller), and gradients have no
component along the axis. The tangent-point kernel of a planar curve in R³ is exactly the planar
kernel (`|T × d|` is the 2-D cross product), so nothing else changes. Why not generic: the paper
and the reference work in R³; the reference's 2-D examples (`scenes/BluenoiseCurves`,
`scenes/PlanarGraphDrawing`) are R³ curves on the plane `y = 0`; Blender data is R³; and the
BVH, the block cluster tree, the CCD and the obstacles would all need two implementations for
no visible benefit. (Path planning, the only n-dimensional use in the paper, is out of scope.)

### 2.3 Input validation

`CurveNetwork::new(positions, edges)` rejects, with the element named in the message:

* `VertexOutOfRange { edge, vertex, n_vertices }`; `SelfLoop { edge, vertex }`;
  `DuplicateEdge { a, b, edges }` (two edges on the same unordered vertex pair);
* `IsolatedVertex { vertex }` (degree 0: no tangent, no mass; Blender strips loose vertices only
  when asked, see §9.3);
* `NonFinitePosition { vertex }`; `ZeroLengthEdge { edge, a, b }`;
* `CoincidentVertices { a, b }` (distinct vertices at the same point: `k` is infinite);
* `EdgesIntersect { a, b, distance }`: two edges in `𝒫` at distance `≤ 1e-12 ℓ̄₀`. The discrete
  energy is finite on intersecting edges whose vertices are distinct (the barrier is only
  approximate, paper Sec. 9), so this is checked explicitly (BVH broad phase, exact segment
  distance);
* `NoInteractingPairs`: `𝒫` is empty (a single edge, a triangle, a star): the energy and the
  metric are zero.

`Exponents::new(α, β)` rejects `α ≤ 1`, `β < α + 2`, `β ≥ 2α + 1`, non-finite values
(`InvalidExponents { alpha, beta, requirement }`; D24). Option structs implement `Default` with
the values of §10 and `validate()`.

## 3. The discrete energy

### 3.1 Kernel and energy (Eq. 3, 17, 18)

```
k(p, q, T) = |T × d|^α / |d|^β,     d = p − q,  |T × d| = |P_T d|,  P_T = I − T Tᵀ  (|T| = 1)   [formula: rc-kernel]
(k̂)_IJ  = ¼ Σ_{i ∈ {i₁,i₂}} Σ_{j ∈ {j₁,j₂}} k(γ_i, γ_j, T_I)                                    [formula: rc-discrete-kernel]
Ê(γ)    = Σ_{(I,J) ∈ 𝒫} (k̂)_IJ ℓ_I ℓ_J                                                         [formula: rc-energy]
```

Both orders `(I, J)` and `(J, I)` are summed; `k̂` is not symmetric (it uses `T_I`). Pairs of
edges with a common vertex are omitted (paper Sec. 5.1). This is the paper's discretization; the
reference uses another one (D1).

### 3.2 Differential

With `n = P_T d`, `ν = |n|`, `n̂ = n/ν`, `ρ = |d|`:

```
∇_p k = α ν^{α−1} n̂ / ρ^β − β k d / ρ²,   ∇_q k = −∇_p k                                       [formula: rc-kernel-gradient]
∇_T k = −α ν^{α−1} (T·d) n̂ / ρ^β        (the part along T^⊥, which is all that unit-T variations see)
```

`ν^{α−1} n̂` is defined as 0 at `ν = 0`; this is the continuous extension because `α > 1`
(the reference instead returns 0 below `ν < 1e-10`, D24). Edge variations:
`δℓ_I = T_I·(δγ_{i₂} − δγ_{i₁})`, `δT_I = P_{T_I}(δγ_{i₂} − δγ_{i₁})/ℓ_I`,
`δx_I = (δγ_{i₁} + δγ_{i₂})/2`. For `e_IJ = k̂_IJ ℓ_I ℓ_J`, `s₁ = −1`, `s₂ = +1`,
`g_T = ¼ Σ_i Σ_j ∇_T k(γ_i, γ_j, T_I)`:

```
∂e_IJ/∂γ_{i_a} = s_a ℓ_J (k̂_IJ T_I + g_T) + ℓ_I ℓ_J ¼ Σ_{j ∈ J} ∇_p k(γ_{i_a}, γ_j, T_I)
∂e_IJ/∂γ_{j_b} = s_b ℓ_I k̂_IJ T_J       + ℓ_I ℓ_J ¼ Σ_{i ∈ I} ∇_q k(γ_i, γ_{j_b}, T_I)          [formula: rc-energy-differential]
```

*(checked: central differences, max relative error 8e-9 for (3,6) and (2,4.5), two components.)*

**Evaluation order.** Per edge `I`, the row `R_I = Σ_{J : (I,J) ∈ 𝒫} (e_IJ + e_JI)` and its
derivative with respect to the two endpoints of `I` only, `G_{I,1}, G_{I,2}`. Then
`∂Ê/∂γ_v = Σ_{I ∋ v} G_{I,a(v)}` summed in increasing edge index, and
`Ê = Σ_I Σ_J e_IJ`. Every (pair, endpoint) derivative is computed exactly once. Rows are
independent (parallel over edges); every reduction has a fixed order, so the result is bitwise
independent of the thread count (D25).

### 3.3 Invariances and scaling

`Ê` is invariant under translations, rotations, reflections and edge re-orientation, and
homogeneous: `Ê(cγ) = c^{2+α−β} Ê(γ)` for `c > 0` [formula: rc-energy-scaling]. Hence
[formula: rc-gradient-identities]

```
Σ_v ∂Ê/∂γ_v = 0,     Σ_v γ_v × ∂Ê/∂γ_v = 0,     Σ_v γ_v · ∂Ê/∂γ_v = (2 + α − β) Ê
```

*(checked: 4e-15, 1e-14, and the Euler identity to 1e-14 relative.)* The exponent `2 + α − β` is
−1 for (3, 6) and −½ for (2, 4.5): the energy falls when the curve grows.

## 4. Inner products

### 4.1 Operators on vertex functions

```
(D u)_I = (u_{i₂} − u_{i₁}) T_I / ℓ_I ∈ R³,  D ∈ R^{3|E|×|V|},  D_I = (1/ℓ_I)[−T_I  T_I]   (Sec. 5.2.1)   [formula: rc-derivative-operator]
(E u)_I = (u_{i₁} + u_{i₂}) / 2,                E ∈ R^{|E|×|V|}                         (B.3.4)
```

### 4.2 High-order term `B` (Eq. 21)

```
uᵀ B v = Σ_{(I,J) ∈ 𝒫} w_IJ ⟨D_I u[I] − D_J u[J], D_I v[I] − D_J v[J]⟩
w_IJ = ¼ ℓ_I ℓ_J Σ_{i ∈ I} Σ_{j ∈ J} |γ_i − γ_j|^{−(2σ+1)}       (trapezoid, paper)
w_IJ =   ℓ_I ℓ_J |x_I − x_J|^{−(2σ+1)}                           (midpoint, reference and Eq. 32; D2)
```

Increments, for every `(I, J) ∈ 𝒫` and `a, b ∈ {1, 2}` (paper, Sec. 5.2.2):
`B_{i_a i_b} += (−1)^{a+b} w_IJ/ℓ_I²`, `B_{i_a j_b} −= (−1)^{a+b} w_IJ ⟨T_I, T_J⟩/(ℓ_I ℓ_J)`,
`B_{j_a j_b} += (−1)^{a+b} w_IJ/ℓ_J²`, `B_{j_a i_b} −= (−1)^{a+b} w_IJ ⟨T_J, T_I⟩/(ℓ_J ℓ_I)`
[formula: rc-high-order-gram]. *(checked against the quadratic form: 4e-15.)*

### 4.3 Low-order term `B⁰` (Eq. 20)

```
uᵀ B⁰ v = Σ_{(I,J) ∈ 𝒫} w⁰_IJ (u_I − u_J)(v_I − v_J),     u_I = (E u)_I
w⁰_IJ = ¼ ℓ_I ℓ_J Σ_{i ∈ I} Σ_{j ∈ J} k²₄(γ_i, γ_j, T_I) / |γ_i − γ_j|^{2σ+1}
      = ¼ ℓ_I ℓ_J Σ_{i ∈ I} Σ_{j ∈ J} k²_{2σ+5}(γ_i, γ_j, T_I)                 (trapezoid, paper)
w⁰_IJ =   ℓ_I ℓ_J k²_{2σ+5}(x_I, x_J, T_I)                                     (midpoint, reference and Eq. 32)
```

Increments: `B⁰_{i_a i_b} += ¼ w⁰_IJ`, `B⁰_{i_a j_b} −= ¼ w⁰_IJ`, `B⁰_{j_a i_b} −= ¼ w⁰_IJ`,
`B⁰_{j_a j_b} += ¼ w⁰_IJ` [formula: rc-low-order-gram]. *(checked: 3e-16.)*

### 4.4 Kernel-matrix form (Eq. 31, 32) and the quadrature (D2)

For `(I, J) ∈ 𝒫` let `K^B_IJ = w_IJ + w_JI` and `K⁰_IJ = w⁰_IJ + w⁰_JI`; both are 0 on the
diagonal and on adjacent pairs. Then [formula: rc-kernel-matrix-form]

```
B  = Dᵀ ((diag(K^B 𝟙) − K^B) ⊗ I₃) D,      B⁰ = Eᵀ (diag(K⁰ 𝟙) − K⁰) E
```

because `uᵀ Mᵀ(diag(K𝟙) − K) M v = ½ Σ_I Σ_J K_IJ ⟨(Mu)_I − (Mu)_J, (Mv)_I − (Mv)_J⟩` for
symmetric `K`. With the midpoint weights, `K⁰_IJ = (k²_{2σ+5}(x_I, x_J, T_I) + k²_{2σ+5}(x_J,
x_I, T_J)) ℓ_I ℓ_J` is exactly Eq. 32, and `K^B_IJ = 2 ℓ_I ℓ_J |x_I − x_J|^{−(2σ+1)} =
(k⁰_{2σ+1}(x_I, x_J, ·) + k⁰_{2σ+1}(x_J, x_I, ·)) ℓ_I ℓ_J` is Eq. 32 with `k⁰_{2σ+1}`, as B.3.4
says. *(checked: Eq. 31–32 reproduce the midpoint increments of §4.2–4.3 to 6e-16.)*

**There is no contradiction between Eq. 32 and Eq. 13/20 in the kernel exponents:**
`k²₄(p, q, T) · |p − q|^{−(2σ+1)} = |T × d|² / |d|^{4 + 2σ + 1} = k²_{2σ+5}(p, q, T)`. What
differs is (a) the quadrature (Eq. 20–21 trapezoid over the four endpoint pairs, Eq. 31–32 one
point at the midpoints) and (b) the symmetrization (Eq. 31 sums each unordered pair once with
`K = w + wᵀ` and the factor ½ of the identity above; Eq. 20 sums ordered pairs). (b) is an
identity. (a) is a real difference: on a 9-vertex test curve the two `B⁰` differ by 67 %
(entrywise, relative to the largest entry), on fine curves by `O(ℓ²/d²)`. The reference uses the
midpoint rule in its dense matrix (`ref: src/sobo_slobo.cpp`, `IntegrateLocalMatrix` with one
point at 0.5) and in its block cluster tree. Decision **D2**: `GramQuadrature::{Midpoint,
Trapezoid}`, default `Midpoint`. The dense and the hierarchical back ends use the same `K`, so
they agree as the admissibility parameters go to 0 (§13.4). The default is confirmed by the
measurement of §14.3.

### 4.5 The metric `Ā` (Eq. 19) and its null space

`A = B + B⁰ ∈ R^{|V|×|V|}`, `Ā = A ⊗ I₃` acting on the interleaved vector
`(γ₁ₓ, γ₁ᵧ, γ₁_z, γ₂ₓ, …)` (D29; the paper displays the block form). Properties
[formula: rc-gram-psd-kernel]:

* `B`, `B⁰`, `A` are symmetric positive semidefinite (sums of `w ≥ 0` times squares).
* `u ∈ ker B ⇔ D_I u = D_J u` for every `(I, J)` with `w_IJ > 0`; `ker B` contains the functions
  that are constant on each component. `u ∈ ker B⁰ ⇔ u_I = u_J` for every `(I, J)` with
  `w⁰_IJ + w⁰_JI > 0`; `ker B⁰` contains the global constants and, since `u_I` is the edge
  average, the alternating ±1 colouring of every bipartite component (even loops, open curves,
  trees; WP2 finding). `ker A = ker B ∩ ker B⁰`.
* For a generic network `ker A = span{𝟙}`: `B⁰` couples all pairs of edges, also across
  components. *(checked: two odd loops, dim ker B = 2, dim ker B⁰ = dim ker A = 1; with even
  loops `ker B⁰` is larger, `ker A = span{𝟙}` still.)* So `Ā` has
  exactly the three translations in its kernel (Sec. 4.2.4); a constraint must remove them
  (§5.2, D5).
* Degenerate: if all edges are collinear, `k² ≡ 0`, `B⁰ = 0` and `ker A` is larger. The solver
  reports `SingularMetric` (the energy and its differential are 0 there anyway).

Scaling [formula: rc-gram-scaling]: `w(cγ) = c^{1−2σ} w(γ)`, `w⁰(cγ) = c^{−(2σ+1)} w⁰(γ)`, hence
`B(cγ) = c^{−(2σ+1)} B(γ)`, `B⁰(cγ) = c^{−(2σ+1)} B⁰(γ)` *(checked: exponents −7/3 and −5/2 for
(3,6) and (2,4.5))*. With `dÊ(cγ) = c^{1+α−β} dÊ(γ)`, the Sobolev gradient scales as
`c^{2σ+2+α−β}` (1/3 for (3,6), 1 for (2,4.5)): rescaling the input changes its length, not its
direction (A.4.2; the printed exponents there are wrong, D3). With constraints whose Jacobian
rows are homogeneous of degree 0 (all of §5.1 except `Surface` with a non-homogeneous `f`), the
projected direction of §6.1 scales the same way [formula: rc-gradient-scaling].

### 4.6 The comparison metrics L², H¹, H² (Supp. 2.1, D23)

```
L²:  M = diag(m_v)
H¹:  L_ij = −1/ℓ_ij for every edge (i, j),  L_ii = Σ_{I ∋ i} 1/ℓ_I            (weak Laplacian of the curve)
H²:  L M⁻¹ L                                                                   [formula: rc-classic-metrics]
```

Sparse; `ker L` and `ker LM⁻¹L` are the per-component constants, so with several components the
constraints must fix every component's translation (`SingularMetric { components }` otherwise).
The Sec. 7 knots have one component.

## 5. Constraints, penalties, length schedules

### 5.1 Constraints (Sec. 8.1)

`C = dΦ(γ)`, rows in the order the constraints are listed. Targets default to the values at
session start ("keep"). Tolerances per row type, used by the projection (§6.2): length rows
`τ_len`, dimensionless rows `τ_len / ℓ̄₀`, surface rows `τ_len · |∇f(γ_v)|`.

| constraint | `Φ` | rows | nonzero Jacobian blocks | row unit |
|---|---|---|---|---|
| `Barycenter { target }` | `(Σ_I ℓ_I x_I)/(Σ_I ℓ_I) − x₀` | 3 | `∂/∂γ_v = (1/L) Σ_{I ∋ v} [ (ℓ_I/2) I₃ + s_{I,v} (x_I − x̄) T_Iᵀ ]` | length |
| `TotalLength { target }` | `L⁰ − Σ_I ℓ_I` | 1 | `∂/∂γ_v = −Σ_{I ∋ v} s_{I,v} T_Iᵀ` | length |
| `EdgeLengths { edges, targets }` | `ℓ⁰_I − ℓ_I` | 1 per edge | `i₁: +T_Iᵀ`, `i₂: −T_Iᵀ` | length |
| `Point { vertex, target }` | `γ_v − x_v` | 3 | `I₃` | length |
| `Surface { vertices, f }` | `f(γ_v)` | 1 per vertex | `∇f(γ_v)ᵀ` | f |
| `Tangent { site, direction X }` | `(e₁·s, e₂·s)` | 2 | `±e_kᵀ P_{T_I}/ℓ_I` per edge in the site | 1 |

`s_{I,v} = +1` if `v = i₂`, `−1` if `v = i₁`; `L = Σ ℓ_I`; `x̄` the current barycenter
[formulas: rc-constraint-barycenter, rc-constraint-total-length, rc-constraint-edge-length,
rc-constraint-point, rc-constraint-surface, rc-constraint-tangent].

* **Barycenter (D6).** The paper's `Σ_I ℓ_I (x_I − x₀)` divided by `L`: same zero set, same
  tangent space at the constraint set, but the violation is in length units like the other rows.
  Exact Jacobian (the reference uses dual-length weights `Σ/deg` and treats them as constant).
* **Tangent (D7).** `site` is `Edge(I)` (`s = T_I`, the paper's `Φ_tangent,I`) or `Vertex(v)`
  for `deg v ∈ {1, 2}` (`s = T` of the edge at an endpoint, oriented away from `v`; at a degree-2
  vertex with neighbours `a` (first incident edge in edge order) and `b`:
  `s = (γ_v − γ_a)/|γ_v − γ_a| + (γ_b − γ_v)/|γ_b − γ_v|`, the reference's vertex tangent up to
  normalization). `e₁, e₂` complete `X/|X|` to an orthonormal frame, fixed at setup (deterministic
  construction, Duff et al. 2017). `(e₁·s, e₂·s) = 0 ⇔ s ∥ X` [formula: rc-tangent-frame]. The
  paper's `T_I − X` has three rows of rank 2 (`P_T T = 0`), which makes `C Cᵀ` singular. The zero
  set also contains `s` antiparallel to `X`: `s·X > 0` is required at setup
  (`InfeasibleConstraint`) and a projection that ends with `s·X ≤ 0` counts as failed.
  `Vertex(v)` with `deg v ≥ 3` is `InvalidTangentSite`.
* **Surface (D8).** `f` is a user object (`ImplicitSurface`: value and gradient), with built-in
  sphere, torus, plane, union (the reference's set, `ref: src/implicit_surface.cpp`). `f` need
  not be a signed distance; the row tolerance is scaled by `|∇f|` so that it means length.
* **Point.** Targets may differ from the current positions only within the projection tolerance
  (all constraints must hold at the start, as in the reference's `scenes/FORMATS.md`;
  `InfeasibleConstraint { constraint, violation }` otherwise).

### 5.2 Translation, rank and redundancy (D5, D9)

Checked when a session starts and whenever the constraint set or the topology changes:

1. **Translation.** `C [t_x t_y t_z]` (the three translation fields; two in planar mode) must
   have full rank, otherwise `Ā` is singular on `ker C`: `TranslationNotFixed`. The reference
   adds a barycenter constraint when a scene has none; here the default constraint list contains
   `Barycenter` explicitly (§10), and nothing is added silently.
2. **Row rank.** `C Cᵀ` is factored (sparse Cholesky, `ddg_core::sparse`, `Llt`). Known
   redundancies are resolved and reported in `dropped_rows`: an `EdgeLengths` row whose two
   endpoints are both `Point`-constrained (its value is then fixed; checked for consistency), and
   the `Barycenter` rows when every vertex of the network is `Point`-constrained. Any other rank
   deficiency is `RedundantConstraints { rows }` naming the constraints. The paper is silent.

### 5.3 Penalties (Sec. 8.1)

Added to the objective and to its differential; the metric is unchanged (paper: lower order).

| penalty | value | gradient |
|---|---|---|
| `TotalLength { w }` | `w Σ_I ℓ_I` | `i₂ += w T_I`, `i₁ −= w T_I` |
| `LengthDifference { w }` | `w Σ_{v : deg v = 2} (ℓ_{I_v} − ℓ_{J_v})²` | `2w(ℓ_I − ℓ_J)(∇ℓ_I − ∇ℓ_J)` |
| `FieldAlignment { X, w }` | `w Σ_I ℓ_I \|T_I × X(x_I)\|²` | below |
| `SurfacePotential { obstacle, w, p }` | below | below |
| `Hencky { w, rest }` | `w Σ_I ℓ⁰_I ln²(ℓ_I/ℓ⁰_I)` | `∂/∂ℓ_I = 2w (ℓ⁰_I/ℓ_I) ln(ℓ_I/ℓ⁰_I)` |

[formulas: rc-penalty-total-length, rc-penalty-length-difference, rc-penalty-field,
rc-penalty-surface, rc-penalty-hencky]

* Weights of `TotalLength` and `LengthDifference` may be negative (the reference's scene format
  allows it); the others must be ≥ 0 (`InvalidOption`).
* **Field alignment.** With `c = T_I × X_I`, `X_I = X(x_I)`, `J_X` the field's Jacobian:
  `∂(ℓ|c|²)/∂γ_{i_a} = s_a |c|² T_I + 2 s_a P_{T_I}(X_I × c) + ℓ_I J_X(x_I)ᵀ(c × T_I)`, times
  `w`. `X` is a user object (`VectorField`: value and Jacobian), built-in constant and circular
  fields. `X` is not normalized (the paper assumes a unit field).
* **Surface potential (D19).** The obstacle is a triangle soup (no manifold requirement), rigid.
  `Ê_M = w Σ_v m_v Σ_t A_t δ(γ_v, t)^{−p}` with `δ` the point–triangle distance (Ericson 2004,
  §5.1.5) and `p = β − α` by default (paper: "a Coulomb potential of the same order"). Near
  field exactly as written; far field lumped per node `N` of a 3-D BVH over the triangles:
  `A_N |γ_v − c_N|^{−p}`, `c_N` the area-weighted centroid, used when `γ_v` lies outside `N`'s
  box and `r_N ≤ θ_M |γ_v − c_N|`. Gradient: `δ` is `C¹` away from the triangle with
  `∇δ = (γ − π_t(γ))/δ` (`π_t` the closest point) [formula: rc-distance-gradient], plus the
  `∂m_v/∂γ` terms. The paper's `∬ |x_M − γ|^{−(β−α)}` sampled at points (the reference: mesh
  vertices with dual areas) has no barrier inside a face of a coarse obstacle; `δ` puts the
  barrier on the faces.
* **Hencky** (Supp. 2.2: "Hencky strain terms … for each edge, weighted by edge lengths of the
  initial configuration") is the soft replacement of `EdgeLengths` for the evaluation. The weight
  is not given; it is calibrated (§14.2, D21).

### 5.4 Target-length schedules and subdivision (Sec. 8.1, D13)

`LengthSchedule { factor, rate }` moves the target of `TotalLength`, or all targets of
`EdgeLengths` (`ScheduleWithoutLengthConstraint` if neither is present). Step `k` uses the total
target `L*_k = L⁰ + sign(factor − 1) · min(k · rate · L⁰, |factor − 1| L⁰)`; edge targets are
`ℓ⁰_I · L*_k / L⁰` (one common factor, ratios kept). Default `rate = 0.01` (the reference grows
the total by `L⁰/100` per step). The reference instead adds the same absolute amount to every
edge target, which changes the ratios [formula: rc-length-schedule].

`Subdivision { max_rounds }` (reference `subdivide_limit`; default 0 = off): when the mean edge
length exceeds twice its value at session start, every edge is split at its midpoint. Numbering:
old vertices keep their indices, the midpoint of edge `I` is vertex `|V| + I`, child `(i₁, m)`
keeps index `I`, child `(m, i₂)` is `|E| + I`. Constraints map as follows: `EdgeLengths` targets
are halved on both children; `TotalLength`, `Barycenter`, `Point` unchanged; `Tangent(Edge(I))`
applies to both children; `Tangent(Vertex(v))` unchanged; `Surface` rows are added for midpoints
of edges whose two endpoints are surface-constrained. All constraints still hold after the split
except `Surface` on the new vertices, which the next re-target projection (§7.1) fixes
[formula: rc-subdivision-constraints].

## 6. Descent direction and projection

### 6.1 Projected Sobolev gradient (Sec. 5.3.1, Eq. 22)

```
[ Ā  Cᵀ ] [ g̃ ]   [ dFᵀ ]
[ C  0  ] [ λ  ] = [ 0   ]                                         [formula: rc-saddle-gradient]
```

`dF` is the differential of `Ê` plus the penalties. If `Ā` is positive definite on `ker C`
(guaranteed by §4.5 and §5.2), the solution exists, is unique, and `g̃` minimizes
`½ ‖g̃ − g‖²_Ā` over `C g̃ = 0` (with `g = Ā⁻¹ dFᵀ` where `Ā` is invertible). Moreover
`dF · g̃ = g̃ᵀ Ā g̃ ≥ 0`, with equality only for `g̃ = 0`, so `−g̃` is a descent direction tangent
to the constraint set [formula: rc-descent-direction]. The report carries `‖g̃‖₂` and the
squared Sobolev norm `dF·g̃` (the reference's `soboDot`).

### 6.2 Constraint projection (Sec. 5.3.2, Eq. 23)

From `γ̃ = γ − τ g̃`, repeat

```
[ Ā  Cᵀ ] [ x ]   [ 0       ]
[ C  0  ] [ µ ] = [ −Φ(γ̃)  ],    γ̃ ← γ̃ + x                      [formula: rc-chord-projection]
```

with `Ā` and `C` frozen at `γ` (a chord method: the paper's "approximation of Newton's method",
which lets the dense back end reuse one factorization for Eq. 22 and every iteration, Sec. 7.4).
Each iteration gives the `Ā`-smallest `x` solving the linearized constraints (proved). Success
when every row is within its tolerance (§5.1) after at least one iteration (a trial already
within the tolerance is projected once too, unless `Φ = 0` exactly; §17l, found by WP12b).
**Departure from the paper (§17l, default):** once the largest relative violation is at most
`FlowOptions::newton_polish` (100) times its tolerance, the remaining iterations are Newton steps
with the Jacobian at the iterate and the lumped-mass metric (`x = M⁻¹Cᵀ(CM⁻¹Cᵀ)⁻¹(−Φ)`, a sparse
factorization each) instead of chord steps with `Ā`; `newton_polish: None` is the paper's chord
method throughout. **Also a departure (default):** a chord iteration whose observed rate cannot
reach its target within the remaining iterations is stopped early (`projection_early_abort`). An optional second-order correction of the trial (`γ − τp + τ²q`,
`C q = −½Φ''[p, p]`, off by default) is described in §17l. Failure after `max_projection_iterations` (3; the
paper: "rarely more than three", the reference: 3 dense, 2 multigrid), when the largest violation
grows by more than 10×, or when a tangent site ends with `s·X ≤ 0`. A failure rejects the trial
step (§7). Convergence of the chord method near a regular point: cited (D10).

### 6.3 Null-space projector (B.3.3, Eq. 30, 33; D4)

```
C† = Cᵀ (C Cᵀ)⁻¹,     P = I − C† C                                [formula: rc-projector]
```

for `C` of full row rank: `C P = 0`, `P² = P`, `Pᵀ = P`, `range P = ker C`. The paper prints
`P := C C†` (which is the `k×k` identity) and `C† := (C Cᵀ)⁻¹ Cᵀ` (whose factors do not even
have compatible shapes). `P` is applied with the sparse Cholesky factor of `C Cᵀ`, once per level
per step.

* Eq. 30: `y` solves `Pᵀ Ā P y = Pᵀ a` ⇔ `x = P y` solves Eq. 29 (with `Ā` positive definite on
  `ker C`; `x` is unique, `y` is not) [formula: rc-projected-system]. *(checked)*
* Eq. 33: with `z = C† b`, `b = −Φ(γ̃)`, and `y` from Eq. 29 with right side `Ā z`, `x = z − y`
  solves Eq. 23 [formula: rc-projection-transform]. *(checked)*

### 6.4 Soft-constraint methods (Supp. 3.2–3.3; evaluation only)

Objective `F = Ê + Hencky (+ other penalties)`, no hard constraints. Translations are removed
from the direction only (Supp.: "gradient … orthogonal (in the L² sense) to the space of
constant vector fields"): the preconditioned gradient `z` solves Eq. 22 with `C` replaced by the
three rows `M_c = [m₁I₃ … m_|V|I₃]` [formula: rc-mean-free-direction]. The metric is that of
the current curve at every iteration ("dynamic", Supp. 3.3).

* **L-BFGS** (Liu & Nocedal 1989), memory 30 (Supp. 3.3), two-loop recursion (Nocedal & Wright
  Alg. 7.4) with the dynamic initial inverse Hessian: the solve above in the middle of the two
  loops. The line search is Armijo only (§7), so a pair with `sᵀy ≤ 1e-12 |s||y|` is not stored
  (D22).
* **NCG**, Polak–Ribière with the same preconditioner, `β = max(0, ⟨r_k, z_k − z_{k−1}⟩ /
  ⟨r_{k−1}, z_{k−1}⟩)` (PR+), `d_k = −z_k + β d_{k−1}`, restart with `d_k = −z_k` when
  `dF·d_k ≥ 0` (D22).
* The L², H¹, H² projected gradients are §6.1 with `Ā` replaced by the metrics of §4.6.

## 7. Time stepping (Sec. 5.4)

### 7.1 One step

0. **Re-target** (only when a schedule moves the targets or a subdivision happened): update the
   targets, project `γ` onto them (§6.2, matrices at `γ`), check the projection path with the
   CCD of §7.4 when the line search is collision-guarded. On failure, halve the target increment
   of this step (up to 5 times), else keep the old targets and report `RetargetStalled`. No
   Armijo test: moving a target may raise the energy. (The reference mixes this into the descent
   step; D13.)
1. Evaluate `F`, `dF` (exact or Barnes–Hut), assemble the metric (dense) or build BVH, BCT and
   the multigrid hierarchy (accelerated).
2. Direction `g̃` (§6.1). Stop test (§7.5).
3. Line search (§7.2 or §7.3). Every trial `τ`: `γ_τ = γ − τ g̃`, projection (§6.2) to `γ̂_τ`,
   Armijo test `F(γ̂_τ) ≤ F(γ) − c τ dF·g̃` on the projected curve [formula: rc-armijo].
4. Accept, or `Stuck` when the trial step underflows.
5. Subdivision check (§5.4).

The step report (§9.1) records every one of these parts with its time.

### 7.2 Collision-guarded line search (the default, D11)

`τ_hi = diam(bbox γ) / max_v |g̃_v|` (no vertex moves farther than the curve is wide).
`τ_max` = first contact time of the linear motion `γ − τ g̃`, `τ ∈ [0, τ_hi]`, by the CCD of §7.4
(∞ if none). Start `τ₀ = ⅔ τ_max` (the paper's heuristic) or `τ_hi` if `τ_max = ∞`; halve until
(i) Armijo holds at `γ̂_τ`, (ii) the projection succeeds, (iii) the projection displacements
(each chord step is a linear motion) are also CCD-free. At most 40 halvings, or `τ < 1e-12 τ₀`:
`Stuck`. Boyd & Vandenberghe Alg. 9.2 with `c = 1e-4`, shrink ½. Because every accepted step
checks the whole path (move and projection), no two edges of `𝒫` and no edge and obstacle ever
cross during the flow: the isotopy class is kept (paper Sec. 5.4; the paper guards only the
move, the reference has no CCD at all).

### 7.3 Normalized Armijo line search (Sec. 5.4, applications; D12)

The paper: "normalize the gradient and perform backtracking starting with τ = 1, until both
(i) the Armijo condition is satisfied and (ii) constraint projection succeeds". Its units are
those of the positions; the reference, which produced Sec. 8, uses
(`ref: src/tpe_flow_sc.cpp`, `LineSearchStep`) `n = ‖g̃‖₂`, `τ₀ = 1/n` if `n > 1` else `1/√n`,
and if the previous accepted step `τ_prev > 1e-5`, `τ₀ = min(1.5 τ_prev, 4 τ₀)`. We use that
rule, then halve until (i) and (ii); same limits as §7.2. No collision guarantee. The Armijo
test is the standard one; the reference's `σ τ ‖g̃‖ ⟨dE, g̃⟩/(‖dE‖‖g̃‖)` (divided by `‖dE‖`) is
not used.

**Default (D11).** `CollisionGuarded`: it is the more reliable strategy (it keeps the isotopy
class; the other can tunnel), and the user's rule makes the more reliable option the default.
Confirmed by the measurement of §14.3: if it is not at least as successful as the normalized
search, or more than 1.5× slower to the 1.1×-reference energy, the numbers go to the user.
The start is a `StartRule` (§17f): `TwoThirdsTauMax` (the paper's) or `GrowPrevious { factor }`,
`τ₀ = min(⅔ τ_max, factor · τ_prev)` (the paper's rule on the first step). Both keep every trial
bounded by `τ_max` and checked by the CCD. **The default is `GrowPrevious { factor: 2 }`** (§17f).

### 7.4 Continuous collision detection (D11)

Primitives: segment–segment for every pair in `𝒫`; with obstacles also point–triangle (curve
vertex, obstacle triangle) and segment–segment (curve edge, triangle edge). Broad phase: a BVH
over the swept boxes of the motion. Narrow phase: conservative advancement (Additive CCD of Li,
Kaufman, Jiang, *Codimensional Incremental Potential Contact*, ACM TOG 40(4), 2021): the distance
between two moving primitives is Lipschitz in `τ` with constant `l` = the sum of the largest
point displacements per unit `τ` on each side, so `dist(τ) ≥ dist(τ_k) − (τ − τ_k) l`
[formula: rc-ccd-advancement]; advance `τ_{k+1} = τ_k + (dist(τ_k) − s d₀)/l` until the distance
reaches `s d₀` (gap `s = 0.1` of the pair's initial distance `d₀`) or `τ = τ_hi`; at most 100
advances, then the current `τ` (still conservative). `τ_max` = minimum over pairs. Distances:
exact segment–segment and point–triangle closest points (Ericson 2004, §5.1.9 and §5.1.5).
Works unchanged in planar mode (it uses distances, not coplanarity polynomials, which degenerate
for curves that always lie in one plane).

### 7.5 Stop rules and status

* `GradientNorm { norm, eps }`, default `‖g̃‖₂ < 1e-4` with the Euclidean norm over all
  coordinates (D14; paper Sec. 5.4 and Supp. 1.3: "the L² norm of the fractional Sobolev
  gradient below ε = 10⁻⁴"; `norm: MassL2` gives `√(Σ_v m_v |g̃_v|²)`. The reference stops when
  `dE·g̃ < 1e-4`.) Both are absolute and scale dependent (§4.5), as in the paper.
* `RelativeGradient { ratio }` (`‖g̃_k‖/‖g̃_0‖`), `EnergyStall { relative, window }`: scale-free
  rules for Blender.
* `max_iterations`, `time_limit`. A stuck line search ends the run with `Stuck` unless a schedule
  is still moving targets (the reference keeps going until the target is reached).

Statuses: `Converged`, `MaxIterations`, `TimeLimit`, `Stuck`, `Cancelled`.

## 8. Acceleration (Sec. 6, App. B)

### 8.1 Bounding volume hierarchy (B.1.1, D15)

Over the edges, each with the tangent-point `p_I = (T_I, x_I) ∈ R⁶` and mass `ℓ_I`. Top down:
the split axis cycles through the six coordinates starting with position x; the split position
along the sorted coordinates minimizes the sum of squared diameters of the two children's 6-D
boxes (paper; computed with prefix and suffix boxes; the reference minimizes squared 1-D widths
along the axis). Leaves hold at most `leaf_size` edges (default 8; reference 1). Node data
[formula: rc-bvh-lumped]: `L_N = Σ ℓ_I`, `x̄_N = Σ ℓ_I x_I / L_N`, `T̂_N = Σ ℓ_I T_I / |Σ ℓ_I T_I|`
(normalized, as the reference; the paper's `T̄_N` is not normalized, D15), the box of the
midpoints (`r_x` = half its diagonal), the box of the tangents (`r_T` = half its diagonal), and
the box of the edge **endpoints** (for the adjacency guard). Rebuilt every step; during one line
search only refitted (topology kept), so every trial energy uses the same approximation (D26).
The code is generic over the dimension and the element (§16).

### 8.2 Barnes–Hut energy and differential (B.1.2, Eq. 26, 27)

For each edge `I`, traverse from the root. An internal node `N` is *admissible* for `I` iff
[formula: rc-bh-admissibility]

```
(guard) N's endpoint box contains neither γ_{i₁} nor γ_{i₂},     r_x^N ≤ θ_x |x_I − x̄_N|,     r_T^N ≤ θ_T
```

Admissible: add `k(x_I, x̄_N, T_I) ℓ_I L_N` (Eq. 26) and stop. Leaves are always evaluated
exactly: `Σ_{J ∈ N, (I,J) ∈ 𝒫} k̂_IJ ℓ_I ℓ_J` (Eq. 18). The guard makes the exclusion of adjacent
pairs exact (any edge sharing a vertex with `I` has that vertex inside the node's endpoint box);
the paper relies on `r_x/d` being small. Differential (B.1.2): for an admissible node,
`∂/∂γ_{i_a} [ℓ_I L_N (k(x_I, x̄_N, T_I) + k(x̄_N, x_I, T̂_N))]` with `L_N, x̄_N, T̂_N` held fixed;
at a leaf, the exact `G_{I,a}` of §3.2 restricted to the leaf's edges
[formula: rc-bh-differential]. As `θ_x, θ_T → 0` no internal node is admissible and both converge
to Eq. 17 and its exact differential (§13.4). The approximation error bound is the paper's
Taylor argument, cited. Since WP3d the differential the session uses is the exact derivative of the
evaluated energy on the frozen interaction lists (the discrete adjoint, §17i); the increment above
is kept as `energy_and_differential_lumped_pairs`.

### 8.3 Block cluster tree and matrix–vector product (B.2, Eq. 28; D16)

Blocks are pairs `(𝒜, ℬ)` of nodes of the edge BVH, starting from `(root, root)`. A block is an
*admissible leaf* iff `𝒜 ≠ ℬ`, the endpoint boxes are disjoint (guard: no adjacent pair inside),
`max(r_x^𝒜, r_x^ℬ) ≤ θ_bct |x̄_𝒜 − x̄_ℬ|` and `max(r_T^𝒜, r_T^ℬ) ≤ θ_T,bct` (Eq. 28)
[formula: rc-bct-admissibility]; an *inadmissible leaf* iff `|𝒜|·|ℬ| ≤ leaf_pairs` (64) or both
are BVH leaves; otherwise it is split into the children pairs (only the non-leaf side splits).

`φ = K ψ` for `K ∈ {K^B, K⁰}`: admissible blocks use the rank-1 form `K̂_𝒜ℬ = ℓ[𝒜] κ(𝒜, ℬ)
ℓ[ℬ]ᵀ` with `κ = 2|x̄_𝒜 − x̄_ℬ|^{−(2σ+1)}` for `K^B` and `κ = k²_{2σ+5}(x̄_𝒜, x̄_ℬ, T̂_𝒜) +
k²_{2σ+5}(x̄_ℬ, x̄_𝒜, T̂_ℬ)` for `K⁰` [formula: rc-bct-lowrank]; inadmissible blocks use the
exact `K_IJ` of §4.4 with the session's quadrature. Implementation as a 0th-order fast multipole
pass: upward sums `ℓ[ℬ]ᵀψ[ℬ]` per node, one scalar per admissible block into `𝒜`, downward
accumulation times `ℓ_I`. `K𝟙` is computed once per tree. Then `Ā u` is `Dᵀ((diag(K^B𝟙) − K^B)
⊗ I₃) D u + Eᵀ(diag(K⁰𝟙) − K⁰) E u` per coordinate: 12 scalar channels (3 coordinates × 3
components of `Du`, plus 3 of `Eu`), done in one tree pass with 12-vectors. As the parameters go
to 0, the product converges to the dense `Ā` with the same quadrature.

### 8.4 Multigrid (B.3; D17)

**Coarsening** [formula: rc-coarsen]. Per component, breadth-first from its smallest vertex index
(the reference's rule, `ref: src/poly_curve_network.cpp`, `Coarsen`): a vertex is black (kept)
iff the BFS proposes it (alternating), or `deg ≠ 2`, or it carries a `Point` or `Tangent(Vertex)`
constraint, or an already explored neighbour is white; loops with fewer than 6 vertices are not
coarsened (all black), so no level gets two edges between the same vertices. Coarse vertices keep
their fine positions (paper "removing white vertices"; the reference uses `J⁺γ`, D17). Coarse
edges: fine edges with two black ends, and `(a, b)` for each white vertex with neighbours `a, b`.
**Prolongation** `J ∈ R^{|V_f|×|V_c|}`: 1 at black vertices, ½ and ½ from the two black
neighbours of a white one (paper; "linear interpolation by edge length did not improve");
`J 𝟙 = 𝟙` [formula: rc-prolongation]. Levels until `|V| ≤ coarsest_vertices` (512), no white
vertex is left, or `max_levels` (20).

**Constraints per level.** Same types on the coarse network: `Barycenter`, `TotalLength`,
`EdgeLengths` on coarse edges (when the fine level constrains all edges; for a subset, on coarse
edges whose fine edges were all constrained), `Point` (the vertex survives), `Tangent(Vertex)`
(the vertex survives, with its coarse neighbours), `Tangent(Edge)` on the coarse edge that
contains the fine edge, `Surface` on surviving constrained vertices (the reference drops surface
rows on coarse levels). Rank deficiency on a coarse level only drops rows there (it is a
preconditioner) and is reported. `P_ℓ = I − C_ℓ† C_ℓ` per level.

**Operators.** `H_ℓ = P_ℓ Ā_ℓ P_ℓ` with `Ā_ℓ` the metric of the coarse curve (rediscretized,
BCT product; not Galerkin), prolongation `Π_ℓ = P_ℓ (J_ℓ ⊗ I₃) P_{ℓ+1}`, restriction `Π_ℓᵀ`
(the reference's `MatrixProjectorOperator`). Coarsest level: dense LU of
`[Ā_L C_Lᵀ; C_L 0]`, which solves `H_L x = P_L b` in `range P_L`. V-cycle: `ν = 3` CG iterations
pre- and post-smoothing on `H_ℓ x = b` (paper: "a standard conjugate gradient smoother"),
residual, restrict, recurse, prolong, correct. Start: full multigrid (restrict the right side to
the coarsest level, solve, prolong with smoothing up to the finest; B.3.1 "Initialization").
Then V-cycles until `‖P(b − Āx)‖ ≤ tol ‖P b‖` (`tol = 1e-3`, the paper; the reference uses
`1e-2`) or 20 cycles. If three cycles in a row reduce the residual by less than 10 %, the
remaining iterations use flexible CG with the V-cycle as preconditioner (Notay 2000), reported.
[formula: rc-multigrid]

**Uses.** Gradient: Eq. 30 with `a = dFᵀ`, `g̃ = P y`. Projection: Eq. 33, one multigrid solve
per chord iteration (no factorization to reuse; Sec. 7.4).

### 8.5 Back ends

`Energy::{Exact, BarnesHut(BhOptions)}` and `Solver::{Dense, Multigrid(MgOptions)}` are chosen
independently (Fig. 20's "Barnes–Hut" curve is BH energy with the dense solve);
`Backend::Auto { dense_max_vertices }` picks exact + dense for networks of at most
`dense_max_vertices` vertices (default 2500, measured in §17l; 1500 in §17k) and BH + multigrid above, at the
session's start and again after each live subdivision. Dense solve (§17l): the Schur-complement
method of §17b, a Cholesky factorization of `A + ρ u uᵀ` (`u = 𝟙/√|V|`) and its inverse, the
`k × k` Schur complement `C (M⁻¹ ⊗ I₃) Cᵀ` (Cholesky) and a `3 × 3` block for the translations,
one step of iterative refinement per solve (`DenseMethod::Schur`, about 30× fewer flops than the
LU); it falls back to one LU with partial pivoting (`faer`, multithreaded) of the `(3|V| + k)`
saddle matrix (`DenseMethod::Lu`) when a Cholesky factor fails or is nearly singular, so singular
saddles are reported as before. The factorization is reused for Eq. 22 and every chord
iteration. The L², H¹, H² metrics use a sparse LU of their saddle matrices (`ddg_core::sparse`).
Every solver runs on the session's pools (`SolverRequest::{pool, dense_pool}`), never on rayon's
global pool. When an iterative solver fails (`MultigridDiverged`, `MultigridNotConverged`), the
step falls back to the dense solver if `|V| ≤ dense_fallback_max_vertices` (§17h).

## 9. Public API

### 9.1 Rust (`ddg_repulsive_curves`)

```rust
pub struct CurveNetwork { /* validated topology: edges, adjacency, components, degrees */ }
impl CurveNetwork {
    pub fn new(n_vertices: usize, edges: Vec<[u32; 2]>) -> Result<Self, RcError>;
    pub fn validate_positions(&self, positions: &[[f64; 3]]) -> Result<(), RcError>; // §2.3 geometric checks
    pub fn components(&self) -> &[Component]; pub fn degree(&self, v: u32) -> usize; /* ... */
}
pub fn read_curve_obj(path: impl AsRef<Path>) -> Result<(CurveNetwork, Vec<[f64; 3]>), RcError>; // v, l (polylines too)
pub fn write_curve_obj(path, network, positions) -> Result<(), RcError>;

pub struct Exponents { alpha: f64, beta: f64 }            // Exponents::new validates (D24); Default (3, 6)
pub enum GramQuadrature { Midpoint, Trapezoid }           // D2
pub struct BhOptions { pub theta_x: f64, pub theta_t: f64, pub leaf_size: usize }
pub struct BctOptions { pub theta_x: f64, pub theta_t: f64, pub leaf_pairs: usize }
pub struct MgOptions { pub bct: BctOptions, pub coarsest_vertices: usize, pub max_levels: usize,
                       pub smoothing_steps: usize, pub tolerance: f64, pub max_cycles: usize,
                       /* … §17g options … */ pub max_threads: Option<usize> }   // default Some(16) (§17l)
pub enum DenseMethod { Schur, Lu }                        // DenseSaddle::with_method; default Schur (§17l)
pub enum EnergyEval { Exact, BarnesHut(BhOptions) }
pub enum Solver { Dense, Multigrid(MgOptions), Custom(Arc<dyn SolverFactory>) }
pub enum Backend { Explicit { energy: EnergyEval, solver: Solver }, Auto { dense_max_vertices: usize } }
pub enum Metric { Sobolev { quadrature: GramQuadrature }, L2, H1, H2 }      // L2/H1/H2: evaluation (§4.6)

pub enum TangentSite { Edge(u32), Vertex(u32) }
pub enum Constraint {
    Barycenter { target: Option<[f64; 3]> },
    TotalLength { target: Option<f64> },
    EdgeLengths { edges: Option<Vec<u32>>, targets: Option<Vec<f64>> },   // None = all edges / current
    Point { vertex: u32, target: Option<[f64; 3]> },
    Surface { vertices: Option<Vec<u32>>, surface: Arc<dyn ImplicitSurface> },
    Tangent { site: TangentSite, direction: Option<[f64; 3]> },
}
pub trait ImplicitSurface: Send + Sync { fn value(&self, x: [f64; 3]) -> f64; fn gradient(&self, x: [f64; 3]) -> [f64; 3]; }
pub trait VectorField: Send + Sync { fn value(&self, x: [f64; 3]) -> [f64; 3]; fn jacobian(&self, x: [f64; 3]) -> [[f64; 3]; 3]; }
pub struct Obstacle { /* triangle soup + 3-D BVH; Obstacle::new(positions, triangles) validates */ }
pub enum Penalty {
    TotalLength { weight: f64 }, LengthDifference { weight: f64 },
    SurfacePotential { obstacle: Arc<Obstacle>, weight: f64, exponent: Option<f64>, theta: f64 },
    FieldAlignment { field: Arc<dyn VectorField>, weight: f64 },
    Hencky { weight: f64, rest_lengths: Option<Vec<f64>> },
}
pub enum LineSearch {
    CollisionGuarded { armijo: f64, shrink: f64, start: StartRule, gap: f64, max_halvings: u32 },
    NormalizedArmijo { armijo: f64, shrink: f64, max_halvings: u32 },
}
pub enum StartRule { TwoThirdsTauMax, GrowPrevious { factor: f64 } }            // guarded search only (§17f)
pub enum Method { ProjectedGradient, Lbfgs { memory: usize }, Ncg }   // Lbfgs/Ncg: soft problems only (§6.4)
pub struct ProjectionOptions { pub tolerance: Tolerance, pub max_iterations: u32 }
pub enum Tolerance { Absolute(f64), MeanEdgeLength(f64) }            // MeanEdgeLength(t) = t · ℓ̄₀
pub enum StopRule { GradientNorm { norm: GradNorm, eps: f64 }, RelativeGradient { ratio: f64 },
                    EnergyStall { relative: f64, window: usize } }
pub struct LengthSchedule { pub factor: f64, pub rate: f64 }
pub struct FlowOptions {
    pub exponents: Exponents, pub metric: Metric, pub backend: Backend, pub method: Method,
    pub constraints: Vec<Constraint>, pub penalties: Vec<Penalty>, pub planar: Option<Axis>,
    pub line_search: LineSearch, pub projection: ProjectionOptions, pub stop: StopRule,
    pub max_iterations: usize, pub time_limit: Option<Duration>,
    pub schedule: Option<LengthSchedule>, pub max_subdivisions: u32, pub threads: Option<usize>,
    pub dense_threads: Option<usize>,
    pub dense_fallback_max_vertices: usize,            // default 5000; 0: never fall back (§17h)
    pub min_edge_length_ratio: Option<f64>,            // default Some(1e-3); collapse guard (§17j)
    pub back_projection: bool,                         // default true; false: Fig. 20's "without projection" (§17l)
    pub second_order_correction: bool,                 // default false; trials γ − τp + τ²q (§17l)
    pub newton_polish: Option<f64>,                    // default Some(100); Newton steps below 100 tol (§17l)
    pub projection_early_abort: bool,                  // default true; rate test of the chord projection (§17l)
}

pub struct EnergyReport { pub energy: f64, pub penalties: Vec<f64>, pub differential: Vec<[f64; 3]>,
                          pub pairs_exact: u64, pub nodes_lumped: u64, pub seconds: f64 }
pub fn energy(network: &CurveNetwork, positions: &[[f64; 3]], exponents: Exponents,
              eval: &EnergyEval, penalties: &[Penalty], threads: Option<usize>) -> Result<EnergyReport, RcError>;
pub fn gram(network, positions, exponents, metric: &Metric) -> Result<DenseGram, RcError>; // B, B0, A (tests, oracle)

/// Step-at-a-time flow with state kept between steps (Blender's live mode).
pub struct Session { /* network, positions, targets, factorizations/trees of the current step, history */ }
impl Session {
    pub fn new(network: CurveNetwork, positions: Vec<[f64; 3]>, options: FlowOptions) -> Result<Session, RcError>;
    pub fn step(&mut self) -> Result<StepReport, RcError>;         // one step of §7.1
    pub fn run(&mut self, progress: Option<&mut dyn FnMut(&StepReport) -> Control>) -> Result<FlowReport, RcError>;
    pub fn positions(&self) -> &[[f64; 3]];  pub fn network(&self) -> &CurveNetwork;   // may change (subdivision)
    pub fn cancel_flag(&self) -> Arc<AtomicBool>;                 // checked between steps and line-search trials
    pub fn set_schedule(&mut self, s: Option<LengthSchedule>) -> Result<(), RcError>;
    pub fn status(&self) -> Status;  pub fn options(&self) -> &FlowOptions;
}
pub fn flow(network: CurveNetwork, positions: Vec<[f64; 3]>, options: FlowOptions,
            progress: Option<&mut dyn FnMut(&StepReport) -> Control>) -> Result<FlowResult, RcError>;

pub struct StepReport {
    pub iteration: usize, pub status: StepStatus,      // Accepted | Stuck | Converged | RetargetStalled
    pub energy_before: f64, pub energy_after: f64, pub penalties_after: Vec<f64>,
    pub gradient_norm: f64, pub sobolev_norm_sq: f64, pub tau: f64, pub tau_max: Option<f64>, pub halvings: u32,
    pub projection_iterations: u32, pub max_violation: f64, pub dropped_rows: Vec<usize>,
    pub linear: LinearReport,                          // iterations, residual, fallback used, levels
    pub retarget: Option<RetargetReport>, pub subdivided: bool,
    pub schedule: Option<ScheduleProgress>,            // factor asked for / imposed / measured (§5.4)
    pub solver_fallback: Option<SolverFallback>,       // the dense solver took over in this step (§17h)
    pub seconds: Timings,                              // energy, differential, metric, solve, line search, projection, ccd
}
pub struct FlowReport { pub status: Status, pub steps: Vec<StepReport>, pub final_energy: f64 }
// `step` returns the error of a failing step (status `Error(e)`); `run` never returns it: it ends
// with `status == Error(e)` (`FlowReport::error()`), the steps before it kept, the curve that of
// the last completed step. `RcError::Cancelled` gives `Status::Cancelled` in both.
```

**Errors** (`RcError`, `thiserror`-free hand-written `Display` like the other crates): the
validation errors of §2.3; `InvalidExponents`, `InvalidOption { name, value, requirement }`,
`NotPlanar { vertex, axis, value }`; constraint errors `TranslationNotFixed`,
`RedundantConstraints { rows }`, `InfeasibleConstraint { constraint, violation }`,
`InvalidTangentSite { vertex, degree }`, `ConstraintOutOfRange { constraint, index }`,
`ScheduleWithoutLengthConstraint`, `MethodNeedsSoftProblem`; solve errors `SingularMetric`,
`SingularSaddle { pivot }`, `MultigridDiverged { residuals }` (the residual grew),
`MultigridNotConverged { residuals, failure_residual }` (it stagnated above the bound),
`NonFiniteEnergy { step }`;
obstacle errors `ObstacleTriangleOutOfRange`, `DegenerateObstacleTriangle`; `Io`, `Parse { line }`;
`Cancelled`. Nothing mutates the caller's data; `Session` owns copies.

### 9.2 Python (`crates/ddg-py`, module `ddg_tools`)

Arrays in and out as elsewhere in `ddg_tools` (buffer protocol, `array.array` results):

* `ddg_tools.repulsive_energy(positions, edges, alpha=3.0, beta=6.0) -> (energy, differential)`.
* `ddg_tools.RepulsiveSession(positions, edges, *, alpha, beta, length="edges"|"total"|"free",
  grow_factor=1.0, grow_rate=0.01, max_subdivisions=0, pinned=(), tangent_vertices=(),
  obstacles=[(positions, triangles)], obstacle_weight, planar_axis=None,
  line_search="collision"|"armijo", backend="auto"|"dense"|"accelerated", threads=None,
  stop_ratio=1e-3, max_steps)` (a `#[pyclass]` holding a `Mutex<Session>`):
  `step(n=1) -> report dict` (**releases the GIL**, so Blender can run it on a worker thread),
  `snapshot() -> (positions, edges, iteration, energy, status)` (does not wait for a running step:
  it reads the last completed state, stored separately), `cancel()`, `reset()` (back to the
  initial curve), `report()`. Keys of the report dict mirror `StepReport`.
* `ddg_tools.repulsive_flow(positions, edges, steps, **options) -> (positions, edges, report)`:
  `flow` in one call.

### 9.3 Blender

Input: the active object, either an **edge-only mesh** (graphs with junctions; a mesh with faces
is an error naming the count; loose vertices are an error unless "Ignore Loose Vertices" is on)
or a **Curve** object (poly and Bézier splines are converted to polylines, written back; Q2). The
flow runs in world space (obstacles are in world space; the energy is not scale invariant, so the
object's scale matters as it would physically) and writes back in object space.

* **`object.ddg_repulsive_live`** (modal): runs the flow in the viewport until stopped. A timer
  (1/30 s) collects the last completed step from `snapshot()` and writes it to the object while a
  Python worker thread runs `session.step()`. **Esc / right click** cancels and restores the
  original object (positions and, after subdivision, topology); **Enter / left click** finishes
  and keeps the result as one undo step. The header shows iteration, energy, step size, status.
  On convergence or `Stuck` it stops stepping and waits for Enter/Esc.
* **`object.ddg_repulsive_flow`** (`REGISTER`, `UNDO`): N steps, or until convergence (`steps =
  0`, capped by `max_steps`). Adjust Last Operation reruns from the original.
* Options (both): α, β (3, 6); Length: Keep Edge Lengths (default) / Keep Total Length / Free;
  Grow/Shrink factor (1 = off) and rate; Subdivide limit; pinned points from a vertex group
  (meshes) or the selected points (curves); keep tangents at a vertex group's vertices;
  obstacles: the mesh objects of a collection, with weight; planar (object-local axis); line
  search; back end; threads. Status line: steps, energy before/after, time per step (native and
  Blender side), largest constraint violation.

Code: `blender/ddg_uv_unwrap/repulsive_operator.py` (+ `repulsive_io.py` for reading/writing
meshes and curves), registration in `__init__.py`, tests `blender/tests/test_repulsive.py`, an
end-to-end table script (§14.5).

## 10. Defaults

| option | default | reason |
|---|---|---|
| `(α, β)` | (3, 6) | paper Sec. 3.2 ("(3, 6) elsewhere"), the reference's default and all its scenes but two; the evaluation uses (2, 4.5) as Figs. 14–19 and Supp. 1.2 |
| quadrature | `Midpoint` | D2: what the reference ran (dense and BCT); 4× cheaper; the far field is midpoint anyway; confirmed in §14.3 |
| constraints | `[Barycenter, EdgeLengths(all, current)]` | Supp. 2.1: barycenter removes the translations, edge lengths make a minimizer exist; `Ê(cγ) = c^{2+α−β} Ê` falls with growth |
| line search | `CollisionGuarded { armijo 1e-4, shrink ½, start GrowPrevious { 2 }, gap 0.1, 40 halvings }` | D11 (more reliable); the start `min(⅔ τ_max, 2 τ_prev)` is measured in §17f (the paper's `⅔ τ_max` alone costs 1.5× Armijo's time per step); `c = 1e-4` the usual lenient Armijo constant (Nocedal & Wright §3.1) |
| `NormalizedArmijo` | same Armijo, reference's start rule | D12 |
| projection | `MeanEdgeLength(1e-4)`, 3 iterations, at least one step per trial; Newton polish below 100 × tolerance; early abort by the rate test; no second-order correction | paper: 1e-4 (absolute) and ≤ 3 iterations; relative to `ℓ̄₀` so that Blender scenes of any scale behave alike; the evaluation uses `Absolute(1e-4)` (D10); polish and correction measured on 24 knots in §17l |
| contact floor | `CONTACT_FLOOR = 1e-6` (of `ℓ̄₀`) | §17l: a guarded step from a curve with two disjoint edges closer than this takes no step (Stuck) |
| stop | `GradientNorm { Euclidean, 1e-4 }`, `max_iterations` 1000 | paper Sec. 5.4, Supp. 1.3 (D14); Blender uses `RelativeGradient(1e-3)` |
| Barnes–Hut | `theta_x 0.125, theta_t 0.25, leaf_size 8` | reference `diag/d < 0.25` ⇔ `r_x/d < 0.125`; tangent bound per Eq. 27; final values by §14.3 (largest θ with differential error ≤ 1e-2 and unchanged iteration counts within 5 %) |
| BCT | `theta_x 0.5, theta_t 0.5, leaf_pairs 64` | reference `sep = 1.0` on `diag/(d − r)` ≈ `r/d ≲ 0.5`; measured as above |
| multigrid | coarsest 512 vertices, 20 levels max, 12+12 CG smoothing (24 in the start), FCG, tol 1e-3, 20 cycles | paper B.3.1 (≤ 6 V-cycles to 1e-3); 512 vertices = a small dense coarsest solve, milliseconds; 12+12 with FCG measured in WP6 |
| multigrid BCT | `theta_x 0.25, theta_t 2` (position-only admissibility), `leaf_pairs 64` | §17l: the reference's rule; `theta_t 0.25` leaves the trees of the perturbed knot nearly dense; direction error against dense 2e-3–8.5e-3 (< 1e-2 of §13.4), 2.3–2.8× faster on knot 96 |
| multigrid threads | `max_threads Some(16)` | §17l: the hierarchy's products and projections stop gaining past 8–16 threads and get slower at 64 (§17k); results do not depend on it |
| dense solver | Schur complement, LU fallback | §17l: same solutions as the LU to 1e-13 relative, residual as small or smaller, bitwise across thread counts, 15× faster at 2048 vertices on one thread |
| back end | `Auto { dense_max_vertices: 2500 }` | measured alone in §17l (Schur dense solve): exact + dense and BH + multigrid cost the same per step at about 2000 vertices (8 threads) and 3000 (64 threads); exact + dense is the reference computation and wins the tie (§17k had 1500 with the LU, §17h 2000 under load) |
| edge-collapse guard | `min_edge_length_ratio Some(1e-3)` | §17j: an edge that shrinks below 1e-3 of the initial mean edge length fails the step with a named error instead of a singular metric later |
| dense fallback | `dense_fallback_max_vertices 5000` | a failed multigrid solve (§17h) is repeated with the dense solver up to this size (the Schur method holds about `4\|V\|² + k²` doubles, 1 GB at 5000; its LU fallback `(3\|V\|+k)²`, 3.2 GB) |
| schedule | none; `rate 0.01` when set | reference: total length moves by 1 % of `L⁰` per step |
| subdivision | 0 rounds | reference default |
| surface potential | `p = β − α`, `θ_M = 0.25`, leaf 8 | paper; reference uses `p = β − α`, `θ = 0.25` |
| L-BFGS | memory 30 | Supp. 3.3 |
| Hencky weight | calibrated (§14.2) | the paper gives none |
| threads | all cores | explicit pool |

## 11. Decisions: the paper as printed, the reference, and what we do

* **D1 — Energy discretization.** Paper Eq. 17–18: edge pairs, trapezoid over the four endpoint
  pairs, pairs with a common vertex omitted. Reference (`ref: src/tpe_energy_sc.cpp`, `tpe_total`,
  and its Barnes–Hut in `src/spatial/tpe_bvh.cpp`): **vertex** pairs with vertex tangents
  (normalized sum of incident edge tangents), vertex masses `Σℓ/deg`, all pairs `i ≠ j` including
  neighbours, and vertices of degree ≥ 3 contribute nothing as the tangent point. Ours: the paper
  (the user listed Eq. 17–18). The reference's discretization is not implemented (Q1); so its
  energies and differentials cannot serve as golden fixtures, its Gram matrices can (§13.7).
* **D2 — Eq. 32 vs Eq. 13/20.** No contradiction in the exponents (`k²₄ |d|^{−(2σ+1)} =
  k²_{2σ+5}`); Eq. 31 is Eq. 20–21 written with `K = w + wᵀ` and ½ (identity, checked). The real
  difference is midpoint (Eq. 32, reference) against trapezoid (Eq. 20–21) quadrature. Option,
  default midpoint, same `K` in dense and BCT (§4.4).
* **D3 — A.4.2 exponents.** Printed: "the term `1/|γ(x)−γ(y)|^{2s+1}` scales by `1/c^{2s+1}` …
  `B_σ` scales by `c²/(c²c^{2s+1})` … `L_σ` scales by `c²/c^{2s+1}`". The denominator of Eq. 12
  and 25 is `|γ(x) − γ(y)|^{2σ+1}`, so the factors are `c^{−(2σ+1)}` for `B_σ` and `B⁰_σ` and
  `c^{1−2σ}` for `L_σ`; the conclusion (same scaling of both terms) is right. Discrete version
  exact and proved (§4.5). Also A.3 writes `W^{s,p}` for `W^{s,α}`.
* **D4 — Projector.** `P := CC†`, `C† := (CCᵀ)⁻¹Cᵀ` (B.3.3) are wrong; `C† = Cᵀ(CCᵀ)⁻¹`,
  `P = I − C†C` (§6.3).
* **D5 — Translations.** `ker Ā` = translations (not per-component ones: `B⁰` couples
  components; checked). A translation-fixing constraint is required (`TranslationNotFixed`); the
  default list contains `Barycenter`. The reference adds one silently when a scene has no
  constraint.
* **D6 — Barycenter.** Length-normalized form of the paper's constraint, exact Jacobian (§5.1).
  Reference: dual-length weights `Σℓ/deg` with a constant-weight Jacobian.
* **D7 — Tangent.** Two rows in a frame orthogonal to `X` plus an orientation check, edge or
  vertex sites. Paper: three rows `T_I − X` (rank 2, singular `CCᵀ`). Reference: vertex tangents,
  three rows, and the vertex's own Jacobian block has the wrong sign
  (`ref: src/flow/constraint_functions.cpp`, `AddTangentTriplets`: `derivTangentCenter = -1 *
  vertex_tangent_wrt_vert(v, v)`, which is already the correct derivative; the blocks then do not
  sum to zero as translation invariance requires).
* **D8 — Surface constraint.** One row per constrained vertex, in order. Reference: the violation
  is written to row `rowStart + vertex index` but the Jacobian to row `rowStart + i`, which
  disagree unless all vertices are constrained (`NegativeSurfaceViolation` vs
  `AddSurfaceTriplets`).
* **D9 — Redundant rows.** Paper silent; §5.2 (two known patterns dropped and reported, others
  are errors).
* **D10 — Projection.** Chord method with `Ā`, `C` at `γ`, ≤ 3 iterations (paper, reference);
  per-row tolerances, default `1e-4 ℓ̄₀` for length rows. Paper: `10⁻⁴` absolute (Knot128 has
  `ℓ̄ ≈ 0.5`, the Freedman–He–Wang curves `ℓ̄ ≈ 0.01–0.04`). The evaluation uses the paper's
  absolute value; §14.3 compares the two.
* **D11 — Collision-guarded line search, default.** CCD also on the projection path, obstacles
  included, conservative advancement with gap 0.1; `τ_hi` when nothing collides (paper silent);
  Armijo on the projected curve. Reference: no CCD at all.
* **D12 — Normalized Armijo.** Start rule from the reference (the paper's "normalize, τ = 1"
  leaves the units open), standard Armijo (the reference divides by `‖dE‖`).
* **D13 — Length schedules.** Common factor on edge targets (reference: equal absolute
  increments), separate re-target sub-step without Armijo, CCD-guarded; subdivision as the
  reference's `subdivide_limit` with the mapping of §5.4.
* **D14 — Stop rule.** Euclidean 2-norm of `g̃` < 1e-4 (paper "L² norm", ambiguous; reference
  `dE·g̃ < 1e-4`).
* **D15 — BVH and Barnes–Hut.** Both admissibility criteria of Eq. 27 plus the endpoint-box
  guard (exact neighbour exclusion); leaves exact; `T̂_N` normalized (reference) rather than the
  paper's `T̄_N`; split by squared 6-D diameters (paper) rather than 1-D widths (reference);
  leaf size 8. Reference: vertex bodies, spatial criterion only (`diag/d < 0.25`), leaf 1.
* **D16 — Block cluster tree.** Eq. 28 both criteria plus the guard; reference: spatial only
  (`sep = 1`), and admissible blocks may contain neighbouring edge pairs (only inadmissible leaves
  skip them, `ref: src/product/block_cluster_tree.cpp`).
* **D17 — Multigrid.** Coarse vertices keep fine positions (paper; reference `J⁺γ`); small loops
  frozen; surface rows kept on coarse levels (reference drops them); tolerance 1e-3 (paper;
  reference 1e-2); stagnation fallback to flexible CG (ours).
* **D18 — Dimension.** R³ with `Planar { axis }` (§2.2); reference: plane constraint `yplane`.
* **D19 — Surface potential.** Point–triangle distance near field, centroid lumping far field
  (§5.3). Reference: unweighted per curve vertex, mesh vertices as point masses, and its energy
  omits the node mass while its force includes it (`ref: src/obstacles/mesh_obstacle.cpp`,
  `bodyEnergy` vs `bodyForce`).
* **D20 — Penalties as printed in Sec. 8.1.** Reference differs: total length `½wL²` with a
  gradient indexing bug (`GetVertex(j)` for the neighbour), length difference with ½, field
  energy without the factor `ℓ_I` while its gradient has it (`ref: src/extra_potentials.cpp`).
* **D21 — Soft problem of Sec. 7.** Hencky as in §5.3; weight calibrated (§14.2); translations
  removed from the direction by `M_c` (Supp. 2.2).
* **D22 — L-BFGS / NCG.** Memory 30, dynamic metric (Supp. 3.3), pairs with `sᵀy ≤ 0` skipped;
  PR+ with restart. Line search as §7.2 (Supp.: same as the main article).
* **D23 — Masses.** `m_v = ½Σ_{I∋v} ℓ_I` (the standard lumped mass); reference `Σℓ/deg`. L² =
  lumped mass, H¹ = cotan-free curve Laplacian `1/ℓ`, H² = `LM⁻¹L` (Supp. 2.1).
* **D24 — Exponents and the kernel at `ν = 0`.** Exponents outside Blatt's range are errors (no
  clamping). `ν^{α−1}n̂ := 0` at `ν = 0` (continuous for `α > 1`); reference cuts at `1e-10`.
* **D25 — Determinism.** Results are bitwise independent of the thread count.
* **D26 — Barnes–Hut in line searches.** Tree topology frozen during one line search (refit).
* **D27 — Input validation** of §2.3 (the reference trusts its input).
* **D28 — Cosmetic typos.** Sec. 5.3.2 "Φ(g) = 0" means `Φ(γ̃) = 0`; Sec. 5.2.2 says the weights
  come from the denominator "in Equation 25", which is that of Eq. 12 as well (same exponent).
* **D29 — Layout.** Interleaved coordinates, `Ā = A ⊗ I₃` (the paper displays blocks).
* **D30 — Fig. 20 setup.** The paper does not say which curve or exponents; §14.1 uses the
  Freedman–He–Wang curve and knot 96 at several resolutions, with both exponent pairs.

## 12. Formula ids

Rust functions carry `[formula: <id>]`; each id is a Lean theorem in
`lean/proofs/DdgProofs/RepulsiveCurves/*.lean` ("proved") or a heading in `lean/proofs/CITED.md`
("cited"). Planned Rust locations refer to the module layout of §15.

| id | statement | Rust (planned) | Lean / CITED |
|---|---|---|---|
| `rc-edge-quantities` | §2.1 | `network` | proved, `Basic.lean` (definitions + `\|T\| = 1`) |
| `rc-dual-length` | `Σ_v m_v = Σ_I ℓ_I` | `network` | proved, `Basic.lean` |
| `rc-kernel` | Eq. 3; `\|T×d\|² = \|d\|² − (T·d)²`; even in `T`; rigid invariance; `k(cp,cq,T) = c^{α−β}k` | `kernel` | proved, `Kernel.lean` |
| `rc-kernel-gradient` | §3.2 gradients for `ν, ρ > 0`; continuity at `ν = 0` for `α > 1` | `kernel` | proved, `KernelGradient.lean` (if the rpow calculus proves heavy: FD + oracle only, no tag) |
| `rc-discrete-kernel` | Eq. 18 | `energy::exact` | proved (definition), `Energy.lean` |
| `rc-energy` | Eq. 17; finite iff no coincident vertices across `𝒫` | `energy::exact` | proved, `Energy.lean` |
| `rc-energy-differential` | §3.2 pair formulas = chain rule | `energy::exact` | proved, `Energy.lean` (same caveat) |
| `rc-energy-scaling` | `Ê(cγ) = c^{2+α−β}Ê` | `energy` (doc) | proved, `Energy.lean` |
| `rc-gradient-identities` | invariance ⇒ `Σg = 0`, `Σγ×g = 0`; homogeneity ⇒ Euler identity | `energy` (doc) | proved, `Energy.lean` (reuses `dev-gradient-invariance`) |
| `rc-exponent-range` | Blatt: finite energy ⇔ `W^{s,α}` for `α > 1, β ∈ [α+2, 2α+1)` | `Exponents::new` | cited |
| `rc-differential-order` | `dℰ` has order `2s` (A.3) | — (docs) | cited |
| `rc-fractional-operator` | Eq. 24/25: `(−Δ)^σ`, `L_σ` of order `2σ` | — (docs) | cited |
| `rc-quadrature-consistency` | omitting adjacent pairs is insignificant under refinement (Sec. 5.1, Fig. 9) | — (docs) | cited |
| `rc-derivative-operator` | `D_I`, `E` | `metric` | proved, `Gram.lean` |
| `rc-high-order-gram` | Eq. 21 increments = quadratic form | `metric::sobolev_dense` | proved, `Gram.lean` |
| `rc-low-order-gram` | Eq. 20 increments = quadratic form | `metric::sobolev_dense` | proved, `Gram.lean` |
| `rc-kernel-matrix-form` | Eq. 31 identity; Eq. 32 = midpoint weights | `metric`, `spatial::bct` | proved, `Gram.lean` |
| `rc-gram-psd-kernel` | PSD; kernels of §4.5; `ker(B+B⁰) = ker B ∩ ker B⁰` | `metric` (doc) | proved, `Gram.lean` |
| `rc-gram-scaling` | `B, B⁰ ∝ c^{−(2σ+1)}` | `metric` (doc) | proved, `Gram.lean` |
| `rc-gradient-scaling` | projected direction scales by `c^{2σ+2+α−β}` | `step` (doc) | proved, `Saddle.lean` |
| `rc-classic-metrics` | §4.6; `ker L` = per-component constants | `metric::sparse_classic` | proved, `Gram.lean` |
| `rc-saddle-gradient` | Eq. 22 ⇔ `Ā`-closest point in `ker C`, unique | `solve` | proved, `Saddle.lean` |
| `rc-descent-direction` | `dF·g̃ = g̃ᵀĀg̃ > 0` unless `g̃ = 0` | `step` | proved, `Saddle.lean` |
| `rc-chord-projection` | one iteration = `Ā`-least-norm solution of the linearized constraints | `solve` | proved, `Saddle.lean` |
| `rc-chord-convergence` | local linear convergence of the chord method | `step::projection` (doc) | cited (Ortega & Rheinboldt 1970, Ch. 10, Ostrowski's theorem; Kelley 1995, Ch. 5, the chord method) |
| `rc-projector` | `C†`, `P` properties | `solve::projector` | proved, `Saddle.lean` |
| `rc-projected-system` | Eq. 30 ⇔ Eq. 29 | `multigrid` | proved, `Saddle.lean` |
| `rc-projection-transform` | Eq. 33 | `multigrid` | proved, `Saddle.lean` |
| `rc-mean-free-direction` | §6.4 direction rows | `optimize` | proved, `Saddle.lean` |
| `rc-constraint-barycenter` | Jacobian of §5.1 | `constraints` | proved, `Constraints.lean` |
| `rc-constraint-total-length` | ″ | `constraints` | proved, `Constraints.lean` |
| `rc-constraint-edge-length` | ″ | `constraints` | proved, `Constraints.lean` |
| `rc-constraint-point` | ″ | `constraints` | proved, `Constraints.lean` |
| `rc-constraint-surface` | ″ | `constraints` | proved, `Constraints.lean` |
| `rc-constraint-tangent` | Jacobian; `P_T T = 0` (rank 2) | `constraints` | proved, `Constraints.lean` |
| `rc-tangent-frame` | `(e₁·s, e₂·s) = 0 ⇔ s ∥ X` | `constraints` | proved, `Constraints.lean` |
| `rc-penalty-total-length` | value/gradient | `penalties` | proved, `Penalties.lean` |
| `rc-penalty-length-difference` | ″ | `penalties` | proved, `Penalties.lean` |
| `rc-penalty-field` | ″ (§5.3 formula) | `penalties` | proved, `Penalties.lean` |
| `rc-penalty-hencky` | ″ | `penalties` | proved, `Penalties.lean` |
| `rc-penalty-surface` | chain rule given `∇δ` | `penalties::surface` | proved, `Penalties.lean` |
| `rc-distance-gradient` | distance to a closed convex set is `C¹` outside it, `∇δ = (x − π(x))/δ` | `obstacle` | cited (Bauschke & Combettes 2017, Ch. 12: `½δ²` of a closed convex set is differentiable with gradient `x − π(x)`; Ericson §5.1.5 for the closest point) |
| `rc-closest-points` | segment–segment and point–triangle closest points (Ericson §5.1.9, §5.1.5) | `ccd`, `obstacle` | cited |
| `rc-length-schedule` | common factor keeps ratios; total target formula | `step::schedule` | proved, `Schedule.lean` |
| `rc-subdivision-constraints` | §5.4 mapping keeps every constraint but new surface rows | `subdivide` | proved, `Schedule.lean` |
| `rc-armijo` | finite termination for a descent direction | `step::line_search` | cited (Nocedal & Wright §3.1 and Alg. 3.1; Boyd & Vandenberghe Alg. 9.2) |
| `rc-ccd-advancement` | Lipschitz lower bound ⇒ no contact before the advanced time | `ccd` | proved, `Ccd.lean` |
| `rc-accd` | additive CCD algorithm | `ccd` | cited (Li, Kaufman, Jiang 2021) |
| `rc-isotopy` | a CCD-free piecewise linear path of embeddings is an isotopy | `step` (doc) | cited (paper Sec. 5.4; standard) |
| `rc-bvh-lumped` | `L_N, x̄_N, T̂_N`; parent = mass-weighted combination of children | `spatial::bvh` | proved, `Hierarchy.lean` |
| `rc-bh-admissibility` | Eq. 27 error bound | `energy::barnes_hut` | cited (paper B.1.2; Barnes & Hut 1986) |
| `rc-bh-differential` | B.1.2 increment = derivative of the lumped term | `energy::barnes_hut` | proved, `Hierarchy.lean` |
| `rc-bh-adjoint` | differential of the evaluated BH energy on frozen lists (chain rule through `L_N`, `x̄_N`) | `energy::barnes_hut` | proved, `BarnesHutAdjoint.lean` |
| `rc-bh-adjoint-shift` | re-centred downward pass of the node sums | `energy::barnes_hut` | proved, `BarnesHutAdjoint.lean` |
| `rc-bct-admissibility` | Eq. 28, low-rank error | `spatial::bct` | cited (Hackbusch 2015, admissibility and low-rank approximation of asymptotically smooth kernels; Greengard & Rokhlin 1997) |
| `rc-bct-lowrank` | rank-1 block = lumped kernel; exact as the clusters shrink | `spatial::bct` | proved, `Hierarchy.lean` (algebraic part) |
| `rc-coarsen` | black/white rule: no two adjacent whites; whites have two black neighbours | `multigrid::coarsen` | proved, `Multigrid.lean` |
| `rc-prolongation` | `J𝟙 = 𝟙`; `Π_ℓ = P_ℓ (J_ℓ ⊗ I₃) P_{ℓ+1}` has its range in `range P_ℓ = ker C_ℓ` | `multigrid::coarsen` | proved, `Multigrid.lean` |
| `rc-multigrid` | V-cycle/FCG convergence for this operator | `multigrid` | cited (Braess & Sarazin 1997; Ainsworth & Glusa 2017; Notay 2000; no proof for `L_σ`-type operators) |
| `rc-lbfgs-two-loop` | two-loop = `H_k r` | `optimize::lbfgs` | cited (reuse `dev-lbfgs-two-loop`) |
| `rc-ncg-polak-ribiere` | PR+ restart keeps descent | `optimize::ncg` | cited (Polak & Ribière 1969; Nocedal & Wright §5.2) |

**Oracle** (§13.6): energy, differential, Gram matrices, one projected step.

## 13. Verification

### 13.1 Layers (CLAUDE.md, adapted)

There is no certifiable minimizer (the problem is nonconvex). Layers: (1) finite differences;
(2) invariances and exact scaling; (3) algebraic properties of the metric and of the solves;
(4) dense against accelerated as the parameters go to 0; (5) constraint and CCD guarantees;
(6) the Lean oracle (differential tests, golden fixtures); (7) the Lean proofs; (8) golden
fixtures from the reference program; (9) end-to-end runs (§14).

### 13.2 Finite differences

Central differences, `h = 1e-6 ℓ̄`, relative error ≤ 1e-6 on: `Ê` for (3,6), (2,4.5), (4,8) on
random closed curves, open curves, a 3-junction graph, two linked loops; every penalty; every
constraint Jacobian; the Barnes–Hut differential against the derivative of its own lumped terms
(fixed tree); `k` near `ν = 0` (no NaN; one-sided limits agree).

### 13.3 Invariances and scaling

`Σ_v ∂Ê/∂γ_v = 0`, `Σ γ_v × ∂Ê/∂γ_v = 0` to 1e-12 relative; Euler identity; `Ê(cγ) =
c^{2+α−β}Ê(γ)` and `A(cγ) = c^{−(2σ+1)}A(γ)` to 1e-13 for `c ∈ {1e-3, 0.7, 1e3}`; the projected
direction for `cγ` equals `c^{2σ+2+α−β}` times that for `γ` (barycenter + edge lengths); energy
unchanged by rigid motions, reflection, edge re-orientation, vertex renumbering.

### 13.4 Metric and solves

`B`, `B⁰` symmetric; increments equal the quadratic forms (random `u, v`); smallest eigenvalue
`≥ −1e-12 λ_max`; dimension of the numerical kernels as §4.5 on one and two components; the
collinear case reports `SingularMetric`. Eq. 31 form equals the increments (midpoint).
**Dense vs accelerated:** with `θ → 0` (θ ∈ {0.5, 0.25, 0.1, 0.03, 0}), Barnes–Hut energy and
differential and the BCT product converge to the dense ones and are equal (≤ 1e-12 relative) at
θ = 0; the error decreases monotonically on average and is ≤ 1e-2 at the defaults (recorded).
Multigrid: `‖P(b − Āx)‖ ≤ tol‖Pb‖` at exit; its `g̃` agrees with the dense one to the expected
order of `tol`; `C g̃ = 0` to 1e-10 relative in both back ends.

### 13.5 Constraints, line search, CCD

Projection: after success every row within tolerance; `C g̃ = 0`; `P² = P`, `CP = 0` to 1e-12.
Rank checks produce the errors and drops of §5.2 on constructed cases. CCD: segment pairs with
known contact times (parallel, crossing, skew, touching at endpoints, planar), `τ_max ≤` the true
contact time always (randomized: 10⁵ pairs against a high-precision bisection); every accepted
collision-guarded step leaves the minimum distance over `𝒫` positive; an adversarial step that
would tunnel through (a large `τ` through a thin gap) is stopped; the normalized search on the
same case is allowed to tunnel and the audit reports it. Determinism: bitwise identical results
for 1, 7 and 64 threads.

### 13.6 Lean oracle (`lean/oracle`, no Mathlib; `DdgOracle/RepulsiveCurves.lean` and
`DdgOracle/RepulsiveCurves/`)

Commands `rc [IN [OUT]]`, `rc-selftest`, `rc-fixtures write|check DIR`; text format v1 with hex
floats as for `dev`. Tasks:

* `energy`: network, positions, `(α, β)` → `Ê`, per-vertex differential.
* `gram`: → `B`, `B⁰`, `A` (dense, both quadratures), `K^B`, `K⁰`.
* `saddle`: + constraints (barycenter, total length, edge lengths, point, tangent edge/vertex,
  surface = sphere) → `C`, `Φ`, `g̃`, `λ` (dense LU of `Dense.lean`).
* `step`: one normalized-Armijo step with given parameters → the trial `τ` sequence, the chord
  iterations and violations, the accepted `τ` and positions.
* `penalties`: values and gradients (an obstacle of a few triangles).

Fixtures `lean/oracle/fixtures/rc/`: circle (analytic symmetry), trefoil (80 vertices), open
curve with tangent sites, 3-junction graph, two linked loops, a planar case. Rust test
`tests/oracle.rs`: relative agreement ≤ 1e-12 (energy, Gram), ≤ 1e-10 (solves), identical trial
sequences.

### 13.7 Golden fixtures from the reference program (once it builds)

Comparable with ours: the dense Gram matrix (`SobolevCurves::SobolevGramMatrix` is our midpoint
`A` exactly: same neighbour exclusion, midpoint weights, `½(k_IJ + k_JI)` low-order kernel), the
edge-length, total-length and point Jacobians, the coarsening and prolongation (if our BFS
order is the reference's, which §8.4 copies), BCT products only against the reference's own
dense matrix (its admissibility differs). Not comparable: energies and differentials (D1),
barycenter, tangent and surface constraints (D6–D8), penalties (D19, D20). Fixtures in
`crates/ddg-repulsive-curves/tests/fixtures/reference/`, produced by a headless driver that
links the reference (kept outside this repository, with the build).

### 13.8 Lean proofs and citations

Proved (`lean/proofs/DdgProofs/RepulsiveCurves/`): every "proved" row of §12. Notable: the
increments of `B` and `B⁰` (finite-sum algebra over ordered pairs with abstract weights), PSD and
kernels, `ker(B+B⁰) = ker B ∩ ker B⁰` for PSD matrices, the Eq. 31 identity, the projector
identities and Eq. 30/33, the saddle optimality and descent property, the scaling exponents
(energy and Gram, as exact homogeneity), invariance ⇒ gradient identities (reusing
`dev-gradient-invariance`), the constraint Jacobians and the tangent frame, the penalties'
gradients, the CCD advancement lemma, the coarsening and prolongation facts, the schedule and
subdivision mappings. Cited (`lean/proofs/CITED.md`, section "Repulsive Curves"): every "cited"
row of §12 — Blatt's lemma, the order of the differential and of `L_σ`, the quadrature
consistency, Barnes–Hut and H-matrix error bounds, multigrid convergence, chord-method
convergence, Armijo termination, L-BFGS and PR+, ACCD, closest-point algorithms, the regularity
of the distance function, the isotopy statement.

## 14. Measurement plan

Single-thread timings and benches run alone; independent runs (per knot, per scene, sweeps) run
in parallel on the 64 cores, one log per job under `runs_tmp/`.

### 14.1 Fig. 20 scaling (D30)

Curves: Freedman–He–Wang 600/1200/2400 edges and uniform refinements to 4800, 9600, 19200;
knot 96 (`Knot128/perturbed/knot0096.obj`, the Supp. multi-resolution knot) refined to
1024–16384. Exponents (3, 6) and (2, 4.5); barycenter + edge lengths. Configurations: Exact
(exact energy, dense solve), Barnes–Hut (BH energy, dense solve), Multigrid (BH + BCT +
multigrid) with and without projection. Metric: mean time per iteration over up to 500
iterations or until convergence (paper). Threads 1 and 64. Paper values read from Fig. 20:
Exact ≈ 6 s at 1000 vertices, ≈ 145 s at 4000; Barnes–Hut ≈ 10 s at 2000, ≈ 60 s at 4000;
Multigrid ≈ 6 s at 2000, ≈ 25 s at 8000, ≈ 57 s at 16000; without projection ≈ 3 s at 2000,
≈ 38 s at 16000. **Gate:** single-thread at or below the paper at every size; report the
64-thread speed-up, memory, the dense/accelerated crossover (sets `Auto`), and the time split
(energy, metric, solve, line search, projection, CCD).
Results: §17k.

### 14.2 Sec. 7 comparison

Data: the 228 knots (Knot128 perturbed + Trefoil100), (α, β) = (2, 4.5), `Absolute(1e-4)`
projection tolerance, stop ‖g̃‖₂ < 1e-4. **Hard** (barycenter + edge lengths, collision-guarded
line search): Hs projected gradient (reference method), L², H¹, H² projected gradients, plus
Hs with the normalized search (for §14.3). **Soft** (`Ê + Hencky`): Hs NCG (reference method),
Hs L-BFGS(30), and for context H¹ NCG and H¹ L-BFGS. Per knot and method: iterations and time to
reach 1.1× the reference method's minimal energy; failures split into *stuck* (line search
collapses) and *nonconvergent* (budget exhausted); knots on which every method is stuck are
excluded and counted, as in Fig. 19 / Supp. Fig. 2, 5. Output: the four scatter plots and the
failure table, next to the paper's numbers (hard: Hs 6.5 % stuck; H¹ 64.5 % failed; H² 92.1 %;
L² ≥ 97.7 %; soft: Hs NCG 3.3 % stuck, Hs L-BFGS 14.8 %).

*Hencky calibration:* the smallest weight `w_H ∈ {10^k}` for which the median, over 10 knots,
of the final largest `|ln(ℓ/ℓ⁰)|` under Hs NCG is ≤ 1e-2, fixed before the full run. "Final" means at the soft minimizer, not after a fixed number of iterations: §17m measures it
through the edge tensions of the hard reference runs (strain ≈ `max_I λ_I/(2w)`), which gives `w_H = 10`.

*Budget (Q4):* phase A — all 228 × methods in parallel, one thread per job, with an iteration
cap (the paper's 24-minute cap is wall time on a 2019 Xeon), recording iterations and energies;
phase B — timed alone, a stratified subset (12 Knot128 + 12 Trefoil100) × all methods, 1 and 64
threads, for the time scatter plots.

### 14.3 Defaults confirmed by measurement

* Line search (D11): collision-guarded vs normalized on the 228 knots (Hs projected gradient)
  and the application scenes: stuck rate, isotopy violations of the normalized search (audit
  CCD on its accepted steps), iterations and time to 1.1× reference energy.
* Quadrature (D2): midpoint vs trapezoid on 40 knots: iterations, stuck rate, time.
* Projection tolerance (D10): `MeanEdgeLength(1e-4)` vs `Absolute(1e-4)`: projection failures,
  backtracks.
* Barnes–Hut and BCT parameters, leaf sizes {1, 4, 8, 16, 32}, multigrid smoothing steps
  {1, 2, 3, 5} and coarsest size {128, 512, 2048}: error vs time curves; pick by the rule in §10.

### 14.4 Applications (examples, tests, showcase)

A reader for the reference's scene format (`scenes/FORMATS.md`, the keywords used by these
scenes) runs the same setups. Per scene: no crossing (CCD audit of every step), constraint
violation within tolerance, energy decrease (except re-target steps), time; renders to
`showcase/repulsive-curves/`.

| application | figure | reference scene |
|---|---|---|
| packing in a domain (surface penalty, growing edges, subdivision) | 1, 23 | `SimpleCurvePacking` (circle, rect), `BunnyCurvePacking`, `RingsPacking` |
| packing on a surface (implicit constraint, growing) | 24 | `ImplicitTorus` (two), `ImplicitSphere`, `ImplicitDoubleTorus` (two) |
| shrinking on a surface (Dehn–Thurston curve) | 25 | `DoubleTorus` (`fix_edgelengths 0.76`) |
| disks with a fixed boundary curve (planar) | 22 | `BluenoiseCurves` (plane `y = 0`, obstacle outline, edges ×1.5) |
| graph drawing, junctions, isometric embedding | 26, 27 | `GraphDrawing` (K5, K3,3), `PlanarGraphDrawing` |
| vascular network | 28 | `CapillaryNetwork` (endpoints pinned, edge lengths, hand obstacle) |
| self-avoiding splines: points, growing length | 30 | `CurveInterpolation` ((4,8), length and length-difference penalties, pins) |
| tangent constraints | 29 | `TangentInterpolation` (special vertices pinned with tangents, length penalty) |
| muscle fibres (partial surface constraint) | 21 | `PartialSurfaceConstraint` (only the endpoints on spheres) |
| untangling | 2, 15, 33 | `FreedmanUnknot`, `CableUntangle` |

### 14.5 Blender end-to-end table

Scenes: Freedman 1200, a Knot128 member, K5 graph (mesh), a packing scene with an obstacle,
a Bézier spline with pinned control points. For `ddg_repulsive_flow` with N = 50 and for 10 s of
`ddg_repulsive_live`: native time per step, Blender-side time (read, write, redraw), steps per
second in live mode, final energy, largest violation. Script pattern of `scripts/e2e_table.py`.

### 14.6 Thread scaling

1, 2, 4, 8, 16, 32, 64 threads for: exact energy + differential, Barnes–Hut, dense Gram
assembly, dense LU, BCT product, multigrid solve, CCD; curves of 2k, 8k and 32k vertices;
efficiency table; bitwise equality of the results across thread counts.

### 14.7 Black-box comparison with the authors' program

Once the reference builds (another agent): a headless driver runs its scenes with dense,
Barnes–Hut and multigrid; its built-in performance log gives the per-iteration split. Compare:
time per iteration on the §14.1 curves (`OMP_NUM_THREADS` 1 and 64, alone), iterations to its
stop rule, and final energies, all evaluated with our Eq. 17 on its output curves (its own
energy is the vertex discretization, D1).

## 15. Work packages

Each package is sized for one agent. "Files" are the only paths it may edit; module files are
created as stubs by WP1 so that parallel packages never touch the same file. `src/lib.rs` (module
list, re-exports) is edited only by WP1 and by the orchestrator between rounds. The crate's
`Cargo.toml` dependencies are set by WP1 (`ddg-core`, `faer` 0.24, `rayon`); `[[bench]]` entries
by WP14. A package that needs a change outside its files reports it instead.

| WP | content | files | depends on | model |
|---|---|---|---|---|
| WP1 | core: network + validation, OBJ curves, kernel, exact energy & differential, errors, module stubs | `Cargo.toml`, `src/{lib,error,network,io,kernel}.rs`, `src/energy/{mod,exact}.rs`, stubs of all other modules, `tests/energy.rs` | — | Sonnet |
| WP2 | metrics (Sobolev dense both quadratures, L²/H¹/H²), constraints + rank checks, dense & sparse saddle solves, projector, `SaddleSolver` trait | `src/metric/{mod,sobolev_dense,sparse_classic}.rs`, `src/constraints/{mod,jacobian,rank}.rs`, `src/solve/{mod,dense,sparse,projector}.rs`, `tests/{metric,constraints}.rs` | WP1 | Opus |
| WP3 | session, step anatomy, both line searches, CCD, schedules, subdivision, stop rules, L-BFGS/NCG | `src/step/{mod,line_search,schedule,stop}.rs`, `src/ccd.rs`, `src/subdivide.rs`, `src/optimize/{mod,lbfgs,ncg}.rs`, `tests/{flow,ccd}.rs` | WP1, WP2 | Opus |
| WP4 | generic BVH, Barnes–Hut energy and differential | `src/spatial/{mod,aabb,bvh}.rs`, `src/energy/barnes_hut.rs`, `tests/barnes_hut.rs` | WP1 | Sonnet |
| WP5 | block cluster tree, hierarchical `Ā` product | `src/spatial/bct.rs`, `src/metric/hierarchical.rs`, `tests/bct.rs` | WP4, WP2 | Opus |
| WP6 | multigrid: coarsening, levels with constraints, V-cycle, FMG, FCG fallback, Eq. 30/33 solver | `src/multigrid/{mod,coarsen,vcycle}.rs`, `tests/multigrid.rs` | WP5, WP2 | Opus |
| WP7 | penalties, obstacles (triangle BVH, point–triangle distance), implicit surfaces, vector fields | `src/penalties/{mod,length,field,surface,hencky}.rs`, `src/obstacle.rs`, `src/surfaces.rs`, `tests/penalties.rs` | WP1, WP4 (BVH) | Sonnet |
| WP8 | Lean proofs + CITED | `lean/proofs/DdgProofs/RepulsiveCurves/*.lean`, `lean/proofs/DdgProofs/RepulsiveCurves.lean`, one import line in `lean/proofs/DdgProofs.lean`, a section in `lean/proofs/CITED.md` | SPEC | Opus |
| WP9 | Lean oracle + fixtures; Rust differential test | `lean/oracle/DdgOracle/RepulsiveCurves.lean`, `lean/oracle/DdgOracle/RepulsiveCurves/*.lean`, dispatch lines in `lean/oracle/Main.lean`, `lean/oracle/fixtures/rc/`, `tests/oracle.rs` | SPEC; Rust test after WP3 | Sonnet |
| WP10 | Python bindings | `crates/ddg-py/src/repulsive.rs`, `crates/ddg-py/src/python/repulsive.rs`, registration lines in `crates/ddg-py/src/{lib,python}.rs`, dependency line in `crates/ddg-py/Cargo.toml`, `ddg_tools.pyi`, README section, `python-tests/test_repulsive.py` | WP3, WP7 | Sonnet |
| WP11 | Blender operators | `blender/ddg_uv_unwrap/{repulsive_operator,repulsive_io}.py`, registration lines in `__init__.py`, `blender/tests/test_repulsive.py`, README section | WP10 | Sonnet |
| WP12 | evaluation (Sec. 7) | `examples/eval_knots.rs`, `scripts/rc-eval.sh`, `scripts/rc_eval_table.py`, results section of this file and README | WP3 (dense suffices) | Sonnet |
| WP13 | applications and showcase | `src/scene.rs`, `examples/scenes.rs`, `tests/applications.rs`, `showcase/repulsive-curves/` | WP3, WP7 (WP6 for the large packing scenes) | Sonnet |
| WP14 | benches: Fig. 20, thread scaling, components | `benches/{fig20,threads,components}.rs`, `[[bench]]` in `Cargo.toml` | WP4–WP6 | Sonnet |
| WP15 | reference parity fixtures | `tests/reference_fixtures.rs`, `tests/fixtures/reference/` (the driver lives with the reference build, outside the repo) | reference build, WP2, WP6 | Sonnet |

Rounds with at most three agents (timed work of WP12 phase B and WP14 must run alone):

1. WP1, WP8, WP9 (Lean side).
2. WP2, WP4, (WP8/WP9 continue).
3. WP3, WP5, WP7.
4. WP6, WP12 phase A, WP9 (Rust side).
5. WP10, WP13, WP15 (if the reference builds).
6. WP11; then, alone: WP14, WP12 phase B, the §14.3 sweeps that time things.

## 16. Shared infrastructure

* **What `ddg-core` gives.** OBJ reading/writing for polygon meshes (`ddg_core::io`): it ignores
  `l` lines, so curves get a small reader here (`v`, and `l` with two or more indices); obstacle
  meshes are read with `ddg_core::io::read_obj` (polygons fan-triangulated by the obstacle loader,
  which needs no manifold check). Sparse matrices and factorizations (`ddg_core::sparse`: `Llt`
  for `C Cᵀ`, `Lu`/`Ldlt` for the sparse saddles of §4.6) and the explicit `ThreadPool`. Dense
  factorizations come from `faer` (the version `ddg-core` uses). There is no BVH, no tree, no
  hierarchical matrix and no multigrid anywhere in the workspace (searched).
* **Reuse by later papers** (`../RepulsiveSurfaces`, `../RepulsiveShells`): the surface version
  needs the same three tools on faces (a BVH over `(x, n) ∈ R⁶`, a block cluster tree for
  fractional kernels, multigrid with constraint projection). Recommendation: write
  `spatial::bvh` generic over the dimension (const generic) and over an element trait (box,
  centroid, mass, payload), `spatial::bct` generic over a cluster-kernel trait, and the V-cycle
  over a `Level` trait (matrix-vector product, smoother, prolongation/restriction, projector,
  coarse solve), with no curve types in those signatures. Keep them in this crate now and move
  them into a shared crate (e.g. `ddg-hierarchy`) when the surfaces paper starts; do not build
  abstractions curves do not need. Coarsening, kernels and constraints stay curve specific.

## 17. Open items

Things the paper and the reference leave open, with how they are settled:

1. Hencky weight (Supp. 2.2 gives none): calibrated by the rule of §14.2.
2. The curve and exponents of Fig. 20: not stated; §14.1 covers two curves and both exponent
   pairs.
3. Multigrid internals of the reference live in `libgmultigrid`, which was not available
   (`deps/` empty): smoothing steps, coarsest size and V-cycle stopping are ours and measured;
   when the library is fetched, check its licence, then compare.
4. The Armijo constants of the Mathematica evaluation (Sec. 7) are not given: we use
   `c = 1e-4`, shrink ½ for all methods.
5. "L² norm of the gradient" in the stop rule (Euclidean vector norm or the function's L² norm):
   Euclidean by default (D14), both available.
6. The 24-minute budget of Sec. 7 is wall time on 2019 hardware with Mathematica: iteration caps
   for the parallel phase, timed subset alone (Q4).
7. The reference's scene files assume all constraints hold initially; so do we (§5.1, error
   otherwise).

## 17b. Noted after WP2 (orchestrator)

* §8.5: faer's parallel triangular solve is not bitwise deterministic across thread counts, so
  the dense solve runs sequentially (D25 holds); its LU stops scaling past 16-32 threads, so the
  dense back end sizes its own pool by measurement (WP14), not 64 by default.
* Possible speed-up, not adopted yet: a Schur-complement solve that uses `Ā = A ⊗ I₃` (Cholesky
  of `A + ρ𝟙𝟙ᵀ`, bordered by the translations, plus a dense `k × k` system; about 6x fewer
  flops). Decided by measurement in WP14 against the multigrid back end.
* §5.2: rank is checked on the row-normalized `C Cᵀ` by inverse iteration (threshold 1e-10),
  since a plain sparse Cholesky can succeed with tiny pivots.

## 17f. Noted after WP3 (orchestrator)

*(The Barnes–Hut noise floor and `Stuck` ending noted below are fixed by §17i: the differential is now the derivative of the evaluated energy and Freedman 600 converges.)*

* **Line-search start rule (D11).** Measured on 24 perturbed Knot128 knots (every 5th: 1, 6, ...,
  116; 1500-2750 vertices), `(α, β) = (3, 6)`, `Absolute(1e-4)`, 3 projection iterations, exact
  energy + dense solver, 100 iterations, one process per knot and variant (8 threads, 8 processes
  at a time; the machine was shared, load average 55-75, so wall-clock ratios are approximate).
  Logs, driver and table: `runs_tmp/rc-wp3b-ls/` (`table.md`; driver `tests/flow.rs`
  `ls_compare_one_knot`, `RC_LS=guarded|grow|armijo`). `E*` is the best energy any variant reached
  on that knot.

  | search | stuck steps | isotopy violations (CCD audit) | trials/step | s/step (median) | paired s/step vs Armijo | reached 1.1 E* | median it, s to 1.1 E* | median final E / E* |
  |---|---|---|---|---|---|---|---|---|
  | guarded, paper start | 0 / 2400 | 0 | 6.0 | 21.0 | 1.54 (1.24-1.93) | 18 / 24 | 71, 1530 | 1.011 |
  | guarded, grow-previous | 0 / 2400 | 0 | 2.0 | 16.7 | 1.14 (1.00-1.29) | 21 / 24 | 66, 1133 | 1.000 |
  | normalized Armijo | 0 / 2400 | not guaranteed; 0 found | 1.0 | 14.4 | 1 | 0 / 24 | none | 32.3 |

  The normalized Armijo start takes very small steps on these perturbed knots: after 100
  iterations its energy is still 6-100 times `E*` (median 32), so none of the 24 reaches 1.1 E*;
  to a common target (`max(1.1 E*, Armijo's final energy)`) both guarded variants need 0.12 of its
  time (median 150 s against 1440 s, 8-10 iterations against 100). The per-step cost is dominated
  by the metric (Gram + dense factorization), not by trials, so 6 trials against 2 cost 1.5×
  against 1.14× per step. Decision: **grow-previous (factor 2) is the default**, because it is
  within 1.5× of Armijo's time per step and is far faster to any energy level; the paper's start
  stays available as `StartRule::TwoThirdsTauMax` (its median per-step ratio 1.54 is just outside
  the D11 limit). The isotopy guarantee is the same for both (every trial is bounded by `τ_max` and
  CCD-checked); the audit (`path_contact`, gap 1e-6, 1e5 advances) found no violation in any of the
  72 runs. (WP3's first measurement on the Knot128 originals, which are at rest for `(3, 6)`, is
  superseded; per §17c it was not a valid set.)
* **Barnes–Hut in the session.** `BarnesHutEnergy`: `prepare` builds the tree at the base curve
  (each step, so also after an accepted step and after a subdivision); trial energies refit it
  (D26); `energy_and_differential` and `energy` also rebuild when the network's vertex or edge
  count differs from the tree's (a live subdivision asks for the energy before the next
  `prepare`; fixed, regression test `barnes_hut_flow_survives_subdivision`). Both sides of every
  Armijo comparison are Barnes–Hut energies (§17d).
* **Barnes–Hut stall near the minimum (found by WP6).** Freedman 600 with Barnes–Hut reached the
  near-circle (step 659) and then took steps of `τ ≈ 1e-14` with `‖g‖ ≈ 4e-3` until the cap.
  Diagnosis (`bh_differential_consistency` in `tests/flow.rs`, directional derivatives along
  random `v` at the final curve): (i) the *evaluated* Barnes–Hut energy is not smooth under the
  refit, because the admissibility tests flip as the curve moves and every flip changes the energy
  by about the approximation error: central differences of `Ê_bh` at `h = 1e-5` were 0.12 off
  against derivatives of order 1e-2 (the differential and the exact differential agreed to 1e-3),
  so Armijo compared noise; (ii) even without flips the Barnes–Hut differential is not the
  derivative of `Ê_bh` (§17d): FD of `Ê_bh` against `⟨dÊ_bh, v⟩` differs by 6e-4 to 3.5e-3, the
  same size as `⟨dÊ_bh − dÊ_exact, v⟩`. Fix applied: the trial energies of a step use the
  interaction lists (lumped nodes and exact leaf edges per edge, `BhTree::interaction_lists`) of
  the step's base curve with refitted node data (`BhTree::energy_with_lists`), so `Ê_bh` is
  smooth within the line search; this removed (i): the run now ends `Stuck` (five failed steps, as
  the reference does) at step 937 with `‖g‖ = 8e-4` after 355 s instead of spinning (the exact
  flow reaches 1e-4 at step 603). The remaining floor is (ii): the principled complete fix is the
  discrete adjoint of `Ê_bh` on the frozen lists (row `I` lumped term `ℓ_I M_N k(x_I, x̄_N, T_I)`
  plus the pushed-down derivative with respect to the node's `x̄_N` and `M_N` through the edges of
  `N`, and the exact pairs with only the `(I, J)` term), which makes the Barnes–Hut flow a descent
  method on the function it evaluates; not done here, left to a separate package. Until then a
  Barnes–Hut flow ends `Stuck` at `‖g‖` ≈ 1e-3, near-circle and energy within 0.4 % of the exact
  flow's (6.1015 against 6.0816 at the same curve, the Barnes–Hut energy 0.33 % high).
* **WP3 API differences from §9.1.** `PenaltyTerm` is a trait (`evaluate`, `subdivided`,
  `obstacle`) and `FlowOptions.penalties: Vec<Arc<dyn PenaltyTerm>>` replaces the `Penalty` enum
  there; `FlowOptions.obstacles` (CCD obstacles, `CollisionObstacle`), `dense_threads` (the dense
  LU pool, default `min(16, threads)`, `DEFAULT_DENSE_THREADS`, from the scaling measured in WP3:
  16 threads is the knee) and `record_path` (each `StepReport` carries the piecewise-linear path
  for audits) are added. `Status` has `Running` and `Error(RcError)`; `Solver` has `Custom` (the
  hook for other solvers; `Solver::Multigrid` since WP3c); `Backend::Auto` is exact +
  dense at every size until WP3c (§17h). The CCD Lipschitz constant is taken on velocities relative to the
  mean velocity of the motion (`max|u − ū| + max|w − ū|`, still the Lean lemma's constant for the
  shifted motion), which is sharper for a nearly rigid move. `StepStatus::RetargetStalled`: the
  length re-target stalled and the descent step did not move; five in a row end the run `Stuck`.
  L-BFGS uses the scaling `H₀ = γ M⁻¹` with `γ = sᵀy / yᵀM⁻¹y` (one extra solve) and caps trials
  at `τ₀ ≤ 1` once curvature pairs exist.
* **WP7.** The surface-potential `θ` default is 0.1 (measured error 6.9e-2 at 0.1, 1.5e-1 at
  0.25); a vertex on an obstacle gives +∞ (a gradient request is an error); Barycenter + Surface
  on a plane is redundant (translation along the normal is fixed twice), so planar-surface scenes
  pin a Point instead; implicit surfaces include the reference's double torus and smooth union
  plus `RoundedBox` and `Metaballs`.

## 17e. Noted after WP5 (orchestrator)

* **BCT default ε = 0.25, not 0.5.** Measured operator error of `Ā` against dense (random
  vector): at 0.5, 2.4e-2 to 5.9e-2 on Freedman 600-2400 and a perturbed Knot128 knot (projected
  gradient 3 % off dense on Freedman 600), above the 1e-2 target of §13.4; at 0.25, 8e-4 to
  6.2e-3 (gradient 0.15 % off). Default set to 0.25 by the reliability rule; §14.3 measures the
  speed cost against multigrid convergence.
* The reference admits a block when `2 r_A / (d − r_B) < 1` (position boxes only), i.e. about
  `r/d < 1/3` for equal clusters, with no tangent criterion; ours is Eq. 28 with both criteria.
* `BctOptions` carries `leaf_size` (BVH leaf, default 8); θ = 0 means no block is admissible;
  κ⁰ is computed as `|T̂ × d|²` (never negative by rounding).
* 0th-order blocks are coarse entry by entry (worst entry per block off by 24 % on average at
  θ = 0.1); the operator error is far smaller because entry errors average out.

## 17d. Noted after WP4 (orchestrator)

* **D15 amended.** The split rule as printed (minimize the sum of squared 6-D diameters of the
  two children over all positions) peels single tangent-points off the root on tangled curves
  (Freedman 1200: depth 1001, mean leaf 1.2 edges, quadratic build, no speed-up). Split positions
  are restricted to leave at least 25 % of the elements on each side
  (`bvh::DEFAULT_MIN_SPLIT_FRACTION`; 0 restores the printed rule): depth 13, build 2 ms for
  2368 edges.
* **Leaf lumping** (`BhTree::set_lump_leaves`, off): lumping admissible leaves too is about 3x
  faster on Freedman 2400 (one thread) with differential error 1e-2 instead of 1.9e-3. Decided
  in §14.3 by measurement (convergence and iterations, not only error).
* **Option names:** `theta_x = 0.125` is the reference's `diameter / distance < 0.25` written as a
  radius; `theta_t = 0.25`.
* The §8.2 differential is the derivative of the lumped `(I, J) + (J, I)` terms, not of the
  evaluated Barnes–Hut energy (which lumps `(I, J)` only): 4 % apart at ε = 1, 1e-3 at 0.5.
  Documented and tested; the line search compares Barnes–Hut energies with each other.

## 17c. The authors' program (built 2026-10-26)

Built outside the repository: `/var/data/discrete-geometry/RepulsiveCurves-build/` (working copy
`src/` with submodules at the pinned commits, binary `build/bin/rcurves_app`, g++-13 + OpenMP;
`rcurves_shared` does not compile and is not needed). A batch mode was added as one patch
(`batch-mode.patch`, no change to energy, solver or line search):
`rcurves_app --batch --steps N --out DIR [--threads T] [--obj-every K | --no-obj] [--no-sobolev]
[--no-backproj] [--no-bh] [--multigrid] [--no-exact-energy] [--ignore-stop] scene.txt`, writing
`steps.csv` (energy exact and Barnes–Hut, gradient norms, step size, per-phase ms), `objs/` and
`summary.txt`. Ready scenes in `RepulsiveCurves-build/scenes/`, runs in `RepulsiveCurves-build/runs/`.

* `libgmultigrid` is MIT too (Copyright (c) 2020 Christopher Yu): it may be read and ported; add
  its notice to `LICENSES/` when anything is ported.
* **Datasets:** the Knot128 `original` curves are already at rest for `(3, 6)` (knot0001: gradient
  norm 6e-7 at the start). The evaluation uses the `perturbed` set (about 2000 vertices) and
  Trefoil100.
* Reference defaults (file:line in the working copy): α, β = 3, 6 (`lws_options.cpp:15`; Freedman
  scenes use 2, 4.5); Barnes–Hut opening θ = 0.25 (`spatial/tpe_bvh.cpp:168`); BCT separation 1.0
  with dense blocks for s₁ + s₂ ≤ 8 (`product/block_cluster_tree.cpp:100`); multigrid tolerance 1e-2
  (relative ∞-norm), at most 20 V-cycles, `max(1, log2|V| − 4)` levels, CG smoother 12 iterations
  per level (24 for the initial guess), direct coarsest solve (`libgmultigrid/multigrid_hierarchy.h`);
  Armijo σ = 0.01 with the Sobolev–L² cosine, start `1/|g|` (or `1/√|g|`) or `min(1.5 τ_prev, 4·start)`,
  halving only, fail below 1e-15 (`tpe_flow_sc.cpp:176-196`); back-projection 1e-4 absolute, up to
  10 halvings × 3 Newton steps; stop when the Sobolev dot < 1e-4 (absolute), or 5 consecutive
  failed steps; no constraints in a scene means Barycenter + EdgeLengths; threads `max/2 + 2`.
* Measured alone (34 threads, s/step): Freedman 600 BH+multigrid 0.94, BH+dense 0.40, exact dense
  1.29; 1200: 2.07, 1.62, 5.21; 2400: 4.19, 8.0, 21.8. Freedman 600 (BH+dense) converges in 505
  steps (203 s), energy 2276 to 6.03, with a long plateau near 1510 for about 150 steps.

## 17g. Multigrid on networks (WP6b)

**Symptom.** The authors' CapillaryNetwork scene (two interleaved networks, junctions of degree
3–5, 297 pinned endpoints, all edge lengths fixed, `(α, β) = (3, 6)`): on the 2146-vertex curve
the cycles stagnated at a relative residual of 0.2–0.4 (20 cycles, returned without error) with a
projected gradient 98 % off the dense one; on the 4296-vertex curve they diverged (residuals 5–13).
The authors' program crashes on this scene with `--multigrid` ("Couldn't orient edge pair"), so
there is nothing to compare with.

**Cause: near-contacts, not junctions.** The scene has 262 (2146 vertices) and 142 (4296) edge
pairs whose midpoints are closer than half an edge length (branches of the two networks that
nearly cross; worst 0.095 ℓ and 0.049 ℓ). Their weights `ℓ_I ℓ_J |x_I − x_J|^{−(2σ+1)}` and
`k²_{2σ+5}` (which grows like `d^{−(2σ+3)}` for a crossing) dwarf their neighbours'. The
rediscretized coarse metric (§8.4) does not reproduce them: the coarse edges' midpoints are
elsewhere. Measured with a smooth coarse vector `x` (an approximate coarse solution):
`(Πx)ᵀ H_ℓ (Πx) / xᵀ H_{ℓ+1} x = 53, 5.5, 3.4` on levels 0→1, 1→2, 2→3 of the 4296-vertex curve,
so the coarse correction overshoots by up to 53×. With the exact Galerkin operators `Πᵀ H Π`
(matrix-free, `tests/multigrid_networks.rs`, `probe_galerkin`) four levels reach 1e-3 in about four
cycles: the coarse spaces (prolongation, coarse constraints, projectors) are adequate, the coarse
operator is not. Ruled out by measurement: the block cluster tree (θ = 0 behaves the same);
coarse `EdgeLengths` rows (dropping them on coarse edges, or at junctions only, does not help);
a Galerkin correction of the pairs of coarse edges that meet at a junction (no effect); the
projectors between levels (`|J x − P J x| / |J x|` ≈ 1.5e-2, energy ratios with and without
`P` equal to 1 %).

**Fix** (`src/multigrid/near.rs`, `MgOptions`):

1. *Galerkin near-pair correction* (`near_field_factor`, default 1): every coarse edge pair that
   contains a pair of the finer level with `|x_I − x_J| < 1 · max(ℓ_I, ℓ_J)`, or a pair corrected
   there, uses `S = Σ Jᵀ E J` over the finer pairs inside it instead of its rediscretized term
   `R` (pair terms are 4 × 4 scalar blocks of Eq. 31; the level metric is `Ā_ℓ + Σ (S − R) ⊗ I`,
   added to the dense matrix on the coarsest level). Checked: the blocks sum to the dense metric
   (1e-12), and with every pair corrected the coarse metric is `Jᵀ A J` exactly (1e-10). The
   smooth-vector ratios above become 1.03, 1.03, 1.13.
2. *Block-preconditioned smoother* (`smoother_factor`, default 1): on levels with near pairs or
   corrections, the CG smoother is preconditioned by `D_rest + Σ E(I, J)` over the near pairs
   (`D_rest` the rest of the metric's diagonal from the BCT row sums; sparse Cholesky, `⊗ I_d`,
   between projections). Near-contacts make the metric stiff (synthetic network, `κ(P Ā P)` on
   `ker C`: 4.7e4 without crossings, 2.5e7 with crossings at 0.1 h, 4.4e8 at 0.05 h; only
   `λ_max` grows); the block solve removes those couplings. Levels without near pairs keep plain
   CG (closed curves unchanged).
3. *Stop at fixed-vertex-dominated levels* (`max_fixed_fraction`, default 0.4): coarsening stops
   when more than 40 % of a level's vertices cannot be removed (degree ≠ 2, pins), if the level has
   at most 8 · `coarsest_vertices` vertices. On the scene, levels below 1073 vertices shrink by
   3–20 % per level (547 junctions and endpoints), cost a product each and raised the cycle count
   from 6 (two levels) to 15 (seven levels).
4. *Energy-norm check* (`energy_check`, default on): flexible CG also requires
   `rᵀ z ≤ tol² bᵀ x` (`z` the V-cycle of `r`) before it stops. The Euclidean residual relative
   to `‖P b‖` is dominated by the near-contact entries of `b`; on a synthetic network it reached
   8e-4 with a direction 20 % off the dense one (6.7e-3 with the check, one extra cycle at most).
5. *No silent failure* (`failure_residual`, default 1e-2): a solve that ends unconverged above it
   is `RcError::MultigridDiverged` (the message says "diverged" also for stagnation).

**Results** (default options, BCT θ = 0.25, 32 threads on a shared machine, so times are
indicative). 2146 vertices: levels 2146 → 1071, 1 cycle, residual 5.4e-4, gradient 2.8e-3 off
dense (dense 22 s; multigrid build 5 s + solve 2 s); 20 steps of the exact-energy flow (scene
obstacle included) with multigrid and with dense: energies agree to 1e-6–3e-5 relative, positions
to 2.6e-4 edge lengths, 1–2 cycles per step. 4296 vertices: levels 4296 → 2146 → 1073, 9 cycles,
residual 5.4e-4, gradient 3.3e-3 off dense (dense 80 s; build 9 s + solve 17 s);
3 flow steps: energies agree to 5e-7–5e-6, positions to 3.5e-4 edge lengths, 9–11 cycles
per step. Freedman 600–4736 (closed): same or fewer cycles than before (2, 3, 2 → 1, 3, 2), same
error against dense. A synthetic network with junctions of degree 3, 4, 5, pinned endpoints and
branches crossing at 0.1 h (`tests/multigrid_networks.rs`): 1–9 cycles, error 1.7e-3 at the
default tolerance and 2e-8 at 1e-8.

**Open.** (a) The product of the block cluster tree is slow on this scene (72 ms per product at
4296 vertices), which dominates the multigrid time; WP5/WP14. (b) Crossings at 0.05 h exactly
over edge midpoints (synthetic) still do not converge in 20 cycles; this is now reported as
`MultigridDiverged`, and the session should fall back to the dense solver (WP3b: `src/step`).
(c) The multigrid `tolerance` keeps its meaning (Euclidean residual); the direction error at
1e-3 is a few 1e-3 on stiff networks.

## 17h. Multigrid in the session, fallback, `Backend::Auto` (WP3c)

**Plumbing.** `SolverRequest` carries the session's two pools (`pool`: Gram assembly, block cluster
tree, hierarchy; `dense_pool`: the `faer` LU). The factories hold no pools any more
(`DirectSolvers::new()`, `MultigridSolvers::new(options)`); a `Custom` factory builds on the pools
of its request, so `FlowOptions::threads` / `dense_threads` govern every parallel loop of the
solvers. (`SurfacePotential::evaluate` still uses rayon's global pool for a per-vertex map; it is
elementwise, hence deterministic.) Tested bitwise at 1, 3 and 8 threads for exact/BH × dense/multigrid
(`tests/session_backends.rs`).

**`Solver::Multigrid(MgOptions)`** is `MultigridSolvers` on the session's pools, with both energy
models. The hierarchy is part of the saddle solver, which is built at every step and every
re-target, so a subdivision (new network, new constraint set) rebuilds it; there is no state to
invalidate. Constraints: the hard constraints of the step are coarsened (§8.4); soft problems use
the mean-free rows. Tests: trefoil with live subdivision (60 → 120 vertices, levels 2 → 3, no
fallback, both energies), an open curve with both ends pinned, a fixed end tangent and all edge
lengths. The multigrid projected gradient differs from the dense one by the block cluster tree's
approximation of the metric (θ = 0.5): 0.9 % in the gradient norm on the open test curve at any
multigrid tolerance (1e-3 … 1e-10), 3e-3 on closed curves (§17g); later steps of two such flows
differ because the line search is path dependent.

**Failure classes.** `MultigridDiverged { residuals }`: the final residual is larger than the
residual after the full-multigrid start (or not finite, or above `‖b‖`): the cycles are not a
contraction. `MultigridNotConverged { residuals, failure_residual }`: it did not grow, but the
budget ended above `failure_residual` (stagnation; before WP3c this was also called "diverged").

**Fallback.** `FallbackSaddle` (`src/step/fallback.rs`) wraps the solver of every step and
re-target of a non-dense session. On either error, and when `|V| ≤
FlowOptions::dense_fallback_max_vertices` (default 5000, 0 = never), it builds the dense saddle
solver of the same step (same metric at the base curve, same constraint rows, same pools),
repeats the failed solve, and answers all later solves of the step with it (the line search's
chord iterations use the same metric). The next step tries multigrid again. `StepReport::
solver_fallback = Some(SolverFallback { from, to: "direct", cause, dense_solves, build_seconds })`
records it; `linear.solver` of that solve names the dense solver. Above the size limit, and for
every other error, the error goes up (`Status::Error`). A fallback step is bitwise the dense
session's step (test). Not exercised on a real network by a natural failure: the synthetic
crossings at 0.05 h of §17g(b) are the intended case.

**Crossover measurement** (`tests/session_backends.rs`, `crossover_one_configuration`, ignored;
logs `runs_tmp/rc-wp3c-x/`). Freedman–He–Wang, `(α, β) = (2, 4.5)`, default options, 8 threads
for the main and the dense pool, 5–8 steps from the start curve, mean seconds per step. The
machine was shared (load 30–40 on 64 cores, 4–6 of these jobs at once), so repeated runs differ
by up to 40 %; **indicative only, WP14 re-measures alone.** 300 = every other vertex of the 600
file; 1600 and 2000 = the 1200 curve with evenly spread edges split at their midpoints (same
shape). 1600 and 2000 were run for exact + dense, exact + MG and BH + MG only.

| vertices | exact + dense | BH + dense | exact + MG | BH + MG |
|---:|---:|---:|---:|---:|
| 300 | 0.13 | 0.13 | 0.19 | 0.19 |
| 600 | 1.0 (0.3 alone) | 0.6 | 1.7 | 1.3 |
| 1200 | 3.6, 2.8, 2.5 | 3.4, 2.7, 2.5 | 6.7, 4.6, 4.4 | 3.7, 4.6, 4.2 |
| 1600 | 6.8 | | 11.3 | 11.2 |
| 2000 | 11.7 | | 11.8 | 11.8 |
| 2400 | 16.9 | 17.2 | 12.4 | 11.7 |

Where the time goes at 2400 (exact + dense): metric/LU 13.3 s, line search 2.6 s, solve 0.6 s,
energy 0.25 s; BH + MG: line search 10 s (the chord iterations of the projections, one multigrid
solve each), metric 0.26 s, solve 1.3 s, energy 0.2 s. At 1200 the dense LU (1.8–2.4 s) is still
cheaper than the multigrid projections (3.4–5 s). Exponents: dense ≈ n^2.3–2.7, BH + MG ≈ n^1.5.
Barnes–Hut with the dense solver saves nothing (the energy is 3–8 % of the step); the
multigrid with the exact energy equals BH + MG. The curves cross at about 2000 vertices.

**Rule (superseded by §17k: 1500, measured alone).** `Backend::Auto { dense_max_vertices: 2000 }`: exact + dense for `|V| ≤ 2000`, BH
(`BhOptions::default()`) + multigrid (`MgOptions::default()`) above. A tie goes to the exact
computation. The rule is applied to the vertex count at session start and after each subdivision
(`StepReport::{energy_backend, solver_backend}` name what ran). The multigrid line search is the
obvious thing to speed up (its projections rebuild nothing but cost a multigrid solve per chord
iteration); a faster projection would move the crossover down.

**`run` and errors.** `Session::step` returns the error of a failing step and sets `Status::
Error`; `Session::run` ends with `Status::Error(e)` (`FlowReport::error()`), keeps the steps
before and the curve of the last completed step, and returns `Ok`. Documented in both and tested.

**Length schedule in the report.** `StepReport::schedule = Some(ScheduleProgress { progress,
span, fraction, factor_target, factor_imposed, factor_measured, finished })`: the final factor
asked for, the factor the length constraints impose after the step, and the total length of the
curve over its total length when the schedule was set.

## 17i. Barnes–Hut adjoint (WP3d)

**Problem** (§17d, §17f). The trial energies of a step are evaluated on the base curve's
interaction lists, so within a line search the Barnes–Hut energy is a smooth function; but the
§8.2 differential (the derivative of the lumped `(I, J) + (J, I)` pair with `L_N, x̄_N, T̂_N`
fixed) is not its derivative: finite differences differ by 6e-4 to 3.5e-3, and Freedman 600 ended
`Stuck` at `‖g‖ ≈ 8e-4`. Fix: differentiate the evaluated energy exactly (the discrete adjoint).

**The function.** With `𝓛(I)` the lumped nodes and `𝓧(I)` the exactly evaluated edges of row `I`
(`BhTree::interaction_lists`, frozen during a step) [formula: rc-bh-adjoint]

```
Ê_bh(γ) = Σ_I [ Σ_{N ∈ 𝓛(I)} ℓ_I L_N k(x_I, x̄_N, T_I) + Σ_{J ∈ 𝓧(I)} e_IJ ],
L_N = Σ_{J∈N} ℓ_J,   x̄_N = Σ_{J∈N} ℓ_J x_J / L_N.
```

`T̂_N`, the radii and the boxes only enter the admissibility test, which the frozen lists replace;
so the node data that move with the curve are `L_N` and `x̄_N`.

**One lumped term** `f = ℓ_I L_N k(x_I, x̄_N, T_I)` (`k, ∇_p k, g_T = ∇_T k` at `(x_I, x̄_N, T_I)`,
`∇_q k = −∇_p k`, `s₁ = −1, s₂ = +1`):

```
δf = L_N [ (k T_I + g_T)·(δγ_{i₂} − δγ_{i₁}) + ½ ℓ_I ∇_p k·(δγ_{i₁} + δγ_{i₂}) ]        (own part)
   + ℓ_I k δL_N + ℓ_I L_N ∇_q k·δx̄_N                                                   (node part)
δL_N = Σ_{J∈N} δℓ_J,   δx̄_N = L_N⁻¹ Σ_{J∈N} (δℓ_J (x_J − x̄_N) + ℓ_J δx_J).
```

With `b = ℓ_I k` and `u = ℓ_I ∇_q k = −ℓ_I ∇_p k` the node part is
`Σ_{J∈N} [(b + u·(x_J − x̄_N)) δℓ_J + ℓ_J u·δx_J]`, and with `δℓ_J = T_J·(δγ_{j₂} − δγ_{j₁})`,
`δx_J = ½(δγ_{j₁} + δγ_{j₂})`:

```
∂f/∂γ_{i_a} = L_N (s_a (k T_I + g_T) + ½ ℓ_I ∇_p k),
∂f/∂γ_{j_b} = s_b λ_J T_J + ½ ℓ_J u,    λ_J = b + u·(x_J − x̄_N)      (J ∈ N).
```

The own part is the §8.2 increment without the `(J, I)` term; `L_N` cancels from the node part
(`u` carries no `L_N`).

**All lumped terms at once.** `λ_J` and `u` are linear in `(b, u)`, so per node
`b_N = Σ_{I : N ∈ 𝓛(I)} ℓ_I k`, `u_N = Σ_I (−ℓ_I ∇_p k)` suffice, and an edge `J` collects
`Σ_{N ∋ J} (b_N + u_N·(x_J − x̄_N))` and `ℓ_J Σ_{N ∋ J} u_N`. The nodes containing `J` are the
path from the root to `J`'s leaf. Downward pass [formula: rc-bh-adjoint-shift]:

```
β_root = b_root, U_root = u_root;    child C of P:  β_C = β_P + U_P·(x̄_C − x̄_P) + b_C,  U_C = U_P + u_C
⇒ Σ_{A on the path to Λ} (b_A + u_A·(x − x̄_A)) = β_Λ + U_Λ·(x − x̄_Λ)
⇒ ∂/∂γ_{j_b} (all lumped terms) = s_b λ_J T_J + ½ ℓ_J U_Λ,   λ_J = β_Λ + U_Λ·(x_J − x̄_Λ).
```

Re-centring at every node avoids the cancellation in `Σb − Σu·x̄ + (Σu)·x` for a curve far from
the origin (checked: `|Σ_v g_v|/‖g‖ = 2e-16` at an offset of `1e4`).

**Exact pairs.** `∂e_IJ/∂γ_I` in row `I`, `∂e_IJ/∂γ_J` (§3.2) for the partner. The lists are not
symmetric (`J ∈ 𝓧(I)` while `I`'s leaf is inside a node that `J` lumps), so the partner part is
the transposed relation. Since a row's exact edges are all non-adjacent edges of the leaves its
traversal opened, the transpose is taken over opened leaves (6–8× fewer entries than pairs), and a
parallel pass over the leaves adds, for each row `I` that opened leaf `Λ` (increasing `I`),
`∂e_IJ/∂γ_J` for every `J ∈ Λ` not adjacent to `I`.

**Lean** (`DdgProofs/RepulsiveCurves/BarnesHutAdjoint.lean`): `hasDerivAt_nodeMass`,
`hasDerivAt_nodeMoment`, `hasDerivAt_nodeCentre` (the derivative of `x̄_N`), `bh_adjoint` (the
whole chain rule of one lumped term along a motion of the curve: own part plus node part,
`ν > 0`, `L_N ≠ 0`), the algebraic identities `adj_own_alg`, `adj_node_alg`, and
`bh_adjoint_shift`, `bh_adjoint_shift_root` (the downward pass as a fold over the path) — ids
`rc-bh-adjoint`, `rc-bh-adjoint-shift`, both proved. The exact pairs are
`hasDerivAt_pairEnergy` (`rc-energy-differential`); summing over the lists is the sum rule (not
restated). Nothing is cited. (§12 and §8.2 should list the two ids and point to this section;
not edited here.)

**Implementation** (`src/energy/barnes_hut.rs`). `BhTree::energy_and_differential` is now the
adjoint on the lists of the current tree and positions (it computes them; same signature, so
`step::BarnesHutEnergy` uses it unchanged); `energy_and_differential_with_lists(lists, pool,
exponents)` takes the lists (error `InvalidOption { name: "lists" }` for lists of another tree);
the WP4 differential stays as `energy_and_differential_lumped_pairs`. The energy returned is
`energy_with_lists` on the same lists bit for bit, and `BhTree::energy` (traversal) now sums the
lumped terms before the exact pairs, so it is too: the base energy of a step and its trial
energies are the same function, summed in the same order. `classify` / `interaction_lists`
follow `set_lump_leaves`. Order of work: rows in parallel in fixed chunks of the tree order
(at most 64 chunks, at least 16 rows each, at most 2²² chunk-local node slots: a function of the
sizes only) with chunk-local `b, u`; node sums over the chunks in chunk order (parallel over
nodes); downward pass (sequential, `O(#nodes)`); transpose of the opened leaves (sequential,
`O(#opened leaves)`); partner pass (parallel over leaves); vertex sums in edge order. Bitwise
independent of the thread count (tested at 1, 3, 8 threads, also on refitted curves with frozen
lists). Kernel evaluations: one per lumped entry (WP4: two), eight per exact pair (as WP4).

**Finite differences** (`tests/barnes_hut_adjoint.rs`; central differences at `h` and `h/2`,
Richardson-extrapolated, `h = 10⁻³ ℓ̄`, `energy_with_lists` on the frozen lists; perturbed trefoil
240, open helix 300, a network of four meridian arcs between two degree-4 junctions plus a
branch at a degree-3 junction (174 vertices); `(3, 6)` and `(2, 4.5)`; `ε = θ_T = 2θ_x ∈ {0.25,
0.5}` and `θ_x = θ_T = 0.5`, leaf size 4). Worst `|FD − ⟨g, v⟩| / (‖g‖‖v‖)` over 4 random `v`, and
`|FD − ‖g‖²| / ‖g‖²` along `v = g`:

| | adjoint, random `v` | adjoint, along `g` | WP4 differential, random `v` | WP4, along `g` |
|---|---|---|---|---|
| trefoil 240 | 8e-12 – 1.4e-9 | 2e-11 – 4e-10 | 2e-5 – 2.8e-3 | 3e-6 – 2.3e-3 |
| helix 300 | 1e-12 – 2.4e-9 | 1e-11 – 3e-10 | 2.5e-6 – 4.4e-4 | 6e-6 – 5e-4 |
| network 174 | 2e-12 – 4.6e-10 | 5e-13 – 2e-12 | 1.5e-5 – 1.7e-3 | 4e-6 – 9e-4 |

(the `1e-9` entries are `(3, 6)`, the FD noise of the larger exponent). Coordinate by coordinate
on a 72-vertex network (`θ = 0.5`, leaf size 2, 891 lumped entries): `‖g − FD‖/‖FD‖ = 3e-11`
(WP4: 0.12). At a moved curve (perturbed by `0.3 ℓ̄`) with the base curve's lists, as in a line
search: 6.5e-11 / 1e-11. `θ = 1e-12` (nothing lumped): equal to the exact differential to
`1.4e-15`. Against the exact differential the adjoint is about as close as the WP4 one
(perturbed trefoil 400, `(3, 6)`: 8.4e-3, 1.5e-3, 2.4e-4 at `ε = 1, 0.5, 0.25` against 6.8e-3,
9.5e-4, 7.5e-5; the existing tests of `tests/barnes_hut.rs` pass unchanged).

**Flow.** Freedman–He–Wang 600, `(2, 4.5)`, Barnes–Hut (default options) with the dense solver,
default flow options otherwise, 8 threads (`freedman_600_barnes_hut_flow_converges`, ignored,
release): **`Converged`** at step 711 with `‖g̃‖ = 9.86e-5` (before: `Stuck` at step 937,
`‖g‖ = 8e-4`); Barnes–Hut energy 6.10156, exact energy of the final curve 6.08164 (the exact flow:
6.0816 at step 603); 314–343 s on the shared machine (two runs, load 30–40).

**Cost** (`timing_old_and_new_differential`, ignored; default options, 8 threads, minimum of 5,
load 30–40 on 64 cores, indicative):

| vertices | WP4 differential | adjoint (computes the lists) | adjoint, lists given | lists | energy on lists | lumped entries | exact pairs |
|---|---|---|---|---|---|---|---|
| 592 (Freedman 600) | 17.4 ms | 21.1 ms | 20.5 ms | 0.5 ms | 6.3 ms | 313 | 345 063 |
| 2368 (Freedman 2400) | 78.3 ms | 76.1 ms | 77.1 ms | 2.8 ms | 29.3 ms | 276 015 | 1 852 338 |
| 9472 (2400 subdivided twice) | 160.3 ms | 172.3 ms | 144.1 ms | 12.7 ms | 65.3 ms | 2 459 754 | 3 066 039 |

A first version that transposed the exact pairs and scattered per-entry node terms sequentially
cost 1.5× the WP4 differential at 9472 vertices (transpose 35–70 ms, scatter 35–65 ms).

**Step.** `BarnesHutEnergy::energy_and_differential` (`src/step/mod.rs`) computes the lists and
then calls `tree.energy_and_differential`, which computes them again (about 3 ms at 2368
vertices, 13 ms at 9472); passing them (`energy_and_differential_with_lists(&lists, …)`) saves
that. Nothing else changes.

## 17j. End-edge collapse and the guard (WP3c)

**Symptom.** The authors' `TangentInterpolation` scene (curve of 132 vertices, ends and two
interior vertices pinned, tangents fixed at both ends, **no** length constraint, a total-length
penalty of weight 0.5, `(α, β) = (3, 6)`, `iteration_limit 65`): in our session an end edge shrinks
by a factor of about 0.6 per step (step 3: 2.7e-2, step 10: 6e-4, step 20: 1.6e-5, step 30: 1e-7;
mean edge 5e-2) and at step 33 the step fails with `SingularMetric`. (The WP13 note said 60 steps;
with the scene's own settings it is 33, and both end edges collapse, the second one later.)

**The reference does not do this.** Run through the authors' program on the same scene
(`rcurves_app --batch --steps 100 --ignore-stop`, private copy, `runs_tmp/rc-wp3c-auth/`): the end
edge dips to 1.4e-2 at step 1 and then grows; the curve is at rest from step 22 (energy 36.65,
minimum edge 3.7e-2, last edge 9.6e-2, length 8.024; step sizes 1e-15 afterwards). The reference's
tangent constraint is the same as ours (the unit tangent of the end edge, three rows, Jacobian of
the normalized difference, `src/flow/constraint_functions.cpp`), so fixing the end tangent does
not by itself allow the collapse. Without the tangent rows our end edge still shrinks (to 1.4e-3)
but stalls there (exploration, `tests/end_edge.rs` history); with them it collapses, because the
constraint is then satisfied by a vanishing edge whose direction is fixed while the next edge is
free to turn: the tangent condition is circumvented, and the discrete energy has no term that
resists it (an edge's pair terms are weighted by its length and vanish with it).

**Cause (hypothesis, not proved).** The reference's energy is the vertex-pair sum
`Σ_{i,j} K(x_i, x_j, T_i) l_i l_j` with dual lengths `l_i` (`tpe_energy_sc.cpp`, `tpe_total`);
ours is the paper's edge-pair energy (D1, decision Q1 of §18: the vertex-pair variant is not
implemented). In the reference an end vertex and its neighbour keep nonzero weight and a vertex
tangent as the end edge shrinks, and the measured behaviour is a local minimum with a long end
edge; in ours the lowest energy is reached by dropping the end edge. We did not port the
vertex-pair variant to test this (it is excluded by Q1), so there is no fix "in the
paper/reference" for the edge-pair energy.

**Guard** (`FlowOptions::min_edge_length_ratio`, default `Some(1e-3)`, `None` disables; not in the
paper or the reference; it changes no energy and no step that is taken). After the line search of a
step and before the positions are taken, if some edge is shorter than `ratio · ℓ̄₀` (`ℓ̄₀` the mean
edge length at the start of the session) **and** shorter than it was before the step, the step
fails with `RcError::EdgeCollapse { edge, length, length_before, threshold, step }`; the curve is
that of the last accepted step and the status is `Error`. On the scene this happens at step 15–20
instead of `SingularMetric` at step 33, and the message names the edge and the remedies (fix the
edge lengths with an `EdgeLengths` constraint, or disable the guard). The check is
`O(|E|)` per step. Tests: `tests/end_edge.rs` (guard fires on the scene, curve unchanged by the
rejected step, `None` gives the old `SingularMetric`, quiet on a trefoil, option validation);
Python `min_edge_length_ratio`. A scene that needs the full 65 steps must fix its edge lengths or
set the option to `None` (then it fails with `SingularMetric` near step 33).

**Open.** Whether the vertex-pair energy of the reference removes the collapse (a differential
test against the authors' program on this scene, WP15, would show it).

## 17k. Benchmarks (WP14, 2026-10-26)

*(After WP14b — Schur dense solve, Newton polish, early abort, position-only multigrid trees, capped multigrid pool — the gate passes everywhere, the knot included, and `Backend::Auto` is 2500: §17l.)*

Benches: `benches/fig20.rs` (one flow configuration per call, one CSV row), `benches/threads.rs`
(every component at 1…64 threads with efficiency and a bitwise check), `benches/components.rs` (one
component, any thread list), shared helpers in `benches/common/`. Plain `harness = false` mains:
`cargo bench --bench fig20 -- --curve fhw|fhw-r|knot --size N --config exact|bh|mg|mg-noproj
--threads T [--dense-threads D --steps 20 --min-steps 5 --cap-s 240]`. Driver scripts, CSVs and
logs: `runs_tmp/rc-wp14/` (`jobs-*.txt`, `runjobs.sh`, `fig20-*.csv`, `threads.csv`, `ref/`).

**Protocol.** The machine was alone for all timing: the comparison jobs of the other agent were
gone at 16:14, the 1-minute load average was 3.96 at the start (16:15) and 1.1–1.5 during the
1-thread blocks. Blocks ran one after another; the `load_start/load_end` columns of the CSVs are
the 1-minute averages around each run. During 34- and 64-thread runs they include the benchmark's
own threads, and they lag at block boundaries (the first 1-thread runs after a 64-thread block
show 9–22 for that reason); spot checks with `ps` during the blocks showed no other process above
12 % CPU. Each
configuration is a *running flow from the start curve*, not a converged one: the mean over the
first 20 steps (the paper: up to 500), stopping earlier when the next step would pass a budget
of 240 s (at least 5 steps, down to 2–3 where one step takes 100–200 s; the `steps` column says
how many). The medians are in the CSVs; they agree with the means within 10 % for the dense configurations
and most multigrid ones, not where single steps vary a lot (the barycenter-only runs, knot 16384). Exponents (2, 4.5) on
Freedman–He–Wang, (3, 6) on the knot; hard constraints barycenter + edge lengths, projection
tolerance `MeanEdgeLength(1e-4)`, 3 chord iterations; multigrid and Barnes–Hut with their
defaults; `dense_fallback_max_vertices = 0` (no silent switch of solver). Main pool `T` threads,
dense pool `min(16, T)`. Sizes: the Freedman files are named 600/1200/2400 but hold 592/1184/2368 vertices (and as many
edges); they are refined by midpoint splitting to 4736/9472/18944 (nominal 4800/9600/19200). The
tables give the actual vertex counts. To compare with the paper's *sizes* the 2368 curve was
also resampled by arclength to exactly 1000, 2000, 4000, 8000, 16000 (`fhw-r`). Knot 96
(`Knot128/perturbed/knot0096.obj`, 2106 vertices) is resampled by arclength to 1024…16384.
"Multigrid without projection": the library has no switch for the back-projection (not changed:
no algorithm work in WP14), so `mg-noproj` runs with the barycenter as the only constraint (no
length rows); that flow differs from the paper's (edge lengths drift, 2–7 cycles per solve) and
its projections are not free (see below). The better proxy is the derived `total − projection`.

### Fig. 20 against the paper, 1 thread (the gate)

Paper values read off Fig. 20 (§14.1). Ours: Freedman–He–Wang resampled to the paper's sizes,
(2, 4.5), 1 thread (main and dense pool), seconds per iteration.

| configuration | n | paper | ours | verdict |
|---|---:|---:|---:|---|
| Exact | 1000 | ≈ 6 | 1.64 | below |
| Exact | 4000 | ≈ 145 | 69.5 (5 steps) | below |
| Barnes–Hut + dense | 2000 | ≈ 10 | 9.01 | below |
| Barnes–Hut + dense | 4000 | ≈ 60 | 63.0 (5 steps) | **5 % above** (figure reading error is of that size) |
| Multigrid | 2000 | ≈ 6 | 3.57 | below |
| Multigrid | 8000 | ≈ 25 | 8.39 | below |
| Multigrid | 16000 | ≈ 57 | 14.8 | below |
| Multigrid w/o projection, derived (total − projection) | 2000 | ≈ 3 | 1.56 | below |
| same | 16000 | ≈ 38 | 4.76 | below |
| same, measured (barycenter only) | 2000 / 16000 | ≈ 3 / ≈ 38 | 3.13 / 11.8 | 4 % above / below |

**Gate verdict.** Passed at every size for Exact and Multigrid (Multigrid by 1.7× at 2000, 3× at 8000, 3.9×
at 16000), failed by 5 % for Barnes–Hut + dense at 4000: that configuration is the dense LU (60.6 of
63.0 s), which runs at about 19 GFlop/s on one core, and the Schur-complement solve of §17b
(about 6× fewer flops) would pass it; `Backend::Auto` never selects Barnes–Hut with the dense
solver. The knot is harder than the Freedman curve (`mg` per iteration on 1 thread: 2.85 s at 1024, 11.7 s at 2048, 41 s at 4096, 127 s at 8192,
216 s at 16384, the paper's multigrid curve being 6 s at 2000, 25 s at 8000, 57 s at 16000):
the paper does not say which curve Fig. 20 uses; the Freedman curve is the one the datasets'
README names for it, and we pass there. On the knot we are above the paper from 2048 vertices
(2× at 2048, 5× at 8192), because the line-search projections of the near-contact knot cost
more (below).

### Fig. 20 on the listed sizes, seconds per iteration (means; steps in parentheses where < 20)

Freedman–He–Wang (2, 4.5). `metric` = Gram assembly + factorization (dense) or hierarchy
build (multigrid); `proj` = the projections of the line search, `f` its objective evaluations.

| n | config | thr | total | energy | metric | solve | line search | proj | ccd | f | peak GB |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 592 | Exact | 1 | 0.451 | 0.089 | 0.250 | 0.006 | 0.070 | 0.031 | 0.003 | 0.036 | 0.1 |
| 1184 | Exact | 1 | 2.58 | 0.356 | 1.73 | 0.021 | 0.328 | 0.092 | 0.007 | 0.229 | 0.3 |
| 2368 | Exact | 1 | 16.1 | 1.42 | 12.9 | 0.070 | 1.22 | 0.309 | 0.016 | 0.900 | 0.9 |
| 4736 | Exact | 1 | 113 (3) | 5.71 | 100 | 0.256 | 4.86 | 1.77 | 0.037 | 3.05 | 3.6 |
| 592 | BH + dense | 1 | 0.463 | 0.097 | 0.251 | 0.006 | 0.072 | 0.030 | 0.003 | 0.039 | 0.1 |
| 1184 | BH + dense | 1 | 2.46 | 0.303 | 1.72 | 0.021 | 0.291 | 0.091 | 0.007 | 0.194 | 0.3 |
| 2368 | BH + dense | 1 | 14.2 (16) | 0.466 | 12.9 | 0.069 | 0.580 | 0.277 | 0.015 | 0.288 | 0.9 |
| 4736 | BH + dense | 1 | 104 (3) | 0.634 | 101 | 0.252 | 2.15 | 1.78 | 0.037 | 0.329 | 3.6 |
| 592 | Multigrid | 1 | 1.08 | 0.096 | 0.053 | 0.121 | 0.767 | 0.726 | 0.003 | 0.038 | 0.1 |
| 1184 | Multigrid | 1 | 2.60 | 0.302 | 0.091 | 0.347 | 1.74 | 1.54 | 0.006 | 0.192 | 0.1 |
| 2368 | Multigrid | 1 | 4.27 | 0.459 | 0.130 | 0.606 | 2.90 | 2.56 | 0.015 | 0.324 | 0.1 |
| 4736 | Multigrid | 1 | 6.51 | 0.566 | 0.179 | 0.951 | 4.60 | 4.14 | 0.043 | 0.409 | 0.1 |
| 9472 | Multigrid | 1 | 11.8 | 0.743 | 0.286 | 1.64 | 8.83 | 8.19 | 0.115 | 0.533 | 0.2 |
| 18944 | Multigrid | 1 | 27.7 (8) | 1.31 | 0.478 | 3.37 | 22.0 | 21.0 | 0.328 | 0.731 | 0.3 |
| 592 | Multigrid, barycenter only | 1 | 2.07 | 0.095 | 0.042 | 0.280 | 1.61 | 1.56 | 0.005 | 0.041 | 0.0 |
| 2368 | same | 1 | 3.63 | 0.450 | 0.092 | 0.606 | 2.30 | 1.91 | 0.029 | 0.363 | 0.1 |
| 4736 | same | 1 | 3.86 | 0.507 | 0.142 | 0.757 | 2.26 | 1.78 | 0.078 | 0.400 | 0.1 |
| 9472 | same | 1 | 7.26 | 0.730 | 0.228 | 1.37 | 4.66 | 3.86 | 0.236 | 0.560 | 0.2 |
| 18944 | same | 1 | 16.1 (14) | 1.36 | 0.397 | 2.31 | 11.6 | 9.73 | 0.776 | 1.04 | 0.3 |
| 592 | Exact | 64 | 0.154 | 0.003 | 0.113 | 0.005 | 0.031 | 0.026 | 0.004 | 0.001 | 0.1 |
| 1184 | Exact | 64 | 0.556 | 0.010 | 0.438 | 0.018 | 0.086 | 0.072 | 0.007 | 0.007 | 0.3 |
| 2368 | Exact | 64 | 2.98 | 0.039 | 2.63 | 0.057 | 0.233 | 0.197 | 0.011 | 0.025 | 1.1 |
| 4736 | Exact | 64 | 24.3 (9) | 0.146 | 22.7 | 0.195 | 1.15 | 1.03 | 0.027 | 0.095 | 3.7 |
| 9472 | Exact | 64 | 170 (2) | 0.576 | 163 | 0.727 | 5.54 | 5.27 | 0.052 | 0.211 | 14.2 |
| 2368 | BH + dense | 64 | 2.68 | 0.024 | 2.39 | 0.054 | 0.203 | 0.181 | 0.011 | 0.010 | 1.1 |
| 4736 | BH + dense | 64 | 21.5 (11) | 0.033 | 20.3 | 0.188 | 0.944 | 0.903 | 0.028 | 0.013 | 3.7 |
| 9472 | BH + dense | 64 | 193 (2) | 0.061 | 186 | 0.725 | 5.37 | 5.31 | 0.049 | 0.010 | 14.2 |
| 592 | Multigrid | 64 | 0.392 | 0.007 | 0.036 | 0.051 | 0.296 | 0.289 | 0.004 | 0.002 | 0.1 |
| 1184 | Multigrid | 64 | 0.814 | 0.019 | 0.061 | 0.132 | 0.597 | 0.579 | 0.008 | 0.010 | 0.1 |
| 2368 | Multigrid | 64 | 1.40 | 0.028 | 0.091 | 0.247 | 1.03 | 1.00 | 0.013 | 0.013 | 0.2 |
| 4736 | Multigrid | 64 | 2.39 | 0.038 | 0.129 | 0.412 | 1.80 | 1.76 | 0.026 | 0.015 | 0.3 |
| 9472 | Multigrid | 64 | 4.89 | 0.057 | 0.205 | 0.784 | 3.83 | 3.76 | 0.054 | 0.019 | 0.4 |
| 18944 | Multigrid | 64 | 12.3 (19) | 0.099 | 0.355 | 1.86 | 9.98 | 9.83 | 0.120 | 0.029 | 0.7 |
| 4736 | Multigrid, barycenter only | 64 | 1.39 | 0.037 | 0.103 | 0.366 | 0.881 | 0.837 | 0.028 | 0.016 | 0.2 |
| 18944 | same | 64 | 6.38 | 0.103 | 0.266 | 1.20 | 4.79 | 4.62 | 0.135 | 0.034 | 0.5 |

Knot 96 (3, 6), same columns, totals only (`mg`, 1 thread / 64 threads; dense `exact` = `bh`
within 3 %): 1024: exact 1.52 / 0.416, mg 2.85 / 0.727; 2048: exact 9.92 / 2.05, mg 11.7 / 2.35;
4096: exact 72.0 (4 steps) / 15.3, mg 41.0 (5) / 6.09; 8192: mg 127 (2) / 16.0; 16384: mg 216 (2)
/ 29.9 (7; the median is 14.7); full rows in `runs_tmp/rc-wp14/fig20-knot.csv`.

**Skipped.** Dense at 9472 on one thread (the LU alone would take about 13 min per step: 28 416
unknowns, 1.5e13 flops at the measured 19 GFlop/s) and dense at 18944 on any thread count: the
saddle matrix is `(3·18944)² · 8 B` = 26 GB, plus the Gram matrices and a copy for the residual
check (about 56 GB peak, extrapolating the measured 3.6 and 14.2 GB with a factor 4 per doubling; 94 GB
machine shared with other sessions) and about 16 min per LU on 16 threads (1.2e14 flops at the
measured 125 GFlop/s). Dense peak memory measured: 3.6 GB at 4736, 14.2 GB at 9472; multigrid 0.3 GB at 18944.

**64-thread speed-up** (1 → 64 threads; dense pool 16): Exact 4736: 113 → 24.3 s (4.7×); BH +
dense 4736: 104 → 21.5 (4.8×); Multigrid 592: 2.8×, 2368: 3.1×, 4736: 2.7×, 9472: 2.4×,
18944: 2.3× (27.7 → 12.3 s). The dense configurations are limited by the LU (below), the
multigrid ones by the projections, whose work is mostly sequential (below).

**Dense pool size** (Exact, 2368, 64 main threads, 8 steps): dense pool 8: 3.71 s, **16: 3.12**,
32: 3.46, 64: 6.33; the LU alone on 24 576 unknowns: 4 threads 159 s, 8: 124, **16: 121**, 32:
126, 64: 136 (2048: 16 threads 1.63 s, 32: 1.79, 64: 4.65). The default `DEFAULT_DENSE_THREADS =
16` is the measured optimum (§17b); unchanged.

### Where the time goes

* **Dense (exact or BH):** the metric (Gram assembly + saddle LU) is 87 % of an Exact iteration at
  4000 vertices (60.6 of 69.5 s) and 96 % of a Barnes–Hut one (60.6 of 63.0 s); the exact energy
  costs 4.1 s (6 %), Barnes–Hut 0.56 s (0.9 %); direction solve 0.3 %, line search 2–4 %. At 64
  threads (dense pool 16) the LU is still 93–96 %.
* **Multigrid: the line search, i.e. its projections, dominates**: projections are
  67 % of the iteration at 592 vertices, 59–60 % at 1184–2368, 64 % at 4736, 69 % at 9472 and
  **76 % at 18944** on one thread (21.0 of 27.7 s); at 64 threads 74 % at 592 rising to **80 %**
  at 18944 (9.83 of 12.3 s). The knot is worse: 83 % at 16384 (179 of 216 s). Per iteration the
  projections cost about 6 direction solves (21.0 s vs 3.4 s at 18944; 0.73 vs 0.12 at 592),
  because each trial of the line search runs a chord projection of 1–3 multigrid solves. The
  rest at 18944 / 1 thread: direction solve 12 %, hierarchy build 1.7 %, energy 4.7 %,
  objective evaluations of the line search 2.6 %, CCD 1.2 %. Even with the barycenter as the
  only (linear) constraint the projection phase takes 80–85 % of the line search (`mg-noproj`): the
  multigrid direction satisfies `C p = 0` only to the solver tolerance, which is consistent with
  every trial needing a chord solve even for a linear constraint (not verified here). Without
  the projections (total − projection) an iteration at 16000 would cost 4.8 s instead of 14.8 s
  (3×): the paper's gain from dropping the projection is 1.5× (57 → 38 s). This is the first
  thing to speed up (cheaper or fewer chord solves, or a convergence test that accepts the
  predicted state); not done here.
* **CCD** is under 3 % of every configuration (1.2 % of the multigrid iteration at 18944 on one
  thread; 4.8 % in the barycenter-only run there).

### Thread scaling (§14.6)

Components on Freedman–He–Wang resampled to 2048, 8192, 32768 vertices (2, 4.5); `best` of 3–20
timed runs after one warm-up, 1/2/4/8/16/32/64 threads, one pool per run. Efficiency
`T₁ / (t · T_t)`. Components: exact energy + differential (`energy-exact`), BH energy +
differential (`energy-bh`, default θ), dense Sobolev Gram assembly (`gram`: `B`, `B⁰`, `A`), dense
saddle LU (`lu`: assemble the `3n + k` matrix and factor; the solve follows sequentially),
block cluster tree build and one product (`bct-build`, `bct-product`), multigrid hierarchy build
and one solve to the default tolerance (`mg-build`, `mg-solve`), CCD of one typical step
(`ccd`, all pairs with a motion of 0.2 mean edge lengths, no obstacle). Dense components are
skipped at 32768 (`gram` 3 × 8.6 GB, `lu` 3·32768 = 98 304 unknowns, 77 GB). `lu` at 8192 was run
from 4 threads (1 thread would be about 13 min; efficiency is relative to 4 threads there).

**Bitwise equality: every row (173) has the same result hash as the first thread count of its series** (hash of
the differential and energy, the Gram matrix, the LU solution, the product / solution vectors,
the CCD report): D25 holds for every parallel component, including the dense LU.

| n | component | T₁ (s) | eff 2 | 4 | 8 | 16 | 32 | 64 | fastest (thr, s, speed-up) |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---|
| 2048 | energy-exact | 1.052 | 0.94 | 0.98 | 0.98 | 0.95 | 0.81 | 0.49 | 64, 0.0338, 31× |
| 2048 | energy-bh | 0.484 | 0.94 | 0.90 | 0.90 | 0.81 | 0.53 | 0.28 | 64, 0.0265, 18× |
| 2048 | gram | 0.387 | 0.84 | 0.74 | 0.53 | 0.35 | 0.19 | 0.09 | 32, 0.0635, 6.1× |
| 2048 | lu | 8.22 | 0.87 | 0.71 | 0.52 | 0.32 | 0.14 | 0.03 | 16, 1.63, 5.1× |
| 2048 | bct-build | 0.0392 | 1.04 | 0.58 | 0.29 | 0.15 | 0.08 | 0.03 | 16, 0.0161, 2.4× |
| 2048 | bct-product | 0.00262 | 0.89 | 0.73 | 0.57 | 0.34 | 0.15 | 0.05 | 16, 0.00049, 5.3× |
| 2048 | mg-build | 0.242 | 0.71 | 0.44 | 0.26 | 0.13 | 0.07 | 0.03 | 32, 0.114, 2.1× |
| 2048 | mg-solve | 0.527 | 0.80 | 0.62 | 0.41 | 0.22 | 0.10 | 0.04 | 16, 0.149, 3.5× |
| 2048 | ccd | 0.00135 | 0.73 | 0.47 | 0.25 | 0.12 | 0.05 | 0.02 | 8, 0.00069, 2.0× |
| 8192 | energy-exact | 17.1 | 0.98 | 0.99 | 0.99 | 0.98 | 0.95 | 0.60 | 64, 0.446, 38× |
| 8192 | energy-bh | 0.719 | 0.96 | 0.90 | 0.81 | 0.69 | 0.44 | 0.24 | 64, 0.0469, 15× |
| 8192 | gram | 6.40 | 0.91 | 0.69 | 0.55 | 0.35 | 0.20 | 0.10 | 64, 1.01, 6.3× |
| 8192 | lu | 159 (at 4 thr) | – | 1.00 | 0.64 | 0.33 | 0.16 | 0.07 | 16, 121, 5.3× vs 1 thr (est.) |
| 8192 | bct-build | 0.0532 | 0.90 | 0.49 | 0.27 | 0.14 | 0.07 | 0.03 | 16, 0.0242, 2.2× |
| 8192 | bct-product | 0.00361 | 0.71 | 0.62 | 0.38 | 0.21 | 0.09 | 0.04 | 16, 0.00106, 3.4× |
| 8192 | mg-build | 0.397 | 0.69 | 0.40 | 0.22 | 0.11 | 0.06 | 0.03 | 32, 0.220, 1.8× |
| 8192 | mg-solve | 1.35 | 0.74 | 0.51 | 0.33 | 0.17 | 0.08 | 0.03 | 16, 0.487, 2.8× |
| 8192 | ccd | 0.00547 | 0.75 | 0.46 | 0.25 | 0.13 | 0.06 | 0.03 | 16, 0.0026, 2.1× |
| 32768 | energy-exact | 271 | 1.03 | 0.98 | 0.98 | 0.98 | 0.98 | 0.61 | 64, 6.93, 39× |
| 32768 | energy-bh | 1.72 | 0.96 | 0.86 | 0.79 | 0.60 | 0.36 | 0.18 | 64, 0.147, 11.7× |
| 32768 | bct-build | 0.140 | 0.81 | 0.51 | 0.28 | 0.14 | 0.07 | 0.03 | 8, 0.0627, 2.2× |
| 32768 | bct-product | 0.00969 | 0.72 | 0.44 | 0.30 | 0.14 | 0.06 | 0.03 | 8, 0.00409, 2.4× |
| 32768 | mg-build | 0.856 | 0.61 | 0.33 | 0.18 | 0.09 | 0.04 | 0.02 | 16, 0.590, 1.4× |
| 32768 | mg-solve | 3.57 | 0.78 | 0.45 | 0.25 | 0.12 | 0.05 | 0.02 | 8, 1.80, 2.0× |
| 32768 | ccd | 0.0192 | 0.63 | 0.37 | 0.20 | 0.11 | 0.05 | 0.03 | 32, 0.0111, 1.7× |

Reading: the exact energy scales to 32 threads at 95–98 % efficiency on 8192 and 32768 vertices
(81 % on 2048; at 64 threads 49–61 %: the machine has 64 hardware threads and the OS shares them).
Barnes–Hut falls off earlier (90 % at 8 threads and 53 % at 32 on 2048; 79 % and 36 % on 32768;
not investigated). The dense Gram assembly and the LU reach 5–6× at best, with the optimum at 16
threads (LU) and 32–64 (Gram); the LU of 2048 vertices is slower at 64 threads than at 4. The hierarchical components
(BCT build and product, multigrid build and solve, CCD) are **sequential-bound at these sizes**:
they stop gaining at 8–16 threads (1.4–3.5×) and the multigrid solve is *slower* at 64 threads
than at 16 on 2048, 8192 and 32768 vertices (0.198 vs 0.149 s; 0.662 vs 0.487 s; 2.51 vs 1.85 s).
One BCT product at 8192 vertices takes 3.6 ms on one thread while one multigrid solve takes 1.35 s
(a solve is about 2 cycles of 6 CG smoothing steps and a residual on each of 3–4 levels, so on
the order of 80 products of decreasing size): the products are a part of the solve, not
obviously all of it, and where the rest goes (coarse levels, the coarsest dense solve, the
projectors) is not identified here. Follow-up: profile `MultigridSaddle::solve`, and cap the
multigrid pool at 8–16 threads (as the dense pool is capped at 16) or parallelize what the
profile shows.

### Auto crossover (§17h re-measured)

Exact + dense against BH + multigrid, mean seconds per iteration of the first 20 steps of the
flow (Freedman resampled to the sizes, (2, 4.5)), machine alone. **Default session threads**
(all 64 cores, dense pool 16) and the 8-thread configuration of §17h (8 and 8):

| n | 64 thr: exact + dense | BH + MG | ratio | 8 thr: exact + dense | BH + MG | ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 600 | 0.148 | 0.348 | 0.43 | | | |
| 800 | 0.253 | 0.492 | 0.52 | 0.283 | 0.401 | 0.71 |
| 1000 | 0.385 | 0.606 | 0.64 | | | |
| 1200 | 0.568 | 0.815 | 0.70 | 0.689 | 0.737 | 0.93 |
| 1400 | 0.825 | 0.980 | 0.84 | | | |
| 1500 | 0.927 | 1.085 | 0.85 | | | |
| 1600 | 1.108 | 1.093 | 1.01 | 1.345 | 0.913 | 1.47 |
| 1700 | 1.301 | 1.130 | 1.15 | | | |
| 1800 | 1.478 | 1.133 | 1.30 | | | |
| 2000 | 1.919 | 1.235 | 1.55 | 2.370 | 1.086 | 2.18 |
| 2500 | 3.478 | 1.536 | 2.26 | 4.361 | 1.304 | 3.35 |
| 3000 | 5.827 | 1.840 | 3.17 | 7.593 | 1.459 | 5.21 |
| 4000 | 14.05 | 2.252 | 6.24 | | | |

(ratio > 1: multigrid is faster.) The crossover is at about **1600 vertices at 64 threads and
1300 at 8 threads**, below §17h's 2000, which was measured under load 30–40 and at 8 threads
with a slower dense side. The 64-thread rows at 1800 and 2000 were measured twice (ratios 1.25 and 1.30; 1.48
and 1.55; the table shows the second pass), which gives the run-to-run spread, 5 %. Scaling
between 2000 and 4000 vertices: exact + dense n^2.9, BH + multigrid n^0.9.

**Rule (replaces §17h's; superseded by §17l: 2500 with the Schur dense solve).** `Backend::Auto { dense_max_vertices: 1500 }`: exact + dense for
`|V| ≤ 1500`, BH + multigrid above; at the worst point inside 1300…1600 either side loses 30 %
(8 threads at 1500: exact + dense about 1.3× slower than multigrid, interpolated; 64 threads at
1500: multigrid 1.2× slower). **Changed from 2000** in `src/step/mod.rs` (`DEFAULT_AUTO_DENSE_MAX_VERTICES`,
rustdoc), the Python default (`ddg_tools.pyi`, `python-tests/test_repulsive.py`; the option's
default itself follows the Rust constant) and the Blender tooltip; the test of the rule
(`tests/session_backends.rs`) now checks 1500/1501. At 2000 vertices the old default was
1.5× (64 threads) and 2.2× (8 threads) slower than the accelerated pair. A faster projection
(above) would move the crossover down.

### The authors' program, alone (§17c re-run)

`rcurves_app --batch --steps 6 (4 at 1 thread) --no-obj --no-exact-energy --ignore-stop`, scenes
`freedman{600,1200,2400}`, (2, 4.5), barycenter + edge lengths, `OMP_NUM_THREADS = --threads =
T`; mean of steps 1…N of `steps.csv` (`seconds_step`), the machine alone (load 1–2 at 1 thread).
Their "Barnes–Hut" is BH energy + dense solve, "multigrid" is BH + multigrid (their
`--multigrid`), "exact" is `--no-bh` (exact energy of their vertex discretization, dense).
Ours: the same curves, same configurations; steps of the two programs differ (line search,
step-size rule), so these are times per step of two flows from the same start, not per equal work.
Seconds per step, **ours / authors'**:

| curve | config | 1 thread | 34 threads (their default) | ours at 64 |
|---|---|---:|---:|---:|
| 592 | exact + dense | 0.451 / 1.362 | 0.151 / 1.310 | 0.154 |
| 592 | BH + dense | 0.463 / 0.718 | 0.154 / 0.421 | 0.155 |
| 592 | BH + multigrid | 1.08 / 2.995 | 0.385 / 1.492 | 0.392 |
| 1184 | exact + dense | 2.58 / 6.289 | 0.576 / 5.226 | 0.556 |
| 1184 | BH + dense | 2.46 / 3.278 | 0.585 / 1.646 | 0.556 |
| 1184 | BH + multigrid | 2.60 / 5.967 | 0.954 / 2.935 | 0.814 |
| 2368 | exact + dense | 16.1 / 30.8 | 3.08 / 21.8 | 2.98 |
| 2368 | BH + dense | 14.2 / 17.7 | 3.04 / 7.54 | 2.68 |
| 2368 | BH + multigrid | 4.27 / 16.1 | 1.58 / 7.23 | 1.40 |

Ours is faster in every cell: 1.25–3.8× at one thread, 2.5–9× at 34 threads. Their dense solve
scales 2.2× from 1 to 34 threads at 2368 vertices (BH + dense: 16.4 → 7.3 s), ours 4.9× (the
metric: 12.9 → 2.6 s). Their split at 2400
(ms/step, 1 thread): BH energy 1004, solve 16 412 (dense) / 5 946 (multigrid), line search 101,
**back-projection 9 028 with multigrid** (56 % of the step; 59 % at 34 threads): the authors'
multigrid program is also dominated by the constraint projection, the same effect as in §17k above.
Their exact energy costs 13.3 s at 2400 against 1.4 s for ours (different discretizations, D1).
Final energies are not compared here.

### Follow-ups found

1. **Multigrid line search**: 60–80 % of a multigrid iteration; about 6 solves worth of work per
   iteration, and it scales worst (2.3× on 64 threads at 18944).
2. **Multigrid solve and hierarchy build do not scale** past 8–16 threads and get slower at 64;
   a cap on the pool, or a profile of what in the solve is sequential.
3. **Dense LU** at 19 GFlop/s single-core, 5× on 16 threads; the Schur solve of §17b
   (Cholesky of `A + ρ𝟙𝟙ᵀ`, about 6× fewer flops) would move the dense side and the crossover up
   and would close the 5 % gate miss of Barnes–Hut + dense at 4000.
4. **Knot 96** is 2–5× above the paper's multigrid curve from 2048 vertices on one thread:
   projections of the near-contact knot (83 % of the iteration at 16384).
5. `Backend::Auto` is untouched in its rule of form; a thread-aware crossover (1300 at 8 threads,
   1600 at 64) would gain at most 20 % near the crossing.

## 17l. Performance after WP14b (2026-10-26)

Logs, job files and CSVs: `runs_tmp/rc-wp14b/` (`projdiag-*.log`, `diag-*.log`, `trial/` and
`trial2/` for the knot study, `fig20-*.csv`, `components.csv`). Benches: `benches/fig20.rs`
(new flags `--trace 1`, `--soc`, `--polish`, `--proj-it`, `--bct-tt`, `--bct-tx`,
`--mg-threads`; `mg-noproj` is now the real switch, the barycenter-only proxy of §17k is
`mg-bary`), `benches/projdiag.rs` (the projection of the trials of one step with several
solvers, violation history per chord iteration), `benches/trialrule.rs` (the knot study).

### 1. Why the projection cost about six direction solves per iteration

Measured with `--trace 1` (every trial: `τ/τ₀`, chord iterations, outcome) and `projdiag`
(Freedman resampled to 2000 / 8000, knot 96 at 2048 / 8192, multigrid, 1 thread):

* **Not the causes suspected in §17k.** The direction satisfies `C p = 0` to rounding (the
  multigrid's `x = P y` uses the exact sparse projector), so the violation of a trial is the
  curvature term `Φ(γ − τp) = τ² ½Φ''[p, p] + O(τ³)`; even the barycenter row is nonlinear
  (length-weighted, D6). The hierarchy is built once per step and reused by every solve. The
  multigrid accuracy is not the problem either: the chord histories with multigrid (tolerance
  1e-3) and with the dense LU agree to three digits (`1.46e3 → 3.97e2 → 2.81e1 → 2.13e0 →
  2.42e-1` both); a looser multigrid tolerance (1e-1) makes the chord iteration converge more
  slowly (8 instead of 4 iterations at `τ₀`) and is not cheaper.
* **The cause: failed trials.** With the grow-previous start the first trial is `2 τ_prev`;
  its chord iteration (rate 0.07–0.27 at that step size) usually needs 4 iterations and fails at
  the limit of 3 (`MaxIterations`), wasting 3 solves; the halved trial then takes 2–3. Typical
  steps: `[1.000:3:MaxIterations 0.500:2:ok]`, `[1.000:3:MaxIterations 0.500:3:MaxIterations
  0.250:2:ok]`; 4.4 projection solves per step on Freedman 2000 (20 steps), 5.3 on knot 2048,
  each a full multigrid solve (full-multigrid start plus 1–2 FCG cycles, 0.47 s at 2000 = the
  direction solve). On the knot, `τ₀ = ⅔ τ_max` binds (the CCD) and that trial failed in every
  step.
* **The knot's solves are 4–15× more expensive** than Freedman's at the same size because its
  block cluster trees are almost dense (item 3 below).
* **The authors' program** (`ref: include/tpe_flow_sc.h`, `LSBackprojectMultigrid`) runs the
  Armijo search on the *unprojected* trials and projects only the accepted step: at most 2
  chord steps per attempt with the multigrid (3 dense), the same threshold `1e-4 ℓ̄`
  (`mg_backproj_threshold`), halving `τ` on failure, up to 10 attempts; each chord step is a
  multigrid solve from scratch (`VCycleSolve` at tolerance 1e-2). It is projection-dominated too
  (56–59 %, §17k). We keep the projected Armijo test (D11: the guarantee is on the projected,
  CCD-checked path).

**Fix 1: the Newton polish (default).** Once the largest relative violation of the projection is
at most `newton_polish = 100` times its tolerance, the remaining iterations are Newton steps
with the Jacobian *at the iterate* and the lumped-mass metric,
`x = M⁻¹Cᵀ(CM⁻¹Cᵀ)⁻¹(−Φ)` (`solve::newton_correction`; one sparse Cholesky factorization of
the `k × k` matrix, `O(k)`, about 4 ms at 2000 vertices against 470 ms for a multigrid solve),
instead of chord steps with `Ā` and `C` frozen at the base curve. Newton converges quadratically
(`6.30e0 → 1.87e-3`, `1.05e2 → 2.25e0 → 2.41e-3` relative violation) where the chord method
converges linearly at a rate proportional to the step. Every correction above 100 × tolerance
(`1e-2 ℓ̄` for the default) is still the `Ā`-smallest chord step of Eq. 23, so the large, smooth
part of every projection is the paper's; the polish only places the last `≤ 1e-2 ℓ̄`, where the
metric does not matter. A rank-deficient Jacobian at the iterate (a collapsing edge) falls back
to chord steps. **This departs from the paper's chord projection**; `newton_polish: None`
restores it (the Lean oracle replays the paper's step and pins it).

**Fix 2 (option, off): the second-order correction.** Trials `γ − τp + τ²q` with `C q = −Φ₂`,
`Φ₂ = ½Φ''[p, p]` by a central second difference of the constraint values (step
`h max|p| = 1e-3 ℓ̄₀`), one solve per step (`solve::second_order_correction`; the second-order
correction of SQP methods, Nocedal & Wright §18.3). The piece from `γ − τp` to the predicted
trial is CCD-checked like a chord step. It cuts the trial violation 18–80× (`1.46e3 → 8.2e1` at
`τ₀`, `3.24e2 → 6.2e0` at `τ₀/2`, trials at `τ₀/4` already within the tolerance), but the step
then grows and the failing first trial moves to larger `τ`: on the knots it saves no time
(below), so it is off by default (`second_order_correction: false`).

**Fix 3 (found by WP12b): every trial gets at least one projection step.** A trial already
within the tolerance used to be accepted unprojected. Tangential steps lengthen the edges at
second order, which lowers `Ê`, so the curve drifted unprojected to `|Φ| = tol`, and every larger
trial was then projected back and paid `λ · Σ δℓ`, far more than the Armijo decrease: on a
40-vertex rough trefoil (exact + dense, (2, 4.5)) `‖g̃‖` stalled at 1.1e-2 with the base
violation at the tolerance and 989 of 1000 accepted steps unprojected; the floor scaled with the
tolerance (4.9e-3 / 1.0e-3 / 9e-5 / 7e-6 for `MeanEdgeLength` 1e-4 / 1e-6 / 1e-8 / 1e-10).
Now `project_constraints` and the line search's projection take at least one step unless
`Φ = 0` exactly (the paper: "update γ̃ ← γ̃ + x and repeat until …"; the reference always calls
`BackprojectConstraints`); the same flow reaches `‖g̃‖ < 1e-4` (test
`hard_flow_of_a_rough_trefoil_reaches_a_small_gradient`; 1.6e-5 at step 1000 in WP12b's run).
The Lean oracle (`chordProject`) follows; fixture `rc-circle-step` changes (its trial now takes
one chord step). With the polish the extra step is a Newton step (cheap).

**Fix 4: early abort of hopeless chord projections (default).** In the knot's steps the first
trial (`2 τ_prev`) still failed after 3 chord solves almost every time
(`[1.000:3:MaxIterations 0.500:3:MaxIterations 0.250:4:ok]` at 8192). A chord iteration
converges linearly, so after one step its rate `θ = h_k/h_{k−1}` predicts how many more steps
reach the target (the polish threshold, or the tolerance); the projection now gives the trial up
(`MaxIterations`, fewer iterations) when `θ ≥ 1` or that number exceeds the remaining budget —
the rate test of simplified-Newton iterations in implicit integrators (Hairer & Wanner, *Solving
ODEs II*, §IV.8, (8.12)). Example (rough trefoil, 3 iterations allowed): `2836 → 591` (`θ = 0.21`,
4.1 more steps needed) stops after 1 solve instead of 3 and ends rejected either way; converging
trials are never stopped (test `early_abort_stops_hopeless_projections_only`).
`FlowOptions::projection_early_abort`, default `true` (below); `false` is the paper's.

**Measurement of the defaults (as §17f).** 24 perturbed Knot128 knots (every 5th, 1500–2750
vertices, all steps on Barnes–Hut + multigrid by `Backend::Auto`), (3, 6), `Absolute(1e-4)`,
3 chord iterations, 100 iterations, CCD audit of every accepted path (`path_contact`, gap 1e-6,
1e5 advances), 8 processes × 4 threads on a shared machine (load 25–35; per-step ratios are
paired per knot). `E*` is the best energy any variant reached on that knot. With fix 3 in place
(`runs_tmp/rc-wp14b/trial2/`, `table.md`; the run before fix 3 in `trial/` gave the same
ranking):

| trial rule | stuck steps | isotopy violations | trials/step | projection solves/step | s/step (median) | paired s/step | reached 1.1 E* | paired time to 1.1 E* (n) | median final E / E* |
|---|---|---|---|---|---|---|---|---|---|
| straight trials, chord projection (the paper) | 0 / 2400 | 0 | 1.98 | 4.25 | 1.67 | 1 | 19 / 24 | 1 | 1.025 |
| second-order correction | 0 / 2400 | 0 | 1.97 | 4.46 | 1.70 | 1.03 (1.01–1.07) | 23 / 24 | 0.97 (0.93–1.04, 19) | 1.012 |
| correction + Newton polish | 0 / 2400 | 0 | 1.98 | 2.43 | 1.20 | 0.75 (0.70–0.84) | 19 / 24 | 0.71 (0.63–0.85, 18) | 1.015 |
| **Newton polish (default)** | 0 / 2400 | 0 | 1.97 | 2.31 | 1.20 | 0.69 (0.66–0.80) | 21 / 24 | 0.69 (0.62–0.81, 18) | 1.000 |

(projection solves/step include the correction's solve; the paired time to 1.5 E* and 2 E*:
polish 0.67 and 0.71, correction 1.01 and 1.02.) Decision: **`newton_polish = Some(100)` by
default** (as reliable, 0.69× the time to the energy target, best final energies);
**`second_order_correction = false`** (no faster alone, slower with the polish).

The early abort on top of the polish (same setup, `trial3/`, machine less loaded):

| trial rule | stuck steps | isotopy violations | projection solves/step | s/step (median) | paired s/step | reached 1.1 E* | paired time, iterations to 1.1 E* (n) | median final E / E* |
|---|---|---|---|---|---|---|---|---|
| Newton polish | 0 / 2400 | 0 | 2.31 | 0.96 | 1 | 24 / 24 | 1, 1 | 1.000 |
| **polish + early abort (default)** | 0 / 2400 | 0 | 1.54 | 0.78 | 0.79 (0.77–0.82) | 24 / 24 | 0.76 (0.74–0.79), 1.00 (24) | 1.000 |

Same iterations to every target (the aborted trials would have failed), 0.74–0.76× the time:
**`projection_early_abort = true` by default.** The Lean oracle replays the paper's step
(straight trials, chord projection to its limit) and pins all three options off.

### 2. The back-projection switch

`FlowOptions::back_projection` (default `true`). `false` is the paper's "Multigrid without
projection": the direction is still the projected gradient of Eq. 22 (tangent to the
constraints), the trials are taken unprojected, nonlinear rows drift by `O(τ²)` per step, and
the re-target sub-step runs only for a moving length schedule. `benches/fig20.rs`'s `mg-noproj`
uses it (§17k used a barycenter-only proxy, now `mg-bary`). Test:
`without_back_projection_no_projection_solve_runs_and_the_lengths_drift`.

### 3. Knot 96: the block cluster trees were almost dense

The multigrid's trees used Eq. 28's tangent test (`theta_t 0.25`). On the perturbed knot every
cluster of more than a few edges has a large tangent spread, so almost no block is admissible:
level 0 stores 1930 exact entries per row at 2048 vertices (94 % of the dense matrix), 1980 at
8192, 2640 at 4096; Freedman: 90–260. The metric's kernel `K^B` does not depend on the tangents
and `K⁰` only through its far coefficient; the reference admits blocks by position only.
`MgOptions::default().bct.theta_t = 2` (`MG_BCT_THETA_T`; every unit-tangent cluster passes,
i.e. position only; `BctOptions::default()` is unchanged). Multigrid direction against the dense
direction, and its solve time (projdiag, after 1 step, machine shared):

| curve | n | `theta_t 0.25`: exact entries/row, error, s | `theta_t 2`: entries/row, error, s |
|---|---:|---|---|
| knot 96 | 1024 | 1016, 7.5e-5, 0.41 | 703, 2.2e-3, 0.32 |
| knot 96 | 2048 | 1850, 9.9e-4, 2.07 | 780, 6.8e-3, 0.89 |
| knot 96 | 4096 | 2641, 2.4e-3, 7.11 | 706, 8.5e-3, 2.49 |
| Freedman | 2000 | 238, 4.4e-3, 0.44 | 122, 6.7e-3, 0.18 |
| Freedman | 4000 | 89, 4.0e-3, 0.72 | 72, 6.8e-3, 0.57 |

All errors stay below the 1e-2 target of §13.4; the multigrid tests (`tests/multigrid.rs`,
`tests/multigrid_networks.rs` with the capillary scene) and the Freedman 600 flows (exact energy
6.0816 within 0.5 %) pass.

### 4. Multigrid thread scaling

Profile of one BCT product (32768 vertices, 1 thread): the two parallel loops (near field per
row, far field per node) are about 54 % of a multigrid solve's time and 83 % of a product; the
upward and downward passes over the nodes and the `ψ`/output loops of
`HierarchicalMetric::apply_unchecked` are sequential (about 17 % of a product), and so are the
projector's sparse triangular solves and the CG vector operations of the smoother. With rayon's
overhead per parallel region these cap the speed-up at 2–3.5×, and at 64 threads the regions'
cost exceeds their gain (§17k). The sequential parts are in `src/spatial` and `src/metric`,
outside this package. Fix here: `MgOptions::max_threads` (default `Some(16)`,
`MG_DEFAULT_MAX_THREADS`): `MultigridSolvers` builds and solves on a pool of its own with at most
16 threads when the session's pool is larger (created once per factory; results bitwise the
same, test `multigrid_thread_cap_does_not_change_the_result`). Measured below.

### 5. Dense solve: the Schur complement (§17b)

`DenseSaddle` factors by the Schur method (module docs of `solve::dense`): `M = A + ρ uuᵀ`
(`u = 𝟙/√n`, `ρ = max_v A_vv`, positive definite since `ker A = span{𝟙}`, and `M u = ρu`),
Cholesky of `M` and its explicit inverse, `S = C(M⁻¹ ⊗ I₃)Cᵀ` (`k × k`, from the sparse rows in
`O(n · nnz C)`), Cholesky of `S`, and the `3 × 3` matrix `W = Eᵀ S⁻¹ E`, `E = CU/ρ`, for the
translations (nonsingular iff the constraints fix them, D5); every solve does one step of
iterative refinement. `k` (≈ n with edge lengths) does not make it as large as the original:
`S` is `k × k` but costs only `k³/3`; total about `1.3 n³` flops against `(2/3)(4n)³ ≈ 43 n³` for
the LU of the `(3n + k)` matrix, memory `≈ 4n² + k²` doubles against `(3n + k)²`. It falls back
to the LU when a Cholesky factor fails or has `min L_ii²/max L_ii² ≤ 1e-13`, so singular saddles
are still decided (and reported) by the LU. Accuracy (`tests/dense_schur.rs`, against the LU on
closed curves, an open helix with pins and an end tangent, a theta network, a planar loop, few
rows, the mean-free rows of a soft problem, Freedman 592; right sides of Eq. 22, of a projection
and both): `x` agrees to 8e-16–3.3e-13 relative, `λ` to 3e-15–1.8e-12, the Schur residual is
as small as the LU's or smaller (4e-15–2.7e-11 against 1.6e-14–2.4e-10). Bitwise identical at
1/3/8/16 threads on the test curves and at 1/7/16/32 threads at 8192 vertices (result hashes).
Time of the factorization on one thread (shared machine): 2048 vertices 0.58 s against 8.69 s
for the LU (15×), 8192 vertices 36 s (the LU about 13 min, §17k).

### 6. Contact floor of the guarded line search (requested after the evaluation's smoke runs)

With metrics blind to the repulsion (H² projected gradient, soft H¹ L-BFGS, the final curve of
the H¹ projected gradient), the CCD bounded every step by a fraction (the gap, 0.1) of the
current distance of the closest pair, which then shrank geometrically until two edges touched
in floating point (1e-15…0): `validate_positions` failed and the path audit reported EdgeEdge
contacts; Hs runs never came close. Now the guarded search first checks the base curve: two
disjoint edges closer than `CONTACT_FLOOR · ℓ̄₀ = 1e-6 ℓ̄₀` count as a contact, `τ_max = 0`, no
step, `StepStatus::Stuck` (`StepReport::contact_floor`); the last accepted curve keeps a distance
of at least `gap · 1e-6 ℓ̄₀ = 1e-7 ℓ̄₀`, 10⁵ times the validation tolerance (`1e-12 ℓ̄`). The
broad phase is a uniform grid (cell = longest edge), `O(|E|)` per step. Choice of `1e-6`: far
above rounding and the validation tolerance, far below anything a Sobolev flow reaches (the
tangent-point energy of a pair at `10⁻⁶ ℓ` is about 10¹⁸ times its value at `ℓ` for (3, 6)), so
Hs flows are bitwise unchanged unless they hit it (they do not: test
`contact_floor_is_never_hit_by_a_sobolev_flow`). Regression test
`contact_floor_stops_a_flow_that_drives_two_edges_together`: an L² flow with a penalty that
overpowers the energy pulls two strands together; with the floor set to 0 the curve fails
`validate_positions` after step 25 (edges 9 and 29 at 1.8e-13); with it the flow ends Stuck at
step 20 at 5.3e-7 ℓ̄₀, every path embedded. Obstacles (vertex–triangle) are not covered by the
floor (no failure seen).

### 7. Timing alone (2026-10-26/27)

Protocol as §17k (`runs_tmp/rc-wp14b/runjobs.sh`, one job after another, `load_start`/`load_end`
columns; a launcher waited for the 1-minute load to stay below 4). No other agent's job was
running (`ps` spot checks); the remaining load (1.1–2.9 at the 1-thread starts) came from an unrelated desktop process
(about 1.5 cores) and from the lag after the benchmark's own multi-threaded runs. Defaults of
this section: Newton polish, early abort, no second-order correction, `theta_t 2`, Schur dense,
multigrid pool capped at 16. CSVs: `fig20-gate2.csv` (the gate), `fig20-long.csv`,
`fig20-items.csv`, `fig20-items2.csv`, `fig20-mgthreads.csv`, `fig20-cross.csv`,
`fig20-cross2.csv`, `components.csv`; the first gate pass `fig20-gate1.csv` (before the early
abort) differs from the second by < 5 % except on the knot (knot 8192: 46.2 → 28.1 s).

**Fig. 20 gate, 1 thread, seconds per iteration** (Freedman resampled, (2, 4.5); mean over the
first 20 steps or the 240 s budget):

| configuration | n | paper | §17k (WP14) | now | verdict |
|---|---:|---:|---:|---:|---|
| Exact | 1000 | ≈ 6 | 1.64 | 0.756 | below |
| Exact | 4000 | ≈ 145 | 69.5 | 13.9 (17 steps) | below |
| Barnes–Hut + dense | 2000 | ≈ 10 | 9.01 | 1.84 | below |
| Barnes–Hut + dense | 4000 | ≈ 60 | 63.0 | 5.90 | **below** (was 5 % above) |
| Multigrid | 2000 | ≈ 6 | 3.57 | 1.74 | below |
| Multigrid | 8000 | ≈ 25 | 8.39 | 3.67 | below |
| Multigrid | 16000 | ≈ 57 | 14.8 | 5.19 | below |
| Multigrid without projection (`back_projection = false`) | 2000 | ≈ 3 | 3.13 (proxy) | 1.33 | below |
| same | 16000 | ≈ 38 | 11.8 (proxy) | 4.09 | below |

Knot 96 ((3, 6), multigrid, 1 thread), against the paper's multigrid curve (6 / 25 / 57 s at
2000 / 8000 / 16000):

| n | §17k | now (240 s budget) | now, 20 steps | verdict |
|---:|---:|---:|---:|---|
| 1024 | 2.85 | 0.845 | | below |
| 2048 | 11.7 | 2.45 | | **below** (was 2× above) |
| 4096 | 41.0 (5 steps) | 6.98 | | below |
| 8192 | 127 (2) | 28.1 (8 steps; median 25.2) | 18.5 (median 18.2) | **below over 20 steps**; the first 8 steps 12 % above |
| 16384 | 216 (2) | 40.4 (5 steps) | 33.0 (median 31.1) | **below** |

(Knot without projection: 1.43 s at 2048, 21.6 at 16384. Exact + dense on the knot: 0.470 /
2.16 / 9.89 s at 1024 / 2048 / 4096, was 1.52 / 9.92 / 72.0.) The paper averages up to 500
iterations; the first steps of the perturbed knot are the most expensive (2–3 failed trials per
step), so the 20-step means are the closer comparison; at 8192 the 240 s budget stops after 8.

**Before / after per item** (1 thread, 10 steps unless noted; `fig20-items*.csv`):

| change | Freedman 8000 | knot 2048 | knot 8192 (5 steps) |
|---|---:|---:|---:|
| WP14 configuration (chord projection, `theta_t 0.25`) | 9.43 | 12.5 | 127 (§17k, 2 steps) |
| + `theta_t 2` only | 6.15 | 5.52 | 43.6 |
| + Newton polish only (`theta_t 0.25`) | 5.90 | 10.5 | 131 (2 steps) |
| `theta_t 2` + polish, no early abort | 5.05 | 4.80 | 50.8 |
| `theta_t 2` + polish + early abort (default) | 4.52 | 3.86 | 35.1 |
| `theta_t 2` + polish + second-order correction | | 4.63 | |

Dense: the metric phase (Gram + factorization) of Barnes–Hut + dense at 4000 is 3.99 s against
60.6 s with the LU. Schur factorization alone (`components.csv`, result hash identical across
each row's thread counts): 2048 vertices 0.51 / 0.27 / 0.21 s at 1 / 4 / 16 threads (LU 8.32 /
2.90 / 1.63), 4096: 3.50 / 1.45 / 1.13 / 1.00 / 1.02 s at 1 / 4 / 8 / 16 / 32, 8192: 22.9 / 7.2 /
6.2 / 5.4 s at 1 / 8 / 16 / 32.

**Multigrid threads.** Components with `theta_t 2` (best of 3): solve at 8192 vertices 0.82 /
0.58 / 0.47 / 0.40 / 0.39 / 0.43 / 0.54 s at 1 / 2 / 4 / 8 / 16 / 32 / 64 threads; at 32768
2.66 / 2.12 / 1.56 / 1.56 / 1.59 / 1.71 / 2.04; build 0.19 → 0.14 s (8192) and 0.61 → 0.51 s
(32768) at best; one BCT product 4.8 → 1.3 ms at 8 threads (8192). Flows at 64 main threads
(`fig20-mgthreads.csv`, 10 steps): Freedman 16000 2.83 s without a cap, 2.36 / 2.37 / 2.51 with
8 / 16 / 32; knot 8192 9.26 → 7.88 (16); Freedman 2000 0.756 → 0.567 (16). Default
`max_threads = Some(16)` (15–25 % faster at 64 threads, the same at ≤ 16).

**Auto crossover** (exact + dense against Barnes–Hut + multigrid, Freedman resampled, 10 steps,
`fig20-cross.csv` before the early abort, `fig20-cross2.csv` after; seconds per iteration):

| n | 64 thr: exact + dense | BH + MG | 8 thr: exact + dense | BH + MG |
|---:|---:|---:|---:|---:|
| 1500 | 0.196 | 0.534 | 0.335 | 0.576 |
| 2000 | 0.379 | 0.609 | 0.620 | 0.636 |
| 2500 | 0.716 | 0.782 | 1.106 | 0.856 |
| 3000 | 1.037 | 1.013 | 1.613 | 1.069 |
| 4000 | 1.62 | 1.07 | 2.78 | 1.15 |
| 6000 | 3.76 | 1.48 | 6.38 | 1.57 |
| 8000 | 7.28 | 1.75 | 12.1 | 1.83 |

(2000–3000 from the second pass; the others from the first, which had no early abort: the
multigrid side would be a few per cent faster.) The curves cross at about 2000 vertices at 8
threads and 3000 at 64 (§17k: 1300 and 1600, before the Schur solve). **`Backend::Auto`'s
default is now 2500** (`DEFAULT_AUTO_DENSE_MAX_VERTICES`, rustdoc, the Python default and
`.pyi`, `python-tests/test_repulsive.py`, the Blender tooltip, `tests/session_backends.rs`
checks 2500/2501): at 2500 exact + dense is 1.29× slower at 8 threads and multigrid 1.09× at
64; within 2000…3000 neither side loses more than 30 %.

**What is still open.** (a) The first steps of the knot at 8192 are above the paper within the
240 s budget (28 s, 12 %), below over 20 steps. (b) The sequential parts of the block cluster
tree product (`src/spatial/bct.rs`: upward/downward passes; `src/metric/hierarchical.rs`: `ψ`
and output loops) cap the multigrid's thread scaling at about 2–3.5×; parallelizing them needs
work outside this package. (c) A far field of `K⁰` that is exact in the tangents (second
moments of the cluster tangents, a few extra numbers per node) would let the trees keep
`theta_t` without the knot's dense trees; also `src/spatial`. (d) A dedicated formula id for the
second-order correction and the Newton polish (both tagged `rc-chord-projection` for now) needs
an entry in `lean/proofs/CITED.md`.

## 17m. Soft Hs and the stop test of the Sec. 7 comparison (WP12b)

Two findings came out of the evaluation probe: trefoil0022, 1857 vertices, (α, β) = (2, 4.5),
`Absolute(1e-4)`, 1 thread, 24 minutes, `Backend::Auto` (Barnes–Hut + multigrid).

* Soft Hs NCG with `w_H = 1e4` reached only Ê 1320. Hard Hs reached 10.7 after 100 iterations.
* Hard Hs never met `‖g̃‖₂ < 1e-4`. It was still at 0.28 after 2316 iterations.

Probes: `tests/wp12b_probe.rs` (`flow_probe`, `multipliers_probe`, `small_hard_tail`). Test:
`tests/soft_problem.rs`. Logs, curves and scripts are in `runs_tmp/rc-wp12b/`:

* `A/`: soft runs at several weights;
* `B/`: back ends compared from the plateau curve;
* `C/`: full hard runs;
* `T/`: tension calibration.

Unless a row says otherwise, runs use 1 thread.

### Finding A: the Hencky weight, not Hs, made soft Hs slow

**What the paper says.**

* Supp. 2.2 replaces the edge lengths "with Hencky strain terms in the objective for each edge,
  weighted by edge lengths of the initial configuration". It gives no further weight and no
  normalization by edge count or by total length.
* Supp. 3.2 uses an inner product `A` in NCG by "evaluating all gradients and inner products in
  the NCG algorithm using the new inner product". The metric is the preconditioner, with nothing
  added for the penalty.
* Supp. 2.2 expects soft methods to do better than hard ones: "allowing trajectories to leave
  the constraint manifold can only result in more direct paths to minimizers".
* Supp. 4.2 reports Hs NCG and Hs L-BFGS "roughly equally well" (stuck 3.3 % and 14.8 %). H¹ NCG
  and H¹ L-BFGS take about twice as long and get stuck more than 50 % of the time.
* The authors' program has no soft path: there is no Hencky term, L-BFGS or NCG in `src/`.

Our objective and our dynamic Hs preconditioner are the paper's. `w_H = 1e4` came from our own
rule (§14.2, D21), not from the paper.

**What we checked.**

* *The Hencky gradient.* Central differences at the plateau curve of trefoil0022 agree to
  1.5e-6 relative. The existing small-curve FD test passes. No bug.
* *Scale (edge tensions).* At a hard minimizer, Eq. 22 gives `dÊ = −Σ λ_I ∇ℓ_I`. The soft
  problem is stationary where `2w ln(ℓ_I/ℓ⁰_I) ≈ λ_I`, so its equilibrium strain is `λ_I / (2w)`.
  * The Euler identity (§3.3) gives `Σ λ_I ℓ_I = ½ Ê` for (2, 4.5), so `λ̄ ≈ Ê/(2L)`.
  * Trefoil0022 plateau (Ê 10.33, L 898): `λ_I` lies in 5.6e-3 … 6.0e-3, which predicts a strain
    of 3.0e-3 / w.
  * Measured after the 24-minute soft Hs NCG runs: 3.12e-3 at w = 1 and 3.1e-4 at w = 10.
  * Test `soft_equilibrium_strain_is_tension_over_two_w` (40-vertex trefoil) checks `Σλℓ = ½Ê`
    and the strain rule to O(1/w), at w = 10, 30, 100.
* *Conditioning (the cause).* The Hencky Hessian is `2w Σ (δℓ_I)²/ℓ⁰_I`, a second-order operator
  on the tangential components, scaled by `2w`. It is not in the Hs metric, which has order
  2s = 3.5.
  * At large `w` the stiff, low-frequency stretching modes dominate the Hs-preconditioned problem,
    and the descent is cut down to the steps those modes allow.
  * H¹ matches the order of a dominant Hencky term. That is why soft H¹ NCG beat soft Hs NCG at
    w = 1e4 (E 220 against 1320). This is an artifact of the weight.
  * On the 40-vertex trefoil, the iteration count of soft Hs NCG near the minimizer depends on the
    numerical path. At w = 100 it was 413, 2702 or 3842, depending on the dense solver and the
    starting curve, so no iteration count is asserted.
* *The weight on trefoil0022* (soft Hs NCG, auto back end, 24 min):

  | w_H | Ê at the end | final max \|ln ℓ/ℓ⁰\| | to 1.1 × 10.334 | to 1.01 × | largest strain in the first 300 s |
  |---|---|---|---|---|---|
  | 1 | 10.318 | 3.1e-3 | 49 s (it 60) | 80 s (it 142) | 1.77 (it 6) |
  | 10 | 10.333 | 3.1e-4 | 197 s (it 320) | 400 s (it 921) | 1.44 (it 8) |
  | 100 | 10.661 | 1.2e-4 | 932 s (it 1724) | not reached | 1.12 (it 4) |
  | 1e4 (eval probe) | 1320 | 1.5e-3 | not reached | not reached | 0.096 (it 24) |

  The transient column is the largest `max_I |ln ℓ/ℓ⁰|` in the first 300 s, measured with
  `flow_probe` on the branch head. Hard Hs PG reaches 1.1 × at 284 s in the eval probe and at
  148 s on the branch head. On the branch head, at w = 10 (`A/head.*`):

  | method | to 1.1 × | to 1.01 × | outcome |
  |---|---|---|---|
  | Hs NCG | 167 s | 356 s | |
  | Hs L-BFGS | 143 s | 195 s | |
  | H¹ NCG | | | `Stuck` at iteration 24, Ê 2999 |
  | H¹ L-BFGS | | | `Stuck` at iteration 29, Ê 2900 |

  At w = 1 (older build): Hs L-BFGS took 90 s and 117 s. H¹ NCG and H¹ L-BFGS stayed above
  Ê 2100, and both final curves failed validation (touching edges, distance 0 and 3e-15; built
  before WP14b's contact floor). With `w ≤ 10` the paper's soft ranking holds: Hs NCG ≈ Hs
  L-BFGS, H¹ slower or stuck, soft at least as fast as hard.
* *The NCG formula and the line search.* PR+ with restarts and the guarded grow-previous search
  reach the hard minimum's energy within 1 % in 80 s at w = 1. They are not the bottleneck.
  Nothing changed.

**Why the old rule picked 1e4.** It measured the strain *after 50 iterations*. At that point the
curve is still untangling: Ê is in the hundreds or thousands, and the strain follows the
current, large tension `Ê/(2L)`. A stiffer weight also slows the descent and keeps Ê high. So
the rule measured transient stiffness and rewarded the weights that stall the flow. The median
Ê after 50 iterations was 91 at w = 10 and 1926 at 1e4.

**Decision (for the evaluation; no library change).**

1. **The calibration measures the strain at the soft minimizer, not after 50 iterations.**
   * Cheap form: the edge tensions of the hard reference runs. Take the smallest decade with
     median over the 10 calibration knots of `max_I λ_I/(2w) ≤ 1e-2`.
   * `λ` comes from `multipliers_probe` at the final hard Hs curve (exact energy, one dense
     Eq. 22 solve per knot).
   * Calibration knots, hard Hs on the branch head, auto back end, 600 s:

     | knot | Ê | λ range | `max λ/2` (strain at w = 1) |
     |---|---|---|---|
     | knot0074 | 54.2 | 2.9e-2 … 3.5e-2 | 1.75e-2 |
     | knot0094 | 55.9 | 2.8e-2 … 3.6e-2 | 1.82e-2 |
     | knot0114 | 65.3 | 1.6e-2 … 9.8e-2 | 4.9e-2 (not converged, dE·g̃ 0.34) |
     | knot0125 | 65.5 | 3.9e-2 … 5.3e-2 | 2.67e-2 |
     | trefoil0013 | 11.0 | 6.9e-3 … 7.6e-3 | 3.8e-3 |
     | trefoil0027 | 10.2 | 5.5e-3 … 5.9e-3 | 2.9e-3 |

   * Median: 1.8e-2 at w = 1, 1.8e-3 at w = 10. The largest is 4.9e-3 at w = 10.
   * knot0023, knot0086, knot0102 and trefoil0046 were stopped for WP14b's timing block. Three
     are Knot128 knots (Ê ≈ 55–65, strain at w = 1 above 1e-2), so they cannot move the median
     below 1e-2 at w = 1.
   * **Result: `w_H = 10`.**
   * The literal reading of Supp. 2.2 is w = 1 (initial-length weights only). It meets the target
     on the trefoils but not on the Knot128 knots (up to 5e-2).
   * The driver's default `--hencky-w 1e3` should become 10.
2. **Record the transient strain** (largest `max_I |ln ℓ/ℓ⁰|` over the run) next to the final
   one.
   * At every `w ≤ 100` the first steps stretch the edges near close contacts by `e^1` or more,
     then relax.
   * So the soft trajectories leave the hard constraint set by far more than the final strain
     shows. Only `w ≈ 1e4` keeps the transient below 0.1, and it makes the Hs-preconditioned
     problem stiff.
   * The soft group is compared only within itself (as in the paper), so this is acceptable.
     Whether the transient matters is a question for the user.
3. **No metric change.** Adding the Hencky stiffness to the Hs preconditioner might make a
   better soft method, but the paper's soft methods do not do it.

### Finding B: why hard Hs never met `‖g̃‖₂ < 1e-4`

Three causes, measured separately, plus the norm question.

**(d) A bug: the flow stalled at the projection-tolerance boundary.**

*Mechanism.*

* `project_constraints` and the line search's `project_with_path_polished` tested the tolerance
  before the first chord step. A trial already within tolerance was not projected at all.
* A step along the tangent of the edge-length constraint lengthens every edge at second order
  (`ℓ(γ − τg̃)² = ℓ² + τ²|Δg̃|²`).
* With positive tensions, stretching lowers Ê. The curve drifted, unprojected and "for free",
  until `|Φ| = tol`.
* From there, every trial that left the tolerance band was projected back and paid
  `≈ Σ λ_I δℓ_I`.
* On trefoil0022 that cost is `5.7e-3 · 1857 · O(1e-4) ≈ 1e-3`. It is far larger than the
  Armijo decrease `c τ dE·g̃` (about 1e-10).
* So the search halved `τ` to about 1e-5. Only steps that stay inside the band were accepted.

*Measurement* (`small_hard_tail`: 40-vertex trefoil, exact + dense, (2, 4.5), 3000 iterations):

| projection tolerance (× ℓ̄₀) | ‖g̃‖ floor | base violation | unprojected accepted steps |
|---|---|---|---|
| 1e-4 | 4.9e-3 (1.1e-2 with the Newton polish of §17l on) | = tol | 2983 / 3000 |
| 1e-6 | 1.0e-3 | = tol | 2892 / 3000 |
| 1e-8 | 9.1e-5 | = tol | 2573 / 3000 |
| 1e-10 | 7.4e-6 | = tol | 2181 / 3000 |
| 1e-4, at least one projection step | 1.6e-5 at iteration 1000, linear decrease, τ ≈ 0.6 | 2e-16 | 0 |

*What the paper and the reference do.*

* The paper projects by "update our guess via γ̃ ← γ̃ + x and repeat until the constraint
  violation is numerically small" (Sec. 5.3.2).
* The reference's `LSBackproject` always calls `BackprojectConstraints` once before it tests
  (`ref: src/tpe_flow_sc.cpp`, l. 306–313).

*Fix.* At least one projection step per trial, unless `Φ = 0` exactly.

* Applied by WP14b in `src/solve/mod.rs` and `src/step/line_search.rs` (§6.2, §17l, commit
  0e7848e).
* Verified first on a patched copy: last row of the table.
* Cost: one extra chord or Newton step for each trial that was within tolerance.

**(a) The Barnes–Hut energy makes a noise floor of `‖g̃‖`; the multigrid does not.** Four back
ends, 150 steps each, from the saved plateau curve of the eval's hs-pg run (`B/`, fresh length
targets):

| back end | ‖g̃‖, step 1 → 150 | F over the run |
|---|---|---|
| exact + dense | 1.32 → 0.024, monotone | 10.334368 → 10.334365, monotone |
| exact + multigrid (tol 1e-3) | 1.32 → 0.025, monotone (within 4 % of dense) | same |
| Barnes–Hut + dense | 0.07 … 0.29, no trend | jitters by ±2e-6 |
| Barnes–Hut + multigrid | 0.08 … 0.28, no trend | jitters by ±2e-6 (BH energy 0.11 % above exact) |

* *Multigrid: no floor.* Its tolerance is relative to the projected right side `P dE`, which
  vanishes at the minimum.
* *Approximate metric (BCT): no floor.* An approximate metric only rescales `g̃`.
* *Barnes–Hut: the floor.*
  * The interaction lists are rebuilt every step (§17f, §17i), so the evaluated function changes
    by about the approximation error from step to step.
  * The Euclidean `‖g̃‖` amplifies that change, because it is dominated by low-frequency modes
    with tiny `Ā` eigenvalues (see (c)).
  * The 0.1–0.3 floor of the eval's run is this one.
* From the start, Barnes–Hut + multigrid with the fix of (d) (patched copy) is at `‖g̃‖₂` 2.0e-2
  after 1680 iterations. Without the fix it is at 5.7e-2 after 3102 iterations.

**(c) A genuinely slow tail in the Euclidean norm.**

* *Rate.* Even exact, dense and fixed, `‖g̃‖₂` falls by only about 1.1 % per step on the plateau
  (`τ` 0.18–0.36, Armijo-limited).
* *Time to 1e-4.* From 0.024 at step 150, it needs about 500 more steps to reach 1e-4: about 25
  minutes at 3 s per exact dense step. Meanwhile the energy is flat to 3e-7 relative.
* *Why so slow.* At the plateau curve, `|dE| = 3.0e-3` and `‖g̃‖₂ = 1.32`, but
  `dE·g̃ = g̃ᵀĀg̃ = 1.9e-5`. The Rayleigh quotient of `g̃` is 1e-5, so its Euclidean norm is
  carried by modes the energy barely sees.
* *Full runs from the start* (trefoil0022, 24 min, `C/`):
  * Every back end first has `dE·g̃ < 1e-4` at iteration 412–420.
  * Barnes–Hut + multigrid gets there at 318–388 s (E 10.3458, its plateau value).
  * Exact + dense gets there at about 1370 s (E 10.3350).
  * No configuration reaches `‖g̃‖₂ < 1e-4` within 24 minutes.

**(b) The norm.**

* *Which vector.* `gradient_norm` is the Euclidean norm of the projected `g̃` itself. `g̃` lies in
  `ker C`, so no constraint-violating component is included (the dense saddle residual is
  1e-12).
* *Which norm.*
  * `gradient_mass_norm` is `√(Σ m_v |g̃_v|²)`: the discrete L² norm of the function,
    independent of resolution.
  * The Euclidean norm grows like `√|V|` for a fixed field.
  * Here the two differ by a factor of 0.71 (0.94 against 1.32), which does not change the
    conclusion.
* *Scale.* `g̃` scales as `c¹` under `γ → cγ` for (2, 4.5) (§4.5). These knots are 30 units wide
  with L ≈ 900, so the threshold 1e-4 is 3e-6 of the curve's size.
* *Supplemental.* It says "the L² norm of the preconditioned gradient is below a fixed tolerance
  ε = 10⁻⁴". Supp. 1.3 says the time cap exists because some methods "would require days to
  satisfy the first condition".
* *Reference.* The authors' program stops on `dE·g̃ < 1e-4` and prints that quantity as
  "Sobolev gradient norm" (`ref: src/tpe_flow_sc.cpp`, `soboNormZero = (soboDot < 1e-4)`).

**Decision (for the evaluation).**

* *Stop rule.* Use **`dE·g̃ < 1e-4`** as the evaluation's convergence test. It is the reference
  program's rule, and the reference prints the quantity as "Sobolev gradient norm".
  * The Sec. 7 measures do not depend on the stop rule, only on the energy history: time and
    iterations to 1.1 × the reference's minimal energy, and stuck or nonconvergent.
  * The stop only ends runs early and, for the reference method, fixes `E_ref`.
  * Hard methods: `StepReport::sobolev_norm_sq`.
  * Soft methods: the same field holds `dF·p` for the NCG or L-BFGS direction, which equals
    `dF·z` at restarts.
  * Keep `‖g̃‖₂` and `‖g̃‖_{L²}` in the per-step history.
* *Effect on trefoil0022.* The rule ends hard Hs at about 318 s, at its plateau energy, instead
  of running the full 24 minutes. That is about 4.5× less compute for the reference method on
  this knot.
* *Library default unchanged.* `‖g̃‖₂ < 1e-4` stays the library default (D14). It is an absolute,
  scale-dependent rule, and unit-size curves reach it: the 40-vertex trefoil does within 1000
  iterations with the fix.
* *Where the evaluation's stop goes.* Either the driver (`examples/eval_knots.rs`: end the loop
  when `r.sobolev_norm_sq < 1e-4` after an accepted step), or a `StopRule` variant in
  `src/step/stop.rs`. Neither was changed here, because both have other owners.
  `Session::step` does not apply `FlowOptions::time_limit` (only `run` does). The driver has its
  own `--max-seconds`.
* *Back end.* Barnes–Hut + multigrid departs from the paper's Sec. 7 setup (exact derivatives,
  dense solves) on these 1500–3000-vertex knots.
  * It moves energies by 0.1 % and adds the floor of (a).
  * It reaches `dE·g̃ < 1e-4` at the same iteration as exact + dense.
  * Exact + dense costs about 3 s per step against 0.4–0.8 s, so Auto stays the budget choice.
  * The driver re-evaluates final energies exactly anyway.

## 18. The user's decisions (2026-10-26)

* Q1: the paper's edge-pair energy only (D1); the reference's vertex-pair variant is not
  implemented. The comparison with the authors' program uses what is comparable (Gram matrix,
  Jacobians, coarsening, final curves); energies are reported side by side, not matched.
* Q2: Blender writes Bézier/NURBS Curve objects back as poly splines (exact), with an optional,
  approximate "Keep Bézier" re-fit; pins on Curve objects are the selected control points.
* Q3: the live operator can subdivide while it runs, as an option that is off by default.
* Q4: evaluation runs all 228 knots in parallel with iteration caps, then wall-clock timings
  alone on a 24-knot subset.
* Q5 (after WP3c): the edge-collapse guard of §17j stays on by default
  (`min_edge_length_ratio = Some(1e-3)`).
* Q6: the timed evaluation (WP12 phase B) records curve frames for every algorithm on the 24-knot
  subset, at equal wall-clock intervals, in a compact format, and renders a side-by-side grid video
  (one column per algorithm), as in the authors' video. Phase A records energies and final curves only.
* Q7 (after WP12 part 1): phase A uses the paper's budget, 24 minutes of wall time per run at one
  thread (60 runs in parallel), no iteration cap below it. Phase B runs the 24-knot subset alone at
  64 threads (reference methods 300 s, the others 120 s, frames every second), plus the two
  reference methods at one thread on 6 of those knots; the other one-thread times come from phase A.
* Q8 (after WP12b, §17m Finding A): the soft problems of the evaluation use the Hencky weight
  `w_H = 10`, from the tension rule (median over knots of the maximum strain `λ/(2w)` at most 1e-2),
  not the earlier 1e4.
* Q9 (after WP12b, §17m Finding B): the evaluation stops a run on the energy rate,
  `dE·g̃ < 1e-4` (`StepReport::sobolev_norm_sq`, the reference's rule), or on its wall budget.
  The library's default stop rule (`‖g̃‖₂ < 1e-4`, D14) is unchanged.
