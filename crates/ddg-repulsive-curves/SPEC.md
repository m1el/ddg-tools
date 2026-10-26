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
  `w⁰_IJ + w⁰_JI > 0`; `ker B⁰` contains the global constants. `ker A = ker B ∩ ker B⁰`.
* For a generic network `ker A = span{𝟙}`: `B⁰` couples all pairs of edges, also across
  components. *(checked: two loops, dim ker B = 2, dim ker B⁰ = dim ker A = 1.)* So `Ā` has
  exactly the three translations in its kernel (Sec. 4.2.4); a constraint must remove them
  (§5.2, D5).
* Degenerate: if all edges are collinear, `k² ≡ 0`, `B⁰ = 0` and `ker A` is larger. The solver
  reports `SingularMetric` (the energy and its differential are 0 there anyway).

Scaling [formula: rc-gram-scaling]: `w(cγ) = c^{1−2σ} w(γ)`, `w⁰(cγ) = c^{1−2σ} w⁰(γ)`, hence
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
when every row is within its tolerance (§5.1). Failure after `max_projection_iterations` (3; the
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
Taylor argument, cited.

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
`Backend::Auto { dense_max_vertices }` picks exact + dense below the crossover measured in §14.1
(placeholder 2000 vertices) and BH + multigrid above. Dense solve: one LU with partial pivoting
(`faer`, multithreaded) of the `(3|V| + k)` saddle matrix, reused for Eq. 22 and every chord
iteration. The L², H¹, H² metrics use a sparse LU of their saddle matrices (`ddg_core::sparse`).

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
                       pub smoothing_steps: usize, pub tolerance: f64, pub max_cycles: usize }
pub enum EnergyEval { Exact, BarnesHut(BhOptions) }
pub enum Solver { Dense, Multigrid(MgOptions) }
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
    CollisionGuarded { armijo: f64, shrink: f64, start_fraction: f64, gap: f64, max_halvings: u32 },
    NormalizedArmijo { armijo: f64, shrink: f64, max_halvings: u32 },
}
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
    pub seconds: Timings,                              // energy, differential, metric, solve, line search, projection, ccd
}
pub struct FlowReport { pub status: Status, pub steps: Vec<StepReport>, pub final_energy: f64 }
```

**Errors** (`RcError`, `thiserror`-free hand-written `Display` like the other crates): the
validation errors of §2.3; `InvalidExponents`, `InvalidOption { name, value, requirement }`,
`NotPlanar { vertex, axis, value }`; constraint errors `TranslationNotFixed`,
`RedundantConstraints { rows }`, `InfeasibleConstraint { constraint, violation }`,
`InvalidTangentSite { vertex, degree }`, `ConstraintOutOfRange { constraint, index }`,
`ScheduleWithoutLengthConstraint`, `MethodNeedsSoftProblem`; solve errors `SingularMetric`,
`SingularSaddle { pivot }`, `MultigridDiverged { residuals }`, `NonFiniteEnergy { step }`;
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
| line search | `CollisionGuarded { armijo 1e-4, shrink ½, start ⅔, gap 0.1, 40 halvings }` | D11 (more reliable); ⅔ is the paper's; `c = 1e-4` the usual lenient Armijo constant (Nocedal & Wright §3.1) |
| `NormalizedArmijo` | same Armijo, reference's start rule | D12 |
| projection | `MeanEdgeLength(1e-4)`, 3 iterations | paper: 1e-4 (absolute) and ≤ 3 iterations; relative to `ℓ̄₀` so that Blender scenes of any scale behave alike; the evaluation uses `Absolute(1e-4)` (D10) |
| stop | `GradientNorm { Euclidean, 1e-4 }`, `max_iterations` 1000 | paper Sec. 5.4, Supp. 1.3 (D14); Blender uses `RelativeGradient(1e-3)` |
| Barnes–Hut | `theta_x 0.125, theta_t 0.25, leaf_size 8` | reference `diag/d < 0.25` ⇔ `r_x/d < 0.125`; tangent bound per Eq. 27; final values by §14.3 (largest θ with differential error ≤ 1e-2 and unchanged iteration counts within 5 %) |
| BCT | `theta_x 0.5, theta_t 0.5, leaf_pairs 64` | reference `sep = 1.0` on `diag/(d − r)` ≈ `r/d ≲ 0.5`; measured as above |
| multigrid | coarsest 512 vertices, 20 levels max, 3+3 CG smoothing, tol 1e-3, 20 cycles | paper B.3.1 (≤ 6 V-cycles to 1e-3); 512 vertices = a 1536-row dense LU, milliseconds |
| back end | `Auto { dense_max_vertices: 2000 }` | placeholder; the crossover of §14.1 |
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
of the final largest `|ln(ℓ/ℓ⁰)|` under Hs NCG is ≤ 1e-2, fixed before the full run.

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

## 18. The user's decisions (2026-10-26)

* Q1: the paper's edge-pair energy only (D1); the reference's vertex-pair variant is not
  implemented. The comparison with the authors' program uses what is comparable (Gram matrix,
  Jacobians, coarsening, final curves); energies are reported side by side, not matched.
* Q2: Blender writes Bézier/NURBS Curve objects back as poly splines (exact), with an optional,
  approximate "Keep Bézier" re-fit; pins on Curve objects are the selected control points.
* Q3: the live operator can subdivide while it runs, as an option that is off by default.
* Q4: evaluation runs all 228 knots in parallel with iteration caps, then wall-clock timings
  alone on a 24-knot subset.
