import DdgOracle.Formulas.Sum
import DdgOracle.Formulas.Threshold

/-!
# Optimality system, residual maps and semismooth differentials

The canonical discrete problem (`CLAUDE.md`), restricted to the `n` interior vertices:

  minimize   `½ Σ wE_i m_i u_i² + Σ t_i |µ_i|`
  subject to `L u = µ − Ω`,  `lo_i ≤ µ_i ≤ hi_i`,          (`t_i = λ·wR_i`)

with `L = L_II` the interior block of the positive semidefinite cotan Laplacian.
With the multiplier `φ` of the constraint `L u − µ + Ω = 0` the Lagrangian is
`½ uᵀ W u + g(µ) + φᵀ (L u − µ + Ω)` (`W = diag(wE_i m_i)`), and its stationarity gives

  `L u = µ − Ω`        (state equation; `u` is the true log scale factor)
  `L φ = −W u`         (adjoint equation)
  `φ_i ∈ ∂g_i(µ_i)`    (inclusion; unbounded: `|φ_i| ≤ t_i`, `φ_i = t_i sign µ_i` where `µ_i ≠ 0`)

This `φ` equals the paper's adjoint; the paper's `u` is the negative of ours, which is why the
paper prints `L u = Ω − µ`, `L φ = M u` (Eq. 17). The paper's Eq. 20 (`L u − Ω + µ`) is this
system in the paper's sign; App. B's `L u − Ω − µ` is a misprint.

Vectors are functions `Nat → α` on the indices `0, …, n−1`, the matrix `L` is a function
`Nat → Nat → α`. Residuals and differentials are given both per block and flattened with the
unknowns ordered `(u, φ, µ)` and the rows ordered `(state, adjoint, cone)`.
-/

namespace Ddg

variable {α : Type} [Scalar α]

/-- Data of the canonical problem on the `n` interior vertices. -/
structure Problem (α : Type) where
  /-- number of interior vertices -/
  n : Nat
  /-- `L_II`, interior block of the PSD cotan Laplacian -/
  L : Nat → Nat → α
  /-- lumped masses `m_i` (normalized by the total area) -/
  mass : Nat → α
  /-- distortion weights `wE_i ≥ 0` -/
  wE : Nat → α
  /-- angle defects `Ω_i` -/
  omega : Nat → α
  /-- thresholds `t_i = λ·wR_i ≥ 0` -/
  t : Nat → α
  /-- lower bounds on the cone angles -/
  lo : Nat → α
  /-- upper bounds on the cone angles -/
  hi : Nat → α

namespace Problem

/-- `(L v)_i`. -/
def applyL (p : Problem α) (v : Nat → α) (i : Nat) : α :=
  sumTo p.n fun j => p.L i j * v j

/-! ## Residuals -/

/-- State residual `(L u)_i − µ_i + Ω_i`. -/
def stateResidual (p : Problem α) (u mu : Nat → α) (i : Nat) : α :=
  p.applyL u i - mu i + p.omega i

/-- Adjoint residual `(L φ)_i + wE_i m_i u_i`. -/
def adjointResidual (p : Problem α) (u phi : Nat → α) (i : Nat) : α :=
  p.applyL phi i + p.wE i * p.mass i * u i

/-- Cone residual `µ_i − P(φ_i + µ_i)` with `P` the (bounded) soft threshold. -/
def coneResidual (p : Problem α) (phi mu : Nat → α) (i : Nat) : α :=
  mu i - boundedThreshold (p.t i) (p.lo i) (p.hi i) (phi i + mu i)

/-- The regularized measure `µ_γ(φ)_i = clamp((m_i/γ) · P_{t_i}(φ_i))` of Stage I.
The paper prints `µ = P_λ(φ)/γ` (Eq. 21), which lacks the mass factor: `µ_i` is an integrated
quantity while the Moreau–Yosida term `(γ/2)‖µ‖²_{L²}` acts on the density `µ_i / m_i`. -/
def regMeasure (p : Problem α) (gamma : α) (phi : Nat → α) (i : Nat) : α :=
  regularizedMeasure (p.mass i / gamma) (p.t i) (p.lo i) (p.hi i) (phi i)

/-- State residual of the regularized problem, `(L u)_i − µ_γ(φ)_i + Ω_i`. -/
def regStateResidual (p : Problem α) (gamma : α) (u phi : Nat → α) (i : Nat) : α :=
  p.applyL u i - p.regMeasure gamma phi i + p.omega i

/-- The residual map `F(u, φ, µ) ∈ α^{3n}` of the exact problem (paper Eq. 20, canonical signs):
rows `[0, n)` state, `[n, 2n)` adjoint, `[2n, 3n)` cone. -/
def F (p : Problem α) (u phi mu : Nat → α) (r : Nat) : α :=
  if r < p.n then p.stateResidual u mu r
  else if r < 2 * p.n then p.adjointResidual u phi (r - p.n)
  else p.coneResidual phi mu (r - 2 * p.n)

/-- The residual map `F_γ(u, φ) ∈ α^{2n}` of the regularized problem (paper Eq. 22, canonical
signs and with the mass factor): rows `[0, n)` state, `[n, 2n)` adjoint. -/
def Freg (p : Problem α) (gamma : α) (u phi : Nat → α) (r : Nat) : α :=
  if r < p.n then p.regStateResidual gamma u phi r
  else p.adjointResidual u phi (r - p.n)

/-! ## Active sets and semismooth differentials -/

/-- State of vertex `i` in Stage II: the piece of `P` at `φ_i + µ_i` (paper: `A^±(φ + µ)`). -/
def state (p : Problem α) (phi mu : Nat → α) (i : Nat) : ConeState :=
  coneState 1 (p.t i) (p.lo i) (p.hi i) (phi i + mu i)

/-- State of vertex `i` in Stage I: the piece of `µ_γ` at `φ_i` (paper: `A^±(φ)`). -/
def regState (p : Problem α) (gamma : α) (phi : Nat → α) (i : Nat) : ConeState :=
  coneState (p.mass i / gamma) (p.t i) (p.lo i) (p.hi i) (phi i)

/-- Entry `(r, c)` of the semismooth differential `DF(u, φ, µ)` (paper Eq. 26, canonical signs),
a `3n × 3n` matrix with columns ordered `(u, φ, µ)`:
```
  [ L     0      −I   ]
  [ W     L       0   ]      W = diag(wE_i m_i),  D = diag(indicator of A⁺ ∪ A⁻)
  [ 0    −D     I − D ]
``` -/
def DF (p : Problem α) (phi mu : Nat → α) (r c : Nat) : α :=
  if r < p.n then
    if c < p.n then p.L r c
    else if c < 2 * p.n then 0
    else if c - 2 * p.n = r then -1 else 0
  else if r < 2 * p.n then
    if c < p.n then (if c = r - p.n then p.wE (r - p.n) * p.mass (r - p.n) else 0)
    else if c < 2 * p.n then p.L (r - p.n) (c - p.n)
    else 0
  else
    if c < p.n then 0
    else if c < 2 * p.n then
      (if c - p.n = r - 2 * p.n then -(p.state phi mu (r - 2 * p.n)).indicator else 0)
    else
      (if c - 2 * p.n = r - 2 * p.n then 1 - (p.state phi mu (r - 2 * p.n)).indicator else 0)

/-- Entry `(r, c)` of the semismooth differential `DF_γ(u, φ)` (paper Eq. 25, canonical signs and
with the mass factor), a `2n × 2n` matrix with columns ordered `(u, φ)`:
```
  [ L    −(1/γ) M D ]
  [ W        L      ]
``` -/
def DFreg (p : Problem α) (gamma : α) (phi : Nat → α) (r c : Nat) : α :=
  if r < p.n then
    if c < p.n then p.L r c
    else if c - p.n = r then -(p.mass r / gamma * (p.regState gamma phi r).indicator) else 0
  else
    if c < p.n then (if c = r - p.n then p.wE (r - p.n) * p.mass (r - p.n) else 0)
    else p.L (r - p.n) (c - p.n)

/-! ## The regularized problem with `µ` kept as an unknown

The paper eliminates `µ = µ_γ(φ)` (Eq. 21) and applies Newton's method to `F_γ(u, φ)`. That form
has a defect once cone angles are bounded: `µ_γ` rises from `0` to `hi` over the interval
`t ≤ φ ≤ t + γ·hi/m`, which shrinks to a point as `γ → 0`, so a vertex jumps between "no cone"
and "cone at the bound" without the linearization ever seeing the values in between, and the
iteration falls into two-cycles. The optimality condition of the regularized problem,
`φ − (γ/m) µ ∈ ∂g(µ)`, is a fixed-point equation of a proximal map for every step size; with
the step size `c = m/(m + γ)` it reads

  `µ = clamp(c · P_t(φ + µ))`,   `c = m / (m + γ) ∈ (0, 1]`,

where the linear piece has the width `hi/c ≥ hi` for every `γ`. On that piece it says
`µ = (m/γ)(φ − t)`, the same relation as `µ = µ_γ(φ)`. For `γ = 0` it is the exact cone
condition, so the exact system `F` is the member `γ = 0` of the family `Fimp γ`.
-/

/-- `c_i = m_i / (m_i + γ)`. -/
def coneScale (p : Problem α) (gamma : α) (i : Nat) : α :=
  p.mass i / (p.mass i + gamma)

/-- Cone residual of the regularized problem, `µ_i − clamp(c_i · P_{t_i}(φ_i + µ_i))`. -/
def regConeResidual (p : Problem α) (gamma : α) (phi mu : Nat → α) (i : Nat) : α :=
  mu i - regularizedMeasure (p.coneScale gamma i) (p.t i) (p.lo i) (p.hi i) (phi i + mu i)

/-- State of vertex `i` for the regularized problem with `µ` kept: the piece of
`x ↦ clamp(c_i · P_{t_i}(x))` at `x = φ_i + µ_i`. -/
def impState (p : Problem α) (gamma : α) (phi mu : Nat → α) (i : Nat) : ConeState :=
  coneState (p.coneScale gamma i) (p.t i) (p.lo i) (p.hi i) (phi i + mu i)

/-- Residual map `F^γ(u, φ, µ) ∈ α^{3n}` of the regularized problem with `µ` kept as an unknown:
rows `[0, n)` state, `[n, 2n)` adjoint, `[2n, 3n)` cone. Its zeros are the zeros of `Freg`
together with `µ = µ_γ(φ)`, and `Fimp 0 = F`. -/
def Fimp (p : Problem α) (gamma : α) (u phi mu : Nat → α) (r : Nat) : α :=
  if r < p.n then p.stateResidual u mu r
  else if r < 2 * p.n then p.adjointResidual u phi (r - p.n)
  else p.regConeResidual gamma phi mu (r - 2 * p.n)

/-! ## The linear system of a given choice of pieces

`Fimp γ` is piecewise affine. Given a piece `σ_i` for every vertex, `Fpiece γ σ` is the affine map
that agrees with `Fimp γ` wherever the point lies on those pieces, and `DFpiece γ σ` is its
(constant) Jacobian. A Newton step with these is the exact solution of the linear optimality
system of the pieces `σ`. With `σ = impState` (the pieces the current point lies on) this is the
semismooth Newton step of the paper; the executable passes the safeguarded states
(`nextState`), which differ from `impState` only when a vertex would skip a piece.
-/

/-- Cone residual on the piece `σ_i`: `µ_i − (value of the piece at φ_i + µ_i)`. -/
def pieceConeResidual (p : Problem α) (gamma : α) (sigma : Nat → ConeState) (phi mu : Nat → α)
    (i : Nat) : α :=
  mu i - (sigma i).value (p.coneScale gamma i) (p.t i) (p.lo i) (p.hi i) (phi i + mu i)

/-- The affine map that agrees with `Fimp γ` on the pieces `σ`. -/
def Fpiece (p : Problem α) (gamma : α) (sigma : Nat → ConeState) (u phi mu : Nat → α)
    (r : Nat) : α :=
  if r < p.n then p.stateResidual u mu r
  else if r < 2 * p.n then p.adjointResidual u phi (r - p.n)
  else p.pieceConeResidual gamma sigma phi mu (r - 2 * p.n)

/-- Entry `(r, c)` of the Jacobian of `Fpiece γ σ`, columns ordered `(u, φ, µ)`:
```
  [ L     0        −I    ]
  [ W     L         0    ]      C = diag(c_i),  D = diag(indicator of σ_i ∈ {pos, neg})
  [ 0   −C D     I − C D ]
``` -/
def DFpiece (p : Problem α) (gamma : α) (sigma : Nat → ConeState) (r c : Nat) : α :=
  if r < p.n then
    if c < p.n then p.L r c
    else if c < 2 * p.n then 0
    else if c - 2 * p.n = r then -1 else 0
  else if r < 2 * p.n then
    if c < p.n then (if c = r - p.n then p.wE (r - p.n) * p.mass (r - p.n) else 0)
    else if c < 2 * p.n then p.L (r - p.n) (c - p.n)
    else 0
  else
    if c < p.n then 0
    else if c < 2 * p.n then
      (if c - p.n = r - 2 * p.n then
        -(p.coneScale gamma (r - 2 * p.n) * (sigma (r - 2 * p.n)).indicator)
       else 0)
    else
      (if c - 2 * p.n = r - 2 * p.n then
        1 - p.coneScale gamma (r - 2 * p.n) * (sigma (r - 2 * p.n)).indicator
       else 0)

/-- The semismooth differential of `Fimp γ`: the Jacobian of the pieces the point lies on.
For `γ = 0` it is `DF` (paper Eq. 26). -/
def DFimp (p : Problem α) (gamma : α) (phi mu : Nat → α) (r c : Nat) : α :=
  p.DFpiece gamma (p.impState gamma phi mu) r c

/-! ## Objective -/

/-- Distortion energy `½ Σ wE_i m_i u_i²`. -/
def distortion (p : Problem α) (u : Nat → α) : α :=
  (sumTo p.n fun i => p.wE i * p.mass i * (u i * u i)) / 2

/-- Cone cost `Σ t_i |µ_i| = λ Σ wR_i |µ_i|`. -/
def coneCost (p : Problem α) (mu : Nat → α) : α :=
  sumTo p.n fun i => p.t i * Scalar.abs (mu i)

/-- Objective of the canonical problem. -/
def objective (p : Problem α) (u mu : Nat → α) : α :=
  p.distortion u + p.coneCost mu

end Problem

end Ddg
