import Mathlib

/-!
# Applications of intrinsic triangulations, Sec. 5.4 and 5.5: closed forms

`crates/ddg-intrinsic-apps/SPEC.md` (decisions N40 to N49), the paper Sec. 5.4 and 5.5.

* **`app-vector-face-holonomy`** — with the frame angles of a face's three halfedges related by
  the signpost update at each corner (`φ_ik = φ_ij + s_i θ_i` up to full turns), the transport
  angles `ρ_ij = φ_ji + π − φ_ij` (Sharp, Soliman, Crane 2019, VHM Eq. 10) of the three
  halfedges add up to `Σ (s_c − 1) θ_c` up to full turns (`transport_around_face`), and at an
  interior vertex with `s = 2π/Θ` the corner terms of all faces around it add up to the angle
  defect `2π − Θ` (`scaled_minus_sum`).
* **`app-vector-fd-connection-laplacian`** — the contribution of one edge with cotan weight
  `w` and transport `r` (`|r| = 1`) to `u* L u`, `w (|u_i|² + |u_j|²) − 2 w Re(conj u_j · r u_i)`,
  is `w |u_j − r u_i|²` (`edge_energy`), so the matrix is positive semidefinite when every
  weight is non-negative; its off-diagonal entries are Hermitian (`edge_entries_hermitian`).
* **`app-amr-map-distortion`** — for a 2×2 matrix `J = [[p, q], [r, s]]` with conformal part
  `a = ((p + s)/2, (r − q)/2)` and anticonformal part `b = ((p − s)/2, (r + q)/2)`,
  `|a|² − |b|² = det J` and `2(|a|² + |b|²) = ‖J‖²_F` (`conformal_parts_det`,
  `conformal_parts_frobenius`); hence `σ₁ = |a| + |b|`, `σ₂ = |a| − |b|` have
  `σ₁ σ₂ = det J` and `σ₁² + σ₂² = ‖J‖²_F` (`singular_values_from_parts`), the two equations
  that determine the singular values.
-/

namespace DDG

namespace Apps

open Real

/-- [formula: app-vector-face-holonomy]
The transport angles around a face add up to the curvature the tangent frames put into it.
At corner `i` the frame angles of the outgoing halfedges `ij`, `ik` differ by `s_i θ_i` (the
scaled corner angle) up to `a` full turns, and likewise at `j` (from `jk` to `ji`) and `k`
(from `ki` to `kj`); the corner angles add up to `π`. -/
theorem transport_around_face (φij φik φji φjk φki φkj si sj sk θi θj θk : ℝ) (a b c : ℤ)
    (hi : φik = φij + si * θi + 2 * π * a)
    (hj : φji = φjk + sj * θj + 2 * π * b)
    (hk : φkj = φki + sk * θk + 2 * π * c)
    (hsum : θi + θj + θk = π) :
    (φji + π - φij) + (φkj + π - φjk) + (φik + π - φki) =
      (si - 1) * θi + (sj - 1) * θj + (sk - 1) * θk + 2 * π * (a + b + c + 2) := by
  rw [hi, hj, hk]
  linear_combination hsum

/-- [formula: app-vector-face-holonomy]
At an interior vertex with angle sum `Θ > 0` and scale `s = 2π/Θ`, the terms `(s − 1) θ` of its
corners (which add up to `Θ`) add up to the angle defect. -/
theorem scaled_minus_sum (Θ : ℝ) (hΘ : Θ ≠ 0) : (2 * π / Θ - 1) * Θ = 2 * π - Θ := by
  field_simp

/-- [formula: app-vector-fd-connection-laplacian]
The energy of one edge in the finite-difference connection Laplacian. -/
theorem edge_energy (w : ℝ) (r ui uj : ℂ) (hr : Complex.normSq r = 1) :
    w * (Complex.normSq ui + Complex.normSq uj) - 2 * w * ((starRingEnd ℂ) uj * (r * ui)).re =
      w * Complex.normSq (uj - r * ui) := by
  have h1 : Complex.normSq (r * ui) = Complex.normSq ui := by
    rw [Complex.normSq_mul, hr, one_mul]
  rw [Complex.normSq_sub, h1]
  simp only [Complex.mul_re, Complex.mul_im, Complex.conj_re, Complex.conj_im]
  ring

/-- [formula: app-vector-fd-connection-laplacian]
The entries `L_ij = −w conj r` and `L_ji = −w r` of one edge are Hermitian. -/
theorem edge_entries_hermitian (w : ℝ) (r : ℂ) :
    (starRingEnd ℂ) (-(w : ℂ) * (starRingEnd ℂ) r) = -(w : ℂ) * r := by
  simp

/-- [formula: app-amr-map-distortion] -/
theorem conformal_parts_det (p q r s : ℝ) :
    (((p + s) / 2) ^ 2 + ((r - q) / 2) ^ 2) - (((p - s) / 2) ^ 2 + ((r + q) / 2) ^ 2) =
      p * s - q * r := by
  ring

/-- [formula: app-amr-map-distortion] -/
theorem conformal_parts_frobenius (p q r s : ℝ) :
    2 * ((((p + s) / 2) ^ 2 + ((r - q) / 2) ^ 2) + (((p - s) / 2) ^ 2 + ((r + q) / 2) ^ 2)) =
      p ^ 2 + q ^ 2 + r ^ 2 + s ^ 2 := by
  ring

/-- [formula: app-amr-map-distortion]
With `A = |a|²` and `B = |b|²`, `σ₁ = √A + √B` and `σ₂ = √A − √B` satisfy `σ₁ σ₂ = det J` and
`σ₁² + σ₂² = ‖J‖²_F`. -/
theorem singular_values_from_parts (p q r s : ℝ) :
    let A := ((p + s) / 2) ^ 2 + ((r - q) / 2) ^ 2
    let B := ((p - s) / 2) ^ 2 + ((r + q) / 2) ^ 2
    (√A + √B) * (√A - √B) = p * s - q * r ∧
      (√A + √B) ^ 2 + (√A - √B) ^ 2 = p ^ 2 + q ^ 2 + r ^ 2 + s ^ 2 := by
  intro A B
  have hA : 0 ≤ A := by positivity
  have hB : 0 ≤ B := by positivity
  have sA := Real.sq_sqrt hA
  have sB := Real.sq_sqrt hB
  constructor
  · have : (√A + √B) * (√A - √B) = √A ^ 2 - √B ^ 2 := by ring
    rw [this, sA, sB]
    exact conformal_parts_det p q r s
  · have : (√A + √B) ^ 2 + (√A - √B) ^ 2 = 2 * (√A ^ 2 + √B ^ 2) := by ring
    rw [this, sA, sB]
    exact conformal_parts_frobenius p q r s

end Apps

end DDG
