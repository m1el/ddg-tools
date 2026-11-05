import DDG.Geometry.Triangle

/-!
# Triangle meshes given by edge lengths

`TriMesh V F`: every face `f` has three vertices `a f`, `b f`, `c f` and the edge lengths are
given by `len`. The per-face quantities (side lengths `la lb lc`, cotangents `cotA`–`cotC`, angles
`angleA`–`angleC`, `area`), the hypothesis `TriMesh.Valid` (three distinct vertices, strict
triangle inequalities) and uniform scaling (`TriMesh.scale`, `angle*_scale`).

(This module was the first part of `DdgProofs.Cotan`; the cotan formula, the assembled cotan
Laplacian and the lumped mass are in `DdgProofs.Laplace.Cotan`.)
-/

noncomputable section

namespace DDG

/-! ## The structure -/

/-- A triangle mesh given intrinsically: every face `f` has three vertices `a f`, `b f`, `c f`,
and `len i j` is the length of the edge between vertices `i` and `j` (only evaluated on pairs of
vertices of a common face). -/
structure TriMesh (V F : Type*) where
  a : F → V
  b : F → V
  c : F → V
  len : V → V → ℝ

namespace TriMesh

variable {V F : Type*} (M : TriMesh V F)

/-- Length of the side of `f` opposite its corner `a`. -/
def la (f : F) : ℝ := M.len (M.b f) (M.c f)
/-- Length of the side of `f` opposite its corner `b`. -/
def lb (f : F) : ℝ := M.len (M.c f) (M.a f)
/-- Length of the side of `f` opposite its corner `c`. -/
def lc (f : F) : ℝ := M.len (M.a f) (M.b f)

/-- Cotangent of the corner of `f` at `a`. -/
def cotA (f : F) : ℝ := cornerCot (M.la f) (M.lb f) (M.lc f)
/-- Cotangent of the corner of `f` at `b`. -/
def cotB (f : F) : ℝ := cornerCot (M.lb f) (M.lc f) (M.la f)
/-- Cotangent of the corner of `f` at `c`. -/
def cotC (f : F) : ℝ := cornerCot (M.lc f) (M.la f) (M.lb f)

/-- Interior angle of `f` at `a`. -/
def angleA (f : F) : ℝ := cornerAngle (M.la f) (M.lb f) (M.lc f)
/-- Interior angle of `f` at `b`. -/
def angleB (f : F) : ℝ := cornerAngle (M.lb f) (M.lc f) (M.la f)
/-- Interior angle of `f` at `c`. -/
def angleC (f : F) : ℝ := cornerAngle (M.lc f) (M.la f) (M.lb f)

/-- Area of `f`. -/
def area (f : F) : ℝ := triArea (M.la f) (M.lb f) (M.lc f)

/-- Hypotheses on the mesh used below: every face has three distinct vertices and its three edge
lengths satisfy the strict triangle inequalities. (No manifoldness, no consistency of the lengths
across faces and no Delaunay condition are needed for the results of this file.) -/
structure Valid : Prop where
  ab : ∀ f, M.a f ≠ M.b f
  bc : ∀ f, M.b f ≠ M.c f
  ca : ∀ f, M.c f ≠ M.a f
  tri : ∀ f, IsTriangle (M.la f) (M.lb f) (M.lc f)

/-! ## Uniform scaling -/

/-- The mesh with all edge lengths multiplied by `s`. -/
def scale (s : ℝ) : TriMesh V F := { M with len := fun i j => s * M.len i j }


/-- The corner angles are invariant under uniform scaling. -/
theorem angleA_scale {s : ℝ} (hs : s ≠ 0) (f : F) : (M.scale s).angleA f = M.angleA f := by
  simp only [angleA, la, lb, lc, scale, cornerAngle_smul hs]

theorem angleB_scale {s : ℝ} (hs : s ≠ 0) (f : F) : (M.scale s).angleB f = M.angleB f := by
  simp only [angleB, la, lb, lc, scale, cornerAngle_smul hs]

theorem angleC_scale {s : ℝ} (hs : s ≠ 0) (f : F) : (M.scale s).angleC f = M.angleC f := by
  simp only [angleC, la, lb, lc, scale, cornerAngle_smul hs]

end TriMesh

end DDG
