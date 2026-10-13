import Std.Data.HashMap
import DdgOracle.Generators
import DdgOracle.Cuts.Format

/-!
# Test meshes and cases for the cut flow
-/

namespace Ddg.Cuts

/-- Icosphere of the given subdivision level on the unit sphere (12, 42, 162, 642, … vertices),
faces oriented outward. -/
def icosphere (level : Nat) : RawMesh := Id.run do
  let t := (1 + Float.sqrt 5) / 2
  let raw : Array (Point Float) :=
    #[⟨-1, t, 0⟩, ⟨1, t, 0⟩, ⟨-1, -t, 0⟩, ⟨1, -t, 0⟩, ⟨0, -1, t⟩, ⟨0, 1, t⟩, ⟨0, -1, -t⟩,
      ⟨0, 1, -t⟩, ⟨t, 0, -1⟩, ⟨t, 0, 1⟩, ⟨-t, 0, -1⟩, ⟨-t, 0, 1⟩]
  let norm := fun (p : Point Float) =>
    let n := Float.sqrt (p.x * p.x + p.y * p.y + p.z * p.z)
    (⟨p.x / n, p.y / n, p.z / n⟩ : Point Float)
  let mut pos := raw.map norm
  let mut faces : Array Face :=
    #[⟨0, 11, 5⟩, ⟨0, 5, 1⟩, ⟨0, 1, 7⟩, ⟨0, 7, 10⟩, ⟨0, 10, 11⟩, ⟨1, 5, 9⟩, ⟨5, 11, 4⟩,
      ⟨11, 10, 2⟩, ⟨10, 7, 6⟩, ⟨7, 1, 8⟩, ⟨3, 9, 4⟩, ⟨3, 4, 2⟩, ⟨3, 2, 6⟩, ⟨3, 6, 8⟩,
      ⟨3, 8, 9⟩, ⟨4, 9, 5⟩, ⟨2, 4, 11⟩, ⟨6, 2, 10⟩, ⟨8, 6, 7⟩, ⟨9, 8, 1⟩]
  for _ in [0:level] do
    let mut cache : Std.HashMap (Nat × Nat) Nat := {}
    let mut next : Array Face := #[]
    for f in faces do
      let mut mids : Array Nat := #[]
      for (x, y) in [(f.a, f.b), (f.b, f.c), (f.c, f.a)] do
        let key := (min x y, max x y)
        match cache.get? key with
        | some m => mids := mids.push m
        | none =>
          let p := pos[x]!
          let q := pos[y]!
          let m := pos.size
          pos := pos.push (norm ⟨(p.x + q.x) / 2, (p.y + q.y) / 2, (p.z + q.z) / 2⟩)
          cache := cache.insert key m
          mids := mids.push m
      let ab := mids[0]!
      let bc := mids[1]!
      let ca := mids[2]!
      next := ((((next.push ⟨f.a, ab, ca⟩).push ⟨f.b, bc, ab⟩).push ⟨f.c, ca, bc⟩).push ⟨ab, bc, ca⟩)
    faces := next
  return { pos, faces }

/-- Two meshes side by side (disjoint union); the second is translated by `shift`. -/
def RawMesh.union (m1 m2 : RawMesh) (shift : Point Float) : RawMesh :=
  let off := m1.pos.size
  { pos := m1.pos ++ m2.pos.map fun p => ⟨p.x + shift.x, p.y + shift.y, p.z + shift.z⟩
    faces := m1.faces ++ m2.faces.map fun f => ⟨f.a + off, f.b + off, f.c + off⟩ }

/-- A case with default parameters (Dirichlet, Gaussian curvature, `α_L = 0.05`, fixed `τ = 0.02`,
`w = (1, 1)`, `β = (0, 0)`, `δ = 10⁻³`, `ρ₀ = 5`, one smoothing pass); the edge lengths are
computed from the positions. `phi` is in the caller's unit. -/
def mkCase (name : String) (m : RawMesh) (phi : FloatArray) : CutCase := Id.run do
  let nv := m.pos.size
  let mut edges : Array (Nat × Nat × Float) := #[]
  match buildTopology nv m.faces (connected := false) with
  | .error _ => pure ()
  | .ok topo =>
    for a in [0:nv] do
      for b in topo.neighbors[a]! do
        if a < b then edges := edges.push (a, b, dist m.pos[a]! m.pos[b]!)
  return { name, nv, faces := m.faces, edges, phi
           energy := .dirichlet, curv := .gaussian, alphaL := 0.05, eta := none
           penalty := .gradient, w := (1.0, 1.0), beta := (0.0, 0.0), step := .fixed 0.02
           delta := 1e-3, rho0 := 5.0, passes := 1, classes := none }

/-- `φ_i = f(position i)`. -/
def phiOf (m : RawMesh) (f : Point Float → Float) : FloatArray :=
  FloatArray.ofFn' m.pos.size fun i => f m.pos[i]!

end Ddg.Cuts
