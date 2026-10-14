import DdgOracle.Intrinsic.Format
import DdgOracle.Cuts.Gen

/-!
# Cases of the intrinsic operations built from generated meshes
-/

namespace Ddg.Intrinsic

open Ddg Ddg.Cuts

/-- A case of the given task on a generated mesh; the edge lengths are computed from the positions. -/
def mkICase (name task : String) (m : RawMesh) : ICase :=
  let c := mkCase name m (zeros m.pos.size)
  { name, task, nv := c.nv, faces := c.faces, edges := c.edges }

/-- The input triangulation of a case (panics-free: `none` on an invalid case). -/
def triOf (c : ICase) : Option Tri :=
  match buildInput c with
  | .ok (t, _) => some t
  | .error _ => none

/-- The first edge `(a, b)` of the mesh (halfedge order) whose flip is refused for the given reason, or
accepted (`none` reason). -/
def firstEdgeWith (c : ICase) (reason : Option Refusal) : Option (Nat × Nat) := Id.run do
  match triOf c with
  | none => return none
  | some t =>
    for h in [0:t.nh] do
      if t.eid[h]! == h then
        let r := flipHalfedge t h
        let hit := match reason, r with
          | none, .ok _ => true
          | some x, .error y => x == y
          | _, _ => false
        if hit then return some (t.tail[h]!, t.head h)
    return none

end Ddg.Intrinsic
