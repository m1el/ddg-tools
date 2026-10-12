import DdgOracle.Selftest

/-!
# Throughput benchmark

Times the complete text pipeline (parse, operators, solve, certificate, print) on meshes of
about 50 and about 200 vertices of every topology.
-/

namespace Ddg

/-- The benchmark meshes of roughly `size` vertices: sphere-like, torus, disk, annulus. -/
def benchMeshes (size : Nat) : Array (String × RawMesh) :=
  if size ≤ 100 then #[
    ("sphere", (uvSphere 7 8 fun p => ⟨1.3 * p.x, p.y, 0.8 * p.z⟩).jitter 0.02 1),          -- 50
    ("torus", (torus 10 5 1.0 0.4 fun a b => 0.2 * Float.sin (2 * a) * Float.cos b).jitter 0.01 2), -- 50
    ("disk", (gridDisk 7 7 fun x y => 0.4 * Float.sin (4 * x) * Float.sin (3 * y)).jitter 0.01 3),  -- 49
    ("annulus", (polar 4 10 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.2 * y * y).jitter 0.01 4)] -- 50
  else #[
    ("sphere", (uvSphere 12 18 fun p => ⟨1.3 * p.x, p.y, 0.8 * p.z⟩).jitter 0.01 5),         -- 200
    ("torus", (torus 20 10 1.0 0.4 fun a b => 0.2 * Float.sin (2 * a) * Float.cos b).jitter 0.005 6), -- 200
    ("disk", (gridDisk 14 14 fun x y => 0.4 * Float.sin (4 * x) * Float.sin (3 * y)).jitter 0.005 7), -- 196
    ("annulus", (polar 9 20 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.2 * y * y).jitter 0.005 8)] -- 200

/-- Run the benchmark and print one line per mesh and task. -/
def bench : IO Unit := do
  for size in #[50, 200] do
    for (kind, mesh) in benchMeshes size do
      let nv := mesh.pos.size
      let interior := (interiorVertices mesh).size
      -- three parameter sets: uniform, weighted with a free site, bounded
      let variants : Array Params := #[
        uniformParams nv 0.05,
        ((uniformParams nv 0.02).withRandomWeights 2.0 9).withFree (interiorVertices mesh |>.extract 0 1),
        (uniformParams nv 0.05).withBounds 0.0 posInf]
      let text := "".intercalate (variants.toList.map fun p => printCase { name := "bench", mesh, params := p })
      let reps := if size ≤ 100 then 4 else 1
      let t0 ← IO.monoMsNow
      for _ in [0:reps] do
        let out ← runText text (some .solve)
        if (out.splitOn "status converged").length != variants.size + 1 then
          IO.println s!"  (not all cases converged for {kind})"
      let t1 ← IO.monoMsNow
      for _ in [0:reps * 5] do
        let _ ← runText text (some .operators)
      let t2 ← IO.monoMsNow
      let solves := reps * variants.size
      let msSolve := Float.ofNat (t1 - t0) / Float.ofNat solves
      let msOps := Float.ofNat (t2 - t1) / Float.ofNat (solves * 5)
      IO.println s!"{kind} nv={nv} interior={interior}: solve {msSolve} ms/case ({1000.0 / msSolve} cases/s), operators {msOps} ms/case ({1000.0 / msOps} cases/s)"

end Ddg
