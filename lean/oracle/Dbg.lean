import DdgOracle
open Ddg

def main (_args : List String) : IO UInt32 := do
  let c : Case := { name := "", mesh := cube, params := (uniformParams 8 0.3).withBounds 0.0 (floatPi / 2) }
  if let some (_, r) := solveCase c then
    IO.println s!"{repr r.status} {repr r.stop} {r.itersOne} {r.itersTwo}"
    for v in [0:8] do
      IO.println s!"{v}: u={r.u.get! v} phi={r.phi.get! v} mu={floatToHex (r.mu.get! v)} {floatToHex (floatPi/2)}"
  IO.println s!"{sci (cornerAngle (1e-7 : Float) 1.0 1.0 - 2 * Float.asin 0.5e-7)}"
  let disk := (polar 4 10 0.0 1.0 fun x y => 0.6 * (1 - x * x - y * y)).jitter 0.01 41
  for l in #[0.002, 0.005, 0.01, 0.02, 0.05, 0.1] do
    if let some (ops, r) := solveCase { name := "", mesh := disk, params := uniformParams disk.pos.size l } then
      let cones := (List.range ops.nv).filter (fun v => r.mu.get! v != 0.0) |>.length
      IO.println s!"lambda={l} cones={cones}"
  return 0
