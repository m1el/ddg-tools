import DdgOracle.Generators
import DdgOracle.PointsAsTori.Format

/-!
# Seeded generators for the `pat` selftest and fixtures
-/

namespace Ddg.PointsAsTori

abbrev G := StateM Rng

def uni (a b : Float) : G Float := modifyGet fun r => r.range a b
def below (n : Nat) : G Nat := modifyGet fun r => r.below n

def twoPi : Float := 6.283185307179586

def mkRng (seed : Nat) : Rng := ⟨seed.toUInt64 * 0x9E3779B97F4A7C15 + 17⟩

/-- Uniform direction. -/
def unitVec : G V3 := do
  let z ← uni (-1.0) 1.0
  let a ← uni 0.0 twoPi
  let r := Float.sqrt (1.0 - z * z)
  return ⟨r * Float.cos a, r * Float.sin a, z⟩

def vecIn (lo hi : Float) : G V3 := do
  return ⟨← uni lo hi, ← uni lo hi, ← uni lo hi⟩

/-- `n` Fibonacci points on a sphere of radius `rad` about `c`, with outward unit normals. -/
def fibonacci (n : Nat) (rad : Float) (c : V3) : Array V3 × Array V3 := Id.run do
  let ga := 2.399963229728653
  let mut ps := #[]
  let mut ns := #[]
  for i in [0:n] do
    let z := 1.0 - (2.0 * i.toFloat + 1.0) / n.toFloat
    let r := Float.sqrt (1.0 - z * z)
    let a := ga * i.toFloat
    let u : V3 := ⟨r * Float.cos a, r * Float.sin a, z⟩
    ps := ps.push (c + u.smul rad)
    ns := ns.push u
  return (ps, ns)

/-- A cloud case skeleton. -/
def cloudCase (name task : String) (ps ns : Array V3) (areas : Option (Array Float) := none) : PatCase :=
  { name, task, pts := ps, nrm := ns, areas := areas.getD (ps.map fun _ => 1.0) }

/-- Random coefficients of moderate size (curvatures of order 0.3 to 1.5, tilts up to 0.3). -/
def randCoef : G (Array Float) := do
  let a00 ← uni (-0.03) 0.03
  let a01 ← uni (-0.3) 0.3
  let a10 ← uni (-0.3) 0.3
  let a11 ← uni (-0.6) 0.6
  let sg1 ← below 2
  let sg2 ← below 2
  let m1 ← uni 0.15 0.8
  let m2 ← uni 0.15 0.8
  let a02 := (if sg1 == 0 then 1.0 else -1.0) * m1
  let a20 := (if sg2 == 0 then 1.0 else -1.0) * m2
  return #[a00, a01, a10, a11, a02, a20]

/-- Random queries in a box. -/
def randQueries (n : Nat) (lo hi : Float) : G (Array V3) := do
  let mut qs := #[]
  for _ in [0:n] do qs := qs.push (← vecIn lo hi)
  return qs

end Ddg.PointsAsTori
