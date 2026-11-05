import DdgOracle.IntegerCoords.Formulas
import DDG.NormalCoords.Formulas

/-!
# Oracle bridges for the integer-coordinate formulas

The library (`DDG/NormalCoords/Formulas.lean`, namespace `DDG.NormalCoords.Formulas`) defines
the pure integer formulas natively. This app module relates them to the oracle's `Ddg.IntegerCoords`
definitions: `rfl` for every definition whose type mentions no oracle inductive, and a conversion map
plus a commutation theorem for those that do (`Step`, `Pt`, `Region`, `Decoded`).
-/

namespace DDG
namespace Oracle
namespace NormalCoords

open DDG.NormalCoords.Formulas (nPlus)

theorem nPlus_eq : Ddg.IntegerCoords.nPlus = DDG.NormalCoords.Formulas.nPlus := rfl
theorem nMinus_eq : Ddg.IntegerCoords.nMinus = DDG.NormalCoords.Formulas.nMinus := rfl
theorem emanating_eq : Ddg.IntegerCoords.emanating = DDG.NormalCoords.Formulas.emanating := rfl
theorem cornerNumerator_eq :
    Ddg.IntegerCoords.cornerNumerator = DDG.NormalCoords.Formulas.cornerNumerator := rfl
theorem cornerCoord_eq : Ddg.IntegerCoords.cornerCoord = DDG.NormalCoords.Formulas.cornerCoord := rfl
theorem cornerCoord?_eq :
    Ddg.IntegerCoords.cornerCoord? = DDG.NormalCoords.Formulas.cornerCoord? := rfl
theorem faceValid_eq : Ddg.IntegerCoords.faceValid = DDG.NormalCoords.Formulas.faceValid := rfl
theorem next_eq : Ddg.IntegerCoords.next = DDG.NormalCoords.Formulas.next := rfl
theorem prev_eq : Ddg.IntegerCoords.prev = DDG.NormalCoords.Formulas.prev := rfl
theorem side_eq : Ddg.IntegerCoords.side = DDG.NormalCoords.Formulas.side := rfl
theorem faceCorner_eq : Ddg.IntegerCoords.faceCorner = DDG.NormalCoords.Formulas.faceCorner := rfl
theorem faceEmanating_eq :
    Ddg.IntegerCoords.faceEmanating = DDG.NormalCoords.Formulas.faceEmanating := rfl
theorem reverseCrossing_eq :
    Ddg.IntegerCoords.reverseCrossing = DDG.NormalCoords.Formulas.reverseCrossing := rfl
theorem iabs_eq : Ddg.IntegerCoords.iabs = DDG.NormalCoords.Formulas.iabs := rfl
theorem flipDoubled_eq : Ddg.IntegerCoords.flipDoubled = DDG.NormalCoords.Formulas.flipDoubled := rfl
theorem flipNormal_eq : Ddg.IntegerCoords.flipNormal = DDG.NormalCoords.Formulas.flipNormal := rfl
theorem roundaboutUpdate_eq :
    Ddg.IntegerCoords.roundaboutUpdate = DDG.NormalCoords.Formulas.roundaboutUpdate := rfl
theorem splitEdgeNormal_eq :
    Ddg.IntegerCoords.splitEdgeNormal = DDG.NormalCoords.Formulas.splitEdgeNormal := rfl
theorem identifyIndex_eq :
    Ddg.IntegerCoords.identifyIndex = DDG.NormalCoords.Formulas.identifyIndex := rfl

/-! ## Types: conversion maps (the oracle's and the native inductives are isomorphic copies) -/

open DDG.NormalCoords.Formulas in
/-- Native `Step` to the oracle's. -/
def stepToOracle : Step → Ddg.IntegerCoords.Step
  | .toIK p => .toIK p
  | .toKJ p => .toKJ p
  | .vertex => .vertex

open DDG.NormalCoords.Formulas in
/-- Native `Pt` to the oracle's. -/
def ptToOracle : Pt → Ddg.IntegerCoords.Pt
  | .corner k => .corner k
  | .side k p => .side k p

open DDG.NormalCoords.Formulas in
/-- Native `Region` to the oracle's. -/
def regionToOracle : Region → Ddg.IntegerCoords.Region
  | .corner k t => .corner k t
  | .middle s => .middle s

open DDG.NormalCoords.Formulas in
/-- Native `Decoded` to the oracle's. -/
def decodedToOracle (d : Decoded) : Ddg.IntegerCoords.Decoded :=
  { n := d.n, c := d.c, kappa := d.kappa, e := d.e }

attribute [local simp] nPlus_eq nMinus_eq emanating_eq cornerNumerator_eq cornerCoord_eq faceValid_eq
  next_eq prev_eq side_eq faceCorner_eq faceEmanating_eq

theorem traceStep_eq (nij nik nkj p : ℤ) :
    Ddg.IntegerCoords.traceStep nij nik nkj p =
      stepToOracle (DDG.NormalCoords.Formulas.traceStep nij nik nkj p) := by
  simp only [Ddg.IntegerCoords.traceStep, DDG.NormalCoords.Formulas.traceStep, nPlus_eq,
    cornerCoord_eq, apply_ite stepToOracle]
  split_ifs <;> simp_all [stepToOracle]

theorem partner_eq (n c : ℕ → ℕ) (k p : ℕ) :
    Ddg.IntegerCoords.partner n c k p = ptToOracle (DDG.NormalCoords.Formulas.partner n c k p) := by
  simp only [Ddg.IntegerCoords.partner, DDG.NormalCoords.Formulas.partner, next_eq, prev_eq]
  by_cases h1 : p < c k
  · simp [h1, ptToOracle]
  · by_cases h2 : p + c (DDG.NormalCoords.Formulas.next k) < n k
    · simp [h1, h2, ptToOracle]
    · simp [h1, h2, ptToOracle]

theorem decode?_eq (n0 n1 n2 : ℤ) :
    Ddg.IntegerCoords.decode? n0 n1 n2 =
      (DDG.NormalCoords.Formulas.decode? n0 n1 n2).map decodedToOracle := by
  simp only [Ddg.IntegerCoords.decode?, DDG.NormalCoords.Formulas.decode?, faceValid_eq, nPlus_eq,
    side_eq, faceCorner_eq, faceEmanating_eq]
  by_cases h : DDG.NormalCoords.Formulas.faceValid n0 n1 n2 = true
  · simp only [h, Bool.not_true, Bool.false_eq_true, ite_false]
    generalize List.find? (fun ke : ℕ × ℕ => decide (ke.2 > 0))
      (List.map (fun k : ℕ => (k, (DDG.NormalCoords.Formulas.faceEmanating n0 n1 n2 k).toNat))
        [0, 1, 2]) = r
    rcases r with _ | ⟨k, e⟩ <;> rfl
  · simp [h]

theorem regions_eq (d : DDG.NormalCoords.Formulas.Decoded) :
    Ddg.IntegerCoords.regions (decodedToOracle d) =
      (DDG.NormalCoords.Formulas.regions d).map regionToOracle := by
  simp only [Ddg.IntegerCoords.regions, DDG.NormalCoords.Formulas.regions, decodedToOracle,
    List.map_append, List.flatMap_cons, List.flatMap_nil, List.map_map, List.append_nil]
  rfl

theorem region_valid_eq (d : DDG.NormalCoords.Formulas.Decoded)
    (R : DDG.NormalCoords.Formulas.Region) :
    Ddg.IntegerCoords.Region.valid (decodedToOracle d) (regionToOracle R) =
      DDG.NormalCoords.Formulas.Region.valid d R := by
  cases R <;> rfl

theorem splitCount_eq (d : DDG.NormalCoords.Formulas.Decoded)
    (R : DDG.NormalCoords.Formulas.Region) (j : ℕ) :
    Ddg.IntegerCoords.splitCount (decodedToOracle d) (regionToOracle R) j =
      DDG.NormalCoords.Formulas.splitCount d R j := by
  cases R <;> rfl

theorem faceSplit?_eq (n0 n1 n2 : ℤ) (R : DDG.NormalCoords.Formulas.Region) :
    Ddg.IntegerCoords.faceSplit? n0 n1 n2 (regionToOracle R) =
      DDG.NormalCoords.Formulas.faceSplit? n0 n1 n2 R := by
  unfold Ddg.IntegerCoords.faceSplit? DDG.NormalCoords.Formulas.faceSplit?
  rw [decode?_eq]
  cases DDG.NormalCoords.Formulas.decode? n0 n1 n2 with
  | none => rfl
  | some d =>
    simp only [Option.map_some, region_valid_eq, splitCount_eq]

end NormalCoords
end Oracle
end DDG
