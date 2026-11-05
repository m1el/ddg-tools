import DDG.NormalCoords.Formulas
import DDG.NormalCoords.FaceModel
import DDG.NormalCoords.FaceSplit
import DDG.NormalCoords.Flip
import DDG.Intrinsic.VertexRemoval
import DDG.NormalCoords.Chords
import DDG.NormalCoords.Recover
import DDG.FEM.Transfer
import DDG.NormalCoords.SplitEdge
import DDG.NormalCoords.SplitAtCrossing
import DDG.Intrinsic.Strip
import DDG.NormalCoords.Complex.Basic
import DDG.NormalCoords.Complex.Trace
import DDG.NormalCoords.Complex.Labels
import DDG.NormalCoords.Complex.QuadFlip
import DDG.NormalCoords.Complex.Flip
import DDG.NormalCoords.Complex.Merge
import DDG.NormalCoords.Complex.SplitAtCrossing
import DDG.NormalCoords.Complex.Extract
import DDG.NormalCoords.Frames.Angle
import DDG.NormalCoords.Frames.Input
import DDG.NormalCoords.Frames.Inserted
import DDG.NormalCoords.Frames.Wedge

/-!
# Integer Coordinates for Intrinsic Geometry Processing (Gillespie, Sharp, Crane 2021)

Tiers A, B, B′ and C (`Complex/`, the complex model; `Frames/`, the input-aligned frames) of `crates/ddg-integer-coords/SPEC.md` §9.2. The face model and how tier C extends
it: `DdgProofs/IntegerCoords/DESIGN.md`. Cited results: `CITED.md`, section "Integer Coordinates".
-/
