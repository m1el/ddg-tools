import DdgProofs.IntegerCoords.Formulas
import DdgProofs.IntegerCoords.FaceModel
import DdgProofs.IntegerCoords.FaceSplit
import DdgProofs.IntegerCoords.Flip
import DDG.Intrinsic.VertexRemoval
import DdgProofs.IntegerCoords.Chords
import DdgProofs.IntegerCoords.Recover
import DdgProofs.IntegerCoords.Transfer
import DdgProofs.IntegerCoords.SplitEdge
import DdgProofs.IntegerCoords.SplitAtCrossing
import DDG.Intrinsic.Strip
import DdgProofs.IntegerCoords.Complex.Basic
import DdgProofs.IntegerCoords.Complex.Trace
import DdgProofs.IntegerCoords.Complex.Labels
import DdgProofs.IntegerCoords.Complex.QuadFlip
import DdgProofs.IntegerCoords.Complex.Flip
import DdgProofs.IntegerCoords.Complex.Merge
import DdgProofs.IntegerCoords.Complex.SplitAtCrossing
import DdgProofs.IntegerCoords.Complex.Extract
import DdgProofs.IntegerCoords.Frames.Angle
import DdgProofs.IntegerCoords.Frames.Input
import DdgProofs.IntegerCoords.Frames.Inserted
import DdgProofs.IntegerCoords.Frames.Wedge

/-!
# Integer Coordinates for Intrinsic Geometry Processing (Gillespie, Sharp, Crane 2021)

Tiers A, B, B′ and C (`Complex/`, the complex model; `Frames/`, the input-aligned frames) of `crates/ddg-integer-coords/SPEC.md` §9.2. The face model and how tier C extends
it: `DdgProofs/IntegerCoords/DESIGN.md`. Cited results: `CITED.md`, section "Integer Coordinates".
-/
