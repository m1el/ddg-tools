import DDG.NormalCoords.Chords
import DDG.NormalCoords.Complex.Basic
import DDG.NormalCoords.Complex.Extract
import DDG.NormalCoords.Complex.Flip
import DDG.NormalCoords.Complex.Labels
import DDG.NormalCoords.Complex.Merge
import DDG.NormalCoords.Complex.Quad
import DDG.NormalCoords.Complex.QuadAdm1
import DDG.NormalCoords.Complex.QuadAdm2
import DDG.NormalCoords.Complex.QuadCount
import DDG.NormalCoords.Complex.QuadFlip
import DDG.NormalCoords.Complex.QuadOuterIL
import DDG.NormalCoords.Complex.QuadOuterJK
import DDG.NormalCoords.Complex.QuadOuterKI
import DDG.NormalCoords.Complex.QuadOuterLJ
import DDG.NormalCoords.Complex.QuadSide
import DDG.NormalCoords.Complex.SplitAtCrossing
import DDG.NormalCoords.Complex.Trace
import DDG.NormalCoords.FaceModel
import DDG.NormalCoords.FaceSplit
import DDG.NormalCoords.Flip
import DDG.NormalCoords.Formulas
import DDG.NormalCoords.Frames.Angle
import DDG.NormalCoords.Frames.Input
import DDG.NormalCoords.Frames.Inserted
import DDG.NormalCoords.Frames.Wedge
import DDG.NormalCoords.Recover
import DDG.NormalCoords.SplitAtCrossing
import DDG.NormalCoords.SplitEdge

/-!
# DDG.NormalCoords: normal (integer) coordinates for intrinsic triangulations

**Scope.** The combinatorial and arithmetic core of integer-coordinate intrinsic triangulations:
a triangulation is stored as, per edge, the number of crossings with the input mesh, with no
floating-point geometry. Statements are over `ℤ` and `ℕ`. The face model and how the complex model
extends it are in `DDG/NormalCoords/DESIGN.md`.

**Main results (tiers A, B, B′, C).**
* `Formulas`: the integer formulas, defined natively (corner coordinates, emanating counts, face
  validity, Eq. 1–4, 6, trace step).
* `FaceModel`, `Chords`, `FaceSplit`, `SplitEdge`, `SplitAtCrossing`, `Recover`: the face model
  (`FaceModel.Valid`, structure theorem `Valid.eq_canon`), chords of a convex polygon cross iff
  their ends interleave, splitting a face or a crossed edge, positions of an inserted vertex.
* `Flip`: the quadrilateral model and Eq. 3 (`flip_count`, `flipNormal_flipNormal`).
* `Complex/*`: the complex model `Cx` (faces glued along sides), the flip of a complex
  (`Complex.Flip.flip_valid`), tracing terminates (`Valid.trace_terminates`), labels and
  roundabouts (Alg. 3, 4), merging the two halves of a split edge, extraction of an input edge.
* `Frames/*`: input-aligned frames at input and inserted vertices.

**Sources.** M. Gillespie, N. Sharp, K. Crane, *Integer Coordinates for Intrinsic Geometry
Processing*, ACM TOG 40(6), 2021, Sec. 2.3–3.7, Eq. 1–6, App. A, D and its supplement, written
from the paper's text; N. Sharp, Y. Soliman, K. Crane, *Navigating Intrinsic Triangulations*,
ACM TOG 2019. Where a reference implementation was consulted for choices the paper leaves open,
it was the authors' MIT-licensed geometry-central (see `crates/ddg-integer-coords/README.md`).

**Depends on:** `DDG.Intrinsic` (and, through it, `DDG.Geometry`).
-/
