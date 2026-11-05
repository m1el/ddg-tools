import DDG.NormalCoords.FaceSplit

/-!
# Integer coordinates: splitting a crossed edge (`ic-split-crossed-edge`, tier B)

`crates/ddg-integer-coords/SPEC.md` §4.5 (I20); the Rust `formulas::split_side`.

A point `m` on side `0` (`a→b`) of a face `abc`, with `before` crossings of the side between `a`
and `m` (`0 ≤ before ≤ n⁺_ab`, `m` not on a crossing), is the gap `g = before` of the face's
boundary (between the points of rank `before` and `before + 1`). The new straight edge `mc`
crosses exactly the arcs that separate `m` from corner `c` (their ends lie on different sides of the
chord `mc`, `ic-chord-interleave`), which is `sepCount n c before 2` (`FaceSplit.lean`). Here:

* `splitSide_eq`: that number is the Rust formula
  `n_mc = c_c + e_a + e_b + max(0, c_a − before) + max(0, before − (n⁺_ab − c_b))`, with `c`, `e`
  of Eq. 1–2 (for side `s` of a face, the same applied to the face rotated by `s`; the face model
  is the canonical one of its counts whatever the numbering);
* `emanatingBefore_eq`: the arcs from `c` that reach `ab` before `m` are
  `clamp(before − c_a, 0, e_c)` (the Rust `emanating_before`).
-/

namespace DDG.IntegerCoords

open DDG.IntegerCoords.Formulas (cornerCoord emanating nPlus faceValid)

/-- **`ic-split-crossed-edge`: the new edge to the opposite corner of a split side.** For a valid
face (admissible corner counts `c`) and a split point after `before ≤ n_0` crossings of side `0`,
the number of arcs separating the split point from corner `2` is
`c₂ + e₀ + e₁ + max(0, c₀ − before) + max(0, before − (n₀ − c₁))`: the arcs around `c`, those from
`a` to `bc` and from `b` to `ca`, the arcs around `a` beyond `m` and those around `b` before `m`.
With `c`, `e` as Eq. 1–2 compute them from the side counts, this is the Rust `split_side`'s
`to_opposite`. [formula: ic-split-crossed-edge] -/
theorem splitSide_eq {n c : Fin 3 → ℕ} (hc : Admissible n c) {before : ℕ} (hb : before ≤ n 0) :
    (sepCount n c before 2 : ℤ) =
      cornerCoord (n 0) (n 2) (n 1) + emanating (n 1) (n 0) (n 2) + emanating (n 2) (n 1) (n 0) +
        max 0 (cornerCoord (n 1) (n 0) (n 2) - before) +
        max 0 (before - (n 0 - cornerCoord (n 2) (n 1) (n 0))) := by
  have hg : before < total n := by unfold total; omega
  have hs : gapSide n before = 0 := by unfold gapSide; simp [hb]
  rw [sepCount_eq hc hg]
  obtain ⟨t0, t1, t2, f0, f1, f2, -⟩ := values_side0 (c := c) hs
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [Fin.sum_univ_three, l1, l2, l3, l4, l5, l6, Fin.isValue, Fin.reduceEq, ite_true,
    ite_false, t0, t1, t2, f0, f1, f2]
  have k0 := hc.coords 0
  have k1 := hc.coords 1
  have k2 := hc.coords 2
  try simp only [l1, l2, l3, l4, l5, l6] at k0 k1 k2
  have e2 := hc.n_eq 0
  have a1 := hc.le 1
  try simp only [l1] at a1
  rw [k0.1, k1.1, k2.1, k0.2, k1.2]
  push_cast
  have := hc.le 0
  try simp only [l1] at this
  omega

/-- The arcs from corner `2` that reach side `0` before the split point: `clamp(before − c₀, 0,
e₂)` (the Rust `emanating_before`, for the roundabout of the new halfedge `c→m`).
[formula: ic-split-crossed-edge] -/
theorem emanatingBefore_eq {n c : Fin 3 → ℕ} (before : ℕ) :
    ((Finset.range (em n c 2)).filter (fun s => c 0 + s < before)).card =
      min (em n c 2) (before - c 0) := by
  have : (Finset.range (em n c 2)).filter (fun s => c 0 + s < before) =
      Finset.range (min (em n c 2) (before - c 0)) := by
    ext s; simp only [Finset.mem_filter, Finset.mem_range]; omega
  rw [this, Finset.card_range]

end DDG.IntegerCoords
