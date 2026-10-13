import DdgProofs

/-!
# Axiom audit of every declaration

`lake env lean DdgProofs/AuditAll.lean` collects the axioms of **every** declaration defined in a
`DdgProofs.*` or `DdgOracle.*` module (theorems, definitions, instances, auxiliary and private
declarations) and reports the union and every declaration that uses an axiom other than `propext`,
`Classical.choice`, `Quot.sound`. `scripts/check-formulas.sh --audit-all` runs it.

The declarations are taken from the module data of this project's modules
(`env.header.moduleData[i].constNames`), not by filtering every constant of the environment: with
Mathlib imported that is about 785 000 constants, and classifying each of them by module took
over six minutes in the interpreter, while the audit of the declarations found (1989 then, 2538
with the cut-flow modules) takes milliseconds (this Lean version records each declaration's
axioms when its module is compiled).
The whole file now runs in a few seconds, nearly all of it loading Mathlib.

This file is not imported by the library.
-/

open Lean Elab Command

/-- Collect the axioms of every declaration defined in a `DdgProofs.*` or `DdgOracle.*` module. -/
elab "#audit_all" : command => do
  let env ← getEnv
  let mut cnt : Nat := 0
  let mut bad : Array (Name × Array Name) := #[]
  let mut all : Std.HashSet Name := {}
  for i in [0:env.header.moduleNames.size] do
    let m := env.header.moduleNames[i]!
    if (`DdgProofs).isPrefixOf m || (`DdgOracle).isPrefixOf m then
      for c in env.header.moduleData[i]!.constNames do
        cnt := cnt + 1
        let axs ← Lean.collectAxioms c
        for a in axs do all := all.insert a
        let extra := axs.filter fun a =>
          a != ``propext && a != ``Classical.choice && a != ``Quot.sound
        if extra.size > 0 then bad := bad.push (c, extra)
  logInfo m!"declarations checked: {cnt}; union of axioms: {all.toList}; offenders: {bad.size}"
  for (c, e) in bad.toList.take 50 do logInfo m!"{c}: {e}"

#audit_all
