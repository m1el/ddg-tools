import DdgProofs

/-!
# Axiom audit of every declaration

`lake env lean DdgProofs/AuditAll.lean` collects the axioms of **every** declaration defined in a
`DdgProofs.*` or `DdgOracle.*` module (theorems, definitions, instances, auxiliary declarations)
and reports the union and every declaration that uses an axiom other than `propext`,
`Classical.choice`, `Quot.sound`. It takes several minutes; `Audit.lean` (the main theorems, by
name) is the quick check that `scripts/check-formulas.sh --build` runs.

This file is not imported by the library.
-/

open Lean Elab Command

/-- Collect the axioms of every declaration defined in a `DdgProofs.*` or `DdgOracle.*` module. -/
elab "#audit_all" : command => do
  let env ← getEnv
  let mut cnt : Nat := 0
  let mut bad : Array (Name × Array Name) := #[]
  let mut all : Std.HashSet Name := {}
  for (c, _) in env.constants.toList do
    match env.getModuleIdxFor? c with
    | none => pure ()
    | some idx =>
      let m := env.header.moduleNames[idx.toNat]!
      if (`DdgProofs).isPrefixOf m || (`DdgOracle).isPrefixOf m then
        cnt := cnt + 1
        let axs ← Lean.collectAxioms c
        for a in axs do all := all.insert a
        let extra := axs.filter fun a =>
          a != ``propext && a != ``Classical.choice && a != ``Quot.sound
        if extra.size > 0 then bad := bad.push (c, extra)
  logInfo m!"declarations checked: {cnt}; union of axioms: {all.toList}; offenders: {bad.size}"
  for (c, e) in bad.toList.take 50 do logInfo m!"{c}: {e}"

#audit_all
