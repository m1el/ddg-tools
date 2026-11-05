/-
Declaration-level cross-module dependency extractor (LIBRARY-PLAN.md, Appendix A).

For every declaration of every module whose name starts with one of the given prefixes, list the
constants of the *other* prefixed modules that its type or value uses. One TSV line per use:

    module <TAB> declaration <TAB> used-module <TAB> used-constant

Usage (run from a Lake package whose environment can import the modules; ~10 s on built oleans):

    cd lean/proofs
    lake env lean --run ../../scripts/lean-uses.lean [OUT.tsv [PREFIX,PREFIX,...]]

  OUT.tsv   default: ../../runs_tmp/uses.tsv
  PREFIXES  comma-separated module-name roots, default: DdgProofs,DdgOracle,DDG

The environment is built by `importModules` on the umbrella modules in the env var ROOTS
(comma-separated, default `DdgProofs`; after S4 use `ROOTS=DDG,DdgProofs`). The search path comes
from LEAN_PATH, which `lake env` sets.

Lazily generated auxiliary declarations (`.eq_N`, `.match_N`, `_proof_N`, `splitter`) are NOT
filtered here; use `scripts/lean-uses-aggregate.py`, which drops them by default.
-/
import Lean
open Lean

def splitCommas (s : String) : List String :=
  (s.splitOn ",").filter (· ≠ "")

def main (args : List String) : IO UInt32 := do
  let out := args.getD 0 "../../runs_tmp/uses.tsv"
  let prefixes := (splitCommas (args.getD 1 "DdgProofs,DdgOracle,DDG")).map String.toName
  let roots := (splitCommas ((← IO.getEnv "ROOTS").getD "DdgProofs")).map String.toName
  initSearchPath (← findSysroot)
  let env ← importModules (roots.map fun r => { module := r }).toArray {} 0
  let header := env.header
  let isOurs (m : Name) : Bool := prefixes.any (·.isPrefixOf m)
  let h ← IO.FS.Handle.mk out .write
  let mut n := 0
  for i in [0:header.moduleNames.size] do
    let modName := header.moduleNames[i]!
    if !isOurs modName then continue
    for c in header.moduleData[i]!.constNames do
      let some ci := env.find? c | continue
      for u in ci.getUsedConstantsAsSet do
        if let some j := env.getModuleIdxFor? u then
          let m2 := header.moduleNames[j.toNat]!
          if m2 != modName && isOurs m2 then
            h.putStrLn s!"{modName}\t{c}\t{m2}\t{u}"
            n := n + 1
  IO.eprintln s!"lean-uses: {n} uses written to {out}"
  return 0
