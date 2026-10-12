import DdgOracle
open Ddg

def usage : String :=
"ddg-oracle: dense Float reference implementation of
  'Optimal Cone Singularities for Conformal Flattening' (Soliman, Slepcev, Crane 2018)

usage:
  ddg-oracle operators [IN [OUT]]   curvatures, masses and cotan Laplacian of every case
  ddg-oracle solve     [IN [OUT]]   solve the canonical problem for every case
  ddg-oracle batch     [IN [OUT]]   like the above; each case names its task in a 'task' line
  ddg-oracle selftest  [SWEEP]      internal consistency checks (SWEEP random cases, default 300)
  ddg-oracle version

IN and OUT default to standard input and output ('-'). All three processing commands read any
number of cases from one stream and flush after every case, so one process can serve a fuzzer.
The text format is specified in README.md."

def openInput (path : String) : IO IO.FS.Stream := do
  if path == "-" then IO.getStdin
  else return IO.FS.Stream.ofHandle (← IO.FS.Handle.mk path .read)

def openOutput (path : String) : IO IO.FS.Stream := do
  if path == "-" then IO.getStdout
  else return IO.FS.Stream.ofHandle (← IO.FS.Handle.mk path .write)

def process (fixed : Option Task) (rest : List String) : IO UInt32 := do
  if rest.length > 2 then
    IO.eprintln usage
    return 2
  let input ← openInput (rest.getD 0 "-")
  let output ← openOutput (rest.getD 1 "-")
  let _ ← runStream input output fixed
  output.flush
  return 0

def main (args : List String) : IO UInt32 := do
  match args with
  | "operators" :: rest => process (some .operators) rest
  | "solve" :: rest => process (some .solve) rest
  | "batch" :: rest => process none rest
  | ["selftest"] => return if (← selftest) == 0 then 0 else 1
  | ["selftest", n] =>
    match n.toNat? with
    | some k => return if (← selftest k) == 0 then 0 else 1
    | none => IO.eprintln usage; return 2
  | ["version"] => IO.println "ddg-oracle 0.1.0 (text format 1)"; return 0
  | _ => IO.eprintln usage; return 2
