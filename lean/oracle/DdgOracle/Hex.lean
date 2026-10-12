/-!
# Bit-exact text encoding of `Float`

A `Float` crosses the process boundary as the 16 hexadecimal digits of its IEEE-754 binary64
bit pattern, most significant digit first (`3ff0000000000000` is `1.0`,
`7ff0000000000000` is `+∞`). Decimal is never used.
-/

namespace Ddg

/-- Lower-case hexadecimal digit of `n < 16`. -/
def hexDigit (n : Nat) : Char :=
  if n < 10 then Char.ofNat (48 + n) else Char.ofNat (87 + n)

/-- Value of a hexadecimal digit (either case). -/
def hexValue? (c : Char) : Option Nat :=
  let n := c.toNat
  if 48 ≤ n ∧ n ≤ 57 then some (n - 48)
  else if 97 ≤ n ∧ n ≤ 102 then some (n - 87)
  else if 65 ≤ n ∧ n ≤ 70 then some (n - 55)
  else none

/-- The 16 hex digits of the bit pattern of `x`. -/
def floatToHex (x : Float) : String := Id.run do
  let bits := x.toBits.toNat
  let mut s := ""
  for k in [0:16] do
    s := s.push (hexDigit ((bits >>> (4 * (15 - k))) % 16))
  return s

/-- Inverse of `floatToHex`; `none` unless the token is exactly 16 hex digits. -/
def hexToFloat? (s : String) : Option Float := Id.run do
  let cs := s.toList
  if cs.length != 16 then return none
  let mut bits : Nat := 0
  for c in cs do
    match hexValue? c with
    | some d => bits := bits * 16 + d
    | none => return none
  return some (Float.ofBits bits.toUInt64)

/-- Whitespace-separated tokens of a line (space, tab, CR, LF). -/
def tokens (s : String) : Array String := Id.run do
  let mut out : Array String := #[]
  let mut cur := ""
  for c in s.toList do
    if c == ' ' || c == '\t' || c == '\r' || c == '\n' then
      if cur != "" then
        out := out.push cur
        cur := ""
    else
      cur := cur.push c
  if cur != "" then out := out.push cur
  return out

/-- `+∞`. -/
def posInf : Float := Float.ofBits 0x7ff0000000000000

/-- `−∞`. -/
def negInf : Float := Float.ofBits 0xfff0000000000000

end Ddg
